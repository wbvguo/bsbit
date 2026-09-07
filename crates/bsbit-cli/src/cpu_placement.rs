//! Best-effort internal CPU placement for alignment pipeline roles.
//!
//! Placement is deliberately not part of the command-line contract.  The
//! inherited CPU set remains the authority (including `taskset`, containers,
//! and job schedulers); this module only partitions that set so mapping and
//! streaming workers do not compete for the same physical cores.

use std::collections::BTreeSet;

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct PhysicalCore {
    package: i64,
    core: i64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct AvailableCpu {
    logical: usize,
    physical: PhysicalCore,
}

/// One automatically derived, process-local role partition.
#[derive(Clone, Debug, Default)]
pub(crate) struct CpuPlacement {
    mapping: Vec<usize>,
    auxiliary: Vec<usize>,
}

impl CpuPlacement {
    /// Detects CPUs allowed by the parent process and assigns mapping workers
    /// to distinct physical cores whenever topology information permits it.
    #[must_use]
    pub(crate) fn detect(mapping_workers: usize) -> Self {
        let available = platform::available_cpus();
        Self::from_available(&available, mapping_workers)
    }

    fn from_available(available: &[AvailableCpu], mapping_workers: usize) -> Self {
        if available.is_empty() || mapping_workers == 0 {
            return Self::default();
        }

        let assigned_cpu_count = mapping_workers.min(available.len());
        let mut represented_cores = BTreeSet::new();
        let mut mapping = Vec::new();
        if mapping.try_reserve_exact(assigned_cpu_count).is_err() {
            return Self::default();
        }
        for cpu in available {
            if represented_cores.insert(cpu.physical) {
                mapping.push(cpu.logical);
                if mapping.len() == assigned_cpu_count {
                    break;
                }
            }
        }

        // When workers outnumber physical cores, consume otherwise unused SMT
        // lanes before sharing a logical CPU.
        if mapping.len() < assigned_cpu_count {
            for cpu in available {
                if !mapping.contains(&cpu.logical) {
                    mapping.push(cpu.logical);
                    if mapping.len() == assigned_cpu_count {
                        break;
                    }
                }
            }
        }

        let mapping_cores = available
            .iter()
            .filter(|cpu| mapping.contains(&cpu.logical))
            .map(|cpu| cpu.physical)
            .collect::<BTreeSet<_>>();
        let mut auxiliary = available
            .iter()
            .filter(|cpu| !mapping_cores.contains(&cpu.physical))
            .map(|cpu| cpu.logical)
            .collect::<Vec<_>>();
        if auxiliary.is_empty() {
            auxiliary.extend(available.iter().map(|cpu| cpu.logical));
        }

        Self { mapping, auxiliary }
    }

    /// Pins one mapping worker to its stable logical CPU when supported.
    pub(crate) fn pin_mapping_worker(&self, ordinal: usize) {
        if let Some(cpu) = self.mapping_cpu(ordinal) {
            platform::set_current_affinity(&[cpu]);
        }
    }

    fn mapping_cpu(&self, ordinal: usize) -> Option<usize> {
        (!self.mapping.is_empty()).then(|| self.mapping[ordinal % self.mapping.len()])
    }

    /// Pins a streaming role (FASTQ or BAM) to the non-mapping CPU pool.
    pub(crate) fn pin_auxiliary_worker(&self) {
        platform::set_current_affinity(&self.auxiliary);
    }

    /// Pins the calling coordinator temporarily and restores its inherited set
    /// when the returned guard is dropped.
    #[must_use]
    pub(crate) fn pin_auxiliary_scoped(&self) -> CurrentAffinityGuard {
        let original = platform::current_affinity();
        platform::set_current_affinity(&self.auxiliary);
        CurrentAffinityGuard { original }
    }
}

/// Restores the coordinator CPU set after one alignment transaction.
pub(crate) struct CurrentAffinityGuard {
    original: Vec<usize>,
}

impl Drop for CurrentAffinityGuard {
    fn drop(&mut self) {
        platform::set_current_affinity(&self.original);
    }
}

#[cfg(target_os = "linux")]
#[allow(unsafe_code)]
mod platform {
    use super::{AvailableCpu, PhysicalCore};
    use std::path::PathBuf;

    const BITS_PER_WORD: usize = core::mem::size_of::<usize>() * 8;
    const INITIAL_AFFINITY_WORDS: usize =
        core::mem::size_of::<libc::cpu_set_t>().div_ceil(core::mem::size_of::<usize>());
    // This is only a guard against an unexpectedly unbounded retry. One MiB
    // represents millions of logical CPUs on the supported 64-bit targets.
    const MAX_AFFINITY_WORDS: usize = (1024 * 1024) / core::mem::size_of::<usize>();

    pub(super) fn available_cpus() -> Vec<AvailableCpu> {
        current_affinity()
            .into_iter()
            .map(|logical| AvailableCpu {
                logical,
                physical: read_physical_core(logical),
            })
            .collect()
    }

    pub(super) fn current_affinity() -> Vec<usize> {
        read_current_affinity_words().map_or_else(Vec::new, |words| decode_affinity_words(&words))
    }

    pub(super) fn set_current_affinity(cpus: &[usize]) {
        if cpus.is_empty() {
            return;
        }
        let Some(current) = read_current_affinity_words() else {
            return;
        };
        let Some(set) = encode_affinity_words(cpus, current.len()) else {
            return;
        };
        let byte_len = set.len() * core::mem::size_of::<usize>();
        // Best effort by design: an inherited container/job affinity policy is
        // authoritative, and alignment correctness never depends on placement.
        // SAFETY: PID 0 selects the calling thread. `set` is suitably aligned,
        // initialized, and readable for the supplied dynamic mask size.
        let _ =
            unsafe { libc::sched_setaffinity(0, byte_len, set.as_ptr().cast::<libc::cpu_set_t>()) };
    }

    fn read_current_affinity_words() -> Option<Vec<usize>> {
        let mut word_count = INITIAL_AFFINITY_WORDS.max(1);
        loop {
            let mut words = vec![0_usize; word_count];
            let byte_len = words.len() * core::mem::size_of::<usize>();
            // Linux accepts a dynamically sized CPU mask. `Vec<usize>` gives
            // the native-word alignment used by the kernel mask, and the
            // supplied byte length prevents access beyond the allocation.
            // SAFETY: PID 0 selects the calling thread; `words` is writable
            // and live for exactly `byte_len` bytes.
            let status = unsafe {
                libc::sched_getaffinity(0, byte_len, words.as_mut_ptr().cast::<libc::cpu_set_t>())
            };
            if status == 0 {
                return Some(words);
            }
            if std::io::Error::last_os_error().raw_os_error() != Some(libc::EINVAL)
                || word_count >= MAX_AFFINITY_WORDS
            {
                return None;
            }
            word_count = word_count.saturating_mul(2).min(MAX_AFFINITY_WORDS);
        }
    }

    fn decode_affinity_words(words: &[usize]) -> Vec<usize> {
        let mut cpus = Vec::new();
        for (word_index, &word) in words.iter().enumerate() {
            let mut remaining = word;
            while remaining != 0 {
                let bit = remaining.trailing_zeros() as usize;
                cpus.push(word_index * BITS_PER_WORD + bit);
                remaining &= remaining - 1;
            }
        }
        cpus
    }

    fn encode_affinity_words(cpus: &[usize], minimum_words: usize) -> Option<Vec<usize>> {
        let required_words = cpus
            .iter()
            .copied()
            .max()
            .map_or(0, |cpu| cpu / BITS_PER_WORD + 1);
        let word_count = minimum_words.max(required_words);
        if word_count > MAX_AFFINITY_WORDS {
            return None;
        }
        let mut words = vec![0_usize; word_count];
        for &cpu in cpus {
            words[cpu / BITS_PER_WORD] |= 1_usize << (cpu % BITS_PER_WORD);
        }
        Some(words)
    }

    fn read_physical_core(logical: usize) -> PhysicalCore {
        let mut root = PathBuf::from("/sys/devices/system/cpu");
        root.push(format!("cpu{logical}"));
        root.push("topology");
        let package = read_integer(root.join("physical_package_id"));
        let core = read_integer(root.join("core_id"));
        match (package, core) {
            (Some(package), Some(core)) => PhysicalCore { package, core },
            _ => PhysicalCore {
                package: -1,
                core: i64::try_from(logical).unwrap_or(i64::MAX),
            },
        }
    }

    fn read_integer(path: PathBuf) -> Option<i64> {
        std::fs::read_to_string(path).ok()?.trim().parse().ok()
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn dynamic_masks_round_trip_large_cpu_ids() {
            let high_cpu = 2_049;
            let cpus = [0, BITS_PER_WORD - 1, BITS_PER_WORD, high_cpu];
            let words = encode_affinity_words(&cpus, INITIAL_AFFINITY_WORDS).unwrap();

            assert_eq!(decode_affinity_words(&words), cpus);
        }
    }
}

#[cfg(not(target_os = "linux"))]
mod platform {
    use super::AvailableCpu;

    pub(super) fn available_cpus() -> Vec<AvailableCpu> {
        Vec::new()
    }

    pub(super) fn current_affinity() -> Vec<usize> {
        Vec::new()
    }

    pub(super) fn set_current_affinity(_cpus: &[usize]) {}
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distinct_mapping_cores_leave_complete_cores_for_streaming() {
        let available = [
            cpu(0, 0),
            cpu(4, 0),
            cpu(1, 1),
            cpu(5, 1),
            cpu(2, 2),
            cpu(6, 2),
        ];
        let placement = CpuPlacement::from_available(&available, 2);
        assert_eq!(placement.mapping, vec![0, 1]);
        assert_eq!(placement.auxiliary, vec![2, 6]);
    }

    #[test]
    fn insufficient_cores_use_smt_lanes_before_sharing() {
        let available = [cpu(2, 0), cpu(3, 0), cpu(6, 1), cpu(7, 1)];
        let placement = CpuPlacement::from_available(&available, 4);
        assert_eq!(placement.mapping, vec![2, 6, 3, 7]);
        assert_eq!(placement.auxiliary, vec![2, 3, 6, 7]);
    }

    #[test]
    fn oversubscribed_workers_cycle_a_bounded_cpu_assignment() {
        let available = [cpu(2, 0), cpu(6, 1)];
        let placement = CpuPlacement::from_available(&available, 1_000_000);
        assert_eq!(placement.mapping, vec![2, 6]);
        assert_eq!(placement.mapping_cpu(0), Some(2));
        assert_eq!(placement.mapping_cpu(1), Some(6));
        assert_eq!(placement.mapping_cpu(2), Some(2));
        assert_eq!(placement.mapping_cpu(999_999), Some(6));
    }

    fn cpu(logical: usize, core: i64) -> AvailableCpu {
        AvailableCpu {
            logical,
            physical: PhysicalCore { package: 0, core },
        }
    }
}
