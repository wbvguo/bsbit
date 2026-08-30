//! Best-effort internal CPU placement for alignment pipeline roles.
//!
//! Placement is deliberately not part of the command-line contract.  The
//! inherited CPU set remains the authority (including `taskset`, containers,
//! and job schedulers); this module only partitions that set so mapping and
//! streaming workers do not compete for the same physical cores.

#![allow(unsafe_code)]

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

        let mut represented_cores = BTreeSet::new();
        let mut mapping = Vec::with_capacity(mapping_workers);
        for cpu in available {
            if represented_cores.insert(cpu.physical) {
                mapping.push(cpu.logical);
                if mapping.len() == mapping_workers {
                    break;
                }
            }
        }

        // When workers outnumber physical cores, consume otherwise unused SMT
        // lanes before sharing a logical CPU.
        if mapping.len() < mapping_workers {
            for cpu in available {
                if !mapping.contains(&cpu.logical) {
                    mapping.push(cpu.logical);
                    if mapping.len() == mapping_workers {
                        break;
                    }
                }
            }
        }
        if mapping.len() < mapping_workers {
            let distinct = mapping.len();
            for ordinal in distinct..mapping_workers {
                mapping.push(mapping[ordinal % distinct]);
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
        if let Some(&cpu) = self.mapping.get(ordinal) {
            platform::set_current_affinity(&[cpu]);
        }
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
mod platform {
    use super::{AvailableCpu, PhysicalCore};
    use std::path::PathBuf;

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
        // SAFETY: an all-zero bit pattern is a valid empty `cpu_set_t`.
        let mut set = unsafe { core::mem::zeroed::<libc::cpu_set_t>() };
        // SAFETY: PID 0 selects the calling thread and `set` is writable for
        // exactly the supplied `cpu_set_t` size.
        let status = unsafe {
            libc::sched_getaffinity(0, core::mem::size_of::<libc::cpu_set_t>(), &raw mut set)
        };
        if status != 0 {
            return Vec::new();
        }
        (0..libc::CPU_SETSIZE as usize)
            .filter(|&cpu| {
                // SAFETY: `cpu` is within `CPU_SETSIZE` and `set` is live.
                unsafe { libc::CPU_ISSET(cpu, &set) }
            })
            .collect()
    }

    pub(super) fn set_current_affinity(cpus: &[usize]) {
        if cpus.is_empty() {
            return;
        }
        // SAFETY: an all-zero bit pattern is a valid empty `cpu_set_t`.
        let mut set = unsafe { core::mem::zeroed::<libc::cpu_set_t>() };
        // SAFETY: `set` is a live, writable `cpu_set_t`.
        unsafe { libc::CPU_ZERO(&mut set) };
        for &cpu in cpus {
            if cpu < libc::CPU_SETSIZE as usize {
                // SAFETY: `cpu` is in range and `set` remains live.
                unsafe { libc::CPU_SET(cpu, &mut set) };
            }
        }
        // Best effort by design: an inherited container/job affinity policy is
        // authoritative, and alignment correctness never depends on placement.
        // SAFETY: PID 0 selects the calling thread and `set` is initialized.
        let _ = unsafe {
            libc::sched_setaffinity(0, core::mem::size_of::<libc::cpu_set_t>(), &raw const set)
        };
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

    fn cpu(logical: usize, core: i64) -> AvailableCpu {
        AvailableCpu {
            logical,
            physical: PhysicalCore { package: 0, core },
        }
    }
}
