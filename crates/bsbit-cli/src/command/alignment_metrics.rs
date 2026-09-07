//! Shared timing and stable alignment-metrics value rendering.

use std::time::Instant;

use bsbit_align::AlignmentOutputPolicy;
use bsbit_align::library::LibraryProfile;
use bsbit_align::paired_end::PairedSearchMode;
use bsbit_hts::AlignmentAuxiliaryMode;

use super::align_options::{Options, ReadOutputMode, SearchMode};

#[derive(Clone, Copy)]
pub(super) struct MetricsTimer(pub(super) Option<Instant>);

impl MetricsTimer {
    pub(super) fn start(enabled: bool) -> Self {
        Self(enabled.then(Instant::now))
    }

    pub(super) fn elapsed_ns(self) -> u128 {
        self.0.map_or(0, |started| started.elapsed().as_nanos())
    }
}

/// One schema-checked metrics row whose names and values are declared together.
pub(super) struct MetricsRow<const N: usize> {
    fields: [(&'static str, String); N],
}

impl<const N: usize> MetricsRow<N> {
    pub(super) const fn new(fields: [(&'static str, String); N]) -> Self {
        Self { fields }
    }

    pub(super) fn header(&self) -> String {
        self.fields
            .iter()
            .map(|(name, _)| *name)
            .collect::<Vec<_>>()
            .join("\t")
    }

    pub(super) fn values(&self) -> String {
        self.fields
            .iter()
            .map(|(_, value)| value.as_str())
            .collect::<Vec<_>>()
            .join("\t")
    }
}

pub(super) const fn output_contract_name(mode: AlignmentAuxiliaryMode) -> &'static str {
    match mode {
        AlignmentAuxiliaryMode::Minimal => "minimal",
        AlignmentAuxiliaryMode::Bismark => "bismark",
    }
}

pub(super) const fn library_profile_name(profile: LibraryProfile) -> &'static str {
    match profile {
        LibraryProfile::Directional => "directional",
        LibraryProfile::NonDirectional => "non-directional",
    }
}

pub(super) fn soft_clip_fallback_name(
    mode: PairedSearchMode,
    policy: &AlignmentOutputPolicy,
) -> &'static str {
    match policy.soft_clip_mode() {
        bsbit_align::SoftClipMode::None => "none",
        bsbit_align::SoftClipMode::Adapter => {
            if policy.adapter().is_some()
                && policy.adapter_maximum_clip_bases() != 0
                && policy.maximum_soft_clip_bases() != 0
            {
                "adapter"
            } else {
                "none"
            }
        }
        bsbit_align::SoftClipMode::Auto => match mode {
            PairedSearchMode::Sensitive if policy.maximum_soft_clip_bases() != 0 => "semi-global",
            PairedSearchMode::Default
                if policy.adapter().is_some()
                    && policy.adapter_maximum_clip_bases() != 0
                    && policy.maximum_soft_clip_bases() != 0 =>
            {
                "adapter"
            }
            PairedSearchMode::Default | PairedSearchMode::Sensitive => "none",
        },
    }
}

pub(super) const fn soft_clip_mode_name(mode: bsbit_align::SoftClipMode) -> &'static str {
    match mode {
        bsbit_align::SoftClipMode::Auto => "auto",
        bsbit_align::SoftClipMode::None => "none",
        bsbit_align::SoftClipMode::Adapter => "adapter",
    }
}

pub(super) fn adapter_name(policy: &AlignmentOutputPolicy) -> String {
    let Some(adapter) = policy.adapter() else {
        return "none".to_owned();
    };
    let default = AlignmentOutputPolicy::default();
    if default.adapter() == Some(adapter) {
        "illumina".to_owned()
    } else {
        String::from_utf8(adapter.to_vec()).expect("validated adapters contain only ASCII")
    }
}

pub(super) const fn mate_rescue_name(mode: PairedSearchMode) -> &'static str {
    match mode {
        PairedSearchMode::Default => "off",
        PairedSearchMode::Sensitive => "windowed",
    }
}

pub(super) const fn search_mode_name(mode: SearchMode) -> &'static str {
    match mode {
        SearchMode::Default => "default",
        SearchMode::Sensitive => "sensitive",
    }
}

pub(super) const fn read_output_name(mode: ReadOutputMode) -> &'static str {
    match mode {
        ReadOutputMode::Complete => "complete",
        ReadOutputMode::MappedOnly => "mapped-only",
    }
}

pub(super) const fn sensitive_mapq_zero_strategy_id() -> &'static str {
    "sensitive-bounded-integrated-mapq0-hash-tie-v1"
}

pub(super) const fn sensitive_read_complete_strategy_id() -> &'static str {
    "sensitive-bounded-integrated-read-complete-hash-tie-v1"
}

pub(super) fn strategy_id(options: &Options) -> &'static str {
    // These identifiers establish the pre-release v1 baseline. After the
    // first public release, change an identifier only when its observable
    // alignment strategy changes; behavior-preserving refactors retain it.
    match (options.search_mode, options.read_output) {
        (SearchMode::Sensitive, ReadOutputMode::Complete) => sensitive_read_complete_strategy_id(),
        (SearchMode::Sensitive, ReadOutputMode::MappedOnly) => sensitive_mapq_zero_strategy_id(),
        (SearchMode::Default, ReadOutputMode::Complete) => {
            "balanced-d5-adapter-recovery-read-complete-hash-tie-v1"
        }
        (SearchMode::Default, ReadOutputMode::MappedOnly) => {
            "balanced-d5-adapter-recovery-mapq0-hash-tie-v1"
        }
    }
}
