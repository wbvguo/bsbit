//! Shared ordered site identity and per-sample counts.

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub(crate) struct SiteKey {
    pub(crate) contig: u32,
    pub(crate) start: u64,
    pub(crate) end: u64,
    pub(crate) strand: u8,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct Counts {
    pub(crate) methylated: u64,
    pub(crate) total: u64,
}
