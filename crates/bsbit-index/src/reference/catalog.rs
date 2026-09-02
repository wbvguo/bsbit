//! Ordered reference-catalog inputs, limits, and validation.

use bsbit_core::sequence::NormalizedSequence;

use super::{
    ReferenceBuildError, ReferenceResource, apply_limit, checked_build_add, physical_to_logical,
};

/// One owned contig supplied to reference construction.
#[derive(Clone, Debug)]
pub struct ContigInput {
    pub(crate) name: Vec<u8>,
    pub(crate) sequence: NormalizedSequence,
}

impl ContigInput {
    /// Creates one owned contig.
    #[must_use]
    pub const fn new(name: Vec<u8>, sequence: NormalizedSequence) -> Self {
        Self { name, sequence }
    }

    /// Returns the exact contig name bytes.
    #[must_use]
    pub fn name(&self) -> &[u8] {
        &self.name
    }

    /// Returns the retained normalized sequence.
    #[must_use]
    pub const fn sequence(&self) -> &NormalizedSequence {
        &self.sequence
    }
}

/// Aggregate dimensions of a validated ordered reference catalog before any
/// search index is constructed.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ReferenceCatalogMetrics {
    contig_count: u64,
    total_name_bytes: u64,
    total_reference_bases: u64,
}

impl ReferenceCatalogMetrics {
    /// Returns the number of ordered contigs.
    #[must_use]
    pub const fn contig_count(self) -> u64 {
        self.contig_count
    }

    /// Returns aggregate exact contig-name bytes.
    #[must_use]
    pub const fn total_name_bytes(self) -> u64 {
        self.total_name_bytes
    }

    /// Returns aggregate normalized bases including `N`.
    #[must_use]
    pub const fn total_reference_bases(self) -> u64 {
        self.total_reference_bases
    }
}

/// Limits for catalog validation that do not model or construct FM lanes.
#[allow(clippy::struct_field_names)]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ReferenceCatalogLimits {
    max_contigs: u64,
    max_total_name_bytes: u64,
    max_total_reference_bases: u64,
}

impl ReferenceCatalogLimits {
    /// Limits admitting every representable catalog dimension.
    pub const MAX: Self = Self {
        max_contigs: u64::MAX,
        max_total_name_bytes: u64::MAX,
        max_total_reference_bases: u64::MAX,
    };

    /// Sets the maximum ordered contig count.
    #[must_use]
    pub const fn with_max_contigs(mut self, value: u64) -> Self {
        self.max_contigs = value;
        self
    }

    /// Sets the maximum aggregate exact name bytes.
    #[must_use]
    pub const fn with_max_total_name_bytes(mut self, value: u64) -> Self {
        self.max_total_name_bytes = value;
        self
    }

    /// Sets the maximum aggregate normalized bases.
    #[must_use]
    pub const fn with_max_total_reference_bases(mut self, value: u64) -> Self {
        self.max_total_reference_bases = value;
        self
    }
}

impl Default for ReferenceCatalogLimits {
    fn default() -> Self {
        Self::MAX
    }
}

/// Validates ordered catalog semantics and dimensions without constructing FM
/// lanes or allocating a reference owner.
///
/// # Errors
///
/// Returns the same catalog-prefix validation errors and priority used by
/// [`ReferenceIndex::build`].
pub fn validate_reference_catalog(
    contigs: &[ContigInput],
    limits: ReferenceCatalogLimits,
) -> Result<ReferenceCatalogMetrics, ReferenceBuildError> {
    validate_catalog_and_measure(contigs, limits)
}

pub(super) fn validate_catalog_and_measure(
    contigs: &[ContigInput],
    limits: ReferenceCatalogLimits,
) -> Result<ReferenceCatalogMetrics, ReferenceBuildError> {
    if contigs.is_empty() {
        return Err(ReferenceBuildError::EmptyReference);
    }

    let contig_count = physical_to_logical(contigs.len(), ReferenceResource::Contigs)?;
    apply_limit(ReferenceResource::Contigs, contig_count, limits.max_contigs)?;

    let mut total_name_bytes = 0_u64;
    for contig in contigs {
        let name_len = physical_to_logical(contig.name.len(), ReferenceResource::TotalNameBytes)?;
        total_name_bytes = checked_build_add(
            ReferenceResource::TotalNameBytes,
            total_name_bytes,
            name_len,
        )?;
    }
    apply_limit(
        ReferenceResource::TotalNameBytes,
        total_name_bytes,
        limits.max_total_name_bytes,
    )?;

    for (duplicate_storage, contig) in contigs.iter().enumerate() {
        let duplicate_ordinal = physical_to_logical(duplicate_storage, ReferenceResource::Contigs)?;
        if contig.name.is_empty() {
            return Err(ReferenceBuildError::EmptyContigName {
                contig_ordinal: duplicate_ordinal,
            });
        }
        for (first_storage, prior) in contigs.iter().take(duplicate_storage).enumerate() {
            if prior.name == contig.name {
                return Err(ReferenceBuildError::DuplicateContigName {
                    first_ordinal: physical_to_logical(first_storage, ReferenceResource::Contigs)?,
                    duplicate_ordinal,
                });
            }
        }
        if contig.sequence.is_empty() {
            return Err(ReferenceBuildError::EmptyContigSequence {
                contig_ordinal: duplicate_ordinal,
            });
        }
    }

    let mut total_reference_bases = 0_u64;
    for contig in contigs {
        total_reference_bases = checked_build_add(
            ReferenceResource::TotalReferenceBases,
            total_reference_bases,
            contig.sequence.len(),
        )?;
    }
    apply_limit(
        ReferenceResource::TotalReferenceBases,
        total_reference_bases,
        limits.max_total_reference_bases,
    )?;

    Ok(ReferenceCatalogMetrics {
        contig_count,
        total_name_bytes,
        total_reference_bases,
    })
}
