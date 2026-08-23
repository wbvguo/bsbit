//! Bounded indexed-FASTA context lookup for methylation calling.

use std::path::Path;

use bsbit_core::reference::{ReferenceSemanticDigest, ReferenceSemanticDigestBuilder};
use bsbit_hts::IndexedFastaReader;

use crate::call_input::BamReference;
use crate::evidence::{ContextClass, CytosineContext, EvidenceStrand};
use crate::region::CallRegion;
use crate::{CallError, CallErrorKind};

pub(crate) struct CallReferenceReader {
    reader: IndexedFastaReader,
    fasta_reference_ids: Vec<u32>,
}

impl CallReferenceReader {
    pub(crate) fn open(path: &Path, references: &[BamReference]) -> Result<Self, CallError> {
        let reader = IndexedFastaReader::open(path).map_err(|error| {
            CallError::with_source(
                CallErrorKind::Input,
                format!("open indexed reference FASTA {}", path.display()),
                error,
            )
        })?;
        let mut fasta_reference_ids = Vec::new();
        fasta_reference_ids
            .try_reserve_exact(references.len())
            .map_err(|error| {
                CallError::with_source(
                    CallErrorKind::Input,
                    "reserve reference dictionary mapping",
                    error,
                )
            })?;
        for (bam_id, bam_reference) in references.iter().enumerate() {
            let matches = reader
                .references()
                .iter()
                .enumerate()
                .filter(|(_, fasta_reference)| fasta_reference.name() == bam_reference.name)
                .collect::<Vec<_>>();
            let [(fasta_id, fasta_reference)] = matches.as_slice() else {
                return Err(CallError::input(if matches.is_empty() {
                    format!(
                        "reference FASTA {} is missing BAM contig {} (`{}`)",
                        path.display(),
                        bam_id,
                        String::from_utf8_lossy(&bam_reference.name)
                    )
                } else {
                    format!(
                        "reference FASTA {} contains duplicate contig `{}`",
                        path.display(),
                        String::from_utf8_lossy(&bam_reference.name)
                    )
                }));
            };
            if fasta_reference.length() != u64::from(bam_reference.length) {
                return Err(CallError::input(format!(
                    "reference FASTA {} contig `{}` length {} differs from BAM dictionary length {}",
                    path.display(),
                    String::from_utf8_lossy(&bam_reference.name),
                    fasta_reference.length(),
                    bam_reference.length
                )));
            }
            fasta_reference_ids.push(
                u32::try_from(*fasta_id).map_err(|_| {
                    CallError::input("reference FASTA dictionary ordinal exceeds u32")
                })?,
            );
        }
        Ok(Self {
            reader,
            fasta_reference_ids,
        })
    }

    pub(crate) fn fetch_context_window(
        &mut self,
        region: CallRegion,
        references: &[BamReference],
    ) -> Result<ReferenceWindow, CallError> {
        let reference = references
            .get(usize::try_from(region.reference).expect("u32 fits usize"))
            .ok_or_else(|| CallError::operation("calling region references a missing contig"))?;
        let fasta_reference_id = *self
            .fasta_reference_ids
            .get(usize::try_from(region.reference).expect("u32 fits usize"))
            .ok_or_else(|| CallError::operation("reference dictionary mapping is incomplete"))?;
        let start = region.start.saturating_sub(2);
        let end = region.end.saturating_add(2).min(reference.length);
        let mut bases = self
            .reader
            .fetch(fasta_reference_id, u64::from(start), u64::from(end))
            .map_err(|error| {
                CallError::with_source(
                    CallErrorKind::Input,
                    format!(
                        "fetch reference FASTA contig `{}` interval {}-{}",
                        String::from_utf8_lossy(&reference.name),
                        start,
                        end
                    ),
                    error,
                )
            })?;
        bases.make_ascii_uppercase();
        Ok(ReferenceWindow {
            reference: region.reference,
            start,
            bases,
        })
    }

    pub(crate) fn validate_semantic_digest(
        &mut self,
        references: &[BamReference],
        expected: ReferenceSemanticDigest,
    ) -> Result<(), CallError> {
        let contig_count = u64::try_from(references.len())
            .map_err(|_| CallError::input("BAM reference count exceeds u64"))?;
        let mut builder = ReferenceSemanticDigestBuilder::new(contig_count);
        for (ordinal, reference) in references.iter().enumerate() {
            let fasta_reference_id = *self.fasta_reference_ids.get(ordinal).ok_or_else(|| {
                CallError::operation("reference dictionary mapping is incomplete")
            })?;
            builder
                .begin_ascii_contig(&reference.name, u64::from(reference.length))
                .map_err(|error| {
                    CallError::with_source(
                        CallErrorKind::Input,
                        format!(
                            "normalize reference FASTA contig `{}` for provenance validation",
                            String::from_utf8_lossy(&reference.name)
                        ),
                        error,
                    )
                })?;
            let mut start = 0_u64;
            let length = u64::from(reference.length);
            while start < length {
                let end = start.saturating_add(8 * 1024 * 1024).min(length);
                let bases = self
                    .reader
                    .fetch(fasta_reference_id, start, end)
                    .map_err(|error| {
                        CallError::with_source(
                            CallErrorKind::Input,
                            format!(
                                "read reference FASTA contig `{}` interval {start}-{end} for provenance validation",
                                String::from_utf8_lossy(&reference.name)
                            ),
                            error,
                        )
                    })?;
                builder.push_ascii_bases(&bases).map_err(|error| {
                    CallError::with_source(
                        CallErrorKind::Input,
                        format!(
                            "normalize reference FASTA contig `{}` for provenance validation",
                            String::from_utf8_lossy(&reference.name)
                        ),
                        error,
                    )
                })?;
                start = end;
            }
            builder.end_ascii_contig().map_err(|error| {
                CallError::with_source(
                    CallErrorKind::Input,
                    format!(
                        "finish reference FASTA contig `{}` provenance validation",
                        String::from_utf8_lossy(&reference.name)
                    ),
                    error,
                )
            })?;
        }
        let observed = builder.finish().map_err(|error| {
            CallError::with_source(
                CallErrorKind::Input,
                "finish reference FASTA provenance validation",
                error,
            )
        })?;
        if observed != expected {
            return Err(CallError::input(format!(
                "reference FASTA semantic digest {observed} differs from BAM provenance {expected}"
            )));
        }
        Ok(())
    }

    pub(crate) fn close(self) -> Result<(), CallError> {
        self.reader.close().map_err(|error| {
            CallError::with_source(CallErrorKind::Input, "close indexed reference FASTA", error)
        })
    }
}

pub(crate) struct ReferenceWindow {
    reference: u32,
    start: u32,
    bases: Vec<u8>,
}

impl ReferenceWindow {
    pub(crate) fn base(&self, reference: u32, position: u32) -> Option<u8> {
        if reference != self.reference {
            return None;
        }
        let offset = usize::try_from(position.checked_sub(self.start)?).ok()?;
        self.bases.get(offset).copied()
    }

    pub(crate) fn context(
        &self,
        reference: u32,
        position: u32,
        strand: EvidenceStrand,
    ) -> Option<CytosineContext> {
        let (first, second) = match strand {
            EvidenceStrand::Top => {
                let first = canonical(self.base(reference, position.checked_add(1)?)?)?;
                if first == b'G' {
                    return Some(CytosineContext {
                        class: ContextClass::Cg,
                        second: first,
                    });
                }
                let second = canonical(self.base(reference, position.checked_add(2)?)?)?;
                (first, second)
            }
            EvidenceStrand::Bottom => {
                let first = complement(self.base(reference, position.checked_sub(1)?)?)?;
                if first == b'G' {
                    return Some(CytosineContext {
                        class: ContextClass::Cg,
                        second: first,
                    });
                }
                let second = complement(self.base(reference, position.checked_sub(2)?)?)?;
                (first, second)
            }
        };
        Some(CytosineContext {
            class: if second == b'G' {
                ContextClass::Chg
            } else {
                ContextClass::Chh
            },
            second: first,
        })
    }
}

const fn canonical(base: u8) -> Option<u8> {
    match base {
        b'A' | b'C' | b'G' | b'T' => Some(base),
        _ => None,
    }
}

const fn complement(base: u8) -> Option<u8> {
    match base {
        b'A' => Some(b'T'),
        b'C' => Some(b'G'),
        b'G' => Some(b'C'),
        b'T' => Some(b'A'),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::ReferenceWindow;
    use crate::evidence::{ContextClass, CytosineContext, EvidenceStrand};

    #[test]
    fn reference_window_resolves_context_on_both_strands_and_edges() {
        let window = ReferenceWindow {
            reference: 0,
            start: 8,
            bases: b"ACGCTGG".to_vec(),
        };
        assert_eq!(
            window.context(0, 9, EvidenceStrand::Top),
            Some(CytosineContext {
                class: ContextClass::Cg,
                second: b'G',
            })
        );
        assert_eq!(
            window.context(0, 13, EvidenceStrand::Bottom),
            Some(CytosineContext {
                class: ContextClass::Chg,
                second: b'A',
            })
        );
        assert_eq!(window.context(0, 8, EvidenceStrand::Bottom), None);
        assert_eq!(window.context(1, 9, EvidenceStrand::Top), None);
    }
}
