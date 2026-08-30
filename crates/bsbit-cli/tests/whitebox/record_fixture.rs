//! Small exhaustive fixture mapper used only by record-format unit tests.

use std::collections::BTreeSet;

use bsbit_align::extension::VerifiedAlignment;
use bsbit_align::materialize::traceback_read_placement;
use bsbit_core::bisulfite::BisulfiteStrand;
use bsbit_core::coordinate::{ReferenceInterval, ReferenceLength};
use bsbit_core::sequence::{NormalizedSequence, normalize_dna};
use bsbit_hts::RecordMappingQuality;
use bsbit_index::reference::ReferenceIndex;

pub(super) struct SingleFixture {
    pub(super) query: NormalizedSequence,
    pub(super) alignment: Option<VerifiedAlignment>,
    pub(super) mapping_quality: RecordMappingQuality,
}

/// Exhaustively enumerates the tiny in-memory references used by format
/// tests. This is deliberately test-only: product single-end mapping always
/// uses `SingleBatchAligner` through `bsbit align`.
pub(super) fn single_fixture(
    reference: &ReferenceIndex,
    raw: &[u8],
    maximum_edit_distance: u64,
) -> SingleFixture {
    let query = normalize_dna(raw).expect("record fixture is normalized DNA");
    let budget = u8::try_from(maximum_edit_distance).expect("record fixture budget fits u8");
    let mut best_distance = None;
    let mut best = Vec::new();
    for contig_ordinal in 0..reference.contig_count() {
        let contig = reference
            .contig_by_ordinal(contig_ordinal)
            .expect("record fixture contig exists");
        let contig_id = reference
            .contig_id(contig_ordinal)
            .expect("record fixture contig id exists");
        let contig_length = contig.sequence().len();
        for start in 0..contig_length {
            for end in start + 1..=contig_length {
                let interval =
                    ReferenceInterval::new(start, end, ReferenceLength::new(contig_length))
                        .expect("enumerated record fixture interval is bounded");
                for strand in [BisulfiteStrand::OT, BisulfiteStrand::OB] {
                    let Ok(alignment) = traceback_read_placement(
                        reference, &query, &contig_id, interval, strand, budget,
                    ) else {
                        continue;
                    };
                    let distance = alignment.distance().get();
                    match best_distance {
                        None => {
                            best_distance = Some(distance);
                            best.push(alignment);
                        }
                        Some(current) if distance < current => {
                            best_distance = Some(distance);
                            best.clear();
                            best.push(alignment);
                        }
                        Some(current) if distance == current => best.push(alignment),
                        Some(_) => {}
                    }
                }
            }
        }
    }

    let placement_count = best
        .iter()
        .map(|alignment| {
            (
                alignment.contig().ordinal(),
                alignment.interval().start(),
                alignment.interval().end(),
                alignment.orientation(),
            )
        })
        .collect::<BTreeSet<_>>()
        .len();
    let mapping_quality = match placement_count {
        0 => RecordMappingQuality::Unmapped,
        1 => RecordMappingQuality::Unavailable,
        _ => RecordMappingQuality::Tied,
    };
    SingleFixture {
        query,
        alignment: best.into_iter().next(),
        mapping_quality,
    }
}
