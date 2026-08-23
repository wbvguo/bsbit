//! Dictionary-aware region planning shared by all calling entry points.

use std::collections::HashMap;
use std::io::{BufRead, BufReader};

use bsbit_hts::{BedInterval, DecodedReader};

use super::{GenomicInterval, RegionSelection};
use crate::call_input::BamReference;
use crate::{CallError, CallErrorKind};

/// One bounded unit of reference-coordinate work.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct CallRegion {
    pub(crate) ordinal: usize,
    pub(crate) reference: u32,
    pub(crate) start: u32,
    pub(crate) end: u32,
}

pub(crate) fn plan_call_regions(
    references: &[BamReference],
    region_bases: u32,
    selection: &RegionSelection,
) -> Result<Vec<CallRegion>, CallError> {
    if region_bases == 0 {
        return Err(CallError::configuration(
            "calling region length must be nonzero",
        ));
    }
    let restricted = !selection.intervals.is_empty() || selection.regions_file.is_some();
    let mut requested = selection.intervals.clone();
    if let Some(path) = selection.regions_file.as_deref() {
        requested.extend(read_regions_file(path)?);
    }
    if restricted && requested.is_empty() {
        return Err(CallError::configuration(
            "explicit calling region selection contains no intervals",
        ));
    }

    let mut reference_by_name = HashMap::with_capacity(references.len());
    for (ordinal, reference) in references.iter().enumerate() {
        if reference_by_name
            .insert(reference.name.as_slice(), ordinal)
            .is_some()
        {
            return Err(CallError::input(format!(
                "BAM dictionary repeats reference `{}`",
                String::from_utf8_lossy(&reference.name)
            )));
        }
    }

    let mut selected = vec![Vec::<(u32, u32)>::new(); references.len()];
    if restricted {
        for interval in requested {
            validate_and_add_interval(&interval, references, &reference_by_name, &mut selected)?;
        }
    } else {
        for (ordinal, reference) in references.iter().enumerate() {
            if reference.length != 0 {
                selected[ordinal].push((0, reference.length));
            }
        }
    }

    let mut regions = Vec::new();
    for (reference_ordinal, intervals) in selected.iter_mut().enumerate() {
        intervals.sort_unstable();
        let merged = merge_intervals(intervals);
        let reference = u32::try_from(reference_ordinal)
            .map_err(|_| CallError::input("BAM reference ordinal exceeds u32"))?;
        for (start, end) in merged {
            let mut chunk_start = start;
            while chunk_start < end {
                let chunk_end = chunk_start.saturating_add(region_bases).min(end);
                regions.push(CallRegion {
                    ordinal: regions.len(),
                    reference,
                    start: chunk_start,
                    end: chunk_end,
                });
                chunk_start = chunk_end;
            }
        }
    }
    if regions.is_empty() {
        return Err(CallError::configuration(
            "calling region selection contains no nonempty reference bases",
        ));
    }
    Ok(regions)
}

fn validate_and_add_interval(
    interval: &GenomicInterval,
    references: &[BamReference],
    reference_by_name: &HashMap<&[u8], usize>,
    selected: &mut [Vec<(u32, u32)>],
) -> Result<(), CallError> {
    if interval.contig.is_empty()
        || interval
            .contig
            .bytes()
            .any(|byte| byte.is_ascii_control() || byte.is_ascii_whitespace())
    {
        return Err(CallError::configuration(
            "calling region contig must be nonempty and contain no whitespace or control bytes",
        ));
    }
    if interval.start >= interval.end {
        return Err(CallError::configuration(format!(
            "calling region {}:{}-{} must be nonempty",
            interval.contig, interval.start, interval.end
        )));
    }
    let Some(&reference_ordinal) = reference_by_name.get(interval.contig.as_bytes()) else {
        return Err(CallError::input(format!(
            "calling region names reference `{}` absent from the BAM dictionary",
            interval.contig
        )));
    };
    let reference_length = u64::from(references[reference_ordinal].length);
    if interval.end > reference_length {
        return Err(CallError::input(format!(
            "calling region {}:{}-{} exceeds BAM reference length {reference_length}",
            interval.contig, interval.start, interval.end
        )));
    }
    let start = u32::try_from(interval.start)
        .map_err(|_| CallError::input("calling region start exceeds u32"))?;
    let end = u32::try_from(interval.end)
        .map_err(|_| CallError::input("calling region end exceeds u32"))?;
    selected[reference_ordinal].push((start, end));
    Ok(())
}

fn merge_intervals(intervals: &[(u32, u32)]) -> Vec<(u32, u32)> {
    let mut merged: Vec<(u32, u32)> = Vec::with_capacity(intervals.len());
    for &(start, end) in intervals {
        if let Some(last) = merged.last_mut()
            && start <= last.1
        {
            last.1 = last.1.max(end);
            continue;
        }
        merged.push((start, end));
    }
    merged
}

fn read_regions_file(path: &std::path::Path) -> Result<Vec<GenomicInterval>, CallError> {
    let decoded = DecodedReader::open(path).map_err(|error| {
        CallError::with_source(
            CallErrorKind::Input,
            format!("open calling regions file {}", path.display()),
            error,
        )
    })?;
    let mut reader = BufReader::new(decoded);
    let parsed = parse_regions_file(&mut reader, path);
    let close = reader.into_inner().close().map_err(|error| {
        CallError::with_source(
            CallErrorKind::Input,
            format!("close calling regions file {}", path.display()),
            error,
        )
    });
    match (parsed, close) {
        (Err(error), _) | (Ok(_), Err(error)) => Err(error),
        (Ok(intervals), Ok(())) => Ok(intervals),
    }
}

fn parse_regions_file(
    reader: &mut impl BufRead,
    path: &std::path::Path,
) -> Result<Vec<GenomicInterval>, CallError> {
    let mut intervals = Vec::new();
    let mut line = Vec::new();
    let mut line_number = 0_u64;
    loop {
        line.clear();
        let read = reader.read_until(b'\n', &mut line).map_err(|error| {
            CallError::with_source(
                CallErrorKind::Input,
                format!("read calling regions file {}", path.display()),
                error,
            )
        })?;
        if read == 0 {
            break;
        }
        line_number += 1;
        while matches!(line.last(), Some(b'\n' | b'\r')) {
            line.pop();
        }
        let Some(interval) = BedInterval::parse_line(&line).map_err(|error| {
            CallError::with_source(
                CallErrorKind::Input,
                format!(
                    "calling regions file {} line {line_number} has invalid BED3+ syntax",
                    path.display()
                ),
                error,
            )
        })?
        else {
            continue;
        };
        let contig = std::str::from_utf8(interval.contig()).map_err(|error| {
            CallError::with_source(
                CallErrorKind::Input,
                format!(
                    "calling regions file {} line {line_number} has a non-UTF-8 contig",
                    path.display()
                ),
                error,
            )
        })?;
        intervals.push(GenomicInterval {
            contig: contig.to_owned(),
            start: interval.start(),
            end: interval.end(),
        });
    }
    Ok(intervals)
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::{GenomicInterval, RegionSelection, plan_call_regions};
    use crate::{CallErrorKind, call_input::BamReference};

    fn unique_directory(label: &str) -> PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!(
            "bsbit-call-regions-{label}-{}-{nonce}",
            std::process::id()
        ))
    }

    fn references() -> Vec<BamReference> {
        vec![
            BamReference {
                name: b"chr1".to_vec(),
                length: 100,
            },
            BamReference {
                name: b"chr2".to_vec(),
                length: 50,
            },
        ]
    }

    #[test]
    fn whole_dictionary_is_chunked_in_dictionary_order() {
        let regions = plan_call_regions(&references(), 64, &RegionSelection::default()).unwrap();
        assert_eq!(
            regions
                .iter()
                .map(|region| (region.ordinal, region.reference, region.start, region.end))
                .collect::<Vec<_>>(),
            vec![(0, 0, 0, 64), (1, 0, 64, 100), (2, 1, 0, 50)]
        );
    }

    #[test]
    fn direct_overlaps_are_merged_and_chunked_once() {
        let selection = RegionSelection {
            intervals: vec![
                GenomicInterval {
                    contig: String::from("chr1"),
                    start: 20,
                    end: 80,
                },
                GenomicInterval {
                    contig: String::from("chr1"),
                    start: 10,
                    end: 30,
                },
                GenomicInterval {
                    contig: String::from("chr2"),
                    start: 5,
                    end: 6,
                },
            ],
            regions_file: None,
        };
        let regions = plan_call_regions(&references(), 40, &selection).unwrap();
        assert_eq!(
            regions
                .iter()
                .map(|region| (region.reference, region.start, region.end))
                .collect::<Vec<_>>(),
            vec![(0, 10, 50), (0, 50, 80), (1, 5, 6)]
        );
    }

    #[test]
    fn unknown_or_out_of_bounds_intervals_fail_closed() {
        for interval in [
            GenomicInterval {
                contig: String::from("chrX"),
                start: 0,
                end: 1,
            },
            GenomicInterval {
                contig: String::from("chr1"),
                start: 99,
                end: 101,
            },
        ] {
            let error = plan_call_regions(
                &references(),
                64,
                &RegionSelection {
                    intervals: vec![interval],
                    regions_file: None,
                },
            )
            .unwrap_err();
            assert!(error.to_string().contains("calling region"));
        }
    }

    #[test]
    fn bed_regions_are_unioned_with_direct_intervals() {
        let directory = unique_directory("union");
        fs::create_dir(&directory).unwrap();
        let bed = directory.join("targets.bed");
        fs::write(
            &bed,
            b"# targets\nchr1\t10\t20\tfirst\nchr1\t18\t30\nchr2\t5\t6\n",
        )
        .unwrap();
        let selection = RegionSelection {
            intervals: vec![GenomicInterval {
                contig: String::from("chr1"),
                start: 25,
                end: 40,
            }],
            regions_file: Some(bed),
        };
        let regions = plan_call_regions(&references(), 64, &selection).unwrap();
        assert_eq!(
            regions
                .iter()
                .map(|region| (region.reference, region.start, region.end))
                .collect::<Vec<_>>(),
            vec![(0, 10, 40), (1, 5, 6)]
        );
        fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn bed_codec_and_text_errors_retain_path_and_line_context() {
        let directory = unique_directory("errors");
        fs::create_dir(&directory).unwrap();
        let bed = directory.join("targets.bed");
        let selection = RegionSelection {
            intervals: Vec::new(),
            regions_file: Some(bed.clone()),
        };

        fs::write(&bed, b"# targets\nchr1\tnot-a-number\t20\n").unwrap();
        let syntax_error = plan_call_regions(&references(), 64, &selection).unwrap_err();
        assert_eq!(syntax_error.kind(), CallErrorKind::Input);
        let message = syntax_error.to_string();
        assert!(message.contains(&bed.display().to_string()));
        assert!(message.contains("line 2"));
        assert!(message.contains("invalid BED3+ syntax"));
        assert!(message.contains("BED column 2"));

        fs::write(&bed, b"# targets\nchr\xff\t10\t20\n").unwrap();
        let text_error = plan_call_regions(&references(), 64, &selection).unwrap_err();
        assert_eq!(text_error.kind(), CallErrorKind::Input);
        let message = text_error.to_string();
        assert!(message.contains(&bed.display().to_string()));
        assert!(message.contains("line 2"));
        assert!(message.contains("non-UTF-8 contig"));

        fs::remove_dir_all(directory).unwrap();
    }
}
