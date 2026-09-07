//! Process-level ground truth for the thin CLI.

use std::ffi::{OsStr, OsString};
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

use bsbit_align::extension::VerifiedAlignment;
use bsbit_align::materialize::traceback_read_placement;
use bsbit_core::bisulfite::{AlignmentOrientation, BisulfiteStrand};
use bsbit_core::coordinate::{ReferenceInterval, ReferenceLength};
use bsbit_core::reference::ReferenceSequenceMd5;
use bsbit_core::sequence::normalize_dna;
use bsbit_hts::{
    AlignmentRecord, AlignmentRecordLimits, BamStagingWriter, DecodedReader, MappedAlignmentRecord,
    RecordMappingQuality, RecordReference, RecordSegment, SamHeader, SamHeaderReference,
    SamSortOrder, build_bam_index_create_new,
};
use bsbit_index::reference::{ContigInput, ReferenceBuildLimits, ReferenceIndex};

const BAM_CIGAR_CODES: &[u8; 10] = b"MIDNSHP=XB";
const BAM_BASES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

fn unique_directory(label: &str) -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock is after epoch")
        .as_nanos();
    std::env::temp_dir().join(format!("bsbit-cli-{label}-{}-{nonce}", std::process::id()))
}

fn decoded_text(path: &Path) -> String {
    let mut reader = DecodedReader::open(path).expect("encoded text output opens");
    let mut decoded = String::new();
    reader
        .read_to_string(&mut decoded)
        .expect("encoded text output decodes");
    decoded
}

fn indexed_call_fixture(directory: &Path) -> (PathBuf, PathBuf) {
    const REFERENCE: &[u8] = b"ACGTTGCACTGATCGATGCTAGCTACGATCGTTCGAGTACCTGACGTA";
    let mut observed = REFERENCE.to_vec();
    observed[0] = b'G';
    let reference = ReferenceIndex::build(
        vec![ContigInput::new(
            b"chr1".to_vec(),
            normalize_dna(REFERENCE).expect("fixture reference is canonical"),
        )],
        ReferenceBuildLimits::MAX,
    )
    .expect("fixture reference builds");
    let read = normalize_dna(&observed).expect("fixture read is canonical");
    let contig_id = reference.contig_id(0).expect("fixture contig id");
    let interval = ReferenceInterval::new(
        0,
        read.len(),
        ReferenceLength::new(
            u64::try_from(REFERENCE.len()).expect("fixture reference length fits u64"),
        ),
    )
    .expect("fixture interval is bounded");
    let mapping = traceback_read_placement(
        &reference,
        &read,
        &contig_id,
        interval,
        BisulfiteStrand::OT,
        1,
    )
    .expect("fixture read placement materializes");
    let header = fixture_sam_header(&reference)
        .with_bsbit_metadata(
            bsbit_hts::BsbitHeaderMetadata::new(
                bsbit_hts::BsbitAlignmentMode::DirectionalPairedEnd,
            ),
            AlignmentRecordLimits::default(),
        )
        .expect("fixture metadata fits")
        .with_sort_order(SamSortOrder::Coordinate);
    let staging = directory.join("fixture.bam.tmp");
    let input = directory.join("fixture.bam");
    let mut writer =
        BamStagingWriter::create_new(&staging, &header, AlignmentRecordLimits::default())
            .expect("fixture BAM opens");
    let qualities = vec![b'I'; observed.len()];
    for ordinal in 0..24 {
        let query_name = format!("read-{ordinal:02}");
        let record = fixture_alignment_record(
            &reference,
            query_name.as_bytes(),
            &read,
            &qualities,
            &mapping,
        );
        writer
            .write_record_as_bam(&record)
            .expect("fixture BAM record writes");
    }
    writer
        .finish()
        .expect("fixture BAM finishes")
        .publish_create_new(&input)
        .expect("fixture BAM publishes");
    build_bam_index_create_new(&input, input.with_extension("bam.bai"), 1)
        .expect("fixture BAI builds");
    let fasta = directory.join("reference.fa");
    let mut fasta_contents = b">chr1\n".to_vec();
    fasta_contents.extend_from_slice(REFERENCE);
    fasta_contents.push(b'\n');
    fs::write(&fasta, fasta_contents).expect("fixture FASTA writes");
    fs::write(
        fasta.with_extension("fa.fai"),
        format!(
            "chr1\t{}\t6\t{}\t{}\n",
            REFERENCE.len(),
            REFERENCE.len(),
            REFERENCE.len() + 1
        ),
    )
    .expect("fixture FAI writes");
    (input, fasta)
}

fn fixture_sam_header(reference: &ReferenceIndex) -> SamHeader {
    let mut entries = Vec::new();
    for ordinal in 0..reference.contig_count() {
        let id = reference.contig_id(ordinal).expect("fixture contig id");
        let contig = reference.resolve_contig(&id).expect("fixture contig");
        entries.push(
            SamHeaderReference::new(ordinal, contig.name(), contig.sequence().len())
                .expect("fixture header entry")
                .with_md5(ReferenceSequenceMd5::from_normalized(
                    contig.sequence().bases(),
                )),
        );
    }
    SamHeader::new(entries, AlignmentRecordLimits::default()).expect("fixture header builds")
}

fn fixture_alignment_record(
    reference: &ReferenceIndex,
    query_name: &[u8],
    read: &bsbit_core::sequence::NormalizedSequence,
    quality: &[u8],
    alignment: &VerifiedAlignment,
) -> AlignmentRecord {
    let contig = reference
        .resolve_contig(alignment.contig())
        .expect("fixture alignment contig");
    let record_reference = RecordReference::new(
        contig.ordinal(),
        contig.name(),
        contig.sequence().len(),
        alignment.interval(),
    )
    .expect("fixture record reference");
    let literal_nm = alignment
        .cached_literal_nm()
        .unwrap_or_else(|| alignment.distance().get());
    let mapping = MappedAlignmentRecord::new(
        record_reference,
        alignment.orientation(),
        alignment.strand(),
        alignment.cytosine_strand(),
        alignment.cigar().clone(),
        read.len(),
        u32::try_from(literal_nm).expect("fixture NM fits u32"),
        None,
        None,
        AlignmentRecordLimits::default(),
    )
    .expect("fixture mapped record");
    let mut sequence = read
        .bases()
        .iter()
        .map(|base| base.as_ascii())
        .collect::<Vec<_>>();
    let mut quality = quality.to_vec();
    if matches!(alignment.orientation(), AlignmentOrientation::Reverse) {
        sequence = read
            .bases()
            .iter()
            .rev()
            .map(|base| base.complement().as_ascii())
            .collect();
        quality.reverse();
    }
    AlignmentRecord::new(
        query_name,
        RecordSegment::Unpaired,
        false,
        RecordMappingQuality::Calibrated(60),
        Some(mapping),
        None,
        0,
        &sequence,
        Some(&quality),
        AlignmentRecordLimits::default(),
    )
    .expect("fixture BAM record builds")
}

fn run(arguments: impl IntoIterator<Item = impl AsRef<OsStr>>) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_bsbit"));
    command.args(arguments);
    command.output().expect("bsbit process starts")
}

fn index(reference: &Path, output: &Path) -> Output {
    run([
        OsString::from("index"),
        OsString::from("--reference"),
        reference.as_os_str().to_owned(),
        OsString::from("--output"),
        output.as_os_str().to_owned(),
    ])
}

fn internal_index_prefix(index: &Path) -> PathBuf {
    const FNV_OFFSET_BASIS: u64 = 0xcbf2_9ce4_8422_2325;
    const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;
    let mut hash = FNV_OFFSET_BASIS;
    for byte in index
        .file_name()
        .expect("index filename")
        .as_encoded_bytes()
    {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(FNV_PRIME);
    }
    index.with_file_name(format!(".bsbit-index-{hash:016x}"))
}

fn align(snapshot: &Path, read1: &Path, read2: Option<&Path>, output_bam: &Path) -> Output {
    let mut arguments = vec![
        OsString::from("align"),
        OsString::from("--index"),
        snapshot.as_os_str().to_owned(),
        OsString::from("-1"),
        read1.as_os_str().to_owned(),
    ];
    if let Some(path) = read2 {
        arguments.push(OsString::from("-2"));
        arguments.push(path.as_os_str().to_owned());
    }
    arguments.extend([
        OsString::from("--output"),
        output_bam.as_os_str().to_owned(),
    ]);
    run(arguments)
}

fn align_single_sensitive(snapshot: &Path, read1: &Path, output_bam: &Path) -> Output {
    run([
        OsString::from("align"),
        OsString::from("--index"),
        snapshot.as_os_str().to_owned(),
        OsString::from("--read1"),
        read1.as_os_str().to_owned(),
        OsString::from("--output"),
        output_bam.as_os_str().to_owned(),
        OsString::from("--total-threads"),
        OsString::from("2"),
        OsString::from("--sensitive"),
    ])
}

fn align_single_nondirectional(snapshot: &Path, read1: &Path, output_bam: &Path) -> Output {
    run([
        OsString::from("align"),
        OsString::from("--index"),
        snapshot.as_os_str().to_owned(),
        OsString::from("--read1"),
        read1.as_os_str().to_owned(),
        OsString::from("--output"),
        output_bam.as_os_str().to_owned(),
        OsString::from("--non-directional"),
    ])
}

fn assert_success(output: &Output) {
    assert!(
        output.status.success(),
        "status={:?}\nstdout={}\nstderr={}",
        output.status.code(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stderr.is_empty());
}

fn assert_human_progress(output: &Output, command: &str) {
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        stdout.contains(&format!(
            "[bsbit::{command}] start architecture={}",
            std::env::consts::ARCH
        )),
        "missing architecture log from {stdout:?}"
    );
    assert!(
        stdout.contains(" backend="),
        "missing backend from {stdout:?}"
    );
    assert!(
        stdout.contains(" instruction_set="),
        "missing instruction set from {stdout:?}"
    );
    assert!(
        stdout.contains(&format!("[bsbit::{command}] completed ")),
        "missing completion log from {stdout:?}"
    );
}

fn assert_metrics_only(output: Output) {
    assert_success(&output);
    let stdout = String::from_utf8(output.stdout).expect("metrics TSV is UTF-8");
    assert!(!stdout.contains("[bsbit::"));
    assert_eq!(stdout.lines().count(), 2);
    assert!(stdout.starts_with("schema\tpairs\t"));
    let mut lines = stdout.lines();
    let header = lines.next().unwrap().split('\t').collect::<Vec<_>>();
    let values = lines.next().unwrap().split('\t').collect::<Vec<_>>();
    assert_eq!(header.len(), values.len());
    let field = |name| {
        let index = header
            .iter()
            .position(|candidate| *candidate == name)
            .unwrap_or_else(|| panic!("missing metrics field {name:?}"));
        values[index]
    };
    assert_eq!(field("architecture"), std::env::consts::ARCH);
    assert!(!field("backend").is_empty());
    assert!(!field("instruction_set").is_empty());
    assert_eq!(field("alignment_policy"), bsbit_align::ALIGNMENT_POLICY_ID);
    assert_eq!(field("mapq_policy"), bsbit_align::MAPQ_POLICY_ID);
    assert!(
        values
            .first()
            .is_some_and(|value| *value == "bsbit-alignment-metrics-paired-end-v1")
    );
}

fn assert_single_metrics_only(output: Output) {
    assert_success(&output);
    let stdout = String::from_utf8(output.stdout).expect("single metrics TSV is UTF-8");
    assert!(!stdout.contains("[bsbit::"));
    assert_eq!(stdout.lines().count(), 2);
    assert!(stdout.starts_with("schema\treads\t"));
    let mut lines = stdout.lines();
    let header = lines.next().unwrap().split('\t').collect::<Vec<_>>();
    let values = lines.next().unwrap().split('\t').collect::<Vec<_>>();
    assert_eq!(header.len(), values.len());
    let field = |name| {
        let index = header
            .iter()
            .position(|candidate| *candidate == name)
            .unwrap_or_else(|| panic!("missing single metrics field {name:?}"));
        values[index]
    };
    assert_eq!(field("schema"), "bsbit-alignment-metrics-single-end-v1");
    assert_eq!(field("max_edit_distance"), "5");
    assert_eq!(field("alignment_policy"), bsbit_align::ALIGNMENT_POLICY_ID);
    assert_eq!(field("mapq_policy"), bsbit_align::MAPQ_POLICY_ID);
    assert_eq!(field("architecture"), std::env::consts::ARCH);
    assert!(!field("backend").is_empty());
    assert!(!field("instruction_set").is_empty());
}

fn take<'a>(bytes: &'a [u8], offset: &mut usize, length: usize) -> &'a [u8] {
    let end = offset
        .checked_add(length)
        .expect("BAM offset does not wrap");
    let value = bytes.get(*offset..end).expect("complete BAM field");
    *offset = end;
    value
}

fn bam_u16(bytes: &[u8], offset: &mut usize) -> u16 {
    u16::from_le_bytes(take(bytes, offset, 2).try_into().expect("two bytes"))
}

fn bam_u32(bytes: &[u8], offset: &mut usize) -> u32 {
    u32::from_le_bytes(take(bytes, offset, 4).try_into().expect("four bytes"))
}

fn bam_i32(bytes: &[u8], offset: &mut usize) -> i32 {
    i32::from_le_bytes(take(bytes, offset, 4).try_into().expect("four bytes"))
}

fn decode_process_bam(path: &Path) -> Vec<Vec<Vec<u8>>> {
    let mut reader = DecodedReader::open(path).expect("process BAM opens");
    let mut bytes = Vec::new();
    reader.read_to_end(&mut bytes).expect("process BAM decodes");
    reader.close().expect("process BAM closes");
    let mut offset = 0;
    assert_eq!(take(&bytes, &mut offset, 4), b"BAM\x01");
    let header_length = usize::try_from(bam_i32(&bytes, &mut offset)).expect("header length");
    take(&bytes, &mut offset, header_length);
    let reference_count = usize::try_from(bam_i32(&bytes, &mut offset)).expect("reference count");
    let mut references = Vec::new();
    for _ in 0..reference_count {
        let name_length = usize::try_from(bam_i32(&bytes, &mut offset)).expect("name length");
        let name = take(&bytes, &mut offset, name_length);
        assert_eq!(name.last(), Some(&0));
        references.push(name[..name.len() - 1].to_vec());
        assert!(bam_i32(&bytes, &mut offset) >= 0);
    }

    let mut records = Vec::new();
    while offset != bytes.len() {
        let block_length = usize::try_from(bam_i32(&bytes, &mut offset)).expect("block length");
        let block = take(&bytes, &mut offset, block_length);
        records.push(decode_bam_record(block, &references));
    }
    records
}

fn decode_process_bam_header(path: &Path) -> Vec<u8> {
    let mut reader = DecodedReader::open(path).expect("process BAM opens");
    let mut bytes = Vec::new();
    reader.read_to_end(&mut bytes).expect("process BAM decodes");
    reader.close().expect("process BAM closes");
    let mut offset = 0;
    assert_eq!(take(&bytes, &mut offset, 4), b"BAM\x01");
    let header_length = usize::try_from(bam_i32(&bytes, &mut offset)).expect("header length");
    take(&bytes, &mut offset, header_length).to_vec()
}

fn decode_bam_record(bytes: &[u8], references: &[Vec<u8>]) -> Vec<Vec<u8>> {
    let mut offset = 0;
    let reference_id = bam_i32(bytes, &mut offset);
    let position = bam_i32(bytes, &mut offset);
    let read_name_length = usize::from(take(bytes, &mut offset, 1)[0]);
    let mapping_quality = take(bytes, &mut offset, 1)[0];
    let _bin = bam_u16(bytes, &mut offset);
    let cigar_count = usize::from(bam_u16(bytes, &mut offset));
    let flag = bam_u16(bytes, &mut offset);
    let sequence_length = usize::try_from(bam_i32(bytes, &mut offset)).expect("sequence length");
    let mate_reference_id = bam_i32(bytes, &mut offset);
    let mate_position = bam_i32(bytes, &mut offset);
    let template_length = bam_i32(bytes, &mut offset);
    let raw_name = take(bytes, &mut offset, read_name_length);
    assert_eq!(raw_name.last(), Some(&0));
    let name = raw_name[..raw_name.len() - 1].to_vec();

    let mut cigar = Vec::new();
    for _ in 0..cigar_count {
        let encoded = bam_u32(bytes, &mut offset);
        cigar.extend_from_slice((encoded >> 4).to_string().as_bytes());
        cigar.push(BAM_CIGAR_CODES[usize::try_from(encoded & 0xf).expect("CIGAR code")]);
    }
    if cigar.is_empty() {
        cigar.push(b'*');
    }

    let packed = take(bytes, &mut offset, sequence_length.div_ceil(2));
    let sequence = (0..sequence_length)
        .map(|index| {
            let value = packed[index / 2];
            let code = if index % 2 == 0 {
                value >> 4
            } else {
                value & 0xf
            };
            BAM_BASES[usize::from(code)]
        })
        .collect::<Vec<_>>();
    let raw_quality = take(bytes, &mut offset, sequence_length);
    let quality = if raw_quality.iter().all(|value| *value == u8::MAX) {
        b"*".to_vec()
    } else {
        raw_quality.iter().map(|value| value + 33).collect()
    };

    let mut fields = vec![
        name,
        flag.to_string().into_bytes(),
        reference_name(reference_id, references),
        if position < 0 {
            b"0".to_vec()
        } else {
            (position + 1).to_string().into_bytes()
        },
        mapping_quality.to_string().into_bytes(),
        cigar,
        mate_reference_name(reference_id, mate_reference_id, references),
        if mate_position < 0 {
            b"0".to_vec()
        } else {
            (mate_position + 1).to_string().into_bytes()
        },
        template_length.to_string().into_bytes(),
        if sequence.is_empty() {
            b"*".to_vec()
        } else {
            sequence
        },
        quality,
    ];
    while offset != bytes.len() {
        fields.push(decode_bam_aux(bytes, &mut offset));
    }
    fields
}

fn reference_name(reference_id: i32, references: &[Vec<u8>]) -> Vec<u8> {
    usize::try_from(reference_id)
        .ok()
        .and_then(|ordinal| references.get(ordinal))
        .cloned()
        .unwrap_or_else(|| b"*".to_vec())
}

fn mate_reference_name(
    reference_id: i32,
    mate_reference_id: i32,
    references: &[Vec<u8>],
) -> Vec<u8> {
    if mate_reference_id >= 0 && mate_reference_id == reference_id {
        b"=".to_vec()
    } else {
        reference_name(mate_reference_id, references)
    }
}

fn decode_bam_aux(bytes: &[u8], offset: &mut usize) -> Vec<u8> {
    let tag = take(bytes, offset, 2);
    let physical_type = take(bytes, offset, 1)[0];
    let (logical_type, value) = match physical_type {
        b'c' => (
            b'i',
            i8::from_le_bytes([take(bytes, offset, 1)[0]])
                .to_string()
                .into_bytes(),
        ),
        b'C' => (b'i', take(bytes, offset, 1)[0].to_string().into_bytes()),
        b's' => (
            b'i',
            i16::from_le_bytes(take(bytes, offset, 2).try_into().expect("i16"))
                .to_string()
                .into_bytes(),
        ),
        b'S' => (b'i', bam_u16(bytes, offset).to_string().into_bytes()),
        b'i' => (b'i', bam_i32(bytes, offset).to_string().into_bytes()),
        b'I' => (b'i', bam_u32(bytes, offset).to_string().into_bytes()),
        b'Z' => {
            let end = bytes[*offset..]
                .iter()
                .position(|byte| *byte == 0)
                .expect("terminated string auxiliary");
            let value = take(bytes, offset, end).to_vec();
            take(bytes, offset, 1);
            (b'Z', value)
        }
        _ => panic!("unsupported test auxiliary type {physical_type}"),
    };
    let mut field = tag.to_vec();
    field.extend_from_slice(&[b':', logical_type, b':']);
    field.extend_from_slice(&value);
    field
}

#[test]
fn general_help_version_and_usage_errors_are_golden() {
    let help = run(["--help"]);
    assert_eq!(help.status.code(), Some(0));
    assert_eq!(help.stdout, bsbit_cli::GENERAL_HELP.as_bytes());
    assert!(help.stderr.is_empty());

    for argument in ["--version", "-v"] {
        let version = run([argument]);
        assert_eq!(version.status.code(), Some(0));
        assert_eq!(
            version.stdout,
            concat!("bsbit ", env!("CARGO_PKG_VERSION"), "\n").as_bytes()
        );
        assert!(version.stderr.is_empty());
    }

    let uppercase_version = run(["-V"]);
    assert_eq!(uppercase_version.status.code(), Some(2));
    assert!(uppercase_version.stdout.is_empty());
    assert_eq!(
        uppercase_version.stderr,
        b"bsbit: unknown command `-V`; run `bsbit --help`\n"
    );

    let missing_call_module = run(["call"]);
    assert_eq!(missing_call_module.status.code(), Some(2));
    assert_eq!(
        missing_call_module.stderr,
        b"bsbit: missing call module; run `bsbit call --help`\n"
    );

    let missing = run(std::iter::empty::<&str>());
    assert_eq!(missing.status.code(), Some(2));
    assert!(missing.stdout.is_empty());
    assert_eq!(
        missing.stderr,
        b"bsbit: missing command; run `bsbit --help`\n"
    );

    {
        let unknown = run(["align", "--pbat", "yes"]);
        assert_eq!(unknown.status.code(), Some(2));
        assert_eq!(unknown.stderr, b"bsbit: unknown option --pbat\n");
    }
}

#[test]
fn cpu_command_reports_independent_features_and_validates_backend_names() {
    let automatic = run(["cpu"]);
    assert_eq!(automatic.status.code(), Some(0));
    assert!(automatic.stderr.is_empty());
    let report = String::from_utf8(automatic.stdout).expect("CPU report is UTF-8");
    let value = |key: &str| {
        report
            .lines()
            .find_map(|line| line.strip_prefix(key))
            .unwrap_or_else(|| panic!("missing {key:?} from {report:?}"))
    };
    assert_eq!(value("architecture="), std::env::consts::ARCH);
    for key in [
        "sse2=",
        "sse4.1=",
        "sse4.2=",
        "popcnt=",
        "avx2=",
        "avx512f=",
        "avx512bw=",
        "neon=",
    ] {
        assert!(
            matches!(value(key), "0" | "1"),
            "invalid {key} in {report:?}"
        );
    }
    #[cfg(target_arch = "x86_64")]
    {
        assert_eq!(value("neon="), "0");
    }
    #[cfg(target_arch = "aarch64")]
    {
        assert_eq!(value("sse2="), "0");
        assert_eq!(value("sse4.1="), "0");
        assert_eq!(value("sse4.2="), "0");
        assert_eq!(value("popcnt="), "0");
        assert_eq!(value("avx2="), "0");
        assert_eq!(value("avx512f="), "0");
        assert_eq!(value("avx512bw="), "0");
        assert_eq!(value("neon="), "1");
    }
    assert_eq!(value("backend="), expected_backend_from_report(&report));

    let forced_name = "scalar";
    let forced = run(["cpu", "--simd-backend", forced_name]);
    assert_eq!(forced.status.code(), Some(0));
    assert!(forced.stderr.is_empty());
    let forced_report = String::from_utf8_lossy(&forced.stdout);
    assert!(
        forced_report
            .lines()
            .any(|line| line == format!("backend={forced_name}")),
        "unexpected forced report: {forced_report:?}"
    );
    assert!(
        forced_report
            .lines()
            .any(|line| line.starts_with("instruction_set="))
    );

    let wrong_architecture = if cfg!(target_arch = "aarch64") {
        "avx2"
    } else {
        "neon"
    };
    let unsupported = run(["cpu", "--simd-backend", wrong_architecture]);
    assert_eq!(unsupported.status.code(), Some(1));
    assert!(unsupported.stdout.is_empty());
    let unsupported_error = String::from_utf8(unsupported.stderr).expect("CPU error is UTF-8");
    assert!(
        unsupported_error.contains("requires")
            && unsupported_error.contains("current architecture is"),
        "unexpected architecture error: {unsupported_error:?}"
    );

    let invalid = run(["cpu", "--simd-backend", "native"]);
    assert_eq!(invalid.status.code(), Some(2));
    assert!(invalid.stdout.is_empty());
    assert_eq!(
        invalid.stderr,
        b"bsbit: invalid --simd-backend `native`: expected auto, scalar, sse2, sse4.2, avx2, avx512, or neon\n"
    );
}

fn expected_backend_from_report(report: &str) -> &'static str {
    let enabled = |key: &str| report.lines().any(|line| line == format!("{key}=1"));
    #[cfg(target_arch = "x86_64")]
    {
        if enabled("avx512f") && enabled("avx512bw") && enabled("avx2") && enabled("popcnt") {
            "avx512"
        } else if enabled("avx2") && enabled("popcnt") {
            "avx2"
        } else if enabled("sse4.1") && enabled("sse4.2") && enabled("popcnt") {
            "sse4.2"
        } else if enabled("sse2") {
            "sse2"
        } else {
            "scalar"
        }
    }
    #[cfg(target_arch = "aarch64")]
    {
        if enabled("neon") { "neon" } else { "scalar" }
    }
    #[cfg(not(any(target_arch = "x86_64", target_arch = "aarch64")))]
    {
        let _ = enabled;
        "scalar"
    }
}

#[test]
fn subcommand_help_and_remote_path_error_are_golden() {
    let remote = run([
        "index",
        "--reference",
        "https://example.invalid/ref.fa",
        "--output",
        "out",
    ]);
    assert_eq!(remote.status.code(), Some(2));
    assert_eq!(
        remote.stderr,
        b"bsbit: unsupported non-local path `https://example.invalid/ref.fa` for `--reference`\n"
    );

    let index_help = run(["index", "--help"]);
    assert_eq!(index_help.status.code(), Some(0));
    assert_eq!(index_help.stdout, bsbit_cli::INDEX_HELP.as_bytes());
    let align_help = run(["align", "--help"]);
    assert_eq!(align_help.status.code(), Some(0));
    assert_eq!(align_help.stdout, bsbit_cli::ALIGN_HELP.as_bytes());

    let combine_help = run(["combine", "--help"]);
    assert_eq!(combine_help.status.code(), Some(0));
    assert_eq!(combine_help.stdout, bsbit_cli::COMBINE_HELP.as_bytes());
}

#[test]
fn option_usage_errors_are_golden() {
    let missing = run(["index", "--output", "out"]);
    assert_eq!(missing.status.code(), Some(2));
    assert_eq!(
        missing.stderr,
        b"bsbit: missing required option `--reference`\n"
    );
    let duplicate = run([
        "index",
        "--reference",
        "one",
        "--reference",
        "two",
        "--output",
        "out",
    ]);
    assert_eq!(duplicate.status.code(), Some(2));
    assert_eq!(duplicate.stderr, b"bsbit: duplicate option `--reference`\n");
    let unknown = run(["align", "--unknown-option"]);
    assert_eq!(unknown.status.code(), Some(2));
    assert_eq!(unknown.stderr, b"bsbit: unknown option --unknown-option\n");
}

#[test]
fn call_help_is_available_only_through_the_umbrella_command() {
    let call_help = run(["call", "--help"]);
    assert_eq!(call_help.status.code(), Some(0));
    assert_eq!(call_help.stdout, bsbit_cli::CALL_HELP.as_bytes());
    assert!(call_help.stderr.is_empty());

    let meth_help = run(["call", "meth", "--help"]);
    assert_eq!(meth_help.status.code(), Some(0));
    assert_eq!(meth_help.stdout, bsbit_cli::CALL_METH_HELP.as_bytes());

    let snp_help = run(["call", "snp", "--help"]);
    assert_eq!(snp_help.status.code(), Some(0));
    assert_eq!(snp_help.stdout, bsbit_cli::CALL_SNP_HELP.as_bytes());

    let joint_help = run(["call", "joint", "--help"]);
    assert_eq!(joint_help.status.code(), Some(0));
    assert_eq!(joint_help.stdout, bsbit_cli::CALL_JOINT_HELP.as_bytes());
}

#[test]
fn umbrella_call_validates_inputs_and_joint_prefix() {
    let input = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../external/htslib/test/range.bam");
    let output_directory = unique_directory("call-validation");
    fs::create_dir(&output_directory).expect("output directory");
    let output = output_directory.join("calls.cgmap.gz");
    let meth = run([
        OsString::from("call"),
        OsString::from("meth"),
        OsString::from("-i"),
        input.as_os_str().to_owned(),
        OsString::from("--reference"),
        OsString::from("missing.fa"),
        OsString::from("-o"),
        output.as_os_str().to_owned(),
        OsString::from("-f"),
        OsString::from("cgmap"),
    ]);
    assert_eq!(meth.status.code(), Some(1));
    assert!(meth.stdout.is_empty());
    assert!(String::from_utf8_lossy(&meth.stderr).contains("missing.fa"));
    assert_eq!(
        fs::metadata(&output).expect("direct output exists").len(),
        0
    );

    let prefix_directory = unique_directory("joint-prefix");
    fs::create_dir(&prefix_directory).expect("prefix directory");
    let prefix = prefix_directory.join("calls");
    let meth_output = prefix.with_extension("CGmap");
    let vcf_output = prefix.with_extension("vcf");
    let joint = run([
        OsString::from("call"),
        OsString::from("joint"),
        OsString::from("-i"),
        input.as_os_str().to_owned(),
        OsString::from("--reference"),
        OsString::from("missing.fa"),
        OsString::from("--prefix"),
        prefix.as_os_str().to_owned(),
        OsString::from("--meth-format"),
        OsString::from("cgmap"),
    ]);
    assert_eq!(joint.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&joint.stderr).contains("missing.fa"));
    assert_eq!(
        fs::metadata(&meth_output)
            .expect("direct methylation output exists")
            .len(),
        0
    );
    assert_eq!(
        fs::metadata(&vcf_output)
            .expect("direct VCF output exists")
            .len(),
        0
    );
    fs::remove_dir_all(output_directory).expect("output cleanup");
    fs::remove_dir_all(prefix_directory).expect("prefix cleanup");
}

#[test]
fn combine_command_builds_a_filtered_named_methylation_matrix() {
    let directory = unique_directory("combine-command");
    fs::create_dir(&directory).expect("fresh directory");
    let first = directory.join("first.bed");
    let second = directory.join("second.bed");
    let matrix = directory.join("matrix.bed");
    fs::write(
        &first,
        concat!(
            "chr1\t0\t1\tm,CG,0\t10\t+\t0\t1\t255,0,0\t10\t70.00\t7\t3\t0\t0\t0\t0\t0\n",
            "chr1\t2\t3\tm,CG,0\t2\t+\t2\t3\t255,0,0\t2\t50.00\t1\t1\t0\t0\t0\t0\t0\n",
        ),
    )
    .expect("first bedMethyl");
    fs::write(
        &second,
        concat!(
            "chr1\t0\t1\tm,CG,0\t5\t+\t0\t1\t255,0,0\t5\t20.00\t1\t4\t0\t0\t0\t0\t0\n",
            "chr1\t1\t2\tm,CHH,0\t6\t+\t1\t2\t255,0,0\t6\t0.00\t0\t6\t0\t0\t0\t0\t0\n",
        ),
    )
    .expect("second bedMethyl");
    let result = run([
        OsString::from("combine"),
        OsString::from("--input"),
        OsString::from(format!("{},{}", first.display(), second.display())),
        OsString::from("--sample-name"),
        OsString::from("case,control"),
        OsString::from("--output"),
        matrix.as_os_str().to_owned(),
        OsString::from("--matrix"),
        OsString::from("count"),
        OsString::from("--min-count"),
        OsString::from("5"),
        OsString::from("--min-prop"),
        OsString::from("1"),
        OsString::from("--threads"),
        OsString::from("2"),
    ]);
    assert_eq!(result.status.code(), Some(0), "{:?}", result.stderr);
    assert!(result.stdout.is_empty());
    assert!(result.stderr.is_empty());
    assert_eq!(
        decoded_text(&matrix),
        concat!(
            "##bsbit_matrix_format=count\n",
            "##bsbit_min_count=5\n",
            "##bsbit_min_prop=1.000000000\n",
            "##bsbit_cg_only=false\n",
            "#chrom\tstart\tend\tmodification\tscore\tstrand",
            "\tcase_meth_count\tcase_total_count",
            "\tcontrol_meth_count\tcontrol_total_count\n",
            "chr1\t0\t1\tm,CG,0\t0\t+\t7\t10\t1\t5\n",
        )
    );
    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[test]
fn combine_command_uses_input_path_when_sample_names_are_omitted() {
    let directory = unique_directory("combine-path-name");
    fs::create_dir(&directory).expect("fresh directory");
    let input = directory.join("sample one.bed");
    let matrix = directory.join("matrix.bed");
    fs::write(
        &input,
        "chr1\t0\t1\tm,CG,0\t4\t+\t0\t1\t255,0,0\t4\t75.00\t3\t1\t0\t0\t0\t0\t0\n",
    )
    .expect("input bedMethyl");
    let result = run([
        OsString::from("combine"),
        OsString::from("--input"),
        input.as_os_str().to_owned(),
        OsString::from("--output"),
        matrix.as_os_str().to_owned(),
    ]);
    assert_eq!(result.status.code(), Some(0), "{:?}", result.stderr);
    assert!(result.stdout.is_empty());
    assert!(result.stderr.is_empty());
    assert_eq!(
        decoded_text(&matrix),
        format!(
            concat!(
                "##bsbit_matrix_format=level\n",
                "##bsbit_min_count=1\n",
                "##bsbit_min_prop=0.000000000\n",
                "##bsbit_cg_only=false\n",
                "#chrom\tstart\tend\tmodification\tscore\tstrand\t{}\n",
                "chr1\t0\t1\tm,CG,0\t0\t+\t0.750000\n",
            ),
            input.display()
        )
    );
    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[test]
fn combine_requires_one_comma_separated_input_and_sample_name_list() {
    let repeated_input = run([
        "combine",
        "--input",
        "one.bed",
        "--input",
        "two.bed",
        "--output",
        "matrix.bed",
    ]);
    assert_eq!(repeated_input.status.code(), Some(2));
    assert_eq!(
        repeated_input.stderr,
        b"bsbit: duplicate option `--input`; provide one comma-separated list\n"
    );

    let repeated_names = run([
        "combine",
        "--input",
        "one.bed,two.bed",
        "--sample-name",
        "one",
        "--sample-name",
        "two",
        "--output",
        "matrix.bed",
    ]);
    assert_eq!(repeated_names.status.code(), Some(2));
    assert_eq!(
        repeated_names.stderr,
        b"bsbit: duplicate option `--sample-name`; provide one comma-separated list\n"
    );
}

fn assert_command_succeeded(output: &Output) {
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn run_meth_call(input: &Path, reference: &Path, output: &Path) {
    let result = run([
        OsString::from("call"),
        OsString::from("meth"),
        OsString::from("-i"),
        input.as_os_str().to_owned(),
        OsString::from("--reference"),
        reference.as_os_str().to_owned(),
        OsString::from("-o"),
        output.as_os_str().to_owned(),
        OsString::from("-f"),
        OsString::from("cgmap"),
        OsString::from("-c"),
        OsString::from("true"),
        OsString::from("-t"),
        OsString::from("2"),
    ]);
    assert_command_succeeded(&result);
}

fn run_snp_call(input: &Path, reference: &Path, output: &Path) {
    let result = run([
        OsString::from("call"),
        OsString::from("snp"),
        OsString::from("-i"),
        input.as_os_str().to_owned(),
        OsString::from("--reference"),
        reference.as_os_str().to_owned(),
        OsString::from("-o"),
        output.as_os_str().to_owned(),
        OsString::from("-c"),
        OsString::from("true"),
        OsString::from("-t"),
        OsString::from("2"),
    ]);
    assert_command_succeeded(&result);
}

fn run_joint_call(input: &Path, reference: &Path, prefix: &Path) {
    let result = run([
        OsString::from("call"),
        OsString::from("joint"),
        OsString::from("-i"),
        input.as_os_str().to_owned(),
        OsString::from("--reference"),
        reference.as_os_str().to_owned(),
        OsString::from("--prefix"),
        prefix.as_os_str().to_owned(),
        OsString::from("--meth-format"),
        OsString::from("cgmap"),
        OsString::from("-c"),
        OsString::from("true"),
        OsString::from("-t"),
        OsString::from("2"),
    ]);
    assert_command_succeeded(&result);
}

#[test]
fn unindexed_uncompressed_fasta_warns_and_remains_usable_by_the_cli() {
    let directory = unique_directory("unindexed-reference-warning");
    fs::create_dir(&directory).expect("call fixture directory");
    let (input, reference) = indexed_call_fixture(&directory);
    fs::remove_file(reference.with_extension("fa.fai")).expect("remove fixture FAI");
    let output = directory.join("meth.cgmap");

    let result = run([
        OsString::from("call"),
        OsString::from("meth"),
        OsString::from("-i"),
        input.as_os_str().to_owned(),
        OsString::from("--reference"),
        reference.as_os_str().to_owned(),
        OsString::from("-o"),
        output.as_os_str().to_owned(),
        OsString::from("-f"),
        OsString::from("cgmap"),
    ]);

    assert_command_succeeded(&result);
    let stderr = String::from_utf8(result.stderr).expect("warning is UTF-8");
    assert!(stderr.contains("bsbit: warning:"));
    assert!(stderr.contains("no adjacent .fai index"));
    assert!(stderr.contains("temporary in-memory line-layout index"));
    assert!(output.exists());
    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[test]
fn umbrella_call_modules_are_consistent_on_real_bam() {
    let directory = unique_directory("call-e2e");
    fs::create_dir(&directory).expect("call fixture directory");
    let (input, reference) = indexed_call_fixture(&directory);
    let meth_output = directory.join("meth.cgmap.gz");
    let snp_output = directory.join("snp.vcf.gz");
    let joint_prefix = directory.join("joint");
    let joint_meth_output = directory.join("joint.CGmap.gz");
    let joint_vcf_output = directory.join("joint.vcf.gz");

    run_meth_call(&input, &reference, &meth_output);
    run_snp_call(&input, &reference, &snp_output);
    run_joint_call(&input, &reference, &joint_prefix);

    assert_eq!(
        fs::read(&meth_output).expect("methylation bytes"),
        fs::read(&joint_meth_output).expect("joint methylation bytes")
    );
    assert_eq!(
        fs::read(&snp_output).expect("variant bytes"),
        fs::read(&joint_vcf_output).expect("joint variant bytes")
    );

    let expected_meth = fs::read(&meth_output).expect("expected methylation bytes");
    let expected_snp = fs::read(&snp_output).expect("expected variant bytes");
    for path in [
        &meth_output,
        &snp_output,
        &joint_meth_output,
        &joint_vcf_output,
    ] {
        fs::write(path, b"existing output\n").expect("existing output fixture");
    }
    run_meth_call(&input, &reference, &meth_output);
    run_snp_call(&input, &reference, &snp_output);
    run_joint_call(&input, &reference, &joint_prefix);
    assert_eq!(
        fs::read(&meth_output).expect("replaced meth"),
        expected_meth
    );
    assert_eq!(fs::read(&snp_output).expect("replaced SNP"), expected_snp);
    assert_eq!(
        fs::read(&joint_meth_output).expect("replaced joint meth"),
        expected_meth
    );
    assert_eq!(
        fs::read(&joint_vcf_output).expect("replaced joint VCF"),
        expected_snp
    );

    let mut decoded = String::new();
    let mut reader = DecodedReader::open(&snp_output).expect("VCF opens");
    reader.read_to_string(&mut decoded).expect("VCF decodes");
    reader.close().expect("VCF closes");
    assert!(decoded.contains("##source=bsbit\n"));
    assert!(decoded.contains("\nchr1\t1\t.\tA\tG\t"));
    fs::remove_dir_all(directory).expect("call fixture cleanup");
}

#[test]
fn index_rejects_a_reference_that_cannot_be_written_to_alignment_bam() {
    let directory = unique_directory("invalid-alignment-reference");
    fs::create_dir(&directory).expect("fresh directory");
    let reference = directory.join("reference.fa");
    let output = directory.join("reference.bsbit");
    let internal = internal_index_prefix(&output);
    fs::write(&reference, b">chr,invalid\nACGT\n").expect("invalid reference fixture");
    fs::write(&output, b"previous index").expect("existing output fixture");

    let failure = index(&reference, &output);
    assert_eq!(failure.status.code(), Some(1));
    let message = String::from_utf8(failure.stderr).expect("index error is UTF-8");
    assert!(message.contains("cannot be represented in alignment BAM output"));
    assert!(message.contains("correct the FASTA and rerun `bsbit index`"));
    assert_eq!(
        fs::metadata(&output).expect("direct output exists").len(),
        0
    );
    assert!(!internal.exists());

    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[test]
fn index_overwrites_existing_public_and_internal_files() {
    let directory = unique_directory("opaque-index");
    fs::create_dir(&directory).expect("fresh directory");
    let reference = directory.join("reference.fa");
    let output = directory.join("reference.bsbit");
    let internal = internal_index_prefix(&output);
    fs::write(&reference, b">chr\nACGTACGT\n").expect("reference fixture");
    fs::write(&internal, b"caller-owned").expect("internal collision fixture");

    assert_success(&index(&reference, &output));
    assert!(output.is_file());
    assert!(
        internal.is_file(),
        "internal search data was built by index"
    );
    assert_ne!(
        fs::read(&internal).expect("replaced internal component"),
        b"caller-owned"
    );
    assert_success(&index(&reference, &output));

    let hidden_subcommand = run([
        OsString::from("index"),
        OsString::from("combined"),
        OsString::from("--snapshot"),
        output.as_os_str().to_owned(),
    ]);
    assert_eq!(hidden_subcommand.status.code(), Some(2));

    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[test]
// Keep these layout and output assertions on one built index; splitting
// them would multiply the most expensive setup in this process-level suite.
#[allow(clippy::too_many_lines)]
fn standard_align_selects_single_or_paired_layout_and_writes_complete_bam() {
    let directory = unique_directory("standard-align-layouts");
    fs::create_dir(&directory).expect("fresh directory");
    let reference = directory.join("reference.fa");
    let index_path = directory.join("reference.bsbit");
    let single_reads = directory.join("single.fq");
    let read1 = directory.join("r1.fq");
    let read2 = directory.join("r2.fq");
    let single_bam = directory.join("single.bam");
    let single_sensitive_bam = directory.join("single-sensitive.bam");
    let single_nondirectional_bam = directory.join("single-nondirectional.bam");
    let single_bismark_bam = directory.join("single-bismark.bam");
    let single_metrics_bam = directory.join("single-metrics.bam");
    let single_mapped_only_bam = directory.join("single-mapped-only.bam");
    let unmapped_reads = directory.join("unmapped.fq");
    let paired_bam = directory.join("paired.bam");
    let metrics_bam = directory.join("metrics.bam");

    fs::write(&reference, b">chr\nAACCGTGATCTAGGCTTACGGAAT\n").expect("reference");
    fs::write(&single_reads, b"@single\nCCGTGA\n+\nIIIIII\n").expect("single reads");
    fs::write(&read1, b"@pair/1\nCCGTGA\n+\nIIIIII\n").expect("R1");
    fs::write(&read2, b"@pair/2\nTCCGTA\n+\nJJJJJJ\n").expect("R2");
    fs::write(&unmapped_reads, b"@unmapped\nNNNNNN\n+\nIIIIII\n").expect("unmapped single read");

    let index_result = index(&reference, &index_path);
    assert_success(&index_result);
    assert_human_progress(&index_result, "index");
    let index_log = String::from_utf8_lossy(&index_result.stdout);
    assert!(index_log.contains("completed bases=24 "));
    assert!(index_log.contains("contigs=1"));

    let single_result = align(&index_path, &single_reads, None, &single_bam);
    assert_success(&single_result);
    assert_human_progress(&single_result, "align");
    assert!(String::from_utf8_lossy(&single_result.stdout).contains("completed reads=1 "));

    let single_sensitive =
        align_single_sensitive(&index_path, &single_reads, &single_sensitive_bam);
    assert_success(&single_sensitive);
    assert_human_progress(&single_sensitive, "align");
    assert!(
        String::from_utf8_lossy(&single_sensitive.stdout)
            .contains("layout=single-end mapping_threads=1 compression_threads=1")
    );
    assert_eq!(decode_process_bam(&single_sensitive_bam).len(), 1);

    let single_nondirectional = run([
        OsString::from("align"),
        OsString::from("-x"),
        index_path.as_os_str().to_owned(),
        OsString::from("-1"),
        single_reads.as_os_str().to_owned(),
        OsString::from("-o"),
        single_nondirectional_bam.as_os_str().to_owned(),
        OsString::from("--non-directional"),
    ]);
    assert_success(&single_nondirectional);
    let nondirectional_header = decode_process_bam_header(&single_nondirectional_bam);
    assert!(
        nondirectional_header
            .windows(b"read-layout=single-end;library-profile=non-directional".len())
            .any(|window| window == b"read-layout=single-end;library-profile=non-directional")
    );

    let single_bismark = run([
        OsString::from("align"),
        OsString::from("--index"),
        index_path.as_os_str().to_owned(),
        OsString::from("--read1"),
        single_reads.as_os_str().to_owned(),
        OsString::from("--output"),
        single_bismark_bam.as_os_str().to_owned(),
        OsString::from("--output-contract"),
        OsString::from("bismark"),
    ]);
    assert_success(&single_bismark);
    assert_eq!(decode_process_bam(&single_bismark_bam).len(), 1);

    let single_metrics = run([
        OsString::from("align"),
        OsString::from("--index"),
        index_path.as_os_str().to_owned(),
        OsString::from("--read1"),
        single_reads.as_os_str().to_owned(),
        OsString::from("--output"),
        single_metrics_bam.as_os_str().to_owned(),
        OsString::from("--metrics"),
    ]);
    assert_single_metrics_only(single_metrics);

    let single_mapped_only = run([
        OsString::from("align"),
        OsString::from("--index"),
        index_path.as_os_str().to_owned(),
        OsString::from("--read1"),
        unmapped_reads.as_os_str().to_owned(),
        OsString::from("--output"),
        single_mapped_only_bam.as_os_str().to_owned(),
        OsString::from("--mapped-only"),
    ]);
    assert_success(&single_mapped_only);
    assert!(decode_process_bam(&single_mapped_only_bam).is_empty());

    let paired_result = align(&index_path, &read1, Some(&read2), &paired_bam);
    assert_success(&paired_result);
    assert_human_progress(&paired_result, "align");
    let paired_log = String::from_utf8_lossy(&paired_result.stdout);
    assert!(paired_log.contains("completed pairs=1 "));
    assert!(paired_log.contains("reads=2"));

    let metrics = run([
        OsString::from("align"),
        OsString::from("--index"),
        index_path.as_os_str().to_owned(),
        OsString::from("--read1"),
        read1.as_os_str().to_owned(),
        OsString::from("--read2"),
        read2.as_os_str().to_owned(),
        OsString::from("--output"),
        metrics_bam.as_os_str().to_owned(),
        OsString::from("--metrics"),
    ]);
    assert_metrics_only(metrics);

    let single = decode_process_bam(&single_bam);
    assert_eq!(single.len(), 1);
    assert_eq!(single[0][0], b"single");
    let single_flag = std::str::from_utf8(&single[0][1])
        .expect("flag UTF-8")
        .parse::<u16>()
        .expect("flag integer");
    assert_eq!(single_flag & 0x1, 0);

    let paired = decode_process_bam(&paired_bam);
    assert_eq!(paired.len(), 2);
    assert_eq!(paired[0][0], b"pair");
    assert_eq!(paired[1][0], b"pair");
    let first_flag = std::str::from_utf8(&paired[0][1])
        .expect("R1 flag UTF-8")
        .parse::<u16>()
        .expect("R1 flag integer");
    let second_flag = std::str::from_utf8(&paired[1][1])
        .expect("R2 flag UTF-8")
        .parse::<u16>()
        .expect("R2 flag integer");
    assert_eq!(first_flag & 0x41, 0x41);
    assert_eq!(second_flag & 0x81, 0x81);

    let single_header = decode_process_bam_header(&single_bam);
    assert!(
        single_header
            .windows(b"read-layout=single-end;library-profile=directional".len())
            .any(|window| window == b"read-layout=single-end;library-profile=directional")
    );
    let paired_header = decode_process_bam_header(&paired_bam);
    assert!(
        paired_header
            .windows(b"read-layout=paired-end;library-profile=directional".len())
            .any(|window| window == b"read-layout=paired-end;library-profile=directional")
    );

    let occupied = directory.join("occupied.bam");
    fs::write(&occupied, b"caller-owned").expect("occupied target");
    assert_success(&align(&index_path, &single_reads, None, &occupied));
    let replaced = decode_process_bam(&occupied);
    assert_eq!(replaced.len(), 1);
    assert_eq!(replaced[0][0], b"single");

    let malformed = directory.join("malformed.fq");
    let partial = directory.join("partial.bam");
    fs::write(&malformed, b"@broken\nACGT\n+\n").expect("malformed reads");
    fs::write(&partial, b"previous BAM").expect("existing output fixture");
    let rejected = align(&index_path, &malformed, None, &partial);
    assert_eq!(rejected.status.code(), Some(1));
    assert!(partial.exists());
    assert_ne!(
        fs::read(&partial).expect("partial BAM remains"),
        b"previous BAM"
    );

    fs::remove_dir_all(directory).expect("fixture cleanup");
}

fn assert_single_adapter_soft_clip(output: &Path, alignment_mode: &[u8]) {
    let records = decode_process_bam(output);
    assert_eq!(records.len(), 1);
    assert_eq!(records[0][0], b"adapter");
    assert_eq!(records[0][5], b"70M13S");
    assert_eq!(records[0][9].len(), 83);
    let header = decode_process_bam_header(output);
    assert!(
        header
            .windows(alignment_mode.len())
            .any(|window| window == alignment_mode)
    );
}

fn assert_single_adapter_metrics(metrics: Output) {
    assert_success(&metrics);
    let stdout = String::from_utf8(metrics.stdout).expect("single metrics TSV is UTF-8");
    let mut lines = stdout.lines();
    let header = lines.next().unwrap().split('\t').collect::<Vec<_>>();
    let values = lines.next().unwrap().split('\t').collect::<Vec<_>>();
    let field = |name| {
        let index = header
            .iter()
            .position(|candidate| *candidate == name)
            .unwrap_or_else(|| panic!("missing adapter metric {name:?}"));
        values[index]
    };
    assert_eq!(field("adapter_unique_reads"), "1");
    assert_eq!(field("adapter_clipped_bases"), "13");
    assert_eq!(field("direct_ungapped_records"), "1");
    assert_eq!(field("traceback_records"), "0");
}

fn assert_single_configurable_clipping(
    index_path: &Path,
    reads: &Path,
    custom_output: &Path,
    unclipped_output: &Path,
) {
    assert_success(&run([
        OsString::from("align"),
        OsString::from("-x"),
        index_path.as_os_str().to_owned(),
        OsString::from("-1"),
        reads.as_os_str().to_owned(),
        OsString::from("-o"),
        custom_output.as_os_str().to_owned(),
        OsString::from("--adapter"),
        OsString::from("agatcggaagagc"),
        OsString::from("--adapter-min-overlap"),
        OsString::from("13"),
        OsString::from("--adapter-max-clip"),
        OsString::from("13"),
        OsString::from("--soft-clip"),
        OsString::from("adapter"),
        OsString::from("--max-soft-clip"),
        OsString::from("13"),
    ]));
    assert_single_adapter_soft_clip(
        custom_output,
        b"read-layout=single-end;library-profile=directional",
    );

    assert_success(&run([
        OsString::from("align"),
        OsString::from("-x"),
        index_path.as_os_str().to_owned(),
        OsString::from("-1"),
        reads.as_os_str().to_owned(),
        OsString::from("-o"),
        unclipped_output.as_os_str().to_owned(),
        OsString::from("--soft-clip"),
        OsString::from("none"),
    ]));
    let unclipped = decode_process_bam(unclipped_output);
    assert_eq!(unclipped.len(), 1);
    assert_eq!(unclipped[0][5], b"*");
}

#[test]
fn single_end_exact_three_prime_adapter_recovery_emits_a_soft_clipped_record() {
    let directory = unique_directory("single-adapter-recovery");
    fs::create_dir(&directory).expect("fresh directory");
    let reference = directory.join("reference.fa");
    let index_path = directory.join("reference.bsbit");
    let reads = directory.join("adapter.fq");
    let output = directory.join("adapter.bam");
    let nondirectional_output = directory.join("adapter-nondirectional.bam");
    let bismark_output = directory.join("adapter-bismark.bam");
    let metrics_output = directory.join("adapter-metrics.bam");
    let custom_output = directory.join("adapter-custom.bam");
    let unclipped_output = directory.join("adapter-unclipped.bam");

    let mut reference_bases = Vec::with_capacity(180);
    let mut state = 0x6a09_e667_f3bc_c909_u64;
    for _ in 0..180 {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        reference_bases.push(b"ACGT"[usize::try_from((state >> 32) & 3).unwrap()]);
    }
    let mut reference_document = b">chr\n".to_vec();
    reference_document.extend_from_slice(&reference_bases);
    reference_document.push(b'\n');
    fs::write(&reference, reference_document).expect("reference fixture");

    let retained = &reference_bases[40..110];
    let adapter = b"AGATCGGAAGAGC";
    let mut adapter_read_document = b"@adapter\n".to_vec();
    adapter_read_document.extend_from_slice(retained);
    adapter_read_document.extend_from_slice(adapter);
    adapter_read_document.extend_from_slice(b"\n+\n");
    adapter_read_document.extend(std::iter::repeat_n(b'I', retained.len() + adapter.len()));
    adapter_read_document.push(b'\n');
    fs::write(&reads, adapter_read_document).expect("adapter read fixture");

    assert_success(&run([
        OsString::from("index"),
        OsString::from("-r"),
        reference.as_os_str().to_owned(),
        OsString::from("-o"),
        index_path.as_os_str().to_owned(),
        OsString::from("--index-speed"),
        OsString::from("fast"),
    ]));
    assert_success(&align(&index_path, &reads, None, &output));
    assert_single_adapter_soft_clip(
        &output,
        b"read-layout=single-end;library-profile=directional",
    );

    assert_single_configurable_clipping(&index_path, &reads, &custom_output, &unclipped_output);

    let nondirectional = align_single_nondirectional(&index_path, &reads, &nondirectional_output);
    assert_success(&nondirectional);
    assert_single_adapter_soft_clip(
        &nondirectional_output,
        b"read-layout=single-end;library-profile=non-directional",
    );

    let bismark = run([
        OsString::from("align"),
        OsString::from("--index"),
        index_path.as_os_str().to_owned(),
        OsString::from("--read1"),
        reads.as_os_str().to_owned(),
        OsString::from("--output"),
        bismark_output.as_os_str().to_owned(),
        OsString::from("--output-contract"),
        OsString::from("bismark"),
    ]);
    assert_success(&bismark);
    let bismark_record = decode_process_bam(&bismark_output).remove(0);
    assert_eq!(bismark_record[5], b"70M13S");
    assert!(
        bismark_record
            .iter()
            .any(|field| field.starts_with(b"MD:Z:"))
    );
    assert!(
        bismark_record
            .iter()
            .any(|field| field.starts_with(b"XM:Z:"))
    );

    assert_single_adapter_metrics(run([
        OsString::from("align"),
        OsString::from("-x"),
        index_path.as_os_str().to_owned(),
        OsString::from("-1"),
        reads.as_os_str().to_owned(),
        OsString::from("-o"),
        metrics_output.as_os_str().to_owned(),
        OsString::from("--metrics"),
    ]));

    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[cfg(target_os = "linux")]
fn create_directory_with_absolute_length(root: &Path, expected_length: usize) -> PathBuf {
    use std::os::unix::ffi::OsStrExt;

    let mut current = root.to_path_buf();
    loop {
        let current_length = current.as_os_str().as_bytes().len();
        if current_length == expected_length {
            return current;
        }
        let remaining = expected_length
            .checked_sub(current_length)
            .expect("requested path length exceeds root length");
        assert!(remaining >= 2, "path needs room for slash and component");
        let component_length = if remaining == 202 {
            199
        } else {
            (remaining - 1).min(200)
        };
        current.push("d".repeat(component_length));
        fs::create_dir(&current).expect("deep path component");
    }
}

#[cfg(target_os = "linux")]
fn assert_complete_path_boundaries(directory: &Path, missing_reference: &Path) {
    use std::os::unix::ffi::OsStrExt;

    let exact_root = directory.join("exact-path");
    fs::create_dir(&exact_root).expect("exact-path root");
    let exact_parent = create_directory_with_absolute_length(&exact_root, 4_093);
    let exact_target = exact_parent.join("x");
    assert_eq!(exact_target.as_os_str().as_bytes().len(), 4_095);
    let exact_failure = index(missing_reference, &exact_target);
    assert_eq!(exact_failure.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&exact_failure.stderr).contains("open reference"));
    assert_eq!(
        fs::metadata(&exact_target)
            .expect("direct boundary output exists")
            .len(),
        0
    );

    let next_root = directory.join("next-path");
    fs::create_dir(&next_root).expect("next-path root");
    let next_parent = create_directory_with_absolute_length(&next_root, 4_094);
    let next_target = next_parent.join("x");
    assert_eq!(next_target.as_os_str().as_bytes().len(), 4_096);
    let next_failure = index(missing_reference, &next_target);
    assert_eq!(next_failure.status.code(), Some(1));
    assert!(
        String::from_utf8_lossy(&next_failure.stderr).contains("open output"),
        "{}",
        String::from_utf8_lossy(&next_failure.stderr)
    );
    assert!(!next_target.exists());
}

#[test]
fn output_component_and_read_limits_have_exact_process_boundaries() {
    let directory = unique_directory("limits");
    fs::create_dir(&directory).expect("fresh directory");
    let reference = directory.join("reference.fa");
    fs::write(&reference, b">chr\nACGT\n").expect("reference fixture");

    let maximum_component_target = directory.join("x".repeat(255));
    assert_success(&index(&reference, &maximum_component_target));
    assert!(maximum_component_target.is_file());

    let overlong_target = directory.join("x".repeat(256));
    let missing_reference = directory.join("missing.fa");
    let path_failure = index(&missing_reference, &overlong_target);
    assert_eq!(path_failure.status.code(), Some(1));
    assert!(!overlong_target.exists());
    assert!(String::from_utf8_lossy(&path_failure.stderr).contains("open output"));

    #[cfg(target_os = "linux")]
    assert_complete_path_boundaries(&directory, &missing_reference);

    let exact_reads = directory.join("exact.fastq");
    let exact_target = directory.join("exact.bam");
    let exact_name = "q".repeat(254);
    let mut exact_fastq = format!("@{exact_name}\n").into_bytes();
    exact_fastq.resize(exact_fastq.len() + 192, b'N');
    exact_fastq.extend_from_slice(b"\n+\n");
    exact_fastq.resize(exact_fastq.len() + 192, b'I');
    exact_fastq.push(b'\n');
    fs::write(&exact_reads, exact_fastq).expect("exact-limit read fixture");
    assert_success(&align(
        &maximum_component_target,
        &exact_reads,
        None,
        &exact_target,
    ));
    let exact_records = decode_process_bam(&exact_target);
    assert_eq!(exact_records.len(), 1);
    assert_eq!(exact_records[0][0].len(), 254);
    assert_eq!(exact_records[0][9].len(), 192);
    assert_eq!(exact_records[0][10].len(), 192);

    let long_name_reads = directory.join("long-name.fastq");
    let long_name_target = directory.join("long-name.bam");
    fs::write(
        &long_name_reads,
        format!("@{}\nAAA\n+\nIII\n", "q".repeat(255)),
    )
    .expect("overlong-name read fixture");
    let name_failure = align(
        &maximum_component_target,
        &long_name_reads,
        None,
        &long_name_target,
    );
    assert_eq!(name_failure.status.code(), Some(1));
    assert!(long_name_target.exists());

    let oversized_reads = directory.join("oversized.fastq");
    let oversized_target = directory.join("oversized.bam");
    let mut oversized = b"@oversized\n".to_vec();
    oversized.resize(oversized.len() + 193, b'A');
    oversized.extend_from_slice(b"\n+\n");
    oversized.resize(oversized.len() + 193, b'I');
    oversized.push(b'\n');
    fs::write(&oversized_reads, oversized).expect("oversized read fixture");
    let read_failure = align(
        &maximum_component_target,
        &oversized_reads,
        None,
        &oversized_target,
    );
    assert_eq!(read_failure.status.code(), Some(1));
    assert!(oversized_target.exists());

    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[cfg(unix)]
#[test]
fn input_and_output_permission_failures_leave_direct_outputs() {
    use std::os::unix::fs::PermissionsExt;

    let directory = unique_directory("permissions");
    fs::create_dir(&directory).expect("fresh directory");
    let reference = directory.join("reference.fa");
    let index_path = directory.join("reference.bsbit");
    let reads = directory.join("reads.fastq");
    fs::write(&reference, b">chr\nTTTACGTAAA\n").expect("reference fixture");
    fs::write(&reads, b"@read\nACGT\n+\nIIII\n").expect("read fixture");

    fs::set_permissions(&reference, fs::Permissions::from_mode(0o000))
        .expect("reference permissions");
    if fs::File::open(&reference).is_ok() {
        fs::set_permissions(&reference, fs::Permissions::from_mode(0o600))
            .expect("restore reference permissions after elevated-user probe");
        fs::remove_dir_all(directory).expect("fixture cleanup after elevated-user probe");
        assert!(
            std::env::var_os("BSBIT_REQUIRE_PERMISSION_DENIAL").is_none(),
            "permission denial was required, but this user can read a mode-000 file"
        );
        eprintln!("permission-denial test skipped: this user can read a mode-000 file");
        return;
    }
    let unreadable_index = directory.join("unreadable-reference.bsbit");
    let reference_failure = index(&reference, &unreadable_index);
    assert_eq!(reference_failure.status.code(), Some(1));
    assert_eq!(
        fs::metadata(&unreadable_index)
            .expect("direct index output exists")
            .len(),
        0
    );
    fs::set_permissions(&reference, fs::Permissions::from_mode(0o600))
        .expect("restore reference permissions");
    assert_success(&index(&reference, &index_path));

    fs::set_permissions(&reads, fs::Permissions::from_mode(0o000)).expect("read permissions");
    let unreadable_output = directory.join("unreadable-reads.bam");
    let read_failure = align(&index_path, &reads, None, &unreadable_output);
    assert_eq!(read_failure.status.code(), Some(1));
    assert_eq!(
        fs::metadata(&unreadable_output)
            .expect("direct BAM output exists")
            .len(),
        0
    );
    fs::set_permissions(&reads, fs::Permissions::from_mode(0o600))
        .expect("restore read permissions");

    let unwritable_target = directory.join("unwritable-existing.bam");
    fs::write(&unwritable_target, b"previous BAM").expect("existing output fixture");
    fs::set_permissions(&unwritable_target, fs::Permissions::from_mode(0o400))
        .expect("existing output permissions");
    let overwrite_failure = align(&index_path, &reads, None, &unwritable_target);
    assert_eq!(overwrite_failure.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&overwrite_failure.stderr).contains("open output"));
    assert_eq!(
        fs::read(&unwritable_target).expect("unwritable output remains unchanged"),
        b"previous BAM"
    );
    fs::set_permissions(&unwritable_target, fs::Permissions::from_mode(0o600))
        .expect("restore existing output permissions");

    let output_directory = directory.join("readonly-output");
    fs::create_dir(&output_directory).expect("output directory");
    fs::set_permissions(&output_directory, fs::Permissions::from_mode(0o555))
        .expect("output permissions");
    let target = output_directory.join("output.bam");
    let output_failure = align(&index_path, &reads, None, &target);
    assert_eq!(output_failure.status.code(), Some(1));
    assert!(!target.exists());
    assert_eq!(
        fs::read_dir(&output_directory)
            .expect("readonly directory can be listed")
            .count(),
        0
    );
    fs::set_permissions(&output_directory, fs::Permissions::from_mode(0o700))
        .expect("restore output permissions");

    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[test]
fn output_access_is_checked_before_inputs_are_opened() {
    let missing_parent = unique_directory("output-access");
    let missing_reference = missing_parent.join("reference.fa");
    let missing_index = missing_parent.join("reference.bsbit");
    let missing_reads = missing_parent.join("reads.fastq");

    let index_output = missing_parent.join("index").join("reference.bsbit");
    let index_failure = index(&missing_reference, &index_output);
    assert_eq!(index_failure.status.code(), Some(1));
    let index_error = String::from_utf8_lossy(&index_failure.stderr);
    assert!(index_error.contains("open output"));
    assert!(index_error.contains(&index_output.to_string_lossy()[..]));

    let align_output = missing_parent.join("align").join("sample.bam");
    let align_failure = align(&missing_index, &missing_reads, None, &align_output);
    assert_eq!(align_failure.status.code(), Some(1));
    let align_error = String::from_utf8_lossy(&align_failure.stderr);
    assert!(align_error.contains("open output"));
    assert!(align_error.contains(&align_output.to_string_lossy()[..]));

    let call_output = missing_parent.join("call").join("calls.CGmap");
    let call_failure = run([
        OsString::from("call"),
        OsString::from("meth"),
        OsString::from("--input"),
        missing_parent.join("reads.bam").into_os_string(),
        OsString::from("--reference"),
        missing_reference.into_os_string(),
        OsString::from("--output"),
        call_output.as_os_str().to_owned(),
        OsString::from("--format"),
        OsString::from("cgmap"),
    ]);
    assert_eq!(call_failure.status.code(), Some(1));
    let call_error = String::from_utf8_lossy(&call_failure.stderr);
    assert!(call_error.contains("open output"));
    assert!(call_error.contains(&call_output.to_string_lossy()[..]));
}

#[test]
#[allow(clippy::too_many_lines)]
fn direct_outputs_never_truncate_inputs_or_hard_link_aliases() {
    let directory = unique_directory("direct-output-input-collision");
    fs::create_dir(&directory).expect("fixture directory");

    let reference = directory.join("reference.fa");
    fs::write(&reference, b">chr1\nACGTACGT\n").expect("reference");
    let same_reference = run([
        OsStr::new("index"),
        OsStr::new("-r"),
        reference.as_os_str(),
        OsStr::new("-o"),
        reference.as_os_str(),
    ]);
    assert_eq!(same_reference.status.code(), Some(1));
    assert_eq!(
        fs::read(&reference).expect("preserved reference"),
        b">chr1\nACGTACGT\n"
    );

    let reads = directory.join("reads.fastq");
    let reads_alias = directory.join("reads-alias.fastq");
    let missing_index = directory.join("missing.bsbit");
    fs::write(&reads, b"@r1\nACGT\n+\nIIII\n").expect("reads");
    fs::hard_link(&reads, &reads_alias).expect("reads hard link");
    let aligned_over_reads = run([
        OsStr::new("align"),
        OsStr::new("-x"),
        missing_index.as_os_str(),
        OsStr::new("-1"),
        reads.as_os_str(),
        OsStr::new("-o"),
        reads_alias.as_os_str(),
    ]);
    assert_eq!(aligned_over_reads.status.code(), Some(1));
    assert_eq!(
        fs::read(&reads).expect("preserved reads"),
        b"@r1\nACGT\n+\nIIII\n"
    );

    let mut hash = 0xcbf2_9ce4_8422_2325_u64;
    for byte in missing_index.file_name().unwrap().as_encoded_bytes() {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    let internal_index = directory.join(format!(".bsbit-index-{hash:016x}"));
    fs::write(&internal_index, b"internal index metadata").expect("internal index placeholder");
    let aligned_over_index_component = run([
        OsStr::new("align"),
        OsStr::new("-x"),
        missing_index.as_os_str(),
        OsStr::new("-1"),
        reads.as_os_str(),
        OsStr::new("-o"),
        internal_index.as_os_str(),
    ]);
    assert_eq!(aligned_over_index_component.status.code(), Some(1));
    assert_eq!(
        fs::read(&internal_index).expect("preserved internal index"),
        b"internal index metadata"
    );

    let bam = directory.join("input.bam");
    fs::write(&bam, b"not-a-bam").expect("BAM placeholder");
    let called_over_bam = run([
        OsStr::new("call"),
        OsStr::new("meth"),
        OsStr::new("-i"),
        bam.as_os_str(),
        OsStr::new("-r"),
        reference.as_os_str(),
        OsStr::new("-o"),
        bam.as_os_str(),
        OsStr::new("-f"),
        OsStr::new("cgmap"),
    ]);
    assert_eq!(called_over_bam.status.code(), Some(1));
    assert_eq!(fs::read(&bam).expect("preserved BAM"), b"not-a-bam");

    let missing_bam = directory.join("missing-input.bam");
    let bam_index = directory.join("missing-input.bam.bai");
    fs::write(&bam_index, b"BAM index placeholder").expect("BAM index placeholder");
    let called_over_bam_index = run([
        OsStr::new("call"),
        OsStr::new("meth"),
        OsStr::new("-i"),
        missing_bam.as_os_str(),
        OsStr::new("-r"),
        reference.as_os_str(),
        OsStr::new("-o"),
        bam_index.as_os_str(),
        OsStr::new("-f"),
        OsStr::new("cgmap"),
    ]);
    assert_eq!(called_over_bam_index.status.code(), Some(1));
    assert_eq!(
        fs::read(&bam_index).expect("preserved BAM index"),
        b"BAM index placeholder"
    );

    let prefix = directory.join("joint");
    let meth_output = directory.join("joint.CGmap");
    let vcf_output = directory.join("joint.vcf");
    fs::write(&meth_output, b"existing joint output").expect("joint output");
    fs::hard_link(&meth_output, &vcf_output).expect("joint output hard link");
    let joint_alias = run([
        OsStr::new("call"),
        OsStr::new("joint"),
        OsStr::new("-i"),
        directory.join("missing.bam").as_os_str(),
        OsStr::new("-r"),
        reference.as_os_str(),
        OsStr::new("-p"),
        prefix.as_os_str(),
        OsStr::new("-f"),
        OsStr::new("cgmap"),
    ]);
    assert_eq!(joint_alias.status.code(), Some(1));
    assert_eq!(
        fs::read(&meth_output).expect("preserved joint output"),
        b"existing joint output"
    );

    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[test]
fn direct_commands_truncate_outputs_before_missing_inputs() {
    let directory = unique_directory("direct-output-filesystem-smoke");
    fs::create_dir(&directory).expect("fixture directory");

    let cases = [
        (
            directory.join("index.bsbit"),
            vec![
                OsString::from("index"),
                OsString::from("-r"),
                directory.join("missing.fa").into_os_string(),
                OsString::from("-o"),
                directory.join("index.bsbit").into_os_string(),
            ],
        ),
        (
            directory.join("align.bam"),
            vec![
                OsString::from("align"),
                OsString::from("-x"),
                directory.join("missing.bsbit").into_os_string(),
                OsString::from("-1"),
                directory.join("missing.fastq").into_os_string(),
                OsString::from("-o"),
                directory.join("align.bam").into_os_string(),
            ],
        ),
        (
            directory.join("call.CGmap"),
            vec![
                OsString::from("call"),
                OsString::from("meth"),
                OsString::from("-i"),
                directory.join("missing.bam").into_os_string(),
                OsString::from("-r"),
                directory.join("missing.fa").into_os_string(),
                OsString::from("-o"),
                directory.join("call.CGmap").into_os_string(),
                OsString::from("-f"),
                OsString::from("cgmap"),
            ],
        ),
    ];
    for (output, arguments) in cases {
        fs::write(&output, b"old output").expect("old output");
        let result = run(arguments);
        assert_eq!(result.status.code(), Some(1));
        assert_eq!(fs::metadata(&output).expect("output metadata").len(), 0);
    }

    fs::remove_dir_all(directory).expect("fixture cleanup");
}
