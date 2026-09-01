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
    let mut digest = bsbit_core::reference::ReferenceSemanticDigestBuilder::new(1);
    digest
        .push_ascii_contig(b"chr1", REFERENCE)
        .expect("fixture semantic digest input");
    let header = fixture_sam_header(&reference)
        .with_bsbit_provenance(
            bsbit_hts::BsbitProgramProvenance::new(
                digest
                    .finish()
                    .expect("fixture semantic digest")
                    .into_bytes(),
                bsbit_hts::BsbitAlignmentMode::CallerCompatibleDirectionalPaired,
            ),
            AlignmentRecordLimits::default(),
        )
        .expect("fixture provenance fits")
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
                .expect("fixture header entry"),
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
        OsString::from("--output-bam"),
        output_bam.as_os_str().to_owned(),
    ]);
    run(arguments)
}

fn align_single_sensitive(snapshot: &Path, read1: &Path, output_bam: &Path) -> Output {
    run([
        OsString::from("align"),
        OsString::from("--index"),
        snapshot.as_os_str().to_owned(),
        OsString::from("-1"),
        read1.as_os_str().to_owned(),
        OsString::from("--output-bam"),
        output_bam.as_os_str().to_owned(),
        OsString::from("--sensitive"),
        OsString::from("--threads"),
        OsString::from("2"),
    ])
}

fn align_single_metrics(snapshot: &Path, read1: &Path, output_bam: &Path) -> Output {
    run([
        OsString::from("align"),
        OsString::from("--index"),
        snapshot.as_os_str().to_owned(),
        OsString::from("-1"),
        read1.as_os_str().to_owned(),
        OsString::from("--output-bam"),
        output_bam.as_os_str().to_owned(),
        OsString::from("--metrics"),
        OsString::from("--threads"),
        OsString::from("2"),
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

    let version = run(["--version"]);
    assert_eq!(version.status.code(), Some(0));
    assert_eq!(
        version.stdout,
        concat!("bsbit ", env!("CARGO_PKG_VERSION"), "\n").as_bytes()
    );
    assert!(version.stderr.is_empty());

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

    let retired = run(["align-general", "--help"]);
    assert_eq!(retired.status.code(), Some(2));
    assert_eq!(
        retired.stderr,
        b"bsbit: unknown command `align-general`; run `bsbit --help`\n"
    );
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
    let retired_output = run([
        "align",
        "--index",
        "index",
        "--read1",
        "reads",
        "--output-format",
        "cram",
    ]);
    assert_eq!(retired_output.status.code(), Some(2));
    assert_eq!(
        retired_output.stderr,
        b"bsbit: unknown option --output-format\n"
    );

    let retired_reads = run([
        "align",
        "--index",
        "index",
        "--reads",
        "reads",
        "--output-bam",
        "out.bam",
    ]);
    assert_eq!(retired_reads.status.code(), Some(2));
    assert_eq!(retired_reads.stderr, b"bsbit: unknown option --reads\n");
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
fn umbrella_call_validates_inputs_and_joint_destinations() {
    let input = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../external/htslib/test/range.bam");
    let output = unique_directory("call-validation").join("calls.cgmap.gz");
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
    assert!(String::from_utf8_lossy(&meth.stderr).contains("@PG ID:bsbit PN:bsbit"));
    assert!(!output.exists());

    let vcf = unique_directory("joint-destinations").join("same.vcf.gz");
    let joint = run([
        OsString::from("call"),
        OsString::from("joint"),
        OsString::from("-i"),
        input.as_os_str().to_owned(),
        OsString::from("--reference"),
        OsString::from("missing.fa"),
        OsString::from("--meth"),
        vcf.as_os_str().to_owned(),
        OsString::from("--meth-format"),
        OsString::from("cgmap"),
        OsString::from("--vcf"),
        vcf.as_os_str().to_owned(),
    ]);
    assert!(!joint.status.success());
    assert!(String::from_utf8_lossy(&joint.stderr).contains("different output files"));
    assert!(!vcf.exists());
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
        fs::read_to_string(&matrix).expect("combined matrix"),
        concat!(
            "##bsbit_matrix_format=count\n",
            "##bsbit_min_count=5\n",
            "##bsbit_min_prop=1.000000000\n",
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
        fs::read_to_string(&matrix).expect("combined matrix"),
        format!(
            concat!(
                "##bsbit_matrix_format=level\n",
                "##bsbit_min_count=1\n",
                "##bsbit_min_prop=0.000000000\n",
                "#chrom\tstart\tend\tmodification\tscore\tstrand\t{}\n",
                "chr1\t0\t1\tm,CG,0\t0\t+\t0.750000\n",
            ),
            input.display()
        )
    );
    fs::remove_dir_all(directory).expect("fixture cleanup");
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

fn run_joint_call(input: &Path, reference: &Path, meth_output: &Path, vcf_output: &Path) {
    let result = run([
        OsString::from("call"),
        OsString::from("joint"),
        OsString::from("-i"),
        input.as_os_str().to_owned(),
        OsString::from("--reference"),
        reference.as_os_str().to_owned(),
        OsString::from("--meth"),
        meth_output.as_os_str().to_owned(),
        OsString::from("--meth-format"),
        OsString::from("cgmap"),
        OsString::from("--vcf"),
        vcf_output.as_os_str().to_owned(),
        OsString::from("-c"),
        OsString::from("true"),
        OsString::from("-t"),
        OsString::from("2"),
    ]);
    assert_command_succeeded(&result);
}

#[test]
fn umbrella_call_modules_are_consistent_on_real_bam() {
    let directory = unique_directory("call-e2e");
    fs::create_dir(&directory).expect("call fixture directory");
    let (input, reference) = indexed_call_fixture(&directory);
    let meth_output = directory.join("meth.cgmap.gz");
    let snp_output = directory.join("snp.vcf.gz");
    let joint_meth_output = directory.join("joint.cgmap.gz");
    let joint_vcf_output = directory.join("joint.vcf.gz");

    run_meth_call(&input, &reference, &meth_output);
    run_snp_call(&input, &reference, &snp_output);
    run_joint_call(&input, &reference, &joint_meth_output, &joint_vcf_output);

    assert_eq!(
        fs::read(&meth_output).expect("methylation bytes"),
        fs::read(&joint_meth_output).expect("joint methylation bytes")
    );
    assert_eq!(
        fs::read(&snp_output).expect("variant bytes"),
        fs::read(&joint_vcf_output).expect("joint variant bytes")
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
fn removed_oracle_options_are_rejected() {
    let retired_command = run(["align-general", "--help"]);
    assert_eq!(retired_command.status.code(), Some(2));
    assert_eq!(
        retired_command.stderr,
        b"bsbit: unknown command `align-general`; run `bsbit --help`\n"
    );

    let backend = run([
        "align",
        "--index",
        "index",
        "--reference-backend",
        "scalar",
        "--read1",
        "reads",
        "--output-bam",
        "out.bam",
    ]);
    assert_eq!(backend.status.code(), Some(2));
    assert_eq!(
        backend.stderr,
        b"bsbit: unknown option --reference-backend\n"
    );

    let removed_cache = run(["cache", "--index", "index"]);
    assert_eq!(removed_cache.status.code(), Some(2));
    assert_eq!(
        removed_cache.stderr,
        b"bsbit: unknown command `cache`; run `bsbit --help`\n"
    );
}

#[test]
fn index_is_one_public_operation_and_rolls_back_an_incomplete_bundle() {
    let directory = unique_directory("opaque-index");
    fs::create_dir(&directory).expect("fresh directory");
    let reference = directory.join("reference.fa");
    let output = directory.join("reference.bsbit");
    let internal = internal_index_prefix(&output);
    fs::write(&reference, b">chr\nACGTACGT\n").expect("reference fixture");
    fs::write(&internal, b"caller-owned").expect("internal collision fixture");

    let failure = index(&reference, &output);
    assert_eq!(failure.status.code(), Some(1));
    assert!(!output.exists(), "logical index must roll back as a unit");
    assert_eq!(
        fs::read(&internal).expect("caller-owned collision retained"),
        b"caller-owned"
    );

    fs::remove_file(&internal).expect("remove collision fixture");
    assert_success(&index(&reference, &output));
    assert!(output.is_file());
    assert!(
        internal.is_file(),
        "internal search data was built by index"
    );

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
fn standard_align_selects_single_or_paired_layout_and_publishes_complete_bam() {
    let directory = unique_directory("standard-align-layouts");
    fs::create_dir(&directory).expect("fresh directory");
    let reference = directory.join("reference.fa");
    let index_path = directory.join("reference.bsbit");
    let single_reads = directory.join("single.fq");
    let read1 = directory.join("r1.fq");
    let read2 = directory.join("r2.fq");
    let single_bam = directory.join("single.bam");
    let paired_bam = directory.join("paired.bam");

    fs::write(&reference, b">chr\nAACCGTGATCTAGGCTTACGGAAT\n").expect("reference");
    fs::write(&single_reads, b"@single\nCCGTGA\n+\nIIIIII\n").expect("single reads");
    fs::write(&read1, b"@pair/1\nCCGTGA\n+\nIIIIII\n").expect("R1");
    fs::write(&read2, b"@pair/2\nTCCGTA\n+\nJJJJJJ\n").expect("R2");

    assert_success(&index(&reference, &index_path));
    assert_success(&align(&index_path, &single_reads, None, &single_bam));
    assert_success(&align(&index_path, &read1, Some(&read2), &paired_bam));

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
            .windows(b"alignment-mode=caller-compatible-directional-single".len())
            .any(|window| window == b"alignment-mode=caller-compatible-directional-single")
    );
    let paired_header = decode_process_bam_header(&paired_bam);
    assert!(
        paired_header
            .windows(b"alignment-mode=caller-compatible-directional-paired".len())
            .any(|window| window == b"alignment-mode=caller-compatible-directional-paired")
    );

    let occupied = directory.join("occupied.bam");
    fs::write(&occupied, b"caller-owned").expect("occupied target");
    let rejected = align(&index_path, &single_reads, None, &occupied);
    assert_eq!(rejected.status.code(), Some(1));
    assert_eq!(
        fs::read(&occupied).expect("target retained"),
        b"caller-owned"
    );

    let malformed = directory.join("malformed.fq");
    let unpublished = directory.join("unpublished.bam");
    fs::write(&malformed, b"@broken\nACGT\n+\n").expect("malformed reads");
    let rejected = align(&index_path, &malformed, None, &unpublished);
    assert_eq!(rejected.status.code(), Some(1));
    assert!(!unpublished.exists());

    fs::remove_dir_all(directory).expect("fixture cleanup");
}

fn assert_single_adapter_recovery_records(
    records: &[Vec<Vec<u8>>],
    supported: &[u8],
    unsupported: &[u8],
    clean: &[u8],
    mapped_supported: &[u8],
    mapped_retained_bases: usize,
) {
    assert_eq!(records.len(), 4);
    assert_eq!(records[0][0], b"supported");
    let supported_flag = std::str::from_utf8(&records[0][1])
        .expect("supported flag UTF-8")
        .parse::<u16>()
        .expect("supported flag integer");
    assert_eq!(
        supported_flag & 0x5,
        0,
        "recovered read is mapped and unpaired"
    );
    let supported_mapq = std::str::from_utf8(&records[0][4])
        .expect("supported MAPQ UTF-8")
        .parse::<u8>()
        .expect("supported MAPQ integer");
    assert!((1..=20).contains(&supported_mapq));
    assert_eq!(
        records[0][5],
        format!("{}M{}S", clean.len(), supported.len() - clean.len()).as_bytes()
    );
    assert_eq!(records[0][9], supported);
    assert_eq!(records[0][10], vec![b'I'; supported.len()]);

    assert_eq!(records[1][0], b"unsupported");
    let unsupported_flag = std::str::from_utf8(&records[1][1])
        .expect("unsupported flag UTF-8")
        .parse::<u16>()
        .expect("unsupported flag integer");
    assert_ne!(
        unsupported_flag & 0x4,
        0,
        "non-exact adapter remains unmapped"
    );
    assert_eq!(records[1][5], b"*");
    assert_eq!(records[1][9], unsupported);

    assert_eq!(records[2][0], b"clean");
    let clean_flag = std::str::from_utf8(&records[2][1])
        .expect("clean flag UTF-8")
        .parse::<u16>()
        .expect("clean flag integer");
    assert_eq!(
        clean_flag & 0x5,
        0,
        "complete read remains mapped and unpaired"
    );
    assert_eq!(records[2][2], records[0][2]);
    assert_eq!(records[2][3], records[0][3]);
    assert_eq!(records[2][5], format!("{}M", clean.len()).as_bytes());
    assert_eq!(records[2][9], clean);

    assert_eq!(records[3][0], b"mapped-supported");
    let mapped_flag = std::str::from_utf8(&records[3][1])
        .expect("mapped adapter flag UTF-8")
        .parse::<u16>()
        .expect("mapped adapter flag integer");
    assert_eq!(mapped_flag & 0x5, 0, "mapped adapter remains mapped");
    assert_eq!(
        records[3][5],
        format!(
            "{}M{}S",
            mapped_retained_bases,
            mapped_supported.len() - mapped_retained_bases
        )
        .as_bytes()
    );
    assert_eq!(records[3][9], mapped_supported);
    assert_eq!(records[3][10], vec![b'I'; mapped_supported.len()]);
}

#[test]
fn single_end_recovers_exact_three_prime_adapter_and_preserves_full_read() {
    const RETAINED: &[u8] = b"ACGTCAGATGCTACGAGTACCGATGACCTAGCATGCATGATCGTACGATCGTAGCTAGCATGCA";
    const MAPPED_RETAINED: &[u8] =
        b"GTCAGTGACCATGCTGACGATCGTACCTGAGTCCAGTACGATGCTAGTCAGGATCGTACGATGC";
    const ADAPTER: &[u8] = b"AGATCGGAAGAGC";
    let directory = unique_directory("single-adapter-recovery");
    fs::create_dir(&directory).expect("fresh directory");
    let reference = directory.join("reference.fa");
    let index_path = directory.join("reference.bsbit");
    let reads = directory.join("single.fq");
    let output_bam = directory.join("single.bam");
    let sensitive_bam = directory.join("single-sensitive.bam");
    let metrics_bam = directory.join("single-metrics.bam");

    let mut reference_bytes = b">chr\n".to_vec();
    reference_bytes.extend(std::iter::repeat_n(b'G', 40));
    reference_bytes.extend_from_slice(RETAINED);
    reference_bytes.extend(std::iter::repeat_n(b'T', 40));
    reference_bytes.extend(std::iter::repeat_n(b'N', 40));
    reference_bytes.extend_from_slice(MAPPED_RETAINED);
    reference_bytes.extend(std::iter::repeat_n(b'A', 8));
    reference_bytes.extend(std::iter::repeat_n(b'G', 40));
    reference_bytes.push(b'\n');
    fs::write(&reference, reference_bytes).expect("reference fixture");

    let retained_read = RETAINED
        .iter()
        .map(|base| if *base == b'C' { b'T' } else { *base })
        .collect::<Vec<_>>();
    let mut supported = retained_read.clone();
    supported.extend_from_slice(ADAPTER);
    let mut unsupported = supported.clone();
    unsupported[RETAINED.len()] = b'C';
    let clean = retained_read.clone();
    let mut mapped_supported = MAPPED_RETAINED
        .iter()
        .map(|base| if *base == b'C' { b'T' } else { *base })
        .collect::<Vec<_>>();
    mapped_supported.extend_from_slice(&ADAPTER[..8]);
    let mut fastq = Vec::new();
    for (name, sequence) in [
        (b"supported".as_slice(), &supported),
        (b"unsupported", &unsupported),
        (b"clean", &clean),
        (b"mapped-supported", &mapped_supported),
    ] {
        fastq.push(b'@');
        fastq.extend_from_slice(name);
        fastq.push(b'\n');
        fastq.extend_from_slice(sequence);
        fastq.extend_from_slice(b"\n+\n");
        fastq.extend(std::iter::repeat_n(b'I', sequence.len()));
        fastq.push(b'\n');
    }
    fs::write(&reads, fastq).expect("single adapter reads");

    assert_success(&index(&reference, &index_path));
    assert_success(&align(&index_path, &reads, None, &output_bam));

    let records = decode_process_bam(&output_bam);
    assert_single_adapter_recovery_records(
        &records,
        &supported,
        &unsupported,
        &clean,
        &mapped_supported,
        MAPPED_RETAINED.len(),
    );

    assert_success(&align_single_sensitive(&index_path, &reads, &sensitive_bam));
    let sensitive = decode_process_bam(&sensitive_bam);
    assert_single_adapter_recovery_records(
        &sensitive,
        &supported,
        &unsupported,
        &clean,
        &mapped_supported,
        MAPPED_RETAINED.len(),
    );

    let metrics = align_single_metrics(&index_path, &reads, &metrics_bam);
    assert_success(&metrics);
    let stdout = std::str::from_utf8(&metrics.stdout).expect("metrics stdout UTF-8");
    let lines = stdout.lines().collect::<Vec<_>>();
    assert_eq!(lines.len(), 2);
    let header = lines[0].split('\t').collect::<Vec<_>>();
    let fields = lines[1].split('\t').collect::<Vec<_>>();
    assert_eq!(header.len(), fields.len());
    let field = |name: &str| {
        fields[header
            .iter()
            .position(|candidate| *candidate == name)
            .expect(name)]
    };
    assert_eq!(field("schema"), "bsbit-single-alignment-metrics-v1");
    assert_eq!(field("reads"), "4");
    assert_eq!(field("unique"), "3");
    assert_eq!(field("unmapped"), "1");
    assert_eq!(field("bam_records"), "4");
    assert_eq!(field("adapter_attempted_reads"), "2");
    assert_eq!(field("adapter_unique_reads"), "2");

    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[test]
fn single_sensitive_completes_repeat_frontier_above_the_default_hit_cap() {
    const MOTIF: &[u8] = b"ACGTCAGATGCTACGAGTACCGATGACCTAGCATGCATGA";
    const COPIES: usize = 1_001;
    let directory = unique_directory("single-sensitive-repeat-frontier");
    fs::create_dir(&directory).expect("fresh directory");
    let reference = directory.join("reference.fa");
    let index_path = directory.join("reference.bsbit");
    let reads = directory.join("reads.fq");
    let default_bam = directory.join("default.bam");
    let sensitive_bam = directory.join("sensitive.bam");

    let mut reference_contents = b">chr\n".to_vec();
    for _ in 0..COPIES {
        reference_contents.extend_from_slice(MOTIF);
        reference_contents.extend_from_slice(b"NNNNNNNNNN");
    }
    reference_contents.push(b'\n');
    fs::write(&reference, reference_contents).expect("repeat reference");
    let mut read_records = Vec::new();
    for ordinal in 0..3 {
        read_records.extend_from_slice(format!("@repeat-{ordinal}\n").as_bytes());
        read_records.extend_from_slice(MOTIF);
        read_records.extend_from_slice(b"\n+\n");
        read_records.extend(std::iter::repeat_n(b'I', MOTIF.len()));
        read_records.push(b'\n');
    }
    fs::write(&reads, read_records).expect("repeat read");

    assert_success(&index(&reference, &index_path));
    assert_success(&align(&index_path, &reads, None, &default_bam));
    assert_success(&align_single_sensitive(&index_path, &reads, &sensitive_bam));

    let default = decode_process_bam(&default_bam);
    let sensitive = decode_process_bam(&sensitive_bam);
    assert_eq!(default.len(), 3);
    assert_eq!(sensitive.len(), 3);
    for (default, sensitive) in default.iter().zip(&sensitive) {
        let default_flag = std::str::from_utf8(&default[1])
            .expect("default flag UTF-8")
            .parse::<u16>()
            .expect("default flag integer");
        let sensitive_flag = std::str::from_utf8(&sensitive[1])
            .expect("sensitive flag UTF-8")
            .parse::<u16>()
            .expect("sensitive flag integer");
        assert_ne!(
            default_flag & 0x4,
            0,
            "default hit cap leaves read unmapped"
        );
        assert_eq!(sensitive_flag & 0x4, 0, "sensitive frontier maps the read");
        assert_eq!(sensitive[4], b"0", "repeat tie must use MAPQ 0");
    }

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
    let exact_parent = create_directory_with_absolute_length(&exact_root, 4_052);
    let exact_target = exact_parent.join("x");
    assert_eq!(exact_target.as_os_str().as_bytes().len(), 4_054);
    let exact_failure = index(missing_reference, &exact_target);
    assert_eq!(exact_failure.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&exact_failure.stderr).contains("open reference"));

    let next_root = directory.join("next-path");
    fs::create_dir(&next_root).expect("next-path root");
    let next_parent = create_directory_with_absolute_length(&next_root, 4_053);
    let next_target = next_parent.join("x");
    assert_eq!(next_target.as_os_str().as_bytes().len(), 4_055);
    let next_failure = index(missing_reference, &next_target);
    assert_eq!(next_failure.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&next_failure.stderr).contains("inspect staging path"));
    assert!(!exact_target.exists());
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
    assert!(String::from_utf8_lossy(&path_failure.stderr).contains("inspect destination"));

    #[cfg(target_os = "linux")]
    assert_complete_path_boundaries(&directory, &missing_reference);

    let exact_reads = directory.join("exact.fastq");
    let exact_target = directory.join("exact.bam");
    let exact_name = "q".repeat(254);
    let mut exact_fastq = format!("@{exact_name}\n").into_bytes();
    exact_fastq.resize(exact_fastq.len() + 1_000_000, b'N');
    exact_fastq.extend_from_slice(b"\n+\n");
    exact_fastq.resize(exact_fastq.len() + 1_000_000, b'I');
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
    assert_eq!(exact_records[0][9].len(), 1_000_000);
    assert_eq!(exact_records[0][10].len(), 1_000_000);

    let long_name_reads = directory.join("long-name.fastq");
    let long_name_target = directory.join("long-name.bam");
    fs::write(&long_name_reads, format!("@{}\nA\n+\nI\n", "q".repeat(255)))
        .expect("overlong-name read fixture");
    let name_failure = align(
        &maximum_component_target,
        &long_name_reads,
        None,
        &long_name_target,
    );
    assert_eq!(name_failure.status.code(), Some(1));
    assert!(!long_name_target.exists());

    let oversized_reads = directory.join("oversized.fastq");
    let oversized_target = directory.join("oversized.bam");
    let mut oversized = b"@oversized\n".to_vec();
    oversized.resize(oversized.len() + 1_000_001, b'A');
    oversized.extend_from_slice(b"\n+\n");
    oversized.resize(oversized.len() + 1_000_001, b'I');
    oversized.push(b'\n');
    fs::write(&oversized_reads, oversized).expect("oversized read fixture");
    let read_failure = align(
        &maximum_component_target,
        &oversized_reads,
        None,
        &oversized_target,
    );
    assert_eq!(read_failure.status.code(), Some(1));
    assert!(!oversized_target.exists());

    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[cfg(unix)]
#[test]
fn input_and_output_permission_failures_publish_nothing() {
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
    assert!(!unreadable_index.exists());
    fs::set_permissions(&reference, fs::Permissions::from_mode(0o600))
        .expect("restore reference permissions");
    assert_success(&index(&reference, &index_path));

    fs::set_permissions(&reads, fs::Permissions::from_mode(0o000)).expect("read permissions");
    let unreadable_output = directory.join("unreadable-reads.bam");
    let read_failure = align(&index_path, &reads, None, &unreadable_output);
    assert_eq!(read_failure.status.code(), Some(1));
    assert!(!unreadable_output.exists());
    fs::set_permissions(&reads, fs::Permissions::from_mode(0o600))
        .expect("restore read permissions");

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
