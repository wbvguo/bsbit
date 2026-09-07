//! Public API integration tests for the `bsbit-call` crate.

use std::error::Error;
use std::fs;
use std::io::Read;
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

use bsbit_core::bisulfite::BisulfiteStrand;
use bsbit_core::reference::ReferenceSequenceMd5;
use bsbit_hts::{
    AlignmentAuxiliaryMode, AlignmentCigarOp, AlignmentCigarRun, AlignmentRecordLimits,
    BamStagingWriter, BorrowedAlignmentRecord, BsbitAlignmentMode, BsbitHeaderMetadata,
    DecodedReader, SamHeader, SamHeaderReference, SamSortOrder, build_bam_index_create_new,
};

use bsbit_call::region::{GenomicInterval, RegionSelection};
use bsbit_call::{CallErrorKind, joint, meth, snp};

const FIXTURE_REFERENCE: &[u8] = b"ACGTTGCACTGATCGATGCTAGCTACGATCGTTCGAGTACCTGACGTA";

fn unique_directory(label: &str) -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock is after epoch")
        .as_nanos();
    std::env::temp_dir().join(format!("bsbit-call-{label}-{}-{nonce}", std::process::id()))
}

fn indexed_bsbit_fixture(directory: &std::path::Path) -> PathBuf {
    indexed_bsbit_fixture_with_mode(directory, AlignmentAuxiliaryMode::Minimal)
}

fn indexed_bsbit_fixture_with_mode(
    directory: &std::path::Path,
    auxiliary_mode: AlignmentAuxiliaryMode,
) -> PathBuf {
    indexed_bsbit_fixture_with_contract(
        directory,
        auxiliary_mode,
        BsbitAlignmentMode::DirectionalPairedEnd,
        true,
    )
}

fn indexed_bsbit_fixture_with_contract(
    directory: &std::path::Path,
    auxiliary_mode: AlignmentAuxiliaryMode,
    alignment_mode: BsbitAlignmentMode,
    include_reference_md5: bool,
) -> PathBuf {
    let mut observed = FIXTURE_REFERENCE.to_vec();
    observed[0] = b'G';
    let limits = AlignmentRecordLimits::default();
    let reference_length = u64::try_from(FIXTURE_REFERENCE.len()).expect("fixture length fits u64");
    let reference =
        SamHeaderReference::new(0, b"chr1", reference_length).expect("fixture dictionary entry");
    let reference = if include_reference_md5 {
        reference.with_md5(ReferenceSequenceMd5::from_ascii(FIXTURE_REFERENCE))
    } else {
        reference
    };
    let header = SamHeader::new(vec![reference], limits)
        .expect("fixture header builds")
        .with_bsbit_metadata(BsbitHeaderMetadata::new(alignment_mode), limits)
        .expect("fixture metadata fits")
        .with_sort_order(SamSortOrder::Coordinate);
    let staging = directory.join("fixture.bam.tmp");
    let input = directory.join("fixture.bam");
    let mut writer =
        BamStagingWriter::create_new(&staging, &header, limits).expect("fixture BAM opens");
    let qualities = vec![b'I'; observed.len()];
    let cigar = [
        AlignmentCigarRun::new(AlignmentCigarOp::Match, reference_length)
            .expect("fixture CIGAR run"),
    ];
    let md = format!("0A{}", reference_length - 1);
    let xm = vec![b'.'; observed.len()];
    let (md, xm) = match auxiliary_mode {
        AlignmentAuxiliaryMode::Minimal => (None, None),
        AlignmentAuxiliaryMode::Bismark => (Some(md.as_bytes()), Some(xm.as_slice())),
    };
    for ordinal in 0..24 {
        let query_name = format!("read-{ordinal:02}");
        let record = BorrowedAlignmentRecord::new(
            query_name.as_bytes(),
            0,
            Some(0),
            1,
            60,
            &cigar,
            None,
            0,
            0,
            &observed,
            Some(&qualities),
            1,
            auxiliary_mode,
            md,
            BisulfiteStrand::OT,
            xm,
            limits,
        )
        .expect("fixture BAM record builds");
        writer
            .write_borrowed_alignment_record(&record)
            .expect("fixture BAM record writes");
    }
    writer
        .finish()
        .expect("fixture BAM finishes")
        .publish_create_new(&input)
        .expect("fixture BAM publishes");
    let index = input.with_extension("bam.bai");
    build_bam_index_create_new(&input, &index, 1).expect("fixture BAI builds");
    let published_index = fs::read(&index).expect("fixture BAI bytes");
    build_bam_index_create_new(&input, &index, 1)
        .expect_err("fixture BAI publication is create-only");
    assert_eq!(
        fs::read(&index).expect("fixture BAI remains"),
        published_index
    );
    input
}

fn indexed_fasta_fixture(directory: &std::path::Path) -> PathBuf {
    let fasta = directory.join("reference.fa");
    let mut contents = b">chr1\n".to_vec();
    contents.extend_from_slice(FIXTURE_REFERENCE);
    contents.push(b'\n');
    fs::write(&fasta, contents).expect("fixture FASTA writes");
    fs::write(
        fasta.with_extension("fa.fai"),
        format!(
            "chr1\t{}\t6\t{}\t{}\n",
            FIXTURE_REFERENCE.len(),
            FIXTURE_REFERENCE.len(),
            FIXTURE_REFERENCE.len() + 1
        ),
    )
    .expect("fixture FAI writes");
    fasta
}

#[test]
#[allow(clippy::too_many_lines)]
fn module_entry_points_validate_before_opening_inputs() {
    let meth_error = meth::call(&meth::Options {
        input: PathBuf::from("missing.bam"),
        reference: PathBuf::from("missing.fa"),
        regions: RegionSelection::default(),
        output: PathBuf::from("missing.cgmap"),
        format: meth::OutputFormat::Cgmap,
        compress: false,
        threads: 0,
        compression_threads: 0,
        parameters: meth::Parameters::default(),
    })
    .expect_err("zero methylation threads must fail before input is opened");
    assert!(meth_error.to_string().contains("thread count"));
    assert_eq!(meth_error.kind(), CallErrorKind::Configuration);
    assert!(meth_error.source().is_none());

    for (compress, compression_threads) in [(false, 1), (true, u32::MAX)] {
        let error = meth::call(&meth::Options {
            input: PathBuf::from("missing.bam"),
            reference: PathBuf::from("missing.fa"),
            regions: RegionSelection::default(),
            output: PathBuf::from("missing.cgmap"),
            format: meth::OutputFormat::Cgmap,
            compress,
            threads: 1,
            compression_threads,
            parameters: meth::Parameters::default(),
        })
        .expect_err("invalid compression workers must fail before input is opened");
        assert_eq!(error.kind(), CallErrorKind::Configuration);
        assert!(error.to_string().contains("compression thread count"));
    }

    let invalid_meth_error = meth::call(&meth::Options {
        input: PathBuf::from("missing.bam"),
        reference: PathBuf::from("missing.fa"),
        regions: RegionSelection::default(),
        output: PathBuf::from("missing.cgmap"),
        format: meth::OutputFormat::Cgmap,
        compress: false,
        threads: 1,
        compression_threads: 0,
        parameters: meth::Parameters {
            minimum_base_quality: 94,
            ..meth::Parameters::default()
        },
    })
    .expect_err("invalid methylation quality must fail before input is opened");
    assert!(invalid_meth_error.to_string().contains("base quality"));
    assert_eq!(invalid_meth_error.kind(), CallErrorKind::Configuration);

    let invalid_parameters = snp::Parameters {
        minimum_depth: 0,
        ..snp::Parameters::default()
    };
    let snp_error = snp::call(&snp::Options {
        input: PathBuf::from("missing.bam"),
        reference: PathBuf::from("missing.fa"),
        sample_name: None,
        regions: RegionSelection::default(),
        output: PathBuf::from("missing.vcf"),
        compress: false,
        threads: 1,
        compression_threads: 0,
        parameters: invalid_parameters,
    })
    .expect_err("invalid SNP parameters must fail before input is opened");
    assert!(snp_error.to_string().contains("must be nonzero"));
    assert_eq!(snp_error.kind(), CallErrorKind::Configuration);

    let invalid_sample = snp::call(&snp::Options {
        input: PathBuf::from("missing.bam"),
        reference: PathBuf::from("missing.fa"),
        sample_name: Some(String::from("invalid sample")),
        regions: RegionSelection::default(),
        output: PathBuf::from("missing.vcf"),
        compress: false,
        threads: 1,
        compression_threads: 0,
        parameters: snp::Parameters::default(),
    })
    .expect_err("invalid sample names fail before input is opened");
    assert_eq!(invalid_sample.kind(), CallErrorKind::Configuration);
    assert!(invalid_sample.to_string().contains("sample name"));

    let joint_error = joint::call(&joint::Options {
        input: PathBuf::from("missing.bam"),
        reference: PathBuf::from("missing.fa"),
        sample_name: None,
        regions: RegionSelection::default(),
        meth_output: PathBuf::from("missing.cgmap"),
        meth_format: meth::OutputFormat::Cgmap,
        vcf_output: PathBuf::from("missing.vcf"),
        compress: false,
        threads: u64::from(u32::MAX) + 1,
        compression_threads: 0,
        cg_only: false,
        parameters: snp::Parameters::default(),
    })
    .expect_err("excess joint threads must fail before input is opened");
    assert!(joint_error.to_string().contains("thread count"));
    assert_eq!(joint_error.kind(), CallErrorKind::Configuration);

    let same_output = PathBuf::from("same-output.gz");
    let joint_path_error = joint::call(&joint::Options {
        input: PathBuf::from("missing.bam"),
        reference: PathBuf::from("missing.fa"),
        sample_name: None,
        regions: RegionSelection::default(),
        meth_output: same_output.clone(),
        meth_format: meth::OutputFormat::Cgmap,
        vcf_output: same_output,
        compress: true,
        threads: 1,
        compression_threads: 0,
        cg_only: false,
        parameters: snp::Parameters::default(),
    })
    .expect_err("colliding joint outputs must fail before input is opened");
    assert!(joint_path_error.to_string().contains("must differ"));
    assert_eq!(joint_path_error.kind(), CallErrorKind::Configuration);
    assert!(joint_path_error.source().is_some());

    let directory = unique_directory("missing-input");
    fs::create_dir(&directory).expect("fixture directory");
    let direct_output = directory.join("unused-output.cgmap");
    let missing_input_error = meth::call(&meth::Options {
        input: PathBuf::from("definitely-missing-input.bam"),
        reference: PathBuf::from("missing.fa"),
        regions: RegionSelection::default(),
        output: direct_output.clone(),
        format: meth::OutputFormat::Cgmap,
        compress: false,
        threads: 1,
        compression_threads: 0,
        parameters: meth::Parameters::default(),
    })
    .expect_err("missing indexed BAM must retain its HTS error source");
    assert_eq!(missing_input_error.kind(), CallErrorKind::Input);
    assert!(missing_input_error.source().is_some());
    assert_eq!(
        fs::metadata(direct_output)
            .expect("direct output exists")
            .len(),
        0
    );
    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[test]
fn reference_backed_call_rejects_same_length_wrong_fasta_even_with_md() {
    let directory = unique_directory("reference-mismatch");
    fs::create_dir(&directory).expect("fixture directory");
    let input = indexed_bsbit_fixture_with_mode(&directory, AlignmentAuxiliaryMode::Bismark);
    let reference = indexed_fasta_fixture(&directory);
    let mut contents = fs::read(&reference).expect("fixture FASTA reads");
    contents[6] = b'T';
    fs::write(&reference, contents).expect("mismatching FASTA writes");
    let output = directory.join("mismatch.cgmap");
    fs::write(&output, b"previous result\n").expect("existing output fixture");

    let error = meth::call(&meth::Options {
        input,
        reference,
        regions: RegionSelection::default(),
        output: output.clone(),
        format: meth::OutputFormat::Cgmap,
        compress: false,
        threads: 2,
        compression_threads: 0,
        parameters: meth::Parameters::default(),
    })
    .expect_err("BAM M5 must reject a same-name, same-length wrong FASTA");

    assert!(error.to_string().contains("M5"));
    assert_eq!(
        fs::metadata(&output).expect("direct output exists").len(),
        0
    );
    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[test]
fn caller_does_not_branch_on_informational_alignment_mode() {
    let directory = unique_directory("informational-alignment-mode");
    fs::create_dir(&directory).expect("fixture directory");
    let input = indexed_bsbit_fixture_with_contract(
        &directory,
        AlignmentAuxiliaryMode::Minimal,
        BsbitAlignmentMode::DirectionalSingleEnd,
        true,
    );
    let reference = indexed_fasta_fixture(&directory);
    let output = directory.join("single.cgmap");

    meth::call(&meth::Options {
        input,
        reference,
        regions: RegionSelection::default(),
        output: output.clone(),
        format: meth::OutputFormat::Cgmap,
        compress: false,
        threads: 1,
        compression_threads: 0,
        parameters: meth::Parameters::default(),
    })
    .expect("caller accepts reference-identical BAM regardless of alignment mode metadata");

    assert!(output.exists());
    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[test]
fn caller_requires_a_standard_m5_on_every_sq_record() {
    let directory = unique_directory("missing-reference-m5");
    fs::create_dir(&directory).expect("fixture directory");
    let input = indexed_bsbit_fixture_with_contract(
        &directory,
        AlignmentAuxiliaryMode::Minimal,
        BsbitAlignmentMode::DirectionalPairedEnd,
        false,
    );
    let reference = indexed_fasta_fixture(&directory);
    let output = directory.join("missing-m5.cgmap");
    fs::write(&output, b"previous result\n").expect("existing output fixture");

    let error = meth::call(&meth::Options {
        input,
        reference,
        regions: RegionSelection::default(),
        output: output.clone(),
        format: meth::OutputFormat::Cgmap,
        compress: false,
        threads: 1,
        compression_threads: 0,
        parameters: meth::Parameters::default(),
    })
    .expect_err("caller must require per-contig M5 reference identity");

    assert!(error.to_string().contains("lacks the standard M5"));
    assert_eq!(
        fs::metadata(&output).expect("direct output exists").len(),
        0
    );
    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[test]
fn unindexed_plain_fasta_succeeds_with_a_reported_warning() {
    let directory = unique_directory("unindexed-plain-fasta-warning");
    fs::create_dir(&directory).expect("fixture directory");
    let input = indexed_bsbit_fixture(&directory);
    let reference = indexed_fasta_fixture(&directory);
    fs::remove_file(reference.with_extension("fa.fai")).expect("remove fixture FAI");
    let output = directory.join("unindexed-reference.cgmap");

    let report = meth::call(&meth::Options {
        input,
        reference,
        regions: RegionSelection::default(),
        output: output.clone(),
        format: meth::OutputFormat::Cgmap,
        compress: false,
        threads: 1,
        compression_threads: 0,
        parameters: meth::Parameters::default(),
    })
    .expect("unindexed plain FASTA remains supported");

    assert!(output.exists());
    assert!(report.warnings().iter().any(|warning| {
        warning.message().contains("no adjacent .fai index")
            && warning.message().contains("temporary in-memory")
    }));
    fs::remove_dir_all(directory).expect("fixture cleanup");
}

#[test]
#[allow(clippy::too_many_lines)]
fn real_meth_snp_and_joint_calls_share_outputs() {
    let directory = unique_directory("public-e2e");
    fs::create_dir(&directory).expect("fixture directory");
    let input = indexed_bsbit_fixture(&directory);
    let reference = indexed_fasta_fixture(&directory);
    let meth_output = directory.join("meth.cgmap.gz");
    let synchronous_meth_output = directory.join("meth-synchronous.cgmap.gz");
    let parallel_meth_output = directory.join("meth-parallel.cgmap.gz");
    let restricted_meth_output = directory.join("meth-targeted.cgmap.gz");
    let snp_output = directory.join("snp.vcf.gz");
    let named_snp_output = directory.join("named.vcf.gz");
    let joint_meth_output = directory.join("joint.cgmap.gz");
    let joint_vcf_output = directory.join("joint.vcf.gz");
    let parameters = snp::Parameters::default();

    meth::call(&meth::Options {
        input: input.clone(),
        reference: reference.clone(),
        regions: RegionSelection::default(),
        output: meth_output.clone(),
        format: meth::OutputFormat::Cgmap,
        compress: true,
        threads: 2,
        compression_threads: 1,
        parameters: meth::Parameters::default(),
    })
    .expect("real methylation call succeeds");
    for (output, compression_threads) in [(&synchronous_meth_output, 0), (&parallel_meth_output, 3)]
    {
        meth::call(&meth::Options {
            input: input.clone(),
            reference: reference.clone(),
            regions: RegionSelection::default(),
            output: output.clone(),
            format: meth::OutputFormat::Cgmap,
            compress: true,
            threads: 2,
            compression_threads,
            parameters: meth::Parameters::default(),
        })
        .expect("methylation output is independent of compression worker count");
    }
    meth::call(&meth::Options {
        input: input.clone(),
        reference: reference.clone(),
        regions: RegionSelection {
            intervals: vec![
                GenomicInterval {
                    contig: String::from("chr1"),
                    start: 0,
                    end: 6,
                },
                GenomicInterval {
                    contig: String::from("chr1"),
                    start: 4,
                    end: 10,
                },
            ],
            regions_file: None,
        },
        output: restricted_meth_output.clone(),
        format: meth::OutputFormat::Cgmap,
        compress: true,
        threads: 2,
        compression_threads: 1,
        parameters: meth::Parameters::default(),
    })
    .expect("targeted methylation call succeeds");
    snp::call(&snp::Options {
        input: input.clone(),
        reference: reference.clone(),
        sample_name: None,
        regions: RegionSelection::default(),
        output: snp_output.clone(),
        compress: true,
        threads: 2,
        compression_threads: 1,
        parameters,
    })
    .expect("real SNP call succeeds");
    snp::call(&snp::Options {
        input: input.clone(),
        reference: reference.clone(),
        sample_name: Some(String::from("donor-A")),
        regions: RegionSelection::default(),
        output: named_snp_output.clone(),
        compress: true,
        threads: 2,
        compression_threads: 1,
        parameters,
    })
    .expect("explicitly named SNP call succeeds");
    joint::call(&joint::Options {
        input: input.clone(),
        reference,
        sample_name: None,
        regions: RegionSelection::default(),
        meth_output: joint_meth_output.clone(),
        meth_format: meth::OutputFormat::Cgmap,
        vcf_output: joint_vcf_output.clone(),
        compress: true,
        threads: 2,
        compression_threads: 1,
        cg_only: false,
        parameters,
    })
    .expect("real joint call succeeds");

    assert_eq!(
        fs::read(&meth_output).expect("standalone methylation output"),
        fs::read(&joint_meth_output).expect("joint methylation output")
    );
    assert_eq!(
        fs::read(&meth_output).expect("one-worker methylation output"),
        fs::read(&synchronous_meth_output).expect("synchronous methylation output")
    );
    assert_eq!(
        fs::read(&meth_output).expect("one-worker methylation output"),
        fs::read(&parallel_meth_output).expect("parallel methylation output")
    );
    assert_eq!(
        fs::read(&snp_output).expect("standalone SNP output"),
        fs::read(&joint_vcf_output).expect("joint SNP output")
    );
    let vcf = fs::read(&snp_output).expect("compressed VCF bytes");
    assert!(vcf.starts_with(&[0x1f, 0x8b]));
    assert!(fs::metadata(&meth_output).expect("CGmap output").len() > 0);
    let mut full_meth = String::new();
    let mut reader = DecodedReader::open(&meth_output).expect("full CGmap opens");
    reader
        .read_to_string(&mut full_meth)
        .expect("full CGmap decodes");
    reader.close().expect("full CGmap closes");
    let expected_targeted = full_meth
        .lines()
        .filter(|line| {
            line.split('\t')
                .nth(2)
                .and_then(|position| position.parse::<u64>().ok())
                .is_some_and(|position| position <= 10)
        })
        .collect::<Vec<_>>()
        .join("\n");
    let mut targeted_meth = String::new();
    let mut reader = DecodedReader::open(&restricted_meth_output).expect("targeted CGmap opens");
    reader
        .read_to_string(&mut targeted_meth)
        .expect("targeted CGmap decodes");
    reader.close().expect("targeted CGmap closes");
    assert!(!expected_targeted.is_empty());
    assert_eq!(targeted_meth.trim_end(), expected_targeted);
    let mut decoded_vcf = String::new();
    let mut reader = DecodedReader::open(&snp_output).expect("VCF opens");
    reader
        .read_to_string(&mut decoded_vcf)
        .expect("VCF decodes");
    reader.close().expect("VCF closes");
    assert!(decoded_vcf.contains("\tFORMAT\tfixture\n"));
    assert!(decoded_vcf.contains("\nchr1\t1\t.\tA\tG\t"));

    let mut named_vcf = String::new();
    let mut reader = DecodedReader::open(&named_snp_output).expect("named VCF opens");
    reader
        .read_to_string(&mut named_vcf)
        .expect("named VCF decodes");
    reader.close().expect("named VCF closes");
    assert!(named_vcf.contains("\tFORMAT\tdonor-A\n"));

    fs::remove_dir_all(directory).expect("fixture cleanup");
}
