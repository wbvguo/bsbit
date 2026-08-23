//! Deterministic VCF header and record rendering for caller results.

use std::io::{self, Write};

use super::result::{
    FILTER_LOW_ALLELE_QUALITY, FILTER_LOW_ALTERNATE_DEPTH, FILTER_LOW_GENOTYPE_QUALITY, SnpConfig,
    VariantCall, has_conversion_confounded_pair, informative_allele_depth, selected_alleles,
    total_allele_depth,
};

#[allow(clippy::too_many_lines)]
pub(crate) fn render_vcf_header(
    writer: &mut (impl Write + ?Sized),
    references: &[(&[u8], u32)],
    config: SnpConfig,
    sample_name: &[u8],
) -> io::Result<()> {
    writer.write_all(b"##fileformat=VCFv4.3\n")?;
    writer.write_all(b"##source=bsbit\n")?;
    render_vcf_model_header(writer)?;
    writeln!(
        writer,
        "##bsbit_min_base_quality={}",
        config.minimum_base_quality
    )?;
    writeln!(
        writer,
        "##bsbit_min_mapping_quality={}",
        config.minimum_mapping_quality
    )?;
    writeln!(writer, "##bsbit_min_depth={}", config.minimum_depth)?;
    writeln!(
        writer,
        "##bsbit_min_alternate_count={}",
        config.minimum_alternate_count
    )?;
    writeln!(
        writer,
        "##bsbit_min_alternate_fraction={:.9}",
        f64::from(config.minimum_alternate_fraction_parts_per_billion) / 1_000_000_000.0
    )?;
    writeln!(
        writer,
        "##bsbit_min_genotype_quality={}",
        config.minimum_genotype_quality
    )?;
    writeln!(
        writer,
        "##bsbit_min_allele_quality={}",
        config.minimum_allele_quality
    )?;
    writeln!(
        writer,
        "##bsbit_heterozygosity={:.9}",
        config.heterozygosity_rate
    )?;
    writeln!(
        writer,
        "##bsbit_underconversion_rate={:.9}",
        config.underconversion_rate
    )?;
    writeln!(
        writer,
        "##bsbit_overconversion_rate={:.9}",
        config.overconversion_rate
    )?;
    writer.write_all(
        b"##bsbit_methylation_marginalization=adaptive-log-concave-simpson-uniform-prior\n",
    )?;
    writer.write_all(
        b"##bsbit_overlap_policy=quality-eligible-then-canonical-base-then-combined-base-mapping-quality-then-r1\n",
    )?;
    for (name, length) in references {
        writer.write_all(b"##contig=<ID=")?;
        writer.write_all(name)?;
        writeln!(writer, ",length={length}>")?;
    }
    writer.write_all(
        b"##FILTER=<ID=LowAD,Description=\"Bisulfite-informative alternate depth is below --min-alt-count\">\n",
    )?;
    writer.write_all(
        b"##FILTER=<ID=LowGQ,Description=\"Factorized ALT-identity and dosage genotype quality is below --min-gq\">\n",
    )?;
    writer.write_all(
        b"##FILTER=<ID=LowAQ,Description=\"Posterior ALT-presence quality is below --min-aq\">\n",
    )?;
    writer.write_all(
        b"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"High-quality fragment depth\">\n",
    )?;
    writer.write_all(
        b"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Prior-free expected ALT fraction conditional on the selected ALT set\">\n",
    )?;
    writer.write_all(
        b"##INFO=<ID=BSI,Number=1,Type=String,Description=\"Bisulfite discrimination uses BOTH strands or ONE unaffected strand\">\n",
    )?;
    writer.write_all(
        b"##INFO=<ID=BS8,Number=8,Type=Integer,Description=\"Top A,C,G,T then bottom A,C,G,T counts after fragment overlap collapse\">\n",
    )?;
    writer
        .write_all(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Unphased genotype\">\n")?;
    writer.write_all(
        b"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Phred-scaled confidence in posterior ALT identity and prior-free conditional dosage\">\n",
    )?;
    writer.write_all(
        b"##FORMAT=<ID=AQ,Number=A,Type=Integer,Description=\"Phred-scaled posterior quality that each selected ALT is present in the genotype\">\n",
    )?;
    writer.write_all(
        b"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"High-quality fragment depth\">\n",
    )?;
    writer.write_all(
        b"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Raw exact-base depths for REF and ALT alleles\">\n",
    )?;
    writer.write_all(
        b"##FORMAT=<ID=IAD,Number=R,Type=Integer,Description=\"Bisulfite-informative exact-base depths for REF and ALT alleles\">\n",
    )?;
    writer.write_all(
        b"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized Phred-scaled genotype likelihoods\">\n",
    )?;
    writer.write_all(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")?;
    writer.write_all(sample_name)?;
    writer.write_all(b"\n")
}

fn render_vcf_model_header(writer: &mut (impl Write + ?Sized)) -> io::Result<()> {
    writer.write_all(b"##bsbit_model=bisulfite-diploid-bayesian-v5\n")?;
    writer.write_all(
        b"##bsbit_genotype_decision=site-posterior-alt-plus-maximum-likelihood-dosage\n",
    )?;
    writer.write_all(b"##bsbit_site_posterior=reference-centered-heterozygosity-prior\n")
}

pub(crate) fn render_vcf_call(
    writer: &mut (impl Write + ?Sized),
    reference_name: &[u8],
    call: &VariantCall,
) -> io::Result<()> {
    let (alternate_storage, alternate_count) = call.alternates();
    let alternates = &alternate_storage[..alternate_count];
    debug_assert!(!alternates.is_empty());
    let (allele_storage, allele_count) = selected_alleles(call.reference, alternates);
    let alleles = &allele_storage[..allele_count];
    let (genotype_left, genotype_right) = call.genotype_indices(alternates);
    let raw_depth_storage =
        allele_storage.map(|allele| total_allele_depth(call.strand_counts, allele));
    let raw_depths = &raw_depth_storage[..allele_count];
    let informative_depth_storage =
        allele_storage.map(|allele| informative_allele_depth(call.strand_counts, allele, alleles));
    let informative_depths = &informative_depth_storage[..allele_count];
    writer.write_all(reference_name)?;
    write!(writer, "\t{}\t.\t", u64::from(call.position) + 1)?;
    writer.write_all(&[call.reference.ascii()])?;
    writer.write_all(b"\t")?;
    for (index, alternate) in alternates.iter().copied().enumerate() {
        if index != 0 {
            writer.write_all(b",")?;
        }
        writer.write_all(&[alternate.ascii()])?;
    }
    write!(writer, "\t{:.2}\t", call.quality)?;
    write_filters(writer, call.filters)?;
    write!(writer, "\tDP={};AF=", call.depth)?;
    write_floats(
        writer,
        &call.conditional_alternate_frequencies[..alternate_count],
    )?;
    let discrimination = if has_conversion_confounded_pair(alleles) {
        "ONE"
    } else {
        "BOTH"
    };
    write!(writer, ";BSI={discrimination};BS8=")?;
    for strand in 0..2 {
        for base in 0..4 {
            if strand != 0 || base != 0 {
                writer.write_all(b",")?;
            }
            write!(writer, "{}", call.strand_counts[strand][base])?;
        }
    }
    write!(
        writer,
        "\tGT:GQ:AQ:DP:AD:IAD:PL\t{genotype_left}/{genotype_right}:{}:",
        call.genotype_quality
    )?;
    write_values(writer, &call.allele_qualities[..alternate_count])?;
    write!(writer, ":{}:", call.depth)?;
    write_values(writer, raw_depths)?;
    writer.write_all(b":")?;
    write_values(writer, informative_depths)?;
    writer.write_all(b":")?;
    let genotype_count = allele_count * (allele_count + 1) / 2;
    write_values(writer, &call.phred_likelihoods[..genotype_count])?;
    writer.write_all(b"\n")
}

fn write_values<T: std::fmt::Display>(
    writer: &mut (impl Write + ?Sized),
    values: &[T],
) -> io::Result<()> {
    for (index, value) in values.iter().enumerate() {
        if index != 0 {
            writer.write_all(b",")?;
        }
        write!(writer, "{value}")?;
    }
    Ok(())
}

fn write_floats(writer: &mut (impl Write + ?Sized), values: &[f64]) -> io::Result<()> {
    for (index, value) in values.iter().enumerate() {
        if index != 0 {
            writer.write_all(b",")?;
        }
        write!(writer, "{value:.6}")?;
    }
    Ok(())
}

fn write_filters(writer: &mut (impl Write + ?Sized), filters: u8) -> io::Result<()> {
    if filters == 0 {
        return writer.write_all(b"PASS");
    }
    let mut wrote_filter = false;
    for (mask, name) in [
        (FILTER_LOW_ALTERNATE_DEPTH, b"LowAD".as_slice()),
        (FILTER_LOW_GENOTYPE_QUALITY, b"LowGQ".as_slice()),
        (FILTER_LOW_ALLELE_QUALITY, b"LowAQ".as_slice()),
    ] {
        if filters & mask != 0 {
            if wrote_filter {
                writer.write_all(b";")?;
            }
            writer.write_all(name)?;
            wrote_filter = true;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::snp::result::{Base, Genotype};

    #[test]
    fn vcf_header_declares_posterior_and_evidence_contract() {
        let mut output = Vec::new();
        render_vcf_header(
            &mut output,
            &[(b"chr1".as_slice(), 100)],
            SnpConfig::default(),
            b"tumor",
        )
        .unwrap();
        let header = String::from_utf8(output).unwrap();
        assert!(header.contains("##bsbit_model=bisulfite-diploid-bayesian-v5\n"));
        assert!(header.contains(
            "##bsbit_genotype_decision=site-posterior-alt-plus-maximum-likelihood-dosage\n"
        ));
        assert!(
            header.contains("##bsbit_site_posterior=reference-centered-heterozygosity-prior\n")
        );
        assert!(header.contains(
            "##bsbit_methylation_marginalization=adaptive-log-concave-simpson-uniform-prior\n"
        ));
        assert!(header.contains("##bsbit_heterozygosity=0.001000000\n"));
        assert!(header.contains("##bsbit_min_alternate_fraction=0.100000000\n"));
        assert!(header.contains("##bsbit_min_genotype_quality=0\n"));
        assert!(header.contains("##bsbit_min_allele_quality=30\n"));
        assert!(header.ends_with("\tFORMAT\ttumor\n"));
        for field in ["AF", "BSI", "BS8", "AQ", "AD", "IAD", "PL"] {
            assert!(
                header.contains(&format!("ID={field},")),
                "{field}: {header}"
            );
        }
    }

    #[test]
    fn vcf_renderer_exposes_eight_strand_specific_counts() {
        let call = VariantCall {
            position: 9,
            reference: Base::A,
            genotype: Genotype {
                left: Base::A,
                right: Base::G,
            },
            depth: 8,
            genotype_quality: 42,
            allele_qualities: [37, 0],
            quality: 50.25,
            conditional_alternate_frequencies: [0.5, 0.0],
            phred_likelihoods: [50, 0, 50, 0, 0, 0],
            strand_counts: [[1, 2, 3, 4], [5, 6, 7, 8]],
            filters: 0,
        };
        let mut output = Vec::new();
        render_vcf_call(&mut output, b"chr1", &call).unwrap();
        assert_eq!(
            output,
            b"chr1\t10\t.\tA\tG\t50.25\tPASS\tDP=8;AF=0.500000;BSI=ONE;BS8=1,2,3,4,5,6,7,8\tGT:GQ:AQ:DP:AD:IAD:PL\t0/1:42:37:8:6,10:1,10:50,0,50\n"
        );
    }
}
