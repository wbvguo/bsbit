//! Exact quality likelihoods and bisulfite-aware diploid SNP decisions.

use super::candidate::CandidateSite;
use super::result::{
    Base, FILTER_LOW_ALLELE_QUALITY, FILTER_LOW_ALTERNATE_DEPTH, FILTER_LOW_GENOTYPE_QUALITY,
    GENOTYPES, Genotype, SnpConfig, VariantCall, alternate_alleles, filtered_observation,
    genotype_index, has_exact_alternate_set, informative_allele_depth, selected_alleles,
    validate_config,
};
use crate::CallError;
use crate::evidence::fragment::combined_observation_error;
use crate::evidence::{EvidenceObservation, EvidenceStrand};

const OBSERVATION_HISTOGRAM_FLUSH: usize = 1_024;
const METHYLATION_MODE_ITERATIONS: usize = 64;
const METHYLATION_INTEGRATION_TOLERANCE: f64 = 1e-12;
const METHYLATION_INTEGRATION_MAX_DEPTH: u8 = 48;
const METHYLATION_INTEGRATION_MAX_EVALUATIONS: usize = 16_384;
const LOG10_SCALE: f64 = 10.0 / std::f64::consts::LN_10;

struct LikelihoodSite {
    reference: Base,
    // Evidence that is independent of methylation is accumulated once per genotype.
    // Conversion-sensitive evidence is retained as compact encoded observations and
    // run-length encoded once it becomes deep. This permits stable adaptive
    // integration without retaining one floating-point array per grid node.
    constant_log_likelihoods: [f64; 10],
    conversion_observations: ObservationHistogram,
    depth: u32,
    strand_counts: [[u32; 4]; 2],
}

impl LikelihoodSite {
    fn new(reference: Base) -> Self {
        Self {
            reference,
            constant_log_likelihoods: [0.0; 10],
            conversion_observations: ObservationHistogram::default(),
            depth: 0,
            strand_counts: [[0; 4]; 2],
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct ObservationBin {
    encoded: u16,
    count: u32,
}

#[derive(Debug, Default)]
struct ObservationHistogram {
    pending: Vec<u16>,
    bins: Vec<ObservationBin>,
}

impl ObservationHistogram {
    fn observe(
        &mut self,
        observed: Base,
        strand: EvidenceStrand,
        base_quality: u8,
        mapping_quality: u8,
    ) -> Result<(), CallError> {
        let encoded = encode_observation(observed, strand, base_quality, mapping_quality);
        if self.bins.is_empty() && self.pending.len() < OBSERVATION_HISTOGRAM_FLUSH {
            self.pending.try_reserve(1).map_err(|error| {
                CallError::with_source(
                    crate::CallErrorKind::Calling,
                    "reserve SNP conversion-observation histogram",
                    error,
                )
            })?;
            self.pending.push(encoded);
            return Ok(());
        }
        if !self.pending.is_empty() {
            self.flush_pending()?;
        }
        self.increment_bin(encoded)
    }

    fn into_bins(mut self) -> Result<Vec<ObservationBin>, CallError> {
        if !self.pending.is_empty() {
            self.flush_pending()?;
        }
        Ok(self.bins)
    }

    fn flush_pending(&mut self) -> Result<(), CallError> {
        self.pending.sort_unstable();
        let pending = std::mem::take(&mut self.pending);
        self.bins.try_reserve(pending.len()).map_err(|error| {
            CallError::with_source(
                crate::CallErrorKind::Calling,
                "reserve compressed SNP conversion-observation histogram",
                error,
            )
        })?;
        for encoded in pending {
            self.increment_bin(encoded)?;
        }
        Ok(())
    }

    fn increment_bin(&mut self, encoded: u16) -> Result<(), CallError> {
        match self
            .bins
            .binary_search_by_key(&encoded, |candidate| candidate.encoded)
        {
            Ok(index) => {
                self.bins[index].count =
                    self.bins[index].count.checked_add(1).ok_or_else(|| {
                        CallError::operation("SNP conversion-observation count overflowed u32")
                    })?;
            }
            Err(index) => {
                self.bins.try_reserve(1).map_err(|error| {
                    CallError::with_source(
                        crate::CallErrorKind::Calling,
                        "grow compressed SNP conversion-observation histogram",
                        error,
                    )
                })?;
                self.bins
                    .insert(index, ObservationBin { encoded, count: 1 });
            }
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct AffineLogFactor {
    intercept: f64,
    slope: f64,
    count: u32,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct LikelihoodModel {
    config: SnpConfig,
    genotype_log_priors: [[f64; 10]; 4],
}

impl LikelihoodModel {
    fn new(config: SnpConfig) -> Self {
        let genotype_log_priors =
            Base::ALL.map(|reference| genotype_log_priors(reference, config.heterozygosity_rate));
        Self {
            config,
            genotype_log_priors,
        }
    }
}

/// Exact quality-aware likelihood state allocated only for first-pass sites.
pub(crate) struct LikelihoodRegion {
    start: u32,
    site_by_offset: Vec<u16>,
    sites: Vec<(u32, LikelihoodSite)>,
    model: LikelihoodModel,
}

pub(crate) const fn likelihood_site_bytes() -> usize {
    std::mem::size_of::<(u32, LikelihoodSite)>()
        + OBSERVATION_HISTOGRAM_FLUSH * std::mem::size_of::<u16>()
}

impl LikelihoodRegion {
    pub(crate) fn new(candidates: &[CandidateSite], config: SnpConfig) -> Result<Self, CallError> {
        validate_config(config)?;
        let model = LikelihoodModel::new(config);
        let Some(first) = candidates.first() else {
            return Ok(Self {
                start: 0,
                site_by_offset: Vec::new(),
                sites: Vec::new(),
                model,
            });
        };
        if candidates.len() >= usize::from(u16::MAX) {
            return Err(CallError::operation(
                "SNP likelihood batch has too many candidate sites",
            ));
        }
        let last = candidates
            .last()
            .expect("nonempty likelihood candidates have a last element");
        let span = last
            .position
            .checked_sub(first.position)
            .and_then(|span| span.checked_add(1))
            .ok_or_else(|| {
                CallError::operation("SNP likelihood candidates are not position-sorted")
            })?;
        let span = usize::try_from(span)
            .map_err(|_| CallError::operation("SNP likelihood span is not addressable"))?;
        let mut site_by_offset = Vec::new();
        site_by_offset.try_reserve_exact(span).map_err(|error| {
            CallError::with_source(
                crate::CallErrorKind::Calling,
                format!("allocate {span} SNP likelihood lookup entries"),
                error,
            )
        })?;
        site_by_offset.resize(span, u16::MAX);
        let mut sites = Vec::new();
        sites.try_reserve_exact(candidates.len()).map_err(|error| {
            CallError::with_source(
                crate::CallErrorKind::Calling,
                format!("allocate {} SNP likelihood sites", candidates.len()),
                error,
            )
        })?;
        for (index, candidate) in candidates.iter().enumerate() {
            if let Some(previous) = index.checked_sub(1).and_then(|index| candidates.get(index))
                && candidate.position < previous.position
            {
                return Err(CallError::operation(
                    "SNP likelihood candidates are not position-sorted",
                ));
            }
            let offset = candidate
                .position
                .checked_sub(first.position)
                .and_then(|offset| usize::try_from(offset).ok())
                .filter(|offset| *offset < site_by_offset.len())
                .ok_or_else(|| {
                    CallError::operation("SNP likelihood candidate escaped its lookup span")
                })?;
            let encoded = u16::try_from(index)
                .map_err(|_| CallError::operation("SNP likelihood site index exceeds u16"))?;
            if site_by_offset[offset] != u16::MAX {
                return Err(CallError::operation(format!(
                    "duplicate SNP candidate at position {}",
                    candidate.position
                )));
            }
            site_by_offset[offset] = encoded;
            sites.push((candidate.position, LikelihoodSite::new(candidate.reference)));
        }
        Ok(Self {
            start: first.position,
            site_by_offset,
            sites,
            model,
        })
    }

    pub(crate) fn observe_fragment(
        &mut self,
        observations: &[EvidenceObservation],
    ) -> Result<(), CallError> {
        for observation in observations {
            let Some(offset) = observation
                .position
                .checked_sub(self.start)
                .and_then(|offset| usize::try_from(offset).ok())
            else {
                continue;
            };
            let Some(&encoded) = self.site_by_offset.get(offset) else {
                continue;
            };
            if encoded == u16::MAX {
                continue;
            }
            let site = &mut self.sites[usize::from(encoded)].1;
            let Some((reference, observed, base_quality, mapping_quality)) =
                filtered_observation(*observation, self.model.config)
            else {
                continue;
            };
            if reference != site.reference {
                return Err(CallError::operation(format!(
                    "SNP likelihood reference changed at position {}",
                    observation.position
                )));
            }
            site.depth = site
                .depth
                .checked_add(1)
                .ok_or_else(|| CallError::operation("SNP likelihood depth overflowed u32"))?;
            let strand = match observation.strand {
                EvidenceStrand::Top => 0,
                EvidenceStrand::Bottom => 1,
            };
            site.strand_counts[strand][observed.index()] = site.strand_counts[strand]
                [observed.index()]
            .checked_add(1)
            .ok_or_else(|| CallError::operation("SNP strand/base count overflowed u32"))?;
            let error = combined_observation_error(base_quality, mapping_quality);
            for (genotype_index, genotype) in GENOTYPES.iter().copied().enumerate() {
                if genotype.is_conversion_sensitive_on(observation.strand) {
                    continue;
                }
                let probability = genotype_observation_probability(
                    genotype,
                    observed,
                    observation.strand,
                    error,
                    0.5,
                );
                site.constant_log_likelihoods[genotype_index] += probability.ln();
            }
            site.conversion_observations.observe(
                observed,
                observation.strand,
                base_quality,
                mapping_quality,
            )?;
        }
        Ok(())
    }

    pub(crate) fn calls(self) -> Result<Vec<VariantCall>, CallError> {
        let mut calls = Vec::new();
        for (position, site) in self.sites {
            let Some(call) = call_site(position, site, &self.model)? else {
                continue;
            };
            calls.try_reserve(1).map_err(|error| {
                CallError::with_source(
                    crate::CallErrorKind::Calling,
                    "reserve SNP variant result",
                    error,
                )
            })?;
            calls.push(call);
        }
        Ok(calls)
    }
}

fn call_site(
    position: u32,
    site: LikelihoodSite,
    model: &LikelihoodModel,
) -> Result<Option<VariantCall>, CallError> {
    if site.depth < model.config.minimum_depth {
        return Ok(None);
    }
    let conversion_observations = site.conversion_observations.into_bins()?;
    let mut likelihoods = site.constant_log_likelihoods;
    for (index, genotype) in GENOTYPES.iter().copied().enumerate() {
        if genotype.is_conversion_sensitive() {
            likelihoods[index] +=
                marginalized_genotype_log_likelihood(genotype, &conversion_observations, model)?;
        }
    }
    let log_posteriors: [f64; 10] = std::array::from_fn(|index| {
        likelihoods[index] + model.genotype_log_priors[site.reference.index()][index]
    });
    let normalizer = log_sum_exp(&log_posteriors);
    let posteriors: [f64; 10] =
        log_posteriors.map(|log_posterior| (log_posterior - normalizer).exp());
    let mut best_posterior_index = 0;
    for index in 1..GENOTYPES.len() {
        if posteriors[index] > posteriors[best_posterior_index] {
            best_posterior_index = index;
        }
    }
    let posterior_genotype = GENOTYPES[best_posterior_index];
    if posterior_genotype.left == site.reference && posterior_genotype.right == site.reference {
        return Ok(None);
    }

    let reference_index = genotype_index(site.reference, site.reference);
    let quality = phred_error_probability(posteriors[reference_index], 999.0);
    let (alternate_storage, alternate_count) =
        alternate_alleles(site.reference, posterior_genotype);
    let alternates = &alternate_storage[..alternate_count];

    let (genotype, dosage_confidence, conditional_alternate_frequencies) =
        maximum_likelihood_dosage(site.reference, posterior_genotype, alternates, &likelihoods);
    let alternate_set_probability = GENOTYPES
        .iter()
        .copied()
        .zip(posteriors)
        .filter(|(candidate, _)| has_exact_alternate_set(site.reference, *candidate, alternates))
        .map(|(_, posterior)| posterior)
        .sum::<f64>();
    let genotype_quality = bounded_rounded_u8(
        phred_error_probability(1.0 - alternate_set_probability * dosage_confidence, 99.0),
        99,
    );

    let mut allele_qualities = [0; 2];
    for (alternate_index, alternate) in alternates.iter().copied().enumerate() {
        let absence_probability = GENOTYPES
            .iter()
            .copied()
            .zip(posteriors)
            .filter(|(candidate, _)| !candidate.contains(alternate))
            .map(|(_, posterior)| posterior)
            .sum();
        allele_qualities[alternate_index] =
            bounded_rounded_u8(phred_error_probability(absence_probability, 99.0), 99);
    }
    let (allele_storage, allele_count) = selected_alleles(site.reference, alternates);
    let phred_likelihoods =
        selected_phred_likelihoods(&allele_storage[..allele_count], &likelihoods);
    let mut filters = 0;
    if alternates.iter().copied().any(|alternate| {
        informative_allele_depth(
            site.strand_counts,
            alternate,
            &allele_storage[..allele_count],
        ) < model.config.minimum_alternate_count
    }) {
        filters |= FILTER_LOW_ALTERNATE_DEPTH;
    }
    if genotype_quality < model.config.minimum_genotype_quality {
        filters |= FILTER_LOW_GENOTYPE_QUALITY;
    }
    if allele_qualities[..alternate_count]
        .iter()
        .any(|quality| *quality < model.config.minimum_allele_quality)
    {
        filters |= FILTER_LOW_ALLELE_QUALITY;
    }
    Ok(Some(VariantCall {
        position,
        reference: site.reference,
        genotype,
        depth: site.depth,
        genotype_quality,
        allele_qualities,
        quality,
        conditional_alternate_frequencies,
        phred_likelihoods,
        strand_counts: site.strand_counts,
        filters,
    }))
}

// Preserve the site posterior's selected ALT set, then estimate copy number
// without the rare-site prior. One ALT compares REF/ALT with ALT/ALT; two
// selected ALTs have one possible diploid dosage.
fn maximum_likelihood_dosage(
    reference: Base,
    posterior_genotype: Genotype,
    alternates: &[Base],
    likelihoods: &[f64; 10],
) -> (Genotype, f64, [f64; 2]) {
    let [alternate] = alternates else {
        debug_assert_eq!(alternates.len(), 2);
        return (posterior_genotype, 1.0, [0.5, 0.5]);
    };
    let heterozygous_index = genotype_index(reference, *alternate);
    let homozygous_alternate_index = genotype_index(*alternate, *alternate);
    let dosage_normalizer = log_sum_exp_pair(
        likelihoods[heterozygous_index],
        likelihoods[homozygous_alternate_index],
    );
    let heterozygous_probability = (likelihoods[heterozygous_index] - dosage_normalizer).exp();
    let homozygous_alternate_probability =
        (likelihoods[homozygous_alternate_index] - dosage_normalizer).exp();
    let (genotype, confidence) = if homozygous_alternate_probability > heterozygous_probability {
        (
            GENOTYPES[homozygous_alternate_index],
            homozygous_alternate_probability,
        )
    } else {
        (GENOTYPES[heterozygous_index], heterozygous_probability)
    };
    (
        genotype,
        confidence,
        [
            heterozygous_probability * 0.5 + homozygous_alternate_probability,
            0.0,
        ],
    )
}

fn genotype_log_priors(reference: Base, heterozygosity: f64) -> [f64; 10] {
    let reference_frequency = (1.0 - heterozygosity).sqrt();
    let alternate_frequency = (1.0 - reference_frequency) / 3.0;
    GENOTYPES.map(|genotype| {
        let left = if genotype.left == reference {
            reference_frequency
        } else {
            alternate_frequency
        };
        let right = if genotype.right == reference {
            reference_frequency
        } else {
            alternate_frequency
        };
        let probability = if genotype.left == genotype.right {
            left * right
        } else {
            2.0 * left * right
        };
        probability.ln()
    })
}

fn marginalized_genotype_log_likelihood(
    genotype: Genotype,
    observations: &[ObservationBin],
    model: &LikelihoodModel,
) -> Result<f64, CallError> {
    let mut factors = Vec::new();
    factors.try_reserve(observations.len()).map_err(|error| {
        CallError::with_source(
            crate::CallErrorKind::Calling,
            "reserve adaptive methylation likelihood factors",
            error,
        )
    })?;
    for observation in observations {
        let (observed, strand, base_quality, mapping_quality) =
            decode_observation(observation.encoded)?;
        if !genotype.is_conversion_sensitive_on(strand) {
            continue;
        }
        let error = combined_observation_error(base_quality, mapping_quality);
        let intercept = genotype_observation_probability(
            genotype,
            observed,
            strand,
            error,
            model.config.underconversion_rate,
        );
        let at_one = genotype_observation_probability(
            genotype,
            observed,
            strand,
            error,
            1.0 - model.config.overconversion_rate,
        );
        factors.push(AffineLogFactor {
            intercept,
            slope: at_one - intercept,
            count: observation.count,
        });
    }
    adaptive_log_integral(&factors)
}

fn adaptive_log_integral(factors: &[AffineLogFactor]) -> Result<f64, CallError> {
    if factors.is_empty() {
        return Ok(0.0);
    }
    let mode = likelihood_mode(factors);
    let maximum = affine_log_likelihood(factors, mode);
    if !maximum.is_finite() {
        return Err(CallError::operation(
            "adaptive methylation likelihood has no finite mode",
        ));
    }
    let mut evaluations = 0_usize;
    let mut evaluate = |methylation: f64| -> Result<f64, CallError> {
        evaluations = evaluations.checked_add(1).ok_or_else(|| {
            CallError::operation("adaptive methylation integration evaluation count overflowed")
        })?;
        if evaluations > METHYLATION_INTEGRATION_MAX_EVALUATIONS {
            return Err(CallError::operation(
                "adaptive methylation integration exceeded its evaluation budget",
            ));
        }
        Ok((affine_log_likelihood(factors, methylation) - maximum).exp())
    };
    let integral = if mode <= f64::EPSILON || mode >= 1.0 - f64::EPSILON {
        adaptive_simpson_interval(0.0, 1.0, METHYLATION_INTEGRATION_TOLERANCE, &mut evaluate)?
    } else {
        adaptive_simpson_interval(
            0.0,
            mode,
            METHYLATION_INTEGRATION_TOLERANCE * 0.5,
            &mut evaluate,
        )? + adaptive_simpson_interval(
            mode,
            1.0,
            METHYLATION_INTEGRATION_TOLERANCE * 0.5,
            &mut evaluate,
        )?
    };
    if integral <= 0.0 || !integral.is_finite() {
        return Err(CallError::operation(
            "adaptive methylation integration produced a nonpositive result",
        ));
    }
    Ok(maximum + integral.ln())
}

fn likelihood_mode(factors: &[AffineLogFactor]) -> f64 {
    let derivative_at_zero = affine_log_likelihood_derivative(factors, 0.0);
    if derivative_at_zero <= 0.0 {
        return 0.0;
    }
    let derivative_at_one = affine_log_likelihood_derivative(factors, 1.0);
    if derivative_at_one >= 0.0 {
        return 1.0;
    }
    let mut lower = 0.0_f64;
    let mut upper = 1.0_f64;
    for _ in 0..METHYLATION_MODE_ITERATIONS {
        let middle = (lower + upper) * 0.5;
        if affine_log_likelihood_derivative(factors, middle) > 0.0 {
            lower = middle;
        } else {
            upper = middle;
        }
    }
    (lower + upper) * 0.5
}

fn affine_log_likelihood(factors: &[AffineLogFactor], methylation: f64) -> f64 {
    factors
        .iter()
        .map(|factor| {
            let probability = factor.intercept + factor.slope * methylation;
            if probability <= 0.0 {
                f64::NEG_INFINITY
            } else {
                f64::from(factor.count) * probability.ln()
            }
        })
        .sum()
}

fn affine_log_likelihood_derivative(factors: &[AffineLogFactor], methylation: f64) -> f64 {
    factors
        .iter()
        .map(|factor| {
            let probability = factor.intercept + factor.slope * methylation;
            if probability <= 0.0 {
                if factor.slope > 0.0 {
                    f64::INFINITY
                } else if factor.slope < 0.0 {
                    f64::NEG_INFINITY
                } else {
                    0.0
                }
            } else {
                f64::from(factor.count) * factor.slope / probability
            }
        })
        .sum()
}

fn adaptive_simpson_interval(
    start: f64,
    end: f64,
    tolerance: f64,
    evaluate: &mut impl FnMut(f64) -> Result<f64, CallError>,
) -> Result<f64, CallError> {
    let middle = (start + end) * 0.5;
    let at_start = evaluate(start)?;
    let at_middle = evaluate(middle)?;
    let at_end = evaluate(end)?;
    let whole = simpson_estimate(start, end, at_start, at_middle, at_end);
    adaptive_simpson_refine(
        start,
        end,
        at_start,
        at_middle,
        at_end,
        whole,
        tolerance,
        METHYLATION_INTEGRATION_MAX_DEPTH,
        evaluate,
    )
}

#[allow(clippy::too_many_arguments)]
fn adaptive_simpson_refine(
    start: f64,
    end: f64,
    at_start: f64,
    at_middle: f64,
    at_end: f64,
    whole: f64,
    tolerance: f64,
    depth: u8,
    evaluate: &mut impl FnMut(f64) -> Result<f64, CallError>,
) -> Result<f64, CallError> {
    let middle = (start + end) * 0.5;
    let left_middle = (start + middle) * 0.5;
    let right_middle = (middle + end) * 0.5;
    let at_left_middle = evaluate(left_middle)?;
    let at_right_middle = evaluate(right_middle)?;
    let left = simpson_estimate(start, middle, at_start, at_left_middle, at_middle);
    let right = simpson_estimate(middle, end, at_middle, at_right_middle, at_end);
    let refined = left + right;
    let correction = refined - whole;
    if correction.abs() <= 15.0 * tolerance {
        return Ok(refined + correction / 15.0);
    }
    if depth == 0 {
        return Err(CallError::operation(
            "adaptive methylation integration did not converge",
        ));
    }
    Ok(adaptive_simpson_refine(
        start,
        middle,
        at_start,
        at_left_middle,
        at_middle,
        left,
        tolerance * 0.5,
        depth - 1,
        evaluate,
    )? + adaptive_simpson_refine(
        middle,
        end,
        at_middle,
        at_right_middle,
        at_end,
        right,
        tolerance * 0.5,
        depth - 1,
        evaluate,
    )?)
}

fn simpson_estimate(start: f64, end: f64, at_start: f64, at_middle: f64, at_end: f64) -> f64 {
    (end - start) * (at_start + 4.0 * at_middle + at_end) / 6.0
}

fn log_sum_exp(values: &[f64]) -> f64 {
    let maximum = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    maximum
        + values
            .iter()
            .map(|value| (*value - maximum).exp())
            .sum::<f64>()
            .ln()
}

fn log_sum_exp_pair(left: f64, right: f64) -> f64 {
    let maximum = left.max(right);
    maximum + ((left - maximum).exp() + (right - maximum).exp()).ln()
}

fn phred_error_probability(probability: f64, maximum: f64) -> f64 {
    if probability <= 0.0 {
        maximum
    } else {
        (-probability.log10() * 10.0).clamp(0.0, maximum)
    }
}

#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn bounded_rounded_u8(value: f64, maximum: u8) -> u8 {
    value.round().clamp(0.0, f64::from(maximum)) as u8
}

#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn bounded_rounded_u16(value: f64, maximum: u16) -> u16 {
    value.round().clamp(0.0, f64::from(maximum)) as u16
}

fn selected_phred_likelihoods(alleles: &[Base], likelihoods: &[f64; 10]) -> [u16; 6] {
    let mut selected = [f64::NEG_INFINITY; 6];
    let mut count = 0;
    for second in 0..alleles.len() {
        for first in 0..=second {
            selected[count] = likelihoods[genotype_index(alleles[first], alleles[second])];
            count += 1;
        }
    }
    let maximum = selected[..count]
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let mut phred = [0; 6];
    for (destination, likelihood) in phred[..count].iter_mut().zip(&selected[..count]) {
        *destination = bounded_rounded_u16((maximum - likelihood) * LOG10_SCALE, 999);
    }
    phred
}

fn encode_observation(
    observed: Base,
    strand: EvidenceStrand,
    base_quality: u8,
    mapping_quality: u8,
) -> u16 {
    let strand = match strand {
        EvidenceStrand::Top => 0_u16,
        EvidenceStrand::Bottom => 1,
    };
    u16::try_from(observed.index()).expect("four canonical bases fit two bits")
        | (strand << 2)
        | (u16::from(base_quality.min(60)) << 3)
        | (u16::from(mapping_quality.min(60)) << 9)
}

fn decode_observation(encoded: u16) -> Result<(Base, EvidenceStrand, u8, u8), CallError> {
    if encoded >> 15 != 0 {
        return Err(CallError::operation(
            "encoded SNP conversion observation exceeds fifteen bits",
        ));
    }
    let observed = Base::ALL
        .get(usize::from(encoded & 0x03))
        .copied()
        .ok_or_else(|| CallError::operation("encoded SNP observation base is invalid"))?;
    let strand = if encoded & (1 << 2) == 0 {
        EvidenceStrand::Top
    } else {
        EvidenceStrand::Bottom
    };
    let base_quality = u8::try_from((encoded >> 3) & 0x3f)
        .map_err(|_| CallError::operation("encoded SNP base quality exceeds u8"))?;
    let mapping_quality = u8::try_from((encoded >> 9) & 0x3f)
        .map_err(|_| CallError::operation("encoded SNP mapping quality exceeds u8"))?;
    if base_quality > 60 || mapping_quality > 60 {
        return Err(CallError::operation(
            "encoded SNP effective quality exceeds the model cap",
        ));
    }
    Ok((observed, strand, base_quality, mapping_quality))
}

fn genotype_observation_probability(
    genotype: Genotype,
    observed: Base,
    strand: EvidenceStrand,
    error: f64,
    retention: f64,
) -> f64 {
    let left = allele_observation_probability(genotype.left, observed, strand, error, retention);
    let right = allele_observation_probability(genotype.right, observed, strand, error, retention);
    ((left + right) * 0.5).max(f64::MIN_POSITIVE)
}

fn allele_observation_probability(
    allele: Base,
    observed: Base,
    strand: EvidenceStrand,
    error: f64,
    retention: f64,
) -> f64 {
    let mut truth = [0.0_f64; 4];
    match (allele, strand) {
        (Base::C, EvidenceStrand::Top) => {
            truth[Base::C.index()] = retention;
            truth[Base::T.index()] = 1.0 - retention;
        }
        (Base::G, EvidenceStrand::Bottom) => {
            truth[Base::G.index()] = retention;
            truth[Base::A.index()] = 1.0 - retention;
        }
        _ => truth[allele.index()] = 1.0,
    }
    truth
        .iter()
        .copied()
        .enumerate()
        .map(|(base, probability)| {
            probability
                * if base == observed.index() {
                    1.0 - error
                } else {
                    error / 3.0
                }
        })
        .sum::<f64>()
        .max(f64::MIN_POSITIVE)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::snp::vcf::render_vcf_call;

    fn observation(
        reference: u8,
        observed: u8,
        strand: EvidenceStrand,
        quality: u8,
    ) -> EvidenceObservation {
        EvidenceObservation {
            reference: 0,
            position: 10,
            reference_base: reference,
            query_base: Some(observed),
            base_quality: Some(quality),
            mapping_quality: 60,
            strand,
            context: None,
        }
    }
    #[test]
    fn quality_likelihood_calls_supported_nonconversion_snv() {
        let config = SnpConfig {
            minimum_depth: 4,
            minimum_alternate_count: 2,
            minimum_genotype_quality: 10,
            ..SnpConfig::default()
        };
        let candidates = [CandidateSite {
            position: 10,
            reference: Base::A,
        }];
        let mut likelihood = LikelihoodRegion::new(&candidates, config).unwrap();
        let observations = (0..24)
            .map(|_| observation(b'A', b'G', EvidenceStrand::Top, 40))
            .collect::<Vec<_>>();
        likelihood.observe_fragment(&observations).unwrap();
        let calls = likelihood.calls().unwrap();
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].reference, Base::A);
        assert_eq!(
            calls[0].genotype,
            Genotype {
                left: Base::G,
                right: Base::G
            }
        );
        assert!(calls[0].genotype_quality >= 10);
    }

    #[test]
    fn methylation_is_marginalized_and_prior_probabilities_normalize() {
        for (depth, methylated) in [(100_u32, 1_u32), (1_000, 330), (10_000, 100)] {
            let factors = [
                AffineLogFactor {
                    intercept: 0.0,
                    slope: 1.0,
                    count: methylated,
                },
                AffineLogFactor {
                    intercept: 1.0,
                    slope: -1.0,
                    count: depth - methylated,
                },
            ];
            let expected = log_factorial(methylated) + log_factorial(depth - methylated)
                - log_factorial(depth + 1);
            let actual = adaptive_log_integral(&factors).unwrap();
            assert!(
                (actual - expected).abs() < 1e-8,
                "depth {depth}, methylated {methylated}: {actual} versus {expected}"
            );
        }

        for reference in Base::ALL {
            let priors = genotype_log_priors(reference, 0.001);
            let total = priors.into_iter().map(f64::exp).sum::<f64>();
            assert!((total - 1.0).abs() < 1e-12);
            let reference_probability = priors[genotype_index(reference, reference)].exp();
            assert!((reference_probability - 0.999).abs() < 1e-12);
        }
    }

    #[test]
    fn adaptive_marginalization_resolves_a_deep_mode_between_fixed_grid_nodes() {
        let depth = 10_000_u32;
        let methylated = 3_300_u32;
        let factors = [
            AffineLogFactor {
                intercept: 0.0,
                slope: 1.0,
                count: methylated,
            },
            AffineLogFactor {
                intercept: 1.0,
                slope: -1.0,
                count: depth - methylated,
            },
        ];
        let expected = log_factorial(methylated) + log_factorial(depth - methylated)
            - log_factorial(depth + 1);
        let maximum = affine_log_likelihood(&factors, 0.33);
        let weighted_sum = (0..=20)
            .map(|index| {
                let methylation = f64::from(index) * 0.05;
                let weight = if index == 0 || index == 20 {
                    1.0
                } else if index % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                weight * (affine_log_likelihood(&factors, methylation) - maximum).exp()
            })
            .sum::<f64>();
        let fixed_grid = maximum + (0.05 / 3.0 * weighted_sum).ln();
        assert!((fixed_grid - expected).abs() > 5.0);
        assert!((adaptive_log_integral(&factors).unwrap() - expected).abs() < 1e-8);
    }

    #[test]
    fn genotype_likelihood_not_site_prior_determines_alt_dosage() {
        let heterozygous = genotype_index(Base::A, Base::G);
        let homozygous_alternate = genotype_index(Base::G, Base::G);

        for heterozygosity_rate in [0.0001, 0.001, 0.1] {
            let mut site = LikelihoodSite::new(Base::A);
            site.depth = 10;
            site.constant_log_likelihoods = [-100.0; GENOTYPES.len()];
            site.constant_log_likelihoods[heterozygous] = 0.0;
            site.constant_log_likelihoods[homozygous_alternate] = 10.0_f64.ln();
            let model = LikelihoodModel::new(SnpConfig {
                heterozygosity_rate,
                ..SnpConfig::default()
            });
            let call = call_site(10, site, &model)
                .unwrap()
                .expect("non-reference call");
            assert_eq!(
                call.genotype,
                Genotype {
                    left: Base::G,
                    right: Base::G,
                }
            );
        }
    }

    fn log_factorial(value: u32) -> f64 {
        (1..=value).map(|item| f64::from(item).ln()).sum()
    }

    #[test]
    fn deep_conversion_only_evidence_is_not_miscalled_as_c_to_t() {
        let config = SnpConfig {
            minimum_depth: 1,
            minimum_alternate_count: 1,
            minimum_genotype_quality: 0,
            ..SnpConfig::default()
        };
        let candidates = [CandidateSite {
            position: 10,
            reference: Base::C,
        }];
        let mut likelihood = LikelihoodRegion::new(&candidates, config).unwrap();
        let mut observations = Vec::with_capacity(1_000);
        observations.extend((0..500).map(|_| observation(b'C', b'C', EvidenceStrand::Top, 40)));
        observations.extend((0..500).map(|_| observation(b'C', b'T', EvidenceStrand::Top, 40)));
        likelihood.observe_fragment(&observations).unwrap();
        assert!(likelihood.calls().unwrap().is_empty());
    }

    #[test]
    fn unaffected_strand_support_resolves_a_c_to_t_variant() {
        let config = SnpConfig {
            minimum_depth: 4,
            minimum_alternate_count: 2,
            minimum_genotype_quality: 0,
            ..SnpConfig::default()
        };
        let candidates = [CandidateSite {
            position: 10,
            reference: Base::C,
        }];
        let mut likelihood = LikelihoodRegion::new(&candidates, config).unwrap();
        let mut observations = Vec::with_capacity(48);
        observations.extend((0..32).map(|_| observation(b'C', b'T', EvidenceStrand::Top, 40)));
        observations.extend((0..8).map(|_| observation(b'C', b'C', EvidenceStrand::Bottom, 40)));
        observations.extend((0..8).map(|_| observation(b'C', b'T', EvidenceStrand::Bottom, 40)));
        likelihood.observe_fragment(&observations).unwrap();
        let calls = likelihood.calls().unwrap();
        assert_eq!(calls.len(), 1);
        assert_eq!(
            calls[0].genotype,
            Genotype {
                left: Base::C,
                right: Base::T,
            }
        );
        assert_eq!(calls[0].filters & FILTER_LOW_ALTERNATE_DEPTH, 0);
    }

    #[test]
    fn genotype_uncertainty_filters_but_does_not_erase_a_variant() {
        let config = SnpConfig {
            minimum_depth: 1,
            minimum_alternate_count: 1,
            minimum_genotype_quality: 99,
            minimum_allele_quality: 20,
            ..SnpConfig::default()
        };
        let candidates = [CandidateSite {
            position: 10,
            reference: Base::A,
        }];
        let mut likelihood = LikelihoodRegion::new(&candidates, config).unwrap();
        let observations = (0..4)
            .map(|_| observation(b'A', b'G', EvidenceStrand::Top, 40))
            .collect::<Vec<_>>();
        likelihood.observe_fragment(&observations).unwrap();
        let calls = likelihood.calls().unwrap();
        assert_eq!(calls.len(), 1);
        assert_ne!(calls[0].filters & FILTER_LOW_GENOTYPE_QUALITY, 0);
        assert_eq!(calls[0].filters & FILTER_LOW_ALTERNATE_DEPTH, 0);
        assert_eq!(calls[0].filters & FILTER_LOW_ALLELE_QUALITY, 0);
        let mut output = Vec::new();
        render_vcf_call(&mut output, b"chr1", &calls[0]).unwrap();
        assert!(String::from_utf8(output).unwrap().contains("\tLowGQ\t"));
    }
}
