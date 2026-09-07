//! Frozen constants and shared arithmetic for the current combined-index image.
//!
//! Keeping the reader and writer on one private layout definition prevents a
//! format change on one side from silently diverging from the other.

pub(crate) const BWT_WORDS_PER_128_ROWS: u64 = 5;
pub(crate) const SA_FLAG_WORDS_PER_256_ROWS: u64 = 5;

pub(crate) const META_BYTES: usize = 120;
pub(crate) const META_BYTES_U32: u32 = 120;
pub(crate) const META_EXTENSION_MAGIC: &[u8; 8] = b"BSBICMB1";
pub(crate) const META_EXTENSION_MAJOR: u16 = 1;
pub(crate) const META_EXTENSION_MINOR: u16 = 0;
pub(crate) const META_EXTENSION_MINOR_SA8: u16 = 1;
pub(crate) const META_EXTENSION_OFFSET: usize = 68;
pub(crate) const META_DIGEST_OFFSET: usize = 84;

pub(crate) const LOOKUP_BASES: usize = 16;
pub(crate) const LOOKUP_KEYS: u64 = 43_046_721;
#[cfg(feature = "index-construction")]
pub(crate) const LOOKUP_KEYS_USIZE: usize = 43_046_721;
pub(crate) const LOOKUP_ENTRIES: u64 = LOOKUP_KEYS + 1;
#[cfg(feature = "index-construction")]
pub(crate) const LOOKUP_GAP_BITS: u32 = 4;
#[cfg(feature = "index-construction")]
pub(crate) const LOOKUP_BOUNDARY_HIGH_MASK: u64 = 0x0fff_ffff;

#[cfg(all(test, feature = "index-construction"))]
pub(crate) const SA_STRIDE: u64 = 16;
pub(crate) const SA_VALUE_BITS: u32 = 30;
pub(crate) const SA_VALUE_MASK: u64 = (1_u64 << SA_VALUE_BITS) - 1;
pub(crate) const OCC_STRIDE: u32 = 64;
pub(crate) const HIGH_OCC_STRIDE: u32 = 128;

/// Returns all three LF boundaries from one validated packed-rank boundary.
///
/// Storage validation remains with the combined image reader; the builder
/// shares this arithmetic so encoding and decoding cannot silently diverge.
#[inline]
#[cfg(feature = "index-construction")]
pub(crate) fn lf_all_boundaries(
    boundary: u64,
    suffix_count: u64,
    sentinel_row: u64,
    first_occurrence: [u64; 4],
    mut bwt_word: impl FnMut(u64) -> u64,
    mut high_occ: impl FnMut(u64) -> u64,
) -> Option<[u64; 3]> {
    if boundary > suffix_count {
        return None;
    }
    let line = boundary - u64::from(boundary > sentinel_row);
    let high_word = (line >> 7).checked_mul(BWT_WORDS_PER_128_ROWS)?;
    let low_block = (line & 127) >> 6;
    let plane_start = high_word.checked_add(1 + (low_block << 1))?;
    let high_occ_block = (line >> 16).checked_mul(2)?;
    let counter_word = bwt_word(high_word);
    let first_plane = bwt_word(plane_start);
    let second_plane = bwt_word(plane_start + 1);
    let first_absolute = high_occ(high_occ_block);
    let second_absolute = high_occ(high_occ_block + 1);
    let counter_shift = low_block << 5;
    let packed = counter_word >> (32 - counter_shift);
    let nonzero = ((packed >> 16) & 0xffff) + (packed & 0xffff);
    let at_block = [
        ((line >> 6) << 6).checked_sub(
            first_absolute
                .checked_add(second_absolute)?
                .checked_add(nonzero)?,
        )?,
        first_absolute + ((counter_word >> (48 - counter_shift)) & 0xffff),
        second_absolute + ((counter_word >> (32 - counter_shift)) & 0xffff),
    ];
    let need = u32::try_from(line & 63).expect("six bits fit u32");
    let within = if need == 0 {
        [0_u64; 3]
    } else {
        let shift = 64 - need;
        [
            u64::from(((!(first_plane | second_plane)) >> shift).count_ones()),
            u64::from((first_plane >> shift).count_ones()),
            u64::from((second_plane >> shift).count_ones()),
        ]
    };
    Some([
        first_occurrence[0]
            .checked_add(at_block[0])?
            .checked_add(within[0])?,
        first_occurrence[1]
            .checked_add(at_block[1])?
            .checked_add(within[1])?,
        first_occurrence[2]
            .checked_add(at_block[2])?
            .checked_add(within[2])?,
    ])
}
