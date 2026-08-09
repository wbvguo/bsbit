//! Bounded-memory construction of the combined three-letter combined-index BWT.
//!
//! The complete projected text is retained at two bits per symbol.  A
//! libsais32 suffix sort initializes the rightmost block; preceding blocks are
//! inserted with exact FM-gap ranks.  Global rows are `usize`, while block
//! offsets and SA16 values remain narrow.  BWT rows and row-ordered SA16
//! quotients are merged in place, so construction never materializes a global
//! 64-bit suffix array or a second complete BWT.

use core::cmp::Ordering;
use core::ffi::c_int;
use core::fmt;
use core::ptr;
use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};
use std::thread;

use crate::reference::ContigInput;
use bsbit_core::alphabet::Base;

use crate::build::libsais::libsais_omp;

const SYMBOL_MASK: u8 = 3;
const SAMPLE_FLAG: u8 = 4;
const SENTINEL_CODE: u8 = 3;
const DEFAULT_SAMPLE_STRIDE: usize = 16;
const LOCAL_RANK_STRIDE: usize = 64;
const SUPER_RANK_STRIDE: usize = 65_536;
const SAMPLE_RANK_STRIDE: usize = 256;
const RADIX_BITS: u32 = 16;
const RADIX_BUCKETS: usize = 1 << RADIX_BITS;
const MIN_PARALLEL_RADIX_ROWS: usize = 1 << 18;
const PROJECTION_MIN_BYTES_PER_WORKER: usize = 1 << 20;
const MIN_MEMORY_HEADROOM: u64 = 256 << 20;
#[cfg(not(test))]
const PARALLEL_MERGE_CHUNK_ROWS: usize = 128 * 1024 * 1024;
// Keep test chunks small enough to exercise high-to-low chunk boundaries,
// including boundaries that split one packed output byte.
#[cfg(test)]
const PARALLEL_MERGE_CHUNK_ROWS: usize = 64;
#[cfg(not(test))]
const PARALLEL_MERGE_MIN_ROWS: usize = 1 << 20;
#[cfg(test)]
const PARALLEL_MERGE_MIN_ROWS: usize = 32;

/// Resource limits for the bounded constructor.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct BoundedBwtConfig {
    memory_mib: u64,
    threads: u32,
    block_bases: Option<usize>,
    sample_stride: usize,
}

impl BoundedBwtConfig {
    pub(crate) fn new(memory_mib: u64, threads: u32) -> Result<Self, BoundedBwtError> {
        if memory_mib == 0 {
            return Err(BoundedBwtError::InvalidConfiguration(
                "memory MiB must be positive",
            ));
        }
        if threads == 0 || c_int::try_from(threads).is_err() {
            return Err(BoundedBwtError::InvalidConfiguration(
                "thread count must fit a positive c_int",
            ));
        }
        Ok(Self {
            memory_mib,
            threads,
            block_bases: None,
            sample_stride: DEFAULT_SAMPLE_STRIDE,
        })
    }

    pub(crate) fn with_sample_stride(
        mut self,
        sample_stride: usize,
    ) -> Result<Self, BoundedBwtError> {
        if !matches!(sample_stride, 8 | 16) {
            return Err(BoundedBwtError::InvalidConfiguration(
                "sample stride must be 8 or 16",
            ));
        }
        self.sample_stride = sample_stride;
        Ok(self)
    }

    #[cfg(test)]
    pub(crate) fn with_block_bases(mut self, block_bases: usize) -> Result<Self, BoundedBwtError> {
        if block_bases == 0 || c_int::try_from(block_bases).is_err() {
            return Err(BoundedBwtError::InvalidConfiguration(
                "block size must fit a positive c_int",
            ));
        }
        self.block_bases = Some(block_bases);
        Ok(self)
    }

    fn block_bases(self, text_len: usize) -> Result<usize, BoundedBwtError> {
        if let Some(block_bases) = self.block_bases {
            return Ok(block_bases.min(text_len));
        }
        let budget = self
            .memory_mib
            .checked_mul(1 << 20)
            .ok_or(BoundedBwtError::SizeOverflow)?;
        let text_bytes =
            u64::try_from(text_len.div_ceil(4)).map_err(|_| BoundedBwtError::SizeOverflow)?;
        let row_bytes =
            u64::try_from((text_len + 1).div_ceil(2)).map_err(|_| BoundedBwtError::SizeOverflow)?;
        let sample_bytes = u64::try_from(text_len / self.sample_stride + 1)
            .map_err(|_| BoundedBwtError::SizeOverflow)?
            .checked_mul(4)
            .ok_or(BoundedBwtError::SizeOverflow)?;
        let persistent = text_bytes
            .checked_add(row_bytes)
            .and_then(|value| value.checked_add(sample_bytes))
            .ok_or(BoundedBwtError::SizeOverflow)?;
        // During ranking, one u64 key per block row overlaps the compact rank
        // index.  During sorting, two u64 key arrays overlap the persistent
        // state.  The latter is the tighter bound; allocator/I/O overhead gets
        // at least 256 MiB or 5% of the requested budget.
        let headroom = (budget / 20).max(MIN_MEMORY_HEADROOM);
        let available = budget
            .checked_sub(persistent)
            .and_then(|value| value.checked_sub(headroom))
            .ok_or(BoundedBwtError::InvalidConfiguration(
                "memory budget cannot hold the packed text and final BWT",
            ))?;
        let radix_rows = available / 16;
        // Local-rank data uses about 5/64 bytes per existing row.  This bound
        // is deliberately rounded up to one byte per eight rows.
        let rank_bytes =
            u64::try_from(text_len.div_ceil(8)).map_err(|_| BoundedBwtError::SizeOverflow)?;
        let rank_rows = budget
            .checked_sub(persistent)
            .and_then(|value| value.checked_sub(rank_bytes))
            .and_then(|value| value.checked_sub(headroom))
            .map_or(0, |value| value / 8);
        let rows = radix_rows.min(rank_rows);
        let maximum = u64::try_from(c_int::MAX).expect("c_int maximum fits u64");
        let rows = rows.min(maximum).min(text_len as u64);
        usize::try_from(rows.max(1)).map_err(|_| BoundedBwtError::SizeOverflow)
    }
}

/// Construction failure for the bounded backend.
#[derive(Debug)]
pub(crate) enum BoundedBwtError {
    InvalidConfiguration(&'static str),
    InvalidInput(&'static str),
    SizeOverflow,
    Allocation(&'static str),
    NativeStatus(i32),
    Invariant(&'static str),
    WorkerPanic,
}

impl fmt::Display for BoundedBwtError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidConfiguration(message) => {
                write!(formatter, "invalid configuration: {message}")
            }
            Self::InvalidInput(message) => write!(formatter, "invalid input: {message}"),
            Self::SizeOverflow => formatter.write_str("bounded combined-index BWT size overflow"),
            Self::Allocation(label) => write!(formatter, "allocation failed for {label}"),
            Self::NativeStatus(status) => write!(formatter, "libsais32 returned status {status}"),
            Self::Invariant(message) => {
                write!(
                    formatter,
                    "bounded combined-index BWT invariant failed: {message}"
                )
            }
            Self::WorkerPanic => formatter.write_str("bounded combined-index BWT worker panicked"),
        }
    }
}

impl std::error::Error for BoundedBwtError {}

#[derive(Clone, Copy)]
struct ReferenceSegment<'a> {
    start: usize,
    end: usize,
    bases: &'a [Base],
}

/// Complete combined projection at two bits per G/T/A symbol.
#[derive(Debug)]
pub(crate) struct PackedProjectedText {
    bytes: Vec<u8>,
    len: usize,
    reference_bases: u64,
}

impl PackedProjectedText {
    #[cfg(test)]
    pub(crate) fn from_digits(digits: &[u8]) -> Result<Self, BoundedBwtError> {
        Self::from_projected_digits(digits.to_vec(), (digits.len() / 2) as u64)
    }

    #[cfg(test)]
    pub(crate) fn from_projected_digits(
        digits: Vec<u8>,
        reference_bases: u64,
    ) -> Result<Self, BoundedBwtError> {
        if digits.is_empty() || digits.iter().any(|&digit| digit > 2) {
            return Err(BoundedBwtError::InvalidInput(
                "packed projected text accepts nonempty G/T/A digits",
            ));
        }
        let len = digits.len();
        let mut bytes = try_vec_filled(len.div_ceil(4), 0_u8, "packed projected text")?;
        for (position, &digit) in digits.iter().enumerate() {
            bytes[position / 4] |= digit << (2 * (position % 4));
        }
        drop(digits);
        Ok(Self {
            bytes,
            len,
            reference_bases,
        })
    }

    #[inline]
    pub(crate) fn get(&self, position: usize) -> u8 {
        debug_assert!(position < self.len);
        (self.bytes[position / 4] >> (2 * (position % 4))) & 3
    }

    pub(crate) const fn len(&self) -> usize {
        self.len
    }

    pub(crate) const fn reference_bases(&self) -> u64 {
        self.reference_bases
    }

    pub(crate) fn decode_range(
        &self,
        start: usize,
        end: usize,
    ) -> Result<Vec<u8>, BoundedBwtError> {
        if start > end || end > self.len {
            return Err(BoundedBwtError::Invariant(
                "decode range exceeds packed text",
            ));
        }
        let mut decoded = try_vec_filled(end - start, 0_u8, "decoded libsais32 block")?;
        for (offset, digit) in decoded.iter_mut().enumerate() {
            *digit = self.get(start + offset);
        }
        Ok(decoded)
    }

    fn compare_suffixes(&self, left: usize, right: usize) -> Ordering {
        if left == right {
            return Ordering::Equal;
        }
        let available = (self.len - left).min(self.len - right);
        let mut offset = 0_usize;
        while offset + 32 <= available {
            let left_word = self.symbol_word(left + offset);
            let right_word = self.symbol_word(right + offset);
            let difference = left_word ^ right_word;
            if difference != 0 {
                let symbol = usize::try_from(difference.trailing_zeros() / 2)
                    .expect("u64 symbol index fits usize");
                return self
                    .get(left + offset + symbol)
                    .cmp(&self.get(right + offset + symbol));
            }
            offset += 32;
        }
        while offset < available {
            match self.get(left + offset).cmp(&self.get(right + offset)) {
                Ordering::Equal => offset += 1,
                ordering => return ordering,
            }
        }
        (self.len - left).cmp(&(self.len - right))
    }

    #[inline]
    fn symbol_word(&self, position: usize) -> u64 {
        debug_assert!(position + 32 <= self.len);
        let byte = position / 4;
        let shift = u32::try_from(2 * (position % 4)).expect("packed shift fits u32");
        // SAFETY: 32 available symbols cover eight complete bytes from this
        // byte. A nonzero intra-byte shift implies that one spill byte also
        // exists in the packed allocation.
        let low = unsafe { ptr::read_unaligned(self.bytes.as_ptr().add(byte).cast::<u64>()) };
        if shift == 0 {
            low
        } else {
            let high = u64::from(self.bytes[byte + 8]);
            (low >> shift) | (high << (64 - shift))
        }
    }
}

/// Projects a consumed catalog directly into the two-bit combined text.
#[allow(clippy::too_many_lines)]
pub(crate) fn project_combined_packed_text(
    contigs: &[ContigInput],
    projection_salt: u64,
    threads: u32,
) -> Result<PackedProjectedText, BoundedBwtError> {
    if threads == 0 {
        return Err(BoundedBwtError::InvalidConfiguration(
            "projection thread count must be positive",
        ));
    }
    let reference_bases = contigs.iter().try_fold(0_usize, |total, contig| {
        let length =
            usize::try_from(contig.sequence().len()).map_err(|_| BoundedBwtError::SizeOverflow)?;
        total
            .checked_add(length)
            .ok_or(BoundedBwtError::SizeOverflow)
    })?;
    if reference_bases == 0 {
        return Err(BoundedBwtError::InvalidInput("reference catalog is empty"));
    }
    let len = reference_bases
        .checked_mul(2)
        .ok_or(BoundedBwtError::SizeOverflow)?;
    let mut segments = Vec::new();
    try_reserve(&mut segments, contigs.len(), "projection segment catalog")?;
    let mut start = 0_usize;
    for contig in contigs {
        let bases = contig.sequence().bases();
        let end = start
            .checked_add(bases.len())
            .ok_or(BoundedBwtError::SizeOverflow)?;
        if !bases.is_empty() {
            segments.push(ReferenceSegment { start, end, bases });
        }
        start = end;
    }
    if start != reference_bases || segments.is_empty() {
        return Err(BoundedBwtError::Invariant(
            "projection segment dimensions disagree",
        ));
    }

    let packed_len = len.div_ceil(4);
    let mut bytes = try_vec_filled(packed_len, 0_u8, "packed combined projection")?;
    let requested_workers = usize::try_from(threads).map_err(|_| BoundedBwtError::SizeOverflow)?;
    let useful_workers = packed_len.div_ceil(PROJECTION_MIN_BYTES_PER_WORKER).max(1);
    let workers = requested_workers.min(useful_workers).min(packed_len);
    let next_worker = AtomicUsize::new(0);
    let output_address = bytes.as_mut_ptr().expose_provenance();
    thread::scope(|scope| {
        let mut handles = Vec::with_capacity(workers);
        for _ in 0..workers {
            handles.push(scope.spawn(|| {
                loop {
                    let worker = next_worker.fetch_add(1, AtomicOrdering::Relaxed);
                    if worker >= workers {
                        break;
                    }
                    let byte_start = packed_len * worker / workers;
                    let byte_end = packed_len * (worker + 1) / workers;
                    let position_start = byte_start * 4;
                    let position_end = (byte_end * 4).min(len);
                    let forward_end = position_end.min(reference_bases);
                    if position_start < forward_end {
                        let mut coordinate = position_start;
                        let mut segment = segments.partition_point(|item| item.end <= coordinate);
                        while coordinate < forward_end {
                            let item = segments.get(segment).ok_or(BoundedBwtError::Invariant(
                                "forward projection segment absent",
                            ))?;
                            let stop = forward_end.min(item.end);
                            while coordinate < stop {
                                let base = item.bases[coordinate - item.start];
                                let canonical =
                                    deterministic_canonical_code(base, coordinate, projection_salt);
                                let digit = [1, 0, 1, 2][usize::from(canonical)];
                                write_packed_digit(output_address, coordinate, digit);
                                coordinate += 1;
                            }
                            segment += 1;
                        }
                    }

                    let reverse_start = position_start.max(reference_bases);
                    if reverse_start < position_end {
                        let mut position = reverse_start;
                        let mut coordinate = len - 1 - position;
                        let mut segment =
                            segments.partition_point(|item| item.start <= coordinate) - 1;
                        while position < position_end {
                            let item = segments.get(segment).ok_or(BoundedBwtError::Invariant(
                                "reverse projection segment absent",
                            ))?;
                            loop {
                                let base = item.bases[coordinate - item.start];
                                let canonical =
                                    deterministic_canonical_code(base, coordinate, projection_salt);
                                let digit = [2, 1, 0, 1][usize::from(canonical)];
                                write_packed_digit(output_address, position, digit);
                                position += 1;
                                if position == position_end {
                                    break;
                                }
                                coordinate -= 1;
                                if coordinate < item.start {
                                    break;
                                }
                            }
                            if position < position_end {
                                segment = segment.checked_sub(1).ok_or(
                                    BoundedBwtError::Invariant("reverse projection underflow"),
                                )?;
                            }
                        }
                    }
                }
                Ok::<_, BoundedBwtError>(())
            }));
        }
        for handle in handles {
            handle.join().map_err(|_| BoundedBwtError::WorkerPanic)??;
        }
        Ok::<_, BoundedBwtError>(())
    })?;

    Ok(PackedProjectedText {
        bytes,
        len,
        reference_bases: u64::try_from(reference_bases)
            .map_err(|_| BoundedBwtError::SizeOverflow)?,
    })
}

#[inline]
fn write_packed_digit(address: usize, position: usize, digit: u8) {
    // SAFETY: projection workers own disjoint complete byte ranges. Every
    // destination starts at zero and receives each of its at most four fields
    // exactly once.
    unsafe {
        let output = ptr::with_exposed_provenance_mut::<u8>(address).add(position / 4);
        *output |= digit << (2 * (position % 4));
    }
}

fn deterministic_canonical_code(base: Base, coordinate: usize, projection_salt: u64) -> u8 {
    match base {
        Base::A => 0,
        Base::C => 1,
        Base::G => 2,
        Base::T => 3,
        _ => {
            let mut value = u64::try_from(coordinate)
                .unwrap_or(u64::MAX)
                .wrapping_add(projection_salt)
                .wrapping_add(0x9e37_79b9_7f4a_7c15);
            value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            u8::try_from((value ^ (value >> 31)) & 3).expect("two bits fit u8")
        }
    }
}

/// Final row-order BWT and SA16 quotient state.
#[derive(Debug)]
pub(crate) struct BoundedBwt {
    packed_rows: Vec<u8>,
    sample_quotients: Vec<u32>,
    rows: usize,
    sentinel_row: usize,
    sample_stride: usize,
}

impl BoundedBwt {
    pub(crate) const fn rows(&self) -> usize {
        self.rows
    }

    pub(crate) const fn text_len(&self) -> usize {
        self.rows - 1
    }

    pub(crate) const fn sentinel_row(&self) -> usize {
        self.sentinel_row
    }

    pub(crate) const fn sample_stride(&self) -> usize {
        self.sample_stride
    }

    #[inline]
    pub(crate) fn nibble(&self, row: usize) -> u8 {
        debug_assert!(row < self.rows);
        (self.packed_rows[row / 2] >> (4 * (row % 2))) & 15
    }

    #[inline]
    pub(crate) fn code(&self, row: usize) -> u8 {
        self.nibble(row) & SYMBOL_MASK
    }

    #[inline]
    pub(crate) fn has_sample(&self, row: usize) -> bool {
        self.nibble(row) & SAMPLE_FLAG != 0
    }

    #[inline]
    pub(crate) fn transformed_digit(&self, line: usize) -> u8 {
        debug_assert!(line < self.text_len());
        let row = line + usize::from(line >= self.sentinel_row);
        let code = self.code(row);
        debug_assert!(code < SENTINEL_CODE);
        code
    }

    pub(crate) fn sample_quotients(&self) -> &[u32] {
        &self.sample_quotients
    }

    pub(crate) fn row_ordered_samples(&self) -> impl Iterator<Item = (u64, u32)> + '_ {
        (0..self.rows)
            .filter(|&row| self.has_sample(row))
            .zip(self.sample_quotients.iter().copied())
            .map(|(row, quotient)| {
                (
                    u64::try_from(row).expect("validated combined row fits u64"),
                    quotient,
                )
            })
    }
}

#[derive(Debug)]
struct BwtState {
    packed_rows: Vec<u8>,
    sample_quotients: Vec<u32>,
    rows: usize,
    text_start: usize,
    text_end: usize,
    start_row: usize,
    counts: [u64; 3],
    rank: Option<RankIndex>,
    sample_stride: usize,
}

impl BwtState {
    fn from_last_block(
        text: &PackedProjectedText,
        start: usize,
        threads: u32,
        sample_stride: usize,
    ) -> Result<Self, BoundedBwtError> {
        let block = text.decode_range(start, text.len())?;
        let length = c_int::try_from(block.len()).map_err(|_| {
            BoundedBwtError::InvalidConfiguration("initial block exceeds libsais32")
        })?;
        let mut suffixes = try_vec_filled(block.len(), 0_i32, "initial libsais32 suffix array")?;
        // SAFETY: the decoded block and exact-length suffix output are live and
        // nonoverlapping for the complete native call.
        let status = unsafe {
            libsais_omp(
                block.as_ptr(),
                suffixes.as_mut_ptr(),
                length,
                0,
                ptr::null_mut(),
                c_int::try_from(threads).expect("validated thread count fits c_int"),
            )
        };
        if status != 0 {
            return Err(BoundedBwtError::NativeStatus(status));
        }

        let final_rows = text
            .len()
            .checked_add(1)
            .ok_or(BoundedBwtError::SizeOverflow)?;
        let rows = block
            .len()
            .checked_add(1)
            .ok_or(BoundedBwtError::SizeOverflow)?;
        let mut packed_rows = Vec::new();
        try_reserve(
            &mut packed_rows,
            final_rows.div_ceil(2),
            "reserved final nibble BWT",
        )?;
        let mut sample_quotients = Vec::new();
        try_reserve(
            &mut sample_quotients,
            text.len() / sample_stride + 1,
            "reserved final SA16 quotients",
        )?;
        let mut written_rows = 0_usize;
        push_row(
            &mut packed_rows,
            &mut sample_quotients,
            written_rows,
            text.get(text.len() - 1),
            sample_quotient(text.len(), sample_stride),
        );
        written_rows += 1;
        let mut start_row = None;
        for &local in &suffixes {
            let local = usize::try_from(local)
                .map_err(|_| BoundedBwtError::Invariant("negative libsais32 suffix"))?;
            if local >= block.len() {
                return Err(BoundedBwtError::Invariant(
                    "libsais32 suffix exceeds initial block",
                ));
            }
            let absolute = start
                .checked_add(local)
                .ok_or(BoundedBwtError::SizeOverflow)?;
            let code = if local == 0 {
                start_row = Some(written_rows);
                SENTINEL_CODE
            } else {
                block[local - 1]
            };
            push_row(
                &mut packed_rows,
                &mut sample_quotients,
                written_rows,
                code,
                sample_quotient(absolute, sample_stride),
            );
            written_rows += 1;
        }
        if written_rows != rows || packed_rows.len() != rows.div_ceil(2) {
            return Err(BoundedBwtError::Invariant(
                "initial packed BWT dimensions disagree",
            ));
        }
        let counts = count_digits(&block)?;
        Ok(Self {
            packed_rows,
            sample_quotients,
            rows,
            text_start: start,
            text_end: text.len(),
            start_row: start_row
                .ok_or(BoundedBwtError::Invariant("initial start suffix is absent"))?,
            counts,
            rank: None,
            sample_stride,
        })
    }

    #[inline]
    fn nibble(&self, row: usize) -> u8 {
        (self.packed_rows[row / 2] >> (4 * (row % 2))) & 15
    }

    #[inline]
    fn code(&self, row: usize) -> u8 {
        self.nibble(row) & SYMBOL_MASK
    }

    fn rebuild_rank_index(&mut self) -> Result<(), BoundedBwtError> {
        self.rank = Some(RankIndex::build(self)?);
        Ok(())
    }

    fn rank(&self) -> Result<&RankIndex, BoundedBwtError> {
        self.rank.as_ref().ok_or(BoundedBwtError::Invariant(
            "rank index is absent during FM ranking",
        ))
    }

    fn backward_rank(&self, code: u8, boundary: usize) -> Result<usize, BoundedBwtError> {
        if code >= SENTINEL_CODE || boundary > self.rows {
            return Err(BoundedBwtError::Invariant(
                "backward-rank input exceeds the FM domain",
            ));
        }
        let first = 1_u64
            + self.counts[..usize::from(code)]
                .iter()
                .copied()
                .sum::<u64>();
        usize::try_from(
            first
                .checked_add(self.rank()?.occ(self, code, boundary)?)
                .ok_or(BoundedBwtError::SizeOverflow)?,
        )
        .map_err(|_| BoundedBwtError::SizeOverflow)
    }

    fn compute_gap_keys(
        &self,
        text: &PackedProjectedText,
        block: &[u8],
        block_start: usize,
        block_end: usize,
        threads: u32,
        layout: KeyLayout,
    ) -> Result<Vec<u64>, BoundedBwtError> {
        if block_end != self.text_start
            || block_start >= block_end
            || block.len() != block_end - block_start
        {
            return Err(BoundedBwtError::Invariant(
                "prepended block is not adjacent to the current tail",
            ));
        }
        let length = block.len();
        let workers = usize::try_from(threads)
            .expect("validated thread count fits usize")
            .min(length);
        let workers = if self.text_end - self.text_start < self.sample_stride {
            1
        } else {
            workers
        };
        let mut boundaries = Vec::with_capacity(workers + 1);
        for worker in 0..=workers {
            boundaries.push(length * worker / workers);
        }
        let initial_ranks = thread::scope(|scope| {
            let mut handles = Vec::with_capacity(workers);
            for worker in 0..workers {
                let end = boundaries[worker + 1];
                handles.push(scope.spawn(move || {
                    if end == length {
                        Ok(self.start_row)
                    } else {
                        self.rank_of_external_suffix(text, block_start + end)
                    }
                }));
            }
            let mut ranks = Vec::with_capacity(workers);
            for handle in handles {
                ranks.push(handle.join().map_err(|_| BoundedBwtError::WorkerPanic)??);
            }
            Ok::<_, BoundedBwtError>(ranks)
        })?;

        let mut keys = try_vec_filled(length, 0_u64, "bounded gap/offset keys")?;
        thread::scope(|scope| {
            let mut remaining = keys.as_mut_slice();
            let mut handles = Vec::with_capacity(workers);
            for worker in 0..workers {
                let begin = boundaries[worker];
                let end = boundaries[worker + 1];
                let (chunk, rest) = remaining.split_at_mut(end - begin);
                remaining = rest;
                let mut rank = initial_ranks[worker];
                handles.push(scope.spawn(move || -> Result<(), BoundedBwtError> {
                    for local in (begin..end).rev() {
                        rank = self.backward_rank(block[local], rank)?;
                        chunk[local - begin] = layout.pack(rank, local)?;
                    }
                    Ok(())
                }));
            }
            for handle in handles {
                handle.join().map_err(|_| BoundedBwtError::WorkerPanic)??;
            }
            Ok::<_, BoundedBwtError>(())
        })?;
        Ok(keys)
    }

    fn rank_of_external_suffix(
        &self,
        text: &PackedProjectedText,
        query: usize,
    ) -> Result<usize, BoundedBwtError> {
        if query >= self.text_start {
            return Err(BoundedBwtError::Invariant(
                "external suffix does not precede the current tail",
            ));
        }
        let mut lower = 0_usize;
        let mut upper = self.rows;
        while lower < upper {
            let middle = lower + (upper - lower) / 2;
            let suffix = self.locate_row(middle)?;
            if text.compare_suffixes(suffix, query) == Ordering::Less {
                lower = middle + 1;
            } else {
                upper = middle;
            }
        }
        Ok(lower)
    }

    fn locate_row(&self, mut row: usize) -> Result<usize, BoundedBwtError> {
        let rank = self.rank()?;
        for steps in 0..=(2 * self.sample_stride - 2) {
            if self.nibble(row) & SAMPLE_FLAG != 0 {
                let quotient = usize::try_from(rank.sample_value(self, row)?)
                    .expect("u32 SA16 quotient fits usize");
                let sample = quotient
                    .checked_mul(self.sample_stride)
                    .ok_or(BoundedBwtError::SizeOverflow)?;
                let local_sample =
                    sample
                        .checked_sub(self.text_start)
                        .ok_or(BoundedBwtError::Invariant(
                            "partial-tail sample precedes its text",
                        ))?;
                let local = (local_sample + steps) % self.rows;
                return self
                    .text_start
                    .checked_add(local)
                    .ok_or(BoundedBwtError::SizeOverflow);
            }
            row = rank.lf(self, row)?;
        }
        Err(BoundedBwtError::Invariant(
            "sampled LF walk exceeded the partial-tail bound",
        ))
    }

    #[allow(clippy::too_many_arguments)]
    fn prepend_block_in_place(
        &mut self,
        text: &PackedProjectedText,
        block_start: usize,
        block_end: usize,
        keys: &[u64],
        layout: KeyLayout,
        added_counts: [u64; 3],
        threads: u32,
    ) -> Result<(), BoundedBwtError> {
        if threads > 1 && self.rows + keys.len() >= PARALLEL_MERGE_MIN_ROWS {
            self.prepend_block_parallel(
                text,
                block_start,
                block_end,
                keys,
                layout,
                added_counts,
                threads,
            )
        } else {
            self.prepend_block_bulk(text, block_start, block_end, keys, layout, added_counts)
        }
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn prepend_block_parallel(
        &mut self,
        text: &PackedProjectedText,
        block_start: usize,
        block_end: usize,
        keys: &[u64],
        layout: KeyLayout,
        added_counts: [u64; 3],
        threads: u32,
    ) -> Result<(), BoundedBwtError> {
        let block_len = block_end - block_start;
        let sample_stride = self.sample_stride;
        if block_end != self.text_start || block_len != keys.len() || self.rank.is_some() {
            return Err(BoundedBwtError::Invariant(
                "parallel merge dimensions or lifecycle disagree",
            ));
        }
        if self.code(self.start_row) != SENTINEL_CODE {
            return Err(BoundedBwtError::Invariant(
                "old start row does not carry the sentinel",
            ));
        }
        let old_rows = self.rows;
        let output_rows = old_rows
            .checked_add(block_len)
            .ok_or(BoundedBwtError::SizeOverflow)?;
        self.packed_rows.resize(output_rows.div_ceil(2), 0);
        let old_sample_count = self.sample_quotients.len();
        let output_sample_count = old_sample_count
            .checked_add(sampled_coordinates(block_start, block_end, sample_stride)?)
            .ok_or(BoundedBwtError::SizeOverflow)?;
        self.sample_quotients.resize(output_sample_count, 0);

        let workers = usize::try_from(threads)
            .expect("validated thread count fits usize")
            .max(1);
        let preceding = text.get(block_end - 1);
        let mut old_rows_buffer = Vec::new();
        try_reserve(
            &mut old_rows_buffer,
            PARALLEL_MERGE_CHUNK_ROWS.div_ceil(2),
            "parallel merge old-row buffer",
        )?;
        let mut old_samples_buffer = Vec::new();
        try_reserve(
            &mut old_samples_buffer,
            PARALLEL_MERGE_CHUNK_ROWS / sample_stride + 1,
            "parallel merge old-sample buffer",
        )?;
        let mut output_end = output_rows;
        let mut old_sample_end = old_sample_count;
        let mut output_sample_end = output_sample_count;
        let mut new_start_row = None;
        while output_end != 0 {
            let output_begin = output_end.saturating_sub(PARALLEL_MERGE_CHUNK_ROWS);
            let new_begin = keys.partition_point(|key| layout.output_position(*key) < output_begin);
            let new_end = keys.partition_point(|key| layout.output_position(*key) < output_end);
            let old_begin =
                output_begin
                    .checked_sub(new_begin)
                    .ok_or(BoundedBwtError::Invariant(
                        "parallel merge old-row prefix underflow",
                    ))?;
            let old_end = output_end
                .checked_sub(new_end)
                .ok_or(BoundedBwtError::Invariant(
                    "parallel merge old-row boundary underflow",
                ))?;
            let old_rows_in_chunk =
                old_end
                    .checked_sub(old_begin)
                    .ok_or(BoundedBwtError::Invariant(
                        "parallel merge old-row range regressed",
                    ))?;
            copy_nibbles_to_zero(
                &self.packed_rows,
                old_begin,
                old_rows_in_chunk,
                &mut old_rows_buffer,
            );

            let chunk_rows = output_end - output_begin;
            let worker_count = workers.min(chunk_rows.div_ceil(2).max(1));
            let mut boundaries = Vec::with_capacity(worker_count + 1);
            boundaries.push(output_begin);
            for worker in 1..worker_count {
                let raw = output_begin + chunk_rows * worker / worker_count;
                let boundary = raw.saturating_add(1) & !1_usize;
                if boundary < output_end
                    && boundary > *boundaries.last().expect("initial merge boundary exists")
                {
                    boundaries.push(boundary);
                }
            }
            if *boundaries.last().expect("initial merge boundary exists") != output_end {
                boundaries.push(output_end);
            }

            let mut parts = Vec::with_capacity(boundaries.len() - 1);
            let mut old_sample_cursor = 0_usize;
            let mut chunk_sample_count = 0_usize;
            for window in boundaries.windows(2) {
                let begin = window[0];
                let end = window[1];
                let key_begin = new_begin
                    + keys[new_begin..new_end]
                        .partition_point(|key| layout.output_position(*key) < begin);
                let key_end = new_begin
                    + keys[new_begin..new_end]
                        .partition_point(|key| layout.output_position(*key) < end);
                let old_local_begin = (begin - output_begin)
                    .checked_sub(key_begin - new_begin)
                    .ok_or(BoundedBwtError::Invariant(
                        "parallel merge partition old prefix underflow",
                    ))?;
                let old_local_end = (end - output_begin)
                    .checked_sub(key_end - new_begin)
                    .ok_or(BoundedBwtError::Invariant(
                        "parallel merge partition old boundary underflow",
                    ))?;
                let old_samples = usize::try_from(count_sample_range(
                    &old_rows_buffer,
                    old_local_begin,
                    old_local_end - old_local_begin,
                ))
                .map_err(|_| BoundedBwtError::SizeOverflow)?;
                let new_samples = count_new_key_samples(
                    &keys[key_begin..key_end],
                    layout,
                    block_start,
                    sample_stride,
                );
                let samples = old_samples
                    .checked_add(new_samples)
                    .ok_or(BoundedBwtError::SizeOverflow)?;
                parts.push(ParallelMergePart {
                    output_begin: begin,
                    output_end: end,
                    key_begin,
                    key_end,
                    old_begin: old_local_begin,
                    old_end: old_local_end,
                    old_sample_begin: old_sample_cursor,
                    old_sample_end: old_sample_cursor + old_samples,
                    output_sample_begin: 0,
                    output_sample_end: samples,
                });
                old_sample_cursor += old_samples;
                chunk_sample_count = chunk_sample_count
                    .checked_add(samples)
                    .ok_or(BoundedBwtError::SizeOverflow)?;
            }
            let old_sample_begin =
                old_sample_end
                    .checked_sub(old_sample_cursor)
                    .ok_or(BoundedBwtError::Invariant(
                        "parallel merge old samples exceed remaining input",
                    ))?;
            old_samples_buffer.clear();
            old_samples_buffer
                .extend_from_slice(&self.sample_quotients[old_sample_begin..old_sample_end]);
            let output_sample_begin = output_sample_end.checked_sub(chunk_sample_count).ok_or(
                BoundedBwtError::Invariant("parallel merge samples exceed remaining output"),
            )?;
            let mut sample_cursor = output_sample_begin;
            for part in &mut parts {
                let count = part.output_sample_end;
                part.output_sample_begin = sample_cursor;
                sample_cursor += count;
                part.output_sample_end = sample_cursor;
            }
            if sample_cursor != output_sample_end || old_sample_cursor != old_samples_buffer.len() {
                return Err(BoundedBwtError::Invariant(
                    "parallel merge sample partition accounting mismatch",
                ));
            }

            let packed_address = self.packed_rows.as_mut_ptr().expose_provenance();
            let sample_address = self.sample_quotients.as_mut_ptr().expose_provenance();
            let old_start_row = self.start_row;
            let old_rows_buffer = &old_rows_buffer;
            let old_samples_buffer = &old_samples_buffer;
            let starts = thread::scope(|scope| {
                let mut handles = Vec::with_capacity(parts.len());
                for part in parts {
                    handles.push(scope.spawn(move || {
                        fill_parallel_merge_part(
                            text,
                            block_start,
                            keys,
                            layout,
                            old_begin,
                            old_start_row,
                            preceding,
                            old_rows_buffer,
                            old_samples_buffer,
                            part,
                            packed_address,
                            sample_address,
                            sample_stride,
                        )
                    }));
                }
                let mut starts = Vec::new();
                for handle in handles {
                    if let Some(row) = handle.join().map_err(|_| BoundedBwtError::WorkerPanic)?? {
                        starts.push(row);
                    }
                }
                Ok::<_, BoundedBwtError>(starts)
            })?;
            for row in starts {
                if new_start_row.replace(row).is_some() {
                    return Err(BoundedBwtError::Invariant(
                        "parallel merge produced multiple new start rows",
                    ));
                }
            }
            output_end = output_begin;
            old_sample_end = old_sample_begin;
            output_sample_end = output_sample_begin;
        }
        if old_sample_end != 0 || output_sample_end != 0 {
            return Err(BoundedBwtError::Invariant(
                "parallel merge did not consume every sample",
            ));
        }
        for (total, added) in self.counts.iter_mut().zip(added_counts) {
            *total = total
                .checked_add(added)
                .ok_or(BoundedBwtError::SizeOverflow)?;
        }
        self.rows = output_rows;
        self.text_start = block_start;
        self.start_row = new_start_row.ok_or(BoundedBwtError::Invariant(
            "parallel merge omitted the prepended start suffix",
        ))?;
        Ok(())
    }

    #[allow(clippy::too_many_lines)]
    fn prepend_block_bulk(
        &mut self,
        text: &PackedProjectedText,
        block_start: usize,
        block_end: usize,
        keys: &[u64],
        layout: KeyLayout,
        added_counts: [u64; 3],
    ) -> Result<(), BoundedBwtError> {
        let block_len = block_end - block_start;
        if block_end != self.text_start || block_len != keys.len() || self.rank.is_some() {
            return Err(BoundedBwtError::Invariant(
                "in-place merge dimensions or lifecycle disagree",
            ));
        }
        let old_rows = self.rows;
        let output_rows = old_rows
            .checked_add(block_len)
            .ok_or(BoundedBwtError::SizeOverflow)?;
        self.packed_rows.resize(output_rows.div_ceil(2), 0);

        let old_sample_count = self.sample_quotients.len();
        let added_samples = sampled_coordinates(block_start, block_end, self.sample_stride)?;
        let output_sample_count = old_sample_count
            .checked_add(added_samples)
            .ok_or(BoundedBwtError::SizeOverflow)?;
        self.sample_quotients.resize(output_sample_count, 0);

        if self.code(self.start_row) != SENTINEL_CODE {
            return Err(BoundedBwtError::Invariant(
                "old start row does not carry the sentinel",
            ));
        }
        let preceding = text.get(block_end - 1);
        let mut output_cursor = output_rows;
        let mut old_cursor = old_rows;
        let mut old_sample = old_sample_count;
        let mut output_sample = output_sample_count;
        let mut new_start_row = None;
        let mut mapped_old_start_row = None;
        for new_ordinal in (0..keys.len()).rev() {
            let output_row = layout.output_position(keys[new_ordinal]);
            if output_row >= output_cursor {
                return Err(BoundedBwtError::Invariant(
                    "new-row positions are not strictly increasing",
                ));
            }
            let old_run = output_cursor - output_row - 1;
            let old_begin = old_cursor
                .checked_sub(old_run)
                .ok_or(BoundedBwtError::Invariant(
                    "old-row run exceeds remaining input",
                ))?;
            move_old_run_toward_end(
                &mut self.packed_rows,
                &mut self.sample_quotients,
                old_begin,
                output_row + 1,
                old_run,
                &mut old_sample,
                &mut output_sample,
                self.start_row,
                &mut mapped_old_start_row,
            )?;
            old_cursor = old_begin;
            output_cursor = output_row;

            let local = layout.local(keys[new_ordinal]);
            if local >= block_len {
                return Err(BoundedBwtError::Invariant(
                    "new suffix offset exceeds its block",
                ));
            }
            let absolute = block_start + local;
            let code = if local == 0 {
                if new_start_row.replace(output_row).is_some() {
                    return Err(BoundedBwtError::Invariant(
                        "multiple new start rows were produced",
                    ));
                }
                SENTINEL_CODE
            } else {
                text.get(absolute - 1)
            };
            let nibble = if let Some(quotient) = sample_quotient(absolute, self.sample_stride) {
                output_sample = output_sample
                    .checked_sub(1)
                    .ok_or(BoundedBwtError::Invariant("new sample underflow"))?;
                self.sample_quotients[output_sample] = quotient;
                code | SAMPLE_FLAG
            } else {
                code
            };
            set_nibble(&mut self.packed_rows, output_row, nibble);
        }
        if output_cursor != old_cursor {
            return Err(BoundedBwtError::Invariant(
                "old-row prefix differs from remaining output prefix",
            ));
        }
        move_old_run_toward_end(
            &mut self.packed_rows,
            &mut self.sample_quotients,
            0,
            0,
            old_cursor,
            &mut old_sample,
            &mut output_sample,
            self.start_row,
            &mut mapped_old_start_row,
        )?;
        let mapped_old_start_row = mapped_old_start_row.ok_or(BoundedBwtError::Invariant(
            "old sentinel row was omitted during merge",
        ))?;
        let old_start_nibble =
            (self.packed_rows[mapped_old_start_row / 2] >> (4 * (mapped_old_start_row % 2))) & 15;
        if old_start_nibble & SYMBOL_MASK != SENTINEL_CODE {
            return Err(BoundedBwtError::Invariant(
                "moved old sentinel row changed symbol",
            ));
        }
        set_nibble(
            &mut self.packed_rows,
            mapped_old_start_row,
            (old_start_nibble & SAMPLE_FLAG) | preceding,
        );
        if old_sample != 0
            || output_sample != 0
            || self.sample_quotients.len() != output_sample_count
        {
            return Err(BoundedBwtError::Invariant(
                "in-place merge did not consume every row and sample",
            ));
        }
        for (total, added) in self.counts.iter_mut().zip(added_counts) {
            *total = total
                .checked_add(added)
                .ok_or(BoundedBwtError::SizeOverflow)?;
        }
        self.rows = output_rows;
        self.text_start = block_start;
        self.start_row = new_start_row.ok_or(BoundedBwtError::Invariant(
            "prepended block start suffix is absent",
        ))?;
        Ok(())
    }

    fn finish(self, expected_text_len: usize) -> Result<BoundedBwt, BoundedBwtError> {
        if self.text_start != 0
            || self.text_end != expected_text_len
            || self.rows != expected_text_len + 1
            || self.code(self.start_row) != SENTINEL_CODE
            || self.sample_quotients.len() != expected_text_len / self.sample_stride + 1
        {
            return Err(BoundedBwtError::Invariant(
                "final BWT dimensions, sentinel, or SA16 count disagree",
            ));
        }
        Ok(BoundedBwt {
            packed_rows: self.packed_rows,
            sample_quotients: self.sample_quotients,
            rows: self.rows,
            sentinel_row: self.start_row,
            sample_stride: self.sample_stride,
        })
    }
}

/// Builds a complete exact BWT while honoring the configured working budget.
pub(crate) fn build_bounded_bwt(
    text: PackedProjectedText,
    config: BoundedBwtConfig,
) -> Result<BoundedBwt, BoundedBwtError> {
    if text.len() == 0 {
        return Err(BoundedBwtError::InvalidInput(
            "projected text must be nonempty",
        ));
    }
    let block_bases = config.block_bases(text.len())?;
    let first_start = text.len() - text.len().min(block_bases);
    let mut state =
        BwtState::from_last_block(&text, first_start, config.threads, config.sample_stride)?;
    if first_start != 0 {
        state.rebuild_rank_index()?;
    }

    let mut block_end = first_start;
    while block_end != 0 {
        let block_start = block_end.saturating_sub(block_bases);
        let block = text.decode_range(block_start, block_end)?;
        let added_counts = count_digits(&block)?;
        let layout = KeyLayout::new(block.len(), state.rows)?;
        let mut keys = state.compute_gap_keys(
            &text,
            &block,
            block_start,
            block_end,
            config.threads,
            layout,
        )?;
        drop(block);
        state.rank = None;

        parallel_radix_sort_gaps(&mut keys, config.threads, layout, state.rows)?;
        sort_equal_gap_suffixes(&mut keys, &text, block_start, config.threads, layout)?;
        gaps_to_output_positions(&mut keys, state.rows, layout)?;

        state.prepend_block_in_place(
            &text,
            block_start,
            block_end,
            &keys,
            layout,
            added_counts,
            config.threads,
        )?;
        drop(keys);
        block_end = block_start;
        if block_end != 0 {
            state.rebuild_rank_index()?;
        }
    }
    let final_state = state.finish(text.len())?;
    drop(text);
    Ok(final_state)
}

#[derive(Clone, Copy, Debug)]
struct KeyLayout {
    local_bits: u32,
    local_mask: u64,
}

impl KeyLayout {
    fn new(block_len: usize, old_rows: usize) -> Result<Self, BoundedBwtError> {
        if block_len == 0 {
            return Err(BoundedBwtError::Invariant("zero-length FM-gap block"));
        }
        let local_bits = (usize::BITS - (block_len - 1).leading_zeros()).max(1);
        let gap_bits = (usize::BITS - old_rows.leading_zeros()).max(1);
        if local_bits + gap_bits > u64::BITS {
            return Err(BoundedBwtError::SizeOverflow);
        }
        let local_mask = (1_u64 << local_bits) - 1;
        Ok(Self {
            local_bits,
            local_mask,
        })
    }

    fn pack(self, gap: usize, local: usize) -> Result<u64, BoundedBwtError> {
        let gap = u64::try_from(gap).map_err(|_| BoundedBwtError::SizeOverflow)?;
        let local = u64::try_from(local).map_err(|_| BoundedBwtError::SizeOverflow)?;
        if local > self.local_mask || gap > (u64::MAX >> self.local_bits) {
            return Err(BoundedBwtError::SizeOverflow);
        }
        Ok((gap << self.local_bits) | local)
    }

    #[inline]
    fn gap(self, key: u64) -> usize {
        usize::try_from(key >> self.local_bits).expect("validated global gap fits usize")
    }

    #[inline]
    fn local(self, key: u64) -> usize {
        usize::try_from(key & self.local_mask).expect("validated block offset fits usize")
    }

    #[inline]
    fn output_position(self, key: u64) -> usize {
        self.gap(key)
    }
}

#[derive(Debug)]
struct RankIndex {
    super_counts: Vec<[u64; 3]>,
    local_12: Vec<u32>,
    sample_counts: Vec<u32>,
}

impl RankIndex {
    fn build(state: &BwtState) -> Result<Self, BoundedBwtError> {
        let mut super_counts = Vec::new();
        try_reserve(
            &mut super_counts,
            state.rows / SUPER_RANK_STRIDE + 1,
            "FM super-rank checkpoints",
        )?;
        let mut local_12 = Vec::new();
        try_reserve(
            &mut local_12,
            state.rows / LOCAL_RANK_STRIDE + 1,
            "FM local-rank checkpoints",
        )?;
        let mut sample_counts = Vec::new();
        try_reserve(
            &mut sample_counts,
            state.rows / SAMPLE_RANK_STRIDE + 1,
            "FM sample-rank checkpoints",
        )?;

        let mut counts = [0_u64; 3];
        let mut super_base = [0_u64; 3];
        let mut samples = 0_u32;
        let mut first = 0_usize;
        loop {
            if first.is_multiple_of(SUPER_RANK_STRIDE) {
                super_base = counts;
                super_counts.push(counts);
            }
            let delta_1 = u16::try_from(counts[1] - super_base[1])
                .map_err(|_| BoundedBwtError::Invariant("digit-1 local rank exceeds u16"))?;
            let delta_2 = u16::try_from(counts[2] - super_base[2])
                .map_err(|_| BoundedBwtError::Invariant("digit-2 local rank exceeds u16"))?;
            local_12.push(u32::from(delta_1) | (u32::from(delta_2) << 16));
            if first.is_multiple_of(SAMPLE_RANK_STRIDE) {
                sample_counts.push(samples);
            }
            if first == state.rows {
                break;
            }
            let rows = (state.rows - first).min(LOCAL_RANK_STRIDE);
            for offset in (0..rows).step_by(16) {
                let width = (rows - offset).min(16);
                let word = packed_word(&state.packed_rows, first + offset, width);
                for code in 0..3_u8 {
                    counts[usize::from(code)] = counts[usize::from(code)]
                        .checked_add(u64::from(count_code_word(word, code, width)))
                        .ok_or(BoundedBwtError::SizeOverflow)?;
                }
                samples = samples
                    .checked_add(count_sample_word(word, width))
                    .ok_or(BoundedBwtError::SizeOverflow)?;
            }
            first = first
                .checked_add(rows)
                .ok_or(BoundedBwtError::SizeOverflow)?;
            if rows < LOCAL_RANK_STRIDE {
                break;
            }
        }
        if counts != state.counts
            || usize::try_from(samples).expect("u32 sample count fits usize")
                != state.sample_quotients.len()
        {
            return Err(BoundedBwtError::Invariant(
                "compact rank accounting disagrees with BWT state",
            ));
        }
        Ok(Self {
            super_counts,
            local_12,
            sample_counts,
        })
    }

    fn occ(&self, state: &BwtState, code: u8, boundary: usize) -> Result<u64, BoundedBwtError> {
        if boundary > state.rows || code >= SENTINEL_CODE {
            return Err(BoundedBwtError::Invariant(
                "compact occurrence query exceeds domain",
            ));
        }
        let local_first = boundary / LOCAL_RANK_STRIDE * LOCAL_RANK_STRIDE;
        let local = *self
            .local_12
            .get(boundary / LOCAL_RANK_STRIDE)
            .ok_or(BoundedBwtError::Invariant("local rank checkpoint absent"))?;
        let super_counts = *self
            .super_counts
            .get(boundary / SUPER_RANK_STRIDE)
            .ok_or(BoundedBwtError::Invariant("super rank checkpoint absent"))?;
        let count_1 = super_counts[1] + u64::from(local & 0xffff);
        let count_2 = super_counts[2] + u64::from(local >> 16);
        let base = match code {
            1 => count_1,
            2 => count_2,
            0 => {
                let sentinel = u64::from(state.start_row < local_first);
                u64::try_from(local_first)
                    .map_err(|_| BoundedBwtError::SizeOverflow)?
                    .checked_sub(sentinel)
                    .and_then(|value| value.checked_sub(count_1))
                    .and_then(|value| value.checked_sub(count_2))
                    .ok_or(BoundedBwtError::Invariant(
                        "inferred digit-0 rank underflow",
                    ))?
            }
            _ => unreachable!(),
        };
        Ok(base
            + u64::from(count_code_range(
                &state.packed_rows,
                local_first,
                boundary - local_first,
                code,
            )))
    }

    fn sample_value(&self, state: &BwtState, row: usize) -> Result<u32, BoundedBwtError> {
        if state.nibble(row) & SAMPLE_FLAG == 0 {
            return Err(BoundedBwtError::Invariant(
                "sample lookup requested for an unmarked row",
            ));
        }
        let checkpoint = row / SAMPLE_RANK_STRIDE;
        let first = checkpoint * SAMPLE_RANK_STRIDE;
        let ordinal = self.sample_counts[checkpoint]
            .checked_add(count_sample_range(&state.packed_rows, first, row - first))
            .ok_or(BoundedBwtError::SizeOverflow)?;
        state
            .sample_quotients
            .get(usize::try_from(ordinal).expect("u32 sample ordinal fits usize"))
            .copied()
            .ok_or(BoundedBwtError::Invariant(
                "sample rank exceeds stored quotients",
            ))
    }

    fn lf(&self, state: &BwtState, row: usize) -> Result<usize, BoundedBwtError> {
        let code = state.code(row);
        if code == SENTINEL_CODE {
            return Ok(0);
        }
        let first = 1_u64
            + state.counts[..usize::from(code)]
                .iter()
                .copied()
                .sum::<u64>();
        usize::try_from(
            first
                .checked_add(self.occ(state, code, row)?)
                .ok_or(BoundedBwtError::SizeOverflow)?,
        )
        .map_err(|_| BoundedBwtError::SizeOverflow)
    }
}

fn parallel_radix_sort_gaps(
    keys: &mut Vec<u64>,
    threads: u32,
    layout: KeyLayout,
    maximum_gap: usize,
) -> Result<(), BoundedBwtError> {
    if keys.len() < MIN_PARALLEL_RADIX_ROWS || threads == 1 {
        keys.sort_unstable_by_key(|value| layout.gap(*value));
        return Ok(());
    }
    let workers = usize::try_from(threads)
        .expect("validated thread count fits usize")
        .min(keys.len());
    let mut scratch = try_vec_filled(keys.len(), 0_u64, "bounded radix scratch")?;
    let gap_bits = (usize::BITS - maximum_gap.leading_zeros()).max(1);
    let passes = gap_bits.div_ceil(RADIX_BITS);
    for pass in 0..passes {
        let shift = layout.local_bits + pass * RADIX_BITS;
        let counts = thread::scope(|scope| {
            let mut handles = Vec::with_capacity(workers);
            for worker in 0..workers {
                let begin = keys.len() * worker / workers;
                let end = keys.len() * (worker + 1) / workers;
                let input = &keys[begin..end];
                handles.push(scope.spawn(move || {
                    let mut local = vec![0_usize; RADIX_BUCKETS];
                    for &value in input {
                        local[radix_bucket(value, shift)] += 1;
                    }
                    local
                }));
            }
            let mut rows = Vec::with_capacity(workers);
            for handle in handles {
                rows.push(handle.join().map_err(|_| BoundedBwtError::WorkerPanic)?);
            }
            Ok::<_, BoundedBwtError>(rows)
        })?;
        let mut positions = counts;
        let mut total = 0_usize;
        for bucket in 0..RADIX_BUCKETS {
            let mut position = total;
            for worker_positions in &mut positions {
                let count = worker_positions[bucket];
                worker_positions[bucket] = position;
                position += count;
            }
            total = position;
        }
        if total != keys.len() {
            return Err(BoundedBwtError::Invariant(
                "parallel radix counts do not cover every key",
            ));
        }
        let output_address = scratch.as_mut_ptr().expose_provenance();
        thread::scope(|scope| {
            let mut handles = Vec::with_capacity(workers);
            for (worker, mut worker_positions) in positions.into_iter().enumerate() {
                let begin = keys.len() * worker / workers;
                let end = keys.len() * (worker + 1) / workers;
                let input = &keys[begin..end];
                handles.push(scope.spawn(move || {
                    let output = ptr::with_exposed_provenance_mut::<u64>(output_address);
                    for &value in input {
                        let bucket = radix_bucket(value, shift);
                        let position = worker_positions[bucket];
                        // SAFETY: worker-specific bucket prefixes are disjoint
                        // and cover the complete initialized destination.
                        unsafe { output.add(position).write(value) };
                        worker_positions[bucket] += 1;
                    }
                }));
            }
            for handle in handles {
                handle.join().map_err(|_| BoundedBwtError::WorkerPanic)?;
            }
            Ok::<_, BoundedBwtError>(())
        })?;
        core::mem::swap(keys, &mut scratch);
    }
    Ok(())
}

fn sort_equal_gap_suffixes(
    keys: &mut [u64],
    text: &PackedProjectedText,
    block_start: usize,
    threads: u32,
    layout: KeyLayout,
) -> Result<(), BoundedBwtError> {
    if keys.len() < 2 {
        return Ok(());
    }
    let workers = usize::try_from(threads)
        .expect("validated thread count fits usize")
        .min(keys.len());
    let mut boundaries = Vec::with_capacity(workers + 1);
    boundaries.push(0);
    for worker in 1..workers {
        let mut boundary = keys.len() * worker / workers;
        while boundary < keys.len() && layout.gap(keys[boundary - 1]) == layout.gap(keys[boundary])
        {
            boundary += 1;
        }
        if boundary > *boundaries.last().expect("initial boundary exists") {
            boundaries.push(boundary);
        }
    }
    if *boundaries.last().expect("initial boundary exists") != keys.len() {
        boundaries.push(keys.len());
    }
    thread::scope(|scope| {
        let mut remaining = keys;
        let mut handles = Vec::with_capacity(boundaries.len() - 1);
        for window in boundaries.windows(2) {
            let (chunk, rest) = remaining.split_at_mut(window[1] - window[0]);
            remaining = rest;
            handles.push(scope.spawn(move || {
                let mut first = 0_usize;
                while first < chunk.len() {
                    let gap = layout.gap(chunk[first]);
                    let last =
                        first + chunk[first..].partition_point(|value| layout.gap(*value) == gap);
                    if last - first > 1 {
                        chunk[first..last].sort_unstable_by(|left, right| {
                            text.compare_suffixes(
                                block_start + layout.local(*left),
                                block_start + layout.local(*right),
                            )
                        });
                    }
                    first = last;
                }
            }));
        }
        for handle in handles {
            handle.join().map_err(|_| BoundedBwtError::WorkerPanic)?;
        }
        Ok::<_, BoundedBwtError>(())
    })
}

fn gaps_to_output_positions(
    keys: &mut [u64],
    old_rows: usize,
    layout: KeyLayout,
) -> Result<(), BoundedBwtError> {
    let output_rows = old_rows
        .checked_add(keys.len())
        .ok_or(BoundedBwtError::SizeOverflow)?;
    let mut previous = None;
    for (ordinal, value) in keys.iter_mut().enumerate() {
        let gap = layout.gap(*value);
        if gap > old_rows {
            return Err(BoundedBwtError::Invariant(
                "sorted gap exceeds the old suffix domain",
            ));
        }
        let position = gap
            .checked_add(ordinal)
            .ok_or(BoundedBwtError::SizeOverflow)?;
        if position >= output_rows || previous.is_some_and(|prior| prior >= position) {
            return Err(BoundedBwtError::Invariant(
                "new output rows are not strictly increasing",
            ));
        }
        *value = layout.pack(position, layout.local(*value))?;
        previous = Some(position);
    }
    Ok(())
}

fn radix_bucket(value: u64, shift: u32) -> usize {
    usize::try_from((value >> shift) & u64::from(u16::MAX)).expect("masked radix bucket fits usize")
}

fn sample_quotient(coordinate: usize, sample_stride: usize) -> Option<u32> {
    coordinate.is_multiple_of(sample_stride).then(|| {
        u32::try_from(coordinate / sample_stride).expect("validated sparse-SA quotient fits u32")
    })
}

fn sampled_coordinates(
    start: usize,
    end: usize,
    sample_stride: usize,
) -> Result<usize, BoundedBwtError> {
    if start >= end {
        return Ok(0);
    }
    let first = start
        .checked_add(sample_stride - 1)
        .ok_or(BoundedBwtError::SizeOverflow)?
        / sample_stride;
    let last = (end - 1) / sample_stride;
    Ok(if first > last { 0 } else { last - first + 1 })
}

fn push_row(
    packed: &mut Vec<u8>,
    samples: &mut Vec<u32>,
    rows: usize,
    code: u8,
    sample: Option<u32>,
) {
    debug_assert!(code <= SENTINEL_CODE);
    let nibble = code | (u8::from(sample.is_some()) * SAMPLE_FLAG);
    if rows.is_multiple_of(2) {
        packed.push(nibble);
    } else {
        *packed.last_mut().expect("odd row has a packed byte") |= nibble << 4;
    }
    if let Some(sample) = sample {
        samples.push(sample);
    }
}

#[inline]
fn set_nibble(packed: &mut [u8], row: usize, nibble: u8) {
    let shift = 4 * (row % 2);
    let mask = 15_u8 << shift;
    packed[row / 2] = (packed[row / 2] & !mask) | ((nibble & 15) << shift);
}

#[derive(Clone, Copy, Debug)]
struct ParallelMergePart {
    output_begin: usize,
    output_end: usize,
    key_begin: usize,
    key_end: usize,
    old_begin: usize,
    old_end: usize,
    old_sample_begin: usize,
    old_sample_end: usize,
    output_sample_begin: usize,
    output_sample_end: usize,
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn fill_parallel_merge_part(
    text: &PackedProjectedText,
    block_start: usize,
    keys: &[u64],
    layout: KeyLayout,
    old_global_begin: usize,
    old_start_row: usize,
    preceding: u8,
    old_rows: &[u8],
    old_samples: &[u32],
    part: ParallelMergePart,
    packed_address: usize,
    sample_address: usize,
    sample_stride: usize,
) -> Result<Option<usize>, BoundedBwtError> {
    let mut key = part.key_begin;
    let mut old = part.old_begin;
    let mut old_sample = part.old_sample_begin;
    let mut output_sample = part.output_sample_begin;
    let mut new_start_row = None;
    for output_row in part.output_begin..part.output_end {
        let is_new = key < part.key_end && layout.output_position(keys[key]) == output_row;
        let nibble = if is_new {
            let local = layout.local(keys[key]);
            let absolute = block_start
                .checked_add(local)
                .ok_or(BoundedBwtError::SizeOverflow)?;
            let code = if local == 0 {
                if new_start_row.replace(output_row).is_some() {
                    return Err(BoundedBwtError::Invariant(
                        "merge partition produced multiple start rows",
                    ));
                }
                SENTINEL_CODE
            } else {
                text.get(absolute - 1)
            };
            key += 1;
            if let Some(quotient) = sample_quotient(absolute, sample_stride) {
                write_parallel_sample(sample_address, output_sample, quotient);
                output_sample += 1;
                code | SAMPLE_FLAG
            } else {
                code
            }
        } else {
            if key < part.key_end && layout.output_position(keys[key]) < output_row {
                return Err(BoundedBwtError::Invariant(
                    "merge partition skipped a new row",
                ));
            }
            if old >= part.old_end {
                return Err(BoundedBwtError::Invariant(
                    "merge partition exhausted old rows early",
                ));
            }
            let mut nibble = (old_rows[old / 2] >> (4 * (old % 2))) & 15;
            let old_global_row = old_global_begin
                .checked_add(old)
                .ok_or(BoundedBwtError::SizeOverflow)?;
            old += 1;
            if nibble & SAMPLE_FLAG != 0 {
                let quotient = *old_samples
                    .get(old_sample)
                    .ok_or(BoundedBwtError::Invariant(
                        "merge partition old sample value is absent",
                    ))?;
                old_sample += 1;
                write_parallel_sample(sample_address, output_sample, quotient);
                output_sample += 1;
            }
            if old_global_row == old_start_row {
                if nibble & SYMBOL_MASK != SENTINEL_CODE {
                    return Err(BoundedBwtError::Invariant(
                        "merge partition old sentinel changed symbol",
                    ));
                }
                nibble = (nibble & SAMPLE_FLAG) | preceding;
            }
            nibble
        };
        write_parallel_nibble(packed_address, output_row, nibble);
    }
    if key != part.key_end
        || old != part.old_end
        || old_sample != part.old_sample_end
        || output_sample != part.output_sample_end
    {
        return Err(BoundedBwtError::Invariant(
            "merge partition did not consume its exact rows and samples",
        ));
    }
    Ok(new_start_row)
}

#[inline]
fn write_parallel_nibble(address: usize, row: usize, nibble: u8) {
    // SAFETY: merge partition boundaries are even, except the sequential outer
    // chunk boundaries. Concurrent workers therefore own disjoint destination
    // bytes and every row lies in the resized packed allocation.
    unsafe {
        let byte = ptr::with_exposed_provenance_mut::<u8>(address).add(row / 2);
        let shift = 4 * (row % 2);
        let mask = 15_u8 << shift;
        *byte = (*byte & !mask) | ((nibble & 15) << shift);
    }
}

#[inline]
fn write_parallel_sample(address: usize, ordinal: usize, quotient: u32) {
    // SAFETY: sample-prefix accounting assigns each partition one disjoint,
    // in-bounds range in the resized quotient vector.
    unsafe {
        ptr::with_exposed_provenance_mut::<u32>(address)
            .add(ordinal)
            .write(quotient);
    }
}

fn count_new_key_samples(
    keys: &[u64],
    layout: KeyLayout,
    block_start: usize,
    sample_stride: usize,
) -> usize {
    keys.iter()
        .filter(|&&key| (block_start + layout.local(key)).is_multiple_of(sample_stride))
        .count()
}

fn copy_nibbles_to_zero(source: &[u8], first_row: usize, rows: usize, output: &mut Vec<u8>) {
    output.clear();
    output.resize(rows.div_ceil(2), 0);
    if rows == 0 {
        return;
    }
    if first_row.is_multiple_of(2) {
        output.copy_from_slice(&source[first_row / 2..(first_row + rows).div_ceil(2)]);
        if !rows.is_multiple_of(2) {
            *output.last_mut().expect("odd row range has a byte") &= 15;
        }
        return;
    }

    let source_byte = first_row / 2;
    let complete_bytes = rows / 2;
    let mut offset = 0_usize;
    while offset + 8 <= complete_bytes {
        // SAFETY: sixteen requested rows beginning at a high nibble cover the
        // eight loaded bytes and the following spill byte. Output has eight
        // initialized destination bytes at this offset.
        unsafe {
            let low = ptr::read_unaligned(source.as_ptr().add(source_byte + offset).cast::<u64>());
            let spill = u64::from(source[source_byte + offset + 8]);
            ptr::write_unaligned(
                output.as_mut_ptr().add(offset).cast::<u64>(),
                (low >> 4) | (spill << 60),
            );
        }
        offset += 8;
    }
    while offset < complete_bytes {
        output[offset] =
            (source[source_byte + offset] >> 4) | ((source[source_byte + offset + 1] & 15) << 4);
        offset += 1;
    }
    if !rows.is_multiple_of(2) {
        output[complete_bytes] = source[source_byte + complete_bytes] >> 4;
    }
}

#[allow(clippy::too_many_arguments)]
fn move_old_run_toward_end(
    packed: &mut [u8],
    samples: &mut [u32],
    source: usize,
    destination: usize,
    rows: usize,
    old_sample_end: &mut usize,
    output_sample_end: &mut usize,
    old_start_row: usize,
    mapped_old_start_row: &mut Option<usize>,
) -> Result<(), BoundedBwtError> {
    if destination < source
        || source.checked_add(rows).is_none()
        || destination.checked_add(rows).is_none()
    {
        return Err(BoundedBwtError::Invariant(
            "old-row run move exceeds coordinate domain",
        ));
    }
    if rows == 0 {
        return Ok(());
    }
    let run_samples = usize::try_from(count_sample_range(packed, source, rows))
        .map_err(|_| BoundedBwtError::SizeOverflow)?;
    let old_sample_begin =
        old_sample_end
            .checked_sub(run_samples)
            .ok_or(BoundedBwtError::Invariant(
                "old sample run exceeds remaining input",
            ))?;
    let output_sample_begin =
        output_sample_end
            .checked_sub(run_samples)
            .ok_or(BoundedBwtError::Invariant(
                "old sample run exceeds remaining output",
            ))?;
    if output_sample_begin < old_sample_begin {
        return Err(BoundedBwtError::Invariant(
            "in-place sample run would overwrite unread input",
        ));
    }
    samples.copy_within(old_sample_begin..*old_sample_end, output_sample_begin);
    *old_sample_end = old_sample_begin;
    *output_sample_end = output_sample_begin;

    if (source..source + rows).contains(&old_start_row) {
        let mapped = destination + old_start_row - source;
        if mapped_old_start_row.replace(mapped).is_some() {
            return Err(BoundedBwtError::Invariant(
                "old sentinel row was mapped more than once",
            ));
        }
    }
    copy_nibbles_toward_end(packed, source, destination, rows);
    Ok(())
}

fn copy_nibbles_toward_end(
    packed: &mut [u8],
    mut source: usize,
    mut destination: usize,
    mut rows: usize,
) {
    debug_assert!(destination >= source);
    debug_assert!((source + rows).div_ceil(2) <= packed.len());
    debug_assert!((destination + rows).div_ceil(2) <= packed.len());
    if rows == 0 || source == destination {
        return;
    }

    // Save an odd destination prefix before any overlapping high-to-low copy.
    // It is written last because its destination nibble may still be an input
    // nibble for the remainder of a one-row-shifted run.
    let prefix = if destination.is_multiple_of(2) {
        None
    } else {
        let value = (packed[source / 2] >> (4 * (source % 2))) & 15;
        source += 1;
        destination += 1;
        rows -= 1;
        Some((destination - 1, value))
    };

    if !rows.is_multiple_of(2) {
        let source_tail = source + rows - 1;
        let destination_tail = destination + rows - 1;
        let value = (packed[source_tail / 2] >> (4 * (source_tail % 2))) & 15;
        set_nibble(packed, destination_tail, value);
        rows -= 1;
    }

    let bytes = rows / 2;
    if bytes != 0 && source.is_multiple_of(2) {
        // SAFETY: both byte ranges are in bounds. `ptr::copy` has memmove
        // semantics and therefore preserves overlapping moves toward higher
        // addresses.
        unsafe {
            ptr::copy(
                packed.as_ptr().add(source / 2),
                packed.as_mut_ptr().add(destination / 2),
                bytes,
            );
        }
    } else if bytes != 0 {
        // The source begins at a high nibble while the destination begins at a
        // low nibble. Eight destination bytes are one shifted u64 source word
        // plus one spill byte. Work high-to-low to retain memmove semantics.
        let source_byte = source / 2;
        let destination_byte = destination / 2;
        let mut remaining = bytes;
        while remaining >= 8 {
            let first = remaining - 8;
            // SAFETY: sixteen source nibbles starting at an odd row occupy the
            // eight loaded bytes plus the validated spill byte. Destination
            // bytes are an in-bounds disjoint chunk of the output range.
            unsafe {
                let low =
                    ptr::read_unaligned(packed.as_ptr().add(source_byte + first).cast::<u64>());
                let spill = u64::from(packed[source_byte + first + 8]);
                let shifted = (low >> 4) | (spill << 60);
                ptr::write_unaligned(
                    packed
                        .as_mut_ptr()
                        .add(destination_byte + first)
                        .cast::<u64>(),
                    shifted,
                );
            }
            remaining = first;
        }
        while remaining != 0 {
            let index = remaining - 1;
            let low = packed[source_byte + index] >> 4;
            let high = (packed[source_byte + index + 1] & 15) << 4;
            packed[destination_byte + index] = low | high;
            remaining = index;
        }
    }

    if let Some((row, value)) = prefix {
        set_nibble(packed, row, value);
    }
}

fn count_digits(digits: &[u8]) -> Result<[u64; 3], BoundedBwtError> {
    let mut counts = [0_u64; 3];
    for &digit in digits {
        let count = counts
            .get_mut(usize::from(digit))
            .ok_or(BoundedBwtError::InvalidInput(
                "projected text contains a digit outside G/T/A",
            ))?;
        *count = count.checked_add(1).ok_or(BoundedBwtError::SizeOverflow)?;
    }
    Ok(counts)
}

fn packed_word(packed: &[u8], first_row: usize, rows: usize) -> u64 {
    debug_assert!(rows <= 16);
    if first_row.is_multiple_of(2) && rows == 16 {
        // SAFETY: sixteen rows occupy exactly eight bytes and the caller only
        // requests complete words inside the packed row extent.
        return unsafe { ptr::read_unaligned(packed.as_ptr().add(first_row / 2).cast::<u64>()) };
    }
    let mut bytes = [0_u8; 8];
    for row in 0..rows {
        let source = first_row + row;
        let nibble = (packed[source / 2] >> (4 * (source % 2))) & 15;
        bytes[row / 2] |= nibble << (4 * (row % 2));
    }
    u64::from_le_bytes(bytes)
}

fn count_code_range(packed: &[u8], first_row: usize, rows: usize, code: u8) -> u32 {
    let mut count = 0_u32;
    for offset in (0..rows).step_by(16) {
        let width = (rows - offset).min(16);
        count += count_code_word(packed_word(packed, first_row + offset, width), code, width);
    }
    count
}

fn count_sample_range(packed: &[u8], first_row: usize, rows: usize) -> u32 {
    let mut count = 0_u32;
    for offset in (0..rows).step_by(16) {
        let width = (rows - offset).min(16);
        count += count_sample_word(packed_word(packed, first_row + offset, width), width);
    }
    count
}

fn count_code_word(word: u64, code: u8, rows: usize) -> u32 {
    let low = 0x1111_1111_1111_1111_u64;
    let used = if rows == 16 {
        low
    } else {
        low & ((1_u64 << (rows * 4)) - 1)
    };
    let symbols = word & 0x3333_3333_3333_3333;
    let mut matches = used;
    for bit in 0..2 {
        let values = (symbols >> bit) & low;
        matches &= if code & (1 << bit) == 0 {
            !values & low
        } else {
            values
        };
    }
    matches.count_ones()
}

fn count_sample_word(word: u64, rows: usize) -> u32 {
    let flags = (word >> 2) & 0x1111_1111_1111_1111;
    let used = if rows == 16 {
        u64::MAX
    } else {
        (1_u64 << (rows * 4)) - 1
    };
    (flags & used).count_ones()
}

fn try_reserve<T>(
    values: &mut Vec<T>,
    additional: usize,
    label: &'static str,
) -> Result<(), BoundedBwtError> {
    values
        .try_reserve_exact(additional)
        .map_err(|_| BoundedBwtError::Allocation(label))
}

fn try_vec_filled<T: Clone>(
    length: usize,
    value: T,
    label: &'static str,
) -> Result<Vec<T>, BoundedBwtError> {
    let mut values = Vec::new();
    try_reserve(&mut values, length, label)?;
    values.resize(length, value);
    Ok(values)
}

#[cfg(test)]
mod tests {
    use super::*;
    use bsbit_core::sequence::normalize_dna;

    fn contig(name: &[u8], bases: &[u8]) -> ContigInput {
        ContigInput::new(name.to_vec(), normalize_dna(bases).expect("valid test DNA"))
    }

    fn naive(digits: &[u8], sample_stride: usize) -> (Vec<u8>, usize, Vec<u32>) {
        let mut suffixes = (0..=digits.len()).collect::<Vec<_>>();
        suffixes.sort_unstable_by(|&left, &right| digits[left..].cmp(&digits[right..]));
        let sentinel = suffixes.iter().position(|&suffix| suffix == 0).unwrap();
        let mut transformed = Vec::with_capacity(digits.len());
        let mut samples = Vec::new();
        for &suffix in &suffixes {
            if suffix != 0 {
                transformed.push(digits[suffix - 1]);
            }
            if let Some(quotient) = sample_quotient(suffix, sample_stride) {
                samples.push(quotient);
            }
        }
        (transformed, sentinel, samples)
    }

    fn assert_builder(digits: &[u8], block: usize, threads: u32, sample_stride: usize) {
        let text = PackedProjectedText::from_digits(digits).unwrap();
        let config = BoundedBwtConfig::new(64, threads)
            .unwrap()
            .with_sample_stride(sample_stride)
            .unwrap()
            .with_block_bases(block)
            .unwrap();
        let actual = build_bounded_bwt(text, config).unwrap();
        let (transformed, sentinel, samples) = naive(digits, sample_stride);
        assert_eq!(
            actual.sentinel_row(),
            sentinel,
            "digits={digits:?}, block={block}"
        );
        assert_eq!(
            (0..digits.len())
                .map(|line| actual.transformed_digit(line))
                .collect::<Vec<_>>(),
            transformed,
            "digits={digits:?}, block={block}"
        );
        assert_eq!(
            actual.sample_quotients(),
            samples,
            "digits={digits:?}, block={block}"
        );
    }

    #[test]
    fn packed_projection_matches_frozen_byte_projection() {
        let catalog = [contig(b"one", b"ACN"), contig(b"two", b"GTNAC")];
        for threads in [1, 2, 8] {
            let packed = project_combined_packed_text(&catalog, 17, threads).unwrap();
            let expected = (0..packed.len())
                .map(|position| packed.get(position))
                .collect::<Vec<_>>();
            assert_eq!(expected.len(), 16);
            assert_eq!(packed.reference_bases(), 8);
            assert!(expected.iter().all(|&digit| digit <= 2));
            if threads > 1 {
                let scalar = project_combined_packed_text(&catalog, 17, 1).unwrap();
                assert_eq!(scalar.bytes, packed.bytes);
            }
        }
    }

    #[test]
    fn bounded_builder_exhaustively_matches_naive_binary_texts() {
        for length in 2..=10 {
            for mut bits in 0_usize..(1_usize << length) {
                let mut digits = vec![0_u8; length];
                for digit in &mut digits {
                    *digit = u8::from(bits & 1 != 0) * 2;
                    bits >>= 1;
                }
                for block in 1..=length.min(4) {
                    for sample_stride in [8, 16] {
                        assert_builder(&digits, block, 1, sample_stride);
                        assert_builder(&digits, block, 2, sample_stride);
                    }
                }
            }
        }
    }

    #[test]
    fn bounded_builder_matches_naive_across_stride_and_rank_boundaries() {
        for length in [15_usize, 16, 17, 63, 64, 65, 255, 256, 257, 1_031] {
            let digits = (0..length)
                .map(|offset| u8::try_from((offset * 11 + offset / 7) % 3).unwrap())
                .collect::<Vec<_>>();
            for block in [3, 17, 61, 127] {
                for sample_stride in [8, 16] {
                    assert_builder(&digits, block.min(length), 4, sample_stride);
                }
            }
        }
    }
}
