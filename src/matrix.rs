//! Activation matrices, stored one row per bitmask: bit `c` of `rows[r]` is
//! `matrix[r][c]`, i.e. "channel `r` feeds channel `c`".

use crate::rng::SplitMix64;

/// Largest channel count supported (rows are `u16` bitmasks).
pub const MAX_CHANNELS: usize = 16;

/// Number of distinct `C x C` 0/1 matrices.
pub fn matrix_count(channels: usize) -> u128 {
    1u128 << (channels * channels)
}

/// Whether every matrix for `channels` can be addressed by a `u64` index.
pub fn enumerable(channels: usize) -> bool {
    channels * channels <= 64
}

/// Decode an index into rows. Bit order matches the original implementation:
/// the most significant bit of the `C*C`-bit number is `matrix[0][0]`, the
/// least significant bit is `matrix[C-1][C-1]` (row-major, MSB first).
pub fn decode<const C: usize>(index: u64) -> [u16; C] {
    debug_assert!(enumerable(C));
    let bits = C * C;
    let mut rows = [0u16; C];
    for (r, row) in rows.iter_mut().enumerate() {
        // Bits of row r occupy positions [(C-1-r)*C, (C-r)*C); column 0 is the top bit.
        let shift = (bits - (r + 1) * C) as u32;
        let chunk = (index >> shift) & ((1u64 << C) - 1);
        // Reverse so that column c maps to bit c.
        *row = (chunk.reverse_bits() >> (64 - C)) as u16;
    }
    rows
}

/// Random rows with independent fair bits.
pub fn random<const C: usize>(rng: &mut SplitMix64) -> [u16; C] {
    let mask = ((1u32 << C) - 1) as u16;
    let mut rows = [0u16; C];
    for row in rows.iter_mut() {
        *row = (rng.next_u64() as u16) & mask;
    }
    rows
}

/// Rows as a nested 0/1 vector, for display and tests.
pub fn to_vec<const C: usize>(rows: &[u16; C]) -> Vec<Vec<u8>> {
    rows.iter()
        .map(|&row| (0..C).map(|c| ((row >> c) & 1) as u8).collect())
        .collect()
}
