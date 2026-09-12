//! One universe: `C/2` elements, each with a positive channel `2e` (its `count`)
//! and a negative channel `2e+1` (its `antiCount`).
//!
//! One step:
//! 1. **activation** – every active channel `r` pushes its matrix row into the
//!    elements: a 1 in column `c` adds +1 to channel `c` if `c` is even
//!    (count) or −1 if `c` is odd (antiCount).
//! 2. **evaluation** – each element collapses to its net value: a positive net
//!    keeps only `count`, a negative net keeps only `antiCount`, zero resets
//!    both. The sign of the net value selects which of the element's channels
//!    is active on the next step.
//!
//! The run stops when no active channel produced any change, when no channel
//! is active, or after `max_steps` steps.

/// Default number of steps before a universe is classified.
pub const DEFAULT_MAX_STEPS: usize = 10;
/// Upper bound on `max_steps` (history is a fixed stack buffer).
pub const MAX_STEPS_CAP: usize = 32;

/// Classification of one channel after a run.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[repr(u8)]
pub enum Symbol {
    /// The run terminated before `max_steps` (dead or converged).
    Terminated = 0,
    /// The final value already occurred earlier in the history (bounded / periodic).
    Bounded = 1,
    /// Positive channel never repeated a value (monotone growth).
    Growth = 2,
    /// Negative channel never repeated a value (monotone decline).
    Decline = 3,
}

impl Symbol {
    /// Text used in signatures, identical to the original implementation.
    pub fn as_str(self) -> &'static str {
        match self {
            Symbol::Terminated => "|",
            Symbol::Bounded => "0",
            Symbol::Growth => "1",
            Symbol::Decline => "-1",
        }
    }

    pub fn from_code(code: u32) -> Symbol {
        match code & 3 {
            0 => Symbol::Terminated,
            1 => Symbol::Bounded,
            2 => Symbol::Growth,
            _ => Symbol::Decline,
        }
    }
}

/// Behaviour options shared by every universe in a run.
#[derive(Clone, Copy, Debug)]
pub struct Rules {
    pub max_steps: usize,
    /// Reproduce the original classifier, which compares the final snapshot
    /// against every earlier snapshot *except* the immediately preceding one.
    /// `false` compares against all earlier snapshots.
    pub skip_previous: bool,
}

impl Default for Rules {
    fn default() -> Self {
        Rules {
            max_steps: DEFAULT_MAX_STEPS,
            skip_previous: true,
        }
    }
}

/// Reusable simulator; call [`Universe::reset`] between runs to avoid reallocating.
#[derive(Clone, Debug)]
pub struct Universe<const C: usize> {
    pub rules: Rules,
    /// Row `r` is the bitmask of channels fed by channel `r`.
    pub rows: [u16; C],
    /// Current value per channel: even = count (≥ 0), odd = antiCount (≤ 0) after evaluation.
    pub values: [i32; C],
    /// Bitmask of channels that fire on the next activation.
    pub active: u16,
    /// One snapshot of `values` per completed step.
    pub history: [[i32; C]; MAX_STEPS_CAP],
    pub steps: usize,
}

impl<const C: usize> Universe<C> {
    pub fn new(rows: [u16; C], initial_active: u16, rules: Rules) -> Self {
        assert!(
            C >= 2 && C.is_multiple_of(2) && C <= 16,
            "channels must be even and at most 16"
        );
        assert!(
            rules.max_steps >= 1 && rules.max_steps <= MAX_STEPS_CAP,
            "max_steps out of range"
        );
        Universe {
            rules,
            rows,
            values: [0; C],
            active: initial_active,
            history: [[0; C]; MAX_STEPS_CAP],
            steps: 0,
        }
    }

    /// Reinitialise for a new matrix without reallocating.
    #[inline]
    pub fn reset(&mut self, rows: [u16; C], initial_active: u16) {
        self.rows = rows;
        self.values = [0; C];
        self.active = initial_active;
        self.steps = 0;
    }

    /// Bitmask of every channel: all fire on the first step.
    pub const ALL: u16 = ((1u32 << C) - 1) as u16;

    /// Push every active row into the channels. Returns `false` if nothing was added.
    #[inline]
    pub fn activation(&mut self) -> bool {
        let mut fed = 0u16;
        let mut active = self.active;
        while active != 0 {
            let r = active.trailing_zeros() as usize;
            active &= active - 1;
            let mut row = self.rows[r];
            fed |= row;
            while row != 0 {
                let c = row.trailing_zeros() as usize;
                row &= row - 1;
                // even channel: +1 to count, odd channel: -1 to antiCount
                self.values[c] += 1 - 2 * (c as i32 & 1);
            }
        }
        fed != 0
    }

    /// Collapse each element to its net value, record the snapshot, compute the next active set.
    #[inline]
    pub fn evaluation(&mut self) {
        let mut next = 0u16;
        for e in 0..C / 2 {
            let (p, n) = (2 * e, 2 * e + 1);
            let net = self.values[p] + self.values[n];
            if net > 0 {
                self.values[p] = net;
                self.values[n] = 0;
                next |= 1 << p;
            } else if net < 0 {
                self.values[p] = 0;
                self.values[n] = net;
                next |= 1 << n;
            } else {
                self.values[p] = 0;
                self.values[n] = 0;
            }
        }
        self.history[self.steps] = self.values;
        self.steps += 1;
        self.active = next;
    }

    /// Run until quiescence or `max_steps`. Returns the number of completed steps.
    pub fn run(&mut self) -> usize {
        while self.active != 0 && self.steps < self.rules.max_steps {
            if !self.activation() {
                break;
            }
            self.evaluation();
        }
        self.steps
    }

    /// Classify one channel after [`Universe::run`].
    pub fn symbol(&self, channel: usize) -> Symbol {
        if self.steps < self.rules.max_steps {
            return Symbol::Terminated;
        }
        let last = self.history[self.steps - 1][channel];
        let compare = if self.rules.skip_previous {
            self.steps - 2
        } else {
            self.steps - 1
        };
        if self.history[..compare]
            .iter()
            .any(|snap| snap[channel] == last)
        {
            return Symbol::Bounded;
        }
        if channel.is_multiple_of(2) {
            Symbol::Growth
        } else {
            Symbol::Decline
        }
    }

    /// Packed signature: two bits per channel, channel `c` at bits `2c..2c+2`.
    pub fn code(&self) -> u32 {
        (0..C).fold(0u32, |code, c| code | ((self.symbol(c) as u32) << (2 * c)))
    }

    /// Human-readable signature, identical to the original output format:
    /// `count,anti` symbols for every element, comma-joined.
    pub fn signature(&self) -> String {
        decode_signature::<C>(self.code())
    }
}

/// Convert a packed code back to the text signature.
pub fn decode_signature<const C: usize>(code: u32) -> String {
    (0..C)
        .map(|c| Symbol::from_code(code >> (2 * c)).as_str())
        .collect::<Vec<_>>()
        .join(",")
}
