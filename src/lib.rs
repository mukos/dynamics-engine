//! Simulation of small "activation matrix" systems.
//!
//! A universe of `N` elements has `C = 2N` channels; a `C x C` 0/1 matrix says
//! which channels a firing channel feeds. The multiverse enumerates (or samples)
//! matrices, runs each universe for a fixed number of steps, and tallies the
//! resulting behaviour signatures. See `README.md` for the model.

pub mod matrix;
pub mod multiverse;
pub mod rng;
pub mod universe;
