// declare all local modules
mod solver;
mod engine;
mod geometry;
mod helpers;

// re-export solver to be visible in main.rs
pub use solver::*;