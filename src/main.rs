//! driver binary crate
//!
//! acts as a main function for the simulation

// import crates and local waveguide library
use anyhow::{Ok, Result};
use waveguide::solver::Solver;

/// main function
///
/// # Arguments
///
/// # Returns
///
/// `Result<()>`
///
/// # Errors
/// - `Solver::new()` errors
/// - `Solver::update()` errors
fn main() -> Result<()> {
    // construct new solver
    let mut solver = Solver::new("./config.toml")?;

    // update solver to a target_time
    solver.update(10e-9)?;

    Ok(())
}
