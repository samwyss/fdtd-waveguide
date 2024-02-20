//! driver binary crate
//!
//! acts as a main function for the simulation

// import crates and local waveguide library
use anyhow::{Ok, Result};
use waveguide::solver::{Solver, Config};

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

    // create a new config from supplied path (see ./src/helpers/mod.rs)
    let config: Config = Config::new("./config.toml")?;

    // construct new solver
    let mut solver = Solver::new(&config)?;

    // update solver to a target_time
    solver.update(&config.end_time)?;

    Ok(())
}
