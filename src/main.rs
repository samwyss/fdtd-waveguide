//! driver binary crate
//!
//! this is compiled into a simulation executable

// import crates and local waveguide library
use anyhow::{Ok, Result};
use waveguide::solver::{Config, Solver};

/// main function
///
/// # Arguments
///
/// # Returns
///
/// `Result<()>`
///
/// # Errors
/// - `config::new()` errors
/// - `Solver::new()` errors
/// - `Solver::update()` errors
fn main() -> Result<()> {
    // create a new config from supplied path (see ./src/helpers/mod.rs)
    let config: Config = Config::new("./config.toml")?;

    // construct new solver
    let mut solver = Solver::new(config)?;

    // update solver to a target_time
    solver.update()?;

    Ok(())
}
