// import crates
use anyhow::{Ok, Result};
use fdtd_waveguide::solver::Solver;

// main driver function
fn main() -> Result<()> {
    // construct new solver
    let mut solver = Solver::new("./config.toml")?;

    // update solver to a target_time
    solver.update(10e-9)?;

    Ok(())
}
