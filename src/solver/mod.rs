//! Solver module
//!
//! Acts as a high level interface for controlling the simulation

// import local crates and modules
use super::{engine::Engine, geometry::Geometry, helpers::Config};
use anyhow::{Ok, Result};

/// Solver struct
///
/// Contains all data & methods needed for a high level simulation interface
#[derive(Debug)]
pub struct Solver {
    geometry: Geometry, // Geometry struct containing all data & code pertaining to the geometry of the simulation (see ./src/geometry/mod.rs)
    engine: Engine, // Engine struct containing all data relevant to the state of the simulation and code needed to evolve said state (see ./src/engine/mod.rs)
}

impl Solver {
    /// Solver constructor
    ///
    /// # Arguments
    ///
    /// `path_str` &str path to desired configuration file
    ///
    /// # Returns
    ///
    /// `Result<Solver>`
    ///
    /// # Errors
    ///
    pub fn new(path_str: &str) -> Result<Solver> {
        // create a new config from supplied path (see ./src/helpers/mod.rs)
        let config: Config = Config::new(path_str)?;

        // create a new geometry for given configuration (see ./src/geometry/mod.rs)
        let geometry = Geometry::new(config)?;

        // create a new engine for the given geometry (see ./src/engine/mod.rs)
        let engine = Engine::new(&geometry)?;

        Ok(Solver { geometry, engine })
    }

    pub fn update(&mut self, target_time: f64) -> Result<()> {
        // update the engine to the target time
        self.engine.update(&self.geometry, target_time)?;

        Ok(())
    }
}
