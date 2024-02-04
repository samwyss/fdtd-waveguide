//! Solver module
//!
//! Acts as a high level interface for controlling the simulation containing relevant data & code to achieve said goal

// import local crates and modules
use super::{geometry::Geometry, helpers::Config, Ok, Result};

/// Solver struct
///
/// Contains all data & methods needed for a high level simulation interface
#[derive(Debug)]
pub struct Solver {
    geometry: Geometry, // Geometry struct containing all data & code pertaining to the geometry of the simulation (see ./src/geometry/mod.rs)
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

        // time increment is set to 0 in constructor
        let dt = 0.0;

        // current time assumed to be 0 in constructor
        let cur_time = 0.0;

        // current time index assumed to be 0 in constructor
        let cur_time_idx: usize = 0;

        // target time assumed to be 0 in constructor
        let targ_time = 0.0;

        // target time index assumed to be 0 in constructor
        let targ_time_idx: usize = 0;

        Ok(Solver { geometry })
    }
}
