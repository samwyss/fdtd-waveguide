//! Solver module
//!
//! Acts as a high level interface for controlling the simulation

// import local crates and modules
use super::{engine::Engine, geometry::Geometry};
use anyhow::{anyhow, Ok, Result};
use serde_derive::Deserialize;
use std::fs::read_to_string;
use std::path::PathBuf;
use std::str::FromStr;
use toml::from_str;

/// Solver struct
///
/// Contains all data & methods needed for a high level simulation interface
#[derive(Debug)]
pub struct Solver {
    geometry: Geometry, // Geometry struct containing all data & code pertaining to the geometry of the simulation (see ./src/geometry/mod.rs)
    engine: Engine, // Engine struct containing all data relevant to the state of the simulation and code needed to evolve said state (see ./src/engine/mod.rs)
    config: Config,
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
    pub fn new(config: Config) -> Result<Solver> {
        // create a new geometry for given configuration (see ./src/geometry/mod.rs)
        let geometry = Geometry::new(&config)?;

        // create a new engine for the given geometry (see ./src/engine/mod.rs)
        let engine = Engine::new(&config, &geometry)?;

        Ok(Solver {
            config,
            geometry,
            engine,
        })
    }

    pub fn update(&mut self) -> Result<()> {
        // update the engine to the target time
        self.engine.update(&self.config, &self.geometry)?;

        Ok(())
    }
}

/// Config Struct
///
/// Contains all data imported from ./config.toml which is used to derive other simulation structures
#[derive(Debug, Deserialize)]
pub struct Config {
    pub max_frequency: f64,               // [Hz] maximum frequency to resolve
    pub voxels_per_min_wavelength: usize, // [] number of voxels per minimum wavelength
    pub voxels_per_min_feature: usize,    // [] number of voxels per minimum feature size
    pub x_len: f64,                       // [m] length of waveguide in x-direction
    pub y_len: f64,                       // [m] length of waveguide in y-direction
    pub z_len: f64, // [m] length of waveguide in z-direction (direction of propagation)
    pub ep_r: f64,  // [] diagonally isotropic relative permittivity inside waveguide
    pub mu_r: f64,  // [] diagonally isotropic relative permeability inside waveguide
    pub sigma: f64, // [S/m] diagonally isotropic conductivity of material
    pub snapshot_steps: usize,
    pub buffered_snapshots: usize,
    pub end_time: f64, // [s] end time of the simulation
    pub frequency: f64,
    pub delay_time: f64,
    pub ramp_time: f64,
}

impl Config {
    /// Config constructor
    ///
    /// # Arguments
    ///
    /// `path_str` &str path to desired configuration file
    ///
    /// # Returns
    ///
    /// `Result<Config>`
    ///
    /// # Errors
    ///
    /// - `path_str` is not able to be converted to PathBuf
    /// - `path_str` is not a valid path on the filesystem
    /// - `path_str` does not correspond to a file
    /// - file at `path_str` cannot be read into a String
    /// - file string from `path_str` cannot be deserialized into toml
    pub fn new(path_str: &str) -> Result<Config> {
        // convert path_str to a cross platform PathBuf
        let path = PathBuf::from_str(path_str)?;

        // ensure path exists on filesystem
        if !path.exists() {
            return Err(anyhow!(
                "provided path does not exist on the host filesystem"
            ));
        }

        // ensure path corresponds to a file on filesystem
        if !path.is_file() {
            return Err(anyhow!(
                "provided path does not correspond to a valid file on the host filesystem"
            ));
        }

        // load in toml string from file
        let file_string = read_to_string(path)?;

        // deserialize file_string into toml
        let config: Config = from_str(&file_string)?;

        Ok(config)
    }
}
