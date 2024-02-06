//! Helpers Module
//!
//! Contains helper methods and structures that are not explicitly tied to the physics of FDTD but useful nonetheless

// import local crates and modules
use anyhow::{anyhow, Ok, Result};
use serde_derive::Deserialize;
use std::fs::read_to_string;
use std::path::PathBuf;
use std::str::FromStr;
use toml::from_str;

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
