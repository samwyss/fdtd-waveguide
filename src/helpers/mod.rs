//! Helpers Module
//!
//! Contains helper methods and structures that are not explicitly tied to the physics of FDTD but useful nonetheless

// import local crates and modules
use super::{anyhow, from_str, read_to_string, Deserialize, FromStr, Ok, PathBuf, Result};

/// Config Struct
///
/// Contains all data imported from ./config.toml which is used to derive other simulation structures
#[derive(Debug, Deserialize)]
pub struct Config {
    max_frequency: f64,               // [Hz] maximium frequency to resolve
    voxels_per_min_wavelength: usize, // [] number of voxels per minimum wavelength
    voxels_per_min_feature: usize,    // [] number of voxels per minimum feature size
    x_len: f64,                       // [m] length of waveguide in x-direction
    y_len: f64,                       // [m] length of waveguide in y-direction
    z_len: f64, // [m] length of waveguide in z-direction (direction of propagation)
    ep_r: f64,  // [] diagonally isotropic relative permittivity inside waveguide
    mu_r: f64,  // [] diagonally isotropic relative permeability inside waveguide
}

impl Config {
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

        // parse data
        let config: Config = from_str(&file_string)?;

        Ok(config)
    }
}
