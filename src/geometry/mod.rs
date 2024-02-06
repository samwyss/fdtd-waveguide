//! Geometry Module
//!
//! Contains all data & code relevant to the geometry of this simulation.
//!
//! This is admittedly overkill for this simulation as the permittivity, permeability, and conductivity distributions are very simple. However this was left in for stylistic purposes as well as for future expandability if desired.

// import local modules and cargo crates
use crate::{helpers::Config, C_0, EP_0, MU_0};
use anyhow::{Ok, Result};

/// Geometry struct
///
/// Contains all data & methods associated with the geometry of this simulation
#[derive(Debug)]
pub struct Geometry {
    pub ep: f64,     // [F/m] diagonally isotropic permittivity of material inside waveguide
    pub mu: f64,     // [H/m] diagonally isotropic permeability of material inside waveguide
    pub sigma: f64,  // [S/m] diagonally isotropic conductivity of material inside waveguide
    pub dx_inv: f64, // [m^-1] inverse spatial increment in x-direction
    pub dy_inv: f64, // [m^-1] inverse spatial increment in y-direction
    pub dz_inv: f64, // [m^-1] inverse spatial increment in z-direction
    pub num_vox_x: usize, // [] number of voxels in x-direction
    pub num_vox_y: usize, // [] number of voxels in y-direction
    pub num_vox_z: usize, // [] number of voxels in z-direction
    pub num_vox: usize, // [] number of voxels in the simulation domain
}

impl Geometry {
    /// Geometry constructor
    ///
    /// # Arguments
    ///
    /// `config` helpers::Config configuration struct
    ///
    /// # Returns
    ///
    /// `Result<Geometry>`
    ///
    /// # Errors
    ///
    pub fn new(config: Config) -> Result<Geometry> {
        // assign x_len
        let x_len: f64 = config.x_len;

        // assign y_len
        let y_len: f64 = config.y_len;

        // assign z_len
        let z_len: f64 = config.z_len;

        // assign ep_r
        let ep_r: f64 = config.ep_r;

        // assign mu_r
        let mu_r: f64 = config.mu_r;

        // assign ep
        let ep = ep_r * EP_0;

        // assign mu
        let mu = mu_r * MU_0;

        // assign sigma
        let sigma = config.sigma;

        // determine the minimum spatial step based on the minimum wavelength (maximum frequency) and the number of voxels desired to resolve a single wavelength
        let minimum_wavelength: f64 = C_0 / (config.max_frequency * (ep_r * mu_r).sqrt());
        let ds_minimum_wavelength: f64 =
            minimum_wavelength / config.voxels_per_min_wavelength as f64;

        // determine the minimum spatial step based on the minimum feature size
        // since the material inside of the waveguide is homogenous, the minimum feature size becomes the smallest side length of the waveguide
        let minimum_feature_size: f64 = vec![x_len, y_len, z_len]
            .iter()
            .fold(f64::INFINITY, |a, &b| a.min(b));

        let ds_minimum_feature_size: f64 =
            minimum_feature_size / config.voxels_per_min_feature as f64;

        // the minimum of these two calculated spatial increments is the smallest step required
        let ds: f64 = vec![ds_minimum_wavelength, ds_minimum_feature_size]
            .iter()
            .fold(f64::INFINITY, |a, &b| a.min(b));

        // this step now needs to be "snapped" to the actual geometry we are simulating
        // assign num_vox_x
        let num_vox_x: usize = (x_len / ds).ceil() as usize;

        // assign num_vox_y
        let num_vox_y: usize = (y_len / ds).ceil() as usize;

        // assign num_vox_z
        let num_vox_z: usize = (z_len / ds).ceil() as usize;

        // assign num_vox
        let num_vox: usize = num_vox_x * num_vox_y * num_vox_z;

        // use the snapped number of voxels to back solve for the spacing in all directions
        // These are stored as inverses as they are only every divided by in the update equations. Thus storing their inverses is a low-skill optimization as divides are computationally expensive in comparison to multiplies especially when done many times in loops
        // assign dx_inv
        let dx_inv: f64 = num_vox_x as f64 / x_len;

        // assign dy_inv
        let dy_inv: f64 = num_vox_y as f64 / y_len;

        // assign dz_inv
        let dz_inv: f64 = num_vox_z as f64 / z_len;

        Ok(Geometry {
            ep,
            mu,
            sigma,
            dx_inv,
            dy_inv,
            dz_inv,
            num_vox_x,
            num_vox_y,
            num_vox_z,
            num_vox,
        })
    }
}
