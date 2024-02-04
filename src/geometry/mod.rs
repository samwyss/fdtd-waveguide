//! Geometry Module
//!
//! Contains all data & code relevant to the geometry of this simulation.
//!
//! This is admittedly overkill for this simulation as the relative permittivity and permeability distributions are very simple. However this was left in for stylistic purposes as well as for future expandability if desired.

// import local crates and modules
use super::{helpers::Config, Ok, Result, C_0};

/// Geometry struct
///
/// Contains all data & methods associated with the geometry of this simulation
#[derive(Debug)]
pub struct Geometry {
    x_len: f64,       // [m] length of waveguide in x-direction
    y_len: f64,       // [m] length of waveguide in y-direction
    z_len: f64,       // [m] length of waveguide in z-direction (direction of propagation)
    ep_r: f64,        // [] diagonally isotropic relative permittivity of material inside waveguide
    mu_r: f64,        // [] diagonally isotropic relative permeability of material inside waveguide
    n: f64,           // [] refractive index of material inside of waveguide
    dx: f64,          // [m] spatial increment in x-direction
    dy: f64,          // [m] spatial increment in y-direction
    dz: f64,          // [m] spatial increment in z-direction
    num_vox_x: usize, // [] number of voxels in x-direction
    num_vox_y: usize, // [] number of voxels in y-direction
    num_vox_z: usize, // [] number of voxels in z-direction
    num_vox: usize,   // [] number of voxels in the simulation domain
    ep_r_inv: f64,    // [] inverse of relative permittivity of material inside waveguide
    mu_r_inv: f64,    // [] inverse of relative permeability of material inside waveguide
    dx_inv: f64,      // [m^-1] inverse of spatial increment in x-direction
    dy_inv: f64,      // [m^-1] inverse of spatial increment in y-direction
    dz_inv: f64,      // [m^-1] inverse of spatial increment in z-direction
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

        // assign n
        let n: f64 = (ep_r * mu_r).sqrt();

        // determine the minimum spatial step based on the minimum wavelength (maximum frequency) and the number of voxels desired to resolve a single wavelength
        let minimum_wavelength: f64 = C_0 / (config.max_frequency * n);
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
        // assign dx
        let dx: f64 = x_len / num_vox_x as f64;

        // assign dy
        let dy: f64 = y_len / num_vox_y as f64;

        // assign dz
        let dz: f64 = z_len / num_vox_z as f64;

        // calculate the inverse of properties that are frequently divided
        // this is done as a low-skill optimization as divides are computationally expensive in comparison to multiplies especially when done many times in loops

        // assign ep_r_inv
        let ep_r_inv: f64 = 1.0 / ep_r;

        // assign mu_r_inv
        let mu_r_inv: f64 = 1.0 / mu_r;

        // assign dx_inv
        let dx_inv: f64 = 1.0 / dx;

        // assign dy_inv
        let dy_inv: f64 = 1.0 / dy;

        // assign dz_inv
        let dz_inv: f64 = 1.0 / dz;

        Ok(Geometry {
            x_len,
            y_len,
            z_len,
            ep_r,
            mu_r,
            n,
            dx,
            dy,
            dz,
            num_vox_x,
            num_vox_y,
            num_vox_z,
            num_vox,
            ep_r_inv,
            mu_r_inv,
            dx_inv,
            dy_inv,
            dz_inv,
        })
    }
}
