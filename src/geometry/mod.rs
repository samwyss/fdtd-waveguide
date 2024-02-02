use super::Error;

trait GeometryInterface {
    fn new(
        x_len: f64,
        y_len: f64,
        z_len: f64,
        ep_r: f64,
        mu_r: f64,
    ) -> Result<Geometry, Box<dyn Error>>; // Geometry struct constructor
}

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
