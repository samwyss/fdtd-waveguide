pub struct Geometry{
    pub ds: f64,                // spatial increment [m]
    pub ds_inv: f64,            // inverse of ds used for optimization [m^-1]
    len: f64,                   // side length of bounding box [m]
    pub num_vox: usize,         // number of voxels in simulation []
    pub num_vox_side: usize,    // number of voxels along one of the sides of bounding box []
    pub ep_r: Vec<f64>,         // numerical dispersion corrected relative permittivity of all points in space []
    pub mu_r: Vec<f64>,         // numerical dispersion corrected relative permeability of all points in space []
    pub ep_r_inv: Vec<f64>,     // inverse of relative permittivity of all points in space []
    pub mu_r_inv: Vec<f64>,     // inverse of relative permeability of all points in space []
    pub n: Vec<f64>,            // refractive index of all points in space []
}