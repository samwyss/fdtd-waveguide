// cargo crate imports
use std::error::Error;

// declare constants
pub const ETA_0: f64 = 376.730313668; // [Ohms] free space impedance
pub const EP_0: f64 = 8.8541878128e-12; // [F/m] free space vacuum permittivity
pub const MU_0: f64 = 1.25663706212e-6; // [H/m] free space vacuum permeability
pub const C_0: f64 = 299792458.0; // [m/s] free space speed of light

// declare all local modules and conditionally expose them to main.rs
mod engine;
mod geometry;
mod helpers;
pub mod solver;

// declare traits (interfaces) to be implemented
trait GeometryInterface {
    fn new(
        x_len: f64,
        y_len: f64,
        z_len: f64,
        ep_r: f64,
        mu_r: f64,
    ) -> Result<geometry::Geometry, Box<dyn Error>>; // Geometry struct constructor
}

trait EngineInterface {
    fn new(geometry: geometry::Geometry) -> Result<engine::Engine, Box<dyn Error>>; // Engine struct constructor
}

trait SolverInterface {
    fn new(path_str: &str) -> Result<solver::Solver, Box<dyn Error>>; // Solver struct constructor
}
