// cargo crate imports
use std::error::Error;

// declare all local modules and conditionally expose them to main.rs
pub mod solver;
mod engine;
mod geometry;
mod helpers;

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