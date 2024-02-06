//! Engine module
//!
//! Contains all data relevant to the state of the simulation and code needed to evolve said state.
//!
//! Operates on a Geometry (see ./src/geometry/mod.rs).

// import local modules and cargo crates
use crate::geometry::Geometry;
use anyhow::{Ok, Result};

#[derive(Debug)]
pub struct Engine {
    dt: f64,
    cur_time: f64,
    cur_time_idx: usize,
    targ_time: f64,
    targ_time_idx: usize,
    ex: Vec<f64>,
    ey: Vec<f64>,
    ez: Vec<f64>,
    hx: Vec<f64>,
    hy: Vec<f64>,
    hz: Vec<f64>,
}

impl Engine {
    pub fn new(geometry: &Geometry) -> Result<Engine> {
        // assign dt to zero as it is not known at compile time
        let dt: f64 = 0.0;

        // assign cur_time to 0 as that is initial state of engine
        let cur_time: f64 = 0.0;

        // assign cur_time_idx to 0 as that is initial state of engine
        let cur_time_idx: usize = 0;

        // assign targ_time to 0 in constructor as it is not known at compile time
        let targ_time: f64 = 0.0;

        // assign targ_time_idx to zero as it is not known at compile time
        let targ_time_idx: usize = 0;

        // extract num_vox
        let num_vox = geometry.num_vox;

        // assign initial ex
        let ex: Vec<f64> = vec![0.0; num_vox];

        // assign initial ey
        let ey: Vec<f64> = vec![0.0; num_vox];

        // assign initial ez
        let ez: Vec<f64> = vec![0.0; num_vox];

        // assign initial hx
        let hx: Vec<f64> = vec![0.0; num_vox];

        // assign initial hy
        let hy: Vec<f64> = vec![0.0; num_vox];

        // assign initial hz
        let hz: Vec<f64> = vec![0.0; num_vox];

        // TODO: initialize constant and curl arrays

        Ok(Engine {
            dt,
            cur_time,
            cur_time_idx,
            targ_time,
            targ_time_idx,
            ex,
            ey,
            ez,
            hx,
            hy,
            hz,
        })
    }

    pub fn update(&mut self, target_time: f64) -> Result<()> {
        println!("updating to {}", target_time);
        Ok(())
    }
}
