//! Engine module
//!
//! Contains all data relevant to the state of the simulation and code needed to evolve said state.
//!
//! Operates on a Geometry (see ./src/geometry/mod.rs).

// constant declarations
const ONE_OVER_TWO: f64 = 1.0 / 2.0;

// import local modules and cargo crates
use crate::{
    geometry::{self, Geometry},
    C_0,
};
use anyhow::{Ok, Result};

#[derive(Debug)]
pub struct Engine {
    cur_time: f64,
    ex: ScalarField,
    ey: ScalarField,
    ez: ScalarField,
    hx: ScalarField,
    hy: ScalarField,
    hz: ScalarField,
}

impl Engine {
    pub fn new(geometry: &Geometry) -> Result<Engine> {
        // assign dt to zero as it is not known at compile time
        let dt: f64 = 0.0;

        // assign cur_time to 0 as that is initial state of engine
        let cur_time: f64 = 0.0;

        // extract num_vox
        let num_vox = geometry.num_vox;

        // assign ex to 0 as that is initial state of engine
        let ex: ScalarField = ScalarField::new(
            0.0,
            geometry.num_vox,
            geometry.num_vox_x,
            geometry.num_vox_y,
        )?;

        // assign ey to 0 as that is initial state of engine
        let ey: ScalarField = ScalarField::new(
            0.0,
            geometry.num_vox,
            geometry.num_vox_x,
            geometry.num_vox_y,
        )?;

        // assign ez to 0 as that is initial state of engine
        let ez: ScalarField = ScalarField::new(
            0.0,
            geometry.num_vox,
            geometry.num_vox_x,
            geometry.num_vox_y,
        )?;

        // assign hx to 0 as that is initial state of engine
        let hx: ScalarField = ScalarField::new(
            0.0,
            geometry.num_vox,
            geometry.num_vox_x,
            geometry.num_vox_y,
        )?;

        // assign hy to 0 as that is initial state of engine
        let hy: ScalarField = ScalarField::new(
            0.0,
            geometry.num_vox,
            geometry.num_vox_x,
            geometry.num_vox_y,
        )?;

        // assign hz to 0 as that is initial state of engine
        let hz: ScalarField = ScalarField::new(
            0.0,
            geometry.num_vox,
            geometry.num_vox_x,
            geometry.num_vox_y,
        )?;

        Ok(Engine {
            cur_time,
            ex,
            ey,
            ez,
            hx,
            hy,
            hz,
        })
    }

    pub fn update(&mut self, geometry: &Geometry, target_time: f64) -> Result<()> {
        // calculate first pass at time step based on Courant–Friedrichs–Lewy stability condition
        let mut dt: f64 = (C_0
            * (geometry.dx_inv.powi(2) + geometry.dy_inv.powi(2) + geometry.dz_inv.powi(2)).sqrt())
        .powi(-1);

        // snap this time step to a specific number of loop iterations using target_time
        let time_steps: usize = (target_time / dt).ceil() as usize;

        // recalculate the time step based on the snapped number of loop iterations
        dt = target_time / time_steps as f64;

        // pre-process loop constants
        let ea: f64 = (geometry.ep / dt + geometry.sigma / 2.0).powi(-1);
        let eb: f64 = geometry.ep / dt - geometry.sigma / 2.0;
        let hax: f64 = dt * geometry.dx_inv / geometry.mu;
        let hay: f64 = dt * geometry.dy_inv / geometry.mu;
        let haz: f64 = dt * geometry.dz_inv / geometry.mu;

        // time loop
        for t in 0..time_steps {
            // update magnetic field
            self.update_h(&geometry, &hax, &hay, &haz)?;

            // update current engine time after magnetic field update
            self.cur_time += ONE_OVER_TWO * dt;

            // update electric field

            // update current engine time after electric field update
            self.cur_time += ONE_OVER_TWO * dt;

            //TODO remove me
            println!("{} of {}", t, time_steps);
        }

        Ok(())
    }

    fn update_h(&mut self, geometry: &Geometry, hax: &f64, hay: &f64, haz: &f64) -> Result<()> {
        // update hx
        self.update_hx(geometry, hay, haz)?;

        // update hy
        self.update_hy(geometry, hax, haz)?;

        // update hz
        self.update_hz(geometry, hax, hay)?;

        Ok(())
    }

    fn update_hx(&mut self, geometry: &Geometry, hay: &f64, haz: &f64) -> Result<()> {
        for k in 0..(geometry.num_vox_z - 1) {
            for j in 0..(geometry.num_vox_y - 1) {
                for i in 0..geometry.num_vox_x {
                    // hx update equation for non j-high and k-high volume
                    *self.hz.idxm(i, j, k) += -hay
                        * (&self.ez.idx(i, j + 1, k) - self.ez.idx(i, j, k))
                        + haz * (self.ey.idx(i, j, k + 1) - self.ey.idx(i, j, k));
                }
            }

            for i in 0..geometry.num_vox_x {
                // hx update equation for j-high plane
                *self.hz.idxm(i, geometry.num_vox_y - 1, k) += -hay
                    * (0.0 - self.ez.idx(i, geometry.num_vox_y - 1, k))
                    + haz
                        * (self.ey.idx(i, geometry.num_vox_y - 1, k + 1)
                            - self.ey.idx(i, geometry.num_vox_y - 1, k));
            }
        }

        for j in 0..(geometry.num_vox_y - 1) {
            for i in 0..geometry.num_vox_x {
                // hx update equation for k-high plane
                *self.hz.idxm(i, j, geometry.num_vox_z - 1) += -hay
                    * (&self.ez.idx(i, j + 1, geometry.num_vox_z - 1)
                        - self.ez.idx(i, j, geometry.num_vox_z - 1))
                    + haz * (0.0 - self.ey.idx(i, j, geometry.num_vox_z - 1));
            }
        }

        for i in 0..geometry.num_vox_x {
            // hx update equation for j-high, k-high line
            *self
                .hz
                .idxm(i, geometry.num_vox_y - 1, geometry.num_vox_z - 1) += -hay
                * (0.0
                    - self
                        .ez
                        .idx(i, geometry.num_vox_y - 1, geometry.num_vox_z - 1))
                + haz
                    * (0.0
                        - self
                            .ey
                            .idx(i, geometry.num_vox_y - 1, geometry.num_vox_z - 1));
        }

        Ok(())
    }

    fn update_hy(&mut self, geometry: &Geometry, hax: &f64, haz: &f64) -> Result<()> {
        for k in 0..(geometry.num_vox_z - 1) {
            for j in 0..geometry.num_vox_y {
                for i in 0..(geometry.num_vox_x - 1) {
                    // hy update for non i-high and k-high volume
                    *self.hy.idxm(i, j, k) += -haz
                        * (self.ex.idx(i, j, k + 1) - self.ex.idx(i, j, k))
                        + hax * (self.ez.idx(i + 1, j, k) - self.ez.idx(i, j, k));
                }
            }

            for j in 0..geometry.num_vox_y {
                // hy update equation for i-high plane
                *self.hy.idxm(geometry.num_vox_x - 1, j, k) += -haz
                    * (self.ex.idx(geometry.num_vox_x - 1, j, k + 1)
                        - self.ex.idx(geometry.num_vox_x, j, k))
                    + hax * (0.0 - self.ez.idx(geometry.num_vox_x - 1, j, k));
            }
        }

        for j in 0..geometry.num_vox_y {
            for i in 0..(geometry.num_vox_x - 1) {
                // hy update equation for k-high plane
                *self.hy.idxm(i, j, geometry.num_vox_z - 1) += -haz
                    * (0.0 - self.ex.idx(i, j, geometry.num_vox_z - 1))
                    + hax
                        * (self.ez.idx(i + 1, j, geometry.num_vox_z - 1)
                            - self.ez.idx(i, j, geometry.num_vox_z - 1));
            }

            // hy update equation for i-high, k-high line
            *self
                .hy
                .idxm(geometry.num_vox_x - 1, j, geometry.num_vox_z - 1) += -haz
                * (0.0
                    - self
                        .ex
                        .idx(geometry.num_vox_x - 1, j, geometry.num_vox_z - 1))
                + hax
                    * (0.0
                        - self
                            .ez
                            .idx(geometry.num_vox_x - 1, j, geometry.num_vox_z - 1));
        }

        Ok(())
    }

    fn update_hz(&mut self, geometry: &Geometry, hax: &f64, hay: &f64) -> Result<()> {
        for k in 0..geometry.num_vox_z {
            for j in 0..(geometry.num_vox_y - 1) {
                for i in 0..(geometry.num_vox_x - 1) {
                    // hz update for non i-high and j-high volume
                    *self.hz.idxm(i, j, k) += -hax
                        * (self.ey.idx(i + 1, j, k) - self.ey.idx(i, j, k))
                        + hay * (self.ex.idx(i, j + 1, k) - self.ex.idx(i, j, k));
                }

                // hz update equation for i-high plane
                *self.hz.idxm(geometry.num_vox_x - 1, j, k) += -hax
                    * (0.0 - self.ey.idx(geometry.num_vox_x - 1, j, k))
                    + hay
                        * (self.ex.idx(geometry.num_vox_x - 1, j + 1, k)
                            - self.ex.idx(geometry.num_vox_x - 1, j, k));
            }

            for i in 0..(geometry.num_vox_x - 1) {
                // hz update equation for j-high plane
                *self.hz.idxm(i, geometry.num_vox_y - 1, k) += -hax
                    * (self.ey.idx(i + 1, geometry.num_vox_y - 1, k)
                        - self.ey.idx(i, geometry.num_vox_y - 1, k))
                    + hay * (0.0 - self.ex.idx(i, geometry.num_vox_y - 1, k));
            }
        }

        Ok(())
    }
}

#[derive(Debug)]
struct ScalarField {
    field: Vec<f64>,
    row_offset: usize,
    column_offset: usize,
}

impl ScalarField {
    pub fn new(
        initial_value: f64,
        size: usize,
        row_offset: usize,
        column_offset: usize,
    ) -> Result<ScalarField> {
        let field = vec![initial_value; size];
        Ok(ScalarField {
            field,
            row_offset,
            column_offset,
        })
    }

    pub fn idx(&self, i: usize, j: usize, k: usize) -> f64 {
        self.field[i + j * self.row_offset + k * self.row_offset * self.column_offset]
    }

    pub fn idxm(&mut self, i: usize, j: usize, k: usize) -> &mut f64 {
        &mut self.field[i + j * self.row_offset + k * self.row_offset * self.column_offset]
    }
}
