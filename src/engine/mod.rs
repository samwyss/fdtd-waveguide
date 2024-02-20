//! Engine module
//!
//! Contains all data relevant to the state of the simulation and code needed to evolve said state.
//!
//! Operates on a Geometry (see ./src/geometry/mod.rs).

// constant declarations
const ONE_OVER_TWO: f64 = 1.0 / 2.0;
const TFSF_SRC_END_OFFSET: usize = 3;

// import local modules and cargo crates
use crate::{geometry::Geometry, helpers::write_buf_vec, solver::Config, C_0};
use anyhow::{Ok, Result};
use csv::WriterBuilder;
use std::f64::consts::{PI, TAU};
use std::fs::{create_dir, remove_dir_all, File, OpenOptions};
use std::io::BufWriter;

#[derive(Debug)]
pub struct Engine {
    cur_time: f64,
    ex: ScalarField,
    ey: ScalarField,
    ez: ScalarField,
    hx: ScalarField,
    hy: ScalarField,
    hz: ScalarField,
    hx_wtr: csv::Writer<BufWriter<File>>,
    hy_wtr: csv::Writer<BufWriter<File>>,
    hz_wtr: csv::Writer<BufWriter<File>>,
    ex_wtr: csv::Writer<BufWriter<File>>,
    ey_wtr: csv::Writer<BufWriter<File>>,
    ez_wtr: csv::Writer<BufWriter<File>>,
    snapshot_steps: usize,
}

impl Engine {
    pub fn new(config: &Config, geometry: &Geometry) -> Result<Engine> {
        // assign cur_time to 0 as that is initial state of engine
        let cur_time: f64 = 0.0;

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

        // remove results from previous simulation runs
        let _ = remove_dir_all("./out");

        // make new out directory for simulation results
        create_dir("./out")?;

        // create paths for output files of all field values
        let hx_path = "./out/hx.csv";
        let hy_path = "./out/hy.csv";
        let hz_path = "./out/hz.csv";
        let ex_path = "./out/ex.csv";
        let ey_path = "./out/ey.csv";
        let ez_path = "./out/ez.csv";

        // create file descriptors for all field values
        let hx_file = OpenOptions::new()
            .append(true)
            .create(true)
            .open(&hx_path)?;
        let hy_file = OpenOptions::new()
            .append(true)
            .create(true)
            .open(&hy_path)?;
        let hz_file = OpenOptions::new()
            .append(true)
            .create(true)
            .open(&hz_path)?;
        let ex_file = OpenOptions::new()
            .append(true)
            .create(true)
            .open(&ex_path)?;
        let ey_file = OpenOptions::new()
            .append(true)
            .create(true)
            .open(&ey_path)?;
        let ez_file = OpenOptions::new()
            .append(true)
            .create(true)
            .open(&ez_path)?;

        // create buffered writers for all field values
        let hx_wtr: csv::Writer<BufWriter<File>> = WriterBuilder::new()
            .has_headers(false)
            .from_writer(BufWriter::new(hx_file));
        let hy_wtr: csv::Writer<BufWriter<File>> = WriterBuilder::new()
            .has_headers(false)
            .from_writer(BufWriter::new(hy_file));
        let hz_wtr: csv::Writer<BufWriter<File>> = WriterBuilder::new()
            .has_headers(false)
            .from_writer(BufWriter::new(hz_file));
        let ex_wtr: csv::Writer<BufWriter<File>> = WriterBuilder::new()
            .has_headers(false)
            .from_writer(BufWriter::new(ex_file));
        let ey_wtr: csv::Writer<BufWriter<File>> = WriterBuilder::new()
            .has_headers(false)
            .from_writer(BufWriter::new(ey_file));
        let ez_wtr: csv::Writer<BufWriter<File>> = WriterBuilder::new()
            .has_headers(false)
            .from_writer(BufWriter::new(ez_file));

        // assign snapshot_steps
        let snapshot_steps = config.snapshot_steps;

        Ok(Engine {
            cur_time,
            ex,
            ey,
            ez,
            hx,
            hy,
            hz,
            hx_wtr,
            hy_wtr,
            hz_wtr,
            ex_wtr,
            ey_wtr,
            ez_wtr,
            snapshot_steps,
        })
    }

    pub fn update(&mut self, geometry: &Geometry, target_time: &f64) -> Result<()> {
        // calculate first pass at time step based on Courant–Friedrichs–Lewy stability condition
        let mut dt: f64 = (C_0
            * (geometry.dx_inv.powi(2) + geometry.dy_inv.powi(2) + geometry.dz_inv.powi(2)).sqrt())
        .powi(-1);

        // snap this time step to a specific number of loop iterations using target_time
        let time_steps: usize = (target_time / dt).ceil() as usize;

        // recalculate the time step based on the snapped number of loop iterations
        dt = target_time / time_steps as f64;

        // calculate the number of steps between each snapshot
        let snapshot_mod_steps: usize;
        if self.snapshot_steps >= time_steps {
            snapshot_mod_steps = 1;
        } else {
            snapshot_mod_steps = (time_steps as f64 / self.snapshot_steps as f64).ceil() as usize;
        }

        // pre-process loop constants
        let ea: f64 = (geometry.ep / dt + geometry.sigma / 2.0).powi(-1);
        let eb: f64 = geometry.ep / dt - geometry.sigma / 2.0;
        let hax: f64 = dt * geometry.dx_inv / geometry.mu;
        let hay: f64 = dt * geometry.dy_inv / geometry.mu;
        let haz: f64 = dt * geometry.dz_inv / geometry.mu;

        // log basic time loop stats to stdout
        println!("time steps: {}, dt: {}[s]", time_steps, dt);
        println!("beginning time loop");

        // time loop
        for t in 0..time_steps {
            // update magnetic field
            self.update_h(&geometry, &hax, &hay, &haz)?;

            // update TFSF source
            self.update_tfsf_source(&geometry, &haz)?;

            // update current engine time after magnetic field update
            self.cur_time += ONE_OVER_TWO * dt;

            // update electric field
            self.update_e(&geometry, &ea, &eb)?;

            // update current engine time after electric field update
            self.cur_time += ONE_OVER_TWO * dt;

            if t % snapshot_mod_steps == 0 {
                self.snapshot_fields()?;
            }
        }

        // flush any remaining data in buffers
        print!("flushing remaining data in buffers");
        self.flush_fields()?;

        Ok(())
    }

    fn snapshot_fields(&mut self) -> Result<()> {
        write_buf_vec(&mut self.hx_wtr, &self.hx.field)?;
        write_buf_vec(&mut self.hy_wtr, &self.hy.field)?;
        write_buf_vec(&mut self.hz_wtr, &self.hz.field)?;
        write_buf_vec(&mut self.ex_wtr, &self.ex.field)?;
        write_buf_vec(&mut self.ey_wtr, &self.ey.field)?;
        write_buf_vec(&mut self.ez_wtr, &self.ez.field)?;

        Ok(())
    }

    fn flush_fields(&mut self) -> Result<()> {
        self.hx_wtr.flush()?;
        self.hx_wtr.flush()?;
        self.hy_wtr.flush()?;
        self.ex_wtr.flush()?;
        self.ey_wtr.flush()?;
        self.ez_wtr.flush()?;

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
                    // hx update equation for non j-high, k-high volume
                    *self.hx.idxm(i, j, k) += -hay
                        * (self.ez.idx(i, j + 1, k) - self.ez.idx(i, j, k))
                        + haz * (self.ey.idx(i, j, k + 1) - self.ey.idx(i, j, k));
                }
            }

            for i in 0..geometry.num_vox_x {
                // hx update equation for j-high plane
                *self.hx.idxm(i, geometry.num_vox_y - 1, k) += -hay
                    * (0.0 - self.ez.idx(i, geometry.num_vox_y - 1, k))
                    + haz
                        * (self.ey.idx(i, geometry.num_vox_y - 1, k + 1)
                            - self.ey.idx(i, geometry.num_vox_y - 1, k));
            }
        }

        for j in 0..(geometry.num_vox_y - 1) {
            for i in 0..geometry.num_vox_x {
                // hx update equation for k-high plane
                *self.hx.idxm(i, j, geometry.num_vox_z - 1) += -hay
                    * (self.ez.idx(i, j + 1, geometry.num_vox_z - 1)
                        - self.ez.idx(i, j, geometry.num_vox_z - 1))
                    + haz * (0.0 - self.ey.idx(i, j, geometry.num_vox_z - 1));
            }
        }

        for i in 0..geometry.num_vox_x {
            // hx update equation for j-high, k-high line
            *self
                .hx
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
                    // hy update for non i-high, k-high volume
                    *self.hy.idxm(i, j, k) += -haz
                        * (self.ex.idx(i, j, k + 1) - self.ex.idx(i, j, k))
                        + hax * (self.ez.idx(i + 1, j, k) - self.ez.idx(i, j, k));
                }
            }

            for j in 0..geometry.num_vox_y {
                // hy update equation for i-high plane
                *self.hy.idxm(geometry.num_vox_x - 1, j, k) += -haz
                    * (self.ex.idx(geometry.num_vox_x - 1, j, k + 1)
                        - self.ex.idx(geometry.num_vox_x - 1, j, k))
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
                    // hz update for non i-high, j-high volume
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

            // hz update equation for i-high, j-high line
            *self
                .hz
                .idxm(geometry.num_vox_x - 1, geometry.num_vox_y - 1, k) += -hax
                * (0.0
                    - self
                        .ey
                        .idx(geometry.num_vox_x - 1, geometry.num_vox_y - 1, k))
                + hay
                    * (0.0
                        - self
                            .ex
                            .idx(geometry.num_vox_x - 1, geometry.num_vox_y - 1, k));
        }

        Ok(())
    }

    fn update_e(&mut self, geometry: &Geometry, ea: &f64, eb: &f64) -> Result<()> {
        // update ex
        self.update_ex(geometry, &ea, &eb)?;

        // update ey
        self.update_ey(geometry, &ea, &eb)?;

        // update ez
        self.update_ez(geometry, &ea, &eb)?;

        Ok(())
    }

    fn update_ex(&mut self, geometry: &Geometry, ea: &f64, eb: &f64) -> Result<()> {
        for k in 1..geometry.num_vox_z {
            for j in 1..geometry.num_vox_y {
                // 0 is inside PEC
                for i in 1..geometry.num_vox_x {
                    // ex update equation for all non j-low, k-low volume
                    *self.ex.idxm(i, j, k) = ea
                        * (eb * self.ex.idx(i, j, k)
                            + geometry.dy_inv * (self.hz.idx(i, j, k) - self.hz.idx(i, j - 1, k))
                            - geometry.dz_inv * (self.hy.idx(i, j, k) - self.hy.idx(i, j, k - 1)));
                }
            }
            /* Removed as these field components stretch into PEC
            for i in 0..geometry.num_vox_x {
                // ex update equation for j-low surface
                *self.ex.idxm(i, 0, k) = ea
                    * (eb * self.ex.idx(i, 0, k) + geometry.dy_inv * (self.hz.idx(i, 0, k) - 0.0)
                        - geometry.dz_inv * (self.hy.idx(i, 0, k) - self.hy.idx(i, 0, k - 1)));
            }
            */
        }
        /* Removed as these field components stretch into PEC
        for j in 1..geometry.num_vox_y {
            for i in 0..geometry.num_vox_x {
                // ex update equation for k-low surface
                *self.ex.idxm(i, j, 0) = ea
                    * (eb * self.ex.idx(i, j, 0)
                        + geometry.dy_inv * (self.hz.idx(i, j, 0) - self.hz.idx(i, j - 1, 0))
                        - geometry.dz_inv * (self.hy.idx(i, j, 0) - 0.0));
            }
        }

        for i in 0..geometry.num_vox_x {
            // ex update for j-low, k-low line
            *self.ex.idxm(i, 0, 0) = ea
                * (eb * self.ex.idx(i, 0, 0) + geometry.dy_inv * (self.hz.idx(i, 0, 0) - 0.0)
                    - geometry.dz_inv * (self.hy.idx(i, 0, 0) - 0.0));
        }
        */

        Ok(())
    }

    fn update_ey(&mut self, geometry: &Geometry, ea: &f64, eb: &f64) -> Result<()> {
        for k in 1..geometry.num_vox_z {
            // 0 is inside PEC
            for j in 1..geometry.num_vox_y {
                /* Removed as these field components stretch into PEC
                // ey update equation for i-low surface
                *self.ey.idxm(0, j, k) = ea
                    * (eb * self.ey.idx(0, j, k)
                        + geometry.dz_inv * (self.hx.idx(0, j, k) - self.hx.idx(0, j, k - 1))
                        - geometry.dx_inv * (self.hz.idx(0, j, k) - 0.0));
                */

                for i in 1..geometry.num_vox_x {
                    // ey update equation for all non i-low, k-low volume
                    *self.ey.idxm(i, j, k) = ea
                        * (eb * self.ey.idx(i, j, k)
                            + geometry.dz_inv * (self.hx.idx(i, j, k) - self.hx.idx(i, j, k - 1))
                            - geometry.dx_inv * (self.hz.idx(i, j, k) - self.hz.idx(i - 1, j, k)));
                }
            }
        }
        /* Removed as these field components stretch into PEC
        for j in 0..geometry.num_vox_y {
            // ey update equation for i-low, k-low line
            *self.ey.idxm(0, j, 0) = ea
                * (eb * self.ey.idx(0, j, 0) + geometry.dz_inv * (self.hx.idx(0, j, 0) - 0.0)
                    - geometry.dx_inv * (self.hz.idx(0, j, 0) - 0.0));

            for i in 1..geometry.num_vox_x {
                // ey update equation for k-low surface
                *self.ey.idxm(i, j, 0) = ea
                    * (eb * self.ey.idx(i, j, 0) + geometry.dz_inv * (self.hx.idx(i, j, 0) - 0.0)
                        - geometry.dx_inv * (self.hz.idx(i, j, 0) - self.hz.idx(i - 1, j, 0)));
            }
        }
        */

        Ok(())
    }

    fn update_ez(&mut self, geometry: &Geometry, ea: &f64, eb: &f64) -> Result<()> {
        // 0 is inside PEC
        for k in 1..geometry.num_vox_z {
            /* Removed as these field components stretch into PEC
            // ez update equation for i-low, j-low line
            *self.ez.idxm(0, 0, k) = ea
                * (eb * self.ez.idx(0, 0, k) + geometry.dx_inv * (self.hy.idx(0, 0, k) - 0.0)
                    - geometry.dy_inv * (self.hx.idx(0, 0, k) - 0.0));

            // ez update equation for j-low surface
            for i in 1..geometry.num_vox_x {
                *self.ez.idxm(i, 0, k) = ea
                    * (eb * self.ez.idx(i, 0, k)
                        + geometry.dx_inv * (self.hy.idx(i, 0, k) - self.hy.idx(i - 1, 0, k))
                        - geometry.dy_inv * (self.hx.idx(i, 0, k) - 0.0));
            }
            */

            for j in 1..geometry.num_vox_y {
                /* Removed as these field components stretch into PEC
                // ez update equation for i-low surface
                *self.ez.idxm(0, j, k) = ea
                    * (eb * self.ez.idx(0, j, k) + geometry.dx_inv * (self.hy.idx(0, j, k) - 0.0)
                        - geometry.dy_inv * (self.hx.idx(0, j, k) - self.hx.idx(0, j - 1, k)));
                */
                for i in 1..geometry.num_vox_x {
                    // ez update equation for all non i-low, j-low volume
                    *self.ez.idxm(i, j, k) = ea
                        * (eb * self.ez.idx(i, j, k)
                            + geometry.dx_inv * (self.hy.idx(i, j, k) - self.hy.idx(i - 1, j, k))
                            - geometry.dy_inv * (self.hx.idx(i, j, k) - self.hx.idx(i, j - 1, k)));
                }
            }
        }

        Ok(())
    }

    fn update_tfsf_source(&mut self, geometry: &Geometry, hay: &f64) -> Result<()> {
        // inject Ey field into Hx, Hz as source field
        for j in 1..(geometry.num_vox_y - 1) {
            for i in 1..(geometry.num_vox_x - 1) {
                // scattered field corrections
                *self.hx.idxm(i, j, geometry.num_vox_z - TFSF_SRC_END_OFFSET) -= hay
                    * (TAU * 1.2e9 * self.cur_time).sin()
                    * (PI * i as f64 * geometry.dx / geometry.x_len).sin();

                *self.hz.idxm(i, j, geometry.num_vox_z - TFSF_SRC_END_OFFSET) -= hay
                    * (TAU * 1.2e9 * self.cur_time).sin()
                    * (PI * i as f64 * geometry.dx / geometry.x_len).sin();

                // total field corrections
                *self
                    .hx
                    .idxm(i, j, geometry.num_vox_z - TFSF_SRC_END_OFFSET - 1) += hay
                    * (TAU * 1.2e9 * self.cur_time).sin()
                    * (PI * i as f64 * geometry.dx / geometry.x_len).sin();

                *self
                    .hz
                    .idxm(i, j, geometry.num_vox_z - TFSF_SRC_END_OFFSET - 1) += hay
                    * (TAU * 1.2e9 * self.cur_time).sin()
                    * (PI * i as f64 * geometry.dx / geometry.x_len).sin();
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
