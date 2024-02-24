//! Engine module
//!
//! Contains all data relevant to the state of the simulation and code needed to evolve said state.
//!
//! Operates on a Geometry (see ./src/geometry/mod.rs) and a Config (see ./src/solver/mod.rs).

// import local modules and cargo crates
use crate::{geometry::Geometry, solver::Config, C_0};
use anyhow::{Ok, Result};
use csv::{Writer, WriterBuilder};
use std::f64::consts::{E, PI, TAU};
use std::fs::{create_dir, remove_dir_all, File, OpenOptions};
use std::io::BufWriter;

// constant declarations
const ONE_OVER_TWO: f64 = 1.0 / 2.0;

/// ScalarField struct
///
/// Provides a no cost abstraction of a Vec<f64> to function as a 3D scalar field using Fortran-style ordering as this intuitively simpler than C-style in my opinion
#[derive(Debug, Clone)]
struct ScalarField {
    field: Vec<f64>, // [implied unit] contains all scalar field data of implied unit (e.g. A/m, V/m, etc.)
    row_offset: usize, // [] row offset into 1D Vec
    column_offset: usize, // [] column offset into 1D Vec
}

impl ScalarField {
    /// ScalarField constructor
    ///
    /// # Arguments
    ///
    /// - `initial_value` [implied unit] initial value for all points in ScalarField
    /// - `size` [] total number of discrete locations in field
    /// - `row_offset` [] row offset into 1D Vec
    /// - `column_offset` [] column offset into 1D Vec
    ///
    /// # Returns
    ///
    /// `Result<ScalarField>`
    ///
    /// # Errors
    ///
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

    /// returns an aliased value at position [i,j,k] in ScalarField
    ///
    /// used to get an immutable value of a ScalarField at position [i,j,k]
    ///
    /// # Arguments
    ///
    /// -`i` [] i-index into ScalarField using Fortran-style ordering
    /// -`j` [] j-index into ScalarField using Fortran-style ordering
    /// -`k` [] k-index into ScalarField using Fortran-style ordering
    ///
    /// # Returns
    ///
    /// `f64`
    ///
    /// # Errors
    ///
    /// panics if [i,j,k] is out of bounds. this was not error handled properly to reduce the number of checks as this function is used very heavily
    pub fn idx(&self, i: usize, j: usize, k: usize) -> f64 {
        self.field[i + j * self.row_offset + k * self.row_offset * self.column_offset]
    }

    /// returns a mutable reference to value at position [i,j,k] in ScalarField
    ///
    /// used to get a mutable reference to a value of a ScalarField at position [i,j,k]
    ///
    /// # Arguments
    ///
    /// -`i` [] i-index into ScalarField using Fortran-style ordering
    /// -`j` [] j-index into ScalarField using Fortran-style ordering
    /// -`k` [] k-index into ScalarField using Fortran-style ordering
    ///
    /// # Returns
    ///
    /// `&mut f64`
    ///
    /// # Errors
    ///
    /// panics if [i,j,k] is out of bounds. this was not error handled properly to reduce the number of checks as this function is used very heavily
    pub fn idxm(&mut self, i: usize, j: usize, k: usize) -> &mut f64 {
        &mut self.field[i + j * self.row_offset + k * self.row_offset * self.column_offset]
    }
}

/// Engine struct
///
/// contains all code and data necessary to evolve field states to a given time
#[derive(Debug)]
pub struct Engine {
    cur_time: f64,                   // [s] the current time of the simulation
    ex: ScalarField,                 // [V/m] Ex ScalarField
    ey: ScalarField,                 // [V/m] Ey ScalarField
    ez: ScalarField,                 // [V/m] Ez ScalarField
    hx: ScalarField,                 // [A/m] Hx ScalarField
    hy: ScalarField,                 // [A/m] Hy ScalarField
    hz: ScalarField,                 // [A/m] Hz ScalarField
    hx_wtr: Writer<BufWriter<File>>, // [] Ex Field Buffered CSV Writer
    hy_wtr: Writer<BufWriter<File>>, // [] Ey Field Buffered CSV Writer
    hz_wtr: Writer<BufWriter<File>>, // [] Ez Field Buffered CSV Writer
    ex_wtr: Writer<BufWriter<File>>, // [] Hx Field Buffered CSV Writer
    ey_wtr: Writer<BufWriter<File>>, // [] Hy Field Buffered CSV Writer
    ez_wtr: Writer<BufWriter<File>>, // [] Hz Field Buffered CSV Writer
    snapshot_steps: usize,           // [] number of steps to snapshot
}

impl Engine {
    /// Engine constructor
    ///
    /// # Arguments
    ///
    /// `config` solver::Config struct
    /// `geometry` geometry::Geometry struct
    ///
    /// # Returns
    ///
    /// `Result<Engine>`
    ///
    /// # Errors
    /// - any `ScalarField` constructors error
    /// -`fs::create_dir()` is unable to create directory "./out"
    /// - creation of any file descriptors fail
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

        // create buffered writers of desired size for all field values
        let hx_wtr: Writer<BufWriter<File>> =
            WriterBuilder::new()
                .has_headers(false)
                .from_writer(BufWriter::with_capacity(
                    8 * geometry.num_vox * config.buffered_snapshots,
                    hx_file,
                ));
        let hy_wtr: Writer<BufWriter<File>> =
            WriterBuilder::new()
                .has_headers(false)
                .from_writer(BufWriter::with_capacity(
                    8 * geometry.num_vox * config.buffered_snapshots,
                    hy_file,
                ));
        let hz_wtr: Writer<BufWriter<File>> =
            WriterBuilder::new()
                .has_headers(false)
                .from_writer(BufWriter::with_capacity(
                    8 * geometry.num_vox * config.buffered_snapshots,
                    hz_file,
                ));
        let ex_wtr: Writer<BufWriter<File>> =
            WriterBuilder::new()
                .has_headers(false)
                .from_writer(BufWriter::with_capacity(
                    8 * geometry.num_vox * config.buffered_snapshots,
                    ex_file,
                ));
        let ey_wtr: Writer<BufWriter<File>> =
            WriterBuilder::new()
                .has_headers(false)
                .from_writer(BufWriter::with_capacity(
                    8 * geometry.num_vox * config.buffered_snapshots,
                    ey_file,
                ));
        let ez_wtr: Writer<BufWriter<File>> =
            WriterBuilder::new()
                .has_headers(false)
                .from_writer(BufWriter::with_capacity(
                    8 * geometry.num_vox * config.buffered_snapshots,
                    ez_file,
                ));

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

    /// updates all fields to time specified in config.end_time
    ///
    /// # Arguments
    ///
    /// `config` solver::Config struct
    /// `geometry` geometry::Geometry struct
    ///
    /// # Returns
    ///
    /// `Result<()>`
    ///
    /// # Errors
    ///
    /// - any `snapshot_fields()` call errors
    /// - `flush_fields()` call errors
    ///
    pub fn update(&mut self, config: &Config, geometry: &Geometry) -> Result<()> {
        // calculate first pass at time step based on Courant–Friedrichs–Lewy stability condition
        let mut dt: f64 = (C_0
            * (geometry.dx_inv.powi(2) + geometry.dy_inv.powi(2) + geometry.dz_inv.powi(2)).sqrt())
        .powi(-1);

        // snap this time step to a specific number of loop iterations using target_time
        let time_steps: usize = (config.end_time / dt).ceil() as usize;

        // recalculate the time step based on the snapped number of loop iterations
        dt = config.end_time / time_steps as f64;

        // calculate the number of steps between each snapshot
        let snapshot_mod_steps: usize;
        if self.snapshot_steps >= time_steps {
            snapshot_mod_steps = 1;
        } else {
            snapshot_mod_steps = (time_steps as f64 / self.snapshot_steps as f64).ceil() as usize;
        }

        // pre-process loop constants
        let ea: f64 = (geometry.ep / dt + geometry.sigma / 2.0).powi(-1); // electric field 'alpha' update coefficient
        let eb: f64 = geometry.ep / dt - geometry.sigma / 2.0; // electric field 'beta' update coefficient
        let hax: f64 = dt * geometry.dx_inv / geometry.mu; // Hx update coefficient
        let hay: f64 = dt * geometry.dy_inv / geometry.mu; // Hy update coefficient
        let haz: f64 = dt * geometry.dz_inv / geometry.mu; // Hz update coefficient

        // pre-process Mur ABC constants
        let c_mur = C_0 / (geometry.mu_r * geometry.ep_r).sqrt(); // [m/s] speed of light in dielectric inside waveguide
        let abc1 = (c_mur * dt - geometry.dz) / (c_mur * dt + geometry.dz); // 1st Mur ABC update coefficient
        let abc2 = (2.0 * geometry.dz) / (c_mur * dt + geometry.dz); // 2nd Mur ABC update coefficient
        let abc3 = ((c_mur * dt).powi(2) * geometry.dz)
            / (2.0 * (geometry.dy).powi(2) * c_mur * dt + geometry.dz); // 3rd Mur ABC update coefficient
        let abc4 = ((c_mur * dt).powi(2) * geometry.dz)
            / (2.0 * (geometry.dx).powi(2) * c_mur * dt + geometry.dz); // 4th Mur ABC update coefficient

        // create clones of self.ey at n and n-1 for use in Mur ABC
        let mut ey_n: ScalarField = self.ey.clone();
        let mut ey_n_1: ScalarField = ey_n.clone();

        // log basic time loop stats to stdout //TODO this could be removed
        println!("time steps: {}, dt: {}[s]", time_steps, dt);
        println!("\nbeginning time loop");

        // time loop
        for t in 0..time_steps {
            // update magnetic field
            self.update_h(&geometry, &hax, &hay, &haz);

            // update TFSF source by correcting the curl of Hx and Hz for an Ey driver source
            self.update_tfsf_source(&config, &geometry, &hay);

            // update current engine time after magnetic field update
            self.cur_time += ONE_OVER_TWO * dt;

            // update electric field
            self.update_e(&geometry, &ea, &eb);

            // update Mur ABC
            self.update_mur_abc(&geometry, &abc1, &abc2, &abc3, &abc4, &ey_n, &ey_n_1);

            // update clones of Ey for Mur ABC
            // TODO this can be improved by only keeping the desired array slice(s) instead of the entire field but this was convenient
            ey_n_1 = ey_n.clone();
            ey_n = self.ey.clone();

            // update current engine time after electric field update
            self.cur_time += ONE_OVER_TWO * dt;

            // conditionally snapshot all field values
            if t % snapshot_mod_steps == 0 {
                println!(
                    "took snapshot at t={}[s], step {}/{}",
                    self.cur_time, t, time_steps
                );
                self.snapshot_fields()?;
            }
        }

        // flush any remaining data in buffers
        print!("flushing remaining data in buffers");
        self.flush_fields()?;

        Ok(())
    }

    /// snapshots all field values
    ///
    /// # Arguments
    ///
    /// `&mut self` mutable reference to Engine struct
    ///
    /// # Returns
    ///
    /// `Result<()>`
    ///
    /// # Errors
    ///
    /// - `write_buf_vec()` calls fail for any field
    ///
    fn snapshot_fields(&mut self) -> Result<()> {
        Engine::write_buf_vec(&mut self.hx_wtr, &self.hx.field)?;
        Engine::write_buf_vec(&mut self.hy_wtr, &self.hy.field)?;
        Engine::write_buf_vec(&mut self.hz_wtr, &self.hz.field)?;
        Engine::write_buf_vec(&mut self.ex_wtr, &self.ex.field)?;
        Engine::write_buf_vec(&mut self.ey_wtr, &self.ey.field)?;
        Engine::write_buf_vec(&mut self.ez_wtr, &self.ez.field)?;

        Ok(())
    }

    /// flushes remaining data in buffers to disk if present
    ///
    /// # Arguments
    ///
    /// `&mut self` mutable reference to Engine struct
    ///
    /// # Returns
    ///
    /// `Result<()>`
    ///
    /// # Errors
    ///
    /// - `wtr.flush()` calls fail for any field
    ///
    fn flush_fields(&mut self) -> Result<()> {
        self.hx_wtr.flush()?;
        self.hy_wtr.flush()?;
        self.hz_wtr.flush()?;
        self.ex_wtr.flush()?;
        self.ey_wtr.flush()?;
        self.ez_wtr.flush()?;

        Ok(())
    }

    /// writes data found in &[f64] into a buffered writer
    ///
    /// # Arguments
    ///
    /// `buffered_writer` &mut Writer<BufWriter<File>> mutable reference to a buffered writer
    /// `data` &[f64] array slice to be written to disk
    ///
    /// # Returns
    ///
    /// `Result<()>`
    ///
    /// # Errors
    ///
    /// - `buffered_writer()` call fails
    ///
    fn write_buf_vec(buffered_writer: &mut Writer<BufWriter<File>>, data: &[f64]) -> Result<()> {
        // write data into buffered_writer
        buffered_writer
            .write_record(data.iter().map(|x| x.to_string()).collect::<Vec<String>>())?;

        Ok(())
    }

    /// updates all magnetic field components using standard 3D Yee scheme
    ///
    /// # Arguments
    ///
    /// `&mut self` &mut Engine mutable reference to Engine struct
    /// `geometry` &Geometry reference to a Geometry struct
    /// `hax` &f64 reference to Hx update coefficient
    /// `hay` &f64 reference to Hy update coefficient
    /// `haz` &f64 reference to Hz update coefficient
    ///
    /// # Returns
    ///
    /// # Errors
    ///
    /// - Will not error, will panic at out of bounds accesses in `ScalarField.idx()` or `ScalarField.idxm()` calls inside of any field update equations
    ///
    fn update_h(&mut self, geometry: &Geometry, hax: &f64, hay: &f64, haz: &f64) {
        // update hx
        self.update_hx(geometry, hay, haz);

        // update hy
        self.update_hy(geometry, hax, haz);

        // update hz
        self.update_hz(geometry, hax, hay);
    }

    /// updates all Hx components using standard 3D Yee scheme
    ///
    /// # Arguments
    ///
    /// `&mut self` &mut Engine mutable reference to Engine struct
    /// `geometry` &Geometry reference to a Geometry struct
    /// `hay` &f64 reference to Hy update coefficient
    /// `haz` &f64 reference to Hz update coefficient
    ///
    /// # Returns
    ///
    /// # Errors
    ///
    /// - Will not error, will panic at out of bounds accesses in `ScalarField.idx()` or `ScalarField.idxm()` calls
    ///
    fn update_hx(&mut self, geometry: &Geometry, hay: &f64, haz: &f64) {
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
    }

    /// updates all Hy components using standard 3D Yee scheme
    ///
    /// # Arguments
    ///
    /// `&mut self` &mut Engine mutable reference to Engine struct
    /// `geometry` &Geometry reference to a Geometry struct
    /// `hax` &f64 reference to Hx update coefficient
    /// `haz` &f64 reference to Hz update coefficient
    ///
    /// # Returns
    ///
    /// # Errors
    ///
    /// - Will not error, will panic at out of bounds accesses in `ScalarField.idx()` or `ScalarField.idxm()` calls
    ///
    fn update_hy(&mut self, geometry: &Geometry, hax: &f64, haz: &f64) {
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
    }

    // updates all Hz components using standard 3D Yee scheme
    ///
    /// # Arguments
    ///
    /// `&mut self` &mut Engine mutable reference to Engine struct
    /// `geometry` &Geometry reference to a Geometry struct
    /// `hax` &f64 reference to Hx update coefficient
    /// `hay` &f64 reference to Hy update coefficient
    ///
    /// # Returns
    ///
    /// # Errors
    ///
    /// - Will not error, will panic at out of bounds accesses in `ScalarField.idx()` or `ScalarField.idxm()` calls
    ///
    fn update_hz(&mut self, geometry: &Geometry, hax: &f64, hay: &f64) {
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
    }

    /// updates all electric field components using standard 3D Yee scheme
    ///
    /// # Arguments
    ///
    /// `&mut self` &mut Engine mutable reference to Engine struct
    /// `geometry` &Geometry reference to a Geometry struct
    /// `ea` &f64 reference electric field 'alpha' update coefficient
    /// `eb` &f64 reference electric field 'beta' update coefficient
    ///
    /// # Returns
    ///
    /// # Errors
    ///
    /// - Will not error, will panic at out of bounds accesses in `ScalarField.idx()` or `ScalarField.idxm()` calls inside of any field update equations
    ///
    fn update_e(&mut self, geometry: &Geometry, ea: &f64, eb: &f64) {
        // update ex
        self.update_ex(geometry, &ea, &eb);

        // update ey
        self.update_ey(geometry, &ea, &eb);

        // update ez
        self.update_ez(geometry, &ea, &eb);
    }

    // updates all Ex components using standard 3D Yee scheme
    ///
    /// # Arguments
    ///
    /// `&mut self` &mut Engine mutable reference to Engine struct
    /// `geometry` &Geometry reference to a Geometry struct
    /// `ea` &f64 reference electric field 'alpha' update coefficient
    /// `eb` &f64 reference electric field 'beta' update coefficient
    ///
    /// # Returns
    ///
    /// # Errors
    ///
    /// - Will not error, will panic at out of bounds accesses in `ScalarField.idx()` or `ScalarField.idxm()` calls
    ///
    fn update_ex(&mut self, geometry: &Geometry, ea: &f64, eb: &f64) {
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
    }

    // updates all Ey components using standard 3D Yee scheme
    ///
    /// # Arguments
    ///
    /// `&mut self` &mut Engine mutable reference to Engine struct
    /// `geometry` &Geometry reference to a Geometry struct
    /// `ea` &f64 reference electric field 'alpha' update coefficient
    /// `eb` &f64 reference electric field 'beta' update coefficient
    ///
    /// # Returns
    ///
    /// # Errors
    ///
    /// - Will not error, will panic at out of bounds accesses in `ScalarField.idx()` or `ScalarField.idxm()` calls
    ///
    fn update_ey(&mut self, geometry: &Geometry, ea: &f64, eb: &f64) {
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
    }

    // updates all Ez components using standard 3D Yee scheme
    ///
    /// # Arguments
    ///
    /// `&mut self` &mut Engine mutable reference to Engine struct
    /// `geometry` &Geometry reference to a Geometry struct
    /// `ea` &f64 reference electric field 'alpha' update coefficient
    /// `eb` &f64 reference electric field 'beta' update coefficient
    ///
    /// # Returns
    ///
    /// # Errors
    ///
    /// - Will not error, will panic at out of bounds accesses in `ScalarField.idx()` or `ScalarField.idxm()` calls
    ///
    fn update_ez(&mut self, geometry: &Geometry, ea: &f64, eb: &f64) {
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
    }

    // injects an Ey source field into Hx, Hz using the TF/SF method with PEC boundary condition for waveguide source at z-high plane
    ///
    /// # Arguments
    ///
    /// `&mut self` &mut Engine mutable reference to Engine struct
    /// `config` &Config reference to a Config struct
    /// `geometry` &Geometry reference to a Geometry struct
    /// `hay` &f64 reference to Hy update coefficient
    ///
    /// # Returns
    ///
    /// # Errors
    ///
    /// - Will not error, will panic at out of bounds accesses in `ScalarField.idx()` or `ScalarField.idxm()` calls
    ///
    fn update_tfsf_source(&mut self, config: &Config, geometry: &Geometry, hay: &f64) {
        // inject Ey field into Hx, Hz as source field
        for j in 1..(geometry.num_vox_y - 1) {
            for i in 1..(geometry.num_vox_x - 1) {
                // scattered field corrections
                *self.hx.idxm(i, j, geometry.num_vox_z - 4) -= hay
                    * self.tapered_sin(config)
                    * (PI * i as f64 * geometry.dx / geometry.x_len).sin();

                // total field corrections
                *self.hx.idxm(i, j, geometry.num_vox_z - 5) += hay
                    * self.tapered_sin(config)
                    * (PI * i as f64 * geometry.dx / geometry.x_len).sin();
            }
        }
    }

    /// computes the value of a tapered sin wave using parameters defined in Config struct
    ///
    /// # Arguments
    ///
    /// `self` &Engine reference to Engine struct
    /// `config` &Config reference to a Config struct
    ///
    /// # Returns
    ///
    /// # Errors
    ///
    fn tapered_sin(&self, config: &Config) -> f64 {
        (1.0 - E.powf(-(self.cur_time - config.delay_time) / (config.ramp_time)))
            * (TAU * config.frequency * self.cur_time).sin()
    }

    // updates Ey field points on z-low plane using Mur ABC scheme
    ///
    /// # Arguments
    ///
    /// `&mut self` &mut Engine mutable reference to Engine structt
    /// `geometry` &Geometry reference to a Geometry struct
    /// `abc1` &f64 reference to 1st Mur ABC update coefficient
    /// `abc2` &f64 reference to 2nd Mur ABC update coefficient
    /// `abc3` &f64 reference to 3rd Mur ABC update coefficient
    /// `abc4` &f64 reference to 4th Mur ABC update coefficient
    /// `ey_n` &ScalarField reference to cloned Ey field at t=n
    /// `ey_n_1` &ScalarField reference to cloned Ey field at t=n-1
    ///
    /// # Returns
    ///
    /// # Errors
    ///
    /// - Will not error, will panic at out of bounds accesses in `ScalarField.idx()` or `ScalarField.idxm()` calls
    ///
    fn update_mur_abc(
        &mut self,
        geometry: &Geometry,
        abc1: &f64,
        abc2: &f64,
        abc3: &f64,
        abc4: &f64,
        ey_n: &ScalarField,
        ey_n_1: &ScalarField,
    ) {
        for j in 1..geometry.num_vox_y {
            for i in 1..geometry.num_vox_x {
                *self.ey.idxm(i, j, 0) = -ey_n_1.idx(i, j, 1)
                    + abc1 * (self.ey.idx(i, j, 1) + ey_n_1.idx(i, j, 0))
                    + abc2 * (ey_n.idx(i, j, 0) + ey_n.idx(i, j, 1))
                    + abc3
                        * (ey_n.idx(i, j + 1, 0) - 2.0 * ey_n.idx(i, j, 0)
                            + ey_n.idx(i, j - 1, 0)
                            + ey_n.idx(i, j + 1, 1)
                            - 2.0 * ey_n.idx(i, j, 1)
                            + ey_n.idx(i, j - 1, 0))
                    + abc4
                        * (ey_n.idx(i + 1, j, 0) - 2.0 * ey_n.idx(i, j, 0)
                            + ey_n.idx(i - 1, j, 0)
                            + ey_n.idx(i + 1, j, 1)
                            - 2.0 * ey_n.idx(i, j, 1)
                            + ey_n.idx(i - 1, j, 1));
            }
        }
    }
}
