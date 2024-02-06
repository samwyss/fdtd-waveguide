//! fdtd_waveguide library
//!
//! Contains all relevant imports, declarations, etc for simulating a 3D waveguide
//!
//! Think of this as something similar-ish to a header file in C/C++

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
