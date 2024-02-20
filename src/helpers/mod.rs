//! Helpers Module
//!
//! Contains helper methods and structures that are not explicitly tied to the physics of FDTD but useful nonetheless

// import local crates and modules
use anyhow::{Ok, Result};
use std::fs::File;
use std::io::BufWriter;

pub fn write_buf_vec(
    buffered_writer: &mut csv::Writer<BufWriter<File>>,
    data: &[f64],
) -> Result<()> {
    // write data into buffered_writer
    buffered_writer.write_record(data.iter().map(|x| x.to_string()).collect::<Vec<String>>())?;

    Ok(())
}
