# fdtd-waveguide
A 3D Finite-Difference Time-Domain solver for ECE 61800 Project 1.
Please reach out if you are interested in the project / course details or the report.

This codebase, as well as others from this class, are not designed to be user-friendly as they were only intended to be ran by myself for the purpose of a class project.
With that said however all code is reasonably documented to the standards outlined in ECE61800.

## Release Builds
To compile a release build, please follow these steps.
1. [Download Rust](https://www.rust-lang.org/)
2. Set the following environment variable in your preferred terminal to generate optimized SIMD instructions for your CPU.
    
    
    RUSTFLAGS = "-C target-cpu=native"

3. Compile with the following command.
    

    cargo build --release

## Running
Prior to running ensure all model parameters in `config.toml` are set appropriately. 
The model will pass these parameters to the engine. 

**WARNING** It is recommended to set `snapshot_steps` to a small number unless you are certain you need a large number of time-steps. 
It is very easy to inadvertently generate Gigabytes of data per run if you are not careful.

Once compiled, the release build can be found under `./target/release/driver.<>` where `.<>` corresponds to the file extension of executables on your platform.

## Data
Field data is stored in the `./out/` directory in `.csv` format. 
For example the x component of the E field is stored in `./out/ex.csv`.
Rows of these csv files correspond to field values at a given time step in Fortran-style column-ordering.

## Documentation
To generate full project documentation in the form of a Doxygen-style webpage, run the following command in your terminal from the project root.


    cargo doc --document-private-items

After running this command, all generated documentation can be viewed by opening either `./target/doc/waveguide/index.html` or `./target/doc/driver/index.html`.

## Example Results
Below are several labeled figures from the report.
1. High-level waveguide diagram
![Model diagram](./figures/model.png)


2. Comparison of source and transmitted tapered sine fields to that of the analytic cutoff frequency.
![Monochromatic Tapered Sine Wave](./figures/monochromatic-source.png)


3. Field profiles of (top color) y component of electric field (top arrows) x and y components of magnetic field, (bottom color) x component of magnetic field, and (bottom arrows) y and z components of electric field.
![Monochromatic Source Field Profiles](./figures/labeled-monochromatic-source-profile.png)

   
4. Comparison of source and transmitted wideband gaussian pulse fields to that of the analytic cutoff frequency.
![Wideband Gaussian Pulse Comparison](./figures/wideband-spectrum.png)


5. Comparison of Beeswax to Beryllia filled cavity resonator with Mur ABC replaced by PEC Wall
![Comparison of Beeswax to Beryllia](./figures/comp.png)
