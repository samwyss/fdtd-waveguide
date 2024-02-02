use super::Error;

trait SolverInterface {
    fn new(path_str: &str) -> Result<Solver, Box<dyn Error>>; // Solver struct constructor
}

pub struct Solver {}
