use super::geometry::Geometry;
use super::Error;

trait EngineInterface {
    fn new(geometry: Geometry) -> Result<Engine, Box<dyn Error>>; // Engine struct constructor
}

#[derive(Debug)]
pub struct Engine {}
