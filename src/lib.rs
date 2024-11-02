#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(iter_map_windows)]

pub mod rng;
pub mod body;
pub mod vector;
pub mod simulation;
pub mod consts;
pub mod presets;
pub mod spatialtree;
pub mod partition;

pub mod prelude {
    pub use crate::vector::VecN;
    pub use crate::body::Body;
    pub use crate::rng::Rng;
    pub use crate::spatialtree::SpatialTree;
    pub use crate::partition::Partition;
    pub use crate::consts::*;
}
