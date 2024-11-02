/// Number of spatial dimensions to simulate
pub const DIMENSIONS: usize = 3;

/// Number of bodies to spawn
pub const N_BODIES: usize = 1024;

/// Maximum number of bodies in a [`crate::spatialtree::SpatialTree`] leaf before it splits.
///
/// Setting it to `1` will result in the classical Barnes-Hut behavior where
/// each leaf holds a single body.
pub const LEAF_CAP: usize = 1;

/// Timestep
pub const DT: f64 = 0.05;

/// Treshold for subdivision of nodes within a SpatialTree.
/// Used to determine whether a node's mass should be used as a whole mass of
/// all of the contained children.
pub const THETA: f64 = 1.0;

/// Shortening factor to avoid division by zero in gravity calculations
pub const EPSILON: f64 = 1.0;
