use crate::vector::VecN;

/// A body (planet, asteroid, star, black hole, etc.)
#[derive(Clone, Copy, Debug)]
pub struct Body<const DIMENSIONS: usize> {
    /// Mass of the body
    pub mass: f64,

    /// Radius of the body
    pub rad: f64,

    /// Position of the body
    pub pos: VecN<DIMENSIONS>,

    /// Velocity of the body
    pub vel: VecN<DIMENSIONS>,

    /// Acceleration of the body
    pub acc: VecN<DIMENSIONS>,
}

impl<const DIMENSIONS: usize> Body<DIMENSIONS> {
    /// Creates a new body with zero acceleration
    pub fn new(mass: f64, rad: f64, pos: VecN<DIMENSIONS>,
               vel: VecN<DIMENSIONS>) -> Self {
        Self {
            mass,
            rad,
            pos,
            vel,
            acc: VecN::<DIMENSIONS>::default(),
        }
    }

    /// Updates the body's attributes by a timestep `dt`
    pub fn update(&mut self, dt: f64) {
        self.vel += self.acc * dt;
        self.pos += self.vel * dt;
    }
}
