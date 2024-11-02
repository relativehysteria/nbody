use crate::prelude::*;

pub struct Simulation<const DIMENSIONS: usize> {
    pub bodies: Vec<Body<DIMENSIONS>>,
    pub tree: SpatialTree<DIMENSIONS>,
}

impl Simulation<DIMENSIONS> {
    pub fn new(bodies: Vec<Body<DIMENSIONS>>) -> Self {
        let tree = SpatialTree::new();
        Self { bodies, tree }
    }

    pub fn step(&mut self) {
        self.update();
        self.attract();
    }

    pub fn update(&mut self) {
        self.bodies.iter_mut().for_each(|b| b.update(DT));
    }

    pub fn attract(&mut self) {
        self.tree.build(&mut self.bodies);

        for idx in 0..self.bodies.len() {
            self.bodies[idx].acc =
                self.tree.accel(self.bodies[idx].pos, &self.bodies);
        }
    }
}
