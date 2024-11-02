// use macroquad::prelude::*;
use nbody::rng::Rng;
use nbody::simulation::Simulation;
use nbody::presets::*;

fn rdtsc() -> usize {
    unsafe { core::arch::x86_64::_rdtsc() as usize }
}

fn main() {
    let seed = rdtsc();
    let mut rng = Rng::new(seed);
    println!("SEED: {seed}");

    let bodies = uniform_disc(&mut rng);
    let mut simulation = Simulation::new(bodies);

    loop { simulation.step(); }
}
