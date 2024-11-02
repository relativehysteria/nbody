use crate::prelude::*;

/// Creates a uniform gravitational disc over a central body
///
/// This spawns bodies only in the first two dimensions, leaving others at 0.
pub fn uniform_disc<const DIMENSIONS: usize>(
        rng: &mut Rng) -> Vec<Body<DIMENSIONS>> {
    // Create the vector to hold the bodies
    let mut bodies = Vec::<Body<DIMENSIONS>>::with_capacity(N_BODIES);

    // Inner radius for the central body, outer radius for the disc
    let inner_radius = 25.0;
    let outer_radius = (N_BODIES as f64).sqrt() * 5.0;

    // Spawn the central body
    let mass = 1e6;
    bodies.push(Body::new(
        mass as f64,
        inner_radius,
        VecN::default(),
        VecN::default()
    ));

    // Generate the rest of the bodies distributed within the disc
    while bodies.len() < bodies.capacity() {
        // Generate an angle for polar coordinates within the circle and its
        // trig identities for positioning
        let angle = (rng.rand() as f64) * core::f64::consts::TAU;
        let (sin, cos) = angle.sin_cos();

        // Calculate the normalized radial distance within the range
        let rad = (inner_radius / outer_radius).powi(2);
        let rad = (rng.rand() as f64) * (1.0 - rad) + rad;

        // TODO: These keep the z coordinate at 0.
        // Read up on Householder transformation and Gram-Schmidt process.

        // Calculate the position for the new body based on polar coords
        let mut pos = VecN::<DIMENSIONS>::default();
        pos[0] = cos;
        pos[1] = sin;
        pos *= outer_radius * rad.sqrt();

        // Tangential velocity perpendicular to the radial position
        let mut vel = VecN::<DIMENSIONS>::default();
        vel[0] = sin;
        vel[1] = -cos;

        // Mass and mass-dependent computed radius
        let mass = 1.0f64;
        let rad  = mass.cbrt();

        // Spawn the body
        bodies.push(Body::new(mass, rad, pos, vel));
    }

    // Sort the bodies based on their distance from the origin, ascending
    bodies
        .sort_by(|a, b| a.pos.magnitude_sq().total_cmp(&b.pos.magnitude_sq()));

    // Accumulate mass and adjust velocities to create circular orbit dynamics
    let mut mass = bodies[0].mass;
    for idx in 1..N_BODIES {
        // Add up the mass for orbital velocity calculation
        mass += bodies[idx].mass;

        // Orbital velocity based on cumulative mass and distance from origin
        let vel = (mass / bodies[idx].pos.magnitude()).sqrt();

        // Scale the velocity to match orbital velocity for stable orbits
        bodies[idx].vel *= vel;
    }

    bodies
}
