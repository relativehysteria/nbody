/// A deterministic PRNG (xorshift)
pub struct Rng(usize);

impl Rng {
    /// Creates a new RNG
    pub fn new(seed: usize) -> Self {
        let mut rng = Self(seed);

        // First couple of runs to avoid bad seeds
        (0..64).for_each(|_| { rng.rand(); });

        rng
    }

    /// Returns a pseudo-random (predetermined) number
    pub fn rand(&mut self) -> usize {
        let ret = self.0;
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 17;
        self.0 ^= self.0 << 43;
        ret
    }

    /// Returns a pseudo-random (predetermined) number within a given range
    pub fn range(&mut self, min: usize, max: usize) -> usize {
        (self.rand() % (max - min + 1)) + min
    }
}
