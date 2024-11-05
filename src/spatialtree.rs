use core::mem::MaybeUninit;
use core::ops::Range;
use crate::prelude::*;


/// Precomputed square of [`THETA`]
const THETA_SQ: f64 = THETA * THETA;

/// Precomputed square of [`EPSILON`]
const EPSILON_SQ: f64 = EPSILON * EPSILON;

/// Precomputed number of subregion splits based on dimensionality
const N_SPLITS: usize = 1 << DIMENSIONS;


/// Representation of an area within a dimensional space
#[derive(Debug, Clone, Copy)]
pub struct Region<const DIMENSIONS: usize> {
    /// Center position of the area
    pub center: VecN<DIMENSIONS>,

    /// Length of one side of the square/cube/...
    pub length: f64,
}

impl<const DIMENSIONS: usize> Region<DIMENSIONS> {
    /// Creates a new space that encompasses all given bodies
    pub fn new_containing(bodies: &[Body<DIMENSIONS>]) -> Self {
        // Initialize minimum and maximum boundaries
        let mut mins = [f64::MAX; DIMENSIONS];
        let mut maxs = [f64::MIN; DIMENSIONS];

        // Update boundaries based on each body's position
        for body in bodies {
            for dimension in 0..DIMENSIONS {
                mins[dimension] = mins[dimension].min(body.pos[dimension]);
                maxs[dimension] = maxs[dimension].max(body.pos[dimension]);
            }
        }

        // Calculate center of `Region` that encompasses all bodies
        let mut center = VecN::<DIMENSIONS>::default();
        mins.iter().zip(maxs.iter())
            .enumerate()
            .map(|(idx, (min, max))| (idx, (min + max)))
            .for_each(|(idx, component)| { center[idx] = component * 0.5 });

        // And calculate the length of each side that encompasses all bodies
        let length = maxs.iter().zip(mins.iter())
            .map(|(max, min)| max - min)
            .fold(0.0f64, |a, b| a.max(b));

        Self { center, length }
    }

    /// Calculates a region by halving the length and adjusting the center
    /// position based on the `region` index
    pub fn into_subregion(mut self, region: usize) -> Self {
        // Halve the size
        self.length *= 0.5;

        // Adjust each dimension's center position based on the corresponding
        // bit in `region`
        for dim in 0..DIMENSIONS {
            // Move the center of a `dimension` based on its `region` bit
            let offset = (region >> dim) & 1 == 1;
            let offset = -0.5 + (offset as usize as f64);
            self.center[dim] += offset * self.length;
        }

        self
    }

    /// Divides the Region into 4 smaller Regions
    pub fn subdivide(&self) -> [Region<DIMENSIONS>; N_SPLITS] {
        // Create the array for regions
        let regions: MaybeUninit<[Region<DIMENSIONS>; N_SPLITS]>
            = MaybeUninit::uninit();
        let mut regions = unsafe { regions.assume_init() };

        // Initialize it and return it
        for i in 0..N_SPLITS { regions[i] = self.into_subregion(i); }
        regions
    }
}

/// Representation of a node within a [`SpatialTree`]
#[derive(Debug, Clone)]
pub struct Node<const DIMENSIONS: usize> {
    /// Index of the first child node, or 0 if this node is a leaf
    pub children: usize,

    /// Index of the next sibling node, 0 if no sibling exists
    pub next: usize,

    /// Position vector for the center of mass
    pub pos: VecN<DIMENSIONS>,

    /// Total mass of bodies within this node
    pub mass: f64,

    /// The region represented by this node
    pub region: Region<DIMENSIONS>,

    /// Range of indices of bodies contained in this node
    pub bodies: Range<usize>,
}

impl<const DIMENSIONS: usize> Node<DIMENSIONS> {
    /// Creates a new `Node`
    pub fn new(next: usize, region: Region<DIMENSIONS>,
               bodies: Range<usize>) -> Self {
        Self {
            children: 0,
            next,
            pos: VecN::default(),
            mass: 0.0,
            region,
            bodies,
        }
    }

    /// Checks whether this node is a leaf (no children)
    pub fn is_leaf(&self) -> bool {
        self.children == 0
    }

    /// Checks whether this node is a branch (has children)
    pub fn is_branch(&self) -> bool {
        self.children != 0
    }

    /// Checks if this node has no mass (is empty)
    pub fn is_empty(&self) -> bool {
        self.mass == 0.0
    }
}

pub struct SpatialTree<const DIMENSIONS: usize> {
    /// List of nodes in the tree
    pub nodes: Vec<Node<DIMENSIONS>>,

    /// Stack to track nodes with children (for mass propagation)
    pub parents: Vec<usize>,
}

impl<const DIMENSIONS: usize> SpatialTree<DIMENSIONS> {
    /// Index for the root node
    const ROOT: usize = 0;

    /// Creates a new `SpatialTree`, preallocating space for [`N_BODIES`].
    pub fn new() -> Self {
        // Estimate maximum nodes based on the number of bodies and leaf cap
        let max_nodes = LEAF_CAP * N_SPLITS;
        let max_nodes = (N_BODIES + LEAF_CAP - 1) / max_nodes;

        Self {
            nodes: Vec::with_capacity(max_nodes),
            parents: Vec::with_capacity(max_nodes),
        }
    }

    /// Clears the nodes and parents within this tree
    pub fn clear(&mut self) {
        self.nodes.clear();
        self.parents.clear();
    }

    /// Subdivides a node if it exceeds the leaf capacity,
    /// redistributing bodies into subregions based on their positions.
    fn subdivide(&mut self, node: usize, bodies: &mut [Body<DIMENSIONS>],
                     range: Range<usize>) {
        // Get the center of current node's region for partitioning
        let center = self.nodes[node].region.center;

        // Loop through each dimension to partition bodies into subregions.
        // For each dimension, it checks each body's position relative to the
        // center of the current node's region, updating the `split` array with
        // the indices where bodies are redistributed based on the partitioning
        // criteria.
        //
        // For example, for the case where `DIMENSIONS == 2`,
        // the code will be evaluated like so:
        //
        // ```
        // let pred = |body: &Body| body.pos[1] < center[1];
        // split[2] = split[0] + bodies[split[0]..split[4]].partition(pred);
        //
        // let pred = |body: &Body| body.pos[0] < center[0];
        // split[1] = split[0] + bodies[split[0]..split[2]].partition(pred);
        // split[3] = split[2] + bodies[split[2]..split[4]].partition(pred);
        // ```

        // Array of the partition start and end indices
        let mut split = [0; N_SPLITS + 1];
        *(split.first_mut().unwrap()) = range.start;
        *(split.last_mut().unwrap())  = range.end;

        for dim in 0..DIMENSIONS {
            // Determine the axis being processed in reverse order
            let d = DIMENSIONS - dim - 1;

            let predicate = |b: &Body<DIMENSIONS>| b.pos[d] < center[d];

            let step = N_SPLITS / 2usize.pow(dim as u32);
            let mut windows = (0..=N_SPLITS).step_by(step)
                .map_windows(|window: &[usize; 2]| (window[0], window[1]));

            for key in (step / 2..N_SPLITS).step_by(step) {
                let (x, y) = windows.next().unwrap();
                split[key] = split[x] + bodies[split[x]..split[y]]
                    .partition(predicate);
            }
        }

        // Record the current node as a parent and set up child nodes
        self.parents.push(node);
        let children = self.nodes.len();
        self.nodes[node].children = children;

        // Calculate the indices for each child node
        let mut nexts = [self.nodes[node].next; N_SPLITS];
        for idx in 0..(N_SPLITS - 1) {
            nexts[idx] = children + idx + 1;
        }

        // Divide the node's region into subregions
        let subregions = self.nodes[node].region.subdivide();

        // Create child nodes for each subregion
        for i in 0..N_SPLITS {
            let body_range = split[i]..split[i + 1];
            self.nodes.push(Node::new(nexts[i], subregions[i], body_range));
        }
    }

    /// Propagates position and mass information from children to parents
    fn propagate(&mut self) {
        // Propagate total mass and positions
        for &node in self.parents.iter().rev() {
            let i = self.nodes[node].children;

            // Sum of the positions of children
            self.nodes[node].pos = (0..N_SPLITS).map(|j| self.nodes[i + j].pos)
                .fold(VecN::<DIMENSIONS>::default(), |a, x| a + x);

            // Sum of the masses of children
            self.nodes[node].mass= (0..N_SPLITS).map(|j| self.nodes[i + j].mass)
                .sum();
        }

        // Normalize the center of mass by dividing by total mass
        for node in &mut self.nodes {
            node.pos /= node.mass.max(f64::MIN_POSITIVE);
        }
    }

    /// Builds the `SpatialTree` from a list of bodies
    pub fn build(&mut self, bodies: &mut [Body<DIMENSIONS>]) {
        // Clear any existing nodes
        self.clear();

        // Initialize the root node to encompass all bodies
        let region = Region::new_containing(bodies);
        self.nodes.push(Node::new(0, region, 0..bodies.len()));

        let mut node = 0;
        while node < self.nodes.len() {
            // Get the range of bodies associated with the current node
            let range = self.nodes[node].bodies.clone();

            // If the number of bodies in the current node exceeds leaf
            // capacity, subdivide. Otherwise accumulate mass and position for
            // each body
            if range.len() > LEAF_CAP {
                self.subdivide(node, bodies, range.clone());
            } else {
                for i in range {
                    // Accumulate weighted position and total mass
                    self.nodes[node].pos  += bodies[i].pos * bodies[i].mass;
                    self.nodes[node].mass += bodies[i].mass;
                }
            }
            // Move to the next node
            node += 1;
        }

        // Propagate the calculated mass and position
        self.propagate();
    }

    /// Calculates the acceleration on a position in space due to gravitational
    /// attraction from `bodies`
    pub fn accel(&self, pos: VecN<DIMENSIONS>,
                 bodies: &[Body<DIMENSIONS>]) -> VecN<DIMENSIONS> {
        // Start with zero acceleration
        let mut acc = VecN::<DIMENSIONS>::default();

        // Helper function to calculate acceleration from mass and position
        let calc_acc = |mass: f64, dist: VecN<DIMENSIONS>|  {
            let dist_sq = dist.magnitude_sq();
            let denom = (dist_sq + EPSILON_SQ) * dist_sq.sqrt();
            dist * (mass / denom).min(f64::MAX)
        };

        // Begin with the root node
        let mut node = Self::ROOT;

        loop {
            let n = self.nodes[node].clone();

            // Distance from node's center of mass to position
            let dist = n.pos - pos;
            let dist_sq = dist.magnitude_sq();

            // Determine if we should use the node's mass as a whole or explore
            // children
            if n.region.length.powi(2) < dist_sq * THETA_SQ {
                // Compute acceleration due to node's mass
                acc += calc_acc(n.mass, dist);

                // Exit if no more nodes to process
                if n.next == 0 { break; }
                node = n.next;
            } else if n.is_leaf() {
                // Compute acceleration due to each body in this leaf node
                acc += n.bodies.map(|i| {
                    let body = &bodies[i];
                    let dist = body.pos - pos;
                    calc_acc(body.mass, dist)
                }).fold(VecN::<DIMENSIONS>::default(), |acc, v| acc + v);

                // Exit of no more nodes to process
                if n.next == 0 { break; }
                node = n.next;
            } else {
                node = n.children;
            }
        }

        acc
    }
}
