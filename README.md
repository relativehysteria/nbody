dimensionally generic n-body gravitational simulation.
uses a modified (more cache-friendly) barnes-hut.

the core logic was mostly done by [_deadlock_](https://github.com/DeadlockCode);
refactoring, generalization and docs by me.

it's `#![no_std]` but requires `alloc::Vec`, there's no visualization or
parallelization, so the simulation can be directly transplanted and executed in
a kernel.

the generalization was done by implementing the simulation over an _n-d tree_
(some call it _hyperoctree_, some _orthree_, etc.) instead of a _quadtree_,
directly translating the logic to be generic over different dimensions.  
this is _not a k-d tree_ mind you; _k-d trees_ split along a dimension whereas
_n-d trees_ (quadtree in 2D, octree in 3D, hyperoctree/n-d tree in further
dimensions) split around a point. they are therefore guaranteed to be "cubical",
having the aspect ratio of 1:1 for each side.
