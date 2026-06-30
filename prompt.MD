# Fix Instructions: Reassembly GA Failure (Shard 6 / Shard 2 Connection)

You are tasked with fixing a bug in our Genetic Algorithm (GA) reassembly pipeline that prevents **Shard 6** from connecting to **Shard 2** correctly on Pot G. The issues have been traced to a feature corruption bug and a missing axis-refinement step during GA iterations.

Please implement the following fixes in the codebase.

---

## Issue 1: Feature Corruption Bug in `main.cpp`

### What is happening:
At the end of each GA iteration in [main.cpp](file:///home/taha/Desktop/structure-from-sherds-pp/main.cpp), the code attempts to recompute matches on the newly reconstructed poses. However, it applies the assembled transformation `T_best[i]` to the shards **before** calling `CalculateFeatureAxisless`:
```cpp
// From main.cpp:
shard = shard_original;
for (int i = 0; i < SHARD_NUMBER; ++i) {
    T_best[i].Output(R, t);
    shard[i].Move(R, t, true); // <--- Moves shard to assembled (tilted/translated) pose
}

for (int i = 0; i < SHARD_NUMBER; ++i) {
    CalculateFeatureAxisless(shard[i], 0); // <--- Computes features on tilted pose
}
```

### Why this is a bug:
`CalculateFeatureAxisless` assumes the shard's symmetry axis is aligned with the global $Z$-axis:
```cpp
// From filter.cpp:
for (int i = 0; i < number_of_points; i++) {
    Vector3d r = { shard.edge_line_.point_(0, i), shard.edge_line_.point_(1, i), 0 };
    Dist(i) = r.norm(); // Norm of X and Y (assumes Z is the axis!)
    Height(i) = shard.edge_line_.point_(2, i);
}
```
If the shard is moved/tilted to its assembled pose, its symmetry axis is no longer the global $Z$-axis. The radial distance and height features become completely corrupted, leading to invalid matching candidates and degrading match counts (e.g., pruned matches drop over iterations).

### How to fix it:
Do **not** compute features on the assembled poses. Features must always be calculated when the shards are in their axis-aligned states. 
* To update features, you must keep the shards in their axis-aligned poses.
* If you have refined the axis of a shard, you must apply the inverse transform to bring the shard back to the aligned pose before recalculating its features.

---

## Issue 2: Missing Axis Refinement (The Chicken-and-Egg Dilemma)

### What is happening:
In the initial LCS run, the estimated symmetry axis of Shard 2 is too inaccurate ("wobbly"), which prevents the LCS algorithm from finding a valid match for the `2-6` edge. In the initial candidate audit, the only matches found for `2-6` have massive tilt/translation errors:
```text
[GT AUDIT] Edge 2-6 (GT=1) Candidates: 2
  Choice 1: inliner= 18 score=0.898 density=9.484 rad=0.4604 trans=13.92 (26 deg tilt)
  Choice 2: inliner= 30 score=0.714 density=17.508 rad=2.2824 trans=7.82 (130 deg tilt)
```
The true high-quality match (LCS size 22, rad error 0.04) **does not exist in the initial candidate set**.

In the original paper's method (SfS++), the assembly is built incrementally, and **global axis refinement (ICP/Ceres solver)** is run at each step. This refinement corrects the wobbly symmetry axis of Shard 2. Once the axis is corrected, the updated descriptors allow the LCS algorithm to discover the hidden, high-quality match between Shard 6 and Shard 2.

Because the GA loop only runs discrete crossover/mutation and lacks continuous axis refinement between generations, the wobbly axes are never corrected. The correct `2-6` match is never generated.

### How to fix it:
To resolve this, integrate intermediate axis refinement into the GA loop between iterations:
1. At the end of a GA iteration, take the current best reconstructed pose configuration.
2. Run the continuous optimizer ([IcpIncGraphAxis](file:///home/taha/Desktop/structure-from-sherds-pp/class/reconstruction.cpp#L1714) / `GraphAxisRefinement`) on this configuration to align and refine the symmetry axes of the assembled shards.
3. Call `T_best[i].InvOut(R, t)` to get the inverse transforms of these refined poses.
4. Move the shards back to their axis-aligned states using these inverse transforms: `shard[i].Move(R, t)`.
5. Recalculate the axis-aligned features using `CalculateFeatureAxisless`.
6. Run `FeatureComp` and `PairwisePruning` to regenerate the candidates list `LCS_out` for the next GA iteration. This will generate the true `2-6` match.
