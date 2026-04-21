# Proposal: Sherd Collision Repair Step

## Context

We have disabled `use_overlap_penalty` because it was unreliable — its parameters were not size-agnostic, and it failed to detect cases where two sherds occupy the same spatial position but are not geometrically overlapping in the traditional sense (e.g. stacked slightly above each other).

We observed a concrete failure in Pot A where two sherds ended up at nearly identical global coordinates because they were similar in shape and both appeared to fit against the same anchor sherd. Only one of them belongs there — the other should connect elsewhere. The current fitness function has no mechanism to detect or resolve this.

---

## Proposed Addition: Position Collision Repair

We want to add a **repair step** (not a penalty) that runs after `GuidedRepair` and `EnsureSherdCoverage` in the generation loop. The idea:

1. **Build global poses** for the current chromosome by running the same BFS-based pose propagation already used in `BuildOutputsFromSelection` — transform each sherd's centroid into global coordinates using the assembled T_to_root poses.

2. **Detect collisions** — for every pair of sherds (i, j) in the assembled chromosome, check if their global centroids are within a threshold distance. The threshold should be size-agnostic, expressed as a fraction of the smaller of the two sherds' radii (e.g. `0.5 * min(radius[i], radius[j])`). Sherd radii are already precomputed in `shard_radius_`.

3. **Resolve collisions** — for each colliding pair (i, j), identify which gene (active match) is responsible for placing the lower-density sherd at that position. "Lower-density" is defined by the same inlier-based density score already used elsewhere (`inliner / (1 + score)`). Force that gene to either:
   - Switch to the next best alternative candidate in its pair group, or
   - Set to inactive (0) if no alternative exists

4. **Repeat** until no collisions remain or a maximum iteration count is reached (to avoid infinite loops in degenerate cases).

---

## Key Differences from the Old Overlap Penalty

- This is a **hard repair**, not a soft penalty. The GA never produces a chromosome with position collisions — it doesn't need to learn to avoid them.
- The threshold is **size-agnostic** — expressed relative to sherd geometry rather than an absolute constant.
- It operates on **actual assembled global positions**, not on radius estimates that may not reflect true physical extent.
- It targets the **specific gene responsible** for the bad placement, rather than penalizing the chromosome as a whole.

---

## Request

Please review this proposal and give your feedback on:
1. Whether the implementation is feasible given the current `genetic_algorithm.h` structure
2. Any edge cases or failure modes you foresee
3. Whether the BFS pose propagation inside the repair step needs special handling (e.g. disconnected components, identity poses for unreached sherds)
4. Any suggestions to improve the approach before we proceed to implementation

**Do not implement anything yet.** This is a design review only.