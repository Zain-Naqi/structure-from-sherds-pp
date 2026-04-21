# Genetic Algorithm Deep Analysis (structure-from-sherds-pp)

This document explains the GA implemented in class/genetic_algorithm.h in plain language.

It focuses on:
- what each stage is doing,
- how the chromosome is encoded,
- how each fitness term works,
- and what likely issues may be causing sub-optimal reconstruction.

---

## 1) High-level goal

The GA is trying to choose one candidate match per shard-pair group (or choose no match), so that all selected matches together form a globally consistent pot assembly.

In simple terms:
- Local match quality is good (high inlier, low score),
- Global graph structure is good (connected, low contradiction),
- Geometry is plausible (low overlap, low residual inconsistency).

Main entry points:
- Run(): class/genetic_algorithm.h:120
- EvaluateFitness(): class/genetic_algorithm.h:406
- BuildOutputsFromSelection(): class/genetic_algorithm.h:1111

---

## 2) Representation (Chromosome)

A chromosome stores genes as integers:
- genes[g] = 0 means this shard-pair group is inactive,
- genes[g] = k (k >= 1) means pick the k-th candidate in that group.

Data structure:
- struct Chromosome: class/genetic_algorithm.h:40

How groups are built:
- BuildPairGroups() groups all LCS candidates by unordered shard pair (min(x,y), max(x,y)).
- So each gene controls exactly one shard pair group.
- BuildPairGroups(): class/genetic_algorithm.h:1277

---

## 3) GA lifecycle step-by-step

### 3.1 Initialization

InitializePopulation(): class/genetic_algorithm.h:352

For each chromosome:
1. Sample each gene with SampleGroupChoice().
2. Run EnsureSherdCoverage() so active shards are not left isolated with degree 0 (where possible).

Sampling logic (SampleGroupChoice(), class/genetic_algorithm.h:1589):
- with probability kInitialPairInactiveRate (0.1): choose inactive (0),
- otherwise choose an active candidate,
- with probability kGreedyPairChoiceRate (0.7): choose best local density,
- else weighted random by density.

### 3.2 Evaluate and sort

EvaluatePopulation(): class/genetic_algorithm.h:373

- Each chromosome gets fitness from EvaluateFitness().
- Population is sorted descending by fitness.

### 3.3 Survivor selection (elitism)

In Run(), top kElitismCount (2) chromosomes are copied directly to next generation.

This preserves the current best solutions.

### 3.4 Parent selection

TournamentSelect(): class/genetic_algorithm.h:939

- Tournament size is 3 (best out of 3 random population picks).
- This is moderate pressure toward fitter parents.

### 3.5 Crossover

Crossover(): class/genetic_algorithm.h:957

- One-point crossover.
- Child takes prefix from parent1 and suffix from parent2.

### 3.6 Mutation

Mutate(): class/genetic_algorithm.h:983

For each gene independently:
- mutate with probability kMutationRate (0.08),
- small chance symmetry-flip branch (kSymmetryFlipRate = 0.15),
- otherwise resample from SampleGroupChoice(), with attempts to avoid no-op.

### 3.7 Local repair

GuidedRepair(): class/genetic_algorithm.h:1057

- Triggered only for one selected mutated gene.
- Tries a few alternatives (kGuidedRepairTrials = 3) and keeps the one with best full fitness.
- This is a tiny local search around mutation.

### 3.8 Coverage repair

EnsureSherdCoverage(): class/genetic_algorithm.h:1404

- If a valid shard has zero selected incident edges, force-add best local connection from incident groups.
- This pushes solutions away from isolated shards.

### 3.9 Diversity enforcement

EnforceDiversity(): class/genetic_algorithm.h:1563

- Compares chromosome pairs by Hamming distance.
- If too similar, replace one using GenerateBiasedReplacement().

GenerateBiasedReplacement(): class/genetic_algorithm.h:1503
- Keep top active genes by inliers from a survivor,
- resample others,
- optional guided repair,
- then coverage repair.

---

## 4) Fitness function in detail

EvaluateFitness(): class/genetic_algorithm.h:406

The fitness is a sum of rewards minus penalties.

General structure:
- Build active edge set from genes,
- compute graph/geometry metrics,
- combine with weights.

### 4.1 Inlier-score reward

Flag: use_inlier_score (default true)

For each active selected candidate:
- density = inlier / (1 + score)
- add density to fitness

Code anchor: class/genetic_algorithm.h:476

Meaning:
- prefers high inlier matches,
- prefers lower score,
- gives a local evidence reward.

### 4.2 Pair-choice penalty

Flag: use_pair_choice_penalty (default true)

For each active gene:
- compare selected density to best possible density in that group,
- compute gap ratio = (best - selected) / best,
- accumulate positive gaps,
- penalty = kPairChoicePenaltyWeight * average_gap.

Code anchors:
- gap accumulation: class/genetic_algorithm.h:485
- final penalty: class/genetic_algorithm.h:551

Meaning:
- discourages selecting a weak candidate when a stronger local alternative exists.

### 4.3 Neighbor penalty

Flag: use_neighbor_penalty (default true)

- Count selected degree of each shard.
- Penalize degree above kMaxNeighbors.

Code anchor: class/genetic_algorithm.h:507

Meaning:
- intended to discourage over-connected shards.

### 4.4 Active pair range penalty

Flag: use_active_pair_range_penalty (default false in this file)

- Compute valid_group_count.
- min_active = ceil(kMinActivePairRatio * valid_group_count)
- max_active = floor(kMaxActivePairRatio * valid_group_count)
- enforce min active also at least spanning_min = valid_shard_count - 1.
- Penalize active_pair_count outside [min_active, max_active].

Code anchors:
- gate and bounds: class/genetic_algorithm.h:524
- constant values: class/genetic_algorithm.h:1663, class/genetic_algorithm.h:1664

Meaning:
- regularizes sparsity/density of selected graph.

### 4.5 Cycle consistency penalty (translation triangle check)

Flag: use_cycle_penalty (default true)

- Build active_transforms map for selected edges.
- For each shard triplet (a,b,c), if all 3 edges exist:
  - compare composed translation (ab * bc) against direct ac translation,
  - if error > threshold, add error.
- total penalty = 2.0 * cycle_penalty.

Code anchors:
- active_transforms: class/genetic_algorithm.h:559
- cycle loop: class/genetic_algorithm.h:587
- apply penalty: class/genetic_algorithm.h:639

Meaning:
- tries to reduce translational inconsistency around loops.

### 4.6 Edge residual penalty (global pose consistency)

Flag: use_edge_residual (default true)

- Build pair_edges with orientation normalized to (a<b).
- Build adjacency pose graph.
- For each root shard:
  - BFS to estimate global poses,
  - for every reachable selected edge, compare predicted edge translation to measured edge translation,
  - penalize squared excess over threshold.
- Average over roots, weight by kEdgeResidualPenalty.

Code anchors:
- pair_edges build: class/genetic_algorithm.h:655
- root loop: class/genetic_algorithm.h:715
- translation residual: class/genetic_algorithm.h:747
- apply weighted penalty: class/genetic_algorithm.h:781

Meaning:
- encourages selected edges to agree with a single global pose assignment.

### 4.7 Rotational residual penalty

Flag: use_rot_residual (default true)

- Same global pose loop as above.
- Compute rotational discrepancy using matrix log norm.
- Penalize excess over kEdgeRotResidualThreshold.
- Weight by kEdgeRotResidualPenalty.

Code anchors:
- rotational residual: class/genetic_algorithm.h:756
- apply weighted penalty: class/genetic_algorithm.h:788

Meaning:
- discourages rotational contradictions between selected edges.

### 4.8 Overlap penalty

Flag: use_overlap_penalty (default true)

- From first valid shard, BFS one connected component in pose graph.
- Transform shard centroids to global frame.
- For each pair of placed shards:
  - if centroid distance < scaled sum of radii, count as overlap violation.
- penalty = overlap_violations * kOverlapPenalty.

Code anchors:
- overlap block: class/genetic_algorithm.h:797
- violation check: class/genetic_algorithm.h:834

Meaning:
- physically discourages impossible occupancy/symmetry collapse.

### 4.9 Connectivity reward and component penalty

Flag: use_connectivity_reward (default true)

- Compute connected components on selected adjacency.
- reward = kConnectivityReward * largest_component
- penalty = kConnectivityComponentPenalty * (num_components - 1) when disconnected.

Code anchors:
- AnalyzeConnectedComponents(): class/genetic_algorithm.h:878
- reward/penalty usage: class/genetic_algorithm.h:850

Meaning:
- pushes assembly toward one large connected structure.

---

## 5) Output transform construction after GA

BuildOutputsFromSelection(): class/genetic_algorithm.h:1111

After best chromosome is chosen:
1. Build graph_ adjacency for selected edges.
2. BFS from first valid shard to assign initial transforms T_to_root.
3. Run translation-only relaxation for kPoseRelaxIterations.
4. Save transforms_ per shard (unreached nodes remain identity).

Important detail:
- relaxation updates only translation, not rotation.

---

## 6) Diagnostics currently in this file

- enable_debug_logging, enable_swap_diagnostics, enable_pose_debug_logging are currently true in this version:
  - class/genetic_algorithm.h:35
  - class/genetic_algorithm.h:36
  - class/genetic_algorithm.h:37

- swap diagnostics run extra full-fitness evaluations for candidate swap tests:
  - class/genetic_algorithm.h:286

These are useful for analysis but expensive.

---

## 7) Suspected issues / possible bugs / logic concerns

This section is intentionally critical and prioritizes what may hurt reconstruction quality.

### 7.1 Debug flags default to true despite "off by default" comment

- Comment says diagnostics should be off by default.
- Actual values are true.

Why it matters:
- can significantly slow GA,
- can change practical search quality if runtime budget is fixed.

Anchors:
- class/genetic_algorithm.h:34-38

### 7.2 Cycle transform orientation likely inconsistent

In cycle penalty:
- key is stored as sorted pair (min,max),
- but transform inserted is raw T_xy without orientation normalization.

If an LCS edge is originally y->x while key is (x,y), cycle composition may use wrong direction.

Why it matters:
- cycle penalty can become numerically wrong,
- GA may reject good chromosomes or prefer bad ones.

Anchor:
- class/genetic_algorithm.h:577-584

### 7.3 Pair-choice penalty depends on selected_density, which depends on use_inlier_score

selected_density is initialized to 0 and only updated if use_inlier_score is true.

If someone disables use_inlier_score but keeps use_pair_choice_penalty on,
- selected_density stays 0,
- gap ratio becomes artificially high,
- penalty behavior becomes unintuitive.

Anchor:
- class/genetic_algorithm.h:474-494

### 7.4 Neighbor penalty threshold may be ineffective

kMaxNeighbors is computed from incident candidate counts and then clipped by num_shards-1.
In practice this can often be so high that neighbor_count almost never exceeds it.

Why it matters:
- term may be effectively dead,
- gives false confidence that over-connection is controlled.

Anchors:
- class/genetic_algorithm.h:412
- class/genetic_algorithm.h:1320-1340
- class/genetic_algorithm.h:507

### 7.5 Overlap penalty checks only one connected component

Overlap BFS starts from the first valid shard and does not iterate all components.
Disconnected components are not overlap-checked against each other.

Why it matters:
- physically implausible multi-component solutions may avoid overlap penalty.

Anchor:
- class/genetic_algorithm.h:797-843

### 7.6 Output pose construction uses first-arrival BFS and translation-only relaxation

BuildOutputsFromSelection:
- rotation follows first BFS arrival path,
- only translation is relaxed afterward.

Why it matters:
- rotational inconsistencies can remain,
- final visualization/result can look suboptimal even if selected edges are good.

Anchor:
- class/genetic_algorithm.h:1167-1229

### 7.7 Fixed RNG seed every Run

srand(42) is called inside Run().
If Run() is invoked repeatedly in outer loops, randomness repeats deterministically each run.

Why it matters:
- lower stochastic diversity across repeated GA runs,
- easier to get stuck in repeated search patterns.

Anchor:
- class/genetic_algorithm.h:122

### 7.8 Active pair range penalty is disabled by default

kMinActivePairRatio and kMaxActivePairRatio are present but inactive unless
use_active_pair_range_penalty is turned on.

Why it matters:
- graph sparsity/density regularization is currently not influencing GA.

Anchors:
- class/genetic_algorithm.h:32
- class/genetic_algorithm.h:524

---

## 8) Recommended future work (practical roadmap)

### Priority A (likely highest impact)

1. Fix cycle orientation consistency.
- Normalize transform direction in active_transforms exactly like pair_edges does.
- Re-validate cycle penalty values before/after.

2. Check all components in overlap penalty.
- Build poses for every connected component, not only from first valid node.
- Penalize overlaps intra-component and optionally inter-component.

3. Improve output pose construction.
- Replace first-arrival BFS pose with weighted spanning tree initialization (max confidence),
- then run joint translation+rotation relaxation.

### Priority B (stability and calibration)

4. Decouple pair-choice penalty from use_inlier_score.
- Always compute selected density for pair-choice term regardless of inlier-score reward toggle,
- or explicitly gate pair-choice term by use_inlier_score.

5. Revisit neighbor penalty threshold.
- Use unique selected neighbor constraints or a fixed domain-informed cap,
- not candidate-incidence-based cap.

6. Turn diagnostics off by default in production runs.
- keep switches for debug sessions.

### Priority C (search quality improvements)

7. Adaptive mutation/diversification.
- Increase mutation when population stagnates,
- reduce when progress is steady.

8. Memetic step (critical feature proposal).
- After crossover/mutation, run a small local optimizer on 1-3 uncertain genes
  using incremental delta-fitness approximations.
- This can significantly improve final quality over pure GA.

9. Better uncertainty-aware edge confidence.
- Build edge weight from inlier, score, geometric residual, and maybe local ambiguity,
- use this confidence consistently in selection, cycle checks, and output graph building.

---

## 9) Quick mental model (simple)

Think of this GA as selecting one glue option per shard pair, then asking:
- Are the selected glues locally strong?
- Do all glues agree globally in loops?
- Does the assembled shape stay connected and physically plausible?

Current code already tries to do all three, but several consistency details (especially cycle orientation and pose construction) can still push it toward suboptimal final assemblies.
