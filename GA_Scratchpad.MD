# 🧬 Genetic Algorithm Loop & Mechanics Cheat Sheet

This scratchpad serves as a practical reference for understanding the core mechanics of our Genetic Algorithm (GA) pipeline for 3D pottery reconstruction. It details how chromosomes are structured, how reproduction and mutation operate, how we maintain diversity to avoid traps, and how you can tune these systems to build a more robust assembler.

---

## 📐 Chromosome Representation

Our chromosome represents a candidate 3D assembly graph in a discrete format:

$$\mathbf{c} = [x_1, x_2, \dots, x_m]$$

*   **Chromosome Length ($m$):** The number of unique undirected pairs of shards that remain after pairwise curve matching and pruning.
*   **Gene Value ($x_i$):**
    *   $x_i = 0$: The connection between this pair of shards is **deactivated** (rejected).
    *   $x_i = j \in [1, S]$: The connection is **active**, using the $j$-th precomputed LCS rigid 3D transformation candidate (where $S$ is the total candidate pool size for this pair).

> [!TIP]
> **Why this representation works:** By representing the choice between two shards as a single gene, the chromosome naturally guarantees that *at most one* transform can be active between any two shards at any time. This mirrors physical reality where two fragments cannot be put together in multiple different alignments simultaneously.

---

## 🚀 Initialization

To start the evolution, we initialize a population of **`kPopulationSize = 100`** chromosomes using a smart, guided approach rather than pure random noise.

For each gene, we sample an initial state using **`SampleGroupChoice`**:
1.  **Stochastic Sparsity Check:** There is a baseline probability (`kInitialPairInactiveRate = 0.05` or **5%**) that a gene starts as `0` (deactivated). This ensures the initial assembly graphs are relatively sparse and don't start with excessive point collisions.
2.  **Density-Based Rank Selection:** If the edge is active, we rank all its precomputed transform candidates using a **density score**:
    $$\text{Density} = \frac{\text{inliers}}{1.0 + \text{score}}$$
    *   **Inlier Count:** The *quantity* of matching points (how much boundary is in contact). More is better!
    *   **RMS Score:** The *quality* of the fit (the residual misalignment error). Lower is better!
3.  **Linear Rank Selection:** The candidates are sorted by density. We assign a probability weight to each candidate linearly based on its sorted rank: the best candidate gets the highest probability weight, and the worst gets the lowest. We then spin a weighted wheel to choose the initial transform.

---

## 🏋️ Fitness Calculation

The fitness function evaluates how globally consistent, spatially closed, and physically valid a candidate 3D assembly (chromosome) is. It balances rewards for good contacts against heavy penalties for geometric drift or physical interpenetration:

$$F(\mathbf{c}) = \sum q_k + F_{\text{conn}} - F_{\text{cycle}} - F_{\text{res}}^{\text{t}} - F_{\text{res}}^{\text{r}} - F_{\text{overlap}}$$

Here is how each component is calculated and why it matters:

### 1. Match Quality Reward ($\sum q_k$)
For every active connection (gene) in the chromosome, we sum their density scores:
$$q_k = \frac{\text{inliers}}{1.0 + \text{score}}$$
*   **Why?** This rewards assemblies that have a large number of contact points (high inliers) with very precise geometric locks (low dissimilarity score).

### 2. Connectivity Reward & Component Penalty ($F_{\text{conn}}$)
Encourages the shards to group into a single, cohesive pot:
$$\text{Reward} = k_{\text{connectivity\_reward}} \times \text{largest\_component\_size} - k_{\text{connectivity\_component\_penalty}} \times (\text{num\_components} - 1)$$
*   **Key Parameters:** `kConnectivityReward = 10.0`, `kConnectivityComponentPenalty = 0.0`.
*   **Why?** It rewards chromosomes that build the largest possible connected subgraph, preventing the GA from producing highly fit but fragmented sub-assemblies.

### 3. Cycle Consistency Penalty ($F_{\text{cycle}}$)
Evaluates spatial loops of three shards (triplets) where all three connections are active:
*   **Composed Relative Transform:** $T_{\text{composed}} = T_{bc} \times T_{ab}$ (transforming from A to B, then B to C)
*   **Direct Relative Transform:** $T_{ac}$ (transforming directly from A to C)
*   **Error:** $\text{error} = \|t_{\text{composed}} - t_{ac}\|$
*   If the composed translation differs from the direct translation by more than `cycle_threshold = 30.0` mm, the excess error is accumulated.
*   **Why?** If three shards form a physical loop, going from $A \rightarrow B \rightarrow C \rightarrow A$ should return to the exact same position. If there is spatial drift, the loop will fail to close. The penalty is doubled and subtracted from the fitness.

### 4. Edge & Rotational Residuals ($F_{\text{res}}$)
Evaluates global relative transform inconsistencies:
*   **Propagation:** The algorithm uses a Breadth-First Search (BFS) starting from a root shard to propagate transforms globally along active edges.
*   **Residual Measurement:** For each active edge, it measures the translational discrepancy $\|t_{\text{pred}} - t_{ab}\|$ and rotational discrepancy $\|\text{log}(R_{\text{pred}} \cdot R_{ab}^T)\|$ against the precomputed pairwise matches.
*   **Why?** If the active edges are geometrically incompatible, propagating poses will yield massive global gaps. The worst-case (MAX-based) errors are penalized to prevent bad edges from hiding in locally good clusters.

### 5. Overlap Penalty ($F_{\text{overlap}}$)
Evaluates physical interpenetration of sherd point-clouds using KD-Tree nearest neighbor distance queries:
*   Points near the fracture boundaries ($< 2.0$ mm) are excluded to allow natural contact, and the remaining surface points are voxel-downsampled.
*   If two shards in the same component overlap, we measure their interpenetration:
    $$\text{Penalty} = k_{\text{overlap\_penalty}} \times \text{mixed\_overlap} \times \text{depth\_term}$$
*   **Key Parameters:** `kOverlapPenalty = 1200.0`.
*   **Why?** This is the ultimate weapon against the **Symmetry Trap** (e.g. mirrored pieces that fit boundaries perfectly but collide physically inside).

---

## 🗳️ Parent Selection

When breeding the next generation, we select parents using **Tournament Selection**:
*   **Tournament Size = 3:** For each parent slot, the algorithm randomly draws 3 individuals from the population and pits them against each other.
*   **The Winner:** The individual with the highest fitness score wins the tournament and is selected as a parent.
*   **Why 3?** It balances *selection pressure*. A larger tournament size (e.g., 5) would make selection too greedy (highly fit parents dominate too quickly, causing premature convergence), while a size of 2 would allow too many weak parents to pass on bad genes.

---

## 🔀 Shard-Partition Crossover

Standard crossover operators randomly swap individual edges, which breaks the physical grouping of parts. To fix this, we use **Shard-Partition Crossover** to preserve local 3D spatial clusters:

1.  **Partitioning:** We randomly assign every physical shard to either `Parent 1` (owner `0`) or `Parent 2` (owner `1`) using a 50/50 coin flip.
2.  **Inheritance Rules:**
    *   **Internal Edges:** If both shards in a gene belong to the same parent partition, the child inherits that parent's gene choice. *(This keeps cohesive 3D clusters intact!)*
    *   **Bridge Edges:** If the shards are split between different parents, the child inherits from either parent with a 50/50 probability.

---

## 👾 Shard-Centric Mutation

Instead of mutating random genes globally, **Shard-Centric Mutation** targets a specific fragment to refine its local connections:

*   **Trigger Rate:** Occurs on an individual chromosome with a probability of `kMutationRate = 0.10` (**10%**).
*   **Localization:** It randomly selects one physical shard, and then picks one of the candidate matching pairs (edges) incident to that shard.
*   **Action Flip (50/50 split):**
    *   **50% chance (Deactivate):** The connection is turned `0` (off) to break potential structural conflicts or overlaps.
    *   **50% chance (Explore):** It runs `SampleGroupChoice` to pick a different candidate transform. If the same choice is drawn, it forces a toggle to guarantee change.

---

## 🛡️ Survivor Selection & Diversity Enforcer

To form the next generation, we combine elitism with active clone protection:

1.  **Elitism (`kElitismCount = 6`):** The absolute top 6 fittest chromosomes are copied directly to the new population without any change.
2.  **Generational Offspring:** The remaining 94 slots are filled by the offspring created via crossover and mutation.
3.  **Premature Convergence Sweep (`EnforceDiversity`):**
    *   We compare all non-elite chromosomes. If two individuals share too many similar gene values (Hamming distance is less than $1/3$ of the chromosome length), they are deemed clones.
    *   **Redundant Clone Eradication:** We replace the clone with a **`GenerateBiasedReplacement`** chromosome.
    *   **Biased Replacement:** It copies the top **40%** (`kBiasInheritRatio = 0.40`) of the clone's active genes that have the highest inlier counts, and re-samples the remaining 60% from scratch. This preserves the strong structural contacts while exploring new alignments for the rest of the pot!

---

## 🏁 Termination & Refining Cycles

The pipeline terminates using two nested layers:

### Layer 1: The Inner Evolutionary Loop (`class/genetic_algorithm.h`)
*   Runs for exactly **`kMaxGenerations = 50`** generations per seed.
*   Runs over **`kNumSeeds = 1`** seed to yield the best overall global 3D assembly.

### Layer 2: The Outer Iterative GA Refinement Loop (`main.cpp`)
Because initial matching is done on scattered sherds, the pipeline runs iteratively:
1.  Runs the inner GA to find a global assembly.
2.  Reposition the sherds in 3D space.
3.  Recomputes matches on this partially-assembled state to uncover "tight-fit" matches.
4.  Runs the GA again.
*   **Patience Termination:** The loop runs up to `kMaxGAIterations = 1` time. If set higher, it stops early if the fitness improvement between successive runs is less than `kConvergenceThreshold = 5.0` for `kMaxPatience = 5` successive times.


1) Make first shard random so we get different bfs trees.
