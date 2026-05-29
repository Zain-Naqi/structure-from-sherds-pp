# Onboarding Agent Guide: 3D Pottery Reconstruction using Genetic Algorithms

Welcome, Agent! This document is designed to immediately onboard you to this repository, explain our objective, walk you through the core Genetic Algorithm (GA) architecture, and highlight the key weaknesses and opportunities for improvement in the GA.

Our main goal is to **improve the Genetic Algorithm** to achieve robust reassembly of complex 3D pottery from scattered sherds, outperforming the legacy greedy beam-search approach on all datasets.

---

## 1. Project Goal & Domain Overview

The goal is to automatically reassemble broken 3D pottery artifacts from unordered geometric scans (sherds). This is a highly challenging combinatorial optimization problem due to:
* **The "Symmetry Trap":** Ceramic vessels are axially symmetric. Fractures are thin and lack distinct features, meaning many sherds can fit together in multiple false-positive orientations (e.g., upside down or mirrored) that yield high geometric overlap but are archaeologically incorrect.
* **Sparse Connectivity:** Valid assemblies represent sparse graphs where only adjacent sherds share an edge.
* **Search Space Complexity:** Even a small pot with 10 pieces has over $3.6 \times 10^6$ combinations, not accounting for 6-DoF continuous transforms.

### The Pipeline Stack (SfS++ Core):
Our codebase builds on top of the **Structure-from-Sherds++ (SfS++)** pipeline, preserving the initial stages but replacing the greedy search:
1. **Feature Extraction:** Estimates the symmetric axis of each sherd, aligns it to the Z-axis, segments the rim/base, and computes 1D curve descriptors.
2. **Pairwise Matching (`FeatureComp`):** Matches curves using Longest Common Subsequence (LCS) to generate a pool of candidate rigid 3D transformations for each pair.
3. **Genetic Algorithm (`class/genetic_algorithm.h`):** Performs global combinatorial search over the candidate pool to identify the true edges and correct transforms.
4. **Post-GA ICP Refinement (`class/reconstruction.cpp`):** Ceres-based non-linear optimization to snap complementary breaklines together.

---

## 2. GA Representation & Mechanics

The GA operates on a discrete chromosome representation of candidate matching pairs.

### 🧬 Genetic Representation
* **Chromosome:** A fixed-length integer vector $\mathbf{c}$ of size $m = |\mathcal{G}|$, where each gene $c_k$ corresponds to a unique undirected pair of shards (**Pair Group**).
* **Gene values:**
  * $c_k = 0$: The connection between this pair of shards is deactivated (rejected).
  * $c_k = j \ge 1$: The $j$-th precomputed LCS transform candidate is selected (activated).

This representation naturally guarantees that at most one transform is selected between any pair of shards, reflecting physical constraints.

### 🏋️ Fitness Function Components
The fitness function evaluates global geometric and topological consistency:
$$F(\mathbf{c}) = \sum q_k + F_{\text{conn}} - F_{\text{cycle}} - F_{\text{res}}^{\text{t}} - F_{\text{res}}^{\text{r}} - F_{\text{overlap}}$$

1. **Match Quality Score ($q_k$):** Proportional to $\frac{\text{inliers}}{1 + \text{score}}$, rewarding high-inlier counts combined with low dissimilarity scores (precision).
2. **Connectivity Reward ($F_{\text{conn}}$):** Rewards the size of the largest connected component and penalizes fragmented component counts.
3. **Cycle Consistency Penalty ($F_{\text{cycle}}$):** Penalizes relative transform drift in spatial loops of three shards (triplets) by composing relative transforms ($T_{bc} \times T_{ab}$) and calculating the distance to the direct transform ($T_{ac}$).
4. **Edge & Rotational Residuals ($F_{\text{res}}$):** After BFS-based pose propagation, it measures global relative transform inconsistencies. Uses a **MAX-based residual** to prevent localized errors from being diluted.
5. **Overlap Penalty ($F_{\text{overlap}}$):** Evaluates physical interpenetration of sherd point-clouds using voxel-downsampled collision clouds and kd-tree distance queries.

### 🔀 Genetic Operators
* **Shard-Partition Crossover:** Randomly assigns shards to one of two parents. 
  * If both shards in a gene belong to the same parent, the child inherits that parent's gene.
  * If they are split, the child inherits from either parent with equal probability.
  This preserves locally consistent spatial clusters.
* **Shard-Centric Mutation (Probability $p_m = 0.10$):** Focuses on a single random shard, selects one of its incident genes, and either deactivates it (50% chance) or samples a new candidate transform (50% chance).

---

## 3. Results and Performance Baseline

Our GA achieves excellent results on moderate or simple pots (A, C, F) but struggles on larger, highly fragmented pots (B, G, H, J) when compared to the original SfS++ BAISER beam search:

### Quantitative Results (Accuracy Comparison)
| Pottery | Sherds | SfS++ SA (%) | SfS++ EA (%) | GA SA (%) | GA EA (%) | Status |
| :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| **Pot A** | 8 | 100.0 | 100.0 | **100.0** | **100.0** | Perfect match |
| **Pot B** | 9 | 100.0 | 100.0 | **100.0** | **60.0** | Lags on edges |
| **Pot C** | 4 | 100.0 | 100.0 | **100.0** | **100.0** | Perfect match |
| **Pot F** | 6 | 100.0 | 100.0 | **100.0** | **83.3** | Missed edges |
| **Pot G** | 7 | 100.0 | 100.0 | **85.7** | **35.0** | Significant lag |
| **Pot H** | 11 | 90.9 | 88.2 | **36.6** | **11.8** | Low accuracy |
| **Pot J** | 11 | 72.7 | 76.2 | **54.5** | **19.0** | Low accuracy |

*Note: SA = Sherd Accuracy (correct placement), EA = Edge Accuracy (correct pairwise connections).*

---

## 4. Identified Weaknesses & GA Bottlenecks

Every agent must understand these core limitations, which are the primary targets of our optimization work:

### 🔍 1. The "Needle in a Haystack" Constraint Problem (Critical)
Reassembling a pot requires a specific, sparse subgraph of edges to be active **at the same time** to create a valid 3D structure.
* **The Problem:** The current GA penalizes an entire chromosome (individual) for having contradictory or floating edges, but **it cannot identify which specific edge (gene) is the bad actor**.
* **Impact:** Random crossover and mutation steps disrupt correct subsets of active genes. We are trying to find a perfect combination of pairs purely by stochastic trial-and-error across a massive search space.

### 📐 2. Absence of a Global Vessel Profile Prior
* **The Problem:** The fitness function only evaluates local pairwise relationships (LCS matches, loop closure, overlaps). It has no concept of the expected overall profile shape or global rotational symmetry of the finished vessel.
* **Impact:** The GA can converge to "false positive" assemblies—solutions with high fitness (strong connectivity, zero cycle penalties, zero overlaps) that are geometrically closed but form archaeologically incorrect shapes.

### 💥 3. Boundary Precision vs. Physical Collision Conflict
* **The Problem:** The current overlap check excludes points within 2.0 mm of the breakline boundaries to prevent penalizing valid contact, and voxel-downsampled the remainder.
* **Impact:** This makes the collision detection too coarse to distinguish between healthy contact at the fracture boundary (good) and slight physical interpenetration (bad). Shards often "hover" or interpenetrate slightly in local space, which halts correct downstream alignment.

### 🔀 4. Lack of Hierarchical and Incremental Building
* **The Problem:** The GA treats reassembly as a single global optimization task, placing all sherds simultaneously.
* **Impact:** This approach scales poorly as the sherd count increases. In contrast, SfS++'s beam-search BAISER builds assemblies incrementally, registering one sherd at a time to filter out incorrect connections early.

---

## 5. Architectural Map of Main Files

* **`class/genetic_algorithm.h`:** **The main target of our work.** Defines `GeneticAssembler`, `Chromosome` structs, crossover, mutation, and the full multi-objective fitness evaluation.
* **`main.cpp`:** Entry point. Orchestrates feature matching, iterative GA loops, post-ICP modes, and keyboard handlers.
* **`class/data_structure.h` / `.cpp`:** Contains `Geom`, `Trans`, and `LCSIndex` definitions. Contains `CountResult` which evaluates the predicted active edges against Ground Truth.
* **`class/feature_matching.h` / `.cpp`:** LCS pairwise feature curve alignment and candidate pruning routines.
* **`class/reconstruction.h` / `.cpp`:** Automated joint-ICP alignment and Ceres solver settings.

---

## 6. Guidelines for Your Work

When designing modifications to improve the GA, prioritize:
1. **Targeted Gene Pruning & Local Fitness:** Developing mechanisms to locate and mutate/prune bad individual edges rather than penalizing the entire chromosome blindly.
2. **Spatial-Cluster Preserving Crossover:** Advanced crossover techniques that preserve entire 3D connected sub-assemblies instead of randomly splitting genes.
3. **Global Vessel Profile Prior:** Integrating axis-distance profiles or global symmetry constraints into the fitness score to eliminate geometrically valid false positives.
4. **Hierarchical / Cooperative Co-evolution GA:** Grouping fragments into stable sub-assemblies (sub-populations) and evolving them collectively to scale to larger pots.

Keep your edits clean, preserve the PCL visualization toggle features, and respect the Ceres safety parameter-block constraints!
