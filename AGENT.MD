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
* **Spanning Tree Constraint ($n-1$ Active Edges):** For a valid reassembly of $n$ active shards, the physical structure must form a connected tree containing exactly $n-1$ active connections. Any candidate chromosome is passed through a **BFS Spanning Tree Repair Operator** (`RepairChromosome()`) that deactivates (sets to $0$) any redundant edges that form cycles, ensuring the final phenotype represents a valid spanning tree of at most $n-1$ edges.

### 🚀 Initialization & Pair Choice Selection
* **Sampling Initial Choices (`SampleGroupChoice`):** For each gene in the initial population of size $N = 350$, we sample a choice:
  * **Sparsity Rate (`kInitialPairInactiveRate = 0.05`):** A $5\%$ chance that a gene starts deactivated ($0$).
  * **Density Ranking:** Transform candidates for active pairs are ranked by a quality density score:
    $$\text{Density} = \frac{\text{inliers}}{1.0 + \text{score}}$$
  * **Linear Rank Selection:** Candidates are sorted by density, and we use a linear-rank probability distribution to select the candidate, biasing the population toward strong initial alignments.

### 🏋️ Fitness Function Components
The fitness function evaluates global geometric and topological consistency. Since the **BFS Spanning Tree Repair Operator** prunes all non-tree edges to prevent loops, cycles are eliminated and residual propagation errors on tree edges are mathematically zero. Therefore, cycle penalties ($F_{\text{cycle}}$) and edge/rotational residuals ($F_{\text{res}}$) are disabled in practice (`use_cycle_penalty = false`, `use_edge_residual = false`, `use_rot_residual = false`). 

The active fitness function simplifies to:
$$F(\mathbf{c}) = \sum q_k + F_{\text{conn}} + F_{\text{consensus}} - F_{\text{overlap}}$$

1. **Match Quality Score ($q_k$):** Proportional to $\frac{\text{inliers}}{1 + \text{score}}$, rewarding high-inlier counts combined with low dissimilarity scores (precision).
2. **Connectivity Reward ($F_{\text{conn}}$):** Rewards the size of the largest connected component and penalizes fragmented component counts.
3. **Consensus Support Reward ($F_{\text{consensus}}$):** Rewards structural consensus. For each placed component, we check every alternative (inactive) pair of shards. If there is a candidate LCS match that aligns them within a translation threshold ($12.0$ mm) and rotation threshold ($0.22$ rad) of their current relative pose, we reward the fitness proportionally to the match's inlier count ($0.25 \times \text{inliers}$).
4. **Overlap Penalty ($F_{\text{overlap}}$):** Evaluates physical interpenetration of sherd point-clouds using voxel-downsampled collision clouds and kd-tree distance queries.

### 🔀 Genetic Operators
* **Shard-Partition Crossover:** Randomly assigns shards to one of two parents. 
  * If both shards in a gene belong to the same parent, the child inherits that parent's gene.
  * If they are split, the child inherits from either parent with equal probability.
  This preserves locally consistent spatial clusters.
* **Shard-Centric Mutation (Probability $p_m = 0.20$):** Focuses on a single random shard, selects one of its incident genes, and either deactivates it (50% chance) or samples a new candidate transform (50% chance). (Tuned from 0.10 to enhance exploration of symmetric choices).

---

## 3. Results and Performance Baseline

Our GA achieves perfect results on moderate or simple pots (A, B, C, F) but struggles on larger, highly fragmented pots (G, H, J) when compared to the original SfS++ BAISER beam search:

### Quantitative Results (Accuracy Comparison)
| Pottery | Sherds | SfS++ SA (%) | SfS++ EA (%) | GA SA (%) | GA EA (%) | Status |
| :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| **Pot A** | 8 | 100.0 | 100.0 | **100.0** | **100.0** | Perfect match |
| **Pot B** | 9 | 100.0 | 100.0 | **100.0** | **100.0** | Perfect match |
| **Pot C** | 4 | 100.0 | 100.0 | **100.0** | **100.0** | Perfect match |
| **Pot F** | 6 | 100.0 | 100.0 | **100.0** | **100.0** | Perfect match |
| **Pot G** | 7 | 100.0 | 100.0 | **85.7** | **35.0** | Significant lag |
| **Pot H** | 11 | 90.9 | 88.2 | **36.6** | **11.8** | Low accuracy |
| **Pot J** | 11 | 72.7 | 76.2 | **54.5** | **19.0** | Low accuracy |

*Note: SA = Sherd Accuracy (correct placement), EA = Edge Accuracy (correct pairwise connections).*

---

## 4. Identified Weaknesses & GA Bottlenecks

### 🔍 Unidentified Ground Truth Matches at Start
* **The Problem:** The pairwise curve matching algorithm (`FeatureComp`) does not identify all valid ground truth matches at the start of the process. This is due to complexities in estimating the initial symmetric axis and the resulting alignment of the curves.
* **Impact:** Because the search space is limited to matches computed at the start, the GA cannot find the correct overall assembly if critical edges are missing. This is the main issue that plagues reconstruction on highly fragmented pots like **Pot G**.
* **Implication:** The pottery needs to be assembled incrementally so that as some sherds are correctly aligned and positioned, the remaining matches can be resolved and uncovered.

---

## 5. Architectural Map of Main Files

* **`class/genetic_algorithm.h`:** **The main target of our work.** Defines `GeneticAssembler`, `Chromosome` structs, crossover, mutation, and the full multi-objective fitness evaluation.
* **`main.cpp`:** Entry point. Orchestrates feature matching, iterative GA loops, post-ICP modes, and keyboard handlers.
* **`class/data_structure.h` / `.cpp`:** Contains `Geom`, `Trans`, and `LCSIndex` definitions. Contains `CountResult` which evaluates the predicted active edges against Ground Truth.
* **`class/feature_matching.h` / `.cpp`:** LCS pairwise feature curve alignment and candidate pruning routines.
* **`class/reconstruction.h` / `.cpp`:** Automated joint-ICP alignment and Ceres solver settings.

---

## 6. Guidelines for Your Work

> [!CAUTION]
> **CRITICAL WARNING FOR AI AGENTS: PRESERVE VISUALIZATION & TRANSFORMATION CODE**
> AI models have a high failure rate when editing coordinate transformations and visualization logic in this codebase.
> * **Transformation Integrity:** Ensure that the absolute pose transformation matrices (`T_best`, `T_ga_eval`, `T_ga_vis`) mapping local coordinates to assembled spaces are kept correct and clean. Any corruption here cascades into breaking PCL viewer positioning, downstream Axis Refinement, and the final edge validation step (`CountResult`).
> * **PCL & Ceres Constraints:** Keep your edits clean, strictly preserve the interactive PCL visualization toggle features, and respect the Ceres solver safety parameter-block constraints!

When designing modifications to improve the GA, prioritize:
1. **Targeted Gene Pruning & Local Fitness:** Developing mechanisms to locate and mutate/prune bad individual edges rather than penalizing the entire chromosome blindly.
2. **Spatial-Cluster Preserving Crossover:** Advanced crossover techniques that preserve entire 3D connected sub-assemblies instead of randomly splitting genes.
3. **Global Vessel Profile Prior:** Integrating axis-distance profiles or global symmetry constraints into the fitness score to eliminate geometrically valid false positives.
4. **Hierarchical / Cooperative Co-evolution GA:** Grouping fragments into stable sub-assemblies (sub-populations) and evolving them collectively to scale to larger pots.
