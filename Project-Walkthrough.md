# Comprehensive Project Structure & History

This document outlines the current state, architecture, and evolutionary history of the Automated 3D Pottery Reconstruction project.

## 1. Project Objective and Pipeline
**Goal:** Automatically reconstruct fractured 3D pottery artifacts from unordered geometric scans (sherds).

### The Reconstruction Pipeline
1. **Data Loading & Preprocessing:** Loads origin point clouds (`.pcd`) and 3D meshes (`.obj`), calculates normals, and identifies symmetric axes for initial alignment (`z-axis`).
2. **Feature Matching ([FeatureComp](file:///home/zain/Zain/Semester%206/CI/Project/SfS/structure-from-sherds-pp/class/feature_matching.cpp#1375-1495))**: Identifies potential matching fracture boundaries between pairs of sherds, outputting candidate match transformations and inlier scores.
3. **Global Assembly Optimization (Genetic Algorithm)**: Traverses the huge search space of potential pairwise match combinations to find the globally consistent structure of the complete pot.
4. **Local Refinement (ICP):** Applies an Iterative Closest Point fine registration on the graph produced by the GA to precisely lock the sherds together.
5. **Evaluation & Visualization:** Computes accuracy scores ([CountResult](file:///home/zain/Zain/Semester%206/CI/Project/SfS/structure-from-sherds-pp/class/data_structure.h#432-435)) against Ground Truth (GT) and saves the transformed `.xyz` and `.obj` outputs.

---

## 2. Evolution of the Assembly Algorithm
The core problem in the pipeline is selecting the correct combinatorial subset of pairwise matches to assemble the full pot out of all potential candidates.

### **Phase 1: Beam Search (Original Implementation)**
The original assembly engine relied on a Beam Search algorithm ([RankingSubgraph](file:///home/zain/Zain/Semester%206/CI/Project/SfS/structure-from-sherds-pp/class/ranking_system.h#65-81)).
- **Mechanism:** Incrementally built graphs by selecting the highest-scoring edges heuristically.
- **Limitation:** Highly prone to early local optima. If an incorrect edge scored highly early in the search, the algorithm struggled to "backtrack" and find the true global shape, particularly for highly degraded artifacts.

### **Phase 2: Permutation-Based Evolutionary Algorithm**
We first replaced Beam Search with an Evolutionary Algorithm (EA).
- **Mechanism:** The chromosome was represented as a continuous permutation string (ordering) of matches. The decoder would walk through the permutation and add matches to the assembly sequentially, skipping invalid overlaps.

### **Phase 3: Binary-Choice Genetic Algorithm (Current Architecture)**
We pivoted to a far more intuitive and effective **Binary-Choice Representation**.
- **Mechanism:** Each candidate pairwise match (from [FeatureComp](file:///home/zain/Zain/Semester%206/CI/Project/SfS/structure-from-sherds-pp/class/feature_matching.cpp#1375-1495)) is assigned a boolean gene.
  - `TRUE` = The match is active and included in the assembly.
  - `FALSE` = The match is deactivated.
- **Guided Repair (`RepairBinaryChromosome`):** Because multiple conflicting matches can exist for the exact same pair of sherds, we implemented a deterministic repair step. After crossover and mutation, the repair mechanism scans the bit-vector. If multiple matches are active for the exact same shard pair, it forcibly retains only the one with the highest inlier count, disabling the rest.

---

## 3. The Architecture of [GeneticAssembler](file:///home/zain/Zain/Semester%206/CI/Project/SfS/structure-from-sherds-pp/class/genetic_algorithm.h#21-1021) ([genetic_algorithm.h](file:///home/zain/Zain/Semester%206/CI/Project/SfS/structure-from-sherds-pp/class/genetic_algorithm.h))
The GA replaces the old Beam Search and is defined entirely within [genetic_algorithm.h](file:///home/zain/Zain/Semester%206/CI/Project/SfS/structure-from-sherds-pp/class/genetic_algorithm.h). We integrated it cleanly into [main.cpp](file:///home/zain/Zain/Semester%206/CI/Project/SfS/structure-from-sherds-pp/main.cpp) so it can intercept the list of pruned matches (`LCS_out`) before ICP.

### **The Fitness Function**
The heart of our Genetic Algorithm is the fitness evaluator ([EvaluateFitness](file:///home/zain/Zain/Semester%206/CI/Project/SfS/structure-from-sherds-pp/class/genetic_algorithm.h#228-540)). Because simple pairwise scores aren't enough to guarantee a physically plausible 3D object, the fitness function incorporates a wide array of geometric constraints:

| Component | Goal | Effect on Fitness |
| :--- | :--- | :--- |
| **Inlier Score Reward** | Reward matches with high structural consistency. | `+ (inliners / (1 + score))` |
| **Connectivity Reward** | Encourage the formation of large contiguous groups (a fully assembled pot). | `+ (Largest Component Size * Const)` |
| **Component Penalty** | Penalize disjointed islands of sherds (fragmented assemblies). | `- (Number of distinct islands * Const)` |
| **Cycle Consistency** | Penalize physically impossible spatial loops (e.g. A->B->C->A creates a spatial gap). | `- (Sum of Loop Discrepancy Norms)` |
| **Edge & Rot Residuals** | Penalize relative positional and rotational errors within the connected graph. | `- (Residual * PenaltyWeight)` |
| **Neighbor Penalty** | Prevent a single sherd from connecting to an unrealistic number of neighbors. | `- (Neighbor Count - Limit)` |

---

## 4. Current Phase: Ablation Studies & Experiments
Because our fitness function is highly complex and multi-objective, our current phase involves diagnosing exactly which components are driving the reconstruction accuracy and which might be creating noise.

### **The Ablation Infrastructure in [main.cpp](file:///home/zain/Zain/Semester%206/CI/Project/SfS/structure-from-sherds-pp/main.cpp)**
We have constructed a massive automated experiment block inside [main.cpp](file:///home/zain/Zain/Semester%206/CI/Project/SfS/structure-from-sherds-pp/main.cpp) to independently toggle fitness components on and off.
- **The Flags:** Six boolean flags inside [GeneticAssembler](file:///home/zain/Zain/Semester%206/CI/Project/SfS/structure-from-sherds-pp/class/genetic_algorithm.h#21-1021): `use_inlier_score`, `use_connectivity_reward`, `use_cycle_penalty`, `use_edge_residual`, `use_rot_residual`, and `use_neighbor_penalty`.
- **The Loop:** We test 9 distinct configurations (e.g., `no_cycle_penalty`, `all_components`, `no_penalties`).
- **ICP Toggle:** Every single configuration is tested twice—once with the final ICP Refinement pass enabled (`ICP_True`) and once relying purely on the raw GA transformation data (`ICP_False`). This helps isolate whether the GA is fundamentally broken or if ICP is just fixing a close-enough EA guess.
- **Statistical Rigidity:** Each test is run across `5` independent random seeds.
- **Outputs:** 
  - The results are grouped dynamically into `Result/ablation/<ICP_State>/<Configuration_Name>/trial_<N>/` folders.
  - Generates comprehensive [.txt](file:///home/zain/Zain/Semester%206/CI/Project/SfS/structure-from-sherds-pp/CMakeLists.txt) summary files outlining Average Sherd Accuracy, Average Edge Accuracy, and Average Fitness for each configuration.
  - Saves the raw spatial `.obj` files for visual verification in MeshLab/Blender.

---

## 5. Summary of Recent Bug Fixes
As part of integrating the Ablation Loop and GA visualization, we solved the following critical issues:
1. **Visualization Distortion (Double-Transforms):** Fixed a bug where viewing the GA-only assembly resulted in a scrambled mess because axis-alignment rotations were being redundantly applied onto already-transformed sherds.
2. **`std::experimental::filesystem` Linker Error:** Replaced failing C++17 filesystem calls with legacy UNIX `system("mkdir -p ...")` commands to ensure compilation stability within the target Linux container.
3. **PCL Typenames:** Corrected a `curvatures`/`normals` missing field error when writing out raw meshes during ablation trials by correctly typing the point clouds to `pcl::PointXYZ` instead of `pcl::PointNormal`.
