# 3D Pottery Reconstruction: Project Walkthrough & Architecture

## 1. Project Objective
The goal of this project is to automatically reassemble fractured 3D pottery artifacts from unordered geometric scans (sherds). This is a complex combinatorial optimization problem where the system must identify which sherds belong together and how they should be oriented in 3D space.

## 2. The Core Challenge: The "Symmetry Trap"
The primary difficulty in pottery reconstruction arises from the geometric symmetry of the artifacts. 
*   **False Positives:** Many fracture boundaries can "fit" together in multiple ways as they can appear on any side of the fracture boundary.
*   **Inlier Bias:** Standard surface-fit metrics (inliers) may reward false positives even more than the correct match because they hit the surface perfectly.
*   **The "Quality" Score:** While Inliers measure the *quantity* of contact, the **Score** represents the *quality* (RMS error). A high inlier count with a high score indicates a "sloppy" or forced fit, whereas a low score indicates a precise geometric lock. Though, still this can produce false positives so it doesn't gurantee that the two sherds actually go together.
*   **Search Space:** With dozens of sherds and hundreds of candidate matches, the number of possible assemblies is astronomically large.

## 3. Legacy Approach: Beam Search
The original implementation (based on existing literature) utilized a **Beam Search** algorithm.
*   **Mechanism:** It built the assembly incrementally by greedily selecting the highest-scoring edges while maintaining a "beam" of the top-$N$ most promising partial graphs.

## 4. Modern Approach: Diagnostic-Driven Genetic Algorithm
We have replaced Beam Search with a **Genetic Algorithm (GA)**.

### **Genetic Representation**
*   **Chromosome:** An array of integers where each index represents a potential "Edge" (connection between two shards) and the value represents the "Choice" (which specific transformation from the candidate list to use).
*   **Choice 0:** Deactivates the connection entirely.

### **Selection Policy (Exploration over Exploitation)**
To prevent the GA from being "tricked" by high-inlier flips, we use the following selection model:
*   **Density-Based Weighting:** Selection is proportional to `inlier / (1 + score)`, ensuring that geometric quality (low error score) can out-compete raw inlier counts.
*   **Increased Population (200):** A larger population ensures that "True" matches (which may rank 2nd or 3rd in inliers) are well-represented in the gene pool.

### **Iterative Refinement Loop**
Instead of a single pass, the system runs in cycles:
1.  **GA Search:** Finds the best global assembly.
2.  **Repositioning:** Shards are moved to their new positions.
3.  **Feature Re-computation:** `FeatureComp` is re-run on the assembled state to find "tight-fit" matches that were invisible in the initial scatter.

---

## 5. Diagnostic Tools & Validation
We integrated advanced diagnostics to "audit" the GA's performance against known reality:

*   **Ground Truth (GT) Candidate Audit:** This tool compares every pre-computed candidate match against the Ground Truth transformation. It accounts for coordinate frame discrepancies (`[Wrapped]` errors) to identify exactly which choices are "True" and which are "Flips."
*   **Collision Detection:** The system detects interpenetrating sherd volumes using a distance-based metric between centroids and radii.

---

## 6. The Fitness Function Components
The GA evaluates solutions using a multi-objective fitness function:

| Component | Purpose |
| :--- | :--- |
| **Inlier Reward** | Rewards matches with high surface consistency. |
| **Connectivity Reward** | Encourages the formation of a single, large contiguous pot. |
| **Cycle Penalty** | Penalizes spatial loops that don't close perfectly (indicating drift). |
| **Residual Penalties** | Penalizes relative positional and rotational errors in the graph (Max-based signal). |
| **Overlap Penalty** | Penalizes physical occupancy to handle symmetry (Currently Centroid-Based). |

---

## 7. Current Problems & Future Work
two core roadblocks remain:

### **1. Signal Dilution in Residuals**
Initially, spatial inconsistencies were averaged across the entire assembly, which diluted the "signal" of single bad edges. We have transitioned to a **MAX-based residual penalty** which ensures that the worst-case geometric error is always penalized, preventing bad edges from hiding in a "locally good" cluster.

### **2. Boundary Precision vs. Physical Collision**
Standard spherical approximations for shards (centroids + radii) are sometimes too coarse to distinguish between a "tight fit" and a "collision."
*   **The Conflict:** If we make collision detection too strict (e.g., using point-cloud proximity), we risk penalizing valid neighbors that naturally touch at the fracture boundary.
*   **The Result:** Currently, sherds that are neighbors often "overlap" slightly in 3D space rather than meeting perfectly at the boundary. This localized overlap prevents them from correctly aligning with *other* neighbors, leading to incomplete assemblies. 
*   **Next Steps:** We need a collision model that can distinguish between "Interpenetration" (bad) and "Contact" (good).

### **2. Basic Crossover Logic**
The current crossover mechanism is a standard edge recombination (inherited from TSP-style problems).
*   **The Problem:** It does not respect **3D spatial clusters**. A crossover might mix the "top half" of one good assembly with the "bottom half" of another, but if they don't share a common coordinate frame, the result is a disjointed mess that the GA cannot easily fix.
*   **The Goal:** We need a "Spatial-Cluster Crossover" that preserves geometrically consistent groups of sherds during reproduction.
