# Agent Handover: Pottery Reconstruction (Genetic Algorithm)

This document serves as a "mental state" transfer for any agent taking over this project. It outlines the current technical status, the methodology that led to our breakthroughs, and the specific "thinking style" required to pair-program with the USER on this codebase.

## 1. Project Context
The goal is the automated 3D reconstruction of pottery from irregular shards. The core engine is a **Genetic Algorithm (GA)** that selects the best pre-computed "Choices" (geometric matches) to form a coherent pot.

### The Core Loop:
1. **FeatureComp**: Pre-computes potential snapped matches (Candidates) between shard pairs.
2. **Genetic Algorithm**: Evolves a population of "Chromosomes" (potential assemblies) to maximize global fitness.
3. **Fitness Function**: Evaluates assemblies based on connectivity, feature alignment, and collision penalties.

---

## 2. Current State of the Engine
We have recently transitioned from a "Greedy/Exploitative" search to a "Diversity/Diagnostic" search.

### Key Implementation Details:
- **Selection Logic (`SampleGroupChoice`)**: 
    - **0.0 Greedy Rate**: We explicitly stopped "snapping" to Choice 1.
    - **Density-Based Weighting**: Weights are calculated as `inliner / (1.0 + score)`. This ensures that geometric quality (`score`) can outweigh a high point count (`inliner`) if the latter is a false positive.
    - **Roulette Wheel Sampling**: Uses `std::mt19937` for thread-safe, high-quality randomness to pick initial genes.
- **Fitness Function**:
    - **MAX-based Residuals**: We use the maximum residual error across active edges rather than an average to prevent "good" edges from masking a single "catastrophic" misalignment.
- **Collision Logic**: 
    - **Centroid-Based Repair**: Currently uses a linear depth-proportional penalty. 
    - **Roadblock**: Aggressive point-cloud collision was rejected due to performance bottlenecks and the risk of penalizing valid boundary contacts.

---

## 3. The "Agent Mindset" (How to perform like me)

### A. Data-Driven Skepticism
Never assume the "Best Choice" (Rank 1) is correct. Our biggest breakthrough came when we implemented a **Ground Truth Audit** in `genetic_algorithm.h`. 
- **Action**: When the GA fails, don't just tweak parameters. Implement a diagnostic that cross-references the GA's current decisions with the Ground Truth (`GT_trans`).
- **Insight**: We discovered that selected choice often has more inliers due to plausible local fits that hurt globally.

### B. Analytical Debugging
When the USER reports a failure (e.g., "6/8 accuracy"), your first step must be to analyze the `ga_diagnostics.log`.
- Look for **rad_err** and **trans_err** in the audit.
- If `rad_err` is ~3.14, the shard is upside down.
- If `rad_err` is low but the shard is deactivated, it likely triggered a **Collision Repair**.

### C. Explaining Complex Concepts
The USER values clear explainations.

---

## 4. Interaction History & Methodology
Our successful interaction pattern has been:
1. **Identify the Roadblock**: (e.g., Shard 4 won't snap correctly).
2. **Implement the "Truth-Check"**: Modify the code to print exactly where the "Truth" is hiding in the candidate list.
3. **Analyze the Weights**: Check why the GA isn't picking the Truth.
4. **Tune the Selection**: Shift the bias from raw numbers (inliers) to geometric density (`inliner / score`).

---

## 5. Active Roadblocks & Next Steps
- **Boundary interpenetration**: Sherds are currently "meeting" but sometimes overlapping slightly in 3D space. 
- **"Stiff" Assemblies**: Because collision penalties are coarse, the assembly can get stuck in a local minimum where one overlap prevents other shards from reaching their edges.
- **Constraint Gating**: We need a way to use high-accuracy collision checks *only* for non-neighboring shards, while being more lenient with verified neighbors.

---

> [!IMPORTANT]
> **To the next agent**: Always keep the `Project-Walkthrough.md` updated. It is the USER's primary map of our technical journey. Do not touch `README.md` unless explicitly asked—it is the project's permanent documentation.
