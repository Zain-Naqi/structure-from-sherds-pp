# Implementation Plan: Surface-Aware Collision Detection

This plan addresses the persistent issue of sherd interpenetration in the Genetic Algorithm assembly. By replacing the coarse centroid-based spherical overlap check with a high-precision surface-to-surface proximity check, we aim to resolve symmetry flips and boundary alignment errors.

## Problem Statement
The current **Centroid-Based** collision logic treats every sherd as a sphere. This leads to:
1.  **False Positives**: Valid connections are penalized because the "spheres" overlap even if the surfaces just touch.
2.  **False Negatives**: Long/thin shards can "stab" through each other without their centroids getting close enough to trigger a penalty.
3.  **Boundary Blindness**: The GA cannot distinguish between a tight fracture snap and a shard being buried several millimeters inside its neighbor.

## Proposed Changes

### [Genetic Algorithm]

#### [MODIFY] [genetic_algorithm.h](file:///home/taha/Desktop/CI_PROJECT/apr2_zain_work_build/structure-from-sherds-pp/class/genetic_algorithm.h)

I will implement a two-tier collision system (Broad-Phase + Narrow-Phase).

**1. Data Structures**
- Add `vector<MatrixXd> collision_clouds_` to the `GeneticAssembler` class to store lightweight surface representations for each sherd.

**2. Initialization (`PrecomputeCollisionClouds`)**
During the constructor or `Run()` initialization:
- **Voxel-Grid Subsampling**: For each sherd, I will iterate through `sur_in_` and `sur_out_` points and pick a representative point every 5mm (voxel-based) to ensure uniform coverage.
- **Edge Filtering**: Any surface point within 2.0mm of the `edge_line_` (the fracture boundary) will be **excluded** from the collision cloud.
- **Result**: A "Collision Cloud" of ~50-100 points that represents the "belly" of the shard but leaves the edges "empty" so they can touch neighbors without penalty.

**3. Fitness Evaluation (`EvaluateFitness`)**
    - Transform the lightweight `collision_clouds_` for both sherds into the global assembly space.
    - Check the distance between these points.
    - **Penalty**: If any point from Shard A is closer than a threshold (e.g., 2.0mm) to a point in Shard B, apply a heavy `kOverlapPenalty`.

**4. Parameter Updates**
- Define `kCollisionPointEpsilon = 2.0` (mm).
- Define `kVoxelSize = 5.0` (mm).

---

## User Review Required

> [!IMPORTANT]
> **Edge Filtering Distance**: I am setting the "Safety Zone" around the fracture edge to 2.0mm. If your pottery has extremely thin fractures, we might need to reduce this to 1.0mm to avoid missing shallow collisions.

---
