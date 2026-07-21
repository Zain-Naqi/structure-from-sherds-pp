#pragma once

// Check # 01: Having both pragma once and this is redundant so the below two lines can be removed apparently
#ifndef _GENETIC_ALGORITHM_H_
#define _GENETIC_ALGORITHM_H_

#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#ifdef _WIN32
#include <direct.h>
#define mkdir(path, mode) _mkdir(path)
#endif
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <random>
#include <queue>
#include <array>
#include <map>
#include <set>
#include <cstdlib>
#include <ctime>
#include <limits>

// Check # 02: #include "data_structure.h" should also work
#include "../class/data_structure.h"
#include "../class/KDTree.h"
#include "../class/reconstruction.h"

// extern: there is a shard_on_off array somewhere in the codebase, I want to use it here, but I'm not defining it
extern bool shard_on_off[];

class GeneticAssembler {
public:
    bool use_inlier_score        = true;
    bool use_connectivity_reward = true;
    bool use_component_penalty   = true;
    bool use_cycle_penalty       = false;  // Always 0 after BFS repair (no cycles in a tree)
    bool use_edge_residual       = false;  // Always 0 after BFS repair (no cycles in a tree)
    bool use_rot_residual        = false;  // Always 0 after BFS repair (no cycles in a tree)
    bool use_overlap_penalty     = true;

    // Diagnostic controls: keep expensive logging/tests off by default.
    bool enable_debug_logging = false;
    bool enable_swap_diagnostics = false;
    bool enable_pose_debug_logging = false;
    int max_swap_diagnostics = 3;

    struct Chromosome {
        vector<int> genes;
        double fitness;
        int root_shard = -1;  // BFS root for repair; valid shard index
    };

    struct FitnessBreakdown {
        double inlier_reward = 0.0;
        double neighbor_penalty = 0.0;
        double active_pair_range_penalty = 0.0;
        double pair_choice_penalty = 0.0;
        double cycle_penalty = 0.0;
        double edge_residual_penalty = 0.0;
        double edge_rot_residual_penalty = 0.0;
        double overlap_penalty = 0.0;
        double connectivity_reward = 0.0;
        double connectivity_component_penalty = 0.0;
        double consensus_reward = 0.0;
        vector<int> sherd_consensus_counts;  // Per-sherd consensus match count
        int active_pair_count = 0;
        int valid_group_count = 0;
        int largest_component = 0;
        int num_components = 0;
        double total_fitness = 0.0;
    };

    struct CollisionCloudStats {
        int sur_in_count = 0;
        int sur_out_count = 0;
        int edge_filtered_count = 0;
        int candidate_count = 0;
        int voxel_count = 0;
        int final_count = 0;
        bool used_fallback = false;
    };

    struct PairOverlapDiagnostic {
        int idx_a = -1;
        int idx_b = -1;
        bool same_component = false;
        bool cloud_a_empty = false;
        bool cloud_b_empty = false;
        bool broad_phase_rejected = false;
        double center_distance = 0.0;
        double broad_radius = 0.0;
        int query_size = 0;
        int target_size = 0;
        int collision_hits = 0;
        double hit_ratio = 0.0;
        double avg_depth = 0.0;
        double min_nn_distance = 0.0;
        double pair_penalty = 0.0;
    };

    //-----------------------------------------------------------------------------------------------------------------//

    GeneticAssembler(const vector<Geom>& shard, list<LCSIndex>& LCS_out, int num_shards, int max_edges = -1)
        : shard_(shard), num_shards_(num_shards), max_edges_(max_edges)
    {

        for (list<LCSIndex>::const_iterator it = LCS_out.begin(); it != LCS_out.end(); ++it) {
            matches_.push_back(*it);
        }
        BuildPairGroups();

        // """
        // BuildPairGroups() populates the pair_groups_ vector like this:
        // pair_groups_[0] → {2, 5, 11}    // all matches between shards (1,3)
        // pair_groups_[1] → {0, 7}        // all matches between shards (2,4)
        // pair_groups_[2] → {3, 6, 9, 14} // all matches between shards (3,7)
        // ...
        // """

        graph_ = MatrixXd::Zero(num_shards_, num_shards_);
        transforms_.resize(num_shards_);
        Matrix3d I = Matrix3d::Identity();
        Vector3d zero = Vector3d::Zero();

        for (int i = 0; i < num_shards_; ++i) {
            transforms_[i].Set(I, zero, i + 1, 1);
        }

        // --- Pre-calculate Sherd Geometry (for Overlap Penalty) ---
        shard_centroids_.resize(num_shards_, Vector3d::Zero());
        shard_radius_.resize(num_shards_, 0.0);

        for (int i = 0; i < num_shards_; ++i) {
            if (!IsShardValidAndOn(i)) continue;

            int total_cols = shard_[i].edge_line_.point_.cols()
               + shard_[i].sur_in_.point_.cols()
               + shard_[i].sur_out_.point_.cols();

            if (total_cols == 0) continue;

            MatrixXd pts(3, total_cols);
            int offset = 0;

            int c1 = shard_[i].edge_line_.point_.cols();
            if (c1 > 0) {
                pts.block(0, offset, 3, c1) = shard_[i].edge_line_.point_;
                offset += c1;
            }
            int c2 = shard_[i].sur_in_.point_.cols();
            if (c2 > 0) {
                pts.block(0, offset, 3, c2) = shard_[i].sur_in_.point_;
                offset += c2;
            }
            int c3 = shard_[i].sur_out_.point_.cols();
            if (c3 > 0) {
                pts.block(0, offset, 3, c3) = shard_[i].sur_out_.point_;
            }

            // Calculate Centroid (Center of Mass)
            Vector3d centroid = Vector3d::Zero();
            for (int p = 0; p < pts.cols(); ++p) {
                centroid += pts.col(p);
            }
            centroid /= static_cast<double>(pts.cols());
            shard_centroids_[i] = centroid;

            // Calculate Adaptive Radius (Mean distance to centroid)
            double avg_dist = 0.0;
            for (int p = 0; p < pts.cols(); ++p) {
                avg_dist += (pts.col(p) - centroid).norm();
            }
            shard_radius_[i] = avg_dist / static_cast<double>(pts.cols());
        }

        PrecomputeCollisionClouds();
        LogCollisionCloudSummary();
    }

    ~GeneticAssembler()
    {
        for (int i = 0; i < num_shards_; ++i) {
            if (i < static_cast<int>(collision_kdtrees_.size()) && collision_kdtrees_[i] != nullptr) {
                kd_free(collision_kdtrees_[i]);
                collision_kdtrees_[i] = nullptr;
            }
        }
    }

    //-----------------------------------------------------------------------------------------------------------------//

    void Run(const MatrixXd& GT_graph, const vector<Trans>& GT_trans, const vector<Trans>& T_axis)
    {
        rng_.seed(rng_seed_);
        if (enable_debug_logging) {
            cout << "[GA DEBUG] RNG seed: " << rng_seed_ << endl;
        }

        if (pair_groups_.empty()) {
            cout << "PAIR GROUPS ARE EMPTY...WHY?" << endl;
            BuildOutputsFromSelection(vector<int>());
            return;
        }

        if (enable_debug_logging && !matches_.empty()) {
            cout << "[GA DEBUG] First Match: shard_x_ = " << matches_[0].shard_x_
            << " shard_y_ = " << matches_[0].shard_y_ << endl;
        }

        if (enable_debug_logging) {
            cout << "[GA DEBUG] Chromosome Length (Total Connections): " << pair_groups_.size() << endl;
        }

        InitializePopulation();
        AuditGroundTruthCandidates(GT_graph, GT_trans, T_axis);

        Chromosome overall_best;
        overall_best.fitness = -1e18;
        int best_seed_idx = -1;
        int best_seed_generations = 0;

        for (int seed_idx = 0; seed_idx < kNumSeeds; ++seed_idx) {
            std::mt19937::result_type seed = static_cast<std::mt19937::result_type>(42 + seed_idx * 7);
            rng_.seed(seed);
            cout << "[GA] === Seed " << seed_idx + 1 << " / " << kNumSeeds << " (seed=" << seed << ") ===" << endl;

            InitializePopulation();
            EvaluatePopulation();
            sort(population_.begin(), population_.end(), [](const Chromosome& a, const Chromosome& b) {
                return a.fitness > b.fitness;
            });
            if (!population_.empty()) {
                FitnessBreakdown best_breakdown;
                EvaluateFitness(population_.front(), &best_breakdown);
                sherd_consensus_counts_ = best_breakdown.sherd_consensus_counts;
            }
            double best_fitness_so_far = -1e9;
            if (!population_.empty()) {
                best_fitness_so_far = population_.front().fitness;
            }
            int generations_since_improvement = 0;
            current_mutation_rate_ = kBaseMutationRate;
            std::uniform_real_distribution<double> real_dist(0.0, 1.0);
            int final_generation = 0;
            for (int generation = 0; ; ++generation) {
                final_generation = generation;

                cout << "[GA] Generation " << generation << " best fitness: " << population_.front().fitness << endl;

                vector<Chromosome> next_population;
                next_population.reserve(kPopulationSize);

                int elitism = min(kElitismCount, static_cast<int>(population_.size()));
                for (int i = 0; i < elitism; ++i) {
                    next_population.push_back(population_[i]);
                }

                while (static_cast<int>(next_population.size()) < kPopulationSize) {
                    const Chromosome& parent1 = TournamentSelect();
                    const Chromosome& parent2 = TournamentSelect();

                    Chromosome child1 = Crossover(parent1, parent2);
                    int mutated_gene1;
                    if (real_dist(rng_) < 0.30) {
                        mutated_gene1 = NeighborhoodMutate(child1);
                    } else {
                        mutated_gene1 = Mutate(child1);
                    }

                    next_population.push_back(child1);

                    if (static_cast<int>(next_population.size()) < kPopulationSize) {
                        Chromosome child2 = Crossover(parent2, parent1);
                        int mutated_gene2;
                        if (real_dist(rng_) < 0.30) {
                            mutated_gene2 = NeighborhoodMutate(child2);
                        } else {
                            mutated_gene2 = Mutate(child2);
                        }
                        next_population.push_back(child2);
                    }
                }

                population_ = next_population;

                EvaluatePopulation();
                sort(population_.begin(), population_.end(), [](const Chromosome& a, const Chromosome& b) {
                    return a.fitness > b.fitness;
                });
                // Cache per-sherd consensus counts from the best individual for guided mutation
                if (!population_.empty()) {
                    FitnessBreakdown best_breakdown;
                    EvaluateFitness(population_.front(), &best_breakdown);
                    sherd_consensus_counts_ = best_breakdown.sherd_consensus_counts;
                }
                // --- Plateau Detection ---
                if (!population_.empty()) {
                    double current_best_fitness = population_.front().fitness;
                    if (current_best_fitness > best_fitness_so_far + 1e-5) {
                        best_fitness_so_far = current_best_fitness;
                        generations_since_improvement = 0;
                        if (current_mutation_rate_ != kBaseMutationRate) {
                            cout << "[GA] Best fitness improved to " << best_fitness_so_far
                                 << ". Resetting mutation rate to " << kBaseMutationRate << endl;
                            current_mutation_rate_ = kBaseMutationRate;
                        }
                    } else {
                        generations_since_improvement++;
                        if (generations_since_improvement >= kStagnationThreshold) {
                            if (current_mutation_rate_ != kHyperMutationRate) {
                                cout << "[GA] Stagnation detected (" << generations_since_improvement
                                     << " generations). Raising mutation rate to " << kHyperMutationRate << endl;
                                current_mutation_rate_ = kHyperMutationRate;
                            }
                        }
                        if (generations_since_improvement >= kEarlyTerminationThreshold) {
                            cout << "[GA] Hard stagnation detected (" << generations_since_improvement
                                 << " generations). Early termination at generation " << generation << endl;
                            if (kEnableSnapshots) {
                                SaveAssemblySnapshot(generation, population_.front());
                            }
                            break;
                        }
                    }
                }
                if (generation % 10 == 0 && !population_.empty()) {
                    LogBestChromosomeBreakdown(generation, population_.front());
                }


                if (seed_idx == 0 && generation == 0) {
                    mkdir(kResultFolder.c_str(), 0777);
                    if (kEnableConvergenceLog) {
                        string filename = kResultFolder + "/ga_convergence.csv";
                        remove(filename.c_str()); // Start fresh
                    }
                }

                if (kEnableConvergenceLog) {
                    LogConvergenceCSV(generation);
                }

                if (kEnableSnapshots && (generation % kSnapshotInterval == 0)) {
                    SaveAssemblySnapshot(generation, population_.front());
                }
            }

            // Population is already evaluated and sorted from the last generation's
            // EvaluatePopulation + sort.

            if (!population_.empty()) {
                double seed_best = population_.front().fitness;
                cout << "[GA] Seed " << seed_idx + 1 << " final best fitness: " << seed_best << endl;
                if (seed_best > overall_best.fitness) {
                    overall_best = population_.front();
                    best_seed_idx = seed_idx;
                    best_seed_generations = final_generation;
                }
            }
        }

        cout << "[GA] Best seed: " << best_seed_idx + 1 << " / " << kNumSeeds
             << " (fitness=" << overall_best.fitness << ")" << endl;

        population_.clear();
        population_.push_back(overall_best);

        LogBestChromosomeBreakdown(best_seed_generations, overall_best);
        AuditChromosomeCollisions(overall_best, "FINAL BEST");
        PrintConsensusDiagnostics(overall_best, "FINAL BEST");

        const Chromosome& best = population_.front();
        const bool run_pair_diagnostics = enable_debug_logging || enable_swap_diagnostics;
        FitnessBreakdown selected_breakdown;
        double selected_total_fitness = 0.0;
        int active_count = 0;
        int swap_tests_ran = 0;

        if (enable_swap_diagnostics) {
            selected_total_fitness = EvaluateFitness(best, &selected_breakdown);
        }

        auto DensityScore = [](const LCSIndex& lcs) {
            return static_cast<double>(lcs.inliner_) / (1.0 + lcs.score_);
        };

        if (run_pair_diagnostics) {
            for (size_t group_idx = 0; group_idx < best.genes.size(); ++group_idx) {
                const vector<size_t>& group = pair_groups_[group_idx];
                if (group.empty()) {
                    continue;
                }

                int pair_x = matches_[group[0]].shard_x_;
                int pair_y = matches_[group[0]].shard_y_;

                int best_local_idx = -1;
                double best_local_density = -1.0;
                for (size_t local = 0; local < group.size(); ++local) {
                    const LCSIndex& candidate = matches_[group[local]];
                    double density = DensityScore(candidate);
                    if (density > best_local_density) {
                        best_local_density = density;
                        best_local_idx = static_cast<int>(local);
                    }
                }

                int choice = best.genes[group_idx];
                if (choice <= 0) {
                    if (enable_debug_logging && best_local_idx >= 0) {
                        const LCSIndex& best_candidate = matches_[group[best_local_idx]];
                        cout << "[GA DEBUG][PAIR " << pair_x << "-" << pair_y << "] "
                            << "inactive"
                            << " | best_idx=" << (best_local_idx + 1)
                            << " inlier=" << best_candidate.inliner_
                            << " score=" << best_candidate.score_
                            << " density=" << best_local_density
                            << endl;
                    }
                    continue;
                }

                int local_idx = choice - 1;
                if (local_idx < 0 || local_idx >= static_cast<int>(group.size())) {
                    if (enable_debug_logging) {
                        cout << "[GA DEBUG][PAIR " << pair_x << "-" << pair_y << "] "
                            << "invalid_choice=" << choice
                            << " (group size=" << group.size() << ")" << endl;
                    }
                    continue;
                }

                const LCSIndex& lcs = matches_[group[local_idx]];
                active_count++;

                if (enable_debug_logging) {
                    cout << "[GA DEBUG] Active match: shard_x_="
                        << lcs.shard_x_
                        << " shard_y_=" << lcs.shard_y_
                        << " inliner_=" << lcs.inliner_
                        << " score_=" << lcs.score_ << endl;

                    cout << "[GA DEBUG][PAIR " << pair_x << "-" << pair_y << "] "
                        << "selected_idx=" << (local_idx + 1)
                        << " / " << group.size()
                        << " density=" << DensityScore(lcs)
                        << " | best_idx=" << (best_local_idx + 1)
                        << " density=" << best_local_density;

                    if (best_local_idx != local_idx) {
                        cout << " (non-best candidate selected)";
                    }
                    cout << endl;
                }

                if (enable_swap_diagnostics && best_local_idx != local_idx && best_local_idx >= 0) {
                    if (swap_tests_ran >= max_swap_diagnostics) {
                        continue;
                    }

                    Chromosome best_candidate_swap = best;
                    best_candidate_swap.genes[group_idx] = best_local_idx + 1;

                    FitnessBreakdown best_candidate_breakdown;
                    double best_candidate_total_fitness = EvaluateFitness(best_candidate_swap, &best_candidate_breakdown);

                    cout << "[GA DEBUG][PAIR " << pair_x << "-" << pair_y << "][SWAP TEST] "
                         << "selected_total=" << selected_total_fitness
                         << " best_candidate_total=" << best_candidate_total_fitness
                         << " delta(best-selected)=" << (best_candidate_total_fitness - selected_total_fitness)
                         << endl;

                    PrintFitnessBreakdown(pair_x, pair_y, "selected_terms", selected_breakdown);
                    PrintFitnessBreakdown(pair_x, pair_y, "best_candidate_terms", best_candidate_breakdown);
                    swap_tests_ran++;
                }
            }

            if (enable_debug_logging) {
                cout << "[GA DEBUG] Total active matches: " << active_count << endl;
            }
        }

        BuildOutputsFromSelection(population_.front().genes, population_.front().root_shard);

    }

    //-----------------------------------------------------------------------------------------------------------------//

    vector<Trans> GetTransforms() const {
        return transforms_;
    }

    //-----------------------------------------------------------------------------------------------------------------//

    MatrixXd GetGraph() const
    {
        return graph_;
    }

    //-----------------------------------------------------------------------------------------------------------------//

    double GetBestFitness() const
    {
        if (population_.empty()) return 0.0;

        double best = population_[0].fitness;

        for (size_t i = 1; i < population_.size(); ++i) {
            if (population_[i].fitness > best) {
                best = population_[i].fitness;
            }
        }

        return best;
    }

    void SeedElites(const vector<Chromosome>& elites)
    {
        seeded_elites_ = elites;
    }

    vector<Chromosome> GetTopChromosomes(int count) const
    {
        int n = min(count, static_cast<int>(population_.size()));
        vector<Chromosome> top(population_.begin(), population_.begin() + n);
        return top;
    }

    //-----------------------------------------------------------------------------------------------------------------//

private:

    void InitializePopulation()
    {
        population_.clear();
        population_.reserve(kPopulationSize);

        for (int i = 0; i < kPopulationSize; ++i) {
            Chromosome chromosome;
            chromosome.genes.resize(pair_groups_.size(), 0);
            chromosome.fitness = 0.0;
            chromosome.root_shard = 0; // Force base sherd to Shard 1 (index 0)

            for (size_t gene_idx = 0; gene_idx < chromosome.genes.size(); ++gene_idx) {
                chromosome.genes[gene_idx] = SampleGroupChoice(gene_idx);
            }
            population_.push_back(chromosome);
        }

        // Apply seeded elites if they exist
        if (!seeded_elites_.empty()) {
            int num_elites = min(static_cast<int>(seeded_elites_.size()), kPopulationSize / 2);
            for (int i = 0; i < num_elites; ++i) {
                int target_idx = kPopulationSize - 1 - i;  // Replace from the back
                population_[target_idx] = seeded_elites_[i];

                // Add exactly one random pair to the chromosome to fill the new edge budget
                vector<int> inactive_genes;
                for (size_t g_idx = 0; g_idx < population_[target_idx].genes.size(); ++g_idx) {
                    if (population_[target_idx].genes[g_idx] == 0 && !pair_groups_[g_idx].empty()) {
                        inactive_genes.push_back(static_cast<int>(g_idx));
                    }
                }
                if (!inactive_genes.empty()) {
                    std::uniform_int_distribution<int> dist(0, static_cast<int>(inactive_genes.size()) - 1);
                    int rand_gene = inactive_genes[dist(rng_)];
                    population_[target_idx].genes[rand_gene] = SampleGroupChoice(rand_gene);
                }
            }
        }
    }

    //-----------------------------------------------------------------------------------------------------------------//

    void EvaluatePopulation()
    {
        const int pop_size = static_cast<int>(population_.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int i = 0; i < pop_size; ++i) {
            RepairChromosome(population_[i]);
            population_[i].fitness = EvaluateFitness(population_[i]);
        }
    }

    //-----------------------------------------------------------------------------------------------------------------//

    void PrintFitnessBreakdown(int pair_x, int pair_y, const string& tag, const FitnessBreakdown& breakdown) const
    {
        cout << "[GA DEBUG][PAIR " << pair_x << "-" << pair_y << "][" << tag << "] "
             << "total=" << breakdown.total_fitness
             << " inlier_reward=" << breakdown.inlier_reward
             << " pair_choice_pen=" << breakdown.pair_choice_penalty
             << " range_pen=" << breakdown.active_pair_range_penalty
             << " neighbor_pen=" << breakdown.neighbor_penalty
             << " cycle_pen=" << breakdown.cycle_penalty
             << " edge_res_pen=" << breakdown.edge_residual_penalty
             << " edge_rot_pen=" << breakdown.edge_rot_residual_penalty
             << " overlap_pen=" << breakdown.overlap_penalty
             << " conn_reward=" << breakdown.connectivity_reward
             << " conn_comp_pen=" << breakdown.connectivity_component_penalty
             << " active_pairs=" << breakdown.active_pair_count
             << " valid_groups=" << breakdown.valid_group_count
             << " largest_comp=" << breakdown.largest_component
             << " num_comp=" << breakdown.num_components
             << endl;
    }

    //-----------------------------------------------------------------------------------------------------------------//

    // Maximum Spanning Forest (MST) Repair Operator (Lamarckian): Prune active edges
    // that form cycles, prioritizing keeping the highest-density matches (based on inliers and scores).
    // This makes the genotype honest about what the phenotype (assembly) actually uses.
    void RepairChromosome(Chromosome& chromosome) const
    {
        struct ActiveEdge {
            int u;
            int v;
            int group_idx;
            double weight;
        };

        vector<ActiveEdge> active_edges;
        active_edges.reserve(chromosome.genes.size());

        // Step 1: Collect active edges and compute their density weights
        for (size_t group_idx = 0; group_idx < chromosome.genes.size(); ++group_idx) {
            int choice = chromosome.genes[group_idx];
            if (choice <= 0) continue;

            const vector<size_t>& group = pair_groups_[group_idx];
            int local_idx = choice - 1;
            if (local_idx < 0 || local_idx >= static_cast<int>(group.size())) continue;

            const LCSIndex& lcs = matches_[group[local_idx]];
            int x = lcs.shard_x_ - 1;
            int y = lcs.shard_y_ - 1;
            if (!IsShardValidAndOn(x) || !IsShardValidAndOn(y)) continue;

            double weight = static_cast<double>(lcs.inliner_) / (1.0 + lcs.score_);
            active_edges.push_back({x, y, static_cast<int>(group_idx), weight});
        }

        // Step 2: Sort active edges in descending order of weight (Maximum Spanning Tree)
        std::sort(active_edges.begin(), active_edges.end(), [](const ActiveEdge& a, const ActiveEdge& b) {
            return a.weight > b.weight;
        });

        // Step 3: Run Kruskal's algorithm using a Disjoint Set Union (DSU) helper
        struct DSU {
            vector<int> parent;
            DSU(int n) {
                parent.resize(n);
                for (int i = 0; i < n; ++i) parent[i] = i;
            }
            int find(int i) {
                if (parent[i] == i)
                    return i;
                return parent[i] = find(parent[i]); // Path compression
            }
            bool unite(int i, int j) {
                int root_i = find(i);
                int root_j = find(j);
                if (root_i != root_j) {
                    parent[root_i] = root_j;
                    return true;
                }
                return false;
            }
        };

        DSU dsu(num_shards_);
        vector<bool> is_tree_edge(pair_groups_.size(), false);
        int edge_count = 0;

        for (const auto& edge : active_edges) {
            if (max_edges_ >= 0 && edge_count >= max_edges_) {
                break;
            }
            if (dsu.unite(edge.u, edge.v)) {
                is_tree_edge[edge.group_idx] = true;
                edge_count++;
            }
        }

        // Step 4: Prune non-tree active edges (permanently modifies chromosome)
        for (size_t group_idx = 0; group_idx < chromosome.genes.size(); ++group_idx) {
            if (chromosome.genes[group_idx] > 0 && !is_tree_edge[group_idx]) {
                chromosome.genes[group_idx] = 0;
            }
        }
    }

    //-----------------------------------------------------------------------------------------------------------------//

    // Check # 13: The fitness function might need a thorough analysis again. Moreover, why are we not using functions in other files (feature_matchings, ranking_system, etc.) at all? We have only used LCS so far?
    double EvaluateFitness(const Chromosome& chromosome, FitnessBreakdown* breakdown = nullptr) const
    {
        if (breakdown != nullptr) {
            *breakdown = FitnessBreakdown();
            breakdown->sherd_consensus_counts.assign(num_shards_, 0);
        }

        const int kMaxNeighbors = max_neighbors_cap_;
        const int valid_shard_count = valid_shard_count_;
        const int valid_group_count = valid_group_count_;

        if (breakdown != nullptr) {
            breakdown->valid_group_count = valid_group_count;
        }

        double fitness = 0.0;
        vector<int> neighbor_count(num_shards_, 0);
        vector<vector<int>> adjacency(num_shards_);
        vector<int> active_groups;
        active_groups.reserve(pair_groups_.size());
        vector<char> group_active(pair_groups_.size(), 0);
        vector<Matrix4d> group_T_ab(pair_groups_.size(), Matrix4d::Identity());
        vector<double> group_weight(pair_groups_.size(), 0.0);

        double pair_choice_gap_sum = 0.0;
        int active_pair_count = 0;

        for (size_t group_idx = 0; group_idx < chromosome.genes.size(); ++group_idx) {
            int choice = chromosome.genes[group_idx];
            if (choice <= 0) {
                continue;
            }

            const vector<size_t>& group = pair_groups_[group_idx];
            int local_idx = choice - 1;
            if (local_idx < 0 || local_idx >= static_cast<int>(group.size())) {
                continue;
            }

            const LCSIndex& lcs = matches_[group[local_idx]];
            int x = lcs.shard_x_ - 1;
            int y = lcs.shard_y_ - 1;
            if (!IsShardValidAndOn(x) || !IsShardValidAndOn(y)) {
                continue;
            }

            // Check # 09: This fitness function might be incorrect. A higher inliner is better while a lower score is better which means that both contradict each other.
            // fitness += static_cast<double>(lcs.inliner_) * lcs.score_;

            // Implementing the above Check # 09:
            // Higher inliers = better, lower score = better
            // Invert score so both terms pull in the same direction
            double selected_density = 0.0;
            if (use_inlier_score) {
                double score_weight = 1.0 / (1.0 + lcs.score_);
                // Apply logarithmic scaling to compress inlier count differences
                selected_density = kInlierScale * std::log(static_cast<double>(lcs.inliner_) + 1.0) * score_weight;
                fitness += selected_density;
                if (breakdown != nullptr) {
                    breakdown->inlier_reward += selected_density;
                }
            }

            active_pair_count++;

            neighbor_count[x]++;
            neighbor_count[y]++;
            adjacency[x].push_back(y);
            adjacency[y].push_back(x);

            if (!group_active[group_idx]) {
                group_active[group_idx] = 1;
                active_groups.push_back(static_cast<int>(group_idx));

                Matrix4d T_xy = Matrix4d::Identity();
                lcs.trans_.Output(T_xy);

                int a = min(x, y);
                int b = max(x, y);
                if (x == a && y == b) {
                    group_T_ab[group_idx] = T_xy;
                } else {
                    group_T_ab[group_idx] = T_xy.inverse();
                }
                group_weight[group_idx] = max(1.0, static_cast<double>(lcs.inliner_));
            }
        }

        // Check # 11: It can be tried to change the threshold or weight (currently -2.0) or re-verify the logic of this cycle consistency code
        if (use_cycle_penalty) {
            double cycle_penalty = 0.0;
            const double cycle_threshold = 30.0;

            for (const array<int, 3>& tri : valid_shard_triplets_) {
                const int a = tri[0];
                const int b = tri[1];
                const int c = tri[2];

                const int group_ab = PairGroupIndex(a, b);
                const int group_bc = PairGroupIndex(b, c);
                const int group_ac = PairGroupIndex(a, c);

                if (group_ab < 0 || group_bc < 0 || group_ac < 0) {
                    continue;
                }
                if (!group_active[group_ab] || !group_active[group_bc] || !group_active[group_ac]) {
                    continue;
                }

                const Matrix4d& T_ab = group_T_ab[group_ab];
                const Matrix4d& T_bc = group_T_ab[group_bc];
                const Matrix4d& T_ac = group_T_ab[group_ac];

                Matrix4d T_composed = T_bc * T_ab;
                Vector3d t_composed = T_composed.block<3, 1>(0, 3);
                Vector3d t_direct = T_ac.block<3, 1>(0, 3);
                double error = (t_composed - t_direct).norm();
                if (error > cycle_threshold) {
                    cycle_penalty += error;
                }
            }

            double penalty = 2.0 * cycle_penalty;
            fitness -= penalty;
            if (breakdown != nullptr) {
                breakdown->cycle_penalty = penalty;
            }
        }

        struct PairEdge {
            int group_idx;
            int a;
            int b;
            Matrix4d T_ab;
            double weight;
        };

        vector<PairEdge> pair_edges;
        pair_edges.reserve(active_groups.size());
        for (int group_idx : active_groups) {
            if (group_idx < 0 || group_idx >= static_cast<int>(pair_groups_.size())) {
                continue;
            }
            if (!group_active[group_idx]) {
                continue;
            }

            if (group_idx >= static_cast<int>(group_rep_x_.size()) || group_idx >= static_cast<int>(group_rep_y_.size())) {
                continue;
            }

            int x = group_rep_x_[group_idx];
            int y = group_rep_y_[group_idx];
            if (!IsShardValidAndOn(x) || !IsShardValidAndOn(y)) {
                continue;
            }

            int a = min(x, y);
            int b = max(x, y);
            pair_edges.push_back({ group_idx, a, b, group_T_ab[group_idx], group_weight[group_idx] });
        }

        vector<vector<pair<int, Matrix4d>>> pose_adj(num_shards_);
        for (const PairEdge& edge : pair_edges) {
            pose_adj[edge.a].push_back(make_pair(edge.b, edge.T_ab.inverse()));
            pose_adj[edge.b].push_back(make_pair(edge.a, edge.T_ab));
        }

        double edge_residual_penalty = 0.0;
        double edge_rot_residual_penalty = 0.0;

        if (use_edge_residual || use_rot_residual) {
            vector<double> max_edge_residuals(pair_groups_.size(), 0.0);
            vector<double> max_rot_residuals(pair_groups_.size(), 0.0);

            for (int root : valid_shard_indices_) {
                vector<bool> pose_visited(num_shards_, false);
                vector<Matrix4d> T_pose(num_shards_, Matrix4d::Identity());
                queue<int> q;
                q.push(root);
                pose_visited[root] = true;

                while (!q.empty()) {
                    int current = q.front();
                    q.pop();

                    for (size_t i = 0; i < pose_adj[current].size(); ++i) {
                        int next = pose_adj[current][i].first;
                        const Matrix4d& T_next_current = pose_adj[current][i].second;
                        if (pose_visited[next]) {
                            continue;
                        }

                        T_pose[next] = T_pose[current] * T_next_current;
                        pose_visited[next] = true;
                        q.push(next);
                    }
                }

                for (const PairEdge& edge : pair_edges) {
                    if (!pose_visited[edge.a] || !pose_visited[edge.b]) {
                        continue;
                    }

                    Matrix4d T_pred_ab = T_pose[edge.b].inverse() * T_pose[edge.a];
                    Vector3d t_pred = T_pred_ab.block<3, 1>(0, 3);
                    Vector3d t_ab = edge.T_ab.block<3, 1>(0, 3);
                    double residual = (t_pred - t_ab).norm();

                    if (residual > kEdgeResidualThreshold) {
                        double excess = residual - kEdgeResidualThreshold;
                        max_edge_residuals[edge.group_idx] = max(max_edge_residuals[edge.group_idx], excess * excess);
                    }

                    Matrix4d T_diff = T_pred_ab * edge.T_ab.inverse();
                    double trace = T_diff(0, 0) + T_diff(1, 1) + T_diff(2, 2);
                    double rot_residual = acos(std::max(-1.0, std::min(1.0, (trace - 1.0) * 0.5)));

                    if (rot_residual > kEdgeRotResidualThreshold) {
                        double excess = rot_residual - kEdgeRotResidualThreshold;
                        max_rot_residuals[edge.group_idx] = max(max_rot_residuals[edge.group_idx], excess);
                    }
                }
            }

            for (int group_idx : active_groups) {
                edge_residual_penalty += max_edge_residuals[group_idx];
                edge_rot_residual_penalty += max_rot_residuals[group_idx];
            }
        }

        if (use_edge_residual) {
            double penalty = kEdgeResidualPenalty * edge_residual_penalty;
            fitness -= penalty;
            if (breakdown != nullptr) {
                breakdown->edge_residual_penalty = penalty;
            }
        }
        if (use_rot_residual) {
            double penalty = kEdgeRotResidualPenalty * edge_rot_residual_penalty;
            fitness -= penalty;
            if (breakdown != nullptr) {
                breakdown->edge_rot_residual_penalty = penalty;
            }
        }

        // Build active tree adjacency lookup for overlap & consensus checking
        vector<vector<bool>> active_adj(num_shards_, vector<bool>(num_shards_, false));
        for (const auto& edge : pair_edges) {
            active_adj[edge.a][edge.b] = true;
            active_adj[edge.b][edge.a] = true;
        }

        // --- Step 1: BFS Component Assembly (Required for both metrics) ---
        vector<bool> vis_global(num_shards_, false);

        struct ComponentAssembly {
            vector<int> placed;
            vector<Matrix4d> T_comp;
        };
        vector<ComponentAssembly> components;

        auto assemble_component_bfs = [&](int start_node) {
            if (vis_global[start_node]) return;

            ComponentAssembly comp;
            comp.T_comp.resize(num_shards_, Matrix4d::Identity());

            queue<int> q;
            q.push(start_node);
            vis_global[start_node] = true;
            comp.placed.push_back(start_node);

            while (!q.empty()) {
                int curr = q.front();
                q.pop();

                for (const auto& edge : pose_adj[curr]) {
                    if (vis_global[edge.first]) continue;

                    comp.T_comp[edge.first] = comp.T_comp[curr] * edge.second;
                    vis_global[edge.first] = true;
                    comp.placed.push_back(edge.first);
                    q.push(edge.first);
                }
            }
            components.push_back(comp);
        };

        // Run component assembly starting at the root shard, then other valid active shards
        assemble_component_bfs(chromosome.root_shard);
        for (int start_node : valid_shard_indices_) {
            assemble_component_bfs(start_node);
        }

        // --- Step 2: Overlap Penalty (Physical occupancy check) ---
        double overlap_penalty_sum = 0.0;
        if (use_overlap_penalty) {
            vector<double> sherd_accumulated_overlap(num_shards_, 0.0);

            for (const auto& comp : components) {
                for (size_t i = 0; i < comp.placed.size(); ++i) {
                    for (size_t j = i + 1; j < comp.placed.size(); ++j) {
                        int idx1 = comp.placed[i];
                        int idx2 = comp.placed[j];

                        if (collision_clouds_[idx1].empty() || collision_clouds_[idx2].empty()) continue;

                        // Skip active adjacent pairs
                        if (active_adj[idx1][idx2]) continue;

                        double hit_ratio = ComputeCloudOverlapRatioIndexed(
                            idx1, idx2, comp.T_comp[idx1], comp.T_comp[idx2]);

                        sherd_accumulated_overlap[idx1] += hit_ratio;
                        sherd_accumulated_overlap[idx2] += hit_ratio;
                    }
                }
            }

            for (int i = 0; i < num_shards_; ++i) {
                if (sherd_accumulated_overlap[i] >= kOverlapMinHitRatio) {
                    overlap_penalty_sum += kOverlapPenalty * sherd_accumulated_overlap[i];
                }
            }
        }
        if (use_overlap_penalty) {
            fitness -= overlap_penalty_sum;
        }
        if (breakdown != nullptr) {
            breakdown->overlap_penalty = overlap_penalty_sum;
        }

        // --- Step 3: Consensus Support Reward (Candidate match agreement) ---
        double consensus_reward = 0.0;
        if (use_consensus_reward) {
            for (const auto& comp : components) {
                for (int idx : comp.placed) {
                    for (int placed_shard : comp.placed) {
                        if (placed_shard == idx) continue;

                        // Symmetrically skip active tree edges
                        if (active_adj[idx][placed_shard]) continue;

                        int g_idx = PairGroupIndex(idx, placed_shard);
                        if (g_idx < 0) continue;

                        const vector<size_t>& group = pair_groups_[g_idx];
                        for (size_t c_idx : group) {
                            const LCSIndex& lcs = matches_[c_idx];
                            int x = lcs.shard_x_ - 1;
                            int y = lcs.shard_y_ - 1;

                            Matrix4d T_candidate = Matrix4d::Identity();
                            lcs.trans_.Output(T_candidate);

                            Matrix4d T_landed;
                            if (idx == y && placed_shard == x) {
                                T_landed = comp.T_comp[placed_shard] * T_candidate.inverse();
                            } else if (idx == x && placed_shard == y) {
                                T_landed = comp.T_comp[placed_shard] * T_candidate;
                            } else {
                                continue;
                            }

                            Matrix4d T_diff = T_landed * comp.T_comp[idx].inverse();

                            Vector3d t_diff(T_diff(0, 3), T_diff(1, 3), T_diff(2, 3));
                            double t_err_sq = t_diff.squaredNorm();
                            if (t_err_sq >= kConsensusTransThreshold * kConsensusTransThreshold) continue;

                            double trace = T_diff(0, 0) + T_diff(1, 1) + T_diff(2, 2);
                            if (trace > 2.0 * kCosConsensusRotThreshold + 1.0) {
                                 consensus_reward += kConsensusWeight; // Flat reward per matching pair
                                 if (breakdown != nullptr && !breakdown->sherd_consensus_counts.empty()) {
                                     breakdown->sherd_consensus_counts[idx]++;
                                 }
                                 break; // Max one supporting candidate per alternative pair
                             }
                        }
                    }
                }
            }
        }
        if (use_consensus_reward) {
            fitness += consensus_reward;
            if (breakdown != nullptr) {
                breakdown->consensus_reward = consensus_reward;
            }
        }

        int largest_component = 0;
        int num_components = 0;

        AnalyzeConnectedComponents(adjacency, largest_component, num_components);
        if (use_connectivity_reward) {
            double reward = kConnectivityReward * static_cast<double>(largest_component);
            fitness += reward;
            if (breakdown != nullptr) {
                breakdown->connectivity_reward = reward;
            }

            int expected_components = (max_edges_ >= 0) ? (valid_shard_count_ - max_edges_) : 1;
            if (use_component_penalty && num_components > expected_components) {
                double penalty = kConnectivityComponentPenalty * static_cast<double>(num_components - expected_components);
                fitness -= penalty;
                if (breakdown != nullptr) {
                    breakdown->connectivity_component_penalty = penalty;
                }
            }
        }

        if (breakdown != nullptr) {
            breakdown->active_pair_count = active_pair_count;
            breakdown->largest_component = largest_component;
            breakdown->num_components = num_components;
            breakdown->total_fitness = fitness;
        }

        return fitness;
    }

    //-----------------------------------------------------------------------------------------------------------------//

    void AnalyzeConnectedComponents(const vector<vector<int>>& adjacency, int& largest_component, int& num_components) const
    {
        largest_component = 0;
        num_components = 0;
        vector<bool> visited(num_shards_, false);

        for (int start = 0; start < num_shards_; ++start) {
            if (!IsShardValidAndOn(start) || visited[start]) {
                continue;
            }

            num_components++;
            queue<int> q;
            q.push(start);
            visited[start] = true;

            int component_size = 0;

            while (!q.empty()) {
                int current = q.front();
                q.pop();
                component_size++;

                for (size_t i = 0; i < adjacency[current].size(); ++i) {
                    int next = adjacency[current][i];
                    if (!IsShardValidAndOn(next) || visited[next]) {
                        continue;
                    }

                    visited[next] = true;
                    q.push(next);
                }
            }

            if (component_size > largest_component) {
                largest_component = component_size;
            }
        }

        // Isolated valid sherds with no selected edges should still count as components.
        // Check # 12: This loop might not be required because the above loop also handles those shards not connected with any other shards.
        for (int node = 0; node < num_shards_; ++node) {
            if (!IsShardValidAndOn(node) || visited[node]) {
                continue;
            }

            num_components++;
            visited[node] = true;

            if (largest_component < 1) {
                largest_component = 1;
            }
        }

        if (num_components == 0) {
            largest_component = 0;
        }
    }

    //-----------------------------------------------------------------------------------------------------------------//

    const Chromosome& TournamentSelect() const
    {
        int size = static_cast<int>(population_.size());
        std::uniform_int_distribution<int> dist(0, size - 1);
        int best_index = dist(rng_);

        // Check # 13: The tournament size (2 here) can be experimented with...
        for (int i = 1; i < 2; ++i) {
            int candidate_index = dist(rng_);
            if (population_[candidate_index].fitness > population_[best_index].fitness) {
                best_index = candidate_index;
            }
        }

        return population_[best_index];
    }

    //-----------------------------------------------------------------------------------------------------------------//

    Chromosome Crossover(const Chromosome& parent1, const Chromosome& parent2) const
    {
        Chromosome child;
        child.fitness = 0.0;
        child.genes.resize(parent1.genes.size(), 0);

        child.root_shard = 0; // Force base sherd to Shard 1 (index 0)

        if (child.genes.empty()) return child;

        // Shard-Partition Crossover: Divide shards between parents to preserve local consensus
        vector<int> shard_owner(num_shards_);
        std::uniform_int_distribution<int> coin(0, 1);
        for (int i = 0; i < num_shards_; ++i) {
            shard_owner[i] = coin(rng_);
        }

        for (size_t gene_idx = 0; gene_idx < child.genes.size(); ++gene_idx) {
            int x = group_rep_x_[gene_idx];
            int y = group_rep_y_[gene_idx];

            if (shard_owner[x] == shard_owner[y]) {
                // Internal edge: inherit from the parent who "owns" these shards
                child.genes[gene_idx] = (shard_owner[x] == 0) ? parent1.genes[gene_idx] : parent2.genes[gene_idx];
            } else {
                // Bridge edge: inherit from a random parent
                child.genes[gene_idx] = (coin(rng_) == 0) ? parent1.genes[gene_idx] : parent2.genes[gene_idx];
            }
        }
        return child;
    }

    //-----------------------------------------------------------------------------------------------------------------//

    // Check # 16: There is a high potential in increasing the chances of mutation as based on this code, there are very less. Multiple genes can be mutated at once as well.
    int Mutate(Chromosome& chromosome) const
    {
        std::uniform_real_distribution<double> real_dist(0.0, 1.0);

        // Root mutation disabled; base sherd remains 0
        chromosome.root_shard = 0;

        if (real_dist(rng_) >= current_mutation_rate_) {
            return -1;
        }

        // --- Weighted Sherd Selection (inverse consensus) ---
        vector<double> weights(num_shards_, 0.0);
        double total_weight = 0.0;
        for (int i = 0; i < num_shards_; ++i) {
            if (sherd_incident_groups_[i].empty()) continue;
            int cc = (i < static_cast<int>(sherd_consensus_counts_.size()))
                     ? sherd_consensus_counts_[i] : 0;
            weights[i] = 1.0 / (1.0 + static_cast<double>(cc));
            total_weight += weights[i];
        }

        if (total_weight <= 0.0) return -1;

        // Roulette wheel selection
        double pick = real_dist(rng_) * total_weight;
        int shard_idx = -1;
        double cumulative = 0.0;
        for (int i = 0; i < num_shards_; ++i) {
            if (sherd_incident_groups_[i].empty()) continue;
            cumulative += weights[i];
            if (cumulative >= pick) {
                shard_idx = i;
                break;
            }
        }
        // Fallback in case of floating-point precision issues:
        if (shard_idx == -1) {
            for (int i = num_shards_ - 1; i >= 0; --i) {
                if (!sherd_incident_groups_[i].empty()) {
                    shard_idx = i;
                    break;
                }
            }
        }

        if (shard_idx == -1 || sherd_incident_groups_[shard_idx].empty()) {
            return -1;
        }

        // --- Adaptive Pivot/Tweak Probability ---
        int cc = (shard_idx < static_cast<int>(sherd_consensus_counts_.size()))
                 ? sherd_consensus_counts_[shard_idx] : 0;
        double pivot_swap_prob;
        if (cc == 0) {
            pivot_swap_prob = 0.70;  // 70% pivot swap, 30% alignment tweak
        } else if (cc == 1) {
            pivot_swap_prob = 0.50;  // 50/50
        } else {
            pivot_swap_prob = 0.30;  // 30% pivot swap, 70% alignment tweak
        }

        if (real_dist(rng_) < pivot_swap_prob) {
            // --- Pivot Swap ---
            vector<int> active_incident_genes;
            vector<int> inactive_incident_genes;
            for (int g_idx : sherd_incident_groups_[shard_idx]) {
                if (chromosome.genes[g_idx] > 0) {
                    active_incident_genes.push_back(g_idx);
                } else {
                    inactive_incident_genes.push_back(g_idx);
                }
            }

            if (!active_incident_genes.empty()) {
                // Pick a random active edge to turn off (deactivate)
                std::uniform_int_distribution<int> active_dist(0, static_cast<int>(active_incident_genes.size()) - 1);
                int deactivate_gene = active_incident_genes[active_dist(rng_)];
                chromosome.genes[deactivate_gene] = 0;

                // Pick a random inactive edge to turn on (activate)
                if (!inactive_incident_genes.empty()) {
                    std::uniform_int_distribution<int> inactive_dist(0, static_cast<int>(inactive_incident_genes.size()) - 1);
                    int activate_gene = inactive_incident_genes[inactive_dist(rng_)];
                    chromosome.genes[activate_gene] = SampleGroupChoice(activate_gene);
                    return activate_gene;
                }
                return deactivate_gene;
            } else {
                // If the shard has no active connections, turn a random inactive one on
                std::uniform_int_distribution<int> gene_dist(0, static_cast<int>(sherd_incident_groups_[shard_idx].size()) - 1);
                int gene_idx = static_cast<int>(sherd_incident_groups_[shard_idx][gene_dist(rng_)]);
                chromosome.genes[gene_idx] = SampleGroupChoice(gene_idx);
                return gene_idx;
            }
        }
        else {
            // --- Alignment Tweak ---
            std::uniform_int_distribution<int> gene_dist(0, static_cast<int>(sherd_incident_groups_[shard_idx].size()) - 1);
            int gene_idx = static_cast<int>(sherd_incident_groups_[shard_idx][gene_dist(rng_)]);

            int old_choice = chromosome.genes[gene_idx];
            int new_choice = SampleGroupChoice(gene_idx);
            if (new_choice == old_choice && !pair_groups_[gene_idx].empty()) {
                new_choice = (old_choice == 0) ? 1 : 0;
            }

            chromosome.genes[gene_idx] = new_choice;
            return gene_idx;
        }
    }

    //-----------------------------------------------------------------------------------------------------------------//

    int NeighborhoodMutate(Chromosome& chromosome) const
    {
        std::uniform_real_distribution<double> real_dist(0.0, 1.0);
        chromosome.root_shard = 0;

        if (real_dist(rng_) >= current_mutation_rate_) {
            return -1;
        }

        // 1. Pick focal sherd using inverse-consensus weighting
        vector<double> weights(num_shards_, 0.0);
        double total_weight = 0.0;
        for (int i = 0; i < num_shards_; ++i) {
            if (sherd_incident_groups_[i].empty()) continue;
            int cc = (i < static_cast<int>(sherd_consensus_counts_.size()))
                     ? sherd_consensus_counts_[i] : 0;
            weights[i] = 1.0 / (1.0 + static_cast<double>(cc));
            total_weight += weights[i];
        }
        if (total_weight <= 0.0) return -1;

        double pick = real_dist(rng_) * total_weight;
        int focal = -1;
        double cumulative = 0.0;
        for (int i = 0; i < num_shards_; ++i) {
            if (sherd_incident_groups_[i].empty()) continue;
            cumulative += weights[i];
            if (cumulative >= pick) {
                focal = i;
                break;
            }
        }
        if (focal == -1) {
            for (int i = num_shards_ - 1; i >= 0; --i) {
                if (!sherd_incident_groups_[i].empty()) {
                    focal = i;
                    break;
                }
            }
        }
        if (focal == -1) return -1;

        // 2. Find all ACTIVE edges incident to focal shard and save originals
        vector<int> incident_genes;
        vector<int> original_choices;
        for (int g_idx : sherd_incident_groups_[focal]) {
            if (chromosome.genes[g_idx] > 0) {
                incident_genes.push_back(g_idx);
                original_choices.push_back(chromosome.genes[g_idx]);
            }
        }
        int degree = static_cast<int>(incident_genes.size());
        if (degree == 0) return -1;

        // 3. Deactivate all incident edges
        for (int g : incident_genes) {
            chromosome.genes[g] = 0;
        }

        // 4. Build DSU of the remaining tree (without focal's edges)
        struct DSU {
            vector<int> parent;
            DSU(int n) {
                parent.resize(n);
                for (int i = 0; i < n; ++i) parent[i] = i;
            }
            int find(int i) {
                if (parent[i] == i) return i;
                return parent[i] = find(parent[i]);
            }
            bool unite(int i, int j) {
                int ri = find(i), rj = find(j);
                if (ri != rj) {
                    parent[ri] = rj;
                    return true;
                }
                return false;
            }
        };

        DSU dsu(num_shards_);
        for (size_t g = 0; g < chromosome.genes.size(); ++g) {
            if (chromosome.genes[g] <= 0) continue;
            int gx = group_rep_x_[g];
            int gy = group_rep_y_[g];
            if (gx == focal || gy == focal) continue;
            if (gx >= 0 && gy >= 0) {
                dsu.unite(gx, gy);
            }
        }

        // 5. Collect ALL candidate pair groups incident to focal
        vector<int> all_candidates;
        for (int g_idx : sherd_incident_groups_[focal]) {
            if (pair_groups_[g_idx].empty()) continue;
            all_candidates.push_back(g_idx);
        }

        // Shuffle candidates for randomness
        for (int i = static_cast<int>(all_candidates.size()) - 1; i > 0; --i) {
            std::uniform_int_distribution<int> swap_dist(0, i);
            std::swap(all_candidates[i], all_candidates[swap_dist(rng_)]);
        }

        // 6. Greedily pick candidates that reconnect distinct components
        int edges_added = 0;
        for (int g_idx : all_candidates) {
            if (edges_added >= degree) break;

            int gx = group_rep_x_[g_idx];
            int gy = group_rep_y_[g_idx];
            int other = (gx == focal) ? gy : gx;
            if (other < 0) continue;

            if (dsu.unite(focal, other)) {
                chromosome.genes[g_idx] = SampleGroupChoice(g_idx);
                edges_added++;
            }
        }

        // 7. Fallback: if we couldn't reconnect all components, restore originals
        if (edges_added < degree) {
            for (size_t g = 0; g < chromosome.genes.size(); ++g) {
                int gx = group_rep_x_[g];
                int gy = group_rep_y_[g];
                if (gx == focal || gy == focal) {
                    chromosome.genes[g] = 0;
                }
            }
            for (int i = 0; i < degree; ++i) {
                chromosome.genes[incident_genes[i]] = original_choices[i];
            }
            return -1;
        }

        return incident_genes[0];
    }

    //-----------------------------------------------------------------------------------------------------------------//

    void BuildOutputsFromSelection(const vector<int>& genes, int root_shard = -1)
    {
        graph_ = MatrixXd::Zero(num_shards_, num_shards_);

        struct AdjEdge {
            int to;
            Matrix4d T_to_current;
            double weight;
            int match_idx;
        };

        vector<vector<AdjEdge>> adjacency(num_shards_);

        for (size_t group_idx = 0; group_idx < genes.size(); ++group_idx) {
            int choice = genes[group_idx];
            if (choice <= 0) {
                continue;
            }

            const vector<size_t>& group = pair_groups_[group_idx];
            int local_idx = choice - 1;
            if (local_idx < 0 || local_idx >= static_cast<int>(group.size())) {
                continue;
            }

            const LCSIndex& lcs = matches_[group[local_idx]];
            int x = lcs.shard_x_ - 1;
            int y = lcs.shard_y_ - 1;
            if (!IsShardValidAndOn(x) || !IsShardValidAndOn(y)) {
                continue;
            }

            graph_(x, y) = 1;
            graph_(y, x) = 1;

            Matrix4d T_xy = Matrix4d::Identity();
            lcs.trans_.Output(T_xy);
            Matrix4d T_yx = T_xy.inverse();
            double edge_weight = max(1.0, static_cast<double>(lcs.inliner_));
            int match_idx = static_cast<int>(group[local_idx]);

            adjacency[x].push_back({ y, T_yx, edge_weight, match_idx });
            adjacency[y].push_back({ x, T_xy, edge_weight, match_idx });
        }

        vector<bool> visited(num_shards_, false);
        vector<Matrix4d> T_to_root(num_shards_, Matrix4d::Identity());

        // Progressive BFS lambda
        auto build_progressive_bfs = [&](int root) {
            if (!IsShardValidAndOn(root) || visited[root]) return;

            visited[root] = true;

            // Set up a temporary ranking subgraph for IcpIncGraphAxis
            std::list<LCSIndex> matches_list(matches_.begin(), matches_.end());
            RankingSubgraph graph_icp(matches_list, num_shards_);
            graph_icp.node_[root] = true;
            graph_icp.root_node_ = root + 1;
            graph_icp.ResetMatchedIndex(shard_);

            vector<int> placed_edges;

            queue<int> bfs_q;
            bfs_q.push(root);

            while (!bfs_q.empty()) {
                int curr = bfs_q.front();
                bfs_q.pop();

                for (size_t i = 0; i < adjacency[curr].size(); ++i) {
                    const AdjEdge& edge = adjacency[curr][i];
                    if (visited[edge.to]) continue;

                    // 1. Mark visited and place initially
                    visited[edge.to] = true;
                    T_to_root[edge.to] = T_to_root[curr] * edge.T_to_current;
                    graph_icp.node_[edge.to] = true;
                    placed_edges.push_back(edge.match_idx);

                    // 2. Generate temporary aligned shard geom copies for optimization
                    vector<Geom> temp_shards = shard_;
                    for (int k = 0; k < num_shards_; ++k) {
                        if (graph_icp.node_[k]) {
                            Matrix3d R_k = T_to_root[k].block<3, 3>(0, 0);
                            Vector3d t_k = T_to_root[k].block<3, 1>(0, 3);
                            temp_shards[k].Move(R_k, t_k, true);
                        }
                    }

                    // 3. Optimize the current subgraph using Ceres ICP
                    vector<Matrix3d> R_step(num_shards_, Matrix3d::Identity());
                    vector<Vector3d> t_step(num_shards_, Vector3d::Zero());
                    vector<RankingSubgraph> dummy_pregraph;
                    int inlier_step = 0;
                    bool rim_restrain = true, axis_restrain = true;

                    // Populate priority_list_ with currently placed edges to prevent out-of-bounds in graph.EdgeOut()
                    graph_icp.priority_list_.clear();
                    Chunk dummy_chunk;
                    dummy_chunk.i_edge = placed_edges;
                    graph_icp.priority_list_.push_back(dummy_chunk);
                    graph_icp.priority_index_ = 0;

                    // IcpIncGraphAxis refines the relative poses of all placed pieces
                    IcpIncGraphAxis(
                        temp_shards,
                        R_step,
                        t_step,
                        graph_icp,
                        dummy_pregraph,
                        inlier_step,
                        rim_restrain,
                        axis_restrain
                    );

                    // 4. Update accumulated transforms with delta optimizations
                    for (int k = 0; k < num_shards_; ++k) {
                        if (graph_icp.node_[k]) {
                            Matrix4d Ti = Matrix4d::Identity();
                            Ti.block<3, 3>(0, 0) = R_step[k];
                            Ti.block<3, 1>(0, 3) = t_step[k];
                            T_to_root[k] = Ti * T_to_root[k];
                        }
                    }

                    bfs_q.push(edge.to);
                }
            }
        };

        // Run progressive BFS
        if (root_shard >= 0 && root_shard < num_shards_) {
            build_progressive_bfs(root_shard);
        }
        for (int root = 0; root < num_shards_; ++root) {
            build_progressive_bfs(root);
        }

        if (enable_pose_debug_logging) {
            int reached = 0;
            for (int i = 0; i < num_shards_; i++) {
                if (visited[i]) {
                    reached++;
                }
            }

            cout << "[GA DEBUG] BFS reached " << reached << " / " << num_shards_ << " sherds" << endl;

            for (int i = 0; i < num_shards_; i++) {
                if (!IsShardValidAndOn(i)) continue;
                cout << "[GA DEBUG] T_to_root[" << i << "] translations: "
                << T_to_root[i](0, 3) << " "
                << T_to_root[i](1, 3) << " "
                << T_to_root[i](2, 3) << endl;
            }
        }

        // Apply final optimized transforms to transforms_
        Matrix3d I = Matrix3d::Identity();
        Vector3d zero = Vector3d::Zero();
        for (int i = 0; i < num_shards_; ++i) {
            transforms_[i].Set(I, zero, i + 1, 1);
            if (!IsShardValidAndOn(i)) continue;
            transforms_[i].Set(T_to_root[i], i + 1, 1);
        }
    }

    //-----------------------------------------------------------------------------------------------------------------//

    bool IsShardValidAndOn(int shard_idx) const
    {
        return (shard_idx >= 0) && (shard_idx < num_shards_) && shard_on_off[shard_idx];
    }

    int PairGroupIndex(int a, int b) const
    {
        if (a < 0 || b < 0 || a >= num_shards_ || b >= num_shards_ || a == b) {
            return -1;
        }
        return pair_group_lookup_[a * num_shards_ + b];
    }

    //-----------------------------------------------------------------------------------------------------------------//

    void BuildPairGroups()
    {
        pair_groups_.clear();
        sherd_incident_groups_.clear();
        sherd_incident_groups_.resize(num_shards_);
        map<pair<int, int>, vector<size_t>> grouped;

        for (size_t i = 0; i < matches_.size(); ++i) {
            int x = matches_[i].shard_x_ - 1;
            int y = matches_[i].shard_y_ - 1;
            if (x < 0 || x >= num_shards_ || y < 0 || y >= num_shards_) {
                continue;
            }
            if (x == y) {
                continue;
            }

            pair<int, int> key = make_pair(min(x, y), max(x, y));
            grouped[key].push_back(i);
        }

        for (map<pair<int, int>, vector<size_t>>::const_iterator it = grouped.begin(); it != grouped.end(); ++it) {
            pair_groups_.push_back(it->second);

            if (it->second.empty()) {
                continue;
            }

            const LCSIndex& representative = matches_[it->second[0]];
            int x = representative.shard_x_ - 1;
            int y = representative.shard_y_ - 1;
            if (x >= 0 && x < num_shards_ && y >= 0 && y < num_shards_) {
                size_t group_idx = pair_groups_.size() - 1;
                sherd_incident_groups_[x].push_back(group_idx);
                sherd_incident_groups_[y].push_back(group_idx);
            }
        }

        PrecomputeStaticCaches();
    }

    //-----------------------------------------------------------------------------------------------------------------//

    void PrecomputeStaticCaches()
    {
        max_neighbors_cap_ = 0;
        valid_shard_indices_.clear();
        valid_shard_triplets_.clear();
        pair_group_lookup_.assign(num_shards_ * num_shards_, -1);

        for (int i = 0; i < num_shards_; ++i) {
            if (IsShardValidAndOn(i)) {
                valid_shard_indices_.push_back(i);
            }
        }
        valid_shard_count_ = static_cast<int>(valid_shard_indices_.size());

        for (size_t i = 0; i < valid_shard_indices_.size(); ++i) {
            for (size_t j = i + 1; j < valid_shard_indices_.size(); ++j) {
                for (size_t k = j + 1; k < valid_shard_indices_.size(); ++k) {
                    valid_shard_triplets_.push_back({
                        valid_shard_indices_[i],
                        valid_shard_indices_[j],
                        valid_shard_indices_[k]
                    });
                }
            }
        }

        vector<int> incident_count(num_shards_, 0);

        for (size_t i = 0; i < matches_.size(); ++i) {
            int x = matches_[i].shard_x_ - 1;
            int y = matches_[i].shard_y_ - 1;
            if (x < 0 || x >= num_shards_ || y < 0 || y >= num_shards_ || x == y) {
                continue;
            }

            incident_count[x]++;
            incident_count[y]++;
        }

        for (int i = 0; i < num_shards_; ++i) {
            max_neighbors_cap_ = max(max_neighbors_cap_, incident_count[i]);
        }
        max_neighbors_cap_ = min(max_neighbors_cap_, max(0, num_shards_ - 1));

        group_best_density_.assign(pair_groups_.size(), 0.0);
        group_rep_x_.assign(pair_groups_.size(), -1);
        group_rep_y_.assign(pair_groups_.size(), -1);
        valid_group_count_ = 0;

        for (size_t group_idx = 0; group_idx < pair_groups_.size(); ++group_idx) {
            const vector<size_t>& group = pair_groups_[group_idx];
            if (group.empty()) {
                continue;
            }

            const LCSIndex& representative = matches_[group[0]];
            group_rep_x_[group_idx] = representative.shard_x_ - 1;
            group_rep_y_[group_idx] = representative.shard_y_ - 1;

            int gx = group_rep_x_[group_idx];
            int gy = group_rep_y_[group_idx];
            if (gx >= 0 && gx < num_shards_ && gy >= 0 && gy < num_shards_) {
                int a = min(gx, gy);
                int b = max(gx, gy);
                pair_group_lookup_[a * num_shards_ + b] = static_cast<int>(group_idx);
                pair_group_lookup_[b * num_shards_ + a] = static_cast<int>(group_idx);
            }
            if (IsShardValidAndOn(gx) && IsShardValidAndOn(gy)) {
                valid_group_count_++;
            }

            double best_density = 0.0;
            for (size_t local_idx = 0; local_idx < group.size(); ++local_idx) {
                const LCSIndex& lcs = matches_[group[local_idx]];
                double density = static_cast<double>(lcs.inliner_) / (1.0 + lcs.score_);
                if (density > best_density) {
                    best_density = density;
                }
            }

            group_best_density_[group_idx] = best_density;
        }
    }

    //-----------------------------------------------------------------------------------------------------------------//

    void PrecomputeCollisionClouds()
    {
        for (size_t i = 0; i < collision_kdtrees_.size(); ++i) {
            if (collision_kdtrees_[i] != nullptr) {
                kd_free(collision_kdtrees_[i]);
                collision_kdtrees_[i] = nullptr;
            }
        }

        collision_clouds_.assign(num_shards_, vector<Vector3d>());
        collision_kdtrees_.assign(num_shards_, nullptr);
        collision_cloud_centroids_.assign(num_shards_, Vector3d::Zero());
        collision_cloud_radii_.assign(num_shards_, 0.0);
        collision_cloud_stats_.assign(num_shards_, CollisionCloudStats());

        const double edge_exclusion_sq = kCollisionEdgeExclusion * kCollisionEdgeExclusion;

        for (int shard_idx = 0; shard_idx < num_shards_; ++shard_idx) {
            vector<Vector3d> candidate_points;
            CollisionCloudStats stats;

            const MatrixXd& sur_in_pts = shard_[shard_idx].sur_in_.point_;
            stats.sur_in_count = sur_in_pts.cols();
            for (int c = 0; c < sur_in_pts.cols(); ++c) {
                Vector3d p = sur_in_pts.col(c);
                if (!IsNearEdge(shard_idx, p, edge_exclusion_sq)) {
                    candidate_points.push_back(p);
                } else {
                    stats.edge_filtered_count++;
                }
            }

            const MatrixXd& sur_out_pts = shard_[shard_idx].sur_out_.point_;
            stats.sur_out_count = sur_out_pts.cols();
            for (int c = 0; c < sur_out_pts.cols(); ++c) {
                Vector3d p = sur_out_pts.col(c);
                if (!IsNearEdge(shard_idx, p, edge_exclusion_sq)) {
                    candidate_points.push_back(p);
                } else {
                    stats.edge_filtered_count++;
                }
            }
            stats.candidate_count = static_cast<int>(candidate_points.size());

            // Safety fallback: if filtering is too strict, fall back to all available surface points.
            if (candidate_points.empty()) {
                stats.used_fallback = true;
                for (int c = 0; c < sur_in_pts.cols(); ++c) {
                    candidate_points.push_back(sur_in_pts.col(c));
                }
                for (int c = 0; c < sur_out_pts.cols(); ++c) {
                    candidate_points.push_back(sur_out_pts.col(c));
                }
                stats.candidate_count = static_cast<int>(candidate_points.size());
            }

            vector<Vector3d> cloud = VoxelDownsample(candidate_points, kCollisionVoxelSize);
            stats.voxel_count = static_cast<int>(cloud.size());
            // CapPointCloudSize(cloud, kCollisionCloudMaxPoints);
            stats.final_count = static_cast<int>(cloud.size());
            collision_cloud_stats_[shard_idx] = stats;

            if (cloud.empty()) {
                collision_clouds_[shard_idx] = cloud;
                collision_cloud_centroids_[shard_idx] = Vector3d::Zero();
                collision_cloud_radii_[shard_idx] = 0.0;
                continue;
            }

            Vector3d centroid = Vector3d::Zero();
            for (size_t i = 0; i < cloud.size(); ++i) {
                centroid += cloud[i];
            }
            centroid /= static_cast<double>(cloud.size());

            double radius = 0.0;
            for (size_t i = 0; i < cloud.size(); ++i) {
                radius = max(radius, (cloud[i] - centroid).norm());
            }

            collision_clouds_[shard_idx] = cloud;
            collision_cloud_centroids_[shard_idx] = centroid;
            collision_cloud_radii_[shard_idx] = radius;

            kdtree* tree = kd_create(3);
            if (tree != nullptr) {
                for (size_t i = 0; i < cloud.size(); ++i) {
                    kd_insert3(tree, cloud[i](0), cloud[i](1), cloud[i](2), nullptr);
                }
                collision_kdtrees_[shard_idx] = tree;
            }
        }
    }

    void LogCollisionCloudSummary() const
    {
        stringstream ss;
        ss << "\n--- COLLISION CLOUD SUMMARY ---" << endl;
        ss << "  Params: voxel=" << kCollisionVoxelSize
           << " edge_exclusion=" << kCollisionEdgeExclusion
           << " epsilon=" << kCollisionPointEpsilon
           << " max_points=" << kCollisionCloudMaxPoints << endl;

        int valid_shards = 0;
        int empty_clouds = 0;
        for (int shard_idx = 0; shard_idx < num_shards_; ++shard_idx) {
            if (!IsShardValidAndOn(shard_idx)) {
                continue;
            }
            valid_shards++;

            const CollisionCloudStats& stats = collision_cloud_stats_[shard_idx];
            if (stats.final_count == 0) {
                empty_clouds++;
            }

            ss << "  Sherd " << (shard_idx + 1)
               << ": in=" << stats.sur_in_count
               << " out=" << stats.sur_out_count
               << " filtered=" << stats.edge_filtered_count
               << " candidate=" << stats.candidate_count
               << " voxel=" << stats.voxel_count
               << " final=" << stats.final_count
               << " radius=" << fixed << setprecision(3) << collision_cloud_radii_[shard_idx];
            if (stats.used_fallback) {
                ss << " [FALLBACK]";
            }
            if (stats.final_count == 0) {
                ss << " [EMPTY]";
            }
            ss << endl;
        }

        ss << "  Valid sherds=" << valid_shards
           << " empty_clouds=" << empty_clouds << endl;
        LogDiagnostic(ss.str());
    }

    bool IsNearEdge(int shard_idx, const Vector3d& p, double edge_exclusion_sq) const
    {
        if (shard_idx < 0 || shard_idx >= num_shards_) {
            return false;
        }

        const MatrixXd& edge_pts = shard_[shard_idx].edge_line_.point_;
        if (edge_pts.cols() == 0) {
            return false;
        }

        for (int c = 0; c < edge_pts.cols(); ++c) {
            if ((p - edge_pts.col(c)).squaredNorm() <= edge_exclusion_sq) {
                return true;
            }
        }

        return false;
    }

    vector<Vector3d> VoxelDownsample(const vector<Vector3d>& input, double voxel_size) const
    {
        if (input.empty()) {
            return vector<Vector3d>();
        }
        if (voxel_size <= 0.0) {
            return input;
        }

        struct VoxelAccum {
            Vector3d sum = Vector3d::Zero();
            int count = 0;
        };

        map<array<int, 3>, VoxelAccum> voxels;
        const double inv_voxel = 1.0 / voxel_size;

        for (size_t i = 0; i < input.size(); ++i) {
            const Vector3d& p = input[i];
            array<int, 3> key = {
                static_cast<int>(floor(p(0) * inv_voxel)),
                static_cast<int>(floor(p(1) * inv_voxel)),
                static_cast<int>(floor(p(2) * inv_voxel))
            };

            VoxelAccum& voxel = voxels[key];
            voxel.sum += p;
            voxel.count++;
        }

        vector<Vector3d> output;
        output.reserve(voxels.size());
        for (const auto& kv : voxels) {
            const VoxelAccum& voxel = kv.second;
            output.push_back(voxel.sum / static_cast<double>(max(1, voxel.count)));
        }
        return output;
    }

    void CapPointCloudSize(vector<Vector3d>& cloud, int max_points) const
    {
        if (max_points <= 0) {
            cloud.clear();
            return;
        }
        if (static_cast<int>(cloud.size()) <= max_points) {
            return;
        }

        vector<Vector3d> reduced;
        reduced.reserve(max_points);

        const double step = static_cast<double>(cloud.size()) / static_cast<double>(max_points);
        for (int i = 0; i < max_points; ++i) {
            size_t idx = static_cast<size_t>(floor(static_cast<double>(i) * step));
            if (idx >= cloud.size()) {
                idx = cloud.size() - 1;
            }
            reduced.push_back(cloud[idx]);
        }

        cloud.swap(reduced);
    }

    double ComputeCloudOverlapPenalty(const vector<Vector3d>& cloud_a, const vector<Vector3d>& cloud_b) const
    {
        return ComputeCloudOverlapDiagnostic(cloud_a, cloud_b).pair_penalty;
    }

    PairOverlapDiagnostic ComputeCloudOverlapDiagnosticIndexed(int idx_a,
                                                               int idx_b,
                                                               const Matrix4d& T_a,
                                                               const Matrix4d& T_b) const
    {
        PairOverlapDiagnostic diag;
        diag.idx_a = idx_a;
        diag.idx_b = idx_b;

        if (idx_a < 0 || idx_a >= num_shards_ || idx_b < 0 || idx_b >= num_shards_) {
            return diag;
        }

        const vector<Vector3d>& cloud_a = collision_clouds_[idx_a];
        const vector<Vector3d>& cloud_b = collision_clouds_[idx_b];
        diag.cloud_a_empty = cloud_a.empty();
        diag.cloud_b_empty = cloud_b.empty();
        diag.query_size = static_cast<int>(cloud_a.size());
        diag.target_size = static_cast<int>(cloud_b.size());
        diag.min_nn_distance = 0.0;

        if (diag.cloud_a_empty || diag.cloud_b_empty ||
            idx_a >= static_cast<int>(collision_kdtrees_.size()) ||
            idx_b >= static_cast<int>(collision_kdtrees_.size()) ||
            collision_kdtrees_[idx_a] == nullptr ||
            collision_kdtrees_[idx_b] == nullptr) {
            return diag;
        }

        // Broad-phase bounding sphere check
        Vector4d c_a_h, c_b_h;
        c_a_h << collision_cloud_centroids_[idx_a], 1.0;
        c_b_h << collision_cloud_centroids_[idx_b], 1.0;
        Vector3d c_a_global = (T_a * c_a_h).head<3>();
        Vector3d c_b_global = (T_b * c_b_h).head<3>();
        double dist = (c_a_global - c_b_global).norm();
        double r_sum = collision_cloud_radii_[idx_a] + collision_cloud_radii_[idx_b];

        diag.center_distance = dist;
        diag.broad_radius = r_sum;
        if (dist > r_sum + kCollisionPointEpsilon) {
            diag.broad_phase_rejected = true;
            diag.pair_penalty = 0.0;
            return diag;
        }
        diag.broad_phase_rejected = false;

        Matrix4d T_b_inv = T_b.inverse();
        Matrix4d T_a_inv = T_a.inverse();

        struct DirectionalStats {
            int query_size = 0;
            int target_size = 0;
            int hits = 0;
            double depth_sum = 0.0;
            double min_nn_distance = 0.0;
        };

        auto QueryDirection = [&](const vector<Vector3d>& query_local,
                                  const Matrix4d& T_query,
                                  const Matrix4d& T_target_inv,
                                  kdtree* target_tree) -> DirectionalStats {
            DirectionalStats stats;
            stats.query_size = static_cast<int>(query_local.size());
            stats.target_size = static_cast<int>(cloud_b.size());
            if (query_local.empty() || target_tree == nullptr) {
                return stats;
            }

            const double epsilon_sq = kCollisionPointEpsilon * kCollisionPointEpsilon;
            double min_best_sq = numeric_limits<double>::max();
            Matrix4d T_combined = T_target_inv * T_query;
            Matrix3d R_combined = T_combined.block<3, 3>(0, 0);
            Vector3d t_combined = T_combined.block<3, 1>(0, 3);

            for (size_t i = 0; i < query_local.size(); ++i) {
                Vector3d q = R_combined * query_local[i] + t_combined;

                kdres* result = kd_nearest3(target_tree, q(0), q(1), q(2));
                if (!result) {
                    continue;
                }

                double np[3] = { 0.0, 0.0, 0.0 };
                kd_res_item(result, np);
                kd_res_free(result);

                double dx = q(0) - np[0];
                double dy = q(1) - np[1];
                double dz = q(2) - np[2];
                double best_sq = dx * dx + dy * dy + dz * dz;

                if (best_sq < min_best_sq) {
                    min_best_sq = best_sq;
                }

                if (best_sq < epsilon_sq) {
                    double depth_ratio = 1.0 - (sqrt(best_sq) / kCollisionPointEpsilon);
                    stats.depth_sum += max(0.0, depth_ratio);
                    stats.hits++;
                }
            }

            if (min_best_sq < numeric_limits<double>::max()) {
                stats.min_nn_distance = sqrt(min_best_sq);
            }
            return stats;
        };

        DirectionalStats ab = QueryDirection(cloud_a, T_a, T_b_inv, collision_kdtrees_[idx_b]);
        DirectionalStats ba = QueryDirection(cloud_b, T_b, T_a_inv, collision_kdtrees_[idx_a]);

        DirectionalStats primary = ab;
        if (ba.hits > ab.hits) {
            primary = ba;
        }

        diag.query_size = max(1, primary.query_size);
        diag.target_size = primary.target_size;
        diag.collision_hits = primary.hits;

        if (ab.min_nn_distance > 0.0 && ba.min_nn_distance > 0.0) {
            diag.min_nn_distance = min(ab.min_nn_distance, ba.min_nn_distance);
        }
        else if (ab.min_nn_distance > 0.0) {
            diag.min_nn_distance = ab.min_nn_distance;
        }
        else {
            diag.min_nn_distance = ba.min_nn_distance;
        }

        const int total_hits = ab.hits + ba.hits;
        if (total_hits <= 0) {
            return diag;
        }

        const double ratio_ab = static_cast<double>(ab.hits) / static_cast<double>(max(1, ab.query_size));
        const double ratio_ba = static_cast<double>(ba.hits) / static_cast<double>(max(1, ba.query_size));
        diag.hit_ratio = max(ratio_ab, ratio_ba);
        diag.avg_depth = (ab.depth_sum + ba.depth_sum) / static_cast<double>(total_hits);

        if (diag.hit_ratio < kOverlapMinHitRatio) {
            diag.pair_penalty = 0.0;
        } else {
            diag.pair_penalty = kOverlapPenalty * diag.hit_ratio;
        }
        return diag;
    }

    PairOverlapDiagnostic ComputeCloudOverlapDiagnostic(const vector<Vector3d>& cloud_a,
                                                        const vector<Vector3d>& cloud_b,
                                                        int idx_a = -1,
                                                        int idx_b = -1) const
    {
        PairOverlapDiagnostic diag;
        diag.idx_a = idx_a;
        diag.idx_b = idx_b;
        diag.cloud_a_empty = cloud_a.empty();
        diag.cloud_b_empty = cloud_b.empty();
        diag.query_size = static_cast<int>(cloud_a.size());
        diag.target_size = static_cast<int>(cloud_b.size());
        diag.min_nn_distance = 0.0;

        if (diag.cloud_a_empty || diag.cloud_b_empty) {
            return diag;
        }

        struct DirectionalStats {
            int query_size = 0;
            int target_size = 0;
            int hits = 0;
            double depth_sum = 0.0;
            double min_nn_distance = 0.0;
        };

        auto ComputeDirectionalStats = [](const vector<Vector3d>& query,
                                          const vector<Vector3d>& target) -> DirectionalStats {
            DirectionalStats stats;
            stats.query_size = static_cast<int>(query.size());
            stats.target_size = static_cast<int>(target.size());
            if (query.empty() || target.empty()) {
                return stats;
            }

            const double epsilon_sq = GeneticAssembler::kCollisionPointEpsilon * GeneticAssembler::kCollisionPointEpsilon;
            double min_best_sq = numeric_limits<double>::max();

            // Build KD-tree on target points for fast nearest-neighbor queries.
            kdtree* tree = kd_create(3);
            if (tree == nullptr) {
                for (size_t i = 0; i < query.size(); ++i) {
                    const Vector3d& q = query[i];
                    double best_sq = numeric_limits<double>::max();

                    for (size_t j = 0; j < target.size(); ++j) {
                        double d_sq = (q - target[j]).squaredNorm();
                        if (d_sq < best_sq) {
                            best_sq = d_sq;
                        }
                    }

                    if (best_sq < min_best_sq) {
                        min_best_sq = best_sq;
                    }

                    if (best_sq < epsilon_sq) {
                        double d = sqrt(best_sq);
                        double depth_ratio = 1.0 - (d / GeneticAssembler::kCollisionPointEpsilon);
                        stats.depth_sum += max(0.0, depth_ratio);
                        stats.hits++;
                    }
                }

                if (min_best_sq < numeric_limits<double>::max()) {
                    stats.min_nn_distance = sqrt(min_best_sq);
                }
                return stats;
            }

            for (size_t j = 0; j < target.size(); ++j) {
                kd_insert3(tree, target[j](0), target[j](1), target[j](2), nullptr);
            }

            for (size_t i = 0; i < query.size(); ++i) {
                kdres* result = kd_nearest3(tree, query[i](0), query[i](1), query[i](2));
                if (result == nullptr) {
                    continue;
                }

                double nearest_pos[3] = { 0.0, 0.0, 0.0 };
                kd_res_item(result, nearest_pos);
                kd_res_free(result);

                const double dx = query[i](0) - nearest_pos[0];
                const double dy = query[i](1) - nearest_pos[1];
                const double dz = query[i](2) - nearest_pos[2];
                const double best_sq = (dx * dx) + (dy * dy) + (dz * dz);

                if (best_sq < min_best_sq) {
                    min_best_sq = best_sq;
                }

                if (best_sq < epsilon_sq) {
                    double d = sqrt(best_sq);
                    double depth_ratio = 1.0 - (d / GeneticAssembler::kCollisionPointEpsilon);
                    stats.depth_sum += max(0.0, depth_ratio);
                    stats.hits++;
                }
            }

            kd_free(tree);

            if (min_best_sq < numeric_limits<double>::max()) {
                stats.min_nn_distance = sqrt(min_best_sq);
            }
            return stats;
        };

        const DirectionalStats ab = ComputeDirectionalStats(cloud_a, cloud_b);
        const DirectionalStats ba = ComputeDirectionalStats(cloud_b, cloud_a);

        DirectionalStats primary = ab;
        if (ba.hits > ab.hits) {
            primary = ba;
        }

        diag.query_size = max(1, primary.query_size);
        diag.target_size = primary.target_size;
        diag.collision_hits = primary.hits;

        if (ab.min_nn_distance > 0.0 && ba.min_nn_distance > 0.0) {
            diag.min_nn_distance = min(ab.min_nn_distance, ba.min_nn_distance);
        }
        else if (ab.min_nn_distance > 0.0) {
            diag.min_nn_distance = ab.min_nn_distance;
        }
        else {
            diag.min_nn_distance = ba.min_nn_distance;
        }

        const int total_hits = ab.hits + ba.hits;
        if (total_hits <= 0) {
            return diag;
        }

        const double ratio_ab = static_cast<double>(ab.hits) / static_cast<double>(max(1, ab.query_size));
        const double ratio_ba = static_cast<double>(ba.hits) / static_cast<double>(max(1, ba.query_size));
        diag.hit_ratio = max(ratio_ab, ratio_ba);
        diag.avg_depth = (ab.depth_sum + ba.depth_sum) / static_cast<double>(total_hits);

        if (diag.hit_ratio < kOverlapMinHitRatio) {
            diag.pair_penalty = 0.0;
        } else {
            diag.pair_penalty = kOverlapPenalty * diag.hit_ratio;
        }
        return diag;
    }

    void CollectChromosomeOverlapDiagnostics(const Chromosome& chromosome,
                                             vector<PairOverlapDiagnostic>& pair_diags,
                                             double& overlap_sum,
                                             int& same_component_pair_count,
                                             int& broad_phase_reject_count,
                                             int& narrow_phase_pair_count,
                                             int& active_collision_pair_count) const
    {
        pair_diags.clear();
        overlap_sum = 0.0;
        same_component_pair_count = 0;
        broad_phase_reject_count = 0;
        narrow_phase_pair_count = 0;
        active_collision_pair_count = 0;

        struct AdjEdge {
            int to;
            Matrix4d T_to_current;
        };

        vector<vector<AdjEdge>> adjacency(num_shards_);
        for (size_t group_idx = 0; group_idx < chromosome.genes.size(); ++group_idx) {
            int choice = chromosome.genes[group_idx];
            if (choice <= 0) {
                continue;
            }

            const vector<size_t>& group = pair_groups_[group_idx];
            int local_idx = choice - 1;
            if (local_idx < 0 || local_idx >= static_cast<int>(group.size())) {
                continue;
            }

            const LCSIndex& lcs = matches_[group[local_idx]];
            int x = lcs.shard_x_ - 1;
            int y = lcs.shard_y_ - 1;
            if (!IsShardValidAndOn(x) || !IsShardValidAndOn(y)) {
                continue;
            }

            Matrix4d T_xy = Matrix4d::Identity();
            lcs.trans_.Output(T_xy);
            adjacency[x].push_back({ y, T_xy.inverse() });
            adjacency[y].push_back({ x, T_xy });
        }

        vector<bool> visited(num_shards_, false);
        vector<Matrix4d> T_to_root(num_shards_, Matrix4d::Identity());
        vector<int> component_id(num_shards_, -1);
        int current_comp_id = 0;

        for (int start_node = 0; start_node < num_shards_; ++start_node) {
            if (!IsShardValidAndOn(start_node) || visited[start_node]) {
                continue;
            }

            queue<int> q;
            q.push(start_node);
            visited[start_node] = true;
            component_id[start_node] = current_comp_id;

            while (!q.empty()) {
                int curr = q.front();
                q.pop();

                for (const auto& edge : adjacency[curr]) {
                    if (visited[edge.to]) {
                        continue;
                    }

                    T_to_root[edge.to] = T_to_root[curr] * edge.T_to_current;
                    visited[edge.to] = true;
                    component_id[edge.to] = current_comp_id;
                    q.push(edge.to);
                }
            }
            current_comp_id++;
        }

        vector<Vector3d> global_centers(num_shards_, Vector3d::Zero());
        for (int i = 0; i < num_shards_; ++i) {
            if (!visited[i] || !IsShardValidAndOn(i)) {
                continue;
            }

            Vector4d c_h;
            c_h << collision_cloud_centroids_[i], 1.0;
            global_centers[i] = (T_to_root[i] * c_h).head<3>();
        }

        for (int i = 0; i < num_shards_; ++i) {
            if (!visited[i] || !IsShardValidAndOn(i)) {
                continue;
            }

            for (int j = i + 1; j < num_shards_; ++j) {
                if (!visited[j] || !IsShardValidAndOn(j)) {
                    continue;
                }

                PairOverlapDiagnostic diag;
                diag.idx_a = i;
                diag.idx_b = j;
                diag.same_component = (component_id[i] >= 0 && component_id[i] == component_id[j]);
                if (!diag.same_component) {
                    pair_diags.push_back(diag);
                    continue;
                }

                // Skip adjacent pairs in the active reconstruction graph to match EvaluateFitness
                bool is_adjacent = false;
                for (const auto& edge : adjacency[i]) {
                    if (edge.to == j) {
                        is_adjacent = true;
                        break;
                    }
                }
                if (is_adjacent) {
                    continue;
                }

                same_component_pair_count++;
                diag.cloud_a_empty = collision_clouds_[i].empty();
                diag.cloud_b_empty = collision_clouds_[j].empty();
                if (diag.cloud_a_empty || diag.cloud_b_empty) {
                    pair_diags.push_back(diag);
                    continue;
                }

                narrow_phase_pair_count++;
                PairOverlapDiagnostic narrow = ComputeCloudOverlapDiagnosticIndexed(i, j, T_to_root[i], T_to_root[j]);
                diag.center_distance = narrow.center_distance;
                diag.broad_radius = narrow.broad_radius;
                diag.broad_phase_rejected = narrow.broad_phase_rejected;

                if (diag.broad_phase_rejected) {
                    broad_phase_reject_count++;
                    pair_diags.push_back(diag);
                    continue;
                }

                diag.query_size = narrow.query_size;
                diag.target_size = narrow.target_size;
                diag.collision_hits = narrow.collision_hits;
                diag.hit_ratio = narrow.hit_ratio;
                diag.avg_depth = narrow.avg_depth;
                diag.min_nn_distance = narrow.min_nn_distance;
                diag.pair_penalty = narrow.pair_penalty;

                pair_diags.push_back(diag);
            }
        }

        // Recompute the accumulated overlap sum using the new sherd-centric logic
        vector<double> sherd_accumulated_overlap(num_shards_, 0.0);
        for (const auto& d : pair_diags) {
            if (d.same_component) {
                // Skip adjacent active edges
                bool is_adjacent = false;
                for (const auto& edge : adjacency[d.idx_a]) {
                    if (edge.to == d.idx_b) {
                        is_adjacent = true;
                        break;
                    }
                }
                if (is_adjacent) continue;

                if (!d.cloud_a_empty && !d.cloud_b_empty && !d.broad_phase_rejected) {
                    sherd_accumulated_overlap[d.idx_a] += d.hit_ratio;
                    sherd_accumulated_overlap[d.idx_b] += d.hit_ratio;
                }
            }
        }

        for (int i = 0; i < num_shards_; ++i) {
            if (sherd_accumulated_overlap[i] >= kOverlapMinHitRatio) {
                overlap_sum += kOverlapPenalty * sherd_accumulated_overlap[i];
            }
        }

        // Update the active collision count for reporting
        for (const auto& d : pair_diags) {
            if (d.same_component && !d.cloud_a_empty && !d.cloud_b_empty && !d.broad_phase_rejected) {
                if (d.hit_ratio > 0.0) {
                    active_collision_pair_count++;
                }
            }
        }
    }

    //-----------------------------------------------------------------------------------------------------------------//

    int BestLocalChoiceForGroup(size_t group_idx) const
    {
        if (group_idx >= pair_groups_.size()) {
            return 0;
        }

        const vector<size_t>& group = pair_groups_[group_idx];
        if (group.empty()) {
            return 0;
        }

        int best_choice = 0;
        double best_density = -1.0;

        for (size_t local_idx = 0; local_idx < group.size(); ++local_idx) {
            const LCSIndex& lcs = matches_[group[local_idx]];
            int x = lcs.shard_x_ - 1;
            int y = lcs.shard_y_ - 1;
            if (!IsShardValidAndOn(x) || !IsShardValidAndOn(y)) {
                continue;
            }

            double density = static_cast<double>(lcs.inliner_) / (1.0 + lcs.score_);
            if (density > best_density) {
                best_density = density;
                best_choice = static_cast<int>(local_idx) + 1;
            }
        }

        return best_choice;
    }

    //-----------------------------------------------------------------------------------------------------------------//

    void EnsureSherdCoverage(Chromosome& chromosome) const
    {
        if (chromosome.genes.size() != pair_groups_.size()) {
            return;
        }

        vector<int> sherd_degree(num_shards_, 0);

        for (size_t group_idx = 0; group_idx < chromosome.genes.size(); ++group_idx) {
            int choice = chromosome.genes[group_idx];
            if (choice <= 0) {
                continue;
            }

            const vector<size_t>& group = pair_groups_[group_idx];
            int local_idx = choice - 1;
            if (local_idx < 0 || local_idx >= static_cast<int>(group.size())) {
                continue;
            }

            const LCSIndex& lcs = matches_[group[local_idx]];
            int x = lcs.shard_x_ - 1;
            int y = lcs.shard_y_ - 1;
            if (!IsShardValidAndOn(x) || !IsShardValidAndOn(y)) {
                continue;
            }

            sherd_degree[x]++;
            sherd_degree[y]++;
        }

        for (int sherd_idx = 0; sherd_idx < num_shards_; ++sherd_idx) {
            if (!IsShardValidAndOn(sherd_idx)) {
                continue;
            }
            if (sherd_degree[sherd_idx] > 0) {
                continue;
            }
            if (sherd_idx >= static_cast<int>(sherd_incident_groups_.size())) {
                continue;
            }

            const vector<size_t>& incident_groups = sherd_incident_groups_[sherd_idx];
            if (incident_groups.empty()) {
                continue;
            }

            size_t best_group_idx = 0;
            int best_choice = 0;
            double best_density = -1.0;
            bool found = false;

            for (size_t i = 0; i < incident_groups.size(); ++i) {
                size_t group_idx = incident_groups[i];
                int candidate_choice = BestLocalChoiceForGroup(group_idx);
                if (candidate_choice <= 0) {
                    continue;
                }

                const vector<size_t>& group = pair_groups_[group_idx];
                int local_idx = candidate_choice - 1;
                if (local_idx < 0 || local_idx >= static_cast<int>(group.size())) {
                    continue;
                }

                const LCSIndex& lcs = matches_[group[local_idx]];
                double density = static_cast<double>(lcs.inliner_) / (1.0 + lcs.score_);
                if (!found || density > best_density) {
                    found = true;
                    best_density = density;
                    best_group_idx = group_idx;
                    best_choice = candidate_choice;
                }
            }

            if (!found) {
                continue;
            }

            int previous_choice = chromosome.genes[best_group_idx];
            chromosome.genes[best_group_idx] = best_choice;

            if (previous_choice <= 0) {
                const vector<size_t>& group = pair_groups_[best_group_idx];
                if (!group.empty()) {
                    const LCSIndex& representative = matches_[group[0]];
                    int x = representative.shard_x_ - 1;
                    int y = representative.shard_y_ - 1;
                    if (IsShardValidAndOn(x) && IsShardValidAndOn(y)) {
                        sherd_degree[x]++;
                        sherd_degree[y]++;
                    }
                }
            }
        }
    }

    //-----------------------------------------------------------------------------------------------------------------//

    void LogDiagnostic(const string& message) const
    {
        static ofstream log_file("ga_diagnostics.log", ios::trunc);
        if (log_file.is_open()) {
            log_file << message << endl;
        }
    }

    void AuditChromosomeCollisions(const Chromosome& chromosome, const string& label) const
    {
        FitnessBreakdown breakdown;
        EvaluateFitness(chromosome, &breakdown);

        vector<PairOverlapDiagnostic> pair_diags;
        double overlap_sum = 0.0;
        int same_component_pairs = 0;
        int broad_phase_rejects = 0;
        int narrow_phase_pairs = 0;
        int active_collision_pairs = 0;

        CollectChromosomeOverlapDiagnostics(chromosome,
                                            pair_diags,
                                            overlap_sum,
                                            same_component_pairs,
                                            broad_phase_rejects,
                                            narrow_phase_pairs,
                                            active_collision_pairs);

        stringstream ss;
        ss << "\n--- COLLISION AUDIT [" << label << "] ---" << endl;
        ss << "  Fitness overlap penalty: " << fixed << setprecision(4) << breakdown.overlap_penalty << endl;
        ss << "  Recomputed indexed overlap sum: " << fixed << setprecision(4) << overlap_sum
           << " (delta=" << fabs(overlap_sum - breakdown.overlap_penalty) << ")" << endl;
        ss << "  Pair counts: total=" << pair_diags.size()
           << " same_comp=" << same_component_pairs
           << " broad_reject=" << broad_phase_rejects
           << " narrow_checked=" << narrow_phase_pairs
           << " active=" << active_collision_pairs << endl;

        bool found_active = false;
        for (size_t i = 0; i < pair_diags.size(); ++i) {
            const PairOverlapDiagnostic& diag = pair_diags[i];
            if (!diag.same_component) {
                continue;
            }

            ss << "  Pair (" << diag.idx_a << "," << diag.idx_b << "):";
            if (diag.cloud_a_empty || diag.cloud_b_empty) {
                ss << " [WARN] EMPTY CLOUD";
                ss << " cloudA=" << (diag.cloud_a_empty ? "EMPTY" : "OK")
                   << " cloudB=" << (diag.cloud_b_empty ? "EMPTY" : "OK") << endl;
                continue;
            }

            ss << " center_dist=" << fixed << setprecision(3) << diag.center_distance
               << " broad_radius=" << diag.broad_radius;

            if (diag.broad_phase_rejected) {
                ss << " [OK] BROAD-PHASE SEPARATED" << endl;
                continue;
            }

            ss << " min_nn=" << diag.min_nn_distance
               << " hits=" << diag.collision_hits << "/" << diag.query_size
               << " hit_ratio=" << diag.hit_ratio
               << " avg_depth=" << diag.avg_depth
               << " pair_pen=" << diag.pair_penalty;

            if (diag.pair_penalty > 1.0e-12) {
                ss << " [!!!] ACTIVE COLLISION UNRESOLVED";
                found_active = true;
            } else {
                ss << " [OK] NARROW-PHASE CLEAR";
            }
            ss << endl;
        }
        if (!found_active) {
            ss << "  No active cloud collisions detected." << endl;
        }
        LogDiagnostic(ss.str());
        cout << ss.str();
    }

    void PrintConsensusDiagnostics(const Chromosome& chromosome, const string& label) const
    {
        stringstream ss;
        ss << "\n=== CONSENSUS SUPPORT DIAGNOSTICS [" << label << "] ===" << endl;

        // BFS to get the active tree and absolute transforms
        vector<bool> vis_global(num_shards_, false);
        vector<Matrix4d> T_comp(num_shards_, Matrix4d::Identity());
        vector<int> parent_shard(num_shards_, -1);
        queue<int> q_comp;
        vector<int> placed;

        int start_node = chromosome.root_shard;
        q_comp.push(start_node);
        vis_global[start_node] = true;
        placed.push_back(start_node);

        struct BFSAdjacent {
            int first;
            Matrix4d second;
        };
        vector<vector<BFSAdjacent>> pose_adj(num_shards_);
        vector<vector<bool>> active_adj(num_shards_, vector<bool>(num_shards_, false));

        for (size_t group_idx = 0; group_idx < chromosome.genes.size(); ++group_idx) {
            int choice = chromosome.genes[group_idx];
            if (choice <= 0) continue;

            const vector<size_t>& group = pair_groups_[group_idx];
            int local_idx = choice - 1;
            if (local_idx < 0 || local_idx >= static_cast<int>(group.size())) continue;

            const LCSIndex& lcs = matches_[group[local_idx]];
            int x = lcs.shard_x_ - 1;
            int y = lcs.shard_y_ - 1;
            if (!IsShardValidAndOn(x) || !IsShardValidAndOn(y)) continue;

            Matrix4d T_xy = Matrix4d::Identity();
            lcs.trans_.Output(T_xy);
            Matrix4d T_yx = T_xy.inverse();

            pose_adj[x].push_back({ y, T_yx });
            pose_adj[y].push_back({ x, T_xy });
            active_adj[x][y] = true;
            active_adj[y][x] = true;
        }

        while (!q_comp.empty()) {
            int curr = q_comp.front();
            q_comp.pop();

            for (const auto& edge : pose_adj[curr]) {
                if (vis_global[edge.first]) {
                    continue;
                }

                T_comp[edge.first] = T_comp[curr] * edge.second;
                parent_shard[edge.first] = curr;
                vis_global[edge.first] = true;
                placed.push_back(edge.first);
                q_comp.push(edge.first);
            }
        }

        ss << "  Active Spanning Tree Component: " << placed.size() << " shards." << endl;
        double total_consensus_reward = 0.0;

        for (int idx : placed) {
            int parent = parent_shard[idx];
            if (parent >= 0) {
                ss << "  Shard " << (idx + 1) << " (placed by parent " << (parent + 1) << "):" << endl;
            } else {
                ss << "  Shard " << (idx + 1) << " (Root Shard):" << endl;
            }

            for (int placed_shard : placed) {
                if (placed_shard == idx) {
                    continue; // Skip self
                }

                bool is_active_edge = active_adj[idx][placed_shard];

                if (is_active_edge) {
                    ss << "    Alternative Pair with Shard " << (placed_shard + 1) << " (ACTIVE IN TREE) [SKIPPED]" << endl;
                    continue; // Skip evaluating active tree edges
                } else {
                    ss << "    Alternative Pair with Shard " << (placed_shard + 1) << " (NOT ACTIVE IN TREE):" << endl;
                }

                int g_idx = PairGroupIndex(idx, placed_shard);
                if (g_idx < 0) {
                    ss << "      [No match candidates found]" << endl;
                    continue;
                }

                const vector<size_t>& group = pair_groups_[g_idx];
                for (size_t c_idx : group) {
                    const LCSIndex& lcs = matches_[c_idx];
                    int x = lcs.shard_x_ - 1;
                    int y = lcs.shard_y_ - 1;

                    Matrix4d T_candidate = Matrix4d::Identity();
                    lcs.trans_.Output(T_candidate);

                    Matrix4d T_landed;
                    if (idx == y && placed_shard == x) {
                        T_landed = T_comp[placed_shard] * T_candidate.inverse();
                    } else if (idx == x && placed_shard == y) {
                        T_landed = T_comp[placed_shard] * T_candidate;
                    } else {
                        continue;
                    }

                    // Measure error
                    Matrix4d T_diff = T_landed * T_comp[idx].inverse();
                    Matrix3d R_diff; Vector3d t_diff, w;
                    for (int r = 0; r < 3; r++) {
                        R_diff.row(r) << T_diff(r, 0), T_diff(r, 1), T_diff(r, 2);
                        t_diff[r] = T_diff(r, 3);
                    }
                    Matrix3d log_R = R_diff.log();
                    w << -log_R(1, 2), log_R(0, 2), -log_R(0, 1);
                    double r_err = w.norm();
                    double t_err = t_diff.norm();

                    bool passes = (r_err < kConsensusRotThreshold && t_err < kConsensusTransThreshold);
                    double added_reward = passes ? kConsensusWeight : 0.0; // Flat reward per matching pair
                    if (passes) {
                        total_consensus_reward += added_reward;
                    }

                    ss << "      Candidate choice=" << (c_idx - group[0] + 1) << " inlier=" << lcs.inliner_
                       << " | RotErr = " << fixed << setprecision(4) << r_err << " rad (" << (r_err * 180.0 / M_PI) << " deg)"
                       << " | TransErr = " << fixed << setprecision(2) << t_err << " mm"
                       << (passes ? " [PASSES]" : " [FAILS]");
                    if (passes) {
                        ss << " -> Reward: +" << fixed << setprecision(2) << added_reward;
                    }
                    ss << endl;
                }
            }
        }
        ss << "  Total Logged Consensus Reward: " << fixed << setprecision(2) << total_consensus_reward << endl;
        ss << "=========================================================" << endl;
        LogDiagnostic(ss.str());
        cout << ss.str();
    }

    void LogBestChromosomeBreakdown(int generation, const Chromosome& best) const
    {
        FitnessBreakdown breakdown;
        EvaluateFitness(best, &breakdown);

        stringstream ss;
        ss << "\n--- GENERATION " << generation << " BEST BREAKDOWN ---" << endl;
        ss << "  Total Fitness: " << fixed << setprecision(2) << best.fitness << endl;
        ss << "  Inlier Reward: " << breakdown.inlier_reward << endl;
        ss << "  Cycle Penalty: " << breakdown.cycle_penalty << endl;
        ss << "  Overlap Penalty: " << breakdown.overlap_penalty << endl;
        ss << "  Edge Residual Penalty: " << breakdown.edge_residual_penalty << endl;
        ss << "  Edge Rot Residual Penalty: " << breakdown.edge_rot_residual_penalty << endl;
        ss << "  Connectivity Reward: " << breakdown.connectivity_reward << endl;
        ss << "  Component Penalty: " << breakdown.connectivity_component_penalty << endl;
        ss << "  Consensus Support Reward: " << breakdown.consensus_reward << endl;
        ss << "  Largest Component: " << breakdown.largest_component << endl;
        ss << "  Num Components: " << breakdown.num_components << endl;
        ss << "  Active Pairs: " << breakdown.active_pair_count << endl;

        vector<PairOverlapDiagnostic> pair_diags;
        double overlap_sum = 0.0;
        int same_component_pairs = 0;
        int broad_phase_rejects = 0;
        int narrow_phase_pairs = 0;
        int active_collision_pairs = 0;
        CollectChromosomeOverlapDiagnostics(best,
                                            pair_diags,
                                            overlap_sum,
                                            same_component_pairs,
                                            broad_phase_rejects,
                                            narrow_phase_pairs,
                                            active_collision_pairs);

        ss << "  OverlapDiag Sum (indexed): " << fixed << setprecision(4) << overlap_sum
           << " (delta_vs_fitness=" << fabs(overlap_sum - breakdown.overlap_penalty) << ")" << endl;
        ss << "  OverlapDiag Pairs: total=" << pair_diags.size()
           << " same_comp=" << same_component_pairs
           << " broad_reject=" << broad_phase_rejects
           << " narrow_checked=" << narrow_phase_pairs
           << " active=" << active_collision_pairs << endl;

        vector<PairOverlapDiagnostic> active_diags;
        for (size_t i = 0; i < pair_diags.size(); ++i) {
            if (pair_diags[i].pair_penalty > 1.0e-12) {
                active_diags.push_back(pair_diags[i]);
            }
        }
        sort(active_diags.begin(), active_diags.end(),
             [](const PairOverlapDiagnostic& a, const PairOverlapDiagnostic& b) {
                 return a.pair_penalty > b.pair_penalty;
             });

        if (active_diags.empty()) {
            ss << "  OverlapDiag Active Pairs: none" << endl;
        } else {
            ss << "  OverlapDiag Top Active Pairs:" << endl;
            int top_k = min(static_cast<int>(active_diags.size()), 5);
            for (int i = 0; i < top_k; ++i) {
                const PairOverlapDiagnostic& diag = active_diags[i];
                ss << "    Pair (" << diag.idx_a << "," << diag.idx_b << ")"
                   << " pen=" << diag.pair_penalty
                   << " hits=" << diag.collision_hits << "/" << diag.query_size
                   << " hit_ratio=" << diag.hit_ratio
                   << " avg_depth=" << diag.avg_depth
                   << " min_nn=" << diag.min_nn_distance << endl;
            }
        }

        ss << "\n  --- ACTIVE GENE MAP ---" << endl;
        for (size_t i = 0; i < best.genes.size(); ++i) {
            int choice = best.genes[i];
            if (choice > 0) {
                const vector<size_t>& group = pair_groups_[i];
                int local_idx = choice - 1;
                if (local_idx >= 0 && local_idx < static_cast<int>(group.size())) {
                    const LCSIndex& lcs = matches_[group[local_idx]];
                    ss << "    Gene " << setw(2) << i << " (Pair " << setw(2) << lcs.shard_x_ << "-" << setw(2) << lcs.shard_y_ << "): choice=" << choice << " (inl=" << lcs.inliner_ << ")" << endl;
                }
            }
        }

        // Specific Debug for Shard 4 (index 3)
        int target_sherd = 3;
        ss << "\n  --- SHARD " << (target_sherd + 1) << " CANDIDATE ANALYSIS ---" << endl;
        for (size_t i = 0; i < pair_groups_.size(); ++i) {
            const vector<size_t>& group = pair_groups_[i];
            if (group.empty()) continue;

            const LCSIndex& rep = matches_[group[0]];
            int x = rep.shard_x_ - 1;
            int y = rep.shard_y_ - 1;

            if (x == target_sherd || y == target_sherd) {
                int selected_choice = best.genes[i];
                ss << "    Pair (" << (x + 1) << "-" << (y + 1) << "): Selected=" << selected_choice << " / " << group.size() << " choices" << endl;
                for (size_t c = 0; c < group.size(); ++c) {
                    const LCSIndex& lcs = matches_[group[c]];
                    ss << "      Choice " << (c + 1) << ": inliners=" << lcs.inliner_ << ( (int(c+1) == selected_choice) ? " [SELECTED]" : "" ) << endl;
                }
            }
        }

        LogDiagnostic(ss.str());
    }

    //-----------------------------------------------------------------------------------------------------------------//

    Chromosome GenerateBiasedReplacement(const Chromosome& survivor)
    {
        Chromosome replacement;
        replacement.genes.resize(survivor.genes.size(), 0);
        replacement.fitness = 0.0;

        struct ActiveGene {
            size_t gene_idx;
            int choice;
            int inliners;
        };
        vector<ActiveGene> active_genes;

        for (size_t i = 0; i < survivor.genes.size(); ++i) {
            int choice = survivor.genes[i];
            if (choice > 0) {
                const vector<size_t>& group = pair_groups_[i];
                int local_idx = choice - 1;
                if (local_idx >= 0 && local_idx < static_cast<int>(group.size())) {
                    int inl = matches_[group[local_idx]].inliner_;
                    active_genes.push_back({i, choice, inl});
                }
            }
        }

        sort(active_genes.begin(), active_genes.end(), [](const ActiveGene& a, const ActiveGene& b) {
            return a.inliners > b.inliners;
        });

        int keep_count = static_cast<int>(active_genes.size() * kBiasInheritRatio);
        vector<bool> inherited(survivor.genes.size(), false);

        for (int i = 0; i < keep_count; ++i) {
            replacement.genes[active_genes[i].gene_idx] = active_genes[i].choice;
            inherited[active_genes[i].gene_idx] = true;
        }

        vector<size_t> random_genes;
        random_genes.reserve(survivor.genes.size() - keep_count);

        for (size_t i = 0; i < survivor.genes.size(); ++i) {
            if (!inherited[i]) {
                replacement.genes[i] = SampleGroupChoice(i);
                random_genes.push_back(i);
            }
        }

        if (!random_genes.empty()) {
            std::uniform_int_distribution<int> dist(0, static_cast<int>(random_genes.size()) - 1);
            size_t repair_idx = random_genes[dist(rng_)];
            // GuidedRepair(replacement, repair_idx);
        }

        // EnsureSherdCoverage(replacement);

        // RepairChromosome(replacement);
        replacement.fitness = EvaluateFitness(replacement);
        return replacement;
    }

    //-----------------------------------------------------------------------------------------------------------------//

    void EnforceDiversity()
    {
        if (population_.empty()) return;
        int chromosome_length = static_cast<int>(population_.front().genes.size());
        if (chromosome_length == 0) return;

        int kDiversityK = chromosome_length / 3;

        for (size_t i = kElitismCount; i < population_.size(); ++i) {
            for (size_t j = i + 1; j < population_.size(); ++j) {
                int hamming_distance = 0;
                for (int g = 0; g < chromosome_length; ++g) {
                    if (population_[i].genes[g] != population_[j].genes[g]) {
                        hamming_distance++;
                        if (hamming_distance >= kDiversityK) {
                            break;
                        }
                    }
                }

                if (hamming_distance < kDiversityK) {
                    population_[j] = GenerateBiasedReplacement(population_[i]);
                }
            }
        }
    }

    //-----------------------------------------------------------------------------------------------------------------//

    int SampleGroupChoice(size_t group_idx) const
    {
        if (group_idx >= pair_groups_.size()) {
            return 0;
        }

        std::uniform_real_distribution<double> dist(0.0, 1.0);

        if (dist(rng_) < kInitialPairInactiveRate) {
            return 0;
        }

        const vector<size_t>& group = pair_groups_[group_idx];
        if (group.empty()) {
            return 0;
        }

        // --- Rank-Based Selection (ranks assigned by density score) ---
        vector<pair<double, int>> ranked_candidates;
        ranked_candidates.reserve(group.size());

        for (size_t i = 0; i < group.size(); ++i) {
            const LCSIndex& lcs = matches_[group[i]];
            // density = inliner / (1 + score)
            double density = static_cast<double>(lcs.inliner_) / (1.0 + lcs.score_);
            ranked_candidates.push_back(make_pair(density, static_cast<int>(i)));
        }

        std::stable_sort(ranked_candidates.begin(), ranked_candidates.end(),
            [](const pair<double, int>& a, const pair<double, int>& b) {
                return a.first > b.first;
            });

        double total_weight = 0.0;
        vector<double> cumulative;
        cumulative.reserve(ranked_candidates.size());
        const int rank_count = static_cast<int>(ranked_candidates.size());
        for (int rank = 0; rank < rank_count; ++rank) {
            // Linear rank weighting: best rank gets the largest weight.
            double w = static_cast<double>(rank_count - rank);
            total_weight += w;
            cumulative.push_back(total_weight);
        }

        if (total_weight <= 0.0) {
            std::uniform_int_distribution<int> choice_dist(1, static_cast<int>(group.size()));
            return choice_dist(rng_);
        }

        double r = dist(rng_) * total_weight;
        for (size_t rank = 0; rank < cumulative.size(); ++rank) {
            if (r <= cumulative[rank]) {
                return ranked_candidates[rank].second + 1;
            }
        }

        return ranked_candidates.front().second + 1;
    }

    void AuditGroundTruthCandidates(const MatrixXd& GT_graph, const vector<Trans>& GT_trans, const vector<Trans>& T_axis)
    {
        cout << "#################### GROUND TRUTH CANDIDATE AUDIT ####################" << endl;

        // Summary tracking for final table
        struct GTEdgeSummary {
            int shard_a, shard_b;
            int num_candidates;
            int best_choice;
            double best_rad_err;
            double best_trans_err;
            bool used_inverse;
            string classification;
        };
        vector<GTEdgeSummary> summaries;

        for (size_t group_idx = 0; group_idx < pair_groups_.size(); ++group_idx) {
            const vector<size_t>& group = pair_groups_[group_idx];
            if (group.empty()) continue;

            const LCSIndex& rep = matches_[group[0]];
            int gx = rep.shard_x_ - 1;
            int gy = rep.shard_y_ - 1;

            int mi = min(gx, gy);
            int ma = max(gx, gy);

            if (GT_graph(mi, ma) == 0) continue;

            cout << "[GT AUDIT] Edge " << (mi + 1) << "-" << (ma + 1) << " (GT=1) Candidates: " << group.size() << endl;

            int best_choice = -1;
            double min_rad_err = 1e9;
            double min_trans_err = 1e9;
            bool best_used_inverse = false;

            for (size_t i = 0; i < group.size(); ++i) {
                const LCSIndex& lcs = matches_[group[i]];
                int x = lcs.shard_x_ - 1;
                int y = lcs.shard_y_ - 1;

                // Build GT relative transform in axis-aligned space (y_aligned -> x_aligned)
                Matrix4d T_ax, T_ay_inv, T_gx_inv, T_gy;
                T_axis[x].Output(T_ax);
                T_axis[y].InvOut(T_ay_inv);
                GT_trans[x].InvOut(T_gx_inv);
                GT_trans[y].Output(T_gy);

                Matrix4d T_gt_wrapped = T_ax * T_gx_inv * T_gy * T_ay_inv;

                Matrix4d T_candidate = Matrix4d::Identity();
                lcs.trans_.Output(T_candidate);

                auto calc_err = [](const Matrix4d& T1, const Matrix4d& T2, double& r_err, double& t_err) {
                    Matrix4d T_diff = T1 * T2.inverse();
                    Matrix3d R_diff; Vector3d t_diff, w;
                    for (int r = 0; r < 3; r++) {
                        R_diff.row(r) << T_diff(r, 0), T_diff(r, 1), T_diff(r, 2);
                        t_diff[r] = T_diff(r, 3);
                    }
                    Matrix3d log_R = R_diff.log();
                    w << -log_R(1, 2), log_R(0, 2), -log_R(0, 1);
                    r_err = w.norm();
                    t_err = t_diff.norm();
                };

                // Try BOTH forward and inverse to handle direction convention mismatches.
                // The LCS transform direction may not match the GT wrapped direction.
                double r_fwd, t_fwd, r_inv, t_inv;
                calc_err(T_gt_wrapped, T_candidate, r_fwd, t_fwd);
                calc_err(T_gt_wrapped, T_candidate.inverse(), r_inv, t_inv);

                // Pick whichever direction gives lower rotation error
                double r_best, t_best;
                bool used_inverse;
                if (r_fwd <= r_inv) {
                    r_best = r_fwd; t_best = t_fwd; used_inverse = false;
                } else {
                    r_best = r_inv; t_best = t_inv; used_inverse = true;
                }

                double density = static_cast<double>(lcs.inliner_) / (1.0 + lcs.score_);

                cout << "  Choice " << (i + 1) << ": inliner=" << setw(3) << lcs.inliner_
                     << " score=" << fixed << setprecision(3) << lcs.score_
                     << " density=" << fixed << setprecision(3) << density
                     << " rad=" << fixed << setprecision(4) << r_best
                     << " trans=" << fixed << setprecision(2) << t_best;
                if (used_inverse) cout << " [INV]";

                // Check if candidate matches original threshold (rad < 0.35 && trans < 50.0)
                if (r_best < 0.35 && t_best < 50.0) {
                    cout << " [VALID GT MATCH]";
                }
                cout << endl;

                if (r_best < min_rad_err) {
                    min_rad_err = r_best;
                    min_trans_err = t_best;
                    best_choice = static_cast<int>(i + 1);
                    best_used_inverse = used_inverse;
                }
            }

            // Classify the best match for this GT edge
            string classification = (min_rad_err < 0.35 && min_trans_err < 50.0) ? "VALID" : "INVALID";

            cout << "  -> Best GT Match: Choice " << best_choice
                 << " (rad=" << fixed << setprecision(4) << min_rad_err
                 << " trans=" << fixed << setprecision(2) << min_trans_err
                 << " " << classification << ")"
                 << (best_used_inverse ? " [INV]" : "") << endl;

            summaries.push_back({mi + 1, ma + 1, static_cast<int>(group.size()),
                                 best_choice, min_rad_err, min_trans_err,
                                 best_used_inverse, classification});
        }

        // Print compact summary table
        cout << "\n--- GT CANDIDATE SUMMARY ---" << endl;
        cout << setw(8) << "Edge" << setw(7) << "Cands" << setw(7) << "Best"
             << setw(10) << "Rad" << setw(10) << "Trans" << setw(6) << "Inv?"
             << setw(12) << "Status" << endl;
        int valid_matches = 0;
        for (size_t s = 0; s < summaries.size(); ++s) {
            const GTEdgeSummary& sm = summaries[s];
            cout << setw(4) << sm.shard_a << "-" << sm.shard_b
                 << setw(7) << sm.num_candidates
                 << setw(7) << sm.best_choice
                 << setw(10) << fixed << setprecision(4) << sm.best_rad_err
                 << setw(10) << fixed << setprecision(2) << sm.best_trans_err
                 << setw(6) << (sm.used_inverse ? "Y" : "N")
                 << setw(12) << sm.classification << endl;
            if (sm.classification == "VALID") {
                valid_matches++;
            }
        }
        cout << "Total GT edges: " << summaries.size()
             << " | Valid GT Matches: " << valid_matches
             << " | Missing/Invalid: " << (summaries.size() - valid_matches) << endl;
        cout << "######################################################################" << endl;

        // Hardcode GA pairs chosen for the Ground Truth audit as requested
        vector<int> gt_genes(pair_groups_.size(), 0);
        int start_node = 0; // Root shard index for BFS component ordering

        auto set_gt_gene = [&](int a, int b, int choice) {
            int g_idx = PairGroupIndex(a - 1, b - 1);
            if (g_idx >= 0) {
                gt_genes[g_idx] = choice;
                cout << "  Spanning Tree Edge: Shard " << a << " - Shard " << b
                     << " (Choice " << choice << ") [HARDCODED]" << endl;
            } else {
                cout << "  [WARNING] Could not find group index for Shard Pair " << a << "-" << b << endl;
            }
        };

        cout << "\n--- SELECTED GT SPANNING TREE EDGES ---" << endl;
        // for pot f
        // set_gt_gene(3, 1, 2);
        // set_gt_gene(6, 1, 1);
        // set_gt_gene(3, 2, 4);
        // set_gt_gene(5, 2, 1);
        // set_gt_gene(4, 5, 1);

        // for pot a
        // set_gt_gene(2, 1, 5);
        // set_gt_gene(3, 1, 1);
        // set_gt_gene(6, 4, 2);
        // set_gt_gene(8, 1, 2);
        // set_gt_gene(4, 1, 3);
        // set_gt_gene(5, 3, 2);
        // set_gt_gene(7, 6, 2);
        // set_gt_gene(9, 6, 1);

        // for pot b
        // set_gt_gene(2, 1, 5);
        // set_gt_gene(3, 1, 1);
        // set_gt_gene(5, 1, 1);
        // set_gt_gene(6, 4, 2);
        // set_gt_gene(9, 1, 2);
        // set_gt_gene(8, 2, 1);
        // set_gt_gene(4, 1, 3);
        // set_gt_gene(9, 7, 1);

        // for pot c
        // set_gt_gene(3, 1, 2);
        // set_gt_gene(3, 2, 11);
        // set_gt_gene(4, 2, 6);

        // for pot g
        // set_gt_gene(1, 3, 2);
        // set_gt_gene(4, 7, 1);
        // set_gt_gene(1, 4, 1);
        // set_gt_gene(2, 5, 1);
        // set_gt_gene(3, 6, 2);
        // set_gt_gene(5, 7, 2);

        // for pot e
        struct ValidGTEdge {
            int a;
            int b;
            int choice;
            int inliers = 0;
        };

        vector<ValidGTEdge> gt_candidates = {
            {1, 2, 4}, {1, 4, 3}, {1, 6, 1}, {1, 14, 1}, {1, 19, 5},
            {2, 4, 1}, {2, 14, 6}, {2, 21, 1},
            {3, 17, 5}, {3, 25, 1}, {3, 26, 1},
            {4, 9, 1}, {4, 18, 3}, {4, 20, 3},
            {5, 8, 5}, {5, 13, 2}, {5, 16, 1}, {5, 22, 3}, {5, 25, 3},
            {6, 8, 5}, {6, 14, 1}, {6, 29, 3},
            {7, 15, 1}, {7, 16, 1}, {7, 22, 6},
            {8, 10, 2}, {8, 16, 2}, {8, 30, 2},
            {9, 11, 1}, {9, 19, 2},
            {10, 15, 1}, {10, 16, 3}, {10, 29, 1},
            {11, 12, 1}, {11, 19, 2}, {11, 27, 2},
            {12, 13, 3}, {12, 15, 4}, {12, 22, 6}, {12, 29, 1},
            {13, 25, 1}, {13, 27, 1}, {13, 28, 1},
            {14, 30, 1},
            {17, 18, 1}, {17, 21, 1},
            {18, 21, 2}, {18, 23, 1},
            {23, 24, 1},
            {25, 31, 1},
            {27, 28, 1}, {28, 31, 1}
        };

        // Fill in actual inlier counts at runtime
        for (auto& edge : gt_candidates) {
            int g_idx = PairGroupIndex(edge.a - 1, edge.b - 1);
            if (g_idx >= 0) {
                const vector<size_t>& group = pair_groups_[g_idx];
                int local_idx = edge.choice - 1;
                if (local_idx >= 0 && local_idx < static_cast<int>(group.size())) {
                    edge.inliers = matches_[group[local_idx]].inliner_;
                }
            }
        }

        // Sort edges by inlier count descending
        sort(gt_candidates.begin(), gt_candidates.end(), [](const ValidGTEdge& e1, const ValidGTEdge& e2) {
            return e1.inliers > e2.inliers;
        });

        // Run Kruskal's to find the Maximum Spanning Tree
        struct DSU {
            vector<int> parent;
            DSU(int n) {
                parent.resize(n);
                iota(parent.begin(), parent.end(), 0);
            }
            int find(int i) {
                if (parent[i] == i) return i;
                return parent[i] = find(parent[i]);
            }
            bool unite(int i, int j) {
                int r_i = find(i);
                int r_j = find(j);
                if (r_i != r_j) {
                    parent[r_i] = r_j;
                    return true;
                }
                return false;
            }
        };

        DSU dsu(num_shards_);
        int mst_edges_built = 0;

        cout << "\n=== COMPUTING MAXIMUM SPANNING TREE FOR POT E GT ===" << endl;
        for (const auto& edge : gt_candidates) {
            if (dsu.unite(edge.a - 1, edge.b - 1)) {
                set_gt_gene(edge.a, edge.b, edge.choice);
                mst_edges_built++;
            }
        }
        cout << "Maximum Spanning Tree completed with " << mst_edges_built << " edges." << endl;
        cout << "====================================================\n" << endl;

        cout << "---------------------------------------" << endl;

        Chromosome gt_chromosome;
        gt_chromosome.genes = gt_genes;
        gt_chromosome.root_shard = start_node;
        gt_chromosome.fitness = 0.0;

        FitnessBreakdown breakdown;
        EvaluateFitness(gt_chromosome, &breakdown);
        gt_chromosome.fitness = breakdown.total_fitness;

        cout << "\n==============================================" << endl;
        cout << "       GROUND TRUTH GA FITNESS AUDIT          " << endl;
        cout << "==============================================" << endl;
        cout << "  Total Fitness Score:         " << breakdown.total_fitness << endl;
        cout << "  Inlier Reward (+):           " << breakdown.inlier_reward << endl;
        cout << "  Connectivity Reward (+):     " << breakdown.connectivity_reward << endl;
        cout << "  Consensus Reward (+):        " << breakdown.consensus_reward << endl;
        cout << "  Cycle Penalty (-):           " << breakdown.cycle_penalty << endl;
        cout << "  Overlap Penalty (-):         " << breakdown.overlap_penalty << endl;
        cout << "  Edge Residual Penalty (-):   " << breakdown.edge_residual_penalty << endl;
        cout << "  Edge Rot Residual Penalty (-): " << breakdown.edge_rot_residual_penalty << endl;
        cout << "  Component Penalty (-):       " << breakdown.connectivity_component_penalty << endl;
        cout << "  Active Pairs (Edges):        " << breakdown.active_pair_count << endl;
        cout << "==============================================" << endl;

        // Audit the collisions specifically for the ground truth chromosome
        AuditChromosomeCollisions(gt_chromosome, "GROUND TRUTH");
        PrintConsensusDiagnostics(gt_chromosome, "GROUND TRUTH");

        // Save the GT chromosome to the GA results so main.cpp visualizes it
        BuildOutputsFromSelection(gt_genes, start_node);
        population_.clear();
        population_.push_back(gt_chromosome);
    }

    //-----------------------------------------------------------------------------------------------------------------//

private:
    double ComputeCloudOverlapRatioIndexed(int idx_a,
                                           int idx_b,
                                           const Matrix4d& T_a,
                                           const Matrix4d& T_b) const
    {
        if (idx_a < 0 || idx_a >= num_shards_ || idx_b < 0 || idx_b >= num_shards_) {
            return 0.0;
        }
        if (idx_a >= static_cast<int>(collision_kdtrees_.size()) ||
            idx_b >= static_cast<int>(collision_kdtrees_.size()) ||
            collision_kdtrees_[idx_a] == nullptr ||
            collision_kdtrees_[idx_b] == nullptr) {
            return 0.0;
        }

        const vector<Vector3d>& cloud_a = collision_clouds_[idx_a];
        const vector<Vector3d>& cloud_b = collision_clouds_[idx_b];
        if (cloud_a.empty() || cloud_b.empty()) {
            return 0.0;
        }

        // Broad-phase bounding sphere check
        Vector4d c_a_h, c_b_h;
        c_a_h << collision_cloud_centroids_[idx_a], 1.0;
        c_b_h << collision_cloud_centroids_[idx_b], 1.0;
        Vector3d c_a_global = (T_a * c_a_h).head<3>();
        Vector3d c_b_global = (T_b * c_b_h).head<3>();
        double dist = (c_a_global - c_b_global).norm();
        double r_sum = collision_cloud_radii_[idx_a] + collision_cloud_radii_[idx_b];
        if (dist > r_sum + kCollisionPointEpsilon) {
            return 0.0;
        }

        Matrix4d T_b_inv = T_b.inverse();
        Matrix4d T_a_inv = T_a.inverse();

        struct DirectionalStats {
            int query_size = 0;
            int hits = 0;
            double depth_sum = 0.0;
            double min_nn_distance = 0.0;
        };

        auto QueryDirection = [&](const vector<Vector3d>& query_local,
                                  const Matrix4d& T_query,
                                  const Matrix4d& T_target_inv,
                                  kdtree* target_tree) -> DirectionalStats {
            DirectionalStats stats;
            stats.query_size = static_cast<int>(query_local.size());
            if (query_local.empty() || target_tree == nullptr) {
                return stats;
            }

            const double epsilon_sq = kCollisionPointEpsilon * kCollisionPointEpsilon;
            double min_best_sq = numeric_limits<double>::max();
            Matrix4d T_combined = T_target_inv * T_query;
            Matrix3d R_combined = T_combined.block<3, 3>(0, 0);
            Vector3d t_combined = T_combined.block<3, 1>(0, 3);

            for (size_t i = 0; i < query_local.size(); ++i) {
                Vector3d q = R_combined * query_local[i] + t_combined;

                kdres* result = kd_nearest3(target_tree, q(0), q(1), q(2));
                if (!result) {
                    continue;
                }

                double np[3] = { 0.0, 0.0, 0.0 };
                kd_res_item(result, np);
                kd_res_free(result);

                double dx = q(0) - np[0];
                double dy = q(1) - np[1];
                double dz = q(2) - np[2];
                double best_sq = dx * dx + dy * dy + dz * dz;

                if (best_sq < min_best_sq) {
                    min_best_sq = best_sq;
                }

                if (best_sq < epsilon_sq) {
                    double depth_ratio = 1.0 - (sqrt(best_sq) / kCollisionPointEpsilon);
                    stats.depth_sum += max(0.0, depth_ratio);
                    stats.hits++;
                }
            }

            if (min_best_sq < numeric_limits<double>::max()) {
                stats.min_nn_distance = sqrt(min_best_sq);
            }
            return stats;
        };

        DirectionalStats ab = QueryDirection(cloud_a, T_a, T_b_inv, collision_kdtrees_[idx_b]);
        DirectionalStats ba = QueryDirection(cloud_b, T_b, T_a_inv, collision_kdtrees_[idx_a]);

        const double ratio_ab = static_cast<double>(ab.hits) / static_cast<double>(max(1, ab.query_size));
        const double ratio_ba = static_cast<double>(ba.hits) / static_cast<double>(max(1, ba.query_size));
        return max(ratio_ab, ratio_ba);
    }

    // static std::mt19937::result_type CreateSeed()
    // {
    //     std::random_device rd;
    //     const unsigned int time_seed = static_cast<unsigned int>(std::time(nullptr));
    //     const unsigned int clock_seed = static_cast<unsigned int>(std::clock());
    //     std::seed_seq seed_seq { rd(), rd(), rd(), rd(), time_seed, clock_seed };
    //     std::array<std::mt19937::result_type, 1> seed_data = { 0u };
    //     seed_seq.generate(seed_data.begin(), seed_data.end());
    //     return seed_data[0];
    // }


    static constexpr int kPopulationSize = 350;
    static constexpr int kMaxGenerations = 100;
    static constexpr int kElitismCount = 1;
    static constexpr int kNumSeeds = 1;
    static constexpr double kBaseMutationRate = 0.15;       // Base mutation rate (allows convergence)
    static constexpr double kHyperMutationRate = 0.40;      // Hyper-mutation rate (escapes local traps)
    static constexpr int kStagnationThreshold = 15;         // Generations without improvement to trigger hyper-mutation
    static constexpr int kEarlyTerminationThreshold = 70;   // Generations without improvement to exit early
    static constexpr double kBiasInheritRatio = 0.4;
    static constexpr double kRootMutationRate = 0.15;
    static constexpr double kInitialPairInactiveRate = 0.05;
    static constexpr double kEdgeResidualThreshold = 5.0;
    static constexpr double kEdgeResidualPenalty = 2.0;
    static constexpr double kEdgeRotResidualThreshold = 0.2;
    static constexpr double kEdgeRotResidualPenalty = 5.0;
    static constexpr double kConnectivityReward = 10.0;
    static constexpr double kConnectivityComponentPenalty = 100;

    // Overlap constants
    static constexpr double kOverlapPenalty = 1200.0;
    static constexpr double kOverlapMinHitRatio = 0.01;
    static constexpr double kCollisionPointEpsilon = 1.0;   // mm
    static constexpr double kCollisionEdgeExclusion = 10.0;  // mm
    static constexpr double kCollisionVoxelSize = 3.0;      // mm

    // Consensus Support parameters
    static constexpr bool use_consensus_reward = true;
    static constexpr double kInlierScale = 6.0;             // Multiplier for log-inlier reward
    static constexpr double kConsensusWeight = 20.0;         // Flat reward per supporting match
    static constexpr double kConsensusRotThreshold = 0.22;   // ~12.6 degrees rotation error limit
    static constexpr double kConsensusTransThreshold = 12.0; // 12.0 mm translation error limit

    static constexpr int kSnapshotInterval = 25;
    static constexpr bool kEnableSnapshots = true;
    static constexpr bool kEnableConvergenceLog = true;
    const string kResultFolder = "result_paper";

    void LogConvergenceCSV(int generation)
    {
        string filename = kResultFolder + "/ga_convergence.csv";
        bool exists = false;
        { ifstream f(filename); exists = f.good(); }

        ofstream out(filename, ios::app);
        if (!out) return;

        if (!exists) {
            out << "generation,best_fitness,avg_fitness,worst_fitness,"
                << "inlier_reward,cycle_penalty,edge_residual_penalty,rot_residual_penalty,"
                << "overlap_penalty,connectivity_reward,active_pairs,largest_component,num_components" << endl;
        }

        double sum_fitness = 0;
        for (const auto& c : population_) sum_fitness += c.fitness;
        double avg_fitness = sum_fitness / population_.size();

        const Chromosome& best = population_.front();
        const Chromosome& worst = population_.back();

        FitnessBreakdown bd;
        EvaluateFitness(best, &bd);

        out << generation << "," << best.fitness << "," << avg_fitness << "," << worst.fitness << ","
            << bd.inlier_reward << "," << bd.cycle_penalty << "," << bd.edge_residual_penalty << "," << bd.edge_rot_residual_penalty << ","
            << bd.overlap_penalty << "," << bd.connectivity_reward << "," << bd.active_pair_count << ","
            << bd.largest_component << "," << bd.num_components << endl;
        out.close();
    }

    void SaveAssemblySnapshot(int generation, const Chromosome& best)
    {
        string filename = kResultFolder + "/ga_snapshot_gen_" + to_string(generation) + ".ply";
        ofstream out(filename);
        if (!out) return;

        BuildOutputsFromSelection(best.genes, best.root_shard);

        int total_points = 0;
        for (int i = 0; i < num_shards_; ++i) {
            total_points += shard_[i].edge_line_.point_.cols();
            total_points += shard_[i].sur_in_.point_.cols();
            total_points += shard_[i].sur_out_.point_.cols();
        }

        out << "ply" << endl;
        out << "format ascii 1.0" << endl;
        out << "element vertex " << total_points << endl;
        out << "property float x" << endl;
        out << "property float y" << endl;
        out << "property float z" << endl;
        out << "property float nx" << endl;
        out << "property float ny" << endl;
        out << "property float nz" << endl;
        out << "property int shard_id" << endl;
        out << "end_header" << endl;

        for (int i = 0; i < num_shards_; ++i) {
            Matrix4d T;
            transforms_[i].Output(T);
            Matrix3d R = T.block<3, 3>(0, 0);

            auto write_points = [&](const MatrixXd& pts, const MatrixXd& nms) {
                if (pts.cols() == 0) return;
                for (int j = 0; j < pts.cols(); ++j) {
                    Vector4d p;
                    p << pts.col(j), 1.0;
                    Vector3d p_w = (T * p).head<3>();
                    Vector3d n_w = (R * nms.col(j)).normalized();
                    out << p_w.x() << " " << p_w.y() << " " << p_w.z() << " "
                        << n_w.x() << " " << n_w.y() << " " << n_w.z() << " " << i << endl;
                }
            };

            write_points(shard_[i].edge_line_.point_, shard_[i].edge_line_.normal_);
            write_points(shard_[i].sur_in_.point_, shard_[i].sur_in_.normal_);
            write_points(shard_[i].sur_out_.point_, shard_[i].sur_out_.normal_);
        }
        out.close();
    }

    static constexpr int kCollisionCloudMaxPoints = 300;

    vector<Geom> shard_;
    vector<LCSIndex> matches_;
    vector<vector<size_t>> pair_groups_;
    vector<vector<size_t>> sherd_incident_groups_;
    int num_shards_;

    vector<Chromosome> population_;
    vector<Chromosome> seeded_elites_;

    vector<Trans> transforms_;
    MatrixXd graph_;

    // Adaptive Geometry for overlap checks
    vector<Vector3d> shard_centroids_;
    vector<double> shard_radius_;
    vector<vector<Vector3d>> collision_clouds_;
    vector<kdtree*> collision_kdtrees_;
    vector<CollisionCloudStats> collision_cloud_stats_;
    vector<Vector3d> collision_cloud_centroids_;
    vector<double> collision_cloud_radii_;

    // Static caches to avoid recomputing identical data in every fitness call.
    int max_neighbors_cap_ = 0;
    vector<double> group_best_density_;
    vector<int> group_rep_x_;
    vector<int> group_rep_y_;
    vector<int> valid_shard_indices_;
    vector<array<int, 3>> valid_shard_triplets_;
    vector<int> pair_group_lookup_;
    int valid_shard_count_ = 0;
    int valid_group_count_ = 0;
    vector<int> sherd_consensus_counts_; // Per-sherd consensus match count from best individual
    double current_mutation_rate_ = kBaseMutationRate;

    const double kCosConsensusRotThreshold = cos(kConsensusRotThreshold);

    int max_edges_ = -1;

    // std::mt19937::result_type rng_seed_ = CreateSeed();
    std::mt19937::result_type rng_seed_ = 42;
    mutable std::mt19937 rng_{rng_seed_};
};

#endif
