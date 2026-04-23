#pragma once

// Check # 01: Having both pragma once and this is redundant so the below two lines can be removed apparently
#ifndef _GENETIC_ALGORITHM_H_   
#define _GENETIC_ALGORITHM_H_   

#include <iostream>
#include <fstream>
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

// extern: there is a shard_on_off array somewhere in the codebase, I want to use it here, but I'm not defining it
extern bool shard_on_off[];

class GeneticAssembler {
public:
    bool use_inlier_score        = true;
    bool use_connectivity_reward = true;
    bool use_cycle_penalty       = true;
    bool use_edge_residual       = true;
    bool use_rot_residual        = true;
    bool use_neighbor_penalty    = false;
    bool use_overlap_penalty     = true;
    bool use_pair_choice_penalty = false;
    bool use_active_pair_range_penalty = false;

    // Diagnostic controls: keep expensive logging/tests off by default.
    bool enable_debug_logging = false;
    bool enable_swap_diagnostics = false;
    bool enable_pose_debug_logging = false;
    int max_swap_diagnostics = 3;

    struct Chromosome {
        vector<int> genes;
        double fitness;
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

    GeneticAssembler(const vector<Geom>& shard, list<LCSIndex>& LCS_out, int num_shards) 
        : shard_(shard), num_shards_(num_shards) 
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

            const MatrixXd& pts = shard_[i].edge_line_.point_;
            if (pts.cols() == 0) continue;

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

    //-----------------------------------------------------------------------------------------------------------------//

    void Run(const MatrixXd& GT_graph, const vector<Trans>& GT_trans, const vector<Trans>& T_axis)
    {
        srand(42);

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

        for (int generation = 0; generation < kMaxGenerations; ++generation) {
            EvaluatePopulation();
            sort(population_.begin(), population_.end(), [](const Chromosome& a, const Chromosome& b) {
                return a.fitness > b.fitness;
            });

            cout << "[GA] Generation " << generation << " best fitness: " << population_.front().fitness << endl;

            vector<Chromosome> next_population;
            next_population.reserve(kPopulationSize);

            int elitism = min(kElitismCount, static_cast<int>(population_.size()));
            for (int i = 0; i < elitism; ++i) {
                next_population.push_back(population_[i]);
            }

            while (static_cast<int>(next_population.size()) < kPopulationSize) {
                Chromosome parent1 = TournamentSelect();
                Chromosome parent2 = TournamentSelect();

                Chromosome child1 = Crossover(parent1, parent2);
                int mutated_gene1 = Mutate(child1);

                if (mutated_gene1 >= 0) {
                    // GuidedRepair(child1, static_cast<size_t>(mutated_gene1));
                }
                EnsureSherdCoverage(child1);
                // RepairPositionCollisions(child1);
                next_population.push_back(child1);

                if (static_cast<int>(next_population.size()) < kPopulationSize) {
                    Chromosome child2 = Crossover(parent2, parent1);
                    int mutated_gene2 = Mutate(child2);
                    if (mutated_gene2 >= 0) {
                        // GuidedRepair(child2, static_cast<size_t>(mutated_gene2));
                    }

                    EnsureSherdCoverage(child2);
                    // RepairPositionCollisions(child2);

                    next_population.push_back(child2);
                }
            }

            population_ = next_population;

            EvaluatePopulation();
            sort(population_.begin(), population_.end(), [](const Chromosome& a, const Chromosome& b) {
                return a.fitness > b.fitness;
            });

            if (generation % 10 == 0 && !population_.empty()) {
                // Diagnostic pass with logging for the current best
                // RepairPositionCollisions(population_.front(), true);
                LogBestChromosomeBreakdown(generation, population_.front());
            }

            EnforceDiversity();

            sort(population_.begin(), population_.end(), [](const Chromosome& a, const Chromosome& b) {
                return a.fitness > b.fitness;
            });
        }

        EvaluatePopulation();
        sort(population_.begin(), population_.end(), [](const Chromosome& a, const Chromosome& b) {
            return a.fitness > b.fitness;
        });

        if (!population_.empty()) {
            LogBestChromosomeBreakdown(kMaxGenerations, population_.front());
            AuditChromosomeCollisions(population_.front(), "FINAL BEST");
        }

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

        BuildOutputsFromSelection(population_.front().genes);
        
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

            for (size_t gene_idx = 0; gene_idx < chromosome.genes.size(); ++gene_idx) {
                chromosome.genes[gene_idx] = SampleGroupChoice(gene_idx);
            }

            EnsureSherdCoverage(chromosome);
            // RepairPositionCollisions(chromosome);
            population_.push_back(chromosome);
        }
    }

    //-----------------------------------------------------------------------------------------------------------------//

    void EvaluatePopulation()
    {
        for (size_t i = 0; i < population_.size(); ++i) {
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

    // Check # 13: The fitness function might need a thorough analysis again. Moreover, why are we not using functions in other files (feature_matchings, ranking_system, etc.) at all? We have only used LCS so far?
    double EvaluateFitness(const Chromosome& chromosome, FitnessBreakdown* breakdown = nullptr) const
    {   
        if (breakdown != nullptr) {
            *breakdown = FitnessBreakdown();
        }

        const int kMaxNeighbors = max_neighbors_cap_;

        int valid_shard_count = 0;
        for (int shard_idx = 0; shard_idx < num_shards_; ++shard_idx) {
            if (IsShardValidAndOn(shard_idx)) {
                valid_shard_count++;
            }
        }

        int valid_group_count = 0;
        for (size_t group_idx = 0; group_idx < pair_groups_.size(); ++group_idx) {
            if (group_idx >= group_rep_x_.size() || group_idx >= group_rep_y_.size()) {
                continue;
            }

            int gx = group_rep_x_[group_idx];
            int gy = group_rep_y_[group_idx];
            if (gx < 0 || gy < 0) {
                continue;
            }

            if (!IsShardValidAndOn(gx) || !IsShardValidAndOn(gy)) {
                continue;
            }

            valid_group_count++;
        }

        if (breakdown != nullptr) {
            breakdown->valid_group_count = valid_group_count;
        }

        double fitness = 0.0;
        vector<int> neighbor_count(num_shards_, 0);
        vector<vector<int>> adjacency(num_shards_);
        double pair_choice_gap_sum = 0.0;
        int active_pair_count = 0;

        for (size_t group_idx = 0; group_idx < chromosome.genes.size(); group_idx++) {
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
                selected_density = static_cast<double>(lcs.inliner_) * score_weight;
                fitness += selected_density;
                if (breakdown != nullptr) {
                    breakdown->inlier_reward += selected_density;
                }
            }

            if (use_pair_choice_penalty) {
                double best_density = 0.0;
                if (group_idx < group_best_density_.size()) {
                    best_density = group_best_density_[group_idx];
                }

                if (best_density > kPairChoiceDensityEps) {
                    double gap_ratio = (best_density - selected_density) / best_density;
                    if (gap_ratio > 0.0) {
                        pair_choice_gap_sum += gap_ratio;
                    }
                }
            }

            active_pair_count++;

            neighbor_count[x]++;
            neighbor_count[y]++;
            adjacency[x].push_back(y);
            adjacency[y].push_back(x);
        }

        if (use_neighbor_penalty) {
            for (int shard_idx = 0; shard_idx < num_shards_; ++shard_idx) {
                if (!IsShardValidAndOn(shard_idx)) {
                    continue;
                }

                // Check # 10: This condition might be totally illogical because neighbor_count[shard_idx] can never exceed kMaxNeighbors as per the above implementation. A more useful function migh help here.
                if (neighbor_count[shard_idx] > kMaxNeighbors) {
                    double penalty = 50.0 * static_cast<double>(neighbor_count[shard_idx] - kMaxNeighbors);
                    fitness -= penalty;
                    if (breakdown != nullptr) {
                        breakdown->neighbor_penalty += penalty;
                    }
                }
            }
        }

        if (use_active_pair_range_penalty && valid_group_count > 0) {
            int min_active = static_cast<int>(std::ceil(kMinActivePairRatio * static_cast<double>(valid_group_count)));
            int max_active = static_cast<int>(std::floor(kMaxActivePairRatio * static_cast<double>(valid_group_count)));

            int spanning_min = max(0, valid_shard_count - 1);
            spanning_min = min(spanning_min, valid_group_count);
            min_active = max(min_active, spanning_min);
            max_active = max(max_active, min_active);

            if (active_pair_count < min_active) {
                double deficit = static_cast<double>(min_active - active_pair_count);
                double penalty = kActivePairRangePenaltyWeight * deficit * deficit;
                fitness -= penalty;
                if (breakdown != nullptr) {
                    breakdown->active_pair_range_penalty += penalty;
                }
            }
            else if (active_pair_count > max_active) {
                double excess = static_cast<double>(active_pair_count - max_active);
                double penalty = (kActivePairRangePenaltyWeight * kExcessActivePairPenaltyScale) * excess * excess;
                fitness -= penalty;
                if (breakdown != nullptr) {
                    breakdown->active_pair_range_penalty += penalty;
                }
            }
        }

        if (use_pair_choice_penalty && active_pair_count > 0) {
            double avg_gap = pair_choice_gap_sum / static_cast<double>(active_pair_count);
            double penalty = kPairChoicePenaltyWeight * avg_gap;
            fitness -= penalty;
            if (breakdown != nullptr) {
                breakdown->pair_choice_penalty = penalty;
            }
        }

        map<pair<int, int>, Matrix4d> active_transforms;

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

            int a = min(x, y);
            int b = max(x, y);
            Matrix4d T_ab = Matrix4d::Identity();
            if (x == a && y == b) {
                T_ab = T_xy;
            } else {
                T_ab = T_xy.inverse();
            }

            auto key = make_pair(a, b);
            if (active_transforms.find(key) == active_transforms.end()) {
                active_transforms[key] = T_ab;
            }
        }

        // Check # 11: It can be tried to change the threshold or weight (currently -2.0) or re-verify the logic of this cycle consistency code
        double cycle_penalty = 0.0;
        const double cycle_threshold = 30.0;

        for (int a = 0; a < num_shards_; ++a) {
            if (!IsShardValidAndOn(a)) {
                continue;
            }

            for (int b = a + 1; b < num_shards_; ++b) {
                if (!IsShardValidAndOn(b)) {
                    continue;
                }

                for (int c = b + 1; c < num_shards_; ++c) {
                    if (!IsShardValidAndOn(c)) {
                        continue;
                    }

                    auto key_ab = make_pair(a, b);
                    auto key_bc = make_pair(b, c);
                    auto key_ac = make_pair(a, c);

                    bool has_ab = active_transforms.count(key_ab) > 0;
                    bool has_bc = active_transforms.count(key_bc) > 0;
                    bool has_ac = active_transforms.count(key_ac) > 0;

                    if (!has_ab || !has_bc || !has_ac) {
                        continue;
                    }

                    Matrix4d T_ab = active_transforms[key_ab];
                    Matrix4d T_bc = active_transforms[key_bc];
                    Matrix4d T_ac = active_transforms[key_ac];

                    Matrix4d T_composed = T_bc * T_ab;

                    Vector3d t_composed = T_composed.block<3, 1>(0, 3);
                    Vector3d t_direct = T_ac.block<3, 1>(0, 3);
                    double error = (t_composed - t_direct).norm();

                    if (error > cycle_threshold) {
                        cycle_penalty += error;
                    }

                }
            }
        }

        if (use_cycle_penalty) {
            double penalty = 2.0 * cycle_penalty;
            fitness -= penalty;
            if (breakdown != nullptr) {
                breakdown->cycle_penalty = penalty;
            }
        }

        struct PairEdge {
            int a;
            int b;
            Matrix4d T_ab;
            double weight;
        };

        map<pair<int, int>, PairEdge> pair_edges;

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

            int a = min(x, y);
            int b = max(x, y);
            Matrix4d T_xy = Matrix4d::Identity();
            lcs.trans_.Output(T_xy);

            Matrix4d T_ab = Matrix4d::Identity();
            if (x == a && y == b) {
                T_ab = T_xy;
            } else {
                T_ab = T_xy.inverse();
            }

            auto key = make_pair(a, b);
            auto it = pair_edges.find(key);
            double w = max(1.0, static_cast<double>(lcs.inliner_));
            if (it == pair_edges.end() || w > it->second.weight) {
                pair_edges[key] = { a, b, T_ab, w };
            }
        }

        vector<vector<pair<int, Matrix4d>>> pose_adj(num_shards_);

        for (const auto& kv : pair_edges) {
            const PairEdge& edge = kv.second;
            pose_adj[edge.a].push_back(make_pair(edge.b, edge.T_ab.inverse()));
            pose_adj[edge.b].push_back(make_pair(edge.a, edge.T_ab));
        }

        map<pair<int, int>, double> max_edge_residuals;
        map<pair<int, int>, double> max_rot_residuals;
        double edge_residual_penalty = 0.0;
        double edge_rot_residual_penalty = 0.0;
        
        for (int root = 0; root < num_shards_; ++root) {
            if (!IsShardValidAndOn(root)) {
                continue;
            }

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

            for (const auto& kv : pair_edges) {
                const PairEdge& edge = kv.second;
                if (!pose_visited[edge.a] || !pose_visited[edge.b]) {
                    continue;
                }

                auto key = make_pair(edge.a, edge.b);

                Matrix4d T_pred_ab = T_pose[edge.b].inverse() * T_pose[edge.a];
                Vector3d t_pred = T_pred_ab.block<3, 1>(0, 3);
                Vector3d t_ab = edge.T_ab.block<3, 1>(0, 3);
                double residual = (t_pred - t_ab).norm();

                if (residual > kEdgeResidualThreshold) {
                    double excess = residual - kEdgeResidualThreshold;
                    max_edge_residuals[key] = max(max_edge_residuals[key], excess * excess);
                }

                Matrix4d T_diff = T_pred_ab * edge.T_ab.inverse();
                Matrix3d R_diff;
                for (int r = 0; r < 3; ++r) {
                    R_diff.row(r) << T_diff(r, 0), T_diff(r, 1), T_diff(r, 2);
                }

                Matrix3d log_R = R_diff.log();
                Vector3d w;
                w << -log_R(1, 2), log_R(0, 2), -log_R(0, 1);
                double rot_residual = w.norm();
                
                if (rot_residual > kEdgeRotResidualThreshold) {
                    double excess = rot_residual - kEdgeRotResidualThreshold;
                    max_rot_residuals[key] = max(max_rot_residuals[key], excess);
                }
            }
        }

        for (auto const& [key, val] : max_edge_residuals) {
            edge_residual_penalty += val;
        }
        for (auto const& [key, val] : max_rot_residuals) {
            edge_rot_residual_penalty += val;
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

        // --- Overlap Penalty (Physical occupancy check to handle symmetry) ---
        if (use_overlap_penalty) {
            vector<bool> vis_global(num_shards_, false);
            double overlap_penalty_sum = 0.0;

            for (int start_node = 0; start_node < num_shards_; ++start_node) {
                if (!IsShardValidAndOn(start_node) || vis_global[start_node]) {
                    continue;
                }

                vector<Matrix4d> T_comp(num_shards_, Matrix4d::Identity());
                queue<int> q_comp;
                vector<int> placed;

                q_comp.push(start_node);
                vis_global[start_node] = true;
                placed.push_back(start_node);

                while (!q_comp.empty()) {
                    int curr = q_comp.front();
                    q_comp.pop();

                    for (const auto& edge : pose_adj[curr]) {
                        if (vis_global[edge.first]) {
                            continue;
                        }

                        T_comp[edge.first] = T_comp[curr] * edge.second;
                        vis_global[edge.first] = true;
                        placed.push_back(edge.first);
                        q_comp.push(edge.first);
                    }
                }

                vector<vector<Vector3d>> cloud_global(num_shards_);
                vector<Vector3d> center_global(num_shards_, Vector3d::Zero());
                for (size_t i = 0; i < placed.size(); ++i) {
                    int shard_idx = placed[i];
                    const vector<Vector3d>& cloud_local = collision_clouds_[shard_idx];
                    const Matrix4d& T = T_comp[shard_idx];

                    vector<Vector3d>& transformed_cloud = cloud_global[shard_idx];
                    transformed_cloud.reserve(cloud_local.size());

                    for (size_t p = 0; p < cloud_local.size(); ++p) {
                        Vector4d h;
                        h << cloud_local[p], 1.0;
                        transformed_cloud.push_back((T * h).head<3>());
                    }

                    Vector4d c_h;
                    c_h << collision_cloud_centroids_[shard_idx], 1.0;
                    center_global[shard_idx] = (T * c_h).head<3>();
                }

                for (size_t i = 0; i < placed.size(); ++i) {
                    for (size_t j = i + 1; j < placed.size(); ++j) {
                        int idx1 = placed[i];
                        int idx2 = placed[j];

                        const vector<Vector3d>& cloud1 = cloud_global[idx1];
                        const vector<Vector3d>& cloud2 = cloud_global[idx2];
                        if (cloud1.empty() || cloud2.empty()) {
                            continue;
                        }

                        double broad_radius = collision_cloud_radii_[idx1] + collision_cloud_radii_[idx2] + kCollisionPointEpsilon;
                        double broad_radius_sq = broad_radius * broad_radius;
                        double center_dist_sq = (center_global[idx1] - center_global[idx2]).squaredNorm();
                        if (center_dist_sq > broad_radius_sq) {
                            continue;
                        }

                        overlap_penalty_sum += ComputeCloudOverlapPenalty(cloud1, cloud2);
                    }
                }
            }

            fitness -= overlap_penalty_sum;
            if (breakdown != nullptr) {
                breakdown->overlap_penalty = overlap_penalty_sum;
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

            if (num_components > 1) {
                double penalty = kConnectivityComponentPenalty * static_cast<double>(num_components - 1);
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

    Chromosome TournamentSelect() const
    {
        int size = static_cast<int>(population_.size());
        int best_index = rand() % size;

        // Check # 13: The tournament size (2 here) can be experimented with...
        for (int i = 1; i < 3; ++i) {
            int candidate_index = rand() % size;
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

        if (child.genes.empty()) {
            return child;
        }

        int point = rand() % static_cast<int>(child.genes.size());

        for (int i = 0; i < point; ++i) {
            child.genes[i] = parent1.genes[i];
        }

        for (size_t i = point; i < child.genes.size(); ++i) {
            child.genes[i] = parent2.genes[i];
        }

        return child;
    }

    //-----------------------------------------------------------------------------------------------------------------//

    // Check # 16: There is a high potential in increasing the chances of mutation as based on this code, there are very less. Multiple genes can be mutated at once as well.
    int Mutate(Chromosome& chromosome) const
    {
        int selected_mutated_gene = -1;
        int mutated_count = 0;

        for (size_t i = 0; i < chromosome.genes.size(); ++i) {

            // Check # 14: This might be a poor randomness technique which has to upgraded with moder C++ standards
            double r = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);

            if (r < kMutationRate) {
                int old_choice = chromosome.genes[i];
                const vector<size_t>& group = pair_groups_[i];
                int num_options = 1 + static_cast<int>(group.size()); // 0 (inactive) + candidates

                if (num_options <= 1) {
                    continue;
                }

                int new_choice = old_choice;

                // --- Symmetry Aware Mutation (Symmetry Flip) ---
                // Keep a small chance of random alternative to preserve exploration.
                double symmetry_flip_prob = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
                if (symmetry_flip_prob < kSymmetryFlipRate) {
                    // Forced Symmetry Flip: Pick any other candidate except the current one or 'inactive'
                    // This forces the GA to explore the alternative high-quality matches in the group
                    int num_candidates = static_cast<int>(group.size());
                    if (num_candidates > 1) {
                        //  int draw = rand() % num_candidates;
                        //  new_choice = draw + 1; 
                        //  if (new_choice == old_choice) {
                        //      new_choice = (draw + 1) % num_candidates + 1;
                        //  }
                        while (new_choice == old_choice) {
                            new_choice = SampleGroupChoice(i);
                        }
                    } else {
                        //  new_choice = SampleGroupChoice(i);
                        if (old_choice == 0){
                            new_choice = 1;
                        }
                        else {
                            new_choice = 0;
                        }
                    }
                } 
                // else {
                //     const int kMaxResampleAttempts = 16;
                //     for (int attempt = 0; attempt < kMaxResampleAttempts; ++attempt) {
                //         new_choice = SampleGroupChoice(i);
                //         if (new_choice != old_choice) {
                //             break;
                //         }
                //     }
                // }

                // if (new_choice == old_choice) {
                //     int clamped_old = old_choice;
                //     if (clamped_old < 0 || clamped_old >= num_options) {
                //         clamped_old = 0;
                //     }

                //     int draw = rand() % (num_options - 1);
                //     if (draw >= clamped_old) {
                //         draw++;
                //     }
                //     new_choice = draw;
                // }

                chromosome.genes[i] = new_choice;
                mutated_count++;

                if ((rand() % mutated_count) == 0) {
                    selected_mutated_gene = static_cast<int>(i);
                }
            }
        }

        return selected_mutated_gene;
    }

    //-----------------------------------------------------------------------------------------------------------------//

    void GuidedRepair(Chromosome& chromosome, size_t gene_idx) const
    {
        if (gene_idx >= chromosome.genes.size() || gene_idx >= pair_groups_.size()) {
            return;
        }

        const vector<size_t>& group = pair_groups_[gene_idx];
        int num_options = 1 + static_cast<int>(group.size());
        
        if (num_options <= 1) {
            return;
        }

        int current_choice = chromosome.genes[gene_idx];
        vector<int> alternatives;
        alternatives.reserve(num_options - 1);

        for (int choice = 0; choice < num_options; ++choice) {
            if (choice == current_choice) {
                continue;
            }
            alternatives.push_back(choice);
        }

        if (alternatives.empty()) {
            return;
        }

        random_shuffle(alternatives.begin(), alternatives.end());
        int trials = min(kGuidedRepairTrials, static_cast<int>(alternatives.size()));

        if (trials <= 0) {
            return;
        }

        double best_fitness = EvaluateFitness(chromosome);
        int best_choice = current_choice;

        for (int trial = 0; trial < trials; ++trial) {
            chromosome.genes[gene_idx] = alternatives[trial];
            double candidate_fitness = EvaluateFitness(chromosome);

            if (candidate_fitness > best_fitness) {
                best_fitness = candidate_fitness;
                best_choice = alternatives[trial];
            }
        }

        chromosome.genes[gene_idx] = best_choice;
        chromosome.fitness = best_fitness;
    }

    //-----------------------------------------------------------------------------------------------------------------//

    void BuildOutputsFromSelection(const vector<int>& genes)
    {
        graph_ = MatrixXd::Zero(num_shards_, num_shards_);
        
        struct AdjEdge {
            int to;
            Matrix4d T_to_current;
            double weight; 
        };

        struct FrontierEdge {
            double weight;
            int from;
            int to;
            Matrix4d T_to_from;
        };

        struct FrontierEdgeCompare {
            bool operator()(const FrontierEdge& lhs, const FrontierEdge& rhs) const
            {
                return lhs.weight < rhs.weight;
            }
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

            adjacency[x].push_back({ y, T_yx, edge_weight });
            adjacency[y].push_back({ x, T_xy, edge_weight });
        }

        vector<bool> visited(num_shards_, false);
        vector<bool> is_component_root(num_shards_, false);
        vector<Matrix4d> T_to_root(num_shards_, Matrix4d::Identity());

        // Build an initial pose per connected component using highest-confidence frontier expansion.
        for (int root = 0; root < num_shards_; ++root) {
            if (!IsShardValidAndOn(root) || visited[root]) {
                continue;
            }

            visited[root] = true;
            is_component_root[root] = true;

            priority_queue<FrontierEdge, vector<FrontierEdge>, FrontierEdgeCompare> frontier;
            for (size_t i = 0; i < adjacency[root].size(); ++i) {
                const AdjEdge& edge = adjacency[root][i];
                if (!visited[edge.to]) {
                    frontier.push({ edge.weight, root, edge.to, edge.T_to_current });
                }
            }

            while (!frontier.empty()) {
                FrontierEdge top = frontier.top();
                frontier.pop();

                if (visited[top.to]) {
                    continue;
                }

                visited[top.to] = true;
                T_to_root[top.to] = T_to_root[top.from] * top.T_to_from;

                for (size_t i = 0; i < adjacency[top.to].size(); ++i) {
                    const AdjEdge& edge = adjacency[top.to][i];
                    if (!visited[edge.to]) {
                        frontier.push({ edge.weight, top.to, edge.to, edge.T_to_current });
                    }
                }
            }
        }

        // Refine component poses by averaging both translation and rotation predictions.
        for (int iter = 0; iter < kPoseRelaxIterations; ++iter) {
            vector<Matrix4d> updated = T_to_root;
            bool changed = false;

            for (int node = 0; node < num_shards_; ++node) {
                if (!visited[node] || !IsShardValidAndOn(node) || is_component_root[node]) {
                    continue;
                }

                Vector3d t_acc = Vector3d::Zero();
                Vector4d q_acc = Vector4d::Zero();
                double w_sum = 0.0;
                bool have_rotation = false;
                Quaterniond q_ref(1.0, 0.0, 0.0, 0.0);

                for (int from = 0; from < num_shards_; ++from) {
                    if (!visited[from] || !IsShardValidAndOn(from)) {
                        continue;
                    }

                    for (size_t k = 0; k < adjacency[from].size(); ++k) {
                        if (adjacency[from][k].to != node) {
                            continue;
                        }

                        Matrix4d pred = T_to_root[from] * adjacency[from][k].T_to_current;
                        Vector3d t_pred = pred.block<3, 1>(0, 3);
                        Matrix3d R_pred = pred.block<3, 3>(0, 0);
                        Quaterniond q_pred(R_pred);

                        if (!have_rotation) {
                            q_ref = q_pred;
                            have_rotation = true;
                        }

                        if (q_ref.dot(q_pred) < 0.0) {
                            q_pred.coeffs() *= -1.0;
                        }

                        double w = adjacency[from][k].weight;
                        t_acc += w * t_pred;
                        q_acc += w * q_pred.coeffs();
                        w_sum += w;
                    }
                }

                if (!have_rotation || w_sum <= 0.0) {
                    continue;
                }

                Vector3d t_old = T_to_root[node].block<3, 1>(0, 3);
                Matrix3d R_old = T_to_root[node].block<3, 3>(0, 0);
                Quaterniond q_old(R_old);

                Vector3d t_est = t_acc / w_sum;
                Quaterniond q_est;
                q_est.coeffs() = q_acc / w_sum;
                if (q_est.norm() <= 1.0e-12) {
                    continue;
                }
                q_est.normalize();

                if (q_old.dot(q_est) < 0.0) {
                    q_est.coeffs() *= -1.0;
                }

                Quaterniond q_new = q_old.slerp(kPoseRelaxAlpha, q_est);
                q_new.normalize();
                Vector3d t_new = ((1.0 - kPoseRelaxAlpha) * t_old) + (kPoseRelaxAlpha * t_est);

                if ((t_new - t_old).norm() > 1.0e-6 || q_old.angularDistance(q_new) > 1.0e-6) {
                    changed = true;
                }

                updated[node].setIdentity();
                updated[node].block<3, 3>(0, 0) = q_new.toRotationMatrix();
                updated[node].block<3, 1>(0, 3) = t_new;
            }

            T_to_root = updated;
            if (!changed) {
                break;
            }
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

        Matrix3d I = Matrix3d::Identity();
        Vector3d zero = Vector3d::Zero();
        for (int i = 0; i < num_shards_; ++i) {
            transforms_[i].Set(I, zero, i + 1, 1);
            if (!IsShardValidAndOn(i)) {
                continue;
            }
            transforms_[i].Set(T_to_root[i], i + 1, 1);
        }   
    }

    //-----------------------------------------------------------------------------------------------------------------//

    bool IsShardValidAndOn(int shard_idx) const
    {
        return (shard_idx >= 0) && (shard_idx < num_shards_) && shard_on_off[shard_idx];
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

        for (size_t group_idx = 0; group_idx < pair_groups_.size(); ++group_idx) {
            const vector<size_t>& group = pair_groups_[group_idx];
            if (group.empty()) {
                continue;
            }

            const LCSIndex& representative = matches_[group[0]];
            group_rep_x_[group_idx] = representative.shard_x_ - 1;
            group_rep_y_[group_idx] = representative.shard_y_ - 1;

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
        collision_clouds_.assign(num_shards_, vector<Vector3d>());
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

        const vector<Vector3d>* query = &cloud_a;
        const vector<Vector3d>* target = &cloud_b;
        if (query->size() > target->size()) {
            query = &cloud_b;
            target = &cloud_a;
        }
        diag.query_size = static_cast<int>(query->size());
        diag.target_size = static_cast<int>(target->size());

        const double epsilon_sq = kCollisionPointEpsilon * kCollisionPointEpsilon;
        double depth_sum = 0.0;
        double min_best_sq = numeric_limits<double>::max();

        for (size_t i = 0; i < query->size(); ++i) {
            const Vector3d& q = (*query)[i];
            double best_sq = numeric_limits<double>::max();

            for (size_t j = 0; j < target->size(); ++j) {
                double d_sq = (q - (*target)[j]).squaredNorm();
                if (d_sq < best_sq) {
                    best_sq = d_sq;
                }
            }

            if (best_sq < min_best_sq) {
                min_best_sq = best_sq;
            }

            if (best_sq < epsilon_sq) {
                double d = sqrt(best_sq);
                double depth_ratio = 1.0 - (d / kCollisionPointEpsilon);
                depth_sum += max(0.0, depth_ratio);
                diag.collision_hits++;
            }
        }

        if (min_best_sq < numeric_limits<double>::max()) {
            diag.min_nn_distance = sqrt(min_best_sq);
        }

        if (diag.collision_hits <= 0) {
            return diag;
        }

        diag.hit_ratio = static_cast<double>(diag.collision_hits) / static_cast<double>(max(1, diag.query_size));
        diag.avg_depth = depth_sum / static_cast<double>(diag.collision_hits);
        diag.pair_penalty = kOverlapPenalty * diag.hit_ratio * diag.avg_depth;
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

        vector<vector<Vector3d>> global_clouds(num_shards_);
        vector<Vector3d> global_centers(num_shards_, Vector3d::Zero());
        for (int i = 0; i < num_shards_; ++i) {
            if (!visited[i] || !IsShardValidAndOn(i)) {
                continue;
            }

            const Matrix4d& T = T_to_root[i];
            const vector<Vector3d>& cloud_local = collision_clouds_[i];
            vector<Vector3d>& cloud_global = global_clouds[i];
            cloud_global.reserve(cloud_local.size());

            for (size_t p = 0; p < cloud_local.size(); ++p) {
                Vector4d h;
                h << cloud_local[p], 1.0;
                cloud_global.push_back((T * h).head<3>());
            }

            Vector4d c_h;
            c_h << collision_cloud_centroids_[i], 1.0;
            global_centers[i] = (T * c_h).head<3>();
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

                same_component_pair_count++;
                const vector<Vector3d>& cloud_i = global_clouds[i];
                const vector<Vector3d>& cloud_j = global_clouds[j];
                diag.cloud_a_empty = cloud_i.empty();
                diag.cloud_b_empty = cloud_j.empty();
                if (diag.cloud_a_empty || diag.cloud_b_empty) {
                    pair_diags.push_back(diag);
                    continue;
                }

                diag.center_distance = (global_centers[i] - global_centers[j]).norm();
                diag.broad_radius = collision_cloud_radii_[i] + collision_cloud_radii_[j] + kCollisionPointEpsilon;
                if (diag.center_distance * diag.center_distance > diag.broad_radius * diag.broad_radius) {
                    diag.broad_phase_rejected = true;
                    broad_phase_reject_count++;
                    pair_diags.push_back(diag);
                    continue;
                }

                narrow_phase_pair_count++;
                PairOverlapDiagnostic narrow = ComputeCloudOverlapDiagnostic(cloud_i, cloud_j, i, j);
                diag.query_size = narrow.query_size;
                diag.target_size = narrow.target_size;
                diag.collision_hits = narrow.collision_hits;
                diag.hit_ratio = narrow.hit_ratio;
                diag.avg_depth = narrow.avg_depth;
                diag.min_nn_distance = narrow.min_nn_distance;
                diag.pair_penalty = narrow.pair_penalty;

                if (diag.pair_penalty > 1.0e-12) {
                    active_collision_pair_count++;
                    overlap_sum += diag.pair_penalty;
                }
                pair_diags.push_back(diag);
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

    void RepairPositionCollisions(Chromosome& chromosome, bool verbose = false) const
    {
        for (int iter = 0; iter < kCollisionRepairMaxIter; ++iter) {
            bool collision_resolved = false;

            // Step 1: Build Adjacency List for Pose Propagation
            struct AdjEdge {
                int to;
                Matrix4d T_to_current;
                double density;
                int gene_idx;
            };

            vector<vector<AdjEdge>> adjacency(num_shards_);
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

                double density = static_cast<double>(lcs.inliner_) / (1.0 + lcs.score_);
                Matrix4d T_xy = Matrix4d::Identity();
                lcs.trans_.Output(T_xy);
                Matrix4d T_yx = T_xy.inverse();

                adjacency[x].push_back({ y, T_yx, density, static_cast<int>(group_idx) });
                adjacency[y].push_back({ x, T_xy, density, static_cast<int>(group_idx) });
            }

            // BFS State
            vector<bool> visited(num_shards_, false);
            vector<Matrix4d> T_to_root(num_shards_, Matrix4d::Identity());
            vector<int> component_id(num_shards_, -1);
            vector<int> incoming_gene_idx(num_shards_, -1);
            vector<double> incoming_density(num_shards_, 0.0);
            int current_comp_id = 0;

            for (int start_node = 0; start_node < num_shards_; ++start_node) {
                if (!IsShardValidAndOn(start_node) || visited[start_node]) continue;

                queue<int> q;
                q.push(start_node);
                visited[start_node] = true;
                component_id[start_node] = current_comp_id;
                incoming_gene_idx[start_node] = -1; // Root
                incoming_density[start_node] = 1e9; // Sentinel high density for root

                while (!q.empty()) {
                    int curr = q.front();
                    q.pop();

                    for (const auto& edge : adjacency[curr]) {
                        if (visited[edge.to]) continue;

                        T_to_root[edge.to] = T_to_root[curr] * edge.T_to_current;
                        visited[edge.to] = true;
                        component_id[edge.to] = current_comp_id;
                        incoming_gene_idx[edge.to] = edge.gene_idx;
                        incoming_density[edge.to] = edge.density;
                        q.push(edge.to);
                    }
                }
                current_comp_id++;
            }

            // Compute Global Centroids
            vector<Vector3d> global_centroids(num_shards_, Vector3d::Zero());
            for (int i = 0; i < num_shards_; ++i) {
                if (visited[i]) {
                    Vector4d c_h;
                    c_h << shard_centroids_[i], 1.0;
                    global_centroids[i] = (T_to_root[i] * c_h).head<3>();
                }
            }

            // Step 2: Detect and Step 3: Resolve Collisions
            bool local_collision_found = false;
            for (int i = 0; i < num_shards_; ++i) {
                if (!visited[i]) continue;
                for (int j = i + 1; j < num_shards_; ++j) {
                    if (!visited[j] || component_id[i] != component_id[j]) continue;

                    // Condition 1: Centroid Distance
                    double dist = (global_centroids[i] - global_centroids[j]).norm();
                    double threshold = kCollisionCentroidScale * (shard_radius_[i] + shard_radius_[j]);
                    if (dist >= threshold) continue;

                    // Condition 2: Rotation Similarity
                    // Rotation check removed to sync with fitness evaluation. 
                    // Any spatial overlap is now considered a collision.

                    // Collision detected!
                    local_collision_found = true;

                    int loser_idx = -1;
                    if (incoming_gene_idx[i] == -1) {
                        loser_idx = j; // i is root, j must lose
                    } else if (incoming_gene_idx[j] == -1) {
                        loser_idx = i; // j is root, i must lose
                    } else {
                        // 50/50 chance to pick loser to allow breaking out of local optima
                        loser_idx = (rand() % 2 == 0) ? i : j;
                    }

                    if (loser_idx != -1 && incoming_gene_idx[loser_idx] != -1) {
                        int bad_gene_idx = incoming_gene_idx[loser_idx];
                        int current_choice = chromosome.genes[bad_gene_idx];
                        
                        const vector<size_t>& group = pair_groups_[bad_gene_idx];
                        vector<int> candidates;
                        vector<double> weights;
                        double total_weight = 0.0;

                        // Collect all OTHER matches in this group and their densities
                        for (size_t c_idx = 0; c_idx < group.size(); ++c_idx) {
                            int candidate = static_cast<int>(c_idx) + 1;
                            if (candidate == current_choice) continue;

                            const LCSIndex& lcs = matches_[group[c_idx]];
                            double w = static_cast<double>(lcs.inliner_) / (1.0 + lcs.score_);
                            if (w > 0) {
                                candidates.push_back(candidate);
                                weights.push_back(w);
                                total_weight += w;
                            }
                        }

                        // Perform weighted random selection, fallback to 0 on last attempt
                        bool force_zero = (iter == kCollisionRepairMaxIter - 1);
                        if (!force_zero && total_weight > 0) {
                            double r = static_cast<double>(rand()) / RAND_MAX * total_weight;
                            double cumulative = 0.0;
                            int selected = 0;
                            for (size_t k = 0; k < candidates.size(); ++k) {
                                cumulative += weights[k];
                                if (r <= cumulative) {
                                    selected = candidates[k];
                                    break;
                                }
                            }
                            chromosome.genes[bad_gene_idx] = selected;
                            
                            if (verbose) {
                                stringstream ss;
                                ss << "[REPAIR] Shard " << loser_idx << " resampled gene " << bad_gene_idx 
                                << " from choice " << current_choice << " to " << selected;
                                LogDiagnostic(ss.str());
                            }
                        } else {
                            chromosome.genes[bad_gene_idx] = 0; 
                            if (verbose) {
                                string reason = (force_zero) ? " (Hard Fallback)" : " (No alternatives)";
                                LogDiagnostic("[REPAIR] Shard " + to_string(loser_idx) + " gene " + to_string(bad_gene_idx) + " set to 0" + reason);
                            }
                        }
                        collision_resolved = true;
                    }

                    // To reflect changes (reachability/poses), break and re-run BFS
                    goto next_repair_iter;
                }
            }

            if (!local_collision_found) break;

        next_repair_iter:
            if (!collision_resolved) break;
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
        ss << "  Recomputed cloud-overlap sum: " << fixed << setprecision(4) << overlap_sum
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
    }

    void LogBestChromosomeBreakdown(int generation, const Chromosome& best) const
    {
        FitnessBreakdown breakdown;
        EvaluateFitness(best, &breakdown);
        
        stringstream ss;
        ss << "\n--- GENERATION " << generation << " BEST BREAKDOWN ---" << endl;
        ss << "  Total Fitness: " << fixed << setprecision(2) << best.fitness << endl;
        ss << "  Inlier Reward: " << breakdown.inlier_reward << endl;
        ss << "  Overlap Penalty: " << breakdown.overlap_penalty << endl;
        ss << "  Edge Residual Penalty: " << breakdown.edge_residual_penalty << endl;
        ss << "  Edge Rot Residual Penalty: " << breakdown.edge_rot_residual_penalty << endl;
        ss << "  Connectivity Reward: " << breakdown.connectivity_reward << endl;
        ss << "  Component Penalty: " << breakdown.connectivity_component_penalty << endl;
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

        ss << "  OverlapDiag Sum: " << fixed << setprecision(4) << overlap_sum
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
            size_t repair_idx = random_genes[rand() % random_genes.size()];
            // GuidedRepair(replacement, repair_idx);
        }

        EnsureSherdCoverage(replacement);

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

        static thread_local std::mt19937 gen(std::random_device{}());
        std::uniform_real_distribution<double> dist(0.0, 1.0);

        if (dist(gen) < kInitialPairInactiveRate) {
            return 0;
        }

        const vector<size_t>& group = pair_groups_[group_idx];
        if (group.empty()) {
            return 0;
        }

        // --- Density-Based Selection ---
        double total_weight = 0.0;
        vector<double> cumulative;
        cumulative.reserve(group.size());

        for (size_t i = 0; i < group.size(); ++i) {
            const LCSIndex& lcs = matches_[group[i]];
            // density = inliner / (1 + score)
            double density = static_cast<double>(lcs.inliner_) / (1.0 + lcs.score_);
            double w = max(1.0e-6, density);
            total_weight += w;
            cumulative.push_back(total_weight);
        }

        if (total_weight <= 0.0) {
            return 1 + (gen() % static_cast<int>(group.size()));
        }

        double r = dist(gen) * total_weight;
        for (size_t i = 0; i < cumulative.size(); ++i) {
            if (r <= cumulative[i]) {
                return static_cast<int>(i + 1);
            }
        }

        return 1;
    }
    
    void AuditGroundTruthCandidates(const MatrixXd& GT_graph, const vector<Trans>& GT_trans, const vector<Trans>& T_axis)
    {
        cout << "#################### GROUND TRUTH CANDIDATE AUDIT ####################" << endl;
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

            for (size_t i = 0; i < group.size(); ++i) {
                const LCSIndex& lcs = matches_[group[i]];
                int x = lcs.shard_x_ - 1;
                int y = lcs.shard_y_ - 1;

                Matrix4d T_ax, T_ay_inv, T_gx_inv, T_gy;
                T_axis[x].Output(T_ax);
                T_axis[y].InvOut(T_ay_inv);
                GT_trans[x].InvOut(T_gx_inv);
                GT_trans[y].Output(T_gy);

                Matrix4d T_gt_wrapped = T_ax * T_gx_inv * T_gy * T_ay_inv;
                Matrix4d T_gt_raw = T_gx_inv * T_gy;

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

                double r_wrapped, t_wrapped, r_raw, t_raw;
                calc_err(T_gt_wrapped, T_candidate, r_wrapped, t_wrapped);
                calc_err(T_gt_raw, T_candidate, r_raw, t_raw);

                double density = static_cast<double>(lcs.inliner_) / (1.0 + lcs.score_);

                cout << "  Choice " << (i + 1) << ": inliner=" << setw(3) << lcs.inliner_ 
                     << " score=" << fixed << setprecision(3) << lcs.score_ 
                     << " density=" << fixed << setprecision(3) << density
                     << " [Wrapped] rad=" << fixed << setprecision(4) << r_wrapped << " trans=" << t_wrapped;
                
                if (r_wrapped < 0.35 && t_wrapped < 50.0) {
                    cout << " [VALID GT MATCH]";
                }
                cout << endl;

                if (r_wrapped < min_rad_err) {
                    min_rad_err = r_wrapped;
                    min_trans_err = t_wrapped;
                    best_choice = static_cast<int>(i + 1);
                }
            }
            cout << "  -> Best Match for GT: Choice " << best_choice << " (rad_err=" << min_rad_err << ")" << endl;
        }
        cout << "######################################################################" << endl;
    }

    //-----------------------------------------------------------------------------------------------------------------//

private:
    static constexpr int kPopulationSize = 200;
    static constexpr int kMaxGenerations = 100;
    static constexpr int kElitismCount = 10;
    static constexpr double kMutationRate = 0.10; // Increased from 0.05 for better symmetry exploration
    static constexpr double kSymmetryFlipRate = 0.15;
    static constexpr double kBiasInheritRatio = 0.4;
    static constexpr int kGuidedRepairTrials = 3;
    static constexpr double kInitialPairInactiveRate = 0.05;
    static constexpr double kGreedyPairChoiceRate = 0.0;
    static constexpr double kPairChoicePenaltyWeight = 35.0;
    static constexpr double kPairChoiceDensityEps = 1.0e-6;
    static constexpr double kMinActivePairRatio = 0.65;
    static constexpr double kMaxActivePairRatio = 0.90;
    static constexpr double kActivePairRangePenaltyWeight = 8.0;
    static constexpr double kExcessActivePairPenaltyScale = 0.5;
    static constexpr double kEdgeResidualThreshold = 80.0;
    static constexpr double kEdgeResidualPenalty = 2.0;
    static constexpr double kEdgeRotResidualThreshold = 1.2;
    static constexpr double kEdgeRotResidualPenalty = 5.0;
    static constexpr double kConnectivityReward = 100.0;
    static constexpr double kConnectivityComponentPenalty = 100;
    static constexpr int kPoseRelaxIterations = 0;
    static constexpr double kPoseRelaxAlpha = 0.5;
    
    // Overlap constants
    static constexpr double kOverlapPenalty = 1500.0;
    static constexpr double kCollisionPointEpsilon = 2.0;   // mm
    static constexpr double kCollisionEdgeExclusion = 2.0;  // mm
    static constexpr double kCollisionVoxelSize = 5.0;      // mm
    static constexpr int kCollisionCloudMaxPoints = 64;
    static constexpr double kCollisionCentroidScale = 0.70;     // fraction of combined radii for centroid threshold
    static constexpr double kCollisionRotAngleThreshold = 0.53; // radians (~90 degrees)
    static constexpr int kCollisionRepairMaxIter = 20;           // max repair iterations per chromosome

    vector<Geom> shard_;
    vector<LCSIndex> matches_;
    vector<vector<size_t>> pair_groups_;
    vector<vector<size_t>> sherd_incident_groups_;
    int num_shards_;

    vector<Chromosome> population_;

    vector<Trans> transforms_;
    MatrixXd graph_;

    // Adaptive Geometry for overlap checks
    vector<Vector3d> shard_centroids_;
    vector<double> shard_radius_;
    vector<vector<Vector3d>> collision_clouds_;
    vector<CollisionCloudStats> collision_cloud_stats_;
    vector<Vector3d> collision_cloud_centroids_;
    vector<double> collision_cloud_radii_;

    // Static caches to avoid recomputing identical data in every fitness call.
    int max_neighbors_cap_ = 0;
    vector<double> group_best_density_;
    vector<int> group_rep_x_;
    vector<int> group_rep_y_;
};

#endif
