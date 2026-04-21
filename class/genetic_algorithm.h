#pragma once

// Check # 01: Having both pragma once and this is redundant so the below two lines can be removed apparently
#ifndef _GENETIC_ALGORITHM_H_   
#define _GENETIC_ALGORITHM_H_   

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <map>
#include <queue>
#include <vector>

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
    bool use_overlap_penalty     = false;
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
    }

    //-----------------------------------------------------------------------------------------------------------------//

    void Run()
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
                RepairPositionCollisions(child1);
                next_population.push_back(child1);

                if (static_cast<int>(next_population.size()) < kPopulationSize) {
                    Chromosome child2 = Crossover(parent2, parent1);
                    int mutated_gene2 = Mutate(child2);
                    if (mutated_gene2 >= 0) {
                        // GuidedRepair(child2, static_cast<size_t>(mutated_gene2));
                    }

                    EnsureSherdCoverage(child2);
                    RepairPositionCollisions(child2);

                    next_population.push_back(child2);
                }
            }

            population_ = next_population;

            EvaluatePopulation();
            sort(population_.begin(), population_.end(), [](const Chromosome& a, const Chromosome& b) {
                return a.fitness > b.fitness;
            });

            EnforceDiversity();

            sort(population_.begin(), population_.end(), [](const Chromosome& a, const Chromosome& b) {
                return a.fitness > b.fitness;
            });
        }

        EvaluatePopulation();
        sort(population_.begin(), population_.end(), [](const Chromosome& a, const Chromosome& b) {
            return a.fitness > b.fitness;
        });

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
            RepairPositionCollisions(chromosome);
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

        double edge_residual_penalty = 0.0;
        double edge_rot_residual_penalty = 0.0;
        int evaluated_roots = 0;
        
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

            bool root_used = false;

            for (const auto& kv : pair_edges) {
                const PairEdge& edge = kv.second;
                if (!pose_visited[edge.a] || !pose_visited[edge.b]) {
                    continue;
                }

                root_used = true;

                Matrix4d T_pred_ab = T_pose[edge.b].inverse() * T_pose[edge.a];
                Vector3d t_pred = T_pred_ab.block<3, 1>(0, 3);
                Vector3d t_ab = edge.T_ab.block<3, 1>(0, 3);
                double residual = (t_pred - t_ab).norm();

                if (residual > kEdgeResidualThreshold) {
                    double excess = residual - kEdgeResidualThreshold;
                    edge_residual_penalty += excess * excess;
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
                    edge_rot_residual_penalty += (rot_residual - kEdgeRotResidualThreshold);
                }
            }

            if (root_used) {
                evaluated_roots++;
            }
        }

        if (evaluated_roots > 0) {
            edge_residual_penalty /= static_cast<double>(evaluated_roots);
            edge_rot_residual_penalty /= static_cast<double>(evaluated_roots);
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
            int overlap_violations = 0;

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

                for (size_t i = 0; i < placed.size(); ++i) {
                    for (size_t j = i + 1; j < placed.size(); ++j) {
                        int idx1 = placed[i];
                        int idx2 = placed[j];
                        Vector4d c1_h;
                        c1_h << shard_centroids_[idx1], 1.0;
                        Vector4d c2_h;
                        c2_h << shard_centroids_[idx2], 1.0;
                        Vector3d gc1 = (T_comp[idx1] * c1_h).head<3>();
                        Vector3d gc2 = (T_comp[idx2] * c2_h).head<3>();
                        double dist = (gc1 - gc2).norm();
                        double min_dist = (shard_radius_[idx1] + shard_radius_[idx2]) * kOverlapThresholdScale;
                        if (dist < min_dist) {
                            overlap_violations++;
                        }
                    }
                }
            }

            double penalty = overlap_violations * kOverlapPenalty;
            fitness -= penalty;
            if (breakdown != nullptr) {
                breakdown->overlap_penalty = penalty;
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
                        while (new_choice != old_choice) {
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

    void RepairPositionCollisions(Chromosome& chromosome) const
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
                    double threshold = kCollisionCentroidScale * min(shard_radius_[i], shard_radius_[j]);
                    if (dist >= threshold) continue;

                    // Condition 2: Rotation Similarity
                    Matrix3d Ri = T_to_root[i].block<3, 3>(0, 0);
                    Matrix3d Rj = T_to_root[j].block<3, 3>(0, 0);
                    double tr = (Ri.transpose() * Rj).trace();
                    double angle = acos(max(-1.0, min(1.0, (tr - 1.0) / 2.0)));

                    if (angle >= kCollisionRotAngleThreshold) continue;

                    // Collision detected!
                    local_collision_found = true;

                    int loser_idx = -1;
                    if (incoming_gene_idx[i] == -1) {
                        loser_idx = j; // i is root, j must lose
                    } else if (incoming_gene_idx[j] == -1) {
                        loser_idx = i; // j is root, i must lose
                    } else if (incoming_density[i] < incoming_density[j]) {
                        loser_idx = i;
                    } else {
                        loser_idx = j;
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

                        // Perform weighted random selection
                        if (total_weight > 0) {
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
                        } else {
                            chromosome.genes[bad_gene_idx] = 0; // No other valid matches for this pair
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

        double inactive_prob = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
        if (inactive_prob < kInitialPairInactiveRate) {
            return 0;
        }

        const vector<size_t>& group = pair_groups_[group_idx];
        if (group.empty()) {
            return 0;
        }

        double greedy_pick_prob = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
        if (greedy_pick_prob < kGreedyPairChoiceRate) {
            int best_idx = -1;
            double best_density = -1.0;
            for (size_t i = 0; i < group.size(); ++i) {
                const LCSIndex& lcs = matches_[group[i]];
                double density = static_cast<double>(lcs.inliner_) / (1.0 + lcs.score_);
                if (density > best_density) {
                    best_density = density;
                    best_idx = static_cast<int>(i);
                }
            }

            if (best_idx >= 0) {
                return best_idx + 1;
            }
        }

        double total_weight = 0.0;
        vector<double> cumulative;
        cumulative.reserve(group.size());

        for (size_t i = 0; i < group.size(); ++i) {
            const LCSIndex& lcs = matches_[group[i]];
            double density = static_cast<double>(lcs.inliner_) / (1.0 + lcs.score_);
            double w = max(1.0e-6, density);
            total_weight += w;
            cumulative.push_back(total_weight);
        }

        if (total_weight <= 0.0) {
            return 1 + (rand() % static_cast<int>(group.size()));
        }

        double r = (static_cast<double>(rand()) / static_cast<double>(RAND_MAX)) * total_weight;
        for (size_t i = 0; i < cumulative.size(); ++i) {
            if (r <= cumulative[i]) {
                return static_cast<int>(i) + 1;
            }
        }

        return static_cast<int>(group.size());
    }

    //-----------------------------------------------------------------------------------------------------------------//

private:
    static constexpr int kPopulationSize = 100;
    static constexpr int kMaxGenerations = 100;
    static constexpr int kElitismCount = 10;
    static constexpr double kMutationRate = 0.10; // Increased from 0.05 for better symmetry exploration
    static constexpr double kSymmetryFlipRate = 0.15;
    static constexpr double kBiasInheritRatio = 0.4;
    static constexpr int kGuidedRepairTrials = 3;
    static constexpr double kInitialPairInactiveRate = 0.05;
    static constexpr double kGreedyPairChoiceRate = 0.5;
    static constexpr double kPairChoicePenaltyWeight = 35.0;
    static constexpr double kPairChoiceDensityEps = 1.0e-6;
    static constexpr double kMinActivePairRatio = 0.65;
    static constexpr double kMaxActivePairRatio = 0.90;
    static constexpr double kActivePairRangePenaltyWeight = 8.0;
    static constexpr double kExcessActivePairPenaltyScale = 0.5;
    static constexpr double kEdgeResidualThreshold = 50.0;
    static constexpr double kEdgeResidualPenalty = 1.0;
    static constexpr double kEdgeRotResidualThreshold = 0.35;
    static constexpr double kEdgeRotResidualPenalty = 1.0;
    static constexpr double kConnectivityReward = 100.0;
    static constexpr double kConnectivityComponentPenalty = 100;
    static constexpr int kPoseRelaxIterations = 0;
    static constexpr double kPoseRelaxAlpha = 0.5;
    
    // Overlap constants
    static constexpr double kOverlapPenalty = 1000.0;
    static constexpr double kOverlapThresholdScale = 0.5; // % of combined radii to trigger overlap penalty
    static constexpr double kCollisionCentroidScale = 0.5;     // fraction of min radius for centroid threshold
    static constexpr double kCollisionRotAngleThreshold = 1.57; // radians (~30 degrees)
    static constexpr int kCollisionRepairMaxIter = 5;           // max repair iterations per chromosome

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

    // Static caches to avoid recomputing identical data in every fitness call.
    int max_neighbors_cap_ = 0;
    vector<double> group_best_density_;
    vector<int> group_rep_x_;
    vector<int> group_rep_y_;
};

#endif

