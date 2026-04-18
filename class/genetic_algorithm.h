#pragma once

// Check # 01: Having both pragma once and this is redundant so the below two lines can be removed apparently
#ifndef _GENETIC_ALGORITHM_H_   
#define _GENETIC_ALGORITHM_H_   

#include <algorithm>
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
    bool use_neighbor_penalty    = true;
    bool use_overlap_penalty     = true;

    struct Chromosome {
        vector<int> genes;
        double fitness;
    };

    //-----------------------------------------------------------------------------------------------------------------//

    GeneticAssembler(const vector<Geom>& shard, list<LCSIndex>& LCS_out, int num_shards) 
        : shard_(shard), num_shards_(num_shards) 
    {   

        for (list<LCSIndex>::const_iterator it = LCS_out.begin(); it != LCS_out.end(); ++it) {

            // Check # 05: We can try to change this threshold of 10.0 everywhere in the file
            if (it->score_ >= 500.0) {
                continue;
            }
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

        if (!matches_.empty()) {
            cout << "[GA DEBUG] First Match: shard_x_ = " << matches_[0].shard_x_ 
            << " shard_y_ = " << matches_[0].shard_y_ << endl;
        }

        cout << "[GA DEBUG] Chromosome Length (Total Connections): " << pair_groups_.size() << endl;

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
                    GuidedRepair(child1, static_cast<size_t>(mutated_gene1));
                }
                next_population.push_back(child1);

                if (static_cast<int>(next_population.size()) < kPopulationSize) {
                    Chromosome child2 = Crossover(parent2, parent1);
                    int mutated_gene2 = Mutate(child2);
                    if (mutated_gene2 >= 0) {
                        GuidedRepair(child2, static_cast<size_t>(mutated_gene2));
                    }

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
        int active_count = 0;
        for (size_t group_idx = 0; group_idx < best.genes.size(); ++group_idx) {
            int choice = best.genes[group_idx];
            if (choice <= 0) {
                continue;
            }

            const vector<size_t>& group = pair_groups_[group_idx];
            int local_idx = choice - 1;
            if (local_idx < 0 || local_idx >= static_cast<int>(group.size())) {
                continue;
            }

            const LCSIndex& lcs = matches_[group[local_idx]];
            active_count++;
            cout << "[GA DEBUG] Active match: shard_x_="
                << lcs.shard_x_
                << " shard_y_=" << lcs.shard_y_
                << " inliner_=" << lcs.inliner_
                << " score_=" << lcs.score_ << endl;
        }

        cout << "[GA DEBUG] Total active matches: " << active_count << endl;

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

    // Check # 13: The fitness function might need a thorough analysis again. Moreover, why are we not using functions in other files (feature_matchings, ranking_system, etc.) at all? We have only used LCS so far?
    double EvaluateFitness(const Chromosome& chromosome) const
    {   
                
        int kMaxNeighbors = 0;
        
        // Check # 07: This logic might be incorrect given that it can count each shard multiple times because shard pairs can repeat in matches_
        for (int i = 0; i < num_shards_; ++i) {
            int count = 0;
            for (size_t j = 0; j < matches_.size(); ++j) {
                if (matches_[j].shard_x_ - 1 == i || matches_[j].shard_y_ - 1 == i) {
                    count++;
                }
            }
            kMaxNeighbors = max(kMaxNeighbors, count);
            // Check # 08: Why is there a need of this max() code above? It can be removed.
        }

        // Cap at num_shards_ - 1 (maximum possible connections)
        kMaxNeighbors = min(kMaxNeighbors, num_shards_ - 1);

        double fitness = 0.0;
        vector<int> neighbor_count(num_shards_, 0);
        vector<vector<int>> adjacency(num_shards_);

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
            if (use_inlier_score) {
                double score_weight = 1.0 / (1.0 + lcs.score_);
                double density_score = static_cast<double>(lcs.inliner_) * score_weight;
                fitness += density_score;
            }

            neighbor_count[x]++;
            neighbor_count[y]++;
            adjacency[x].push_back(y);
            adjacency[y].push_back(x);
        }

        for (int shard_idx = 0; shard_idx < num_shards_; ++shard_idx) {
            if (!IsShardValidAndOn(shard_idx)) {
                continue;
            }
            // Check # 10: This condition might be totally illogical because neighbor_count[shard_idx] can never exceed kMaxNeighbors as per the above implementation. A more useful function migh help here.
            if (neighbor_count[shard_idx] > kMaxNeighbors) {
                if (use_neighbor_penalty) {
                    fitness -= 50.0 * static_cast<double>(neighbor_count[shard_idx] - kMaxNeighbors);
                }
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

            auto key = make_pair(min(x, y), max(x, y));
            if (active_transforms.find(key) == active_transforms.end()) {
                active_transforms[key] = T_xy;
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

                    Matrix4d T_composed = T_ab * T_bc;

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
            fitness -= 2.0 * cycle_penalty;
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

                Matrix4d T_pred_ab = T_pose[edge.a].inverse() * T_pose[edge.b];
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
            fitness -= kEdgeResidualPenalty * edge_residual_penalty;
        }
        if (use_rot_residual) {
            fitness -= kEdgeRotResidualPenalty * edge_rot_residual_penalty;
        }

        // --- Overlap Penalty (Physical occupancy check to handle symmetry) ---
        if (use_overlap_penalty) {
            int start_node = -1;
            for (int i = 0; i < num_shards_; ++i) {
                if (IsShardValidAndOn(i)) { start_node = i; break; }
            }

            if (start_node != -1) {
                vector<bool> vis(num_shards_, false);
                vector<Matrix4d> T_comp(num_shards_, Matrix4d::Identity());
                queue<int> q_comp;
                q_comp.push(start_node);
                vis[start_node] = true;
                vector<int> placed;
                placed.push_back(start_node);

                while (!q_comp.empty()) {
                    int curr = q_comp.front(); q_comp.pop();
                    for (const auto& edge : pose_adj[curr]) {
                        if (!vis[edge.first]) {
                            T_comp[edge.first] = T_comp[curr] * edge.second;
                            vis[edge.first] = true;
                            placed.push_back(edge.first);
                            q_comp.push(edge.first);
                        }
                    }
                }

                int overlap_violations = 0;
                for (size_t i = 0; i < placed.size(); ++i) {
                    for (size_t j = i + 1; j < placed.size(); ++j) {
                        int idx1 = placed[i];
                        int idx2 = placed[j];
                        Vector4d c1_h; c1_h << shard_centroids_[idx1], 1.0;
                        Vector4d c2_h; c2_h << shard_centroids_[idx2], 1.0;
                        Vector3d gc1 = (T_comp[idx1] * c1_h).head<3>();
                        Vector3d gc2 = (T_comp[idx2] * c2_h).head<3>();
                        double dist = (gc1 - gc2).norm();
                        double min_dist = (shard_radius_[idx1] + shard_radius_[idx2]) * kOverlapThresholdScale;
                        if (dist < min_dist) overlap_violations++;
                    }
                }
                fitness -= overlap_violations * kOverlapPenalty;
            }
        }

        int largest_component = 0;
        int num_components = 0;

        AnalyzeConnectedComponents(adjacency, largest_component, num_components);
        if (use_connectivity_reward) {
            fitness += kConnectivityReward * static_cast<double>(largest_component);

            if (num_components > 1) {
                fitness -= kConnectivityComponentPenalty * static_cast<double>(num_components - 1);
            }
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
                // 50% chance to perform a "Symmetry Flip" (pick a non-weighted alternative)
                // 50% chance to follow standard inlier-weighted sampling
                if ((rand() % 100) < 50) {
                    // Forced Symmetry Flip: Pick any other candidate except the current one or 'inactive'
                    // This forces the GA to explore the alternative high-quality matches in the group
                    int num_candidates = static_cast<int>(group.size());
                    if (num_candidates > 1) {
                         int draw = rand() % num_candidates;
                         new_choice = draw + 1; 
                         if (new_choice == old_choice) {
                             new_choice = (draw + 1) % num_candidates + 1;
                         }
                    } else {
                         new_choice = SampleGroupChoice(i);
                    }
                } else {
                    const int kMaxResampleAttempts = 16;
                    for (int attempt = 0; attempt < kMaxResampleAttempts; ++attempt) {
                        new_choice = SampleGroupChoice(i);
                        if (new_choice != old_choice) {
                            break;
                        }
                    }
                }

                if (new_choice == old_choice) {
                    int clamped_old = old_choice;
                    if (clamped_old < 0 || clamped_old >= num_options) {
                        clamped_old = 0;
                    }

                    int draw = rand() % (num_options - 1);
                    if (draw >= clamped_old) {
                        draw++;
                    }
                    new_choice = draw;
                }

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

            if (lcs.score_ >= 100.0) {
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
        vector<Matrix4d> T_to_root(num_shards_, Matrix4d::Identity());

        int bfs_root = -1;
        for (int i = 0; i < num_shards_; i++) {
            if (IsShardValidAndOn(i)) {
                bfs_root = i;
                break;
            }
        }

        if (bfs_root >= 0) {
            queue<int> q;
            q.push(bfs_root);
            visited[bfs_root] = true;

            while (!q.empty()) {
                int current = q.front();
                q.pop();

                for (size_t i = 0; i < adjacency[current].size(); ++i) {
                    int next = adjacency[current][i].to;
                    const Matrix4d& T_next_current = adjacency[current][i].T_to_current; 
                    if (visited[next]) {
                        continue;
                    }

                    T_to_root[next] = T_to_root[current] * T_next_current;
                    visited[next] = true;
                    q.push(next);
                }
            }

            for (int iter = 0; iter < kPoseRelaxIterations; ++iter) {
                vector<Matrix4d> updated = T_to_root;
                bool changed = false;

                for (int node = 0; node < num_shards_; ++node) {
                    if (!visited[node] || node == bfs_root || !IsShardValidAndOn(node)) {
                        continue;
                    }

                    Vector3d t_acc = Vector3d::Zero();
                    double w_sum = 0.0;

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
                            double w = adjacency[from][k].weight;
                            t_acc += w * t_pred;
                            w_sum += w;
                        }
                    }

                    if (w_sum <= 0.0) {
                        continue;
                    }

                    Vector3d t_old = T_to_root[node].block<3, 1>(0, 3);
                    Vector3d t_est = t_acc / w_sum;
                    Vector3d t_new = ((1.0 - kPoseRelaxAlpha) * t_old) +  (kPoseRelaxAlpha * t_est);

                    if ((t_new - t_old).norm() > 1.0e-6) {
                        changed = true;
                    }

                    updated[node].block<3, 1>(0, 3) = t_new;
                }

                T_to_root = updated;
                if (!changed) {
                    break;
                }
            }
        }

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
            GuidedRepair(replacement, repair_idx);
        }

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

        double activate_prob = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
        if (activate_prob >= kInitialPairActivationRate) {
            return 0;
        }

        const vector<size_t>& group = pair_groups_[group_idx];
        if (group.empty()) {
            return 0;
        }

        double total_weight = 0.0;
        vector<double> cumulative;
        cumulative.reserve(group.size());

        for (size_t i = 0; i < group.size(); ++i) {
            const LCSIndex& lcs = matches_[group[i]];
            double w = max(1.0, static_cast<double>(lcs.inliner_));
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
    static constexpr int kElitismCount = 2;
    static constexpr double kMutationRate = 0.08; // Increased from 0.05 for better symmetry exploration
    static constexpr double kBiasInheritRatio = 0.4;
    static constexpr int kGuidedRepairTrials = 3;
    static constexpr double kInitialPairActivationRate = 0.3;
    static constexpr double kEdgeResidualThreshold = 50.0;
    static constexpr double kEdgeResidualPenalty = 0.05;
    static constexpr double kEdgeRotResidualThreshold = 0.35;
    static constexpr double kEdgeRotResidualPenalty = 1.0;
    static constexpr double kConnectivityReward = 10.0;
    static constexpr double kConnectivityComponentPenalty = 8.0;
    static constexpr int kPoseRelaxIterations = 5;
    static constexpr double kPoseRelaxAlpha = 0.5;
    
    // Overlap constants
    static constexpr double kOverlapPenalty = 1000.0;
    static constexpr double kOverlapThresholdScale = 0.4; // % of combined radii to trigger overlap penalty

    vector<Geom> shard_;
    vector<LCSIndex> matches_;
    vector<vector<size_t>> pair_groups_;
    int num_shards_;

    vector<Chromosome> population_;

    vector<Trans> transforms_;
    MatrixXd graph_;

    // Adaptive Geometry for overlap checks
    vector<Vector3d> shard_centroids_;
    vector<double> shard_radius_;
};

#endif

