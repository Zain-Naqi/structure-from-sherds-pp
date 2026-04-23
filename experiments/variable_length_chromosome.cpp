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
    struct Gene {
        size_t group_idx;
        int choice; // 1-indexed choice within the group
    };

    struct Chromosome {
        vector<Gene> genes;
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
        // Check # 04: The referenced relative to shard index is 1 here, which might need to be changed based on feedback from Check # 03
    }

    //-----------------------------------------------------------------------------------------------------------------//

    void Run()
    {
        srand(42);

        if (pair_groups_.empty()) {
            cout << "PAIR GROUPS ARE EMPTY...WHY?" << endl;
            BuildOutputsFromSelection(Chromosome());
            return;
        }

        if (!matches_.empty()) {
            cout << "[GA DEBUG] First Match: shard_x_ = " << matches_[0].shard_x_ 
            << " shard_y_ = " << matches_[0].shard_y_ << endl;
        }

        InitializePopulation();

        for (int generation = 0; generation < kMaxGenerations; ++generation) {
            EvaluatePopulation();
            sort(population_.begin(), population_.end(), [](const Chromosome& a, const Chromosome& b) {
                return a.fitness > b.fitness;
            });

            cout << "[GA] Generation " << generation 
                 << " best fitness: " << population_.front().fitness 
                 << " length: " << population_.front().genes.size() 
                 << " / " << pair_groups_.size() << endl;

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
                Mutate(child1);
                next_population.push_back(child1);

                if (static_cast<int>(next_population.size()) < kPopulationSize) {
                    Chromosome child2 = Crossover(parent2, parent1);
                    Mutate(child2);
                    next_population.push_back(child2);
                }
            }

            population_ = next_population;

        }

        EvaluatePopulation();
        sort(population_.begin(), population_.end(), [](const Chromosome& a, const Chromosome& b) {
            return a.fitness > b.fitness;
        });

        const Chromosome& best = population_.front();
        int active_count = static_cast<int>(best.genes.size());
        for (const auto& gene : best.genes) {
            const vector<size_t>& group = pair_groups_[gene.group_idx];
            const LCSIndex& lcs = matches_[group[gene.choice - 1]];
            cout << "[GA DEBUG] Active match: shard_x_="
                << lcs.shard_x_
                << " shard_y_=" << lcs.shard_y_
                << " inliner_=" << lcs.inliner_
                << " score_=" << lcs.score_ << endl;
        }

        cout << "[GA DEBUG] Total active matches: " << active_count << endl;

        BuildOutputsFromSelection(population_.front());
        
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
            chromosome.fitness = 0.0;

            for (size_t gene_idx = 0; gene_idx < pair_groups_.size(); ++gene_idx) {
                // Biased activation for each pair group
                double r = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
                if (r < kInitialPairActivationRate) {
                    int choice = SelectWeightedCandidate(gene_idx);
                    if (choice > 0) {
                        chromosome.genes.push_back({gene_idx, choice});
                    }
                }
            }
            Repair(chromosome);
            population_.push_back(chromosome);
        }
    }

    //-----------------------------------------------------------------------------------------------------------------//

    void Repair(Chromosome& chromosome) const
    {
        if (chromosome.genes.empty()) return;

        // 1. Remove out-of-bounds gene indices
        vector<Gene> valid_genes;
        for (const auto& g : chromosome.genes) {
            if (g.group_idx < pair_groups_.size()) {
                const vector<size_t>& group = pair_groups_[g.group_idx];
                if (g.choice >= 1 && g.choice <= static_cast<int>(group.size())) {
                    valid_genes.push_back(g);
                }
            }
        }

        // 2. Resolve duplicate group indices (keep the one with best inliner count)
        map<size_t, Gene> group_map;
        for (const auto& g : valid_genes) {
            auto it = group_map.find(g.group_idx);
            if (it == group_map.end()) {
                group_map[g.group_idx] = g;
            } else {
                // Compare inliner counts
                const LCSIndex& current_lcs = matches_[pair_groups_[g.group_idx][g.choice - 1]];
                const LCSIndex& existing_lcs = matches_[pair_groups_[it->second.group_idx][it->second.choice - 1]];
                if (current_lcs.inliner_ > existing_lcs.inliner_) {
                    it->second = g;
                }
            }
        }

        chromosome.genes.clear();
        for (auto const& [idx, gene] : group_map) {
            chromosome.genes.push_back(gene);
        }
    }

    //-----------------------------------------------------------------------------------------------------------------//

    int SelectWeightedCandidate(size_t group_idx) const
    {
        const vector<size_t>& group = pair_groups_[group_idx];
        if (group.empty()) return 0;

        double total_weight = 0.0;
        vector<double> cumulative;
        cumulative.reserve(group.size());

        for (size_t i = 0; i < group.size(); ++i) {
            const LCSIndex& lcs = matches_[group[i]];
            double w = max(1.0, static_cast<double>(lcs.inliner_));
            total_weight += w;
            cumulative.push_back(total_weight);
        }

        if (total_weight <= 0.0) return 1 + (rand() % static_cast<int>(group.size()));

        double r = (static_cast<double>(rand()) / static_cast<double>(RAND_MAX)) * total_weight;
        for (size_t i = 0; i < cumulative.size(); ++i) {
            if (r <= cumulative[i]) return static_cast<int>(i) + 1;
        }
        return static_cast<int>(group.size());
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

        for (const auto& gene : chromosome.genes) {
            const vector<size_t>& group = pair_groups_[gene.group_idx];
            const LCSIndex& lcs = matches_[group[gene.choice - 1]];
            int x = lcs.shard_x_ - 1;
            int y = lcs.shard_y_ - 1;
            if (!IsShardValidAndOn(x) || !IsShardValidAndOn(y)) {
                continue;
            }
            
            double score_weight = 1.0 / (1.0 + lcs.score_);
            double density_score = static_cast<double>(lcs.inliner_) * score_weight;
            fitness += density_score;

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
                fitness -= 50.0 * static_cast<double>(neighbor_count[shard_idx] - kMaxNeighbors);
            }
        }

        map<pair<int, int>, Matrix4d> active_transforms;

        for (const auto& gene : chromosome.genes) {
            const vector<size_t>& group = pair_groups_[gene.group_idx];
            const LCSIndex& lcs = matches_[group[gene.choice - 1]];
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

        fitness -= 2.0 * cycle_penalty;

        struct PairEdge {
            int a;
            int b;
            Matrix4d T_ab;
            double weight;
        };

        map<pair<int, int>, PairEdge> pair_edges;

        for (const auto& gene : chromosome.genes) {
            const vector<size_t>& group = pair_groups_[gene.group_idx];
            const LCSIndex& lcs = matches_[group[gene.choice - 1]];
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

        fitness -= kEdgeResidualPenalty * edge_residual_penalty;
        fitness -= kEdgeRotResidualPenalty * edge_rot_residual_penalty;

        int largest_component = 0;
        int num_components = 0;

        AnalyzeConnectedComponents(adjacency, largest_component, num_components);
        fitness += kConnectivityReward * static_cast<double>(largest_component);

        if (num_components > 1) {
            fitness -= kConnectivityComponentPenalty * static_cast<double>(num_components - 1);
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

        map<size_t, int> map1, map2;
        for (const auto& g : parent1.genes) map1[g.group_idx] = g.choice;
        for (const auto& g : parent2.genes) map2[g.group_idx] = g.choice;

        set<size_t> all_groups;
        for (const auto& pair : map1) all_groups.insert(pair.first);
        for (const auto& pair : map2) all_groups.insert(pair.first);

        for (size_t group_idx : all_groups) {
            bool in1 = map1.count(group_idx);
            bool in2 = map2.count(group_idx);

            if (in1 && in2) {
                if (map1[group_idx] == map2[group_idx]) {
                    // Agree: preserve
                    child.genes.push_back({group_idx, map1[group_idx]});
                } else {
                    // Disagree: split randomly
                    int choice = (rand() % 2 == 0) ? map1[group_idx] : map2[group_idx];
                    child.genes.push_back({group_idx, choice});
                }
            } else {
                // Exclusive: inherit with 50% probability
                if (rand() % 2 == 0) {
                    int choice = in1 ? map1[group_idx] : map2[group_idx];
                    child.genes.push_back({group_idx, choice});
                }
            }
        }

        Repair(child);
        return child;
    }

    //-----------------------------------------------------------------------------------------------------------------//

    void Mutate(Chromosome& chromosome) const
    {
        // Probability map for current genes
        map<size_t, bool> active_map;
        for (const auto& g : chromosome.genes) active_map[g.group_idx] = true;

        // 1. Modify existing genes
        for (auto& gene : chromosome.genes) {
            if ((static_cast<double>(rand()) / RAND_MAX) < kMutationRate) {
                const vector<size_t>& group = pair_groups_[gene.group_idx];
                if (group.size() > 1) {
                    int old_choice = gene.choice;
                    int new_choice = SelectWeightedCandidate(gene.group_idx);
                    if (new_choice == old_choice) {
                        new_choice = (old_choice % static_cast<int>(group.size())) + 1;
                    }
                    gene.choice = new_choice;
                }
            }
        }

        // 2. Delete existing genes
        vector<Gene> filtered;
        for (const auto& gene : chromosome.genes) {
            if ((static_cast<double>(rand()) / RAND_MAX) >= kMutationRate) {
                filtered.push_back(gene);
            } else {
                active_map[gene.group_idx] = false;
            }
        }
        chromosome.genes = filtered;

        // 3. Add new genes for unused groups
        for (size_t i = 0; i < pair_groups_.size(); ++i) {
            if (!active_map[i]) {
                if ((static_cast<double>(rand()) / RAND_MAX) < kMutationRate) {
                    int choice = SelectWeightedCandidate(i);
                    if (choice > 0) {
                        chromosome.genes.push_back({i, choice});
                    }
                }
            }
        }

        Repair(chromosome);
    }

    //-----------------------------------------------------------------------------------------------------------------//



    //-----------------------------------------------------------------------------------------------------------------//

    void BuildOutputsFromSelection(const Chromosome& chromosome)
    {
        graph_ = MatrixXd::Zero(num_shards_, num_shards_);
        
        struct AdjEdge {
            int to;
            Matrix4d T_to_current;
            double weight; 
        };
        vector<vector<AdjEdge>> adjacency(num_shards_);
        
        for (const auto& gene : chromosome.genes) {
            const vector<size_t>& group = pair_groups_[gene.group_idx];
            const LCSIndex& lcs = matches_[group[gene.choice - 1]];
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
    static constexpr double kMutationRate = 0.08; 
    static constexpr double kInitialPairActivationRate = 0.3;
    static constexpr double kEdgeResidualThreshold = 50.0;
    static constexpr double kEdgeResidualPenalty = 0.05;
    static constexpr double kEdgeRotResidualThreshold = 0.35;
    static constexpr double kEdgeRotResidualPenalty = 1.0;
    static constexpr double kConnectivityReward = 10.0;
    static constexpr double kConnectivityComponentPenalty = 8.0;
    static constexpr int kPoseRelaxIterations = 5;
    static constexpr double kPoseRelaxAlpha = 0.5;

    vector<Geom> shard_;
    vector<LCSIndex> matches_;
    vector<vector<size_t>> pair_groups_;
    int num_shards_;

    vector<Chromosome> population_;

    vector<Trans> transforms_;
    MatrixXd graph_;
};

#endif

