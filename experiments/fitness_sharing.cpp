// #pragma once

// // Check # 01: Having both pragma once and this is redundant so the below two lines can be removed apparently
// #ifndef _GENETIC_ALGORITHM_H_   
// #define _GENETIC_ALGORITHM_H_   

// #include <algorithm>
// #include <cstdlib>
// #include <ctime>
// #include <iostream>
// #include <map>
// #include <queue>
// #include <vector>

// // Check # 02: #include "data_structure.h" should also work
// #include "../class/data_structure.h"

// // extern: there is a shard_on_off array somewhere in the codebase, I want to use it here, but I'm not defining it
// extern bool shard_on_off[];

// class GeneticAssembler {
// public:
//     struct Chromosome {
//         vector<int> genes;
//         double fitness;
//         double raw_fitness;  // fitness before sharing penalty
//     };

//     //-----------------------------------------------------------------------------------------------------------------//

//     GeneticAssembler(const vector<Geom>& shard, list<LCSIndex>& LCS_out, int num_shards) 
//         : shard_(shard), num_shards_(num_shards) 
//     {   

//         for (list<LCSIndex>::const_iterator it = LCS_out.begin(); it != LCS_out.end(); ++it) {

//             // Check # 05: We can try to change this threshold of 10.0 everywhere in the file
//             if (it->score_ >= 500.0) {
//                 continue;
//             }
//             matches_.push_back(*it);
//         }
//         BuildPairGroups();

//         graph_ = MatrixXd::Zero(num_shards_, num_shards_);
//         transforms_.resize(num_shards_);
//         Matrix3d I = Matrix3d::Identity();
//         Vector3d zero = Vector3d::Zero();
        
//         for (int i = 0; i < num_shards_; ++i) {
//             transforms_[i].Set(I, zero, i + 1, 1);
//         }
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

//     void Run()
//     {
//         srand(42);

//         if (pair_groups_.empty()) {
//             cout << "PAIR GROUPS ARE EMPTY...WHY?" << endl;
//             BuildOutputsFromSelection(vector<int>());
//             return;
//         }

//         if (!matches_.empty()) {
//             cout << "[GA DEBUG] First Match: shard_x_ = " << matches_[0].shard_x_ 
//             << " shard_y_ = " << matches_[0].shard_y_ << endl;
//         }

//         if (kUseFitnessSharing) {
//             cout << "[GA] Fitness sharing: ON (sigma=" << kSigmaShare << ")" << endl;
//         } else {
//             cout << "[GA] Fitness sharing: OFF" << endl;
//         }

//         InitializePopulation();

//         for (int generation = 0; generation < kMaxGenerations; ++generation) {
//             EvaluatePopulation();
            
//             if (kUseFitnessSharing) {
//                 // Apply sharing to promote diversity
//                 ApplyFitnessSharing();
//             }

//             sort(population_.begin(), population_.end(), [](const Chromosome& a, const Chromosome& b) {
//                 return a.fitness > b.fitness;
//             });

//             int active_genes = 0;
//             if (!population_.empty()) {
//                 const auto& best = population_.front();
//                 for (int g : best.genes) if (g > 0) active_genes++;
//             }

//             int viable_count = 0;
//             for (const auto& c : population_) if (c.raw_fitness > 0.0) viable_count++;

//             cout << "[GA Gen " << generation << "] "
//                  << "Best raw fitness: " << population_.front().raw_fitness
//                  << " Best shared fitness: " << population_.front().fitness
//                  << " Viable: " << viable_count << "/" << population_.size()
//                  << " Active genes: " << active_genes << endl;

//             vector<Chromosome> next_population;
//             next_population.reserve(kPopulationSize);

//             int elitism = min(kElitismCount, static_cast<int>(population_.size()));
//             for (int i = 0; i < elitism; ++i) {
//                 next_population.push_back(population_[i]);
//             }

//             while (static_cast<int>(next_population.size()) < kPopulationSize) {
//                 Chromosome parent1 = TournamentSelect();
//                 Chromosome parent2 = TournamentSelect();

//                 Chromosome child1 = Crossover(parent1, parent2);
//                 int mutated_gene1 = Mutate(child1);

//                 if (mutated_gene1 >= 0) {
//                     GuidedRepair(child1, static_cast<size_t>(mutated_gene1));
//                 }
//                 next_population.push_back(child1);

//                 if (static_cast<int>(next_population.size()) < kPopulationSize) {
//                     Chromosome child2 = Crossover(parent2, parent1);
//                     int mutated_gene2 = Mutate(child2);
//                     if (mutated_gene2 >= 0) {
//                         GuidedRepair(child2, static_cast<size_t>(mutated_gene2));
//                     }

//                     next_population.push_back(child2);
//                 }
//             }

//             if (kUseFitnessSharing) {
//                 // Restore raw fitness for the next generation's baseline evaluation
//                 for (auto& c : population_) c.fitness = c.raw_fitness;
//             }

//             population_ = next_population;

//         }

//         EvaluatePopulation();
//         sort(population_.begin(), population_.end(), [](const Chromosome& a, const Chromosome& b) {
//             return a.fitness > b.fitness;
//         });

//         const Chromosome& best = population_.front();
//         int active_count = 0;
//         for (size_t group_idx = 0; group_idx < best.genes.size(); ++group_idx) {
//             int choice = best.genes[group_idx];
//             if (choice <= 0) {
//                 continue;
//             }

//             const vector<size_t>& group = pair_groups_[group_idx];
//             int local_idx = choice - 1;
//             if (local_idx < 0 || local_idx >= static_cast<int>(group.size())) {
//                 continue;
//             }

//             const LCSIndex& lcs = matches_[group[local_idx]];
//             active_count++;
//             cout << "[GA DEBUG] Active match: shard_x_="
//                 << lcs.shard_x_
//                 << " shard_y_=" << lcs.shard_y_
//                 << " inliner_=" << lcs.inliner_
//                 << " score_=" << lcs.score_ << endl;
//         }

//         cout << "[GA DEBUG] Total active matches: " << active_count << endl;

//         BuildOutputsFromSelection(population_.front().genes);
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

//     vector<Trans> GetTransforms() const {
//         return transforms_;
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

//     MatrixXd GetGraph() const
//     {
//         return graph_;
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

//     double GetBestFitness() const
//     {
//         if (population_.empty()) return 0.0;

//         double best = population_[0].fitness;

//         for (size_t i = 1; i < population_.size(); ++i) {
//             if (population_[i].fitness > best) {
//                 best = population_[i].fitness;
//             }
//         }

//         return best;
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

//     vector<Chromosome> GetTopKDiverseChromosomes(int k)
//     {
//         vector<Chromosome> diverse;
//         if (population_.empty()) return diverse;

//         // Ensure sorted by raw fitness
//         vector<Chromosome> sorted = population_;
//         sort(sorted.begin(), sorted.end(), [](const Chromosome& a, const Chromosome& b) {
//             return a.raw_fitness > b.raw_fitness;
//         });

//         double niche_radius = kSigmaShare * pair_groups_.size();

//         for (const auto& c : sorted) {
//             bool is_diverse = true;
//             for (const auto& selected : diverse) {
//                 if (ComputeHammingDistance(c, selected) < niche_radius) {
//                     is_diverse = false;
//                     break;
//                 }
//             }
//             if (is_diverse) {
//                 diverse.push_back(c);
//             }
//             if (static_cast<int>(diverse.size()) >= k) break;
//         }
//         return diverse;
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

//     void BuildOutputsFromSelection(const vector<int>& genes)
//     {
//         graph_ = MatrixXd::Zero(num_shards_, num_shards_);
        
//         struct AdjEdge {
//             int to;
//             Matrix4d T_to_current;
//             double weight; 
//         };
//         vector<vector<AdjEdge>> adjacency(num_shards_);
        
//         for (size_t group_idx = 0; group_idx < genes.size(); ++group_idx) {
//             int choice = genes[group_idx];
//             if (choice <= 0) {
//                 continue;
//             }

//             const vector<size_t>& group = pair_groups_[group_idx];
//             int local_idx = choice - 1;
//             if (local_idx < 0 || local_idx >= static_cast<int>(group.size())) {
//                 continue;
//             }

//             const LCSIndex& lcs = matches_[group[local_idx]];
//             int x = lcs.shard_x_ - 1;
//             int y = lcs.shard_y_ - 1;
//             if (!IsShardValidAndOn(x) || !IsShardValidAndOn(y)) {
//                 continue;
//             }

//             if (lcs.score_ >= 100.0) {
//                 continue;
//             }

//             graph_(x, y) = 1;
//             graph_(y, x) = 1;

//             Matrix4d T_xy = Matrix4d::Identity();
//             lcs.trans_.Output(T_xy);
//             Matrix4d T_yx = T_xy.inverse();
//             double edge_weight = max(1.0, static_cast<double>(lcs.inliner_));

//             adjacency[x].push_back({ y, T_yx, edge_weight });
//             adjacency[y].push_back({ x, T_xy, edge_weight });
//         }

//         vector<bool> visited(num_shards_, false);
//         vector<Matrix4d> T_to_root(num_shards_, Matrix4d::Identity());

//         for (int i = 0; i < num_shards_; ++i) {
//             if (visited[i] || !IsShardValidAndOn(i)) {
//                 continue;
//             }

//             queue<int> q;
//             q.push(i);
//             visited[i] = true;

//             while (!q.empty()) {
//                 int current = q.front();
//                 q.pop();

//                 for (size_t j = 0; j < adjacency[current].size(); ++j) {
//                     int next = adjacency[current][j].to;
//                     if (visited[next]) {
//                         continue;
//                     }

//                     T_to_root[next] = T_to_root[current] * adjacency[current][j].T_to_current;
//                     visited[next] = true;
//                     q.push(next);
//                 }
//             }
//         }

//         for (int i = 0; i < num_shards_; ++i) {
//             transforms_[i].Input(T_to_root[i]);
//         }
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

// private:

//     void InitializePopulation()
//     {
//         population_.clear();
//         population_.reserve(kPopulationSize);

//         for (int i = 0; i < kPopulationSize; ++i) {
//             Chromosome chromosome;
//             chromosome.genes.resize(pair_groups_.size(), 0);
//             chromosome.fitness = 0.0;

//             for (size_t gene_idx = 0; gene_idx < chromosome.genes.size(); ++gene_idx) {
//                 chromosome.genes[gene_idx] = SampleGroupChoice(gene_idx);
//             }
//             population_.push_back(chromosome);
//         }
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

//     void EvaluatePopulation()
//     {
//         for (size_t i = 0; i < population_.size(); ++i) {
//             double f = EvaluateFitness(population_[i]);
//             population_[i].fitness = f;
//             population_[i].raw_fitness = f;
//         }
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

//     double EvaluateFitness(const Chromosome& chromosome) const
//     {   
//         int kMaxNeighbors = 0;
//         for (int i = 0; i < num_shards_; ++i) {
//             int count = 0;
//             for (size_t j = 0; j < matches_.size(); ++j) {
//                 if (matches_[j].shard_x_ - 1 == i || matches_[j].shard_y_ - 1 == i) {
//                     count++;
//                 }
//             }
//             kMaxNeighbors = max(kMaxNeighbors, count);
//         }

//         kMaxNeighbors = min(kMaxNeighbors, num_shards_ - 1);

//         double fitness = 0.0;
//         vector<int> neighbor_count(num_shards_, 0);
//         vector<vector<int>> adjacency(num_shards_);

//         for (size_t group_idx = 0; group_idx < chromosome.genes.size(); group_idx++) {
//             int choice = chromosome.genes[group_idx];
//             if (choice <= 0) {
//                 continue;
//             }

//             const vector<size_t>& group = pair_groups_[group_idx];
//             int local_idx = choice - 1;
//             if (local_idx < 0 || local_idx >= static_cast<int>(group.size())) {
//                 continue;
//             }

//             const LCSIndex& lcs = matches_[group[local_idx]];
//             int x = lcs.shard_x_ - 1;
//             int y = lcs.shard_y_ - 1;
//             if (!IsShardValidAndOn(x) || !IsShardValidAndOn(y)) {
//                 continue;
//             }
            
//             double score_weight = 1.0 / (1.0 + lcs.score_);
//             double density_score = static_cast<double>(lcs.inliner_) * score_weight;
//             fitness += density_score;

//             neighbor_count[x]++;
//             neighbor_count[y]++;
//             adjacency[x].push_back(y);
//             adjacency[y].push_back(x);
//         }

//         for (int shard_idx = 0; shard_idx < num_shards_; ++shard_idx) {
//             if (!IsShardValidAndOn(shard_idx)) {
//                 continue;
//             }
//             if (neighbor_count[shard_idx] > kMaxNeighbors) {
//                 fitness -= 50.0 * static_cast<double>(neighbor_count[shard_idx] - kMaxNeighbors);
//             }
//         }

//         map<pair<int, int>, Matrix4d> active_transforms;

//         for (size_t group_idx = 0; group_idx < chromosome.genes.size(); ++group_idx) {
//             int choice = chromosome.genes[group_idx];
//             if (choice <= 0) {
//                 continue;
//             }

//             const vector<size_t>& group = pair_groups_[group_idx];
//             int local_idx = choice - 1;
//             const LCSIndex& lcs = matches_[group[local_idx]];
//             int x = lcs.shard_x_ - 1;
//             int y = lcs.shard_y_ - 1;
//             if (!IsShardValidAndOn(x) || !IsShardValidAndOn(y)) {
//                 continue;
//             }

//             Matrix4d T_xy = Matrix4d::Identity();
//             lcs.trans_.Output(T_xy);

//             auto key = make_pair(min(x, y), max(x, y));
//             if (active_transforms.find(key) == active_transforms.end()) {
//                 active_transforms[key] = T_xy;
//             }
//         }

//         double cycle_penalty = 0.0;
//         const double cycle_threshold = 30.0;

//         for (int a = 0; a < num_shards_; ++a) {
//             if (!IsShardValidAndOn(a)) continue;
//             for (int b = a + 1; b < num_shards_; ++b) {
//                 if (!IsShardValidAndOn(b)) continue;
//                 for (int c = b + 1; c < num_shards_; ++c) {
//                     if (!IsShardValidAndOn(c)) continue;

//                     auto key_ab = make_pair(a, b);
//                     auto key_bc = make_pair(b, c);
//                     auto key_ac = make_pair(a, c);

//                     if (!active_transforms.count(key_ab) || !active_transforms.count(key_bc) || !active_transforms.count(key_ac)) {
//                         continue;
//                     }

//                     Matrix4d T_ab = active_transforms[key_ab];
//                     Matrix4d T_bc = active_transforms[key_bc];
//                     Matrix4d T_ac = active_transforms[key_ac];

//                     Matrix4d T_composed = T_ab * T_bc;
//                     Vector3d t_composed = T_composed.block<3, 1>(0, 3);
//                     Vector3d t_direct = T_ac.block<3, 1>(0, 3);
//                     double error = (t_composed - t_direct).norm();

//                     if (error > cycle_threshold) {
//                         cycle_penalty += error;
//                     }
//                 }
//             }
//         }

//         fitness -= 2.0 * cycle_penalty;

//         struct PairEdge {
//             int a;
//             int b;
//             Matrix4d T_ab;
//             double weight;
//         };

//         map<pair<int, int>, PairEdge> pair_edges;

//         for (size_t group_idx = 0; group_idx < chromosome.genes.size(); ++group_idx) {
//             int choice = chromosome.genes[group_idx];
//             if (choice <= 0) continue;

//             const vector<size_t>& group = pair_groups_[group_idx];
//             const LCSIndex& lcs = matches_[group[choice - 1]];
//             int x = lcs.shard_x_ - 1;
//             int y = lcs.shard_y_ - 1;
//             if (!IsShardValidAndOn(x) || !IsShardValidAndOn(y)) continue;

//             int a = min(x, y);
//             int b = max(x, y);
//             Matrix4d T_xy = Matrix4d::Identity();
//             lcs.trans_.Output(T_xy);

//             Matrix4d T_ab = (x == a && y == b) ? T_xy : T_xy.inverse();

//             auto key = make_pair(a, b);
//             auto it = pair_edges.find(key);
//             double w = max(1.0, static_cast<double>(lcs.inliner_));
//             if (it == pair_edges.end() || w > it->second.weight) {
//                 pair_edges[key] = { a, b, T_ab, w };
//             }
//         }

//         vector<vector<pair<int, Matrix4d>>> pose_adj(num_shards_);
//         for (const auto& kv : pair_edges) {
//             const PairEdge& edge = kv.second;
//             pose_adj[edge.a].push_back(make_pair(edge.b, edge.T_ab.inverse()));
//             pose_adj[edge.b].push_back(make_pair(edge.a, edge.T_ab));
//         }

//         double edge_residual_penalty = 0.0;
//         double edge_rot_residual_penalty = 0.0;
//         int evaluated_roots = 0;
        
//         for (int root = 0; root < num_shards_; ++root) {
//             if (!IsShardValidAndOn(root)) continue;

//             vector<bool> pose_visited(num_shards_, false);
//             vector<Matrix4d> T_pose(num_shards_, Matrix4d::Identity());
//             queue<int> q;
//             q.push(root);
//             pose_visited[root] = true;

//             while (!q.empty()) {
//                 int current = q.front();
//                 q.pop();

//                 for (size_t i = 0; i < pose_adj[current].size(); ++i) {
//                     int next = pose_adj[current][i].first;
//                     if (pose_visited[next]) continue;
//                     T_pose[next] = T_pose[current] * pose_adj[current][i].second;
//                     pose_visited[next] = true;
//                     q.push(next);
//                 }
//             }

//             bool root_used = false;
//             for (const auto& kv : pair_edges) {
//                 const PairEdge& edge = kv.second;
//                 if (!pose_visited[edge.a] || !pose_visited[edge.b]) continue;
//                 root_used = true;
//                 Matrix4d T_pred_ab = T_pose[edge.a].inverse() * T_pose[edge.b];
//                 double residual = (T_pred_ab.block<3, 1>(0, 3) - edge.T_ab.block<3, 1>(0, 3)).norm();
//                 if (residual > kEdgeResidualThreshold) {
//                     edge_residual_penalty += pow(residual - kEdgeResidualThreshold, 2);
//                 }
//                 Matrix4d T_diff = T_pred_ab * edge.T_ab.inverse();
//                 Matrix3d R_diff = T_diff.block<3, 3>(0, 0);
//                 Matrix3d log_R = R_diff.log();
//                 double rot_residual = Vector3d(-log_R(1, 2), log_R(0, 2), -log_R(0, 1)).norm();
//                 if (rot_residual > kEdgeRotResidualThreshold) {
//                     edge_rot_residual_penalty += (rot_residual - kEdgeRotResidualThreshold);
//                 }
//             }
//             if (root_used) evaluated_roots++;
//         }

//         if (evaluated_roots > 0) {
//             edge_residual_penalty /= static_cast<double>(evaluated_roots);
//             edge_rot_residual_penalty /= static_cast<double>(evaluated_roots);
//         }

//         fitness -= kEdgeResidualPenalty * edge_residual_penalty;
//         fitness -= kEdgeRotResidualPenalty * edge_rot_residual_penalty;

//         int largest_component = 0;
//         int num_components = 0;
//         AnalyzeConnectedComponents(adjacency, largest_component, num_components);
//         fitness += kConnectivityReward * static_cast<double>(largest_component);
//         if (num_components > 1) fitness -= kConnectivityComponentPenalty * static_cast<double>(num_components - 1);
        
//         return fitness;
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

//     void AnalyzeConnectedComponents(const vector<vector<int>>& adjacency, int& largest_component, int& num_components) const
//     {
//         largest_component = 0;
//         num_components = 0;
//         vector<bool> visited(num_shards_, false);

//         for (int start = 0; start < num_shards_; ++start) {
//             if (!IsShardValidAndOn(start) || visited[start]) continue;
//             num_components++;
//             queue<int> q;
//             q.push(start);
//             visited[start] = true;
//             int component_size = 0;
//             while (!q.empty()) {
//                 int current = q.front(); q.pop();
//                 component_size++;
//                 for (int next : adjacency[current]) {
//                     if (!IsShardValidAndOn(next) || visited[next]) continue;
//                     visited[next] = true;
//                     q.push(next);
//                 }
//             }
//             if (component_size > largest_component) largest_component = component_size;
//         }

//         for (int node = 0; node < num_shards_; ++node) {
//             if (!IsShardValidAndOn(node) || visited[node]) continue;
//             num_components++;
//             visited[node] = true;
//             if (largest_component < 1) largest_component = 1;
//         }
//         if (num_components == 0) largest_component = 0;
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

//     Chromosome TournamentSelect() const
//     {
//         int size = static_cast<int>(population_.size());
//         int best_index = rand() % size;
//         for (int i = 1; i < 3; ++i) {
//             int candidate_index = rand() % size;
//             if (population_[candidate_index].fitness > population_[best_index].fitness) {
//                 best_index = candidate_index;
//             }
//         }
//         return population_[best_index];
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

//     Chromosome Crossover(const Chromosome& parent1, const Chromosome& parent2) const
//     {
//         Chromosome child;
//         child.fitness = 0.0;
//         child.genes.resize(parent1.genes.size(), 0);
//         if (child.genes.empty()) return child;
//         int point = rand() % static_cast<int>(child.genes.size());
//         for (int i = 0; i < point; ++i) child.genes[i] = parent1.genes[i];
//         for (size_t i = point; i < child.genes.size(); ++i) child.genes[i] = parent2.genes[i];
//         return child;
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

//     int Mutate(Chromosome& chromosome) const
//     {
//         int selected_mutated_gene = -1;
//         int mutated_count = 0;
//         for (size_t i = 0; i < chromosome.genes.size(); ++i) {
//             if ((static_cast<double>(rand()) / static_cast<double>(RAND_MAX)) < kMutationRate) {
//                 chromosome.genes[i] = SampleGroupChoice(i);
//                 mutated_count++;
//                 if ((rand() % mutated_count) == 0) selected_mutated_gene = static_cast<int>(i);
//             }
//         }
//         return selected_mutated_gene;
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

//     void GuidedRepair(Chromosome& chromosome, size_t gene_idx) const
//     {
//         if (gene_idx >= chromosome.genes.size()) return;
//         const vector<size_t>& group = pair_groups_[gene_idx];
//         int num_options = 1 + static_cast<int>(group.size());
//         if (num_options <= 1) return;

//         int current_choice = chromosome.genes[gene_idx];
//         vector<int> alternatives;
//         for (int choice = 0; choice < num_options; ++choice) {
//             if (choice != current_choice) alternatives.push_back(choice);
//         }
//         random_shuffle(alternatives.begin(), alternatives.end());
//         int trials = min(kGuidedRepairTrials, static_cast<int>(alternatives.size()));
//         double best_fitness = EvaluateFitness(chromosome);
//         int best_choice = current_choice;

//         for (int trial = 0; trial < trials; ++trial) {
//             chromosome.genes[gene_idx] = alternatives[trial];
//             double candidate_fitness = EvaluateFitness(chromosome);
//             if (candidate_fitness > best_fitness) {
//                 best_fitness = candidate_fitness;
//                 best_choice = alternatives[trial];
//             }
//         }
//         chromosome.genes[gene_idx] = best_choice;
//         chromosome.fitness = best_fitness;
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

//     bool IsShardValidAndOn(int shard_idx) const
//     {
//         return (shard_idx >= 0) && (shard_idx < num_shards_) && shard_on_off[shard_idx];
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

//     void BuildPairGroups()
//     {
//         pair_groups_.clear();
//         map<pair<int, int>, vector<size_t>> grouped;
//         for (size_t i = 0; i < matches_.size(); ++i) {
//             int x = matches_[i].shard_x_ - 1;
//             int y = matches_[i].shard_y_ - 1;
//             if (x < 0 || x >= num_shards_ || y < 0 || y >= num_shards_ || x == y) continue;
//             grouped[make_pair(min(x, y), max(x, y))].push_back(i);
//         }
//         for (const auto& kv : grouped) pair_groups_.push_back(kv.second);
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

//     int SampleGroupChoice(size_t group_idx) const
//     {
//         if (group_idx >= pair_groups_.size()) return 0;
//         if ((static_cast<double>(rand()) / static_cast<double>(RAND_MAX)) >= kInitialPairActivationRate) return 0;
//         const vector<size_t>& group = pair_groups_[group_idx];
//         if (group.empty()) return 0;

//         double total_weight = 0.0;
//         vector<double> cumulative;
//         for (size_t i = 0; i < group.size(); ++i) {
//             total_weight += max(1.0, static_cast<double>(matches_[group[i]].inliner_));
//             cumulative.push_back(total_weight);
//         }
//         if (total_weight <= 0.0) return 1 + (rand() % static_cast<int>(group.size()));
//         double r = (static_cast<double>(rand()) / static_cast<double>(RAND_MAX)) * total_weight;
//         for (size_t i = 0; i < cumulative.size(); ++i) {
//             if (r <= cumulative[i]) return static_cast<int>(i) + 1;
//         }
//         return static_cast<int>(group.size());
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

//     void ApplyFitnessSharing()
//     {
//         if (!kUseFitnessSharing || population_.empty()) return;
//         double niche_radius = kSigmaShare * static_cast<double>(pair_groups_.size());
//         size_t N = population_.size();
//         vector<double> sharing_factors(N, 0.0);
//         for (size_t i = 0; i < N; ++i) {
//             double sum_sh = 0.0;
//             for (size_t j = 0; j < N; ++j) {
//                 double d = static_cast<double>(ComputeHammingDistance(population_[i], population_[j]));
//                 if (d < niche_radius) sum_sh += (1.0 - pow(d / niche_radius, kFitnessSharingAlpha));
//             }
//             sharing_factors[i] = max(1.0, sum_sh);
//         }
//         for (size_t i = 0; i < N; ++i) population_[i].fitness /= sharing_factors[i];
//     }

//     int ComputeHammingDistance(const Chromosome& c1, const Chromosome& c2) const
//     {
//         int distance = 0;
//         size_t len = min(c1.genes.size(), c2.genes.size());
//         for (size_t i = 0; i < len; ++i) if (c1.genes[i] != c2.genes[i]) distance++;
//         distance += abs(static_cast<int>(c1.genes.size()) - static_cast<int>(c2.genes.size()));
//         return distance;
//     }

// private:
//     static constexpr int kPopulationSize = 200;
//     static constexpr int kMaxGenerations = 200;
//     static constexpr int kElitismCount = 7;
//     static constexpr double kMutationRate = 0.25;
//     static constexpr int kGuidedRepairTrials = 3;
//     static constexpr double kInitialPairActivationRate = 0.3;
//     static constexpr double kEdgeResidualThreshold = 50.0;
//     static constexpr double kEdgeResidualPenalty = 0.05;
//     static constexpr double kEdgeRotResidualThreshold = 0.35;
//     static constexpr double kEdgeRotResidualPenalty = 1.0;
//     static constexpr double kConnectivityReward = 10.0;
//     static constexpr double kConnectivityComponentPenalty = 8.0;

//     static constexpr double kSigmaShare = 0.25;
//     static constexpr double kFitnessSharingAlpha = 1.0;
//     static constexpr bool kUseFitnessSharing = true;

//     vector<Geom> shard_;
//     vector<LCSIndex> matches_;
//     vector<vector<size_t>> pair_groups_;
//     int num_shards_;

//     vector<Chromosome> population_;

//     vector<Trans> transforms_;
//     MatrixXd graph_;
// };

// #endif
