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
// #include <random>

// // Check # 02: #include "data_structure.h" should also work
// #include "../class/data_structure.h"

// // extern: there is a shard_on_off array somewhere in the codebase, I want to use it here, but I'm not defining it
// extern bool shard_on_off[];

// class GeneticAssembler {
// public:
//     struct Chromosome {
//         vector<int> genes;
//         double fitness;
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

//         // """
//         // BuildPairGroups() populates the pair_groups_ vector like this:
//         // pair_groups_[0] → {2, 5, 11}    // all matches between shards (1,3)
//         // pair_groups_[1] → {0, 7}        // all matches between shards (2,4)
//         // pair_groups_[2] → {3, 6, 9, 14} // all matches between shards (3,7)
//         // ...
//         // """

//         graph_ = MatrixXd::Zero(num_shards_, num_shards_);
//         transforms_.resize(num_shards_);
//         Matrix3d I = Matrix3d::Identity();
//         Vector3d zero = Vector3d::Zero();
        
//         for (int i = 0; i < num_shards_; ++i) {
//             transforms_[i].Set(I, zero, i + 1, 1);
//         }
//         // Check # 04: The referenced relative to shard index is 1 here, which might need to be changed based on feedback from Check # 03
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

//     void Run()
//     {
//         if (pair_groups_.empty()) {
//             cout << "PAIR GROUPS ARE EMPTY...WHY?" << endl;
//             BuildOutputsFromSelection(vector<int>());
//             return;
//         }

//         if (!matches_.empty()) {
//             cout << "[GA DEBUG] First Match: shard_x_ = " << matches_[0].shard_x_ 
//             << " shard_y_ = " << matches_[0].shard_y_ << endl;
//         }

//         cout << "[GA] Mode: " << (kUseIslandModel ? "Island Model" : "Single Population") << endl;

//         if (kUseIslandModel) {
//             cout << "[GA] Islands: " << kNumIslands 
//                  << " Migration interval: " << kMigrationInterval
//                  << " Migration size: " << kMigrationSize << endl;

//             RunIslandModel();
//         } else {
//             RunSinglePopulation();
//         }
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
//         vector<Chromosome> candidates;

//         if (kUseIslandModel) {
//             // Collect best from each island
//             for (int i = 0; i < kNumIslands; ++i) {
//                 if (!island_populations_.empty() && !island_populations_[i].empty()) {
//                     // Sort island to ensure best is at front
//                     sort(island_populations_[i].begin(), island_populations_[i].end(),
//                          [](const Chromosome& a, const Chromosome& b) { return a.fitness > b.fitness; });
//                     // Take top few from each island
//                     int take = min(k, static_cast<int>(island_populations_[i].size()));
//                     for (int j = 0; j < take; ++j) {
//                         candidates.push_back(island_populations_[i][j]);
//                     }
//                 }
//             }
//         } else {
//             candidates = population_;
//         }

//         if (candidates.empty()) return candidates;

//         // Sort by fitness descending
//         sort(candidates.begin(), candidates.end(),
//              [](const Chromosome& a, const Chromosome& b) { return a.fitness > b.fitness; });

//         // Greedy diverse selection
//         double niche_radius = 0.25 * static_cast<double>(pair_groups_.size());
//         vector<Chromosome> diverse;
//         for (const auto& c : candidates) {
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

// private:

//     // ============ Thread-safe RNG ============

//     static mt19937& GetRNG()
//     {
//         thread_local mt19937 rng(random_device{}());
//         return rng;
//     }

//     static double RandDouble()
//     {
//         thread_local uniform_real_distribution<double> dist(0.0, 1.0);
//         return dist(GetRNG());
//     }

//     static int RandInt(int n)
//     {
//         if (n <= 0) return 0;
//         return static_cast<int>(GetRNG()() % static_cast<unsigned int>(n));
//     }

//     // ============ Single Population Path (original behaviour) ============

//     void RunSinglePopulation()
//     {
//         InitializePopulation(population_);

//         EvolveIsland(population_, kMaxGenerations, 0);

//         FinalizeResult(population_);
//     }

//     // ============ Island Model Path ============

//     void RunIslandModel()
//     {
//         // Initialize N independent populations
//         island_populations_.resize(kNumIslands);
//         for (int i = 0; i < kNumIslands; ++i) {
//             InitializePopulation(island_populations_[i]);
//         }

//         int num_epochs = kMaxGenerations / kMigrationInterval;
//         int remainder_gens = kMaxGenerations % kMigrationInterval;

//         for (int epoch = 0; epoch < num_epochs; ++epoch) {
//             cout << "[GA] === Epoch " << epoch + 1 << " / " << num_epochs << " ===" << endl;

//             // Evolve each island independently (parallelized with OpenMP)
//             #pragma omp parallel for num_threads(kNumIslands)
//             for (int i = 0; i < kNumIslands; ++i) {
//                 EvolveIsland(island_populations_[i], kMigrationInterval, i);
//             }

//             // Migration (must happen after all islands finish)
//             Migrate(island_populations_);

//             // Log diversity after migration
//             LogMigrationDiversity(island_populations_);
//         }

//         // Run remaining generations if kMaxGenerations is not evenly divisible
//         if (remainder_gens > 0) {
//             #pragma omp parallel for num_threads(kNumIslands)
//             for (int i = 0; i < kNumIslands; ++i) {
//                 EvolveIsland(island_populations_[i], remainder_gens, i);
//             }
//         }

//         // Find best chromosome across all islands
//         int best_island_id = 0;
//         double best_fitness = -1e18;
//         for (int i = 0; i < kNumIslands; ++i) {
//             if (!island_populations_[i].empty()) {
//                 sort(island_populations_[i].begin(), island_populations_[i].end(),
//                      [](const Chromosome& a, const Chromosome& b) { return a.fitness > b.fitness; });
//                 if (island_populations_[i][0].fitness > best_fitness) {
//                     best_fitness = island_populations_[i][0].fitness;
//                     best_island_id = i;
//                 }
//             }
//         }

//         cout << "[GA] Final best from Island " << best_island_id 
//              << " with fitness " << best_fitness << endl;

//         // Copy best island's population as the main population (for GetBestFitness compatibility)
//         population_ = island_populations_[best_island_id];

//         FinalizeResult(population_);
//     }

//     // ============ Core Evolution Loop ============

//     void EvolveIsland(vector<Chromosome>& population, int generations, int island_id)
//     {
//         // Evaluate initial population
//         EvaluatePopulation(population);

//         for (int generation = 0; generation < generations; ++generation) {
//             sort(population.begin(), population.end(), [](const Chromosome& a, const Chromosome& b) {
//                 return a.fitness > b.fitness;
//             });

//             if (kUseIslandModel) {
//                 cout << "[GA Island " << island_id << " Gen " << generation 
//                      << "] Best fitness: " << population.front().fitness << endl;
//             } else {
//                 cout << "[GA] Generation " << generation << " best fitness: " << population.front().fitness << endl;
//             }

//             vector<Chromosome> next_population;
//             next_population.reserve(kPopulationSize);

//             int elitism = min(kElitismCount, static_cast<int>(population.size()));
//             for (int i = 0; i < elitism; ++i) {
//                 next_population.push_back(population[i]);
//             }

//             while (static_cast<int>(next_population.size()) < kPopulationSize) {
//                 Chromosome parent1 = TournamentSelect(population);
//                 Chromosome parent2 = TournamentSelect(population);

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

//             population = next_population;

//             // Evaluate the new population for the next generation's sort
//             EvaluatePopulation(population);
//         }
//     }

//     // ============ Migration ============

//     void Migrate(vector<vector<Chromosome>>& island_populations)
//     {
//         // Collect migration buffers from each island (synchronous — all prepare before any receive)
//         vector<vector<Chromosome>> migration_buffers(kNumIslands);

//         for (int i = 0; i < kNumIslands; ++i) {
//             // Sort to ensure best are at front
//             sort(island_populations[i].begin(), island_populations[i].end(),
//                  [](const Chromosome& a, const Chromosome& b) { return a.fitness > b.fitness; });

//             int send_count = min(kMigrationSize, static_cast<int>(island_populations[i].size()));
//             for (int j = 0; j < send_count; ++j) {
//                 migration_buffers[i].push_back(island_populations[i][j]);
//             }
//         }

//         // Distribute: island i sends to island (i+1) % N
//         for (int i = 0; i < kNumIslands; ++i) {
//             int receiver = (i + 1) % kNumIslands;

//             // Replace worst chromosomes on receiver
//             int replace_count = min(kMigrationSize, static_cast<int>(island_populations[receiver].size()));
//             replace_count = min(replace_count, static_cast<int>(migration_buffers[i].size()));

//             for (int j = 0; j < replace_count; ++j) {
//                 // Replace from the back (worst)
//                 int back_idx = static_cast<int>(island_populations[receiver].size()) - 1 - j;
//                 if (back_idx >= 0) {
//                     island_populations[receiver][back_idx] = migration_buffers[i][j];
//                 }
//             }

//             // Re-sort receiver after insertion
//             sort(island_populations[receiver].begin(), island_populations[receiver].end(),
//                  [](const Chromosome& a, const Chromosome& b) { return a.fitness > b.fitness; });
//         }
//     }

//     // ============ Migration Logging ============

//     void LogMigrationDiversity(const vector<vector<Chromosome>>& island_populations) const
//     {
//         cout << "[GA Migration] Island fitnesses: ";
//         for (int i = 0; i < kNumIslands; ++i) {
//             if (!island_populations[i].empty()) {
//                 cout << "I" << i << "=" << island_populations[i][0].fitness << " ";
//             }
//         }
//         cout << endl;

//         // Pairwise diversity between best chromosomes
//         for (int i = 0; i < kNumIslands; ++i) {
//             for (int j = i + 1; j < kNumIslands; ++j) {
//                 if (!island_populations[i].empty() && !island_populations[j].empty()) {
//                     int dist = ComputeHammingDistance(island_populations[i][0], island_populations[j][0]);
//                     cout << "[GA Migration] Diversity I" << i << "-I" << j 
//                          << ": " << dist << " genes differ" << endl;
//                 }
//             }
//         }
//     }

//     // ============ Finalize Result ============

//     void FinalizeResult(vector<Chromosome>& population)
//     {
//         sort(population.begin(), population.end(), [](const Chromosome& a, const Chromosome& b) {
//             return a.fitness > b.fitness;
//         });

//         const Chromosome& best = population.front();
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

//         BuildOutputsFromSelection(population.front().genes);
//     }

//     // ============ Population Initialization ============

//     void InitializePopulation(vector<Chromosome>& population)
//     {
//         population.clear();
//         population.reserve(kPopulationSize);

//         for (int i = 0; i < kPopulationSize; ++i) {
//             Chromosome chromosome;
//             chromosome.genes.resize(pair_groups_.size(), 0);
//             chromosome.fitness = 0.0;

//             for (size_t gene_idx = 0; gene_idx < chromosome.genes.size(); ++gene_idx) {
//                 chromosome.genes[gene_idx] = SampleGroupChoice(gene_idx);
//             }
//             population.push_back(chromosome);
//         }
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

//     void EvaluatePopulation(vector<Chromosome>& population)
//     {
//         for (size_t i = 0; i < population.size(); ++i) {
//             population[i].fitness = EvaluateFitness(population[i]);
//         }
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

//     // Check # 13: The fitness function might need a thorough analysis again. Moreover, why are we not using functions in other files (feature_matchings, ranking_system, etc.) at all? We have only used LCS so far?
//     double EvaluateFitness(const Chromosome& chromosome) const
//     {   
                
//         int kMaxNeighbors = 0;
        
//         // Check # 07: This logic might be incorrect given that it can count each shard multiple times because shard pairs can repeat in matches_
//         for (int i = 0; i < num_shards_; ++i) {
//             int count = 0;
//             for (size_t j = 0; j < matches_.size(); ++j) {
//                 if (matches_[j].shard_x_ - 1 == i || matches_[j].shard_y_ - 1 == i) {
//                     count++;
//                 }
//             }
//             kMaxNeighbors = max(kMaxNeighbors, count);
//             // Check # 08: Why is there a need of this max() code above? It can be removed.
//         }

//         // Cap at num_shards_ - 1 (maximum possible connections)
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
            
//             // Check # 09: This fitness function might be incorrect. A higher inliner is better while a lower score is better which means that both contradict each other.
//             // fitness += static_cast<double>(lcs.inliner_) * lcs.score_;

//             // Implementing the above Check # 09:
//             // Higher inliers = better, lower score = better
//             // Invert score so both terms pull in the same direction
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
//             // Check # 10: This condition might be totally illogical because neighbor_count[shard_idx] can never exceed kMaxNeighbors as per the above implementation. A more useful function migh help here.
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
//             if (local_idx < 0 || local_idx >= static_cast<int>(group.size())) {
//                 continue;
//             }

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

//         // Check # 11: It can be tried to change the threshold or weight (currently -2.0) or re-verify the logic of this cycle consistency code
//         double cycle_penalty = 0.0;
//         const double cycle_threshold = 30.0;

//         for (int a = 0; a < num_shards_; ++a) {
//             if (!IsShardValidAndOn(a)) {
//                 continue;
//             }

//             for (int b = a + 1; b < num_shards_; ++b) {
//                 if (!IsShardValidAndOn(b)) {
//                     continue;
//                 }

//                 for (int c = b + 1; c < num_shards_; ++c) {
//                     if (!IsShardValidAndOn(c)) {
//                         continue;
//                     }

//                     auto key_ab = make_pair(a, b);
//                     auto key_bc = make_pair(b, c);
//                     auto key_ac = make_pair(a, c);

//                     bool has_ab = active_transforms.count(key_ab) > 0;
//                     bool has_bc = active_transforms.count(key_bc) > 0;
//                     bool has_ac = active_transforms.count(key_ac) > 0;

//                     if (!has_ab || !has_bc || !has_ac) {
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

//             int a = min(x, y);
//             int b = max(x, y);
//             Matrix4d T_xy = Matrix4d::Identity();
//             lcs.trans_.Output(T_xy);

//             Matrix4d T_ab = Matrix4d::Identity();
//             if (x == a && y == b) {
//                 T_ab = T_xy;
//             } else {
//                 T_ab = T_xy.inverse();
//             }

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
//             if (!IsShardValidAndOn(root)) {
//                 continue;
//             }

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
//                     const Matrix4d& T_next_current = pose_adj[current][i].second;
//                     if (pose_visited[next]) {
//                         continue;
//                     }

//                     T_pose[next] = T_pose[current] * T_next_current;
//                     pose_visited[next] = true;
//                     q.push(next);
//                 }
//             }

//             bool root_used = false;

//             for (const auto& kv : pair_edges) {
//                 const PairEdge& edge = kv.second;
//                 if (!pose_visited[edge.a] || !pose_visited[edge.b]) {
//                     continue;
//                 }

//                 root_used = true;

//                 Matrix4d T_pred_ab = T_pose[edge.a].inverse() * T_pose[edge.b];
//                 Vector3d t_pred = T_pred_ab.block<3, 1>(0, 3);
//                 Vector3d t_ab = edge.T_ab.block<3, 1>(0, 3);
//                 double residual = (t_pred - t_ab).norm();

//                 if (residual > kEdgeResidualThreshold) {
//                     double excess = residual - kEdgeResidualThreshold;
//                     edge_residual_penalty += excess * excess;
//                 }

//                 Matrix4d T_diff = T_pred_ab * edge.T_ab.inverse();
//                 Matrix3d R_diff;
//                 for (int r = 0; r < 3; ++r) {
//                     R_diff.row(r) << T_diff(r, 0), T_diff(r, 1), T_diff(r, 2);
//                 }

//                 Matrix3d log_R = R_diff.log();
//                 Vector3d w;
//                 w << -log_R(1, 2), log_R(0, 2), -log_R(0, 1);
//                 double rot_residual = w.norm();
                
//                 if (rot_residual > kEdgeRotResidualThreshold) {
//                     edge_rot_residual_penalty += (rot_residual - kEdgeRotResidualThreshold);
//                 }
//             }

//             if (root_used) {
//                 evaluated_roots++;
//             }
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

//         if (num_components > 1) {
//             fitness -= kConnectivityComponentPenalty * static_cast<double>(num_components - 1);
//         }
        
//         return fitness;
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

//     void AnalyzeConnectedComponents(const vector<vector<int>>& adjacency, int& largest_component, int& num_components) const
//     {
//         largest_component = 0;
//         num_components = 0;
//         vector<bool> visited(num_shards_, false);

//         for (int start = 0; start < num_shards_; ++start) {
//             if (!IsShardValidAndOn(start) || visited[start]) {
//                 continue;
//             }

//             num_components++;
//             queue<int> q;
//             q.push(start);
//             visited[start] = true;

//             int component_size = 0;

//             while (!q.empty()) {
//                 int current = q.front();
//                 q.pop();
//                 component_size++;

//                 for (size_t i = 0; i < adjacency[current].size(); ++i) {
//                     int next = adjacency[current][i];
//                     if (!IsShardValidAndOn(next) || visited[next]) {
//                         continue;
//                     }

//                     visited[next] = true;
//                     q.push(next);
//                 }
//             }
            
//             if (component_size > largest_component) {
//                 largest_component = component_size;
//             }
//         }

//         // Isolated valid sherds with no selected edges should still count as components.
//         // Check # 12: This loop might not be required because the above loop also handles those shards not connected with any other shards.
//         for (int node = 0; node < num_shards_; ++node) {
//             if (!IsShardValidAndOn(node) || visited[node]) {
//                 continue;
//             }

//             num_components++;
//             visited[node] = true;
            
//             if (largest_component < 1) {
//                 largest_component = 1;
//             }
//         }

//         if (num_components == 0) {
//             largest_component = 0;
//         }
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

//     Chromosome TournamentSelect(const vector<Chromosome>& population) const
//     {
//         int size = static_cast<int>(population.size());
//         int best_index = RandInt(size);

//         // Check # 13: The tournament size (2 here) can be experimented with...
//         for (int i = 1; i < 3; ++i) {
//             int candidate_index = RandInt(size);
//             if (population[candidate_index].fitness > population[best_index].fitness) {
//                 best_index = candidate_index;
//             }
//         }

//         return population[best_index];
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

//     Chromosome Crossover(const Chromosome& parent1, const Chromosome& parent2) const
//     {
//         Chromosome child;
//         child.fitness = 0.0;
//         child.genes.resize(parent1.genes.size(), 0);

//         if (child.genes.empty()) {
//             return child;
//         }

//         int point = RandInt(static_cast<int>(child.genes.size()));

//         for (int i = 0; i < point; ++i) {
//             child.genes[i] = parent1.genes[i];
//         }

//         for (size_t i = point; i < child.genes.size(); ++i) {
//             child.genes[i] = parent2.genes[i];
//         }

//         return child;
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

//     // Check # 16: There is a high potential in increasing the chances of mutation as based on this code, there are very less. Multiple genes can be mutated at once as well.
//     int Mutate(Chromosome& chromosome) const
//     {
//         int selected_mutated_gene = -1;
//         int mutated_count = 0;

//         for (size_t i = 0; i < chromosome.genes.size(); ++i) {

//             // Check # 14: This might be a poor randomness technique which has to upgraded with moder C++ standards
//             double r = RandDouble();

//             if (r < kMutationRate) {
//                 int old_choice = chromosome.genes[i];
//                 const vector<size_t>& group = pair_groups_[i];
//                 int num_options = 1 + static_cast<int>(group.size()); // 0 (inactive) + candidates

//                 if (num_options <= 1) {
//                     continue;
//                 }

//                 int new_choice = old_choice;
//                 const int kMaxResampleAttempts = 16;

//                 for (int attempt = 0; attempt < kMaxResampleAttempts; ++attempt) {
//                     new_choice = SampleGroupChoice(i);
//                     if (new_choice != old_choice) {
//                         break;
//                     }
//                 }

//                 if (new_choice == old_choice) {
//                     int clamped_old = old_choice;
//                     if (clamped_old < 0 || clamped_old >= num_options) {
//                         clamped_old = 0;
//                     }

//                     int draw = RandInt(num_options - 1);
//                     if (draw >= clamped_old) {
//                         draw++;
//                     }
//                     new_choice = draw;
//                 }

//                 chromosome.genes[i] = new_choice;
//                 mutated_count++;

//                 if (RandInt(mutated_count) == 0) {
//                     selected_mutated_gene = static_cast<int>(i);
//                 }
//             }
//         }

//         return selected_mutated_gene;
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

//     void GuidedRepair(Chromosome& chromosome, size_t gene_idx) const
//     {
//         if (gene_idx >= chromosome.genes.size() || gene_idx >= pair_groups_.size()) {
//             return;
//         }

//         const vector<size_t>& group = pair_groups_[gene_idx];
//         int num_options = 1 + static_cast<int>(group.size());
        
//         if (num_options <= 1) {
//             return;
//         }

//         int current_choice = chromosome.genes[gene_idx];
//         vector<int> alternatives;
//         alternatives.reserve(num_options - 1);

//         for (int choice = 0; choice < num_options; ++choice) {
//             if (choice == current_choice) {
//                 continue;
//             }
//             alternatives.push_back(choice);
//         }

//         if (alternatives.empty()) {
//             return;
//         }

//         // Shuffle using thread-safe RNG
//         for (int i = static_cast<int>(alternatives.size()) - 1; i > 0; --i) {
//             int j = RandInt(i + 1);
//             swap(alternatives[i], alternatives[j]);
//         }

//         int trials = min(kGuidedRepairTrials, static_cast<int>(alternatives.size()));

//         if (trials <= 0) {
//             return;
//         }

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

//         int bfs_root = -1;
//         for (int i = 0; i < num_shards_; i++) {
//             if (IsShardValidAndOn(i)) {
//                 bfs_root = i;
//                 break;
//             }
//         }

//         if (bfs_root >= 0) {
//             queue<int> q;
//             q.push(bfs_root);
//             visited[bfs_root] = true;

//             while (!q.empty()) {
//                 int current = q.front();
//                 q.pop();

//                 for (size_t i = 0; i < adjacency[current].size(); ++i) {
//                     int next = adjacency[current][i].to;
//                     const Matrix4d& T_next_current = adjacency[current][i].T_to_current; 
//                     if (visited[next]) {
//                         continue;
//                     }

//                     T_to_root[next] = T_to_root[current] * T_next_current;
//                     visited[next] = true;
//                     q.push(next);
//                 }
//             }

//             for (int iter = 0; iter < kPoseRelaxIterations; ++iter) {
//                 vector<Matrix4d> updated = T_to_root;
//                 bool changed = false;

//                 for (int node = 0; node < num_shards_; ++node) {
//                     if (!visited[node] || node == bfs_root || !IsShardValidAndOn(node)) {
//                         continue;
//                     }

//                     Vector3d t_acc = Vector3d::Zero();
//                     double w_sum = 0.0;

//                     for (int from = 0; from < num_shards_; ++from) {
//                         if (!visited[from] || !IsShardValidAndOn(from)) {
//                             continue;
//                         }

//                         for (size_t k = 0; k < adjacency[from].size(); ++k) {
//                             if (adjacency[from][k].to != node) {
//                                 continue;
//                             }

//                             Matrix4d pred = T_to_root[from] * adjacency[from][k].T_to_current;
//                             Vector3d t_pred = pred.block<3, 1>(0, 3);
//                             double w = adjacency[from][k].weight;
//                             t_acc += w * t_pred;
//                             w_sum += w;
//                         }
//                     }

//                     if (w_sum <= 0.0) {
//                         continue;
//                     }

//                     Vector3d t_old = T_to_root[node].block<3, 1>(0, 3);
//                     Vector3d t_est = t_acc / w_sum;
//                     Vector3d t_new = ((1.0 - kPoseRelaxAlpha) * t_old) +  (kPoseRelaxAlpha * t_est);

//                     if ((t_new - t_old).norm() > 1.0e-6) {
//                         changed = true;
//                     }

//                     updated[node].block<3, 1>(0, 3) = t_new;
//                 }

//                 T_to_root = updated;
//                 if (!changed) {
//                     break;
//                 }
//             }
//         }

//         int reached = 0;
//         for (int i = 0; i < num_shards_; i++) {
//             if (visited[i]) {
//                 reached++;
//             }
//         }

//         cout << "[GA DEBUG] BFS reached " << reached << " / " << num_shards_ << " sherds" << endl;

//         for (int i = 0; i < num_shards_; i++) {
//             if (!IsShardValidAndOn(i)) continue;
//             cout << "[GA DEBUG] T_to_root[" << i << "] translations: " 
//             << T_to_root[i](0, 3) << " " 
//             << T_to_root[i](1, 3) << " "
//             << T_to_root[i](2, 3) << endl;
//         }

//         Matrix3d I = Matrix3d::Identity();
//         Vector3d zero = Vector3d::Zero();
//         for (int i = 0; i < num_shards_; ++i) {
//             transforms_[i].Set(I, zero, i + 1, 1);
//             if (!IsShardValidAndOn(i)) {
//                 continue;
//             }
//             transforms_[i].Set(T_to_root[i], i + 1, 1);
//         }   
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
//             if (x < 0 || x >= num_shards_ || y < 0 || y >= num_shards_) {
//                 continue;
//             }
//             if (x == y) {
//                 continue;
//             }

//             pair<int, int> key = make_pair(min(x, y), max(x, y));
//             grouped[key].push_back(i);
//         }

//         for (map<pair<int, int>, vector<size_t>>::const_iterator it = grouped.begin(); it != grouped.end(); ++it) {
//             pair_groups_.push_back(it->second);
//         }
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

//     int SampleGroupChoice(size_t group_idx) const
//     {
//         if (group_idx >= pair_groups_.size()) {
//             return 0;
//         }

//         double activate_prob = RandDouble();
//         if (activate_prob >= kInitialPairActivationRate) {
//             return 0;
//         }

//         const vector<size_t>& group = pair_groups_[group_idx];
//         if (group.empty()) {
//             return 0;
//         }

//         double total_weight = 0.0;
//         vector<double> cumulative;
//         cumulative.reserve(group.size());

//         for (size_t i = 0; i < group.size(); ++i) {
//             const LCSIndex& lcs = matches_[group[i]];
//             double w = max(1.0, static_cast<double>(lcs.inliner_));
//             total_weight += w;
//             cumulative.push_back(total_weight);
//         }

//         if (total_weight <= 0.0) {
//             return 1 + RandInt(static_cast<int>(group.size()));
//         }

//         double r = RandDouble() * total_weight;
//         for (size_t i = 0; i < cumulative.size(); ++i) {
//             if (r <= cumulative[i]) {
//                 return static_cast<int>(i) + 1;
//             }
//         }

//         return static_cast<int>(group.size());
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

//     int ComputeHammingDistance(const Chromosome& c1, const Chromosome& c2) const
//     {
//         int distance = 0;
//         size_t len = min(c1.genes.size(), c2.genes.size());
//         for (size_t i = 0; i < len; ++i) {
//             if (c1.genes[i] != c2.genes[i]) {
//                 distance++;
//             }
//         }
//         distance += abs(static_cast<int>(c1.genes.size()) - static_cast<int>(c2.genes.size()));
//         return distance;
//     }

//     //-----------------------------------------------------------------------------------------------------------------//

// private:
//     // Original GA parameters
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
//     static constexpr int kPoseRelaxIterations = 5;
//     static constexpr double kPoseRelaxAlpha = 0.5;

//     // Island Model parameters
//     static constexpr int kNumIslands = 4;
//     static constexpr int kMigrationInterval = 20;   // migrate every N generations
//     static constexpr int kMigrationSize = 3;         // number of chromosomes exchanged
//     static constexpr bool kUseIslandModel = true;    // flag for comparison with single population

//     vector<Geom> shard_;
//     vector<LCSIndex> matches_;
//     vector<vector<size_t>> pair_groups_;
//     int num_shards_;

//     vector<Chromosome> population_;
//     vector<vector<Chromosome>> island_populations_;  // one population per island

//     vector<Trans> transforms_;
//     MatrixXd graph_;
// };

// #endif

