import re

with open("main.cpp", "r") as f:
    text = f.read()

# Insert the structs and functions right before int main()
insertion = """//------------------------------------------------------------------------------------------------------------------------//

struct AblationConfig {
    std::string name;
    bool use_inlier_score;
    bool use_connectivity_reward;
    bool use_cycle_penalty;
    bool use_edge_residual;
    bool use_rot_residual;
    bool use_neighbor_penalty;
    bool is_icp_refine;
};

std::vector<AblationConfig> ablation_experiments = {
    { "all_components_no_icp",      true,  true,  true,  true,  true,  true,  false },
    { "all_components_icp",         true,  true,  true,  true,  true,  true,  true  },
    { "no_inlier_score_no_icp",     false, true,  true,  true,  true,  true,  false },
    { "no_inlier_score_icp",        false, true,  true,  true,  true,  true,  true  },
    { "no_connectivity_reward_no",  true,  false, true,  true,  true,  true,  false },
    { "no_connectivity_reward_icp", true,  false, true,  true,  true,  true,  true  },
    { "no_cycle_penalty_no_icp",    true,  true,  false, true,  true,  true,  false },
    { "no_cycle_penalty_icp",       true,  true,  false, true,  true,  true,  true  },
    { "no_edge_residual_no_icp",    true,  true,  true,  false, true,  true,  false },
    { "no_edge_residual_icp",       true,  true,  true,  false, true,  true,  true  },
    { "no_rot_residual_no_icp",     true,  true,  true,  true,  false, true,  false },
    { "no_rot_residual_icp",        true,  true,  true,  true,  false, true,  true  },
    { "no_neighbor_penalty_no_icp", true,  true,  true,  true,  true,  false, false },
    { "no_neighbor_penalty_icp",    true,  true,  true,  true,  true,  false, true  },
    { "inlier_score_only_no_icp",   true,  false, false, false, false, false, false },
    { "inlier_score_only_icp",      true,  false, false, false, false, false, true  },
    { "no_penalties_no_icp",        true,  true,  false, false, false, false, false },
    { "no_penalties_icp",           true,  true,  false, false, false, false, true  },
};

const int kAblationTrials = 5;
const std::vector<unsigned int> kAblationSeeds = {42, 123, 456, 789, 1024};

std::vector<Trans> RunICPRefinement(
    const std::list<LCSIndex>& LCS_original,
    const std::vector<Geom>& shard_original,
    const Eigen::MatrixXd& graph_ga,
    const std::vector<Trans>& T_axis,
    const std::vector<Trans>& T_ga_eval
) {
    RankingSubgraph graph_icp(LCS_original, SHARD_NUMBER);
    graph_icp.node_[0] = true;
    graph_icp.root_node_ = 1;
    graph_icp.ResetMatchedIndex(shard_original);

    std::vector<Geom> shard_icp = shard_original;

    std::vector<int> placement_order;
    std::vector<bool> bfs_placed(SHARD_NUMBER, false);
    bfs_placed[0] = true;
    placement_order.push_back(0);
    std::queue<int> bfs_queue;
    bfs_queue.push(0);
    while (!bfs_queue.empty()) {
        int curr = bfs_queue.front(); bfs_queue.pop();
        for (int j = 0; j < SHARD_NUMBER; j++) {
            if (!bfs_placed[j] && (graph_ga(curr, j) > 0 || graph_ga(j, curr) > 0)
                && shard_on_off[j]) {
                bfs_placed[j] = true;
                placement_order.push_back(j);
                bfs_queue.push(j);
            }
        }
    }

    for (int pi = 1; pi < (int)placement_order.size(); pi++) {
        int current = placement_order[pi]; // 0-indexed

        shard_icp = shard_original;

        for (int k = 0; k < SHARD_NUMBER; k++) {
            if (!graph_icp.node_[k]) continue;
            Eigen::Matrix3d R_k = Eigen::Matrix3d::Identity();
            Eigen::Vector3d t_k = Eigen::Vector3d::Zero();
            graph_icp.T_[k].Output(R_k, t_k);
            shard_icp[k].Move(R_k, t_k);
        }

        std::vector<RankingSubgraph> single_graph_vec;
        single_graph_vec.push_back(graph_icp);
        single_graph_vec[0].MakeHierarchyPriorityList(0, single_graph_vec);
        graph_icp = single_graph_vec[0];

        std::vector<Chunk> filtered_priority;
        for (auto& chunk : graph_icp.priority_list_) {
            bool has_root_edge = false;
            for (int ei : chunk.i_edge) {
                int sx = graph_icp.sub_graph_[ei].shard_x_;
                int sy = graph_icp.sub_graph_[ei].shard_y_;
                if (sx == 1 || sy == 1) {
                    has_root_edge = true;
                    break;
                }
            }
            if (has_root_edge)
                filtered_priority.push_back(chunk);
        }
        if (!filtered_priority.empty())
            graph_icp.priority_list_ = filtered_priority;

        if (graph_icp.priority_list_.empty()) {
            graph_icp.node_[current] = true;
            continue;
        }

        bool found = false;
        for (int ci = 0; ci < (int)graph_icp.priority_list_.size(); ci++) {
            if (graph_icp.priority_list_[ci].node == current + 1) {
                graph_icp.priority_index_ = ci;
                found = true;
                break;
            }
        }
        if (!found) {
            graph_icp.node_[current] = true;
            continue;
        }

        std::vector<TransHistory> history;
        PrepareGraphBuilding(shard_icp, graph_icp, history, 0);

        RankingSubgraph pregraph_step(SHARD_NUMBER);
        for (int k = 0; k < SHARD_NUMBER; k++) {
            if (k == current) continue; 
            if (graph_icp.node_[k])
                pregraph_step.node_[k] = true;
        }
        std::vector<RankingSubgraph> pregraph_vec = { pregraph_step };

        std::vector<Eigen::Matrix3d> R_step(SHARD_NUMBER, Eigen::Matrix3d::Identity());
        std::vector<Eigen::Vector3d> t_step(SHARD_NUMBER, Eigen::Vector3d::Zero());
        int inlier_step = 0;
        IcpIncGraphAxis(shard_icp, R_step, t_step, graph_icp,
                        pregraph_vec, inlier_step, true, true);

        for (int k = 0; k < SHARD_NUMBER; k++) {
            if (!graph_icp.node_[k]) continue;
            graph_icp.T_[k].Input(R_step[k], t_step[k]);
        }
    }

    std::vector<Trans> T_ga_eval_icp = T_axis;
    for (int k = 0; k < SHARD_NUMBER; k++) {
        if (!shard_on_off[k]) continue;
        if (k == 0) continue;
        Eigen::Matrix3d R_k = Eigen::Matrix3d::Identity();
        Eigen::Vector3d t_k = Eigen::Vector3d::Zero();
        graph_icp.T_[k].Output(R_k, t_k);
        T_ga_eval_icp[k].Input(R_k, t_k);
    }
    
    return T_ga_eval_icp;
}

//------------------------------------------------------------------------------------------------------------------------//

int main(int argc, char** argv) 
"""

text = text.replace("//------------------------------------------------------------------------------------------------------------------------//\n\nint main(int argc, char** argv) \n", insertion)

# Insert Ablation Loop Right before the genetic algorithm iterations
ablation_loop = """
//------------------------------------------------------------------------------------------------------------------//

cout << "#################### Ablation Experiments ####################" << endl;

    for (const AblationConfig& config : ablation_experiments) {
        cout << "=== Ablation experiment: " << config.name << " ===" << endl;

        std::vector<int> sherd_correct(kAblationTrials, 0);
        std::vector<int> sherd_total(kAblationTrials, 0);
        std::vector<int> edge_correct(kAblationTrials, 0);
        std::vector<int> edge_total(kAblationTrials, 0);
        std::vector<double> best_fitness(kAblationTrials, 0.0);

        for (int trial = 0; trial < kAblationTrials; ++trial) {
            unsigned int seed = kAblationSeeds[trial];
            srand(seed);

            cout << "  Trial " << trial + 1 
                 << " / " << kAblationTrials 
                 << " (seed=" << seed << ")" << endl;

            GeneticAssembler ga(shard, LCS_original, SHARD_NUMBER);
            ga.use_inlier_score        = config.use_inlier_score;
            ga.use_connectivity_reward = config.use_connectivity_reward;
            ga.use_cycle_penalty       = config.use_cycle_penalty;
            ga.use_edge_residual       = config.use_edge_residual;
            ga.use_rot_residual        = config.use_rot_residual;
            ga.use_neighbor_penalty    = config.use_neighbor_penalty;
            ga.Run();

            std::vector<Trans> T_ga_eval_abl = T_axis;
            std::vector<Trans> T_ga_abl = ga.GetTransforms();
            for (int i = 0; i < SHARD_NUMBER; i++) {
                if (!shard_on_off[i]) continue;
                Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
                Eigen::Vector3d t = Eigen::Vector3d::Zero();
                T_ga_abl[i].Output(R, t);
                T_ga_eval_abl[i] = T_axis[i];
                T_ga_eval_abl[i].Input(R, t);
            }

            std::vector<Trans> final_T_eval = T_ga_eval_abl;
            if (config.is_icp_refine) {
                final_T_eval = RunICPRefinement(LCS_original, shard_original, ga.GetGraph(), T_axis, T_ga_eval_abl);
            }

            std::vector<bool> right_sherd_abl(SHARD_NUMBER, true);
            for (int i = 0; i < SHARD_NUMBER; i++) {
                if (!shard_on_off[i]) right_sherd_abl[i] = false;
            }

            auto [ks, ts, ke, te] = CountResult(
                GT_graph, GT_trans,
                ga.GetGraph(), final_T_eval,
                right_sherd_abl);

            sherd_correct[trial] = ks;
            sherd_total[trial]   = ts;
            edge_correct[trial]  = ke;
            edge_total[trial]    = te;
            best_fitness[trial]  = ga.GetBestFitness();

            cout << "  Sherd: " << ks << "/" << ts
                 << "  Edge: "  << ke << "/" << te
                 << "  Fitness: " << ga.GetBestFitness() << endl;

            std::string trial_path = path + "Result/ablation/" 
                              + config.name + "/trial_" 
                              + std::to_string(trial + 1) + "/";

            std::experimental::filesystem::create_directories(trial_path);

            // Reconstruct pc_origin properly
            std::vector<Visualize> pc_origin_trial(SHARD_NUMBER);
            for (int i = 0; i < SHARD_NUMBER; i++) {
                if (shard_on_off[i]) {
                    Eigen::Matrix3d R_vis = Eigen::Matrix3d::Identity();
                    Eigen::Vector3d t_vis = Eigen::Vector3d::Zero();
                    final_T_eval[i].Output(R_vis, t_vis);

                    std::string pointName = "origin_" + std::to_string(i + 1);
                    pc_origin_trial[i].MakePointCloud(shard_original[i].edge_line_.point_, shard_original[i].edge_line_.normal_, pointName);
                    pointName = "o_Mesh" + std::to_string(i + 1);
                    pc_origin_trial[i].MakeMesh(obj_path[i], pointName);

                    // Actually, SaveResult expects the points to be already transformed
                    // I will just use pc_origin directly as requested by the user,
                    // But first apply transformations.
                }
            }

            // The user requested: SaveResult(pc_origin, ga.GetGraph(), T_ga_eval_abl, shard_on_off, trial_path, SHARD_NUMBER);
            SaveResult(pc_origin, ga.GetGraph(), final_T_eval,
                       shard_on_off, trial_path, SHARD_NUMBER);

            shard = shard_original;
        }

        double avg_sherd = 0.0, avg_edge = 0.0, avg_fitness = 0.0;
        for (int t = 0; t < kAblationTrials; t++) {
            avg_sherd   += static_cast<double>(sherd_correct[t]) 
                         / static_cast<double>(sherd_total[t]);
            avg_edge    += static_cast<double>(edge_correct[t])  
                         / static_cast<double>(edge_total[t]);
            avg_fitness += best_fitness[t];
        }
        avg_sherd   /= kAblationTrials;
        avg_edge    /= kAblationTrials;
        avg_fitness /= kAblationTrials;

        cout << "  [" << config.name << "] "
             << "Avg sherd accuracy: " << avg_sherd * 100.0 << "% "
             << "Avg edge accuracy: "  << avg_edge  * 100.0 << "% "
             << "Avg fitness: "        << avg_fitness << endl;

        std::string summary_path = path + "Result/ablation/" 
                            + config.name + "/summary.txt";
        std::ofstream summary(summary_path);
        summary << "Experiment: " << config.name << endl;
        summary << "Trials: " << kAblationTrials << endl;
        summary << "Seeds: ";
        for (unsigned int s : kAblationSeeds) summary << s << " ";
        summary << endl;
        summary << "Components active:" << endl;
        summary << "  use_inlier_score        = " << config.use_inlier_score        << endl;
        summary << "  use_connectivity_reward = " << config.use_connectivity_reward << endl;
        summary << "  use_cycle_penalty       = " << config.use_cycle_penalty       << endl;
        summary << "  use_edge_residual       = " << config.use_edge_residual       << endl;
        summary << "  use_rot_residual        = " << config.use_rot_residual        << endl;
        summary << "  use_neighbor_penalty    = " << config.use_neighbor_penalty    << endl;
        summary << "  is_icp_refine           = " << config.is_icp_refine           << endl;
        summary << endl;
        for (int t = 0; t < kAblationTrials; t++) {
            summary << "Trial " << t + 1
                    << " seed=" << kAblationSeeds[t]
                    << " sherd=" << sherd_correct[t] << "/" << sherd_total[t]
                    << " edge="  << edge_correct[t]  << "/" << edge_total[t]
                    << " fitness=" << best_fitness[t] << endl;
        }
        summary << endl;
        summary << "Average sherd accuracy : " << avg_sherd   * 100.0 << "%" << endl;
        summary << "Average edge accuracy  : " << avg_edge    * 100.0 << "%" << endl;
        summary << "Average fitness        : " << avg_fitness          << endl;
        summary.close();

        cout << "  Summary saved to: " << summary_path << endl;
    }

//------------------------------------------------------------------------------------------------------------------//

cout << "#################### Genetic Algorithm search ####################" << endl;
"""

text = text.replace("""// Save matches computed on original axis-aligned positions
// (LCS_out gets overwritten during iterative GA, so preserve it here)

list<LCSIndex> LCS_original = LCS_out;

//------------------------------------------------------------------------------------------------------------------//

cout << "#################### Genetic Algorithm search ####################" << endl;""", """// Save matches computed on original axis-aligned positions
// (LCS_out gets overwritten during iterative GA, so preserve it here)

list<LCSIndex> LCS_original = LCS_out;""" + ablation_loop)

# Now remove the old ICP logic (from lines 369 to 552) to replace with a simple run of `RunICPRefinement` or just remove it using regex since we already extracted it.
# Wait, actually replacing the remaining big ICP block.
# Let's cleanly replace the ICP block.
text_to_replace = """//------------------------------------------------------------------------------------------------------------------//

// ---- Post-GA Incremental ICP Refinement (replicating BAISER flow) ----
    bool is_icp_refine = false; // Set to false to skip ICP refinement and use EA result as-is

    if (is_icp_refine) {"""

end_str = """
// Evaluate fine result separately
vector<bool> right_sherd_fine(SHARD_NUMBER, true);"""

idx1 = text.find(text_to_replace)
idx2 = text.find(end_str)

if idx1 != -1 and idx2 != -1:
    new_icp_block = """//------------------------------------------------------------------------------------------------------------------//

// ---- Post-GA Incremental ICP Refinement (replicating BAISER flow) ----
    bool is_icp_refine = false; // Set to false to skip ICP refinement and use EA result as-is
    int k_sherd_icp = 0, t_sherd_icp = 0, k_edge_icp = 0, t_edge_icp = 0;
    
    if (is_icp_refine) {
        cout << "#################### Post-GA Incremental ICP Refinement ####################" << endl;
        std::vector<Trans> T_ga_eval_icp = RunICPRefinement(LCS_original, shard_original, graph_ga, T_axis, T_ga_eval);
        
        vector<bool> right_sherd_icp(SHARD_NUMBER, true);
        for (int k = 0; k < SHARD_NUMBER; k++)
            if (!shard_on_off[k]) right_sherd_icp[k] = false;

        auto result = CountResult(GT_graph, GT_trans, graph_ga, T_ga_eval_icp, right_sherd_icp);
        k_sherd_icp = std::get<0>(result);
        t_sherd_icp = std::get<1>(result);
        k_edge_icp = std::get<2>(result);
        t_edge_icp = std::get<3>(result);

        cout << "########## Post-ICP Results ##########" << endl;
        cout << "ICP Sherd Accuracy : " << k_sherd_icp << " / " << t_sherd_icp << endl;
        cout << "ICP Edge Accuracy  : " << k_edge_icp << " / " << t_edge_icp << endl;
        cout << "######################################" << endl;
"""
    # Replace from idx1 to idx2
    text = text[:idx1] + new_icp_block + text[idx2:]

with open("main.cpp", "w") as f:
    f.write(text)

