#include <iostream>
#include "glog/logging.h"
#include <time.h>
#include <vector>
#include <fstream>
#include <algorithm>
#include <queue>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <boost/thread/thread.hpp>
#include <pcl/common/common_headers.h>
#include <pcl/features/normal_3d.h> 
#include <pcl/io/pcd_io.h>
#include <pcl/io/obj_io.h>
#include <pcl/io/vtk_lib_io.h>
#include <pcl/io/impl/vtk_lib_io.hpp>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/console/parse.h>
#include <pcl/common/transforms.h>
#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING
#include <experimental/filesystem>
#include "ceres/ceres.h"
#include "class/data_path.h"
#include "class/data_structure.h"
#include "class/visualize.h"
#include "class/reconstruction.h"
#include "class/feature_matching.h"			
#include "class/ranking_system.h"
#include "class/genetic_algorithm.h"

//------------------------------------------------------------------------------------------------------------------------//

//#define NO_RIM_INFO
#define NO_BASE_INFO

using namespace std;
using namespace Eigen;

vector<Geom> shard(SHARD_NUMBER);
vector<Trans> GT_trans(SHARD_NUMBER);
MatrixXd GT_graph(SHARD_NUMBER, SHARD_NUMBER);

pcl::visualization::PCLVisualizer::Ptr viewer(new pcl::visualization::PCLVisualizer("Pot reconstruction"));

VisSwitchVariables vis; 

//------------------------------------------------------------------------------------------------------------------------//

void keyboardEventOccurred(const pcl::visualization::KeyboardEvent& event, void* nothing)
{
	pcl::visualization::PCLVisualizer* viewer = static_cast<pcl::visualization::PCLVisualizer*> (nothing);
	std::string key_string = event.getKeySym();
	bool key_down = event.keyDown();
	vis.KeyEvent(key_string, key_down);
}

//------------------------------------------------------------------------------------------------------------------------//

struct AblationConfig {
    std::string name;
    bool use_inlier_score;
    bool use_connectivity_reward;
    bool use_cycle_penalty;
    bool use_edge_residual;
    bool use_rot_residual;
    bool use_neighbor_penalty;
    bool use_overlap_penalty;
    bool is_icp_refine;
};

std::vector<AblationConfig> ablation_experiments = {
    { "all_components_no_icp",      true,  true,  true,  true,  true,  true,  true,  false },
    { "all_components_icp",         true,  true,  true,  true,  true,  true,  true,  true  },
    { "no_inlier_score_no_icp",     false, true,  true,  true,  true,  true,  true,  false },
    { "no_inlier_score_icp",        false, true,  true,  true,  true,  true,  true,  true  },
    { "no_connectivity_reward_no",  true,  false, true,  true,  true,  true,  true,  false },
    { "no_connectivity_reward_icp", true,  false, true,  true,  true,  true,  true,  true  },
    { "no_cycle_penalty_no_icp",    true,  true,  false, true,  true,  true,  true,  false },
    { "no_cycle_penalty_icp",       true,  true,  false, true,  true,  true,  true,  true  },
    { "no_edge_residual_no_icp",    true,  true,  true,  false, true,  true,  true,  false },
    { "no_edge_residual_icp",       true,  true,  true,  false, true,  true,  true,  true  },
    { "no_rot_residual_no_icp",     true,  true,  true,  true,  false, true,  true,  false },
    { "no_rot_residual_icp",        true,  true,  true,  true,  false, true,  true,  true  },
    { "no_neighbor_penalty_no_icp", true,  true,  true,  true,  true,  false, true,  false },
    { "no_neighbor_penalty_icp",    true,  true,  true,  true,  true,  false, true,  true  },
    { "inlier_score_only_no_icp",   true,  false, false, false, false, false, true,  false },
    { "inlier_score_only_icp",      true,  false, false, false, false, false, true,  true  },
    { "no_penalties_no_icp",        true,  true,  false, false, false, false, true,  false },
    { "no_penalties_icp",           true,  true,  false, false, false, false, true,  true  },
    { "no_overlap_penalty_icp",     true,  true,  true,  true,  true,  true,  false, true  },
};

const int kAblationTrials = 5;
const std::vector<unsigned int> kAblationSeeds = {42, 123, 456, 789, 1024};

// Modern helper for Global ICP refinement in ablation trials
vector<Trans> RunGlobalICP_Refine_Ablation(
    const vector<Geom>& shard_original,
    const MatrixXd& graph_ga,
    const vector<Trans>& T_ga_eval, // Input GT-relative
    const vector<Trans>& T_ga_vis   // Input Axis-assembled relative
) {
    vector<Geom> shard_fine = shard_original;
    vector<Matrix3d> R_fine(SHARD_NUMBER, Matrix3d::Identity());
    vector<Vector3d> t_fine(SHARD_NUMBER, Vector3d::Zero());
    vector<bool> true_node_ga(SHARD_NUMBER, false);

    for (int i = 0; i < SHARD_NUMBER; i++) {
        if (!shard_on_off[i]) continue;
        true_node_ga[i] = true;
        Matrix3d R = Matrix3d::Identity();
        Vector3d t = Vector3d::Zero();
        T_ga_vis[i].Output(R, t);
        shard_fine[i].Move(R, t, true);
    }

    MatrixXd graph_ga_fine = graph_ga;
    IcpFine(shard_fine, R_fine, t_fine, true_node_ga, graph_ga_fine);

    vector<Trans> T_ga_eval_fine = T_ga_eval;
    for (int i = 0; i < SHARD_NUMBER; i++) {
        if (!true_node_ga[i]) continue;
        T_ga_eval_fine[i].Input(R_fine[i], t_fine[i]);
    }
    return T_ga_eval_fine;
}

//------------------------------------------------------------------------------------------------------------------------//

int main(int argc, char** argv) 
{
	//-------------------------------------------------------------------------------------------------------------------//

	double calculation_time(0);
	int s_time(0), e_time(0);

	viewer->setBackgroundColor(0, 0, 0);
	viewer->addCoordinateSystem(1.0);
	viewer->initCameraParameters();

	viewer->registerKeyboardCallback(&keyboardEventOccurred, (void*)viewer.get());

	//------------------------------------------------------------------------------------------------------------------//

	cout << "#################### Pottery Data load ####################" << endl;

	int max_breakline_points(0);

	for (int i = 0; i < SHARD_NUMBER; i++) {
		shard[i].edge_line_.ReadAxis(axis_path[i]);

		if (shard[i].edge_line_.axis_point_.empty()) {
			shard_on_off[i] = false;
			continue;
		}

		// Consider multi-axis shards at the same time
		if(shard_on_off[i])	{
			shard[i].edge_line_.ReadPCDFileWithInfo(file_path[i]);

			if (shard[i].edge_line_.point_.cols() < 50) {
				shard_on_off[i] = false;
				shard[i].edge_line_.Remove();
				continue;
			}

			shard[i].edge_line_.CalculateLineNormal();
			int breakline_points = shard[i].edge_line_.point_.cols();
			max_breakline_points = max(max_breakline_points, breakline_points);
			shard[i].LoadSurface(surface_in[i], surface_out[i], surface_fr[i]);
			shard[i].is_matching_ = true;
			shard[i].sur_frac_.CalculateLineNormal();
		}
	}

	//------------------------------------------------------------------------------------------------------------------//

#ifdef NO_RIM_INFO 
	for (int i = 0; i < SHARD_NUMBER; i++) {
		shard[i].edge_line_.is_seg_rim_ = false;
	}
#endif

#ifdef NO_BASE_INFO
	for (int i = 0; i < SHARD_NUMBER; i++) {
		shard[i].edge_line_.is_seg_base_ = false;
	}
#endif

	//------------------------------------------------------------------------------------------------------------------//

	cout << "#################### Ground Truth data load ####################" << endl;

	for (int i = 0; i < SHARD_NUMBER; i++) {
		for (int j = 0; j < SHARD_NUMBER; j++)
			GT_graph(i, j) = 0;
	}

	int start_index(0);

	for (int i = 0; i < NUM_MIXED_SHERD; i++) {
		MatrixXd single_graph;
		ifstream myfile(gt_graph_path[i]);
		string str;
		vector<string> fileContents;
		stringstream ss;

		while (getline(myfile, str)) {
			fileContents.push_back(str);
		}

		int num_raw = fileContents.size();
		single_graph.resize(num_raw, num_raw);

		for (int j = 0; j < num_raw; j++) {
			ss << fileContents[j];

			for (int k = 0; k < num_raw; k++) {
				ss >> single_graph(j, k);
			}
			ss.clear();
		}

		for (int j = start_index; j < num_raw + start_index; j++) {
			GT_trans[j].Read(gt_T_path[j]);

			for (int k = start_index; k < num_raw + start_index; k++) {
				GT_graph(j, k) = single_graph(j - start_index, k - start_index);
			}
		}

		start_index += num_raw;
	}

	//########## Remove excluded sherd information 
	for (int i = 0; i < SHARD_NUMBER; i++) {
		if (!shard_on_off[i]) {
			for (int j = 0; j < SHARD_NUMBER; j++) {
				GT_graph(i, j) = 0;
				GT_graph(j, i) = 0;
			}
		}
	}
	cout << GT_graph << endl;

	//------------------------------------------------------------------------------------------------------------------//

	cout << "#################### Save initial state ####################" << endl;

	vector<Visualize> pc_origin(SHARD_NUMBER); 

	for (int i = 0; i < SHARD_NUMBER; i++) {
		if (shard_on_off[i]) {
			std::string pointName = "origin_" + std::to_string(i + 1);
			pc_origin[i].MakePointCloud(shard[i].edge_line_.point_, shard[i].edge_line_.normal_, pointName);
			pointName = "o_Mesh" + std::to_string(i + 1);
			pc_origin[i].MakeMesh(obj_path[i], pointName);
		}
	}

	vector<Geom> shard_original = shard;

	//------------------------------------------------------------------------------------------------------------------//

	s_time = clock();

	//------------------------------------------------------------------------------------------------------------------//

	cout << "#################### Change Axis symmetrix to z axis ####################" << endl;

	vector<Trans> T_axis(SHARD_NUMBER);

	for (int i = 0; i < SHARD_NUMBER; i++) {

		if (shard[i].is_matching_) {
			Matrix3d R_d = Matrix3d::Identity();
			Vector3d t_d = { 0, 0, 0 };
			T_axis[i].Set(R_d, t_d, i + 1, i + 1);

			// Align symmetric axis to z-axis
			AxisAlignment(shard[i].edge_line_, R_d, t_d);	
			shard[i].SurMove(R_d, t_d, true);

			pc_origin[i].UpdateData(viewer, shard[i].edge_line_.point_, shard[i].edge_line_.normal_);
			pc_origin[i].AddPointCloud(viewer);
			pc_origin[i].MeshTransform(R_d, t_d, viewer);
			T_axis[i].Input(R_d, t_d);		// Save transformation matrix to z-axis

			CalculateFeatureAxisless(shard[i]);

			//######## Multi axis
			int num_axis = shard[i].edge_line_.axis_norm_.size();

			if (num_axis > 1) {

				for (int j = 1; j < num_axis; j++) {
					Matrix3d R_a, R_i;
					Vector3d t_a, t_i;
					AxisAlignment(shard[i].edge_line_, R_a, t_a, j);
					CalculateFeatureAxisless(shard[i], j);
					R_i = R_a.inverse();
					t_i = -R_i * t_a;
					shard[i].MoveWOSurface(R_i, t_i);
				}
			}
		}
	}

	//------------------------------------------------------------------------------------------------------------------//

	cout << "#################### Feature matching ####################" << endl;

	list<LCSIndex> LCS_out;
	FeatureComp(shard, LCS_out, 25, MINIMUM_NUMBER, 0);
	cout << "Total Number : " << LCS_out.size() << endl;

	//------------------------------------------------------------------------------------------------------------------//

	cout << "#################### Pairwise pruning ####################" << endl;

	PairwisePruning(shard, LCS_out);

	list<LCSIndex> LCS_original = LCS_out;

	//------------------------------------------------------------------------------------------------------------------//

	cout << "#################### Ablation Experiments ####################" << endl;

	for (const AblationConfig& config : ablation_experiments) {
		cout << "=== Ablation experiment: " << config.name << " ===" << endl;

		vector<int> sherd_correct(kAblationTrials, 0);
		vector<int> sherd_total(kAblationTrials, 0);
		vector<int> edge_correct(kAblationTrials, 0);
		vector<int> edge_total(kAblationTrials, 0);
		vector<double> best_fitness(kAblationTrials, 0.0);

		for (int trial = 0; trial < kAblationTrials; ++trial) {
			unsigned int seed = kAblationSeeds[trial];
			srand(seed);

			cout << "  Trial " << trial + 1 << " / " << kAblationTrials << " (seed=" << seed << ")" << endl;

			GeneticAssembler ga(shard, LCS_original, SHARD_NUMBER);
			ga.use_inlier_score        = config.use_inlier_score;
			ga.use_connectivity_reward = config.use_connectivity_reward;
			ga.use_cycle_penalty       = config.use_cycle_penalty;
			ga.use_edge_residual       = config.use_edge_residual;
			ga.use_rot_residual        = config.use_rot_residual;
			ga.use_neighbor_penalty    = config.use_neighbor_penalty;
			ga.use_overlap_penalty     = config.use_overlap_penalty;
			ga.Run();

			vector<Trans> T_ga_vis_abl = ga.GetTransforms();
            vector<Trans> T_ga_eval_abl = T_axis;
            for(int i=0; i<SHARD_NUMBER; i++) {
                if(!shard_on_off[i]) continue;
                Matrix3d R; Vector3d t;
                T_ga_vis_abl[i].Output(R, t);
                T_ga_eval_abl[i].Input(R, t);
            }
            
            if (config.is_icp_refine) {
                T_ga_eval_abl = RunGlobalICP_Refine_Ablation(shard_original, ga.GetGraph(), T_ga_eval_abl, T_ga_vis_abl);
            }

			vector<bool> right_sherd_abl(SHARD_NUMBER, true);
			for (int i = 0; i < SHARD_NUMBER; i++) {
				if (!shard_on_off[i]) right_sherd_abl[i] = false;
			}

			auto [ks, ts, ke, te] = CountResult(GT_graph, GT_trans, ga.GetGraph(), T_ga_eval_abl, right_sherd_abl);

			sherd_correct[trial] = ks;
			sherd_total[trial]   = ts;
			edge_correct[trial]  = ke;
			edge_total[trial]    = te;
			best_fitness[trial]  = ga.GetBestFitness();

			cout << "  Sherd: " << ks << "/" << ts << "  Edge: "  << ke << "/" << te << endl;

            // Save result to disk
            string icp_folder = config.is_icp_refine ? "ICP_True/" : "ICP_False/";
			string trial_path = path + "Result/ablation/" + icp_folder + config.name + "/trial_" + to_string(trial + 1) + "/";
            std::system(("mkdir -p \"" + trial_path + "\"").c_str());

            struct tm* curr_tm;
            time_t curr_time = time(nullptr);
            curr_tm = localtime(&curr_time);
            string com_path = trial_path + to_string(curr_tm->tm_year+1900) + "_" + to_string(curr_tm->tm_mon+1) + "_" + to_string(curr_tm->tm_mday) + "_Ablation_";
            
            for (int i = 0; i < SHARD_NUMBER; i++) {
                if (shard_on_off[i]) {
                    Matrix3d R; Vector3d t;
                    T_ga_eval_abl[i].Output(R, t); // This is Raw->Axis->Assembled[+ICP]
                    BreakLine edge_line_transformed = shard_original[i].edge_line_; // Needs axis transform done separately? 
                    // No, T_ga_eval_abl includes the Axis transform. So we use shard_original (Raw).
                    EdgeLineMove(edge_line_transformed, R, t);
                    
                    Matrix4d T_4d = Matrix4d::Identity();
                    for(int r = 0; r < 3; r++) { T_4d.row(r) << R(r,0), R(r,1), R(r,2), t[r]; }
                    pcl::PolygonMesh save_mesh;
                    pc_origin[i].OutMesh(save_mesh);
                    pcl::PointCloud<pcl::PointNormal> mesh_cloudp;
                    pcl::fromPCLPointCloud2(save_mesh.cloud, mesh_cloudp);
                    pcl::transformPointCloud(mesh_cloudp, mesh_cloudp, T_4d);
                    pcl::toPCLPointCloud2(mesh_cloudp, save_mesh.cloud);

                    SaveEdgeLine(edge_line_transformed, com_path + "Edgeline_" + to_string(i + 1) + ".xyz");
                    pcl::io::saveOBJFile(com_path + "OBJ_" + to_string(i + 1) + ".obj", save_mesh);
                }
            }
		}

		double avg_sherd = 0.0, avg_edge = 0.0, avg_fitness = 0.0;
		for (int t = 0; t < kAblationTrials; t++) {
			avg_sherd += (double)sherd_correct[t] / sherd_total[t];
			avg_edge += (double)edge_correct[t] / edge_total[t];
			avg_fitness += best_fitness[t];
		}
		avg_sherd /= kAblationTrials; avg_edge /= kAblationTrials; avg_fitness /= kAblationTrials;

		cout << "  [" << config.name << "] Avg sherd: " << avg_sherd * 100.0 << "% Avg edge: "  << avg_edge * 100.0 << "%" << endl;
        
        string icp_folder = config.is_icp_refine ? "ICP_True/" : "ICP_False/";
		ofstream summary(path + "Result/ablation/" + icp_folder + config.name + "/summary.txt");
		summary << "Experiment: " << config.name << "\nAvg sherd accuracy: " << avg_sherd*100.0 << "%\nAvg edge accuracy: " << avg_edge*100.0 << "%\n" << endl;
		summary.close();
	}

	cout << "#################### Genetic Algorithm search ####################" << endl;

	const int kMaxGAIterations = 5;
	const double kConvergenceThreshold = 5.0;
	const int kMaxPatience = 5;

	vector<Trans> T_ga;
	MatrixXd graph_ga;
	vector<Trans> T_ga_eval = T_axis;
	vector<Trans> T_ga_vis(SHARD_NUMBER); 
	vector<Trans> T_live(SHARD_NUMBER);
	vector<Trans> T_best(SHARD_NUMBER);
	double prev_best_fitness = -1e9;
	int ga_iteration = 0, patience_counter = 0;

	vector<bool> right_sherd_ga(SHARD_NUMBER, true);
	for (int i = 0; i < SHARD_NUMBER; ++i) {
		if (!shard_on_off[i]) {
			right_sherd_ga[i] = false;
		}
	}


	for (ga_iteration = 0; ga_iteration < kMaxGAIterations; ++ga_iteration) {
		cout << "=== GA Iteration " << ga_iteration + 1 << " / " << kMaxGAIterations << " ===" << endl;
		GeneticAssembler ga_iter(shard, LCS_out, SHARD_NUMBER);
		ga_iter.Run();
		T_ga = ga_iter.GetTransforms();

		for (int i = 0; i < SHARD_NUMBER; ++i) {
			if (!shard_on_off[i]) continue;
			Matrix3d R_d; Vector3d t_d; T_ga[i].Output(R_d, t_d);
			T_live[i].Input(R_d, t_d);
		}

		double current_fitness = ga_iter.GetBestFitness();
		double improvement = current_fitness - prev_best_fitness;

		if (current_fitness >= prev_best_fitness || ga_iteration == 0) {
            graph_ga = ga_iter.GetGraph();
			T_best = T_live;
			for (int i = 0; i < SHARD_NUMBER; ++i) {
				if (!shard_on_off[i]) continue;
				Matrix3d R_b; Vector3d t_b; T_best[i].Output(R_b, t_b);
				T_ga_eval[i] = T_axis[i]; T_ga_eval[i].Input(R_b, t_b);
				T_ga_vis[i] = T_best[i];
			}
		}

		if (ga_iteration > 0) {
			if (improvement < kConvergenceThreshold) {
				if (++patience_counter >= kMaxPatience) break;
			} else patience_counter = 0;
		}
		prev_best_fitness = current_fitness;
		if (ga_iteration == kMaxGAIterations - 1) break;

		shard = shard_original;
		for (int i = 0; i < SHARD_NUMBER; ++i) {
			if (!shard_on_off[i]) continue;
			Matrix3d R; Vector3d t; T_live[i].Output(R, t);
			shard[i].Move(R, t, true);
		}
		LCS_out.clear();
		FeatureComp(shard, LCS_out, 25, MINIMUM_NUMBER, 0);
		PairwisePruning(shard, LCS_out);
	}

	auto [k_sherd, t_sherd, k_edge, t_edge] = CountResult(GT_graph, GT_trans, graph_ga, T_ga_eval, right_sherd_ga);

	//------------------------------------------------------------------------------------------------------------------//
	cout << "#################### Post-GA Global ICP Refinement ####################" << endl;


	vector<Geom> shard_fine = shard_original;
	vector<Matrix3d> R_fine(SHARD_NUMBER, Matrix3d::Identity());
	vector<Vector3d> t_fine(SHARD_NUMBER, Vector3d::Zero());
	vector<bool> true_node_ga(SHARD_NUMBER, false);

	for (int i = 0; i < SHARD_NUMBER; i++) {
		if (!shard_on_off[i]) continue;
		true_node_ga[i] = true;
		Matrix3d R; Vector3d t; T_ga_vis[i].Output(R, t);
		shard_fine[i].Move(R, t, true);
	}

	MatrixXd graph_ga_fine = graph_ga;
	IcpFine(shard_fine, R_fine, t_fine, true_node_ga, graph_ga_fine);

	vector<Trans> T_ga_eval_fine = T_ga_eval;
	vector<Trans> T_ga_vis_fine = T_ga_vis;
	for (int i = 0; i < SHARD_NUMBER; i++) {
		if (!true_node_ga[i]) continue;
		T_ga_eval_fine[i].Input(R_fine[i], t_fine[i]);
		T_ga_vis_fine[i].Input(R_fine[i], t_fine[i]);
	}

	vector<bool> right_sherd_fine = right_sherd_ga;
	auto [k_s_f, t_s_f, k_e_f, t_e_f] = CountResult(GT_graph, GT_trans, graph_ga_fine, T_ga_eval_fine, right_sherd_fine);
		
	cout << "########## Final Refined Results ##########\nFinal Sherd Accuracy : " << k_s_f << " / " << t_s_f << "\nFinal Edge Accuracy  : " << k_e_f << " / " << t_e_f << "\n###########################################" << endl;

	for (int i = 0; i < SHARD_NUMBER; i++) pc_origin[i].TurnOffData(viewer);
	for (int i = 0; i < SHARD_NUMBER; i++) {
		if (!shard_on_off[i]) continue;
		Matrix3d R; Vector3d t; T_ga_vis_fine[i].Output(R, t);
		pc_origin[i].UpdateData(viewer, shard_original[i].edge_line_.point_, shard_original[i].edge_line_.normal_);
		pc_origin[i].MeshTransform(R, t, viewer);
		pc_origin[i].AddPointCloud(viewer); pc_origin[i].AddMesh(viewer);
	}

	viewer->resetCamera();
	cout << "Showing Refined assembly result. Press Q to quit." << endl;
	while (!viewer->wasStopped()) viewer->spinOnce(100);

	return 0;
}