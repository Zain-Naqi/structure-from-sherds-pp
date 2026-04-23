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

	// Save matches computed on original axis-aligned positions
	// (LCS_out gets overwritten during iterative GA, so preserve it here)

	list<LCSIndex> LCS_original = LCS_out;

	//------------------------------------------------------------------------------------------------------------------//

	cout << "#################### Genetic Algorithm search ####################" << endl;

	// Iterative GA Parameters
	const int kMaxGAIterations = 5;
	const double kConvergenceThreshold = 5.0;	// minimum fitness improvement to continue
	const int kMaxPatience = 5;

	vector<Trans> T_ga;
	MatrixXd graph_ga;
	vector<Trans> T_ga_eval = T_axis;
	vector<Trans> T_ga_vis(SHARD_NUMBER); // Tracks overall best GA movement (axis->assembled)
	vector<Trans> T_live(SHARD_NUMBER);   // The currently exploring accumulated transform (axis->assembled)
	vector<Trans> T_best(SHARD_NUMBER);   // The best accumulated transform
	double prev_best_fitness =	 -1e9;
	int ga_iteration = 0;
	int patience_counter = 0;

	// Keep a copy of original axis-aligned shard positions
	// so we can reset between iterations cleanly
	vector<Geom> shard_original = shard;

	for (ga_iteration = 0; ga_iteration < kMaxGAIterations; ++ga_iteration) {

		cout << "=== GA Iteration " << ga_iteration + 1 << " / " << kMaxGAIterations << " ===" << endl;
		
		// Run GA on current match list
		GeneticAssembler ga_iter(shard, LCS_out, SHARD_NUMBER);
		ga_iter.Run(GT_graph, GT_trans, T_axis);
		T_ga = ga_iter.GetTransforms();

		// Accumulate live transformation by composing delta (T_ga) on top of current T_live
		// Trans::Input essentially executes T_live_new = T_delta * T_live_old
		for (int i = 0; i < SHARD_NUMBER; ++i) {
			if (!shard_on_off[i]) continue;
			
			Matrix3d R_delta;
			Vector3d t_delta;
			T_ga[i].Output(R_delta, t_delta);
			T_live[i].Input(R_delta, t_delta);
		}

		// Get best fitness from this run
		double current_fitness = ga_iter.GetBestFitness();
		double improvement = current_fitness - prev_best_fitness;

		// Only accumulate if this iteration improved fitness
		if (current_fitness >= prev_best_fitness || ga_iteration == 0) {
            
            graph_ga = ga_iter.GetGraph();
			T_best = T_live; // Snapshot the best axis-aligned transform

			for (int i = 0; i < SHARD_NUMBER; ++i) {
				if (!shard_on_off[i]) {
					continue;
				}

				Matrix3d R_b;
				Vector3d t_b;
				T_best[i].Output(R_b, t_b);
				
                // T_ga_eval tracks Raw -> Axis -> Assembled
				T_ga_eval[i] = T_axis[i];
				T_ga_eval[i].Input(R_b, t_b);
				T_ga_vis[i] = T_best[i];
			}
		}

		cout << "[GA Iter " << ga_iteration + 1 << "] "
			 << "Best fitness: " << current_fitness
			 << " (improvement: " << improvement << ")" << endl;
			
		// Check convergence
		if (ga_iteration > 0) {
			if (improvement < kConvergenceThreshold) {
				patience_counter++;
				if (patience_counter >= kMaxPatience) {
					cout << "[GA Iter " << ga_iteration + 1 << "] Converged after " << kMaxPatience << " successive non-improving iterations. Stopping." << endl;
					break;
				} else {
					cout << "[GA Iter " << ga_iteration + 1 << "] Did not improve. Patience: " << patience_counter << " / " << kMaxPatience << endl;
				}
			} else {
				patience_counter = 0;
			}
		}

		prev_best_fitness = current_fitness;

		// If this is the last iteration, don't recompute matches
		if (ga_iteration == kMaxGAIterations - 1) {
			break;
		}

		// Apply GA transforms to sherds to get new assembled positions
		// Reset to original positions first, then apply new transforms
		shard = shard_original;

		for (int i = 0; i < SHARD_NUMBER; ++i) {
			if (!shard_on_off[i]) {
				continue;
			}

			Matrix3d R;
			Vector3d t;
			T_live[i].Output(R, t);
			shard[i].Move(R, t, true);
		}

		// Recompute matches on newly assembled positions
		cout << "[GA Iter " << ga_iteration + 1 << "] Recomputing matches on assembled positions..." << endl;

		LCS_out.clear();

		FeatureComp(shard, LCS_out, 25, MINIMUM_NUMBER, 0);

		cout << "[GA Iter " << ga_iteration + 1 << "] New match count: " << LCS_out.size() << endl;

		PairwisePruning(shard, LCS_out);

		cout << "[GA Iter " << ga_iteration + 1 << "] Pruned match count: " << LCS_out.size() << endl;
	}

	cout << "GA converged after " << ga_iteration + 1 << " iteration(s)." << endl;

	pair<int, int> sherd_acc, edge_acc;
	vector<bool> right_sherd_ga(SHARD_NUMBER, true);

	for (int i = 0; i < SHARD_NUMBER; ++i) {
		if (!shard_on_off[i]) {
			right_sherd_ga[i] = false;
		}
	}

	auto [k_sherd, t_sherd, k_edge, t_edge] = CountResult(GT_graph, GT_trans, graph_ga, T_ga_eval, right_sherd_ga);

	//------------------------------------------------------------------------------------------------------------------//
	cout << "#################### Post-GA Global ICP Refinement ####################" << endl;

	// Apply GA transforms to get assembled positions for IcpFine
	vector<Geom> shard_fine = shard_original;
	vector<Matrix3d> R_fine(SHARD_NUMBER, Matrix3d::Identity());
	vector<Vector3d> t_fine(SHARD_NUMBER, Vector3d::Zero());
	vector<bool> true_node_ga(SHARD_NUMBER, false);

	for (int i = 0; i < SHARD_NUMBER; i++) {
		if (!shard_on_off[i]) continue;
		true_node_ga[i] = true;
		Matrix3d R = Matrix3d::Identity();
		Vector3d t = Vector3d::Zero();
		T_ga_vis[i].Output(R, t); // Use T_ga_vis (axis->assembled) to prevent double T_axis application
		shard_fine[i].Move(R, t, true);
		// R_fine and t_fine stay as Identity/Zero
		// IcpFine will output the delta refinement only
	}

	MatrixXd graph_ga_fine = graph_ga;
	IcpFine(shard_fine, R_fine, t_fine, true_node_ga, graph_ga_fine);

	// Compose delta on top of T_ga_eval and T_ga_vis
	vector<Trans> T_ga_eval_fine = T_ga_eval;
	vector<Trans> T_ga_vis_fine = T_ga_vis;
	for (int i = 0; i < SHARD_NUMBER; i++) {
		if (!true_node_ga[i]) continue;
		T_ga_eval_fine[i].Input(R_fine[i], t_fine[i]);
		T_ga_vis_fine[i].Input(R_fine[i], t_fine[i]);
	}

	// Evaluate fine result separately
	vector<bool> right_sherd_fine(SHARD_NUMBER, true);
	for (int i = 0; i < SHARD_NUMBER; i++) {
		if (!shard_on_off[i]) right_sherd_fine[i] = false;
	}
	
	auto [k_sherd_fine, t_sherd_fine, k_edge_fine, t_edge_fine] = CountResult(
		GT_graph, GT_trans, graph_ga_fine, T_ga_eval_fine, right_sherd_fine);
		
	// Final result summary for refined assembly
	sherd_acc = { k_sherd_fine, t_sherd_fine };
	edge_acc = { k_edge_fine, t_edge_fine };
	double time_total = (clock() - s_time) / 1000.0;

	cout << "########## Final Refined Results ##########" << endl;
	cout << "Final Sherd Accuracy : " << k_sherd_fine << " / " << t_sherd_fine << endl;
	cout << "Final Edge Accuracy  : " << k_edge_fine << " / " << t_edge_fine << endl;
	cout << "Total Runtime        : " << time_total << " sec" << endl;
	cout << "###########################################" << endl;

	string path_result_refined = path + "Result/Refined_";
	string refined_mkdir_cmd = "mkdir -p \"" + path_result_refined + "\"";
	std::system(refined_mkdir_cmd.c_str());
	SaveAcc(path_result_refined, sherd_acc, edge_acc, time_total);

	for (int i = 0; i < SHARD_NUMBER; i++) {
		pc_origin[i].TurnOffData(viewer);
	}

	// Apply globally ICP-refined relative transforms and show result.
	for (int i = 0; i < SHARD_NUMBER; i++) {
		if (!shard_on_off[i]) continue;
		Matrix3d R_vis = Matrix3d::Identity();
		Vector3d t_vis = Vector3d::Zero();
		T_ga_vis_fine[i].Output(R_vis, t_vis); // Showcase the globally refined pot exactly from axis-aligned origin
		pc_origin[i].UpdateData(viewer,
			shard_original[i].edge_line_.point_,
			shard_original[i].edge_line_.normal_);
		pc_origin[i].MeshTransform(R_vis, t_vis, viewer);
		pc_origin[i].AddPointCloud(viewer);
		pc_origin[i].AddMesh(viewer);
	}

	viewer->resetCamera();

	//------------------------------------------------------------------------------------------------------------------//

	// Simple viewer loop — press q to quit
	cout << "Showing Refined assembly result. Press Q to quit." << endl;
	while (!viewer->wasStopped()) {
		viewer->spinOnce(100);
	}

	//------------------------------------------------------------------------------------------------------------------//

	return 0;
}