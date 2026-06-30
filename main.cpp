#include <iostream>
#include "glog/logging.h"
#include <chrono>
#include <vector>
#include <fstream>
#include <algorithm>
#include <tuple>
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
	//#################### PCL viewer setting ####################//
	double calculation_time(0);
	double time_ga(0);

	viewer->setBackgroundColor(0, 0, 0);
	viewer->addCoordinateSystem(1.0);
	viewer->initCameraParameters();

	viewer->registerKeyboardCallback(&keyboardEventOccurred, (void*)viewer.get());

	enum class PostGaIcpMode {
		GlobalFine,
		Incremental,
		Compare
	};

	PostGaIcpMode post_icp_mode = PostGaIcpMode::GlobalFine;

	for (int arg_i = 1; arg_i < argc; ++arg_i) {
		string arg = argv[arg_i];
		if (arg == "--icp-global") {
			post_icp_mode = PostGaIcpMode::GlobalFine;
		}
		else if (arg == "--icp-incremental") {
			post_icp_mode = PostGaIcpMode::Incremental;
		}
		else if (arg == "--icp-compare") {
			post_icp_mode = PostGaIcpMode::Compare;
		}
		else if (arg.rfind("--icp-mode=", 0) == 0) {
			string mode = arg.substr(11);
			if (mode == "global") {
				post_icp_mode = PostGaIcpMode::GlobalFine;
			}
			else if (mode == "incremental") {
				post_icp_mode = PostGaIcpMode::Incremental;
			}
			else if (mode == "compare") {
				post_icp_mode = PostGaIcpMode::Compare;
			}
			else {
				cout << "[WARN] Unknown --icp-mode value: " << mode
					 << " (expected global|incremental|compare). Using default global." << endl;
			}
		}
	}

	auto IcpModeName = [](PostGaIcpMode mode) -> string {
		switch (mode) {
		case PostGaIcpMode::GlobalFine:
			return "global";
		case PostGaIcpMode::Incremental:
			return "incremental";
		case PostGaIcpMode::Compare:
			return "compare";
		default:
			return "unknown";
		}
	};

	cout << "[ICP MODE] " << IcpModeName(post_icp_mode)
		 << " (flags: --icp-mode=global|incremental|compare)" << endl;

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

	auto start_time_total = std::chrono::high_resolution_clock::now();

	//------------------------------------------------------------------------------------------------------------------//

	// Keep a copy of original RAW shard positions
	vector<Geom> shard_raw_backup = shard;

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
			pc_origin[i].AddMesh(viewer);
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

	cout << "Pruned match count: " << LCS_out.size() << endl;


	// Save matches computed on original axis-aligned positions
	// (LCS_out gets overwritten during iterative GA, so preserve it here)

	list<LCSIndex> LCS_original = LCS_out;

	//------------------------------------------------------------------------------------------------------------------//

	cout << "#################### Genetic Algorithm search ####################" << endl;

	// Iterative GA Parameters
	int active_shard_count = 0;
	for (int i = 0; i < SHARD_NUMBER; ++i) {
		if (shard_on_off[i]) active_shard_count++;
	}
	const int kMaxGAIterations = std::max(1, active_shard_count);
	const double kConvergenceThreshold = 5.0;	// minimum fitness improvement to continue
	const int kMaxPatience = 5;

	vector<Trans> T_ga;
	MatrixXd graph_ga = MatrixXd::Zero(SHARD_NUMBER, SHARD_NUMBER);
	vector<Trans> T_ga_eval = T_axis;
	vector<Trans> T_ga_vis(SHARD_NUMBER); // Tracks overall best GA movement (axis->assembled)
	vector<Trans> T_live(SHARD_NUMBER);   // The currently exploring accumulated transform (axis->assembled)
	vector<Trans> T_best(SHARD_NUMBER);   // The best accumulated transform
	vector<Matrix3d> R_preview_applied(SHARD_NUMBER, Matrix3d::Identity());
	vector<Vector3d> t_preview_applied(SHARD_NUMBER, Vector3d::Zero());
	double best_fitness_so_far = -1e9;
	int ga_iteration = 0;
	int patience_counter = 0;

	// Keep a copy of original axis-aligned shard positions
	// so we can reset between iterations cleanly
	vector<Geom> shard_original = shard;
	vector<Vector3d> sherd_label_anchor_local(SHARD_NUMBER, Vector3d::Zero());
	vector<bool> sherd_label_valid(SHARD_NUMBER, false);

	for (int i = 0; i < SHARD_NUMBER; ++i) {
		if (!shard_on_off[i]) {
			continue;
		}

		const MatrixXd& pts = shard_original[i].edge_line_.point_;
		if (pts.cols() <= 0) {
			continue;
		}

		Vector3d centroid = Vector3d::Zero();
		for (int c = 0; c < pts.cols(); ++c) {
			centroid += pts.col(c);
		}
		centroid /= static_cast<double>(pts.cols());
		sherd_label_anchor_local[i] = centroid;
		sherd_label_valid[i] = true;
	}

	auto ClearSherdLabels = [&]() {
		for (int i = 0; i < SHARD_NUMBER; ++i) {
			if (!shard_on_off[i]) {
				continue;
			}
			viewer->removeShape("sherd_label_" + to_string(i + 1));
		}
	};

	auto ShowSherdLabels = [&](const vector<Trans>& transforms, double r, double g, double b) {
		ClearSherdLabels();
		const double kLabelScale = 5.0;
		for (int i = 0; i < SHARD_NUMBER; ++i) {
			if (!shard_on_off[i] || !sherd_label_valid[i]) {
				continue;
			}

			Matrix3d R = Matrix3d::Identity();
			Vector3d t = Vector3d::Zero();
			transforms[i].Output(R, t);
			Vector3d p = R * sherd_label_anchor_local[i] + t;
			pcl::PointXYZ label_pt(static_cast<float>(p[0]), static_cast<float>(p[1]), static_cast<float>(p[2]));
			viewer->addText3D(to_string(i + 1),
				label_pt,
				kLabelScale,
				r,
				g,
				b,
				"sherd_label_" + to_string(i + 1));
		}
	};

	auto start_time_ga = std::chrono::high_resolution_clock::now();

	for (ga_iteration = 0; ga_iteration < kMaxGAIterations; ++ga_iteration) {

		cout << "=== GA Iteration " << ga_iteration + 1 << " / " << kMaxGAIterations << " ===" << endl;

		// Run GA on current match list, limiting the number of active edges to the current iteration index + 1
		int target_edges = std::min(active_shard_count - 1, ga_iteration + 1);
		GeneticAssembler ga_iter(shard, LCS_out, SHARD_NUMBER, target_edges);
		ga_iter.Run(GT_graph, GT_trans, T_axis);
		T_ga = ga_iter.GetTransforms();

		// Update live transformation directly to the absolute pose T_ga
		for (int i = 0; i < SHARD_NUMBER; ++i) {
			if (!shard_on_off[i]) continue;

			Matrix3d R_delta;
			Vector3d t_delta;
			T_ga[i].Output(R_delta, t_delta);
			T_live[i].Set(R_delta, t_delta);
		}

		// Get best fitness from this run
		double current_fitness = ga_iter.GetBestFitness();
		double improvement = current_fitness - best_fitness_so_far;
		bool improved = (current_fitness > best_fitness_so_far);

		// Only keep state if it strictly improves the global best fitness.
		if (improved) {
            best_fitness_so_far = current_fitness;

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

		// // Interactive preview: visualize the current iteration's best individual,
		// // then wait for space key before moving to the next GA iteration.
		// for (int i = 0; i < SHARD_NUMBER; ++i) {
		// 	if (!shard_on_off[i]) {
		// 		continue;
		// 	}

		// 	Matrix3d R_target = Matrix3d::Identity();
		// 	Vector3d t_target = Vector3d::Zero();
		// 	T_live[i].Output(R_target, t_target);

		// 	Matrix3d R_delta = R_target * R_preview_applied[i].transpose();
		// 	Vector3d t_delta = t_target - R_delta * t_preview_applied[i];
		// 	pc_origin[i].Transform(R_delta, t_delta, viewer);
		// 	R_preview_applied[i] = R_target;
		// 	t_preview_applied[i] = t_target;
		// }

		// ShowSherdLabels(T_live, 1.0, 1.0, 0.2);

		// cout << "[GA Iter " << ga_iteration + 1
		// 	 << "] Previewing current best individual. Press SPACE to continue." << endl;
		// vis.first_ = false;
		// while (!viewer->wasStopped() && !vis.first_) {
		// 	viewer->spinOnce(50);
		// }
		// vis.first_ = false;
		// if (viewer->wasStopped()) {
		// 	return 0;
		// }

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

		// If this is the last iteration, don't recompute matches
		if (ga_iteration == kMaxGAIterations - 1) {
			break;
		}

		// Recompute matches from best-so-far pose to avoid drifting into
		// weaker candidate sets after non-improving exploratory moves.
		// 1. Move active shards to their current best assembled poses
		shard = shard_original;
		for (int i = 0; i < SHARD_NUMBER; ++i) {
			if (shard_on_off[i]) {
				Matrix3d R;
				Vector3d t;
				T_best[i].Output(R, t);
				shard[i].Move(R, t, true);
			}
		}

		// 2. Identify active shards in the largest connected component of graph_ga
		vector<int> largest_component;
		{
			vector<vector<int>> components;
			vector<bool> visited(SHARD_NUMBER, false);
			for (int i = 0; i < SHARD_NUMBER; ++i) {
				if (!shard_on_off[i] || visited[i]) continue;
				vector<int> comp;
				queue<int> q;
				q.push(i);
				visited[i] = true;
				while (!q.empty()) {
					int curr = q.front();
					q.pop();
					comp.push_back(curr);
					for (int next = 0; next < SHARD_NUMBER; ++next) {
						if (shard_on_off[next] && !visited[next] && (graph_ga(curr, next) > 0 || graph_ga(next, curr) > 0)) {
							visited[next] = true;
							q.push(next);
						}
					}
				}
				components.push_back(comp);
			}
			int max_size = 0;
			for (const auto& comp : components) {
				if (comp.size() > max_size) {
					max_size = comp.size();
					largest_component = comp;
				}
			}
		}

		vector<Geom*> active_geom_ptrs;
		for (int idx : largest_component) {
			active_geom_ptrs.push_back(&shard[idx]);
		}
		int root_idx = largest_component.empty() ? -1 : largest_component[0];

		if (root_idx != -1 && active_geom_ptrs.size() > 1) {
			cout << "[GA Iter " << ga_iteration + 1 << "] Running joint axis refinement on "
				 << active_geom_ptrs.size() << " shards in the largest component..." << endl;

			// Initialize joint cylinder axis direction and point to root shard's axis in Z-aligned assembly
			Vector3d axis_point = { 0.0, 0.0, 0.0 };
			Vector3d axis_normal = { 0.0, 0.0, 1.0 };

			// Run the Ceres axis optimization solver
			RefineAxis(active_geom_ptrs, axis_point, axis_normal, 100, NUMBER_OF_THREAD, 0.5, true);

			// 3. Align this refined joint cylinder axis to the global Z-axis
			shard[root_idx].edge_line_.axis_point_[0] = axis_point;
			shard[root_idx].edge_line_.axis_norm_[0] = axis_normal;

			Matrix3d R_a;
			Vector3d t_a;
			AxisAlignment(shard[root_idx].edge_line_, R_a, t_a, 0);
			shard[root_idx].SurMove(R_a, t_a, true);

			for (int idx : largest_component) {
				if (idx != root_idx) {
					shard[idx].Move(R_a, t_a, true);
				}
			}

			// 4. Update axis-aligning transforms and recalculate local centered descriptors for component shards
			for (int i : largest_component) {
				// Accumulated transform from raw -> refined Z-aligned assembly
				// T_assembly_i = T_align * T_best_i * T_axis_i
				Matrix3d R_best_old;
				Vector3d t_best_old;
				T_best[i].Output(R_best_old, t_best_old);

				Matrix3d R_axis_old;
				Vector3d t_axis_old;
				T_axis[i].Output(R_axis_old, t_axis_old);

				Matrix3d R_combined = R_best_old * R_axis_old;
				Vector3d t_combined = R_best_old * t_axis_old + t_best_old;

				Matrix3d R_assembly = R_a * R_combined;
				Vector3d t_assembly = R_a * t_combined + t_a;

				// Extrapolate the Z-axis of the assembly back to raw coordinates to get the refined axis
				Vector3d v_refined = R_assembly.transpose() * Vector3d(0, 0, 1);
				Vector3d p_refined = R_assembly.transpose() * (-t_assembly);

				// Prepare the raw shard to be aligned to this refined axis
				shard_original[i] = shard_raw_backup[i];
				shard_original[i].edge_line_.axis_norm_[0] = v_refined;
				shard_original[i].edge_line_.axis_point_[0] = p_refined;

				Matrix3d R_align;
				Vector3d t_align;
				AxisAlignment(shard_original[i].edge_line_, R_align, t_align, 0);
				shard_original[i].SurMove(R_align, t_align, true);

				// Update T_axis[i] with the new raw -> axis alignment transform
				T_axis[i].Set(R_align, t_align, i + 1, i + 1);

				// Recalculate 1D descriptors on the refined Z-aligned local pose
				CalculateFeatureAxisless(shard_original[i], 0);

				// Compute the new T_best[i] relative to the new axis-aligned pose:
				// T_best_new_i = T_assembly_i * T_axis_new_i.inverse()
				Matrix3d R_axis_new = R_align;
				Vector3d t_axis_new = t_align;
				Matrix3d R_axis_new_inv = R_axis_new.transpose();
				Vector3d t_axis_new_inv = -R_axis_new_inv * t_axis_new;

				Matrix3d R_best_new = R_assembly * R_axis_new_inv;
				Vector3d t_best_new = R_assembly * t_axis_new_inv + t_assembly;

				T_best[i].Set(R_best_new, t_best_new, i + 1, i + 1);
			}

			// Update T_live to reflect the new axis-aligned starting state for the next GA run
			T_live = T_best;
		}

		// Reset shard to the updated axis-aligned positions
		shard = shard_original;

		// Recompute matches on newly assembled positions
		cout << "[GA Iter " << ga_iteration + 1 << "] Recomputing matches on assembled positions..." << endl;

		LCS_out.clear();

		FeatureComp(shard, LCS_out, 25, MINIMUM_NUMBER, 0);

		cout << "[GA Iter " << ga_iteration + 1 << "] New match count: " << LCS_out.size() << endl;

		PairwisePruning(shard, LCS_out);

		cout << "[GA Iter " << ga_iteration + 1 << "] Pruned match count: " << LCS_out.size() << endl;
	}

	auto end_time_ga = std::chrono::high_resolution_clock::now();
	time_ga = std::chrono::duration<double>(end_time_ga - start_time_ga).count();

	cout << "GA converged after " << ga_iteration + 1 << " iteration(s)." << endl;

	pair<int, int> sherd_acc, edge_acc;
	vector<bool> right_sherd_ga(SHARD_NUMBER, true);

	for (int i = 0; i < SHARD_NUMBER; ++i) {
		if (!shard_on_off[i]) {
			right_sherd_ga[i] = false;
		}
	}

	auto [k_sherd, t_sherd, k_edge, t_edge] = CountResult(
		GT_graph, GT_trans, graph_ga, T_ga_eval, right_sherd_ga);

	auto PrintStage = [](const string& label,
		int k_sherd_stage,
		int t_sherd_stage,
		int k_edge_stage,
		int t_edge_stage)
		{
			cout << "########## " << label << " ##########" << endl;
			cout << "Sherd Accuracy : " << k_sherd_stage << " / " << t_sherd_stage << endl;
			cout << "Edge Accuracy  : " << k_edge_stage << " / " << t_edge_stage << endl;
			cout << "########################################" << endl;
		};

	PrintStage("Pre-ICP GA Results", k_sherd, t_sherd, k_edge, t_edge);

	int k_sherd_global = k_sherd, t_sherd_global = t_sherd;
	int k_edge_global = k_edge, t_edge_global = t_edge;
	vector<Trans> T_ga_eval_global = T_ga_eval;
	vector<Trans> T_ga_vis_global = T_ga_vis;
	MatrixXd graph_ga_global = graph_ga;
	bool has_global_result = false;

	int k_sherd_incremental = k_sherd, t_sherd_incremental = t_sherd;
	int k_edge_incremental = k_edge, t_edge_incremental = t_edge;
	vector<Trans> T_ga_eval_incremental = T_ga_eval;
	vector<Trans> T_ga_vis_incremental = T_ga_vis;
	MatrixXd graph_ga_incremental = graph_ga;
	bool has_incremental_result = false;

	//------------------------------------------------------------------------------------------------------------------//
	if (post_icp_mode == PostGaIcpMode::GlobalFine || post_icp_mode == PostGaIcpMode::Compare) {
		cout << "#################### Post-GA Global ICP Refinement ####################" << endl;

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

		graph_ga_global = graph_ga;
		IcpFine(shard_fine, R_fine, t_fine, true_node_ga, graph_ga_global);

		T_ga_eval_global = T_ga_eval;
		T_ga_vis_global = T_ga_vis;
		for (int i = 0; i < SHARD_NUMBER; i++) {
			if (!true_node_ga[i]) continue;
			T_ga_eval_global[i].Input(R_fine[i], t_fine[i]);
			T_ga_vis_global[i].Input(R_fine[i], t_fine[i]);
		}

		vector<bool> right_sherd_fine(SHARD_NUMBER, true);
		for (int i = 0; i < SHARD_NUMBER; i++) {
			if (!shard_on_off[i]) right_sherd_fine[i] = false;
		}

		tie(k_sherd_global, t_sherd_global, k_edge_global, t_edge_global) = CountResult(
			GT_graph, GT_trans, graph_ga_global, T_ga_eval_global, right_sherd_fine);

		has_global_result = true;
		PrintStage("Post-ICP Global Results", k_sherd_global, t_sherd_global, k_edge_global, t_edge_global);
	}

	if (post_icp_mode == PostGaIcpMode::Incremental || post_icp_mode == PostGaIcpMode::Compare) {
		cout << "#################### Post-GA Incremental ICP Refinement ####################" << endl;

		int root_idx = -1;
		for (int i = 0; i < SHARD_NUMBER; ++i) {
			if (shard_on_off[i]) {
				root_idx = i;
				break;
			}
		}

		if (root_idx >= 0) {
			RankingSubgraph graph_icp(LCS_original, SHARD_NUMBER);
			graph_icp.node_[root_idx] = true;
			graph_icp.root_node_ = root_idx + 1;
			graph_icp.ResetMatchedIndex(shard_original);

			vector<Geom> shard_icp = shard_original;

			vector<int> placement_order;
			vector<bool> bfs_placed(SHARD_NUMBER, false);
			queue<int> bfs_queue;
			bfs_placed[root_idx] = true;
			placement_order.push_back(root_idx);
			bfs_queue.push(root_idx);

			while (!bfs_queue.empty()) {
				int curr = bfs_queue.front();
				bfs_queue.pop();

				for (int j = 0; j < SHARD_NUMBER; ++j) {
					if (bfs_placed[j] || !shard_on_off[j]) {
						continue;
					}
					if (graph_ga(curr, j) > 0 || graph_ga(j, curr) > 0) {
						bfs_placed[j] = true;
						placement_order.push_back(j);
						bfs_queue.push(j);
					}
				}
			}

			for (int pi = 1; pi < static_cast<int>(placement_order.size()); ++pi) {
				int current = placement_order[pi];

				shard_icp = shard_original;

				for (int k = 0; k < SHARD_NUMBER; ++k) {
					if (!graph_icp.node_[k]) continue;
					Matrix3d R_k = Matrix3d::Identity();
					Vector3d t_k = Vector3d::Zero();
					graph_icp.T_[k].Output(R_k, t_k);
					shard_icp[k].Move(R_k, t_k);
				}

				vector<RankingSubgraph> single_graph_vec;
				single_graph_vec.push_back(graph_icp);
				single_graph_vec[0].MakeHierarchyPriorityList(0, single_graph_vec);
				graph_icp = single_graph_vec[0];

				vector<Chunk> filtered_priority;
				for (auto& chunk : graph_icp.priority_list_) {
					bool has_root_edge = false;
					for (int edge_index : chunk.i_edge) {
						int sx = graph_icp.sub_graph_[edge_index].shard_x_;
						int sy = graph_icp.sub_graph_[edge_index].shard_y_;
						if (sx == (root_idx + 1) || sy == (root_idx + 1)) {
							has_root_edge = true;
							break;
						}
					}
					if (has_root_edge) {
						filtered_priority.push_back(chunk);
					}
				}

				if (!filtered_priority.empty()) {
					graph_icp.priority_list_ = filtered_priority;
				}

				if (graph_icp.priority_list_.empty()) {
					graph_icp.node_[current] = true;
					continue;
				}

				bool found = false;
				for (int ci = 0; ci < static_cast<int>(graph_icp.priority_list_.size()); ++ci) {
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

				vector<TransHistory> history;
				PrepareGraphBuilding(shard_icp, graph_icp, history, 0);

				RankingSubgraph pregraph_step(SHARD_NUMBER);
				for (int k = 0; k < SHARD_NUMBER; ++k) {
					if (k == current) continue;
					if (graph_icp.node_[k]) {
						pregraph_step.node_[k] = true;
					}
				}

				vector<RankingSubgraph> pregraph_vec = { pregraph_step };

				vector<Matrix3d> R_step(SHARD_NUMBER, Matrix3d::Identity());
				vector<Vector3d> t_step(SHARD_NUMBER, Vector3d::Zero());
				int inlier_step = 0;
				IcpIncGraphAxis(shard_icp,
					R_step,
					t_step,
					graph_icp,
					pregraph_vec,
					inlier_step,
					true,
					true);

				for (int k = 0; k < SHARD_NUMBER; ++k) {
					if (!graph_icp.node_[k]) continue;
					graph_icp.T_[k].Input(R_step[k], t_step[k]);
				}
			}

			T_ga_eval_incremental = T_axis;
			T_ga_vis_incremental = graph_icp.T_;
			for (int k = 0; k < SHARD_NUMBER; ++k) {
				if (!shard_on_off[k] || k == root_idx) continue;
				Matrix3d R_k = Matrix3d::Identity();
				Vector3d t_k = Vector3d::Zero();
				graph_icp.T_[k].Output(R_k, t_k);
				T_ga_eval_incremental[k].Input(R_k, t_k);
			}

			vector<bool> right_sherd_icp(SHARD_NUMBER, true);
			for (int k = 0; k < SHARD_NUMBER; ++k) {
				if (!shard_on_off[k]) right_sherd_icp[k] = false;
			}

			graph_ga_incremental = graph_ga;
			tie(k_sherd_incremental, t_sherd_incremental, k_edge_incremental, t_edge_incremental) = CountResult(
				GT_graph, GT_trans, graph_ga_incremental, T_ga_eval_incremental, right_sherd_icp);

			has_incremental_result = true;
			PrintStage("Post-ICP Incremental Results",
				k_sherd_incremental,
				t_sherd_incremental,
				k_edge_incremental,
				t_edge_incremental);
		}
		else {
			cout << "[WARN] Incremental ICP skipped: no active shard found." << endl;
		}
	}

	vector<Trans> T_final_vis = T_ga_vis;
	string final_stage_label = "Pre-ICP GA";
	int k_sherd_final = k_sherd, t_sherd_final = t_sherd;
	int k_edge_final = k_edge, t_edge_final = t_edge;

	auto IsResultBetter = [](int sherd_a, int edge_a, int sherd_b, int edge_b) {
		if (sherd_a != sherd_b) {
			return sherd_a > sherd_b;
		}
		return edge_a > edge_b;
	};

	if (post_icp_mode == PostGaIcpMode::GlobalFine && has_global_result) {
		final_stage_label = "Post-ICP Global";
		k_sherd_final = k_sherd_global;
		t_sherd_final = t_sherd_global;
		k_edge_final = k_edge_global;
		t_edge_final = t_edge_global;
		T_final_vis = T_ga_vis_global;
	}
	else if (post_icp_mode == PostGaIcpMode::Incremental && has_incremental_result) {
		final_stage_label = "Post-ICP Incremental";
		k_sherd_final = k_sherd_incremental;
		t_sherd_final = t_sherd_incremental;
		k_edge_final = k_edge_incremental;
		t_edge_final = t_edge_incremental;
		T_final_vis = T_ga_vis_incremental;
	}
	else if (post_icp_mode == PostGaIcpMode::Compare) {
		bool selected_non_pre = false;

		if (has_global_result &&
			IsResultBetter(k_sherd_global, k_edge_global, k_sherd_final, k_edge_final)) {
			final_stage_label = "Post-ICP Global (selected in compare)";
			k_sherd_final = k_sherd_global;
			t_sherd_final = t_sherd_global;
			k_edge_final = k_edge_global;
			t_edge_final = t_edge_global;
			T_final_vis = T_ga_vis_global;
			selected_non_pre = true;
		}

		if (has_incremental_result &&
			IsResultBetter(k_sherd_incremental, k_edge_incremental, k_sherd_final, k_edge_final)) {
			final_stage_label = "Post-ICP Incremental (selected in compare)";
			k_sherd_final = k_sherd_incremental;
			t_sherd_final = t_sherd_incremental;
			k_edge_final = k_edge_incremental;
			t_edge_final = t_edge_incremental;
			T_final_vis = T_ga_vis_incremental;
			selected_non_pre = true;
		}

		if (!selected_non_pre) {
			final_stage_label = "Pre-ICP GA (selected in compare)";
		}
	}

	T_final_vis = T_ga_vis; // FORCE PRE-ICP GA POSES FOR DEBUGGING

	// Final result summary for refined assembly
	sherd_acc = { k_sherd_final, t_sherd_final };
	edge_acc = { k_edge_final, t_edge_final };
	auto end_time_total = std::chrono::high_resolution_clock::now();
	double time_total = std::chrono::duration<double>(end_time_total - start_time_total).count();

	cout << "########## Final Refined Results ##########" << endl;
	cout << "Selected Stage      : " << final_stage_label << endl;
	cout << "Final Sherd Accuracy : " << k_sherd_final << " / " << t_sherd_final << endl;
	cout << "Final Edge Accuracy  : " << k_edge_final << " / " << t_edge_final << endl;
	cout << "Total Runtime        : " << time_total << " sec" << endl;
	cout << "GA Runtime           : " << time_ga << " sec" << endl;
	cout << "###########################################" << endl;

	string path_result_refined = path + "Result/Refined_";
	string refined_mkdir_cmd = "mkdir -p \"" + path_result_refined + "\"";
	std::system(refined_mkdir_cmd.c_str());
	SaveAcc(path_result_refined, sherd_acc, edge_acc, time_total);

	for (int i = 0; i < SHARD_NUMBER; i++) {
		pc_origin[i].TurnOffData(viewer);
	}
	ClearSherdLabels();

	// Build Ground Truth assembly for comparison
	vector<Visualize> pc_gt(SHARD_NUMBER);
	for (int i = 0; i < SHARD_NUMBER; i++) {
		if (!shard_on_off[i]) continue;
		string pointName = "gt_origin_" + to_string(i + 1);
		pc_gt[i].MakePointCloud(shard_original[i].edge_line_.point_, shard_original[i].edge_line_.normal_, pointName);
		pointName = "gt_Mesh" + to_string(i + 1);
		pc_gt[i].MakeMesh(obj_path[i], pointName);

		Matrix3d R_gt;
		Vector3d t_gt;
		GT_trans[i].Output(R_gt, t_gt);
		pc_gt[i].MeshTransform(R_gt, t_gt, viewer);
		pc_gt[i].TurnOffData(viewer);
	}

	// Apply selected final transforms and show initial result.
	for (int i = 0; i < SHARD_NUMBER; i++) {
		if (!shard_on_off[i]) continue;
		Matrix3d R_vis = Matrix3d::Identity();
		Vector3d t_vis = Vector3d::Zero();
		T_final_vis[i].Output(R_vis, t_vis);
		pc_origin[i].UpdateData(viewer,
			shard_original[i].edge_line_.point_,
			shard_original[i].edge_line_.normal_);
		pc_origin[i].Transform(R_vis, t_vis, viewer);
		pc_origin[i].AddPointCloud(viewer);
		pc_origin[i].AddMesh(viewer);
	}
	ShowSherdLabels(T_final_vis, 1.0, 1.0, 0.2);

	viewer->resetCamera();

	// Interactive toggle loop — press 'G' to switch between GA and GT
	bool last_ground = false; // Monitor toggle state
	cout << "Showing Refined assembly. Press 'G' to toggle Ground Truth comparison! Press 'Q' to quit." << endl;

    while (!viewer->wasStopped()) {
		if (vis.ground_ != last_ground) {
			for (int i = 0; i < SHARD_NUMBER; i++) {
				if (!shard_on_off[i]) continue;
				if (vis.ground_) {
					// SHOW GROUND TRUTH (Green tint)
					pc_origin[i].TurnOffData(viewer);
					pc_gt[i].AddMesh(viewer, 0.0, 1.0, 0.5); // Green-ish Cyan
				} else {
					// SHOW GA ASSEMBLY (Default)
					pc_gt[i].TurnOffData(viewer);
					pc_origin[i].AddMesh(viewer, 0.8, 0.8, 0.8); // Standard Grey
				}
			}

			if (vis.ground_) {
				ShowSherdLabels(GT_trans, 0.0, 1.0, 0.5);
			}
			else {
				ShowSherdLabels(T_final_vis, 1.0, 1.0, 0.2);
			}

			last_ground = vis.ground_;
		}
		viewer->spinOnce(100);
	}

	//------------------------------------------------------------------------------------------------------------------//

	return 0;
}
