////////////////////////////////////////////////////////////////////////////////
#include "State.h"
#include <cellogram/convex_hull.h>
#include <cellogram/extrude_mesh.h>
#include <cellogram/image_reader.h>
#include <cellogram/laplace_energy.h>
#include <cellogram/MeshUtils.h>
#include <cellogram/point_source_detection.h>
#include <cellogram/PolygonUtils.h>
#include <cellogram/region_grow.h>
#include <cellogram/remesh_adaptive.h>
#include <cellogram/tri2hex.h>
#include <cellogram/voronoi.h>
#include <polyfem/MeshUtils.hpp>
#include <polyfem/InterpolatedFunction.hpp>
#include <igl/bounding_box.h>
#include <igl/bounding_box_diagonal.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/list_to_matrix.h>
#include <igl/point_in_poly.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/slice.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/enumerable_thread_specific.h>
#include <fstream>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

	// -----------------------------------------------------------------------------

	namespace {
		json default_detection_settings = R"({
			"energy_variation_from_mean": 2.0,
			"lloyd_iterations": 20,
			"perm_possibilities": 15,
			"sigma": 2.0,
			"otsu_multiplier": 1.0
     		})"_json;

		json default_analysis_settings = R"({
			"scaling": 0.2,
			"E": 13.578,
			"I": 0.5,
			"L": 3.0,
			"eps": 0.32967033982276917,
			"formulation": "LinearElasticity",
			"image_from_pillars": false,
			"nu": 0.49,
			"padding_size": 25.0,
			"power": 2.0,
			"uniform_mesh_size": 0.5,
			"adaptive_mesh_size": [
				0.2,
				1.0
			],
			"thickness": 30.0
     		})"_json;
	}


	State::State() {
	}


	State &State::state() {
		static State instance;
		return instance;
	}



	int State::find_region_by_interior_vertex(const int index)
	{
		for(int i = 0; i < regions.size(); i++)
		{
			for (int j = 0; j < regions[i].region_interior.rows(); j++)
			{
				if (regions[i].region_interior(j) == index)
				{
					int region_index = i;
					return region_index;
				}
			}
		}
		return -1;
	}

	void State::find_region_by_boundary_vertex(const int index, std::vector<int> & found_regions)
	{
		for (int i = 0; i < regions.size(); i++)
		{
			for (int j = 0; j < regions[i].region_boundary.rows(); j++)
			{
				if (regions[i].region_boundary(j) == index)
				{
					found_regions.push_back(i);
				}
			}
		}
	}

	void State::split_region(const Eigen::Vector2i & split_end_points)
	{
		// find region that contains for points in boundary and call splitting method of region
		std::vector<int> regions1, regions2;
		find_region_by_boundary_vertex(split_end_points(0), regions1);
		find_region_by_boundary_vertex(split_end_points(1), regions2);

		if (regions1.empty() || regions2.empty())
			return;

		int index;
		for (int i = 0; i < regions1.size(); i++)
		{
			if (std::find(regions2.begin(), regions2.end(), regions1[i]) != regions2.end())
			{
				index = regions1[i];
				break;
			}
		}

		std::cout << "splitting region: " << index << std::endl;
		Region r1, r2;
		regions[index].split_region(mesh, split_end_points,r1,r2);

		regions.erase(regions.begin() + index);
		if (r1.size() > 0)
			regions.push_back(r1);
		if (r2.size() > 0)
			regions.push_back(r2);

		for(int i = 0; i < regions.size();++i)
		{
			regions[i].find_triangles(mesh.triangles);
		}
	}

	// -----------------------------------------------------------------------------
	void State::load_detection_settings(json args) {
		json settings = default_detection_settings;
		settings.merge_patch(args);
		lloyd_iterations = settings["lloyd_iterations"];
		energy_variation_from_mean = settings["energy_variation_from_mean"];
		perm_possibilities = settings["perm_possibilities"];
		sigma = settings["sigma"];
		otsu_multiplier = settings["otsu_multiplier"];
	}
	void State::load_analysis_settings(json args) {
		json settings = default_analysis_settings;
		settings.merge_patch(args);
		std::vector<float> tmp = settings["adaptive_mesh_size"];
		adaptive_mesh_size[0] = tmp[0];
		adaptive_mesh_size[1] = tmp[1];
		power = settings["power"];
		scaling = settings["scaling"];
		padding_size = settings["padding_size"];
		thickness = settings["thickness"];
		uniform_mesh_size = settings["uniform_mesh_size"];
		E = settings["E"];
		nu = settings["nu"];
		eps = settings["eps"];
		I = settings["I"];
		L = settings["L"];
		formulation = settings["formulation"].get<std::string>();
		image_from_pillars = settings["image_from_pillars"];
	}

	void State::load_settings(json args) {
		load_detection_settings(args.value("settings", json::object()));
		load_analysis_settings(args.value("analysis_settings", json::object()));
	}

	void State::load_settings(const std::string &path) {
		json settings = json::object();
		if (!path.empty()) {
			std::ifstream in(path);
			in >> settings;
		}
		load_settings(settings);
	}

	bool State::load(const std::string & path)
	{
		using json = nlohmann::json;

		bool ok = false;
		// clear previous
		hull_vertices.resize(0, 0); //needed for lloyd
		hull_faces.resize(0, 0);
		hull_polygon.resize(0, 0);
		regions.clear();



		// load settings

		std::ifstream json_in(path + "/all.json");

		json unique;
		json_in >> unique;
		load_settings(unique);
		/*if (!unique["settings"].empty())
		{
			json settings = unique["settings"];
			load_settings(settings);
		}*/
		// Load points
		// ok = mesh.load(path);
		ok = mesh.load(unique["mesh"]);

		if (!unique["analysis"].empty() && !image_from_pillars)
			ok = mesh3d.load(unique["analysis"]);

		//if (!unique["analysis_setting"].empty())
		//{
		//	json analysis_settings = unique["analysis_settings"];
		//	load_analysis_settings(analysis_settings);
		//}

		compute_hull();

		return ok;
	}

	bool State::is_data_available(const std::string &path)
	{
// #ifdef WIN32
// 		std::string save_data = path + "\\cellogram\\moved.vert";
// #else
// 		std::string save_data = path + "/cellogram/moved.vert";
// #endif

#ifdef WIN32
		std::string save_data = path + "\\all.json";
#else
		std::string save_data = path + "/all.json";
#endif
		std::ifstream f(save_data.c_str());
		bool ok = f.good();
		return f.good();
	}

	bool State::load_image(const std::string fname)
	{
		bool ok = read_image(fname, img);

		//std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << img << std::endl;

		hull_vertices.resize(0, 0); //needed for lloyd
		hull_faces.resize(0, 0);
		hull_polygon.resize(0, 0);
		regions.clear();
		mesh.clear();
		mesh3d.clear();

		return ok;
	}

	//bool State::load_param(const std::string & path)
	//{
	//	return mesh.load_params(path);
	//}

	bool State::save(const std::string & path, const bool full_path)
	{
		using json = nlohmann::json;

		json unique;

		json json_data;
		json_data["lloyd_iterations"] = lloyd_iterations;
		json_data["energy_variation_from_mean"] = energy_variation_from_mean;
		json_data["perm_possibilities"] = perm_possibilities;
		json_data["sigma"] = sigma;
		json_data["otsu_multiplier"] = otsu_multiplier;

		unique["settings"] = json_data;

		unique["mesh"] = json::object();
		mesh.save(unique["mesh"]);

		unique["analysis"] = json::object();
		if (mesh3d.empty() == false)
		{
			mesh3d.save_mesh(unique["analysis"]);
		}
		if (mesh3d.analysed())
		{
			mesh3d.save_traction(unique["analysis"]);

			json_data.clear();
			json_data["formulation"] = formulation;
			json_data["image_from_pillars"] = image_from_pillars;
			json_data["adaptive_mesh_size"] = adaptive_mesh_size;
			json_data["power"] = power;
			json_data["padding_size"] = padding_size;
			json_data["thickness"] = thickness;
			json_data["uniform_mesh_size"] = uniform_mesh_size;
			json_data["E"] = E;
			json_data["nu"] = nu;
			json_data["eps"] = eps;
			json_data["I"] = I;
			json_data["L"] = L;
			json_data["scaling"] = scaling;

			unique["analysis_settings"] = json_data;
		}

		json hull;
		hull["vertices"] = json::object();
		hull["triangles"] = json::object();
		write_json_mat(hull_vertices, hull["vertices"]);
		write_json_mat(hull_faces, hull["triangles"]);

		unique["hull"] = hull;

		{
			if (full_path)
			{
				// todo: check if .json is already at end of string
				std::ofstream json_out(path + ".json");
				json_out << unique.dump(4) << std::endl;
				json_out.close();
			}
			else
			{
				std::ofstream json_out(path + "/all.json");
				json_out << unique.dump(4) << std::endl;
				json_out.close();
			}
		}

		return true;
	}

	void State::compute_hull()
	{
		//convex_hull(points, boundary);
		if (mesh.moved.rows() < 3)
			return;

		loose_convex_hull(mesh.moved, mesh.boundary, 6);
		int dims = (int)mesh.moved.cols();

		// Compute polygon of the convex hull
		Eigen::VectorXi I = mesh.boundary;
		Eigen::VectorXi J = Eigen::VectorXi::LinSpaced(dims, 0, dims - 1);

		igl::slice(mesh.moved, I, J, hull_polygon);

		// Offset by epsilon
		// offset_polygon(hull_polygon, hull_polygon, 1);

		// Draw filled polygon
		//triangulate_convex_polygon(hull_polygon, hull_vertices, hull_faces);
		triangulate_polygon(hull_polygon, hull_vertices, hull_faces);
	}

	void State::clean_hull()
	{
		// Find and delete the vertices outside of the hull
		//std::vector<int> ind;
		//for (int i = 0; i < points.rows(); i++)
		//{
		//	//bool is_on_hull = false;
		//	//for (int j = 0; j < hull_index.rows(); j++)
		//	//{
		//	//	if (i == hull_index(j))
		//	//	{
		//	//		is_on_hull = true;
		//	//	}
		//	//}
		//	//if (!is_on_hull)
		//	//{
		//		if (!is_inside(hull_vertices, mesh.detected(i, 0), mesh.detected(i, 1)))
		//		{
		//			std::cout << mesh.detected(i, 0) << " " << mesh.detected(i, 1) << "\n";
		//			ind.push_back(i);
		//		}
		//	//}
		//}
		//for (int i = ind.size(); i > 0; i--)
		//{
		//	delete_vertex(ind[i-1]);
		//}
	}

	Eigen::VectorXi State::increase_boundary(const Eigen::VectorXi &boundary)
	{
		std::vector<int> new_boundary;
		for (int i = 0; i < boundary.size(); i++)
		{
			new_boundary.push_back(boundary(i));
			for (int j = 0; j < mesh.adj[boundary(i)].size(); j++)
			{
				new_boundary.push_back(mesh.adj[boundary(i)][j]);
			}
		}

		std::sort(new_boundary.begin(), new_boundary.end());
		new_boundary.erase(std::unique(new_boundary.begin(), new_boundary.end()), new_boundary.end());

		Eigen::VectorXi bboundary(new_boundary.size());
		for (int i = 0; i < new_boundary.size(); i++)
		{
			bboundary(i) = new_boundary[i];
		}
		return bboundary;
	}

	void State::untangle()
	{
		regions.clear();
		mesh.untangle();
	}


	class LocalThreadStorage
	{
	public:
		Eigen::MatrixXd V;
		DetectionParams params;

		LocalThreadStorage(int x)
		{ }
	};


	void State::detect_vertices()
	{
		// clear previous
		hull_vertices.resize(0, 0); //needed for lloyd
		hull_faces.resize(0, 0);
		hull_polygon.resize(0, 0);
		regions.clear();

		mesh3d.clear();

		Eigen::MatrixXd V;
		DetectionParams params;

#ifdef asdUSE_TBB

		typedef tbb::enumerable_thread_specific< LocalThreadStorage > LocalStorage;
		LocalStorage storages(LocalThreadStorage(1));

		const int n_threads = tbb::task_scheduler_init::default_num_threads();
		const float sigma_step = 0.05;
		const float minimum_fraction = 0.4; // higher than 50% removes duplicates
		const float minimum_pixel_spacing = 3;

		Eigen::VectorXf sigmas = Eigen::VectorXf::LinSpaced(n_threads, std::max(1.0, sigma - sigma_step*n_threads/2.), sigma + sigma_step * n_threads / 2);

		tbb::parallel_for(tbb::blocked_range<int>(0, n_threads), [&](const tbb::blocked_range<int> &r) {
			LocalStorage::reference loc_storage = storages.local();
			for (int e = r.begin(); e != r.end(); ++e) {
				point_source_detection(img, sigmas(e), loc_storage.V, loc_storage.params);
			}
		}
		);

		// Concatenate all found points into one vector, and their params into a struct
		std::vector<std::pair<double, double>> V_all;
		DetectionParams ptmp;

		for (LocalStorage::iterator it = storages.begin(); it != storages.end(); ++it)
		{
			const LocalThreadStorage &storage = *it;

			for (int row = 0; row < it->V.rows(); row++)
			{
				double x = storage.V(row, 0);
				double y = storage.V(row, 1);

				V_all.push_back(std::pair<double, double>(x, y));
				ptmp.push_back_params(it->params.get_index(row));
			}
		}

		Eigen::MatrixXd tmpV(V_all.size(), 2);
		for (int i = 0; i < V_all.size(); i++)
		{
			tmpV(i, 0) = V_all[i].first;
			tmpV(i, 1) = V_all[i].second;
		}

		// Vi: indices of selected Vertices, Vj: Index of vertex replacing them
		Eigen::MatrixXi Vi, Vj;
		igl::remove_duplicate_vertices(tmpV, 1.0, V, Vi, Vj);

		// Average V values and count entries
		Eigen::VectorXi count(V.rows());
		V.setZero(V.rows(),2);
		count.setZero(V.rows());

		DetectionParams P_final, P_tmp;
		P_final.setZero(V.rows());
		P_tmp.setZero(1);

		for (int i = 0; i < Vj.size(); i++)
		{
			V(Vj(i), 0) += tmpV(i, 0);
			V(Vj(i), 1) += tmpV(i, 1);

			P_tmp = P_final.get_index(Vj(i));

			P_tmp.sum(ptmp.get_index(i));
			P_final.set_from(P_tmp, Vj(i));

			count(Vj(i))++;

		}

		// Use only vertices that have been detected minimum_fraction * n_threads times
		Eigen::MatrixXd V_final(V.rows(),2);
		int index = 0;
		for (int i = 0; i < V.rows(); i++)
		{
			if (count(i) > minimum_fraction*n_threads)
			{
				V_final.row(index) = V.row(i)/count(i);
				P_tmp = P_final.get_index(i);
				P_tmp.divide(count(i));

				P_final.set_from(P_tmp, index);
				index++;
			}
		}

		P_final.conservative_resize(index);
		V_final.conservativeResize(index, 2);

		V = V_final;
		params = P_final;
#else

		point_source_detection(img, sigma, otsu_multiplier, V, params);

#endif
		if (V.cols() != 3)
		{
			V.conservativeResize(V.rows(), 3);
			V.col(2).setConstant(0);
		}
		mesh.clear();
		mesh.detect_vertices(V, params);
	}

	void State::relax_with_lloyd()
	{
		regions.clear();
		mesh.relax_with_lloyd(lloyd_iterations, hull_vertices, hull_faces, fix_regular_regions);
	}

	void State::detect_bad_regions()
	{
		//fill the region std::vector, no loop

		regions.clear();

		// Calculate the laplacian energy with respect to the original positions
		Eigen::VectorXd energy;
		laplace_energy(mesh.detected, mesh.triangles, energy); //any other energy?
		//laplace_energy(points, triangles, energy);

		// Find the degree of each vertex
		Eigen::VectorXi degree;
		mesh.vertex_degree(degree);

		// Determine whether each vertex passes the criterium for bad mesh
		double avg = energy.mean();
		Eigen::Matrix<bool, 1, Eigen::Dynamic> crit_pass(energy.rows());
		crit_pass.setConstant(false);
		for (int i = 0; i < energy.rows(); i++)
		{
			if (energy(i) > energy_variation_from_mean*avg || degree(i) != 6)
			{
				crit_pass(i) = true;
			}
		}

		// boundaries of the image must not be considered bad, but instead used to enclose regions inside the image
		Eigen::VectorXi boundary = increase_boundary(mesh.boundary);
		boundary = increase_boundary(boundary);

		for (int i = 0; i < boundary.size(); i++)
		{
			crit_pass(boundary(i)) = false;
		}

		// Don't find regions that have been solved
		for (int i = 0; i < mesh.solved_vertex.size(); i++)
		{
			if (mesh.solved_vertex(i))
				crit_pass(i) = false;
		}

		// vertices manually picked as good must be removed from bad
		for (int i = 0; i < mesh.vertex_status_fixed.size(); i++)
		{
			if(mesh.vertex_status_fixed(i) == 1)
				crit_pass(i) = false;
			if (mesh.vertex_status_fixed(i) == -1)
				crit_pass(i) = true;
		}

		// Find connected regions where the criterium was not passed

		//maybe avoid regions id, fill directly region_vertices and use in the next method
		Eigen::VectorXi regions_id;
		region_grow(mesh.adj, crit_pass, regions_id);
		int n_regions = regions_id.maxCoeff();

		regions.resize(n_regions);
		for (int i = 0; i < n_regions; ++i)
		{
			regions[i].find_points(regions_id, i + 1);
			regions[i].find_triangles(mesh.triangles, regions_id, i + 1);
		}

		int index = 0;
		for (auto & r : regions)
		{
			r.bounding(mesh.triangles, mesh.points);

			bool repeat = false;
			do {
				repeat = r.fix_missing_points(mesh.triangles);
			} while (repeat);

			repeat = false;
			do {
				repeat = r.add_orphaned_triangles(mesh);
			} while (repeat);

			// fix pinched regions by closing area at pinch
			repeat = false;
			do {
				repeat = r.fix_pinched_region(mesh);
			} while (repeat);

			repeat = false;
			do {
				repeat = r.fix_missing_points(mesh.triangles);
			} while (repeat);

			++index;

			//r.find_triangles(mesh.triangles);
		}

		erase_small_regions();
	}

	void State::erase_small_regions()
	{
		counter_small_region = 0;
		//first remove/grow the ones with too small
		auto it = regions.begin();
		while (it != regions.end())
		{
			int nPolygon = it->region_boundary.size(); // length of current edge
			if (nPolygon < 7)
			{
				it = regions.erase(it);
			}
			else
			{
				++it;
			}
		}


		//grow the ones with wrong valency
		//split the big ones
	}

	void State::check_regions()
	{
		for (auto & r : regions)
		{
			r.check_region(mesh.detected, mesh.points, mesh.triangles, mesh.adj);
			//std::cout << r.status << std::endl;
		}
	}

	void State::check_region(const int index)
	{
		regions[index].check_region(mesh.detected, mesh.points, mesh.triangles, mesh.adj);
	}

	// void State::fix_regions()
	// {
	// 	bool must_continue = true;

	// 	std::vector<int> to_remove;
	// 	// loop through remaining regions and act depending on status
	// 	for (auto & r : regions)
	// 	{
	// 		// If too many vertices are in the region, find the least likely to be correct and remove
	// 		if (r.status == Region::TOO_MANY_VERTICES)
	// 		{
	// 			int n = r.region_interior.rows();
	// 			Eigen::VectorXd x_pstd(n), y_pstd(n), pval_Ar(n);
	// 			for (int i = 0; i < r.region_interior.size(); i++)
	// 			{
	// 				x_pstd(i) = mesh.params.std_x(r.region_interior(i));
	// 				y_pstd(i) = mesh.params.std_y(r.region_interior(i));
	// 				pval_Ar(i) = mesh.params.pval_Ar(r.region_interior(i));
	// 			}
	// 			int xMax, yMax, pMax;
	// 			x_pstd.maxCoeff(&xMax);
	// 			y_pstd.maxCoeff(&yMax);
	// 			pval_Ar.maxCoeff(&pMax);

	// 			if (xMax == yMax && xMax == pMax)
	// 			{
	// 				int ind = r.region_interior(xMax);
	// 				to_remove.push_back(ind);

	// 				// pop and rebuild current region
	// 				// solve rebuilt region
	// 			}
	// 		}
	// 	}


	// 	std::sort(to_remove.begin(), to_remove.end());

	// 	for (int i = to_remove.size() - 1; i >= 0; --i)
	// 	{
	// 		delete_vertex(to_remove[i]);
	// 	}

	// 	for (auto & r : regions)
	// 		r.find_triangles(mesh.triangles);

	// 	check_regions();
	// 	//resolve_regions();
	// }


	void State::grow_region(const int index)
	{
		regions[index].grow(mesh.points, mesh.triangles);
	}

	void State::grow_regions()
	{
		Eigen::Matrix<bool, 1, Eigen::Dynamic> crit_pass(mesh.points.rows());
		crit_pass.setConstant(false);
		for (int j = 0; j < regions.size(); j++)
		{
			auto &r = regions[j];
			for (size_t i = 0; i < r.region_interior.size(); i++)
			{
				crit_pass(r.region_interior(i)) = true;
			}
			if (j > 0)
			{
				for (size_t i = 0; i < r.region_boundary.size(); i++)
				{
					crit_pass(r.region_boundary(i)) = true;
				}
			}
		}

		// Find connected regions where the criterium was not passed

		//maybe avoid regions id, fill directly region_vertices and use in the next method
		Eigen::VectorXi regions_id;
		region_grow(mesh.adj, crit_pass, regions_id);
		int n_regions = regions_id.maxCoeff();

		regions.resize(n_regions);
		for (int i = 0; i < n_regions; ++i)
		{
			regions[i].find_points(regions_id, i + 1);
			regions[i].find_triangles(mesh.triangles, regions_id, i + 1);
		}

		for (auto & r : regions)
		{
			r.bounding(mesh.triangles, mesh.points);
		}
	}

	void State::resolve_region(const int index)
	{
		Eigen::MatrixXd new_points;
		Eigen::MatrixXi new_triangles;

		counter_invalid_neigh = 0;
		counter_infeasible_region = 0;
		counter_solved_regions = 0;

		//remove from regions all the solved one
		//grow the others and regurobi the regions max 4 time

		auto &region = regions[index];

		double gurobi_max_time = gurobi_time_limit_long;
		region.resolve(mesh.detected, mesh.points, perm_possibilities, gurobi_max_time, new_points, new_triangles, true);
		//gurobi r
		if (region.status == Region::NOT_PROPERLY_CLOSED)
		{
			counter_invalid_neigh++;
			return;
		}
		else if (region.status == Region::NO_SOLUTION)
		{
			counter_infeasible_region++;
			return;
		}

		// Map q.T back to global indices
		assert(region.region_boundary.size() + region.region_interior.size() == new_points.rows());

		//add this into the region
		Eigen::VectorXi local_to_global;
		region.local_to_global(local_to_global);

		mesh.local_update(local_to_global, new_points, new_triangles, region.region_triangles);

		counter_solved_regions++;

		// Mark region as good in mesh
		mesh.mark_vertex_as_solved(region.region_interior);

		regions.erase(regions.begin() + index);

		for (auto &r : regions)
			r.find_triangles(mesh.triangles);
	}

	void State::resolve_regions()
	{
		Eigen::MatrixXd new_points;
		Eigen::MatrixXi new_triangles;

		counter_invalid_neigh = 0;
		counter_infeasible_region = 0;
		counter_solved_regions = 0;

		// loop through all the boundaries and solve individually. Possibly skip first one, as it's the boundary of the image
		auto it = regions.begin();
		while (it != regions.end())
		{
			//remove from regions all the solved one
			//grow the others and regurobi the regions max 4 time

			auto &region = *it;
			region.find_triangles(mesh.triangles);
			double gurobi_max_time = gurobi_time_limit_short;
			region.resolve(mesh.detected, mesh.points, perm_possibilities, gurobi_max_time, new_points, new_triangles);
			//gurobi r
			if (region.status == Region::NOT_PROPERLY_CLOSED)
			{
				counter_invalid_neigh++;
				++it;
				continue;
			}
			else if (region.status == Region::NO_SOLUTION)
			{
				counter_infeasible_region++;
				++it;
				continue;
			}
			else if (region.status != 0)
			{
				++it;
				continue;
			}

			if (region.region_boundary.size() + region.region_interior.size() != new_points.rows())
			{
				counter_infeasible_region++; //maybe new name (point missing or too much)
				++it;
				continue;
			}


			// Map q.T back to global indices
			assert(region.region_boundary.size() + region.region_interior.size() == new_points.rows());


			Eigen::VectorXi local_to_global;
			region.local_to_global(local_to_global);

			mesh.local_update(local_to_global, new_points, new_triangles, region.region_triangles);

			// Mark region as good in mesh
			mesh.mark_vertex_as_solved(region.region_interior);

			it = regions.erase(it);
			counter_solved_regions++;

		}

		for (auto &r : regions)
			r.find_triangles(mesh.triangles);

	}

	void State::final_relax()
	{
		//increase boundary by one row for fixation
		Eigen::VectorXi boundary = increase_boundary(mesh.boundary);
		mesh.final_relax(boundary);
	}

	void State::delete_vertex(const int index)
	{

		if (regions.empty())
		{
			mesh.delete_vertex(index, true);
			//maybe recompute the hull
		}
		else
		{
			//find the region
			int region_ind = find_region_by_interior_vertex(index);
			if (region_ind == -1)
			{
				regions.clear();
				delete_vertex(index);
				return;
			}
			//and call region.delete_vertex(index) which calls the new local_update of mesh
			//the region delete vertex shold get the local id of the vertex compute the local2global without that vertex and retriangulate the region
			regions[region_ind].delete_vertex(mesh,index);

			detect_bad_regions();
		}

		for (auto &r : regions)
			r.find_triangles(mesh.triangles);

		// maybe clear regions and recompute hull
	}

	void State::add_vertex(Eigen::Vector3d new_point)
	{
		mesh.add_vertex(new_point);

	}

	// void State::init_3d_mesh()
	// {
	// 	if(image_from_pillars)
	// 		mesh3d.init_pillars(mesh, eps, I, L);
	// 	else
	// 		mesh3d.init_nano_dots(mesh, padding_size, thickness, E, nu, formulation);
	// }

	void State::mesh_2d_adaptive() {
		// Create background mesh with a scalar field = norm of the displacements of the dots
		Eigen::MatrixXd V;
		Eigen::MatrixXi F;
		Eigen::VectorXd S;
		mesh.get_background_mesh(scaling, V, F, S, padding_size);

		// Rescale displacement field
		double vmin = V.minCoeff();
		double vmax = V.maxCoeff();
		S = (S.array() - S.minCoeff()) / std::max(1e-9, (S.maxCoeff() - S.minCoeff()));
		S = (1.0 - S.array()).pow(power) * (adaptive_mesh_size[1] - adaptive_mesh_size[0]) + adaptive_mesh_size[0];
		std::cout << S.transpose() << std::endl;
		S *= median_edge_length(mesh.detected, mesh.triangles) * scaling / std::max(1e-9, vmax - vmin);
		std::cout << S.transpose() << std::endl;
		std::cout << median_edge_length(mesh.detected, mesh.triangles) * scaling / std::max(1e-9, vmax - vmin) << std::endl;

		// Remesh triangle mesh
		V = V.array() / std::max(1e-9, vmax - vmin);
 		remesh_adaptive_2d(V, F, S, V, F, mmg_options);
		V = V.array() * std::max(1e-9, vmax - vmin);

		// Mesh volume adaptively based on background mesh
		mesh3d.V.resize(V.rows(), 3);
		mesh3d.V.leftCols<2>() = V.leftCols<2>();
		mesh3d.V.col(2).setZero();

		mesh3d.F = F;
		mesh3d.T = Eigen::MatrixXi(0, 4);
	}

	void State::extrude_2d_mesh() {
		if (mesh3d.V.size() == 0) { mesh_2d_adaptive(); }

		double zmin = mesh3d.V.col(2).minCoeff();
		double zmax = mesh3d.V.col(2).maxCoeff();

		// Remesh 2d surface in repeated calls
		if (std::abs(zmin - zmax) > 1e-6) {
			mesh_2d_adaptive();
		}

		// Extrude into a 3d mesh
		extrude_mesh(mesh3d.V, mesh3d.F, -thickness, mesh3d.V, mesh3d.F, mesh3d.T);
	}

	void State::mesh_3d_uniform() {
		extrude_2d_mesh(); // too lazy to recode this

		Eigen::MatrixXd BV;
		Eigen::MatrixXi BF;
		igl::bounding_box(mesh3d.V, BV, BF);

		// Target tetrahedron volume
		double a = uniform_mesh_size * median_edge_length(mesh.detected, mesh.triangles) * scaling;
		double target_volume = std::sqrt(2.0) / 12.0 * a * a * a;

		std::stringstream buf;
		buf.precision(100);
		buf.setf(std::ios::fixed, std::ios::floatfield);
		buf << "Qpq1.414a" << target_volume;

		// Mesh volume
		igl::copyleft::tetgen::tetrahedralize(BV, BF, buf.str(), mesh3d.V, mesh3d.T, mesh3d.F);
		polyfem::orient_closed_surface(mesh3d.V, mesh3d.F);
	}

	void State::remesh_3d_adaptive() {
		// Generate background mesh if needed
		if (mesh3d.V.size() == 0 || mesh3d.T.rows() == 0) { mesh_3d_uniform(); }
		double zmin = mesh3d.V.col(2).minCoeff();
		double zmax = mesh3d.V.col(2).maxCoeff();
		if (std::abs(zmin - zmax) > 1e-6) { mesh_3d_uniform(); }

		// Scalar field to interpolate
		Eigen::MatrixXd V;
		Eigen::MatrixXi F;
		Eigen::VectorXd S;
		mesh.get_background_mesh(scaling, V, F, S, padding_size);

		// Interpolate
		Eigen::MatrixXd VP = mesh3d.V.leftCols<2>();
		polyfem::InterpolatedFunction2d aux(S, V, F);
		S = aux.interpolate(VP);

		// Rescale displacement field
		double vmin = mesh3d.V.minCoeff();
		double vmax = mesh3d.V.maxCoeff();
		S = (S.array() - S.minCoeff()) / std::max(1e-9, (S.maxCoeff() - S.minCoeff()));
		S = (1.0 - S.array()).pow(power) * (adaptive_mesh_size[1] - adaptive_mesh_size[0]) + adaptive_mesh_size[0];
		S *= median_edge_length(mesh.detected, mesh.triangles) * scaling / std::max(1e-9, vmax - vmin);

		// Remesh volume mesh
		mesh3d.V = mesh3d.V.array() / std::max(1e-9, vmax - vmin);
 		remesh_adaptive_3d(mesh3d.V, mesh3d.T, S, mesh3d.V, mesh3d.F, mesh3d.T, mmg_options);
		mesh3d.V = mesh3d.V.array() * std::max(1e-9, vmax - vmin);
	}

	void State::analyze_3d_mesh() {
		if(image_from_pillars)
		{
			mesh3d.init_pillars(mesh, eps, I, L, scaling);
		}
		else
		{
			mesh3d.init_nano_dots(mesh, padding_size, thickness, E, nu, scaling, formulation);
		}
	}

	void State::reset_state()
	{
		mesh.reset();
		mesh3d.clear();
		regions.clear();
		compute_hull();
	}


}// namespace cellogram


