////////////////////////////////////////////////////////////////////////////////
#include "State.h"

#include <cellogram/point_source_detection.h>
#include <cellogram/image_reader.h>
#include <cellogram/convex_hull.h>

#include <cellogram/voronoi.h>

#include <cellogram/laplace_energy.h>
#include <cellogram/tri2hex.h>
#include <cellogram/region_grow.h>
#include <cellogram/PolygonUtils.h>

#include <igl/slice.h>
#include <igl/list_to_matrix.h>
#include <igl/point_in_poly.h>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

	// -----------------------------------------------------------------------------

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

		if (size(regions1) == 0 || size(regions2) == 0)
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

		regions[index].split_region(split_end_points);

	}

	// -----------------------------------------------------------------------------

	bool State::load(const std::string & path)
	{
		bool ok = false;
		// clear previous
		hull_vertices.resize(0, 0); //needed for lloyd
		hull_faces.resize(0,0);
		hull_polygon.resize(0, 0);
		regions.clear();

		// Load points
		ok = mesh.load(path);

		compute_hull();

		// load settings
		using json = nlohmann::json;
		std::ifstream json_in(path + "/cellogram/settings.json");

		json json_data;
		json_in >> json_data;

		std::vector<double> tmp;
		lloyd_iterations = json_data["lloyd_iterations"];
		energy_variation_from_mean = json_data["energy_variation_from_mean"];
		perm_possibilities = json_data["perm_possibilities"];
		sigma = json_data["sigma"];

		return ok;

	}

	bool State::load_image(const std::string fname)
	{
		bool ok = read_image(fname, img);

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

	bool State::save(const std::string & path)
	{
		using json = nlohmann::json;


		json json_data;
		json_data["lloyd_iterations"] = lloyd_iterations;
		json_data["energy_variation_from_mean"] = energy_variation_from_mean;
		json_data["perm_possibilities"] = perm_possibilities;
		json_data["sigma"] = sigma;

		{
			std::ofstream json_out(path + "/settings.json");
			json_out << json_data.dump(4) << std::endl;
			json_out.close();
		}

		mesh.save(path);

		{
			std::ofstream hull_path(path + "/hull.vert");
			hull_path << hull_vertices << std::endl;
			hull_path.close();
		}

		{
			std::ofstream hull_path(path + "/hull.tri");
			hull_path << hull_faces << std::endl;
			hull_path.close();
		}

		//std::vector<Region> regions;
		//std::vector<int> fixed_as_good;



		return true;
	}

	void State::compute_hull()
	{
		//convex_hull(points, boundary);
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

#ifdef WITH_UNTANGLER
	void State::untangle()
	{
		mesh.untangle();
	}
#endif

	void State::detect_vertices()
	{
		// clear previous
		hull_vertices.resize(0, 0); //needed for lloyd
		hull_faces.resize(0, 0);
		hull_polygon.resize(0, 0);
		regions.clear();
		
		Eigen::MatrixXd V;
		DetectionParams params;

		point_source_detection(img, sigma, V, params);

		if (V.cols() != 3)
		{
			V.conservativeResize(V.rows(), 3);
			V.col(2).setConstant(0);
		}

		mesh.detect_vertices(V, params);
	}

	void State::relax_with_lloyd()
	{
		//reset the state
		mesh.relax_with_lloyd(lloyd_iterations, hull_vertices, hull_faces);
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

			++index;
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
			r.status = r.check_region(mesh.detected, mesh.points, mesh.triangles, mesh.adj);
			//std::cout << r.status << std::endl;
		}
	}

	void State::check_region(const int index)
	{
		regions[index].status = regions[index].check_region(mesh.detected, mesh.points, mesh.triangles, mesh.adj);
	}

	void State::fix_regions()
	{
		bool must_continue = true;

		std::vector<int> to_remove;
		// loop through remaining regions and act depending on status
		for (auto & r : regions)
		{
			// If too many vertices are in the region, find the least likely to be correct and remove
			if (r.status == Region::TOO_MANY_VERTICES)
			{
				int n = r.region_interior.rows();
				Eigen::VectorXd x_pstd(n), y_pstd(n), pval_Ar(n);
				for (int i = 0; i < r.region_interior.size(); i++)
				{
					x_pstd(i) = mesh.params.std_x(r.region_interior(i));
					y_pstd(i) = mesh.params.std_y(r.region_interior(i));
					pval_Ar(i) = mesh.params.pval_Ar(r.region_interior(i));
				}
				int xMax, yMax, pMax;
				x_pstd.maxCoeff(&xMax);
				y_pstd.maxCoeff(&yMax);
				pval_Ar.maxCoeff(&pMax);

				if (xMax == yMax && xMax == pMax)
				{
					int ind = r.region_interior(xMax);
					to_remove.push_back(ind);

					// pop and rebuild current region
					// solve rebuilt region
				}
			}
		}


		std::sort(to_remove.begin(), to_remove.end());

		for (int i = to_remove.size() - 1; i >= 0; --i)
		{
			delete_vertex(to_remove[i]);
		}

		for (auto & r : regions)
			r.find_triangles(mesh.triangles);

		check_regions();
		//resolve_regions();
	}


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
		auto status = region.resolve(mesh.detected, mesh.points, perm_possibilities, gurobi_max_time, new_points, new_triangles, true);
		region.status = status;
		//gurobi r
		if (status == Region::NOT_PROPERLY_CLOSED)
		{
			counter_invalid_neigh++;
			return;
		}
		else if (status == Region::NO_SOLUTION)
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
			auto status = region.resolve(mesh.detected, mesh.points, perm_possibilities, gurobi_max_time, new_points, new_triangles);
			region.status = status;
			//gurobi r
			if (status == Region::NOT_PROPERLY_CLOSED)
			{
				counter_invalid_neigh++;
				++it;
				continue;
			}
			else if (status == Region::NO_SOLUTION)
			{
				counter_infeasible_region++;
				++it;
				continue;
			}
			else if (status != 0)
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
		mesh.final_relax();
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

	void State::init_3d_mesh()
	{
		mesh3d.init(mesh, padding_size, thickness);
	}

	void State::reset_state()
	{
		mesh.reset();
		mesh3d.clear();
		regions.clear();
		compute_hull();
	}


}// namespace cellogram


