////////////////////////////////////////////////////////////////////////////////
#include "State.h"

#include <cellogram/load_points.h>
#include <cellogram/load_points.h>
#include <cellogram/convex_hull.h>
#include <cellogram/delaunay.h>
#include <cellogram/voronoi.h>

#include <cellogram/laplace_energy.h>
#include <cellogram/tri2hex.h>
#include <cellogram/vertex_degree.h>
#include <cellogram/region_grow.h>

#include <igl/slice.h>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

	namespace
	{
		void removeRow(Eigen::MatrixXi& matrix, unsigned int rowToRemove)
		{
			unsigned int numRows = matrix.rows() - 1;
			unsigned int numCols = matrix.cols();

			if (rowToRemove < numRows)
				matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.bottomRows(numRows - rowToRemove);

			matrix.conservativeResize(numRows, numCols);
		}

		void replaceTriangles(const Eigen::MatrixXi tNew, const Eigen::MatrixXi &to_remove, Eigen::MatrixXi &triangles)
		{
			std::vector<int> removeIdx(to_remove.data(), to_remove.data() + to_remove.rows() * to_remove.cols());
			std::sort(removeIdx.begin(), removeIdx.end());
			// remove all the triangles that connect to the internal of ROI
			for (int i = removeIdx.size() - 1; i >= 0; i--)
			{
				removeRow(triangles, removeIdx[i]);
			}

			// add new rows at the end of triangles
			Eigen::MatrixXi tmp = triangles;
			triangles.resize(triangles.rows() + tNew.rows(), triangles.cols());

			triangles << tmp, tNew;
		}
	}

	// -----------------------------------------------------------------------------

	State &State::state() {
		static State instance;
		return instance;
	}

	// -----------------------------------------------------------------------------

	bool State::load(const std::string & path)
	{
		if (path.empty()) { return false; }
		// clear previous
		detected.resize(0,0); // detected (unmoved) point positions
		points.resize(0, 0); // relaxed point positions
		triangles.resize(0, 0); // triangular mesh
		adj.clear(); // adjaceny list of triangluar mesh
		boundary.resize(0); // list of vertices on the boundary
		hull_vertices.resize(0, 0); //needed for lloyd
		hull_faces.resize(0,0);
		hull_polygon.resize(0, 0);
		regions.clear();

		// Load points
		load_points(path, detected);
		points = detected;
		compute_triangulation();

		return true;
	}

	bool State::save(const std::string & path)
	{
		return false;
	}

	void State::compute_hull()
	{
		//convex_hull(points, boundary);
		loose_convex_hull(detected, boundary);
		int dims = (int)detected.cols();


		// Compute polygon of the convex hull
		Eigen::VectorXi I = boundary;
		Eigen::VectorXi J = Eigen::VectorXi::LinSpaced(dims, 0, dims - 1);

		igl::slice(detected, I, J, hull_polygon);

		// Offset by epsilon
		// offset_polygon(hull_polygon, hull_polygon, 1);

		// Draw filled polygon
		// triangulate_convex_polygon(hull_polygon, hull_vertices, hull_faces);
		triangulate_polygon(hull_polygon, hull_vertices, hull_faces);
	}

	void State::compute_triangulation()
	{
		delaunay_triangulation(points, triangles);

		// Calculate the graph adjancency
		adjacency_list(triangles, adj);
	}

	void State::relax_with_lloyd()
	{
		//reset the state
		points = detected;
		compute_triangulation();

		lloyd_relaxation(points, boundary, lloyd_iterations, hull_vertices, hull_faces);
		compute_triangulation();
	}

	void State::detect_bad_regions()
	{
		//fill the region std::vector, no loop

		regions.clear();

		// Calculate the laplacian energy with respect to the original positions
		Eigen::VectorXd energy;
		laplace_energy(detected, triangles, energy); //any other energy?
		//laplace_energy(points, triangles, energy);

		// Find the degree of each vertex
		Eigen::VectorXi degree;
		vertex_degree(triangles, degree);

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
		for (int i = 0; i < boundary.size(); i++)
		{
			crit_pass(boundary(i)) = false;
		}

		// Find connected regions where the criterium was not passed

		//maybe avoid regions id, fill directly region_vertices and use in the next method
		Eigen::VectorXi regions_id;
		region_grow(adj, crit_pass, regions_id);
		int n_regions = regions_id.maxCoeff();

		regions.resize(n_regions);
		for (int i = 0; i < n_regions; ++i)
		{
			regions[i].find_points(regions_id, i + 1);
			regions[i].find_triangles(triangles, regions_id, i + 1);
		}

		for (auto & r : regions)
		{
			r.bounding(triangles, points);

			r.fix_missing_points(triangles);
		}



	}

	void State::fix_regions()
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


		//grow the ones with wong valency
		//split the big ones
	}


	void State::grow_region(const int index)
	{
		regions[index].grow(points, triangles);
	}

	void State::grow_regions()
	{
		Eigen::Matrix<bool, 1, Eigen::Dynamic> crit_pass(points.rows());
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
		region_grow(adj, crit_pass, regions_id);
		int n_regions = regions_id.maxCoeff();

		regions.resize(n_regions);
		for (int i = 0; i < n_regions; ++i)
		{
			regions[i].find_points(regions_id, i + 1);
			regions[i].find_triangles(triangles, regions_id, i + 1);
		}

		for (auto & r : regions)
		{
			r.bounding(triangles, points);
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

		auto state = region.resolve(detected, points, triangles, adj, perm_possibilities, new_points, new_triangles);
		//gurobi r
		if (state == Region::NOT_PROPERLY_CLOSED)
		{
			counter_invalid_neigh++;
			return;
		}
		else if (state == Region::NO_SOLUTION)
		{
			counter_infeasible_region++;
			return;
		}

		// Map q.T back to global indices
		assert(region.region_boundary.size() + region.region_interior.size() == new_points.rows());
		std::cout << new_points << std::endl;
		const int n_region_boundary = region.region_boundary.size();
		const int n_region_interior = region.region_interior.size();

		Eigen::VectorXi local_to_global(n_region_boundary + n_region_interior);

		for (int i = 0; i < n_region_boundary; i++)
		{
			const int global_index = region.region_boundary[i];
			local_to_global(i) = global_index;
			points.row(global_index) = new_points.row(i);

		}
		for (int i = 0; i < n_region_interior; i++)
		{
			const int global_index = region.region_interior[i];
			local_to_global(i + n_region_boundary) = global_index;
			points.row(global_index) = new_points.row(i + n_region_boundary);
		}


		Eigen::MatrixXi tGlobal = Eigen::MatrixXi(new_triangles.rows(), 3);
		for (int i = 0; i < new_triangles.rows(); i++)
		{
			for (int j = 0; j < 3; j++)
			{
				tGlobal(i, j) = local_to_global(new_triangles(i, j));
			}

		}

		replaceTriangles(tGlobal, region.region_triangles, triangles);

		counter_solved_regions++;

		adjacency_list(triangles, adj);

		regions.erase(regions.begin() + index);

		for (auto &r : regions)
			r.find_triangles(triangles);
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
		++it; //skip the first one, Tobi knows why
		while (it != regions.end())
		{
			//remove from regions all the solved one
			//grow the others and regurobi the regions max 4 time

			auto &region = *it;

			auto state = region.resolve(detected, points, triangles, adj, perm_possibilities, new_points, new_triangles);
			//gurobi r
			if (state == Region::NOT_PROPERLY_CLOSED)
			{
				counter_invalid_neigh++;
				++it;
				continue;
			}
			else if (state == Region::NO_SOLUTION)
			{
				counter_infeasible_region++;
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
			const int n_region_boundary = region.region_boundary.size();
			const int n_region_interior = region.region_interior.size();

			Eigen::VectorXi local_to_global(n_region_boundary + n_region_interior);

			for (int i = 0; i < n_region_boundary; i++)
			{
				const int global_index = region.region_boundary[i];
				local_to_global(i) = global_index;
				points.row(global_index) = new_points.row(i);

			}
			for (int i = 0; i < n_region_interior; i++)
			{
				const int global_index = region.region_interior[i];
				local_to_global(i + n_region_boundary) = global_index;
				points.row(global_index) = new_points.row(i + n_region_boundary);
			}


			Eigen::MatrixXi tGlobal = Eigen::MatrixXi(new_triangles.rows(), 3);
			for (int i = 0; i < new_triangles.rows(); i++)
			{
				for (int j = 0; j < 3; j++)
				{
					tGlobal(i, j) = local_to_global(new_triangles(i, j));
				}

			}

			replaceTriangles(tGlobal, region.region_triangles, triangles);


			it = regions.erase(it);
			counter_solved_regions++;
		}

		adjacency_list(triangles, adj);

		for (auto &r : regions)
			r.find_triangles(triangles);
	}

	void State::delete_vertex(const int index)
	{
		Eigen::MatrixXd tmp(detected.rows() - 1, detected.cols());
		tmp.block(0, 0, index, detected.cols()) = detected.block(0, 0, index, detected.cols());
		tmp.block(index, 0, tmp.rows() - index, detected.cols()) = detected.block(index+1, 0, tmp.rows() - index, detected.cols());

		detected = tmp;
		points = detected;
		compute_triangulation();
		compute_hull();
	}

	void State::add_vertex(Eigen::Vector3d new_point)
	{
		Eigen::MatrixXd tmp(detected.rows() + 1, detected.cols());
		tmp.block(0, 0, detected.rows(), detected.cols()) = detected;
		//tmp.block(detected.rows(), 0, 1, detected.cols()) = new_point.transpose();
		tmp.row(detected.rows()) = new_point.transpose();

		std::cout << tmp << std::endl;
		std::cout << "\nNew point:\n" << new_point.transpose() << std::endl;

		//detected.resizeLike(tmp);
		//points.resizeLike(tmp);
		detected = tmp;
		points = detected;
		compute_triangulation();
		compute_hull();
	}


}// namespace cellogram


