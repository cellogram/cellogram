////////////////////////////////////////////////////////////////////////////////
#include "cellogram/mesh_solver.h"
#include <cellogram/laplace_energy.h>
#include <cellogram/tri2hex.h>
#include <cellogram/vertex_degree.h>
#include <cellogram/region_grow.h>
#include <gurobi_solver/state.h>
#include <gurobi_solver/generateQ.h>
#include <gurobi_solver/gurobiModel.h>

#include <cellogram/boundary_loop.h>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {
	namespace
	{
		bool is_neigh_valid(const VectorXi &neigh)
		{
			// check if neigh actually closes a loop by tacing
			Eigen::VectorXi turns = neigh.array() - 4;

			Eigen::MatrixXi dirs_even(6, 2);
			Eigen::MatrixXi dirs_odd(6, 2);
			dirs_even << 1, 0,
				0, 1,
				-1, 1,
				-1, 0,
				-1, -1,
				0, -1;
			dirs_odd << 1, 0,
				1, 1,
				0, 1,
				-1, 0,
				0, -1,
				1, -1;

			int dir = 0;

			std::vector<Eigen::Vector2i> pairs;

			Eigen::Vector2i current_pair;
			current_pair << 0, 0;
			pairs.push_back(current_pair);

			// Follow the boundary
			for (unsigned i = 0; i < turns.size(); ++i)
			{
				// move in the current direction
				current_pair += (current_pair(1) % 2) == 0 ? dirs_even.row(dir).transpose() : dirs_odd.row(dir).transpose();

				// save the current position
				pairs.push_back(current_pair);

				//update the direction
				dir = (dir + dirs_even.rows() + turns((i + 1) % turns.size())) % dirs_even.rows();
			}

			Eigen::MatrixXi ret;

			if (pairs.front() == pairs.back() && 0 == dir)
			{
				return true;
			}
			else
			{
				return false;
			}
		}

		Eigen::VectorXi find_n_neighbor(const Eigen::MatrixXi &F2, const Eigen::VectorXi &current_edge_vertices, std::vector<std::vector<int>> &adj) {
			// This function needs to find the number of correct neighbors to the outside. 
			// This can not be done using the inside of the region, because that region may be faulty
			std::vector<int> internalTri(current_edge_vertices.rows(), 0);
			Eigen::VectorXi neigh(current_edge_vertices.rows());
			for (int j = 0; j < current_edge_vertices.rows(); j++)
			{
				for (int f = 0; f < F2.rows(); f++)
				{
					for (int lv = 0; lv < 3; lv++)
					{
						// count how many triangles in the region this vertex belongs to
						//std::cout << F2(f, lv) << " -- " << current_edges(j) << std::endl;
						if (F2(f, lv) == current_edge_vertices(j)) {
							internalTri[j]++;
						}
					}

				}
				neigh(j) = adj[current_edge_vertices(j)].size() - (internalTri[j] - 1);
			}
			return neigh;
		}
	}

	void Region::compute_edges(const Eigen::MatrixXd &V, Eigen::MatrixXd &bad_P1, Eigen::MatrixXd &bad_P2) {
		// Line boundary

		const int n = region_boundary.size();

		bad_P1.resize(n, 3);
		bad_P2.resize(n, 3);

		for (int j = 0; j < n; j++)
		{
			bad_P1.row(j) = V.row(region_boundary(j));
			bad_P2.row(j) = V.row(region_boundary((j + 1) % n));
		}
	}

	void Region::grow()
	{
		//TODO
	}



	void Region::find_points(const Eigen::VectorXi &region_ids, const int id)
	{
		std::vector<int> vertices;
		for (int j = 0; j < region_ids.rows(); j++)
		{
			// only use vertices with region(j) > 0, because those ==0 are considered good
			if (region_ids(j) == id) {
				vertices.push_back(j);
			}
		}

		region_interior.resize(vertices.size());

		for (int j = 0; j < vertices.size(); j++)
		{
			region_interior(j) = vertices[j];
		}
	}


	void Region::find_triangles(const Eigen::MatrixXi &F, const Eigen::VectorXi &region_ids, const int id)
	{
		std::vector<int> region_faces;
		for (int j = 0; j < region_ids.rows(); j++)
		{
			// only use vertices with region(j) > 0, because those ==0 are considered good
			if (region_ids(j) == id)
			{
				// for this vertex find the all the faces
				for (int f = 0; f < F.rows(); ++f) {
					for (int lv = 0; lv < F.cols(); ++lv) {
						if (j == F(f, lv))
						{
							region_faces.push_back(f);
							break;
						}
					}
				}
			}
		}
		
		std::sort(region_faces.begin(), region_faces.end());
		region_faces.erase(std::unique(region_faces.begin(), region_faces.end()), region_faces.end());

		region_triangles.resize(region_faces.size());

		for (int j = 0; j < region_faces.size(); j++)
		{
			region_triangles(j) = region_faces[j];
		}
	}


	void Region::bounding(const Eigen::MatrixXi &F) {
		boundary_loop(get_triangulation(F), region_boundary);
	}


	int Region::resolve(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, std::vector<std::vector<int>> &adj, const int perm_possibilities, Eigen::MatrixXd  &new_points, Eigen::MatrixXi &new_triangles)
	{
		int nPolygon = region_boundary.size(); // length of current edge

		// 2.1
		// Save current edge points into single vector

		// for each region the belonging triangles need to be extracted and duplicates removed
		Eigen::MatrixXi F2 = get_triangulation(F);

		// 2.3 find number of neighbors

		// 2.2
		// Determine number of connections into the cluster to determine "neigh"
		Eigen::VectorXi neigh = find_n_neighbor(F2, region_boundary, adj);
		if (!is_neigh_valid(neigh)) {
			// this mesh is not properly close
			return NOT_PROPERLY_CLOSED;
		}



		// 3 Pepare for Gurobi
		// vB contains the boundary vertices, vI the internal vertices
		Eigen::MatrixXd vB(nPolygon, 2);
		for (int i = 0; i < nPolygon; i++)
		{
			vB(i, 0) = V(region_boundary(i), 0);
			vB(i, 1) = V(region_boundary(i), 1);
		}
		Eigen::MatrixXd vI(region_interior.size(), 2);
		for (int i = 0; i < region_interior.size(); i++)
		{
			vI(i, 0) = V(region_interior[i], 0);
			vI(i, 1) = V(region_interior[i], 1);
		}
		// Generate perfect mesh in ROI
		State s;
		s.init(vB, vI, neigh);
		s.fill_hole();

		// Check whether vertices inside region is equal to the ones expected
		if (verbose && s.Vdeformed.rows() != s.Vdeformed.rows()) {
			std::cout << "\nDeformed\n" << s.Vdeformed.transpose();
			std::cout << "\nPerfect\n" << s.Vperfect.transpose();
		}

		// Generate adjacency matrix and the laplacian
		generateQ q;
		gurobiModel g;
		q.adjacencyMatrix(s.F);
		q.laplacianMatrix();


		// Deriving Q and constraints for optimization
		q.QforOptimization(s.Vperfect, s.Vdeformed, perm_possibilities);
		q.optimizationConstraints(s.V_boundary.rows());


		// Generate and solve model for gurobi
		g.model(q.Q, q.Aeq);
		if (g.resultX(1) == -1) {
			// no solution found
			return NO_SOLUTION;
		}

		// Map back to indices of coordinates
		q.mapBack(g.resultX);
		//std::cout << g.resultX << std::endl;

		std::vector<Eigen::Triplet<double> > permutation_triplets;
		const int n_points = s.Vperfect.rows();
		for (int i = 0; i < n_points; ++i)
		{
			const int index = i * perm_possibilities;

			for (int j = 0; j < perm_possibilities; ++j)
			{
				if (g.resultX(index + j) == 1)
				{
					const int dest_vertex = q.IDX(j, i);
					permutation_triplets.emplace_back(i, dest_vertex, 1);
					continue;
				}
			}
		}

		// Replace vertex
		Eigen::SparseMatrix<double> permutation(n_points, n_points);
		permutation.setFromTriplets(permutation_triplets.begin(), permutation_triplets.end());
		
		new_points = permutation.transpose() * s.Vperfect;
		new_triangles = q.T;

		return OK;
	}

	Eigen::MatrixXi Region::get_triangulation(const Eigen::MatrixXi &F)
	{
		Eigen::MatrixXi F2(region_triangles.size(), 3);
		for (int j = 0; j < region_triangles.size(); j++)
		{
			F2.row(j) = F.row(region_triangles(j));
		}
		return F2;
	}

	
} // namespace cellogram

// -----------------------------------------------------------------------------


