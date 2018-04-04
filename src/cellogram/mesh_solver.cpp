////////////////////////////////////////////////////////////////////////////////
#include "cellogram/mesh_solver.h"
#include <cellogram/laplace_energy.h>
#include <cellogram/tri2hex.h>
#include <cellogram/vertex_degree.h>
#include <cellogram/region_grow.h>
#include <gurobi_solver/state.h>
#include <gurobi_solver/generateQ.h>
#include <gurobi_solver/gurobiModel.h>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

// -----------------------------------------------------------------------------
	typedef std::vector<int> Path;
	State s;
	generateQ q;
	gurobiModel g;
	
	void removeRow(Eigen::MatrixXi& matrix, unsigned int rowToRemove)
	{
		unsigned int numRows = matrix.rows() - 1;
		unsigned int numCols = matrix.cols();

		if (rowToRemove < numRows)
			matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.bottomRows(numRows - rowToRemove);

		matrix.conservativeResize(numRows, numCols);
	}

	void replaceTriangles(const Eigen::MatrixXi tNew, Eigen::MatrixXi &triangles, Path roiInternal)
	{
		// find any triangles that connect to the internal of ROI
		std::vector<int> removeIdx;

		for (size_t i = 0; i < triangles.rows(); i++)
		{
			for (size_t j = 0; j < roiInternal.size(); j++)
			{
				if ((triangles(i, 0) == roiInternal[j]) || (triangles(i, 1) == roiInternal[j]) || (triangles(i, 2) == roiInternal[j]))
				{
					removeIdx.push_back(i);
					std::cout << i << "\n";
				}
			}
		}

		// using default comparison:
		std::vector<int>::iterator it;
		it = std::unique(removeIdx.begin(), removeIdx.end());
		removeIdx.resize(std::distance(removeIdx.begin(), it));



		// remove all the triangles that connect to the internal of ROI
		for (int i = removeIdx.size() - 1; i >= 0; i--)
		{
			removeRow(triangles, removeIdx[i]);
		}

		// add new rows at the end of triangles
		Eigen::MatrixXi tmp = triangles;
		triangles.resize(triangles.rows() + tNew.rows(), 3);

		triangles << tmp, tNew;
	}

	void mesh_solver::find_bad_regions(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) {
		region_edges.clear();
		// Calculate the graph adjancency
		std::vector<std::vector<int>> adjacency_list; // adjaceny list of triangluar mesh
		tri2hex(F, adjacency_list);

		// Calculate the laplacian energy with respect to the original positions
		Eigen::VectorXd current_laplace_energy;
		laplace_energy(V, F, current_laplace_energy);

		// Find the degree of each vertex
		Eigen::VectorXi degree;
		vertex_degree(F, degree);

		// Determine whether each vertex passes the criterium for bad mesh
		double avg = current_laplace_energy.mean();
		std::vector<bool> crit_pass(V.rows(), false);
		for (int i = 0; i < V.rows(); i++)
		{
			if (current_laplace_energy(i) > 1.8*avg || degree(i) != 6)
				//if (degree(i) != 6)
			{
				crit_pass[i] = true;
			}
		}

		// Find connected regions where the criterium was not passed
		region_grow(adjacency_list, crit_pass, region);

		// Find edges of connected region
		region_bounding(F, region, region_edges);
	}

	void mesh_solver::compute_regions_edges(const Eigen::MatrixXd &V, Eigen::MatrixXd &bad_P1, Eigen::MatrixXd &bad_P2) {
		// Line boundary
		bad_P1.resize(V.rows(), 3);
		bad_P2.resize(V.rows(), 3);
		int k = 0;
		for (int i = 0; i < region_edges.size(); i++)
		{
			int n = (int)region_edges[i].size();
			for (int j = 0; j < region_edges[i].size(); j++) //
			{
				bad_P1.row(k) = V.row(region_edges[i][j]);
				bad_P2.row(k) = V.row(region_edges[i][(j + 1) % n]);
				k++;
			}
		}

		bad_P1.conservativeResize(k, 3);
		bad_P2.conservativeResize(k, 3);
	}


	std::vector<Path> find_points_in_region(const Eigen::VectorXi &region)
	{
		std::vector<Path> region_points(region.rows());
		for (int j = 0; j < region.rows(); j++)
		{
			// only use vertices with region(j) > 0, because those ==0 are considered good
			if (region(j) > 0) {
				region_points[region(j) - 1].push_back(j);
			}
		}
		return region_points;
	}


	std::vector<Path> find_triangles_in_region(const Eigen::MatrixXi &F, const Eigen::VectorXi &region)
	{
		std::vector<Path> region_F(region.rows());
		for (int j = 0; j < region.rows(); j++)
		{
			// only use vertices with region(j) > 0, because those ==0 are considered good
			if (region[j] > 0)
			{
				// for this vertex find the all the faces
				for (int f = 0; f < F.rows(); ++f) {
					for (int lv = 0; lv < F.cols(); ++lv) {
						if (j == F(f, lv))
						{
							region_F[region(j) - 1].push_back(f);
						}
					}
				}
			}
		}
		return region_F;
	}

	Eigen::VectorXi find_n_neighbor(const Eigen::MatrixXi &F2,const Eigen::VectorXi &current_edge_vertices, std::vector<Path> &adj) {
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

	Path find_internal_vertices(const Eigen::MatrixXi &F2, const Eigen::VectorXi &current_edges) {
		// add all vertices in F2....
		int nVertices = (int)F2.maxCoeff() + 1;
		Path internal_vertex;
		for (int i = 0; i < F2.rows(); i++)
		{
			for (int j = 0; j < 3; j++)
			{
				internal_vertex.push_back(F2(i, j));
			}
		}
		//... remove duplicates...
		std::sort(internal_vertex.begin(), internal_vertex.end());
		internal_vertex.erase(std::unique(internal_vertex.begin(), internal_vertex.end()), internal_vertex.end());

		//...and the vertices in current_edges
		for (int i = 0; i < internal_vertex.size(); i++)
		{
			for (int j = 0; j < current_edges.rows(); j++)
			{
				if (internal_vertex[i] == current_edges(j))
				{
					internal_vertex[i] = -1;
				}
			}

		}
		// remove all the -1 values
		std::vector<int>::iterator it = internal_vertex.begin();
		while (it != internal_vertex.end())
		{
			if (*it == -1)
			{
				it = internal_vertex.erase(it);
			}
			else
				it++;
		}
		return internal_vertex;
	}

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

	void mesh_solver::solve_regions(const Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
		int nRegions = (int)region_edges.size();
		int counter_invalid_neigh = 0;
		int counter_infeasible_region = 0;
		// 0. Calculate adj for later
		std::vector<Path> adj;
		tri2hex(F,adj);

		// 1.
		// create a 2d vector containing all the vertices belonging to each region
		// as well as one containing the triangles in that region
		
		//std::vector<Path> region_points(nRegions);
		
		region_edges.resize(nRegions);

		// 1.1 find the points and triangles inside
		std::vector<Path> region_points = find_points_in_region(region);
		std::vector<Path> region_F = find_triangles_in_region(F, region);
		

		// 2.
		// loop through all the boundaries and solve individually. Possibly skip first one, as it's the boundary of the image
		for (int i = 1; i < region_edges.size(); i++)
		{
			int nPolygon = (int)region_edges[i].size(); // length of current edge
			
			// 2.1
			// Save current edge points into single vector
			Eigen::VectorXi current_edge_vertices = Eigen::VectorXi::Zero(nPolygon);
			Eigen::VectorXi current_triangles = Eigen::VectorXi::Zero(nPolygon);
			for (int j = 0; j < nPolygon; j++)
			{
				current_edge_vertices(j) = region_edges[i][j];
			}

			//std::cout << "\nCurrent_edge \n" << current_edges.transpose() << std::endl;
			
			// 2.2
			// Determine number of connections into the cluster to determine "neigh"
			Eigen::VectorXi neigh(nPolygon);

			// for each region the belonging triangles need to be extracted and duplicates removed
			std::sort(region_F[i].begin(), region_F[i].end());
			region_F[i].erase(unique(region_F[i].begin(), region_F[i].end()), region_F[i].end());

			Eigen::MatrixXi F2(region_F[i].size(), 3);
			for (int j = 0; j < region_F[i].size(); j++)
			{
				F2.row(j) = F.row(region_F[i][j]);
			}
			
			// 2.3 find number of neighbors
			
			neigh = find_n_neighbor(F2, current_edge_vertices, adj);
			if (!is_neigh_valid(neigh)) {
				// this mesh is not properly close
				//std::cout << "Invalid neighbors\n" << neigh << std::endl;
				counter_invalid_neigh++;
				break;
			}

			

			// 3 Pepare for Gurobi
			// vB contains the boundary vertices, vI the internal vertices
			Eigen::MatrixXd vB(nPolygon, 2);
			for (int i = 0; i < nPolygon; i++)
			{
				vB(i, 0) = V(current_edge_vertices[i], 0);
				vB(i, 1) = V(current_edge_vertices[i], 1);
			}
			Path current_internal = find_internal_vertices(F2, current_edge_vertices);
			Eigen::MatrixXd vI(current_internal.size(), 2);
			for (int i = 0; i < current_internal.size(); i++)
			{
				vI(i, 0) = V(current_internal[i], 0);
				vI(i, 1) = V(current_internal[i], 1);
			}
			// Generate perfect mesh in ROI
			s.init(vB, vI, neigh);
			
			s.fill_hole();

			// Check whether vertices inside region is equal to the ones expected
			if (s.Vdeformed.rows() != s.Vdeformed.rows()) {
				std::cout << "\nDeformed\n" << s.Vdeformed.transpose();
				std::cout << "\nPerfect\n" << s.Vperfect.transpose();
			}

			// Generate adjacency matrix and the laplacian
			q.adjacencyMatrix(s.F);
			q.laplacianMatrix();

			
			// Deriving Q and constraints for optimization
			q.QforOptimization(s.Vperfect, s.Vdeformed, 8);
			q.optimizationConstraints(s.V_boundary.rows());

			
			// Generate and solve model for gurobi
			g.model(q.Q, q.Aeq);
			if (g.resultX(1) == -1) {
				// no solution found
				counter_infeasible_region++;
				break;
			}
			
			// Map back to indices of coordinates
			q.mapBack(g.resultX);

			
			// Map q.T back to global indices
			Eigen::VectorXi vGlobalInd = Eigen::VectorXi::Zero(current_edge_vertices.size() + current_internal.size());
			for (int i = 0; i < current_edge_vertices.size(); i++)
			{
				vGlobalInd(i) = current_edge_vertices[i];
			}
			for (int i = 0; i < current_internal.size(); i++)
			{
				vGlobalInd(i + current_edge_vertices.size()) = current_internal[i];
			}

			Eigen::MatrixXi tGlobal = Eigen::MatrixXi(q.T.rows(), 3);
			for (int i = 0; i < q.T.rows(); i++)
			{
				for (int j = 0; j < 3; j++)
				{
					tGlobal(i, j) = vGlobalInd(q.T(i, j));
				}

			}
			
			replaceTriangles(tGlobal, F, current_internal);


		}
		int nRegions_true = nRegions - 1; // because the first one is skipped, as it is the boundary of the image
		std::cout << "\n\nOptmization complete" << std::endl;
		std::cout << "Total regions:\t" << nRegions_true << std::endl;
		std::cout << "Solved:\t\t" << nRegions_true - counter_invalid_neigh - counter_infeasible_region << std::endl;
		std::cout << "Bad loop:\t" << counter_invalid_neigh << std::endl;
		std::cout << "Infeasible:\t" << counter_infeasible_region << std::endl;
	}

} // namespace cellogram

// -----------------------------------------------------------------------------


