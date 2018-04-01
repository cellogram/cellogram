////////////////////////////////////////////////////////////////////////////////
#include <cellogram/mesh_solver.h>
#include <cellogram/laplace_energy.h>
#include <cellogram/tri2hex.h>
#include <cellogram/vertex_degree.h>
#include <cellogram/region_grow.h>
#include <gurobi_solver/state.h>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

// -----------------------------------------------------------------------------
	typedef std::vector<int> Path;
	State s;

	void mesh_solver::find_bad_regions(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) {
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
			if (current_laplace_energy(i) > 0.8*avg || degree(i) != 6)
				//if (degree(i) != 6)
			{
				crit_pass[i] = true;
			}
		}

		// Find connected regions where the criterium was not passed
		region_grow(adjacency_list, crit_pass, region);

		// Find edges of connected region
		region_bounding(V, F, region, region_edges);

		// Line boundary
		bad_P1 = Eigen::MatrixXd(V.rows(), 3);
		bad_P2 = Eigen::MatrixXd(V.rows(), 3);
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

	Eigen::VectorXi find_n_neighbor(const Eigen::MatrixXi &F2,const Eigen::VectorXi &current_edges) {
		std::vector<int> internalTri(current_edges.rows(), 0);
		Eigen::VectorXi neigh;
		for (int j = 0; j < current_edges.rows(); j++)
		{
			for (int f = 0; f < F2.rows(); f++)
			{
				for (int lv = 0; lv < 3; lv++)
				{
					// count how many triangles in the region this vertex belongs to
					//std::cout << F2(f, lv) << " -- " << current_edges(j) << std::endl;
					if (F2(f, lv) == current_edges(j)) {
						internalTri[j]++;
					}
				}

			}
			neigh(j) = 6 - (internalTri[j] - 1);
		}
		return neigh;
	}

	Eigen::VectorXi find_internal_vertices(const Eigen::MatrixXi &F2, const Eigen::VectorXi &current_edges) {
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
		for (int i = 0; i < nVertices; i++)
		{
			for (int j = 0; j < current_edges.rows(); j++)
			{
				if (internal_vertex[i] == current_edges(j))
				{
					internal_vertex[i] = -1;
				}
			}

		}
		// save correct values into vI
		int k = 0;
		Eigen::VectorXi internal_vertex_eigen = Eigen::VectorXi::Zero(nVertices);
		for (int i = 0; i < nVertices; i++)
		{
			if (internal_vertex[i] != -1)
			{
				internal_vertex_eigen(k) = internal_vertex[i];
				k++;
			}
		}
		return internal_vertex_eigen;
	}

	void mesh_solver::solve_regions(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) {
		int nRegions = (int)region_edges.size();
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
			Eigen::VectorXi current_edges = Eigen::VectorXi::Zero(nPolygon);
			Eigen::VectorXi current_triangles = Eigen::VectorXi::Zero(nPolygon);
			for (int j = 0; j < nPolygon; j++)
			{
				current_edges(j) = region_edges[i][j];
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
			neigh = find_n_neighbor(F2, current_edges);
			
			// 3 Pepare for Gurobi
			// vB contains the boundary vertices, vI the internal vertices
			Eigen::MatrixXd vB(nPolygon, 2);
			for (int i = 0; i < nPolygon; i++)
			{
				vB.row(i) = V.row(current_edges[i]);
			}
			Eigen::VectorXi current_internal = find_internal_vertices(F2, current_edges);
			Eigen::MatrixXd vI(current_internal.size(), 2);
			for (int i = 0; i < current_internal.size(); i++)
			{
				vI.row(i) = V.row(current_internal[i]);
			}
			/*
			// Generate perfect mesh in ROI
			s.init(vB, vI, neigh);
			/*
			s.fill_hole();
			/*
			// Generate adjacency matrix and the laplacian
			q.adjacencyMatrix(s.F);
			q.laplacianMatrix();

			// Deriving Q and constraints for optimization
			q.QforOptimization(s.Vperfect, s.Vdeformed, 8);
			q.optimizationConstraints(s.V_boundary.rows());

			// Generate and solve model for gurobi
			g.model(q.Q, q.Aeq);

			// Map back to indices of coordinates
			q.mapBack(g.resultX);

			// Map q.T back to global indices
			VectorXi vGlobalInd = VectorXi::Zero(vBoundaryInd.size() + vInternalInd.size());
			for (size_t i = 0; i < vBoundaryInd.size(); i++)
			{
				vGlobalInd(i) = vBoundaryInd[i];
			}
			for (size_t i = 0; i < vInternalInd.size(); i++)
			{
				vGlobalInd(i + vBoundaryInd.size()) = vInternalInd[i];
			}

			tGlobal = MatrixXi(q.T.rows(), 3);
			for (size_t i = 0; i < q.T.rows(); i++)
			{
				for (size_t j = 0; j < 3; j++)
				{
					tGlobal(i, j) = vGlobalInd(q.T(i, j));
				}

			}
			*/
		}
	}

} // namespace cellogram

// -----------------------------------------------------------------------------


