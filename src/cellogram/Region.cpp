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


#include <igl/opengl/glfw/Viewer.h>

#include <igl/opengl/glfw/imgui/ImGuiMenu.h>

#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {
	namespace
	{
		bool is_clockwise(const Eigen::MatrixXd &poly) {
			double s = 0;
			for (int i = 0; i < poly.rows(); ++i) {
				int j = (i + 1) % poly.rows();
				double xi = poly(i, 0), yi = poly(i, 1);
				double xj = poly(j, 0), yj = poly(j, 1);
				s += (xj - xi) * (yj + yi);
			}
			return (s > 0);
		}

		bool is_neigh_valid(const VectorXi &neigh)
		{
			// check if number of neigh is smaller than 7
			if (neigh.maxCoeff() > 6)
			{
				return false;
			}

			// check if neigh actually closes a loop by tacing
			Eigen::VectorXi turns = neigh.array() - 4;
			std::cout << neigh << std::endl;
			//std::cout << turns << std::endl;

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

		Eigen::VectorXi find_n_neighbor(const Eigen::MatrixXi &F2, const Eigen::VectorXi &current_edge_vertices, const std::vector<std::vector<int>> &adj) {
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
				//if (neigh(j) > 6) {
				//	std::cout << "\nCheck neigh:\n" << current_edge_vertices(j) << " - " << adj[current_edge_vertices(j)].size() << " " << internalTri[j];
				//	std::cout << "\n\nF2: " << F2.transpose();
				//	std::cout << "\n\nEdge: " << current_edge_vertices.transpose();
				//}
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

	void Region::find_triangles(const Eigen::MatrixXi &F)
	{
		std::vector<int> region_faces;
		for (int j = 0; j < region_interior.rows(); j++)
		{
			int current_vertex = region_interior(j);
			{
				// for this vertex find the all the faces
				for (int f = 0; f < F.rows(); ++f) {
					for (int lv = 0; lv < F.cols(); ++lv) {
						if (current_vertex == F(f, lv))
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

	void Region::bounding(const Eigen::MatrixXi &F, const Eigen::MatrixXd &V) {
		boundary_loop(get_triangulation(F), region_boundary);

		Eigen::MatrixXd points(region_boundary.size(), 2);

		for (int i = 0; i < region_boundary.size(); ++i)
		{
			points(i, 0) = V(region_boundary(i), 0);
			points(i, 1) = V(region_boundary(i), 1);
		}

		if (is_clockwise(points))
			region_boundary.reverseInPlace();
	}


	int Region::resolve(const Eigen::MatrixXd &V_detected, const Eigen::MatrixXd &V_current, const Eigen::MatrixXi &F, const std::vector<std::vector<int>> &adj, const int perm_possibilities, Eigen::MatrixXd  &new_points, Eigen::MatrixXi &new_triangles)
	{
		int nPolygon = region_boundary.size(); // length of current edge

		// update triangle indices in region_triangles
		std::cout << "\n1 - Tri: \n" << region_triangles.transpose();
		find_triangles(F);
		std::cout << "\n2 - Tri: \n" << region_triangles.transpose();
		// for each region the belonging triangles need to be extracted and duplicates removed
		
		Eigen::MatrixXi F2 = get_triangulation(F);
		
		// Determine number of connections into the cluster to determine "neigh"
		Eigen::VectorXi neigh = find_n_neighbor(F2, region_boundary, adj);
		std::cout << "\n3 - Neigh: \n" << neigh.transpose();
		//PLEASE USE ME 
		//igl::opengl::glfw::Viewer viewer;
		//igl::opengl::glfw::imgui::ImGuiMenu menu;
		//viewer.plugins.push_back(&menu);

		//viewer.data().set_mesh(V, F2);
		//for (int i = 0; i < region_boundary.size(); ++i) {
		//	viewer.data().add_label(V.row(region_boundary(i)), std::to_string(neigh(i)));
		//}

		//viewer.launch();



		if (!is_neigh_valid(neigh)) {
			// this mesh is not properly close
			return NOT_PROPERLY_CLOSED;
		}



		// 3 Pepare for Gurobi
		// vB contains the boundary vertices, vI the internal vertices
		Eigen::MatrixXd vB(nPolygon, 2);
		Eigen::MatrixXd v_output(nPolygon,3);
		v_output.col(2).setZero();
		for (int i = 0; i < nPolygon; i++)
		{
			vB(i, 0) = V_detected(region_boundary(i), 0);
			vB(i, 1) = V_detected(region_boundary(i), 1);

			v_output(i, 0) = V_current(region_boundary(i), 0);
			v_output(i, 1) = V_current(region_boundary(i), 1);
		}
		Eigen::MatrixXd vI(region_interior.size(), 2);
		for (int i = 0; i < region_interior.size(); i++)
		{
			vI(i, 0) = V_detected(region_interior[i], 0);
			vI(i, 1) = V_detected(region_interior[i], 1);
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
		if (g.resultX(0) == -1) {
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
		
		std::cout << "\nDef: \n" << s.Vdeformed;
		std::cout << "\nPerf: \n" << s.Vperfect;
		std::cout << "\Vi: \n" << vI;

		new_points = permutation.transpose() * s.Vperfect;

		const Eigen::Vector3d bary_detected = new_points.block(0, 0, nPolygon, 3).colwise().mean();
		const Eigen::Vector3d bary_points = v_output.colwise().mean();
		const Eigen::Vector3d traslation = bary_points - bary_detected;

		for (long i = 0; i < new_points.rows(); ++i) {
			new_points.row(i) += traslation;
		}

		new_points.block(0, 0, nPolygon, 3) = v_output;

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


