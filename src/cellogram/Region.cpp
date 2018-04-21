////////////////////////////////////////////////////////////////////////////////
#include "cellogram/mesh_solver.h"
#include <cellogram/laplace_energy.h>
#include <cellogram/tri2hex.h>
#include <cellogram/vertex_degree.h>
#include <cellogram/region_grow.h>
#include <cellogram/Mesh.h>
#include <gurobi_solver/state.h>
#include <gurobi_solver/generateQ.h>
#include <gurobi_solver/gurobiModel.h>

#include <cellogram/navigation.h>
#include <cellogram/boundary_loop.h>
#include <cellogram/delaunay.h>
#include <cellogram/convex_hull.h>
#include <igl/opengl/glfw/Viewer.h>

#include <igl/triangle/triangulate.h>

#include <igl/opengl/glfw/imgui/ImGuiMenu.h>

#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/Timer.h>
#include <imgui/imgui.h>


#include <igl/opengl/glfw/Viewer.h>
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

		void removeRow(Eigen::VectorXi& matrix, unsigned int rowToRemove)
		{
			unsigned int numRows = matrix.rows() - 1;
			unsigned int numCols = matrix.cols();

			if (rowToRemove < numRows)
				matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.bottomRows(numRows - rowToRemove);

			matrix.conservativeResize(numRows, numCols);
		}

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
				neigh(j) = 6 - (internalTri[j] - 1);
				//neigh(j) = adj[current_edge_vertices(j)].size() - (internalTri[j] - 1);
				//if (neigh(j) > 6) {
				//	std::cout << "\nCheck neigh:\n" << current_edge_vertices(j) << " - " << adj[current_edge_vertices(j)].size() << " " << internalTri[j];
				//	std::cout << "\n\nF2: " << F2.transpose();
				//	std::cout << "\n\nEdge: " << current_edge_vertices.transpose();
				//}
			}
			return neigh;
		}
	}



	std::string Region::pretty_status(const int status)
	{
		switch (status)
		{
		case NOT_CHECKED: return "Not checked";
		case TOO_FEW_VERTICES: return "Too few vertices";
		case TOO_MANY_VERTICES: return "Too many vertices";
		case NOT_PROPERLY_CLOSED: return "Not properly closed";
		case NO_SOLUTION: return "No solution";
		case OK: return "Ok";
		default: return "";
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

	void Region::grow(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
	{
		// Turn boundary vertices to internal vertices and recalculate boundary
		// This function is called on individual regions an will not merge overlapping regions

		Eigen::VectorXi regions_id(V.rows());
		regions_id.setZero();

		//for any pt in region_interior or region_boundary set regions_id to 1
		for (size_t i = 0; i < region_interior.size(); i++)
		{
			regions_id(region_interior(i)) = 1;
		}
		for (size_t i = 0; i < region_boundary.size(); i++)
		{
			regions_id(region_boundary(i)) = 1;
		}

		find_points(regions_id, 1);
		find_triangles(F, regions_id, 1);

		bounding(F, V);
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
		assert(region_triangles.size() > 0);

		for (int j = 0; j < region_faces.size(); j++)
		{
			region_triangles(j) = region_faces[j];
		}
	}

	void Region::find_triangles(const Eigen::MatrixXi &F)
	{
		std::vector<std::vector<int>> vertex_to_tri(9999);

		for (int f = 0; f < F.rows(); ++f) {
			for (int lv = 0; lv < F.cols(); ++lv)
			{
				vertex_to_tri[F(f, lv)].push_back(f);
			}
		}

		std::vector<int> region_faces;
		for (int j = 0; j < region_interior.rows(); j++)
		{
			int current_vertex = region_interior(j);
			region_faces.insert(region_faces.end(), vertex_to_tri[current_vertex].begin(), vertex_to_tri[current_vertex].end());
			//// for this vertex find the all the faces
			//for (int f = 0; f < F.rows(); ++f) {
			//	for (int lv = 0; lv < F.cols(); ++lv) {
			//		if (current_vertex == F(f, lv))
			//		{
			//			region_faces.push_back(f);
			//			break;
			//		}
			//	}
			//}
		}

		std::sort(region_faces.begin(), region_faces.end());
		region_faces.erase(std::unique(region_faces.begin(), region_faces.end()), region_faces.end());

		region_triangles.resize(region_faces.size());
		assert(region_triangles.size() > 0);
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


	int Region::check_region(const Eigen::MatrixXd &V_detected, const Eigen::MatrixXd &V_current, const Eigen::MatrixXi &F, const std::vector<std::vector<int>> &adj)
	{
		int nPolygon = region_boundary.size(); // length of current edge

		// update triangle indices in region_triangles
		//find_triangles(F);
		// for each region the belonging triangles need to be extracted and duplicates removed

		Eigen::MatrixXi F2 = get_triangulation(F);

		// Determine number of connections into the cluster to determine "neigh"
		Eigen::VectorXi neigh = find_n_neighbor(F2, region_boundary, adj);

		////PLEASE USE ME 
		//igl::opengl::glfw::Viewer viewer;
		////igl::opengl::glfw::imgui::ImGuiMenu menu;
		////viewer.plugins.push_back(&menu);

		//viewer.data().set_mesh(V_detected, F2);
		//for (int i = 0; i < region_boundary.size(); ++i) {
		//	//viewer.data().add_label(V_detected.row(region_boundary(i)), std::to_string(neigh(i)));
		//}

		//viewer.launch();



		if (!is_neigh_valid(neigh)) {
			// this mesh is not properly close
			return NOT_PROPERLY_CLOSED;
		}



		// 3 Pepare for Gurobi
		// vB contains the boundary vertices, vI the internal vertices
		Eigen::MatrixXd vB(nPolygon, 2);
		Eigen::MatrixXd v_output(nPolygon, 3);
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
		s.init(vB, vI, neigh);
		s.fill_hole();

		//// compare sizes
		//std::cout << "\n\nregion_boundary: \n" << region_boundary.transpose() << "\n" <<region_boundary.rows();
		//std::cout << "\n\nregion_interior: \n" << region_interior.transpose() << "\n" << region_interior.rows();
		//std::cout << "\n\ns.Vdeformed: \n" << s.Vdeformed.transpose() << "\n" << s.Vdeformed.rows();
		//std::cout << "\n\ns.Vperfect: \n" << s.Vperfect.transpose() << "\n" << s.Vperfect.rows();
		////

		// Check whether vertices inside region is equal to the ones expected
		if (s.Vperfect.rows() > s.Vdeformed.rows()) {
			return TOO_FEW_VERTICES;
		}
		else if (s.Vperfect.rows() < s.Vdeformed.rows()) {
			return TOO_MANY_VERTICES;
		}
		return 0;
	}

	int Region::resolve(const Eigen::MatrixXd &V_detected, const Eigen::MatrixXd &V_current, const int perm_possibilities, Eigen::MatrixXd  &new_points, Eigen::MatrixXi &new_triangles)
	{
		if (status != 0)
			return status;

		int nPolygon = region_boundary.size(); // length of current edge

		Eigen::MatrixXd v_output(nPolygon, 3);
		v_output.col(2).setZero();
		for (int i = 0; i < nPolygon; i++)
		{
			v_output(i, 0) = V_current(region_boundary(i), 0);
			v_output(i, 1) = V_current(region_boundary(i), 1);
		}


		// Generate adjacency matrix and the laplacian
		generateQ q;
		gurobiModel g;
		q.adjacencyMatrix(s.F);
		q.laplacianMatrix();


		// Deriving Q and constraints for optimization
		int K = std::min(perm_possibilities,q.iNrV);
		q.QforOptimization(s.Vperfect, s.Vdeformed, K);
		q.optimizationConstraints(s.V_boundary.rows());


		// Generate and solve model for gurobi
		g.model(q.Q, q.Aeq);
		if (g.resultX(0) == -1) {
			// no solution found
			return NO_SOLUTION;
		}

		//std::cout << "\ng.resultX\n" << g.resultX << std::endl;
		//std::cout << "\nq.IDX\n" << q.IDX << std::endl;

		// Map back to indices of coordinates
		q.mapBack(g.resultX);
		//std::cout << g.resultX << std::endl;

		std::vector<Eigen::Triplet<double> > permutation_triplets;
		const int n_points = s.Vperfect.rows();
		for (int i = 0; i < n_points; ++i)
		{
			const int index = i * K;

			for (int j = 0; j < K; ++j)
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

		//std::cout << Eigen::MatrixXd(permutation) << std::endl;
		
		//std::cout << "\nDef: \n" << s.Vdeformed;
		//std::cout << "\nPerf: \n" << s.Vperfect;
		//std::cout << "\Vi: \n" << vI;

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

	bool Region::fix_missing_points(const Eigen::MatrixXi &F)
	{
		igl::Timer timer;  timer.start();
		std::vector<int> internal(region_interior.data(), region_interior.data() + region_interior.size());
		std::vector<int> boundary(region_boundary.data(), region_boundary.data() + region_boundary.size());

		for (int j = 0; j < region_triangles.size(); j++)
		{
			for (int i = 0; i < 3; ++i)
			{
				const int vId = F(region_triangles(j), i);
				//if (std::find(boundary.begin(), boundary.end(), vId) == boundary.end())
				if (std::count(boundary.begin(), boundary.end(), vId) == 0)
				{
					internal.push_back(vId);
				}
			}
		}


		std::sort(internal.begin(), internal.end());
		internal.erase(std::unique(internal.begin(), internal.end()), internal.end());

		if (internal.size() == region_interior.size())
			return false;

		region_interior.resize(internal.size());

		for (int j = 0; j < internal.size(); j++)
		{
			region_interior(j) = internal[j];
		}

		timer.stop();
		std::cout << " find missing took " << timer.getElapsedTime() << "s" << std::endl;

		timer.start();
		find_triangles(F);

		timer.stop();
		std::cout << " find_triangles took " << timer.getElapsedTime() << "s" << std::endl;
		return true;
	}

	void Region::delete_vertex(Mesh &mesh, const int index)
	{
		int n = region_interior.rows();
		Eigen::VectorXi V_new(n - 1);
		int i;
		for (i = 0; i < n; i++)
		{
			if (region_interior(i) == index)
			{
				removeRow(region_interior, i);
				break;
			}
		}


		const int n_poly_vertices = mesh.adj[index].size();
		Eigen::VectorXi local2global(n_poly_vertices);

		// retriangulate region
		Eigen::MatrixXd boundary_V(n_poly_vertices, mesh.points.cols());

		int f = mesh.vertex_to_tri[index][0];
		int lv;
		for (lv = 0; lv < 3; ++lv)
		{
			if (mesh.triangles(f, lv) == index)
				break;
		}
		assert(ln < 3);

		NavigationData data(mesh.triangles);
		NavigationIndex n_index = index_from_face(mesh.triangles, data, f, lv);

		for (int i = 0; i <n_poly_vertices; i++)
		{
			const int vertex_id = switch_vertex(data, n_index).vertex;
			local2global(i) = vertex_id;
			boundary_V.row(i) = mesh.points.row(vertex_id);
			n_index = next_around_vertex(data, n_index);
		}

		Eigen::MatrixXd dummy;
		Eigen::MatrixXi new_triangles;

		//triangulate_polygon(boundary_V, dummy, new_triangles);
		// TODO fill boundary_V 

		Eigen::MatrixXi E, WE;
		Eigen::VectorXi J;

		int nn = (int)boundary_V.rows();
		E.resize(nn, 2);
		for (int i = 0; i < nn; ++i) {
			E.row(i) << i, (i + 1) % nn;
		}
		Eigen::MatrixXd PV = boundary_V.leftCols<2>();


		//VpqYYS0
		igl::triangle::triangulate(PV, E, MatrixXd(0, 2), "QpqYYS0", dummy, new_triangles);
		//asd.conservativeResize(asd.rows(), 3);
		//asd.col(2).setZero();

		//std::cout << boundary_V << std::endl;
		//std::cout << new_triangles << std::endl;
		//std::cout << PV << std::endl;

		//igl::opengl::glfw::Viewer viewer;
		//viewer.data().set_mesh(boundary_V, new_triangles);
		//viewer.data().add_points(asd, eigen::rowvector3d(1, 0, 0));
		//viewer.data().add_points(pv, eigen::rowvector3d(1, 1, 0));

		//viewer.launch();


		mesh.local_update(local2global, index, new_triangles);
		find_triangles(mesh.triangles);
	}

	void Region::triangluate_region()
	{
	}

	void Region::local_to_global(Eigen::VectorXi & local2global)
	{
		const int n_region_boundary = region_boundary.size();
		const int n_region_interior = region_interior.size();
		//add this into the region
		local2global.resize(n_region_boundary + n_region_interior);

		for (int i = 0; i < n_region_boundary; i++)
		{
			const int global_index = region_boundary[i];
			local2global(i) = global_index;
		}
		for (int i = 0; i < n_region_interior; i++)
		{
			const int global_index = region_interior[i];
			local2global(i + n_region_boundary) = global_index;
		}

	}

	bool Region::add_orphaned_triangles(const Mesh &mesh)
	{
		Eigen::VectorXi n_internal_neighbors(region_boundary.size());
		n_internal_neighbors.setZero();
		for (int i = 0; i < region_boundary.size(); i++)
		{
			for (int neigh : mesh.adj[region_boundary(i)])
			{
				for (int j = 0; j < region_interior.size(); j++)
				{
					if (region_interior(j) == neigh)
					{
						n_internal_neighbors(i)++;
						break;
					}
				}
			}
		}

		std::vector<int> new_points;
		for (int i = 0; i < n_internal_neighbors.size(); i++)
		{
			if (n_internal_neighbors(i) == 4) {
				bool must_add = true;
				for (int j = 0; j < mesh.boundary.size(); j++) {
					if (region_boundary(i) == mesh.boundary(j))
					{
						must_add = false;
						break;
					}
				}
				if(must_add)
					new_points.push_back(region_boundary(i));
			}
		}
		if (new_points.empty())
			return false;


		for (int i = 0; i < region_interior.size(); ++i)
			new_points.push_back(region_interior(i));

		std::sort(new_points.begin(), new_points.end());
		new_points.erase(std::unique(new_points.begin(), new_points.end()), new_points.end());

		if (new_points.size() == region_interior.size())
			return false;

		region_interior.resize(new_points.size());
		for (int i = 0; i < region_interior.size(); ++i)
			region_interior(i) = new_points[i];


/*		for (int i : mesh.adj[660])
			std::cout << i << std::endl;
		for (int i = 0; i < n_internal_neighbors.size(); ++i)
			std::cout << region_boundary(i) << " " << n_internal_neighbors(i) << std::endl;

		std::cout << region_interior << std::endl;*/

		find_triangles(mesh.triangles);

		region_boundary.resize(0);
		bounding(mesh.triangles, mesh.points);

		return true;
	}

	
} // namespace cellogram

// -----------------------------------------------------------------------------


