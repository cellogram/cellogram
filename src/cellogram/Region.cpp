////////////////////////////////////////////////////////////////////////////////
#include "Region.h"
#include <cellogram/mesh_solver.h>
#include <cellogram/boundary_loop.h>
#include <cellogram/convex_hull.h>
#include <cellogram/delaunay.h>
#include <cellogram/dijkstra.h>
#include <cellogram/laplace_energy.h>
#include <cellogram/Mesh.h>
#include <cellogram/navigation.h>
#include <cellogram/region_grow.h>
#include <cellogram/State.h>
#include <cellogram/tri2hex.h>
#include <cellogram/vertex_degree.h>
#ifdef CELLOGRAM_WITH_GUROBI
#include <gurobi_solver/generateQ.h>
#include <gurobi_solver/gurobiModel.h>
#include <gurobi_solver/state.h>
#endif
#include <igl/dijkstra.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/sort.h>
#include <igl/Timer.h>
#include <igl/triangle/triangulate.h>
#include <imgui/imgui.h>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {
	namespace
	{
		double polygon_area(Eigen::VectorXd X, Eigen::VectorXd Y) {
			assert(X.size() == Y.size());
			double area = 0;
			int j = X.size()-1;
			for (int i = 0; i < X.size(); i++)
			{
				area += (X(i) + X(j))*(Y(i) - Y(j));
				j = i;
			}
			return 0.5 * std::abs(area);
		}

		double inline det(const Eigen::Vector2d &u, const Eigen::Vector2d &v) {
			Eigen::Matrix2d M;
			M.col(0) = u;
			M.col(1) = v;
			return M.determinant();
		}
		// Return true iff [a,b] intersects [c,d], and store the intersection in ans
		bool intersect_segment(const Eigen::Vector2d &a, const Eigen::Vector2d &b, const Eigen::Vector2d &c, const Eigen::Vector2d &d) {
			const double eps = 1e-10; // small epsilon for numerical precision
			double x = det(c - a, d - c);
			double y = det(b - a, a - c);
			double z = det(b - a, d - c);
			// ab and cd are parallel ||
			if (std::abs(z) < eps || x * z < 0 || x * z > z*z || y * z < 0 || y * z > z*z) return false;
			return true;
		}

		bool is_inside(const Eigen::MatrixXd &poly, const Eigen::Vector2d &query) {
			Eigen::Vector2d outside(-1000, -1000);
			size_t n = poly.rows();
			bool tmp, ans = false;
			for (size_t i = 0; i < poly.rows(); ++i) {
				const size_t ip = (i + 1) % n;
				Eigen::Vector2d p0(poly(i, 0), poly(i, 1));
				Eigen::Vector2d p1(poly(ip, 0), poly(ip, 1));
				tmp = intersect_segment(query, outside, p0, p1);
				ans = (ans != tmp);
			}
			return ans;
		}

		void flip_triangles(const Eigen::MatrixXd &pts, const Eigen::MatrixXi &tri_original,  std::vector<std::vector<int>>& adj)
		{
			int lv = 0;

			NavigationData data(tri_original);

			for (int f = 0; f < tri_original.rows(); f++)
			{
				NavigationIndex index = index_from_face(tri_original, data, f, lv);

				for (int i = 0; i < 3; i++)
				{
					if (switch_face(data, index).face < 0)
						continue;
					const int v0 = index.vertex;
					const int v1 = switch_vertex(data, index).vertex;
					const int v2 = switch_vertex(data, switch_edge(data, index)).vertex;
					const int v3 = switch_vertex(data, switch_edge(data, switch_face(data, index))).vertex;
					const double d0 = (pts.row(v0) - pts.row(v1)).norm();
					const double d2 = (pts.row(v2) - pts.row(v3)).norm();

					if (std::abs(d0 - d2) / d0 < 0.5)
						adj[v2].push_back(v3);

					index = next_around_face(data, index);
				}
			}
		}

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

		bool is_neigh_valid(const Eigen::VectorXi &neigh)
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

		//struct edge {
		//	int x, y, w;
		//	edge(void) : x(), y(), w() {}
		//	edge(int a, int b, int c) : x(a), y(b), w(c) {}
		//	bool operator< (const edge &e) const {
		//		return w > e.w; // Extract min-cost edges first
		//	}
		//};
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
		case REGION_TOO_LARGE: return "Region too large";
		case OK: return "Ok";
		default: return "";
		}
	}


	int Region::size()
	{
		return region_interior.size();
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
		assert(region_triangles.size() > 0 && "Region size is zero, check for redundant vertices");

		for (int j = 0; j < region_faces.size(); j++)
		{
			region_triangles(j) = region_faces[j];
		}
	}

	void Region::find_triangles(const Eigen::MatrixXi &F, bool force_add_boundaries)
	{
		std::vector<std::vector<int>> vertex_to_tri(9999);

		std::vector<int> region_faces;

		for (int f = 0; f < F.rows(); ++f) {
			int counter = 0;
			for (int lv = 0; lv < F.cols(); ++lv)
			{
				const int vid = F(f, lv);
				if(force_add_boundaries && lv == counter)
					counter += (region_boundary.array() - vid).abs().minCoeff() == 0 ? 1 : 0;
				vertex_to_tri[vid].push_back(f);
			}

			if (counter == F.cols())
				region_faces.push_back(f);
		}

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


	void Region::check_region(const Eigen::MatrixXd &V_detected, const Eigen::MatrixXd &V_current, const Eigen::MatrixXi &F, const std::vector<std::vector<int>> &adj)
	{

		// update triangle indices in region_triangles
		//find_triangles(F);
		// for each region the belonging triangles need to be extracted and duplicates removed

		Eigen::MatrixXi F2 = get_triangulation(F);
		check_region_from_local_tris(V_detected, V_current, F2, adj);
	}

	void Region::check_region_from_local_tris(const Eigen::MatrixXd &V_detected, const Eigen::MatrixXd &V_current, const Eigen::MatrixXi &F2, const std::vector<std::vector<int>> &adj)
	{
#ifdef CELLOGRAM_WITH_GUROBI
		int nPolygon = region_boundary.size(); // length of current edge


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
			status = NOT_PROPERLY_CLOSED;
			return;
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
		//if (s.Vdeformed.rows() == 63)
		//{
		//	std::cout << "\n\nregion_triangles: \n" << region_triangles.transpose() << std::endl;
		//	std::cout << "\n\nregion_boundary: \n" << region_boundary.transpose() << "\n\n#:" << region_boundary.rows() << std::endl;
		//	std::cout << "\n\nregion_boundary: \n" << region_boundary.transpose() << "\n\n#:" <<region_boundary.rows() << std::endl;
		//	std::cout << "\n\nregion_interior: \n" << region_interior.transpose() << "\n\n#:" << region_interior.rows() << std::endl;
		//	std::cout << "\n\ns.Vdeformed: \n" << s.Vdeformed.transpose() << "\n\n#:" << s.Vdeformed.rows() << std::endl;
		//	std::cout << "\n\ns.Vperfect: \n" << s.Vperfect.transpose() << "\n\n#:" << s.Vperfect.rows() << std::endl;
		//}
		////

		// Check whether vertices inside region is equal to the ones expected
		if (s.Vperfect.rows() > s.Vdeformed.rows()) {
			points_delta = s.Vdeformed.rows() - s.Vperfect.rows();
			status = TOO_FEW_VERTICES;
			return;
		}
		else if (s.Vperfect.rows() < s.Vdeformed.rows()) {
			points_delta = s.Vdeformed.rows() - s.Vperfect.rows();
			status = TOO_MANY_VERTICES;
			return;
		}

		if (region_interior.size() > cellogram::State::max_region_vertices)
		{
			status = REGION_TOO_LARGE;
			return;
		}

		//// temporary to save vB,vI,neigh for comparing solvers
		//if (vI.rows() > 10)
		//{
		//	int nError;
		//	std::string path = "C:/Users/Tobias/Dropbox/NY/test_regions_" + std::to_string(std::rand() % 100);
		//	std::wstring widestr = std::wstring(path.begin(), path.end());
		//	nError = _wmkdir(widestr.c_str()); // can be used on Windows

		//	{
		//		std::ofstream vB_path(path + "/vB.vert");
		//		vB_path << vB << std::endl;
		//		vB_path.close();
		//	}
		//	{
		//		std::ofstream vI_path(path + "/vI.vert");
		//		vI_path << vI << std::endl;
		//		vI_path.close();
		//	}
		//	{
		//		std::ofstream neigh_path(path + "/neigh.txt");
		//		neigh_path << neigh << std::endl;
		//		neigh_path.close();
		//	}
		//}
		status = 0;
#endif
	}

	void Region::resolve(const Eigen::MatrixXd &V_detected, const Eigen::MatrixXd &V_current, const int perm_possibilities, const double gurobi_time_limit, Eigen::MatrixXd  &new_points, Eigen::MatrixXi &new_triangles, bool force_solve)
	{
		if (status != OK)
		{
			if(!force_solve || status != REGION_TOO_LARGE)
				return;
		}

#ifdef CELLOGRAM_WITH_GUROBI

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
		g.model(q.Q, q.Aeq, gurobi_time_limit);
		if (g.resultX(0) == -1) {
			// no solution found
			status = NO_SOLUTION;
			return;
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
#endif
	}

	bool Region::split_region(Mesh &mesh, const Eigen::Vector2i & split_end_points, Region &r1, Region  &r2)
	{
		Eigen::VectorXi path;
		find_split_path(mesh, split_end_points, path);

		//Eigen::VectorXi boundary1, boundary2;
		find_split_boundaries(split_end_points, path, r1.region_boundary, r2.region_boundary);

		find_interior_V(mesh, r1.region_boundary, r1.region_interior);
		find_interior_V(mesh, r2.region_boundary, r2.region_interior);

		Eigen::MatrixXi new_triangles_r1, new_triangles_r2;
		r1.triangluate_region(mesh, new_triangles_r1);
		r2.triangluate_region(mesh, new_triangles_r2);

		mesh.update_triangles_from_split(region_triangles, new_triangles_r1, new_triangles_r2);

		find_triangles(mesh.triangles);

		//Eigen::MatrixXi region_triangles_copy = region_triangles;
		bool ok = false;
		do {
			ok = r1.clean_up_boundary(mesh.triangles, region_triangles);
		} while (ok);
		do {
			ok = r2.clean_up_boundary(mesh.triangles, region_triangles);
		} while (ok);

		find_interior_V(mesh, r1.region_boundary, r1.region_interior);
		find_interior_V(mesh, r2.region_boundary, r2.region_interior);

		std::vector<std::vector<int>> adj;

		if (r1.size() > 0) {
			r1.find_triangles(mesh.triangles, true);
			auto tmp = r1.get_triangulation(mesh.triangles);
			adjacency_list(tmp, adj);
			r1.check_region_from_local_tris(mesh.detected, mesh.points, tmp, adj);
		}

		if (r2.size() > 0) {
			r2.find_triangles(mesh.triangles, true);
			auto tmp = r2.get_triangulation(mesh.triangles);
			adjacency_list(tmp, adj);
			r2.check_region_from_local_tris(mesh.detected, mesh.points, tmp, adj);
		}

		//if (something)
		//	return false;


		//r1.find_triangles(mesh.triangles);
		//std::cout << r1.region_triangles << std::endl;

		//r2.find_triangles(mesh.triangles);

		//
		//igl::opengl::glfw::Viewer viewer;
		//int mmax = mesh.points.maxCoeff();
		//Eigen::MatrixXd V(4, 3);
		//V <<
		//	0, 0, 0,
		//	mmax, 0, 0,
		//	mmax, mmax, 0,
		//	0, mmax, 0;

		//viewer.core.align_camera_center(V);
		//viewer.data().show_overlay_depth = false;

		//Eigen::MatrixXd asd1, asd2;
		//compute_edges(mesh.points, asd1, asd2);
		//viewer.data().point_size = 10;
		//viewer.data().add_points(asd1, Eigen::RowVector3d(0, 0, 1));

		//viewer.data().point_size = 8;
		//r1.compute_edges(mesh.points, asd1, asd2);
		//viewer.data().add_edges(asd1, asd2, Eigen::RowVector3d(1, 0, 0));
		//viewer.data().add_points(asd1, Eigen::RowVector3d(1, 0, 0));
		//r2.compute_edges(mesh.points, asd1, asd2);
		//viewer.data().add_edges(asd1, asd2, Eigen::RowVector3d(0, 1, 0));
		//viewer.data().add_points(asd1, Eigen::RowVector3d(0, 1, 0));
		//Eigen::MatrixXi asd3(new_triangles_r1.rows(), 3);// +new_triangles_r2.rows(), 3);
		//asd3 << new_triangles_r1;//, new_triangles_r2;
		//viewer.data().set_mesh(mesh.points, asd3);
		//viewer.data().add_points(asd1.row(0), Eigen::RowVector3d(0, 1, 1));

		//viewer.data().point_size = 10;
		//for (int i = 0; i < region_interior.size(); ++i) {
		//	viewer.data().add_points(mesh.points.row(region_interior(i)), Eigen::RowVector3d(0, 0, 1));
		//}

		//viewer.data().point_size = 8;
		//for (int i = 0; i < r1.region_interior.size(); ++i) {
		//	//viewer.data().add_points(mesh.points.row(r1.region_interior(i)), Eigen::RowVector3d(1, 0, 0));
		//}
		//for (int i = 0; i < r2.region_interior.size(); ++i) {
		//	viewer.data().add_points(mesh.points.row(r2.region_interior(i)), Eigen::RowVector3d(1, 1, 0));
		//}
		////viewer.data().set_mesh(boundary_V, new_triangles);
		////viewer.data().add_points(asd, eigen::rowvector3d(1, 0, 0));
		////viewer.data().add_points(pv, eigen::rowvector3d(1, 1, 0));
		//viewer.launch();
		//

		return true;
	}

	void Region::find_split_path(const Mesh &mesh, const Eigen::Vector2i & split_end_points, Eigen::VectorXi & path)
	{
		Eigen::MatrixXi tri(region_triangles.rows(),3);
		for (int i = 0; i < region_triangles.rows(); i++)
		{
			tri.row(i) = mesh.triangles.row(region_triangles(i));
		}
		//Eigen::MatrixXi flipped_region_triangles;

		std::vector<std::vector<int>> adj;
		adjacency_list(tri, adj);
		flip_triangles(mesh.points, tri, adj);

		std::set<int> target;
		target.insert(split_end_points(1));

		//Eigen::VectorXd min_distance;
		//Eigen::VectorXi previous;

		//igl::dijkstra_compute_paths(split_end_points(0), target, adj, min_distance, previous);

		//std::vector<int> path_tmp;
		//igl::dijkstra_get_shortest_path_to(split_end_points(1), previous, path_tmp);

		std::vector<std::vector<edge>> graph;
		graph.resize(adj.size());
		for (int i = 0; i < adj.size(); i++)
		{
			for (int j = 0; j < adj[i].size(); j++)
			{
				double dx = mesh.points(i, 0) - mesh.points(adj[i][j], 0);
				double dy = mesh.points(i, 1) - mesh.points(adj[i][j], 1);
				graph[i].push_back(edge(i, adj[i][j], std::sqrt(dx * dx + dy * dy)));
			}
		}

		std::vector<int> prev, dist;
		dijkstra(graph, prev, dist, split_end_points(0));

		std::vector<int> path_tmp;
		path_tmp.push_back(split_end_points(1));
		do {
			path_tmp.push_back(prev[path_tmp.back()]);
		} while (path_tmp.back() != prev[path_tmp.back()]);

		path.resize(path_tmp.size());
		for (int i = 0; i < path_tmp.size(); i++)
		{
			path(i) = path_tmp[i];
		}
	}

	void Region::find_split_boundaries(const Eigen::Vector2i & split_end_points, const Eigen::VectorXi & path, Eigen::VectorXi & boundary1, Eigen::VectorXi & boundary2)
	{
		//find position of split in boundary
		Eigen::Vector2i split_ind;
		for (int i = 0; i < region_boundary.size(); i++)
		{
			if (region_boundary(i) == split_end_points(0))
				split_ind(0) = i;
			else if (region_boundary(i) == split_end_points(1))
				split_ind(1) = i;
		}

		//make sure that split0 comes before split1
		int dir = 1;
		if (split_ind(0) > split_ind(1))
		{
			int tmp = split_ind(0);
			split_ind(0) = split_ind(1);
			split_ind(1) = tmp;
			dir = -1;
		}

		int l1 = split_ind(1) - split_ind(0) - 1;
		int l2 = region_boundary.size() - (split_ind(1) - split_ind(0) + 1);
		boundary1.resize(l1 + path.size());
		boundary2.resize(l2 + path.size());

		for (int i = 0; i < l1; i++)
		{
			boundary1(i) = region_boundary(split_ind(0) + i + 1);
		}
		for (int i = 0; i < l2; i++)
		{
			boundary2(i) = region_boundary((split_ind(1) + i + 1) % region_boundary.size());
		}
		for (int i = 0; i < path.size(); i++)
		{
			if (dir == 1)
			{
				boundary1(l1 + i) = path(i);
				boundary2(l2 + i) = path(path.size() - 1 - i);
			}
			else if (dir == -1)
			{
				boundary1(l1 + i) = path(path.size() - 1 - i);
				boundary2(l2 + i) = path(i);
			}
		}

	}

	void Region::find_interior_V(const Mesh & mesh, const Eigen::VectorXi & boundary, Eigen::VectorXi & interior_V)
	{
		if (boundary.size() <= 3)
			return;

		Eigen::MatrixXd poly(boundary.size(), 2);


		for (int i = 0; i < boundary.size(); i++)
		{
			poly(i, 0) = mesh.points(boundary(i), 0);
			poly(i, 1) = mesh.points(boundary(i), 1);
		}


		std::vector<int> tmp;
		for (int i = 0; i < region_interior.size(); i++)
		{
			if ((boundary.array() - region_interior(i)).abs().minCoeff() == 0)
				continue;
			Eigen::Vector2d point;
			point(0) = mesh.points(region_interior(i), 0);
			point(1) = mesh.points(region_interior(i), 1);
			if (is_inside(poly, point))
			{
				tmp.push_back(region_interior(i));
			}
		}

		interior_V.resize(tmp.size());
		for (int i = 0; i < tmp.size(); i++)
		{
			interior_V(i) = tmp[i];
		}
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
		//igl::Timer timer;  timer.start();
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

		//timer.stop();
		//std::cout << " find missing took " << timer.getElapsedTime() << "s" << std::endl;

		//timer.start();
		find_triangles(F);

		//timer.stop();
		//std::cout << " find_triangles took " << timer.getElapsedTime() << "s" << std::endl;
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
		assert(lv < 3);

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
		igl::triangle::triangulate(PV, E, Eigen::MatrixXd(0, 2), "QpqYYS0", dummy, new_triangles);
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

	void Region::triangluate_region(const Mesh mesh, Eigen::MatrixXi &new_triangles)
	{
		const int nn = (int)region_boundary.size();

		if (nn <= 3)
			return;

		Eigen::VectorXi local2global;
		local_to_global(local2global);

		//std::cout << local2global << std::endl;

		Eigen::MatrixXd pts(local2global.size(), 2);
		for (int i = 0; i < local2global.size(); i++)
		{
			const int global_index = local2global[i];
			pts(i, 0) = mesh.points(global_index, 0);
			pts(i, 1) = mesh.points(global_index, 1);
		}

		Eigen::MatrixXi E;
		E.resize(nn, 2);
		for (int i = 0; i < nn; ++i) {
			E.row(i) << i, (i + 1) % nn;
		}


		Eigen::VectorXi VM(local2global.size());
		VM.setZero();
		VM.block(0, 0, nn, 1).setConstant(1);

		//std::cout << VM << std::endl;
		Eigen::VectorXi EM;
		Eigen::MatrixXd dummy;
		Eigen::VectorXi dummy1, dummy2;
		Eigen::MatrixXi tmp;
		igl::triangle::triangulate(pts, E, Eigen::MatrixXd(0, 2), VM, EM, "QpqYYS0", dummy, tmp, dummy1, dummy2);

		new_triangles.resizeLike(tmp);

		for (int i = 0; i < tmp.rows(); i++)
		{
			for (int j = 0; j < tmp.cols(); j++)
			{
				new_triangles(i, j) = local2global(tmp(i, j));
			}
		}
	}

	bool Region::clean_up_boundary(const Eigen::MatrixXi & tri,const Eigen::MatrixXi & tri_list)
	{
		// getting region triangles and sorting them row-wise
		Eigen::MatrixXi tri_tmp(tri_list.size(),3), region_tri_sorted, dummy;
		for (int i = 0; i < tri_list.size(); i++)
		{
			tri_tmp.row(i) = tri.row(tri_list(i));
		}

		// sorting necessary
		igl::sort(tri_tmp, 2, 1, region_tri_sorted, dummy);


		const int n = region_boundary.size();
		for (int split_ind = 0; split_ind < n; split_ind++)
		{
			//find adjacent points
			Eigen::MatrixXi adj(1, 3), adj_sorted(1, 3);

			adj = Eigen::RowVector3i(region_boundary((split_ind - 1 + n) % n), region_boundary(split_ind), region_boundary((split_ind + 1) % n));
			igl::sort(adj, 2, 1, adj_sorted, dummy);

			//check if adjacent points create a triangle in tri
			std::vector<int> ind_del;
			for (int i = 0; i < region_tri_sorted.rows(); i++)
			{
				if ((adj_sorted - region_tri_sorted.row(i)).cwiseAbs().sum() == 0)
				{
					ind_del.push_back(i);
				}
			}

			if (ind_del.empty())
				continue;

			//std::cout << "\n\nBefore: \n" << tri_list << std::endl;
			for (int i = ind_del.size() - 1; i >= 0; i--)
			{
				// remove from tri list, so that it does not get removed from mesh.triangles later
				//removeRow(tri_list, ind_del[i]);

				// remove from region_boundary
				removeRow(region_boundary, split_ind);
			}
			//std::cout << "\n\nAfter: \n" << tri_list << std::endl;
			return true;
		}

		return false;
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

	bool Region::fix_pinched_region(const Mesh & mesh)
	{
		if (region_boundary.size() < 7)
			return false;

		int p1 = 0, p2 = 0;
		// loop through current boundary and check for triangles that have 
		// two vertices on the boundary that are not consecutive in the boundary
		for (int i = 1; i < region_boundary.size(); i++)
		{
			auto current_tri = mesh.vertex_to_tri[region_boundary(i)];
			
			// for each current_tri check whether it is somewhere else in the boundary
			for (int j = 0; j < region_boundary.size(); j++)
			{
				if (std::abs(i - j) <= 1)
					continue;
				if (std::abs(i - j) == region_boundary.size()-1)
					continue;

				auto next_tri = mesh.vertex_to_tri[region_boundary(j)];
				for (int m = 0; m < current_tri.size(); m++)
				{
					for (int n = 0; n < next_tri.size(); n++)
					{
						if ((current_tri[m] - next_tri[n]) == 0)
						{
							// if this is reached, it means that a triangle has been found
							// with two boundary vertices in non-consective positions in the boundary
							p1 = i;
							p2 = j;
							break;
						}
					}
				}
				
			}
		}
		if (p1 == p2)
			return false;

		if (p1 > p2)
		{
			int tmp = p1;
			p1 = p2;
			p2 = tmp;
		}

		// split the boundary
		Eigen::VectorXi b1(p1 - p2 + 1 + region_boundary.size()), b2(p2 - p1 + 1);

		b1.segment(0, p1+1) = region_boundary.segment(0, p1+1);
		b1.segment(p1 + 1, region_boundary.size() - p2) = region_boundary.segment(p2, region_boundary.size() - p2);

		b2.segment(0, p2 - p1 + 1) = region_boundary.segment(p1, p2 - p1 + 1);
		std::cout << "pinch found:\n" << std::endl;
		////std::cout << "points\n" << mesh.points.transpose() << std::endl;
		//std::cout << "boundary\n" << region_boundary.transpose() << std::endl;
		//std::cout << "b1\n" << b1.transpose() << std::endl;
		//std::cout << "b2\n" << b2.transpose() << std::endl;
		//std::cout << "\n" << std::endl;

		// determine which one to keep based on area
		Eigen::VectorXd X1(b1.size()), Y1(b1.size());
		for (int i = 0; i < X1.size(); i++)
		{
			X1(i) = mesh.points(b1(i),0);
			Y1(i) = mesh.points(b1(i),1);
		}
		double A1 = polygon_area(X1,Y1);

		Eigen::VectorXd X2(b2.size()), Y2(b2.size());
		for (int i = 0; i < X2.size(); i++)
		{
			X2(i) = mesh.points(b2(i), 0);
			Y2(i) = mesh.points(b2(i), 1);
		}
		double A2 = polygon_area(X2,Y2);
		
		std::cout << A1 << " " << A2 << std::endl;

		Eigen::VectorXi new_boundary;
		if (A1 > A2)
			new_boundary = b1;
		else
			new_boundary = b2;

		// overwrite current region
		region_boundary = new_boundary;

		Eigen::VectorXi tmp;
		find_interior_V(mesh, new_boundary, tmp);
		region_interior = tmp;
		if(tmp.size() > 0)
		{
			find_triangles(mesh.triangles, false);
		}

		return true;
	}


} // namespace cellogram

// -----------------------------------------------------------------------------


