////////////////////////////////////////////////////////////////////////////////
#include "Mesh.h"
#include <cellogram/State.h>
#include <cellogram/delaunay.h>
#include <cellogram/load_points.h>
#include <cellogram/convex_hull.h>
#include <cellogram/delaunay.h>
#include <cellogram/tri2hex.h>
#include <cellogram/voronoi.h>
#include <cellogram/vertex_degree.h>
#include <cellogram/laplace_energy.h>
#include <points_untangler/points_untangler.h>
#include <igl/list_to_matrix.h>
#include <igl/bounding_box.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Sparse>
#include <fstream>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

	namespace
	{
		template<typename T>
		void removeRow(T& matrix, unsigned int rowToRemove)
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

			//std::cout << tNew << "\n\n" << std::endl;

			// add new rows at the end of triangles
			Eigen::MatrixXi tmp = triangles;
			triangles.resize(triangles.rows() + tNew.rows(), triangles.cols());

			triangles << tmp, tNew;
		}

		template<typename Matrix>
		bool load_data(const std::string & path, Matrix & data_matrix)
		{
			typedef typename Matrix::Scalar Scalar;
			std::fstream file;
			file.open(path);

			if (!file.good())
			{
				std::cerr << "Failed to open file : " << path << std::endl;
				file.close();
				return false;
			}

			std::string s;
			std::vector<std::vector<Scalar>> matrix;

			while (getline(file, s))
			{
				std::stringstream input(s);
				double temp;
				matrix.emplace_back();

				auto &currentLine = matrix.back();

				while (input >> temp)
					currentLine.push_back(temp);
			}

			if (!igl::list_to_matrix(matrix, data_matrix))
			{
				std::cerr << "list to matrix error" << std::endl;
				file.close();
				return false;
			}

			return true;
		}
	}


	bool Mesh::load(const nlohmann::json &data)
	{
		params = {};
		detected.resize(0, 0); // detected (unmoved) point positions
		moved.resize(0, 0); // manually moved
		points.resize(0, 0); // relaxed point positions
		triangles.resize(0, 0); // triangular mesh
		adj.clear(); // adjaceny list of triangluar mesh
		vertex_to_tri.clear();
		boundary.resize(0); // list of vertices on the boundary
		vertex_status_fixed.resize(0);

		read_json_mat(data["detected"], detected);
		read_json_mat(data["points"], points);
		read_json_mat(data["moved"], moved);
		read_json_mat(data["triangles"], triangles);

		params.load(data["params"]);

		solved_vertex.resize(points.rows(), 1);
		solved_vertex.setConstant(false);

		adjacency_list(triangles, adj);
		generate_vertex_to_tri();

		vertex_status_fixed.resize(points.rows(), 1);
		vertex_status_fixed.setZero();

		return true;
	}


	// bool Mesh::load(const std::string & path)
	// {
	// 	if (path.empty()) { return false; }
	// 	// clear previous
	// 	params = {};
	// 	detected.resize(0, 0); // detected (unmoved) point positions
	// 	moved.resize(0, 0); // manually moved
	// 	points.resize(0, 0); // relaxed point positions
	// 	triangles.resize(0, 0); // triangular mesh
	// 	adj.clear(); // adjaceny list of triangluar mesh
	// 	vertex_to_tri.clear();
	// 	boundary.resize(0); // list of vertices on the boundary
	// 	vertex_status_fixed.resize(0);

	// 	// Load data
	// 	load_data(path + "/cellogram/detected.vert", detected);
	// 	load_data(path + "/cellogram/points.vert", points);
	// 	load_data(path + "/cellogram/moved.vert", moved);
	// 	load_data(path + "/cellogram/mesh.tri", triangles);

	// 	// change when moved is also saved
	// 	//moved = detected;

	// 	params.load(path + "/cellogram/params.json");

	// 	solved_vertex.resize(points.rows(), 1);
	// 	solved_vertex.setConstant(false);

	// 	adjacency_list(triangles, adj);
	// 	generate_vertex_to_tri();

	// 	vertex_status_fixed.resize(points.rows(), 1);
	// 	vertex_status_fixed.setZero();

	// 	return true;
	// }

	void Mesh::relax_with_lloyd(const int lloyd_iterations, const Eigen::MatrixXd &hull_vertices,const Eigen::MatrixXi &hull_faces, const bool fix_regular_regions)
	{
		//reset the state
		//points = detected;
		points = moved;
		compute_triangulation();

		Eigen::VectorXi fixed_V;
		if (fix_regular_regions)
		{
			Eigen::VectorXd energy;
			laplace_energy(moved, triangles, energy);
			// Determine whether each vertex passes the criterium for bad mesh
			double avg = energy.mean();

			// Find the degree of each vertex
			Eigen::VectorXi degree;
			vertex_degree(degree);

			Eigen::Matrix<bool, 1, Eigen::Dynamic> low_energy(moved.rows());
			low_energy.setConstant(false);
			int count = 0;
			for (int i = 0; i < moved.rows(); i++)
			{
				if (energy(i) < 0.6*avg && degree(i) == 6)
				{
					low_energy(i) = true;
					count++;
				}
			}
			for (int i = 0; i < boundary.size(); i++)
			{
				low_energy(boundary(i)) = true;
				count++;
			}
			fixed_V.resize(low_energy.size());
			count = 0;
			for (int i = 0; i < low_energy.size(); i++)
			{
				if (low_energy(i))
				{
					fixed_V(count) = i;
					count++;
				}
				//check if the vertices were manually moved
				else if((moved.row(i) - detected.row(i)).squaredNorm() > 1.)
				{
					fixed_V(count) = i;
					count++;
				}
			}
			fixed_V.conservativeResize(count);
		}
		else
			fixed_V = boundary;

		lloyd_relaxation(points, fixed_V, lloyd_iterations, hull_vertices, hull_faces);
		compute_triangulation();
	}

	void Mesh::vertex_degree(Eigen::VectorXi & degree)
	{
		cellogram::vertex_degree(triangles, degree);
	}

	void Mesh::detect_vertices(const Eigen::MatrixXd &V, const DetectionParams &params)
	{
		detected.resize(0, 0); // detected (unmoved) point positions
		points.resize(0, 0); // relaxed point positions
		triangles.resize(0, 0); // triangular mesh
		adj.clear(); // adjaceny list of triangluar mesh
		vertex_to_tri.clear();
		boundary.resize(0); // list of vertices on the boundary

		if (V.size() == 0)
			return;

		this->params = params;
		detected = V;

		moved = detected;
		points = detected;
		vertex_status_fixed.resize(points.rows(), 1);
		vertex_status_fixed.setZero();
		solved_vertex.resize(points.rows(), 1);
		solved_vertex.setConstant(false);


		if (V.rows() < 3)
			return;

		//loose_convex_hull(moved, boundary, 6); moved to compute_triangulation 
		compute_triangulation();

		// automatically load params if available

	}

	void Mesh::delete_vertex(const int index, bool recompute_triangulation)
	{
		// Delete vertex
		removeRow(detected, index);
		removeRow(moved, index);
		removeRow(solved_vertex, index);
		removeRow(vertex_status_fixed, index);

		// Delete entry in params
		if (params.A.size() > 0)
		{
			params.remove_index(index);
		}

		for (int i = 0; i < boundary.size(); ++i)
		{
			if (boundary(i) > index)
				--boundary(i);
		}

		if (recompute_triangulation)
		{
			reset();
		}
		else
		{
			// delete from points
			removeRow(points, index);

			// delete triangles
			std::vector<int> ind;
			ind = vertex_to_tri[index];
			std::sort(ind.begin(), ind.end());

			for (int i = ind.size()-1; i >= 0; i--)
			{
				removeRow(triangles,ind[i]);
			}
			// lower index of values higher than index
			for (int i = 0; i < triangles.rows(); i++)
			{
				for (int j = 0; j < triangles.cols(); j++)
				{
					if (triangles(i, j) > index)
						triangles(i, j) = triangles(i, j) - 1;
				}
			}
			adjacency_list(triangles, adj);
			generate_vertex_to_tri();
		}
	}

	void Mesh::add_vertex(Eigen::Vector3d & new_point, bool must_reset)
	{
		// Add vertex
		{
			Eigen::MatrixXd tmp(detected.rows() + 1, detected.cols());
			tmp.block(0, 0, detected.rows(), detected.cols()) = detected;
			tmp.row(detected.rows()) = new_point.transpose();
			detected = tmp;
		}
		// Add vertex to moved
		{
			Eigen::MatrixXd tmp(moved.rows() + 1, moved.cols());
			tmp.block(0, 0, moved.rows(), moved.cols()) = moved;
			tmp.row(moved.rows()) = new_point.transpose();
			moved = tmp;
		}
		// Add new entry to solved_vertex
		Eigen::Matrix<bool, Eigen::Dynamic,1> tmp_bool;
		tmp_bool.resize(solved_vertex.rows()+1);
		tmp_bool.block(0, 0, solved_vertex.rows(), solved_vertex.cols()) = solved_vertex;
		tmp_bool(solved_vertex.rows()) = false;
		solved_vertex = tmp_bool;

		// Add new entry to vertex_status_fixed
		Eigen::MatrixXi tmp_i;
		tmp_i.resize(vertex_status_fixed.rows() + 1, vertex_status_fixed.cols());
		tmp_i.block(0, 0, vertex_status_fixed.rows(), vertex_status_fixed.cols()) = vertex_status_fixed;
		tmp_i(vertex_status_fixed.rows()) = 0;
		vertex_status_fixed = tmp_i;

		// Add zero row to params
		params.push_back_const(0);

		if(must_reset)
			reset();
	}

	void Mesh::local_update(Eigen::VectorXi &local2global, const int global_to_remove, Eigen::MatrixXi & new_triangles)
	{
		//decide if to remove is local or global (better local)
		//remove vertex to_remove from points detected params (ie fix the existing method) delete_vertex(local2global(local_to_remove), false); and update triangulation and the local2global important
		delete_vertex(global_to_remove, false);
		Eigen::MatrixXi tGlobal = Eigen::MatrixXi(new_triangles.rows(), 3);
		for (int i = 0; i < new_triangles.rows(); i++)
		{
			for (int j = 0; j < 3; j++)
			{
				const int gindex = local2global(new_triangles(i, j));
				assert(gindex != global_to_remove);
				tGlobal(i, j) = gindex -(gindex >= global_to_remove ? 1 : 0);
			}

		}

		Eigen::MatrixXi tmp = triangles;
		triangles.resize(triangles.rows() + tGlobal.rows(), triangles.cols());

		triangles << tmp, tGlobal;

		adjacency_list(triangles, adj);
		generate_vertex_to_tri();
	}

	void Mesh::update_triangles_from_split(const Eigen::VectorXi & t_ind_old, const  Eigen::MatrixXi & t1, const Eigen::MatrixXi & t2)
	{
		// first remove all the old triangles
		for (int i = t_ind_old.size()-1; i >= 0; i--)
		{
			removeRow(triangles, t_ind_old(i));
		}

		int nT = triangles.rows();

		// append the new triangles and return their indices in triangle
		int index = triangles.rows();

		triangles.conservativeResize(triangles.rows() + t1.rows() + t2.rows(), 3);
		if (t1.size() > 0)
		{
			triangles.block(index, 0, t1.rows(), 3) = t1;
			index += t1.rows();
		}
		if (t2.size() > 0)
		{
			triangles.block(index, 0, t2.rows(), 3) = t2;
			index += t2.rows();
		}
		assert(index == triangles.rows());
		adjacency_list(triangles, adj);
		generate_vertex_to_tri();
	}

	void Mesh::mark_vertex_as_solved(const Eigen::VectorXi & region_interior)
	{
		for (int i = 0; i < region_interior.size(); i++)
		{
			solved_vertex(region_interior(i)) = true;
		}
	}

	void Mesh::get_physical_bounding_box(double scaling, Eigen::Vector2d & min, Eigen::Vector2d & max) const
	{
		Eigen::RowVector3d maxVal = points.colwise().maxCoeff();
		Eigen::RowVector3d minVal = points.colwise().minCoeff();

		max = maxVal.block<1,2>(0, 0)*scaling;
		min = minVal.block<1,2>(0, 0)*scaling;

		// std::cout << "max:\n" << max << "\n\nmin:\n" << min << std::endl;
	}

	void Mesh::get_background_mesh(double scaling, Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXd &S, double padding) const {
		Eigen::MatrixXd BV;
		Eigen::MatrixXi BF;
		V = points.leftCols<2>();
		igl::bounding_box(V, padding / scaling, BV, BF);
		assert(BV.rows() == 4);

		V.resize(points.rows() + BV.rows(), 2);
		V.topRows(points.rows()) = points.leftCols<2>();
		V.bottomRows(BV.rows()) = BV;

		delaunay_triangulation(V, F);
		V *= scaling;

		S.resize(V.rows());
		S.head(points.rows()) = (detected - points).rowwise().norm();
		S.tail(BV.rows()).setZero();
	}

	void Mesh::local_update(Eigen::VectorXi & local2global, Eigen::MatrixXd & new_points, Eigen::MatrixXi & new_triangles, Eigen::VectorXi & old_triangles)
	{
		for (int i = 0; i < local2global.size(); i++)
		{
			const int global_index = local2global(i);
			points.row(global_index) = new_points.row(i);

		}

		Eigen::MatrixXi tGlobal = Eigen::MatrixXi(new_triangles.rows(), 3);
		for (int i = 0; i < new_triangles.rows(); i++)
		{
			for (int j = 0; j < 3; j++)
			{
				tGlobal(i, j) = local2global(new_triangles(i, j));
			}

		}

		replaceTriangles(tGlobal, old_triangles, triangles);

		adjacency_list(triangles, adj);
		generate_vertex_to_tri();
	}

	void Mesh::final_relax(const Eigen::VectorXi & expanded_boundary)
	{
		int n = points.rows();
		////Find vertices that have lower connectivity than 6

		Eigen::VectorXi neighCount;
		vertex_degree(neighCount);
		
		// Determine fixed vertices based on connectivity and bounding box
		Eigen::VectorXi indFixed = Eigen::VectorXi::Zero(n);

		for (size_t i = 0; i < n; i++)
		{
			if (neighCount(i) != 6)
			{
				indFixed(i) = 1;
			}
		}
		
		// fix vertices on the boundary and the ones connected to the boundary
		for (size_t i = 0; i < expanded_boundary.rows(); i++)
		{
			indFixed(expanded_boundary(i)) = 1;
		}

		//////PLEASE USE ME
		//igl::opengl::glfw::Viewer viewer;

		//viewer.data().set_mesh(points, triangles);
		//Eigen::MatrixXd cc(n, 3);
		//for (int i = 0; i < indFixed.rows(); ++i) {
		//	if(indFixed(i) == 1)
		//		cc.row(i) = Eigen::RowVector3d(1, 0, 0);
		//	else
		//		cc.row(i) = Eigen::RowVector3d(0, 0, 0);
		//}
		//viewer.data().set_points(points, cc);
		//viewer.data().point_size = float(8);
		//viewer.core.set_rotation_type(igl::opengl::ViewerCore::RotationType::ROTATION_TYPE_NO_ROTATION);
		//viewer.core.orthographic = true;
		//viewer.core.is_animating = true;
		//viewer.launch();


		// Rearrange coordinates and triangle indices to have free vertices first...
		//MatrixXd xRearranged = MatrixXd::Zero(n, 1);
		Eigen::SparseVector<double> xRearranged(n);
		Eigen::SparseVector<double> yRearranged(n);
		Eigen::VectorXi indicesMapping = Eigen::VectorXi::Zero(n);
		Eigen::MatrixXi trianglesRearranged = Eigen::MatrixXi::Zero(triangles.rows(), 3);


		int c = 0;
		for (size_t i = 0; i < n; i++)
		{
			if (indFixed(i) == 0)
			{
				indicesMapping(i) = c;
				xRearranged.fill(c) = points(i, 0);
				yRearranged.fill(c) = points(i, 1);

				c++;
			}
		}
		int nFree = c;

		// and then the fixed vertices
		for (size_t i = 0; i < n; i++)
		{
			if (indFixed(i) == 1)
			{
				indicesMapping(i) = c;
				xRearranged.fill(c) = points(i, 0);
				yRearranged.fill(c) = points(i, 1);

				c++;
			}
		}

		// rearrange triangles such that it is congruent with the newly arranged coordinates
		for (size_t j = 0; j < triangles.rows(); j++)
		{
			for (size_t k = 0; k < 3; k++)
			{
				trianglesRearranged(j, k) = indicesMapping(triangles(j, k));
			}
		}

		// Generate Laplacian
		Eigen::VectorXd diag = Eigen::VectorXd::Zero(n);
		Eigen::SparseMatrix<double> L(n, n);
		typedef Eigen::Triplet<int> Trip;
		std::vector< Trip > tripletList;
		tripletList.reserve(n * 7);

		for (size_t i = 0; i < trianglesRearranged.rows(); i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				tripletList.push_back(Trip(trianglesRearranged(i, j), trianglesRearranged(i, (j + 1) % 3), -1));
				tripletList.push_back(Trip(trianglesRearranged(i, (j + 1) % 3), trianglesRearranged(i, j), -1));
				diag(trianglesRearranged(i, j))++;
			}
		}

		for (size_t i = 0; i < diag.rows(); i++)
		{
			tripletList.push_back(Trip(i, i, diag(i)));
		}

		L.setFromTriplets(tripletList.begin(), tripletList.end());

		// Force all non-zeros to be one
		for (int k = 0; k<L.outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(L, k); it; ++it)
			{
				if (it.value() < 0)
				{
					L.coeffRef(it.row(), it.col()) = -1;
				}
			}
		}

		// Solve for xNew and yNew
		Eigen::SparseMatrix<double>  Li = L.block(0, 0, nFree, nFree);
		Eigen::SparseMatrix<double>  Lb = L.block(0, nFree, nFree, L.rows() - nFree);

		//SparseMatrix<int>  Lii = Li.inverse();
		Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
		solver.compute(Li);
		if (solver.info() != Eigen::Success) {
			// decomposition failed
			return;
		}

		Eigen::SparseVector<double> xB(L.rows() - nFree);
		Eigen::SparseVector<double> yB(L.rows() - nFree);
		for (size_t i = 0; i < L.rows() - nFree; i++)
		{
			xB.fill(i) = xRearranged.coeffRef(nFree + i);
			yB.fill(i) = yRearranged.coeffRef(nFree + i);
		}
		Eigen::SparseVector<double> tmp = -Lb * xB;

		Eigen::VectorXd xNew = solver.solve(Eigen::VectorXd(tmp));

		tmp = -Lb * yB;
		Eigen::VectorXd yNew = solver.solve(Eigen::VectorXd(tmp));

		// Overwrite the free vertices of xRearranged and yRearranged
		for (size_t i = 0; i < nFree; i++)
		{
			xRearranged.coeffRef(i) = xNew(i);
			yRearranged.coeffRef(i) = yNew(i);
		}

		// use indices mapping to overwrite "points" with the newly calculated coordinates
		for (size_t i = 0; i < indicesMapping.rows(); i++)
		{
			points(i,0) = xRearranged.coeffRef(indicesMapping(i));
			points(i,1) = yRearranged.coeffRef(indicesMapping(i));
		}
	}

	void Mesh::generate_vertex_to_tri()
	{
		vertex_to_tri.clear();
		vertex_to_tri.resize(points.rows());

		for (int f = 0; f < triangles.rows(); ++f) {
			for (int lv = 0; lv < triangles.cols(); ++lv)
			{
				vertex_to_tri[triangles(f, lv)].push_back(f);
			}
		}
	}

	void Mesh::untangle()
	{
		Eigen::MatrixXd newPts;
		std::vector<int> dropped;
		// std::cout<<points<<std::endl;
		cellogram::PointsUntangler::pointsUntangler(moved, triangles, dropped, newPts);
		// assert(moved.rows() - dropped.size() + newPts.rows() == triangles.maxCoeff() - 1);

		for (int i = 0; i < newPts.rows(); ++i) {
			Eigen::Vector3d tmp; tmp.setZero();
			for (int j = 0; j < newPts.cols(); ++j)
				tmp(j) = newPts(i, j);
			add_vertex(tmp, false);
		}
		if(newPts.rows()>0)
			std::cout << "New points added (" << newPts.rows() << ")" << std::endl;

		for(int gid : dropped)
			delete_vertex(gid, false);

		points = moved;
		solved_vertex.setConstant(false);

		assert(triangles.maxCoeff() < points.rows());

		adjacency_list(triangles, adj);
		generate_vertex_to_tri();
	}

	void Mesh::clear()
	{
		detected.resize(0, 0); // detected (unmoved) point positions
		moved.resize(0, 0); // relaxed point positions
		points.resize(0, 0); // relaxed point positions
		triangles.resize(0, 0); // triangular mesh
		adj.clear(); // adjaceny list of triangluar mesh
		vertex_to_tri.clear();
		boundary.resize(0); // list of vertices on the boundary

	}


	void Mesh::save(nlohmann::json &data)
	{
		data["detected"] = json::object();
		write_json_mat(detected, data["detected"]);

		data["moved"] = json::object();
		write_json_mat(moved, data["moved"]);

		data["points"] = json::object();
		write_json_mat(points, data["points"]);

		data["displacement"] = json::object();
		write_json_mat((points-detected).eval(), data["displacement"]);

		data["triangles"] = json::object();
		write_json_mat(triangles, data["triangles"]);

		data["boundary"] = json::object();
		write_json_mat(boundary, data["boundary"]);

		data["params"] = json::object();

		params.save(data["params"]);
	}

	// void Mesh::save(const std::string & path)
	// {
	// 	{
	// 		std::ofstream out(path + "/detected.vert");
	// 		out << detected << std::endl;
	// 		out.close();
	// 	}

	// 	{
	// 		std::ofstream out(path + "/moved.vert");
	// 		out << moved << std::endl;
	// 		out.close();
	// 	}

	// 	{
	// 		std::ofstream out(path + "/points.vert");
	// 		out << points << std::endl;
	// 		out.close();
	// 	}

	// 	{
	// 		std::ofstream out(path + "/displacement.txt");
	// 		out << (points-detected) << std::endl;
	// 		out.close();
	// 	}

	// 	{
	// 		std::ofstream out(path + "/mesh.tri");
	// 		out << triangles << std::endl;
	// 		out.close();
	// 	}

	// 	{
	// 		std::ofstream out(path + "/boundary.txt");
	// 		out << boundary << std::endl;
	// 		out.close();
	// 	}

	// 	params.save(path);
	// }

	void Mesh::reset()
	{
		//points = detected;
		points = moved;
		solved_vertex.setConstant(false);
		
		//recompute boundary

		compute_triangulation();
	}

	void Mesh::compute_triangulation()
	{
		if (points.size() == 0)
			return;
		// delaunay_triangulation(points, triangles);

		loose_convex_hull(moved, boundary, 6);
		constrained_delaunay_triangulation(points, boundary, triangles);

		// Calculate the graph adjancency
		adjacency_list(triangles, adj);
		generate_vertex_to_tri();

	}

}// namespace cellogram


