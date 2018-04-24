////////////////////////////////////////////////////////////////////////////////
#include "Mesh.h"
#include <cellogram/load_points.h>
#include <cellogram/delaunay.h>
#include <cellogram/tri2hex.h>
#include <cellogram/voronoi.h>
#include <cellogram/vertex_degree.h>

#include <igl/list_to_matrix.h>

#include <fstream>
#include <Eigen/Sparse>

#include <igl/opengl/glfw/Viewer.h>

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

			std::cout << tNew << "\n\n" << std::endl;

			// add new rows at the end of triangles
			Eigen::MatrixXi tmp = triangles;
			triangles.resize(triangles.rows() + tNew.rows(), triangles.cols());

			triangles << tmp, tNew;
		}

		
	}


	bool Mesh::load(const std::string & path)
	{
		if (path.empty()) { return false; }
		// clear previous
		detected.resize(0, 0); // detected (unmoved) point positions
		points.resize(0, 0); // relaxed point positions
		triangles.resize(0, 0); // triangular mesh
		adj.clear(); // adjaceny list of triangluar mesh
		vertex_to_tri.clear();
		boundary.resize(0); // list of vertices on the boundary

		// Load points
		load_points(path, detected);
		points = detected;
		solved_vertex.resize(points.rows(), 1);
		solved_vertex.setConstant(false);

		compute_triangulation();

		// automatically load params if available

		return true;
	}

	bool Mesh::load_params(const std::string & path)
	{
		params.resize(0, 0);
		std::fstream file;
		file.open(path.c_str());

		if (!file.good())
		{
			std::cerr << "Failed to open file : " << path << std::endl;
			file.close();
			return false;
		}


		std::string s;
		std::vector<std::vector<double>> matrix;

		while (getline(file, s))
		{
			std::stringstream input(s);
			double temp;
			matrix.emplace_back();

			std::vector<double> &currentLine = matrix.back();

			while (input >> temp)
				currentLine.push_back(temp);
		}

		if (!igl::list_to_matrix(matrix, params))
		{
			std::cerr << "list to matrix error" << std::endl;
			file.close();
			return false;
		}
		assert(detected.rows() == params.rows());
		return true;
	}

	void Mesh::relax_with_lloyd(const int lloyd_iterations, const Eigen::MatrixXd &hull_vertices,const Eigen::MatrixXi &hull_faces)
	{
		//reset the state
		points = detected;
		compute_triangulation();

		lloyd_relaxation(points, boundary, lloyd_iterations, hull_vertices, hull_faces);
		compute_triangulation();
	}

	void Mesh::vertex_degree(Eigen::VectorXi & degree)
	{
		cellogram::vertex_degree(triangles, degree);
	}

	void Mesh::delete_vertex(const int index, bool recompute_triangulation)
	{
		// Delete vertex
		removeRow(detected, index);
		removeRow(solved_vertex, index);

		// Delete entry in params
		if (params.rows() > 0)
		{
			removeRow(params, index);
		}

		for (int i = 0; i < boundary.size(); ++i)
		{
			if (boundary(i) > index)
				--boundary(i);
		}

		if (recompute_triangulation)
		{
			points = detected;
			compute_triangulation();
		}
		else
		{
			// delete from points
			removeRow(points, index);

			// delete triangles
			std::vector<int> ind, ind2;
			ind = vertex_to_tri[index];

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

	void Mesh::add_vertex(Eigen::Vector3d & new_point)
	{
		// Add vertex
		Eigen::MatrixXd tmp(detected.rows() + 1, detected.cols());
		tmp.block(0, 0, detected.rows(), detected.cols()) = detected;
		tmp.row(detected.rows()) = new_point.transpose();
		detected = tmp;

		// Add new entry to solved_vertex
		Eigen::Matrix<bool, Eigen::Dynamic,1> tmp_bool;
		tmp_bool.resize(solved_vertex.rows()+1);
		tmp_bool.block(0, 0, solved_vertex.rows(), solved_vertex.cols()) = solved_vertex;
		tmp_bool(solved_vertex.rows()) = false;
		solved_vertex = tmp_bool;

		// Add zero row to params
		if (params.rows() > 0)
		{
			Eigen::MatrixXd tmp(params.rows() + 1, params.cols());
			tmp.block(0, 0, params.rows(), params.cols()) = params;
			tmp.row(params.rows()) = Eigen::RowVectorXd::Zero(params.cols());
			params = tmp;
		}
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

	void Mesh::mark_vertex_as_solved(const Eigen::VectorXi & region_interior)
	{
		for (int i = 0; i < region_interior.size(); i++)
		{
			solved_vertex(region_interior(i)) = true;
		}
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

	void Mesh::final_relax()
	{
		int n = points.rows();
		////Find vertices that have lower connectivity than 6
		//Eigen::VectorXi neighCount = Eigen::VectorXi::Zero(n);

		//for (size_t i = 0; i < triangles.rows(); i++)
		//{
		//	for (size_t j = 0; j < 3; j++)
		//	{
		//		neighCount(triangles(i, j))++;
		//	}

		//}

		Eigen::VectorXi neighCount;
		vertex_degree(neighCount);

		//std::cout << neighCount.transpose() << std::endl;
		//std::cout << degree.transpose() << std::endl;

		// Determine fixed vertices based on connectivity and bounding box
		Eigen::VectorXi indFixed = Eigen::VectorXi::Zero(n);

		for (size_t i = 0; i < n; i++)
		{
			if (neighCount(i) != 6)
			{
				indFixed(i) = 1;
			}
		}
		for (size_t i = 0; i < boundary.rows(); i++)
		{
			indFixed(boundary(i)) = 1;
		}

		////PLEASE USE ME 
		//igl::opengl::glfw::Viewer viewer;

		//viewer.data().set_mesh(points, triangles);
		//Eigen::MatrixXd C(n, 3);
		//for (int i = 0; i < indFixed.rows(); ++i) {
		//	if(indFixed(i) == 1)
		//		C.row(i) = Eigen::RowVector3d(1, 0, 0);
		//	else
		//		C.row(i) = Eigen::RowVector3d(0, 0, 0);
		//}
		//viewer.data().set_points(points, C);
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

	void Mesh::reset()
	{
		points = detected;
		compute_triangulation();
	}

	void Mesh::compute_triangulation()
	{
		delaunay_triangulation(points, triangles);

		// Calculate the graph adjancency
		adjacency_list(triangles, adj);
		generate_vertex_to_tri();

	}

}// namespace cellogram

