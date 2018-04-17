////////////////////////////////////////////////////////////////////////////////
#include "Mesh.h"
#include <cellogram/load_points.h>
#include <cellogram/delaunay.h>
#include <cellogram/tri2hex.h>
#include <cellogram/voronoi.h>
#include <cellogram/vertex_degree.h>

#include <igl/list_to_matrix.h>

#include <fstream>
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

		void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
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


	bool Mesh::load(const std::string & path)
	{
		if (path.empty()) { return false; }
		// clear previous
		detected.resize(0, 0); // detected (unmoved) point positions
		points.resize(0, 0); // relaxed point positions
		triangles.resize(0, 0); // triangular mesh
		adj.clear(); // adjaceny list of triangluar mesh
		boundary.resize(0); // list of vertices on the boundary

		// Load points
		load_points(path, detected);
		points = detected;
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
			std::vector<int> ind;
			for (int i = 0; i < triangles.rows(); i++)
			{
				if (triangles(i, 0) == index || triangles(i, 1) == index || triangles(i, 2) == index)
					ind.push_back(i);
			}
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
		}
	}

	void Mesh::add_vertex(Eigen::Vector3d & new_point)
	{
		// Add vertex
		Eigen::MatrixXd tmp(detected.rows() + 1, detected.cols());
		tmp.block(0, 0, detected.rows(), detected.cols()) = detected;
		tmp.row(detected.rows()) = new_point.transpose();
		detected = tmp;

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
	}

}// namespace cellogram


