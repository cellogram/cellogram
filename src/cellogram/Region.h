#pragma once

////////////////////////////////////////////////////////////////////////////////
#include "common.h"
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {
	class Region
	{
	public:
		static const int NOT_PROPERLY_CLOSED = -2;
		static const int NO_SOLUTION = -1;
		static const int OK = 0;

		bool verbose = true;
		Eigen::VectorXi region_boundary;
		Eigen::VectorXi region_interior;
		Eigen::VectorXi region_triangles;

		void compute_edges(const Eigen::MatrixXd &V, Eigen::MatrixXd &bad_P1, Eigen::MatrixXd &bad_P2);

		void grow();

		void find_points(const Eigen::VectorXi &region_ids, const int id);
		void find_triangles(const Eigen::MatrixXi &F, const Eigen::VectorXi &region_ids, const int id);
		void find_triangles(const Eigen::MatrixXi & F);
		void bounding(const Eigen::MatrixXi &F, const Eigen::MatrixXd &V);

		int resolve(const Eigen::MatrixXd &V_detected, const Eigen::MatrixXd &V_current, const Eigen::MatrixXi &F, const std::vector<std::vector<int>> &adj, const int perm_possibilities, Eigen::MatrixXd  &new_points, Eigen::MatrixXi &new_triangles);

	private:
		Eigen::MatrixXi get_triangulation(const Eigen::MatrixXi &F);
	};
} // namespace cellogram
