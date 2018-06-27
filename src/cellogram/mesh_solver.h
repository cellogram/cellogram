#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <cellogram/common.h>
#include <vector>
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {
	class mesh_solver
	{
	public:
		std::vector<std::vector<int>> region_edges;
		Eigen::VectorXi region;


	public:
		///
		/// @brief      { Finds the regions in the mesh that don't pass the low energy requirement}
		///
		/// @param[in]  V		  { #V x dims input point positions }
		/// @param[in]  F		  { #F x 3 input triangle indices }
		/// @param[out] regions   { bad regions found by indices in V}
		//void find_bad_regions(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);
		///
		/// @brief      { Calls the gurobi solver for each region }
		///
		/// @param[in]  V		  { #V x dims input point positions }
		/// @param[in]  F		  { #F x 3 input triangle indices }
		void solve_regions(Eigen::MatrixXd & pts,  Eigen::MatrixXd &V, Eigen::MatrixXi &F);


		//void compute_regions_edges(const Eigen::MatrixXd &V, Eigen::MatrixXd &bad_P1, Eigen::MatrixXd &bad_P2);
	};
} // namespace cellogram
