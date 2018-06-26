#pragma once

////////////////////////////////////////////////////////////////////////////////
#include "common.h"
#include <vector>
#include <geogram/basic/geometry.h>
#include <geogram/mesh/mesh.h>
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

void lloyd_relaxation(std::vector<GEO::vec2> &points, const std::vector<bool> &fixed, int num_iter,
	const GEO::Mesh *domain = nullptr);

///
/// @brief         { Perform Lloyd relaxation on the vertices }
///
/// @param[in,out] P         { #P x dims vertices positions }
/// @param[in]     fixed     { List of indices of fixed vertices }
/// @param[in]     num_iter  { Number of iterations }
/// @param[in]     V         { #V x dims vertex position of the domain }
/// @param[in]     F         { #F x 3 list triangle indices of the domain }
///
void lloyd_relaxation(Eigen::MatrixXd &P, const Eigen::VectorXi &fixed, int num_iter,
	const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);

} // namespace cellogram
