#pragma once

////////////////////////////////////////////////////////////////////////////////
#include "common.h"
#include <geogram/mesh/mesh.h>
#include <Eigen/Dense>
#include <vector>
#include <memory>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

///
/// @brief      { Converts a triangle mesh to a Geogram mesh }
///
/// @param[in]  V     { #V x 3 input mesh vertices }
/// @param[in]  F     { #F x 3 input mesh surface }
/// @param[out] M     { Output Geogram mesh }
///
void to_geogram_mesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, GEO::Mesh &M);

///
/// @brief      { Extract simplices from a Geogram mesh }
///
/// @param[in]  M     { Input Geogram mesh }
/// @param[out] V     { #V x 3 output mesh vertices }
/// @param[out] F     { #F x 3 output mesh faces }
/// @param[out] T     { #T x 4 output mesh tets }
///
void from_geogram_mesh(const GEO::Mesh &M, Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXi &T);

} // namespace poly_fem
