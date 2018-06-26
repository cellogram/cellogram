#pragma once

////////////////////////////////////////////////////////////////////////////////
#include "common.h"
#include <vector>
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

///
/// @brief      { Adds or removes vertices from dataset}
///
/// @param[in]  V     { #V x dims input point positions }
/// @param[out] V     { #V x dims output point positions }
///
void add_vertex(const Eigen::MatrixXd &V);
void delete_vertex(const Eigen::MatrixXd &V);

} // namespace cellogram
