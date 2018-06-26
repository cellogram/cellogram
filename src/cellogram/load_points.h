#pragma once

////////////////////////////////////////////////////////////////////////////////
#include "common.h"
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

///
/// @brief      Loads points from a files.
///
/// @param[in]  filename  { Input filename }
/// @param[out] V         { #V x 3 matrix of point positions }
///
void load_points(const std::string &filename, Eigen::MatrixXd &V);

} // namespace cellogram
