////////////////////////////////////////////////////////////////////////////////
#include "common.h"
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

///
/// @brief      { Computes the connectivity graph from the trimesh}
///
/// @param[in]  F     { #F x 3 input triangle indices }
/// @param[out] Graph     { #m x dims connectivity }
///
void tri2hex(const Eigen::MatrixXi &F, std::vector<std::vector<int>> &Graph);

} // namespace cellogram
