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
void adjacency_list(const Eigen::MatrixXi &F, std::vector<std::vector<int>> &Graph);

void triangle_region_list(const Eigen::VectorXi &vertex_region_id, const Eigen::MatrixXi &F, Eigen::VectorXd &face_region_id);

} // namespace cellogram
