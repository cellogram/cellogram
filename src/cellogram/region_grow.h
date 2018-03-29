////////////////////////////////////////////////////////////////////////////////
#include "common.h"
#include <Eigen/Dense>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

///
/// @brief      { Takes in graph connectivity matrix, a value vector and a criterium and checks which points that fulfill the criterium are connected}
///
/// @param[in]  Graph			{ #m x dims input graph connectivity }
/// @param[in]  Vertex_Value	{ #1 x dims value per vertex }
/// @param[in]  criterium		{ #single criterium value for growing region where vertices pass }
/// @param[out] region			{ #connected vertices that pass criterium }
///
void region_grow(std::vector<std::vector<int>> &Graph, const std::vector<bool> &crit_pass, Eigen::VectorXi &region);

///
/// @brief      { Takes in the regions from "region_grow" and calculates their edge}
///
/// @param[in]  Graph			{ #m x dims input graph connectivity }
/// @param[in]  Vertex_Value	{ #1 x dims value per vertex }
/// @param[in]  criterium		{ #single criterium value for growing region where vertices pass }
/// @param[out] region			{ #connected vertices that pass criterium }
///
void region_bounding(const Eigen::MatrixXd &points, const Eigen::MatrixXi &triangles, const Eigen::VectorXi &region, std::vector<std::vector<int>> &region_edges);

// this function may not be called from external
//void check_crit(const Eigen::MatrixXi &Graph, const Eigen::VectorXi &Vertex_Value, Eigen::VectorXi &region, Eigen::VectorXi &visited, const double criterium, int group, int ind);
} // namespace cellogram
