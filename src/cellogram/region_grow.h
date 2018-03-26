////////////////////////////////////////////////////////////////////////////////
#include "common.h"
#include <Eigen/Dense>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

///
/// @brief      { Takes in graph connectivity matrix, a value vector and a criterium}
///
/// @param[in]  Graph			{ #m x dims input graph connectivity }
/// @param[in]  Vertex_Value	{ #1 x dims value per vertex }
/// @param[in]  criterium		{ #single criterium value for growing region where vertices pass }
/// @param[out] region			{ #connected vertices that pass criterium }
///
void region_grow(const Eigen::MatrixXi &Graph, const Eigen::VectorXi &Vertex_Value, const double criterium, Eigen::VectorXi &region);

// this function may not be called from external
void check_crit(const Eigen::MatrixXi &Graph, const Eigen::VectorXi &Vertex_Value, Eigen::VectorXi &region, Eigen::VectorXi &visited, const double criterium, int group, int ind);
} // namespace cellogram
