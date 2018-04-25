////////////////////////////////////////////////////////////////////////////////
#include "common.h"
#include <vector>
#include <geogram/basic/geometry.h>
#include <geogram/mesh/mesh.h>
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {


void point_source_detection(const Eigen::MatrixXd &img, const double sigma, Eigen::MatrixXd &V, Eigen::MatrixXd &V_std, Eigen::VectorXd &pval_Ar);
} // namespace cellogram
