#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>

#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

typedef struct OptimPara_t {
    double epsilon;  // LBFGS stops when (gradient_norm) < max(x_norm, 1.0) * epsilon
    int maxIt;  // LBFGS max iteration

    OptimPara_t() : epsilon(1e-4), maxIt(50) {}
} OptimPara_t;

////////////////////////////////////////////////
// optim [only has static member functions]

class optim {

private:


public:
    static void Optim_FixDepth(const OptimPara_t &optimPara, const bspline &bsp, const Eigen::MatrixXd &CI, Eigen::MatrixXd &CO, bool invertColor=false);
    static bool Optim_FixDepth(const OptimPara_t &optimPara, const bspline &bsp, const Eigen::VectorXd &CI, Eigen::VectorXd &CO, double &EO, int &IterO, bool invertColor=false);
    static void Optim_FixDepth(const OptimPara_t &optimPara, const bspline &bsp, const Eigen::MatrixXd &CI, Eigen::MatrixXd &CO, Eigen::VectorXd &EO, Eigen::VectorXi &IterO, bool invertColor=false);
    /// Optimize cylinder(s) "CI" using 3D image "bsp"
    /// The depth (z) coordinate is fixed
    /// [NOTE] Image needs to be normalized to [0, 1]
    ///
    /// @param[in]   bsp         { a B-spline solver with image size registered }
    /// @param[in]   CI          { [#cylinder x 4] matrix of input x, y, z, r }
    /// @param[out]  CO          { [#cylinder x 4] matrix of output x, y, z, r }
    /// @param[out]  EO          { [#cylinder] vector of optimized energy }
    /// @param[out]  IterO       { [#cylinder] vector of optimization iteration }
    /// @param[in]   invertColor { [#cylinder] vector of boolean indicating whether treating the color as inverted }
    ///

    static void Optim_WithDepth(const OptimPara_t &optimPara, const bspline &bsp, const int zNum, const double zGap, const Eigen::MatrixXd &CI, Eigen::MatrixXd &CO, bool invertColor=false);
    /// Optimize cylinder(s) "CI" using 3D image "bsp"
    /// Also find the optimal depth. The search range is [z - zNum*zGap, z + zNum*zGap]
    ///
};

}  // namespace zebrafish
