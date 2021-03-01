#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>

#include <Eigen/Dense>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

typedef struct OptimPara_t {
    double epsilon;  // LBFGS stops when (gradient_norm) < max(x_norm, 1.0) * epsilon
    int maxIt;  // LBFGS max iteration
    double zSearchMaxXYDisp;  // max displacement in XY during depth search
    // customized linear search method TODO

    OptimPara_t() : epsilon(1e-4), maxIt(50), zSearchMaxXYDisp(3.0) {}
} OptimPara_t;

typedef struct OptimDepthInfo_t {
    Eigen::MatrixXd C;  // x, y, z, r
    Eigen::VectorXd energy;
    Eigen::VectorXi iter;

    Eigen::MatrixXd ToMat() {  // debug
        Eigen::MatrixXd mat(energy.size(), 6);
        mat << C, energy, iter.cast<double>();
        return mat;
    }
} OptimDepthInfo_t;

typedef enum DepthSearchFlag_t {
    SecondDerivative = -2,
    InvalidEnergy = -1, 
    Unknown = 0,
    Success = 1
} DepthSearchFlag_t;

////////////////////////////////////////////////
// optim [only has static member functions]

class optim {

private:


public:
    static void Optim_FixDepth(const OptimPara_t &optimPara, const bspline &bsp, const Eigen::MatrixXd &CI, Eigen::MatrixXd &CO, const bool invertColor=false);
    static bool Optim_FixDepth(const OptimPara_t &optimPara, const bspline &bsp, const Eigen::VectorXd &CI, Eigen::VectorXd &CO, double &EO, int &IterO, const bool invertColor=false);
    static void Optim_FixDepth(const OptimPara_t &optimPara, const bspline &bsp, const Eigen::MatrixXd &CI, Eigen::MatrixXd &CO, Eigen::VectorXd &EO, Eigen::VectorXi &IterO, const bool invertColor=false);
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

    static void Optim_WithDepth(const OptimPara_t &optimPara, const bspline &bsp, const int zNum, const double zGap, const Eigen::VectorXd &CI, OptimDepthInfo_t &C_depth_info, const bool invertColor=false);
    static void Optim_WithDepth(const OptimPara_t &optimPara, const bspline &bsp, const int zNum, const double zGap, const Eigen::MatrixXd &CI, std::vector<OptimDepthInfo_t> &C_depth_info, const bool invertColor=false);
    /// Optimize cylinder(s) "CI" using 3D image "bsp"
    /// Also find the optimal depth. The search range is [z - zNum*zGap, z + zNum*zGap]
    ///
    /// @param[in]   bsp         { a B-spline solver with image size registered }
    /// @param[in]   zNum        { number of depths to search in each direction (above and below) }
    /// @param[in]   zGap        { search gap }
    /// @param[in]   CI          { [#cylinder x 4] matrix of input x, y, z, r }
    /// @param[out]  CO          { vector{[(2*zNum+1) x 4]} of output x, y, z, r for each row in CI}
    /// @param[in]   invertColor { [#cylinder] vector of boolean indicating whether treating the color as inverted }
    ///

    static void DepthSelection(const OptimPara_t &optimPara, const Eigen::VectorXd &CI, const OptimDepthInfo_t &C_depth_info, Eigen::VectorXd &CO, DepthSearchFlag_t &flag);
    static void DepthSelection(const OptimPara_t &optimPara, const Eigen::MatrixXd &CI, const std::vector<OptimDepthInfo_t> &C_depth_info, Eigen::MatrixXd &CO, std::vector<DepthSearchFlag_t> &flag);
    /// Find the optimal depth
    ///
    /// @param[in]   CI          { [#cylinder x 4] matrix of input x, y, z, r }
    /// @param[in]   C_depth_info{ #cylinder vector of depth search results returned from "Optim_WithDepth" }
    /// @param[out]  CO          { [#cylinder x 4] matrix of resultant x, y, z, r }
    /// @param[out]  flag        { #cylinder vector of success flags }
    ///
};

}  // namespace zebrafish
