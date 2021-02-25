#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>

#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

typedef struct cylinder_t {
// uniquely defines a cylinder
    double x, y, z, r, h;  // bottom{ x(row), y(col), z } + radius + height
} cylinder_t;


typedef struct samplePoints_t {
// The sample points coordinates [x, y, z] and their weights for one cylinder
    Eigen::Matrix<DScalar, Eigen::Dynamic, 2> innerPoints;  // [#points x 2] point positions
    Eigen::Matrix<DScalar, Eigen::Dynamic, 2> outerPoints;
    Eigen::Matrix<DScalar, Eigen::Dynamic, 1> weights;  // [#points] (Gaussian quadrature weight)*(energy function constants)*(subtraction function)
    Eigen::Matrix<double, Eigen::Dynamic, 1> zArray;   // [#depth] depth along z-axis
} samplePoints_t;

////////////////////////////////////////////////
// cylinder [only has static member functions]

class cylinder {

private:
    template<typename T>
    static void EnergyHelper(const bspline &bsp, const Eigen::Matrix<double, Eigen::Dynamic, 1> &zArray,
                      const T &r, const T &x, const T &y, T &resT);

public:
    static double alpha;  // in [0, 1], how severe we should penalize the peripheral area
    static double K;  // extended radius / cylinder radius
    static double H;  // cylinder height (fixed, do not optimize)

    template<typename T>
    static bool IsValid(const bspline &bsp, const T &x, const T &y, const double z, const T &r, const double h);
    /// Check if the cylinder [x, y, z, r, h] is a valid cylinder.
    /// A valid cylinder must reside inside the image cube and has
    /// a reasonable radius.
    ///
    /// @param[in]   bsp       { a B-spline solver with image size registered }
    /// @param[in]   x, y, z   { coordinate of the cylinder bottom center }
    /// @param[in]   r         { cylinder radius. }
    /// @param[in]   h         { cylinder height. }
    /// @return      { whether the cylidner is valid }

    template <typename T>
    static void EvaluateCylinder(const bspline &bsp, T x, T y, double z, T r, double h, T &res, bool reverse = false);
    /// Calculate sample points for the given cylinder and evaluate the energy.
    /// This function does not check whether the cylinder is valid.
    ///
    /// @param[in]   bsp       { a B-spline solver with control points calculated }
    /// @param[in]   x, y, z   { coordinate of the cylinder bottom center }
    /// @param[in]   r         { cylinder radius. Must be positive }
    /// @param[in]   h         { cylinder height. Must be positive }
    /// @param[out]  res       { evaluated energy }
    ///

    static void SubtractionHelper(const Eigen::MatrixXd &points, const Eigen::VectorXd &weight, Eigen::VectorXd &resWeight);
    /// [ Deprecated ]
    /// Multiply the weight array with a subtraction array such that the integral result
    /// measures the subtraction of two areas.
    /// The subtraction array can be hard or soft (sigmoid).
    ///
    /// @param[in]   points     { [#points x 2] sample locations used to calculate disk quadrature }
    /// @param[in]   weight     { [#points] weights corresponding to the sample points }
    /// @param[out]  resWeight  { [#points] returned weight array }
    ///

    // maintenance methods
    cylinder();
    ~cylinder();
};

}  // namespace zebrafish
