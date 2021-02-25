#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <zebrafish/Common.h>

#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

///////////////////////////////////////
// Quadrature

typedef struct quadrature {

    int numDiskPts;    // size of sample points
    int heightLayers;  // number of layer along z-axis
    // hardcoded disk quadrature weights and locations
    Eigen::Matrix<double, Eigen::Dynamic, 2> xyArray;
    Eigen::Matrix<double, Eigen::Dynamic, 1> weightArray;
    // hardcoded 1D gaussian quadrature weights and locations
    Eigen::Matrix<double, Eigen::Dynamic, 1> heightLocArray;
    Eigen::Matrix<double, Eigen::Dynamic, 1> heightWeightArray;

    void LoadDiskParas(int diskQuadMethod);
    /// Select a disk quadrature method and update the inner storage xy & weight array
    /// This can be changed on-the-fly

    // maintenance methods
    quadrature(int diskQuadMethod = 6);
    ~quadrature();
} quadrature;

}  // namespace zebrafish
