#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <zebrafish/autodiff.h>

#include <cassert>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

enum MOUSE_TYPE {
    MOUSEDOWN = 1, 
    MOUSEUP   = 2,
    MOUSEMOVE = 3
};

enum REJECT_MODE {
    REJECT_SINGLE = 0,
    REJECT_AREA   = 1
};

enum COMPRESS_METHOD {
    COMPRESS_MAX = 0,
    COMPRESS_AVG = 1
};

///////////
// types //
///////////
typedef std::vector<Eigen::MatrixXd> image_t;  // a 3D double matrix
typedef std::vector<image_t> imageData_t;  // an array of 3D image (different frames)
typedef Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> texture_t;
typedef std::vector<texture_t> textureArray_t;
typedef Eigen::Vector3d gradient_t;  // gradient
typedef DScalar1<double, gradient_t> DScalar;

}  // namespace zebrafish
