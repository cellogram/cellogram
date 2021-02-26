#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <zebrafish/Common.h>

#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

double QuantileImage(const zebrafish::image_t &image, double q, int layerBegin = -1, int layerEnd = -1);
/// Calculate quantile(image(:), q) where 0 < q < 1
/// Note: for performance consideration, q should be larger than 0.95

void NormalizeImage(image_t &image, double thres);
/// Normalize image to 0-1 and trim values larger than "thres"

}  // namespace zebrafish
