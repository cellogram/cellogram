#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

	bool read_tif_image(const std::string &path, std::vector<Eigen::MatrixXd> & img3D);
	bool read_png_image(const std::string &path, Eigen::MatrixXd &img);

	bool read_image(const std::string &path, std::vector<Eigen::MatrixXd> & img3D);
    bool WriteTif(const std::string &path, const std::vector<Eigen::MatrixXd> &image, int sliceBegin, int sliceEnd);
} // namespace cellogram
