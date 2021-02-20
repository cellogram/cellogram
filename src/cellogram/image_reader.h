#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

	bool read_tif_image(const std::string &path, std::vector<Eigen::MatrixXd> & img3D);
	bool read_png_image(const std::string &path, Eigen::MatrixXd &img);

	bool read_image(const std::string &path, std::vector<Eigen::MatrixXd> & img3D);
} // namespace cellogram
