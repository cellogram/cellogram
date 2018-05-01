#pragma once

////////////////////////////////////////////////////////////////////////////////
#include "common.h"
#include <vector>
#include <geogram/basic/geometry.h>
#include <geogram/mesh/mesh.h>
#include <Eigen/Dense>

#include "fitGaussian2D.h"
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

	class DetectionParams
	{
	public:
		Eigen::VectorXd A;
		Eigen::VectorXd sigma;
		Eigen::VectorXd C;

		Eigen::VectorXd std_x, std_y;
		Eigen::VectorXd std_A;
		Eigen::VectorXd std_sigma;
		Eigen::VectorXd std_C;

		Eigen::VectorXd mean;
		Eigen::VectorXd std;
		Eigen::VectorXd RSS;

		Eigen::VectorXd pval_Ar;

		void resize(const int size);
		void conservative_resize(const int size);

		void set_from(internal::Params &params, const int index);

		void remap(const Eigen::VectorXi &map);

		void remove_index(const int index);

		void push_back(const double value);
	};

	void point_source_detection(const Eigen::MatrixXd &img, const double sigma, Eigen::MatrixXd &V, DetectionParams &params);
} // namespace cellogram
