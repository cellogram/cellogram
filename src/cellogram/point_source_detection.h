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

		void sum(const DetectionParams value);
		void divide(const double divisor);

		DetectionParams get_index(const int index);

		int size() const;

		void resize(const int size);
		void setZero(const int size);

		void conservative_resize(const int size);

		void set_from(internal::Params &params, const int index);
		void set_from(DetectionParams &params, const int index);

		void remap(const Eigen::VectorXi &map);

		void remove_index(const int index);

		void push_back_const(const double value);

		void push_back_params(const DetectionParams &value);

		void save(const std::string &path);

		void load(const std::string &path);

		void print();
	};

	void point_source_detection(const Eigen::MatrixXd &img, const double sigma, Eigen::MatrixXd &V, DetectionParams &params);
} // namespace cellogram
