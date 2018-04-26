#pragma once

#include <Eigen/Dense>

namespace cellogram
{
	namespace internal
	{
		struct Params
		{
			double A;
			double sigma;
			double C;

			double std_x, std_y;
			double std_A;
			double std_sigma;
			double std_C;

			double mean;
			double std;
			double RSS;
		};
	}
	void fitGaussian2D(const Eigen::MatrixXd &window, double x0, double y0, double A0, double sigma0, double C0, Eigen::Vector2d &xy, internal::Params &params);
}