#include <Eigen/Dense>

namespace cellogram
{
	void fitGaussian2D(const Eigen::MatrixXd &window, double x0, double y0, double A0, double sigma0, double C0);
}