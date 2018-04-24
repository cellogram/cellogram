#pragma once

#include <Eigen/Dense>

namespace cellogram
{
	namespace PointsUntangler
	{
		void pointsUntangler(const Eigen::MatrixXd &detected, Eigen::MatrixXi &tris, Eigen::MatrixXd &newPoints);
	}
}

