#pragma once

#include "mesh.h"
#include "grid.h"

#include <Eigen/Dense>


namespace cellogram
{
	namespace PointsUntangler
	{
		void pointsUntangler(Mesh &m, Grid &g, const std::string &outputPath = "");
		void pointsUntangler(const Eigen::MatrixXd &detected, Eigen::MatrixXi &tris, Eigen::MatrixXd &newPoints);
	}
}

