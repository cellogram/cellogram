////////////////////////////////////////////////////////////////////////////////
#pragma once
#include "Mesh.h"

#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

	class Mesh3d {

	public:
		Mesh3d() { }

		Eigen::MatrixXd V;
		Eigen::MatrixXi F;
		Eigen::MatrixXd sol;
		Eigen::MatrixXd traction_forces;

		void init(const Mesh &mesh, float padding_size, float thickness, float lambda, float mu, const std::string &formulation);
		void clear();
	};

} // namespace cellogram
