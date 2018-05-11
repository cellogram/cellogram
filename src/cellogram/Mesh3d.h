////////////////////////////////////////////////////////////////////////////////
#pragma once
#include "Mesh.h"
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

	class Mesh3d {

	public:
		Mesh3d() { }

		Eigen::MatrixXd V;
		Eigen::MatrixXi F;
		Eigen::MatrixXd sol;

		void init(const Mesh &mesh, float padding_size, float thickness);
		void clear();
	};

} // namespace cellogram
