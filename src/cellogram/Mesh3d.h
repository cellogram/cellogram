////////////////////////////////////////////////////////////////////////////////
#pragma once
#include <cellogram/Mesh.h>
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

	class Mesh3d {

	public:
		Mesh3d() { }

		Eigen::MatrixXd V;
		Eigen::MatrixXi F;
		Eigen::MatrixXi T;
		Eigen::MatrixXd displacement;
		Eigen::MatrixXd traction_forces;
		Eigen::VectorXd sizing;

		bool empty();
		bool analysed();

		void init_nano_dots(const Mesh &mesh, float padding_size, float thickness, float E, float nu, double  scaling, const std::string &formulation);
		void init_pillars(const Mesh &mesh, float eps, float I, float L, double scaling);
		void clear();

		bool load(const nlohmann::json &data);
		void save_mesh(nlohmann::json &data);
		void save_traction(nlohmann::json &data);
	};

} // namespace cellogram
