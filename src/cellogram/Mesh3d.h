////////////////////////////////////////////////////////////////////////////////
#pragma once
#include <cellogram/Mesh.h>
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

	class Mesh3d {

	public:
		Mesh3d() { }

		nlohmann::json simulation_out = nlohmann::json({});

		Eigen::MatrixXd V;
		Eigen::MatrixXi F;
		Eigen::MatrixXi T;
		Eigen::MatrixXd displacement;
		Eigen::MatrixXd traction_forces;
		// Eigen::VectorXd sizing;

		bool empty();
		bool analysed();

		void init_nano_dots(const Mesh &mesh, float padding_size, float thickness, float E, float nu, double  scaling, double zscaling, const std::string &formulation, const std::string &save_dir);
		void init_pillars(const Mesh &mesh, float eps, float I, float L, double scaling, double zscaling);
		void clear();

		bool load(const nlohmann::json &data);
		void save_mesh(nlohmann::json &data);
		void save_traction(nlohmann::json &data);

        static void GetMarker3D(const Eigen::MatrixXd &marker4d, Eigen::MatrixXd &marker3d);
	};

} // namespace cellogram
