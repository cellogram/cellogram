////////////////////////////////////////////////////////////////////////////////
#include "Mesh3d.h"
#include <cellogram/Mesh.h>
#include <cellogram/State.h>
#include <cellogram/remesh_adaptive.h>
#include <polyfem/State.hpp>
#include <polyfem/Mesh3D.hpp>
#include <polyfem/MeshUtils.hpp>
#include <igl/writeOBJ.h>
#include <polyfem/PointBasedProblem.hpp>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <geogram/mesh/mesh_io.h>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

	namespace
	{
		nlohmann::json compute_analysis(const Eigen::MatrixXd &vertices, const Eigen::MatrixXi &faces, const Eigen::MatrixXi &tets, const Mesh &mesh,
			float thickness, float E, float nu, const std::string &formulation, double scaling,
			Eigen::MatrixXd &vals, Eigen::MatrixXd &traction_forces)
		{
			//TODO
			//const std::string rbf_function = "gaussian";
			const std::string rbf_function = "thin-plate";
			// const std::string rbf_function = "cubic";
			const double eps = 4;

			static const bool export_data = false;
			assert(tets.cols() == 4);
			assert(vertices.cols() == 3);

			GEO::Mesh M;
			M.vertices.create_vertices((int) vertices.rows());
			for (int i = 0; i < (int) M.vertices.nb(); ++i) {
				GEO::vec3 &p = M.vertices.point(i);
				p[0] = vertices(i, 0);
				p[1] = vertices(i, 1);
				p[2] = vertices(i, 2);
			}

			M.cells.create_tets((int) tets.rows());

			for (int c = 0; c < (int) M.cells.nb(); ++c) {
				for (int lv = 0; lv < tets.cols(); ++lv) {
					M.cells.set_vertex(c, lv, tets(c, lv));
				}
			}
			M.cells.connect();

			if(export_data){
				GEO::mesh_save(M, "mesh.mesh");
				GEO::mesh_load("mesh.mesh", M);
			}


			json j_args = {
				{"problem", "PointBasedTensor"},
				{"normalize_mesh", false},

				{"tensor_formulation", formulation},

				{"discr_order", 1},
				{"vismesh_rel_area", 1000},

				{"solver_params", {"conv_tol", 1e-7, "max_iter", 2000}},



				{"nl_solver_rhs_steps", 4},
				{"n_boundary_samples", 3},

				{"params", {
					{"E", E},
					{"nu", nu},
				}},
			};

			polyfem::State state;
			std::string log_file = "";
			bool is_quiet = false;
			int log_level = 1;

			state.init_logger(log_file, log_level, is_quiet);
			state.init(j_args);

			state.load_mesh(M, [&](const polyfem::RowVectorNd &bary){
				// top, Id = 1
				if(std::abs(bary(2)) < 1e-6){
					return 1;
				}

				//Bottom, Id = 3
				if(std::abs(bary(2)- -thickness) < 1e-8){
					return 3;
				}

				//any other
				return 2;
			});

			// state.compute_mesh_stats();

			polyfem::PointBasedTensorProblem &problem = *dynamic_cast<polyfem::PointBasedTensorProblem *>(state.problem.get());

			Eigen::MatrixXd disp = (mesh.detected - mesh.points) * scaling;  // TODO: make this 3D
			Eigen::MatrixXd pts = mesh.points * scaling;

			if(export_data){
				GEO::mesh_save(M, "mesh.mesh");
				{
					std::ofstream out("problem.json");
					out.precision(100);
					out << j_args.dump(4) << std::endl;
					out.close();
				}
				{
					std::ofstream out("fun.pts");
					out.precision(100);
					out << pts << std::endl;
					out.close();
				}

				{
					std::ofstream out("fun.tri");
					out.precision(100);
					out << mesh.triangles << std::endl;
					out.close();
				}

				{
					std::ofstream out("fun.txt");
					out.precision(100);
					out << disp << std::endl;
					out.close();
				}

				{
					std::ofstream out("tags.txt");

					for(int i = 0; i < state.mesh->n_faces(); ++i)
						out << state.mesh->get_boundary_id(i) << std::endl;
					out.close();
				}
			}


			//Id = 1, func, mesh, coord =2, means skip z for the interpolation
			Eigen::Matrix<bool, 3, 1> dirichet_dims;
			dirichet_dims(0) = dirichet_dims(1) = true;
			dirichet_dims(2) = false; // z is not dirichet
			problem.add_function(1, disp, pts, rbf_function, eps, 2, dirichet_dims);

			//Id = 3, zero Dirichelt
			problem.add_constant(3, Eigen::Vector3d(0,0,0));

			// state.compute_mesh_stats();

			state.build_basis();
			// state.build_polygonal_basis();


			state.assemble_rhs();
			state.assemble_stiffness_mat();

			state.solve_problem();

			//true = compute average insteat of just integral

			// state.interpolate_boundary_function(vertices, faces, state.sol, true, vals);
			// state.interpolate_boundary_tensor_function(vertices, faces, state.sol, true, traction_forces);

			state.interpolate_boundary_function_at_vertices(vertices, faces, state.sol, vals);
			state.interpolate_boundary_tensor_function(vertices, faces, state.sol, vals, true, traction_forces);


			// vals = Eigen::Map<Eigen::MatrixXd>(state.sol.data(), 3, vertices.rows());
			// vals = Eigen::Map<Eigen::MatrixXd>(state.rhs.data(), 3, vertices.rows());
			// vals = vals.transpose().eval();
			// std::cout<<vals<<std::endl;


			if(export_data)
			{
				std::ofstream out("sol.txt");
				out.precision(100);
				out << state.sol << std::endl;
				out.close();
			}

			if(export_data)
			{
				Eigen::MatrixXd nodes(state.n_bases, state.mesh->dimension());
				for(const auto &eb : state.bases)
				{
					for(const auto &b : eb.bases)
					{
						for(const auto &lg : b.global())
						{
							nodes.row(lg.index) = lg.node;
						}
					}
				}
				std::ofstream out("nodes.txt");
				out.precision(100);
				out << nodes;
				out.close();
			}

			nlohmann::json json;
			state.save_json(json);
			return json;


			// auto &tmp_mesh = *dynamic_cast<polyfem::Mesh3D *>(state.mesh.get());
			// igl::opengl::glfw::Viewer viewer;
			// Eigen::MatrixXi asdT;
			// Eigen::MatrixXd asdP, asdC;
			// std::vector<int> asdR;
			// tmp_mesh.triangulate_faces(asdT,asdP,asdR);

			// asdC.resize(asdT.rows(), 3);
			// asdC.setZero();

			// for(std::size_t i = 0; i < tmp_mesh.n_faces(); ++i)
			// {
			// 	const int a =tmp_mesh.get_boundary_id(i);
			// 	const auto bb = tmp_mesh.face_barycenter(i);

			// 	if(a == 1)
			// 		viewer.data().add_points(bb, Eigen::RowVector3d(1, 0, 0));
			// 	else if(a == 3)
			// 		viewer.data().add_points(bb, Eigen::RowVector3d(0, 1, 0));
			// 	else if(a == 2)
			// 		viewer.data().add_points(bb, Eigen::RowVector3d(0, 0, 1));
			// }

			// viewer.data().set_mesh(asdP, asdT);
			// viewer.launch();
		}
	}

	void Mesh3d::init_pillars(const Mesh &mesh, float eps, float I, float L, double scaling)
	{
		clear();

		displacement = (mesh.detected - mesh.points) * scaling;
		V = mesh.points * scaling;

		const float bending_force = 3*eps * I /(L*L*L);
		traction_forces = bending_force * displacement;
	}

	bool Mesh3d::empty()
	{
		return V.size() == 0;
	}

	bool Mesh3d::analysed()
	{
		return traction_forces.size() > 0;
	}

	void Mesh3d::init_nano_dots(const Mesh &mesh, float padding_size, const float thickness, float E, float nu, double scaling, const std::string &formulation)
	{
		//Uncomment to used not adaptive tetgen mesher
// 		clear();

// 		Eigen::Vector2d max_dim;
// 		Eigen::Vector2d min_dim;
// 		mesh.get_physical_bounding_box(min_dim, max_dim);

// 		double xMin = min_dim(0) - padding_size;
// 		double xMax = max_dim(0) + padding_size;
// 		double yMin = min_dim(1) - padding_size;
// 		double yMax = max_dim(1) + padding_size;
// 		double zMin = -thickness;
// 		double zMax = 0;

// 		V.resize(8, 3);
// 		V << xMin, yMin, zMin,
// 		xMin, yMax, zMin,
// 		xMax, yMax, zMin,
// 		xMax, yMin, zMin,

// 			//4
// 		xMin, yMin, zMax,
// 		xMin, yMax, zMax,
// 		xMax, yMax, zMax,
// 		xMax, yMin, zMax;

// 		F.resize(12, 3);
// 		F << 1, 2, 0,
// 		0, 2, 3,

// 		5, 4, 6,
// 		4, 7, 6,

// 		1, 0, 4,
// 		1, 4, 5,

// 		2, 1, 5,
// 		2, 5, 6,

// 		3, 2, 6,
// 		3, 6, 7,

// 		0, 3, 7,
// 		0, 7, 4;

// 		Eigen::MatrixXd TV;
// 		Eigen::MatrixXi TT;
// 		Eigen::MatrixXi TF;
// #ifdef NDEBUG
// 		igl::copyleft::tetgen::tetrahedralize(V, F, "Qpq1.414a100", TV, TT, TF);
// #else
// 		igl::copyleft::tetgen::tetrahedralize(V, F, "Qpq1.414a1000", TV, TT, TF);
// #endif
// 		V = TV;
// 		F = TF;
// 		T = TT;

		simulation_out = compute_analysis(V, F, T, mesh, thickness, E, nu, formulation, scaling, displacement, traction_forces);
	}

	void Mesh3d::clear()
	{
		F.resize(0, 0);
		V.resize(0, 0);

		T.resize(0, 0);
		displacement.resize(0, 0);
		traction_forces.resize(0, 0);

		simulation_out = nlohmann::json({});
	}

	bool Mesh3d::load(const nlohmann::json & data)
	{
		read_json_mat(data["V"], V);
		read_json_mat(data["F"], F);
		//read_json_mat(data["T"], T); // maybe the Tets are unnecessary

		read_json_mat(data["displacement"], displacement);
		read_json_mat(data["traction_forces"], traction_forces);

		return true;
	}

	void Mesh3d::save_mesh(nlohmann::json & data)
	{
		data["V"] = json::object();
		write_json_mat(V, data["V"]);

		data["F"] = json::object();
		write_json_mat(F, data["F"]);

		//data["T"] = json::object();
		//write_json_mat(T, data["T"]);

		data["displacement"] = json::object();
		write_json_mat(displacement, data["displacement"]);
	}
	void Mesh3d::save_traction(nlohmann::json & data)
	{
		data["traction_forces"] = json::object();
		write_json_mat(traction_forces, data["traction_forces"]);
	}

}// namespace cellogram


