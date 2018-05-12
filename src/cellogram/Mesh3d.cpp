////////////////////////////////////////////////////////////////////////////////
#include "Mesh3d.h"
#include "Mesh.h"
#include <Eigen/Dense>

#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <geogram/mesh/mesh_io.h>
#include <State.hpp>
#include <CustomProblem.hpp>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

	namespace
	{
		void compute_analysis(const Eigen::MatrixXd &vertices, const Eigen::MatrixXi &tets, const Mesh &mesh, float thickness, Eigen::MatrixXd &vals)
		{
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

			json j_args = {
				{"problem", "Custom"},
				{"normalize_mesh", false},

				{"tensor_formulation", "LinearElasticity"},

				{"discr_order", 1},

				{"params", {
					{"lambda", 0.75},
					{"mu", 0.375},
				}},
			};

			poly_fem::State &state = poly_fem::State::state();
			state.init(j_args);


			state.load_mesh(M, [thickness](const poly_fem::RowVectorNd &bary){
				//TOP
				if(bary(2) < 1e-8)
					return 1;

				//Bottom
				if(std::abs(bary(2)- -thickness) < 1e-8)
					return 3;

				//any other
				return 2;
			});
			// state.compute_mesh_stats();

			//TODO flags
			poly_fem::CustomProblem &problem = *dynamic_cast<poly_fem::CustomProblem *>(state.problem.get());
			problem.init({1, 3});


			Eigen::MatrixXd disp = (mesh.detected - mesh.points) * mesh.scaling;
			Eigen::MatrixXd pts = mesh.points * mesh.scaling;

			problem.set_function(0, disp, pts, mesh.triangles); //0 is id 1
			problem.set_constant(1, Eigen::Vector3d(0,0,0)); //1 is id 3

			state.build_basis();
			state.build_polygonal_basis();


			state.assemble_rhs();
			state.assemble_stiffness_mat();

			state.solve_problem();

			// state.interpolate_function(vertices.rows(), state.sol, vals);
			//vals = state.sol;
			// vals = Eigen::Map<Eigen::MatrixXd>(state.sol.data(), 3, vertices.rows());
			vals = Eigen::Map<Eigen::MatrixXd>(state.rhs.data(), 3, vertices.rows());
			vals = vals.transpose().eval();
		}
	}

	void Mesh3d::init(const Mesh &mesh, float padding_size, float thickness)
	{
		Eigen::Vector2d max_dim;
		Eigen::Vector2d min_dim;
		mesh.get_physical_bounding_box(min_dim, max_dim);

		double xMin = min_dim(0) - padding_size;
		double xMax = max_dim(0) + padding_size;
		double yMin = min_dim(1) - padding_size;
		double yMax = max_dim(1) + padding_size;
		double zMin = -thickness;
		double zMax = 0;

		V.resize(8, 3);
		V << xMin, yMin, zMin,
		xMin, yMax, zMin,
		xMax, yMax, zMin,
		xMax, yMin, zMin,

			//4
		xMin, yMin, zMax,
		xMin, yMax, zMax,
		xMax, yMax, zMax,
		xMax, yMin, zMax;

		F.resize(12, 3);
		F << 1, 2, 0,
		0, 2, 3,

		5, 4, 6,
		4, 7, 6,

		1, 0, 4,
		1, 4, 5,

		2, 1, 5,
		2, 5, 6,

		3, 2, 6,
		3, 6, 7,

		0, 3, 7,
		0, 7, 4;

		Eigen::MatrixXd TV;
		Eigen::MatrixXi TT;
		Eigen::MatrixXi TF;
		igl::copyleft::tetgen::tetrahedralize(V, F, "Qpq1.414a100000", TV, TT, TF);

		std::cout<<"n tets: "<<TF.rows()<<std::endl;

		compute_analysis(TV, TT, mesh, thickness, sol);

		F = TF;
		V = TV;

	}

	void Mesh3d::clear()
	{
		F.resize(0, 0);
		V.resize(0, 0);
	}

}// namespace cellogram


