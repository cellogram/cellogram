////////////////////////////////////////////////////////////////////////////////
#include "Mesh3d.h"
#include "Mesh.h"
#include <Eigen/Dense>

#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <geogram/mesh/mesh_io.h>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

	namespace
	{
		
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
		igl::copyleft::tetgen::tetrahedralize(V, F, "Qpq1.414a100", TV, TT, TF);

		std::cout<<"n tets: "<<TF.rows()<<std::endl;

		GEO::Mesh M;
		M.vertices.create_vertices((int) TV.rows());
		for (int i = 0; i < (int) M.vertices.nb(); ++i) {
			GEO::vec3 &p = M.vertices.point(i);
			p[0] = TV(i, 0);
			p[1] = TV(i, 1);
			p[2] = TV(i, 2);
		}

		M.cells.create_tets((int) TT.rows());

		for (int c = 0; c < (int) M.cells.nb(); ++c) {
			for (int lv = 0; lv < TT.cols(); ++lv) {
					M.cells.set_vertex(c, lv, TT(c, lv));
			}
		}
		GEO::mesh_save(M, "mesh.mesh");

		F = TF;
		V = TV;

	}

}// namespace cellogram


