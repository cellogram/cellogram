#include <cellogram/remesh_adaptive.h>
#include <igl/triangle/triangulate.h>
#include <igl/write_triangle_mesh.h>
#include <Eigen/Dense>
#include <iostream>

void isotropic_quad(Eigen::MatrixXd &OV, Eigen::MatrixXi &OF, double area) {
	std::stringstream buf;
	buf.precision(100);
	buf.setf(std::ios::fixed, std::ios::floatfield);

	buf << "Qqa" << area;

	Eigen::MatrixXd V(4,2); V <<
		-1,-1,
		-1,1,
		1,1,
		1,-1;

	Eigen::MatrixXi E(4,2); E <<
		0,1,
		1,2,
		2,3,
		3,0;

	Eigen::MatrixXd H(0,2);
	igl::triangle::triangulate(V, E, H, buf.str(), OV, OF);
}

int main(int argc, char** argv) {
	Eigen::MatrixXd V1, V2;
	Eigen::MatrixXi F1, F2;
	Eigen::VectorXd S;

	isotropic_quad(V1, F1, 0.01);

	S.resize(V1.rows());
	for (int v = 0; v < V1.rows(); ++v) {
		Eigen::RowVector2d p = V1.row(v);
		S(v) = 0.001 + 0.1 * p.squaredNorm();
	}

	cellogram::remesh_adaptive_2d(V1, F1, S, V2, F2);
	igl::write_triangle_mesh("output.obj", V2, F2);

	return 0;
}
