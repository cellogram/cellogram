////////////////////////////////////////////////////////////////////////////////
#include "common.h"
#include <vector>
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

struct State {
public:
	static State &state();
	~State() = default;

public:
	Eigen::MatrixXd detected; // detected (unmoved) point positions
	Eigen::MatrixXd points; // relaxed point positions

	Eigen::MatrixXi triangles; // triangular mesh
	std::vector<std::vector<int>> adjacency_list; // adjaceny list of triangluar mesh

	Eigen::VectorXi boundary; // list of vertices on the boundary

	Eigen::MatrixXd hull_vertices;
	Eigen::MatrixXi hull_faces;

	Eigen::VectorXd original_laplace_energy; // Laplace energy at vertices for the ORIGINAL positions with CURRENT mesh
	Eigen::VectorXd current_laplace_energy; // Laplace energy at vertices for the CURRENT positions with CURRENT mesh

public:
	// void foo();
};

} // namespace cellogram
