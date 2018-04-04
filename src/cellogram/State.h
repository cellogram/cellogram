////////////////////////////////////////////////////////////////////////////////
#include "common.h"
#include "mesh_solver.h"
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

	Eigen::VectorXi boundary; // list of vertices on the boundary

	Eigen::MatrixXd hull_vertices;
	Eigen::MatrixXi hull_faces;

	mesh_solver regions;

public:
	// void foo();
};

} // namespace cellogram
