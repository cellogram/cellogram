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

	Eigen::MatrixXi triangles;

	Eigen::VectorXi boundary; // list of vertices on the boundary

	Eigen::MatrixXd hull_vertices;
	Eigen::MatrixXi hull_faces;

public:
	// void foo();
};

} // namespace cellogram
