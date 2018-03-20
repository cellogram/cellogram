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
	Eigen::MatrixXd points; // displaced points

	Eigen::MatrixXd hull;

public:
	// void foo();
};

} // namespace cellogram