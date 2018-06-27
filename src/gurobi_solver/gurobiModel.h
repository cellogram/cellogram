#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <gurobi_c++.h>

class gurobiModel
{
public:

    // Output
	Eigen::VectorXi resultX;

	// Generate and solve Model
	void model(const Eigen::SparseMatrix<double> &Q, const Eigen::SparseMatrix<int> &Aeq, double time_limit);

};

