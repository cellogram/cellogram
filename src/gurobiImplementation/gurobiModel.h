#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <gurobi_c++.h>

using namespace Eigen;

class gurobiModel
{
public:

    // Output
	VectorXd resultX;

	// Generate and solve Model
	void model(const SparseMatrix<double> &Q, const SparseMatrix<int> &Aeq);

};

