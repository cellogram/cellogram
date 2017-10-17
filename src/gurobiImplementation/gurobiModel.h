#include <Eigen/Dense>
#include <Eigen/Sparse>
#ifdef WIN32
#include "C:/gurobi751/win64/include/gurobi_c++.h"
#else
#include "gurobi_c++.h"
#endif

using namespace Eigen;

class gurobiModel
{
public:

    // Input
	//MatrixXd V;

	// Generate Model
	void model(const SparseMatrix<double> &Q, const SparseMatrix<int> &Aeq);

	// Solve Model
	void solve(const SparseMatrix<double> &Q, const SparseMatrix<int> &Aeq);

};

