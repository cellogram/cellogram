////////////////////////////////////////////////////////////////////////////////
#include "delaunay.h"
#include <geogram/mesh/mesh.h>
#include <geogram/delaunay/delaunay.h>
#include <igl/adjacency_matrix.h>
#include <igl/sum.h>
#include <igl/diag.h>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

// -----------------------------------------------------------------------------

void laplace_energy(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::VectorXd &E){
	assert(V.cols() == 2 || V.cols() == 3);
	assert(F.cols() == 3);
	int n = (int)V.rows();

	// Compute the uniform laplacian
	Eigen::SparseMatrix<double> A;
	igl::adjacency_matrix(F, A);
	// sum each row
	Eigen::SparseVector<double> Asum;
	igl::sum(A, 1, Asum);
	// Convert row sums into diagonal of sparse matrix
	Eigen::SparseMatrix<double> Adiag;
	igl::diag(Asum, Adiag);
	// Build uniform laplacian
	Eigen::SparseMatrix<double> L;
	L = A - Adiag;

	{
		Eigen::VectorXd Ex = L * V.col(0);
		Eigen::VectorXd Ey = L * V.col(1);
		E = Ex.array().abs() + Ey.array().abs();
	}

	double Emean = E.array().mean();

	// Normalize Laplacian Energy from 0 to 2 * average, and set everything above to 1
	Eigen::VectorXd Enorm = Eigen::VectorXd::Zero(n);
	for (size_t i = 0; i < n; i++)
	{
		if (E(i) < (2 * Emean))
			Enorm(i) = E(i) / (2 * Emean);
		else
			Enorm(i) = 1;
	}

	E = Enorm;

}

// -----------------------------------------------------------------------------

} // namespace cellogram
