#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <set>

using namespace Eigen;

class generateQ
{
public:

    // Input
	MatrixXd V; // The first vertices are the boundary (sorted as in V_boundary), then the internals (unsorted)
	MatrixXi F;
	
	// Calculations
	int iNrT;
	int iNrV;

	MatrixXi A; // Adjacency matrix
	MatrixXi L; // Laplacian matrix

	int K;
	MatrixXi IDX; //Matrix containing indices of closest neighbors

	SparseMatrix<double> Q;
	SparseMatrix<int> Aeq;

	// Generate Laplacian
	void adjacencyMatrix(const MatrixXi &triangles);
	void laplacianMatrix();

	// Generate Q
	void QforOptimization(const MatrixXd &Vperfect, const MatrixXd &Vdeformed, int K);

	// Compile constraints for opimization
	void optimizationConstraints(int nrBoundaryV);
};

