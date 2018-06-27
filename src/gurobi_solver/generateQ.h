#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <set>

class generateQ
{
public:

    // Input
	Eigen::MatrixXd V; // The first vertices are the boundary (sorted as in V_boundary), then the internals (unsorted)
	Eigen::MatrixXi F;

	// Calculations
	int iNrT;
	int iNrV;

	Eigen::MatrixXi A; // Adjacency matrix
	Eigen::MatrixXi L; // Laplacian matrix

	int K;
	Eigen::MatrixXi IDX; //Matrix containing indices of closest neighbors

	Eigen::SparseMatrix<double> Q;
	Eigen::SparseMatrix<int> Aeq;

	Eigen::VectorXi vMapping; // Vector containing mapping from F (perfect triangles) to input indices (V)
	Eigen::MatrixXi T; // triangular connectivity of V

	// Generate Laplacian
	void adjacencyMatrix(const Eigen::MatrixXi &triangles);
	void laplacianMatrix();

	// Generate Q
	void QforOptimization(const Eigen::MatrixXd &Vperfect, const Eigen::MatrixXd &Vdeformed, int K);

	// Compile constraints for opimization
	void optimizationConstraints(int nrBoundaryV);

	// Map the results from gurobi back to the indices of the coordinates
	void mapBack(const Eigen::VectorXi &x);
};

