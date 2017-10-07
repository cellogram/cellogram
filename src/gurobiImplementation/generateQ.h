#include <Eigen/Dense>
#include <set>

using namespace Eigen;

// Important convenction to store a hex grid in a matrix (see scheme.png)

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


	// Generate Laplacian
	void adjacencyMatrix(const MatrixXi &triangles);
	void laplacianMatrix();

	// Generate Q
	void QforOptimization(const MatrixXd &V0, int K);
};

