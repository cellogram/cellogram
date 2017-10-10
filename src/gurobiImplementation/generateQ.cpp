
#include "generateQ.h"
#include <iostream>
#include <Eigen/Dense>
#include <geogram/points/kd_tree.h>
#include <igl/repmat.h>

using namespace std;
using namespace Eigen;



void generateQ::adjacencyMatrix(const MatrixXi &triangles)
{
	F = triangles;
	iNrT = F.rows();
	iNrV = F.maxCoeff() + 1;

	//cout << "iNrV: " << iNrV << endl;

	A = MatrixXi::Zero(iNrV, iNrV);

	for (int i = 0; i < iNrT; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			A(F(i, j), F(i, 0)) = 1;
			A(F(i, j), F(i, 1)) = 1;
			A(F(i, j), F(i, 2)) = 1;
			A(F(i, 0), F(i, j)) = 1;
			A(F(i, 1), F(i, j)) = 1;
			A(F(i, 2), F(i, j)) = 1;
		}

	}
	for (int i = 0; i < iNrV; i++)
	{
		A(i, i) = 0;
	}

	//cout << "A: " << A << endl;
}

void generateQ::laplacianMatrix()
{
	L = MatrixXi::Zero(iNrV, iNrV);
	for (int i = 0; i < iNrV; i++)
	{
		L(i, i) = A.row(i).sum();
	}
	L = L - A;

	// cout << "L: " << L << endl;
}

void generateQ::QforOptimization(const MatrixXd &Vperfect, const MatrixXd &Vdeformed, int K)
{
	// first check the number of closest neighbors taken into account for meshing and adjust if larger than the available number
	if (iNrV < K)
	{
		K = iNrV;
	}

	// generate a kdtree object with the original (deformed) coordinates
	GEO::NearestNeighborSearch_var nnsearch =
		GEO::NearestNeighborSearch::create((GEO::coord_index_t) 2, "BNN");

	MatrixXd B = Vdeformed.transpose();
	cout << "B: " << endl << B << endl;
	cout << Vperfect.transpose() << endl;
	nnsearch->set_points((GEO::index_t)Vdeformed.rows(), B.data());

	// Query
	vector<GEO::index_t> nearest(K);
	vector<double> dist(K);

	VectorXi Qind1(K*iNrV), Qind2(K*iNrV);
	VectorXd Vx(K*iNrV), Vy(K*iNrV);
	IDX = MatrixXi(K, iNrV);

	for (int i = 0; i < iNrV; i++)
	{
		Vector2d x(Vperfect(i,0), Vperfect(i,1));

		cout << "x: " << x.transpose() << endl;

		nnsearch->get_nearest_neighbors((GEO::index_t)K, x.data(), nearest.data(), dist.data());
		// Save query results
		for (int j = 0; j < K; j++)
		{
			Qind1(i*K + j) = nearest[j];
			Qind2(i*K + j) = i;
			Vx(i*K + j) = Vdeformed(nearest[j], 0);
			Vy(i*K + j) = Vdeformed(nearest[j], 1);
			IDX(j,i) = nearest[j];
		}

	}
	cout << "IDX: " << endl << IDX << endl;

	//VK_x = sparse(Q1, Q2, Vx);
	//VK_y = sparse(Q1, Q2, Vy);
	//Q = VK_x'*(L')*L*VK_x + VK_y'*(L')*L*VK_y;

}

/*

  //
  VectorXd Qind1(K*iNrV), Qind2(K*iNrV), Vx(K*iNrV), Vy(K*iNrV);
  MatrixXd IDX(K,iNrV);

  for (int i = 0; i < iNrV; i++)
  {
	  // K_closest = knnsearch(KDTree_iV,[V(i,1),V(i,2)],'k',K);

  }
  Triplet<double> T;
  vector<T> tripletListX, tripletListY;

  for (int i = 0; i < Qind1.cols(); i++)
  {
	  tripletListX.push_back(T(Qind1(i), Qind2(i), Vx(i)));
	  tripletListY.push_back(T(Qind1(i), Qind2(i), Vy(i)));
  }


  SparseMatrix<double> VK_x.setFromTriplet(tripletListX.begin(), tripletListX.end());
  SparseMatrix<double> VK_y.setFromTriplet(tripletListY.begin(), tripletListY.end());
  */
