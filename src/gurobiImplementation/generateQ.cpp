
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

void generateQ::QforOptimization(const MatrixXd &Vperfect, const MatrixXd &Vdeformed, int nrNeigh)
{
	K = nrNeigh;
	// first check the number of closest neighbors taken into account for meshing and adjust if larger than the available number
	if (iNrV < K)
	{
		K = iNrV;
	}

	// generate a kdtree object with the original (deformed) coordinates
	GEO::NearestNeighborSearch_var nnsearch =
		GEO::NearestNeighborSearch::create((GEO::coord_index_t) 2, "BNN");

	MatrixXd B = Vdeformed.transpose();
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

		nnsearch->get_nearest_neighbors((GEO::index_t)K, x.data(), nearest.data(), dist.data());
		// Save query results
		for (int j = 0; j < K; j++)
		{
			Qind1(i*K + j) = i;
			Qind2(i*K + j) = i*K + j;
			Vx(i*K + j) = Vdeformed(nearest[j], 0);
			Vy(i*K + j) = Vdeformed(nearest[j], 1);
			IDX(j,i) = nearest[j];
		}

	}
	
	/* Checked with working code in matlab:
	Qind1,Qind2,Vx,Vy,IDX are all correctly assembled */

	// generate triplets for sparse matrices
	typedef Eigen::Triplet<double> T;
	std::vector<T> VxList,VyList;
	VxList.reserve(Qind1.rows());
	VyList.reserve(Qind1.rows());

	for (int i = 0; i < Qind1.rows(); i++)
	{
		VxList.push_back(T(Qind1(i), Qind2(i), Vx(i)));
		VyList.push_back(T(Qind1(i), Qind2(i), Vy(i)));
	}

	// generate sparse matrices
	SparseMatrix<double> VK_x(iNrV, Qind1.rows());
	SparseMatrix<double> VK_y(iNrV, Qind1.rows());
	SparseMatrix<double> Ls(L.rows(), L.cols());
	VK_x.setFromTriplets(VxList.begin(), VxList.end());
	VK_y.setFromTriplets(VyList.begin(), VyList.end());
	
	{
		MatrixXd B = L.cast<double>();
		Ls = B.sparseView();
	}
	
	/* Checked with working code in matlab:
	VK_x and VK_y are correct */
	
	// Calculate Q for the optimization:
	// Q = VK_x'*(L')*L*VK_x + VK_y'*(L')*L*VK_y
	Q = VK_x.transpose()*Ls.transpose()*Ls*VK_x + VK_y.transpose()*Ls.transpose()*Ls*VK_y;
	
	/* Checked with working code in matlab:
	Q ist mostly correct with differences of <<1%  */

}

void generateQ::optimizationConstraints(int nrBoundaryV)
{
	// Aeq x = beq
	// constraint 1: only one entry per row of M
	// constraint 2: only one entry per column of M

	// generate triplets for sparse matrices
	typedef Eigen::Triplet<double> T;
	std::vector<T> AeqList;
	AeqList.reserve(2*K*iNrV + nrBoundaryV);


	for (int i = 0; i < iNrV; i++)
	{
		for (int j = 0; j < K; j++)
		{
			AeqList.push_back(T(i, i*K + j, 1));
			AeqList.push_back(T(IDX(j, i) + iNrV,i*K+j, 1));
		}
	}

	// constraint 3 : vertices on cluster edge remain unmoved
	for (int i = 0; i < nrBoundaryV; i++)
	{
		AeqList.push_back(T(i + 2*iNrV, i*K, 1));
	}

	/* AeqList is correct*/
	//for (int i = 0; i < AeqList.size(); i++)
	//{
	//	cout << AeqList[i].row() << "," << AeqList[i].col() << "," << AeqList[i].value() << endl;
	//}

	//Concatenate all Aeq into single sparse matrix
	Aeq = SparseMatrix<int>(2*iNrV + nrBoundaryV, K*iNrV);

	Aeq.setFromTriplets(AeqList.begin(), AeqList.end());

	/*constraints are correct, checked with working code in matlab*/
	//cout << "Aeq: " << endl << Aeq << endl;
}