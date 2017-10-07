
#include <string>
#include <iostream>
#include <fstream>
#include "state.h"
#include "generateQ.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

State s;
generateQ L;

int main()
{
	
  string path = "C:\\Users\\letobias\\Documents\\cellogram\\src\\gurobiImplementation\\Data\\1";

  cout << path << endl;

  // Load the data from text files
  s.load(string(string(path)));

  // Generate perfect mesh in void
  s.fill_hole();

  s.save(string(string(path)));

  cout << s.F << endl;

  // Generate laplacian
  L.adjacencyMatrix(s.F);
  L.laplacianMatrix();


  /*
  // Deriving Q for optimization
  int K = 20; 
  if (iNrV < K)
  {
	  K = iNrV;
  }

  // create KDtree from original coordinates
  
  //
  VectorXd Qind1(K*iNrV), Qind2(K*iNrV), Vx(K*iNrV), Vy(K*iNrV);
  MatrixXd IDX(K,iNrV);
  
  for (int i = 0; i < iNrV; i++)
  {
	  // K_closest = knnsearch(KDTree_iV,[V(i,1),V(i,2)],'k',K);
	  // Qind1.segment(i*K, (i + 1)*K - 1) = K_closest;
	  // Vx.segment(i*K, (i + 1)*K - 1) = s.V(K_closest, 1);
	  // Vy.segment(i*K, (i + 1)*K - 1) = s.V(K_closest, 2);
	  // IDX.cols(i) = K_closest;
  }
  Triplet<double> T;
  vector<T> tripletListX, tripletListY;

  for (int i = 0; i < Qind1.cols(); i++)
  {
	  tripletListX.push_back(T(Qind1(i), Qind2(i), Vx(i)));
	  tripletListY.push_back(T(Qind1(i), Qind2(i), Vy(i)));
  }


  //SparseMatrix<double> VK_x.setFromTriplet(tripletListX.begin(), tripletListX.end());
  SparseMatrix<double> VK_y.setFromTriplet(tripletListY.begin(), tripletListY.end());
  */
  system("pause");
}
