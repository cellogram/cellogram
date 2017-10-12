
#include <string>
#include <iostream>
#include <fstream>
#include "state.h"
#include "generateQ.h"
#include <Eigen/Dense>
#ifdef WIN32
#include "C:/gurobi751/win64/include/gurobi_c++.h"
#endif
#include <geogram/basic/common.h>

using namespace std;
using namespace Eigen;

State s;
generateQ Q;

int main()
{
  GEO::initialize();
  string path = DATA_DIR "1";

  cout << path << endl;

  // Load the data from text files
  s.load(path);

  // Generate perfect mesh in void
  s.fill_hole();

  s.save(path);

  cout << "perfect: " << endl << s.Vperfect << endl;
  cout << "deformed: " << endl << s.Vdeformed << endl;

  // Generate adjacency matrix and the laplacian
  Q.adjacencyMatrix(s.F);
  Q.laplacianMatrix();

  // Deriving Q for optimization
  Q.QforOptimization(s.Vperfect, s.Vdeformed, 15);
  
  #ifdef WIN32
  cin.get();
  #endif
}
