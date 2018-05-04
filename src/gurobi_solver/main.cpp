
#include <string>
#include <iostream>
#include <fstream>
#include "state.h"
#include "generateQ.h"
#include "gurobiModel.h"
#include <geogram/basic/common.h>

using namespace std;
using namespace Eigen;

State s;
generateQ q;
gurobiModel g;

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
  q.adjacencyMatrix(s.F);
  q.laplacianMatrix();

  // Deriving Q and constraints for optimization
  q.QforOptimization(s.Vperfect, s.Vdeformed, 8);
  q.optimizationConstraints(s.V_boundary.rows());

  // Generate and solve model for gurobi
  g.model(q.Q, q.Aeq, 10);

  // Map back to indices of coordinates
  q.mapBack(g.resultX);

  #ifdef WIN32
  cin.get();
  #endif
}
