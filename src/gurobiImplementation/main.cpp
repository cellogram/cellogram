
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
  s.load(string(string(path)));

  // Merge V_boundary and V_internal to contain all original coordinates in a matrix
  MatrixXd Vdeformed(s.V_boundary.rows() + s.V_internal.rows(), 2);
  Vdeformed << s.V_boundary, s.V_internal;

  // Generate perfect mesh in void
  s.fill_hole();

  s.save(string(string(path)));

  cout << "T: " << endl << s.F << endl;

  // Generate laplacian
  Q.adjacencyMatrix(s.F);
  Q.laplacianMatrix();

  // Deriving Q for optimization
  Q.QforOptimization(s.V, Vdeformed, 15);

  // Gurobi optimization test
  /*try {
	  GRBEnv env = GRBEnv();

	  GRBModel model = GRBModel(env);

	  // Create variables

	  GRBVar x = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x");
	  GRBVar y = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y");
	  GRBVar z = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "z");

	  // Set objective: maximize x + y + 2 z

	  model.setObjective(x + y + 2 * z, GRB_MAXIMIZE);

	  // Add constraint: x + 2 y + 3 z <= 4

	  model.addConstr(x + 2 * y + 3 * z <= 4, "c0");

	  // Add constraint: x + y >= 1

	  model.addConstr(x + y >= 1, "c1");

	  // Optimize model

	  model.optimize();

	  cout << x.get(GRB_StringAttr_VarName) << " "
		  << x.get(GRB_DoubleAttr_X) << endl;
	  cout << y.get(GRB_StringAttr_VarName) << " "
		  << y.get(GRB_DoubleAttr_X) << endl;
	  cout << z.get(GRB_StringAttr_VarName) << " "
		  << z.get(GRB_DoubleAttr_X) << endl;

	  cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

  }
  catch (GRBException e) {
	  cout << "Error code = " << e.getErrorCode() << endl;
	  cout << e.getMessage() << endl;
  }
  catch (...) {
	  cout << "Exception during optimization" << endl;
  }
  */
  cin.get();
}
