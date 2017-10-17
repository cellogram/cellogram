
#include "gurobiModel.h"
#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;



void gurobiModel::model(const SparseMatrix<double> &Q, const SparseMatrix<int> &Aeq)
{
	/* Mixed Integer Programming
	Objective: 	    minimize xT Q x + qT x
    Constraints: 	A x = b (linear constraints)
	*/
	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);

	int nVars = Aeq.cols();
	int nConstraints = Aeq.rows();

	// determine variables
	GRBVar *x = model.addVars(nVars, GRB_BINARY);
	model.update();

	// objective function x'Qx
	model.set(GRB_IntAttr_ModelSense, 1);
	for (int i = 0; i < nVars; i++)
	{
		GRBQuadExpr Obj = 0;
		for (int j = 0; j < nVars; j++)
		{
			Obj += x[i]*Q.row(i).col(j)*x[j];
		}
		model.setObjective(Obj);
	}

	// constraint declaration Aeq x = 1
	//cout << endl << endl << Aeq.row(0).col(1);
	for (int i = 0; i < nConstraints; i++)
	{
		GRBLinExpr LHS = 0;
		for (int j = 0; j < nVars; j++)
		{
			LHS += Aeq.row(i).col(j)*x[j];
		}
		model.addConstr(LHS,GRB_EQUAL, 1);
	}

}

void gurobiModel::solve(const SparseMatrix<double> &Q, const SparseMatrix<int> &Aeq)
{

}
