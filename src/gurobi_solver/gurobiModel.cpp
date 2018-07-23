
#include "gurobiModel.h"
#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;



void gurobiModel::model(const SparseMatrix<double> &Q, const SparseMatrix<int> &Aeq, double time_limit)
{
	/* Mixed Integer Programming
	Objective: 	    minimize xT Q x + qT x
	Constraints: 	A x = b (linear constraints)
	*/
	bool solution_found = false;

	try
	{
		GRBEnv env = GRBEnv();
		GRBModel model = GRBModel(env);

		model.set(GRB_DoubleParam_TimeLimit, time_limit);

		int nVars = Aeq.cols();
		int nConstraints = Aeq.rows();

		// determine variables
		GRBVar *x = model.addVars(nVars, GRB_BINARY);
		model.update();

		// objective function x'Qx
		GRBQuadExpr Obj = 0;
		for (int k = 0; k < Q.outerSize(); ++k)
		{
			for (SparseMatrix<double>::InnerIterator it(Q, k); it; ++it)
			{
				Obj += it.value() * x[it.row()] * x[it.col()];
			}
		}
		model.setObjective(Obj);
		//model.set(GRB_IntAttr_ModelSense, 1);

		//cout << endl << Q.nonZeros() << " " << Q.rows() << " " << Q.cols();

		// constraint declaration Aeq x = 1
		//cout << endl << endl << Aeq.row(0).col(1);
		SparseMatrix<int> AeqT = Aeq.transpose();
		for (int k = 0; k < AeqT.outerSize(); ++k)
		{
			GRBLinExpr LHS = 0;
			for (SparseMatrix<int>::InnerIterator it(AeqT, k); it; ++it)
			{
				LHS += it.value() * x[it.row()];
			}
			model.addConstr(LHS, GRB_EQUAL, 1);
		}

		// optimize model

		model.optimize();
		int optimstatus = model.get(GRB_IntAttr_Status);
		cout << "Optimization complete" << endl;
		double objval = 0;
		if (optimstatus == GRB_OPTIMAL) {
			objval = model.get(GRB_DoubleAttr_ObjVal);
			solution_found = true;
		}
		else if (optimstatus == GRB_INF_OR_UNBD) {
			cout << "Model is infeasible or unbounded" << endl;
		}
		else if (optimstatus == GRB_INFEASIBLE) {
			cout << "Model is infeasible" << endl;
		}
		else if (optimstatus == GRB_UNBOUNDED) {
			cout << "Model is unbounded" << endl;
		}
		else {
			cout << "Optimization was stopped with status = "
				<< optimstatus << endl;

			//// Time limit has been reached, but solution may be good
			//if (optimstatus == GRB_TIME_LIMIT)
			//{
			//	objval = model.get(GRB_DoubleAttr_ObjVal);
			//	solution_found = true;
			//}
		}



		// Iterate over the solutions and compute the objectives
		if (solution_found) {
			resultX.resize(Q.rows());
			for (int i = 0; i < Q.rows(); i++) {
				resultX(i) = std::round(x[i].get(GRB_DoubleAttr_X));
			}
		}
		else {
			resultX.resize(1);
			resultX(0) = -1;
		}

	}
	catch (GRBException e)
	{
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...)
	{
		cout << "Exception during optimization" << endl;
	}

	//cout << resultX;
}
