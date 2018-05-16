#include "cellogram/point_source_detection.h"
#include <gurobi_solver/state.h>
#include <gurobi_solver/generateQ.h>
#include <gurobi_solver/gurobiModel.h>


#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <igl/list_to_matrix.h>
#include <igl/Timer.h>

bool load_data(const std::string & path, Eigen::MatrixXd &mat)
{
	std::fstream file;
	file.open(path.c_str());

	if (!file.good())
	{
		std::cerr << "Failed to open file : " << path << std::endl;
		file.close();
		return false;
	}


	std::string s;
	std::vector<std::vector<double>> matrix;

	while (getline(file, s))
	{
		std::stringstream input(s);
		double temp;
		matrix.emplace_back();

		std::vector<double> &currentLine = matrix.back();

		while (input >> temp)
			currentLine.push_back(temp);
	}

	if (!igl::list_to_matrix(matrix, mat))
	{
		std::cerr << "list to matrix error" << std::endl;
		file.close();
		return false;
	}
	return true;
}


int main(int argc, char** argv) {
	Eigen::MatrixXd vB;
	Eigen::MatrixXd vI;
	Eigen::MatrixXd tmp;
	
	const std::string root = "C:/Users/Tobias/Dropbox/NY/test_regions/";
	const std::string file = "1";
	load_data(root + file + "/vB.vert", vB);
	load_data(root + file + "/vI.vert", vI);
	load_data(root + file + "/neigh.txt", tmp);

	std::cout << "\nvB:\n" << vB.transpose() << std::endl;
	std::cout << "\nvI:\n" << vI.transpose() << std::endl;
	std::cout << "\nneigh:\n" << tmp.transpose() << std::endl;

	Eigen::MatrixXi tmpi = tmp.cast<int>();

	Eigen::VectorXi neigh(Eigen::Map<Eigen::VectorXi>(tmpi.data(), tmpi.cols()*tmpi.rows()));

	// Perparing data
	State s;
	s.init(vB, vI, neigh);
	s.fill_hole();

	// Generate adjacency matrix and the laplacian
	generateQ q;
	gurobiModel g;
	q.adjacencyMatrix(s.F);
	q.laplacianMatrix();

	int perm_possibilities = 10;
	// Deriving Q and constraints for optimization
	int K = std::min(perm_possibilities, q.iNrV);
	q.QforOptimization(s.Vperfect, s.Vdeformed, K);
	q.optimizationConstraints(s.V_boundary.rows());

	int gurobi_time_limit = 60;

	// Generate and solve model for gurobi
	g.model(q.Q, q.Aeq, gurobi_time_limit);
	if (g.resultX(0) == -1) {
		// no solution found
		return -1;
	}

	// Map back to indices of coordinates
	q.mapBack(g.resultX);

	// Google solver

	return 0;
}
