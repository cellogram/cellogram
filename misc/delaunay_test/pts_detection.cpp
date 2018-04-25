#include "cellogram/point_source_detection.h"

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <igl/list_to_matrix.h>

bool load_img(const std::string & path, Eigen::MatrixXd &res)
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

	if (!igl::list_to_matrix(matrix, res))
	{
		std::cerr << "list to matrix error" << std::endl;
		file.close();
		return false;
	}
	return true;
}


int main(int argc, char** argv) {
	Eigen::MatrixXd img;
	const double sigma = 2;
	Eigen::MatrixXd V;
	Eigen::MatrixXd V_std;
	Eigen::VectorXd pval_Ar;

	const std::string root = DATA_DIR;

	load_img(root + "1-1.txt", img);

	double min = img.minCoeff();
	double max = img.maxCoeff();

	Eigen::MatrixXd imgNorm;

	imgNorm = (img.array() - min) / (max - min);

	cellogram::point_source_detection(imgNorm, sigma, V, V_std, pval_Ar);


	return 0;
}
