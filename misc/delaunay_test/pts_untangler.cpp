#include "pointsUntangler/PointsUntangler.h"

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <igl/list_to_matrix.h>
#include <igl/write_triangle_mesh.h>

bool load(const std::string & path, Eigen::MatrixXd &res)
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
	Eigen::MatrixXd pts;
	Eigen::MatrixXi F;
	Eigen::MatrixXd newPts;

	const std::string root = DATA_DIR;

	load("/Users/teseo/GDrive/Cellogram/Images/img12/vDetected.xyz", pts);
	if(pts.cols() == 2){
		pts.conservativeResize(pts.rows(), 3);
		pts.col(2).setZero();
	}

	cellogram::PointsUntangler::pointsUntangler(pts, F, newPts);

	Eigen::MatrixXd total(pts.rows()+newPts.rows(), pts.cols());
	// total.setZero();

	total.block(0, 0, pts.rows(), pts.cols()) = pts;
	total.block(pts.rows(), 0, newPts.rows(), newPts.cols()) = newPts;

	igl::write_triangle_mesh("test.obj", total, F);


	return 0;
}
