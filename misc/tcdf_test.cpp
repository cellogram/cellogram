#include "tcdf/tcdf.h"

#include <Eigen/Dense>
#include <iostream>

int main(int argc, char** argv) {
	Eigen::MatrixXd x(5,6);
	Eigen::MatrixXd v(5,6);
	Eigen::MatrixXd res;

	x.setRandom();
	v.setRandom();

	matlabautogen::tcdf(x, v, res);

	std::cout<<x<<std::endl;
	std::cout<<v<<std::endl;
	std::cout<<res<<std::endl;


	return 0;
}
