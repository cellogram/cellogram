#include "cellogram/image_reader.h"

#include <Eigen/Dense>
#include <iostream>

int main(int argc, char** argv)
{
	Eigen::MatrixXd img;
	cellogram::read_image(".../detection_test.tif", img);

	std::cout<<img<<std::endl;

	return 0;
}
