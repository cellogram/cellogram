
#include <cellogram/image_reader.h>
#include <Eigen/Dense>
#include <iostream>

int main(int argc, char** argv) {
	const std::string root = DATA_DIR;
	const std::string filename = root + "cell.tif";


	Eigen::MatrixXd img;
	cellogram::read_tif_image(filename, img);
	std::cout<<img<<std::endl;


	return 0;
}
