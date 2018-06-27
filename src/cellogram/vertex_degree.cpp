////////////////////////////////////////////////////////////////////////////////
#include "vertex_degree.h"
#include <igl/edges.h>
#include <Eigen/Dense>
#include <iostream>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

// -----------------------------------------------------------------------------

void vertex_degree(const Eigen::MatrixXi &F, Eigen::VectorXi &degree)
{
	degree = Eigen::VectorXi::Zero(F.maxCoeff()+1);
	/*for (int i = 0; i < F.rows(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			degree(F(i, j))++;
			degree(F(i, (j + 1) % 3))++;
		}
	}*/

	Eigen::MatrixXi E;
	igl::edges(F, E);

	//std::cout << E.transpose() << std::endl;

	for (int i = 0; i < E.rows(); i++)
	{
		for (int j = 0; j < 2; j++)
		{
			degree(E(i, j))++;
			degree(E(i, (j + 1) % 2))++;
		}

	}
	degree /= 2;
}

// -----------------------------------------------------------------------------

} // namespace cellogram
