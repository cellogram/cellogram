////////////////////////////////////////////////////////////////////////////////
#include <vector>
#include <Eigen/Dense>
#include <iostream>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

// -----------------------------------------------------------------------------
	void check_crit(const Eigen::MatrixXi &Graph, const Eigen::VectorXi &Vertex_Value, Eigen::VectorXi &region, Eigen::VectorXi &visited, const double criterium, int group, int ind)
	{
		if (visited(ind) == 0 && Vertex_Value(ind) > criterium)
		{
			// This vertex has now been visited and is assigned to the current group
			visited(ind) = 1;
			region(ind) = group;

			// Start a new loop going through all the vertex connected directly to this one
			int j = 0;
			while (Graph(ind,j)!=0)
			{
				check_crit(Graph, Vertex_Value, region, visited, criterium, group, Graph(ind, j));
				j++;
			}
			
		}
	}


	void region_grow(const Eigen::MatrixXi &Graph, const Eigen::VectorXi &Vertex_Value, const double criterium, Eigen::VectorXi &region)
	{
		int n = Graph.rows();
		int group = 1;
		Eigen::VectorXi visited = Eigen::VectorXi::Zero(n);
		
		int k = 0;
		while (k < n)
		{
			// check if point has been visited and if criterium is met
			check_crit(Graph, Vertex_Value, region, visited, criterium, group, k);

			// when function returns all vertices indirectly connected to this one are assigned to same group
			group++;
			k++;
		}
	}
// -----------------------------------------------------------------------------

} // namespace cellogram
