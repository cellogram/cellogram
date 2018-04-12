////////////////////////////////////////////////////////////////////////////////
#include "region_grow.h"

#include <queue>
#include <Eigen/Dense>
#include <iostream>
#include <algorithm>
#include "navigation.h"
#include <igl/slice.h>
#include <cellogram/boundary_loop.h>
#include <igl/opengl/glfw/Viewer.h>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

// -----------------------------------------------------------------------------
	void region_grow(std::vector<std::vector<int>> &Graph, const Eigen::Matrix<bool, 1, Eigen::Dynamic> &crit_pass, Eigen::VectorXi &region)
	{
		int n = crit_pass.size();
		int cId = 1;
		std::queue<int> Q;
		Eigen::VectorXi visited = Eigen::VectorXi::Zero(n);
		std::vector<int> id(n, 0);

		for (int i = 0; i < n; i++)
		{
 			Q.push(i);
			while (!Q.empty())
			{
				int u = Q.front();
				Q.pop();
				
				if (crit_pass(u) && visited(u) == 0)
				{
					//add connections to queue if not already visited
					for (std::vector<int>::iterator it = Graph[u].begin(); it != Graph[u].end(); ++it) {
						if (visited(*it) == 0)
						{
							Q.push(*it);
						}
					}
					id[u] = cId;
				}
				// set vertex to visited
				visited(u) = 1;

				//std::cout << i << " - " << crit_pass[i] << " - " << id[u] << " - " << u << std::endl;
			}
			// While loop left which means that the next vertex found is from a new group
			cId++;
		}
		
		// Loop through region to clean up IDs
		std::vector<int> map = id;
		std::sort(map.begin(), map.end());
		auto newEnd = std::unique(map.begin(), map.end());
		map.erase(newEnd, map.end());

		for (int i = 0; i < map.size(); i++)
		{
			std::replace(id.begin(), id.end(), map[i], i);
		}

		// Save into eigen format
		region = Eigen::VectorXi::Zero(n);
		for (int i = 0; i < id.size(); i++)
		{
			region(i) = id[i];
		}
	}


// -----------------------------------------------------------------------------
} // namespace cellogram
