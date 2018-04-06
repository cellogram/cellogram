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

	/*void region_bounding(const Eigen::MatrixXi &F, const Eigen::VectorXi &region, std::vector<Region> &regions)
	{
		int nRegions = region.maxCoeff();
		
		// create a 2d vector containing all the vertices belonging to one region
		// as well as one containing the triangles
		std::vector<std::vector<int>> region_points(nRegions);
		std::vector<std::vector<int>> region_F(nRegions);
		regions.resize(nRegions);

		for (int i = 0; i <region.size(); i++)
		{
			if (region(i) > 0) {
				region_points[region(i)-1].push_back(i);
				// for this vertex find the all the faces
				for (int f = 0; f < F.rows(); ++f) {
					for (int lv = 0; lv < F.cols(); ++lv) {
						if (i == F(f,lv))
						{
							region_F[region(i) - 1].push_back(f);
						}
					}
				}
			}
		}

		// per region find the bounding vertices
		for (int i = 0; i < nRegions; i++)
		{
			// for each region the belonging triangles need to be extracted
			std::sort(region_F[i].begin(), region_F[i].end());
			region_F[i].erase(unique(region_F[i].begin(), region_F[i].end()), region_F[i].end());

			Eigen::MatrixXi F2(region_F[i].size(), 3);
			for (int j = 0; j < region_F[i].size(); j++)
			{
				F2.row(j) = F.row(region_F[i][j]);
			}
			//Eigen::VectorXi L;
			//std::vector<std::vector<int>> Ls;
			// boundary_loop(F2, L);
			Region &current_region = regions[i];
			boundary_loop(F2, current_region.region_boundary);
			//igl::boundary_loop(F2, Ls);

			//igl::opengl::glfw::Viewer v;
			//v.data().set_mesh(points, F2);
			/*for (auto L : Ls) {
				Eigen::MatrixXd pts(L.size(), points.cols());
				for (int i = 0; i < L.size(); ++i) {
					pts.row(i) = points.row(L[i]);
				}
				
				v.data().add_points(pts, Eigen::RowVector3d::Random());
			}*
			//Eigen::MatrixXd pts(L.size(), points.cols());
			//igl::slice(points, L, 1, pts);
			//v.data().add_points(pts, Eigen::RowVector3d::Random());
			//v.launch();

			// Checked with MATLAB - F2 is correct regarding the not closing of the boundary
			//std::cout << "\nTriangles\n" << F2.transpose();
			//std::cout << "\nRegion\n" << L.transpose();

			// save boundary to vector
			//std::cout << std::endl << i << std::endl;

			//Region &current_region = regions[i];
			//current_region.region_boundary = L;
			/*for (int j = 0; j < L.rows(); j++)
			{
				region_edges[i].push_back(L(j));
				//std::cout << region_edges[i][j] << ",";
			}*
		}
		//std::cout << "\nPoints\n" << points.transpose();
	}*/

// -----------------------------------------------------------------------------
} // namespace cellogram
