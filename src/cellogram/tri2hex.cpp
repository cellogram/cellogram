////////////////////////////////////////////////////////////////////////////////
#include <vector>
#include <Eigen/Dense>
#include <iostream>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

// -----------------------------------------------------------------------------

	void adjacency_list(const Eigen::MatrixXi &F, std::vector<std::vector<int>> &adj) {
		int n = F.maxCoeff() + 1;
		adj.resize(n);

		for (int i = 0; i < F.rows(); ++i) {
			for (int j = 0; j < 3; ++j) {
				int v1 = F(i, j);
				int v2 = F(i, (j + 1) % 3);
				adj[v1].push_back(v2);
				adj[v2].push_back(v1);
			}
		}

		for (int i = 0; i < n; i++)
		{
			std::sort(adj[i].begin(), adj[i].end());
			auto it = std::unique(adj[i].begin(), adj[i].end());
			adj[i].resize(std::distance(adj[i].begin(), it));
		}

	}

	void triangle_region_list(const Eigen::VectorXi &vertex_region_id, const Eigen::MatrixXi &F, Eigen::VectorXd &face_region_id) {
		face_region_id = Eigen::VectorXd::Zero(F.rows());
		for (size_t i = 0; i < vertex_region_id.size(); i++)
		{
			for (size_t j = 0; j < F.rows(); j++)
			{
				if (F(j, 0) == i || F(j, 1) == i || F(j, 2) == i)
				{
					face_region_id(j) = (double)vertex_region_id(i);
				}
			}
		}
	}
// -----------------------------------------------------------------------------

} // namespace cellogram
