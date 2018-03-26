////////////////////////////////////////////////////////////////////////////////
#include <vector>
#include <Eigen/Dense>
#include <iostream>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

// -----------------------------------------------------------------------------

	void tri2hex(const Eigen::MatrixXi &F, Eigen::MatrixXi &Graph) {
		int n = F.maxCoeff();
		std::vector<std::vector<int>> adj(n);
		for (int i = 0; i < F.rows(); ++i) {
			for (int j = 0; j < 3; ++j) {
				int v1 = F(i, j);
				int v2 = F(i, (j + 1) % 3);
				adj[v1].push_back(v2);
				adj[v2].push_back(v1);
			}
		}
		int maxValence = 0;
		int dim;
		for (int i = 0; i < n; i++)
		{
			std::sort(adj[i].begin(), adj[i].end());
			auto it = std::unique(adj[i].begin(), adj[i].end());
			dim = std::distance(adj[i].begin(), it);
			adj[i].resize(std::distance(adj[i].begin(), it));
			if (dim > maxValence)
			{
				maxValence = dim;
			}
		}
		std::cout << maxValence;
		//Graph = Eigen::MatrixXi::Zero(n, maxValence);
		for (size_t i = 0; i < n; i++)
		{
			int k = 0;
			for (auto it = adj[i].begin(); it != adj[i].end(); ++it) {
				Graph(i, k) = *it;
				k++;
			}
		}
	}
// -----------------------------------------------------------------------------

} // namespace cellogram
