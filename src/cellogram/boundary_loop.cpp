////////////////////////////////////////////////////////////////////////////////
#include "boundary_loop.h"
#include <cellogram/tri2hex.h>
#include <cellogram/navigation.h>
#include <igl/is_border_vertex.h>
#include <igl/edges.h>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

// -----------------------------------------------------------------------------
	typedef std::vector<int> Path;

Path get_longest_path(const std::vector<Path> &paths) {
	if (paths.empty()) return Path();
	int maxi = 0;
	for (unsigned int i = 1; i < paths.size(); i++) {
		if (paths[maxi].size() < paths[i].size()) {
			maxi = i;
		}
	}
	return paths[maxi];
}

Eigen::VectorXi path_to_vecXi(const Path &p) {
	Eigen::VectorXi res(p.size());
	for (size_t i = 0; i < p.size(); i++) res(i)=p[i];
	return res;
}

void boundary_graph(const Eigen::MatrixXi &F, std::vector<std::vector<int>> &adj) {
	int num_vertices = F.maxCoeff() + 1;
	NavigationData data(F);
	adj.clear();
	adj.resize(num_vertices);
	for (int f = 0; f < F.rows(); ++f) {
		NavigationIndex idx = index_from_face(F, data, f, 0);
		for (int lv = 0; lv < F.cols(); ++lv) {
			if (switch_face(data, idx).face < 0) {
				int v1 = idx.vertex;
				int v2 = switch_vertex(data, idx).vertex;
				adj[v1].push_back(v2);
				adj[v2].push_back(v1);
			}
			idx = next_around_face(data, idx);
		}
	}
}

void boundary_loop(const Eigen::MatrixXi &F, Eigen::VectorXi &longest_path) {

	int n = F.maxCoeff() + 1;

	Eigen::VectorXd Vdummy(F.maxCoeff() + 1, 1);
	std::vector<bool> border_vertex = igl::is_border_vertex(Vdummy, F);

	std::vector<std::vector<int>> adj;
	//tri2hex(F, adj);
	boundary_graph(F, adj);

	std::vector< Path > paths; // paths that we detect

	int num_vertices = adj.size();
	std::vector<bool> possible_starting_point(num_vertices, true);
	for (int v0 = 0; v0 < num_vertices; ++v0) {
		if (adj[v0].size() == 2 && possible_starting_point[v0]) {
			std::vector<bool> in_the_path(num_vertices, false);
			std::vector<bool> has_been_removed(num_vertices, false);
			std::vector<int> prev(num_vertices, -1);
			int x = v0;
			Path path;
			path.push_back(v0);
			in_the_path[v0] = true;
			prev[v0] = v0;
			for (int i = 0; i < num_vertices; ++i) {
				// Check neighbors of x
				bool found_v0 = false;
				for (int y : adj[x]) {
					if (in_the_path[y]) {
						if (y == v0) {
							found_v0 = true;
						} else if (prev[x] != y) {
							int z = x;
							while (z != y) {
								path.pop_back();
								has_been_removed[z] = true;
								in_the_path[z] = false;
								int old_z = z;
								z = prev[z];
								prev[old_z] = -1;
							}
							has_been_removed[y] = false;
							x = y;
							found_v0 = false;
							break;
						}
					} else if (!has_been_removed[y]) {
						prev[y] = x;
						x = y;
						found_v0 = false;
						path.push_back(x);
						in_the_path[x] = true;
						break;
					}
				}

				if (x == v0 || found_v0) {
					break;
				}
			}

			// Cleanup
			for (int x : path) {
				possible_starting_point[x] = false;
			}
			if (!path.empty()) {
				if (path.size() == 3) {
					//std::cout << "boundary_loop.cpp, path size is only three " << v0 << std::endl << "Path not added\n";
				}
				else{
					//std::cout << path.size() << std::endl;
					paths.push_back(path);
				}
			}
		}
	}


	longest_path = path_to_vecXi( get_longest_path(paths) );

}

// -----------------------------------------------------------------------------

} // namespace cellogram
