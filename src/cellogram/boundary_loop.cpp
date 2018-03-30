////////////////////////////////////////////////////////////////////////////////
#include <cellogram/tri2hex.h>
#include <cellogram/navigation.h>
#include <igl/is_border_vertex.h>
#include <igl/edges.h>
#include <Eigen/dense>
#include <vector>
#include <iostream>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

// -----------------------------------------------------------------------------

int find_next_starting_point(const std::vector<bool> & needs_visit ) {

	for (unsigned int i = 0; i < needs_visit.size(); i++) {
		if (needs_visit[i] ) return i;
	}
	return -1;
}

typedef std::vector<int> Path;

int get_next(
	int starting_node, 
	const std::vector<std::vector<int>> &adj, 
	std::vector<int> &pos_in_path, 
	const std::vector<bool> &needs_visit
) 
{
	std::cout << "who is next from " << starting_node << "?\n";

	int time_now = pos_in_path.back();
	for (int neigh : adj[starting_node] ) {
		std::cout << "maybe " << neigh << "?\n";
		if (needs_visit[neigh]) std::cout << " he needs a visit!\n";
		else {
			std::cout << " he does not needs a visit!\n";
			std::cout << " but " << pos_in_path[neigh] << "<" << (time_now - 1) << "\n";
		}
		if (needs_visit[neigh] || (pos_in_path[neigh] < time_now-1)) return neigh;
	}
	assert(false);
	std::cout << "NO NEXT MOVE FOUND ????\n\n\n";
	return 0;
}

void trace_loop(int starting_pos, std::vector<Path> &paths, std::vector<bool> &needs_visit, const std::vector<std::vector<int>> &adj) {
	
	std::cout << "\n\nstarting loop from: " << starting_pos << "\n";

	int n = (int)adj.size();
	const int BIG = 1000000000;
	std::vector<int> pos_in_path(n, +BIG);

	Path p; 
	int i = starting_pos;
	while (true) {
		// visit node i
		needs_visit[i] = false;
		pos_in_path[i] = (int)p.size();
		p.push_back(i);		
		
		std::cout << "visiting node: " << i << "\n";
		i = get_next( i, adj, pos_in_path , needs_visit );
		if (i == starting_pos) break;

		if (pos_in_path[i] != BIG) {
			// pacman thing
			Path subloop( p.begin() + pos_in_path[i], p.end() );

			p.resize(pos_in_path[i]);
			paths.push_back(subloop);
		}
	}
	paths.push_back(p);
}

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
					std::cout << v0 << std::endl;
				}
				//std::cout << path.size() << std::endl;
				paths.push_back(path);
			}
		}
	}

	//while (true) {
	//	int starting_pos = find_next_starting_point(border_vertex);
	//	if (starting_pos == -1) break;

	//	trace_loop(starting_pos, paths, border_vertex, adj); // adds at least one (maybe more) Paths to paths
	//}
	
	longest_path = path_to_vecXi( get_longest_path(paths) );
	
}

// -----------------------------------------------------------------------------

} // namespace cellogram
