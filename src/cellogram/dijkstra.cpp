////////////////////////////////////////////////////////////////////////////////
#include "dijkstra.h"
#include <igl/edges.h>
#include <vector>
#include <queue>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

// -----------------------------------------------------------------------------

// Dijkstra. One-to-all shortest path in positive weighted graph. O(m*log(n))
void dijkstra(std::vector<std::vector<edge> > &graph, std::vector<int> &prev, std::vector<int> &dist, int source) {
	std::priority_queue<edge> file;
	std::vector<bool> mark(graph.size(), false);
	prev.assign(graph.size(), -1);
	dist.assign(graph.size(), INT_MAX);
	file.push(edge(source, source, 0));
	while (not file.empty()) {
		edge e = file.top(); file.pop();
		if (mark[e.y]) continue;
		prev[e.y] = e.x;
		dist[e.y] = e.w; // Distance form source is now correct
		mark[e.y] = true;
		for (int i = 0; i < (int) graph[e.y].size(); i++) {
			edge f = graph[e.y][i];
			f.w += e.w;
			if (dist[f.y] > f.w) file.push(f);
		}
	}
}

// -----------------------------------------------------------------------------

std::vector<std::vector<int>> adjacency_graph(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) {
	Eigen::MatrixXi E;
	igl::edges(F, E);
	std::vector<std::vector<int>> adj(V.rows());
	for (int e = 0; e < E.rows(); ++e) {
		int x = E(e, 0);
		int y = E(e, 1);
		adj[x].push_back(y);
		adj[y].push_back(x);
	}
	return adj;
}

// -----------------------------------------------------------------------------

void dijkstra_grading(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::VectorXd &S, double grading, const std::vector<int> &sources) {
	typedef std::pair<int, int> Segment;
	typedef std::pair<double, Segment> WeightedSegment;
	std::priority_queue<WeightedSegment> q;
	size_t n = V.rows();
	std::vector<bool> marked(n, false);
	for (int x : sources) {
		q.emplace(0, Segment(x, x));
	}
	auto adj = adjacency_graph(V, F);
	while (!q.empty()) {
		auto kv = q.top(); q.pop();
		int x = kv.second.first;
		int y = kv.second.second;
		double dist = -kv.first;
		if (marked[y]) { continue; }
		marked[y] = true;
		S[y] = S[x] + dist * grading;
		for (int z : adj[y]) {
			dist += (V.row(y) - V.row(z)).norm();
			q.emplace(-dist, Segment(x, z));
		}
	}
}

} // namespace cellogram
