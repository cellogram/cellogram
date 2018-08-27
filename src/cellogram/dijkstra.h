#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <cellogram/common.h>
#include <geogram/mesh/mesh.h>
#include <Eigen/Dense>
#include <vector>
#include <queue>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

	struct edge {
		int x, y;
		double w;
		edge(void) : x(), y(), w() {}
		edge(int a, int b, double c) : x(a), y(b), w(c) {}
		bool operator< (const edge &e) const {
			return w > e.w; // Extract min-cost edges first
		}
	};

	void dijkstra(std::vector<std::vector<edge> > &graph, std::vector<int> &prev, std::vector<int> &dist, int source);

	// Build an undirected adjacency graph for the vertices of a triangle mesh
	void build_adjacency_graph(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, std::vector<std::vector<int>> &ajd);

	// Propagate a scalar from the given sources, ensuring a given grading
	void dijkstra_grading(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::VectorXd &S, double grading, const std::vector<int> &sources);
}
