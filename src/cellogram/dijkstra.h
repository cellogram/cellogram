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
}
