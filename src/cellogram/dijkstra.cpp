////////////////////////////////////////////////////////////////////////////////
#include "dijkstra.h"
#include <vector>
#include <queue>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

// -----------------------------------------------------------------------------

/********************************************************************************
* Dijkstra. One-to-all shortest path in positive weighted graph. O(m*log(n))
********************************************************************************/



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
			for (int i = 0; i < graph[e.y].size(); i++) {
				edge f = graph[e.y][i];
				f.w += e.w;
				if (dist[f.y] > f.w) file.push(f);
			}
		}
	}

// -----------------------------------------------------------------------------

} // namespace cellogram
