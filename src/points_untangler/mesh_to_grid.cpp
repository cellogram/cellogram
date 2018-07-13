#include <vector>
#include "grid.h"
#include "mesh.h"
#include "my_assert.h"

namespace cellogram
{
	namespace PointsUntangler
	{

		extern Mesh m;


		struct ExpansionMove {

			int fi; // starting mesh-face
			int ei; // mesh-edge to be traversed (orthogonally, during expansion)

			int pos; // pos on the grid (linearized index)
			int dir; // dir on the grid (0..5)

			scalar score; // priority!

			int time; // for cosmetics only, NO USE

			bool operator < (const ExpansionMove&b) const {
				if (score<b.score) return true;
				if (score>b.score) return false;
				if (ei<b.ei) return true;
				if (ei>b.ei) return false;
				if (pos<b.pos) return true;
				if (pos>b.pos) return false;
				return (dir>b.dir);
			}
		};

		std::vector<ExpansionMove> moves;

		void heapUp(int i) {
			if (i >= 1) {
				int j = i / 2;
				if (moves[j]<moves[i]) {
					std::swap(moves[j], moves[i]);
					heapUp(i);
				}
			}
		}

		void heapDown(int i) {
			int j = i * 2;
			if (j >= (int)moves.size()) return;
			int k = j + 1;
			if (k<(int)moves.size()) {
				if (moves[j]<moves[k]) j = k;
			}
			if (moves[i]<moves[j]) {
				std::swap(moves[i], moves[j]);
				heapDown(j);
			}
		}

		void heapPush(const ExpansionMove &m) {
			moves.push_back(m);
			heapUp(moves.size() - 1);
		}

		void heapPop() {
			moves[0] = moves.back();

			moves.pop_back();
			heapDown(0);
		}

		ExpansionMove heapTop() {
			return moves[0];
		}

		/* FLOOD FILL MODES:
		* 0: follow mesh, build new grid, flip-edges-on-mesh
		* 1: follow mesh + existing grid, flip-edges-on-mesh
		* 2: follow mesh, build new grid, don't change mesh, travel on "fixed" faces only (for demo purposes)
		*/
		void floodFill(Mesh &m, Grid &g, int floodfillMode) {

			std::vector<bool> visited(m.E.size() * 2, false);

			/*
			auto assessFaceRegularity = [&m,&g](int fi)->float{
			const Face& f (m.F[fi]);
			return m.V[ f[0] ].distToIrr + m.V[ f[1] ].distToIrr + m.V[ f[2] ].distToIrr
			+ m.goodTriangle(f[0],f[1],f[2]);
			};*/

			int latestTime = 0;

			auto addMove = [&m, &g, floodfillMode](int ei, int fi, int pos, int dir, int time) {
				ExpansionMove move;
				const Edge& e(m.E[ei]);
				int fj = e.fi[0]; // dest face
				if (fi == fi) fj = e.fi[1];

				if (fj == -1) return;
				const Face& Fj(m.F[fj]);
				const Face& Fi(m.F[fi]);

				switch (floodfillMode) {
				case 0: move.score = -m.parallelogramError(ei); break;
				case 1: move.score = -m.parallelogramError(ei); break;
				case 2: move.score = (Fi.fixed) ? 999 : 0; break; // conquer fixed faces first
																  // other potential scores:
																  //Fi.regularity
																  //assessFaceRegularity(fi) + assessFaceRegularity(fj)
				}


				//int vi = m.F[fb].oppositeVertOfEdge( m.E[ei] );
				//move.score -= m.V[vi].disputed*100;

				//std::cout<<m.F[fa].regularity <<"+"<< m.F[fb].regularity<<"\n";
				move.ei = ei;
				move.fi = fi;
				move.dir = dir;
				move.pos = pos;

				move.time = time;
				//static int tmp = 1000; move.score += (tmp--)/1000.0;

				heapPush(move);
			};


			auto fillGrid = [&m, &g](int gi, int vi) {
				//if (m.V[vi].val!=6) {std::cout<<"WARNING: "<<vi<<" was irregular; D:"<<m.V[vi].distToIrr<<"\n"; return;}
				//if (m.V[vi].dontcare) {std::cout<<"WARNING: dint' care bout "<<vi<<" D:"<<m.V[vi].distToIrr<<"\n"; return; }
				g.assign(gi, vi);
				m.V[vi].timeReached = 0.0;

			};

			int a = 90000;

			auto playMoveFix = [&m, &g, addMove, &visited, &latestTime, &a, floodfillMode](ExpansionMove move) {
				//std::cout<<"play move\n";

				if (floodfillMode == 2) if (move.score<800) return;

				int f0 = m.E[move.ei].fi[0];
				int f1 = m.E[move.ei].fi[1];
				int kk = 0;
				if (move.fi == f0) { std::swap(f0, f1); kk = 1; }

				// expanding from fj to fi

				/*scalar avgError = ( m.V[ m.F[f1].vi[0] ].discrepancy
				+ m.V[ m.F[f1].vi[1] ].discrepancy
				+ m.V[ m.F[f1].vi[2] ].discrepancy) / 3;*/

				//std::cout<<"On to face "<<fi<<"\n";
				const Face &f(m.F[f0]);

				int w0 = f.cornerOfEdge(move.ei);
				int w1 = (w0 + 1) % 3;
				int w2 = (w0 + 2) % 3;

				int gi = g.shiftPos(move.pos, move.dir);
				int vi = f[w2];

				int gj = g.posInGrid[vi];
				int vj = g.grid[gi];

				myAssert((gj == gi) == (vj == vi), "not reciprocal??");

				if ((vj != -1) && (vj != vi)) {
					m.V[vi].disputed++;
					m.V[vj].disputed++;
					//std::cout<<"A dispute! Vert already allocated somewhere else\n";
				}
				else if ((gj != -1) && (gj != gi)) {
					m.V[vi].disputed++;
					//std::cout<<"A dispute! Square already occupied \n";
				}
				if (g.vdesired[vi] == -1) g.vdesired[vi] = gi;

				if (visited[move.ei * 2 + kk]) return; // been there!

				visited[move.ei * 2 + kk] = true;

				if ((gj == -1) && (vj == -1)) {
					g.assign(gi, vi);
					m.V[vi].timeReached = move.time;
					latestTime = std::max(latestTime, move.time);
				}

				// enlist new moves
				addMove(f.ei[w1], f0, move.pos, (move.dir + 5) % 6, move.time + 1);
				addMove(f.ei[w2], f0, gi, (move.dir + 1) % 6, move.time + 1);

				//std::cout<<"placed "<<f[w0]<<","<<f[w1]<<" -> "<<f[w2]<<"\n";

			};

			auto playMove = [&m, &g, addMove, &visited, &latestTime, &a](ExpansionMove move) {
				//std::cout<<"play move\n";
				const Edge &E(m.E[move.ei]);
				int f0 = E.fi[0];
				int f1 = E.fi[1];
				if (move.fi == f0) { std::swap(f0, f1); }
				// expanding from f1 to f0

				// vertex back
				int vB = m.F[f1].oppositeVertOfEdge(move.ei);

				int w0 = m.F[f0].cornerOfEdge(move.ei);
				int w1 = (w0 + 1) % 3;
				int w2 = (w0 + 2) % 3;

				int gi = g.shiftPos(move.pos, move.dir); // grid square to conquer
				int vG = g.grid[gi]; // vertex currently in that square (if any)

									 // GO STRAIGHT?
				int vC = m.F[f0].vi[w2];
				scalar badC = m.parallelogramError(vB, E.vi[0], E.vi[1], vC);
				if ((vG != -1) && (vC != vG)) badC += 1000; // some other vert is here already
				if ((g.posInGrid[vC] != -1) && (g.posInGrid[vC] != gi)) badC += 1000; // this vert is somewhere else already

																					  // GO LEFT?
				int eL = m.F[f0].ei[w1];
				int vL = m.viOnOtherSide(f0, eL);
				scalar badL = 1e30;
				if (vL != -1) {
					badL = m.parallelogramError(vB, E.vi[0], E.vi[1], vL);
					if ((vG != -1) && (vL != vG)) badL += 1000; // some other vert is here already
					if ((g.posInGrid[vL] != -1) && (g.posInGrid[vL] != gi)) badL += 1000; // this vert is somewhere else already
				}

				// GO RIGHT?
				int eR = m.F[f0].ei[w2];
				int vR = m.viOnOtherSide(f0, eR);
				scalar badR = 1e30;
				if (vR != -1) {
					badR = m.parallelogramError(vB, E.vi[0], E.vi[1], vR);
					if ((vG != -1) && (vR != vG)) badR += 1000; // some other vert is here already
					if ((g.posInGrid[vR] != -1) && (g.posInGrid[vR] != gi)) badR += 1000; // this vert is somewhere else already
				}

				int choice = 0;
				scalar minBad = badC;

				if (m.canFlip(eL) && minBad>badL) { choice = 1; minBad = badL; }
				if (m.canFlip(eR) && minBad>badR) { choice = 2; minBad = badR; }

				if (choice != 0) {
					m.applyFlip((choice == 1) ? eL : eR);
					f0 = E.fi[0];
					f1 = E.fi[1];
					if (move.fi == f0) { std::swap(f0, f1); }

					// will redo this move
					addMove(move.ei, move.fi, move.pos, move.dir, move.time);
				}
				else {
					// conquer!
					int vi = m.F[f0].oppositeVertOfEdge(move.ei);

					if (g.vdesired[vi] == -1) g.vdesired[vi] = gi;

					if (visited[move.ei]) return; // been there!
					visited[move.ei] = true;

					if ((g.posInGrid[vi] != -1) && (g.posInGrid[vi] != gi)) return;
					if ((g.grid[gi] != -1) && (g.grid[gi] != vi)) return;
					if (g.grid[gi] == -1) m.V[vi].timeReached = move.time;
					latestTime = std::max(latestTime, move.time);
					g.assign(gi, vi);
					//std::cout<<"placed "<<f[w0]<<","<<f[w1]<<" -> "<<f[w2]<<"\n";

					// enlist new moves
					addMove(m.F[f0].ei[w1], f0, move.pos, (move.dir + 5) % 6, move.time + 1);
					addMove(m.F[f0].ei[w2], f0, gi, (move.dir + 1) % 6, move.time + 1);
				}

			};

			if (floodfillMode != 1) {
				g.clear();
				g.create(200, 200);
				g.createVertices(m.V.size());
			}
			int x = 100, y = 100;
			m.setDistanceToIrr();
			//m.setFaceRegularity();

			// copy geometry
			for (uint i = 0; i<m.V.size();i++) g.vert[i] = m.V[i].p;

			for (Vert& v : m.V) { v.disputed = 0; v.timeReached = -1; }
			//for (Vert& v:m.V) v.discrepancy = 9e99;

			int fi = m.bestFace();
			const Face &f(m.F[fi]);
			//std::cout<<"start from face "<<fi<<" (with "<<"score="<< f.regularity<<"="<<m.V[f[0]].distToIrr<<"+"<<m.V[f[1]].distToIrr<<"+"<<m.V[f[2]].distToIrr<<")\n";

			fillGrid(g.indexOf(x + 0, y + 0), f[0]);
			fillGrid(g.indexOf(x + 1, y + 0), f[1]);
			fillGrid(g.indexOf(x + 1, y + 1), f[2]);
			//m.V[ f[0] ].discrepancy = m.V[ f[1] ].discrepancy = m.V[ f[2] ].discrepancy =  0.0;

			addMove(f.ei[0], fi, g.indexOf(x + 0, y + 0), 0, 1);
			addMove(f.ei[1], fi, g.indexOf(x + 1, y + 0), 2, 1);
			addMove(f.ei[2], fi, g.indexOf(x + 1, y + 1), 4, 1);

			//g.enlargeToInclude(0,3);
			//int k=0;
			while (!moves.empty()) {
				//if (++k>18170) break;

				if (floodfillMode == 2) playMoveFix(heapTop());
				else playMove(heapTop());

				heapPop();
			}

			for (Vert& v : m.V) {
				if (v.timeReached >= 0) v.timeReached /= latestTime; else v.timeReached = 2.0;
				v.disputed /= 6.0;
			}

			g.updatePosInGrid();
			m.updateAverageEdge();
			g.edgeLen = m.avgEdge;
			//g.trimBorders();

		}

		void meshToGrid(Mesh &m, Grid &g) {
			std::cout << "MESH TO GRID ... \n";
			floodFill(m, g, 0);
		}

		void gridToMesh(Grid &g, Mesh &m) {
			std::cout << "GRID TO MESH... \n";
			floodFill(m, g, 1);
			//m.flipAs(g);
			floodFill(m, g, 1);

		}


	}
} // namespaces
