#include <vector>
#include "grid.h"
#include "mesh.h"
#include "my_assert.h"

namespace cellogram
{
namespace PointsUntangler
{

struct ExpansionMove{
    int ei; // edge to be traversed
    int fi; // starting face
    int dir; // on the grid
    scalar score;
    int pos; // grid index
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

void heapUp(int i){
    if (i>=1) {
        int j = i/2;
        if (moves[j]<moves[i]) {
            std::swap(moves[j],moves[i]);
            heapUp(i);
        }
    }
}

void heapDown(int i){
    int j = i*2;
    if (j>=(int)moves.size()) return;
    int k = j++;
    if (k<(int)moves.size()) {
        if (moves[j]<moves[k]) j = k;
    }
    if (moves[i]<moves[j]) {
        std::swap(moves[i],moves[j]);
        heapDown(j);
    }
}

void heapPush(const ExpansionMove &m){
    moves.push_back(m);
    heapUp( moves.size()-1);
}

void heapPop(){
    moves[0] = moves.back();

    moves.pop_back();
    heapDown(0);
}

ExpansionMove heapTop(){
    return moves[0];
}



void meshToGrid(Mesh &m, Grid &g){

    std::vector<bool> visited(m.E.size(),false);

    int time = 0;

    auto addMove = [&m,&g](int ei, int fi, int pos, int dir){
        ExpansionMove move;
        int fa = m.E[ei].fi[0];
        int fb = m.E[ei].fi[1];
        if (fa==-1) return;
        if (fb==-1) return;
        move.score = m.F[fi].regularity; // + m.F[fb].regularity * 10;

        //int vi = m.F[fb].oppositeVertOfEdge( m.E[ei] );
        //move.score -= m.V[vi].disputed*100;

        //std::cout<<m.F[fa].regularity <<"+"<< m.F[fb].regularity<<"\n";
        move.ei = ei;
        move.fi = fi;
        move.dir = dir;
        move.pos = pos;

        //static int tmp = 1000; move.score += (tmp--)/1000.0;

        heapPush(move);
    };


    auto fillGrid = [&m,&g](int gi, int vi) {
        //if (m.V[vi].val!=6) {std::cout<<"WARNING: "<<vi<<" was irregular; D:"<<m.V[vi].distToIrr<<"\n"; return;}
        //if (m.V[vi].dontcare) {std::cout<<"WARNING: dint' care bout "<<vi<<" D:"<<m.V[vi].distToIrr<<"\n"; return; }
        g.assign(gi,vi);
    };

    auto playMove = [&m,&g,addMove,&visited,&time](ExpansionMove move){
        //std::cout<<"play move\n";
        int fi = m.E[move.ei].fi[0];
        if (move.fi==fi) fi = m.E[move.ei].fi[1];

        if (visited[move.ei]) {
            //std::cout<<"been there\n";
            return;
        }
        //std::cout<<"On to face "<<fi<<"\n";
        const Face &f(m.F[fi]);
        int w0=-1;
        for (int ww=0; ww<3; ww++) if (f.ei[ww]==move.ei) w0 = ww;
        myAssert(w0!=-1,"wrong move!");
        int w1 = (w0+1)%3;
        int w2 = (w0+2)%3;

        visited[move.ei] = true;

        int gi = g.shiftPos(move.pos, move.dir);
        int vi = f[w2];


        int gj = g.posInGrid[vi];
        int vj = g.grid[gi];

        myAssert((gj==gi)==(vj==vi),"not reciprocal??");
        if ((vj!=-1) && (vj!=vi)) { m.V[vi].disputed++; m.V[vj].disputed++; }
        else if ((gj!=-1) && (gj!=gi)) {
            m.V[vi].disputed++;
        }
        if (g.vdesired[ vi ]==-1) g.vdesired[ vi ] = gi;

        //if (m.V[f[w2]].dontcare) return;

        // place vertex
        // enlist new moves

        //bool ok = false;
        if ((gj==-1)&&(vj==-1)) {
            g.assign( gi, vi );
            m.V[vi].timeReached = time++;
            //ok = true;
        }
        addMove( f.ei[w1] ,fi, move.pos, (move.dir+5)%6);
        addMove( f.ei[w2] ,fi,   gi    , (move.dir+1)%6);


        //std::cout<<"placed "<<f[w0]<<","<<f[w1]<<" -> "<<f[w2]<<"\n";

    };

    g.clear();
    m.setDistanceToIrr();
    m.setFaceRegularity();

    int x=0,y=0;
    g.create(200,200);x = 100; y = 100;
    g.createVertices( m.V.size() );

    // copy geometry
    for (uint i=0; i<m.V.size();i++) g.vert[i]=m.V[i].p;

    for (Vert& v:m.V) v.disputed = 0;

    int fi = m.mostRegularFace();
    const Face &f(m.F[ fi ]);
    //std::cout<<"start from face "<<fi<<" (with "<<"score="<< f.regularity<<"="<<m.V[f[0]].distToIrr<<"+"<<m.V[f[1]].distToIrr<<"+"<<m.V[f[2]].distToIrr<<")\n";

    std::cout<<"MESH TO GRID ... \n";

    fillGrid( g.indexOf(x+0,y+0) , f[0] );
    fillGrid( g.indexOf(x+1,y+0) , f[1] );
    fillGrid( g.indexOf(x+1,y+1) , f[2] );

    addMove( f.ei[0], fi, g.indexOf(x+0,y+0), 0 );
    addMove( f.ei[1], fi, g.indexOf(x+1,y+0), 2 );
    addMove( f.ei[2], fi, g.indexOf(x+1,y+1), 4 );

    //g.enlargeToInclude(0,3);
    //int k=0;
    while (!moves.empty()) {
        //if (++k>18170) break;
        playMove( heapTop() );
        heapPop();
    }

    for (Vert& v:m.V) {
        v.timeReached /= scalar(time);
        //v.disputed = v.nextDisputed;
    }

    g.updatePosInGrid();
    m.updateAverageEdge();
    g.edgeLen = m.avgEdge;
    g.trimBorders();
}

}
}
