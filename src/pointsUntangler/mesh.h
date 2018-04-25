#pragma once

#include <vector>
#include <string>
#include <set>
#include "vec2.h"
#include <Eigen/Dense>

namespace cellogram
{
namespace PointsUntangler
{

struct Edge;

struct Face{
    int vi[3];
    int ei[3];
    int operator[](int i) const {return vi[i];}
    int& operator[](int i){return vi[i];}
    bool dontcare;

    scalar regularity; // the higer, the safest this face

    bool has(int e) const{ return (ei[0]==e) || (ei[1]==e)  || (ei[2]==e);}
    int oppositeVertOfEdge(const Edge& e) const;
};

struct Edge{
    int vi[2];
    int fi[2];
    int operator[](int i) const {return vi[i];}
    int& operator[](int i){return vi[i];}
    Edge(int a, int b){
        vi[0] = std::max(a,b);
        vi[1] = std::min(a,b);
        fi[0] = fi[1] = -1;
    }
    bool operator < (const Edge&b) const { if(vi[0]<b.vi[0]) return true; if(vi[0]>b.vi[0]) return false; return (vi[1]<b.vi[1]); }
    void substitute(int fa, int fb);
    bool has(int fa) const{ return (fi[0]==fa) || (fi[1]==fa);}
};

struct Vert{
    vec2 p;
    int val;
    bool dontcare = false;
    int distToIrr;

    float timeReached; // for visuzliazion purposes only
    scalar disputed = 0;

    scalar price() const; // how painful is having a irregular vertex here (default 1)

};


struct FlipScore{
    scalar valReduction;
    scalar lenReduction;

    FlipScore(scalar a, scalar b):valReduction(a),lenReduction(b){}

    static FlipScore zero(){ return FlipScore(0,0); } // doing this flip has no effect on energy

    bool operator < (const FlipScore&b) const
    {
        if(valReduction<b.valReduction) return true;
        if(valReduction>b.valReduction) return false;
        return lenReduction<b.lenReduction;
    }
};

class Grid;
class Mesh;
void meshToGrid(Mesh &m, Grid &g);

class Mesh
{

public:
    void fromEigen(const Eigen::MatrixXd &mat);

    bool importXYZ(const std::string& filename );
    bool importXYZv2(const std::string& filename );
    bool exportOBJ(const std::string& filename );
    bool exportOFF(const std::string& filename );
    bool exportPLY(const std::string& filename, bool colorByFloodfill = false);

    void regularizeByFlips(int howDeep, int howWide);
    void delaunay();

    void fuckUp(); // to test heuristic

    std::vector<Vert> V;
    std::vector<Edge> E;
    std::vector<Face> F;

    void smoothDisputed();

friend void meshToGrid(Mesh &m, Grid &g);

private:
    //vec2 avgEdgeDir;
    float avgEdge;

    std::set<int> irregulars;

    std::vector<Face> bestConf;

    int nVal;
    int bestVal = 1000000;

    //void updateAverageVecDir();
    void updateAverageEdge();
    void updateValencies();
    void updateIndices();
    void dontCareAboutBoundaries();
    void propagateDontcareToFaces();
    void propagateDontcareToVerts();
    bool canFlip(int ei);
    scalar deltaEng(int v0, int v1, int delta);
    scalar goodTriangle(int v0,int v1, int v2) const;
    FlipScore evaluateFlip(int ei, bool desperate);
    int bestFlip(bool force );
    bool sanityCheck();
    void updateValence(int vi, int delta);
    void applyFlip(int ei);
    void removeDontcare();
    bool checkIfBest();
    void setDistanceToIrr();
    void setFaceRegularity();
    int mostRegularFace() const;
};

}
}

