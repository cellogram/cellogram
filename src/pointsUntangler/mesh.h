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
    bool fixed;

    int operator[](int i) const {return vi[i];}
    int& operator[](int i){return vi[i];}
    bool dontcare;

    scalar regularity; // the higer, the safest this face

    bool has(int e) const{ return (ei[0]==e) || (ei[1]==e)  || (ei[2]==e);}
    int oppositeVertOfEdge(const Edge& e) const;

    int cornerOfEdge( int ei ) const { for (int w=0; w<3; w++) if (this->ei[w]==ei) return(w); assert(0); return -1; }
    int cornerOfVert( int vi ) const { for (int w=0; w<3; w++) if (this->vi[w]==vi) return(w); assert(0); return -1; }

    int oppositeVertOfEdge( int ei ) const {
        for (int w=0; w<3; w++) if (this->ei[w]==ei) return(vi[(w+2)%3]); assert(0); return -1;
    }
};

struct Edge{
    int vi[2];
    int fi[2];
    bool fixed;
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
    int distToIrr; // distance to closest regular point
    //float discrepancy; // distance from parallelogram rule: how smooth went my inclusion to the grid

    float timeReached; // for visuzliazion purposes only
    scalar disputed = 0;

    scalar price() const; // how painful is having a irregular vertex here (default 1)

};


struct FlipScore{
    scalar valReduction;
    scalar lenReduction;

    FlipScore(scalar a, scalar b):valReduction(a),lenReduction(b){}
    FlipScore():valReduction(0),lenReduction(0){}

    static FlipScore zero(){ return FlipScore(0,0); } // doing this flip has no effect on energy

    bool isPos() const {
        if(valReduction>0) return true;
        return lenReduction > 0.0002;
    }

    void operator += (const FlipScore&b) {
        valReduction+=b.valReduction;
        lenReduction+=b.lenReduction;
    }

    bool operator < (const FlipScore&b) const
    {

        //return (valReduction*1+lenReduction)<(b.valReduction*1+b.lenReduction);

        if(valReduction<b.valReduction) return true;
        if(valReduction>b.valReduction) return false;
        return lenReduction<b.lenReduction;
    }
};

class Grid;
class Mesh;
void meshToGrid(Mesh &m, Grid &g);

typedef enum { BY_VAL , BY_FLOOD , BY_DISPUTED} ColMode;

class Mesh
{

public:
    void fromEigen(const Eigen::MatrixXd &mat);

    bool importXYZ(const std::string& filename );
    bool importXYZv2(const std::string& filename );
    bool importXYZv3(const std::string& filename );
    bool importFVFix(const std::string& fnV , const std::string& fnF, const std::string& fnFix);
    bool exportOBJ(const std::string& filename );
    bool exportOFF(const std::string& filename );
    bool exportPLY(const std::string& filename, ColMode colMode);

    bool exportEdgesPLY(const std::string& filename);
    bool exportEdgesPLY(const std::string& filename, const Grid& g);

    void greedyFlips(int howDeep, Grid &g);
    void flipAs(const Grid& g);
    void delaunay();

    void fuckUp(); // to test heuristic

    std::vector<Vert> V;
    std::vector<Edge> E;
    std::vector<Face> F;

    void smoothDisputed();

    friend void floodFill(Mesh &m, Grid &g, int floodfillMode);
    friend void test01();
    friend void test02();

    void swapVert(int vi, int vj);

    void buildEdgesFromFaces();
    void propagateFixedF2E();
    void updateValencies();

private:
    //vec2 avgEdgeDir;
    float avgEdge;

    std::set<int> irregulars;

    struct{
        std::vector<Face> F;
        std::vector<Edge> E;
        FlipScore score;
    } bestEver;

    int nVal;
    int bestVal = 1000000;

    //void updateAverageVecDir();
    void updateAverageEdge();
    
    void sanityCheckValencies();

    void dontCareAboutBoundaries();
    void propagateDontcareV2F();
    void propagateDontcareF2V();

    bool canFlip(int ei);
    scalar deltaEng(int v0, int v1, int delta);
    scalar goodTriangle(int v0,int v1, int v2) const;
    FlipScore evaluateFlip(int ei, Grid &g);
    int bestFlip(FlipScore &score, Grid &g);
    bool sanityCheck();
    void updateValence(int vi, int delta);
    void applyFlip(int ei);
    void removeDontcare();
    void setDistanceToIrr();
    void setFaceRegularity();
    int mostRegularFace() const;

    vec2 parellelogramRule( int fi, int ei ) const;
    scalar parallelogramError( int ei ) const;
    scalar parallelogramError( int vB,int vL,int vR,int vF ) const;

    int viOnOtherSide(int fa, int ei) const;

};

}
}

