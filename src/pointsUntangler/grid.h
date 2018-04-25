#ifndef GRID_H
#define GRID_H


#include<vector>
#include<string>
#include<set>
#include"vec2.h"
#include <Eigen/Dense>


namespace cellogram
{
namespace PointsUntangler
{

class Grid;
class Mesh;
void meshToGrid(Mesh &m, Grid &g);

class Grid{
public:
    Grid();
    std::vector<vec2> vert;
    std::vector<int> posInGrid;

    std::vector<int> grid;
    std::vector<bool> madeUpVert;

    std::vector<int> vdesired; // this vert desires this grid position

    int sx,sy;
    void create(int sx, int sy);

    void initIndicesOnGrid(int nx,int ny);
    void initVertOnGrid(int nx,int ny);

    int assignUnassignedNiceWay();
    int assignUnassignedHardWay();

    void printf() const;

    void updatePosInGrid();

    scalar computeAvgEdgeLen() const;

    scalar energyTotal() const;

    void exportEigen(Eigen::MatrixXi &tris, Eigen::MatrixXd &newPoints) const;
    bool importXYZ(const std::string& filename );
    bool exportOBJ(const std::string& filename ) const;
    bool exportPLY(const std::string& filename ) const;
    bool exportPLYtartan(const std::string& filename ) const;

    void createVertices(int nv);

    int tryAllSwaps();
    int tryAllSwapsBordersIncluded();

    int fillGapsMakingPtsUp();

friend void meshToGrid(Mesh &m, Grid &g);

private:

    void sanityCheck( int x, int y );
    void sanityCheck();


    int tryAllBiSwaps();
    int tryAllTriSwaps();
    int tryAllQuadriSwaps();
    int tryAllSwapsAround(int gi);

    bool fixUnassignedVertexNiceWay(int vj);
    bool fixUnassignedVertexDijkstra(int vj);

    std::vector<int> distFromBorder;
    void updateDistFromBorder();

    scalar edgeLen = 1.0;

    void trimBorders();

    void assign(int gi, int vi);
    void unassign(int gi);
    int indexOf(int x,int y) const {return x+y*sx;}

    int shiftPos(int gi, int dir) const;

    vec2 baryAround(int gi) const;
    vec2 baryAroundOfExisting(int gi) const;
    vec2 avgPos(int gi) const; // only of assigned

    bool isBoundary(int gi) const;
    void swapTwo(int gi, int gj);

    //scalar energyDeltaForSwap(int gi, int gj); // positive if swap profitable
    scalar energyAround(int gi) const;
    scalar energyAround(int gi, int gj) const;
    scalar energyBetween(int grid, int vj) const;
    scalar energyAroundIf(int gi, int vi) const; // if... grid[gi] is assigned to vi

    scalar energyAroundExcept1(int gi) const;

    scalar energyAround2(int gi) const;
    scalar energyAround2(int gi, int gj) const;
    scalar energyBetween2(int grid, int vj) const;

    int friendsAround(int gi) const;

    bool testAndDoSwap( int gi, int gj);
    bool testAndDoTriSwap( int gi, int gj, int gk);
    bool testAndDoQuadriSwap( int gi, int gj, int gk, int gh);
    bool testAndDoSwapBordersIncluded( int gi, int gj);

    void clear();

    int numberOfAdj(int gi) const;
    void enlargeToInclude(int gi, int buffer);
    void enlargeGrid(int dxMin, int dxMax, int dyMin, int dyMax);
    int hopDistance(int gi, int gj) const;

    //int bestPositionFor(int vj, std::vector<int> except = {});

    int safeGiMin, safeGiMax; // first and last grid elem with a neigh
    int safeGiMinS2, safeGiMaxS2; // first and last grid elem with a 2 star
    int safeGiMinS3, safeGiMaxS3; // first and last grid elem with a 3 star
    int neigh[6];

    scalar triangleQuality(int vi, int vj, int vk ) const;

};

}
}
#endif // GRID_H
