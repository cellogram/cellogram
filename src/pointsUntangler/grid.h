#ifndef GRID_H
#define GRID_H


#include<vector>
#include<string>
#include<set>
#include"vec2.h"
#include"spatial_index.h"
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
    std::vector<bool> fixed;
    std::vector<bool> madeUpVert;


    std::vector<int> vdesired; // this vert desires this grid position


    int sx,sy;
    void create(int sx, int sy);

    void initIndicesOnGrid(int nx,int ny);
    void initVertOnGrid(int nx,int ny);

    //void exportToEigen( Eigen::MatrixXi &t , Eigen::MatrixXd &newv );
    //void inputFromEigen( const Eigen::MatrixXd &v );

    int assignUnassignedNiceWay();
    int assignUnassignedHardWay();

    void setBoundaryAsFixed();
    void setNullAsFixed();

    void printf() const;
    void printf( const std::vector<int> &diff ) const;

    void shuffle(int howMany);

    void heuristicConstruct();
    void heuristic();
    void heuristicSimple();
    void sanityCheck( int x, int y );

    bool updateHappiness();
    void updatePosInGrid();

    scalar computeAvgEdgeLen() const;

    scalar energyTotal() const;

    void exportEigen(Eigen::MatrixXi &tris, Eigen::MatrixXd &newPoints) const;
    bool importXYZ(const std::string& filename );
    bool exportOBJ(const std::string& filename ) const;
    bool exportPLY(const std::string& filename ) const;
    bool exportPLYtartan(const std::string& filename ) const;

    void initSpatialIndex();

    void createVertices(int nv);

    int tryAllSwaps();
    int tryAllSwapsBordersIncluded();

    int fillGapsMakingPtsUp();

friend void meshToGrid(Mesh &m, Grid &g);

private:

    int tryAllBiSwaps();
    int tryAllTriSwaps();
    int tryAllQuadriSwaps();
    int tryAllSwapsAround(int gi);

    bool fixUnassignedVertexNiceWay(int vj);
    bool fixUnassignedVertexHardWay(int vj);
    bool fixUnassignedVertexDijkstra(int vj);

    std::vector<int> distFromBorder;
    void updateDistFromBorder();

    scalar edgeLen = 1.0;


    void assign(int gi, int vi);
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

    SpatialIndex spatialIndex;

    // heuristic stuff
    bool findABetterSuitorFor(int gi , scalar budget, int depth , std::vector<int> chain);
    int nSwapsAttempted; // bookkeeping
    int maxDepth;

    void clear();

    std::vector<bool> isHappy;

    int randomNonFixedPos() const;

    // for construction heuristic
    std::set<int> boundary;
    std::set<int> suitors;
    void conquer(int gi, int vi);
    int numberOfAdj(int gi) const;
    void enlargeToInclude(int gi, int buffer);
    void enlargeGrid(int dxMin, int dxMax, int dyMin, int dyMax);
    int closestPosToInBoundary(std::vector<int> p) const;
    int closestPosTo(std::vector<int> p, std::vector<int> except) const;
    int hopDistance(int gi, int gj) const;


    bool fillOneBoundary();
    bool suitOneSuitor();

    int bestPositionFor(int vj, std::vector<int> except = {});
    int bestApproachingPositionFor(int vj) const;

    int safeGiMin, safeGiMax; // first and last grid elem with a neigh
    int safeGiMinS2, safeGiMaxS2; // first and last grid elem with a 2 star
    int safeGiMinS3, safeGiMaxS3; // first and last grid elem with a 3 star
    int neigh[6];
    int nDone;
};

}
}
#endif // GRID_H
