#ifndef SPATIAL_GRID_H
#define SPATIAL_GRID_H

#include<vector>
#include<string>
#include<set>
#include"vec2.h"

namespace cellogram
{
namespace PointsUntangler
{

struct SpatialIndexIterator{

    vec2 target; // in grid pos
    std::vector< std::pair<scalar,int> > candidates;
    int k;
    float wallDistance[4];
    int minx,miny,maxx,maxy;
};

struct SpatialIndex{
    void init(const std::vector<vec2> &vert, const std::vector<bool> &exceptThese);
    void init(const std::vector<vec2> &vert );

    void setTarget(vec2, SpatialIndexIterator& ite) const;
    int nextClosest(SpatialIndexIterator& ite) const;

    int closestTo(vec2) const;

    int midPoint() const;

private:

    std::vector< std::pair<vec2,int> >  v;

    std::vector< std::vector<int> > lattice;
    int resx,resy;

    vec2 boxmin,boxmax;
    float oneOverCellSize;
    void updateBox();
    void populate();
    void findCellOf(vec2 p, int& cx, int &cy) const;
    vec2 toGridCoords(vec2) const;
    void addLatticeCell( int cx,int cy, SpatialIndexIterator& ite ) const;
    void addLatticeCol(int x, SpatialIndexIterator& ite) const;
    void addLatticeRaw(int y, SpatialIndexIterator& ite) const;
    void debugCheckBounds(int cx, int cy) const;
    void updateWallDistances( SpatialIndexIterator& ite ) const;

};

}
}
#endif // SPATIAL_GRID_H
