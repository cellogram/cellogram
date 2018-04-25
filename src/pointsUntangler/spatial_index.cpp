#include <iostream>
#include <algorithm>
#include "spatial_index.h"

namespace cellogram
{
namespace PointsUntangler
{

void SpatialIndex::init(const std::vector<vec2> &vert){
    v.clear();
    for(uint i=0; i<vert.size(); i++) v.push_back(std::pair<vec2,int>(vert[i],i) );
    updateBox();
    lattice.clear();
    lattice.resize( resx*resy );
    populate();
}

int SpatialIndex::closestTo(vec2 p) const{
    SpatialIndexIterator ite;
    setTarget( p, ite );
    return nextClosest(ite);
}

int SpatialIndex::midPoint() const{
    return closestTo( (boxmin+boxmax)/2 );
}

void SpatialIndex::init(const std::vector<vec2> &vert,const std::vector<bool> &exceptThese){
    //myAssert(vert.size()==exceptThese.size(),"SIZE MISMATCH "<<" "<<vert.size()<<"!="<<exceptThese.size());
    v.clear();
    for(uint i=0; i<vert.size(); i++) if (!exceptThese[i]) v.push_back(std::pair<vec2,int>(vert[i],i) );

    updateBox();
    lattice.clear();
    lattice.resize( resx*resy );
    populate();
}

void SpatialIndex::updateBox(){
    boxmin=boxmax=v[0].first;
    for (const std::pair<vec2,int>& p:v){
        boxmin.x = std::min(boxmin.x,p.first.x);
        boxmin.y = std::min(boxmin.y,p.first.y);
        boxmax.x = std::max(boxmax.x,p.first.x);
        boxmax.y = std::max(boxmax.y,p.first.y);
    }
    vec2 boxsize = boxmax-boxmin;
    float safety = boxsize.x/1000.0;
    boxmin-=vec2(safety,safety);
    boxmax+=vec2(safety,safety);
    boxsize = boxmax-boxmin;

    float boxArea = (boxsize.x*boxsize.y);
    float cellArea = boxArea / v.size(); // approx
    oneOverCellSize = sqrt(1/cellArea);
    resx = 1+boxsize.x * oneOverCellSize;
    resy = 1+boxsize.y * oneOverCellSize;

    for (std::pair<vec2,int>& p:v) p.first = toGridCoords(p.first);
}

void SpatialIndex::populate(){

    for (uint i=0; i<v.size(); i++){
        int x,y;
        findCellOf(v[i].first,x,y);
        lattice[x+y*resx].push_back(i);
    }
}

void SpatialIndex::debugCheckBounds(int cx, int cy) const{
    if ((cx<0) || (cy<0) || (cx>=resx) || (cy>=resy)){
        std::cout<<"ERROR BOUND\n";
    }
}

vec2 SpatialIndex::toGridCoords(vec2 p) const{
    return (p-boxmin)*oneOverCellSize;
}

void SpatialIndex::findCellOf(vec2 p, int &cx, int &cy) const {
    cx = (int)floor(p.x);
    cy = (int)floor(p.y);
    if (cx<0) cx=0;
    if (cy<0) cy=0;
    if (cx>=resx) cx=resx-1;
    if (cy>=resy) cy=resy-1;
    //debugCheckBounds(cx,cy);
}

void SpatialIndex::addLatticeCol(int cx,SpatialIndexIterator &ite) const{
    for (int y=ite.miny; y<ite.maxy; y++) addLatticeCell(cx,y,ite);
}
void SpatialIndex::addLatticeRaw(int cy,SpatialIndexIterator &ite) const{
    for (int x=ite.minx; x<ite.maxx; x++) addLatticeCell(x,cy,ite);
}

void SpatialIndex::addLatticeCell(int cx, int cy, SpatialIndexIterator &ite) const{
    //debugCheckBounds(cx,cy);
    for (const int& i:lattice[cx+cy*resx]) {
        ite.candidates.push_back( std::pair<scalar,int>(distance(ite.target,v[i].first),v[i].second) );
        // bubble sort the new element deeper
        int k = ite.candidates.size()-1;
        while (k>0 && ite.candidates[k]<ite.candidates[k-1]) {
            std::swap(ite.candidates[k],ite.candidates[k-1]);
            k--;
        }
    }
}

// BRUTAL BUT SAFE WAY:
/*
void SpatialIndex::setTarget(vec2 t,SpatialIndexIterator& ite) const{
    ite.candidates.resize(v.size());
    for (uint i=0; i<ite.candidates.size(); i++) {
        ite.candidates[i].second = i;
        ite.candidates[i].first = squaredDistance(v[i],t);
    }
    std::sort(ite.candidates.begin(), ite.candidates.end());
    ite.k=0;
}

int SpatialIndex::nextClosest(SpatialIndexIterator& ite) const{
    if (ite.k==(int)ite.candidates.size()) return -1;
    return ite.candidates[ite.k++].second;
}
*/
static const scalar VERY_FAR = 1000000000.0;

void SpatialIndex::updateWallDistances(SpatialIndexIterator &ite) const{
    ite.wallDistance[0] = (ite.minx==0   )?VERY_FAR:(ite.target.x-ite.minx);
    ite.wallDistance[1] = (ite.maxx==resx)?VERY_FAR:(ite.maxx-ite.target.x);
    ite.wallDistance[2] = (ite.miny==0   )?VERY_FAR:(ite.target.y-ite.miny);
    ite.wallDistance[3] = (ite.maxy==resy)?VERY_FAR:(ite.maxy-ite.target.y);
}

void SpatialIndex::setTarget(vec2 t,SpatialIndexIterator& ite) const{
    ite.target = toGridCoords(t);

    int cx,cy;
    findCellOf(ite.target,cx,cy);
    ite.minx=cx; ite.maxx=cx+1;
    ite.miny=cy; ite.maxy=cy+1;

    addLatticeCell(cx,cy, ite);
    ite.k=0;
    updateWallDistances(ite);
}

int SpatialIndex::nextClosest(SpatialIndexIterator& ite) const{
    //if (!ite.candidates.size()) return -1;

    //std::cout<<"CLOSEST...\n";
    while (1) {
        int a = (ite.wallDistance[0]<ite.wallDistance[1])?0:1;
        int b = (ite.wallDistance[2]<ite.wallDistance[3])?2:3;
        int closestWallI = (ite.wallDistance[a]<ite.wallDistance[b])?a:b;
        float closestWallDist = ite.wallDistance[closestWallI];


        if (ite.k<(int)ite.candidates.size()) {
            auto closest = ite.candidates[ite.k];
            if (closest.first<closestWallDist) {
                ite.k++;
                //std::cout<<"for the "<<ite.k<<"th: tested "<<ite.candidates.size()<<" candidates on a "<<(ite.maxx-ite.minx)<<"x"<<(ite.maxy-ite.miny)<<"\n";
                return closest.second;
            }
        }

        // expand one wall
        //std::cout<<"EXPANDING WALL "<<closestWallI<<"\n";
        if (closestWallDist==VERY_FAR) return -1;
        switch(closestWallI){
        case 0: ite.minx--;addLatticeCol(ite.minx  ,ite); break;
        case 1: ite.maxx++;addLatticeCol(ite.maxx-1,ite); break;
        case 2: ite.miny--;addLatticeRaw(ite.miny  ,ite); break;
        case 3: ite.maxy++;addLatticeRaw(ite.maxy-1,ite); break;
        }
        updateWallDistances(ite);
    }
}

}
}