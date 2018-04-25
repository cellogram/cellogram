#include <fstream>
#include <iostream>
#include <iomanip>      // std::setfill, std::setw
#include <algorithm>

#include "grid.h"

#include "my_assert.h"

namespace cellogram
{
namespace PointsUntangler
{

void Grid::clear(){
    grid.clear();
    posInGrid.clear();
    vdesired.clear();
}

void Grid::createVertices(int nv){
    vert.resize(nv);
    posInGrid.resize(nv,-1);
    vdesired.resize(nv,-1);
    madeUpVert.resize(nv);
    madeUpVert.assign(nv,false);

}

bool Grid::importXYZ(const std::string& filename ){
    clear();
    std::fstream f;
    f.open(filename.data(),std::fstream::in);
    if (!f.is_open()) {
        std::cout<<"Cannot open \""<<filename<<"\"\n";
        return false;
    }
    int n;
    f>>n;
    vert.resize(n);

    for (int i=0; i<n; i++) {
        float dummy;
        f>>vert[i].x>>vert[i].y>>dummy;
        myAssert(dummy==0,"Error reading???");
    }
    f.close();

    posInGrid.resize( vert.size(),-1);
    std::cout<<"Done reading "<<n<<" verts\n";
    return true;
}

bool Grid::exportOBJ(const std::string& filename ) const{
    std::ofstream f;
    f.open(filename.data(),std::fstream::out);
    if (!f.is_open()) return false;
    for (vec2 v:vert) f<<"v "<<v.x<<" "<<v.y<<" 0\n";

    int nf=0;

    for (int y=0; y<sy-1; y++)
    for (int x=0; x<sx-1; x++) {
        int i,j,k;
        i=grid[indexOf(x,y)]; j=grid[indexOf(x,y+1)]; k=grid[indexOf(x+1,y+1)];
        if ((i!=-1) && (j!=-1) && (k!=-1)) { f<<"f "<<j+1<<" "<<i+1<<" "<<k+1<<"\n"; nf++;}
        i=grid[indexOf(x+1,y+1)]; j=grid[indexOf(x+1,y)]; k=grid[indexOf(x,y)];
        if ((i!=-1) && (j!=-1) && (k!=-1)) { f<<"f "<<j+1<<" "<<i+1<<" "<<k+1<<"\n"; nf++;}
    }
    f.close();
    std::cout<<"Done writing OBJ ("<<vert.size()<<" verts, "<<nf<<" faces)\n";
    return true;
}

void tartanColor(int gi, int sx, int &r,int &g,int &b ){
    int x = gi%sx;
    int y = gi/sx;
    int i = x%4;
    int j = y%4;
    int k = (400000+x-y)%4;
    r = 255;
    g = 255;
    b = 255;
    if (i/2%2==0) {r-=10;g-=80;b-=80;}
    if (j/2%2==0) {r-=80;g-=20;b-=80;}
    if (k/2%2==0) {r-=80;g-=70;b-=10;}

}

bool Grid::exportPLYtartan(const std::string &filename) const{

    int nUnassigned  = 0;
    for (uint vi=0; vi<posInGrid.size(); vi++)
        if ((posInGrid[vi]==-1)&& (vdesired[vi]!=-1)) nUnassigned++;

    int nv = vert.size();
    int nf = 0;
    for (int y=0; y<sy-1; y++)
    for (int x=0; x<sx-1; x++) {
        int i,j,k;
        i=grid[indexOf(x,y)]; j=grid[indexOf(x,y+1)]; k=grid[indexOf(x+1,y+1)];
        if ((i!=-1) && (j!=-1) && (k!=-1)) { nf++;}
        i=grid[indexOf(x+1,y+1)]; j=grid[indexOf(x+1,y)]; k=grid[indexOf(x,y)];
        if ((i!=-1) && (j!=-1) && (k!=-1)) { nf++;}
    }

    nf += nUnassigned*4;
    nv += nUnassigned*4;

    std::ofstream f;
    f.open(filename.data(),std::fstream::out);
    if (!f.is_open()) return false;
    f
    <<"ply\n"
    <<"format ascii 1.0\n"
    <<"element vertex "<< nv <<"\n"
    <<"property float x\n"
    <<"property float y\n"
    <<"property float z\n"
    <<"property uchar red\n"
    <<"property uchar green\n"
    <<"property uchar blue\n"
    <<"property uchar alpha\n"
    <<"element face "<< nf<<"\n"
    <<"property list uchar int vertex_indices\n"
    <<"end_header\n";

    for (uint vi=0; vi<vert.size(); vi++) {
        int r=255,g=255,b=255;
        if (posInGrid[vi]==-1) { r=g=b=0; }
        //else if (madeUpVert[i]) {r/=4; b/=4; g=200;}
        else tartanColor( posInGrid[vi], sx, r,g,b);
        f<<vert[vi].x<<" "<<vert[vi].y<<" 0 "<<r<<" "<<g<<" "<<b<<" "<<" 255\n";
    }
    for (uint vi=0; vi<vert.size(); vi++) {
        int r=255,g=255,b=255;
        if ((posInGrid[vi]==-1)&& (vdesired[vi]!=-1)) {
            g/=2; b/=2;
            for  (int w=0; w<4; w++) {
                vec2 p = vert[vi];
                if (w==0) p.y -= edgeLen/3;
                if (w==1) p.x += edgeLen/3;
                if (w==2) p.y += edgeLen/3;
                if (w==3) p.x -= edgeLen/3;
                f<<p.x<<" "<<p.y<<" 0.2 "<<r<<" "<<g<<" "<<b<<" "<<" 255\n";
            }
        }
    }

    for (int y=0; y<sy-1; y++)
    for (int x=0; x<sx-1; x++) {
        int i,j,k;
        i=grid[indexOf(x,y)]; j=grid[indexOf(x,y+1)]; k=grid[indexOf(x+1,y+1)];
        if ((i!=-1) && (j!=-1) && (k!=-1)) { f<<"3 "<<i<<" "<<k<<" "<<j<<"\n"; }
        i=grid[indexOf(x+1,y+1)]; j=grid[indexOf(x+1,y)]; k=grid[indexOf(x,y)];
        if ((i!=-1) && (j!=-1) && (k!=-1)) { f<<"3 "<<i<<" "<<k<<" "<<j<<"\n"; }
    }

    int c = vert.size();
    for (uint vi=0; vi<vert.size(); vi++) {
        if ((posInGrid[vi]==-1)&& (vdesired[vi]!=-1)) {
            int aa,bb,cc,dd;
            aa=c++; bb=c++; cc=c++; dd=c++;
            f<<"3 "<<aa<<" "<<bb<<" "<<vi<<"\n";
            f<<"3 "<<bb<<" "<<cc<<" "<<vi<<"\n";
            f<<"3 "<<cc<<" "<<dd<<" "<<vi<<"\n";
            f<<"3 "<<dd<<" "<<aa<<" "<<vi<<"\n";

        }
    }
    f.close();

    return true;
}

bool Grid::exportPLY(const std::string& filename ) const{
    int nUnassigned = 0;
    for (uint vi=0; vi<posInGrid.size(); vi++)
        if ((posInGrid[vi]==-1)&&(vdesired[vi]!=-1)) nUnassigned++;

    int nv = vert.size();
    int nf = 0;
    for (int y=0; y<sy-1; y++)
    for (int x=0; x<sx-1; x++) {
        int i,j,k;
        i=grid[indexOf(x,y)]; j=grid[indexOf(x,y+1)]; k=grid[indexOf(x+1,y+1)];
        if ((i!=-1) && (j!=-1) && (k!=-1)) { nf++;}
        i=grid[indexOf(x+1,y+1)]; j=grid[indexOf(x+1,y)]; k=grid[indexOf(x,y)];
        if ((i!=-1) && (j!=-1) && (k!=-1)) { nf++;}
    }

    nf += nUnassigned*4;
    nv += nUnassigned*4;

    std::ofstream f;
    f.open(filename.data(),std::fstream::out);
    if (!f.is_open()) return false;
    f
    <<"ply\n"
    <<"format ascii 1.0\n"
    <<"element vertex "<< nv <<"\n"
    <<"property float x\n"
    <<"property float y\n"
    <<"property float z\n"
    <<"property uchar red\n"
    <<"property uchar green\n"
    <<"property uchar blue\n"
    <<"property uchar alpha\n"
    <<"element face "<< nf<<"\n"
    <<"property list uchar int vertex_indices\n"
    <<"end_header\n";

    for (uint vi=0; vi<vert.size(); vi++) {
        int r=255,g=255,b=255;
        if (posInGrid[vi]==-1) {g/=2; b/=2;}
        else if (madeUpVert[vi]) {r/=4; b/=4; g=200;}
        f<<vert[vi].x<<" "<<vert[vi].y<<" 0 "<<r<<" "<<g<<" "<<b<<" "<<" 255\n";
    }
    for (uint vi=0; vi<vert.size(); vi++) {
        int r=255,g=255,b=255;
        if ((posInGrid[vi]==-1)&&(vdesired[vi]!=-1)) {
            g/=2; b/=2;
            for  (int w=0; w<4; w++) {
                vec2 p = vert[vi];
                if (w==0) p.y -= edgeLen/3;
                if (w==1) p.x += edgeLen/3;
                if (w==2) p.y += edgeLen/3;
                if (w==3) p.x -= edgeLen/3;
                f<<p.x<<" "<<p.y<<" 0.2 "<<r<<" "<<g<<" "<<b<<" "<<" 255\n";
            }
        }
    }

    for (int y=0; y<sy-1; y++)
    for (int x=0; x<sx-1; x++) {
        int i,j,k;
        i=grid[indexOf(x,y)]; j=grid[indexOf(x,y+1)]; k=grid[indexOf(x+1,y+1)];
        if ((i!=-1) && (j!=-1) && (k!=-1)) { f<<"3 "<<i<<" "<<k<<" "<<j<<"\n"; }
        i=grid[indexOf(x+1,y+1)]; j=grid[indexOf(x+1,y)]; k=grid[indexOf(x,y)];
        if ((i!=-1) && (j!=-1) && (k!=-1)) { f<<"3 "<<i<<" "<<k<<" "<<j<<"\n"; }
    }

    int c = vert.size();
    for (uint vi=0; vi<vert.size(); vi++) {
        if ((posInGrid[vi]==-1)&&(vdesired[vi]!=-1)) {
            int aa,bb,cc,dd;
            aa=c++; bb=c++; cc=c++; dd=c++;
            f<<"3 "<<aa<<" "<<bb<<" "<<vi<<"\n";
            f<<"3 "<<bb<<" "<<cc<<" "<<vi<<"\n";
            f<<"3 "<<cc<<" "<<dd<<" "<<vi<<"\n";
            f<<"3 "<<dd<<" "<<aa<<" "<<vi<<"\n";

        }
    }
    f.close();

    return true;
}

void Grid::exportEigen(Eigen::MatrixXi &tris, Eigen::MatrixXd &newPoints) const
{
    int nMadeUp  = 0;
    for (uint i=0; i<vert.size(); i++) if (madeUpVert[i]) nMadeUp++;

    // posInGrid[i]==-1 this are the points unused

    int nf = 0;
    for (int y=0; y<sy-1; y++)
    for (int x=0; x<sx-1; x++) {
        int i,j,k;
        i=grid[indexOf(x,y)]; j=grid[indexOf(x,y+1)]; k=grid[indexOf(x+1,y+1)];
        if ((i!=-1) && (j!=-1) && (k!=-1)) { nf++;}
        i=grid[indexOf(x+1,y+1)]; j=grid[indexOf(x+1,y)]; k=grid[indexOf(x,y)];
        if ((i!=-1) && (j!=-1) && (k!=-1)) { nf++;}
    }

    newPoints.resize(nMadeUp, 3);
    int index = 0;
    for (uint i=0; i<vert.size(); i++) {
        if (madeUpVert[i]) {
            const vec2 &p = vert[i];
            newPoints(index, 0) = p.x;
            newPoints(index, 1) = p.y;
            newPoints(index, 2) = 0;
            index++;
        }
    }

    tris.resize(nf, 3);
    index = 0;
    for (int y=0; y<sy-1; y++)
    for (int x=0; x<sx-1; x++) {
        int i,j,k;
        i=grid[indexOf(x,y)]; j=grid[indexOf(x,y+1)]; k=grid[indexOf(x+1,y+1)];
        if ((i!=-1) && (j!=-1) && (k!=-1)) {
            tris(index, 0) = i;
            tris(index, 1) = k;
            tris(index, 2) = j;
            ++index;
        }
        i=grid[indexOf(x+1,y+1)]; j=grid[indexOf(x+1,y)]; k=grid[indexOf(x,y)];
        if ((i!=-1) && (j!=-1) && (k!=-1)) {
            // f<<"3 "<<i<<" "<<j<<" "<<k<<"\n";
            tris(index, 0) = i;
            tris(index, 1) = k;
            tris(index, 2) = j;
            ++index;
        }
    }
}


Grid::Grid(){
    sx = sy = 0;
}

void Grid::updateDistFromBorder(){
    distFromBorder.resize(sx*sy);
    for (int gi=0; gi<sx*sy; gi++) distFromBorder[gi] = (grid[gi]==-1)?0:1000;
    while (1) {
        bool over = true;
        for (int gi=safeGiMin; gi<safeGiMax; gi++) {
            int d = distFromBorder[gi]+1;
            for (int n=0; n<6; n++) {
                if (distFromBorder[gi+neigh[n]]>d) {
                    over = false;
                    distFromBorder[gi+neigh[n]] = d;
                }
            }
        }
        if (over) break;
    }
}



bool Grid::isBoundary( int gi) const{
    int j;
    j = grid[ gi-sx-1 ]; if (j==-1) return true;
    j = grid[ gi-sx   ]; if (j==-1) return true;
    j = grid[ gi-1    ]; if (j==-1) return true;
    j = grid[ gi+1    ]; if (j==-1) return true;
    j = grid[ gi+sx   ]; if (j==-1) return true;
    j = grid[ gi+sx+1 ]; if (j==-1) return true;
    return false;
}

void Grid::swapTwo(int gi, int gj){
    std::swap( grid[gi], grid[gj] );
    int vi = grid[gi];
    int vj = grid[gj];
    if (vi!=-1) posInGrid[vi] = gi;
    if (vj!=-1) posInGrid[vj] = gj;
    //std::swap( posInGrid[grid[gi]], posInGrid[grid[gj]]);
}

vec2 Grid::baryAroundOfExisting(int gi) const{
    vec2 res(0,0);
    int count = 0;
    int j;
    j = grid[ gi-sx-1 ]; if (j!=-1) { res+=vert[j];count++;}
    j = grid[ gi-sx   ]; if (j!=-1) { res+=vert[j];count++;}
    j = grid[ gi-1    ]; if (j!=-1) { res+=vert[j];count++;}
    j = grid[ gi+1    ]; if (j!=-1) { res+=vert[j];count++;}
    j = grid[ gi+sx   ]; if (j!=-1) { res+=vert[j];count++;}
    j = grid[ gi+sx+1 ]; if (j!=-1) { res+=vert[j];count++;}
    return res / count;

}

vec2 Grid::baryAround(int gi) const{
    vec2 res(0,0);
    vec2 c = vert[ grid[gi] ];
    int j;
    j = grid[ gi-sx-1 ]; res+=(j==-1)?c:vert[j];
    j = grid[ gi-sx   ]; res+=(j==-1)?c:vert[j];
    j = grid[ gi-1    ]; res+=(j==-1)?c:vert[j];
    j = grid[ gi+1    ]; res+=(j==-1)?c:vert[j];
    j = grid[ gi+sx   ]; res+=(j==-1)?c:vert[j];
    j = grid[ gi+sx+1 ]; res+=(j==-1)?c:vert[j];
    return res / 6;
}

vec2 Grid::avgPos(int gi) const{
    vec2 res(0,0);
    int div=0;
    int j;
    j = grid[ gi-sx-1 ]; if (j!=-1) { res+= vert[j]; div++;}
    j = grid[ gi-sx   ]; if (j!=-1) { res+= vert[j]; div++;}
    j = grid[ gi-1    ]; if (j!=-1) { res+= vert[j]; div++;}
    j = grid[ gi+1    ]; if (j!=-1) { res+= vert[j]; div++;}
    j = grid[ gi+sx   ]; if (j!=-1) { res+= vert[j]; div++;}
    j = grid[ gi+sx+1 ]; if (j!=-1) { res+= vert[j]; div++;}
    if (div==0) return vec2(0,0);

    res/=div;
    return res;
}

scalar Grid::energyAround2(int gi) const{
    scalar res = 0, div = -0.5;
    scalar e;
    e = energyBetween2(grid[gi],grid[gi-sx-1]); if (e>0) {res+=e; div++;}
    e = energyBetween2(grid[gi],grid[gi-sx]);   if (e>0) {res+=e; div++;}
    e = energyBetween2(grid[gi],grid[gi-1]);    if (e>0) {res+=e; div++;}
    e = energyBetween2(grid[gi],grid[gi+1]);    if (e>0) {res+=e; div++;}
    e = energyBetween2(grid[gi],grid[gi+sx]);   if (e>0) {res+=e; div++;}
    e = energyBetween2(grid[gi],grid[gi+sx+1]); if (e>0) {res+=e; div++;}
    if (div<=0) return edgeLen*edgeLen*40.0;
    return res/div;
}
scalar Grid::energyAround2(int gi, int gj) const{
    return energyAround2(gi)+energyAround2(gj);
}

scalar Grid::energyAroundExcept1(int gi) const{
    scalar e[6];
    e[0] = energyBetween(grid[gi],grid[gi-sx-1]);
    e[1] = energyBetween(grid[gi],grid[gi-sx]);
    e[2] = energyBetween(grid[gi],grid[gi-1]);
    e[3] = energyBetween(grid[gi],grid[gi+1]);
    e[4] = energyBetween(grid[gi],grid[gi+sx]);
    e[5] = energyBetween(grid[gi],grid[gi+sx+1]);
    int max = 0;
    for (int i=1; i<6; i++) if (e[i]>e[max]) max = i;
    e[max] = 0;
    return e[0]+e[1]+e[2]+e[3]+e[4]+e[5];
}

scalar Grid::energyAround(int gi) const{
    return energyBetween(grid[gi],grid[gi-sx-1])+
           energyBetween(grid[gi],grid[gi-sx])+
           energyBetween(grid[gi],grid[gi-1])+
           energyBetween(grid[gi],grid[gi+1])+
           energyBetween(grid[gi],grid[gi+sx])+
           energyBetween(grid[gi],grid[gi+sx+1]);
}
scalar Grid::energyAroundIf(int gi, int vi) const{
    return energyBetween(vi,grid[gi-sx-1])+
           energyBetween(vi,grid[gi-sx])+
           energyBetween(vi,grid[gi-1])+
           energyBetween(vi,grid[gi+1])+
           energyBetween(vi,grid[gi+sx])+
           energyBetween(vi,grid[gi+sx+1]);
}


int Grid::friendsAround(int gi) const{
    return ((grid[gi-sx-1]>=0)?1:0)+
           ((grid[gi-sx]>=0)?1:0)+
           ((grid[gi-1]>=0)?1:0)+
           ((grid[gi+1]>=0)?1:0)+
           ((grid[gi+sx]>=0)?1:0)+
           ((grid[gi+sx+1]>=0)?1:0);
}

scalar Grid::energyAround(int gi, int gj) const{
    return energyAround(gi)+energyAround(gj);
}

scalar Grid::energyTotal() const{
    scalar res = 0.0;
    for (int y=1; y<sy; y++)
    for (int x=0; x<sx; x++){
        int gi=indexOf(x,y);
        res+=energyBetween(grid[gi],grid[gi-sx-1])+
             energyBetween(grid[gi],grid[gi-sx])+
             energyBetween(grid[gi],grid[gi-1]);
    }
    return res;
}

scalar Grid::triangleQuality(int vi, int vj, int vk) const{
    if (vi<0) return -10;
    if (vj<0) return -10;
    if (vk<0) return -10;
    vec2 e0 = (vert[vi]-vert[vj]);
    vec2 e1 = (vert[vj]-vert[vk]);
    vec2 e2 = (vert[vk]-vert[vi]);
    scalar a0 = e0.norm();
    scalar a1 = e1.norm();
    scalar a2 = e2.norm();

    if (a0<a1) std::swap(a0,a1);
    if (a1<a2) std::swap(a1,a2);
    if (a0<a1) std::swap(a0,a1);
    return 1.0 - (a2-a0)/edgeLen;

}
void Grid::trimBorders(){
    int tri[7];
    tri[0]=-sx-1;
    tri[1]=-sx;
    tri[2]=1;
    tri[3]=sx+1;
    tri[4]=sx;
    tri[5]=-1;
    tri[6]=tri[0];

    int count = 0;
    for (int gi=safeGiMin; gi<safeGiMax; gi++) {
        if (grid[gi]==-1) continue;
        int nt = 0, ng = 0; // num tri, num good tri
        for (int n=0; n<6; n++) {
            scalar q = triangleQuality(grid[gi],grid[gi+tri[n]],grid[gi+tri[n+1]]);
            if ( q > 0.95) ng++;
            if ( q > -1) nt++;
            if (nt>2) break;
        }
        if ((nt<2) && (ng<2) ) { unassign(gi); count++;}

    }
    std::cout<<"trimmed "<<count<<" verts at boundaries\n";
}

void Grid::initIndicesOnGrid(int nx,int ny){

    for (int y=0,k=0; y<ny; y++)
    for (int x=0; x<nx; x++,k++)
    {
        int gi = indexOf(x+1+y/2, y+1);
        grid[ gi ] = k;
        posInGrid[ k ] = gi;
    }
}

void Grid::sanityCheck(){
    for (int gi=0; gi<(int)grid.size(); gi++) if (grid[gi]!=-1)
        myAssert( posInGrid[ grid[gi] ] == gi, "Wrong at gi:"<<gi<<" vi:"<<grid[gi]);
}

void Grid::updatePosInGrid(){
    posInGrid.assign(vert.size(),-1);
    for (uint gi=0; gi<grid.size(); gi++) if (grid[gi]!=-1) posInGrid[ grid[gi] ] = gi;
}

void Grid::create(int _sx, int _sy){

    sx = _sx;
    sy = _sy;
    int n = sx*sy;
    grid.resize(n);
    grid.assign(n,-1);

    neigh[0]=sx+1;
    neigh[1]=sx;
    neigh[2]=1;
    neigh[3]=-1;
    neigh[4]=-sx;
    neigh[5]=-sx-1;

    safeGiMin = (sx+1);
    safeGiMax = (sx*sy)-(sx+1);

    safeGiMinS2 = (sx+1)*2;
    safeGiMaxS2 = (sx*sy)-(sx+1)*2;
    safeGiMinS3 = (sx+1)*3;
    safeGiMaxS3 = (sx*sy)-(sx+1)*3;
    //std::cout<<"N = "<<n<<"\n";
}


int Grid::shiftPos(int gi, int dir) const{
    switch(dir){
    case 0: return gi-sx;
    case 1: return gi+1;
    case 2: return gi+sx+1;
    case 3: return gi+sx;
    case 4: return gi-1;
    case 5: return gi-sx-1;
    }

    assert(false);
    return 0;
}


scalar Grid::energyBetween(int vi, int vj) const{
    if (vi==-1 || vj==-1) return edgeLen*edgeLen*2.0;
    return squaredDistance(vert[vi],vert[vj]);
}

scalar Grid::energyBetween2(int vi, int vj) const{
    if (vi==-1 || vj==-1) return -1;
    return squaredDistance(vert[vi],vert[vj]);
}

void Grid::printf() const{
    for (int y=0,k=0; y<sy; y++) {
        for (int i=0; i<(sy-y-1); i++) std::cout<<"  ";
        for (int x=0; x<sx; x++,k++){
            if (grid[k]==-1) std::cout<<"*** ";
            else std::cout << std::setfill('0') << std::setw(3) << grid[k] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "Eng = "<<energyTotal()<<"\n";
}


bool Grid::testAndDoSwapBordersIncluded(int gi, int gj){
    if(grid[gi]==-1 && grid[gj]==-1) return false;
    scalar before = energyAround2(gi,gj);
    swapTwo( gi, gj );
    scalar after = energyAround2(gi,gj);
    if (after+0.00001<before) {
        return true;
    }
    swapTwo( gi, gj );
    return false;
}
//scalar lastGain;
bool Grid::testAndDoSwap(int gi, int gj){
    if((grid[gi]==-1) || (grid[gj]==-1)) return false;
    scalar before = energyAround(gi,gj);
    swapTwo( gi, gj );
    scalar after = energyAround(gi,gj);
    if (after+0.00001<before) {
        //lastGain =
        return true;
    }
    swapTwo( gi, gj );
    return false;
}

bool Grid::testAndDoTriSwap(int gi, int gj, int gk){
    if((grid[gi]==-1) || (grid[gj]==-1) || (grid[gk]==-1)) return false;
    scalar before = energyAround(gi)+energyAround(gj)+energyAround(gk);
    swapTwo( gi, gj );
    swapTwo( gj, gk );
    scalar after = energyAround(gi)+energyAround(gj)+energyAround(gk);
    if (after+0.00001<before) {
        return true;
    }
    swapTwo( gj, gk );
    swapTwo( gi, gj );
    return false;
}

bool Grid::testAndDoQuadriSwap(int gi, int gj, int gk, int gh){
    if((grid[gi]==-1) || (grid[gj]==-1) || (grid[gk]==-1)|| (grid[gh]==-1)) return false;
    scalar before = energyAround(gi)+energyAround(gj)+energyAround(gk)+energyAround(gh);
    swapTwo( gi, gj );
    swapTwo( gj, gk );
    swapTwo( gk, gh );
    scalar after = energyAround(gi)+energyAround(gj)+energyAround(gk)+energyAround(gh);
    if (after+0.00001<before) {
        return true;
    }
    swapTwo( gk, gh );
    swapTwo( gj, gk );
    swapTwo( gi, gj );
    return false;
}


int Grid::fillGapsMakingPtsUp(){

    int count = 0;
    while (1) {
        bool over = true;
        for (int gi=safeGiMin; gi<safeGiMax; gi++) if (grid[gi]==-1){
            int nf = friendsAround(gi);
            //std::cout<<nf<<" ";
            if ( (nf==5)|| (nf==6)) {
                int vi = (int)vert.size();

                vert.push_back( baryAroundOfExisting(gi) );
                posInGrid.push_back(gi);
                madeUpVert.push_back(true);

                myAssert(posInGrid.size()==vert.size(),"Size mismatch!");
                grid[gi] = vi;
                count++;
                over = false;
            }
        }
        break;
        if (over) break;
    }
    std::cout<<"Filled with "<<count<<" made up points\n";
    return count;

}

int Grid::tryAllSwapsAround(int gi){
    int res = 0;
    for(int n=0; n<6; n++)
    if (testAndDoSwap(gi,gi+neigh[n])) res++;
    return res;
}

int Grid::tryAllSwaps(){
    return tryAllBiSwaps()+tryAllTriSwaps()+tryAllQuadriSwaps();
}

int Grid::tryAllBiSwaps(){

    //edgeLen = computeAvgEdgeLen();

    int count = 0;
    while (1) {
        int pass = 0;
        for (int i=safeGiMinS2; i<safeGiMaxS2; i++) {
            if (testAndDoSwap(i,i+1)) pass++;
            if (testAndDoSwap(i,i+sx+1)) pass++;
            if (testAndDoSwap(i,i+sx)) pass++;
        }
        if (pass==0) break;
        count+=pass;
    }
    std::cout<<"Done "<<count<<" greedy swaps;\n";
    return count;
}

int Grid::tryAllTriSwaps(){
    int count = 0;
    while (1) {
        int pass = 0;
        for (int i=safeGiMinS2; i<safeGiMaxS2; i++) {
            if (testAndDoTriSwap(i,i+sx,i+sx+1)) pass++;
            if (testAndDoTriSwap(i,i+sx+1,i+sx)) pass++;
            if (testAndDoTriSwap(i,i+1,i+sx+1)) pass++;
            if (testAndDoTriSwap(i,i+sx+1,i+1)) pass++;
        }
        if (pass==0) break;
        count+=pass;
    }
    std::cout<<"Done "<<count<<" greedy three-swaps;\n";
    return count;
}

int Grid::tryAllQuadriSwaps(){
    int count = 0;
    while (1) {
        int pass = 0;
        for (int i=safeGiMinS3; i<safeGiMaxS3; i++) {
            int a,b,c,d;
            a = i; b = i+1; c=i+sx; d=i+sx+1;
            if (testAndDoQuadriSwap(a,b,c,d)) pass++;
            if (testAndDoQuadriSwap(a,b,d,c)) pass++;
            if (testAndDoQuadriSwap(a,c,b,d)) pass++;
            if (testAndDoQuadriSwap(a,c,d,b)) pass++;
            if (testAndDoQuadriSwap(a,d,b,c)) pass++;
            if (testAndDoQuadriSwap(a,d,c,b)) pass++;
            a = i; b = i+1; c=i+sx+1; d=i+sx+2;
            if (testAndDoQuadriSwap(a,b,c,d)) pass++;
            if (testAndDoQuadriSwap(a,b,d,c)) pass++;
            if (testAndDoQuadriSwap(a,c,b,d)) pass++;
            if (testAndDoQuadriSwap(a,c,d,b)) pass++;
            if (testAndDoQuadriSwap(a,d,b,c)) pass++;
            if (testAndDoQuadriSwap(a,d,c,b)) pass++;
            a = i; b = i+1; c=i-sx; d=i+sx+1;
            if (testAndDoQuadriSwap(a,b,c,d)) pass++;
            if (testAndDoQuadriSwap(a,b,d,c)) pass++;
            if (testAndDoQuadriSwap(a,c,b,d)) pass++;
            if (testAndDoQuadriSwap(a,c,d,b)) pass++;
            if (testAndDoQuadriSwap(a,d,b,c)) pass++;
            if (testAndDoQuadriSwap(a,d,c,b)) pass++;
        }
        if (pass==0) break;
        count+=pass;
    }
    std::cout<<"Done "<<count<<" greedy quadri-swaps;\n";
    return count;
}


int Grid::tryAllSwapsBordersIncluded(){

    //edgeLen = computeAvgEdgeLen();

    int count = 0;
    while (1) {
        int pass = 0;
        for (int i=safeGiMinS2; i<safeGiMaxS2; i++) {
            if (testAndDoSwapBordersIncluded(i,i+1)) pass++;
            if (testAndDoSwapBordersIncluded(i,i+sx+1)) pass++;
            if (testAndDoSwapBordersIncluded(i,i+sx)) pass++;
        }
        if (pass==0) break;
        count+=pass;
    }
    std::cout<<"Done "<<count<<" greedy swaps;\n";
    return count;
}

void Grid::initVertOnGrid(int sx, int sy){
    vert.clear();
    for (int y=0; y<sy; y++){
         for (int x=0; x<sx; x++){
            vert.push_back(vec2(  (x - ((y%2)?0.5f:0.0f))*edgeLen ,
                                  y*(edgeLen*sqrt(3)/2.0)
                               ));
        }
    }
    posInGrid.resize(sx*sy);
}


void Grid::sanityCheck( int testX, int testY ){

    int gj=indexOf(testX,testY);
    vec2 p = vert[ grid[gj] ];
    std::cout << "At ["<<testX<<","<<testY<<"]: vert "<<grid[gj]<<" ("<<p.x<<","<<p.y<<")\n";
    vec2 q = baryAround(gj);
    std::cout << "Baricenter around it: ("<<q.x<<","<<q.y<<")\n";
    /*spatialIndex.init();
    spatialIndex.setTarget(q);
    int k = spatialIndex.nextClosest();
    int h = posInGrid[k];
    std::cout << "Closest to that bary: vert "<<k<<" ["<<h%sx<<","<<h/sx<<"]\n";
    std::cout << "coords: ("<<vert[k].x<<","<<vert[k].y<<")\n";*/
}


void Grid::assign(int gi, int vi){
    if (gi!=-1) grid[gi] = vi;
    if (vi!=-1) posInGrid[vi] = gi;
}
void Grid::unassign(int gi){
    int vi = grid[gi];
    grid[gi] = -1;
    if (vi!=-1) {
        posInGrid[vi] = -1;
        vdesired[vi] = -1;
    }
}




int Grid::numberOfAdj(int gi) const{
    int res=0;
    for (int i=0; i<6; i++) {
        if (grid[ gi+neigh[i] ]!=-1) res++;
    }
    return res;
}

void Grid::enlargeToInclude(int gi,int buffer){
    int x = gi%sx;
    int y = gi/sx;
    if (x<2) { enlargeGrid(buffer-x,0,0,0); x=buffer; }
    if (y<2) { enlargeGrid(0,0,buffer-y,0); y=buffer; }
    if (x>sx-3) { enlargeGrid(0,x-(sx-buffer-1),0,0); x=sx-buffer-1; }
    if (y>sy-3) { enlargeGrid(0,0,0,y-(sy-buffer-1)); y=sy-buffer-1; }
}

void Grid::enlargeGrid(int dxMin, int dxMax, int dyMin, int dyMax){
    auto backup = grid;
    int oldSx = sx;
    int oldSy = sy;
    auto oldToNew = [oldSx,dxMin,dyMin,this](int gi)->int{
        if (gi==-1) return -1;
        return indexOf( (gi%oldSx)+dxMin, gi/oldSx+dyMin);
    };
    sx += dxMin+dxMax;
    sy += dyMin+dyMax;
    create(sx,sy);
    for (int i=0; i<oldSx*oldSy; i++) {
        grid[ oldToNew(i) ] = backup[ i ];
    }
    for (int &p:posInGrid) p = oldToNew(p);

}

scalar squaredOf(scalar x){return x*x;}


int Grid::hopDistance(int gi, int gj) const{
    int xi = gi%sx;
    int yi = gi/sx;
    int xj = gj%sx;
    int yj = gj/sx;
    int xjb = xj-(yj-yi);
    return std::min(
                std::abs(xi-xj)+std::abs(yi-yj),
                std::abs(xi-xjb)+std::abs(yi-yj)
                ); // TODO! -std::max(0,std::)
}

static bool contains(const std::vector<int> &except, int i){
    for (int j:except) if(i==j) return true;
    return false;
}


bool Grid::fixUnassignedVertexNiceWay(int vi){

    while (1){
        int gi = vdesired[vi];
        if (gi==-1) return false;
        int vj = grid[gi];
        if (vj==-1) {
            // found an empty spot!
            assign(gi,vi);
            return true;
        } else {
            scalar before = energyAroundExcept1(gi);
            assign(gi,vi);
            scalar after =  energyAroundExcept1(gi);
            if (before<after) {
                // before was better: undo swap
                assign(gi,vj);
                posInGrid[vi] = -1;
                continue; // try again
                //return false;
            } else {
                return fixUnassignedVertexNiceWay(vj);
            }
        }
    }
}

bool Grid::fixUnassignedVertexDijkstra(int vi){

    scalar maxCost = edgeLen*edgeLen*6*5;
    std::vector<scalar> cost(grid.size(),maxCost);
    std::vector<int> prevStep(grid.size(),-2);
    std::set<int> boundary;
    std::set<int> visited;
    int dest;

    int gi = vdesired[vi]; //bestPositionFor(vi);
    if (gi==-1) {
        std::cout<<"Assign "<<vi<<": I can't even.\n";
        return false;
    }
    prevStep[gi] = -1;
    boundary.insert(gi);
    cost[gi] = 0;
    bool pathFound = false;
    while (!pathFound) {
        if (boundary.empty()) {
            std::cout<<"Assign "<<vi<<": NOT doing it (too expensive).\n";
            return false;
        }

        // find least expensive node on boundary
        int gi = -1; scalar best = 10000000;
        for (int gj:boundary) {
            if (cost[gj]<best)  { best = cost[gj]; gi = gj;}
        }
        boundary.erase(gi);
        visited.insert(gi);

        for (int n=0; n<6; n++) {
            int gj = gi+neigh[n];
            if (visited.find(gj)!=visited.end()) continue;
            if (grid[gj]==-1) {
                prevStep[gj] = gi;
                pathFound = true;
                dest = gj;
                break;
            }
            scalar newCost = cost[gi] + energyAroundIf(gj, grid[gi] ) - energyAround(gj);
            if (cost[gj] > newCost) {
                cost[gj] = newCost;
                boundary.insert(gj);
                prevStep[gj] = gi;
            }
        }
    }
    std::cout<<"Assign "<<vi<<": path found. Applyting it: ";
    {
    //posInGrid[vi]=vdesired[vi];
    int gi = dest;
    while (1) {
        myAssert( grid[gi]==-1 ,"not carry the -1");
        int gj = prevStep[gi];

        //std::cout<<gi<<"<--"<<gj<<"  ";
        if (gj==-1) {
            myAssert((gi!=-1)&&(vi!=-1),"Wrong wrong");
            assign(gi,vi);
            std::cout<<"*\n";
            break;
        }

        std::cout<<"-";

        swapTwo( gj, gi );
        tryAllSwapsAround(gj);
        tryAllSwapsAround(gi);
        gi = gj;
    }
    //std::cout<<"done "<<vi<<" (grid Pos:"<<posInGrid[vi]<<")\n";
    return true;
    }
}


int Grid::assignUnassignedNiceWay(){
    int countYes = 0;
    int countNo=0;
    for (int vi=0; vi<(int)vert.size(); vi++) {
        if ((posInGrid[vi]==-1) && (vdesired[vi]!=-1)) {
            if (fixUnassignedVertexNiceWay(vi)) {
                countYes++;
                std::cout<<"Fix N."<<countYes <<" Nice...\n";
            } else countNo++;
        }
    }
    std::cout<<"Assigned "<<countYes<<" new pts ("<<countNo<<" left)\n";
    assert(false);
    return 0;
}


int Grid::assignUnassignedHardWay(){

    int countYes = 0;
    std::vector<int> toFix;
    for (int vi=0; vi<(int)vert.size(); vi++) {
        if ((posInGrid[vi]==-1) && (vdesired[vi]!=-1)) toFix.push_back(vi);
    }

    updateDistFromBorder();

    std::cout << "There are "<<toFix.size()<<" verts to fix!\n";

    while (toFix.size()>0) {
        int maxDist = -1;
        int win = -1;
        for (int i=0; i<(int)toFix.size(); i++) {
            int vi = toFix[i];
            int dist = distFromBorder[ vdesired[vi] ]; //bestPositionFor(vi) ];
            if (dist>maxDist) { maxDist = dist; win=i;}
        }
        myAssert(win!=-1,"I cannot program sort\n");
        int vi = toFix[ win ];
        std::swap(toFix[win],toFix.back()); toFix.pop_back();

        if (fixUnassignedVertexDijkstra( vi )) {
           countYes++;
           //std::cout<<"Fix N."<<countYes <<" Hard...\n";
           updateDistFromBorder();
        }

    }
    std::cout<<"Fixed "<<countYes<<" unassigned points ("<<toFix.size()<<" left)\n";
    return countYes;

}

}}
