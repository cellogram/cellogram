#include<igl/delaunay_triangulation.h>
#include"mesh.h"


namespace cellogram
{
namespace PointsUntangler
{

short sign(double d){
    if (d<0) return -1; else if (d>0) return +1; else return 0;
}

short orient2D(
    const scalar pa[2],
    const scalar pb[2],
    const scalar pc[2]){

    return sign( cross( vec2(pa)-vec2(pc) , vec2(pb)-vec2(pc) ) );
}

inline vec2 circleCenter(vec2 A, vec2 B, vec2 C) {
    vec2 ab = A + B;
    vec2 bc = B + C;

    vec2 ba = B - A;
    vec2 cb = C - B;
    vec2 ac = A - C;

    vec2 center;

    center.x = (ba.x*cb.y*ab.x - ba.y*bc.x*cb.x + ba.y*cb.y*ac.y)/(2*cross(ba,cb));
    center.y = ab.y/2 - (center.x - ab.x/2)*ba.x/ba.y;

    return center;
}

IGL_INLINE short incircle(
    const scalar pa[2],
    const scalar pb[2],
    const scalar pc[2],
    const scalar pd[2])
{
    vec2 c = circleCenter(pa,pb,pc);
    return sign( squaredDistance(pa,c) - squaredDistance(pd,c));
}

void Mesh::delaunay(){
    Eigen::MatrixXd VV;
    Eigen::MatrixXi FF;

    VV.resize(V.size(),2);
    for (int i=0; i<V.size(); i++) { VV(i,0)=V[i].p.x; VV(i,1)=V[i].p.y; }
    igl::delaunay_triangulation(VV, orient2D, incircle,  FF);

    F.resize(FF.rows());
    for (int i=0; i<FF.rows(); i++) for (int w=0; w<3; w++) { F[i].vi[w] = FF(i,w); }

    updateIndices();
    updateValencies();
}

}
}
