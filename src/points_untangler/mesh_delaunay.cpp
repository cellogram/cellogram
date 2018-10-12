#include<igl/delaunay_triangulation.h>
// #include<igl/incircle.h>
// #include<igl/orient2D.h>
#include<iostream>
#include"mesh.h"


namespace cellogram
{
namespace PointsUntangler
{

short sign(scalar d){
    if (d<0) return -1; else if (d>0) return +1; else return 0;
}

inline short orient2D( const scalar A[2], const scalar B[2], const scalar C[2])
{
    scalar ACx = A[0]-C[0], ACy = A[1]-C[1];
    scalar CBx = C[0]-B[0], CBy = C[1]-B[1];
    return sign( CBx*ACy - CBy*ACx );
}

inline short incircle( const scalar A[2], const scalar B[2], const scalar C[2], const scalar D[2])
{

    scalar BAx = B[0]-A[0], BAy = B[1]-A[1];
    scalar ACx = A[0]-C[0], ACy = A[1]-C[1];
    scalar CBx = C[0]-B[0], CBy = C[1]-B[1];
    scalar BDx = B[0]-D[0], BDy = B[1]-D[1];
    scalar ADx = A[0]-D[0], ADy = A[1]-D[1];

    scalar Txx = BAx*CBx;
    scalar Tyy = BAy*CBy;
    scalar Txy = BAx*CBy;
    scalar Tyx = BAy*CBx;

    return sign( (Tyy*ACy + Tyx*ACx - (Tyx-Txy)*BDx )*(Tyx-Txy)*ADx -
                 (Txx*ACx + Txy*ACy + (Tyx-Txy)*BDy )*(Tyx-Txy)*ADy );

}

void Mesh::initWithDelaunay(){
    Eigen::MatrixXd VV;
    Eigen::MatrixXi FF;

    VV.resize(V.size(),2);
    for (uint i=0; i<V.size(); i++) { VV(i,0)=V[i].p.x; VV(i,1)=V[i].p.y; }
    igl::delaunay_triangulation(VV, orient2D, incircle,  FF);

    F.resize(FF.rows());
    for (int i=0; i<FF.rows(); i++) for (int w=0; w<3; w++) { F[i].vi[w] = FF(i,w); }

    buildEdgesFromFaces();
    updateValencies();
    if(verbose)
    std::cout<<"Delaunay done (" << nVal <<" irregulars)\n";
}

}
}
