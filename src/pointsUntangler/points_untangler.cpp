#include "points_untangler.h"

#include<iostream>
#include"mesh.h"
#include"grid.h"

namespace cellogram
{
namespace PointsUntangler
{

// std::string outputPath = "C:/projects/sandbox/tarini/cell_forces/cell_forces/";

void pointsUntangler(const Eigen::MatrixXd &detected, Eigen::MatrixXi &tris, Eigen::MatrixXd &newPoints)
{
    Mesh m;
    Grid g;
    m.fromEigen(detected);

    //m.fuckUp();

    std::string outputPath = "./";

    m.delaunay();
    // m.exportPLY(outputPath+"1_0Delonay.ply");

    for (char c='A'; c<='D'; c++) {
        m.regularizeByFlips( 1 , 1 );  // ( 10 , 30 );
        // m.exportPLY(outputPath+"1_1Flips"+c+".ply");
        meshToGrid(m,g);
        // m.exportPLY(outputPath+"1_1Floodfill"+c+".ply", true);
        // g.exportPLY(outputPath+"1_1Grid"+c+".ply");
        m.smoothDisputed();
        m.smoothDisputed();
    }

    for (int i=0; i<1; i++) {
        int nDone = 0;
        nDone += g.tryAllSwaps();
        if (nDone!=0) break;
    }
    // g.exportPLY(outputPath+"1_3AssignSoft.ply");

	g.assignUnassignedNiceWay();
    for (int i=0; i<3; i++) {
        int nDone = 0;
		
        nDone += g.assignUnassignedHardWay();
        nDone += g.tryAllSwaps();
        if (nDone==0) break;
    }
    // g.exportPLY(outputPath+"1_4AssignHard.ply");

    g.tryAllSwapsBordersIncluded();
    g.fillGapsMakingPtsUp();
    g.tryAllSwaps();

    g.fillGapsMakingPtsUp();

    g.exportEigen(tris, newPoints);
    g.exportPLY(outputPath+"1_5Final.ply");
    // g.exportPLYtartan(outputPath+"1_5FinalT.ply");

}

}
}
