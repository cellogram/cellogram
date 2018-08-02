#include "points_untangler.h"

#include<iostream>
#include"mesh.h"
#include"grid.h"


std::string inputPath = "C:/projects/github/cellogram/data/debug/macro_bad.txt";
std::string outputPath = "C:/projects/github/cellogram/data/debug/";

namespace cellogram
{
namespace PointsUntangler
{

Grid g;
Mesh m;

void test02(){
    //m.importFVFix(inputPath+"V.vert", inputPath+"F.tri", inputPath+"fixed.txt");
    m.importXYZv3(inputPath);
    pointsUntangler(m, g, outputPath);
    g.exportPLY(outputPath+"final0.ply");
    g.fillGapsMakingPtsUp();
    g.exportOBJ(outputPath+"final1.obj");

}

// void test01(){
//     Mesh m;
//     m.importFVFix(inputPath+"V.vert", inputPath+"F.tri", inputPath+"fixed.txt");
//     m.exportPLY(outputPath+"4_0.ply", ColMode::BY_VAL);

//     int D = 10;
//     m.greedyFlips(D);
//     m.exportPLY(outputPath+"4_1.ply", ColMode::BY_VAL);

//     meshToGrid(m,g);
//     //m.flipAs(g);
//     m.greedyFlips(D);
//     m.exportPLY(outputPath+"4_2.ply", ColMode::BY_VAL);
//     m.exportEdgesPLY(outputPath+"4_2_edges.ply",g);
//     m.exportEdgesPLY(outputPath+"4_2_edges_f.ply");
//     g.exportPLY(outputPath+"4_2_clean.ply");

//     meshToGrid(m,g);
//     //m.flipAs(g);
//     m.greedyFlips(D);
//     m.exportPLY(outputPath+"4_3.ply", ColMode::BY_VAL);
//     m.exportEdgesPLY(outputPath+"4_3_edges.ply",g);
//     g.exportPLY(outputPath+"4_3_clean.ply");

//     meshToGrid(m,g);
//     //m.flipAs(g);
//     m.greedyFlips(D);
//     m.exportPLY(outputPath+"4_4.ply", ColMode::BY_VAL);
//     m.exportEdgesPLY(outputPath+"4_4_edges.ply",g);
//     g.exportPLY(outputPath+"4_4_clean.ply");

//     meshToGrid(m,g);
//     //m.flipAs(g);
//     m.greedyFlips(D);
//     m.exportPLY(outputPath+"4_5.ply", ColMode::BY_VAL);
//     m.exportEdgesPLY(outputPath+"4_5_edges.ply",g);
//     g.exportPLY(outputPath+"4_5_clean.ply");

//     meshToGrid(m,g);
//     //m.flipAs(g);
//     m.greedyFlips(D);
//     m.exportPLY(outputPath+"4_6.ply", ColMode::BY_VAL);
//     m.exportPLY(outputPath+"4_6_FF.ply", ColMode::BY_FLOOD);
//     m.exportEdgesPLY(outputPath+"4_6_edges.ply",g);
//     g.exportPLY(outputPath+"4_6_clean.ply");

//     meshToGrid(m,g);
//     //m.flipAs(g);
//     m.greedyFlips(D);
//     m.exportPLY(outputPath+"4_7.ply", ColMode::BY_VAL);
//     m.exportPLY(outputPath+"4_7_FF.ply", ColMode::BY_FLOOD);
//     m.exportEdgesPLY(outputPath+"4_7_edges.ply",g);
//     g.exportPLY(outputPath+"4_7_clean.ply");

//     meshToGrid(m,g);
//     m.greedyFlips(D);
//     m.exportPLY(outputPath+"4_8.ply", ColMode::BY_VAL);
//     m.exportPLY(outputPath+"4_8_FF.ply", ColMode::BY_FLOOD);
//     m.exportEdgesPLY(outputPath+"4_8_edges.ply",g);
//     g.exportPLY(outputPath+"4_8_clean.ply");

//     /*
//     m.smoothDisputed();
//     m.smoothDisputed();

//     m.regularizeByFlips(5,5);

//     meshToGrid(m,g);
//     m.regularizeByFlips(5,5);
//     m.exportPLY(outputPath+"41_ff.ply",ColMode::BY_FLOOD);
//     g.exportPLY(outputPath+"41_out.ply");
//     m.exportPLY(outputPath+"41_di0.ply",ColMode::BY_DISPUTED);
//     m.exportPLY(outputPath+"41_flips.ply",ColMode::BY_VAL);
//     m.exportEdgesPLY(outputPath+"41_edges.ply");*/

//     //m.smoothDisputed();
//     //m.smoothDisputed();
//     //m.smoothDisputed();
//     //m.exportPLY(outputPath+"40_di1.ply",ColMode::BY_DISPUTED);


// }

// void test00(){
//     Mesh m;
//     Grid g;
//     m.importXYZv2(outputPath+"8.xyz");

//     //m.fuckUp();

//     m.delaunay();
//     //m.exportPLY(outputPath+"1_0Delonay.ply");

//     for (char c='A'; c<='D'; c++) {
//         m.greedyFlips( 1 );  // ( 10 , 30 );
//         m.exportPLY(outputPath+"1_1Flips"+c+".ply",ColMode::BY_VAL);
//         meshToGrid(m,g);
//         m.exportPLY(outputPath+"1_1Floodfill"+c+".ply", ColMode::BY_FLOOD);
//         g.exportPLY(outputPath+"1_1Grid"+c+".ply");
//         m.smoothDisputed();
//         m.smoothDisputed();
//     }

//     for (int i=0; i<1; i++) {
//         int nDone = 0;
//         //nDone += g.assignUnassignedNiceWay();
//         nDone += g.greedySwaps();
//         if (nDone) break;
//     }
//     g.exportPLY(outputPath+"1_3AssignSoft.ply");

//     for (int i=0; i<3; i++) {
//         int nDone = 0;
//         nDone += g.greedyAssignUnassigned();
//         nDone += g.greedySwaps();
//         if (nDone==0) break;
//     }
//     g.exportPLY(outputPath+"1_4AssignHard.ply");

//     g.tryAllSwapsBordersIncluded();
//     g.fillGapsMakingPtsUp();
//     g.greedySwaps();

//     g.fillGapsMakingPtsUp();

//     g.exportPLY(outputPath+"1_5Final.ply");
//     g.exportPLYtartan(outputPath+"1_5FinalT.ply");
// }

}
}

int main(){

    cellogram::PointsUntangler::test02();

}
