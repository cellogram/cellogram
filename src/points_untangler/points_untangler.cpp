#include "points_untangler.h"

#include<iostream>
#include"mesh.h"


namespace cellogram
{
    namespace PointsUntangler
    {

        namespace
        {
            std::string file_name(const std::string &outputPath, const char* label, int phase){
                static int lastPhase = 0;
                static int subPhase = 0;
                if (lastPhase!=phase) subPhase = 0;
                char str[8];
                sprintf(str,"%02d_%02d",phase , subPhase );
                subPhase++;
                lastPhase = phase;
                return outputPath + str + "_" + label;
            }
        }

        void pointsUntangler(const Eigen::MatrixXd &detected, Eigen::MatrixXi &tris, Eigen::MatrixXd &newPoints)
        {
            Mesh m;
            Grid g;
            m.fromEigen(detected);
            pointsUntangler(m, g);
            g.exportEigen(tris, newPoints);
        }

        void pointsUntangler(Mesh &m, Grid &g, const std::string &outputPath)
        {
            m.delaunay();
            m.buildEdgesFromFaces();
            m.propagateFixedF2E();
            m.updateValencies();

            if(!outputPath.empty())
                m.exportPLY( file_name(outputPath, "mesh",0) , ColMode::BY_VAL);

            for (int i=1; i<=5; i++) {

                m.greedyFlips(10, g);
                if(!outputPath.empty())
                    m.exportPLY( file_name(outputPath, "mesh_opt",i), ColMode::BY_VAL);

                meshToGrid(m,g);
                //if(!outputPath.empty())
                    //m.exportEdgesPLY(file_name(outputPath, "mesh_edges",i),g);
                g.computeMatrices();
                g.smoothMatrices(20);

                if(!outputPath.empty())
                {
                    m.exportPLY( file_name(outputPath, "mesh_grid",i), ColMode::BY_VAL);
                    g.exportPLY(file_name(outputPath, "grid_mesh",i) );
                }
                //return;

                {
                    if (i==5) g.greedyAssignUnassigned();
                    g.greedySwaps();
                    g.greedySwaps();
                    if (i==5) g.greedySwaps();
                    if (i==5) g.greedySwaps();

                    if(!outputPath.empty())
                        g.exportPLY(file_name(outputPath, "grid_opt",i) );
                }
                if(!outputPath.empty())
                    m.exportPLY( file_name(outputPath, "mesh_gopt",i), ColMode::BY_VAL);
                m.greedyFlips(10, g);
                m.flipAs(g);
                //m.greedyFlips(10);

                if(!outputPath.empty())
                    m.exportPLY( file_name(outputPath, "mesh_opt",i), ColMode::BY_VAL);

            }
        }

    }
}
