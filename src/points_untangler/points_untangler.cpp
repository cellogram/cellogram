#include "points_untangler.h"

#include<iostream>
#include <fstream>
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
                sprintf(str,"%02d_%c",phase , 'A' + subPhase );
                //sprintf(str,"%c_%02d" , 'A' + subPhase ,phase);
                subPhase++;
                lastPhase = phase;
                return outputPath + str + "_" + label;
            }
        }

        void pointsUntangler(const Eigen::MatrixXd &detected, Eigen::MatrixXi &tris, std::vector<int> &droppedPoints, Eigen::MatrixXd &newPoints)
        {
            bool verbose = true;
            Mesh m;
            Grid g;

            if(verbose){
                std::ofstream file("marco.txt");
                if(file.good())
                    file << detected <<std::endl;
                file.close();
            }

            m.fromEigen(detected);
            pointsUntangler(m, g, verbose ? "./" : "");
            g.fillGapsMakingPtsUp();
            if(verbose)
                g.exportOBJ("./output.obj");
            g.exportEigen(tris, droppedPoints, newPoints);
        }

        void pointsUntangler(Mesh &m, Grid &g, const std::string &outputPath)
        {
            bool verbose = !outputPath.empty();

            m.verbose = verbose;
            g.verbose = verbose;

            m.initWithDelaunay();
            if (verbose) m.exportPLY( file_name(outputPath, "mesh",0) , ColMode::BY_VAL);

            for (int i=1; i<=5; i++) {

                m.greedyFlips(10, g);

                if (verbose) m.exportPLY( file_name(outputPath, "mesh_opt",i), ColMode::BY_VAL);

                meshToGrid(m,g, verbose); // also modifies m

                if (verbose) m.exportPLY( file_name(outputPath, "mesh_grid",i), ColMode::BY_VAL);
                if (verbose) g.exportPLY(file_name(outputPath, "grid_mesh",i) );

                g.greedyOps();
                //m.greedyFlips( 10, g );

                if (verbose) g.exportPLY(file_name(outputPath, "grid_opt",i) );

                gridToMesh(g,m, verbose);

                m.updateValencies();
                if (verbose) m.exportPLY( file_name(outputPath, "mesh_gopt",i), ColMode::BY_VAL);

            }
        }

    }
}
