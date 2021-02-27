// optim test
#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/Optimization.h>
#include <zebrafish/Quantile.h>
#include <zebrafish/Logger.hpp>

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>
#include <catch.hpp>

#include <math.h>
#include <stdlib.h>
#include <random>
#include <queue>
#include <vector>

using namespace std;
using namespace Eigen;
using namespace zebrafish;
using namespace Catch::literals;


const auto func_const = [](auto x, auto y, auto z) {return 1;};
const auto func_quadratic = [](auto x, auto y, auto z) {
    if ((x-14.5)*(x-14.5) + (y-14.5)*(y-14.5) <= 4*4) {
        return 0;
    } else {
        return 1;
    }
};


const auto GenImage = [](image_t &image, auto func) {
	int sizeX, sizeY, sizeZ, i;
    double x, y, z;

    sizeX = 30; // 0, 1, ..., 29
    sizeY = 30;
    sizeZ = 30;

    // generate sample grid (3D)
    double maxPixel = 0;
    for (z=0; z<sizeZ; z++) {
        
        MatrixXd layer(sizeX, sizeY);
        for (x=0; x<sizeX; x++)
            for (y=0; y<sizeY; y++) {
                layer(x, y) = func(x, y, z);
            }

        image.push_back(layer);
        if (layer.maxCoeff() > maxPixel) maxPixel = layer.maxCoeff();
    }
};

////////////////////////////////////////////////////////////////////////////

TEST_CASE("optim_zone", "[OptimTest]") {

	image_t image; // 30 * 30 * 30
    GenImage(image, func_quadratic);
    double thres = QuantileImage(image, 1.0);
    NormalizeImage(image, thres);

    // prepare B-spline
    spdlog::set_level(spdlog::level::warn);
    const int bsplineDegree = 2;
    bspline bsplineSolver;
    bsplineSolver.CalcControlPts(image, 0.7, 0.7, 0.7, bsplineDegree);

    // TBB
    const int num_threads = 3;
    const size_t MB = 1024 * 1024;
    const size_t stack_size = 64 * MB;
    tbb::task_scheduler_init scheduler(num_threads, stack_size);

    // optim
    OptimPara_t optimPara;
    Eigen::MatrixXd CI, CO;
    Eigen::VectorXd EO;
    Eigen::VectorXi IterO;
    const auto GridCI = [&CI](int m, double gap) {
        CI.resize((2*m+1)*(2*m+1), 4);
        int cnt = 0;
        for (int i=-m; i<=m; i++)
            for (int j=-m; j<=m; j++) {
                
                CI.row(cnt) << 14.5+i*gap, 14.5+j*gap, 10, 5.0;
                cnt++;
            }
    };
    GridCI(3, 1.0);

    optim::Optim_FixDepth(optimPara, bsplineSolver, CI, CO, EO, IterO, false);

    const int N = CI.rows();
    Eigen::MatrixXd output(N, 10);
    output.leftCols(4) = CI;
    output.block(0, 4, N, 4) = CO;
    output.col(8) = EO;
    output.col(9) = IterO.cast<double>();

    cout << output << endl;
}
