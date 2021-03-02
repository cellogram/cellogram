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


const auto func_const = [](auto x, auto y, auto z) -> double {return 1;};
const auto func_quadratic = [](auto x, auto y, auto z) -> double {
    if ((x-14.25)*(x-14.25) + (y-14.75)*(y-14.75) <= 4*4) {
        return 0;
    } else {
        return 1;
    }
};
const auto func_smooth = [](auto x, auto y, auto z) -> double {
    double dist = (x-14.25)*(x-14.25) + (y-14.75)*(y-14.75);
    return 1.0 / (1.0 + std::exp(-(dist-6)));
};
const auto func_cubic = [](auto x, auto y, auto z) -> double {
    // (x^2 + y^2)^(3/2)
    return pow((x-14.25)*(x-14.25) + (y-14.75)*(y-14.75), 1.5);
};


const auto GenImage = [](image_t &image, auto func) {
	int sizeX, sizeY, sizeZ, i;
    int x, y, z;

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

TEST_CASE("optim_ideal", "[OptimTest]") {

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

    // log
    const int N = CI.rows();
    int good=0, bad=0, failure=0;
    for (int i=0; i<N; i++) {
        if (IterO(i) >= optimPara.maxIt || EO(i) == 1.0) 
            failure++;
        else {
            if (std::pow(14.25 - CO(i, 0), 2.0) + std::pow(14.75 - CO(i, 1), 2.0) > 0.01*0.01) {  // xy precision require 0.01
                bad++;
                cout << "bad idx = " << i << endl;
            } else {
                good++;
            }
        }
    }

    cout << good << endl << bad << endl << failure << endl;

    Eigen::MatrixXd output(N, 10);
    output.leftCols(4) = CI;
    output.block(0, 4, N, 4) = CO;
    output.col(8) = EO;
    output.col(9) = IterO.cast<double>();

    cout << output << endl;

    REQUIRE(good >= 40);
    REQUIRE(bad <= 1);  // TODO: Debug the failure here idx=42
    REQUIRE(failure <= 9);  // due to discontinuity?
}


TEST_CASE("optim_smooth", "[OptimTest]") {

	image_t image; // 30 * 30 * 30
    GenImage(image, func_smooth);
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
                
                CI.row(cnt) << 14.25+i*gap, 14.5+j*gap, 10, 5.0;
                cnt++;
            }
    };
    GridCI(3, 1.0);

    optim::Optim_FixDepth(optimPara, bsplineSolver, CI, CO, EO, IterO, false);

    // log
    const int N = CI.rows();
    int good=0, bad=0, failure=0;
    for (int i=0; i<N; i++) {
        if (IterO(i) >= optimPara.maxIt || EO(i) == 1.0) 
            failure++;
        else {
            if (std::pow(14.2728 - CO(i, 0), 2.0) + std::pow(14.7272 - CO(i, 1), 2.0) > 0.01*0.01) {  // xy precision require 0.01
                bad++;
            } else {
                good++;
            }
        }
    }

    REQUIRE(good >= 45);
    REQUIRE(bad == 0);
    REQUIRE(failure <= 4);

    cout << good << endl << bad << endl << failure << endl;

    Eigen::MatrixXd output(N, 10);
    output.leftCols(4) = CI;
    output.block(0, 4, N, 4) = CO;
    output.col(8) = EO;
    output.col(9) = IterO.cast<double>();

    cout << output << endl;
}


TEST_CASE("optim_debug", "[OptimTest]") {

	image_t image; // 30 * 30 * 30
    GenImage(image, func_smooth);
    // GenImage(image, func_quadratic);
    double thres = QuantileImage(image, 1.0);
    NormalizeImage(image, thres);

    // prepare B-spline
    spdlog::set_level(spdlog::level::warn);
    const int bsplineDegree = 2;
    bspline bsplineSolver;
    bsplineSolver.CalcControlPts(image, 0.7, 0.7, 0.7, bsplineDegree);

    // TBB
    const int num_threads = 1;  // sequential
    const size_t MB = 1024 * 1024;
    const size_t stack_size = 64 * MB;
    tbb::task_scheduler_init scheduler(num_threads, stack_size);

    // optim
    OptimPara_t optimPara;
    Eigen::MatrixXd CI(1, 4), CO;
    Eigen::VectorXd EO;
    Eigen::VectorXi IterO;
    CI << 11.5, 15.5, 10, 5;

    optim::Optim_FixDepth(optimPara, bsplineSolver, CI, CO, EO, IterO, false);
}


/*
1:  12.494 15.4466 4.90499
   0.0213679 0.00302428    0.11407
2: 12.8994 15.4096 4.44327
  0.00448133 0.00669058   0.141446
3: 12.8994 15.4096 4.44328
  0.00448148 0.00669063   0.141446
4: 12.8994 15.4096 4.44328
  0.00448152 0.00669076   0.141445
5: 12.8994 15.4096 4.44328
  0.00448153 0.00669104   0.141445
6: 12.8994 15.4096 4.44329
  0.00448152 0.00669213   0.141444
7:  12.899 15.4097 4.44343
  0.00448117 0.00671875   0.141411
8: 12.8989 15.4098 4.44347
  0.00448126  0.0067241   0.141404
9: 12.8989 15.4098 4.44347
  0.00448134 0.00672512   0.141403
10: 12.8989 15.4098 4.44348
  0.00448143 0.00672562   0.141403
11: 12.8989 15.4098 4.44348
  0.00448154 0.00672614   0.141402
12: 12.8989 15.4098 4.44349
  0.0044817 0.0067271  0.141401
13: 12.8988 15.4099 4.44353
  0.00448138  0.0067341   0.141393
14: 12.8988 15.4099 4.44354
  0.00448168 0.00673524   0.141392
15: 12.8988 15.4099 4.44357
  0.00448219 0.00673712    0.14139
16: 12.8988   15.41 4.44362
   0.0044833 0.00674127   0.141386
17: 12.8987 15.4101 4.44376
   0.0044865 0.00675322   0.141375
18: 12.8982 15.4111 4.44472
  0.00450772 0.00683241   0.141301
19: -1.18781e+07   2.5141e+07  2.35813e+07  // why did this happen????
  0 0 0
*/