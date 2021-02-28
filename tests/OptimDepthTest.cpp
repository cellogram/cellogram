// optim z test
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

const auto func_const = [](auto x, auto y, auto z) -> double { return 1; };
const auto func_quadratic = [](auto x, auto y, auto z) -> double {
    if ((x - 14.25) * (x - 14.25) + (y - 14.75) * (y - 14.75) <= 4 * 4) {
        return 0;
    } else {
        return 1;
    }
};
const auto func_smooth = [](auto x, auto y, auto z) -> double {
    double dist = (x - 14.25) * (x - 14.25) + (y - 14.75) * (y - 14.75);
    return 1.0 / (1.0 + std::exp(-(dist - 6)));
};
const auto func_cubic = [](auto x, auto y, auto z) -> double {
    // (x^2 + y^2)^(3/2)
    return pow((x - 14.25) * (x - 14.25) + (y - 14.75) * (y - 14.75), 1.5);
};

const auto GenImage = [](image_t &image, auto func) {
    int sizeX, sizeY, sizeZ, i;
    int x, y, z;

    sizeX = 30; // 0, 1, ..., 29
    sizeY = 30;
    sizeZ = 30;

    // generate sample grid (3D)
    double maxPixel = 0;
    for (z = 0; z < sizeZ; z++) {

        MatrixXd layer(sizeX, sizeY);
        for (x = 0; x < sizeX; x++)
            for (y = 0; y < sizeY; y++) {
                layer(x, y) = func(x, y, z);
            }

        image.push_back(layer);
        if (layer.maxCoeff() > maxPixel)
            maxPixel = layer.maxCoeff();
    }
};

////////////////////////////////////////////////////////////////////////////

TEST_CASE("optimz_ideal", "[OptimTest]") {

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
    Eigen::MatrixXd CI(3, 4), CO;
    CI << 14, 14, 10, 5, 
          15, 15, 10, 5, 
          16, 15, 10, 5;
    optim::Optim_WithDepth(optimPara, bsplineSolver, 3, 1.0, CI, CO);

    /*
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
						} else {
								good++;
						}
				}
		}

		REQUIRE(good >= 40);
		REQUIRE(bad == 0);
		REQUIRE(failure <= 9);  // due to discontinuity?

		cout << good << endl << bad << endl << failure << endl;

		Eigen::MatrixXd output(N, 10);
		output.leftCols(4) = CI;
		output.block(0, 4, N, 4) = CO;
		output.col(8) = EO;
		output.col(9) = IterO.cast<double>();

		cout << output << endl;
		*/
}

TEST_CASE("optimz_debug", "[OptimTest]") {

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
    const int num_threads = 1; // sequential
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
