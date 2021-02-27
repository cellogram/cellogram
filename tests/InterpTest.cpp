// Interp TEST
#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/Logger.hpp>

#include <catch.hpp>

#include <math.h>
#include <stdlib.h>
#include <random>
#include <queue>
#include <vector>

using namespace std;
using namespace Eigen;
using namespace zebrafish;

double func_linear(double x, double y, double z) {
    return x + y;
}

double func_quad(double x, double y, double z) {

    // x^2 + y^2
    return (x-14.5)*(x-14.5) + (y-14.5)*(y-14.5);

    // (x^2 + y^2)^(3/2)
    // return pow((x-14.5)*(x-14.5) + (y-14.5)*(y-14.5), 1.5);
    // return (x-14.5)*(x-14.5)*(x-14.5) + (y-14.5)*(y-14.5)*(y-14.5);

    // x^4 + y^4 + 2 * x^2 * y^2
    // return (x-14.5)*(x-14.5)*(x-14.5)*(x-14.5) + (y-14.5)*(y-14.5)*(y-14.5)*(y-14.5) +
    //     2 * (y-14.5)*(y-14.5)* (x-14.5)*(x-14.5);

    // x^5 + y^5
    // return (x-14.5)*(x-14.5)*(x-14.5)*(x-14.5)*(x-14.5) + (y-14.5)*(y-14.5)*(y-14.5)*(y-14.5)*(y-14.5);

    // (x^2 + y^2)^(5/2)
    // return pow((x-14.5)*(x-14.5) + (y-14.5)*(y-14.5), 2.5);
}

template<typename Func>
void interpolate(const Func &func, int &degree, double &mean_error, double &median_error, double &min_error, double &max_error) {

    image_t image; // 30 * 30 * 10
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

    // prepare B-spline
    const int bsplineDegree = degree;
    spdlog::set_level(spdlog::level::warn);
    bspline bsplineSolver;
    bsplineSolver.SetResolution(0.325, 0.325, 0.5);
    bsplineSolver.CalcControlPts_um(image, 0.6, 0.6, 0.6, bsplineDegree);

    // random
    srand(time(NULL));
    uniform_real_distribution<double> unif(0, sizeX - 1);
    default_random_engine re(0);

    // moving median
    double tt, l_;
    priority_queue<double, vector<double>, greater<double>> u;
    priority_queue<double, vector<double>, less<double>> l;

    // interp test
    double err, sumerr = 0, minerr = 1.0, maxerr = 0.0;
    const int trialNum = 1e6;
    Matrix<double, Dynamic, 2> sampleInput;
    Matrix<double, Dynamic, 1> sampleOutput;
    sampleInput.resize(trialNum, 2);
    sampleOutput.resize(trialNum, 1);
    for (i = 0; i<trialNum; i++) {

        x = unif(re);
        y = unif(re);
        sampleInput(i, 0) = x;
        sampleInput(i, 1) = y;
    }
    z = 5;

    // Interp 
    for (i = 0; i<trialNum; i++) {
        
        sampleOutput(i) = bsplineSolver.Interp3D(sampleInput(i, 0), sampleInput(i, 1), z);
    }

    // calculate theoretical output
    for (i = 0; i<trialNum; i++) {

        err = func(sampleInput(i, 0), sampleInput(i, 1), z) - sampleOutput(i);
        err = fabs(err);

        // stats
        sumerr += err;
        if (err > maxerr) maxerr = err;
        if (err < minerr) minerr = err;

        // moving median
        if (l.empty()) {l.push(err); continue;}
        l_ = l.top();
        if (err >= l_) {u.push(err);} else {l.push(err);}
        if (u.size() > l.size()) {
            tt = u.top();
            u.pop();
            l.push(tt);
        }
        if (l.size() > u.size() + 1) {
            tt = l.top();
            l.pop();
            u.push(tt);
        }
    }
}

/////////////////////////////////////////////////////

TEST_CASE("linear_interp", "[InterpTest]") {

    int degree = 2;
    double mean_error;
    double median_error;
    double min_error;
    double max_error;

    interpolate(func_linear, degree, mean_error, median_error, min_error, max_error);

    REQUIRE(mean_error < 1e-10);
    REQUIRE(median_error < 1e-10);
    REQUIRE(min_error < 1e-10);
    REQUIRE(max_error < 1e-10);
    /*
    cout << "Degree = " << degree << endl;
    cout << "Mean error = " << mean_error << endl;
    cout << "Median error = " << median_error << endl;
    cout << "Min  error = " << min_error << endl;
    cout << "Max  error = " << max_error << endl;
    */
}


TEST_CASE("quadra_interp", "[InterpTest]") {

    int degree = 2;
    double mean_error;
    double median_error;
    double min_error;
    double max_error;

    interpolate(func_quad, degree, mean_error, median_error, min_error, max_error);

    REQUIRE(mean_error < 1e-10);
    REQUIRE(median_error < 1e-10);
    REQUIRE(min_error < 1e-10);
    REQUIRE(max_error < 1e-10);
    /*
    cout << "Degree = " << degree << endl;
    cout << "Mean error = " << mean_error << endl;
    cout << "Median error = " << median_error << endl;
    cout << "Min  error = " << min_error << endl;
    cout << "Max  error = " << max_error << endl;
    */
}
