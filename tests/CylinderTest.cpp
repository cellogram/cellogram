// cylinder test
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
using namespace Catch::literals;


auto func_const = [](auto x, auto y, auto z) {return 1;};
auto func_ideal = [](auto x, auto y, auto z) {
    if ((x-14.5)*(x-14.5) + (y-14.5)*(y-14.5) <= 4*4) {
        return 0;
    } else {
        return 1;
    }
};
auto func_cubic = [](auto x, auto y, auto z) {
    return (x-14.5)*(x-14.5)*(x-14.5) + (y-14.5)*(y-14.5)*(y-14.5);
};

auto GenImage = [](image_t &image, auto func){
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

TEST_CASE("cylinder_const", "[CylinderTest]") {

	image_t image; // 30 * 30 * 30
    GenImage(image, func_const);

    // prepare B-spline
    spdlog::set_level(spdlog::level::warn);
    const int bsplineDegree = 2;
    bspline bsplineSolver;
    bsplineSolver.CalcControlPts(image, 0.7, 0.7, 0.7, bsplineDegree);

	// evaluate cylinder
    double x = 14.5;
    double y = 14.5;
    double z = 12;
    double r = 4.0;
    double res;
	bool valid = cylinder::IsValid(bsplineSolver, x, y, z, r, cylinder::H);
    cylinder::EvaluateCylinder(bsplineSolver, x, y, z, r, cylinder::H, res, false);

    REQUIRE(valid == true);
    REQUIRE(res == 0);
}


TEST_CASE("cylinder_ideal", "[CylinderTest]") {

	image_t image; // 30 * 30 * 30
    GenImage(image, func_ideal);

    // prepare B-spline
    spdlog::set_level(spdlog::level::warn);
    const int bsplineDegree = 2;
    bspline bsplineSolver;
    bsplineSolver.CalcControlPts(image, 0.7, 0.7, 0.7, bsplineDegree);

	// evaluate cylinder
    double x = 14.5;
    double y = 14.5;
    double z = 12;
    double r = 4.0;
    double res;
	bool valid = cylinder::IsValid(bsplineSolver, x, y, z, r, cylinder::H);
    cylinder::EvaluateCylinder(bsplineSolver, x, y, z, r, cylinder::H, res, false);

    
    REQUIRE(valid == true);
    REQUIRE(res == Approx(-0.89222).margin(0.0001));
    // cout << valid << endl;
    // cout << res << endl;
}


TEST_CASE("cylinder_cubic", "[CylinderTest]") {

	image_t image; // 30 * 30 * 30
    GenImage(image, func_cubic);

    // prepare B-spline
    spdlog::set_level(spdlog::level::warn);
    const int bsplineDegree = 2;
    bspline bsplineSolver;
    bsplineSolver.CalcControlPts(image, 0.7, 0.7, 0.7, bsplineDegree);

	// evaluate cylinder
    double x = 14.5;
    double y = 14.5;
    double z = 12;
    double r = 4.0;
    double res;
	bool valid = cylinder::IsValid(bsplineSolver, x, y, z, r, cylinder::H);
    cylinder::EvaluateCylinder(bsplineSolver, x, y, z, r, cylinder::H, res, false);

    
    REQUIRE(valid == true);
    REQUIRE(res == Approx(0).margin(0.0001));
    // cout << valid << endl;
    // cout << res << endl;
}
