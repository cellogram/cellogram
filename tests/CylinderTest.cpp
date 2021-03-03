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


auto func_const = [](auto x, auto y, auto z) -> double {return 1;};
auto func_ideal = [](auto x, auto y, auto z) -> double {
    if ((x-14.5)*(x-14.5) + (y-14.5)*(y-14.5) <= 4*4) {
        return 0;
    } else {
        return 1;
    }
};
auto func_cubic = [](auto x, auto y, auto z) -> double {
    return (x-14.5)*(x-14.5)*(x-14.5) + (y-14.5)*(y-14.5)*(y-14.5);
};
auto func_alpha = [](auto x, auto y, auto z) -> double {
    if ((x-14.5)*(x-14.5) + (y-14.5)*(y-14.5) <= 4*4) {
        return 0.2;
    } else {
        if ((x-14.5)*(x-14.5) + (y-14.5)*(y-14.5) <= 4*4*2)
            return 0.6;
        else 
            return 1.0;
    }
};


auto GenImage = [](image_t &image, auto func){
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


TEST_CASE("cylinder_invert", "[CylinderTest]") {

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


TEST_CASE("cylinder_alpha", "[CylinderTest]") {

	image_t image; // 30 * 30 * 30
    GenImage(image, func_alpha);

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

    cylinder::alpha = 0.5;
	bool valid = cylinder::IsValid(bsplineSolver, x, y, z, r, cylinder::H);
    cylinder::EvaluateCylinder(bsplineSolver, x, y, z, r, cylinder::H, res, false);
    REQUIRE(valid == true);
    REQUIRE(res == Approx(-0.4).margin(0.03));

    cylinder::alpha = 0.2;
	valid = cylinder::IsValid(bsplineSolver, x, y, z, r, cylinder::H);
    cylinder::EvaluateCylinder(bsplineSolver, x, y, z, r, cylinder::H, res, false);
    REQUIRE(valid == true);
    REQUIRE(res == Approx(-0.52).margin(0.02));
    cout << "alpha=0.2 " << res << endl;

    cylinder::alpha = 0.7;
	valid = cylinder::IsValid(bsplineSolver, x, y, z, r, cylinder::H);
    cylinder::EvaluateCylinder(bsplineSolver, x, y, z, r, cylinder::H, res, false);
    REQUIRE(valid == true);
    REQUIRE(res == Approx(-0.32).margin(0.05));  // huge error caused by interp
    cout << "alpha=0.7 " << res << endl;

    cylinder::alpha = 0.5;
}
