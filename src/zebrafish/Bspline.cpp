#include <zebrafish/Bspline.h>
#include <zebrafish/Common.h>
#include <zebrafish/autodiff.h>
#include <zebrafish/Logger.hpp>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include <Eigen/IterativeLinearSolvers>
#include <polysolve/LinearSolver.hpp>
#include <string>
#include <cmath>
#include <math.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <array>


namespace zebrafish {

namespace {

// Get the value the variable
template<typename T>
class GetVal {
public:
    double operator()(const T x) const {
        return x.getValue();
    }
};

template<>
class GetVal<double> {
public:
    double operator()(const double x) const {
        return x;
    }
};

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////////

void bspline::SetResolution(const double resX, const double resY, const double resZ) {

    resolutionX = resX;
    resolutionY = resY;
    resolutionZ = resZ;
}


void bspline::CalcControlPts_um(const image_t &image, const double distX, const double distY, const double distZ, const int degree_) {

    assert(degree_ == 2 || degree_ == 3);

    double xratio, yratio, zratio;
    // dimension of sample points
    Nz = image.size();
    Nx = image[0].rows();
    Ny = image[0].cols();
    /// Note: Why the calculation of ratio for degree 2 & 3 B-spline are different?
    ///       Clamped B-spline => different #multiplicity
    if (degree_ == 3) {
        numX = round(resolutionX / distX * (double(Nx) - 1) + 3.0);
        numY = round(resolutionY / distY * (double(Ny) - 1) + 3.0);
        numZ = round(resolutionZ / distZ * (double(Nz) - 1) + 3.0);
    } else {
        numX = round(resolutionX / distX * (double(Nx) - 1) + 2.0);
        numY = round(resolutionY / distY * (double(Ny) - 1) + 2.0);
        numZ = round(resolutionZ / distZ * (double(Nz) - 1) + 2.0);
    }
    xratio = (double(numX) - 0.5) / double(Nx);
    yratio = (double(numY) - 0.5) / double(Ny);
    zratio = (double(numZ) - 0.5) / double(Nz);

    assert(xratio > 0 && yratio > 0 && zratio > 0);
    assert(xratio <= 1 && yratio <= 1 && zratio <= 1);

    CalcControlPts(image, xratio, yratio, zratio, degree_);
}


void bspline::CalcControlPts(const image_t &image, const double xratio, const double yratio, const double zratio, const int degree_) {
// Note: this function will only be called once for each 3D image

    assert(xratio <= 1 && yratio <= 1 && zratio <= 1);
    assert(degree_ == 2 || degree_ == 3);

    // register B-spline degree to class member variable
    degree = degree_;
        logger().info("====================================================");
        logger().info("B-Spline degree = {}", degree);
        logger().info("Start calculating control points...");
        // Warning if ratio too small
        /// A small control point ratio will introduce error to interpolation 
        /// and energy evaluation.
        const double ratioThres = 0.5;
        if (xratio < ratioThres || yratio < ratioThres || zratio < ratioThres)
            logger().warn("Control point size too small: xratio={}, yratio={} and zratio={}.", xratio, yratio, zratio);

    int N, num;
    double centerX, centerY, centerZ, t;
    Eigen::VectorXd inputPts;

    // dimension of sample points
    Nz = image.size();
    Nx = image[0].rows();
    Ny = image[0].cols();
    N = Nx*Ny*Nz;
        logger().debug("Sample points:  N= {} Nx= {} Ny= {} Nz= {}", N, Nx, Ny, Nz);

    // dimension of control points
    numX = ceil(Nx * xratio);
    numY = ceil(Ny * yratio);
    numZ = ceil(Nz * zratio);
    num = numX * numY * numZ;
        logger().debug("Control points: num= {} numX= {} numY= {} numZ={}", num, numX, numY, numZ);

    // A[i, j]: the evaluation of basis j for the i-th sample point
    Eigen::SparseMatrix<double, Eigen::RowMajor> A(N, num);
    Eigen::SparseMatrix<double, Eigen::ColMajor> Atranspose(num, N);
    // Reserve space for least square matrix
    // At most [ (degree+1)^3 ] non-zero elements per row
    // This step is critical for efficiency
    A.reserve(Eigen::VectorXi::Constant(A.rows(), (degree==2 ? 27 : 64)));
    Atranspose.reserve(Eigen::VectorXi::Constant(Atranspose.cols(), (degree==2 ? 27 : 64)));

    // the interval between two control points along a direction
    if (degree == 3) {
        gapX = double(Nx-1) / double(numX-1-2);  // "-2" due to clamped B-spline
        gapY = double(Ny-1) / double(numY-1-2);
        gapZ = double(Nz-1) / double(numZ-1-2);
    } else if (degree == 2) {
        gapX = double(Nx-1) / double(numX-1-1);  // "-1" due to clamped B-spline
        gapY = double(Ny-1) / double(numY-1-1);
        gapZ = double(Nz-1) / double(numZ-1-1);
    }
        logger().debug("xratio = {}  yratio = {}  zratio = {}", xratio, yratio, zratio);
        logger().debug("gapX = {}px  gapY = {}px  gapZ = {}px", gapX, gapY, gapZ);
        logger().debug("gapX = {}um  gapY = {}um  gapZ = {}um", gapX*resolutionX, gapY*resolutionY, gapZ*resolutionZ);

    // map 3D "image" to 1D "inputPts"
        // order matters! z -> y -> x
    inputPts.resize(N);
    t = 0;
    for (auto it=image.begin(); it!=image.end(); it++) {
        inputPts.segment(t*Nx*Ny, Nx*Ny) << Eigen::Map<const Eigen::VectorXd>(it->data(), Nx*Ny);
        t++;
    }

    // Calculate basis function lookup table
    // DScalar version
    CalcBasisFunc<DScalar>(basisX, numX, gapX);
    CalcBasisFunc<DScalar>(basisY, numY, gapY);
      /// CalcBasisFunc<DScalar>(basisZ, numZ, gapZ);
      /// This is useless because depth-axis will not be affected by [x, y, r] and thus has gradient zero
    // double version
    CalcBasisFunc<double>(basisXd, numX, gapX);
    CalcBasisFunc<double>(basisYd, numY, gapY);
    CalcBasisFunc<double>(basisZd, numZ, gapZ);

    // Calculate & fill the least square matrix A
        logger().info("Calculating & filling least square matrix...");
    CalcLeastSquareMat(A, Atranspose);
        logger().info("Least square matrix A size = {} * {}", A.rows(), A.cols());
        logger().info("A: # non-zero elements = {}/{}", A.nonZeros(), A.rows()*A.cols());

    // Solve the least square problem
    // controlPts = inverse(A'*A) * (A'*inputPts)
    SolveLeastSquare(A, Atranspose, inputPts);

        /////////////// DEBUG ONLY ///////////
        // std::cout << ">>>>> control points >>>>>" << std::endl;
        // for (int i=0; i<controlPoints.size(); i++)
        //     printf("%13.10f\n", controlPoints(i));
        // printf("\n");
        // // std::cout << controlPoints << std::endl;
        // std::cout << "<<<<< control points <<<<<" << std::endl;
        // std::cout << ">>>>> input points >>>>>" << std::endl;
        // std::cout << inputPts << std::endl;
        // std::cout << "<<<<< input points <<<<<" << std::endl;
        // std::cout << ">>>>> A >>>>>" << std::endl;
        // for (int i=0; i<A.rows(); i++) {
        //     for (int j=0; j<A.cols(); j++) {
        //         printf("%13.10f  ", A.coeff(i, j));
        //     }
        //     printf("\n");
        // }
        // std::cout << "<<<<< A <<<<<" << std::endl;
        /////////////// DEBUG ONLY ///////////

    if (controlPointsCacheFlag) {
        // optional optimization. By default disabled.
        logger().info("Creating control points cache...");
        CreateControlPtsCache();
    }
        logger().info("Control points calculated.");
        logger().info("====================================================");
}


void bspline::CalcLeastSquareMat(Eigen::SparseMatrix<double, Eigen::RowMajor> &A, Eigen::SparseMatrix<double, Eigen::ColMajor> &Atranspose) {

    int i, j, ix, iy, iz, jx, jy, jz, refIdx_x, refIdx_y, refIdx_z;
    int jxMin, jxMax, jyMin, jyMax, jzMin, jzMax;
    double t, t1, t2, t3;

    // order matters!
    // for column (j): z -> x -> y
    // for row (i): z -> y -> x
    for (iz=0; iz<=Nz-1; iz++)
        for (iy=0; iy<=Ny-1; iy++)
            for (ix=0; ix<=Nx-1; ix++) {

                // iterating sample points
                i = iz * Nx * Ny + iy * Nx + ix;

                // progress bar
                // if (i % int(Nx*Ny*Nz*0.2) == 0)
                //    logger().trace("{} / {}", i, Nx*Ny*Nz);

                refIdx_x = floor(ix / gapX);
                refIdx_y = floor(iy / gapY);
                refIdx_z = floor(iz / gapZ);
                // if the query point lies exactly at one edge of the input grid
                if (refIdx_x == numX-1-(degree-1)) refIdx_x--;
                if (refIdx_y == numY-1-(degree-1)) refIdx_y--;
                if (refIdx_z == numZ-1-(degree-1)) refIdx_z--;

                for (jz=refIdx_z; jz<=refIdx_z+degree; jz++)
                    for (jx=refIdx_x; jx<=refIdx_x+degree; jx++)
                        for (jy=refIdx_y; jy<=refIdx_y+degree; jy++) {

                            // iterating control points' basis function
                            j = jz * numX * numY + jx * numY + jy;
                            t1 = basisXd(jx-refIdx_x, refIdx_x)(ix);
                            t2 = basisYd(jy-refIdx_y, refIdx_y)(iy);
                            t3 = basisZd(jz-refIdx_z, refIdx_z)(iz);
                            t = t1 * t2 * t3;

                            if (fabs(t) > 2e-15) {
                                A.insert(i, j) = t;  // O(1) insertion (space has been reserved)
                                Atranspose.insert(j, i) = t;
                            }
                        }
                ///////////// DEBUG ONLY
                // if (fabs(rowSum-1.0)> 1e-8)  // if (rowSum != 1.0)
                //     logger().error("RowSum error: ix={}, iy={}, iz={}, rowSum={}", ix, iy, iz, rowSum);
                ///////////// DEBUG ONLY
            }
}


void bspline::SolveLeastSquare(const Eigen::SparseMatrix<double, Eigen::RowMajor> &A, const Eigen::SparseMatrix<double, Eigen::ColMajor> &Atranspose, const Eigen::VectorXd &inputPts) {

    const int num = numX * numY * numZ;
    const int N = Nx * Ny * Nz;
    Eigen::VectorXd vectorY;
    Eigen::SparseMatrix<double, Eigen::ColMajor> AtransposeA(num, num);

    switch (leastSquareMethod) {
    case 1: {
        logger().info("Least square using method 1: Normal Function + Hypre Solver");
        logger().error("This method has been disabled");
        // calculate control points based on least square
            logger().info("Calculating A'*A and A'*y...");
        AtransposeA = (Atranspose * A).pruned();  // A' * A
        vectorY.resize(N);
        vectorY = Atranspose * inputPts;  // A' * y
            logger().debug("A' * A size = {} * {}", AtransposeA.rows(), AtransposeA.cols());
            logger().debug("A'A: # non-zero elements = {}", AtransposeA.nonZeros());

        // solve linear system for control points
            // Deprecated eigen solver:
            // Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > chol(AtransposeA);
            // controlPoints = chol.solve(vectorY);
        const std::string solverName = "Hypre";
        auto solver = polysolve::LinearSolver::create(solverName, "");
        const nlohmann::json params = {
            {"max_iter", solverMaxIt},
            {"conv_tol", solverConvTol},
            {"tolerance", solverTol}
        };
        solver->setParameters(params);
            logger().info("Analyzing matrix pattern...");
        solver->analyzePattern(AtransposeA, AtransposeA.rows());
            logger().info("Factorizing matrix...");
        solver->factorize(AtransposeA);
            logger().info("Solving linear system...");
        controlPoints.resize(num, 1);
        solver->solve(vectorY, controlPoints);
        break;
    }
    case 2: {
            logger().info("Least square using method 2: Conjugate Gradient");
        Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > lscg;
        lscg.setTolerance(solverTol); 
            logger().info("LS conjugate gradient computing...");
        lscg.compute(A);

            logger().info("LSCG Estimating initial guess...");
        Eigen::VectorXd initialGuess;
        CalcInitialGuess(inputPts, initialGuess);
        
            logger().info("Solving linear system...");
        controlPoints.resize(num, 1);
        controlPoints = lscg.solveWithGuess(inputPts, initialGuess);

        // log result
        initialGuess = (initialGuess - controlPoints).array().abs();
            logger().debug("guess mean_err = {}", initialGuess.mean());
            logger().debug("guess max_err = {}", initialGuess.maxCoeff());
            logger().info("LSCG #iteration = {}", lscg.iterations());
            logger().info("LSCG error = {}", lscg.error());
        break;
    }
    case 3: {
        assert(false);
        break;
        /*
        // Deprecated QR least square solver
        // This is too slow
        logger().info("Prepare solver");
        Eigen::SparseQR <Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
        A.makeCompressed();
        logger().info("A compressed");
        solver.compute(A);
        logger().info("A factorized");
        controlPoints = solver.solve(inputPts);
        */
    }
    default: {
        assert(false);
        break;
    }
    }
}


void bspline::CalcInitialGuess(const Eigen::VectorXd &inputPts, Eigen::VectorXd &initialGuess) {
// This function calculates the estimated "controlPts" based on "inputPts"
// The guess is used as the starting point of LS conjugate gradient

    initialGuess.resize(numX*numY*numZ, 1);

    int i, j, idx_x, idx_y, idx_z;
    const double delta = (degree == 2) ? 0.5 : 2.0;  // why delta? clamped B-spline

    for (int iz=0; iz<numZ; iz++)
        for (int ix=0; ix<numX; ix++)
            for (int iy=0; iy<numY; iy++) {

                i = iz * numX * numY + ix * numY + iy;
                
                idx_x = round((ix-delta) * gapX);
                idx_y = round((iy-delta) * gapY);
                idx_z = round((iz-delta) * gapZ);
                if (idx_x < 0) idx_x = 0;
                if (idx_y < 0) idx_y = 0;
                if (idx_z < 0) idx_z = 0;
                if (idx_x > Nx-1) idx_x = Nx-1;
                if (idx_y > Ny-1) idx_y = Ny-1;
                if (idx_z > Nz-1) idx_z = Nz-1;

                j = idx_z * Nx * Ny + idx_y * Nx + idx_x;

                initialGuess(i) = inputPts(j);
            }
    
    // TEST only 
    // initialGuess = Eigen::VectorXd::Zero(numX*numY*numZ);
}


template <typename T>
void bspline::CalcBasisFunc(Eigen::Matrix< std::function<T (T) >, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> &basisT, const int numT, const double gapT) {
// NOTE: "basisT" can be "basisX", "basisY" or "basisZ"

    if (degree == 3) {  // cubic basis

        assert(numT >= 8);

        basisT.resize(4, numT-1-2);
        // first 3 columns
            // cubic [11112]
        basisT(0, 0) = [gapT](T x) {x = x/gapT; return - (x-1.0)*(x-1.0)*(x-1.0);};
            // cubic [11123]
        basisT(1, 0) = [gapT](T x) {x = x/gapT; return (7.0/4.0)*x*x*x - (9.0/2.0)*x*x + 3.0*x;};
        basisT(0, 1) = [gapT](T x) {x = x/gapT; return -(1.0/4.0)*x*x*x + (3.0/2.0)*x*x - 3.0*x + 2.0;};
            // cubic [11234]
        basisT(2, 0) = [gapT](T x) {x = x/gapT; return -(11.0/12.0)*x*x*x + (3.0/2.0)*x*x;};
        basisT(1, 1) = [gapT](T x) {x = x/gapT; return (7.0/12.0)*x*x*x - 3.0*x*x + 9.0/2.0*x - (3.0/2.0);};
        basisT(0, 2) = [gapT](T x) {x = x/gapT; return -(1.0/6.0)*x*x*x + (3.0/2.0)*x*x - (9.0/2.0)*x + (9.0/2.0);};
        // last 3 columns
            // cubic [11112]
        basisT(3, numT-4) = [gapT, numT](T x) {x = numT-3-x/gapT; return - (x-1.0)*(x-1.0)*(x-1.0);};
            // cubic [11123]
        basisT(2, numT-4) = [gapT, numT](T x) {x = numT-3-x/gapT; return (7.0/4.0)*x*x*x - (9.0/2.0)*x*x + 3.0*x;};
        basisT(3, numT-5) = [gapT, numT](T x) {x = numT-3-x/gapT; return -(1.0/4.0)*x*x*x + (3.0/2.0)*x*x - 3.0*x + 2.0;};
            // cubic [11234]
        basisT(1, numT-4) = [gapT, numT](T x) {x = numT-3-x/gapT; return -(11.0/12.0)*x*x*x + (3.0/2.0)*x*x;};
        basisT(2, numT-5) = [gapT, numT](T x) {x = numT-3-x/gapT; return (7.0/12.0)*x*x*x - 3.0*x*x + 9.0/2.0*x - (3.0/2.0);};
        basisT(3, numT-6) = [gapT, numT](T x) {x = numT-3-x/gapT; return -(1.0/6.0)*x*x*x + (3.0/2.0)*x*x - (9.0/2.0)*x + (9.0/2.0);};
        // middle columns
            // cubic [12345]
        for (int i=2; i<=numT-5; i++) {
            basisT(3, i-2) = [gapT, i](T x) {x = x/gapT-i; return (1.0/6.0)*x*x*x + x*x + 2.0*x + (4.0/3.0);};
            basisT(2, i-1) = [gapT, i](T x) {x = x/gapT-i; return -(1.0/2.0)*x*x*x - x*x + (2.0/3.0);};
            basisT(1, i  ) = [gapT, i](T x) {x = x/gapT-i; return (1.0/2.0)*x*x*x - x*x + (2.0/3.0);};
            basisT(0, i+1) = [gapT, i](T x) {x = x/gapT-i; return -(1.0/6.0)*x*x*x + x*x - 2.0*x + (4.0/3.0);};
        }
    } else if (degree == 2) {  // quadratic basis

        assert(numT >= 5);

        basisT.resize(3, numT-1-1);
        // first 2 columns
            // quadratic [1112]
        basisT(0, 0) = [gapT](T x) {x = x/gapT; return (x-1.0)*(x-1.0);};
            // quadratic [1123]
        // Why minus one? the basis function is from -1, see notes
        basisT(1, 0) = [gapT](T x) {x = x/gapT-1; return -(3.0/2.0)*x*x - x + (1.0/2.0);};
        basisT(0, 1) = [gapT](T x) {x = x/gapT-1; return (1.0/2.0)*x*x - x + (1.0/2.0);};
        // last 2 columns
            // quadratic [1112]
        basisT(2, numT-3) = [numT, gapT](T x) {x = numT-2-x/gapT; return (x-1.0)*(x-1.0);};
            // quadratic [1123]
        basisT(1, numT-3) = [numT, gapT](T x) {x = numT-2-x/gapT-1; return -(3.0/2.0)*x*x - x + (1.0/2.0);};
        basisT(2, numT-4) = [numT, gapT](T x) {x = numT-2-x/gapT-1; return (1.0/2.0)*x*x - x + (1.0/2.0);};
        // middle columns
            // quadratic [1234]
        for (int i=1; i<=numT-4; i++) {
            basisT(2, i-1) = [gapT, i](T x) {x = x/gapT-i; return (1.0/2.0)*x*x + x + (1.0/2.0);};
            basisT(1, i  ) = [gapT, i](T x) {x = x/gapT-i; return -x*x + x + (1.0/2.0);};
            basisT(0, i+1) = [gapT, i](T x) {x = x/gapT-i; return (1.0/2.0)*x*x - 2.0*x + 2.0;};
        }
    }
}
// explicit template instantiation
template void bspline::CalcBasisFunc(Eigen::Matrix< std::function<DScalar (DScalar) >, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> &basisT, const int numT, const double gapT);
template void bspline::CalcBasisFunc(Eigen::Matrix< std::function<double  (double)  >, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> &basisT, const int numT, const double gapT);


void bspline::CreateControlPtsCache() {
// Optional optimization. Uses significant amount of memory.
// By default disabled.

    assert(degree == 2);

    const int N = numX * numY * numZ;
    int ix, iy, iz, i, jx, jy, jz, j, idx;
    logger().debug("Cached control point size: {}Mb = N * 27.0 * 8 / 1e6", N * 27.0 * 8 / 1e6);
    controlPointsCache.resize(27, N);

    for (iz=0; iz<numZ-2; iz++)
        for (ix=0; ix<numX-2; ix++)
            for (iy=0; iy<numY-2; iy++) {

                i = iz * numX * numY + ix * numY + iy;

                for (jz=0; jz<=2; jz++)
                    for (jx=0; jx<=2; jx++)
                        for (jy=0; jy<=2; jy++) {

                            j = jz * 9 + jx * 3 + jy;
                            idx = i + jz*numX*numY + jx*numY + jy;

                            controlPointsCache(j, i) = controlPoints(idx);
                        }
            }
}


////////////////////////////////////////////////////////////////////////////////////////
// Interp3D & InterpDisk


template <typename T>
T bspline::Interp3D(const T &x, const T &y, const T &z) const {
// This function is only used for test/debug purpose. Do NOT optimize this.
// Warning: This function can only used for test/debug. It is wrong to use
// Interp3D before InterpDisk.

    quad.LoadDiskParas(-2);

    Eigen::Matrix<T, Eigen::Dynamic, 1> resArr;
    const GetVal<T> getVal;

    InterpDisk(x, y, getVal(z), T(1.0), resArr);  // dummy radius
    return resArr(0);
}
// explicit template instantiation
template DScalar bspline::Interp3D(const DScalar &x, const DScalar &y, const DScalar &z) const;
template double  bspline::Interp3D(const double  &x, const double  &y, const double  &z) const;


template <typename basisT, typename T>
void bspline::InterpDiskHelper(const T x, const T y, const double z, const T r, Eigen::Matrix<T, Eigen::Dynamic, 1> &res, 
                             const basisT &basisX_, const basisT &basisY_, const basisd_t &basisZ_) const {

    assert(degree == 2 || degree == 3);
    assert(numX != 0 && numY != 0 && numZ != 0);

    static const GetVal<T> getVal;
    int i, ix, iy, iz, idx;
    std::array<T, 4> basisX_v, basisY_v, basisZ_v;
    int refIdx_x, refIdx_y, refIdx_z = floor(z / gapZ);
    T xT, yT;

    res.resize(quad.numDiskPts, 1);
    for (i=0; i<quad.numDiskPts; i++) {

        // Calculate the transformed location
        xT = quad.xyArray(i, 0) * r + x, 
        yT = quad.xyArray(i, 1) * r + y;
        // Get reference index
        refIdx_x = floor( getVal(xT) / gapX);
        refIdx_y = floor( getVal(yT) / gapY);

        /// A special case: if the query point lies exactly at the end of the edge
        /// NOTE: in practice, we will never query those points
        /// Uncomment the following three lines to fix this at the cost of ~5% more time
        // if (refIdx_x == numX-1-(degree-1)) refIdx_x--;
        // if (refIdx_y == numY-1-(degree-1)) refIdx_y--;
        // if (refIdx_z == numZ-1-(degree-1)) refIdx_z--;

        assert(refIdx_x >= 0 && refIdx_y >= 0 && refIdx_z >= 0);
        assert(refIdx_x <= numX-(degree+1) && refIdx_y <= numY-(degree+1) && refIdx_z <= numZ-(degree+1));

        // Evaluate
        res(i) = 0.0;
        // calculate the 9 basis function evaluated values
        basisX_v[0] = basisX_(0, refIdx_x)(xT);
        basisX_v[1] = basisX_(1, refIdx_x)(xT);
        basisX_v[2] = basisX_(2, refIdx_x)(xT);
        basisY_v[0] = basisY_(0, refIdx_y)(yT);
        basisY_v[1] = basisY_(1, refIdx_y)(yT);
        basisY_v[2] = basisY_(2, refIdx_y)(yT);
        basisZ_v[0] = basisZ_(0, refIdx_z)(z);
        basisZ_v[1] = basisZ_(1, refIdx_z)(z);
        basisZ_v[2] = basisZ_(2, refIdx_z)(z);
        if (degree == 3) {
            basisX_v[3] = basisX_(3, refIdx_x)(xT);
            basisY_v[3] = basisY_(3, refIdx_y)(yT);
            basisZ_v[3] = basisZ_(3, refIdx_z)(z);
        }

        // loop to calculate summation
        for (iz=0; iz<=degree; iz++)
            for (ix=0; ix<=degree; ix++)
                for (iy=0; iy<=degree; iy++) {

                    res(i) += controlPoints(numX*numY*(refIdx_z+iz) + numY*(refIdx_x+ix) + (refIdx_y+iy)) *
                              basisX_v[ix] * basisY_v[iy] * basisZ_v[iz];
                }
    }
}
// explicit template instantiation
template void bspline::InterpDiskHelper(const DScalar x, const DScalar y, const double z, const DScalar r, Eigen::Matrix<DScalar, Eigen::Dynamic, 1> &res, 
            const basis_t  &basisX_, const basis_t  &basisY_, const basisd_t &basisZ_) const;
template void bspline::InterpDiskHelper(const double  x, const double  y, const double z, const double  r, Eigen::Matrix<double,  Eigen::Dynamic, 1> &res, 
            const basisd_t &basisX_, const basisd_t &basisY_, const basisd_t &basisZ_) const;


void bspline::InterpDisk(const DScalar x, const DScalar y, const double z, const DScalar r, Eigen::Matrix<DScalar, Eigen::Dynamic, 1> &res) const {
// NOTE: This interpolation function cannot query points lying exactly on the surface
//       of the 3D sample grid.

    InterpDiskHelper<basis_t, DScalar>(x, y, z, r, res, basisX, basisY, basisZd);
}


void bspline::InterpDisk(const double  x, const double  y, const double z, const double  r, Eigen::Matrix<double, Eigen::Dynamic, 1> &res) const {
// NOTE: This interpolation function cannot query points lying exactly on the surface
//       of the 3D sample grid.

    InterpDiskHelper<basisd_t, double>(x, y, z, r, res, basisXd, basisYd, basisZd);

        // use control point cache
        /*
        idx = numX*numY*refIdx_z + numY*refIdx_x + refIdx_y;
        res(i) =
            controlPointsCache(0 , idx) * basisX_t[0] * basisY_t[0] * basisZ_t[0] +
            controlPointsCache(1 , idx) * basisX_t[0] * basisY_t[1] * basisZ_t[0] +
            controlPointsCache(2 , idx) * basisX_t[0] * basisY_t[2] * basisZ_t[0] +
            controlPointsCache(3 , idx) * basisX_t[1] * basisY_t[0] * basisZ_t[0] +
            controlPointsCache(4 , idx) * basisX_t[1] * basisY_t[1] * basisZ_t[0] +
            controlPointsCache(5 , idx) * basisX_t[1] * basisY_t[2] * basisZ_t[0] +
            controlPointsCache(6 , idx) * basisX_t[2] * basisY_t[0] * basisZ_t[0] +
            controlPointsCache(7 , idx) * basisX_t[2] * basisY_t[1] * basisZ_t[0] +
            controlPointsCache(8 , idx) * basisX_t[2] * basisY_t[2] * basisZ_t[0] +
            controlPointsCache(9 , idx) * basisX_t[0] * basisY_t[0] * basisZ_t[1] +
            controlPointsCache(10, idx) * basisX_t[0] * basisY_t[1] * basisZ_t[1] +
            controlPointsCache(11, idx) * basisX_t[0] * basisY_t[2] * basisZ_t[1] +
            controlPointsCache(12, idx) * basisX_t[1] * basisY_t[0] * basisZ_t[1] +
            controlPointsCache(13, idx) * basisX_t[1] * basisY_t[1] * basisZ_t[1] +
            controlPointsCache(14, idx) * basisX_t[1] * basisY_t[2] * basisZ_t[1] +
            controlPointsCache(15, idx) * basisX_t[2] * basisY_t[0] * basisZ_t[1] +
            controlPointsCache(16, idx) * basisX_t[2] * basisY_t[1] * basisZ_t[1] +
            controlPointsCache(17, idx) * basisX_t[2] * basisY_t[2] * basisZ_t[1] +
            controlPointsCache(18, idx) * basisX_t[0] * basisY_t[0] * basisZ_t[2] +
            controlPointsCache(19, idx) * basisX_t[0] * basisY_t[1] * basisZ_t[2] +
            controlPointsCache(20, idx) * basisX_t[0] * basisY_t[2] * basisZ_t[2] +
            controlPointsCache(21, idx) * basisX_t[1] * basisY_t[0] * basisZ_t[2] +
            controlPointsCache(22, idx) * basisX_t[1] * basisY_t[1] * basisZ_t[2] +
            controlPointsCache(23, idx) * basisX_t[1] * basisY_t[2] * basisZ_t[2] +
            controlPointsCache(24, idx) * basisX_t[2] * basisY_t[0] * basisZ_t[2] +
            controlPointsCache(25, idx) * basisX_t[2] * basisY_t[1] * basisZ_t[2] +
            controlPointsCache(26, idx) * basisX_t[2] * basisY_t[2] * basisZ_t[2];
        */
}


////////////////////////////////////////////////////////////////////////////////////////
// maintenance methods

struct quadrature quad_;

bspline::bspline() : quad(quad_) {

    DiffScalarBase::setVariableCount(3);  // x, y, r

    degree = 0;
    controlPoints.setZero();
    controlPointsCache.setZero();
    numX = 0;
    numY = 0;
    numZ = 0;
    // By default use Eigen::LeastSquaresConjugateGradient
    leastSquareMethod = 2;
    // solver by default high accuracy
    solverMaxIt = 10000;    // 1000
    solverConvTol = 1e-8;  // 1e-8
    solverTol = 1e-8;      // 1e-8
    // whether enable control points cache
    controlPointsCacheFlag = false;
}


bspline::~bspline() {

}

}  // namespace zebrafish
