#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <zebrafish/Common.h>
#include <zebrafish/Quadrature.h>
#include <zebrafish/autodiff.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

// A matrix of pre-calculated basis functions (in the form of lambda functions)
typedef Eigen::Matrix<std::function<double (double)>,  Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> basisd_t;
typedef Eigen::Matrix<std::function<DScalar(DScalar)>, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> basis_t;

///////////////////////////////////////
// bspline

class bspline {

private:
    int degree;  // B-spline degree
                 // "2" for quadratic clamped B-spline
                 // "3" for cubic clamped B-spline
    Eigen::Matrix<double, 27, Eigen::Dynamic, Eigen::ColMajor> controlPointsCache;  // optional optimization. By default disabled.
    Eigen::VectorXd controlPoints;  // [#points] control points value

    int Nx, Ny, Nz;         // the dimension of sample points (Nx * Ny * Nz == #pixels)
                            // this is also the dimension of the image
    int numX, numY, numZ;   // the dimension of control points (numX * numY * numZ == #control points)
    double gapX, gapY, gapZ;  // the interval between two control points along a direction
    double resolutionX, resolutionY, resolutionZ;  // The distance between two pixels (in micrometers)

    int leastSquareMethod;            // Which method to solve the least square problem
    int solverMaxIt;                  // Hypre solver max iterations
    double solverConvTol, solverTol;  // Hypre solver "convergence tolerance" and "tolerance"
    bool controlPointsCacheFlag;  // whether enable control points cache. By default inactive.

    // Pre-calculated lambda basis functions (double & DScalar)
    basisd_t basisXd, basisYd, basisZd;  // double  version
    basis_t  basisX,  basisY;            // DScalar version

    void CalcLeastSquareMat(Eigen::SparseMatrix<double, Eigen::RowMajor> &A, Eigen::SparseMatrix<double, Eigen::ColMajor> &Atranspose);
    void SolveLeastSquare(const Eigen::SparseMatrix<double, Eigen::RowMajor> &A, const Eigen::SparseMatrix<double, Eigen::ColMajor> &Atranspose, const Eigen::VectorXd &inputPts);
    void CalcInitialGuess(const Eigen::VectorXd &inputPts, Eigen::VectorXd &initialGuess);
    template <typename T>
    void CalcBasisFunc(Eigen::Matrix< std::function<T(T)>, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> &basisT, const int numT, const double gapT);
    void CreateControlPtsCache(); //Optional
    template <typename basisT, typename T>
    void InterpDiskHelper(const T x, const T y, const double z, const T r, Eigen::Matrix<T, Eigen::Dynamic, 1> &res, const basisT &basisX_, const basisT &basisY_, const basisd_t &basisZ_) const;

public:
    struct quadrature &quad;  // quadrature struct

    void SetResolution(const double resX, const double resY, const double resZ);
    /// Set the microscope resolution in X, Y and Z direction
    /// Unit in micrometers (um)

    void CalcControlPts(   const image_t &image, const double xratio, const double yratio, const double zratio, const int degree);
    void CalcControlPts_um(const image_t &image, const double distX,  const double distY,  const double distZ,  const int degree);
    /// Use least square to calculate an array of control points and store the result
    /// in private variables. This function must be called before any evaluation.
    ///
    /// @param[in]   xratio     { the ratio of (#control points) to (#sample points) along x-axis }
    /// @param[in]   yratio     { the ratio of (#control points) to (#sample points) along y-axis }
    /// @param[in]   zratio     { the ratio of (#control points) to (#sample points) along z-axis }
    /// @param[in]   distX      { the distance between two control points in X-axis. Unit: micrometer }
    /// @param[in]   distY      { the distance between two control points in Y-axis. Unit: micrometer }
    /// @param[in]   distZ      { the distance between two control points in Z-axis. Unit: micrometer }
    /// @param[in]   degree     { the degree of B-spline. Can be 2 or 3. }
    ///

    template <typename T>
    T Interp3D(const T &x, const T &y, const T &z) const;
    /// Calculate the interpolated B-spline result for a given point.
    /// This interface is only used for debug / test. Do not optimize this function.
    /// Warning: It is wrong to use this function before InterpDisk.

    void InterpDisk(const DScalar x, const DScalar y, const double z, const DScalar r, Eigen::Matrix<DScalar, Eigen::Dynamic, 1> &res) const;
    void InterpDisk(const double  x, const double  y, const double z, const double  r, Eigen::Matrix<double , Eigen::Dynamic, 1> &res) const;
    /// Calculate the interpolated B-spline result for one disk layer.
    /// Note: this function does NOT check for input validity.
    ///
    /// @param[in]   x, y, z   { coordinate of the cylinder bottom center }
    /// @param[in]   r         { cylinder radius. }
    /// @param[out]  res       { #DiskSample x 1 output interpolated results }
    ///

    // getter & setter
    bool Get_controlPointsCacheFlag() const { return controlPointsCacheFlag; }
    void Set_controlPointsCacheFlag(bool x) { controlPointsCacheFlag = x; }
    int  Get_degree() const { return degree; }
    void Set_degree(int x) { degree = x; }
    int Get_Nx() const { return Nx;}
    int Get_Ny() const { return Ny;}
    int Get_Nz() const { return Nz;}
    void Set_leastSquareMethod(int x) { leastSquareMethod = x; }
    void Set_solverTol(double x) { solverTol = x; }

    // maintenance methods
    bspline();
    ~bspline();
};

}  // namespace zebrafish
