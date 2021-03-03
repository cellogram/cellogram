#include <zebrafish/Common.h>
#include <zebrafish/Cylinder.h>
#include <zebrafish/autodiff.h>

#ifdef WIN32
# define M_PI           3.14159265358979323846  /* pi */
#endif

#include <cmath>

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
// class cylinder

template<typename T>
bool cylinder::IsValid(const bspline &bsp, const T &x_, const T &y_, const double z, const T &r_, const double h) {

    assert(h > 0);  // we do not optimize h. It should be manually set to be a positive number

    const int xmax = bsp.Get_Nx();
    const int ymax = bsp.Get_Ny();
    const int zmax = bsp.Get_Nz();
    assert(xmax > 0 && ymax > 0 && zmax > 0);  // must update boundary before eval

    static const GetVal<T> getVal;
    const double x = getVal(x_);
    const double y = getVal(y_);
    const double r = getVal(r_);

    // The extended cylinder (union of inner & outer) has radius K*r
    // Also avoid interpolating points lying on the surface of the sample grid
    if (x - cylinder::K*r - 0.5 < 0 || x + cylinder::K*r + 0.5 > xmax - 1) return false;  // x-axis
    if (y - cylinder::K*r - 0.5 < 0 || y + cylinder::K*r + 0.5 > ymax - 1) return false;  // y-axis
    if (z < 0 || z+h >= zmax - 1) return false;  // depth-axis
    if (r < cylinder::Rmin) return false;  // radius

    return true;
}
// explicit template instantiation
template bool cylinder::IsValid(const bspline &bsp, const DScalar &x_, const DScalar &y_, const double z, const DScalar &r_, const double h);
template bool cylinder::IsValid(const bspline &bsp, const double  &x_, const double  &y_, const double z, const double  &r_, const double h);


template <typename T>
void cylinder::EnergyHelper(const bspline &bsp, const Eigen::Matrix<double, Eigen::Dynamic, 1> &zArray,
                            const T &r, const T &x, const T &y, T &resT) {
// This function calculates an intermediate value used by energy evaluation function
// Note: "weightScalar" not multiplied here

    Eigen::Matrix<T, Eigen::Dynamic, 1> interpRes;  // store the results of interpolation
    int i, depth;
    T layerSum;

    //////////////////////////////////////////////////////////////////////////////////
    /// Note about the rationale behind "InterpDisk":
    /// Originally we did the following and call "Interp3D" to evaluate the points:
    /// > Eigen::Matrix<T, Eigen::Dynamic, 2> points;
    /// > points.resize(numPts, 2);
    /// > for (i=0; i<numPts; i++) {
    /// >     // * r + [x, y]
    /// >     points(i, 0) = xyArray(i, 0) * r + x;
    /// >     points(i, 1) = xyArray(i, 1) * r + y;
    /// > }
    /// But this requires some extra memory (points matrix). For efficiency consideration
    /// we used "InterpDisk".
    //////////////////////////////////////////////////////////////////////////////////
    resT = T(0.0);
    for (depth=0; depth<bsp.quad.heightLayers; depth++) {

        // Interpolate all the points in "depth"-th layer
        bsp.InterpDisk(x, y, zArray(depth), r, interpRes);

        assert(bsp.quad.numDiskPts == interpRes.size());
        layerSum = T(0.0);
        for (i=0; i<bsp.quad.numDiskPts; i++) {
            layerSum = layerSum + interpRes(i) * bsp.quad.weightArray(i);  // disk weight
        }
        resT = resT + layerSum * bsp.quad.heightWeightArray(depth);  // height weight
    }
}
// explicit template instantiation
template void cylinder::EnergyHelper(const bspline &bsp, const Eigen::Matrix<double, Eigen::Dynamic, 1> &zArray, const DScalar &r, const DScalar &x, const DScalar &y, DScalar &resT);
template void cylinder::EnergyHelper(const bspline &bsp, const Eigen::Matrix<double, Eigen::Dynamic, 1> &zArray, const double  &r, const double  &x, const double  &y, double  &resT);


template <typename T>
void cylinder::EvaluateCylinder(const bspline &bsp, T x, T y, double z, T r, double h, T &res, bool invert/*=false*/) {

    // kill the staic, zArray should be an input
    Eigen::Matrix<double, Eigen::Dynamic, 1> zArray;  // store the array of depths
    T resInner, resExt;
    assert(bsp.quad.numDiskPts > 0);  // Quadrature method has been chosen

    // depth array
    /// Note: This is simply a very small array with a few doubles
    zArray.resize(bsp.quad.heightLayers, 1);
    zArray = bsp.quad.heightLocArray.array() * (h/2.0) + ((2.0*z+h) / 2.0);

    // weight correction term
    //////////////////////////////////////////////////////////////////////////////////
    /// Notes about the quadrature weight correction term "scalar":
    /// scalar = [disk quadrature jacobian] / [cylinder volumn] * [depth gaussian correction term]
    /// "quad.heightWeightArray" has been multiplied in energy helper
    /// For the inner cylinder with equidistant depth layer - 
    ///                r^2               H
    ///     scalar = --------------- * -----
    ///                pi * r^2 * H      2
    /// For the extended cylinder with equidistant depth layer - 
    ///                    (sqrt(2)*r)^2          H
    ///     scalar = ------------------------ * -----
    ///               pi * (sqrt(2)*r)^2 * H      2
    /// Thus they share the same correction term.
    //////////////////////////////////////////////////////////////////////////////////
    /// Difference between "quad.heightWeightArray" and "depth gaussian correction term":
    /// The height 1D Legendre-Gauss locations and weights are designed for interval [-1, 1]
    /// To integrate an arbitary interval [a, b], we need to calculate
    ///    b - a                   b - a           a + b
    ///   ------- * SUM{ w_i * F( ------- * x_i + ------- ) }
    ///      2                       2               2
    /// "quad.heightWeightArray" is "w_i"
    /// "depth gaussian correction term" is the outer "(b-a)/2"
    //////////////////////////////////////////////////////////////////////////////////
    // Note: Theoretically "weightScalar" should be DScalar, but the radius got
    // cancelled in the formula so it is OK to use double
    double weightScalar = 1.0 / (M_PI * 2.0);

    // Inner area
    EnergyHelper<T>(bsp, zArray, r, x, y, resInner);

    // Extended area
    EnergyHelper<T>(bsp, zArray, r * cylinder::K, x, y, resExt);

    ////////////////////////////////////////////////////
    /// Notes about enerygy function evaluation:
    /// Energy = (Inner density) - (Periperal density)
    ///        = S_in / V_in - S_peri / V_peri
    /// [note: we enforced V_in == V_peri]
    ///        = (1 / V_in) * (S_in + S_in - S_in - S_peri)
    ///        = (1 / V_in) * (2*S_in - S_ext)
    ///        = 2 * (Inner density - Extended density)
    ////////////////////////////////////////////////////
    // double quadrature solution of the energy
    // res = 2.0*(resInner - resExt) * weightScalar;  // old
    res = (resInner - 2*cylinder::alpha * resExt) * weightScalar;  // assume K=sqrt(2)
    res = 2.0 * (res + cylinder::alpha - 0.5);  // scale from [-alpha, 1-alpha] to [-1, 1]
    if (invert) res *= -1;  // invert color

    // DEBUG ONLY
    // std::cerr << "resInner = " << resInner * weightScalar << std::endl;
    // std::cerr << "resExt = " << resExt * weightScalar << std::endl;
}
// explicit template instantiation
template void cylinder::EvaluateCylinder(const bspline &bsp, DScalar x, DScalar y, double z, DScalar r, double h, DScalar &res, bool invert);
template void cylinder::EvaluateCylinder(const bspline &bsp, double  x, double  y, double z, double  r, double h, double  &res, bool invert);


void cylinder::SubtractionHelper(const Eigen::MatrixXd &points, const Eigen::VectorXd &weight, Eigen::VectorXd &resWeight) {

    // this function can be used to implement sigmoid subtraction function
    // [deprecated]
}

}  // namespace zebrafish
