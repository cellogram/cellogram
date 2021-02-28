#include <zebrafish/Common.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/Optimization.h>
#include <zebrafish/Cylinder.h>
#include <zebrafish/autodiff.h>

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>

#include <LBFGS.h>
#include <cmath>
#include <limits>

namespace zebrafish {

namespace {

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// class optimization

void optim::Optim_FixDepth(const OptimPara_t &optimPara, const bspline &bsp, const Eigen::MatrixXd &CI, Eigen::MatrixXd &CO, const bool invertColor/*=false*/) {

    Eigen::VectorXd EO;
    Eigen::VectorXi IterO;

    optim::Optim_FixDepth(optimPara, bsp, CI, CO, EO, IterO, invertColor);
}


bool optim::Optim_FixDepth(const OptimPara_t &optimPara, const bspline &bsp, const Eigen::VectorXd &CI, Eigen::VectorXd &CO, double &EO, int &IterO, const bool invertColor/*=false*/) {

    Eigen::MatrixXd CI_mat(1, 4);
    Eigen::MatrixXd CO_mat;
    Eigen::VectorXd EO_vec;
    Eigen::VectorXi IterO_vec;

    CI_mat.row(0) = CI;
    optim::Optim_FixDepth(optimPara, bsp, CI_mat, CO_mat, EO_vec, IterO_vec, invertColor);
    CO = CO_mat.row(0);
    EO = EO_vec(0);
    IterO = IterO_vec(0);

    return (IterO < optimPara.maxIt && EO < 1.0);
}


void optim::Optim_FixDepth(
    const OptimPara_t &optimPara, 
    const bspline &bsp, 
    const Eigen::MatrixXd &CI, 
          Eigen::MatrixXd &CO, 
          Eigen::VectorXd &EO, 
          Eigen::VectorXi &IterO, 
    const bool invertColor/*=false*/) {

    // prepare LBFGS
    LBFGSpp::LBFGSParam<double> param;
    param.epsilon = optimPara.epsilon;
    param.max_iterations = optimPara.maxIt;

    // prepare sampleInput & sampleOutput for Newtons
    const int N = CI.rows();
    CO.resize(N, 4);  // output cylinder x, y, z, r
    EO.resize(N);  // output energy array
    IterO.resize(N);  // output iteration array

    // Optimization
    tbb::parallel_for( tbb::blocked_range<int>(0, N),
        //////////////////////////////////////
        // lambda function for parallel_for //
        //////////////////////////////////////
        [&bsp, &CI, &CO, &EO, &IterO, &param, invertColor]
        (const tbb::blocked_range<int> &r) {

        // NOTE: LBFGSSolver is NOT thread safe. This must be instantiated for every thread
        LBFGSpp::LBFGSSolver<double> solver(param);

        // NOTE: the "variable count" used by "Autodiff" will be stored in 
        //       thread-local memory, so this must be set for every thread
        DiffScalarBase::setVariableCount(3);

        for (int ii = r.begin(); ii != r.end(); ++ii) {    

                Eigen::VectorXd vec(3, 1);
                vec(0) = CI(ii, 0);  // x
                vec(1) = CI(ii, 1);  // y
                vec(2) = CI(ii, 3);  // r
                double res;

                ///////////////////////////////////
                // lambda function for optimizer //
                ///////////////////////////////////
                auto func = [&bsp, &CI, invertColor, ii]
                (const Eigen::VectorXd& x, Eigen::VectorXd& grad) {

                        DScalar ans;

                        if (!cylinder::IsValid(bsp, x(0), x(1), CI(ii, 2), x(2), cylinder::H)) {
                            grad.setZero();
                            return 1.0;
                        }
                        cylinder::EvaluateCylinder(bsp, DScalar(0, x(0)), DScalar(1, x(1)), CI(ii, 2), DScalar(2, x(2)), cylinder::H, ans, invertColor);
                        grad = ans.getGradient();
                        return ans.getValue();
                    };
                // NOTE: the template of "solver.minimize" does not accept a temprary variable (due to non-const argument)
                //       so we define a "func" and pass it in
                int it = solver.minimize(func, vec, res);
                ///////////////////////////////////

                CO(ii, 0) = vec(0);     // x
                CO(ii, 1) = vec(1);     // y
                CO(ii, 2) = CI(ii, 2);  // z
                CO(ii, 3) = vec(2);     // r
                EO(ii) = res;        // energy
                IterO(ii) = it;      // iteration
        }
    });
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void optim::Optim_WithDepth(
    const OptimPara_t &optimPara, 
    const bspline &bsp, 
    const int zNum, 
    const double zGap, 
    const Eigen::VectorXd &CI, 
          OptimDepthInfo_t &C_depth_info,
    const bool invertColor/*=false*/) {

    Eigen::MatrixXd CI_mat(1, 4);
    CI_mat.row(0) = CI;
    std::vector<OptimDepthInfo_t> C_depth_info_vec;
    optim::Optim_WithDepth(optimPara, bsp, zNum, zGap, CI_mat, C_depth_info_vec, invertColor);
    C_depth_info = C_depth_info_vec.at(0);
}


void optim::Optim_WithDepth(
    const OptimPara_t &optimPara, 
    const bspline &bsp, 
    const int zNum, 
    const double zGap, 
    const Eigen::MatrixXd &CI, 
          std::vector<OptimDepthInfo_t> &C_depth_info,
    const bool invertColor/*=false*/) {

    const int N = CI.rows();
    const int M = zNum * 2 + 1;
    C_depth_info.resize(N);
    if (N == 0) return;

    // prepare for parallel optimization
    Eigen::VectorXd zDeltaArray = Eigen::VectorXd::LinSpaced(2*zNum+1, -zNum*zGap, zNum*zGap);
    Eigen::MatrixXd CI_withZ(N*M, 4);
    for (int i=0; i<N; i++) {
        CI_withZ.block(i*M, 0, M, 4) = CI.row(i).replicate(M, 1);
        CI_withZ.block(i*M, 2, M, 1) += zDeltaArray;
    }
    Eigen::MatrixXd CO_withZ;
    Eigen::VectorXd EO;
    Eigen::VectorXi IterO;

    // optimization
    optim::Optim_FixDepth(optimPara, bsp, CI_withZ, CO_withZ, EO, IterO, invertColor);

    // output
    for (int i=0; i<N; i++) {
        C_depth_info[i].C = CO_withZ.middleRows(i*M, M);
        C_depth_info[i].energy = EO.segment(i*M, M);
        C_depth_info[i].iter = IterO.segment(i*M, M);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void optim::DepthSelection(
    const OptimPara_t &optimPara, 
    const Eigen::VectorXd &CI, 
    const OptimDepthInfo_t &C_depth_info, 
          Eigen::VectorXd &CO, 
          DepthSearchFlag_t &flag) {
    
    Eigen::MatrixXd CI_mat(1, 4);
    CI_mat.row(0) = CI;
    std::vector<OptimDepthInfo_t> C_depth_info_vec;
    C_depth_info_vec.push_back(C_depth_info);
    Eigen::MatrixXd CO_mat;
    std::vector<DepthSearchFlag_t> flag_vec;
    optim::DepthSelection(optimPara, CI_mat, C_depth_info_vec, CO_mat, flag_vec);
    CO = CO_mat.row(0);
    flag = flag_vec.at(0);
}


void optim::DepthSelection(
    const OptimPara_t &optimPara, 
    const Eigen::MatrixXd &CI, 
    const std::vector<OptimDepthInfo_t> &C_depth_info, 
          Eigen::MatrixXd &CO, 
          std::vector<DepthSearchFlag_t> &flag) {

    const int N = CI.rows();
    CO.resize(N, 4);
    flag.resize(N);
    if (N == 0) return;
    const int M = C_depth_info[0].C.rows();  // depth trials

    Eigen::VectorXd energy_smooth;
    Eigen::VectorXd second_derivative;
    for (int i=0; i<N; i++) {

        const Eigen::VectorXd &E_raw = C_depth_info[i].energy;
        // smooth curve
        // for i-th index, find 5 adjacent numbers, calc mean after discarding the min and max in those 5 numbers
        energy_smooth = Eigen::VectorXd::NullaryExpr(M, [&E_raw](Eigen::Index i) {
            const int l = 2;
            int left = std::max(int(i-l), 0);
            int right = std::min(int(i+l), int(E_raw.size()-1));

            Eigen::VectorXd E_raw_segment = E_raw.segment(left, right-left+1);
            double sum = E_raw_segment.sum();
            double minE = E_raw_segment.minCoeff();
            double maxE = E_raw_segment.maxCoeff();
            return (sum - minE - maxE) / double(right-left-1);
        });
        // 1D discrete laplacian operator
        second_derivative = Eigen::VectorXd::NullaryExpr(M, [&energy_smooth](Eigen::Index i) {
            if (i == 0 || i == int(energy_smooth.size()-1))
                return 0.0;
            return energy_smooth(i-1) - 2 * energy_smooth(i) + energy_smooth(i+1);
        });

        //////////////////////////////////////////////////////////
        double minE = std::numeric_limits<double>::max();  // or simply 1.0
        int minIdx = -1;
        for (int j=0; j<M; j++) {
            // if this valid? (xy does not move to another location)
            Eigen::Vector2d xyDisp;
            xyDisp << C_depth_info[i].C(j, 0) - CI(i, 0), C_depth_info[i].C(j, 1) - CI(i, 1);
            if (xyDisp.norm() > optimPara.zSearchMaxXYDisp) {
                energy_smooth(j) = 1.2;  // a little bit larger than 1.0 to distinguish (for debug)
                continue;
            }
            // smaller energy?
            if (energy_smooth(j) < minE) {
                minE = energy_smooth(j);
                minIdx = j;
            }
        }

        // minE == 1.0
        if (minE == 1.0 || minIdx == -1) {
            char errorMsg[100];
            std::sprintf(errorMsg, "> [warning] No valid energy. Too close to the boundary or other failures. Marker index %d at [%.2f, %.2f, %.2f].", i, CI(i, 0), CI(i, 1), CI(i, 2));
            logger().warn(errorMsg);
            std::cerr << errorMsg << std::endl;
            flag[i] = InvalidEnergy;
            continue;
        }
        // minIdx at end point
        if (M>0 && (minIdx == 0 || minIdx == M-1)) {
            char warnMsg[100];
            std::sprintf(warnMsg, "> [warning] Abnormal second derivative. Marker index %d at [%.2f, %.2f, %.2f].", i, CI(i, 0), CI(i, 1), CI(i, 2));
            logger().debug(warnMsg);
            std::cerr << warnMsg << std::endl;
            flag[i] = SecondDerivative;

            // std::cerr << E_raw.transpose() << std::endl;  // DEBUG PURPOSE
            continue;
        }
        // in case the smoothing shifts the actual min index
        for (int j=std::max(0, minIdx-1); j<=std::min(M-1, minIdx+1); j++) {
            if (E_raw(j) < minE) {
                minE = E_raw(j);
                minIdx = j;
            }
        }
        // save result to CO
        flag[i] = Success;
        CO.row(i) << C_depth_info[i].C.row(minIdx);
    }
}

}  // namespace zebrafish
