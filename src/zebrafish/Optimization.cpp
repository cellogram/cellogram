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

namespace zebrafish {

namespace {

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// class optimization

void optim::Optim_FixDepth(const OptimPara_t &optimPara, const bspline &bsp, const Eigen::MatrixXd &CI, Eigen::MatrixXd &CO, bool invertColor) {

    Eigen::VectorXd EO;
    Eigen::VectorXi IterO;

    optim::Optim_FixDepth(optimPara, bsp, CI, CO, EO, IterO, invertColor);
}


bool optim::Optim_FixDepth(const OptimPara_t &optimPara, const bspline &bsp, const Eigen::VectorXd &CI, Eigen::VectorXd &CO, double &EO, int &IterO, bool invertColor) {

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


void optim::Optim_FixDepth(const OptimPara_t &optimPara, const bspline &bsp, const Eigen::MatrixXd &CI, Eigen::MatrixXd &CO, Eigen::VectorXd &EO, Eigen::VectorXi &IterO, bool invertColor) {

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

}  // namespace zebrafish
