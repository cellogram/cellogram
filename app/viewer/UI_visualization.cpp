////////////////////////////////////////////////////////////////////////////////
#include "UIState.h"
#include <zebrafish/Optimization.h>
#include <zebrafish/Logger.hpp>


////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

using zebrafish::logger;

namespace {

template <typename T>
std::string ToStringWithPrecision(const T a_value, const int n = 0) {
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////
// Special module about visualizing some effects

void UIState::DrawAxisDots() {

    if (!show_axisPoints) return;

    const auto Linspace = [](double s, double e, double gap) -> Eigen::VectorXd {
        int n = std::ceil((e - s) / gap) + 1;
        Eigen::VectorXd res(n);
        for (int i=0; i<n; i++)
            res(i) = i * gap;
        res(n-1) = e;
        return res;
    };

    const auto ToMat = [](const Eigen::VectorXd &v1, const Eigen::VectorXd &v2, const Eigen::VectorXd &v3) -> Eigen::MatrixXd {
        Eigen::MatrixXd res(v1.rows(), 3);
        res << v1, v2, v3;
        return res;
    };

    const int imgRows = state.img3D[0].rows();
    const int imgCols = state.img3D[0].cols();
    const int layerPerImg = state.img3D.size();
    static double gapr = (imgRows > 200) ? double(imgRows)/40.0 : 5;
    static double gapc = (imgCols > 200) ? double(imgCols)/40.0 : 5;
    static double gapz = (layerPerImg > 60) ? double(layerPerImg)/20.0 : 3;
    static int nr = std::ceil(double(imgRows) / gapr)+1;
    static int nc = std::ceil(double(imgCols) / gapc)+1;
    static int nz = std::ceil(double(layerPerImg) / gapz)+1;
    static Eigen::MatrixXd loc_r;
    static Eigen::MatrixXd label_loc_r;
    static Eigen::MatrixXd loc_c;
    static Eigen::MatrixXd label_loc_c;
    static Eigen::MatrixXd loc_z;
    static Eigen::MatrixXd label_loc_z;
    static float visual_z_mult_cache = -1;
    if (visual_z_mult_cache != imgViewer.visual_z_mult) {
        visual_z_mult_cache = imgViewer.visual_z_mult;
        loc_r = ToMat(Eigen::VectorXd::Zero(nr), Linspace(0, imgRows, gapr), Eigen::VectorXd::Zero(nr));
        label_loc_r = ToMat(Eigen::VectorXd::Ones(nr).array() * (-10.0), Linspace(0, imgRows, gapr).array() + 0.5, Eigen::VectorXd::Zero(nr));
        loc_c = ToMat(Linspace(0, imgCols, gapc), Eigen::VectorXd::Ones(nc).array() * imgRows, Eigen::VectorXd::Zero(nc));
        label_loc_c = ToMat(Linspace(0, imgCols, gapc).array() - 0.5, Eigen::VectorXd::Ones(nc).array() * (imgRows+6.0), Eigen::VectorXd::Zero(nc));
        loc_z = ToMat(Eigen::VectorXd::Ones(nz).array() * imgCols, Eigen::VectorXd::Ones(nz).array() * imgRows, Linspace(0, layerPerImg*imgViewer.visual_z_mult, gapz*imgViewer.visual_z_mult));
        label_loc_z = ToMat(Eigen::VectorXd::Ones(nz).array() * (imgCols+2.0), Eigen::VectorXd::Ones(nz).array() * imgRows, Linspace(0, layerPerImg*imgViewer.visual_z_mult, gapz*imgViewer.visual_z_mult).array() + 1.2);
    }

    // add points
    visual_data().add_points(loc_r, colorUI.ref_point_color);
    visual_data().add_points(loc_c, colorUI.ref_point_color);
    visual_data().add_points(loc_z, colorUI.ref_point_color);

    // add labels
    for (int i=0; i<nr; i++)
        visual_data().add_label(label_loc_r.row(i), ToStringWithPrecision(imgRows - loc_r(i, 1)));
    for (int i=0; i<nc; i++)
        visual_data().add_label(label_loc_c.row(i), ToStringWithPrecision(loc_c(i, 0)));
    for (int i=0; i<nz; i++)
        visual_data().add_label(label_loc_z.row(i), ToStringWithPrecision(loc_z(i, 2) / imgViewer.visual_z_mult));
}


void UIState::DrawAllMarkerIdx() {

    if (!show_allIndex) return;
    const int N = state.mesh.marker_4D.rows();
    if (N == 0) return;

    for (int i=0; i<N; i++) {
        Eigen::VectorXd p(3);
        p << state.mesh.marker_4D(i, 0), state.mesh.marker_4D(i, 1), state.mesh.marker_4D(i, 2) * imgViewer.visual_z_mult + 0.05;
        visual_data().add_label(p, std::to_string(i));
    }
}


void UIState::DrawWarnViewer() {

    const auto &dsFlag = state.mesh.dsFlag;
    const int N = dsFlag.size();
    Eigen::MatrixXd pts = state.mesh.marker_4D.leftCols(3);
    pts.col(2) *= imgViewer.visual_z_mult;
    pts.col(2).array() += 0.2;

    warn_data().point_size = 10;
    warn_data().show_faces = false;
	warn_data().show_lines = false;

    for (int i=0; i<N; i++) {
        if (dsFlag[i] == zebrafish::DepthSearchFlag_t::SecondDerivative) {
            warn_data().add_points(pts.row(i), colorUI.warn_snd_derivative);
            visual_data().add_label(pts.row(i), std::to_string(i));
        }
        if (dsFlag[i] == zebrafish::DepthSearchFlag_t::InvalidEnergy) {
            warn_data().add_points(pts.row(i), colorUI.warn_invalid_energy);
            visual_data().add_label(pts.row(i), std::to_string(i));
        }
    }
}

}  // namespace cellogram
