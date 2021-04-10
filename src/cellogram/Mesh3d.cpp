////////////////////////////////////////////////////////////////////////////////
#include "Mesh3d.h"
#include <cellogram/Mesh.h>
#include <cellogram/State.h>
#include <cellogram/remesh_adaptive.h>
#include <geogram/mesh/mesh_io.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/writeOBJ.h>
#include <igl/write_triangle_mesh.h>
#include <igl/slice.h>
#include <igl/remove_unreferenced.h>
#include <polyfem/Mesh3D.hpp>
#include <polyfem/MeshUtils.hpp>
#include <polyfem/PointBasedProblem.hpp>
#include <polyfem/State.hpp>
#include <polyfem/VTUWriter.hpp>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

namespace {

void ToPhysicalUnit(Eigen::MatrixXd &mat, double scaling, double zscaling) {
    mat.leftCols(2) *= scaling;
    mat.col(2) *= zscaling;
}

nlohmann::json compute_analysis(const Eigen::MatrixXd &vertices_, const Eigen::MatrixXi &faces, const Eigen::MatrixXi &tets, const Mesh &mesh,
                                float thickness, float E, float nu, const std::string &formulation, double scaling, double zscaling,
                                Eigen::MatrixXd &displacement, Eigen::MatrixXd &traction_forces, const std::string &save_dir) {
    //TODO
    //const std::string rbf_function = "gaussian";
    const std::string rbf_function = "thin-plate";
    // const std::string rbf_function = "cubic";
    const double eps = 4;

    static const bool export_data = false;
    assert(tets.cols() == 4);
    assert(vertices_.cols() == 3);

    Eigen::MatrixXd vertices = vertices_;
    // ToPhysicalUnit(vertices, scaling, zscaling);
    // xy already in physical unit

    // create a GEOGRAM tet mesh
    GEO::Mesh M;
    M.vertices.create_vertices((int)vertices.rows());
    for (int i = 0; i < (int)M.vertices.nb(); ++i) {
        GEO::vec3 &p = M.vertices.point(i);
        p[0] = vertices(i, 0);
        p[1] = vertices(i, 1);
        p[2] = vertices(i, 2);
    }
    M.cells.create_tets((int)tets.rows());
    for (int c = 0; c < (int)M.cells.nb(); ++c) {
        for (int lv = 0; lv < tets.cols(); ++lv) {
            M.cells.set_vertex(c, lv, tets(c, lv));
        }
    }
    M.cells.connect();

    // igl::write_triangle_mesh("/Users/ziyizhang/Projects/tmp/tetmesh.obj", vertices, faces);
    // std::cerr << "/Users/ziyizhang/Projects/tmp/tetmesh.obj" << std::endl;

    if (export_data) {
        GEO::mesh_save(M, "mesh.mesh");
        GEO::mesh_load("mesh.mesh", M);
    }

    json j_args = {
        {"problem", "PointBasedTensor"},
        {"normalize_mesh", false},

        {"tensor_formulation", formulation},

        {"discr_order", 1},
        {"vismesh_rel_area", 1000},

        {"solver_params", {"conv_tol", 1e-7, "max_iter", 2000}},

        {"nl_solver_rhs_steps", 4},
        {"n_boundary_samples", 3},

        {"params", {
                       {"E", E},
                       {"nu", nu},
                   }},
    };

    polyfem::State state;
    std::string log_file = "";
    bool is_quiet = false;
    int log_level = 1;

    state.init_logger(log_file, log_level, is_quiet);
    state.init(j_args);

    state.load_mesh(M, [&](const polyfem::RowVectorNd &bary) {
        // top, Id = 1
        if (std::abs(bary(2)) < 1e-6) {
            return 1;
        }

        // Bottom, Id = 3
        if (std::abs(bary(2) - -thickness) < 1e-8) {
            return 3;
        }

        //any other
        return 2;
    });

    // state.compute_mesh_stats();

    polyfem::PointBasedTensorProblem &problem = *dynamic_cast<polyfem::PointBasedTensorProblem *>(state.problem.get());

    // Zebrafish depth information
    Eigen::MatrixXd marker_3D;
    Mesh3d::GetMarker3D(mesh.marker_4D, marker_3D);

    // pixel -> um
    Eigen::MatrixXd disp = (marker_3D - mesh.points);
    Eigen::MatrixXd pts = mesh.points;
    ToPhysicalUnit(disp, scaling, zscaling);
    ToPhysicalUnit(pts, scaling, zscaling);

    // igl::write_triangle_mesh("/Users/ziyizhang/Projects/tmp/pts.obj", pts, mesh.triangles);
    // std::cerr << "/Users/ziyizhang/Projects/tmp/pts.obj" << std::endl;
    // std::cerr << "disp" << std::endl << disp << std::endl;
    // std::cerr << "pts" << std::endl << pts << std::endl;
    // igl::write_triangle_mesh("/Users/ziyizhang/Projects/tmp/marker_3d.obj", marker_3D, mesh.triangles);
    // igl::write_triangle_mesh("/Users/ziyizhang/Projects/tmp/mesh_points.obj", mesh.points, mesh.triangles);

    //Id = 1, func, mesh, coord =2, means skip z for the interpolation
    Eigen::Matrix<bool, 3, 1> dirichet_dims;
    dirichet_dims(0) = dirichet_dims(1) = true;
    dirichet_dims(2) = true; // z is not dirichet
    problem.add_function(1, disp, pts, rbf_function, eps, -1, dirichet_dims);

    //Id = 3, zero Dirichelt
    problem.add_constant(3, Eigen::Vector3d(0, 0, 0));

    // state.compute_mesh_stats();

    state.build_basis();
    // state.build_polygonal_basis();

    state.assemble_rhs();
    state.assemble_stiffness_mat();

    state.solve_problem();

    //true = compute average insteat of just integral

    // state.interpolate_boundary_function(vertices, faces, state.sol, true, vals);
    // state.interpolate_boundary_tensor_function(vertices, faces, state.sol, true, traction_forces);

    // zebrafish: this part has been moved to "OutputHelper"
    // state.interpolate_boundary_function_at_vertices(vertices, faces, state.sol, vals);
    // Eigen::MatrixXd stresses, mises;
    // state.interpolate_boundary_tensor_function(vertices, faces, state.sol, vals, true, traction_forces, stresses, mises);

    // vals = Eigen::Map<Eigen::MatrixXd>(state.sol.data(), 3, vertices.rows());
    // vals = Eigen::Map<Eigen::MatrixXd>(state.rhs.data(), 3, vertices.rows());
    // vals = vals.transpose().eval();
    // std::cout<<vals<<std::endl;


    // zebrafish export
    // output
    const auto OutputHelper = [&state, &save_dir, &displacement, &traction_forces](const Eigen::MatrixXd &mesh_v, const MatrixXi &mesh_f) {
        Eigen::MatrixXd stress;
        Eigen::MatrixXd mises;

        state.interpolate_boundary_function_at_vertices(mesh_v, mesh_f, state.sol, displacement);
        state.interpolate_boundary_tensor_function(mesh_v, mesh_f, state.sol, displacement, true, traction_forces, stress, mises);

        // per-triangle data -> per-vertex data by taking average
        const auto ToVertexData = [&mesh_v, &mesh_f](
                                    const Eigen::MatrixXd &traction_forces_tri,
                                    const Eigen::MatrixXd &stress_tri,
                                    const Eigen::MatrixXd &mises_tri,
                                    Eigen::MatrixXd &traction_forces_ver,
                                    Eigen::MatrixXd &stress_ver,
                                    Eigen::MatrixXd &mises_ver) {
            traction_forces_ver = Eigen::MatrixXd::Zero(mesh_v.rows(), 3);
            stress_ver = Eigen::MatrixXd::Zero(mesh_v.rows(), 9); // 3x3 matrix
            mises_ver = Eigen::MatrixXd::Zero(mesh_v.rows(), 1);

            Eigen::VectorXd area, vertex_area(mesh_v.rows());
            vertex_area.setZero();
            igl::doublearea(mesh_v, mesh_f, area);

            for (int f = 0; f < mesh_f.rows(); ++f) {
                for (int d = 0; d < 3; ++d) {
                    int vid = mesh_f(f, d);
                    traction_forces_ver.row(vid) += traction_forces_tri.row(f) * area(f);
                    stress_ver.row(vid) += stress_tri.row(f) * area(f);
                    mises_ver(vid) += mises_tri(f) * area(f);
                    vertex_area(vid) += area(f);
                }
            }

            for (int d = 0; d < 3; ++d) {
                traction_forces_ver.col(d).array() /= vertex_area.array();
            }
            for (int d = 0; d < 9; ++d) {
                stress_ver.col(d).array() /= vertex_area.array();
            }
            mises_ver.array() /= vertex_area.array();
        };

        // triangle to vertex
        Eigen::MatrixXd traction_force_v;
        Eigen::MatrixXd stress_v;
        Eigen::MatrixXd mises_v;

        ToVertexData(traction_forces, stress, mises, traction_force_v, stress_v, mises_v);

        // extract the upper surface only
        const double thres = -1e-9;
        std::vector<int> vid_upper_vec;
        for (int i=0; i<mesh_v.rows(); i++)
            if (mesh_v(i, 2) > thres) vid_upper_vec.push_back(i);
        std::vector<int> fid_upper_vec;
        for (int i=0; i<mesh_f.rows(); i++)
            if (mesh_v(mesh_f(i, 0), 2) > thres && mesh_v(mesh_f(i, 1), 2) > thres && mesh_v(mesh_f(i, 2), 2) > thres) fid_upper_vec.push_back(i);
        Eigen::VectorXi vid_upper = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(vid_upper_vec.data(), vid_upper_vec.size());
        Eigen::VectorXi fid_upper = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(fid_upper_vec.data(), fid_upper_vec.size());

        // write to VTU
        polyfem::VTUWriter VTUwriter;
        Eigen::MatrixXd displacement_upper = igl::slice(displacement, vid_upper, 1);
        VTUwriter.add_field("displacement", displacement_upper);
        Eigen::MatrixXd traction_force_v_upper = igl::slice(traction_force_v, vid_upper, 1);
        VTUwriter.add_field("traction_forces", traction_force_v_upper);
        Eigen::MatrixXd stress_v_upper = igl::slice(stress_v, vid_upper, 1);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                const std::string index = "_" + std::to_string(i) + std::to_string(j);
                VTUwriter.add_field("stress" + index, stress_v_upper.col(i * 3 + j));
            }
        }
        Eigen::MatrixXd mises_v_upper = igl::slice(mises_v, vid_upper, 1);
        VTUwriter.add_field("von_mises", mises_v_upper);

        Eigen::MatrixXd NV;
        Eigen::MatrixXi NF, I;
        igl::remove_unreferenced(mesh_v, igl::slice(mesh_f, fid_upper, 1), NV, NF, I);
        if (VTUwriter.write_mesh(save_dir + ".vtu", NV, NF)) {
            std::cout << "Result saved to " + save_dir + ".vtu" << std::endl;
        } else {
            std::cerr << "Export failed: path not found" << std::endl;
            std::cerr << save_dir + ".vtu" << std::endl;
        }
    };  // OutputHelper lambda function

    OutputHelper(vertices, faces);
    state.save_vtu(save_dir + ".all.vtu", 0);
    // state.save_surface(save_dir + ".surf.vtu");

    if (export_data) {
        std::ofstream out("sol.txt");
        out.precision(100);
        out << state.sol << std::endl;
        out.close();
    }

    if (export_data) {
        Eigen::MatrixXd nodes(state.n_bases, state.mesh->dimension());
        for (const auto &eb : state.bases) {
            for (const auto &b : eb.bases) {
                for (const auto &lg : b.global()) {
                    nodes.row(lg.index) = lg.node;
                }
            }
        }
        std::ofstream out("nodes.txt");
        out.precision(100);
        out << nodes;
        out.close();
    }

    nlohmann::json json;
    state.save_json(json);
    return json;

    // auto &tmp_mesh = *dynamic_cast<polyfem::Mesh3D *>(state.mesh.get());
    // igl::opengl::glfw::Viewer viewer;
    // Eigen::MatrixXi asdT;
    // Eigen::MatrixXd asdP, asdC;
    // std::vector<int> asdR;
    // tmp_mesh.triangulate_faces(asdT,asdP,asdR);

    // asdC.resize(asdT.rows(), 3);
    // asdC.setZero();

    // for(std::size_t i = 0; i < tmp_mesh.n_faces(); ++i)
    // {
    // 	const int a =tmp_mesh.get_boundary_id(i);
    // 	const auto bb = tmp_mesh.face_barycenter(i);

    // 	if(a == 1)
    // 		viewer.data().add_points(bb, Eigen::RowVector3d(1, 0, 0));
    // 	else if(a == 3)
    // 		viewer.data().add_points(bb, Eigen::RowVector3d(0, 1, 0));
    // 	else if(a == 2)
    // 		viewer.data().add_points(bb, Eigen::RowVector3d(0, 0, 1));
    // }

    // viewer.data().set_mesh(asdP, asdT);
    // viewer.launch();
}
} // namespace

void Mesh3d::init_pillars(const Mesh &mesh, float eps, float I, float L, double scaling, double zscaling) {
    clear();

    displacement = (mesh.detected - mesh.points) * scaling;
    V = mesh.points * scaling;
    ToPhysicalUnit(displacement, scaling, zscaling);
    ToPhysicalUnit(V, scaling, zscaling);

    const float bending_force = 3 * eps * I / (L * L * L);
    traction_forces = bending_force * displacement;
}

bool Mesh3d::empty() {
    return V.size() == 0;
}

bool Mesh3d::analysed() {
    return traction_forces.size() > 0;
}

void Mesh3d::init_nano_dots(const Mesh &mesh, float padding_size, const float thickness, float E, float nu, double scaling, double zscaling, const std::string &formulation, const std::string &save_dir) {
    //Uncomment to used not adaptive tetgen mesher
    // 		clear();

    // 		Eigen::Vector2d max_dim;
    // 		Eigen::Vector2d min_dim;
    // 		mesh.get_physical_bounding_box(min_dim, max_dim);

    // 		double xMin = min_dim(0) - padding_size;
    // 		double xMax = max_dim(0) + padding_size;
    // 		double yMin = min_dim(1) - padding_size;
    // 		double yMax = max_dim(1) + padding_size;
    // 		double zMin = -thickness;
    // 		double zMax = 0;

    // 		V.resize(8, 3);
    // 		V << xMin, yMin, zMin,
    // 		xMin, yMax, zMin,
    // 		xMax, yMax, zMin,
    // 		xMax, yMin, zMin,

    // 			//4
    // 		xMin, yMin, zMax,
    // 		xMin, yMax, zMax,
    // 		xMax, yMax, zMax,
    // 		xMax, yMin, zMax;

    // 		F.resize(12, 3);
    // 		F << 1, 2, 0,
    // 		0, 2, 3,

    // 		5, 4, 6,
    // 		4, 7, 6,

    // 		1, 0, 4,
    // 		1, 4, 5,

    // 		2, 1, 5,
    // 		2, 5, 6,

    // 		3, 2, 6,
    // 		3, 6, 7,

    // 		0, 3, 7,
    // 		0, 7, 4;

    // 		Eigen::MatrixXd TV;
    // 		Eigen::MatrixXi TT;
    // 		Eigen::MatrixXi TF;
    // #ifdef NDEBUG
    // 		igl::copyleft::tetgen::tetrahedralize(V, F, "Qpq1.414a100", TV, TT, TF);
    // #else
    // 		igl::copyleft::tetgen::tetrahedralize(V, F, "Qpq1.414a1000", TV, TT, TF);
    // #endif
    // 		V = TV;
    // 		F = TF;
    // 		T = TT;

    simulation_out = compute_analysis(V, F, T, mesh, thickness, E, nu, formulation, scaling, zscaling, displacement, traction_forces, save_dir);
}

void Mesh3d::clear() {
    F.resize(0, 0);
    V.resize(0, 0);

    T.resize(0, 0);
    displacement.resize(0, 0);
    traction_forces.resize(0, 0);

    simulation_out = nlohmann::json({});
}

bool Mesh3d::load(const nlohmann::json &data) {
    read_json_mat(data["V"], V);
    read_json_mat(data["F"], F);
    //read_json_mat(data["T"], T); // maybe the Tets are unnecessary

    read_json_mat(data["displacement"], displacement);
    read_json_mat(data["traction_forces"], traction_forces);

    return true;
}

void Mesh3d::save_mesh(nlohmann::json &data) {
    data["V"] = json::object();
    write_json_mat(V, data["V"]);

    data["F"] = json::object();
    write_json_mat(F, data["F"]);

    //data["T"] = json::object();
    //write_json_mat(T, data["T"]);

    data["displacement"] = json::object();
    write_json_mat(displacement, data["displacement"]);
}
void Mesh3d::save_traction(nlohmann::json &data) {
    data["traction_forces"] = json::object();
    write_json_mat(traction_forces, data["traction_forces"]);
}

void Mesh3d::GetMarker3D(const Eigen::MatrixXd &marker4d, Eigen::MatrixXd &marker3d) {

    marker3d = marker4d.leftCols(3);
    // try to translate so most points have z=0 (they are background and shouldn't move too much)
    marker3d.col(2).array() -= marker3d.col(2).mean();
    double sum_z = 0.0;
    int cnt_z = 0;
    for (int i=0; i<marker3d.rows(); i++) {
        if (std::fabs(marker3d(i, 2)) < 0.5) {
            cnt_z++;
            sum_z += marker3d(i, 2);
        }
    }
    marker3d.col(2).array() -= sum_z / double(cnt_z);
}

} // namespace cellogram
