////////////////////////////////////////////////////////////////////////////////
#pragma once
#include <cellogram/common.h>
#include <cellogram/Mesh.h>
#include <cellogram/Mesh3d.h>
#include <cellogram/Region.h>
#include <cellogram/remesh_adaptive.h>
#include <Eigen/Dense>
#include <vector>
#include <string>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

	template<class Mat>
	void write_json_mat(const Mat &mat, nlohmann::json &data)
	{
		data["rows"] = mat.rows();
		data["cols"] = mat.cols();
		data["data"] = std::vector<typename Mat::Scalar>(mat.data(), mat.data() + mat.rows() * mat.cols());
	}

	template<class Mat>
	void read_json_mat(const nlohmann::json &data, Mat &mat)
	{
		const int rows = data["rows"];
		const int cols = data["cols"];
		mat.resize(rows, cols);
		std::vector<typename Mat::Scalar> tmp = data["data"];
		mat = Eigen::Map<Mat>(&tmp[0], rows, cols);
	}

struct State {
public:
	static State &state();
	~State() = default;

private:
	State();

public:
	int lloyd_iterations;
	float energy_variation_from_mean;
	int perm_possibilities;
	float sigma;

	bool image_from_pillars = false;

	MmgOptions mmg_options;

	float target_mesh_size[2] = {0.001f, 0.1f};
	float power = 2.0;
	float padding_size = 25;
	float thickness = 30;

#ifdef NDEBUG
	float target_volume = 0.25;
#else
	float target_volume = 0.5;
#endif

	float scaling; // [um/px]
	float E = 13.58;
	float nu = 0.49;

	float eps = 0.32967032967032966;
	float I = 0.5;
	float L = 3;

	std::string formulation = "LinearElasticity"; //"LinearElasticity"; // NeoHookean

	static const int max_region_vertices = 45;
	static const int gurobi_time_limit_short = 60;
	static const int gurobi_time_limit_long = 300;
	bool fix_regular_regions = false;

	int counter_invalid_neigh;
	int counter_small_region;
	int counter_infeasible_region;
	int counter_solved_regions;

	Mesh mesh;
	Mesh3d mesh3d;

	Eigen::MatrixXd hull_vertices; //needed for lloyd
	Eigen::MatrixXi hull_faces;
	Eigen::MatrixXd hull_polygon;

	Eigen::MatrixXd img;

	std::vector<Region> regions;

	bool load(const std::string &path);
	void load_settings(const std::string &path);
	void load_settings(json args);
	void load_detection_settings(json args);
	void load_analysis_settings(json args);
	bool is_data_available(const std::string &path);
	bool load_image(const std::string fname);
	//bool load_param(const std::string &path);
	bool save(const std::string &path);
	void compute_hull();
	void clean_hull();

	Eigen::VectorXi increase_boundary(const Eigen::VectorXi &boundary);

	void untangle();
	void detect_vertices();
	void relax_with_lloyd();
	void detect_bad_regions();
	void erase_small_regions();
	void check_regions();
	void check_region(const int index);
	void fix_regions();
	void grow_region(const int index);
	void grow_regions();
	void resolve_region(const int index);
	void resolve_regions();
	void final_relax();

	int find_region_by_interior_vertex(const int index);
	void find_region_by_boundary_vertex(const int index, std::vector<int> & regions);

	void split_region(const Eigen::Vector2i& split_end_points);

	void delete_vertex(const int index);
	void add_vertex(Eigen::Vector3d new_point);

	// FEM part
	void mesh_2d_adaptive();
	void extrude_2d_mesh();
	void mesh_3d_uniform();
	void remesh_3d_adaptive();
	void analyze_3d_mesh();

	void reset_state();

public:
	// void foo();
};

} // namespace cellogram
