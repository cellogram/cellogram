////////////////////////////////////////////////////////////////////////////////
#pragma once

#include "Mesh.h"
#include "Mesh3d.h"
#include "common.h"
#include "Region.h"
#include <vector>
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

struct State {
public:
	static State &state();
	~State() = default;

private:
	State() { }

public:
	int lloyd_iterations = 20;
	float energy_variation_from_mean = 1.7;
	int perm_possibilities = 12;
	float sigma = 2.2;

	float padding_size = 25;
	float thickness = 30;

	static const int max_region_vertices = 50;
	static const int gurobi_time_limit_short = 60;
	static const int gurobi_time_limit_long = 300;

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
	bool is_data_available(const std::string &path);
	bool load_image(const std::string fname);
	//bool load_param(const std::string &path);
	bool save(const std::string &path);
	void compute_hull();
	void clean_hull();

	Eigen::VectorXi increase_boundary(const Eigen::VectorXi &boundary);

#ifdef WITH_UNTANGLER
	void untangle();
#endif
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

	void init_3d_mesh();
	
	void reset_state();

public:
	// void foo();
};

} // namespace cellogram
