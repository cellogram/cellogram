////////////////////////////////////////////////////////////////////////////////
#pragma once

#include "Mesh.h"
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
	int lloyd_iterations = 12;
	double energy_variation_from_mean = 1.8;
	int perm_possibilities = 12;


	int counter_invalid_neigh;
	int counter_small_region;
	int counter_infeasible_region;
	int counter_solved_regions;

	Mesh mesh;

	Eigen::MatrixXd hull_vertices; //needed for lloyd
	Eigen::MatrixXi hull_faces;
	Eigen::MatrixXd hull_polygon;

	std::vector<Region> regions;

	std::vector<int> fixed_as_good;

	bool load(const std::string &path);
	bool load_param(const std::string &path);
	bool save(const std::string &path);
	void compute_hull();
	void clean_hull();

	void relax_with_lloyd();
	void detect_bad_regions();
	void erase_small_regions();
	void fix_regions();
	void grow_region(const int index);
	void grow_regions();
	void resolve_region(const int index);
	void resolve_regions();
	void final_relax();

	int find_region_by_vertex(const int index);

	void delete_vertex(const int index);
	void add_vertex(Eigen::Vector3d new_point);

	void reset_state();

public:
	// void foo();
};

} // namespace cellogram
