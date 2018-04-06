////////////////////////////////////////////////////////////////////////////////
#pragma once

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
	int perm_possibilities = 8;


	int counter_invalid_neigh;
	int counter_small_region;
	int counter_infeasible_region;
	int counter_solved_regions;


	Eigen::MatrixXd detected; // detected (unmoved) point positions
	Eigen::MatrixXd points; // relaxed point positions

	Eigen::MatrixXi triangles; // triangular mesh
	std::vector<std::vector<int>> adj; // adjaceny list of triangluar mesh

	Eigen::VectorXi boundary; // list of vertices on the boundary

	Eigen::MatrixXd hull_vertices; //needed for lloyd
	Eigen::MatrixXi hull_faces;

	Eigen::MatrixXd hull_polygon;

	std::vector<Region> regions;


	bool load(const std::string &path);
	bool save(const std::string &path);
	void compute_hull();
	void compute_triangulation();

	void relax_with_lloyd();
	void detect_bad_regions();
	void fix_regions();
	void grow_region(const int index);
	void resolve_regions();

public:
	// void foo();
};

} // namespace cellogram
