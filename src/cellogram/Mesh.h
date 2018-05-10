////////////////////////////////////////////////////////////////////////////////
#pragma once

#include <cellogram/point_source_detection.h>

#include "common.h"
#include <vector>
#include <Eigen/Dense>

////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

class Mesh {

public:
	Mesh() { }

	Eigen::MatrixXd detected; // detected (unmoved) point positions
	DetectionParams params; // parameters of detected
	Eigen::MatrixXd moved; // manually moved point from detected positions
	Eigen::MatrixXd points; // relaxed point positions
	Eigen::VectorXi boundary; // list of vertices on the boundary

	Eigen::MatrixXi triangles; // triangular mesh
	std::vector<std::vector<int>> adj; // adjaceny list of triangluar mesh
	std::vector<std::vector<int>> vertex_to_tri;

	Eigen::Matrix<bool, Eigen::Dynamic, 1> solved_vertex;
	Eigen::VectorXi vertex_status_fixed;

	float scaling = 1; // [um/px]

	bool load(const std::string &path);

	void relax_with_lloyd(const int lloyd_iterations, const Eigen::MatrixXd &hull_vertices, const Eigen::MatrixXi &hull_faces);
	void vertex_degree(Eigen::VectorXi &degree);

	void detect_vertices(const Eigen::MatrixXd &V, const DetectionParams &params);
	void delete_vertex(const int index, bool recompute_triangulation = true);
	void add_vertex(Eigen::Vector3d &new_point);

	void local_update(Eigen::VectorXi &local2global, Eigen::MatrixXd &new_points, Eigen::MatrixXi &new_triangles, Eigen::VectorXi & old_triangles);
	void local_update(Eigen::VectorXi & local2global, const int global_to_remove, Eigen::MatrixXi & new_triangles);

	void mark_vertex_as_solved(const Eigen::VectorXi & region_interior);

	void get_physical_bounding_box(Eigen::Vector2d &min, Eigen::Vector2d &max) const;

	void final_relax();
	void generate_vertex_to_tri();
	void reset();

#ifdef WITH_UNTANGLER
	void untangle();
#endif

	void clear();

	void save(const std::string &path);

private:
	void compute_triangulation();

};

} // namespace cellogram
