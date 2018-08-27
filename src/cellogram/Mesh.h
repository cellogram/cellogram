////////////////////////////////////////////////////////////////////////////////
#pragma once
#include <cellogram/common.h>
#include <cellogram/point_source_detection.h>
#include <Eigen/Dense>
#include <vector>
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
	Eigen::VectorXd sizing; // size map for adaptive remeshing

	Eigen::MatrixXi triangles; // triangular mesh
	std::vector<std::vector<int>> adj; // adjacency list of triangular mesh
	std::vector<std::vector<int>> vertex_to_tri;

	Eigen::Matrix<bool, Eigen::Dynamic, 1> solved_vertex;
	Eigen::VectorXi vertex_status_fixed;
	Eigen::VectorXi added_by_untangler;
	Eigen::MatrixXd deleted_by_untangler;

	// bool load(const std::string &path);
	bool load(const nlohmann::json &data);

	void relax_with_lloyd(const int lloyd_iterations, const Eigen::MatrixXd &hull_vertices, const Eigen::MatrixXi &hull_faces, const bool fix_regular_regions);
	void vertex_degree(Eigen::VectorXi &degree);

	void detect_vertices(const Eigen::MatrixXd &V, const DetectionParams &params);
	void delete_vertex(const int index, bool recompute_triangulation = true);
	void add_vertex(Eigen::Vector3d &new_point, bool reset = true);

	void local_update(Eigen::VectorXi &local2global, Eigen::MatrixXd &new_points, Eigen::MatrixXi &new_triangles, Eigen::VectorXi & old_triangles);
	void local_update(Eigen::VectorXi & local2global, const int global_to_remove, Eigen::MatrixXi & new_triangles);
	void update_triangles_from_split(const Eigen::VectorXi & t_old, const Eigen::MatrixXi & t1, const Eigen::MatrixXi & t2);


	void mark_vertex_as_solved(const Eigen::VectorXi & region_interior);

	void get_physical_bounding_box(double scaling, Eigen::Vector2d &min, Eigen::Vector2d &max) const;

	// Return a triangle mesh embedded in the physical bounding box of the input mesh,
	// where each point is associated a scalar field = norm of the points displacements
	void get_background_mesh(double scaling, Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXd &S, double padding = 0) const;

	void final_relax(const Eigen::VectorXi & expanded_boundary);
	void generate_vertex_to_tri();
	void reset();

	bool untangle();

	void clear();

	// void save(const std::string &path);
	void save(nlohmann::json &data);

private:
	void compute_triangulation();

};

} // namespace cellogram
