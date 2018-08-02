#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <cellogram/common.h>
#include <cellogram/Mesh.h>
#include <Eigen/Dense>
#ifdef CELLOGRAM_WITH_GUROBI
#include <gurobi_solver/state.h>
#endif
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {
	class Region
	{
	public:
		static const int NOT_CHECKED = -6;
		static const int REGION_TOO_LARGE = -5;
		static const int TOO_FEW_VERTICES = -4;
		static const int TOO_MANY_VERTICES = -3;
		static const int NOT_PROPERLY_CLOSED = -2;
		static const int NO_SOLUTION = -1;
		static const int OK = 0;

		static std::string pretty_status(const int status);

		bool verbose = true;


		int status = NOT_CHECKED;
		int points_delta = 0; // difference between points in region an expected (<0: missing; >0: too many)

		Eigen::VectorXi region_boundary;
		Eigen::VectorXi region_interior;
		Eigen::VectorXi region_triangles;

		Region() { }

		int size();

		void compute_edges(const Eigen::MatrixXd &V, Eigen::MatrixXd &bad_P1, Eigen::MatrixXd &bad_P2);

		void grow(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const std::vector<std::vector<int>> &vertex_to_tri);

		void find_points(const Eigen::VectorXi &region_ids, const int id);
		void find_triangles(const std::vector<std::vector<int>> &vertex_to_tri, const Eigen::VectorXi &region_ids, const int id);
		void find_triangles(const Eigen::MatrixXi & F, bool force_add_boundaries = false);
		void bounding(const Eigen::MatrixXi &F, const Eigen::MatrixXd &V);

		void check_region(const Eigen::MatrixXd &V_detected, const Eigen::MatrixXd &V_current, const Eigen::MatrixXi &F, const std::vector<std::vector<int>> &adj);
		void check_region_from_local_tris(const Eigen::MatrixXd &V_detected, const Eigen::MatrixXd &V_current, const Eigen::MatrixXi &F2, const std::vector<std::vector<int>> &adj);
		void resolve(const Eigen::MatrixXd &V_detected, const Eigen::MatrixXd &V_current, const int perm_possibilities, const double gurobi_time_limit, Eigen::MatrixXd  &new_points, Eigen::MatrixXi &new_triangles, bool force_solve = false);

		bool fix_missing_points(const Eigen::MatrixXi & F);

		void delete_vertex(cellogram::Mesh &mesh, const int index);


		void local_to_global(Eigen::VectorXi &local2global);
		bool add_orphaned_triangles(const Mesh &mesh);
		bool fix_pinched_region(const Mesh &mesh);
		//todo
		//void reduce_index

		bool split_region(Mesh &mesh, const Eigen::Vector2i & split_end_points, Region &r1, Region  &r2);

	private:
		void find_split_path(const Mesh &mesh, const Eigen::Vector2i & split_end_points, Eigen::VectorXi & path);
		void find_split_boundaries(const Eigen::Vector2i & split_end_points, const Eigen::VectorXi & path, Eigen::VectorXi & boundary1, Eigen::VectorXi & boundary2);
		void find_interior_V(const Mesh & mesh, const Eigen::VectorXi & boundary, Eigen::VectorXi & interior);
		void triangluate_region(const Mesh mesh, Eigen::MatrixXi &new_triangles);
		bool clean_up_boundary(const Eigen::MatrixXi &tri, const Eigen::MatrixXi & tri_list);
		Eigen::MatrixXi get_triangulation(const Eigen::MatrixXi &F);

#ifdef CELLOGRAM_WITH_GUROBI
		gurobi::State s;
#endif
	};
} // namespace cellogram
