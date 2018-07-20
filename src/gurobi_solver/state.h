#ifndef GM_STATE
#define GM_STATE

#include <Eigen/Dense>
#include <set>

// Important convenction to store a hex grid in a matrix (see scheme.png)
namespace cellogram {
	namespace gurobi {

		class State
		{
		public:

			// Input
			Eigen::MatrixXd V_boundary;
			Eigen::MatrixXd V_internal;
			Eigen::VectorXi neigh;
			int starting_direction = 0; // does not matter for the algorithm, debug only
			bool oddeven = false; // does not matter for the algorithm, debug only
			bool snap_border = true; // snap the border first
			bool check_contained = true; // never add triangles with other dots inside
			bool avoid_connected_components = true; // never allow to disconnected the unfrozed part

			bool view_vertex_ids = false;

			// added by tobi to be used later and not change data for any subsequent steps by Daniele
			Eigen::MatrixXd Vdeformed;
			Eigen::MatrixXd Vperfect;


			// IO
			void load(std::string path);
			void init(const Eigen::MatrixXd &vB, const Eigen::MatrixXd &vI, const Eigen::VectorXi &n);
			void save(std::string path);


			// Phase 1 Data
			Eigen::MatrixXd V; // The first vertices are the boundary (sorted as in V_boundary), then the internals (unsorted)
			Eigen::MatrixXi F;

			// Phase 2 Data
			Eigen::VectorXi fixed;      // Fixed vertices in V
			Eigen::VectorXi fixed_dotid;      // Fixed vertices pointer to the corresponding id in fixed_dots (-1 is boundary, -2 free)
			Eigen::VectorXi fixed_dots; // Fixed(used) vertices in V, -1 is free, elsewhere is the id in V


			// Phase 3 Data


			// Phase 4 Data
			std::set<int> to_snap;


			// Phase 1
			void reverse_input();

			void create_initial_triangulation();

			Eigen::MatrixXi find_pairs_by_tracing(const Eigen::VectorXi &neigh);

			Eigen::MatrixXi rasterize_pairs(const Eigen::MatrixXi &pairs);

			Eigen::MatrixXi fill(const Eigen::MatrixXi &raster);

			void extract_mesh_from_raster(const Eigen::MatrixXi &raster, Eigen::MatrixXd &V, Eigen::MatrixXi &F);

			// Phase 2
			void fit_triangulation();

			Eigen::MatrixXd compute_new_positions(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXi &fixed);

			// Phase 3
			void snap_closest();
			bool is_mesh_valid(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::VectorXi& fixed, const Eigen::MatrixXd& P, const Eigen::VectorXi& fixed_P);

			void unfreeze_overlapping();


			// All together
			void fill_hole();

		};

	} // namespace gurobi
} // namespace cellogram

#endif
