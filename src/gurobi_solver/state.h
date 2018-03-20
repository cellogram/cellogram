#ifndef GM_STATE
#define GM_STATE

#include <Eigen/Dense>
#include <set>

using namespace Eigen;

// Important convenction to store a hex grid in a matrix (see scheme.png)

class State
{
public:

    // Input
    MatrixXd V_boundary;
    MatrixXd V_internal;
    VectorXi neigh;
    int starting_direction = 0; // does not matter for the algorithm, debug only
    bool oddeven = false; // does not matter for the algorithm, debug only
    bool snap_border = true; // snap the border first
    bool check_contained = true; // never add triangles with other dots inside
    bool avoid_connected_components = true; // never allow to disconnected the unfrozed part

    bool view_vertex_ids = false;

	// added by tobi to be used later and not change data for any subsequent steps by Daniele
	MatrixXd Vdeformed;
	MatrixXd Vperfect;


    // IO
    void load(std::string path);
	void init(const MatrixXd &vB, const MatrixXd &vI, const VectorXi &n);
    void save(std::string path);


    // Phase 1 Data
    MatrixXd V; // The first vertices are the boundary (sorted as in V_boundary), then the internals (unsorted)
    MatrixXi F;

    // Phase 2 Data
    VectorXi fixed;      // Fixed vertices in V
    VectorXi fixed_dotid;      // Fixed vertices pointer to the corresponding id in fixed_dots (-1 is boundary, -2 free)
    VectorXi fixed_dots; // Fixed(used) vertices in V, -1 is free, elsewhere is the id in V


    // Phase 3 Data


    // Phase 4 Data
    std::set<int> to_snap;


    // Phase 1
    void reverse_input();

    void create_initial_triangulation();

    MatrixXi find_pairs_by_tracing(const VectorXi &neigh);

    MatrixXi rasterize_pairs(const MatrixXi &pairs);

    MatrixXi fill(const MatrixXi &raster);

    void extract_mesh_from_raster(const MatrixXi &raster, MatrixXd &V, MatrixXi &F);

    // Phase 2
    void fit_triangulation();

    MatrixXd compute_new_positions(const MatrixXd &V, const MatrixXi &F, const MatrixXi &fixed);

    // Phase 3
	void snap_closest();
    bool is_mesh_valid(const MatrixXd& V, const MatrixXi& F, const VectorXi& fixed, const MatrixXd& P, const VectorXi& fixed_P);

    void unfreeze_overlapping();


    // All together
    void fill_hole();

};

#endif
