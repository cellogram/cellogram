/*
 *  Copyright (c) 2012-2014, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine,
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX
 *     FRANCE
 *
 */

#include <geogram_gfx/basic/GLSL.h>
#include <geogram_gfx/basic/GL.h>
#include <geogram_gfx/GLUP/GLUP.h>
#include <geogram_gfx/glup_viewer/glup_viewer.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/logger.h>
#include <geogram/numerics/predicates.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <fstream>

#include <gurobi_solver/state.h>
#include <gurobi_solver/generateQ.h>
#include <gurobi_solver/gurobiModel.h>
#include <igl/unique.h>

/*********************************************************************/

#include <algorithm>
#include <numeric>



struct Compare {
	const std::vector<GEO::vec2> &points;
	int leftmost; // Leftmost point of the poly

	Compare(const std::vector<GEO::vec2> &points_) : points(points_) { }

	bool operator ()(int i1, int i2) {
		if (i2 == leftmost) { return false; }
		if (i1 == leftmost) { return true; }
		GEO::vec2 u = points[i1] - points[leftmost];
		GEO::vec2 v = points[i2] - points[leftmost];
		double d = det(u, v);
		return (d < 0 || (d == 0 && u.length2() < v.length2()));
	}
};

bool inline salientAngle(const GEO::vec2 &a, const GEO::vec2 &b, const GEO::vec2 &c) {
	return (det(b - a, c - a) >= 0);
}

std::vector<int> convex_hull(const std::vector<GEO::vec2> &points) {
	Compare order(points);
	order.leftmost = 0;
	for(int i = 1; i < (int) points.size(); ++i) {
		if (points[i][0] < points[order.leftmost][0]) {
			order.leftmost = i;
		} else if (points[i][0] == points[order.leftmost][0]
			&& points[i][1] < points[order.leftmost][1])
		{
			order.leftmost = i;
		}
	}
	std::vector<int> index(points.size());
	std::iota(index.begin(), index.end(), 0);
	std::sort(index.begin(), index.end(), order);
	std::vector<int> hull;
	for (int i : index) {
		hull.push_back(i);
		while (hull.size() > 3u && salientAngle(points[hull.end()[-3]],
			points[hull.end()[-2]], points[hull.end()[-1]]))
		{
			hull.end()[-2] = hull.back();
			hull.pop_back(); // Pop inner vertices
		}
	}
	return hull;
}

void mark_fixed(const std::vector<GEO::vec2> &points, std::vector<bool> &fixed, std::vector<GEO::vec2> &border) {
	auto hull = convex_hull(points);
	border.clear();
	for (int i : hull) {
		fixed[i] = true;
		border.push_back(points[i]);
	}
	std::reverse(border.begin(), border.end());
}

/*********************************************************************/

namespace {

	using namespace GEO;

	MatrixXi tGlobal;
	MatrixXi triangles;
	std::vector<std::string> filenames;
	typedef std::complex<double> Point;
	std::vector<bool> fixed;
	std::vector<bool> roi;
	std::vector<bool> roiInternal;
	std::vector<int> vBoundaryInd;
	std::vector<vec2> points;
	std::vector<vec2> vInternal;
	std::vector<vec2> detected;  // same as reference below, i.e. the distorted mesh.
	std::vector<vec2> reference; // same as input, i.e. the distorted mesh.
	bool showCellogramDelaunay = true;
	Delaunay_var delaunay;
	std::vector<int> degree;
	typedef std::vector<vec2> Polygon;
	Polygon border;
	Polygon vBoundary;
	float interpolation = 0.0;

	void Lloyd_relaxation();
	void solve_ROI();
	void relaxLaplacian();

	void convex_clip_polygon(
		const Polygon& P, const Polygon& clip, Polygon& result
	);

	vec2 centroid(const Polygon& P);

	void get_bbox(const std::vector<vec2> &pos, vec2 &xy_min, vec2 &xy_max) {
		xy_min = pos[0];
		xy_max = pos[0];
		for (index_t i = 0; i < pos.size(); ++i) {
			xy_min[0] = std::min(pos[i][0], xy_min[0]);
			xy_min[1] = std::min(pos[i][1], xy_min[1]);
			xy_max[0] = std::max(pos[i][0], xy_max[0]);
			xy_max[1] = std::max(pos[i][1], xy_max[1]);
		}
	}

	/**
	 * \brief Updates the Delaunay triangulation with the current
	 *  vector of points.
	 */
	void update_Delaunay() {
		delaunay->set_vertices((int) points.size(), &points.data()->x);
		degree.resize(points.size());
		for(index_t i=0; i<degree.size(); ++i) {
			degree[i] = 0;
		}
		for(index_t c=0; c<delaunay->nb_cells(); ++c) {
			const signed_index_t* cell = delaunay->cell_to_v() + 3*c;
			for(index_t e=0; e<3; ++e) {
				signed_index_t v1 = cell[e];
				signed_index_t v2 = cell[(e+1)%3];
				if (v1 < v2) {
					degree[v1]++;
					degree[v2]++;
				}
			}
		}
	}

	/**
	 * \brief Creates random points.
	 * \details Points are created uniformly in the [0,1]x[0,1]
	 *  square
	 * \param[in] nb the number of points to create.
	 */
	void create_random_points(index_t nb) {
		for(index_t i=0; i<nb; ++i) {
			fixed.push_back(false);
			points.push_back(
				vec2(0.25, 0.25) +
				vec2(Numeric::random_float64()*0.5, Numeric::random_float64()*0.5)
			);
		}
		update_Delaunay();
	}

	void load_points(const std::string &filename, std::vector<vec2> &target) {
		GEO::Mesh M;
		if(!GEO::mesh_load(filename, M)) {
			return;
		}
		target.clear();
		for (index_t i = 0; i < M.vertices.nb(); ++i) {
			GEO::vec3 p = M.vertices.point(i);
			target.push_back(vec2(p[0], p[1]));
		}
		fixed.resize(target.size());
		roi.resize(target.size());
		roiInternal.resize(target.size());
		update_Delaunay();
	}

	/**
	 * \brief Gets the index of the point that matches a given
	 *  point.
	 * \details The current GLUP point size is taken into account
	 *  to determine the tolerance.
	 * \param[in] p a const reference coordinate to a world-space
	 *  point.
	 * \param[in] get_nearest if true, gets the nearest point, else
	 *  only gets the point that matches p (up to current value of
	 *  point size when reprojected on the screen)
	 * \return the index of the point that corresponds to p if it
	 *  exists, or index_t(-1) if no such point exists.
	 */
	index_t get_picked_point(const vec2& p, bool get_nearest = false) {
		if(points.size() == 0) {
			return index_t(-1);
		}
		double dist_so_far = Numeric::max_float64();
		index_t nearest = index_t(-1);
		for(index_t i=0; i<points.size(); ++i) {
			double dist = distance2(p, points[i]);
			if(dist < dist_so_far) {
				nearest = i;
				dist_so_far = dist;
			}
		}
		if(!get_nearest) {
			vec3 q(points[nearest].x, points[nearest].y, 0.0);
			vec3 p_scr;
			vec3 q_scr;
			glup_viewer_project(p.x, p.y, 0, &p_scr.x, &p_scr.y, &p_scr.z);
			glup_viewer_project(q.x, q.y, 0, &q_scr.x, &q_scr.y, &q_scr.z);
			if(distance2(p_scr,q_scr) > 10.0 * double(glupGetPointSize())) {
				nearest = index_t(-1);
			}
		}
		return nearest;
	}

	/**
	 * \brief The size of all the displayed points.
	 */
	GLint point_size = 10;

	void set_border_as_polygon(index_t nb_sides) {
		border.clear();
		for(index_t i=0; i<nb_sides; ++i) {
			double alpha = double(i) * 2.0 * M_PI / double(nb_sides);
			double s = sin(alpha);
			double c = cos(alpha);
			border.push_back( vec2(0.5*(c + 1.0), 0.5*(s + 1.0)) );
		}
	}

	void set_border_shape(int shape) {
		return;
		static int current_shape=-1;
		if (current_shape != shape) {
			current_shape = shape;
		}
		switch(shape) {
		case 1:
			set_border_as_polygon(3);
			break;
		case 2:
			set_border_as_polygon(5);
			break;
		case 3:
			set_border_as_polygon(100);
			break;
		case 4:
			{
				vec2 xyz_min, xyz_max;
				get_bbox(points, xyz_min, xyz_max);
				border.clear();
				border.push_back(vec2(xyz_min[0], xyz_min[1]));
				border.push_back(vec2(xyz_max[0], xyz_min[1]));
				border.push_back(vec2(xyz_max[0], xyz_max[1]));
				border.push_back(vec2(xyz_min[0], xyz_max[1]));
			}
			break;
		default:
			border.clear();
			border.push_back(vec2(0.0, 0.0));
			border.push_back(vec2(1.0, 0.0));
			border.push_back(vec2(1.0, 1.0));
			border.push_back(vec2(0.0, 1.0));
			break;
		}
	}

	/**
	 * \brief Initializes OpenGL objects.
	 * \details Specifed as glup_viewer_set_init_func() callback.
	 */
	void init() {
		GEO::Graphics::initialize();

		glup_viewer_set_background_color(1.0, 1.0, 1.0);
		glup_viewer_add_toggle(
			'T', glup_viewer_is_enabled_ptr(GLUP_VIEWER_TWEAKBARS),
			"Toggle tweakbars"
		);

		glup_viewer_disable(GLUP_VIEWER_BACKGROUND);
		glup_viewer_disable(GLUP_VIEWER_3D);

		delaunay = Delaunay::create(2,"BDEL2d");
		create_random_points(3);

		load_points(filenames[0], points);
		load_points(filenames[0], detected);
		load_points(filenames[1], reference);
		if (filenames.size() > 2) {
			std::ifstream in(filenames[2]);
			int	x;
			while (in >> x) {
				if (x >= 0 && x < (int) fixed.size()) {
					fixed[x] = true;
				}
			}
		} else {
			mark_fixed(points, fixed, border);
		}

		vec2 xyz_min, xyz_max;
		get_bbox(points, xyz_min, xyz_max);
		glup_viewer_set_region_of_interest(
			(float) xyz_min[0], (float) xyz_min[1], 0.f,
			(float) xyz_max[0], (float) xyz_max[1], 1.f
		);
	}

	/**
	 * \brief Displays the border of the domain.
	 */
	void display_border() {
		glupSetColor3f(GLUP_FRONT_AND_BACK_COLOR, 0.0, 0.0, 0.0);
		glupSetMeshWidth(4);
		glupBegin(GLUP_LINES);
		for(index_t i=0; i<border.size(); ++i) {
			glupVertex(border[i]);
			glupVertex(border[(i+1)%border.size()]);
		}
		glupEnd();
	}

	/**
	 * \brief Displays the points.
	 */
	void display_points(const std::vector<vec2> _points, bool ref = 0) {
		glupEnable(GLUP_LIGHTING);
		glupSetPointSize(GLfloat(point_size));
		glupEnable(GLUP_VERTEX_COLORS);
		glupBegin(GLUP_POINTS);
		for(index_t i=0; i<_points.size(); ++i) {
			if (ref) {
				glupColor3f(1.0f, 1.0f, 0.0f);
			} else if (i < degree.size() && roiInternal[i]) {
				glupColor3f(0.0f, 0.0f, 1.0f);
			} else if (i < roi.size() && roi[i]) {
				glupColor3f(1.0f, 1.0f, 0.0f);
			} else if (i < fixed.size() && fixed[i]) {
				glupColor3f(0.0f, 1.0f, 0.0f);
			} else if (i < degree.size() && degree[i] != 6) {
				glupColor3f(1.0f, 0.0f, 0.0f);
			} else {
				glupColor3f(0.0f, 1.0f, 1.0f);
			}
			if (ref == false) {
				vec2 p = _points[i] * (1.f - interpolation) + interpolation * reference[i];
				glupVertex(p);
			} else {
				glupVertex(_points[i]);
			}
		}
		glupEnd();
		glupColor3f(0.0f, 0.0f, 0.0f);
		glupDisable(GLUP_VERTEX_COLORS);
		glupDisable(GLUP_LIGHTING);
	}

	/**
	 * \brief Displays the expected matching.
	 */
	void display_matching(const std::vector<vec2> from, const std::vector<vec2> to, bool color) {
		glupSetMeshWidth(3);
		if (color) {
			glupSetColor3f(GLUP_FRONT_AND_BACK_COLOR, 1.0f, 0.1f, 0.0f);
		} else {
			glupSetColor3f(GLUP_FRONT_AND_BACK_COLOR, 0.1f, 0.1f, 0.0f);
		}
		glupBegin(GLUP_LINES);
		for(index_t i=0; i<points.size(); ++i) {
			glupVertex(from[i]);
			glupVertex(to[i]);
		}
		glupEnd();
	}

	/**
	 * \brief Displays the Delaunay triangles.
	 */
	void display_Delaunay_triangles() {
		glupSetColor3f(GLUP_FRONT_AND_BACK_COLOR, 0.7f, 0.7f, 0.7f);
		glupSetMeshWidth(1);
		glupBegin(GLUP_LINES);
		if (showCellogramDelaunay)
		{
			for (index_t c = 0; c<delaunay->nb_cells(); ++c) {
				const signed_index_t* cell = delaunay->cell_to_v() + 3 * c;
				for (index_t e = 0; e<3; ++e) {
					signed_index_t v1 = cell[e];
					signed_index_t v2 = cell[(e + 1) % 3];
					vec2 p1 = points[v1] * (1.f - interpolation) + interpolation * reference[v1];
					vec2 p2 = points[v2] * (1.f - interpolation) + interpolation * reference[v2];
					glupVertex(p1);
					glupVertex(p2);
				}
			}
		}
		else
		{
			for (index_t i = 0; i<triangles.rows(); ++i) {
				for (index_t j = 0; j<3; ++j) {
					signed_index_t v1 = triangles(i, j);
					signed_index_t v2 = triangles(i, (j + 1) % 3);
					vec2 p1 = points[v1] * (1.f - interpolation) + interpolation * reference[v1];
					vec2 p2 = points[v2] * (1.f - interpolation) + interpolation * reference[v2];
					glupVertex(p1);
					glupVertex(p2);
				}
			}
		}
		glupEnd();
	}

	/**
	 * \brief Gets the circumcenter of a triangle.
	 * \param[in] t the index of the triangle, in 0..delaunay->nb_cells()-1
	 * \return the circumcenter of triangle \p t
	 */
	vec2 circumcenter(index_t t) {
		signed_index_t v1 = delaunay->cell_to_v()[3*t];
		signed_index_t v2 = delaunay->cell_to_v()[3*t+1];
		signed_index_t v3 = delaunay->cell_to_v()[3*t+2];
		vec2 p1(delaunay->vertex_ptr(index_t(v1)));
		vec2 p2(delaunay->vertex_ptr(index_t(v2)));
		vec2 p3(delaunay->vertex_ptr(index_t(v3)));
		return Geom::triangle_circumcenter(p1,p2,p3);
	}

	/**
	 * \brief Gets an infinite vertex in the direction normal to an
	 *  edge on the boundary of the convex hull.
	 * \param[in] t the index of the triangle, in 0..delaunay->nb_cells()-1
	 * \param[in] e the local index of the edge, in {0,1,2}
	 * \return a point located far away along the direction normal to the
	 *  edge \p e of triangle \p t
	 */
	vec2 infinite_vertex(index_t t, index_t e) {
		index_t lv1 = (e+1)%3;
		index_t lv2 = (e+2)%3;
		index_t v1 = index_t(delaunay->cell_to_v()[3*t+lv1]);
		index_t v2 = index_t(delaunay->cell_to_v()[3*t+lv2]);
		vec2 p1(delaunay->vertex_ptr(v1));
		vec2 p2(delaunay->vertex_ptr(v2));
		vec2 n = normalize(p2-p1);
		n = vec2(n.y, -n.x);
		return 0.5*(p1+p2)+100000.0*n;
	}

	/**
	 * \brief Displays the Voronoi edges.
	 */
	void display_Voronoi_edges() {
		glupSetColor3f(GLUP_FRONT_AND_BACK_COLOR, 0.3f, 0.3f, 0.3f);
		glupSetMeshWidth(2);
		glupBegin(GLUP_LINES);
		for(index_t t=0; t<delaunay->nb_cells(); ++t) {
			vec2 cc = circumcenter(t);
			for(index_t e=0; e<3; ++e) {
				signed_index_t t2 = delaunay->cell_to_cell()[3*t+e];
				if(t2 == -1) {
					glupVertex(cc);
					glupVertex(infinite_vertex(t,e));
				} else if(t2 >signed_index_t(t)) {
					glupVertex(cc);
					glupVertex(circumcenter(index_t(t2)));
				}
			}
		}
		glupEnd();
	}

	/**
	 * \brief Finds the local index of a vertex in a triangle.
	 * \details Throws an assertion failure if the triangle \p t is
	 *  not incident to vertex \p v
	 * \param[in] t the triangle, in 0..delaunay->nb_cells()-1
	 * \param[in] v the vertex, in 0..delaunay->nb_vertices()-1
	 * \return the local index of v in t, in {0,1,2}
	 */
	index_t find_vertex(index_t t, index_t v) {
		for(index_t lv=0; lv<3; ++lv) {
			if(index_t(delaunay->cell_to_v()[3*t+lv]) == v) {
				return lv;
			}
		}
		geo_assert_not_reached;
	}

	/**
	 * \brief Gets a Voronoi cell of a vertex
	 * \details The vertex is specified by a triangle and a local index in
	 *  the triangle
	 * \param[in] t0 the triangle
	 * \param[in] lv the local index of the vertex in triangle \p t0
	 * \param[out] cell a reference to the Voronoi cell
	 */
	void get_Voronoi_cell(index_t t0, index_t lv, Polygon& cell) {
		cell.resize(0);
		index_t v = index_t(delaunay->cell_to_v()[3*t0+lv]);
		bool on_border = false;
		index_t t = t0;

		// First, we turn around the vertex v. To do that, we compute
		// lv, the local index of v in the current triangle. Following
		// the standard numerotation of a triangle, edge number lv is
		// not incident to vertex v. The two other edges (lv+1)%3 and
		// (lv+2)%3 of the triangle are incident to vertex v. By always
		// traversing (lv+1)%3, we turn around vertex v.
		do {
			index_t e = (lv+1)%3;
			signed_index_t neigh_t = delaunay->cell_to_cell()[3*t+e];
			if(neigh_t == -1) {
				on_border = true;
				break;
			}
			cell.push_back(circumcenter(t));
			t = index_t(neigh_t);
			lv = find_vertex(t,v);
		} while(t != t0);


		// If one traversed edge is on the border of the convex hull, then
		// we empty the cell, and start turning around the vertex in the other
		// direction, i.e. by traversing this time edge (lv+2)%3 until we
		// reach the other edge on the border of the convex hull that is
		// incident to v.
		if(on_border) {
			cell.resize(0);
			cell.push_back(infinite_vertex(t,(lv + 1)%3));
			for(;;) {
				cell.push_back(circumcenter(t));
				index_t e = (lv+2)%3;
				signed_index_t neigh_t = delaunay->cell_to_cell()[3*t+e];
				if(neigh_t == -1) {
					cell.push_back(infinite_vertex(t, e));
					break;
				}
				t = index_t(neigh_t);
				lv = find_vertex(t,v);
			}
		}

		Polygon clipped;
		convex_clip_polygon(cell, border, clipped);
		cell.swap(clipped);
	}


	/**
	 * \brief Displays the Voronoi cells as filled polygons with
	 *  random colors.
	 */
	void display_Voronoi_cells() {
		std::vector<bool> v_visited(delaunay->nb_vertices());
		Polygon cell;
		glupEnable(GLUP_VERTEX_COLORS);
		glupBegin(GLUP_TRIANGLES);
		for(index_t t=0; t<delaunay->nb_cells(); ++t) {
			for(index_t lv=0; lv<3; ++lv) {
				index_t v = index_t(delaunay->cell_to_v()[3*t+lv]);
				if(!v_visited[v]) {
					glup_viewer_random_color_from_index(int(v));
					v_visited[v] = true;
					get_Voronoi_cell(t,lv,cell);
					for(index_t i=1; i+1<cell.size(); ++i) {
						glupVertex(cell[0]);
						glupVertex(cell[i]);
						glupVertex(cell[i+1]);
					}
				}
			}
		}
		glupEnd();
		glupDisable(GLUP_VERTEX_COLORS);
	}


	bool show_Voronoi_cells = false;
	bool show_Delaunay_triangles = true;
	bool show_Voronoi_edges = false;
	bool show_points = true;
	bool show_reference = false;
	bool show_matching = false;
	bool show_border = true;
	bool move_vertices = false;
	bool select_fixed = false;
	bool select_roi = false;
	/**
	 * \brief Draws all the elements of the Delaunay triangulation /
	 *  Voronoi diagram.
	 * \details Specified as glup_viewer_set_display_func() callback.
	 */
	void display() {
		if(show_Voronoi_cells) {
			display_Voronoi_cells();
		}
		if(show_Delaunay_triangles) {
			display_Delaunay_triangles();
		}
		if(show_Voronoi_edges) {
			display_Voronoi_edges();
		}
		if(show_points) {
			display_points(points);
		}
		if(show_reference) {
			display_points(reference, true);
		}
		if(show_matching) {
			display_matching(points, reference, false);
			display_matching(points, detected, true);
		}
		if(show_border) {
			display_border();
		}
		if(glup_viewer_is_enabled(GLUP_VIEWER_IDLE_REDRAW)) {
			Lloyd_relaxation();
		}
	}

	/**
	 * \brief Displays and manages the GUI.
	 */
	void overlay() {
		static int border_shape = 4;

		ImGui::SetNextWindowPos(ImVec2(20, 20), ImGuiSetCond_Once);
		ImGui::SetNextWindowSize(ImVec2(180, 380), ImGuiSetCond_Once);

		ImGui::Begin("Delaunay", nullptr, ImGuiWindowFlags_AlwaysAutoResize);

		ImGui::Text("Display options");
		ImGui::Checkbox("Voronoi cells", &show_Voronoi_cells);
		ImGui::Checkbox("Voronoi edges", &show_Voronoi_edges);
		ImGui::Checkbox("Delaunay triangles", &show_Delaunay_triangles);
		ImGui::Checkbox("points", &show_points);
		ImGui::Checkbox("reference", &show_reference);
		ImGui::Checkbox("matching", &show_matching);
		ImGui::Checkbox("border", &show_border);

		ImGui::Separator();
		ImGui::Text("Lloyd");
		if(ImGui::Button("one iteration")) {
			Lloyd_relaxation();
		}
		ImGui::Checkbox(
			"animate",
			(bool*)glup_viewer_is_enabled_ptr(GLUP_VIEWER_IDLE_REDRAW)
		);
		ImGui::Checkbox("move points", &move_vertices);
		ImGui::Checkbox("select fixed", &select_fixed);

		//
		ImGui::PushItemWidth(-60);
		ImGui::DragFloat("advect", &interpolation, 0.01f, 0.f, 1.f);
		ImGui::PopItemWidth();

		ImGui::Separator();
		ImGui::Text("Gurobi");
		ImGui::PushItemWidth(-60);
		//ImGui::Combo("", &border_shape, "square\0triangle\0pentagon\0circle\0bbox\0\0");
		// TOBI
		ImGui::Checkbox("select ROI", &select_roi);
		if (ImGui::Button("Solve ROI")) {
			solve_ROI();
		}
		if (ImGui::Button("relax Laplacian")) {
			relaxLaplacian();
		}
		ImGui::PopItemWidth();

		ImGui::Separator();
		if(ImGui::Button("reset points")) {
			load_points(filenames[0], points);
		}

		set_border_shape(border_shape);

		ImGui::End();
	}


	/**
	 * \brief The callback function called for each mouse event
	 * \param[in] x , y the screen-space coordinates of the mouse pointer
	 * \param[in] button the button that changed state
	 * \param[in] event one of GLUP_VIEWER_UP, GLUP_VIEWER_DOWN, GLUP_VIEWER_MOVE
	 */
	GLboolean mouse(float x, float y, int button, enum GlupViewerEvent event) {
		const index_t NO_POINT = index_t(-1);
		static index_t picked_point = NO_POINT;
		static index_t last_button = index_t(-1);

		GEO::geo_argused(x);
		GEO::geo_argused(y);

		GLdouble xyz[3];
		GLboolean bkgnd;
		glup_viewer_get_picked_point(xyz, &bkgnd);
		vec2 p(xyz[0], xyz[1]);

		if(glup_viewer_is_enabled(GLUP_VIEWER_IDLE_REDRAW)) {
			switch(event) {
			case GLUP_VIEWER_DOWN:
				last_button = index_t(button);
				break;
			case GLUP_VIEWER_UP:
				last_button = index_t(-1);
				break;
			case GLUP_VIEWER_MOVE:
				break;
			default:
				break;
			}
			(void) last_button; // ignore

			if (!move_vertices) { return GL_FALSE; }
			// switch(last_button) {
			// case 0: {
			//     points.push_back(p);
			//     picked_point = points.size() - 1;
			// } break;
			// case 1: {
			//     picked_point = get_picked_point(p,true);
			//     if(points.size() > 3) {
			//         points.erase(points.begin() + int(picked_point));
			//     }
			// } break;
			// }
			// update_Delaunay();
			if(event == GLUP_VIEWER_DOWN) {
				picked_point = get_picked_point(p);
			}
			if(event == GLUP_VIEWER_MOVE && move_vertices) {
				if (picked_point != NO_POINT) {
					points[picked_point] = p;
					update_Delaunay();
				}
				return GL_TRUE;
			}
			if(event == GLUP_VIEWER_MOVE && select_fixed) {
				picked_point = get_picked_point(p);
				if (picked_point != NO_POINT) {
					fixed[picked_point] = true;
				}
				return GL_TRUE;
			}
			// TOBI
			if (event == GLUP_VIEWER_DOWN && select_roi) {
				picked_point = get_picked_point(p);
				if (picked_point != NO_POINT) {
					roi[picked_point] = true;
					//prevent duplicates:
					if (vBoundaryInd.size() == 0)
					{
						vBoundary.push_back(vec2(points[picked_point][0], points[picked_point][1]));
						vBoundaryInd.push_back(picked_point);
					}
					else if (vBoundaryInd.back() != picked_point)
					{
						vBoundary.push_back(vec2(points[picked_point][0], points[picked_point][1]));
						vBoundaryInd.push_back(picked_point);
					}
				}
				return GL_TRUE;
			}

			return GL_FALSE;

		} else {

			if(event == GLUP_VIEWER_DOWN) {
				picked_point = get_picked_point(p);
				if (!move_vertices) { return GL_FALSE; }
				// switch(button) {
				// case 0: {
				//     if(picked_point == NO_POINT) {
				//         points.push_back(p);
				//         picked_point = points.size() - 1;
				//     }
				// } break;
				// case 1: {
				//     if(points.size() > 3) {
				//         if(picked_point != NO_POINT) {
				//             points.erase(points.begin() + int(picked_point));
				//         }
				//     }
				//     picked_point = NO_POINT;
				// } break;
				// }
				// update_Delaunay();
				// return GL_TRUE;
			}
			if(event == GLUP_VIEWER_MOVE && move_vertices) {
				if (picked_point != NO_POINT) {
					points[picked_point] = p;
					update_Delaunay();
				}
				return GL_TRUE;
			}
			if(event == GLUP_VIEWER_MOVE && select_fixed) {
				picked_point = get_picked_point(p);
				if (picked_point != NO_POINT) {
					fixed[picked_point] = true;
				}
				return GL_TRUE;
			}
			// TOBI
			if (event == GLUP_VIEWER_MOVE && select_roi) {
				picked_point = get_picked_point(p);
				if (picked_point != NO_POINT) {
					roi[picked_point] = true;
					//prevent duplicates:
					if (vBoundaryInd.size() == 0)
					{
						vBoundary.push_back(vec2(points[picked_point][0], points[picked_point][1]));
						vBoundaryInd.push_back(picked_point);
					}
					else if (vBoundaryInd.back() != picked_point)
					{
						vBoundary.push_back(vec2(points[picked_point][0], points[picked_point][1]));
						vBoundaryInd.push_back(picked_point);
					}
				}
				return GL_TRUE;
			}

			return GL_FALSE;
		}
	}

	void Lloyd_relaxation() {
		std::vector<bool> v_visited(delaunay->nb_vertices());
		std::vector<vec2> new_sites(points.size());
		Polygon cell;
		for(index_t t=0; t<delaunay->nb_cells(); ++t) {
			for(index_t lv=0; lv<3; ++lv) {
				index_t v = index_t(delaunay->cell_to_v()[3*t+lv]);
				if(!v_visited[v]) {
					v_visited[v] = true;
					get_Voronoi_cell(t,lv,cell);
					if(cell.size() > 0 && (v < fixed.size() && !fixed[v])) {
						new_sites[v] = centroid(cell);
					} else {
						new_sites[v] = points[v];
					}
				}
			}
		}
		for(index_t v=0; v<points.size(); ++v) {
			points[v] = new_sites[v];
		}
		update_Delaunay();
	}


	double inline det(const Point &u, const Point &v) {
		return imag(conj(u) * v);
	}
	// Return true iff [a,b] intersects [c,d], and store the intersection in ans
	bool intersect_segment(const Point &a, const Point &b, const Point &c, const Point &d, Point &ans) {
		const double eps = 1e-10; // small epsilon for numerical precision
		double x = det(c - a, d - c);
		double y = det(b - a, a - c);
		double z = det(b - a, d - c);
		// ab and cd are parallel ||
		if (std::abs(z) < eps || x*z < 0 || x*z > z*z || y*z < 0 || y*z > z*z) return false;
		ans = c + (d - c) * y / z;
		return true;
	}

	bool is_inside(const Polygon &poly, const Point &query) {
		Point outside(-1000, -1000);
		size_t n = poly.size();
		bool tmp, ans = false;
		for (size_t i = 0; i < poly.size(); ++i) {
			Point m; // Coordinates of intersection point
			Point p0(poly[i][0], poly[i][1]);
			Point p1(poly[(i + 1) % n][0], poly[(i + 1) % n][1]);
			tmp = intersect_segment(query, outside, p0, p1, m);
			ans = (ans != tmp);
		}
		return ans;
	}

	void sort(int &a, int &b, int &c) {
		if (a > b) {
			std::swap(a, b);
		}
		if (a > c) {
			std::swap(a, c);
		}
		if (b > c) {
			std::swap(b, c);
		}
	}

	void createTriangles()
	{
		triangles = MatrixXi::Zero(delaunay->nb_cells(), 3);
		for (index_t c = 0; c<delaunay->nb_cells(); ++c)
		{
			const signed_index_t* cell = delaunay->cell_to_v() + 3 * c;

			signed_index_t v1 = cell[0];
			signed_index_t v2 = cell[1];
			signed_index_t v3 = cell[2];

			sort(v1, v2, v3);

			triangles(c, 0) = v1;
			triangles(c, 1) = v2;
			triangles(c, 2) = v3;

		}
	}

	void removeRow(MatrixXi& matrix, unsigned int rowToRemove)
	{
		unsigned int numRows = matrix.rows()-1;
		unsigned int numCols = matrix.cols();

		if (rowToRemove < numRows)
			matrix.block(rowToRemove, 0, numRows - rowToRemove , numCols) = matrix.bottomRows(numRows - rowToRemove);

		matrix.conservativeResize(numRows, numCols);
	}

	void replaceTriangles(const MatrixXi tNew)
	{
		// find any triangles that connect to the internal of ROI
		std::vector<int> removeIdx;

		for (size_t i = 0; i < triangles.rows(); i++)
		{
			for (size_t j = 0; j < roiInternal.size(); j++)
			{
				if (roiInternal[j])
				{
					if ((triangles(i, 0) == j) || (triangles(i, 1) == j) || (triangles(i, 2) == j))
					{
						removeIdx.push_back(i);
					}
				}
			}
		}

		// using default comparison:
		std::vector<int>::iterator it;
		it = std::unique(removeIdx.begin(), removeIdx.end());
		removeIdx.resize(std::distance(removeIdx.begin(), it));

		// remove all the triangles that connect to the internal of ROI
		for (int i = removeIdx.size()-1; i >= 0; i--)
		{
			removeRow(triangles, removeIdx[i]);
		}

		// add new rows at the end of triangles
		MatrixXi tmp = triangles;
		triangles.resize(triangles.rows() + tNew.rows(),3);

		triangles << tmp, tNew;
	}

	void solve_ROI() {
		/*
		1. check if roi contains at least 6 points -> vBoundary
		2. find all vertices in polygon -> vInternal
		3. from these find the outermost and connected vertices
		   along with the connectivity -> neigh
		4. Call gurobi solver with vBoundary,vInternal,neigh
		*/

		/* 0 At this point we assume that lloyd's has finished
		 It is ok to produce the matrix triangles, if it does
		 not already exist */
		State s;
		generateQ q;
		gurobiModel g;
		if (triangles.size() == 0)
		{
			createTriangles();
		}

		// 1
		size_t nPolygon = vBoundary.size();
		nPolygon = vBoundary.size();
		std::vector<int> vInternalInd;

		if (nPolygon < 6)
		{
			return;
		}

		// 2 Find vertices in polygon
		for (int i = 0; i < points.size(); i++)
		{
			Point query(points[i][0], points[i][1]);
			bool candidate = true;
			for (size_t j = 0; j < nPolygon; j++)
			{
				if (vBoundaryInd[j] == i) { candidate = false; }

			}
			if (is_inside(vBoundary, query) && candidate)
			{
				roiInternal[i] = true;
				vInternal.push_back(points[i]);
				vInternalInd.push_back(i);
			}
		}


		// 3 Determine number of connections into the cluster to determine "neigh"
		std::vector<int> internalNeigh;
		VectorXi neigh;
		internalNeigh.resize(nPolygon);
		neigh.resize(nPolygon);

		for (int i = 0; i < vBoundaryInd.size(); i++)
		{
			internalNeigh[i] = 0;
			for (index_t c = 0; c<delaunay->nb_cells(); ++c) {
				const signed_index_t* cell = delaunay->cell_to_v() + 3 * c;
				for (index_t e = 0; e<3; ++e) {
					signed_index_t v1 = cell[e];
					signed_index_t v2 = cell[(e + 1) % 3];
					for (int j = 0; j < vInternal.size(); j++)
					{
						if ((int(v1) == vBoundaryInd[i]) && (int(v2) == vInternalInd[j]))
						{
							internalNeigh[i]++;
						}
					}

				}
			}
			neigh[i] = 6 - internalNeigh[i];
		}


		// 4 Call gubori solver
		MatrixXd vB(nPolygon,2);
		MatrixXd vI(vInternalInd.size(),2);
		for (size_t i = 0; i < nPolygon; i++)
		{
			vB(i, 0) = detected[vBoundaryInd[i]][0];
			vB(i, 1) = detected[vBoundaryInd[i]][1];
		}
		for (size_t i = 0; i < vInternalInd.size(); i++)
		{
			vI(i, 0) = detected[vInternalInd[i]][0];
			vI(i, 1) = detected[vInternalInd[i]][1];
		}

		// Generate perfect mesh in ROI
		s.init(vB, vI, neigh);
		s.fill_hole();

		// Generate adjacency matrix and the laplacian
		q.adjacencyMatrix(s.F);
		q.laplacianMatrix();

		// Deriving Q and constraints for optimization
		q.QforOptimization(s.Vperfect, s.Vdeformed, 8);
		q.optimizationConstraints(s.V_boundary.rows());

		// Generate and solve model for gurobi
		g.model(q.Q, q.Aeq);

		// Map back to indices of coordinates
		q.mapBack(g.resultX);

		// Map q.T back to global indices
		VectorXi vGlobalInd = VectorXi::Zero(vBoundaryInd.size() + vInternalInd.size());
		for (size_t i = 0; i < vBoundaryInd.size(); i++)
		{
			vGlobalInd(i) = vBoundaryInd[i];
		}
		for (size_t i = 0; i < vInternalInd.size(); i++)
		{
			vGlobalInd(i + vBoundaryInd.size()) = vInternalInd[i];
		}

		tGlobal = MatrixXi(q.T.rows(), 3);
		for (size_t i = 0; i < q.T.rows(); i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				tGlobal(i, j) = vGlobalInd(q.T(i, j));
			}

		}

		// Replace tGlobal in triangles
		replaceTriangles(tGlobal);


		// Reset ROI
		/*
		vBoundaryInd.clear();
		vInternalInd.clear();
		for (size_t i = 0; i < roi.size(); i++)
		{
			roi[i] = false;
			roiInternal[i] = false;
		}
		*/
	}

	void relaxLaplacian()
	{
		// Find vertices that have lower connectivity than 6
		VectorXi neighCount = VectorXi::Zero(points.size());

		for (size_t i = 0; i < triangles.rows(); i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				neighCount(triangles(i, j))++;
			}

		}

		// Determine fixed vertices based on connectivity and bounding box
		VectorXi indFixed = VectorXi::Zero(points.size());

		for (size_t i = 0; i < points.size(); i++)
		{
			if (neighCount(i) != 6 || fixed[i])
			{
				indFixed(i) = 1;
			}
		}

		// Rearrange coordinates and triangle indices to have free vertices first...
		//MatrixXd xRearranged = MatrixXd::Zero(points.size(), 1);
		SparseVector<double> xRearranged(points.size());
		SparseVector<double> yRearranged(points.size());
		VectorXi indicesMapping = VectorXi::Zero(points.size());
		MatrixXi trianglesRearranged = MatrixXi::Zero(triangles.rows(), 3);


		int c = 0;
		for (size_t i = 0; i < points.size(); i++)
		{
			if (indFixed(i) == 0)
			{
				indicesMapping(i) = c;
				xRearranged.fill(c) = points[i][0];
				yRearranged.fill(c) = points[i][1];

				c++;
			}
		}
		int nFree = c;

		// and then the fixed vertices
		for (size_t i = 0; i < points.size(); i++)
		{
			if (indFixed(i) == 1)
			{
				indicesMapping(i) = c;
				xRearranged.fill(c) = points[i][0];
				yRearranged.fill(c) = points[i][1];

				c++;
			}
		}

		// rearrange triangles such that it is congruent with the newly arranged coordinates
		for (size_t j = 0; j < triangles.rows(); j++)
		{
			for (size_t k = 0; k < 3; k++)
			{
				trianglesRearranged(j, k) = indicesMapping(triangles(j, k));
			}
		}

		// Generate Laplacian
		VectorXd diag = VectorXd::Zero(points.size());
		SparseMatrix<double> L(points.size(), points.size());
		typedef Triplet<int> Trip;
		std::vector< Trip > tripletList;
		tripletList.reserve(points.size() * 7);

		for (size_t i = 0; i < trianglesRearranged.rows(); i++)
		{
			for (size_t j = 0; j < 3; j++)
			{
				tripletList.push_back(Trip(trianglesRearranged(i, j), trianglesRearranged(i, (j + 1) % 3), -1));
				tripletList.push_back(Trip(trianglesRearranged(i, (j + 1) % 3), trianglesRearranged(i, j), -1));
				diag(trianglesRearranged(i, j))++;
			}
		}

		for (size_t i = 0; i < diag.rows(); i++)
		{
			tripletList.push_back(Trip(i,i, diag(i)));
		}

		L.setFromTriplets(tripletList.begin(), tripletList.end());

		// Force all non-zeros to be one
		for (int k = 0; k<L.outerSize(); ++k)
		{
			for (SparseMatrix<double>::InnerIterator it(L, k); it; ++it)
			{
				if (it.value() < 0)
				{
					L.coeffRef(it.row(), it.col()) = -1;
				}
			}
		}

		// Solve for xNew and yNew
		SparseMatrix<double>  Li = L.block(0, 0, nFree, nFree);
		SparseMatrix<double>  Lb = L.block(0, nFree, nFree, L.rows() - nFree);

		//SparseMatrix<int>  Lii = Li.inverse();
		SimplicialLDLT<SparseMatrix<double> > solver;
		solver.compute(Li);
		if (solver.info() != Success) {
			// decomposition failed
			return;
		}

		SparseVector<double> xB(L.rows() - nFree);
		SparseVector<double> yB(L.rows() - nFree);
		for (size_t i = 0; i < L.rows() - nFree; i++)
		{
			xB.fill(i) = xRearranged.coeffRef(nFree + i);
			yB.fill(i) = yRearranged.coeffRef(nFree + i);
		}
		SparseVector<double> tmp = -Lb*xB;

		VectorXd xNew = solver.solve(VectorXd(tmp));

		tmp = -Lb*yB;
		VectorXd yNew = solver.solve(VectorXd(tmp));

		// Overwrite the free vertices of xRearranged and yRearranged
		for (size_t i = 0; i < nFree; i++)
		{
			xRearranged.coeffRef(i) = xNew(i);
			yRearranged.coeffRef(i) = yNew(i);
		}

		// use indices mapping to overwrite "points" with the newly calculated coordinates
		for (size_t i = 0; i < indicesMapping.rows(); i++)
		{
			points[i][0] = xRearranged.coeffRef(indicesMapping(i));
			points[i][1] = yRearranged.coeffRef(indicesMapping(i));
		}

		// switch to show newly calculated mesh
		showCellogramDelaunay = false;
	}

}

/*********************************************************************/

namespace {
	using namespace GEO;

	// http://astronomy.swin.edu.au/~pbourke/geometry/polyarea/
	double signed_area(const Polygon& P) {
		double result = 0 ;
		for(unsigned int i=0; i<P.size(); i++) {
			unsigned int j = (i+1) % (int) P.size() ;
			const vec2& t1 = P[i] ;
			const vec2& t2 = P[j] ;
			result += t1.x * t2.y - t2.x * t1.y ;
		}
		result /= 2.0 ;
		return result ;
	}

	// http://astronomy.swin.edu.au/~pbourke/geometry/polyarea/
	vec2 centroid(const Polygon& P) {
		geo_assert(P.size() > 0) ;

		double A = signed_area(P) ;

		if(::fabs(A) < 1e-30) {
			return P[0] ;
		}

		double x = 0.0 ;
		double y = 0.0 ;
		for(unsigned int i=0; i<P.size(); i++) {
			unsigned int j = (i+1) % (int) P.size() ;
			const vec2& t1 = P[i] ;
			const vec2& t2 = P[j] ;
			double d = (t1.x * t2.y - t2.x * t1.y) ;
			x += (t1.x + t2.x) * d ;
			y += (t1.y + t2.y) * d ;
		}

		return vec2(
			x / (6.0 * A),
			y / (6.0 * A)
		) ;
	}

	static inline Sign point_is_in_half_plane(
		const vec2& p, const vec2& q1, const vec2& q2
	) {
		return PCK::orient_2d(q1, q2, p);
	}

	static inline bool intersect_segments(
		const vec2& p1, const vec2& p2,
		const vec2& q1, const vec2& q2,
		vec2& result
	) {

		vec2 Vp = p2 - p1;
		vec2 Vq = q2 - q1;
		vec2 pq = q1 - p1;

		double a =  Vp.x;
		double b = -Vq.x;
		double c =  Vp.y;
		double d = -Vq.y;

		double delta = a*d-b*c;
		if(delta == 0.0) {
			return false ;
		}

		double tp = (d * pq.x -b * pq.y) / delta;

		result = vec2(
			(1.0 - tp) * p1.x + tp * p2.x,
			(1.0 - tp) * p1.y + tp * p2.y
		);

		return true;
	}

	void clip_polygon_by_half_plane(
		const Polygon& P,
		const vec2& q1,
		const vec2& q2,
		Polygon& result
	) {
		result.clear() ;

		if(P.size() == 0) {
			return ;
		}

		if(P.size() == 1) {
			if(point_is_in_half_plane(P[0], q1, q2)) {
				result.push_back(P[0]) ;
			}
			return ;
		}

		vec2 prev_p = P[P.size() - 1] ;
		Sign prev_status = point_is_in_half_plane(
			prev_p, q1, q2
		);

		for(unsigned int i=0; i<P.size(); i++) {
			vec2 p = P[i] ;
			Sign status = point_is_in_half_plane(
				p, q1, q2
			);
			if(
				status != prev_status &&
				status != ZERO &&
				prev_status != ZERO
			) {
				vec2 intersect ;
				if(intersect_segments(prev_p, p, q1, q2, intersect)) {
					result.push_back(intersect) ;
				}
			}

			switch(status) {
			case NEGATIVE:
				break ;
			case ZERO:
				result.push_back(p) ;
				break ;
			case POSITIVE:
				result.push_back(p) ;
				break ;
			default:
				break;
			}

			prev_p = p ;
			prev_status = status ;
		}
	}

	void convex_clip_polygon(
		const Polygon& P, const Polygon& clip, Polygon& result
	) {
		Polygon tmp1 = P ;
		Polygon tmp2 ;
		Polygon* src = &tmp1 ;
		Polygon* dst = &tmp2 ;
		for(unsigned int i=0; i<clip.size(); i++) {
			unsigned int j = ((i+1) % (int) clip.size()) ;
			const vec2& p1 = clip[i] ;
			const vec2& p2 = clip[j] ;
			clip_polygon_by_half_plane(*src, p1, p2, *dst);
			geo_swap(src, dst) ;
		}
		result = *src ;
	}
}


/*********************************************************************/

int main(int argc, char** argv) {
	GEO::initialize();
	GEO::Logger::instance()->set_quiet(false);

	GEO::CmdLine::import_arg_group("standard");
	GEO::CmdLine::import_arg_group("algo");
	GEO::CmdLine::import_arg_group("gfx");

	GEO::CmdLine::set_arg("sys:assert","abort");

	// Parse command line options and filenames
	if (!GEO::CmdLine::parse(argc, argv, filenames, "<detected.xyz> <reference.xyz> <fixed>")) {
		return 1;
	}

	if (filenames.size() == 1) {
		filenames.push_back(filenames.front());
	} else if (filenames.empty()) {
		std::string input_filename = DATA_DIR "1.xyz";  // 4 does not work
		filenames.push_back(input_filename);
		filenames.push_back(input_filename);
	}

	char title[] = "[Float] Geogram Delaunay2d";
	glup_viewer_set_window_title(title);

	glup_viewer_set_screen_size(1024,800);

	glup_viewer_set_init_func(init);
	glup_viewer_set_display_func(display);
	glup_viewer_set_overlay_func(overlay);
	glup_viewer_set_mouse_func(mouse);
	glup_viewer_add_key_func(
		'k', Lloyd_relaxation, "One iteration of Lloyd relaxation"
	);
	glup_viewer_add_toggle(
		'a', glup_viewer_is_enabled_ptr(GLUP_VIEWER_IDLE_REDRAW), "Animation"
	);
	glup_viewer_add_toggle(
		'm', &move_vertices, "Toggle move selected vertices"
	);

	if(GEO::CmdLine::get_arg_bool("gfx:full_screen")) {
	   glup_viewer_enable(GLUP_VIEWER_FULL_SCREEN);
	}

	glup_viewer_main_loop(argc, argv);

	return 0;
}
