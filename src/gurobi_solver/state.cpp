#include "state.h"
#include <vector>
#include <iostream>
#include <igl/adjacency_matrix.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/sum.h>
#include <igl/diag.h>
//#include <igl/timer.h>
#include <igl/writeOBJ.h>

using namespace std;

// IO

MatrixXd readVertices(string s)
{
  vector<double> temp;
  ifstream f(s);
  double t;
  while (f >> t)
    temp.push_back(t);
  f.close();

  MatrixXd r(temp.size() / 2, 2);
  for (unsigned i = 0; i < r.rows(); ++i)
    r.row(i) << temp[i * 2 + 0], temp[i * 2 + 1];

  return r;
}

VectorXi readIntegers(string s)
{
  vector<int> temp;
  ifstream f(s);
  int t;
  while (f >> t)
    temp.push_back(t);
  f.close();

  VectorXi r(temp.size(), 1);
  for (unsigned i = 0; i < r.rows(); ++i)
    r.row(i) << temp[i];

  return r;
}

void State::load(std::string path)
{
  V_boundary = readVertices(path + "/V_border.txt");
  V_internal = readVertices(path + "/V_internals.txt");
  neigh = readIntegers(path + "/neigh.txt");

  // Merge V_boundary and V_internal to contain all original coordinates in a matrix
  Vdeformed = MatrixXd(V_boundary.rows() + V_internal.rows(), 2);
  Vdeformed << V_boundary, V_internal;
}

void State::init(const MatrixXd &vB, const MatrixXd &vI, const VectorXi &n)
{
	V_boundary = vB;
	V_internal = vI;
	neigh = n;

	// Merge V_boundary and V_internal to contain all original coordinates in a matrix
	Vdeformed = MatrixXd(V_boundary.rows() + V_internal.rows(), 2);
	Vdeformed << V_boundary, V_internal;
	/*
	
	V_boundary = vB;
	V_internal = vI;
	neigh = n;

	// Merge V_boundary and V_internal to contain all original coordinates in a matrix
	Vdeformed = MatrixXd::Zero(V_boundary.rows() + V_internal.rows(), 3);

	int k = 0;
	for (int i = 0; i < vB.rows(); i++)
	{
		Vdeformed(i, 0) = vB(i, 0);
		Vdeformed(i, 1) = vB(i, 1);
	}
	for (int i = 0; i < vI.rows(); i++)
	{
		Vdeformed(vB.rows() + i, 0) = vI(i, 0);
		Vdeformed(vB.rows() + i, 1) = vI(i, 1);
	}*/
}

void State::save(std::string path)
{
  igl::writeOBJ(path + "/out.obj",V,F);
}

// Phase 1

void State::reverse_input()
{
  VectorXi neigh2 = neigh;
  MatrixXd V_boundary2 = V_boundary;

  for (unsigned i=0;i<neigh.size();++i)
  {
    neigh(i) = neigh2(neigh.size() - i - 1);
    V_boundary.row(i) = V_boundary2.row(neigh.size() - i - 1);
  }
}

void State::create_initial_triangulation()
{
  MatrixXi pairs = find_pairs_by_tracing(neigh);
  MatrixXi raster = rasterize_pairs(pairs);

//  cerr << "Raster:" << endl;
//  cerr << raster << endl;

  MatrixXi raster_filled = fill(raster);

//  cerr << "Raster_filled:" << endl;
//  cerr << raster_filled << endl;

   extract_mesh_from_raster(raster_filled, V, F);
}

MatrixXi State::find_pairs_by_tracing(const VectorXi &neigh)
{
  VectorXi turns = neigh.array() - 4;

  MatrixXi dirs_even(6, 2);
  MatrixXi dirs_odd(6, 2);
  dirs_even << 1, 0,
               0, 1,
               -1, 1,
               -1, 0,
               -1, -1,
               0, -1;
  dirs_odd << 1, 0,
              1, 1,
              0, 1,
              -1, 0,
              0, -1,
              1, -1;

  int dir = starting_direction;

  vector<Vector2i> pairs;

  Vector2i current_pair;
  current_pair << 0, 0;
  pairs.push_back(current_pair);

  // Follow the boundary
  for (unsigned i = 0; i < turns.size(); ++i)
  {
    // move in the current direction
    current_pair += (current_pair(1) % 2) == 0 ? dirs_even.row(dir).transpose() : dirs_odd.row(dir).transpose();

    // save the current position
    pairs.push_back(current_pair);

    //update the direction
    dir = (dir + dirs_even.rows() + turns((i + 1) % turns.size())) % dirs_even.rows();
  }

  MatrixXi ret;

  if (pairs.front() == pairs.back() && starting_direction == dir)
  {
    // Pack it before returning the pairs
    ret.resize(pairs.size(), 2);
    for (unsigned i = 0; i < ret.rows(); ++i)
      ret.row(i) = pairs[i];
  }

  return ret;
}

MatrixXi State::rasterize_pairs(const MatrixXi &pairs)
{
  // find the size of the grid
  Vector2i min = pairs.colwise().minCoeff();
  Vector2i max = pairs.colwise().maxCoeff();

  min = min.array() - 1;
  max = max.array() + 1;

  // very important, rasterize in a grid whose even and odd rows correspond with the convenction in diagram.png
  if (((min(1)+100000) % 2) == 1)
    min(1) = min(1) - 1;
  assert(((min(1)+100000) % 2) == 0);

  MatrixXi raster = MatrixXi::Constant(max(0) - min(0) + 1, max(1) - min(1) + 1, -2);

  // rasterize the pairs
  for (unsigned i = 0; i < pairs.rows() - 1; ++i)
    raster(pairs(i, 0) - min(0), pairs(i, 1) - min(1)) = i;

  return raster;
}


bool fillifvalid(MatrixXi& raster, int i, int j)
{
  if (i>=0 && i<raster.rows() && j >=0 && j <raster.cols() && raster(i,j) == -2)
  {
    raster(i, j) = -3;
    return true;
  }
  return false;
}

MatrixXi State::fill(const MatrixXi &raster)
{

//  MatrixXi raster_filled = raster;

  // first give a different color to the ones outside
  bool done = false;
  MatrixXi raster_filled = raster;

  raster_filled(0,0) = -3; // first pixel is outside

  while(!done)
  {
    done = true;
    for (unsigned i = 0; i < raster.rows(); ++i)
    {
      for (unsigned j = 0; j < raster.cols(); ++j)
      {
        if (raster_filled(i,j) == -3)
        {
          if ((j%2) == oddeven ? 0 : 1)
          {
            if (
                    fillifvalid(raster_filled,i,j-1) ||
                    fillifvalid(raster_filled,i+1,j-1)
                    ||
                    fillifvalid(raster_filled,i-1,j) ||
                    fillifvalid(raster_filled,i+1,j)
                    ||
                    fillifvalid(raster_filled,i,j+1) ||
                    fillifvalid(raster_filled,i+1,j+1)
                    )
              done = false;

          }
          else
          {
            if (
                    fillifvalid(raster_filled,i-1,j-1) ||
                    fillifvalid(raster_filled,i,j-1)
                    ||
                    fillifvalid(raster_filled,i-1,j) ||
                    fillifvalid(raster_filled,i+1,j)
                    ||
                    fillifvalid(raster_filled,i-1,j+1) ||
                    fillifvalid(raster_filled,i,j+1)
                    )
              done = false;
          }
        }
      }
    }

  }

  int count = V_boundary.rows();

  for (unsigned i = 0; i < raster_filled.rows(); ++i)
  {
    for (unsigned j = 0; j < raster_filled.cols(); ++j)
    {
      if (raster_filled(i,j) == -2)
        raster_filled(i,j) = count++;
    }
  }


  return raster_filled;
}

void State::extract_mesh_from_raster(const MatrixXi &raster, MatrixXd &V, MatrixXi &F)
{

  // Count the vertices
  int count = raster.maxCoeff() + 1;
  // Create the vertices
  //cout << V << endl;
  V.resize(count, 3);

  count = 0;
  for (unsigned i = 0; i < raster.rows(); i++)
    for (unsigned j = 0; j < raster.cols(); j++)
      if (raster(i, j) != -3)
        V.row(raster(i, j)) << (i % 2 == 0 ? i * 2 : i * 2 + 1), j * 2, 0;

  vector<Vector3i> Fv;

  for (unsigned i = 0; i < raster.rows() - 1; i++)
  {
    for (unsigned j = 0; j < raster.cols() - 1; j++)
    {
      if ((j%2) == (oddeven ? 0 : 1))
      {
        // odd
        if ((raster(i, j) != -3) && (raster(i + 1, j+1) != -3) && (raster(i, j + 1) != -3))
        {
          Vector3i v;
          v << raster(i, j), raster(i + 1, j+1), raster(i, j + 1);
          Fv.push_back(v);
        }
        if ((raster(i, j) != -3) && (raster(i + 1,j) != -3) && (raster(i + 1, j + 1) != -3))
        {
          Vector3i v;
          v << raster(i, j), raster(i + 1, j), raster(i + 1, j + 1);
          Fv.push_back(v);
        }
      }
      else
      {
        if (i > 0 && (raster(i, j) != -3) && (raster(i, j+1) != -3) && (raster(i-1, j + 1) != -3))
        {
          Vector3i v;
          v << raster(i, j), raster(i, j+1), raster(i-1 , j + 1);
          Fv.push_back(v);
        }
        if ((raster(i, j) != -3) && (raster(i + 1,j) != -3) && (raster(i, j + 1) != -3))
        {
          Vector3i v;
          v << raster(i, j), raster(i + 1, j), raster(i, j + 1);
          Fv.push_back(v);
        }
      }

    }
  }

  // Copy in F;
  F.resize(Fv.size(), 3);
  for (unsigned i = 0; i < Fv.size(); ++i)
    F.row(i) << Fv[i].transpose();

}

// Phase 2
void State::fit_triangulation()
{
  // Allocate fixed
  fixed = VectorXi::Constant(V.rows(), 0);
  fixed_dotid = VectorXi::Constant(V.rows(), -2);

  // Allocate fixed_vertices
  fixed_dots = VectorXi::Constant(V_internal.rows(), -1);

  // Copy the vertex positions and fix the boundary
  for (unsigned i = 0; i < V_boundary.rows(); ++i)
  {
    V.row(i) << V_boundary.row(i), 0;
    fixed(i) = 1;
    fixed_dotid(i) = -1;
  }

  // Compute the new positions
  V = compute_new_positions(V, F, fixed);
  Vperfect = V;
}

MatrixXd State::compute_new_positions(const MatrixXd &V, const MatrixXi &F, const MatrixXi &fixed)
{
  // Compute the uniform laplacian

  Eigen::SparseMatrix<double> A;
  igl::adjacency_matrix(F, A);
  // sum each row
  SparseVector<double> Asum;
  igl::sum(A, 1, Asum);
  // Convert row sums into diagonal of sparse matrix
  SparseMatrix<double> Adiag;
  igl::diag(Asum, Adiag);
  // Build uniform laplacian
  SparseMatrix<double> L;
  L = A - Adiag;

  VectorXi fixed_index = VectorXi::Constant(fixed.sum(), 0);
  MatrixXd fixed_value = MatrixXd::Constant(fixed.sum(), 3, 0);

  int count = 0;
  for (unsigned i = 0; i < fixed.size(); ++i)
  {
    if (fixed(i) == 1)
    {
      fixed_index(count) = i;
      fixed_value.row(count++) = V.row(i);
    }
  }

  // If everything is constrained, dont solve
  if (fixed_index.size() == V.rows())
    return V;

  // Solve with boundary fixed
  igl::min_quad_with_fixed_data<double> data;
  MatrixXd NV;

  // Linear term is 0
  VectorXd B = VectorXd::Zero(V.rows(), 1);
  // Empty constraints
  VectorXd Beq;
  SparseMatrix<double> Aeq;

  //cerr << fixed_index << endl;

  // solve
  igl::min_quad_with_fixed_precompute(L, fixed_index, Aeq, false, data);
  igl::min_quad_with_fixed_solve(data, B, fixed_value, Beq, NV);

  // Copy back the result
  return NV;
}

// Phase 3
void State::snap_closest()
{
  // For every internal point, try to snap it to all possible dots
  // and create a list ordered by distance

  //cerr << to_snap.size() << endl;

  if (to_snap.size() == 0)
  {
    // add all free vertices to the snap list
    if (!snap_border)
    {
      for (unsigned i = 0; i < fixed.size(); ++i)
        if (fixed(i) == 0)
          to_snap.insert(i);
    }
    else
    {
      // add only the ones in triangles where two vertices are already fixed
      VectorXi border = VectorXi::Constant(fixed.size(), 0);
      for (unsigned i = 0; i < F.rows(); ++i)
      {
        if (fixed(F(i, 0)) + fixed(F(i, 1)) + fixed(F(i, 2)) >= 2)
        {
          border(F(i, 0)) = 1;
          border(F(i, 1)) = 1;
          border(F(i, 2)) = 1;
        }
      }

      // use a temporary vector
      vector<int> t;
      for (unsigned i = 0; i < fixed.size(); ++i)
        if (fixed(i) == 0 && border(i) == 1)
          to_snap.insert(i);
    }
  }

  // collect ids of internal vertices
  VectorXi free_internal_ids(to_snap.size());
  int count = 0;
  for(auto it = to_snap.begin(); it != to_snap.end(); ++it)
    free_internal_ids(count++) = *it;

  // collect ids of unused dots
  count = 0;
  for (unsigned i=0;i<fixed_dots.size();++i)
    if (fixed_dots(i) == -1)
      ++count;

  VectorXi free_dots_ids(count);
  count = 0;

  for (unsigned i=0;i<fixed_dots.size();++i)
    if (fixed_dots(i) == -1)
      free_dots_ids(count++) = i;

  // if there are no vertices, give up
  if (free_internal_ids.size() == 0)
    return;


  // produce all the possible pairs with their distances
  vector<Vector3d> q; // Every item contains: distance, id vertex mesh, id dots

  for (unsigned i=0; i<free_internal_ids.size(); ++i)
  {
    for (unsigned j=0; j<free_dots_ids.size(); ++j)
    {
      Vector3d t;
      Vector3d p_mesh = V.row(free_internal_ids(i));

      Vector3d p_dot; p_dot << V_internal.row(free_dots_ids(j)).transpose(),0;

      t << (p_mesh - p_dot).squaredNorm(), free_internal_ids(i), free_dots_ids(j);
      q.push_back(t);
    }
  }

  // sort q by distance
  std::sort(q.begin(), q.end(), [](const Vector3d& a, const Vector3d& b)
  {
      return a(0) < b(0);
  });

  if (check_contained)
  {
    int current = 0;
    MatrixXd V_temp;
    while (current < q.size())
    {
      //cerr << current << "/" << q.size() << endl;
      // round one
      Vector3d selected = q[current];
      int dot_id = round(selected(2));
      int internal_id = round(selected(1));

      assert(fixed_dots(dot_id) == -1);

      V.row(internal_id) << V_internal.row(dot_id), 0;
      fixed(internal_id) = 1;
      fixed_dotid(internal_id) = dot_id;
      fixed_dots(dot_id) = internal_id;

      // Compute the new positions
      V_temp = compute_new_positions(V, F, fixed);

      // if it is valid go ahead
      if (is_mesh_valid(V_temp,F,fixed,V_internal,fixed_dots))
      {
        //cerr << "Check passed, mesh is valid!" << endl;
        // remove it from the to_snap and update the temp data
        to_snap.erase(internal_id);
        V = V_temp;
        // done, return
        return;
      }
      else
      {
        //cerr << "Check failed, mesh is invalid!" << endl;
        fixed(internal_id) = 0;
        fixed_dotid(internal_id) = -2;
        fixed_dots(dot_id) = -1;
      }
      current++;

    }
  }

  // if it arrives here, no smart selection is possible, round the first one!

  // round the first one (for now blindly)
  Vector3d selected = q[0];
  int dot_id = round(selected(2));
  int internal_id = round(selected(1));

  // remove it from the to_snap
  to_snap.erase(internal_id);

  assert(fixed_dots(dot_id) == -1);
  fixed_dots(dot_id) = internal_id;

  V.row(internal_id) << V_internal.row(dot_id), 0;
  fixed(internal_id) = 1;
  fixed_dotid(internal_id) = dot_id;

  cerr << "V.rows(): " << V.rows() << endl;
  cerr << "fixed.size(): " << fixed.size() << endl;
  cerr << "fixed: " << fixed << endl;

  // Compute the new positions
  V = compute_new_positions(V, F, fixed);

}

//bool in_triangle(const Vector3d& a, const Vector3d& b, const Vector3d& c, const Vector3d& p)
//{
//  return (b-a).cross(p-a).squaredNorm() > 0 && (c-b).cross(p-b).squaredNorm() > 0 && (a-c).cross(p-c).squaredNorm() > 0;
//}

float sign (const Vector2d& p1, const Vector3d& p2, const Vector3d& p3)
{
  return (p1(0) - p3(0)) * (p2(1) - p3(1)) - (p2(0) - p3(0)) * (p1(1) - p3(1));
}

bool in_triangle (const Vector3d& v1, const Vector3d& v2, const Vector3d& v3, const Vector2d& pt)
{
  bool b1, b2, b3;

  b1 = sign(pt, v1, v2) < 0.0f;
  b2 = sign(pt, v2, v3) < 0.0f;
  b3 = sign(pt, v3, v1) < 0.0f;

  return ((b1 == b2) && (b2 == b3));
}

bool State::is_mesh_valid(const MatrixXd& V, const MatrixXi& F, const VectorXi& fixed, const MatrixXd& P, const VectorXi& fixed_P)
{
  // check for overlaps
  for (unsigned i=0;i<F.rows();++i)
  {
    if (fixed(F(i,0)) == 1 && fixed(F(i,1)) == 1 && fixed(F(i,2)) == 1)
    {
      //cerr << "checking face: " << i << endl;
      for (unsigned j=0; j<P.rows(); ++j)
      {
        //cerr << "checking vertex: " << j << endl;

        if (fixed_P(j) == -1 && in_triangle(V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)),P.row(j)))
          return false;
      }
    }
  }

  // check for multiple connected components
  if (avoid_connected_components)
  {
    // copy fixed
    VectorXi fixedc = fixed;

    int seed = -1;
    // find an unfixed vertex
    for (unsigned i=0; i<fixedc.size(); ++i)
    {
      if (fixedc(i) == 0)
      {
        seed = i;
        break;
      }
    }

    if (seed != -1)
    {
      fixedc(seed) = -1;
      // expand until possible
      bool done = false;
      while (!done)
      {
        done = true;
        for (unsigned j = 0; j < F.rows(); ++j)
        {
          if ((fixedc(F(j, 0)) == -1) || (fixedc(F(j, 1)) == -1) || (fixedc(F(j, 2)) == -1))
          {
            if (fixedc(F(j, 0)) == 0)
            {
              fixedc(F(j, 0)) = -1;
              done = false;
            }
            if (fixedc(F(j, 1)) == 0)
            {
              fixedc(F(j, 1)) = -1;
              done = false;
            }
            if (fixedc(F(j, 2)) == 0)
            {
              fixedc(F(j, 2)) = -1;
              done = false;
            }
          }
        }
      }

      // check of there are zeros left
      for (unsigned i = 0; i < fixedc.size(); ++i)
      {
        if (fixedc(i) == 0)
        {
          cerr << "connected component failed!" << endl;
          return false;
        }
      }
    }
  }

  return true;
}

void State::unfreeze_overlapping()
{
  // find all vertices overlapping with any triangle (excluding its own vertices)
  vector<int> overlapping_vertices;

  for (unsigned i=0;i<F.rows();++i)
  {
    for (unsigned j=0; j<V.rows(); ++j)
    {
      if (in_triangle(V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)),V.row(j).segment(0,2)))
      {
        // check if it is one of the corners
        // check if it is one of the corners
        if (!((F(i,0) == j) || (F(i,1) == j) || (F(i,2) == j)))
          overlapping_vertices.push_back(j);
      }
    }
  }

  to_snap.clear();

//  overlapping_vertices.push_back(197);
//  overlapping_vertices.push_back(203);
//  overlapping_vertices.push_back(104);
//  overlapping_vertices.push_back(114);
//  overlapping_vertices.push_back(113);

  // Expand for n-rings (for now 1)
  VectorXi mark = VectorXi::Constant(V.rows(),0);
  for (unsigned i=0;i<overlapping_vertices.size();++i)
    mark(overlapping_vertices[i]) = 1;

  VectorXi markcopy = mark;

  for(unsigned i=0; i<1; ++i)
  {
    for(unsigned j=0;j<F.rows();++j)
    {
      if ((mark(F(j,0)) == 1) || (mark(F(j,1)) == 1) || (mark(F(j,2)) == 1))
      {
        markcopy(F(j,0)) = 1;
        markcopy(F(j,1)) = 1;
        markcopy(F(j,2)) = 1;
      }
    }
  }

  mark = markcopy;

  overlapping_vertices.clear();
  for (unsigned i=0;i<mark.size();++i)
    if (mark(i) == 1)
      overlapping_vertices.push_back(i);

//  cerr << "Overlapping vertices: " << endl;

  for (unsigned i=0; i<overlapping_vertices.size(); ++i)
  {
    int internal_id = overlapping_vertices[i];
    int dot_id = fixed_dotid(internal_id);

    if ((dot_id == -1) || (dot_id == -2)) // boundary, do not unfreeze
      continue;

    fixed(internal_id) = 0;
    fixed_dotid(internal_id) = -2;
    fixed_dots(dot_id) = -1;

//    cerr << overlapping_vertices[i] << " ";

  }
//  cerr << endl;

  // Compute the new positions
  V = compute_new_positions(V, F, fixed);



}

// All together
void State::fill_hole()
{
//  igl::Timer timer;
//  timer.start();

  create_initial_triangulation();
  
  fit_triangulation();

  /* commented by Tobi - snapping is pointless if solved by gurobi
  for (unsigned i = 0; i < V_internal.rows(); ++i)

    snap_closest();
	*/
//  timer.stop();

//  cerr << "Total time (ms):" << timer.getElapsedTimeInMilliSec() << endl;

}
