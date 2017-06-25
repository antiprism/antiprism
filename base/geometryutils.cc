/*
   Copyright (c) 2010-2016, Adrian Rossiter

   Antiprism - http://www.antiprism.com

   Permission is hereby granted, free of charge, to any person obtaining a
   copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation
   the rights to use, copy, modify, merge, publish, distribute, sublicense,
   and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions:

      The above copyright notice and this permission notice shall be included
      in all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
  IN THE SOFTWARE.
*/

/*!\file geometryutils.h
   \brief utilities for geometries
*/

#include <stdio.h>
#include <stdlib.h>

#include <vector>

#include "geometry.h"
#include "geometryinfo.h"
#include "geometryutils.h"
#include "mathutils.h"
#include "private_geodesic.h"
#include "private_misc.h"

using std::vector;

namespace anti {

bool make_geodesic_planar(Geometry &geom, const Geometry &base, int m, int n)
{
  if (m < 0 || n < 0 || (m == 0 && n == 0))
    return false; // invalid pattern
  Geodesic geod(base, m, n, 'p');
  geod.make_geo(geom);
  return true; // valid pattern
}

bool make_geodesic_sphere(Geometry &geom, const Geometry &base, int m, int n,
                          Vec3d cent)
{
  if (m < 0 || n < 0 || (m == 0 && n == 0))
    return false; // invalid pattern
  Geodesic geod(base, m, n, 's', cent);
  geod.make_geo(geom);
  return true; // valid pattern
}

void project_onto_sphere(Geometry &geom, Vec3d centre, double radius)
{
  for (Vec3d &v : geom.raw_verts())
    v = (v - centre).with_len(radius);
}

// RK - test points versus hull functions

bool test_points_vs_hull(const vector<Vec3d> &P, const Geometry &hull,
                         const bool inside, const bool surface,
                         const bool outside, const double &eps)
{
  const vector<Vec3d> &verts = hull.verts();
  const vector<vector<int>> &faces = hull.faces();

  Vec3d C = centroid(verts);

  bool answer = true;
  for (const auto &face : faces) {
    Vec3d n = face_norm(verts, face).unit();
    double D = vdot(verts[face[0]] - C, n);
    // if(D < 0)
    if (double_compare(D, 0, eps) < 0) { // Make sure the normal points outwards
      D = -D;
      n = -n;
    }

    for (const auto &j : P) {
      double t = vdot(j - C, n);
      // if (t < D-eps && !inside)
      if ((double_compare(t, D, eps) < 0) && !inside)
        answer = false;
      else
          // if (t > D+eps && !outside)
          if ((double_compare(t, D, eps) > 0) && !outside)
        answer = false;
      else if (!surface)
        answer = false;

      if (!answer)
        break;
    }

    if (!answer)
      break;
  }

  return answer;
}

bool are_points_in_hull(const vector<Vec3d> &points, const Geometry &hull,
                        unsigned int inclusion_test, const double &eps)
{
  if (inclusion_test % 8 == 0 ||
      (inclusion_test & INCLUSION_IN && inclusion_test & INCLUSION_OUT))
    return false;
  else
    return test_points_vs_hull(points, hull, inclusion_test & INCLUSION_IN,
                               inclusion_test & INCLUSION_ON,
                               inclusion_test & INCLUSION_OUT, eps);
}

// RK - Various find functions for geom

int find_vert_by_coords(const Geometry &geom, const Vec3d &coords, double eps)
{
  const vector<Vec3d> &verts = geom.verts();
  int v_idx = -1;
  for (unsigned int i = 0; i < verts.size(); i++) {
    if (!compare(verts[i], coords, eps)) {
      v_idx = i;
      break;
    }
  }
  return v_idx;
}

// elem could be face or another edge
bool edge_exists_in_elem(const vector<int> &elem, const vector<int> &edge)
{
  vector<int> edge1 = make_edge(edge[0], edge[1]);

  bool found = false;

  int sz = elem.size();
  for (int i = 0; i < sz; i++) {
    vector<int> edge2 = make_edge(elem[i], elem[(i + 1) % sz]);

    if (edge1 == edge2) {
      found = true;
      break;
    }
  }

  return found;
}

bool edge_exists_in_face(const vector<int> &face, const vector<int> &edge)
{
  return edge_exists_in_elem(face, edge);
}

vector<int> find_faces_with_edge(const vector<vector<int>> &faces,
                                 const vector<int> &edge)
{
  vector<int> face_idxs;
  for (unsigned int i = 0; i < faces.size(); i++) {
    if (edge_exists_in_elem(faces[i], edge))
      face_idxs.push_back(i);
  }

  return face_idxs;
}

inline bool vertex_exists_in_elem(const vector<int> &face, int v_idx)
{
  return std::find(face.begin(), face.end(), v_idx) != face.end();
}

bool vertex_exists_in_face(const vector<int> &face, const int v_idx)
{
  return vertex_exists_in_elem(face, v_idx);
}

bool vertex_exists_in_edge(const vector<int> &edge, int v_idx)
{
  return vertex_exists_in_elem(edge, v_idx);
}

vector<int> find_elems_with_vertex(const vector<vector<int>> &elems, int v_idx)
{
  vector<int> elem_idxs;
  for (unsigned int i = 0; i < elems.size(); i++) {
    if (vertex_exists_in_elem(elems[i], v_idx))
      elem_idxs.push_back(i);
  }

  return elem_idxs;
}

vector<int> find_faces_with_vertex(const vector<vector<int>> &faces, int v_idx)
{
  return find_elems_with_vertex(faces, v_idx);
}

vector<int> find_edges_with_vertex(const vector<vector<int>> &edges, int v_idx)
{
  return find_elems_with_vertex(edges, v_idx);
}

// find index of an edge in an edge list. -1 if not found
int find_edge_in_edge_list(const vector<vector<int>> &edges,
                           const vector<int> &edge)
{
  int found = -1;
  for (unsigned int i = 0; i < edges.size(); i++) {
    if (edge_exists_in_elem(edges[i], edge)) {
      found = i;
      break;
    }
  }

  return found;
}

vector<vector<int>> find_unmatched_edges(const Geometry &geom)
{
  const vector<vector<int>> &faces = geom.faces();
  vector<vector<int>> edges;

  // can't use get_impl_edges() here because we need to know the duplicates
  for (auto face : faces) {
    int sz = face.size();
    for (int j = 0; j < sz; j++)
      edges.push_back(make_edge(face[j], face[(j + 1) % sz]));
  }

  sort(edges.begin(), edges.end());

  vector<vector<int>> unmatched_edges;

  int sz = edges.size();
  for (int i = 1; i < sz;) {
    if (edges[i - 1] == edges[i])
      i += 2;
    else {
      unmatched_edges.push_back(edges[i - 1]);
      if (i == sz - 1)
        unmatched_edges.push_back(edges[i]);
      i++;
    }
  }

  return unmatched_edges;
}

} // namespace anti
