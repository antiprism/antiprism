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

bool make_geodesic_planar(Geometry *geom, const Geometry &base, int m, int n)
{
  if (m < 0 || n < 0 || (m == 0 && n == 0))
    return false; // invalid pattern
  Geodesic geod(base, m, n, 'p');
  geod.make_geo(*geom);
  return true; // valid pattern
}

bool make_geodesic_sphere(Geometry *geom, const Geometry &base, int m, int n,
                          Vec3d cent)
{
  if (m < 0 || n < 0 || (m == 0 && n == 0))
    return false; // invalid pattern
  Geodesic geod(base, m, n, 's', cent);
  geod.make_geo(*geom);
  return true; // valid pattern
}

void project_onto_sphere(Geometry *geom, Vec3d centre, double radius)
{
  for (Vec3d &v : geom->raw_verts())
    v = (v - centre).with_len(radius);
}

// RK - test points versus hull functions

bool test_points_vs_hull(const vector<Vec3d> &P, const Geometry &hull,
                         const bool &inside, const bool &surface,
                         const bool &outside, const double &eps)
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

// RK - functions for winding number

// put faces numbers in face_idxs into fgeom
Geometry faces_to_geom(const Geometry &geom, const vector<int> &face_idxs)
{
  Geometry fgeom;
  fgeom.add_verts(geom.verts());
  fgeom.colors(VERTS) = geom.colors(VERTS);
  for (int j : face_idxs) {
    fgeom.add_face(geom.faces()[j], geom.colors(FACES).get(j));
  }
  fgeom.del(VERTS, fgeom.get_info().get_free_verts());
  return fgeom;
}

// Copyright 2001, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.

//    a Point is defined by its coordinates {int x, y;}
//===================================================================

// isLeft(): tests if a point is Left|On|Right of an infinite line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 and P1
//            =0 for P2 on the line
//            <0 for P2 right of the line
//    See: the January 2001 Algorithm "Area of 2D and 3D Triangles and Polygons"

// inline int
// isLeft( Point P0, Point P1, Point P2 )

// RK - this works with wn_PnPoly better if a double is passed back

static double isLeft(const Vec3d &P0, const Vec3d &P1, const Vec3d &P2,
                     const int &idx)
{
  int idx1 = (idx + 1) % 3;
  int idx2 = (idx + 2) % 3;

  double P0_x = P0[idx1];
  double P0_y = P0[idx2];

  double P1_x = P1[idx1];
  double P1_y = P1[idx2];

  double P2_x = P2[idx1];
  double P2_y = P2[idx2];

  return ((P1_x - P0_x) * (P2_y - P0_y) - (P2_x - P0_x) * (P1_y - P0_y));
}

//===================================================================

// cn_PnPoly(): crossing number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  0 = outside, 1 = inside
// This code is patterned after [Franklin, 2000]

// bool
// cn_PnPoly( Point P, Point* V, int n )

// RK - try to use same variable names as original pnpoly

/* RK - Not currently used. none of the tests consider epsilon; if used
recommend using functions

bool cn_PnPoly(const Geometry &polygon, const Vec3d &P, const int &idx, int
&crossing_number, double eps)
{
   const vector<Vec3d> &verts = polygon.verts();

   int idx1 = (idx+1)%3;
   int idx2 = (idx+2)%3;

   double testx = P[idx1]; // P.x
   double testy = P[idx2]; // P.y

   const vector<int> &face = polygon.faces()[0];
   int n = face.size();

   int cn = 0;    // the crossing number counter

   // loop through all edges of the polygon
   for (int i=0; i<n; i++) {    // edge from face[i] to face[i+1]
      int j = (i+1)%n;

      double vertx_i = verts[face[i]][idx1]; // V[i].x
      double verty_i = verts[face[i]][idx2]; // V[i].y
      double vertx_j = verts[face[j]][idx1]; // V[i+1].x
      double verty_j = verts[face[j]][idx2]; // V[i+1].y

      if (((verty_i <= testy) && (verty_j > testy))          // an upward
crossing
       || ((verty_i > testy) && (verty_j <= testy))) {       // a downward
crossing
         // compute the actual edge-ray intersect x-coordinate
         double vt = (testy - verty_i) / (verty_j - verty_i);
         if (testx < vertx_i + vt * (vertx_j - vertx_i))     // testx <
intersect
            ++cn;   // a valid crossing of y=testy right of testx
       }
   }

   crossing_number = cn;

   //return (cn&1);    // 0 if even (out), and 1 if odd (in)
   return (cn&1 ? true : false); // false if even (out), and true if odd (in)
}
*/
//===================================================================

// wn_PnPoly(): winding number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  wn = the winding number (=0 only if P is outside V[])

// bool
// wn_PnPoly( Point P, Point* V, int n )

// RK - try to use same variable names as original pnpoly
// RK - there is proof that there are mistakes when epsilon is not considered
// RK - wn is n time too large. winding_number = wn/n

static bool wn_PnPoly(const Geometry &polygon, const Vec3d &P, const int &idx,
                      int &winding_number, double eps)
{
  winding_number = 0;
  if (!P.is_set() || !polygon.faces().size())
    return false;

  const vector<Vec3d> &verts = polygon.verts();

  int idx2 = (idx + 2) % 3;

  double testy = P[idx2]; // P.y

  const vector<int> &face = polygon.faces()[0];
  int n = face.size();

  int wn = 0; // the winding number counter

  // loop through all edges of the polygon
  for (int i = 0; i < n; i++) { // edge from face[i] to face[i+1]
    int j = (i + 1) % n;

    Vec3d vert_i = verts[face[i]]; // V[i]
    Vec3d vert_j = verts[face[j]]; // V[[i+1]

    double verty_i = vert_i[idx2]; // V[i].y
    double verty_j = vert_j[idx2]; // V[i+1].y

    for (int i = 0; i < n; i++) { // edge from face[i] to face[i+1]
      // if (verty_i <= testy) {   // start y <= P.y
      if (double_le(verty_i, testy, eps)) {
        // if (verty_j > testy)     // an upward crossing
        if (double_gt(verty_j, testy, eps))
          // if (isLeft( vert_i, vert_j, P, idx ) > eps)  // P left of edge
          if (double_gt(isLeft(vert_i, vert_j, P, idx), 0, eps))
            ++wn; // have a valid up intersect
      }
      else { // start y > P.y (no test needed)
        // if (verty_j <= testy)              // a downward crossing
        if (double_le(verty_j, testy, eps))
          // if (isLeft( vert_i, vert_j, P, idx ) < -eps)  // P right of edge
          if (double_lt(isLeft(vert_i, vert_j, P, idx), 0, eps))
            --wn; // have a valid down intersect
      }
    }
  }

  if (n)
    winding_number = wn / n;

  return (wn ? true : false);
}

int get_winding_number(const Geometry &polygon, const Vec3d &point, double eps)
{
  int winding_number;
  wn_PnPoly(polygon, point, 2, winding_number, eps);
  return winding_number;
}

// the winding number is correct when the faces are placed on the xy-plane
// facing z-positive
// geom contains one polygon
// if points is empty, all serial triads of points in a face are tested and the
// highest magnitude winding number is returned
// if complex polygon has both -W and +W winding, +W is returned
int get_winding_number_polygon(const Geometry &polygon,
                               const vector<Vec3d> &points,
                               const Normal &face_normal, bool find_direction,
                               double eps)
{
  // make copies
  Geometry pgon = polygon;
  vector<Vec3d> pts = points;

  int psz = (int)pts.size();
  if (!psz) {
    const vector<Vec3d> &verts = pgon.verts();
    const vector<int> &face = pgon.faces()[0];

    vector<Vec3d> triangle;
    int fsz = (int)face.size();
    for (int i = 0; i < fsz; i++) {
      triangle.push_back(verts[face[i]]);
      triangle.push_back(verts[face[(i + 1) % fsz]]);
      triangle.push_back(verts[face[(i + 2) % fsz]]);
      pts.push_back(centroid(triangle));
      triangle.clear();
    }
    psz = (int)pts.size();
  }

  Vec3d norm = face_normal.outward().unit();
  double z = 1;
  if (find_direction) {
    double test = 0.0;
    for (int i = 2; i >= 0; i--) {
      test = norm[i];
      if (!double_eq(test, 0.0, eps))
        break;
    }
    z = (test > 0.0) ? 1 : -1;
  }

  // rotate polygon to face forward
  Trans3d trans = Trans3d::rot(norm + pgon.face_cent(0), Vec3d(0, 0, z));
  pgon.transform(trans);

  // also rotate points per matrix
  Geometry vgeom;
  ;
  for (int i = 0; i < psz; i++)
    vgeom.add_vert(pts[i]);
  vgeom.transform(trans);
  pts.clear();
  for (int i = 0; i < psz; i++)
    pts.push_back(vgeom.verts()[i]);
  vgeom.clear_all();

  int winding_number = 0;
  for (int i = 0; i < psz; i++) {
    // projection index = 2;
    int wn = 0;
    wn_PnPoly(pgon, pts[i], 2, wn, eps);
    if ((abs(wn) > abs(winding_number)) ||
        ((wn > 0) && (wn == -winding_number)))
      winding_number = wn;
  }
  return winding_number;
}

int find_polygon_denominator_signed(const Geometry &geom, int face_idx,
                                    double eps)
{
  const vector<int> &face = geom.faces()[face_idx];

  Normal face_normal(geom, face_idx, centroid(geom.verts()), eps);

  vector<int> sface_idxs;
  sface_idxs.push_back(face_idx);
  Geometry polygon = faces_to_geom(geom, sface_idxs);

  vector<Vec3d> points; // empty points to trigger test triangles
  int d = get_winding_number_polygon(polygon, points, face_normal, false, eps);

  int fsz = (int)face.size();
  if (d < 0)
    d += fsz;

  // check for condition we missed the polygon
  if (!d)
    d = 1;

  // hemispherical faces always are treated as positively wound so return them
  // signed
  if (face_normal.is_hemispherical() && d > fsz / 2)
    d = fsz - d;

  return d;
}

} // namespace anti
