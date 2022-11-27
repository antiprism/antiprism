/*
   Copyright (c) 2003-2022, Adrian Rossiter, Roger Kaufman

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

/*
   Name: canonical_common.cc
   Description: canonical code shared in /src
   Project: Antiprism - http://www.antiprism.com
*/

#include "canonical_common.h"
#include "../base/antiprism.h"

#include <cstdio>
#include <string>
#include <vector>

using std::string;
using std::vector;

using namespace anti;

// RK - find nearpoints radius, sets range minimum and maximum
double edge_nearpoints_radius(const Geometry &geom, double &min, double &max,
                              Vec3d &center)
{
  min = std::numeric_limits<double>::max();
  max = std::numeric_limits<double>::min();

  vector<vector<int>> edges;
  geom.get_impl_edges(edges);

  vector<Vec3d> near_pts;

  double nearpt_radius = 0;
  for (auto &edge : edges) {
    Vec3d P = geom.edge_nearpt(edge, Vec3d(0, 0, 0));
    near_pts.push_back(P);

    double l = P.len();
    nearpt_radius += l;
    if (l < min)
      min = l;
    if (l > max)
      max = l;
  }

  center = centroid(near_pts);

  return nearpt_radius / double(edges.size());
}

// RK - wrapper
double edge_nearpoints_radius(const Geometry &geom)
{
  double min = 0;
  double max = 0;
  Vec3d center;
  return edge_nearpoints_radius(geom, min, max, center);
}

// sets radius of geom to average of edge near points radius
void unitize_nearpoints_radius(Geometry &geom)
{
  double avg = edge_nearpoints_radius(geom);
  geom.transform(Trans3d::scale(1 / avg));
}

// return true if maximum vertex radius is radius_range_percent (0.0 to ...)
// greater than minimum vertex radius (visible for canonical.cc)
bool canonical_radius_range_test(const Geometry &geom,
                                 const double radius_range_percent)
{
  GeometryInfo rep(geom);
  rep.set_center(geom.centroid());

  double min = rep.vert_dist_lims().min;
  double max = rep.vert_dist_lims().max;

  // min and max should always be positive, max should always be larger
  return (((max - min) / ((max + min) / 2.0)) > radius_range_percent) ? true
                                                                      : false;
}

// Addition to algorithm by Adrian Rossiter
// Finds the edge near points centroid
Vec3d edge_nearpoints_centroid(Geometry &geom, const Vec3d cent)
{
  vector<vector<int>> edges;
  geom.get_impl_edges(edges);
  Vec3d e_cent(0, 0, 0);
  for (auto &edge : edges)
    e_cent += geom.edge_nearpt(edge, cent);
  return e_cent / double(edges.size());
}

// reciprocalN() is from the Hart's Conway Notation web page
// make array of vertices reciprocal to given planes (face normals)
// RK - save of verbatim port code
/*
vector<Vec3d> reciprocalN_old(const Geometry &geom)
{
  const vector<vector<int>> &faces = geom.faces();
  const vector<Vec3d> &verts = geom.verts();

  vector<Vec3d> normals;
  for (const auto &face : faces) {
    Vec3d centroid(0, 0, 0);
    Vec3d normal(0, 0, 0);
    double avgEdgeDist = 0;

    int v1 = face.at(face.size() - 2);
    int v2 = face.at(face.size() - 1);
    for (int v3 : face) {
      centroid += verts[v3];
      // orthogonal() was from the Hart's Conway Notation web page. replacement
      // normal += orthogonal(verts[v1], verts[v2], verts[v3]);
      normal += vcross(verts[v3] - verts[v2], verts[v2] - verts[v1]);
      // tangentPoint() was from Hart's Conway Notation web page. replacement
      // avgEdgeDist += tangentPoint(verts[v1], verts[v2]).len();
      Vec3d d = verts[v2] - verts[v1];
      // prevent division by zero
      // avgEdgeDist += (verts[v1] - ((vdot(d,verts[v1])/d.len2()) * d)).len();
      double vdt;
      if (d[0] == 0 && d[1] == 0 && d[2] == 0)
        vdt = 0;
      else
        vdt = vdot(d, verts[v1]) / d.len2();
      avgEdgeDist += (verts[v1] - (vdt * d)).len(); // tangentPoint without call
      v1 = v2;
      v2 = v3;
    }
    centroid *= 1.0 / face.size();
    normal.to_unit();
    avgEdgeDist /= face.size();

    // reciprocal call replace below:
    // prevent division by zero
    // Vec3d ans = reciprocal(normal * vdot(centroid,normal));
    Vec3d v = normal * vdot(centroid, normal);
    Vec3d ans;
    if (v[0] == 0 && v[1] == 0 && v[2] == 0)
      ans = v;
    else {
      ans = v * 1.0 / v.len2();
      ans *= (1 + avgEdgeDist) / 2;
    }
    normals.push_back(ans);
  }

  return normals;
}
*/

/* RK - save of simplified tangent code
      Vec3d d = verts[v2] - verts[v1];
      double vdt = 0;
      // prevent division by zero
      if (d[0] != 0 || d[1] != 0 || d[2] != 0)
        vdt = vdot(d, verts[v1]) / d.len2();
      avgEdgeDist += (verts[v1] - (vdt * d)).len();
*/

/*
// RK - Gives the same answer as the built in face_norm().unit()
// even when nonplanar and measuring all edges
Vec3d face_norm_newell(const Geometry &geom, vector<int> &face)
{
  const vector<Vec3d> &v = geom.verts();

  Vec3d face_normal(0, 0, 0);

  unsigned int sz = face.size();
  for (unsigned int i = 0; i < sz; i++) {
    int v1 = face[i];
    int v2 = face[(i + 1) % sz];

    face_normal[0] += (v[v1][1] - v[v2][1]) * (v[v1][2] + v[v2][2]);
    face_normal[1] += (v[v1][2] - v[v2][2]) * (v[v1][0] + v[v2][0]);
    face_normal[2] += (v[v1][0] - v[v2][0]) * (v[v1][1] + v[v2][1]);
  }

  return face_normal.to_unit();
}

// return the unit normal of all perimeter triangles
Vec3d face_norm_newell(const Geometry &geom, const int f_idx)
{
  vector<int> face = geom.faces(f_idx);
  return face_norm_newell(geom, face);
}
*/

// reciprocalN() is from the Hart's Conway Notation web page
// make array of vertices reciprocal to given planes (face normals)
// RK - has accuracy issues and will have trouble with -l 16
vector<Vec3d> reciprocalN(const Geometry &geom)
{
  vector<Vec3d> normals;

  for (const auto &face : geom.faces()) {
    // RK - the algoritm was written to use triangles for measuring
    // non-planar faces. Now method can be chosen
    Vec3d face_normal = face_norm(geom.verts(), face).unit();
    Vec3d face_centroid = anti::centroid(geom.verts(), face);
    // make sure face_normal points outward
    if (vdot(face_normal, face_centroid) < 0)
      face_normal *= -1.0;

    // RK - find the average lenth of the edge near points
    unsigned int sz = face.size();
    double avgEdgeDist = 0;
    for (unsigned int j = 0; j < sz; j++) {
      int v1 = face[j];
      int v2 = face[(j + 1) % sz];

      avgEdgeDist += geom.edge_nearpt(make_edge(v1, v2), Vec3d(0, 0, 0)).len2();
    }

    // RK - sqrt of length squared here
    avgEdgeDist = sqrt(avgEdgeDist / sz);

    // the face normal height set to intersect face at v
    Vec3d v = face_normal * vdot(face_centroid, face_normal);

    // adjust v to the reciprocal value
    Vec3d ans = v;
    // prevent division by zero
    if (v[0] != 0 || v[1] != 0 || v[2] != 0)
      ans = v * 1.0 / v.len2();

    // edge correction (of v based on all edges of the face)
    ans *= (1 + avgEdgeDist) / 2;

    normals.push_back(ans);
  }

  return normals;
}

// reciprocate on face centers dividing by magnitude squared
vector<Vec3d> reciprocalC_len2(const Geometry &geom)
{
  vector<Vec3d> centers;
  geom.face_cents(centers);
  for (auto &center : centers)
    center /= center.len2();
  return centers;
}

// Implementation of George Hart's planarization and canonicalization algorithms
// http://www.georgehart.com/virtual-polyhedra/conway_notation.html
bool canonicalize_bd(Geometry &base, IterationControl it_ctrl,
                     double radius_range_percent, const bool planarize_only)
{
  bool completed = false;
  it_ctrl.set_finished(false);

  Geometry dual;
  // the dual's initial vertex locations are immediately overwritten
  get_dual(dual, base, 1);
  dual.clear_cols();

  double test_val = it_ctrl.get_test_val();
  double max_diff2 = 0;

  for (it_ctrl.start_iter(); !it_ctrl.is_done(); it_ctrl.next_iter()) {
    vector<Vec3d> base_verts_last = base.verts();

    if (planarize_only) {
      // adjust vertices with side effect of planarization. len2() version
      dual.raw_verts() = reciprocalC_len2(base);
      base.raw_verts() = reciprocalC_len2(dual);
      // move centroid to origin for balance
      base.transform(Trans3d::translate(-centroid(base.verts())));
    }
    else {
      // base/dual canonicalize method
      dual.raw_verts() = reciprocalN(base);
      base.raw_verts() = reciprocalN(dual);
      // re-center for drift
      base.transform(
          Trans3d::translate(-edge_nearpoints_centroid(base, Vec3d(0, 0, 0))));
    }

    // reduces size imbalance problem with this algorithm
    unitize_nearpoints_radius(base);

    string finish_msg;
    if (it_ctrl.is_status_check_iter()) {
      // len2() for difference value to minimize internal sqrt() calls
      max_diff2 = 0;
      for (unsigned int i = 0; i < base.verts().size(); i++) {
        double diff2 = (base.verts(i) - base_verts_last[i]).len2();
        if (diff2 > max_diff2)
          max_diff2 = diff2;
      }

      // break out if volume goes to zero, restore last good vertices
      if (std::isnan(base.verts(0)[0])) {
        it_ctrl.set_finished();
        finish_msg = "error, volume went to zero";
        base.raw_verts() = base_verts_last;
        radius_range_percent = 0; // bypass range test
      }
      else if (sqrt(max_diff2) < test_val) {
        completed = true;
        it_ctrl.set_finished();
        finish_msg = "solved, test value achieved";
      }
      else if (it_ctrl.is_last_iter()) {
        // reached last iteration without solving
        it_ctrl.set_finished();
        finish_msg = "not solved, test value not achieved";
      }

      // check if radius is expanding or contracting unreasonably,
      // but only for the purpose of finishing early
      // if minimum and maximum radius are differing, the polyhedron is
      // crumpling
      if (radius_range_percent &&
          canonical_radius_range_test(base, radius_range_percent)) {
        if (!it_ctrl.is_finished())
          it_ctrl.set_finished();
        finish_msg = "breaking out: radius range detected. try increasing -d";
      }
    }

    if (it_ctrl.is_status_report_iter()) {
      if (it_ctrl.is_finished())
        it_ctrl.print("Final iteration (%s):\n", finish_msg.c_str());

      it_ctrl.print("%-12u max_diff:%17.15e\n", it_ctrl.get_current_iter(),
                    sqrt(max_diff2));
    }
  }

  return completed;
}

// RK - wrapper for basic planarization with base/dual algorithm
// meant to be called with finite num_iters (not -1)
bool planarize_bd(Geometry &geom, IterationControl it_ctrl)
{
  double radius_range_percent = 0;
  bool planarize_only = true;
  return canonicalize_bd(geom, it_ctrl, radius_range_percent, planarize_only);
}

/*
// plane a single face aligned to the z axis
void plane_face(Geometry &polygon)
{
  Vec3d face_normal = polygon.face_norm(0).unit();
  Vec3d face_centroid = polygon.face_cent(0);
  // make sure face_normal points outward
  if (vdot(face_normal, face_centroid) < 0)
    face_normal *= -1.0;

  // this gives the same results (from the mathematica algorithm)
  // for (auto &vert : polygon.raw_verts())
  //  vert += vdot(face_normal, face_centroid - vert) * face_normal;
  // return;

  // rotate face to z axis
  Trans3d trans = Trans3d::rotate(face_normal, Vec3d(0, 0, 1));
  polygon.transform(trans);

  // refresh face centroid
  face_centroid = polygon.face_cent(0);

  // set z of all vertices to height of face centroid
  for (auto &vert : polygon.raw_verts())
    vert[2] = face_centroid[2];

  // rotate face back to original position
  polygon.transform(trans.inverse());
}
*/

// color edges by dihedral angle
void color_edges_by_dihedral(Geometry &geom, const ColorMapMulti &edge_map,
                             double eps)
{
  GeometryInfo info(geom);

  // if negative volume, orientation reversed, dihedral angles opposite
  bool reverse = (GeometryInfo(geom).volume() < 0) ? true : false;

  // color all elements convex
  geom.add_missing_impl_edges();
  Coloring(&geom).e_one_col(edge_map.get_col(0));

  vector<double> dihedrals = info.get_edge_dihedrals();
  // assuming edge order is intact
  for (unsigned int i = 0; i < dihedrals.size(); i++) {
    bool convexity = true;       // assume true;
    bool coplanar_found = false; // assume false;
    double angle = dihedrals[i];
    if (double_eq(angle, M_PI, eps))
      coplanar_found = true;
    else if ((!reverse && (angle > M_PI)) || (reverse && (angle < M_PI)))
      convexity = false;

    if (convexity && !coplanar_found)
      continue;

    // if not coplanar then it was nonconvex
    Color ec = (coplanar_found) ? edge_map.get_col(1) : edge_map.get_col(2);

    geom.colors(EDGES).set(i, ec);
  }
}

// color faces by convexity compared to other faces
void color_faces_by_convexity(Geometry &geom, const ColorMapMulti &face_map,
                              double eps)
{
  GeometryInfo info(geom);

  // if negative volume, orientation reversed, dihedral angles opposite
  bool reverse = (GeometryInfo(geom).volume() < 0) ? true : false;

  // color all elements convex
  geom.add_missing_impl_edges();
  Coloring(&geom).f_one_col(face_map.get_col(0));

  vector<double> dihedrals = info.get_edge_dihedrals();
  // assuming edge order is intact
  for (unsigned int i = 0; i < geom.faces().size(); i++) {
    bool convexity = true;       // assume true;
    bool coplanar_found = false; // assume false;
    vector<int> face = geom.faces(i);
    unsigned int fsz = face.size();
    for (unsigned int j = 0; j < fsz; j++) {
      int v1 = face[j];
      int v2 = face[(j + 1) % fsz];
      vector<int> edge = make_edge(v1, v2);
      double angle = dihedrals[find_edge_in_edge_list(geom.edges(), edge)];
      if (double_eq(angle, M_PI, eps))
        coplanar_found = true;
      else if ((!reverse && (angle > M_PI)) || (reverse && (angle < M_PI)))
        convexity = false;

      // give priority to coplanar
      if (coplanar_found)
        break;
    }

    // if's fail then its convex
    if (coplanar_found)
      geom.colors(FACES).set(i, face_map.get_col(1));
    else if (!convexity)
      geom.colors(FACES).set(i, face_map.get_col(2));
  }
}
