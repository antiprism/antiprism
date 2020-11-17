/*
   Copyright (c) 2003-2016, Adrian Rossiter, Roger Kaufman
   Includes ideas and algorithms by George W. Hart, http://www.georgehart.com

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
   Name: canonic.cc
   Description: canonicalize a polyhedron
                Implementation of George Hart's canonicalization algorithm
                http://library.wolfram.com/infocenter/Articles/2012/
   Project: Antiprism - http://www.antiprism.com
*/

#include "boundbox.h"
#include "geometry.h"
#include "geometryinfo.h"
#include "planar.h"

#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <string>
#include <vector>

using std::map;
using std::set;
using std::string;
using std::vector;

namespace anti {

// RK - find nearpoints radius, sets range minimum and maximum
double edge_nearpoints_radius(const Geometry &geom, double &min, double &max,
                              Vec3d &center)
{
  min = DBL_MAX;
  max = DBL_MIN;

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

// Implementation of George Hart's canonicalization algorithm
// http://library.wolfram.com/infocenter/Articles/2012/
// RK - the model will possibly become non-convex early in the loops.
// if it contorts too badly, the model will implode. Having the input
// model at a radius of near 1 minimizes this problem
bool canonicalize_mm(Geometry &geom, IterationControl it_ctrl,
                     const double edge_factor, const double plane_factor,
                     const double radius_range_percent, const char normal_type,
                     const bool alternate_loop, const bool planarize_only)
{
  bool completed = false;
  it_ctrl.set_finished(false);

  vector<Vec3d> &verts = geom.raw_verts();

  vector<vector<int>> edges;
  geom.get_impl_edges(edges);

  double test_val = it_ctrl.get_test_val();
  double max_diff2 = 0;

  for (it_ctrl.start_iter(); !it_ctrl.is_done(); it_ctrl.next_iter()) {
    vector<Vec3d> verts_last = verts;

    if (!planarize_only) {
      vector<Vec3d> near_pts;
      if (!alternate_loop) {
        for (auto &edge : edges) {
          Vec3d P = geom.edge_nearpt(edge, Vec3d(0, 0, 0));
          near_pts.push_back(P);
          Vec3d offset = edge_factor * (P.len() - 1) * P;
          verts[edge[0]] -= offset;
          verts[edge[1]] -= offset;
        }
      }
      // RK - alternate form causes the near points to be applied in a 2nd loop
      // most often not needed unless the model is off balance
      else {
        for (auto &edge : edges) {
          Vec3d P = geom.edge_nearpt(edge, Vec3d(0, 0, 0));
          near_pts.push_back(P);
          // RK - these 4 lines cause the near points to be applied in a 2nd
          // loop
        }
        int p_cnt = 0;
        for (auto &edge : edges) {
          Vec3d P = near_pts[p_cnt++];
          Vec3d offset = edge_factor * (P.len() - 1) * P;
          verts[edge[0]] -= offset;
          verts[edge[1]] -= offset;
        }
      }
      /*
            // RK - revolving loop. didn't solve the imbalance problem
            else {
              for (unsigned int ee = cnt; ee < edges.size() + cnt; ee++) {
                unsigned int e = ee % edges.size();
                Vec3d P = geom.edge_nearpt(edges[e], Vec3d(0, 0, 0));
                near_pts.push_back(P);
                Vec3d offset = edge_factor * (P.len() - 1) * P;
                verts[edges[e][0]] -= offset;
                verts[edges[e][1]] -= offset;
              }
            }
      */

      // re-center for drift
      Vec3d cent_near_pts = centroid(near_pts);
      for (unsigned int i = 0; i < verts.size(); i++)
        verts[i] -= cent_near_pts;
    }

    // Make a copy of verts into vs and zero out
    // Accumulate vertex changes instead of altering vertices in place
    // This can help relieve when a vertex is pushed towards one plane
    // and away from another
    vector<Vec3d> vs = verts;
    for (auto &v : vs)
      v = Vec3d(0, 0, 0);

    for (unsigned int f = 0; f < geom.faces().size(); f++) {
      if (geom.faces(f).size() == 3)
        continue;
      Vec3d face_normal = face_normal_by_type(geom, f, normal_type).unit();
      Vec3d face_centroid = geom.face_cent(f);
      // make sure face_normal points outward
      if (vdot(face_normal, face_centroid) < 0)
        face_normal *= -1.0;
      // place a planar vertex over or under verts[v]
      // adds or subtracts it to get to the planar verts[v]
      for (int v : geom.faces(f))
        vs[v] += vdot(plane_factor * face_normal, face_centroid - verts[v]) *
                 face_normal;
    }

    // adjust vertices post-loop
    for (unsigned int i = 0; i < vs.size(); i++)
      verts[i] += vs[i];

    string finish_msg;
    if (it_ctrl.is_status_check_iter()) {
      // len2() for difference value to minimize internal sqrt() calls
      max_diff2 = 0;
      for (unsigned int i = 0; i < verts.size(); i++) {
        double diff2 = (verts[i] - verts_last[i]).len2();
        if (diff2 > max_diff2)
          max_diff2 = diff2;
      }

      if (sqrt(max_diff2) < test_val) {
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
          canonical_radius_range_test(geom, radius_range_percent)) {
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

// return the normal of all perimeter triangles
Vec3d face_norm_nonplanar_triangles(const Geometry &geom,
                                    const vector<int> &face)
{
  Vec3d face_normal(0, 0, 0);

  unsigned int sz = face.size();
  for (unsigned int i = 0; i < sz; i++) {
    int v0 = face[i];
    int v1 = face[(i + 1) % sz];
    int v2 = face[(i + 2) % sz];

    face_normal += vcross(geom.verts()[v0] - geom.verts()[v1],
                          geom.verts()[v1] - geom.verts()[v2]);
  }

  return face_normal;
}

// return the normal of all perimeter triangles
Vec3d face_norm_nonplanar_triangles(const Geometry &geom, const int f_idx)
{
  vector<int> face = geom.faces(f_idx);
  return face_norm_nonplanar_triangles(geom, face);
}

// return the normal of all quads in polygon
Vec3d face_norm_nonplanar_quads(const Geometry &geom, const vector<int> &face)
{
  Vec3d face_normal(0, 0, 0);

  unsigned int sz = face.size();
  for (unsigned int i = 0; i < sz; i++) {
    int v0 = face[i];
    int v1 = face[(i + 1) % sz];
    int v2 = face[(i + 2) % sz];
    int v3 = face[(i + 3) % sz];

    face_normal += vcross(geom.verts()[v0] - geom.verts()[v2],
                          geom.verts()[v1] - geom.verts()[v3]);
  }

  return face_normal;
}

// return the normal of quads in polygon
Vec3d face_norm_nonplanar_quads(const Geometry &geom, const int f_idx)
{
  vector<int> face = geom.faces(f_idx);
  return face_norm_nonplanar_quads(geom, face);
}

// select normal by type. Newell, triangles, or quads
// 'n' is the default, if normal_type is not given or is wrong
Vec3d face_normal_by_type(const Geometry &geom, const vector<int> &face,
                          const char normal_type)
{
  Vec3d face_normal;

  if (normal_type == 't')
    face_normal = face_norm_nonplanar_triangles(geom, face);
  else if (normal_type == 'q')
    face_normal = face_norm_nonplanar_quads(geom, face);
  else // if (normal_type == 'n')
    face_normal = geom.face_norm(face);
  return face_normal;
}

// select normal by type. Newell, triangles, or quads
Vec3d face_normal_by_type(const Geometry &geom, const int f_idx,
                          const char normal_type)
{
  vector<int> face = geom.faces(f_idx);
  return face_normal_by_type(geom, face, normal_type);
}

// reciprocalN() is from the Hart's Conway Notation web page
// make array of vertices reciprocal to given planes (face normals)
// RK - has accuracy issues and will have trouble with -l 16
vector<Vec3d> reciprocalN(const Geometry &geom, const char normal_type)
{
  vector<Vec3d> normals;

  for (const auto &face : geom.faces()) {
    // RK - the algoritm was written to use triangles for measuring
    // non-planar faces. Now method can be chosen
    Vec3d face_normal = face_normal_by_type(geom, face, normal_type).unit();
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

// Implementation of George Hart's planarization and canonicalization algorithms
// http://www.georgehart.com/virtual-polyhedra/conway_notation.html
bool canonicalize_bd(Geometry &base, IterationControl it_ctrl,
                     const char canonical_method,
                     const double radius_range_percent, const char centering,
                     const char normal_type)
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

    switch (canonical_method) {
    // base/dual canonicalize method
    case 'b': {
      dual.raw_verts() = reciprocalN(base, normal_type);
      base.raw_verts() = reciprocalN(dual, normal_type);
      // re-center for drift
      if (centering == 'e')
        base.transform(Trans3d::translate(
            -edge_nearpoints_centroid(base, Vec3d(0, 0, 0))));
      else if (centering == 'v')
        base.transform(Trans3d::translate(-centroid(base.verts())));
      break;
    }

    // adjust vertices with side effect of planarization. len2() version
    case 'q':
      // move centroid to origin for balance
      dual.raw_verts() = reciprocalC_len2(base);
      base.raw_verts() = reciprocalC_len2(dual);
      base.transform(Trans3d::translate(-centroid(base.verts())));
      break;
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

      if (sqrt(max_diff2) < test_val) {
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
// for an internal call from conway
bool planarize_bd(Geometry &geom, IterationControl it_ctrl)
{
  char canonical_method = 'q';
  double radius_range_percent = 0;
  char centering = 'x';
  char normal_type = 'n';
  return canonicalize_bd(geom, it_ctrl, canonical_method, radius_range_percent,
                         centering, normal_type);
}

//------------------------------------------------------------------
// Unit edge regular polygon

Status make_regular_faces2(Geometry &base_geom, IterationControl it_ctrl,
                           double shorten_factor, double plane_factor,
                           double radius_factor, const string &sym_str)
{
  Symmetry sym(Symmetry::C1);
  Status stat;
  if (sym_str.size()) {
    // if (!(stat = init_sym(base_geom, sym_str.c_str(), sym)))
    //  return stat;
  }
  bool using_symmetry = (sym.get_sym_type() != Symmetry::C1);

  // No further processing if no faces, but not an error
  if (base_geom.faces().size() == 0)
    return Status::ok();

  double test_val = it_ctrl.get_test_val();
  const double divergence_test2 = 1e30; // test vertex dist^2 for divergence
  // Scale to get edges close to 1
  GeometryInfo info(base_geom);
  double scale = info.iedge_length_lims().sum / info.num_iedges();
  if (scale)
    base_geom.transform(Trans3d::scale(1 / scale));

  SymmetricUpdater sym_updater((using_symmetry) ? base_geom : Geometry(), sym);
  const Geometry &geom =
      (using_symmetry) ? sym_updater.get_geom_working() : base_geom;

  const vector<Vec3d> &verts = geom.verts();
  const vector<vector<int>> &faces = geom.faces();

  double max_diff2 = 0;
  Vec3d origin(0, 0, 0);
  vector<double> rads(faces.size());
  for (unsigned int f = 0; f < faces.size(); f++) {
    int N = faces[f].size();
    int D = std::abs(find_polygon_denominator_signed(geom, f, epsilon));
    if (!D)
      D = 1;
    rads[f] = 0.5 / sin(M_PI * D / N); // circumradius of regular polygon
  }

  // Get a list of the faces that contain a principal vertex of any type
  vector<int> principal_verts;
  vector<int> verts_to_update;
  vector<int> faces_to_process;
  vector<int> edges_to_process;
  if (using_symmetry) {
    vector<set<int>> vert_faces(geom.verts().size());
    for (int f_idx = 0; f_idx < (int)geom.faces().size(); f_idx++)
      for (int v_idx : faces[f_idx])
        vert_faces[v_idx].insert(f_idx);

    vector<set<int>> vert_edges(info.get_impl_edges().size());
    for (int e_idx = 0; e_idx < (int)info.get_impl_edges().size(); e_idx++)
      for (int v_idx : info.get_impl_edges()[e_idx])
        vert_edges[v_idx].insert(e_idx);

    for (auto &v_orbit : sym_updater.get_equiv_sets(VERTS)) {
      int principal_idx = *v_orbit.begin();
      principal_verts.push_back(principal_idx);
      for (int f_idx : vert_faces[principal_idx]) {
        faces_to_process.push_back(f_idx);
        for (int v_idx : faces[f_idx])
          verts_to_update.push_back(v_idx);
      }
      for (int e_idx : vert_edges[principal_idx])
        edges_to_process.push_back(e_idx);
      // implicit edges, so vertices already addded to verts_to_update
    }
    sort(verts_to_update.begin(), verts_to_update.end());
    verts_to_update.erase(
        unique(verts_to_update.begin(), verts_to_update.end()),
        verts_to_update.end());
    sort(faces_to_process.begin(), faces_to_process.end());
    faces_to_process.erase(
        unique(faces_to_process.begin(), faces_to_process.end()),
        faces_to_process.end());
  }
  else { // not using_symmetry
    // all vertices are proncipal vertices
    principal_verts.resize(verts.size());
    std::iota(principal_verts.begin(), principal_verts.end(), 0);
    // all face numbers should be be processed
    faces_to_process.resize(faces.size());
    std::iota(faces_to_process.begin(), faces_to_process.end(), 0);
    // all edges should be be processed
    edges_to_process.resize(info.get_impl_edges().size());
    std::iota(edges_to_process.begin(), edges_to_process.end(), 0);
  }

  vector<Vec3d> offsets(verts.size()); // Vertex adjustments

  for (it_ctrl.start_iter(); !it_ctrl.is_done(); it_ctrl.next_iter()) {
    std::fill(offsets.begin(), offsets.end(), Vec3d::zero);

    if (using_symmetry) {
      // Ensure that the vertices used in adjustment are up to date
      for (int v_idx : verts_to_update)
        sym_updater.update_from_principal_vertex(v_idx);
    }

    for (int f_idx : faces_to_process) {
      const vector<int> &face = faces[f_idx];

      Vec3d norm = geom.face_norm(f_idx).unit();
      Vec3d f_cent = geom.face_cent(f_idx);
      if (vdot(norm, f_cent) < 0)
        norm *= -1.0;

      for (int i = 0; i < (int)face.size(); i++) {
        const int v_idx = face[i];
        // offset for planarity
        offsets[v_idx] +=
            vdot(plane_factor * norm, f_cent - verts[v_idx]) * norm;

        // offset for polygon radius
        Vec3d rad_vec = (verts[v_idx] - f_cent);
        offsets[v_idx] +=
            (rads[f_idx] - rad_vec.len()) * radius_factor * rad_vec;
      }
    }

    // offsets for unit edges
    for (int e_idx : edges_to_process) {
      const auto &edge = info.get_impl_edges()[e_idx];
      Vec3d offset =
          (1 - geom.edge_len(edge)) * shorten_factor * geom.edge_vec(edge);
      offsets[edge[0]] -= 2 * offset;
      offsets[edge[1]] += 2 * offset;
    }

    // adjust vertices post-loop
    if (using_symmetry) {
      // adjust principal vertices
      for (int v_idx : principal_verts)
        sym_updater.update_principal_vertex(v_idx,
                                            verts[v_idx] + offsets[v_idx]);
    }
    else { // not using_symmetry
      // adjust all vertices
      for (unsigned int i = 0; i < verts.size(); i++)
        base_geom.raw_verts()[i] += offsets[i];
    }

    string finish_msg;
    if (it_ctrl.is_status_check_iter()) {

      max_diff2 = 0;
      for (auto &offset : offsets) {
        double diff2 = offset.len2();
        if (diff2 > max_diff2)
          max_diff2 = diff2;
      }

      double width = BoundBox(verts).max_width();
      if (sqrt(max_diff2) / width < test_val) {
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
      if (!it_ctrl.is_finished() && divergence_test2 > 0) {
        for (auto &offset : offsets)
          if (offset.len2() > divergence_test2) {
            it_ctrl.set_finished();
            finish_msg = "not solved, quit early as probably diverging";
          }
      }
    }

    if (it_ctrl.is_status_report_iter()) {
      if (it_ctrl.is_finished())
        it_ctrl.print("Final iteration (%s):\n", finish_msg.c_str());

      it_ctrl.print("%-12u max_diff:%17.15e\n", it_ctrl.get_current_iter(),
                    sqrt(max_diff2));
    }
  }

  if (using_symmetry)
    base_geom = sym_updater.get_geom_final();

  return Status::ok();
}

//------------------------------------------------------------------
// Make faces planar

inline Vec3d nearpoint_on_plane(const Vec3d &P, const Vec3d &point_on_plane,
                                const Vec3d &unit_norm)
{
  return P + vdot(point_on_plane - P, unit_norm) * unit_norm;
}

Status make_planar2(Geometry &base_geom, IterationControl it_ctrl,
                    double plane_factor)
{
  // chosen by experiment
  const double intersect_test_val = 1e-5; // test for coplanar faces
  const double diff2_test_val = 10;       // limit for using plane intersection
  const double readjustment = 1.01;       // to adjust adjustment factor
  const double plane_factor_max = 1.1;    // maximum value for adjustment factor

  auto &geom = base_geom;
  // No further processing if no faces, but not an error
  if (geom.faces().size() == 0)
    return Status::ok();

  double test_val = it_ctrl.get_test_val();
  const vector<Vec3d> &verts = geom.verts();
  const vector<vector<int>> &faces = geom.faces();

  vector<vector<int>> vert_faces(geom.verts().size());
  for (int f_idx = 0; f_idx < (int)geom.faces().size(); f_idx++)
    for (int v_idx : faces[f_idx])
      vert_faces[v_idx].push_back(f_idx);

  vector<Vec3d> offsets(verts.size()); // Vertex adjustments

  double last_max_diff2 = 0.0;
  for (it_ctrl.start_iter(); !it_ctrl.is_done(); it_ctrl.next_iter()) {
    std::fill(offsets.begin(), offsets.end(), Vec3d::zero);

    GeometryInfo info(geom);
    vector<Vec3d> norms;
    geom.face_norms(norms);
    vector<Vec3d> cents;
    geom.face_cents(cents);

    int cnt_proj = 0;
    int cnt_int = 0;
    double max_diff2 = 0.0;
    for (unsigned int v_idx = 0; v_idx < geom.verts().size(); v_idx++) {
      int intersect_cnt = 0;
      const auto &vfaces = vert_faces[v_idx];
      const int vf_sz = vfaces.size();
      bool good_intersections = (vf_sz >= 3);
      for (int f0 = 0; f0 < vf_sz - 2 && good_intersections; f0++) {
        int f0_idx = vfaces[f0];
        for (int f1 = f0 + 1; f1 < vf_sz - 1 && good_intersections; f1++) {
          int f1_idx = vfaces[f1];
          for (int f2 = f1 + 1; f2 < vf_sz && good_intersections; f2++) {
            int f2_idx = vfaces[f2];
            Vec3d intersection;
            good_intersections = three_plane_intersect(
                cents[f0_idx], norms[f0_idx], cents[f1_idx], norms[f1_idx],
                cents[f2_idx], norms[f2_idx], intersection, intersect_test_val);
            offsets[v_idx] += intersection;
            intersect_cnt++;
          }
        }
      }

      if (good_intersections)
        offsets[v_idx] =
            (offsets[v_idx] / intersect_cnt - verts[v_idx]) * plane_factor;

      // no good 3 plane intersections OR
      // moving to much
      if (!good_intersections ||
          offsets[v_idx].len2() / last_max_diff2 > diff2_test_val) {
        // target vertex is centroid of projection of vertex onto planes
        offsets[v_idx] = Vec3d::zero;
        for (int f0 = 0; f0 < vf_sz; f0++) {
          int f0_idx = vfaces[f0];
          offsets[v_idx] +=
              nearpoint_on_plane(verts[v_idx], cents[f0_idx], norms[f0_idx]);
        }
        offsets[v_idx] = (offsets[v_idx] / vf_sz - verts[v_idx]) * plane_factor;
        cnt_proj++;
      }
      else
        cnt_int++;

      auto diff2 = offsets[v_idx].len2();
      if (diff2 > max_diff2)
        max_diff2 = diff2;
    }

    for (unsigned int v_idx = 0; v_idx < verts.size(); v_idx++)
      geom.verts(v_idx) += offsets[v_idx];

    // adjust plane factor
    if (max_diff2 < last_max_diff2)
      plane_factor *= readjustment;
    else
      plane_factor /= readjustment;
    if (plane_factor > plane_factor_max)
      plane_factor = plane_factor_max;
    last_max_diff2 = max_diff2;

    string finish_msg;
    if (it_ctrl.is_status_check_iter()) {
      double width = BoundBox(verts).max_width();
      if (sqrt(max_diff2) / width < test_val) {
        it_ctrl.set_finished();
        finish_msg = "solved, test value achieved";
      }
      else if (it_ctrl.is_last_iter()) {
        // reached last iteration without solving
        it_ctrl.set_finished();
        finish_msg = "not solved, test value not achieved";
      }
    }

    if (it_ctrl.is_status_report_iter()) {
      if (it_ctrl.is_finished())
        it_ctrl.print("Final iteration (%s):\n", finish_msg.c_str());

      it_ctrl.print("%-12u max_diff:%17.15e  -f %-10.5f i/p: %7d/%-7d\n",
                    it_ctrl.get_current_iter(), sqrt(max_diff2),
                    100 * plane_factor, cnt_int, cnt_proj);
    }
  }

  return Status::ok();
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

// P and Q are modified
void move_line_to_point(Vec3d &P, Vec3d &Q, const Vec3d &X)
{
  Vec3d Y = X + (Q - P);
  Vec3d V = P + X;
  Vec3d P2 = lines_intersection(P, V, X, Y, 0);
  if (P2.is_set()) {
    Q += P2 - P;
    P = P2;
  }
}

// RK - edge near points of base seek 1
bool canonicalize_unit(Geometry &geom, IterationControl it_ctrl,
                       const double radius_range_percent, const char centering,
                       const char normal_type, const bool planarize_only)
{
  bool completed = false;
  it_ctrl.set_finished(false);

  vector<vector<int>> edges;
  geom.get_impl_edges(edges);

  vector<Vec3d> &verts = geom.raw_verts();

  double test_val = it_ctrl.get_test_val();
  double max_diff2 = 0;

  for (it_ctrl.start_iter(); !it_ctrl.is_done(); it_ctrl.next_iter()) {
    vector<Vec3d> verts_last = verts;

    if (!planarize_only) {
      for (auto &edge : edges) {
        // unit near point
        Vec3d P = geom.edge_nearpt(edge, Vec3d(0, 0, 0)).unit();
        move_line_to_point(verts[edge[0]], verts[edge[1]], P);
      }

      // re-center for drift
      if (centering == 'e')
        geom.transform(Trans3d::translate(
            -edge_nearpoints_centroid(geom, Vec3d(0, 0, 0))));
      else if (centering == 'v')
        geom.transform(Trans3d::translate(-centroid(geom.verts())));
    }

    for (unsigned int f = 0; f < geom.faces().size(); f++) {
      /*
            // give polygon its own geom. face index needs to reside in vector
            vector<int> face_idxs(1);
            face_idxs[0] = f;
            Geometry polygon = faces_to_geom(geom, face_idxs);
            plane_face(polygon);

            // map vertices back into original geom
            // the numerical order of vertex list in polygon geom is preserved
            vector<int> v_idx;
            for (int v : geom.faces(f))
              v_idx.push_back(v);
            sort(v_idx.begin(), v_idx.end());
            int j = 0;
            for (int v : v_idx)
              verts[v] = polygon.verts(j++);
      */
      // RK - this does formulaically what the above does by brute force
      Vec3d face_normal = face_normal_by_type(geom, f, normal_type).unit();
      Vec3d face_centroid = geom.face_cent(f);
      // make sure face_normal points outward
      if (vdot(face_normal, face_centroid) < 0)
        face_normal *= -1.0;
      // place a planar vertex over or under verts[v]
      // adds or subtracts it to get to the planar verts[v]
      for (int v : geom.faces(f))
        verts[v] += vdot(face_normal, face_centroid - verts[v]) * face_normal;
    }

    string finish_msg;
    if (it_ctrl.is_status_check_iter()) {
      // len2() for difference value to minimize internal sqrt() calls
      max_diff2 = 0;
      for (unsigned int i = 0; i < verts.size(); i++) {
        double diff2 = (verts[i] - verts_last[i]).len2();
        if (diff2 > max_diff2)
          max_diff2 = diff2;
      }

      if (sqrt(max_diff2) < test_val) {
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
          canonical_radius_range_test(geom, radius_range_percent)) {
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

// RK - wrapper for basic planarization with unit algorithm
// meant to be called with finite num_iters (not -1)
// for an internal call from conway
bool planarize_unit(Geometry &geom, IterationControl it_ctrl)
{
  double radius_range_percent = 0;
  char centering = 'x';
  char normal_type = 'n';
  bool planarize_only = true;
  return canonicalize_unit(geom, it_ctrl, radius_range_percent, centering,
                           normal_type, planarize_only);
}

} // namespace anti
