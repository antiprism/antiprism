/*
   Copyright (c) 2003-2022, Adrian Rossiter, Roger Kaufman
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
   Name: canonical.cc
   Description: canonicalize a polyhedron
                Implementation of George Hart's canonicalization algorithm
                http://library.wolfram.com/infocenter/Articles/2012/
   Project: Antiprism - http://www.antiprism.com
*/

#include "boundbox.h"
#include "geometryinfo.h"
#include "utils.h"

using std::string;
using std::vector;

namespace anti {

Geometry base_to_ambo(const Geometry &base, char point_type)
{
  Geometry ambo = base;
  truncate_verts(ambo, 0.5);
  if (point_type == 'n') {
    vector<vector<int>> impl_edges;
    base.get_impl_edges(impl_edges);
    int i = 0;
    for (const auto &edge : impl_edges)
      ambo.verts(i++) = nearest_point(Vec3d(0, 0, 0), base.verts(), edge);
  }
  ambo.transform(Trans3d::translate(-ambo.centroid()));
  return ambo;
}

void update_base_from_ambo(Geometry &base, const Geometry &ambo)
{
  Geometry ambo_d;
  get_pol_recip_verts(ambo_d, ambo, 1, Vec3d::zero);
  auto v_sz = base.verts().size();
  base.raw_verts().clear();
  std::copy(ambo_d.verts().begin(), ambo_d.verts().begin() + v_sz,
            std::back_inserter(base.raw_verts()));
  return;
}

Status make_planar_unit(Geometry &base_geom, IterationControl it_ctrl,
                        double factor, double factor_max, Symmetry sym)
{
  // chosen by experiment
  const double readjust_up = 1.01;    // to adjust adjustment factor up
  const double readjust_down = 0.995; // to adjust adjustment factor down
  const double unit_mult = 1;         // factor for unit adjustment
  const double orth_mult = 0.5;       // extra multiplier for orthogonality adj
  const double overlap_mult = 0.5;    // extra multiplier for overlap adj

  Status stat;

  base_geom.orient(1); // positive orientation

  // transform to centre on centroid, adjust to_std to match
  auto initial_centroid = base_geom.centroid();
  base_geom.transform(Trans3d::translate(base_geom.centroid()));
  sym.set_to_std(sym.get_to_std() * Trans3d::translate(initial_centroid));

  bool using_symmetry = (sym.get_sym_type() > Symmetry::C1);
  auto fixed_subspace = sym.get_fixed_subspace(); // used to find centroid
  SymmetricUpdater sym_updater((using_symmetry) ? base_geom : Geometry(), sym);
  const Geometry &geom =
      (using_symmetry) ? sym_updater.get_geom_working() : base_geom;

  // No further processing if no faces, but not an error
  if (geom.faces().size() == 0)
    return stat;

  const vector<Vec3d> &verts = geom.verts();
  const vector<vector<int>> &faces = geom.faces();

  // List of faces that a vertex is part of
  auto vert_faces = geom.get_info().get_dual().faces();
  vector<vector<int>> vert_figs;
  {
    auto vfigs = geom.get_info().get_vert_figs();
    for (const auto &vfig : vfigs) {
      if (vfig.size() == 0) {
        stat.set_error(
            "intermediate (ambo) model construction error, input polyhedron "
            "may be invalid or degenerate (e.g. contain digons or order-2 "
            "vertices)");
        return stat;
      }
      if (vfig[0].size() != 4) {
        stat.set_error(msg_str(
            "intermediate (ambo) model has vertex with order %d instead of 4",
            (int)vfig[0].size()));
        return stat;
      }
      vert_figs.push_back(vfig[0]);
    }
  }

  // Get a list of the faces that contain a principal vertex of any type
  vector<int> principal_verts;  // first vertex in each vertex orbit
  vector<int> verts_to_update;  // vertices accessed during iteration
  vector<int> faces_to_process; // faces processed during iteration
  if (using_symmetry) {
    principal_verts = sym_updater.get_principal(VERTS);
    faces_to_process = sym_updater.get_associated_elems(vert_faces);
    verts_to_update =
        SymmetricUpdater::get_included_verts(faces_to_process, geom.faces());
    auto vfig_verts_to_update =
        SymmetricUpdater::get_included_verts(principal_verts, vert_figs);
    verts_to_update.insert(verts_to_update.end(), vfig_verts_to_update.begin(),
                           vfig_verts_to_update.end());
    SymmetricUpdater::to_unique_index_list(verts_to_update);
  }
  else { // not using_symmetry
    principal_verts =
        SymmetricUpdater::sequential_index_list(geom.verts().size());
    faces_to_process =
        SymmetricUpdater::sequential_index_list(geom.faces().size());
  }

  // Use oversized arrays to avoid mapping
  vector<Vec3d> offsets(verts.size()); // Vertex adjustments
  vector<Vec3d> norms(faces.size());   // Face normals
  vector<Vec3d> cents(faces.size());   // Face centroids

  double test_val = it_ctrl.get_test_val();
  double last_max_diff2 = 0.0;
  for (it_ctrl.start_iter(); !it_ctrl.is_done(); it_ctrl.next_iter()) {
    std::fill(offsets.begin(), offsets.end(), Vec3d::zero);

    if (using_symmetry) {
      // Ensure that the vertices used in adjustment are up to date
      for (auto v_idx : verts_to_update)
        sym_updater.update_from_principal_vertex(v_idx);
    }

    // Initialize face data for just the necessary faces
    for (auto f_idx : faces_to_process) {
      norms[f_idx] = geom.face_norm(f_idx).unit();
      cents[f_idx] = geom.face_cent(f_idx);
    }

    Vec3d centroid = Vec3d::zero;
    if (using_symmetry) {
      // For each orbit, project a vertex weighted by the orbit size onto
      // the fixed subspace
      const auto &vorbits = sym_updater.get_equiv_sets(VERTS);
      for (const auto &vorbit : vorbits)
        centroid += fixed_subspace.nearest_point(verts[*vorbit.begin()]) *
                    vorbit.size();
      centroid /= verts.size(); // centroid of weighted projected vertices
    }
    else {
      centroid = geom.centroid();
    }

    double max_diff2 = 0.0;
    for (auto v_idx : principal_verts) {
      const auto &vfaces = vert_faces[v_idx];
      const int vf_sz = vfaces.size();
      // target vertex is centroid of projection of vertex onto planes
      for (int f0 = 0; f0 < vf_sz; f0++) {
        int f0_idx = vfaces[f0];
        offsets[v_idx] +=
            nearpoint_on_plane(verts[v_idx], cents[f0_idx], norms[f0_idx]);
      }
      offsets[v_idx] = (offsets[v_idx] / vf_sz - verts[v_idx]) * factor;

      // adjust for centroid
      offsets[v_idx] -= centroid;

      // adjust for orthogonality
      for (int i = 0; i < 2; i++) {
        auto n = vcross(norms[vfaces[i + 2]], norms[vfaces[i]]).unit();
        const auto v_ideal = nearpoint_on_plane(verts[v_idx], Vec3d::zero, n);
        const auto offset = (v_ideal - verts[v_idx]) * factor * orth_mult;
        offsets[v_idx] += offset;
      }

      // adjust for non-overlap
      const auto &vfig = vert_figs[v_idx];
      for (int i = 0; i < 4; i++) {
        if (vtriple(verts[v_idx], verts[vfig[i]], verts[vfig[(i + 1) % 4]]) >
            0) {
          auto v_ideal = anti::centroid(
              {verts[vfig[0]], verts[vfig[1]], verts[vfig[2]], verts[vfig[3]]});
          offsets[v_idx] += (v_ideal - verts[v_idx]) * overlap_mult;
          break;
        }
      }

      auto diff2 = offsets[v_idx].len2();
      if (diff2 > max_diff2)
        max_diff2 = diff2;
    }

    // adjust vertices post-loop
    if (using_symmetry) {
      // adjust principal vertices
      for (int v_idx : principal_verts) {
        auto new_v = verts[v_idx] + offsets[v_idx];
        double new_v_len = new_v.len();
        new_v *= 1 + (1 / new_v_len - 1) * unit_mult;
        sym_updater.update_principal_vertex(v_idx, new_v);
      }
    }
    else { // not using_symmetry
      // adjust all vertices
      for (unsigned int i = 0; i < verts.size(); i++) {
        auto new_v = verts[i] + offsets[i];
        double new_v_len = new_v.len();
        new_v *= 1 + (1 / new_v_len - 1) * unit_mult;
        base_geom.raw_verts()[i] = new_v;
      }
    }

    // adjust plane factor
    if (max_diff2 < last_max_diff2)
      factor *= readjust_up;
    else
      factor *= readjust_down;
    if (factor > factor_max)
      factor = factor_max;
    last_max_diff2 = max_diff2;

    string finish_msg;
    if (it_ctrl.is_status_check_iter()) {
      sym_updater.update_all();
      double width = BoundBox(verts).max_width();
      if (sqrt(max_diff2) / width < test_val) {
        it_ctrl.set_finished();
        finish_msg = "solved, test value achieved";
        stat.set_ok(finish_msg);
      }
      else if (it_ctrl.is_last_iter()) {
        // reached last iteration without solving
        it_ctrl.set_finished();
        finish_msg = "not solved, test value not achieved";
        stat.set_warning(finish_msg);
      }
    }

    if (it_ctrl.is_status_report_iter()) {
      if (it_ctrl.is_finished())
        it_ctrl.print("Final iteration (%s):\n", finish_msg.c_str());

      it_ctrl.print("%-12u max_diff:%17.15e  -f %-10.5f\n",
                    it_ctrl.get_current_iter(), sqrt(max_diff2), 100 * factor);
    }
  }

  if (using_symmetry)
    base_geom = sym_updater.get_geom_final();

  return stat;
}

// was make_canonical_enp in src/canonical.cc
Status make_canonical(Geometry &geom, IterationControl it_ctrl, double factor,
                      double factor_max, char initial_point_type, Symmetry sym)
{
  Geometry ambo = base_to_ambo(geom, initial_point_type);
  Status stat = make_planar_unit(ambo, it_ctrl, factor, factor_max, sym);
  if (!stat.is_error())
    update_base_from_ambo(geom, ambo);
  return stat; // ok - completed, warning - not completed, error - error
}

//------------------------------------------------------------------
// Make faces planar
Status make_planar(Geometry &base_geom, IterationControl it_ctrl,
                   double plane_factor, const Symmetry &sym)
{
  // chosen by experiment
  const double intersect_test_val = 1e-5; // test for coplanar faces
  const double diff2_test_val = 10;       // limit for using plane intersection
  const double readjustment = 1.01;       // to adjust adjustment factor
  const double plane_factor_max = 1.1;    // maximum value for adjustment factor
  bool using_symmetry = (sym.get_sym_type() > Symmetry::C1);

  Status stat;

  SymmetricUpdater sym_updater((using_symmetry) ? base_geom : Geometry(), sym);
  const Geometry &geom =
      (using_symmetry) ? sym_updater.get_geom_working() : base_geom;

  // No further processing if no faces, but not an error
  if (geom.faces().size() == 0)
    return stat;

  const vector<Vec3d> &verts = geom.verts();
  const vector<vector<int>> &faces = geom.faces();

  // List of faces that a vertex is part of
  GeometryInfo info(geom);
  const auto &vert_faces = info.get_vert_faces();

  // Get a list of the faces that contain a principal vertex of any type
  vector<int> principal_verts;  // first vertex in each vertex orbit
  vector<int> verts_to_update;  // vertices accessed during iteration
  vector<int> faces_to_process; // faces processed during iteration
  if (using_symmetry) {
    principal_verts = sym_updater.get_principal(VERTS);
    faces_to_process = sym_updater.get_associated_elems(vert_faces);
    verts_to_update =
        SymmetricUpdater::get_included_verts(faces_to_process, geom.faces());
  }
  else { // not using_symmetry
    principal_verts =
        SymmetricUpdater::sequential_index_list(geom.verts().size());
    faces_to_process =
        SymmetricUpdater::sequential_index_list(geom.faces().size());
  }

  // Use oversized arrays to avoid mapping
  vector<Vec3d> offsets(verts.size()); // Vertex adjustments
  vector<Vec3d> norms(faces.size());   // Face normals
  vector<Vec3d> cents(faces.size());   // Face centroids

  double test_val = it_ctrl.get_test_val();
  double last_max_diff2 = 0.0;
  for (it_ctrl.start_iter(); !it_ctrl.is_done(); it_ctrl.next_iter()) {
    std::fill(offsets.begin(), offsets.end(), Vec3d::zero);

    if (using_symmetry) {
      // Ensure that the vertices used in adjustment are up to date
      for (auto v_idx : verts_to_update)
        sym_updater.update_from_principal_vertex(v_idx);
    }

    // Initialize face data for just the necessary faces
    for (auto f_idx : faces_to_process) {
      norms[f_idx] = geom.face_norm(f_idx).unit();
      cents[f_idx] = geom.face_cent(f_idx);
    }

    int cnt_proj = 0;
    int cnt_int = 0;
    double max_diff2 = 0.0;
    for (auto v_idx : principal_verts) {
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
        stat.set_ok(finish_msg);
      }
      else if (it_ctrl.is_last_iter()) {
        // reached last iteration without solving
        it_ctrl.set_finished();
        finish_msg = "not solved, test value not achieved";
        stat.set_warning(finish_msg);
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

  if (using_symmetry)
    base_geom = sym_updater.get_geom_final();

  return stat;
}

} // namespace anti
