/*
   Copyright (c) 2003-2016, Adrian Rossiter

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
   Name: info_print.cc
   Description: information from OFF file - print functions
   Project: Antiprism - http://www.antiprism.com
*/

#include "rep_print.h"

#include <cstdio>
#include <cstring>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

using std::map;
using std::pair;
using std::set;
using std::string;
using std::vector;

using namespace anti;

Status rep_printer::set_sub_symmetry(const string &sub_sym)
{
  Symmetry sub;
  Status stat = get_symmetry().get_sub_sym(sub_sym, &sub);
  if (stat.is_ok())
    sub_sym_str = sub_sym;
  return stat;
}

string rep_printer::idx2s(int idx, int extra_sz)
{
  if (idx < extra_sz)
    return msg_str("%d", idx);
  else
    return msg_str("x%d", idx - extra_sz);
}

std::string rep_printer::col2s(Color col)
{
  string ret;
  if (col.is_index())
    ret = msg_str("%d", col.get_index());
  else if (col.is_value()) {
    ret = msg_str("%d %d %d", col[0], col[1], col[2]);
    if (col[3] < 255)
      ret += msg_str(" %d", col[3]);
  }

  return ret;
}

void rep_printer::general_sec()
{
  fprintf(ofile, "[general]\n");
  fprintf(ofile, "num_verts = %d\n", num_verts());
  fprintf(ofile, "num_faces = %d\n", num_faces());
  fprintf(ofile, "num_edges = %d\n", num_edges());
  fprintf(ofile, "vertex_centroid = (%s)\n", v2s(geom.centroid()).c_str());
  fprintf(ofile, "volume_centroid = ");
  if (geom.is_oriented() && is_polyhedron())
    fprintf(ofile, "(%s)", v2s(volume_centroid()).c_str());
  else
    fprintf(ofile, "n/a (calculated value: (%s))",
            v2s(volume_centroid()).c_str());
  fprintf(ofile, "\n");

  fprintf(ofile, "oriented = ");
  if (is_known_connectivity())
    fprintf(ofile, "%s", is_oriented() ? "yes" : "no");
  else
    fprintf(ofile, "n/a (calculated value: %s)", is_oriented() ? "yes" : "no");
  fprintf(ofile, "\n");

  fprintf(ofile, "orientable = ");
  if (is_known_connectivity())
    fprintf(ofile, "%s", is_orientable() ? "yes" : "no");
  else
    fprintf(ofile, "n/a (calculated value: %s)",
            is_orientable() ? "yes" : "no");
  fprintf(ofile, "\n");

  fprintf(ofile, "connectivity = %spolyhedron, %sclosed, %seven, %sknown\n",
          is_polyhedron() ? "" : "not ", is_closed() ? "" : "not ",
          is_even_connectivity() ? "" : "not ",
          is_known_connectivity() ? "" : "not ");

  fprintf(ofile, "num_parts = ");
  if (is_known_connectivity())
    fprintf(ofile, "%d", num_parts());
  else
    fprintf(ofile, "n/a (calculated value: %d)", num_parts());
  fprintf(ofile, "\n");

  fprintf(ofile, "genus = ");
  if (is_known_genus())
    fprintf(ofile, "%d", genus());
  else
    fprintf(ofile, "n/a");
  fprintf(ofile, "\n");

  fprintf(ofile, "area = %s\n", d2s(face_areas().sum).c_str());

  fprintf(ofile, "volume = ");
  if (geom.is_oriented() && is_polyhedron())
    fprintf(ofile, "%s", d2s(volume()).c_str());
  else
    fprintf(ofile, "n/a (calculated value: %s)", d2s(volume()).c_str());
  fprintf(ofile, "\n");

  fprintf(ofile, "\n");
}

void rep_printer::faces_sec()
{
  int sz = num_faces();
  fprintf(ofile, "[faces]\n");
  fprintf(ofile, "num_faces = %d\n", sz);
  fprintf(ofile, "area = %s\n", d2s((sz) ? face_areas().sum : 0.0).c_str());
  fprintf(ofile, "face_area_max = %s (%s)\n",
          d2s((sz) ? face_areas().max : 0.0).c_str(),
          (sz) ? d2s(face_areas().idx[ElementLimits::IDX_MAX]).c_str() : "n/a");
  fprintf(ofile, "face_area_min = %s (%s)\n",
          d2s((sz) ? face_areas().min : 0.0).c_str(),
          (sz) ? d2s(face_areas().idx[ElementLimits::IDX_MIN]).c_str() : "n/a");
  fprintf(ofile, "face_area_avg = %s\n",
          d2s((sz) ? face_areas().sum / sz : 0.0).c_str());
  fprintf(ofile, "volume = ");
  if (sz && is_closed())
    fprintf(ofile, "%s\n", d2s(volume()).c_str());
  else
    fprintf(ofile, "n/a (not closed)\n");
  fprintf(ofile, "isoperimetric_quotient = ");
  if (sz && is_closed())
    fprintf(ofile, "%s\n", d2s(isoperimetric_quotient()).c_str());
  else
    fprintf(ofile, "n/a (not closed)\n");

  double max_nonplanar = 0.0;
  double sum_nonplanar = 0.0;
  for (int i = 0; i < sz; i++) {
    double nonplanar = get_f_max_nonplanars()[i];
    sum_nonplanar += nonplanar;
    if (nonplanar > max_nonplanar)
      max_nonplanar = nonplanar;
  }
  fprintf(ofile, "maximum_nonplanarity = %s\n", d2s(max_nonplanar).c_str());
  fprintf(ofile, "average_nonplanarity = %s\n",
          d2s((sz) ? sum_nonplanar / sz : 0.0).c_str());
  fprintf(ofile, "\n");
}

void rep_printer::angles_sec()
{
  fprintf(ofile, "[angles]\n");
  fprintf(ofile, "angle_max = %s\n", d2s(rad2deg(angle_lims().max)).c_str());
  fprintf(ofile, "angle_min = %s\n", d2s(rad2deg(angle_lims().min)).c_str());
  fprintf(ofile, "angle_avg = %s\n",
          d2s(rad2deg(angle_lims().sum / num_angles())).c_str());
  fprintf(ofile, "angle_defect = %s\n", d2s(rad2deg(angle_defect())).c_str());
  fprintf(ofile, "\n");
}

void rep_printer::solid_angles_sec()
{
  fprintf(ofile, "[solid_angles]\n");
  fprintf(ofile, "dihed_angle_max = %s (%d, %d)\n",
          d2s(rad2deg(dihed_angle_lims().max)).c_str(),
          dihed_angle_lims().idx[ElementLimits::IDX_MAX],
          dihed_angle_lims().idx[ElementLimits::IDX_MAX2]);
  fprintf(ofile, "dihed_angle_min = %s (%d, %d)\n",
          d2s(rad2deg(dihed_angle_lims().min)).c_str(),
          dihed_angle_lims().idx[ElementLimits::IDX_MIN],
          dihed_angle_lims().idx[ElementLimits::IDX_MIN2]);
  fprintf(ofile, "dihed_angle_flattest = %s (%d, %d)\n",
          d2s(rad2deg(dihed_angle_lims().zero)).c_str(),
          dihed_angle_lims().idx[ElementLimits::IDX_ZERO],
          dihed_angle_lims().idx[ElementLimits::IDX_ZERO2]);
  fprintf(ofile, "solid_angle_max = %s [%sx4PI] (%d)\n",
          d2s(solid_angle_lims().max).c_str(),
          d2s(solid_angle_lims().max / (4 * M_PI)).c_str(),
          solid_angle_lims().idx[ElementLimits::IDX_MAX]);
  fprintf(ofile, "solid_angle_min = %s [%sx4PI] (%d)\n",
          d2s(solid_angle_lims().min).c_str(),
          d2s(solid_angle_lims().min / (4 * M_PI)).c_str(),
          solid_angle_lims().idx[ElementLimits::IDX_MIN]);
  fprintf(ofile, "\n");
}

void rep_printer::edges_sec()
{
  fprintf(ofile, "[edges]\n");
  fprintf(ofile, "num_edges = %d\n", num_edges());
  fprintf(ofile, "perimeter = %s\n", d2s(edge_length_lims().sum).c_str());
  fprintf(ofile, "edge_length_max = %s (%d,%d)\n",
          d2s(edge_length_lims().max).c_str(),
          edge_length_lims().idx[ElementLimits::IDX_MAX],
          edge_length_lims().idx[ElementLimits::IDX_MAX2]);
  fprintf(ofile, "edge_length_min = %s (%d,%d)\n",
          d2s(edge_length_lims().min).c_str(),
          edge_length_lims().idx[ElementLimits::IDX_MIN],
          edge_length_lims().idx[ElementLimits::IDX_MIN2]);
  fprintf(ofile, "edge_length_avg = %s\n",
          d2s(edge_length_lims().sum / num_edges()).c_str());
  fprintf(ofile, "\n");
}

void rep_printer::distances_sec()
{
  fprintf(ofile, "[distances]\n");
  fprintf(ofile, "given_center = (%s)\n", v2s(get_center()).c_str());
  fprintf(ofile, "vert_min = %s (%d)\n", d2s(vert_dist_lims().min).c_str(),
          vert_dist_lims().idx[ElementLimits::IDX_MIN]);
  fprintf(ofile, "vert_max = %s (%d)\n", d2s(vert_dist_lims().max).c_str(),
          vert_dist_lims().idx[ElementLimits::IDX_MAX]);
  fprintf(ofile, "vert_avg = %s\n",
          d2s(vert_dist_lims().sum / num_verts()).c_str());
  fprintf(ofile, "face_min = %s (%d)\n", d2s(face_dist_lims().min).c_str(),
          face_dist_lims().idx[ElementLimits::IDX_MIN]);
  fprintf(ofile, "face_max = %s (%d)\n", d2s(face_dist_lims().max).c_str(),
          face_dist_lims().idx[ElementLimits::IDX_MAX]);
  fprintf(ofile, "face_avg = %s\n",
          d2s(face_dist_lims().sum / num_faces()).c_str());
  fprintf(ofile, "edge_min = %s (%d,%d)\n", d2s(edge_dist_lims().min).c_str(),
          edge_dist_lims().idx[ElementLimits::IDX_MIN],
          edge_dist_lims().idx[ElementLimits::IDX_MIN2]);
  fprintf(ofile, "edge_max = %s (%d,%d)\n", d2s(edge_dist_lims().max).c_str(),
          edge_dist_lims().idx[ElementLimits::IDX_MAX],
          edge_dist_lims().idx[ElementLimits::IDX_MAX2]);
  fprintf(ofile, "edge_avg = %s\n",
          d2s(edge_dist_lims().sum / num_edges()).c_str());
  fprintf(ofile, "\n");
}

void rep_printer::symmetry()
{
  fprintf(ofile, "[symmetry]\n");
  fprintf(ofile, "type = %s\n", get_symmetry_type_name().c_str());

  fprintf(ofile, "realignments = %u fixed",
          (unsigned int)get_symmetry_autos().get_fixed().size());
  int free_rots = get_symmetry_autos().num_free_rots();
  if (free_rots == 1)
    fprintf(ofile, " x axial rotation");
  else if (free_rots == 3)
    fprintf(ofile, " x full rotation");
  int free_transls = get_symmetry_autos().num_free_transls();
  if (free_transls == 1)
    fprintf(ofile, " x axial translation");
  else if (free_transls == 2)
    fprintf(ofile, " x plane translation");
  else if (free_transls == 3)
    fprintf(ofile, " x space translation");
  fprintf(ofile, "\n");

  fprintf(ofile, "alignment_to_std =");
  Trans3d m = get_symmetry_alignment_to_std();

  for (int i = 0; i < 12; i++)
    fprintf(ofile, "%c%s", (i == 0) ? ' ' : ',', d2s(m[i]).c_str());
  fprintf(ofile, "\n\n");
  fprintf(ofile, "[symmetry_axes]\n");
  map<string, int> ax_cnts;
  const set<SymmetryAxis> &axes = get_symmetry_axes(); // get_symmetry_axes();
  set<SymmetryAxis>::const_iterator ai;
  for (ai = axes.begin(); ai != axes.end(); ++ai)
    ax_cnts[Symmetry(ai->get_sym_type(), ai->get_nfold()).get_symbol()]++;

  map<string, int>::const_iterator mi;
  for (mi = ax_cnts.begin(); mi != ax_cnts.end(); ++mi)
    fprintf(ofile, "%s,%d\n", mi->first.c_str(), mi->second);
  fprintf(ofile, "\n");

  fprintf(ofile, "[symmetry_subgroups_nonconjugate]\n");
  const set<Symmetry> &subs = get_symmetry_subgroups();
  set<Symmetry>::const_iterator si;
  map<string, int> subsym_cnts;
  for (si = subs.begin(); si != subs.end(); ++si)
    subsym_cnts[si->get_symbol()]++;

  for (mi = subsym_cnts.begin(); mi != subsym_cnts.end(); ++mi)
    fprintf(ofile, "%s,%d\n", mi->first.c_str(), mi->second);
  fprintf(ofile, "\n");
}

void rep_printer::vert_heights_cnts()
{
  fprintf(ofile, "[vert_heights_cnts]\n");
  map<double, double_range_cnt, AngleLess>::iterator mi;
  map<double, double_range_cnt, AngleLess> v_heights;
  const vector<Vec3d> &verts = geom.verts();
  unsigned int sz = verts.size();
  for (unsigned int i = 0; i < sz; i++) {
    double val = verts[i][2];
    mi = v_heights.find(val);
    if (mi == v_heights.end()) // new range
      v_heights[val] = double_range_cnt().update(val);
    else // existing range
      mi->second.update(val);
  }

  for (mi = v_heights.begin(); mi != v_heights.end(); ++mi)
    fprintf(ofile, "%s = %d\t(range +/- %s)\n", d2s(mi->second.mid()).c_str(),
            mi->second.cnt, d2s(mi->second.rad()).c_str());
  fprintf(ofile, "\n");
}

void rep_printer::edge_lengths_cnts()
{
  fprintf(ofile, "[edge_lengths_cnts]\n");
  const map<double, double_range_cnt, AngleLess> &edge_lengths =
      get_edge_lengths_by_size();
  map<double, double_range_cnt, AngleLess>::const_iterator ei;
  for (ei = edge_lengths.begin(); ei != edge_lengths.end(); ++ei)
    fprintf(ofile, "%s = %d\t(range +/- %s)\n", d2s(ei->second.mid()).c_str(),
            ei->second.cnt, d2s(ei->second.rad()).c_str());
  fprintf(ofile, "\n");
}

void rep_printer::dihedral_angles_cnts()
{
  fprintf(ofile, "[dihedral_angles_cnts]\n");
  const map<double, double_range_cnt, AngleLess> &dihedrals =
      get_dihedral_angles_by_size();
  map<double, double_range_cnt, AngleLess>::const_iterator di;
  for (di = dihedrals.begin(); di != dihedrals.end(); ++di)
    fprintf(ofile, "%s = %d\t(range +/- %s)\n",
            d2s(rad2deg(di->second.mid())).c_str(), di->second.cnt,
            d2s(rad2deg(di->second.rad())).c_str());
  fprintf(ofile, "\n");
}

void rep_printer::edge_faces_cnts()
{
  fprintf(ofile, "[edge_faces_cnts]\n");

  auto ef_prs = geom.get_edge_face_pairs(false);
  map<int, int> edge_faces_cnts;
  for (auto ei : geom.edges()) {
    auto e2fs = ef_prs.find(ei);
    if (e2fs != ef_prs.end())
      edge_faces_cnts[e2fs->second.size()]++;
  }

  for (auto sz2cnt : edge_faces_cnts)
    fprintf(ofile, "%d = %d\n", sz2cnt.first, sz2cnt.second);
  fprintf(ofile, "\n");
}

void rep_printer::solid_angles_cnts()
{
  fprintf(ofile, "[solid_angles_cnts]\n");
  const map<double, double_range_cnt, AngleLess> &solid_angs =
      get_solid_angles_by_size();
  map<double, double_range_cnt, AngleLess>::const_iterator si;
  for (si = solid_angs.begin(); si != solid_angs.end(); ++si)
    fprintf(ofile, "%s = %d\t(range +/- %s)\n", d2s(si->second.mid()).c_str(),
            si->second.cnt, d2s(si->second.rad()).c_str());
  fprintf(ofile, "\n");
}

void rep_printer::vert_order_cnts()
{
  fprintf(ofile, "[vert_order_cnts]\n");
  map<int, int>::iterator mi;
  map<int, int> cnts;
  const vector<vector<int>> &v_cons = get_vert_cons();
  for (const auto &v_con : v_cons) {
    mi = cnts.find(v_con.size());
    if (mi == cnts.end())
      cnts[v_con.size()] = 1;
    else
      mi->second += 1;
  }
  for (mi = cnts.begin(); mi != cnts.end(); ++mi)
    fprintf(ofile, "%d = %d\n", mi->first, mi->second);
  fprintf(ofile, "\n");
}

void rep_printer::face_sides_cnts()
{
  fprintf(ofile, "[face_sides_cnts]\n");
  map<int, int>::iterator mi;
  map<int, int> cnts;
  for (unsigned int i = 0; i < geom.faces().size(); i++) {
    mi = cnts.find(geom.faces(i).size());
    if (mi == cnts.end())
      cnts[geom.faces(i).size()] = 1;
    else
      mi->second += 1;
  }
  for (mi = cnts.begin(); mi != cnts.end(); ++mi)
    fprintf(ofile, "%d = %d\n", mi->first, mi->second);
  fprintf(ofile, "\n");
}

void rep_printer::face_angles_cnts()
{
  fprintf(ofile, "[face_angles_cnts]\n");
  const map<vector<double>, int, AngleVectLess> &face_angs =
      get_plane_angles_by_size();
  map<vector<double>, int, AngleVectLess>::const_iterator fi;
  for (fi = face_angs.begin(); fi != face_angs.end(); ++fi) {
    for (unsigned int i = 0; i < fi->first.size(); i++)
      fprintf(ofile, "%s%s", d2s(rad2deg(fi->first[i])).c_str(),
              (i < fi->first.size() - 1) ? "," : "");
    fprintf(ofile, " = %d\n", fi->second);
  }
  fprintf(ofile, "\n");
}

void rep_printer::face_winding_cnts(const vector<int> winding_numbers,
                                    const bool signing)
{
  map<pair<int, int>, int>::iterator mi;
  map<pair<int, int>, int> cnts;
  for (unsigned int i = 0; i < geom.faces().size(); i++) {
    pair<int, int> key;
    key.first = geom.faces(i).size();
    int wn = winding_numbers[i];
    if (!signing && (wn > key.first / 2))
      wn = key.first - wn;
    key.second = wn;
    mi = cnts.find(key);
    if (mi == cnts.end())
      cnts[key] = 1;
    else
      mi->second += 1;
  }

  fprintf(ofile, "[face_windings_cnts %s]\n",
          (signing ? "signed" : "unsigned"));

  for (mi = cnts.begin(); mi != cnts.end(); ++mi) {
    pair<int, int> key = mi->first;
    fprintf(ofile, "%d/%d = %d\n", key.first, key.second, mi->second);
  }
  fprintf(ofile, "\n");
}

void rep_printer::vertex_figure_winding_cnts()
{
  Geometry polygon;
  polygon.add_verts(geom.verts());

  map<pair<int, int>, int>::iterator mi;
  map<pair<int, int>, int> cnts;
  for (unsigned int i = 0; i < geom.verts().size(); i++) {
    const vector<vector<int>> &vfigs = get_vert_figs()[i];
    bool mult_vf = (vfigs.size() > 1) ? true : false;

    int fsz_total = 0;
    int winding = 0;
    for (const auto &vfig : vfigs) {
      polygon.add_face(vfig);
      int d = find_polygon_denominator_signed(polygon, 0, epsilon);
      polygon.clear(FACES);

      int fsz = vfig.size();
      // if multiple vertex figures, make negative when signed to arrive a
      // correct total winding
      if (mult_vf) {
        if (d > fsz / 2)
          d -= fsz;
      }
      fsz_total += fsz;
      winding += d;
    }

    // in case this ever happens
    if (mult_vf && winding < 0)
      winding += fsz_total;

    pair<int, int> key;
    key.first = fsz_total;
    key.second = winding;
    mi = cnts.find(key);
    if (mi == cnts.end())
      cnts[key] = 1;
    else
      mi->second += 1;
  }

  fprintf(ofile, "[vertex_figure_windings_cnts signed]\n");

  for (mi = cnts.begin(); mi != cnts.end(); ++mi) {
    pair<int, int> key = mi->first;
    fprintf(ofile, "%d/%d = %d\n", key.first, key.second, mi->second);
  }
  fprintf(ofile, "\n");
}

void rep_printer::windings()
{
  // get signed winding number
  vector<int> winding_numbers;
  for (unsigned int i = 0; i < geom.faces().size(); i++)
    winding_numbers.push_back(
        find_polygon_denominator_signed(geom, i, epsilon));

  face_winding_cnts(winding_numbers, false);
  face_winding_cnts(winding_numbers, true);

  vertex_figure_winding_cnts();
}

void rep_printer::sym_orbit_cnts()
{
  fprintf(ofile, "[sym_orbit_cnts]\n");
  // Get symmetry if necessary
  Symmetry sym;
  if (sub_sym_str == "") { // use full symmetry
    sym = get_symmetry();
  }
  else { // use subsymmetry
    Status stat = get_symmetry().get_sub_sym(sub_sym_str, &sym);
    if (stat.is_error()) {
      fprintf(ofile, "invalid subsymmetry '%s': %s\n\n", sub_sym_str.c_str(),
              stat.c_msg());
      return;
    }
  }

  const char *elems[] = {"verts", "edges", "faces"};
  vector<vector<set<int>>> sym_equivs;
  get_equiv_elems(geom, sym.get_trans(), &sym_equivs);
  for (int i = 0; i < 3; i++) {
    string cnt_list;
    vector<set<int>>::const_iterator vi;
    for (vi = sym_equivs[i].begin(); vi != sym_equivs[i].end(); ++vi)
      cnt_list += msg_str("%d, ", vi->size());
    fprintf(ofile, "%s: %u ", elems[i], (unsigned int)sym_equivs[i].size());
    if (cnt_list.size())
      fprintf(ofile, " (%s)", cnt_list.substr(0, cnt_list.size() - 2).c_str());
    fprintf(ofile, "\n");
  }
  fprintf(ofile, "\n");
}

void rep_printer::v_index(int v_idx)
{
  fprintf(ofile, "%s", vidx2s(v_idx).c_str());
}

void rep_printer::v_coords(int v_idx)
{
  fprintf(ofile, "%s", v2s(geom.verts(v_idx)).c_str());
}

void rep_printer::v_neighbours(int v_idx)
{
  const vector<int> &vcons = get_vert_cons()[v_idx];
  for (unsigned int i = 0; i < vcons.size(); i++)
    fprintf(ofile, "%s%s", vidx2s(vcons[i]).c_str(),
            (i < vcons.size() - 1) ? " " : "");
}

void rep_printer::v_figure(int v_idx)
{
  const vector<vector<int>> &vfigs = get_vert_figs()[v_idx];
  for (unsigned int i = 0; i < vfigs.size(); i++) {
    if (i > 0) // print circuit separator
      fprintf(ofile, ":");
    for (unsigned int j = 0; j < vfigs[i].size(); j++)
      fprintf(ofile, "%s%s", vidx2s(vfigs[i][j]).c_str(),
              (j < vfigs[i].size() - 1) ? " " : "");
  }
}

void rep_printer::v_face_idxs(int v_idx)
{
  const vector<int> &v_faces = get_vert_faces()[v_idx];
  for (size_t i = 0; i < v_faces.size(); i++)
    fprintf(ofile, "%s%s", fidx2s(v_faces[i]).c_str(),
            (i < v_faces.size() - 1) ? " " : "");
}

void rep_printer::v_solid_angle(int v_idx)
{
  fprintf(ofile, "%s", d2s(get_vert_solid_angles()[v_idx]).c_str());
}

void rep_printer::v_order(int v_idx)
{
  fprintf(ofile, "%lu", (unsigned long)get_vert_cons()[v_idx].size());
}

void rep_printer::v_distance(int v_idx)
{
  fprintf(ofile, "%s", d2s((geom.verts(v_idx) - get_center()).len()).c_str());
}

void rep_printer::v_angles(int v_idx)
{
  const Geometry &dual = get_dual();
  const vector<int> &fcons = dual.faces(v_idx);
  for (unsigned int i = 0; i < fcons.size(); i++) {
    pair<int, int> vf_pr(v_idx, fcons[i]);
    fprintf(ofile, "%s%s", d2s(rad2deg(get_plane_angles().at(vf_pr))).c_str(),
            (i < fcons.size() - 1) ? " " : "");
  }
}

void rep_printer::v_color(int v_idx)
{
  fprintf(ofile, "%s", col2s(geom.colors(VERTS).get(v_idx)).c_str());
}

void rep_printer::e_index(int e_idx)
{
  fprintf(ofile, "%s", eidx2s(e_idx).c_str());
}

void rep_printer::e_vert_idxs(int e_idx)
{
  vector<int> edge = geom.edges(e_idx);
  fprintf(ofile, "%s %s", vidx2s(edge[0]).c_str(), vidx2s(edge[1]).c_str());
}

void rep_printer::e_face_idxs(int e_idx)
{
  vector<int> edge = geom.edges(e_idx);
  auto ei = get_edge_face_pairs().find(edge);
  if (ei != get_edge_face_pairs().end()) {
    auto fidxs = ei->second;
    for (unsigned int i = 0; i < fidxs.size(); i++) {
      if (fidxs[i] >= 0)
        fprintf(ofile, "%s%s", fidx2s(fidxs[i]).c_str(),
                (i < fidxs.size() - 1) ? " " : "");
    }
  }
}

void rep_printer::e_dihedral_angle(int e_idx)
{
  fprintf(ofile, "%s", d2s(rad2deg(get_edge_dihedrals()[e_idx])).c_str());
}

void rep_printer::e_central_angle(int e_idx)
{
  Vec3d v0 = geom.verts(geom.edges(e_idx, 0)) - get_center();
  Vec3d v1 = geom.verts(geom.edges(e_idx, 1)) - get_center();
  fprintf(
      ofile, "%s",
      d2s(rad2deg(acos(safe_for_trig(vdot(v0.unit(), v1.unit()))))).c_str());
}

void rep_printer::e_distance(int e_idx)
{
  Vec3d v0 = geom.verts(geom.edges(e_idx, 0)) - get_center();
  Vec3d v1 = geom.verts(geom.edges(e_idx, 1)) - get_center();
  double dist = (nearest_point(get_center(), v0, v1) - get_center()).len();
  fprintf(ofile, "%s", d2s(dist).c_str());
}

void rep_printer::e_centroid(int e_idx)
{
  fprintf(ofile, "%s", v2s(geom.edge_cent(e_idx)).c_str());
}

void rep_printer::e_direction(int e_idx)
{
  Vec3d v0 = geom.verts(geom.edges(e_idx, 0));
  Vec3d v1 = geom.verts(geom.edges(e_idx, 1));
  fprintf(ofile, "%s", v2s((v1 - v0).unit()).c_str());
}

void rep_printer::e_length(int e_idx)
{
  Vec3d v0 = geom.verts(geom.edges(e_idx, 0));
  Vec3d v1 = geom.verts(geom.edges(e_idx, 1));
  fprintf(ofile, "%s", d2s((v1 - v0).len()).c_str());
}

void rep_printer::e_color(int e_idx)
{
  fprintf(ofile, "%s", col2s(geom.colors(EDGES).get(e_idx)).c_str());
}

void rep_printer::f_index(int f_idx)
{
  fprintf(ofile, "%s", fidx2s(f_idx).c_str());
}

void rep_printer::f_vert_idxs(int f_idx)
{
  const vector<int> &face = geom.faces(f_idx);
  for (unsigned int i = 0; i < face.size(); i++)
    fprintf(ofile, "%s%s", vidx2s(face[i]).c_str(),
            (i < face.size() - 1) ? " " : "");
}

void rep_printer::f_neighbours(int f_idx)
{
  map<vector<int>, vector<int>>::const_iterator ei;
  vector<int> edge(2);
  const vector<int> &face = geom.faces(f_idx);
  unsigned int sz = face.size();
  for (unsigned int i = 0; i < sz; i++) {
    ei = get_edge_face_pairs().find(make_edge(face[i], face[(i + 1) % sz]));
    int neigh = (ei->second[0] != f_idx) ? ei->second[0] : ei->second[1];
    fprintf(ofile, "%s%s", fidx2s(neigh).c_str(), (i < sz - 1) ? " " : "");
  }
}

void rep_printer::f_normal(int f_idx)
{
  fprintf(ofile, "%s", v2s(geom.face_norm(f_idx).unit()).c_str());
}

void rep_printer::f_angles(int f_idx)
{
  vector<double> angs;
  geom.face_angles_lengths(f_idx, &angs);
  for (unsigned int i = 0; i < angs.size(); i++)
    fprintf(ofile, "%s%s", d2s(rad2deg(angs[i])).c_str(),
            (i < angs.size() - 1) ? " " : "");
}

void rep_printer::f_sides(int f_idx)
{
  fprintf(ofile, "%lu", (unsigned long)geom.faces(f_idx).size());
}

void rep_printer::f_distance(int f_idx)
{
  double dist = (nearest_point(get_center(), geom.verts(), geom.faces(f_idx)) -
                 get_center())
                    .len();
  fprintf(ofile, "%s", d2s(dist).c_str());
}

void rep_printer::f_area(int f_idx)
{
  fprintf(ofile, "%s", d2s(get_f_areas()[f_idx]).c_str());
}

void rep_printer::f_perimeter(int f_idx)
{
  fprintf(ofile, "%s", d2s(get_f_perimeters()[f_idx]).c_str());
}

void rep_printer::f_centroid(int f_idx)
{
  fprintf(ofile, "%s", v2s(geom.face_cent(f_idx)).c_str());
}

void rep_printer::f_lengths(int f_idx)
{
  const vector<Vec3d> &verts = geom.verts();
  const vector<int> &face = geom.faces(f_idx);
  unsigned int sz = face.size();
  for (unsigned int i = 0; i < sz; i++)
    fprintf(ofile, "%s%s",
            d2s((verts[face[i]] - verts[face[(i + 1) % sz]]).len()).c_str(),
            (i < face.size() - 1) ? " " : "");
}

void rep_printer::f_max_nonplanar(int f_idx)
{
  fprintf(ofile, "%s", d2s(get_f_max_nonplanars()[f_idx]).c_str());
}

void rep_printer::f_color(int f_idx)
{
  fprintf(ofile, "%s", col2s(geom.colors(FACES).get(f_idx)).c_str());
}
