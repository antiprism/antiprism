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

/* \file Coloring.cc
 * \brief classes to color all elements of a type.
 */

#include <algorithm>
#include <limits.h>
#include <limits>
#include <map>
#include <set>
#include <string.h>
#include <string>
#include <vector>

#include "boundbox.h"
#include "coloring.h"
#include "geometryinfo.h"
#include "mathutils.h"
#include "private_prop_col.h"
#include "random.h"
#include "utils.h"

using std::map;
using std::set;
using std::string;
using std::vector;

namespace anti {

Coloring::Coloring(Geometry *geo) : geom(geo), cycle_msecs(0) {}

Coloring::~Coloring() = default;

Coloring::Coloring(const Coloring & /*clrng*/) = default;

Coloring &Coloring::operator=(const Coloring &clrng)
{
  if (this != &clrng) {
    ColorMapMulti::operator=(clrng);
    geom = clrng.geom;
    cycle_msecs = clrng.cycle_msecs;
  }

  return *this;
}

void Coloring::cycle_map_cols() { set_shift(get_shift() + 1); }

void Coloring::set_all_idx_to_val(map<int, Color> &cols)
{
  map<int, Color>::iterator mi;
  for (mi = cols.begin(); mi != cols.end(); mi++)
    if (mi->second.is_index())
      mi->second = get_col(mi->second.get_index());
}

inline double fract(double rng[], double frac)
{
  return fmod(rng[0] + (rng[1] - rng[0]) * frac, 1 + epsilon);
}

Vec3d get_unsigned(Vec3d v)
{
  if (double_lt(v[2], 0))
    v = -v;
  else if (double_eq(v[2], 0)) {
    if (double_lt(v[1], 0))
      v = -v;
    else if (double_eq(v[1], 0)) {
      if (double_lt(v[0], 0))
        v = -v;
      else if (double_eq(v[0], 0))
        v = Vec3d::Z;
    }
  }
  return v;
}

int Coloring::z_gradient(Vec3d vec, Vec3d cent, double height, int def_sz)
{
  int sz = def_sz;
  const vector<ColorMap *> &cmaps = get_cmaps();
  if (cmaps.size() > 0) {
    if (cmaps[0]->effective_size() > 0)
      sz = cmaps[0]->effective_size();
    else
      sz = INT_MAX;
  }

  return (int)floor(sz * (0.5 * height + (vec - cent)[2]) / (height + epsilon));
}

void Coloring::setup_lights(Geometry &lts)
{
  if (lts.verts().size() == 0) {
    lts.add_vert(Vec3d(1, 0, 0), Vec3d(1, 0, 0));
    lts.add_vert(Vec3d(0, 1, 0), Vec3d(0, 1, 0));
    lts.add_vert(Vec3d(0, 0, 1), Vec3d(0, 0, 1));
    lts.add_vert(Vec3d(-1, 0, 0), Vec3d(0, 1, 1));
    lts.add_vert(Vec3d(0, -1, 0), Vec3d(1, 0, 1));
    lts.add_vert(Vec3d(0, 0, -1), Vec3d(1, 1, 0));
  }
  else
    for (unsigned int l = 0; l < lts.verts().size(); l++)
      lts.verts(l).to_unit();
}

Color Coloring::light(Vec3d vec, Geometry &lts)
{
  Vec3d col_sum(0, 0, 0);
  vec.to_unit();
  double dot;
  for (unsigned int l = 0; l < lts.verts().size(); l++)
    if ((dot = vdot(vec, lts.verts(l))) > 0)
      col_sum += dot * lts.colors(VERTS).get(l).get_vec3d();

  for (int j = 0; j < 3; j++)
    if (col_sum[j] > 1)
      col_sum[j] = 1;

  return Color(col_sum);
}

void Coloring::v_apply_cmap()
{
  set_all_idx_to_val(get_geom()->colors(VERTS).get_properties());
}

void Coloring::v_one_col(Color col)
{
  for (unsigned int i = 0; i < get_geom()->verts().size(); i++)
    get_geom()->colors(VERTS).set(i, col);
}

void Coloring::v_unique(bool apply_map)
{
  for (unsigned int i = 0; i < get_geom()->verts().size(); i++) {
    if (apply_map)
      get_geom()->colors(VERTS).set(i, get_col(i));
    else
      get_geom()->colors(VERTS).set(i, i);
  }
}

void Coloring::v_sets(const vector<set<int>> &equivs, bool apply_map)
{
  for (unsigned int i = 0; i < equivs.size(); i++) {
    for (auto si = equivs[i].begin(); si != equivs[i].end(); ++si) {
      if (apply_map)
        get_geom()->colors(VERTS).set(*si, get_col(i));
      else
        get_geom()->colors(VERTS).set(*si, i);
    }
  }
}

void Coloring::v_proper(bool apply_map)
{
  ProperColor prop(get_geom()->verts().size());

  for (unsigned int i = 0; i < get_geom()->faces().size(); ++i)
    for (unsigned int j = 0; j < get_geom()->faces(i).size(); ++j)
      prop.set_adj(get_geom()->faces(i, j), get_geom()->faces_mod(i, j + 1));

  prop.find_colors();
  for (unsigned int i = 0; i < get_geom()->verts().size(); i++)
    if (apply_map)
      get_geom()->colors(VERTS).set(i, get_col(prop.get_color(i)));
    else
      get_geom()->colors(VERTS).set(i, prop.get_color(i));
}

void Coloring::v_order(bool apply_map)
{
  vector<int> f_cnt(get_geom()->verts().size());
  for (unsigned int i = 0; i < get_geom()->faces().size(); i++)
    for (unsigned int j = 0; j < get_geom()->faces(i).size(); j++)
      f_cnt[get_geom()->faces(i, j)]++;
  for (unsigned int i = 0; i < get_geom()->verts().size(); i++)
    if (apply_map)
      get_geom()->colors(VERTS).set(i, get_col(f_cnt[i]));
    else
      get_geom()->colors(VERTS).set(i, f_cnt[i]);
}

void Coloring::v_position(bool apply_map)
{
  BoundBox bb(get_geom()->verts());
  Vec3d cent = bb.get_centre();
  double height = bb.get_max()[2] - bb.get_min()[2];
  for (unsigned int i = 0; i < get_geom()->verts().size(); i++) {
    int idx = z_gradient(get_geom()->verts(i), cent, height);
    if (apply_map)
      get_geom()->colors(VERTS).set(i, get_col(idx));
    else
      get_geom()->colors(VERTS).set(i, idx);
  }
}

void Coloring::v_avg_angle(bool apply_map)
{
  GeometryInfo info(*get_geom());
  Geometry tmp;
  tmp.raw_verts() = get_geom()->verts();
  for (int i = 0; i < info.num_verts(); i++)
    tmp.add_face(info.get_vert_figs()[i].size() ? info.get_vert_figs()[i][0]
                                                : vector<int>());
  Geometry *orig_geom = get_geom();
  set_geom(&tmp);
  f_avg_angle(apply_map);
  orig_geom->colors(VERTS) = tmp.colors(FACES);
  set_geom(orig_geom);
}

void Coloring::v_lights(Geometry lts)
{
  setup_lights(lts);
  Vec3d cent = get_geom()->centroid();
  for (unsigned int i = 0; i < get_geom()->verts().size(); i++)
    get_geom()->colors(VERTS).set(i, light(get_geom()->verts(i) - cent, lts));
}

void Coloring::face_edge_color(const vector<vector<int>> &elems,
                               const ElemProps<Color> &cols)
{
  vector<vector<int>> v_elems(get_geom()->verts().size());
  for (unsigned int i = 0; i < elems.size(); ++i)
    for (unsigned int j = 0; j < elems[i].size(); ++j)
      v_elems[elems[i][j]].push_back(i);

  for (unsigned int i = 0; i < v_elems.size(); ++i) {
    Vec4d col(0, 0, 0, 0);
    int val_cnt = 0;
    int first_idx = -1;
    for (int j : v_elems[i]) {
      Color ecol = cols.get(j);
      if (ecol.is_value()) {
        col += ecol.get_vec4d();
        val_cnt++;
      }
      else if (first_idx == -1 && ecol.is_index())
        first_idx = ecol.get_index();
    }
    if (val_cnt)
      get_geom()->colors(VERTS).set(i, Color(col / double(val_cnt)));
    else if (first_idx != -1)
      get_geom()->colors(VERTS).set(i, Color(first_idx));
  }
}

void Coloring::v_face_color()
{
  face_edge_color(get_geom()->faces(), get_geom()->colors(FACES));
}

void Coloring::v_edge_color()
{
  face_edge_color(get_geom()->edges(), get_geom()->colors(EDGES));
}

void Coloring::f_apply_cmap()
{
  set_all_idx_to_val(get_geom()->colors(FACES).get_properties());
}

void Coloring::f_one_col(Color col)
{
  for (unsigned int i = 0; i < get_geom()->faces().size(); i++)
    get_geom()->colors(FACES).set(i, col);
}

void Coloring::f_sets(const vector<set<int>> &equivs, bool apply_map)
{
  for (unsigned int i = 0; i < equivs.size(); i++) {
    Color col(i);
    for (auto si = equivs[i].begin(); si != equivs[i].end(); ++si) {
      if (apply_map)
        get_geom()->colors(FACES).set(*si, get_col(i));
      else
        get_geom()->colors(FACES).set(*si, i);
    }
  }
}

void Coloring::f_unique(bool apply_map)
{
  for (unsigned int i = 0; i < get_geom()->faces().size(); i++) {
    if (apply_map)
      get_geom()->colors(FACES).set(i, get_col(i));
    else
      get_geom()->colors(FACES).set(i, i);
  }
}

void Coloring::f_proper(bool apply_map)
{
  ProperColor prop(get_geom()->faces().size());

  auto ef_prs = get_geom()->get_edge_face_pairs(false);
  for (const auto &kp : ef_prs) {
    const auto &f_idxs = kp.second;
    for (unsigned int i = 0; i < f_idxs.size(); ++i)
      for (unsigned int j = i + 1; j < f_idxs.size(); ++j)
        prop.set_adj(f_idxs[i], f_idxs[j]);
  }

  prop.find_colors();
  for (unsigned int i = 0; i < get_geom()->faces().size(); i++)
    if (apply_map)
      get_geom()->colors(FACES).set(i, get_col(prop.get_color(i)));
    else
      get_geom()->colors(FACES).set(i, prop.get_color(i));
}

void Coloring::f_sides(bool apply_map)
{
  for (unsigned int i = 0; i < get_geom()->faces().size(); i++)
    if (apply_map)
      get_geom()->colors(FACES).set(i, get_col(get_geom()->faces(i).size()));
    else
      get_geom()->colors(FACES).set(i, get_geom()->faces(i).size());
}

void Coloring::f_avg_angle(bool apply_map)
{
  int faces_sz = get_geom()->faces().size();
  for (int i = 0; i < faces_sz; i++) {
    vector<double> f_angs;
    get_geom()->face_angles_lengths(i, &f_angs);
    double ang_sum = 0;
    for (double f_ang : f_angs)
      ang_sum += f_ang;
    // set invalid faces to have impossible angle 181
    int idx =
        f_angs.size() ? (int)(rad2deg(ang_sum / f_angs.size()) + 0.5) : 181;
    if (apply_map)
      get_geom()->colors(FACES).set(i, get_col(idx));
    else
      get_geom()->colors(FACES).set(i, idx);
  }
}

void Coloring::f_parts(bool apply_map)
{
  vector<vector<int>> parts;
  Geometry gtmp = *get_geom();
  gtmp.orient(&parts);
  for (unsigned int i = 0; i < parts.size(); i++)
    for (unsigned int j = 0; j < parts[i].size(); j++)
      if (apply_map)
        get_geom()->colors(FACES).set(parts[i][j], get_col(i));
      else
        get_geom()->colors(FACES).set(parts[i][j], i);
}

void Coloring::f_normal(bool apply_map)
{
  for (unsigned int i = 0; i < get_geom()->faces().size(); i++) {
    int idx = z_gradient(get_geom()->face_norm(i).unit());
    if (apply_map)
      get_geom()->colors(FACES).set(i, get_col(idx));
    else
      get_geom()->colors(FACES).set(i, idx);
  }
}

void Coloring::f_centroid(bool apply_map)
{
  BoundBox bb(get_geom()->verts());
  Vec3d cent = bb.get_centre();
  double height = bb.get_max()[2] - bb.get_min()[2];
  for (unsigned int i = 0; i < get_geom()->faces().size(); i++) {
    int idx = z_gradient(get_geom()->face_cent(i), cent, height);
    if (apply_map)
      get_geom()->colors(FACES).set(i, get_col(idx));
    else
      get_geom()->colors(FACES).set(i, idx);
  }
}

void Coloring::f_lights(Geometry lts)
{
  setup_lights(lts);
  for (unsigned int i = 0; i < get_geom()->faces().size(); i++)
    get_geom()->colors(FACES).set(i, light(get_geom()->face_norm(i), lts));
}

void Coloring::f_lights2(Geometry lts)
{
  setup_lights(lts);
  for (unsigned int i = 0; i < get_geom()->faces().size(); i++)
    get_geom()->colors(FACES).set(i, light(get_geom()->face_cent(i), lts));
}

void Coloring::e_apply_cmap()
{
  set_all_idx_to_val(get_geom()->colors(EDGES).get_properties());
}

void Coloring::e_one_col(Color col)
{
  for (unsigned int i = 0; i < get_geom()->edges().size(); i++)
    get_geom()->colors(EDGES).set(i, col);
}

void Coloring::e_sets(const vector<set<int>> &equivs, bool apply_map)
{
  for (unsigned int i = 0; i < equivs.size(); i++) {
    for (auto si = equivs[i].begin(); si != equivs[i].end(); ++si) {
      if (apply_map)
        get_geom()->colors(EDGES).set(*si, get_col(i));
      else
        get_geom()->colors(EDGES).set(*si, i);
    }
  }
}

void Coloring::e_unique(bool apply_map)
{
  for (unsigned int i = 0; i < get_geom()->edges().size(); i++) {
    if (apply_map)
      get_geom()->colors(EDGES).set(i, get_col(i));
    else
      get_geom()->colors(EDGES).set(i, i);
  }
}

void Coloring::e_proper(bool apply_map)
{
  auto edges = get_geom()->get_edge_face_pairs(false);
  map<vector<int>, vector<int>>::iterator mi;
  int idx = 0;
  for (mi = edges.begin(); mi != edges.end(); ++mi)
    mi->second[0] = idx++; // assign index numbers to the impl edges
  ProperColor prop(idx);   // idx is total number of impl edges

  vector<int> e_next(2);
  for (unsigned int i = 0; i < get_geom()->faces().size(); ++i) {
    for (unsigned int j = 0; j < get_geom()->faces(i).size(); ++j) {
      vector<int> e0 =
          make_edge(get_geom()->faces(i, j), get_geom()->faces_mod(i, j + 1));
      vector<int> e1 = make_edge(get_geom()->faces_mod(i, j + 1),
                                 get_geom()->faces_mod(i, j + 2));
      // An edge is adjacent to the edge that follows it on a face
      prop.set_adj(edges[e0][0], edges[e1][0]);
    }
  }

  prop.find_colors();
  for (mi = edges.begin(); mi != edges.end(); ++mi) {
    int e_idx = edges[mi->first][0];
    int col_idx = prop.get_color(e_idx);
    Color col = (apply_map) ? get_col(col_idx) : Color(col_idx);
    get_geom()->add_edge(mi->first[0], mi->first[1], col);
  }
}

void Coloring::e_face_color()
{
  const vector<vector<int>> &faces = get_geom()->faces();
  vector<vector<int>> efaces(get_geom()->edges().size());
  for (unsigned int i = 0; i < faces.size(); ++i) {
    for (unsigned int j = 0; j < faces[i].size(); ++j) {
      vector<int> edge =
          make_edge(faces[i][j], faces[i][(j + 1) % faces[i].size()]);
      auto ei = get_geom()->edges().begin();
      while ((ei = find(ei, get_geom()->edges().end(), edge)) !=
             get_geom()->edges().end()) {
        efaces[ei - get_geom()->edges().begin()].push_back(i);
        ++ei;
      }
    }
  }

  for (unsigned int i = 0; i < efaces.size(); ++i) {
    Vec4d col(0, 0, 0, 0);
    int val_cnt = 0;
    int first_idx = -1;
    for (int j : efaces[i]) {
      Color fcol = get_geom()->colors(FACES).get(j);
      if (fcol.is_value()) {
        col += fcol.get_vec4d();
        val_cnt++;
      }
      else if (first_idx == -1 && fcol.is_index())
        first_idx = fcol.get_index();
    }
    if (val_cnt)
      get_geom()->colors(EDGES).set(i, Color(col / double(val_cnt)));
    else if (first_idx != -1)
      get_geom()->colors(EDGES).set(i, Color(first_idx));
  }
}

void Coloring::e_parts(bool apply_map)
{
  GeometryInfo info(*get_geom());
  const vector<vector<int>> parts = info.get_edge_parts();
  for (unsigned int part_idx = 0; part_idx < parts.size(); part_idx++)
    for (unsigned int j = 0; j < parts[part_idx].size(); j++) {
      const int edge_idx = parts[part_idx][j];
      if (apply_map)
        get_geom()->colors(EDGES).set(edge_idx, get_col(part_idx));
      else
        get_geom()->colors(EDGES).set(edge_idx, part_idx);
    }
}

void Coloring::e_direction(bool apply_map)
{
  for (unsigned int i = 0; i < get_geom()->edges().size(); i++) {
    Vec3d v = (get_geom()->edge_vec(i)).with_len(2.0);
    if (v[2] < 0)
      v = -v;
    v -= Vec3d::Z; // put v[2] in the range -1.0 to 1.0;
    int idx = z_gradient(v);
    if (apply_map)
      get_geom()->colors(EDGES).set(i, get_col(idx));
    else
      get_geom()->colors(EDGES).set(i, idx);
  }
}

void Coloring::e_mid_point(bool apply_map)
{
  BoundBox bb(get_geom()->verts());
  Vec3d cent = bb.get_centre();
  double height = bb.get_max()[2] - bb.get_min()[2];
  for (unsigned int i = 0; i < get_geom()->edges().size(); i++) {
    int idx = z_gradient(get_geom()->edge_cent(i), cent, height);
    if (apply_map)
      get_geom()->colors(EDGES).set(i, get_col(idx));
    else
      get_geom()->colors(EDGES).set(i, idx);
  }
}

struct vec_less {
  bool operator()(const Vec3d &v1, const Vec3d &v2) const
  {
    return compare(v1, v2, epsilon) == -1;
  }
};

void Coloring::e_vector(bool apply_map)
{
  map<Vec3d, vector<int>, vec_less> dirs;
  for (unsigned int i = 0; i < get_geom()->edges().size(); i++)
    dirs[get_unsigned(get_geom()->edge_vec(i)).unit()].push_back(i);

  int idx = 0;
  map<Vec3d, vector<int>, vec_less>::const_iterator mi;
  for (mi = dirs.begin(); mi != dirs.end(); ++mi) {
    for (int i : mi->second) {
      if (apply_map)
        get_geom()->colors(EDGES).set(i, get_col(idx));
      else
        get_geom()->colors(EDGES).set(i, idx);
    }
    idx++;
  }
}

void Coloring::e_lengths(double min_diff, bool apply_map)
{
  if (get_geom()->edges().size() == 0)
    return;

  const Geometry &geom = *get_geom();
  vector<std::pair<double, int>> lens(geom.edges().size());
  for (unsigned int i = 0; i < geom.edges().size(); i++)
    lens[i] = {geom.edge_len(i), i};
  std::sort(lens.begin(), lens.end());

  vector<set<int>> equivs;
  // large negative will test different to first edge length
  double seq_first_len = -std::numeric_limits<double>::max();
  for (auto kp : lens) {
    if (!double_eq(kp.first, seq_first_len, min_diff)) {
      seq_first_len = kp.first;
      equivs.push_back(set<int>());
    }
    equivs.back().insert(kp.second);
  }

  e_sets(equivs, apply_map);
}

void Coloring::e_lights(Geometry lts)
{
  setup_lights(lts);
  Vec3d cent = get_geom()->centroid();
  for (unsigned int i = 0; i < get_geom()->edges().size(); i++)
    get_geom()->colors(EDGES).set(
        i, light(get_geom()->edge_nearpt(i, cent) - cent, lts));
}

void Coloring::e_dir_lights(Geometry lts)
{
  if (lts.verts().size() == 0) {
    lts.add_vert(Vec3d(0.57, 0.0, 0.3), Vec3d(1, 0, 0));
    lts.add_vert(Vec3d(-0.57, 0.0, -0.3), Vec3d(1, 0, 0));
    lts.add_vert(Vec3d(-0.28, 0.5, 0.3), Vec3d(0, 1, 0));
    lts.add_vert(Vec3d(0.28, 0.5, -0.3), Vec3d(0, 1, 0));
    lts.add_vert(Vec3d(-0.28, -0.5, 0.3), Vec3d(0, 0, 1));
    lts.add_vert(Vec3d(0.28, -0.5, -0.3), Vec3d(0, 0, 1));
  }
  setup_lights(lts);

  for (unsigned int i = 0; i < get_geom()->edges().size(); i++)
    get_geom()->colors(EDGES).set(
        i, light(get_unsigned(get_geom()->edge_vec(i)), lts));
}

void Coloring::vef_one_col(Color vert_col, Color edge_col, Color face_col)
{
  if (vert_col.is_set())
    v_one_col(vert_col);

  if (edge_col.is_set()) {
    geom->add_missing_impl_edges();
    e_one_col(edge_col);
  }

  if (face_col.is_set())
    f_one_col(face_col);
}

static bool get_cycle_rate(const char *str, double *cps)
{
  size_t len = strlen(str);
  if (len >= 3 && str[len - 2] == 'h' && str[len - 1] == 'z') {
    char str_copy[MSG_SZ];
    strncpy(str_copy, str, len - 2);
    str_copy[len - 2] = '\0';
    double cycs;
    if (read_double(str_copy, &cycs) && cycs >= 0.0) {
      *cps = cycs;
      return true;
    }
  }

  return false;
}

Status read_colorings(Coloring clrngs[], const char *line)
{
  char line_copy[MSG_SZ];
  strncpy(line_copy, line, MSG_SZ);
  line_copy[MSG_SZ - 1] = '\0';

  vector<char *> parts;
  int parts_sz = split_line(line_copy, parts, ",");

  vector<char *> map_names;
  Coloring clrng;
  unsigned int conv_elems = 0;

  Status stat;
  for (int i = 0; i < parts_sz; i++) {
    Status stat2;
    ColorMap *col_map = colormap_from_name(parts[i], &stat2);
    double cps;
    if (get_cycle_rate(parts[i], &cps)) {
      clrng.set_cycle_msecs((int)(1000 / cps));
      if (col_map)
        stat.set_warning(msg_str("cycle_rate '%s' is also a valid "
                                 "colour map name",
                                 parts[i]));
    }
    else if (strspn(parts[i], "vef") == strlen(parts[i])) {
      conv_elems |= 4 * (strchr(parts[i], 'f') != nullptr) +
                    2 * (strchr(parts[i], 'e') != nullptr) +
                    1 * (strchr(parts[i], 'v') != nullptr);
      if (col_map)
        stat.set_warning(msg_str("conversion elements '%s' is also a "
                                 "valid colour map name",
                                 parts[i]));
    }
    else if (col_map) {
      clrng.add_cmap(col_map);
      if (stat2.is_warning())
        stat.set_warning(msg_str("map %d: %s", i + 1, stat2.c_msg()));
    }
    else {
      return stat2; // from init ColorMap
    }
  }

  if (!conv_elems)
    conv_elems = 7; // no elements specified, use "fev"
  for (int i = 0; i < 3; i++) {
    if ((conv_elems & (1 << i)))
      clrngs[i] = clrng;
  }

  return stat;
}

} // namespace anti
