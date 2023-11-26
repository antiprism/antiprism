/*
   Copyright (c) 2003-2023, Adrian Rossiter, Roger Kaufman

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
   Name: color_common.h
   Description: color code shared in /src
   Project: Antiprism - http://www.antiprism.com
*/

#include "color_common.h"
#include "../base/antiprism.h"

#include <cstdio>
#include <set>
#include <string>
#include <vector>

using std::pair;
using std::set;
using std::string;
using std::vector;

using namespace anti;

void OffColor::set_v_col_op(char v)
{
  if (upper)
    v = toupper(v);
  v_col_op = v;
}
void OffColor::set_v_col(const anti::Color v)
{
  v_col = v;
  v_col_op = val_op;
}
void OffColor::set_v_sub_sym(const string v) { v_sub_sym = v; }

void OffColor::set_e_col_op(char e)
{
  if (upper)
    e = toupper(e);
  e_col_op = e;
}
void OffColor::set_e_col(const anti::Color e)
{
  e_col = e;
  e_col_op = val_op;
}
void OffColor::set_e_sub_sym(const string e) { e_sub_sym = e; }
void OffColor::set_e_min_len_diff(const double e) { e_min_len_diff = e; }

void OffColor::set_f_col_op(char f)
{
  if (upper)
    f = toupper(f);
  f_col_op = f;
}
void OffColor::set_f_col(const anti::Color f)
{
  f_col = f;
  f_col_op = val_op;
}
void OffColor::set_f_sub_sym(const string f) { f_sub_sym = f; }

char OffColor::get_v_col_op() { return v_col_op; }
string OffColor::get_v_sub_sym() { return v_sub_sym; }
anti::Color OffColor::get_v_col() { return v_col; }

char OffColor::get_e_col_op() { return e_col_op; }
anti::Color OffColor::get_e_col() { return e_col; }
string OffColor::get_e_sub_sym() { return e_sub_sym; }
double OffColor::get_e_min_len_diff() { return e_min_len_diff; }

char OffColor::get_f_col_op() { return f_col_op; }
anti::Color OffColor::get_f_col() { return f_col; }
string OffColor::get_f_sub_sym() { return f_sub_sym; }

/* not sure about this code. direct access for now
void OffColor::set_clrng(anti::Coloring &clrng, int n) { clrngs[n] = clrng; }

anti::Coloring OffColor::get_clrng(int n) { return clrngs[n]; }
*/

bool OffColor::v_op_check(char *v_col_op, const char *op_str)
{
  return ((strlen(v_col_op) == 1) && (strspn(v_col_op, op_str)));
}

bool OffColor::e_op_check(char *e_col_op, const char *op_str)
{
  return ((strlen(e_col_op) == 1) && (strspn(e_col_op, op_str)));
}

bool OffColor::f_op_check(char *f_col_op, const char *op_str)
{
  return ((strlen(f_col_op) == 1) && (strspn(f_col_op, op_str)));
}

Status OffColor::off_color_main(Geometry &geom)
{
  Status stat;
  Geometry lights; // unset

  // only simulating -e
  if (e_col_op)
    geom.add_missing_impl_edges();

  // Get symmetry if necessary
  Symmetry sym;
  vector<vector<set<int>>> sym_equivs;
  if ((f_col_op && strchr("sS", f_col_op)) ||
      (e_col_op && strchr("sS", e_col_op)) ||
      (v_col_op && strchr("sS", v_col_op))) {
    sym.init(geom, &sym_equivs);
    v_equivs = sym_equivs[0];
    e_equivs = sym_equivs[1];
    f_equivs = sym_equivs[2];

    Symmetry sub;
    if (v_col_op && strchr("sS", v_col_op)) {
      if (!(stat = sym.get_sub_sym(v_sub_sym, &sub)))
        return stat.set_error(msg_str("get_sub_sym 'v': %s", stat.c_msg()));
      get_equiv_elems(geom, sub.get_trans(), &sym_equivs);
      v_equivs = sym_equivs[0];
    }
    if (e_col_op && strchr("sS", e_col_op)) {
      if (!(stat = sym.get_sub_sym(e_sub_sym, &sub)))
        return stat.set_error(msg_str("get_sub_sym 'e': %s", stat.c_msg()));
      get_equiv_elems(geom, sub.get_trans(), &sym_equivs);
      e_equivs = sym_equivs[1];
    }
    if (f_col_op && strchr("sS", f_col_op)) {
      if (!(stat = sym.get_sub_sym(f_sub_sym, &sub)))
        return stat.set_error(msg_str("get_sub_sym 'f': %s", stat.c_msg()));
      get_equiv_elems(geom, sub.get_trans(), &sym_equivs);
      f_equivs = sym_equivs[2];
    }
  }

  Coloring &fc = clrngs[FACES];
  fc.set_geom(&geom);
  if (f_col_op) {
    char op = f_col_op;
    ColorMap *cmap = nullptr;
    if (fc.get_cmaps().size() == 0) {
      if (strchr("GgCc", op))
        cmap = colormap_from_name("range");
      else
        cmap = colormap_from_name("spread");
    }
    if (cmap)
      fc.add_cmap(cmap);
    if (op == val_op)
      fc.f_one_col(f_col);
    else if (strchr("uU", op))
      fc.f_unique(op == 'U');
    else if (strchr("pP", op))
      fc.f_proper(op == 'P');
    else if (strchr("sS", op))
      fc.f_sets(sym_equivs[2], op == 'S');
    else if (strchr("nN", op))
      fc.f_sides(op == 'N');
    else if (strchr("aA", op))
      fc.f_avg_angle(op == 'A');
    else if (strchr("kK", op))
      fc.f_parts(op == 'K');
    else if (strchr("Gg", op))
      fc.f_normal(op == 'G');
    else if (strchr("Cc", op))
      fc.f_centroid(op == 'C');
    else if (strchr("L", op))
      fc.f_lights(lights);
    else if (strchr("l", op))
      fc.f_lights2(lights);
    else if (strchr("M", op))
      fc.f_apply_cmap();
  }

  Coloring &ec = clrngs[EDGES];
  ec.set_geom(&geom);
  if (e_col_op) {
    char op = e_col_op;
    ColorMap *cmap = nullptr;
    if (ec.get_cmaps().size() == 0) {
      if (strchr("GgCc", op))
        cmap = colormap_from_name("range");
      else
        cmap = colormap_from_name("spread");
    }
    if (cmap)
      ec.add_cmap(cmap);
    if (op == val_op)
      ec.e_one_col(e_col);
    else if (strchr("uU", op))
      ec.e_unique(op == 'U');
    else if (strchr("pP", op))
      ec.e_proper(op == 'P');
    else if (strchr("sS", op))
      ec.e_sets(sym_equivs[1], op == 'S');
    else if (strchr("nN", op))
      ec.e_order(op == 'N');
    else if (strchr("jJ", op))
      ec.e_lengths(get_e_min_len_diff(), op == 'J');
    else if (strchr("kK", op))
      ec.e_parts(op == 'K');
    else if (strchr("Gg", op))
      ec.e_direction(op == 'G');
    else if (strchr("Cc", op))
      ec.e_mid_point(op == 'C');
    else if (strchr("Dd", op))
      ec.e_vector(op == 'D');
    else if (strchr("L", op))
      ec.e_lights(lights);
    else if (strchr("l", op))
      ec.e_dir_lights(lights);
    else if (strchr("M", op))
      ec.e_apply_cmap();
  }

  Coloring &vc = clrngs[VERTS];
  vc.set_geom(&geom);
  if (v_col_op) {
    char op = v_col_op;
    ColorMap *cmap = nullptr;
    if (vc.get_cmaps().size() == 0) {
      if (strchr("Cc", op))
        cmap = colormap_from_name("range");
      else
        cmap = colormap_from_name("spread");
    }
    if (cmap)
      vc.add_cmap(cmap);
    if (op == val_op)
      vc.v_one_col(v_col);
    else if (strchr("uU", op))
      vc.v_unique(op == 'U');
    else if (strchr("pP", op))
      vc.v_proper(op == 'P');
    else if (strchr("sS", op))
      vc.v_sets(sym_equivs[0], op == 'S');
    else if (strchr("nN", op))
      vc.v_order(op == 'N');
    else if (strchr("aA", op))
      vc.v_avg_angle(op == 'A');
    else if (strchr("cC", op))
      vc.v_position(op == 'C');
    else if (strchr("L", op))
      vc.v_lights(lights);
    else if (strchr("M", op))
      vc.v_apply_cmap();
  }

  // Average colour values from adjoining elements after converting
  // index numbers
  if (f_col_op == 'V')
    fc.f_from_adjacent(VERTS);
  else if (f_col_op == 'E')
    fc.f_from_adjacent(EDGES);

  if (e_col_op == 'F')
    ec.e_from_adjacent(FACES);
  else if (e_col_op == 'V')
    ec.e_from_adjacent(VERTS);

  if (v_col_op == 'E')
    vc.v_from_adjacent(EDGES);
  else if (v_col_op == 'F')
    vc.v_from_adjacent(FACES);

  return Status::ok();
}

Status apply_transparency(Geometry &geom, const int opacity, const int elem)
{
  return (Coloring(&geom).apply_transparency(opacity, elem));
}

void apply_transparencies(Geometry &geom, const int (&opacity)[3])
{
  Status stat;

  for (int i = 0; i < 3; i++) {
    if (opacity[i] == -1)
      continue;
    stat = apply_transparency(geom, opacity[i], i);
    if (!stat.is_ok()) {
      string elem_str;
      if (i == VERTS)
        elem_str = "vertices";
      else if (i == EDGES)
        elem_str = "edges";
      else if (i == FACES)
        elem_str = "faces";
      fprintf(stderr, "apply transparency (%s): %s\n", elem_str.c_str(),
              stat.c_msg());
    }
  }
}

ColorMapMap *alloc_no_colorMap()
{
  auto *col_map = new ColorMapMap;

  // index INT_MAX as unset edge color
  col_map->set_col(0, Color::maximum_index);
  // this doesn't work because 'no color' is not allowed in a map
  // col_map->set_col(0, Color());  // unset edge color

  col_map->set_wrap();

  return col_map;
}

// color edges by dihedral angle
void color_edges_by_dihedral(Geometry &geom, Coloring &clrng, bool apply_map,
                             double eps)
{
  GeometryInfo info(geom);

  // if negative volume, orientation reversed, dihedral angles opposite
  bool reverse = (GeometryInfo(geom).volume() < 0) ? true : false;

  ColorMapMulti edge_map = clrng;

  // color all elements convex
  geom.add_missing_impl_edges();
  Coloring(&geom).e_one_col(apply_map ? edge_map.get_col(0) : 0);

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
    if (coplanar_found)
      geom.colors(EDGES).set(i, (apply_map ? edge_map.get_col(1) : 1));
    else if (!convexity)
      geom.colors(EDGES).set(i, (apply_map ? edge_map.get_col(2) : 2));
  }
}

// color faces by convexity compared to other faces
void color_faces_by_convexity(Geometry &geom, Coloring &clrng, bool apply_map,
                              double eps)
{
  GeometryInfo info(geom);

  // if negative volume, orientation reversed, dihedral angles opposite
  bool reverse = (GeometryInfo(geom).volume() < 0) ? true : false;

  // color all elements convex
  geom.add_missing_impl_edges();
  Coloring(&geom).f_one_col(apply_map ? clrng.get_col(0) : 0);

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
      geom.colors(FACES).set(i, (apply_map ? clrng.get_col(1) : 1));
    else if (!convexity)
      geom.colors(FACES).set(i, (apply_map ? clrng.get_col(2) : 2));
  }
}

// uses kis and alters geom
void color_faces_by_connection(Geometry &geom, Coloring &clrng, bool apply_map)
{
  // need edges or they'll turn invisible
  geom.add_missing_impl_edges();

  // do kis operation
  Geometry kis;
  wythoff_make_tiling(kis, geom, "k", true, false);
  // remove digons
  vector<int> dels;
  for (unsigned int i = 0; i < kis.faces().size(); i++) {
    if (kis.faces(i).size() < 3)
      dels.push_back((int)i);
  }
  kis.del(FACES, dels);
  kis.orient(1); // positive orientation

  // color edges and vertices invisible
  kis.add_missing_impl_edges();
  Coloring(&kis).vef_one_col(Color::invisible, Color::invisible, Color());

  // transfer edge and vertex colors from geom
  // if elements were invisible, mark them maximum
  // only non-counted element will remain invisible
  for (unsigned int i = 0; i < kis.edges().size(); i++) {
    unsigned int v_idx[2];
    Vec3d v[2];
    for (unsigned int j = 0; j < 2; j++) {
      v_idx[j] = kis.edges(i)[j];
      v[j] = kis.verts(v_idx[j]);
    }
    int geom_edge_no = find_edge_by_coords(geom, v[0], v[1], anti::epsilon);
    Color col;
    if (geom_edge_no > -1) {
      col = geom.colors(EDGES).get(geom_edge_no);
      if (col.is_invisible())
        col = Color::maximum_index;
      else if (!col.is_set())
        col = Color();
      kis.colors(EDGES).set(i, col);
    }
    for (unsigned int j = 0; j < 2; j++) {
      int ev = kis.edges(i)[j];
      int geom_v_idx =
          find_vert_by_coords(geom, kis.verts()[ev], anti::epsilon);
      if (geom_v_idx > -1) {
        col = geom.colors(VERTS).get(geom_v_idx);
        if (col.is_invisible())
          col = Color::maximum_index;
        else if (!col.is_set())
          col = Color();
        kis.colors(VERTS).set(ev, col);
      }
    }
  }

  // color kis faces by connection
  for (unsigned int i = 0; i < kis.faces().size(); i++) {
    vector<int> face = kis.faces()[i];
    unsigned int fsz = face.size();
    // face to face
    // faces connected by invisible edges are ignored
    int connections = 0;
    for (unsigned int j = 0; j < fsz; j++) {
      int v1 = face[j];
      int v2 = face[(j + 1) % fsz];
      vector<int> edge = make_edge(v1, v2);
      vector<int> face_idx = find_faces_with_edge(kis.faces(), edge);
      int edge_no = find_edge_in_edge_list(kis.edges(), edge);
      if (!(kis.colors(EDGES).get(edge_no)).is_invisible())
        connections += face_idx.size();
    }
    kis.colors(FACES).set(
        i, (apply_map ? clrng.get_col(connections) : connections));
  }

  // any elements marked maximum were invisible, set them back
  for (unsigned int i = 0; i < kis.edges().size(); i++) {
    Color col = kis.colors(EDGES).get(i);
    if (col.is_maximum_index()) {
      col = Color::invisible;
      kis.colors(EDGES).set(i, col);
    }
    for (unsigned int j = 0; j < 2; j++) {
      unsigned int v_idx = kis.edges(i)[j];
      col = kis.colors(VERTS).get(v_idx);
      if (col.is_maximum_index()) {
        col = Color::invisible;
        kis.colors(VERTS).set(v_idx, col);
      }
    }
  }

  geom = kis;
}

// duplicate code from stellate and miller
// color_faces_by_connection(
Status color_faces_by_connection_vef(Geometry &geom, OffColor &off_color)
{
  char op = off_color.get_f_col_op();

  // geom is built with face colors from the diagram, if 'q' do nothing
  Geometry geom_save;
  if (op && strchr("hH", op)) {
    color_faces_by_connection(geom, off_color.clrngs[FACES], (op == 'H'));
    // save geom to reassert invisible elements from kis operation
    geom_save = geom;
  }

  // any other color options done by class
  // will color edges invisible if so set
  Status stat;
  if (!(stat = off_color.off_color_main(geom)))
    return Status::error(stat.msg());

  // reassert invisible elements from kis operation
  if (op && strchr("hH", op)) {
    for (unsigned int i = 0; i < geom.edges().size(); i++) {
      unsigned int v_idx[2];
      Vec3d v[2];
      for (unsigned int j = 0; j < 2; j++) {
        v_idx[j] = geom.edges(i)[j];
        v[j] = geom.verts(v_idx[j]);
      }
      Color col;
      int save_edge_no =
          find_edge_by_coords(geom_save, v[0], v[1], anti::epsilon);
      if (save_edge_no > -1) {
        col = geom_save.colors(EDGES).get(save_edge_no);
        if (col.is_invisible())
          geom.colors(EDGES).set(i, Color::invisible);
      }
      for (unsigned int j = 0; j < 2; j++) {
        int ev = geom.edges(i)[j];
        int save_idx =
            find_vert_by_coords(geom_save, geom.verts()[ev], anti::epsilon);
        col = geom_save.colors(VERTS).get(save_idx);
        if (col.is_invisible())
          geom.colors(VERTS).set(ev, col);
      }
    }
  }

  return Status::ok();
}

// for lat_util.cc, bravais.cc and waterman.cc
// Rotational octahedral by Adrian Rossiter
Vec3d sort_Vec3d_chiral(const Vec3d &v, const double eps)
{
  Vec3d c = v;
  // Rotate into positive octant
  if (c[0] < 0) {
    c[0] = -c[0];
    c[2] = -c[2];
  }
  if (c[1] < 0) {
    c[1] = -c[1];
    c[2] = -c[2];
  }
  if (c[2] < 0) {
    std::swap(c[0], c[1]);
    c[2] = -c[2];
  }

  // if c[1] is maximum rotate to first place: 1,2,0
  if (c[1] > c[0] && c[1] > c[2] - eps)
    c = Vec3d(c[1], c[2], c[0]);
  else
    // if c[2] is maximum rotate to first place: 2,0,1
    if (c[2] > c[0] && c[2] > c[1] + eps)
      c = Vec3d(c[2], c[0], c[1]);
    else
      // if c[0] is maximum do nothing
      c = Vec3d(c[0], c[1], c[2]);

  // Check whether c is near negative triangle external boundary, and
  // rotate to corresponding positive triangle boundary if so.
  if (double_eq(c[0], c[1], eps))
    c = Vec3d(c[1], c[0], c[2]);
  if (double_eq(c[2], 0, eps))
    c = Vec3d(c[0], -c[2], c[1]);

  return c;
}

// sort an absolute value of Vec3d without altering original Vec3d
Vec3d sort_Vec3d(Vec3d &v)
{
  vector<double> c;
  c.push_back(fabs(v[0]));
  c.push_back(fabs(v[1]));
  c.push_back(fabs(v[2]));

  sort(c.begin(), c.end());

  return (Vec3d(c[0], c[1], c[2]));
}

void color_by_symmetry_normals(Geometry &geom, const char color_method,
                               const int face_opacity, const double eps)
{
  const vector<vector<int>> &faces = geom.faces();
  const vector<Vec3d> &verts = geom.verts();

  string map_string = "rnd";
  if (face_opacity > -1)
    map_string += msg_str("_A%g", (double)face_opacity / 255);
  // std::unique_ptr<ColorMap> cmap(colormap_from_name(map_name.c_str()));
  ColorMap *cmap = colormap_from_name(map_string.c_str());

  for (unsigned int i = 0; i < faces.size(); i++) {
    Vec3d norm = face_norm(verts, faces[i]).unit();
    if (strchr("yY", color_method))
      norm = sort_Vec3d(norm);
    else if (strchr("zZ", color_method))
      norm = sort_Vec3d_chiral(norm, eps);

    long idx = (long)(norm[0] * 1000000) + (long)(norm[1] * 10000) +
               (long)norm[2] * 100;
    if (strchr("YZ", color_method))
      geom.colors(FACES).set(i, cmap->get_col(idx));
    else
      geom.colors(FACES).set(i, idx);
  }
}

void color_edges_by_sqrt(Geometry &geom, const char color_method)
{
  geom.add_missing_impl_edges();

  string map_string = "rnd";
  // std::unique_ptr<ColorMap> cmap(colormap_from_name("rnd"));
  ColorMap *cmap = colormap_from_name(map_string.c_str());
  // e_coloring clrg(&geom);
  for (unsigned int i = 0; i < geom.edges().size(); i++) {
    // geom.colors(EDGES).set(i, int(floor(pow(geom.edge_len(i),2)+0.5)));
    int idx = int(floor(pow(geom.edge_len(i), 2) + 0.5));
    if (color_method == 'T')
      geom.colors(EDGES).set(i, cmap->get_col(idx));
    // geom.colors(EDGES).set(i,clrg.idx_to_rand_val(idx));
    else
      geom.colors(EDGES).set(i, idx);
  }
}
