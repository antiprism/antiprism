/*
   Copyright (c) 2003-2016, Adrian Rossiter

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
   Name: info.cc
   Description: information from OFF file
   Project: Antiprism - http://www.antiprism.com
*/

#include "geometryinfo.h"
#include "geometryutils.h"
#include "mathutils.h"
#include "private_misc.h"

#include <algorithm>
#include <cstring>
#include <limits>
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

namespace anti {

//---------------------------------------------------------------------
// Comparison

static int cmp_angles(const double &a, const double &b, double min_diff = 1e-6)
{
  if (fabs(a - b) > min_diff) {
    if (a < b)
      return -1;
    else
      return 1;
  }
  return 0;
}

static int cmp_face_angles(const vector<double> &f1, const vector<double> &f2,
                           double min_diff = 1e-6)
{
  if (f1.size() != f2.size())
    return 0;

  for (unsigned int i = 0; i < f1.size(); i++) {
    if (int cmp = cmp_angles(f1[i], f2[i], min_diff))
      return cmp;
  }
  return 0;
}

bool AngleLess::operator()(const double &ang1, const double &ang2) const
{
  return (cmp_angles(ang1, ang2) < 0);
}

bool AngleVectLess::operator()(const vector<double> &angs1,
                               const vector<double> &angs2) const
{
  if (angs1.size() != angs2.size())
    return (angs1.size() < angs2.size());
  return (cmp_face_angles(angs1, angs2) < 0);
}

//---------------------------------------------------------------------
// ElementLimits

ElementLimits::ElementLimits() { init(); }

void ElementLimits::init()
{
  idx[0] = -1;
  max = -1e100;
  min = 1e100;
  zero = 1e100;
  sum = 0;
}

bool ElementLimits::is_set() const { return idx[0] != -1; }

//---------------------------------------------------------------------
// GeometryInfo

GeometryInfo::GeometryInfo(const Geometry &geo, Vec3d center)
    : cent(center), geom(geo)
{
  reset();
}

void GeometryInfo::reset()
{
  oriented = -1;
  orientable = -1;
  found_connectivity = false;
  genus_val = std::numeric_limits<int>::max();
  dual.clear_all();
  sym = Symmetry();
  efpairs.clear();
  edge_parts.clear();
  face_angles.clear();
  vert_dihed.clear();
  dihedral_angles.clear();
  edge_index_numbers.clear();
  e_lengths.clear();
  sol_angles.clear();
  vertex_angles.clear();
  f_areas.clear();
  f_perimeters.clear();
  vert_impl_edges.clear();
  vert_faces.clear();
  vert_cons.clear();
  vert_cons_orig.clear();
  face_cons.clear();
  vert_norms.clear();
  found_free_verts = false;
  free_verts.clear();
  set_center(cent);
}

void GeometryInfo::set_center(Vec3d center)
{
  cent = center;
  v_dists.init();
  e_dists.init();
  ie_dists.init();
  f_dists.init();
}

const Geometry &GeometryInfo::get_geom() const { return geom; }

Vec3d GeometryInfo::get_center() const { return cent; }

void GeometryInfo::find_impl_edges() { geom.get_impl_edges(impl_edges); }

bool GeometryInfo::is_closed()
{
  if (!found_connectivity)
    find_connectivity();
  return closed;
}

bool GeometryInfo::is_polyhedron()
{
  if (!found_connectivity)
    find_connectivity();
  return polyhedron;
}

bool GeometryInfo::is_even_connectivity()
{
  if (!found_connectivity)
    find_connectivity();
  return even_connectivity;
}

bool GeometryInfo::is_known_connectivity()
{
  if (!found_connectivity)
    find_connectivity();
  return known_connectivity;
}

bool GeometryInfo::is_known_genus()
{
  return genus() < std::numeric_limits<int>::max() - 1;
}

// elements
int GeometryInfo::num_verts() const { return geom.verts().size(); }

int GeometryInfo::num_edges() const { return geom.edges().size(); }

int GeometryInfo::num_iedges() { return get_impl_edges().size(); }

int GeometryInfo::num_faces() const { return geom.faces().size(); }

int GeometryInfo::num_parts()
{
  is_orientable();
  return number_parts;
}

ElementLimits GeometryInfo::face_areas()
{
  if (f_areas.size() == 0)
    find_f_areas();
  return area;
}

double GeometryInfo::volume()
{
  if (f_areas.size() == 0)
    find_f_areas();
  return vol;
}

Vec3d GeometryInfo::volume_centroid()
{
  if (f_areas.size() == 0)
    find_f_areas();
  return vol_cent;
}

double GeometryInfo::isoperimetric_quotient()
{
  return 36.0 * M_PI * pow(volume(), 2) / pow(face_areas().sum, 3);
}

ElementLimits GeometryInfo::edge_length_lims()
{
  get_edge_lengths_by_size();
  return edge_len;
}

ElementLimits GeometryInfo::iedge_length_lims()
{
  get_iedge_lengths_by_size();
  return iedge_len;
}

ElementLimits GeometryInfo::dihed_angle_lims()
{
  if (dihedral_angles.size() == 0)
    find_dihedral_angles();
  return dih_angles;
}

ElementLimits GeometryInfo::solid_angle_lims()
{
  if (sol_angles.size() == 0)
    find_solid_angles();
  return so_angles;
}

ElementLimits &GeometryInfo::vert_dist_lims()
{
  if (!v_dists.is_set())
    find_v_dist_lims();
  return v_dists;
}

ElementLimits &GeometryInfo::edge_dist_lims()
{
  if (!e_dists.is_set())
    find_e_dist_lims();
  return e_dists;
}

ElementLimits &GeometryInfo::iedge_dist_lims()
{
  if (!ie_dists.is_set())
    find_ie_dist_lims();
  return ie_dists;
}

ElementLimits &GeometryInfo::face_dist_lims()
{
  if (!f_dists.is_set())
    find_f_dist_lims();
  return f_dists;
}

ElementLimits GeometryInfo::angle_lims()
{
  if (!plane_angles.size())
    find_face_angles();
  return ang;
}

int GeometryInfo::num_angles()
{
  if (!plane_angles.size())
    find_face_angles();
  return num_angs;
}

double GeometryInfo::angle_defect()
{
  return num_verts() * 2 * M_PI - angle_lims().sum;
}

// verts
const vector<Vec3d> &GeometryInfo::get_vert_norms(bool local_orient)
{
  if (!vert_norms.size() || local_orient != vert_norms_local_orient)
    find_vert_norms(local_orient);
  return vert_norms;
}

const vector<vector<int>> &GeometryInfo::get_vert_cons()
{
  if (!vert_cons.size())
    find_vert_cons();
  return vert_cons;
}

const vector<vector<int>> &GeometryInfo::get_vert_faces()
{
  if (!vert_faces.size())
    find_vert_elems(geom.faces(), vert_faces);
  return vert_faces;
}

const vector<vector<int>> &GeometryInfo::get_vert_impl_edges()
{
  if (!vert_impl_edges.size())
    find_vert_elems(get_impl_edges(), vert_impl_edges);
  return vert_impl_edges;
}

const vector<vector<vector<int>>> &GeometryInfo::get_face_cons()
{
  if (!face_cons.size())
    find_face_cons();
  return face_cons;
}

const vector<vector<vector<int>>> &GeometryInfo::get_vert_figs()
{
  if (!vert_figs.size())
    find_vert_figs();
  return vert_figs;
}

const vector<double> &GeometryInfo::get_vert_solid_angles()
{
  if (!vertex_angles.size())
    GeometryInfo::find_solid_angles();
  return vertex_angles;
}

const map<double, double_range_cnt, AngleLess> &
GeometryInfo::get_solid_angles_by_size()
{
  if (!sol_angles.size())
    find_solid_angles();
  return sol_angles;
}

const map<pair<int, int>, double> &GeometryInfo::get_plane_angles()
{
  if (!plane_angles.size())
    find_face_angles();
  return vf_plane_angles;
}

const vector<int> &GeometryInfo::get_free_verts()
{
  if (!found_free_verts)
    find_free_verts();
  return free_verts;
}

// edges

int GeometryInfo::get_edge_index(int v_idx0, int v_idx1)
{
  if (!edge_index_numbers.size())
    find_edge_index_numbers();
  auto e_it = edge_index_numbers.find(make_edge(v_idx0, v_idx1));
  return (e_it != edge_index_numbers.end()) ? e_it->second : -1;
}

const map<vector<int>, vector<int>> &GeometryInfo::get_edge_face_pairs()
{
  if (!efpairs.size())
    find_edge_face_pairs();
  return efpairs;
}

const vector<double> &GeometryInfo::get_edge_dihedrals()
{
  if (!dihedral_angles.size())
    find_dihedral_angles();
  return edge_dihedrals;
}

const vector<vector<int>> &GeometryInfo::get_edge_parts()
{
  if (!edge_parts.size())
    find_edge_parts();
  return edge_parts;
}

const map<double, double_range_cnt, AngleLess> &
GeometryInfo::get_dihedral_angles_by_size()
{
  if (!dihedral_angles.size())
    find_dihedral_angles();
  return dihedral_angles;
}

const map<double, double_range_cnt, AngleLess> &
GeometryInfo::get_edge_lengths_by_size()
{
  if (!e_lengths.size())
    find_e_lengths(e_lengths, geom.edges(), edge_len);
  return e_lengths;
}

// implicit edges
const vector<vector<int>> &GeometryInfo::get_impl_edges()
{
  if (!impl_edges.size())
    find_impl_edges();
  return impl_edges;
}

const map<double, double_range_cnt, AngleLess> &
GeometryInfo::get_iedge_lengths_by_size()
{
  if (!ie_lengths.size())
    find_e_lengths(ie_lengths, get_impl_edges(), iedge_len);
  return ie_lengths;
}

// faces
const map<vector<double>, int, AngleVectLess> &
GeometryInfo::get_plane_angles_by_size()
{
  if (!face_angles.size())
    find_face_angles();
  return face_angles;
}

const vector<double> &GeometryInfo::get_f_areas()
{
  if (!f_areas.size())
    find_f_areas();
  return f_areas;
}

const vector<double> &GeometryInfo::get_f_perimeters()
{
  if (!f_perimeters.size())
    find_f_perimeters();
  return f_perimeters;
}

const vector<double> &GeometryInfo::get_f_max_nonplanars()
{
  if (!f_max_nonplanars.size())
    find_f_max_nonplanars();
  return f_max_nonplanars;
}

const Geometry &GeometryInfo::get_dual()
{
  if (!dual.faces().size())
    anti::get_dual(dual, geom);
  return dual;
}

// Symmetry
const Symmetry &GeometryInfo::get_symmetry()
{
  if (!sym.is_set())
    find_symmetry();
  return sym;
}

string GeometryInfo::get_symmetry_type_name()
{
  if (!sym.is_set())
    find_symmetry();
  return sym.get_symbol();
}

const set<SymmetryAxis> &GeometryInfo::get_symmetry_axes()
{
  if (!sym.is_set())
    find_symmetry();
  return sym.get_axes();
}

const set<Symmetry> &GeometryInfo::get_symmetry_subgroups()
{
  if (!sym.is_set())
    find_symmetry();
  return sym.get_sub_syms();
}

const SymmetryAutos &GeometryInfo::get_symmetry_autos()
{
  if (!sym.is_set())
    find_symmetry();
  return sym.get_autos();
}

Trans3d GeometryInfo::get_symmetry_alignment_to_std()
{
  if (!sym.is_set())
    find_symmetry();
  return sym.get_to_std();
}

bool GeometryInfo::is_oriented()
{
  if (oriented < 0)
    oriented = geom.is_oriented();
  return oriented;
}

bool GeometryInfo::is_orientable()
{
  if (orientable < 0) {
    Geometry geom2 = geom;
    number_parts = geom2.orient();
    orientable = geom2.is_oriented();
  }
  return orientable;
}

void GeometryInfo::find_edge_face_pairs()
{
  if (is_oriented())
    efpairs = geom.get_edge_face_pairs(true);
  else
    efpairs = geom.get_edge_face_pairs(false);
}

void GeometryInfo::find_connectivity()
{
  // get a copy of edge face pairs with all faces around an edge
  map<vector<int>, vector<int>> tmp_efpairs;
  if (is_oriented())
    tmp_efpairs = geom.get_edge_face_pairs(false);
  else if (efpairs.size() == 0)
    find_edge_face_pairs();

  const map<vector<int>, vector<int>> &pairs =
      is_oriented() ? tmp_efpairs : efpairs;

  known_connectivity = true;
  even_connectivity = true;
  polyhedron = true;
  closed = true;
  map<vector<int>, vector<int>>::const_iterator ei;
  for (ei = pairs.begin(); ei != pairs.end(); ++ei) {
    if (ei->second.size() == 1) // One faces at an edge
      closed = false;
    if (ei->second.size() != 2) // Edge not met be exactly 2 faces
      polyhedron = false;
    if (ei->second.size() % 2) // Odd number of faces at an edge
      even_connectivity = false;
    if (ei->second.size() > 2) // More than two faces at an edge
      known_connectivity = false;
  }

  found_connectivity = true;
}

static double face_vol(const Geometry &geom, int f_no, Vec3d *face_vol_cent)
{
  double f_vol = 0;
  Vec3d f_vol_cent = Vec3d(0, 0, 0);
  const vector<Vec3d> &verts = geom.verts();
  const vector<int> &face = geom.faces(f_no);
  Vec3d V = verts[0];
  Vec3d v0 = verts[face[0]];
  for (unsigned int i = 1; i < face.size() - 1; i++) {
    const Vec3d &v1 = verts[face[i]];
    const Vec3d &v2 = verts[face[i + 1]];
    double tet_vol = vtriple(v0 - V, v1 - V, v2 - V);
    f_vol_cent += tet_vol * (V + v0 + v1 + v2); // /4 deferred
    f_vol += tet_vol;                           // /6 deferred
  }
  if (!double_eq(f_vol, 0))
    *face_vol_cent = (f_vol_cent / 4) / f_vol;
  else
    *face_vol_cent = Vec3d(0, 0, 0);

  return f_vol / 6;
}

void GeometryInfo::find_f_areas()
{
  int fsz = geom.faces().size();
  f_areas.resize(fsz);
  area.init();
  vol = 0;
  vol_cent = Vec3d(0, 0, 0);
  for (int i = 0; i < fsz; i++) {
    f_areas[i] = geom.face_norm(i, true).len();
    if (f_areas[i] < area.min) {
      area.min = f_areas[i];
      area.idx[ElementLimits::IDX_MIN] = i;
    }
    if (f_areas[i] > area.max) {
      area.max = f_areas[i];
      area.idx[ElementLimits::IDX_MAX] = i;
    }
    area.sum += f_areas[i];
    Vec3d f_vol_cent;
    double f_vol = face_vol(geom, i, &f_vol_cent);
    vol += f_vol;
    vol_cent += f_vol_cent * f_vol;
  }
  if (!double_eq(vol, 0))
    vol_cent /= vol;
  else
    vol_cent.unset();
}

void GeometryInfo::find_f_perimeters()
{
  f_perimeters.resize(num_faces());
  for (unsigned int i = 0; i < geom.faces().size(); i++) {
    double perim = 0.0;
    for (unsigned int j = 0; j < geom.faces(i).size(); j++)
      perim += geom.edge_vec(geom.faces(i, j), geom.faces_mod(i, j + 1)).len();
    f_perimeters[i] = perim;
  }
}

void GeometryInfo::find_face_angles()
{
  ang.init();
  num_angs = 0;
  map<vector<double>, int, AngleVectLess>::iterator fi;
  map<double, double_range_cnt, AngleLess>::iterator ai;

  for (unsigned int i = 0; i < geom.faces().size(); i++) {
    unsigned int fsz = geom.faces(i).size();
    vector<double> f1(fsz), f2(fsz), f_min(fsz);
    geom.face_angles_lengths(i, &f1);
    f_min = f1;
    for (unsigned int offset = 0; offset < fsz; offset++) {
      pair<int, int> vf_pr(geom.faces(i, offset), i);
      const double f_ang = f1[offset];
      vf_plane_angles[vf_pr] = f_ang;
      ai = plane_angles.find(f_ang);
      if (ai == plane_angles.end())
        plane_angles[f_ang] = double_range_cnt().update(f_ang);
      else
        ai->second.update(f_ang);
      if (ang.max < f_ang)
        ang.max = f_ang;
      if (ang.min > f_ang)
        ang.min = f_ang;
      ang.sum += f_ang;
      num_angs++;

      for (unsigned int k = 0; k < fsz; k++)
        f2[(k + offset) % fsz] = f1[k];
      if (cmp_face_angles(f2, f_min) < 0)
        f_min = f2;
      reverse(f2.begin(), f2.end());
      if (cmp_face_angles(f2, f_min) < 0)
        f_min = f2;
    }
    fi = face_angles.find(f_min);
    if (fi == face_angles.end())
      face_angles[f_min] = 1;
    else
      fi->second += 1;
  }
}

void GeometryInfo::find_dihedral_angles()
{
  if (efpairs.size() == 0)
    find_edge_face_pairs();
  edge_dihedrals.resize(efpairs.size());

  dih_angles.init();
  map<vector<int>, vector<int>>::iterator ei;
  map<double, double_range_cnt, AngleLess>::iterator di;
  double cos_a = 1, sign = 1;
  int e_idx = -1;
  for (ei = efpairs.begin(); ei != efpairs.end(); ++ei) {
    e_idx++;
    if (ei->second.size() == 2 && ei->second[0] >= 0 &&
        ei->second[1] >= 0) { // pair of faces
      Vec3d n0;
      Vec3d n1;
      if (is_oriented()) {
        n0 = geom.face_norm(ei->second[0]).unit();
        n1 = geom.face_norm(ei->second[1]).unit();
        Vec3d e_dir = geom.verts(ei->first[1]) - geom.verts(ei->first[0]);
        sign = vdot(e_dir, vcross(n0, n1));
      }
      else {
        vector<int> f0 = geom.faces(ei->second[0]);
        vector<int> f1 = geom.faces(ei->second[1]);
        orient_face(f0, ei->first[0], ei->first[1]);
        orient_face(f1, ei->first[1], ei->first[0]);
        n0 = face_norm(geom.verts(), f0).unit();
        n1 = face_norm(geom.verts(), f1).unit();
        sign = 1;
      }
      cos_a = -vdot(n0, n1);
    }
    else
      cos_a = 1; // one face on an edge
    double ang = acos(safe_for_trig(cos_a));
    if (sign < 0) // in oriented polyhedron
      ang = 2 * M_PI - ang;

    edge_dihedrals[e_idx] = ang;

    if (ang > dih_angles.max) {
      dih_angles.max = ang;
      dih_angles.idx[ElementLimits::IDX_MAX] = ei->first[0];
      dih_angles.idx[ElementLimits::IDX_MAX2] = ei->first[1];
    }
    if (ang < dih_angles.min) {
      dih_angles.min = ang;
      dih_angles.idx[ElementLimits::IDX_MIN] = ei->first[0];
      dih_angles.idx[ElementLimits::IDX_MIN2] = ei->first[1];
    }
    if (fabs(ang - M_PI) < fabs(dih_angles.zero)) {
      dih_angles.zero = ang;
      dih_angles.idx[ElementLimits::IDX_ZERO] = ei->first[0];
      dih_angles.idx[ElementLimits::IDX_ZERO2] = ei->first[1];
    }

    di = dihedral_angles.find(ang);
    if (di == dihedral_angles.end())
      dihedral_angles[ang] = double_range_cnt().update(ang);
    else
      di->second.update(ang);
  }
}

void GeometryInfo::find_vert_cons_orig()
{
  vert_cons_orig.resize(num_verts(), vector<int>());
  if (!num_faces())
    return;
  Geometry dual = get_dual();

  vector<int> del_faces;
  const vector<vector<int>> &dfaces = dual.faces();

  // this needs checking
  for (unsigned int i = 0; i < dfaces.size(); i++)
    if (dfaces[i].size() < 2)
      del_faces.push_back(i);
  // dual.delete_faces(del_faces);

  if (del_faces.size()) { // set up an unordered list
    vector<vector<int>> es = geom.edges();
    geom.get_impl_edges(es);
    for (auto &e : es) {
      vert_cons_orig[e[0]].push_back(e[1]);
      vert_cons_orig[e[1]].push_back(e[0]);
    }
    return;
  }

  dual.orient();
  const vector<vector<int>> &faces = geom.faces();
  for (unsigned int i = 0; i < dfaces.size(); i++) {
    for (unsigned int j = 0; j < dfaces[i].size(); j++) {
      int f = dfaces[i][j];
      int sz = faces[f].size();
      for (int v = 0; v < sz; v++) {
        if (faces[f][v] == (int)i) { // found vertex in surrounding face
          if (j == 0) {              // first vertex, must be in following face
            int idx = faces[f][(v + 1) % sz];
            const vector<int> &nface =
                faces[dfaces[i][(j + 1) % dfaces[i].size()]];
            if (find(nface.begin(), nface.end(), idx) != nface.end())
              vert_cons_orig[i].push_back(faces[f][(v + 1) % sz]);
            else
              vert_cons_orig[i].push_back(faces[f][(v + sz - 1) % sz]);
          }
          else {
            if (faces[f][(v + 1) % sz] == vert_cons_orig[i][j - 1])
              vert_cons_orig[i].push_back(faces[f][(v + sz - 1) % sz]);
            else
              vert_cons_orig[i].push_back(faces[f][(v + 1) % sz]);
          }
          break;
        }
      }
    }
  }

  // first dual vertex on dual face 0 is a face containing vertex 0
  int f0 = dfaces[0][0];
  // find vertices before and after 0
  vector<int> edge(2);
  int sz = faces[f0].size();
  for (int i = 0; i < sz; i++)
    if (faces[f0][i] == 0) {
      edge[0] = faces[f0][(i - 1 + sz) % sz];
      edge[1] = faces[f0][(i + 1) % sz];
      break;
    }

  for (unsigned int j = 0; j < vert_cons_orig[0].size(); j++) {
    if (vert_cons_orig[0][j] == edge[0]) {
      if (vert_cons_orig[0][(j + 1) % vert_cons_orig[0].size()] != edge[1]) {
        for (auto &v : vert_cons_orig)
          reverse(v.begin(), v.end());
      }
      break;
    }
  }
}

void GeometryInfo::find_vert_cons()
{
  vert_cons.resize(num_verts(), vector<int>());
  vector<vector<int>> es = geom.edges();
  geom.get_impl_edges(es);
  for (auto &e : es) {
    vert_cons[e[0]].push_back(e[1]);
    vert_cons[e[1]].push_back(e[0]);
  }
}

void GeometryInfo::find_vert_elems(const vector<vector<int>> &elems,
                                   vector<vector<int>> &vert_elems)
{
  vert_elems.resize(num_verts(), vector<int>());
  for (size_t elem_idx = 0; elem_idx < elems.size(); elem_idx++)
    for (int v_idx : elems[elem_idx])
      vert_elems[v_idx].push_back(elem_idx);

  for (auto v_idxs : vert_elems) {
    sort(v_idxs.begin(), v_idxs.end());
    v_idxs.erase(unique(v_idxs.begin(), v_idxs.end()), v_idxs.end());
  }
}

void GeometryInfo::find_face_cons()
{
  face_cons.resize(num_faces(), vector<vector<int>>());
  auto e2fs = geom.get_edge_face_pairs(false);
  for (unsigned int f_idx = 0; f_idx < geom.faces().size(); f_idx++) {
    for (unsigned int v = 0; v < geom.faces(f_idx).size(); v++) {
      face_cons[f_idx].push_back(vector<int>());
      auto ei = e2fs.find(
          make_edge(geom.faces(f_idx, v), geom.faces_mod(f_idx, v + 1)));
      if (ei != e2fs.end()) { // test, but shouldn't fail
        for (int i : ei->second) {
          if (i != (int)f_idx)
            face_cons[f_idx][v].push_back(i);
        }
      }
    }
  }
}

void GeometryInfo::find_vert_figs()
{
  vert_figs.resize(num_verts());
  get_vert_cons();
  auto ef_pairs = geom.get_edge_face_pairs(false);

  // find set of faces that each vertex belongs to
  const int v_sz = geom.verts().size();
  vector<set<int>> v_faces(v_sz);
  for (unsigned int f_idx = 0; f_idx < geom.faces().size(); f_idx++)
    for (unsigned int v = 0; v < geom.faces(f_idx).size(); v++)
      v_faces[geom.faces(f_idx, v)].insert(f_idx);

  // copy of vertices to be used for creating the sets of triangles
  Geometry g_fig;
  g_fig.raw_verts() = geom.verts();
  map<vector<int>, int> circuit_edge_cnts;

  for (int i = 0; i < v_sz; i++) {
    g_fig.clear(FACES);
    circuit_edge_cnts.clear();
    bool figure_good = true;
    for (auto si = v_faces[i].begin(); si != v_faces[i].end(); ++si) {
      const int f = *si;
      vector<int> tri(3);
      for (unsigned int n = 0; n < geom.faces(f).size(); n++) {
        if (geom.faces(f, n) == i) {
          tri[0] = geom.faces_mod(f, n - 1);
          tri[1] = geom.faces(f, n);
          tri[2] = geom.faces_mod(f, n + 1);
          if (ef_pairs[make_edge(tri[0], tri[1])].size() != 2 ||
              ef_pairs[make_edge(tri[1], tri[2])].size() != 2) {
            figure_good = false;
            break; // finish processing this face from set
          }
          circuit_edge_cnts[make_edge(tri[2], tri[0])]++;
          g_fig.add_face(tri);
        }
      }
      if (!figure_good)
        break; // finish processing this vertex
    }
    if (figure_good) {
      unsigned int num_tris = g_fig.faces().size();
      close_poly_basic(g_fig);
      for (unsigned int f_idx = num_tris; f_idx < g_fig.faces().size(); f_idx++)
        vert_figs[i].push_back(g_fig.faces(f_idx));
      map<vector<int>, int>::const_iterator ei;
      for (ei = circuit_edge_cnts.begin(); ei != circuit_edge_cnts.end(); ++ei)
        if (ei->second == 2) // closed edge must be digonal figure
          vert_figs[i].push_back(ei->first);
    }
  }
}

// Calculate vertex normals as average of surrounding face normals,
// using existing face orientation (an optimisation when it is known
// that a model is oriented)
static void get_vert_norms_raw_orientation(const Geometry &geom,
                                           vector<Vec3d> &v_norms)
{
  const unsigned int verts_sz = geom.verts().size();
  v_norms.assign(verts_sz, Vec3d::zero);
  vector<int> vert_face_cnt(geom.verts().size(), 0);
  for (unsigned int i = 0; i < geom.faces().size(); i++) {
    Vec3d norm = geom.face_norm(i).unit();
    for (unsigned int j = 0; j < geom.faces(i).size(); j++) {
      const int v_idx = geom.faces(i, j);
      v_norms[v_idx] += norm;
      vert_face_cnt[v_idx]++;
    }
  }

  for (unsigned int i = 0; i < verts_sz; i++) {
    if (int cnt = vert_face_cnt[i])
      v_norms[i] /= cnt;
    else
      v_norms[i].unset();
  }
}

// Calculate vertex normals as average of surrounding face normals,
// but locally orienting faces consistenly around a vertex before calculating
// the normals to average
static void get_vert_norms_general(const Geometry &geom, vector<Vec3d> &v_norms)
{

  // find set of faces that each vertex belongs to
  const int v_sz = geom.verts().size();
  vector<set<int>> v_faces(v_sz);
  for (unsigned int f_idx = 0; f_idx < geom.faces().size(); f_idx++)
    for (unsigned int v = 0; v < geom.faces(f_idx).size(); v++)
      v_faces[geom.faces(f_idx, v)].insert(f_idx);

  // copy of vertices to be used for orientating the sets of faces
  Geometry g_orient;
  g_orient.raw_verts() = geom.verts();

  v_norms.assign(v_sz, Vec3d::zero); // initialise normals to zero
  for (int i = 0; i < v_sz; i++) {
    for (int si : v_faces[i])
      g_orient.add_face(geom.faces(si));
    g_orient.orient();
    vector<Vec3d> f_norms;
    g_orient.face_norms(f_norms);
    for (auto &f_norm : f_norms)
      v_norms[i] += f_norm;
    v_norms[i] /= f_norms.size();
    g_orient.clear(FACES);
  }
}

void GeometryInfo::find_vert_norms(bool local_orient)
{
  vert_norms_local_orient = local_orient;
  if (local_orient)
    get_vert_norms_general(get_geom(), vert_norms);
  else
    get_vert_norms_raw_orientation(get_geom(), vert_norms);
}

void GeometryInfo::find_free_verts()
{
  free_verts.clear();
  const Geometry &geom = get_geom();
  vector<int> cnt(geom.verts().size());
  for (unsigned int i = 0; i < geom.faces().size(); i++)
    for (unsigned int j = 0; j < geom.faces(i).size(); j++)
      cnt[geom.faces(i, j)]++;
  for (unsigned int i = 0; i < geom.edges().size(); i++)
    for (unsigned int j = 0; j < geom.edges(i).size(); j++)
      cnt[geom.edges(i, j)]++;

  for (int i = 0; i < (int)cnt.size(); i++)
    if (cnt[i] == 0)
      free_verts.push_back(i);

  found_free_verts = true;
}

static double sph_tri_area(Vec3d u0, Vec3d u1, Vec3d u2)
{
  double sign = 1 - 2 * (vtriple(u0, u1, u2) > 0);
  Vec3d u[3] = {u0, u1, u2};
  double ang, area = 0;
  for (int i = 0; i < 3; i++) {
    ang = acos(safe_for_trig(vdot(vcross(u[i], u[(i - 1 + 3) % 3]).unit(),
                                  vcross(u[i], u[(i + 1) % 3]).unit())));
    if (sign < 0)
      ang = 2 * M_PI - ang;
    area += ang;
  }
  return area - M_PI;
}

void GeometryInfo::find_solid_angles()
{
  if (!vert_cons_orig.size())
    find_vert_cons_orig();

  vertex_angles = vector<double>(num_verts(), 0);

  so_angles.init();
  map<double, double_range_cnt, AngleLess>::iterator si;
  for (unsigned int i = 0; i < vert_cons_orig.size(); i++) {
    vector<Vec3d> dirs(vert_cons_orig[i].size());
    for (unsigned int j = 0; j < vert_cons_orig[i].size(); j++)
      dirs[j] = geom.verts(i) - geom.verts(vert_cons_orig[i][j]);

    for (unsigned int j = 1; j < vert_cons_orig[i].size() - 1; j++)
      vertex_angles[i] += sph_tri_area(dirs[0], dirs[j], dirs[j + 1]);

    // if(!is_oriented()) {
    //   fmod(vertex_angles[i], 4*M_PI);
    //   if(vertex_angles[i]>2*M_PI)
    //      vertex_angles[i] = 4*M_PI - vertex_angles[i];
    // }

    if (vertex_angles[i] > so_angles.max) {
      so_angles.max = vertex_angles[i];
      so_angles.idx[ElementLimits::IDX_MAX] = i;
    }
    if (vertex_angles[i] < so_angles.min) {
      so_angles.min = vertex_angles[i];
      so_angles.idx[ElementLimits::IDX_MIN] = i;
    }
    if (fabs(vertex_angles[i]) < fabs(so_angles.zero)) {
      so_angles.zero = vertex_angles[i];
      so_angles.idx[ElementLimits::IDX_ZERO] = i;
    }

    double s_ang = vertex_angles[i];
    si = sol_angles.find(s_ang);
    if (si == sol_angles.end())
      sol_angles[s_ang] = double_range_cnt().update(s_ang);
    else
      si->second.update(s_ang);
  }
}

void GeometryInfo::find_edge_index_numbers()
{
  edge_index_numbers.clear();
  for (int i = 0; i < (int)geom.edges().size(); i++)
    edge_index_numbers[geom.edges(i)] = i;
}

void GeometryInfo::find_e_lengths(
    map<double, double_range_cnt, AngleLess> &e_lens,
    const vector<vector<int>> &edges, ElementLimits &lens)
{
  lens.init();
  map<double, double_range_cnt, AngleLess>::iterator ei;
  for (const auto &edge : edges) {
    double dist = geom.edge_len(edge);
    if (dist > lens.max) {
      lens.max = dist;
      lens.idx[ElementLimits::IDX_MAX] = edge[0];
      lens.idx[ElementLimits::IDX_MAX2] = edge[1];
    }
    if (dist < lens.min) {
      lens.min = dist;
      lens.idx[ElementLimits::IDX_MIN] = edge[0];
      lens.idx[ElementLimits::IDX_MIN2] = edge[1];
    }
    lens.sum += dist;
    ei = e_lens.find(dist);
    if (ei == e_lens.end())
      e_lens[dist] = double_range_cnt().update(dist);
    else
      ei->second.update(dist);
  }
}

// get indexes of edge parts
static void get_edge_part(vector<int> &edge_part, const Geometry &geom,
                          const int idx, const vector<vector<int>> &vcons,
                          vector<bool> &seen)
{
  if (seen[idx])
    return;
  else
    seen[idx] = true;

  for (unsigned int i = 0; i < vcons[idx].size(); i++) {
    int next_idx = vcons[idx][i];
    if (idx == next_idx)
      continue;
    vector<int> edge(2);
    edge[0] = idx;
    edge[1] = next_idx;
    edge_part.push_back(find_edge_in_edge_list(geom.edges(), edge));
    get_edge_part(edge_part, geom, next_idx, vcons, seen);
  }
}

// get list of edge parts indexes
void GeometryInfo::find_edge_parts()
{
  vector<vector<int>> vcons(geom.verts().size(), vector<int>());

  const vector<vector<int>> &edges = geom.edges();
  for (const auto &edge : edges) {
    vcons[edge[0]].push_back(edge[1]);
    vcons[edge[1]].push_back(edge[0]);
  }

  vector<bool> seen(vcons.size(), false);
  for (unsigned int i = 0; i < vcons.size(); i++)
    if (!seen[i]) {
      vector<int> edge_part;

      // find the edge parts by index list
      get_edge_part(edge_part, geom, i, vcons, seen);

      // if no edge part is produced it will be empty
      if (edge_part.size()) {
        // the lists may contain duplicates and must be unique
        sort(edge_part.begin(), edge_part.end());
        auto ep = unique(edge_part.begin(), edge_part.end());
        edge_part.resize(ep - edge_part.begin());

        edge_parts.push_back(edge_part);
      }
    }
}

void GeometryInfo::find_f_max_nonplanars()
{
  f_max_nonplanars.resize(geom.faces().size());
  for (unsigned int f = 0; f < geom.faces().size(); f++) {
    if (geom.faces(f).size() == 3) {
      f_max_nonplanars[f] = 0;
      continue;
    }
    Vec3d norm = geom.face_norm(f).unit();
    Vec3d f_cent = geom.face_cent(f);
    double max = 0;
    for (unsigned int v = 0; v < geom.faces(f).size(); v++) {
      double dist = fabs(vdot(norm, f_cent - geom.verts(geom.faces(f, v))));
      if (dist > max)
        max = dist;
    }
    f_max_nonplanars[f] = max;
  }
}

void GeometryInfo::find_symmetry() { sym.init(geom); }

int GeometryInfo::genus()
{
  if (genus_val == std::numeric_limits<int>::max()) {
    genus_val = std::numeric_limits<int>::max() - 1; // 'not known' value
    if (num_parts() == 1 && is_known_connectivity()) {
      int euler_char = num_verts() - num_iedges() + num_faces();
      Geometry geom2 = geom;
      if (close_poly_basic(geom2)) {
        int num_boundaries = geom2.faces().size() - num_faces();
        genus_val = 2 - num_boundaries - euler_char;
        if (is_orientable())
          genus_val /= 2;
        else
          genus_val *= -1; // to indicate non-orientable genus (demigenus)
      }
    }
  }
  return genus_val;
}

//--------------------------------------------------------------------
// Distances

static void set_edge_dists(const Geometry &geom, Vec3d cent,
                           const vector<vector<int>> &edges, ElementLimits &lim)
{
  lim.init();
  for (const auto &edge : edges) {
    double dist = (geom.edge_nearpt(edge, cent) - cent).len();
    if (dist < lim.min) {
      lim.min = dist;
      lim.idx[ElementLimits::IDX_MIN] = edge[0];
      lim.idx[ElementLimits::IDX_MIN2] = edge[1];
    }
    if (dist > lim.max) {
      lim.max = dist;
      lim.idx[ElementLimits::IDX_MAX] = edge[0];
      lim.idx[ElementLimits::IDX_MAX2] = edge[1];
    }
    lim.sum += dist;
  }
}

void GeometryInfo::find_v_dist_lims()
{
  ElementLimits &lim = v_dists;
  lim.init();
  for (unsigned int i = 0; i < geom.verts().size(); i++) {
    double dist = (geom.verts(i) - cent).len();
    if (dist < lim.min) {
      lim.min = dist;
      lim.idx[ElementLimits::IDX_MIN] = i;
    }
    if (dist > lim.max) {
      lim.max = dist;
      lim.idx[ElementLimits::IDX_MAX] = i;
    }
    lim.sum += dist;
  }
}

void GeometryInfo::find_e_dist_lims()
{
  set_edge_dists(geom, cent, geom.edges(), e_dists);
}

void GeometryInfo::find_ie_dist_lims()
{
  vector<vector<int>> edges;
  geom.get_impl_edges(edges);
  set_edge_dists(geom, cent, edges, ie_dists);
}

void GeometryInfo::find_f_dist_lims()
{
  ElementLimits &lim = f_dists;
  lim.init();
  for (unsigned int i = 0; i < geom.faces().size(); i++) {
    double dist = geom.face_nearpt(i, cent).len();
    if (dist < lim.min) {
      lim.min = dist;
      lim.idx[ElementLimits::IDX_MIN] = i;
    }
    if (dist > lim.max) {
      lim.max = dist;
      lim.idx[ElementLimits::IDX_MAX] = i;
    }
    lim.sum += dist;
  }
}

} // namespace anti
