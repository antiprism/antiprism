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
   Name: dual.cc
   Description: creates duals
   Project: Antiprism - http://www.antiprism.com
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "coloring.h"
#include "geometryinfo.h"
#include "mathutils.h"
#include "symmetry.h"
#include "utils.h"

using std::vector;
using std::map;
using std::pair;
using std::swap;
using std::string;

namespace anti {

void sym_repeat(Geometry *geom, const Geometry &part, const Transformations &ts,
                char col_part_elems, Coloring *clrngs)
{
  Coloring tmp_clrngs[3];
  if (!clrngs)
    clrngs = tmp_clrngs;

  Geometry symmetry_unit = part;
  geom->clear_all();
  Transformations::const_iterator si;
  int idx = 0;
  for (si = ts.begin(); si != ts.end(); si++, idx++) {
    Geometry sym_unit = symmetry_unit;
    sym_unit.transform(*si);
    for (int i = 0; i < 3; i++)
      clrngs[i].set_geom(&sym_unit);

    if (col_part_elems & ELEM_VERTS)
      clrngs[VERTS].v_one_col(clrngs[0].get_col(idx));
    if (col_part_elems & ELEM_EDGES) {
      sym_unit.add_missing_impl_edges();
      clrngs[1].e_one_col(clrngs[EDGES].get_col(idx));
    }
    if (col_part_elems & ELEM_FACES)
      clrngs[2].f_one_col(clrngs[FACES].get_col(idx));
    geom->append(sym_unit);
  }
}

bool sym_repeat(Geometry *geom, const Geometry &part, const Symmetry &sym,
                char col_part_elems, Coloring *clrngs)
{
  Transformations ts;
  sym.get_trans(ts);
  if (!ts.is_set())
    return false;
  sym_repeat(geom, part, ts, col_part_elems, clrngs);
  return true;
}

void orient_face(vector<int> &face, int v0, int v1)
{
  for (unsigned int i = 0; i < face.size(); i++) {
    if (face[i] == v0) {
      if (face[(i + 1) % face.size()] != v1)
        reverse(face.begin(), face.end());
      break;
    }
  }
}

int orient_geom(Geometry &geom, vector<vector<int>> *parts)
{
  int part_num = 0;
  const int done = -1;
  auto edges = geom.get_edge_face_pairs(false);
  vector<int> cur_idx(geom.faces().size(), 0);
  vector<int> prev_face(geom.faces().size(), 0);
  vector<int> e_verts(2);
  vector<int> orig_e_verts(2);
  vector<int> e_faces(2);
  for (unsigned int i = 0; i < geom.faces().size(); i++) {
    if (geom.faces(i).size() < 3)
      cur_idx[i] = done; // don't process degenerate faces
    if (cur_idx[i] == done)
      continue;
    int cur_fidx = i;
    if (parts)
      parts->push_back(vector<int>(1, i));
    while (cur_idx[i] != done) {
      int idx = cur_idx[cur_fidx];
      if (idx == done) {
        cur_fidx = prev_face[cur_fidx];
        continue;
      }

      // read off the next edge
      const vector<int> &face = geom.faces(cur_fidx);
      orig_e_verts[0] = face[idx];
      idx = (idx + 1) % face.size();
      orig_e_verts[1] = face[idx];
      cur_idx[cur_fidx] = idx ? idx : done; // set to next idx, or mark done

      e_verts = orig_e_verts;
      if (e_verts[0] > e_verts[1])
        swap(e_verts[0], e_verts[1]);
      e_faces = edges.find(e_verts)->second;

      int next_face = (e_faces[0] != cur_fidx) ? e_faces[0] : e_faces[1];
      if (next_face >= 0 && cur_idx[next_face] == 0) { // face not looked at yet
        orient_face(geom.raw_faces()[next_face], orig_e_verts[1],
                    orig_e_verts[0]);
        if (parts)
          parts->back().push_back(next_face);
        prev_face[next_face] = cur_fidx;
        cur_fidx = next_face;
      }
    }
    part_num++;
  }
  return part_num;
}

// From planar.cc (Roger Kaufman)
// volume to positive=1, negative=2, reverse=3
// or flip=4 which reverse the orientation of the input model\n"
bool orient_geom(Geometry &geom, int option, char *errmsg)
{
  if (errmsg)
    *errmsg = '\0';

  GeometryInfo info(geom);
  bool is_orientable = info.is_orientable();
  if (!is_orientable)
    if (errmsg)
      strcpy_msg(errmsg, "input file contains a non-orientable geometry, "
                         "use 'flip' to reverse face orientations");
  // if model is not oriented, don't do a pre-orientation if we just want
  // orient_reverse
  if (!info.is_oriented() && option != 4)
    geom.orient();
  info.reset();
  double vol = info.volume();
  if (is_orientable && vol == 0 && (option == 1 || option == 2))
    if (errmsg)
      strcpy_msg(errmsg, "volume is zero, use 'reverse' for alternative "
                         "orientation");
  if ((vol < 0 && option == 1) || (vol > 0 && option == 2) || option == 3 ||
      option == 4)
    geom.orient_reverse();

  return option == 4 || is_orientable;
}

void orient_reverse(Geometry &geom)
{
  for (unsigned int i = 0; i < geom.faces().size(); ++i)
    reverse(geom.raw_faces()[i].begin(), geom.raw_faces()[i].end());
}

void get_pol_recip_verts(Geometry *dual, const Geometry &geom, double recip_rad,
                         Vec3d centre, double inf)
{
  const double min_lim = 1e-15;
  // const int r_sign = 1 - 2*(recip_rad<0); // -ve rad will reflect in centre
  const int r_sign = 1; // no reflection in centre

  Vec3d vert;
  const vector<Vec3d> &verts = geom.verts();
  const vector<vector<int>> &faces = geom.faces();
  map<pair<int, int>, pair<int, int>> edges;
  map<pair<int, int>, pair<int, int>>::iterator mi;
  pair<int, int> edge;
  dual->clear(VERTS);
  for (const auto &face : faces) {
    if (recip_rad == 0) { // all dual vertices = centre
      dual->add_vert(centre);
      continue;
    }

    for (unsigned int v1 = 0; v1 + 2 < face.size(); ++v1) {
      for (unsigned int v2 = v1 + 1; v2 + 1 < face.size(); ++v2) {
        for (unsigned int v3 = v2 + 1; v3 < face.size(); ++v3) {
          vert = verts[face[v2]] - verts[face[v1]];
          vert = vcross(vert, verts[face[v3]] - verts[face[v2]]);
          if (vert.len() > min_lim) {
            double face_dist = vdot(vert, verts[face[0]] - centre);
            double dist = r_sign * recip_rad * recip_rad / face_dist;
            if (fabs(face_dist) < min_lim || fabs(dist) > inf)
              dist = r_sign * inf * (1 - 2 * (face_dist < 0)) / vert.len();
            dual->add_vert(vert * dist + centre);
            v1 = v2 = v3 = face.size(); // to move on to next face
            break;
          }
        }
      }
    }
  }
}

class zero_sz {
public:
  bool operator()(const vector<int> &f) const { return f.size() == 0; }
};

void get_dual(Geometry *dual, const Geometry &geom, double recip_rad,
              Vec3d centre, double inf)
{
  get_pol_recip_verts(dual, geom, recip_rad, centre, inf);
  vector<vector<int>> d_faces(geom.verts().size());

  const vector<vector<int>> &faces = geom.faces();
  map<pair<int, int>, pair<int, int>> edges;
  map<pair<int, int>, pair<int, int>>::iterator mi;
  pair<int, int> edge;
  for (unsigned int i = 0; i < faces.size(); ++i) {
    for (unsigned int j = 0; j < faces[i].size(); ++j) {
      edge.first = faces[i][j];
      edge.second = faces[i][(j + 1) % faces[i].size()];
      bool swap_idxs = (edge.first > edge.second);
      if (swap_idxs)
        swap(edge.first, edge.second);
      mi = edges.find(edge);
      if (mi != edges.end()) {
        mi->second.second = i;
        if (swap_idxs)
          swap(mi->second.first, mi->second.second);
      }
      else
        edges[edge].first = i;
    }
  }

  for (mi = edges.begin(); mi != edges.end(); mi++) {
    d_faces[mi->first.first].push_back(mi->second.first);
    d_faces[mi->first.first].push_back(mi->second.second);
    d_faces[mi->first.second].push_back(mi->second.second);
    d_faces[mi->first.second].push_back(mi->second.first);
  }

  vector<int>::iterator vi;
  for (auto &d_face : d_faces) {
    for (unsigned int j = 0; j + 2 < d_face.size(); j += 2) {
      vi = find(d_face.begin() + j + 3, d_face.end(), d_face[j + 1]);
      if (vi == d_face.end())
        continue;

      if (vi - (d_face.begin() + j) != 3)
        swap(*vi, d_face[j + 2]);
      if (is_even(vi - d_face.begin()))
        swap(*(vi + 1), d_face[j + 3]);
      else
        swap(*(vi - 1), d_face[j + 3]);
    }
    vi = unique(d_face.begin(), d_face.end());
    if (vi != d_face.begin())
      d_face.erase(vi - 1, d_face.end());
  }

  dual->clear(EDGES);
  const vector<vector<int>> &g_edges = geom.edges();
  dual->colors(FACES) = geom.colors(VERTS);
  dual->colors(VERTS) = geom.colors(FACES);
  vector<int> g_edge(2);
  vector<int> d_edge(2);
  for (mi = edges.begin(); mi != edges.end(); mi++) {
    g_edge[0] = mi->first.first;
    g_edge[1] = mi->first.second;
    d_edge[0] = mi->second.first;
    d_edge[1] = mi->second.second;
    int gidx, didx;
    vector<vector<int>>::const_iterator ei;
    ei = find(g_edges.begin(), g_edges.end(), g_edge);
    if (ei != g_edges.end()) {
      gidx = ei - g_edges.begin();
      didx = dual->add_edge(d_edge);
      dual->colors(EDGES).set(didx, geom.colors(EDGES).get(gidx));
    }
  }

  for (unsigned int i = 0; i < d_faces.size(); ++i)
    if (d_faces[i].size() >= 3)
      dual->add_face(d_faces[i]);
}

void add_extra_ideal_elems(Geometry *geom, Vec3d centre, double inf)
{
  map<int, int> ideals;
  double inf2 = inf * inf;
  int sz = geom->verts().size();
  for (int i = 0; i < sz; i++) {
    Vec3d v = geom->verts(i) - centre;
    if (v.len2() > inf2) {
      v = v.with_len(inf);
      geom->verts(i) = centre + v;
      ideals[i] = geom->verts().size();
      int idx = geom->add_vert(centre - v);
      geom->colors(VERTS).set(idx, geom->colors(VERTS).get(i));
    }
  }

  map<int, int>::const_iterator mi;
  sz = geom->faces().size();
  for (int i = 0; i < sz; i++) {
    vector<vector<int>> alt_faces;
    alt_faces.push_back(geom->faces(i));
    for (unsigned int j = 0; j < geom->faces(i).size(); j++) {
      mi = ideals.find(geom->faces(i, j));
      if (mi != ideals.end()) {
        unsigned int af_sz = alt_faces.size();
        for (unsigned int k = 0; k < af_sz; k++) {
          alt_faces.push_back(alt_faces[k]);
          alt_faces.back()[j] = mi->second;
        }
      }
    }
    for (unsigned int k = 1; k < alt_faces.size(); k++)
      geom->raw_faces()[i].insert(geom->raw_faces()[i].end(),
                                  alt_faces[k].begin(), alt_faces[k].end());
  }

  sz = geom->edges().size();
  for (int i = 0; i < sz; i++) {
    vector<vector<int>> alt_edges;
    alt_edges.push_back(geom->edges(i));
    for (unsigned int j = 0; j < geom->edges(i).size(); j++) {
      mi = ideals.find(geom->edges(i, j));
      if (mi != ideals.end()) {
        unsigned int ae_sz = alt_edges.size();
        for (unsigned int k = 0; k < ae_sz; k++) {
          alt_edges.push_back(alt_edges[k]);
          alt_edges.back()[j] = mi->second;
        }
      }
    }
    for (unsigned int k = 1; k < alt_edges.size(); k++) {
      int idx = geom->add_edge(alt_edges[k]);
      geom->colors(EDGES).set(idx, geom->colors(EDGES).get(i));
    }
  }
}

void transform_and_repeat(Geometry *geom, string sym_to, string sym_from,
                          Trans3d pos)
{
  Transformations ts;
  ts.min_set(Symmetry(sym_to).get_trans(), Symmetry(sym_from).get_trans(), pos);
  geom->transform(pos);
  sym_repeat(geom, *geom, ts, ELEM_VERTS | ELEM_EDGES | ELEM_FACES);
}

} // namespace anti
