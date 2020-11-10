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

/*\file geom.cc
 * \brief Classes to represent a geometry
 */

#include "geometry.h"
#include "coloring.h"
#include "geometryinfo.h"
#include "private_misc.h"
#include "private_off_file.h"
#include "private_std_polys.h"
#include "symmetry.h"
#include "utils.h"

#include <algorithm>

using std::map;
using std::pair;
using std::string;
using std::swap;
using std::vector;

namespace anti {

// Geometry implementation

int Geometry::add_vert(Vec3d vert, Color col)
{
  int idx = verts().size();
  raw_verts().push_back(vert);
  if (col.is_set())
    colors(VERTS).set(idx, col);
  return idx;
}

int Geometry::add_verts(const vector<Vec3d> &vrts)
{
  for (const auto &v : vrts)
    add_vert(v);
  return verts().size() - 1;
}

int Geometry::add_edge_raw(const vector<int> &edge, Color col)
{
  int idx = edges().size();
  raw_edges().push_back(edge);
  if (col.is_set())
    colors(EDGES).set(idx, col);
  return idx;
}

int Geometry::add_edges_raw(const vector<vector<int>> &edgs)
{
  for (const auto &e : edgs)
    add_edge_raw(e);
  return edges().size() - 1;
}

int Geometry::add_edge(std::vector<int> edge, Color col)
{
  int idx;
  if (edge[0] > edge[1])
    swap(edge[0], edge[1]);
  auto ei = find(edges().begin(), edges().end(), edge);
  if (ei != edges().end()) {
    idx = ei - edges().begin();
    colors(EDGES).set(idx, col);
  }
  else {
    idx = edges().size();
    raw_edges().push_back(edge);
    if (col.is_set())
      colors(EDGES).set(idx, col);
  }
  return idx;
}

int Geometry::add_edge(int v_idx1, int v_idx2, Color col)
{
  return add_edge(make_edge(v_idx1, v_idx2), col);
}

int Geometry::add_edges(const vector<vector<int>> &edgs)
{
  for (const auto &e : edgs)
    add_edge(e);
  return edges().size() - 1;
}

int Geometry::add_face(const std::vector<int> &face, Color col)
{
  int idx = faces().size();
  raw_faces().push_back(face);
  if (col.is_set())
    colors(FACES).set(idx, col);
  return idx;
}

int Geometry::add_faces(const vector<vector<int>> &fces)
{
  for (const auto &f : fces)
    add_face(f);
  return faces().size() - 1;
}

static void delete_verts(Geometry *geom, const vector<int> &v_nos,
                         map<int, int> *vert_map)
{
  vector<int> dels = v_nos;
  map<int, int> tmp;
  map<int, int> *v_map = (vert_map) ? vert_map : &tmp;
  v_map->clear();
  if (!dels.size())
    return;
  sort(dels.begin(), dels.end());
  unsigned int del_verts_cnt = 0;
  int map_to;
  for (unsigned int i = 0; i < geom->verts().size(); i++) {
    if (del_verts_cnt < dels.size() && (int)i == dels[del_verts_cnt]) {
      del_verts_cnt++;
      map_to = -1;
    }
    else {
      map_to = i - del_verts_cnt;
      geom->verts(map_to) = geom->verts(i);
    }
    (*v_map)[i] = map_to;
  }
  geom->raw_verts().resize(geom->verts().size() - del_verts_cnt);

  vector<int> del_faces;
  for (unsigned int i = 0; i < geom->faces().size(); i++) {
    int curr_idx = 0;
    for (unsigned int j = 0; j < geom->faces(i).size(); j++) {
      int new_idx = (*v_map)[geom->faces(i, j)];
      if (new_idx >= 0)
        geom->raw_faces()[i][curr_idx++] = new_idx;
    }
    geom->raw_faces()[i].resize(curr_idx);
    if (curr_idx < 3)
      del_faces.push_back(i);
  }
  geom->del(FACES, del_faces);

  vector<int> del_edges;
  for (unsigned int i = 0; i < geom->edges().size(); i++)
    for (unsigned int j = 0; j < 2; j++) {
      int new_idx = (*v_map)[geom->edges(i, j)];
      if (new_idx >= 0)
        geom->raw_edges()[i][j] = new_idx;
      else {
        del_edges.push_back(i);
        break;
      }
    }
  geom->del(EDGES, del_edges);
}

static void delete_faces(Geometry *geom, const vector<int> &f_nos,
                         map<int, int> *face_map)
{
  vector<int> dels = f_nos;
  if (face_map)
    face_map->clear();
  if (!dels.size())
    return;
  sort(dels.begin(), dels.end());
  unsigned int del_faces_cnt = 0;
  int map_to;
  for (unsigned int i = 0; i < geom->faces().size(); i++) {
    if (del_faces_cnt < dels.size() && (int)i == dels[del_faces_cnt]) {
      del_faces_cnt++;
      map_to = -1;
    }
    else {
      map_to = i - del_faces_cnt;
      geom->raw_faces()[map_to] = geom->faces(i);
    }
    if (face_map)
      (*face_map)[i] = map_to;
  }
  geom->raw_faces().resize(geom->faces().size() - del_faces_cnt);
}

static void delete_edges(Geometry *geom, const vector<int> &e_nos,
                         map<int, int> *edge_map)
{
  vector<int> dels = e_nos;
  if (edge_map)
    edge_map->clear();
  if (!dels.size())
    return;
  sort(dels.begin(), dels.end());
  unsigned int del_edges_cnt = 0;
  int map_to;
  for (unsigned int i = 0; i < geom->edges().size(); i++) {
    if (del_edges_cnt < dels.size() && (int)i == dels[del_edges_cnt]) {
      del_edges_cnt++;
      map_to = -1;
    }
    else {
      map_to = i - del_edges_cnt;
      geom->raw_edges()[map_to] = geom->edges(i);
    }
    if (edge_map)
      (*edge_map)[i] = map_to;
  }
  geom->raw_edges().resize(geom->edges().size() - del_edges_cnt);
}

void Geometry::del(int type, const vector<int> &idxs, map<int, int> *elem_map)
{
  map<int, int> tmp;
  map<int, int> *elm_map = (elem_map) ? elem_map : &tmp;
  if (type == VERTS)
    delete_verts(this, idxs, elm_map);
  else if (type == EDGES)
    delete_edges(this, idxs, elm_map);
  else if (type == FACES)
    delete_faces(this, idxs, elm_map);
  colors(type).remap(*elm_map);
}

void Geometry::del(int type, int idx, map<int, int> *elem_map)
{
  vector<int> del_elems;
  del_elems.push_back(idx);
  del(type, del_elems, elem_map);
}

void remap_shift(vector<vector<int>> &elems, int offset)
{
  for (auto &elem : elems)
    for (int &j : elem)
      j += offset;
}

void Geometry::append(const Geometry &geom)
{
  cols.append(geom.get_cols(), verts().size(), edges().size(), faces().size());
  vector<Vec3d> g_verts = geom.verts();
  vector<vector<int>> g_faces = geom.faces();
  vector<vector<int>> g_edges = geom.edges();

  remap_shift(g_edges, verts().size());
  remap_shift(g_faces, verts().size());

  raw_verts().insert(raw_verts().end(), g_verts.begin(), g_verts.end());
  raw_edges().insert(raw_edges().end(), g_edges.begin(), g_edges.end());
  raw_faces().insert(raw_faces().end(), g_faces.begin(), g_faces.end());
}

void Geometry::clear(int type)
{
  if (type == VERTS)
    raw_verts().clear();
  else if (type == EDGES)
    raw_edges().clear();
  else if (type == FACES)
    raw_faces().clear();

  colors(type).clear();
}

void Geometry::clear_all()
{
  clear(VERTS);
  clear(EDGES);
  clear(FACES);
}

void Geometry::get_impl_edges(vector<vector<int>> &edgs) const
{
  // edgs hasn't been cleared
  for (unsigned int i = 0; i < faces().size(); ++i)
    for (unsigned int j = 0; j < faces(i).size(); ++j)
      edgs.push_back(make_edge(faces(i, j), faces_mod(i, j + 1)));

  // Clear duplicate edges (each edge appears in two faces)
  sort(edgs.begin(), edgs.end());
  auto vit = unique(edgs.begin(), edgs.end());
  edgs.erase(vit, edgs.end());
}

Status Geometry::add_hull(string qh_args, int *dim)
{
  return anti::add_hull(*this, qh_args, dim);
}

Status Geometry::set_hull(string qh_args, int *dim)
{
  return anti::set_hull(*this, qh_args, dim);
}

int Geometry::orient(vector<vector<int>> *parts)
{
  return orient_geom(*this, parts);
}

Status Geometry::orient(int type)
{
  Status stat;
  char errmsg[MSG_SZ];
  if (!orient_geom(*this, type, errmsg))
    stat.set_error(errmsg);
  else if (*errmsg)
    stat.set_warning(errmsg);
  return stat;
}

void Geometry::orient_reverse() { return ::orient_reverse(*this); }

void Geometry::sym_align() { transform(Symmetry(*this).get_to_std()); }

void Geometry::triangulate(Color col, unsigned int winding, vector<int> *fmap)
{
  anti::triangulate(*this, col, winding, fmap);
}

Status Geometry::read(string file_name)
{
  Status stat;
  char errmsg[MSG_SZ];
  if (!off_file_read(file_name, *this, errmsg))
    stat.set_error(errmsg);
  else if (*errmsg)
    stat.set_warning(errmsg);
  return stat;
}

Status Geometry::read(FILE *file)
{
  Status stat;
  char errmsg[MSG_SZ];
  if (!off_file_read(file, *this, errmsg))
    stat.set_error(errmsg);
  else if (*errmsg)
    stat.set_warning(errmsg);
  return stat;
}

Status Geometry::read_resource(string res_name)
{
  Status stat;
  char errmsg[MSG_SZ];
  if (!make_resource_geom(*this, res_name, errmsg))
    stat.set_error(errmsg);
  else if (*errmsg)
    stat.set_warning(errmsg);
  return stat;
}

Status Geometry::write(string file_name, int sig_dgts) const
{
  Status stat;
  char errmsg[MSG_SZ];
  if (!off_file_write(file_name, *this, errmsg, sig_dgts))
    stat.set_error(errmsg);
  else if (*errmsg)
    stat.set_warning(errmsg);
  return stat;
}

void Geometry::write(FILE *file, int sig_dgts) const
{
  off_file_write(file, *this, sig_dgts);
}

Status Geometry::write_crds(string file_name, const char *sep,
                            int sig_dgts) const
{
  Status stat;
  char errmsg[MSG_SZ];
  if (!crds_write(file_name, *this, errmsg, sep, sig_dgts))
    stat.set_error(errmsg);
  else if (*errmsg)
    stat.set_warning(errmsg);
  return stat;
}

void Geometry::write_crds(FILE *file, const char *sep, int sig_dgts) const
{
  crds_write(file, *this, sep, sig_dgts);
}

Status Geometry::write_obj(string file_name, string mtl_file, const char *sep,
                           int sig_dgts) const
{
  Status stat;
  char errmsg[MSG_SZ];
  if (!obj_write(file_name, mtl_file, *this, errmsg, sep, sig_dgts))
    stat.set_error(errmsg);
  else if (*errmsg)
    stat.set_warning(errmsg);
  return stat;
}

void Geometry::write_obj(FILE *file, FILE *mfile, string mtl_file,
                         const char *sep, int sig_dgts) const
{
  obj_write(file, mfile, mtl_file, *this, sep, sig_dgts);
}

GeometryInfo Geometry::get_info() const { return GeometryInfo(*this); }

vector<int> make_edge(int v_idx1, int v_idx2)
{
  vector<int> edge(2);
  edge[0] = v_idx1;
  edge[1] = v_idx2;
  if (edge[0] > edge[1])
    swap(edge[0], edge[1]);
  return edge;
}

bool Geometry::is_oriented() const
{
  std::set<pair<int, int>> edgs;
  pair<int, int> edge;
  for (unsigned int f = 0; f < faces().size(); f++)
    for (unsigned int i = 0; i < faces(f).size(); i++) {
      edge = pair<int, int>(faces(f, i), faces_mod(f, i + 1));
      if (edgs.find(edge) != edgs.end())
        return false;
      else
        edgs.insert(edge);
    }
  return true;
}

std::map<std::vector<int>, std::vector<int>>
Geometry::get_edge_face_pairs(bool oriented) const
{
  map<vector<int>, vector<int>> edge2facepr;
  vector<int> vrts(2);
  for (unsigned int i = 0; i < faces().size(); ++i) {
    for (unsigned int j = 0; j < faces(i).size(); ++j) {
      vrts[0] = faces(i, j);
      vrts[1] = faces_mod(i, j + 1);
      int face_pos = 0;
      if (vrts[0] > vrts[1]) {
        swap(vrts[0], vrts[1]);
        face_pos = 1;
      }
      if (oriented) {
        if (edge2facepr.find(vrts) == edge2facepr.end()) {
          edge2facepr[vrts].resize(2);
          edge2facepr[vrts][(face_pos + 1) % 2] = -1;
        }
        edge2facepr[vrts][face_pos] = i;
      }
      else
        edge2facepr[vrts].push_back(i);
    }
  }
  return edge2facepr;
}

void Geometry::verts_merge(map<int, int> &vmap)
{
  map<int, int>::iterator vmi;
  vector<int> del_faces;
  for (unsigned int i = 0; i < faces().size(); i++) {
    for (unsigned int j = 0; j < faces(i).size(); j++) {
      int idx = faces(i, j);
      vmi = vmap.find(idx);
      if (vmi != vmap.end())
        raw_faces()[i][j] = vmap[idx];
    }
    if (raw_faces()[i].size()) {
      auto vi = unique(raw_faces()[i].begin(), raw_faces()[i].end());
      if (faces(i, 0) == *(vi - 1))
        vi--;
      raw_faces()[i].resize(vi - faces(i).begin());
    }
    if (faces(i).size() < 3)
      del_faces.push_back(i);
  }
  del(FACES, del_faces);

  vector<int> del_edges;
  for (unsigned int i = 0; i < edges().size(); i++) {
    for (unsigned int j = 0; j < edges(i).size(); j++) {
      int idx = edges(i, j);
      vmi = vmap.find(idx);
      if (vmi != vmap.end())
        raw_edges()[i][j] = vmap[idx];
    }
    if (edges(i, 0) == edges(i, 1))
      del_edges.push_back(i);
  }
  del(EDGES, del_edges);

  vector<int> del_verts;
  for (vmi = vmap.begin(); vmi != vmap.end(); vmi++)
    del_verts.push_back(vmi->first);
  del(VERTS, del_verts);
}

void Geometry::face_angles_lengths(int f_idx, vector<double> *angles,
                                   vector<double> *lengths) const
{
  vector<double> l;
  vector<double> &lens = (lengths) ? *lengths : l;
  vector<double> &angs = *angles;
  unsigned int fsz = faces(f_idx).size();
  angs.resize(fsz);
  lens.resize(fsz);
  Vec3d norm;
  double ang_sum = 0;
  for (unsigned int j = 0; j < fsz; j++) {
    Vec3d v0 = face_v(f_idx, j) - face_v(f_idx, (j + fsz - 1) % fsz);
    Vec3d v1 = face_v(f_idx, j) - face_v(f_idx, (j + 1) % fsz);
    if (j == 0)
      lens[fsz - 1] = v0.len();
    lens[j] = v1.len();
    if (!norm.is_set())
      norm = vcross(v0, v1);
    angs[j] = acos(
        safe_for_trig(vdot(v0, v1) / (lens[(j + fsz - 1) % fsz] * lens[j])));
    double sign = vdot(norm, vcross(v0, v1));
    if (sign < 0)
      angs[j] = 2 * M_PI - angs[j];
    ang_sum += angs[j];
  }

  if (ang_sum / fsz > M_PI)
    for (unsigned int j = 0; j < fsz; j++)
      angs[j] = 2 * M_PI - angs[j];
}

void Geometry::add_missing_impl_edges(const Color &col)
{
  // save original edges and colours
  vector<vector<int>> e_edges = edges();
  ElemProps<Color> cols = colors(EDGES);

  clear(EDGES);
  get_impl_edges(raw_edges());
  if (col.is_set()) {
    for (unsigned int i = 0; i < edges().size(); i++)
      colors(EDGES).set(i, col);
  }

  // restore original edges and colours
  for (unsigned int e = 0; e < e_edges.size(); ++e)
    add_edge(e_edges[e], cols.get(e));
}

} // namespace anti
