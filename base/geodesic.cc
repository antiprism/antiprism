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
   Name: geodesics.cc
   Description: Generate geodesic spheres and polyhedra
   Project: Antiprism - http://www.antiprism.com
*/

#include "geometryutils.h"
#include "mathutils.h"
#include "private_geodesic.h"
#include "private_misc.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <vector>

using std::map;
using std::swap;
using std::vector;

using namespace anti;

std::string IJPos::dump()
{
  const char *txt[] = {"face", "e2", "e0", "v0", "e1", "v2", "v1", "out"};
  return txt[pos];
}

void min_idx_first(vector<int> &tri)
{
  int first;
  if (tri[0] < tri[1] && tri[0] < tri[2])
    first = 0;
  else if (tri[1] < tri[2])
    first = 1;
  else
    first = 2;

  rotate(tri.begin(), tri.begin() + first, tri.end());
}

inline int_pr Geodesic::normal_crds(int_pr crds)
{
  int_pr n_crds;
  if (crds.second < 0)
    n_crds = crds;
  else if (crds.first < 0)
    n_crds = rot_f(crds);
  else
    n_crds = rot_f(rot_f(crds));

  return rot_e0(n_crds);
}

inline int Geodesic::get_edge_index(int v0, int v1)
{
  auto mi = edge_idx.find(make_edge(v0, v1));
  return (mi != edge_idx.end()) ? mi->second : -1;
}

vector<int> Geodesic::make_face_indexes(int i, const vector<int> &face)
{
  vector<int> indx(7);
  indx[0] = face[0];
  indx[1] = face[1];
  indx[2] = face[2];
  indx[3] = get_edge_index(face[0], face[1]);
  indx[4] = get_edge_index(face[1], face[2]);
  indx[5] = get_edge_index(face[2], face[0]);
  indx[6] = i;
  return indx;
}

inline int Geodesic::index_map(int i, int j, const vector<int> &indx, int p_idx)
{
  IJPos pos = get_pos(i, j);
  // polyhedron vertex
  if (pos.is_vert()) {
    if (pos == IJPos::v0)
      return indx[0];
    else if (pos == IJPos::v1)
      return indx[1];
    else
      return indx[2];
  }

  int x = grid_x(i, j);
  int y = grid_y(i, j);

  // polyhedron edge interior
  int V_sz = base.verts().size();
  if (pos.is_edge()) {
    if (pos == IJPos::e0) {
      if (indx[0] < indx[1])
        return V_sz + (F - 1) * indx[3] + F * x / freq - 1;
      else
        return V_sz + (F - 1) * indx[3] + F - F * x / freq - 1;
    }
    else if (pos == IJPos::e1) {
      if (indx[1] < indx[2])
        return V_sz + (F - 1) * indx[4] + F * y / freq - 1;
      else
        return V_sz + (F - 1) * indx[4] + F - F * y / freq - 1;
    }
    else {
      if (indx[2] < indx[0])
        return V_sz + (F - 1) * indx[5] + F - F * y / freq - 1;
      else
        return V_sz + (F - 1) * indx[5] + F * y / freq - 1;
    }
  }

  // polyhedron face interior
  if (pos.is_face()) {
    if (p_idx == noindex) {
      auto mi = grid_idxs.find(int_pr(i, j));
      if (mi == grid_idxs.end()) {
        // fprintf(stderr, "not found in grid_idxs i=%d, j=%d\n", i, j);
        return noindex;
      }
      p_idx = mi->second;
    }
    int indx_no = V_sz + (F - 1) * base.edges().size() +
                  (F * F * (m * m + m * n + n * n) - F * 3 + 2) / 2 * indx[6] +
                  p_idx;
    return indx_no;
  }

  if (pos.is_out()) {
    vector<int> edge(2);
    if (y < 0) {
      edge[0] = indx[0];
      edge[1] = indx[1];
    }
    else if (x + y > freq) {
      edge[0] = indx[1];
      edge[1] = indx[2];
    }
    else { // x<0
      edge[0] = indx[2];
      edge[1] = indx[0];
    }

    int nf_pos = 1;
    if (edge[0] > edge[1]) {
      swap(edge[0], edge[1]);
      nf_pos = 0;
    }

    auto mi_ef = edge_faces.find(edge);
    int nf_idx = mi_ef->second[nf_pos]; // index of neighbouring face
    if (nf_idx == -1)                   // no neighbouring face
      return noindex;

    const vector<int> &nface = base.faces()[nf_idx];
    int_pr n_crds = normal_crds(mk_int_pr(x, y));
    int new_origin;
    if (y < 0)
      new_origin = indx[1];
    else if (x < 0)
      new_origin = indx[0];
    else
      new_origin = indx[2];

    if (new_origin == nface[1])
      n_crds = rot_f(n_crds);
    else if (new_origin == nface[2])
      n_crds = rot_f(rot_f(n_crds));

    // if(n_crds.second < n/2)
    //   return noindex;

    vector<int> n_indx = make_face_indexes(nf_idx, nface);
    return index_map(coord_i(n_crds), coord_j(n_crds), n_indx);
  }

  return noindex; // should never get here!
}

Geodesic::Geodesic(const Geometry &base_poly, int mm, int nn, char mthd,
                   Vec3d cen)
    : base(base_poly), m(mm), n(nn), method(mthd), centre(cen)
{
  init();
}

void Geodesic::init()
{
  // "Normalise" the pattern"
  int fact = gcd(m, n);
  m /= fact;
  n /= fact;
  freq = fact * (m * m + m * n + n * n);

  // triangulate_basic(base, true, 0);
  triangulate_basic(base, false, 0);
  base.add_missing_impl_edges();
  for (unsigned int i = 0; i < base.faces().size(); i++) {
    min_idx_first(base.faces(i));
  }

  for (unsigned int i = 0; i < base.edges().size(); i++)
    edge_idx[base.edges(i)] = i;

  // fprintf(stderr, "edges.size()=%d\n", edges.size());
  edge_faces = base.get_edge_face_pairs();
  map<vector<int>, vector<int>>::iterator mi;

  F = freq / (m * m + m * n + n * n);
  make_grid_idxs();
}

void clear_def_edges(Geometry &geom)
{
  vector<int> del_edges;
  int e_sz = geom.edges().size();
  for (int i = 0; i < e_sz; i++) {
    if (!geom.colors(EDGES).get(i).is_set())
      del_edges.push_back(i);
  }
  geom.del(EDGES, del_edges);
}

void Geodesic::make_geo(Geometry &geo)
{
  vector<Vec3d> &gverts = geo.raw_verts();

  gverts = base.verts();
  geo.colors(VERTS) = base.colors(VERTS);

  if (method == 's')
    project_onto_sphere(geo, centre);

  int num_verts =
      base.faces().size() * (F * F * (m * m + m * n + n * n) - F * 3 + 2) / 2 +
      base.edges().size() * (F - 1) + base.verts().size();
  // fprintf(stderr, "num_verts = %d\n", num_verts);
  // fprintf(stderr, "num_verts = %d, fs=%d, es=%d, vs=%d, F=%d, m=%d, n=%d,
  // freq=%d\n", num_verts, faces.size(), edges.size(), verts.size(), F, m, n,
  // freq);

  gverts.resize(num_verts);

  vector<vector<int>> orig_edges;
  for (unsigned int i = 0; i < base.faces().size(); i++) {
    vector<int> indx = make_face_indexes(i, base.faces(i));
    grid_to_points(indx, gverts);
    vector<vector<int>> new_face_tris;
    grid_to_tris(indx, new_face_tris, orig_edges);
    vector<vector<int>> &gfaces = geo.raw_faces();
    gfaces.insert(gfaces.end(), new_face_tris.begin(), new_face_tris.end());
    Color f_col = base.colors(FACES).get(int(i));
    for (unsigned int f = gfaces.size() - new_face_tris.size();
         f < gfaces.size(); f++)
      geo.colors(FACES).set(int(f), f_col);
  }

  geo.add_missing_impl_edges();
  // fprintf(stderr, "orig_edges.size()= %d\n", orig_edges.size());
  for (auto &orig_edge : orig_edges) {
    vector<int> edge(2);
    edge[0] = orig_edge[0];
    edge[1] = orig_edge[1];
    const vector<vector<int>> &gedges = geo.edges();
    int e_idx = find(gedges.begin(), gedges.end(), edge) - gedges.begin();
    Color e_col = base.colors(EDGES).get(orig_edge[2]);
    geo.colors(EDGES).set(e_idx, e_col);
    // char buf[MSG_SZ];
    // vtostr(buf, e_col.get_Vec4d());
    // fprintf(stderr, "adding edge (%d, %d) col=%s\n", edge[0], edge[1],
    // buf);
  }
  clear_def_edges(geo);
}

void Geodesic::grid_to_points(vector<int> indx, vector<Vec3d> &gverts)
{
  const vector<int> face = base.faces(indx[6]);
  // fprintf(stderr, "\n+++++++\t\t\t\tface %d = (%d, %d, %d)\n", indx[6],
  // face[0], face[1], face[2]);

  vector<vector<Vec3d>> v(3);
  for (int vtx = 0; vtx < 3; vtx++) {
    Vec3d A = gverts[face[vtx]];
    Vec3d B = gverts[face[(vtx + 1) % 3]];
    Vec3d edge_vec = B - A;
    if (method == 'p') {
      v[vtx].push_back(Vec3d(0, 0, 0));
      for (int i = 1; i < freq + 1; i++)
        v[vtx].push_back(edge_vec * (double(i) / freq));
    }
    else if (method == 's') {
      // find X, nearest point to origin, O, on line AB
      Vec3d X = (A + B) / 2.0;
      double Xmag = X.len();
      double ang_XOA = -acos(safe_for_trig(vdot(A, X) / Xmag));
      double ang_XOB = acos(safe_for_trig(vdot(B, X) / Xmag));
      double ang_AOB = ang_XOB - ang_XOA;
      // fprintf(stderr, "\nk=%f\nXOA=%f\nXOB=%f\nAOB = %f\nlen edgevec=%f\n",
      // k, ang_XOA, ang_XOB, ang_AOB, edge_vec.len());
      Vec3d unit_edge_vec = edge_vec;
      unit_edge_vec.to_unit();
      for (int i = 0; i < freq + 1; i++) {
        double len = tan(ang_XOA + ang_AOB * double(i) / freq) * Xmag;
        // fprintf(stderr, "\n\nlen = %f\n", len);
        // v[vtx].push_back(unit_edge_vec*len + X - A);
        v[vtx].push_back((X + unit_edge_vec * len).unit());
        // v[vtx].back().dump("v[vtx]");
      }
      // fprintf(stderr,"\n");
    }
  }

  double sign =
      1 - 2 * (vdot(v[1][0], vcross(v[0][0] - v[1][0], v[1][0] - v[2][0])) > 0);
  int test_val = 2 * freq / (m + n);
  // i and j at twice the corner angle (which lies inside the axes)
  for (int i = 0; i < test_val; i++)
    for (int j = 0; j < test_val; j++) {
      IJPos pos = get_pos(i, j);
      if (pos.is_out() || pos.is_vert())
        continue;

      int x = grid_x(i, j);
      int y = grid_y(i, j);
      int n[] = {x, y, freq - x - y};
      Vec3d pt;
      if (method == 'p') {
        Vec3d v_delta = v[0][n[0]] + v[(0 - 1 + 3) % 3][freq - n[(0 + 1) % 3]] -
                        v[(0 - 1 + 3) % 3][freq];
        pt = gverts[face[0]] + v_delta;
      }
      else if (method == 's') {
        Vec3d lnorms[3];
        for (int k = 0; k < 3; k++) {
          Vec3d p1 = v[k][n[k]];
          Vec3d p2 = v[(k + 1) % 3][freq - n[k]];
          lnorms[k] = vcross(p1, p2);
        }
        pt = Vec3d(0, 0, 0);
        for (int k = 0; k < 3; k++)
          pt += vcross(sign * lnorms[k], lnorms[(k - 1 + 3) % 3]).unit();
        pt.to_unit();
      }

      int idx = index_map(i, j, indx);
      // fprintf(stderr, "gverts.size()=%d, idx=%d\n", gverts.size(), idx);
      gverts[idx] = pt;
      // pt.dump("pt");
    }
}

void Geodesic::make_grid_idxs()
{
  int test_val = 2 * freq / (m + n);
  int idx = 0;
  // i and j at twice the corner angle (which lies inside the axes)
  for (int i = 0; i < test_val - 1; i++)
    for (int j = 0; j < test_val - 1; j++)
      if (get_pos(i, j).is_face())
        grid_idxs[int_pr(i, j)] = idx++;
}

inline vector<int> make_tri(int v0, int v1, int v2)
{
  vector<int> tri(3);
  tri[0] = v0;
  tri[1] = v1;
  tri[2] = v2;
  return tri;
}

int orig_edge(IJPos p0_pos, int p0_idx, IJPos p1_pos, int p1_idx,
              vector<int> &indx, vector<int> &e_col)
{
  const int e_to_indx[] = {0, 5, 3, 1, 4, 3, 2, 7};
  int e_no = 0;
  if (p0_pos.is_edge() && (p1_pos.is_vert() || p0_pos == p1_pos))
    e_no = p0_pos.to_int();
  else if (p1_pos.is_edge() && p0_pos.is_vert())
    e_no = p1_pos.to_int();
  else if (p0_pos.is_vert() && p1_pos.is_vert()) { // F1 Class I !
    int vpos_sum = p0_pos.to_int() + p1_pos.to_int();
    if (vpos_sum == 9)
      e_no = 2;
    else if (vpos_sum == 11)
      e_no = 4;
    else
      e_no = 1;
  }

  // fprintf(stderr, "    p0=%s, p1=%s, (%d, %d), e_no=%d, %d\n",
  // p0_pos.dump().c_str(), p1_pos.dump().c_str(), p0_idx, p1_idx, e_no,
  // indx[e_to_indx[e_no]]);
  if (e_no) {
    e_col.clear();
    e_col.push_back(p0_idx);
    e_col.push_back(p1_idx);
    if (e_col[0] > e_col[1])
      swap(e_col[0], e_col[1]);
    e_col.push_back(indx[e_to_indx[e_no]]);
  }

  return e_no;
}

void add_orig_edges(IJPos p0_pos, int p0_idx, IJPos p1_pos, int p1_idx,
                    IJPos p2_pos, int p2_idx, vector<int> &indx,
                    vector<vector<int>> &e_cols)
{
  vector<int> e_col;
  if (orig_edge(p0_pos, p0_idx, p1_pos, p1_idx, indx, e_col))
    e_cols.push_back(e_col);
  if (orig_edge(p1_pos, p1_idx, p2_pos, p2_idx, indx, e_col))
    e_cols.push_back(e_col);
  if (orig_edge(p2_pos, p2_idx, p0_pos, p0_idx, indx, e_col))
    e_cols.push_back(e_col);
}

inline bool Geodesic::tri_test(int i, int j, int di, int dj)
{
  IJPos p0 = get_pos(i, j);
  IJPos p1 = get_pos(i + 1, j + 1);
  IJPos p2 = get_pos(i + di, j + dj);

  int out_cnt = p0.is_out() + p1.is_out() + p2.is_out();
  if (out_cnt > 1) // needs at least 2 points in to be counted
    return false;

  if (out_cnt == 0) // 3 points in so part of face
    return true;

  // face has 1 point out, and 2 points in face or on edges

  int is_face_cnt = p0.is_face() + p1.is_face() + p2.is_face();

  // 2 verts in face so can only belong to this face
  if (is_face_cnt == 2)
    return true;

  // face now has one point out, and zero or one point on face
  if (is_face_cnt == 0)
    return false;
  else
    return dj;
}

void Geodesic::grid_to_tris(vector<int> indx, vector<vector<int>> &new_tris,
                            vector<vector<int>> &orig_edges)
{
  int p0_idx, p1_idx, p2_idx, p3_idx;

  int test_val = 2 * freq / (m + n);
  for (int i = 0; i < test_val - 1; i++)
    for (int j = 0; j < test_val - 1; j++) {
      IJPos p0_pos = get_pos(i, j);
      IJPos p1_pos = get_pos(i + 1, j + 1);
      IJPos p2_pos = get_pos(i + 1, j);
      IJPos p3_pos = get_pos(i, j + 1);

      p0_idx = index_map(i, j, indx);
      p1_idx = index_map(i + 1, j + 1, indx);
      if (p0_idx != noindex && p1_idx != noindex) {
        if (tri_test(i, j, 1, 0)) {
          p2_idx = index_map(i + 1, j, indx);
          if (p2_idx != noindex) {
            new_tris.push_back(make_tri(p0_idx, p1_idx, p2_idx));
            if (m * n == 0)
              add_orig_edges(p0_pos, p0_idx, p1_pos, p1_idx, p2_pos, p2_idx,
                             indx, orig_edges);
          }
        }

        if (tri_test(i, j, 0, 1)) {
          p3_idx = index_map(i, j + 1, indx);
          if (p3_idx != noindex) {
            new_tris.push_back(make_tri(p1_idx, p0_idx, p3_idx));
            if (m * n == 0)
              add_orig_edges(p1_pos, p1_idx, p0_pos, p0_idx, p3_pos, p3_idx,
                             indx, orig_edges);
          }
        }
      }

      // fprintf(stderr, "0=%d %s, 1=%d %s, 2=%d %s, 3=%d %s\n", p0_idx,
      // p0_pos.dump().c_str(), p1_idx, p1_pos.dump().c_str(), p2_idx,
      // p2_pos.dump().c_str(), p3_idx, p3_pos.dump().c_str());
    }
}
