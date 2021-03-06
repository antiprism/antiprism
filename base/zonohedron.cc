/*
   Copyright (c) 2003-2021, Adrian Rossiter

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
   Name: zono.cc
   Description: creating zonohedra
   Project: Antiprism - http://www.antiprism.com
*/

#include "geometry.h"
#include "geometryinfo.h"
#include "polygon.h"
#include "private_misc.h"
#include "utils.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <vector>

using std::logical_not;
using std::map;
using std::pair;
using std::set;
using std::string;
using std::swap;
using std::vector;

namespace anti {

struct vec_less {
  bool operator()(const Vec3d &v1, const Vec3d &v2) const
  {
    return compare(v1, v2, anti::epsilon) ==
           -1; // may need to be larger for Qhull
  }
};

// Normalised direction of the line of a vector
Vec3d normalised_dir(const Vec3d v)
{
  if (v[0] < -anti::epsilon)
    return -v;
  else if (v[0] < anti::epsilon) {
    if (v[1] < -anti::epsilon)
      return -v;
    else if (v[1] < anti::epsilon) {
      if (v[2] < -anti::epsilon)
        return -v;
    }
  }
  return v;
}

vector<Vec3d> get_star(const Geometry &geom, char type, Vec3d centre)
{
  vector<Vec3d> star;
  const vector<Vec3d> &verts = geom.verts();
  switch (type) {
  case 'v': // vertices
    for (const auto &vert : verts) {
      Vec3d v = vert - centre;
      if (v.len2() > anti::epsilon * anti::epsilon)
        star.push_back(v);
    }
    break;
  case 'a': // any vector between two vertices
    for (unsigned int i = 0; i < verts.size(); i++)
      for (unsigned int j = i + 1; j < verts.size(); j++)
        star.push_back(geom.edge_vec(i, j));
    break;
  case 'e': // explicit edge vectors
    for (unsigned int i = 0; i < geom.edges().size(); i++)
      star.push_back(geom.edge_vec(i));
    break;
  case 'i': { // implict edge vectors
    GeometryInfo info(geom);
    const vector<vector<int>> edges = info.get_impl_edges();
    for (const auto &edge : edges)
      star.push_back(geom.edge_vec(edge));
    break;
  }
  }
  return star;
}

void make_zonohedron_1d(Geometry &zono, const vector<Vec3d> &star)
{
  zono.clear_all();
  if (!star.size())
    return;
  Vec3d pos = Vec3d(0, 0, 0);
  Vec3d neg = Vec3d(0, 0, 0);
  for (unsigned int i = 0; i < star.size(); i++) {
    double cos_a = vdot(star[0], star[i]);
    if (cos_a > 0)
      pos += star[i];
    else
      neg += star[i];
  }
  zono.add_vert(pos);
  zono.add_vert(neg);
}

void make_zonohedron_2d(Geometry &zono, const vector<Vec3d> &star, Vec3d fnorm)
{
  zono.clear_all();
  map<Vec3d, set<int>, vec_less> stars;
  for (unsigned int i = 0; i < star.size(); i++)
    stars[normalised_dir(star[i]).unit()].insert(i);

  map<Vec3d, set<int>>::iterator mi;
  for (mi = stars.begin(); mi != stars.end(); ++mi) {
    vector<Vec3d> z_face_star;
    Vec3d pos = Vec3d(0, 0, 0);
    Vec3d neg = Vec3d(0, 0, 0);
    const Vec3d &norm = vcross(mi->first, fnorm).unit();
    for (int i = 0; i < (int)star.size(); i++) {
      if (mi->second.find(i) != mi->second.end())
        z_face_star.push_back(star[i]);
      else {
        double cos_a = vdot(norm, star[i]);
        if (cos_a > 0)
          pos += star[i];
        else
          neg += star[i];
      }
    }
    Geometry zono_face;
    make_zonohedron_1d(zono_face, z_face_star);
    for (unsigned int k = 0; k < zono_face.verts().size(); k++) {
      zono.add_vert(pos + zono_face.verts(k));
      zono.add_vert(neg + zono_face.verts(k));
    }
  }
}

Status make_zonohedron(Geometry &geom, const vector<Vec3d> &star)
{
  geom.clear_all();

  vector<double> star_mags(star.size());
  for (unsigned int i = 0; i < star.size(); i++)
    star_mags[i] = star[i].len();

  map<Vec3d, set<int>, vec_less> stars2d;
  if (star.size()) {
    for (unsigned int i = 0; i < star.size() - 1; i++)
      for (unsigned int j = i + 1; j < star.size(); j++) {
        Vec3d n = vcross(star[i], star[j]) / (star_mags[i] * star_mags[j]);
        if (n.len2() > anti::epsilon * anti::epsilon) {
          Vec3d norm = normalised_dir(n).unit();
          stars2d[norm].insert(i);
          stars2d[norm].insert(j);
        }
      }
  }

  if (stars2d.size()) {
    map<Vec3d, set<int>>::iterator mi;
    for (mi = stars2d.begin(); mi != stars2d.end(); ++mi) {
      vector<Vec3d> z_face_star;
      Vec3d pos = Vec3d(0, 0, 0);
      Vec3d neg = Vec3d(0, 0, 0);
      const Vec3d &norm = mi->first;
      for (int i = 0; i < (int)star.size(); i++) {
        if (mi->second.find(i) != mi->second.end())
          z_face_star.push_back(star[i]);
        else {
          double cos_a = vdot(norm, star[i]);
          if (cos_a > 0)
            pos += star[i];
          else
            neg += star[i];
        }
      }

      Geometry zono_face;
      make_zonohedron_2d(zono_face, z_face_star, norm);
      for (unsigned int k = 0; k < zono_face.verts().size(); k++) {
        geom.add_vert(pos + zono_face.verts(k));
        geom.add_vert(neg + zono_face.verts(k));
      }
    }
  }
  else {
    make_zonohedron_1d(geom, star);
  }

  if (geom.is_set())
    return geom.set_hull("A0.9999999");
  else {
    geom.add_vert(Vec3d::zero);
    return Status::ok();
  }
}

Status make_zonohedrified_polyhedron(Geometry &geom, const Geometry &seed,
                                     const vector<Vec3d> &star, Color col)
{
  geom.clear_all();
  geom.append(seed);

  // Store original face colours by normal
  std::map<Vec3d, Color, vec_less> orig_cols;
  for (unsigned int i = 0; i < geom.faces().size(); i++)
    orig_cols[geom.face_norm(i).to_unit()] = geom.colors(FACES).get(i);

  for (const auto &i : star) {
    int v_sz = geom.verts().size();
    geom.raw_verts().resize(v_sz * 2);
    for (int j = 0; j < v_sz; j++) {
      geom.raw_verts()[j + v_sz] = geom.verts(j) + i;
    }
    Status stat = geom.set_hull("");
    if (!stat)
      return stat;
  }

  // Restore original face colours by normal
  for (unsigned int i = 0; i < geom.faces().size(); i++) {
    auto mi = orig_cols.find(geom.face_norm(i).to_unit());
    geom.colors(FACES).set(i, (mi != orig_cols.end()) ? mi->second : col);
  }

  return Status::ok();
}

static inline int pos_mod(int a, int b) { return (a % b + b) % b; }

static int get_idx(int P, int s, int s_step, int num_spirals, int i, int j,
                   int V)
{
  // each spiral is a rectangular array of quads: P-s_step x s_step
  //  with 0 <= i < P-s_tep and 0 <= y < s_step
  //
  // fprintf(stderr, "   (%d, %d) -> ", i, j);
  int idx = -999999;
  if (i < 0) {
    if (j == 0) {
      idx = -2; // first point
      // fprintf(stderr, "idx_a=%d\n", idx);
    }
    else if (j > P - s_step) {
      // Move to next spiral and go down i points
      idx = ((pos_mod(s - 1, num_spirals) + 1) * (P - s_step) * s_step) + j - P;
      // fprintf(stderr, "idx_b0=%d\n", idx);
    }
    else {
      // Move to previous spiral and go down j points
      idx = (pos_mod(s - 1, num_spirals) * (P - s_step) * s_step) +
            (j - 1) * s_step;
      // fprintf(stderr, "idx_b1=%d (%d->%d : %d->%d)\n", idx, pos_mod(s - 1,
      // num_spirals), pos_mod(s - 1, num_spirals) * (P - s_step) * s_step, j,
      // (j
      // - 1)*s_step);
    }
  }
  else if (j == s_step) {
    if (i == P - s_step - 1) {
      idx = -1;
      // fprintf(stderr, "idx_c=%d\n", idx);
    }
    else if (i >= P - 2 * s_step) {
      // Move to next spiral and go down i points
      idx = (pos_mod(s - 1, num_spirals) * (P - s_step) * s_step) +
            (P - s_step) * s_step + i - (P - s_step - 1);
      // fprintf(stderr, "idx_d=%d\n", idx);
    }
    else {
      // Move to next spiral and go down i points
      idx = (pos_mod(s - 1, num_spirals) * (P - s_step) * s_step) +
            (i + j) * s_step;
      // fprintf(stderr, "idx_e=%d\n", idx);
    }
  }
  else {
    idx = (s * (P - s_step) * s_step) + (i * s_step) + j;
    // fprintf(stderr, "idx_f=%d\n", idx);
  }

  // if(idx > P*(P-s_step)-1)
  // if(idx>47)
  //  idx=-1;

  return V + idx + 2;
}

static inline vector<int> get_face(int P, int s, int s_step, int num_spirals,
                                   int i, int j, int V)
{
  vector<int> face = {get_idx(P, s, s_step, num_spirals, i - 1, j, V),
                      get_idx(P, s, s_step, num_spirals, i, j, V),
                      get_idx(P, s, s_step, num_spirals, i, j + 1, V),
                      get_idx(P, s, s_step, num_spirals, i - 1, j + 1, V)};
  // fprintf(stderr, "(%d, %d) -> %d, %d, %d, %d\n", i, j, face[0], face[1],
  // face[2], face[3]);
  return face;
}

Status make_polar_zonohedron(Geometry &geom, const vector<Vec3d> &star,
                             int step, int spiral_step)
{
  geom.clear_all();
  int N = star.size();       // number of vectors
  int D = step;              // step between vectors
  int num_parts = gcd(N, D); // number of parts in final model

  int P = N / num_parts; // number of vectors in part
  int P_spiral_step =
      (spiral_step == 0) ? 1 : pos_mod(spiral_step / num_parts, P);
  int P_num_spirals = P / gcd(P, P_spiral_step);

  for (int p = 0; p < num_parts; p++) {
    int V = geom.verts().size();
    vector<Vec3d> star_part;
    for (int i = 0; i < P; i++)
      star_part.push_back(star[(p * P / P_num_spirals + i * D) % N]);

    geom.add_vert(Vec3d(0, 0, 0)); // initial point
    geom.add_vert(Vec3d(0, 0, 0)); // final point, set later
    Vec3d A;                       // points along this spiral
    Vec3d B;                       // points along following spiral
    for (int s = 0; s < P_num_spirals; s++) {
      A = Vec3d(0, 0, 0);
      for (int i = 0; i < P - P_spiral_step; i++) {
        int i_idx = (s * P_spiral_step + i) % P;
        A += star_part[i_idx]; // next point on this spiral

        B = Vec3d(0, 0, 0);
        for (int j = 0; j < P_spiral_step; j++) {
          geom.add_vert(A + B);
          // index of next star vector on following spiral
          int j_idx = pos_mod((s - 1) * P_spiral_step + j, P);
          B += star_part[j_idx]; // next point on following spiral
          geom.add_face(get_face(P, s, P_spiral_step, P_num_spirals, i, j, V),
                        Color(p * P_num_spirals + s));
        }
        // index of next star vector on this spiral
      }
    }
    geom.verts(V + 1) = A + B; // set final point
  }
  return Status::ok();
}

Status make_translation_surface(Geometry &geom, const Geometry &star0,
                                const Geometry &star1, const string open_flgs)
{
  // set open_type to two values from d, o, c for default (detect open
  // or close), close (add a close vector), open (don't add a close vector)
  vector<char> open_type(2, open_flgs.empty() ? 'd' : open_flgs[0]);
  if (open_flgs.size() > 1)
    open_type[1] = open_flgs[1];

  // Vectors that will be used, and whether they are a loop (i.e. is the
  // last point reached by the chain identified with the initial point)
  vector<vector<Vec3d>> chains(2, {Vec3d::zero});
  vector<bool> is_loop(2);

  const Geometry *stars[] = {&star0, &star1};
  for (int s = 0; s < 2; s++) {
    Vec3d sum = Vec3d::zero;
    for (size_t i = 0; i < stars[s]->verts().size(); i++)
      sum += stars[s]->verts(i);

    chains[s].insert(chains[s].end(), stars[s]->verts().begin(),
                     stars[s]->verts().end());
    // detect whether vectors form a loop
    is_loop[s] = compare(sum / stars[s]->verts().size(), Vec3d::zero) == 0;
    if (open_type[s] == 'd') {
      if (is_loop[s]) {
        // the last vertex is the same as the first and so can be deleted
        chains[s].resize(chains[s].size() - 1);
      }
    }
    else
      is_loop[s] = (open_type[s] == 'c');
  }

  const int v_sz[2] = {(int)chains[0].size(), (int)chains[1].size()};
  Vec3d cur0 = Vec3d::zero;
  for (int i = 0; i < v_sz[0]; i++) {
    cur0 += chains[0][i];
    Vec3d cur1 = Vec3d::zero;
    for (int j = 0; j < v_sz[1]; j++) {
      cur1 += chains[1][j];
      geom.add_vert(cur0 + cur1);
      vector<int> face = {i * v_sz[1] + j, i * v_sz[1] + (j + 1),
                          (i + 1) * v_sz[1] + (j + 1), (i + 1) * v_sz[1] + j};
      if (i + 1 < v_sz[0] && j + 1 < v_sz[1]) {
        geom.add_face(face);
        geom.add_edge({face[0], face[3]}, stars[0]->colors(VERTS).get(i));
        geom.add_edge({face[0], face[1]}, stars[1]->colors(VERTS).get(j));
      }
      else if (i + 1 < v_sz[0])
        geom.add_edge({face[0], face[3]}, stars[0]->colors(VERTS).get(i));
      else if (j + 1 < v_sz[1])
        geom.add_edge({face[0], face[1]}, stars[1]->colors(VERTS).get(j));
      bool add_close_i_edge = (is_loop[0] && i + 1 == v_sz[0]);
      bool add_close_j_edge = (is_loop[1] && j + 1 == v_sz[1]);
      if (add_close_i_edge || add_close_j_edge) {
        const int i_next = (i + 1) % v_sz[0];
        const int j_next = (j + 1) % v_sz[1];
        vector<int> face2 = {i * v_sz[1] + j, i * v_sz[1] + j_next,
                             i_next * v_sz[1] + j_next, i_next * v_sz[1] + j};
        if (i + 1 != v_sz[0] || j + 1 != v_sz[1] || (is_loop[0] && is_loop[1]))
          geom.add_face(face2);

        if (is_loop[0])
          geom.add_edge({face2[0], face2[3]}, stars[0]->colors(VERTS).get(i));
        if (is_loop[1])
          geom.add_edge({face2[0], face2[1]}, stars[1]->colors(VERTS).get(j));
      }
    }
  }
  return Status::ok();
}

} // namespace anti
