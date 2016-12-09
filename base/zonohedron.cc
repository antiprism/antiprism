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
   Name: zono.cc
   Description: creating zonohedra
   Project: Antiprism - http://www.antiprism.com
*/

#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "geometry.h"
#include "geometryinfo.h"
#include "polygon.h"
#include "private_misc.h"
#include "utils.h"

using std::vector;
using std::set;
using std::map;
using std::pair;
using std::swap;
using std::logical_not;
using std::string;

namespace anti {

struct vec_less {
  bool operator()(const Vec3d &v1, const Vec3d &v2)
  {
    return compare(v1, v2, epsilon) == -1; // may need to be larger for Qhull
  }
};

// Normalised direction of the line of a vector
Vec3d normalised_dir(const Vec3d v)
{
  if (v[0] < -epsilon)
    return -v;
  else if (v[0] < epsilon) {
    if (v[1] < -epsilon)
      return -v;
    else if (v[1] < epsilon) {
      if (v[2] < -epsilon)
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
      if (v.len2() > epsilon * epsilon)
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

void make_zonohedron_1d(Geometry *zono, const vector<Vec3d> &star)
{
  zono->clear_all();
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
  zono->add_vert(pos);
  zono->add_vert(neg);
}

void make_zonohedron_2d(Geometry *zono, const vector<Vec3d> &star, Vec3d fnorm)
{
  zono->clear_all();
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
    make_zonohedron_1d(&zono_face, z_face_star);
    for (unsigned int k = 0; k < zono_face.verts().size(); k++) {
      zono->add_vert(pos + zono_face.verts(k));
      zono->add_vert(neg + zono_face.verts(k));
    }
  }
}

Status make_zonohedron(Geometry *geom, const vector<Vec3d> &star)
{
  geom->clear_all();

  vector<double> star_mags(star.size());
  for (unsigned int i = 0; i < star.size(); i++)
    star_mags[i] = star[i].len();

  map<Vec3d, set<int>, vec_less> stars2d;
  if (star.size()) {
    for (unsigned int i = 0; i < star.size() - 1; i++)
      for (unsigned int j = i + 1; j < star.size(); j++) {
        Vec3d n = vcross(star[i], star[j]) / (star_mags[i] * star_mags[j]);
        if (n.len2() > epsilon * epsilon) {
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
      make_zonohedron_2d(&zono_face, z_face_star, norm);
      for (unsigned int k = 0; k < zono_face.verts().size(); k++) {
        geom->add_vert(pos + zono_face.verts(k));
        geom->add_vert(neg + zono_face.verts(k));
      }
    }
  }
  else {
    make_zonohedron_1d(geom, star);
  }

  if (geom->is_set())
    return geom->set_hull("A0.9999999");
  else {
    geom->add_vert(Vec3d::zero);
    return Status::ok();
  }
}

Status make_zonohedrified_polyhedron(Geometry *geom, const Geometry &seed,
                                     const vector<Vec3d> &star, Color col)
{
  geom->clear_all();
  geom->append(seed);

  // Store original face colours by normal
  std::map<Vec3d, Color, vec_less> orig_cols;
  for (unsigned int i = 0; i < geom->faces().size(); i++)
    orig_cols[geom->face_norm(i).to_unit()] = geom->colors(FACES).get(i);

  for (const auto &i : star) {
    int v_sz = geom->verts().size();
    geom->raw_verts().resize(v_sz * 2);
    for (int j = 0; j < v_sz; j++) {
      geom->raw_verts()[j + v_sz] = geom->verts(j) + i;
    }
    Status stat = geom->set_hull("");
    if (!stat)
      return stat;
  }

  // Restore original face colours by normal
  for (unsigned int i = 0; i < geom->faces().size(); i++) {
    auto mi = orig_cols.find(geom->face_norm(i).to_unit());
    geom->colors(FACES).set(i, (mi != orig_cols.end()) ? mi->second : col);
  }

  return Status::ok();
}

Status make_polar_zonohedron(Geometry *geom, const vector<Vec3d> &star,
                             int step)
{
  geom->clear_all();
  int N = star.size();
  int D = step;
  int num_parts = gcd(N, D);
  int P = N / num_parts;
  for (int p = 0; p < num_parts; p++) {
    int V = geom->verts().size();
    vector<Vec3d> star_part;
    for (int i = 0; i < P; i++)
      star_part.push_back(star[(p + i * D) % N]);

    geom->add_verts(star_part);
    for (int i = 1; i < P - 1; i++)
      for (int j = 0; j < P; j++)
        geom->add_vert(geom->verts(V + (i - 1) * P + j) +
                       star_part[(i + j) % P]);
    geom->add_vert(geom->verts(V + P * (P - 2)) + star_part[P - 1]);
    geom->add_vert(Vec3d(0, 0, 0));
    for (int j = 0; j < P; j++) {
      geom->add_face(V + P * (P - 1) + 1, V + j, V + j + P, V + (j + 1) % P,
                     -1);
      if (P > 2)
        geom->add_face(V + P * (P - 1), V + P * (P - 2) + j,
                       V + P * (P - 3) + (j + 1) % P,
                       V + P * (P - 2) + (j + 1) % P, -1);
    }
    for (int i = 0; i < P - 3; i++)
      for (int j = 0; j < P; j++) {
        geom->add_face(V + i * P + (j + 1) % P, V + (i + 1) * P + j,
                       V + (i + 2) * P + j, V + (i + 1) * P + (j + 1) % P, -1);
      }
  }
  return Status::ok();
}

} // namespace anti
