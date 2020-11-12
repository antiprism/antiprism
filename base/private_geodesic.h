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
   Name: private_geodesic.h
   Description: Generate geodesic spheres and polyhedra
   Project: Antiprism - http://www.antiprism.com
*/

#ifndef GEODESIC_H
#define GEODESIC_H

#include "geometry.h"
#include "geometryutils.h"

#include <map>
#include <string>
#include <vector>

typedef std::pair<int, int> int_pr;
inline int_pr mk_int_pr(int x, int y)
{
  int_pr pr;
  pr.first = x;
  pr.second = y;
  return pr;
}

class IJPos {
private:
  int pos;

public:
  enum { face = 0, v0 = 3, v1 = 6, v2 = 5, e0 = 2, e1 = 4, e2 = 1, out = 7 };
  IJPos(int p) : pos(p) {}
  bool operator==(int pos2) { return pos == pos2; }
  bool operator==(IJPos p2) { return pos == p2.pos; }
  bool is_vert() { return (pos == v0 || pos == v1 || pos == v2); }
  bool is_edge() { return (pos == e0 || pos == e1 || pos == e2); }
  bool is_face() { return (pos == face); }
  bool is_out() { return (pos == out); }
  int to_int() { return pos; }
  std::string dump();
};

class Geodesic {
private:
  enum { noindex = -1 };

  anti::Geometry base;
  int freq;
  int F;
  int m;
  int n;
  char method;
  anti::Vec3d centre;

  std::map<std::vector<int>, int> edge_idx;
  std::map<std::vector<int>, std::vector<int>> edge_faces;
  // std::map<std::vector<int>, int> face_idx;
  std::map<int_pr, int> grid_idxs;

  void init();
  void sphere_projection(anti::Geometry &geom);
  void make_grid_idxs();
  int grid_x(int i, int j) { return i * (-m) + j * (m + n); }
  int grid_y(int i, int j) { return i * (m + n) + j * (-n); }
  int_pr rot_e0(int_pr crds) // half-rot about centre e0
  {
    return mk_int_pr(freq - crds.first, -crds.second);
  }
  int_pr rot_f(int_pr crds) // third-rot about centre f
  {
    return mk_int_pr(freq - crds.first - crds.second, crds.first);
  }
  int_pr normal_crds(int_pr crds);

  int coord_i(int_pr crds)
  {
    return (n * crds.first + (m + n) * crds.second) / (m * m + m * n + n * n);
  }
  int coord_j(int_pr crds)
  {
    return ((m + n) * crds.first + m * crds.second) / (m * m + m * n + n * n);
  }

  void grid_to_points(std::vector<int> indx, std::vector<anti::Vec3d> &gverts);
  bool tri_test(int i, int j, int di, int dj);
  void grid_to_tris(std::vector<int> indx,
                    std::vector<std::vector<int>> &new_tris,
                    std::vector<std::vector<int>> &orig_edges);
  std::vector<int> make_face_indexes(int i, const std::vector<int> &face);
  int index_map(int i, int j, const std::vector<int> &indx,
                int p_idx = noindex);

  int get_edge_index(int v0, int v1);
  std::vector<int> get_face_indexes(std::vector<int> face);

  IJPos get_pos_xy(int x, int y)
  {
    if (x < 0 || y < 0 || x + y > freq)
      return IJPos::out;
    else
      return IJPos((x == 0) + 2 * (y == 0) + 4 * (x + y == freq));
  }
  IJPos get_pos(int i, int j) { return get_pos_xy(grid_x(i, j), grid_y(i, j)); }

public:
  enum { err_not_tri = 1 };
  Geodesic(const anti::Geometry &base_poly, int mm, int nn = 0, char mthd = 's',
           anti::Vec3d cen = anti::Vec3d(0, 0, 0));
  void make_geo(anti::Geometry &geo);
};

#endif // GEODESIC_H
