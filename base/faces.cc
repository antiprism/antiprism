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
   Name: faces.cc
   Description: face manipulations
   Project: Antiprism - http://www.antiprism.com
*/

#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "geometry.h"
#include "geometryutils.h"

using std::map;
using std::pair;
using std::string;
using std::swap;
using std::vector;

namespace anti {

// close a polyhedron, no more than two open edges per vertex
bool close_poly_basic(Geometry &geom, const Color &col)
{
  map<vector<int>, int> edges;
  map<vector<int>, int>::iterator mi, mi2;
  vector<int> edge(2);
  vector<int> r_edge(2);
  for (unsigned int i = 0; i < geom.faces().size(); ++i)
    for (unsigned int j = 0; j < geom.faces(i).size(); ++j) {
      edge[0] = geom.faces(i, j);
      edge[1] = geom.faces_mod(i, j + 1);
      mi = edges.find(edge);
      if (mi != edges.end())
        mi->second++;
      else {
        r_edge[0] = edge[1];
        r_edge[1] = edge[0];
        mi2 = edges.find(r_edge);
        if (mi2 != edges.end())
          mi2->second++;
        else
          edges[r_edge] = 1; // edge is oriented for missing face
      }
    }

  map<int, vector<int>> neighbours;
  map<int, vector<int>>::iterator ni;
  for (mi = edges.begin(); mi != edges.end(); mi++) {
    if (mi->second == 1) {
      for (int i = 0; i < 2; i++) {
        const int from = mi->first[i];
        const int to = mi->first[(i + 1) % 2];
        ni = neighbours.find(from);
        if (ni == neighbours.end()) {
          vector<int> idxs(2, -1);
          idxs[i] = to;
          neighbours[from] = idxs;
        }
        else {
          if (ni->second[0] < 0)
            ni->second[0] = to;
          else if (ni->second[1] < 0)
            ni->second[1] = to;
          else // three open edges at a vertex
            return false;
        }
      }
    }
  }

  while (neighbours.size()) {
    vector<int> face;
    ni = neighbours.begin();
    int prev = -1;
    int curr = ni->first;
    face.push_back(curr);
    while (true) {
      int new_curr = (ni->second[0] != prev) ? ni->second[0] : ni->second[1];
      prev = curr;
      curr = new_curr;
      neighbours.erase(ni);
      ni = neighbours.find(curr);
      if (ni == neighbours.end()) // completed face
        break;
      face.push_back(curr);
    }
    geom.add_face(face, col);
  }
  return true;
}

bool face_bond(Geometry &base, Geometry &brick, int base_f_idx, int brick_f_idx,
               int off, bool merge, bool flip)
{
  int flip_sign = flip ? 1 : -1;
  vector<Vec3d> pts, bpts;
  for (unsigned int i = 0; i < base.faces(base_f_idx).size() && i < 3; i++) {
    pts.push_back(base.face_v_mod(base_f_idx, i));
    bpts.push_back(brick.face_v_mod(brick_f_idx, off + flip_sign * i));
  }

  brick.transform(Trans3d::align(bpts, pts));

  if (merge) {
    map<int, int> vmap;
    for (unsigned int i = 0; i < brick.faces(brick_f_idx).size(); i++)
      vmap[brick.faces_mod(brick_f_idx, off + flip_sign * i) +
           base.verts().size()] = base.faces_mod(base_f_idx, i);

    vector<int> del_faces(1);
    del_faces[0] = base_f_idx;
    base.del(FACES, del_faces);
    del_faces[0] = brick_f_idx;
    brick.del(FACES, del_faces);
    base.append(brick);
    base.verts_merge(vmap);
  }

  return true;
}

int combine_faces(vector<int> &base, vector<int> brick, const vector<int> &edge)
{
  int v0 = edge[0];
  int v1 = edge[1];
  vector<int> *faces[2] = {&base, &brick};
  // write faces so brick is {edge[1], ...., edge[0]} and
  // base is {edge[0], ... interior vertices ..., edge[1]}
  // then the interior vertices of brick can be appended to base
  // swap edge vertices if they appear in reverse order in base
  int found = 0;
  for (int f = 0; f < 2; f++) {
    vector<int> &face = *faces[f];
    int f_sz = face.size();
    for (int i = 0; i < f_sz; i++) {
      if (face[i] == v0) {
        if (face[(i + 1) % f_sz] == v1) {
          found++;
          std::rotate(face.begin(), face.begin() + (i + 1) % f_sz, face.end());
          if (f == 1)
            std::reverse(face.begin(), face.end());
          break;
        }
        else if (face[(i - 1 + f_sz) % f_sz] == v1) {
          found++;
          std::rotate(face.begin(), face.begin() + i, face.end());
          if (f == 0)
            std::swap(v0, v1);
          break;
        }
      }
    }
  }

  if (found != 2) // edge not found in a face
    return 0;

  if (brick.size() > 2)
    base.insert(base.end(), brick.begin() + 1, brick.end() - 1);
  return 1;
}

// Combine face circuits into a single face by bridging between circuits
// with edges
vector<int> make_face_from_contours(const vector<vector<int>> &contours,
                                    vector<vector<int>> *new_edges)
{
  vector<int> face;
  if (new_edges)
    new_edges->clear();
  vector<int> fs; // faces with two or more sides
  for (unsigned int i = 0; i < contours.size(); i++)
    if (contours[i].size() > 1)
      fs.push_back(i);

  for (unsigned int i = 0; i < fs.size(); i++) {
    const int idx = fs[i];
    if (i > 0) {
      const int prev_idx = fs[i - 1];
      face.push_back(contours[prev_idx][0]); // first vert on prev circuit
      if (new_edges)
        new_edges->push_back(
            make_edge(contours[prev_idx][0], contours[idx][0]));
    }
    for (unsigned int j = 0; j < contours[idx].size(); j++)
      face.push_back(contours[i][j]);
  }
  for (int i = (int)fs.size() - 1; i > 0; i--)
    face.push_back(contours[fs[i]][0]);

  return face;
}

} // namespace anti
