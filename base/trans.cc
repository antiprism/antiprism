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
   Name: trans.cc
   Description: polyhedron transformations
   Project: Antiprism - http://www.antiprism.com
*/

#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "geometry.h"
#include "geometryinfo.h"
#include "geometryutils.h"

using std::map;
using std::pair;
using std::string;
using std::swap;
using std::vector;

namespace anti {
// truncate vertices specified by index number
void truncate_verts(Geometry &geom, vector<int> &v_idxs, double ratio,
                    GeometryInfo *info)
{
  GeometryInfo tmp_info(geom);
  if (!info)
    info = &tmp_info;

  Geometry trunc_geom;

  // Convert v_idxs to sorted and without duplicates
  std::sort(v_idxs.begin(), v_idxs.end());
  v_idxs.erase(unique(v_idxs.begin(), v_idxs.end()), v_idxs.end());

  // Add the non truncated vertices to the geometry. Map old to new
  // vertex index numbers
  map<int, int> old2new_idxs;
  { // new scope for tempory objects
    vector<int> v_range(geom.verts().size());
    for (int i = 0; i < (int)geom.verts().size(); i++)
      v_range[i] = i;
    vector<int> non_trunc_idxs;
    std::set_difference(v_range.begin(), v_range.end(), v_idxs.begin(),
                        v_idxs.end(),
                        std::inserter(non_trunc_idxs, non_trunc_idxs.begin()));
    for (int i = 0; i < (int)non_trunc_idxs.size(); i++) {
      int old_idx = non_trunc_idxs[i];
      old2new_idxs[old_idx] = i;
      trunc_geom.add_vert(geom.verts(old_idx), geom.colors(VERTS).get(old_idx));
    }
  }

  // Each old edge is mapped, in each order, to a new vertex index associated
  // with the truncation of the first edge vertex. New vertex is an original
  // vertex location or a truncation vertex lying on the edge.
  map<vector<int>, int> old_es2new_idxs;
  for (int idx0 = 0; idx0 < (int)info->get_vert_cons().size(); idx0++) {
    auto non_trunc_it = old2new_idxs.find(idx0);
    bool truncating = non_trunc_it == old2new_idxs.end();
    for (int j = 0; j < (int)info->get_vert_cons()[idx0].size(); j++) {
      int idx1 = info->get_vert_cons()[idx0][j];
      if (truncating) {
        // is it a merged new vertex - both ends of the edge are trunc'd to 0.5
        bool merged_vertex =
            ratio == 0.5 && old2new_idxs.find(idx1) == old2new_idxs.end();
        // Add the new vertex unless it is the second copy of a merged vertex
        if (!merged_vertex || idx0 < idx1) {
          old_es2new_idxs[{idx0, idx1}] = (int)trunc_geom.verts().size();
          trunc_geom.add_vert(geom.verts(idx0) +
                              ratio * (geom.verts(idx1) - geom.verts(idx0)));
        }
        // Add the second copy of a merged vertex
        if (merged_vertex && idx0 < idx1)
          old_es2new_idxs[{idx1, idx0}] = (int)trunc_geom.verts().size() - 1;
      }
      else
        old_es2new_idxs[{idx0, idx1}] = non_trunc_it->second;
    }
  }

  // Add faces from truncated vertices
  for (int idx0 : v_idxs) {
    for (const auto &cons : info->get_vert_figs()[idx0]) {
      vector<int> face;
      for (int idx1 : cons)
        face.push_back(old_es2new_idxs[{idx0, idx1}]);
      if (face.size() >= 2)
        trunc_geom.add_face(face);
    }
  }

  // Add original faces truncated at the vertices
  for (int f_idx = 0; f_idx < (int)geom.faces().size(); f_idx++) {
    vector<int> trunc_face;
    for (int i = 0; i < (int)geom.faces(f_idx).size(); i++) {
      int idx0 = geom.faces(f_idx, i);
      int idx1 = geom.faces_mod(f_idx, i + 1);
      int new_idx0 = old_es2new_idxs[{idx0, idx1}];
      int new_idx1 = old_es2new_idxs[{idx1, idx0}];
      // Add the first vertex if it is not the same as the previous vertex
      if (!trunc_face.size() || new_idx0 != trunc_face.back())
        trunc_face.push_back(new_idx0);
      // Add the second vertex if it is not the same as the previous vertex
      if (new_idx1 != trunc_face.back())
        trunc_face.push_back(new_idx1);
    }
    // Clear any repeated index of the first and last entry
    if (trunc_face.size() > 1 && trunc_face.front() == trunc_face.back())
      trunc_face.resize(trunc_face.size() - 1);
    trunc_geom.add_face(trunc_face, geom.colors(FACES).get(f_idx));
  }

  geom = trunc_geom;
}

// truncate vertices with vertex order (default order = 0, truncate all)
void truncate_verts(Geometry &geom, double ratio, int order, GeometryInfo *info)
{
  GeometryInfo tmp_info(geom);
  if (!info)
    info = &tmp_info;

  vector<int> v_idxs;
  const vector<vector<int>> &vcons = info->get_vert_cons();
  for (unsigned int i = 0; i < vcons.size(); i++)
    if (order == 0 || (int)vcons[i].size() == order)
      v_idxs.push_back(i);
  truncate_verts(geom, v_idxs, ratio, info);
}

/*
/// Make a geometry with a face for each edge in the original
/ **Like the Conway 'join' operation.
  * \param geom geometry with edges to convert into faces.* /
void make_edges_to_faces(Geometry &geom);

void make_edges_to_faces(Geometry &geom)
{
  int orig_v_sz = geom.verts().size();
  for (auto &f : geom.faces())
    geom.add_vert(geom.face_cent(f));

  auto edges = geom.get_edge_face_pairs(false);
  geom.clear(FACES);

  for (auto &mv : edges) {
    if (mv.second.size() == 2) {
      vector<int> face(4);
      face[0] = mv.first[0];
      face[1] = mv.second[0] + orig_v_sz;
      face[2] = mv.first[1];
      face[3] = mv.second[1] + orig_v_sz;
      geom.add_face(face);
    }
  }

  geom.del(VERTS, geom.get_info().get_free_verts());
}
*/

} // namespace anti
