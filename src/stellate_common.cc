/*
   Copyright (c) 2003-2022, Roger Kaufman

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
   Name: stellate_common.cc
   Description: stellate code shared in /src
   Project: Antiprism - http://www.antiprism.com
*/

#include "stellate_common.h"
#include "../base/antiprism.h"

#include <cstdio>
#include <set>
#include <string>
#include <vector>

using std::set;
using std::string;
using std::vector;

using namespace anti;

// use anti::epsilon so it is not effected by opt epsilon
void color_stellation(Geometry &geom, const char face_coloring_method,
                      const char edge_coloring_method,
                      const char vertex_coloring_method,
                      const Color &face_color, const Color &edge_color,
                      const Color &vertex_color, const int face_opacity,
                      const string &map_string, const string caller)
{
  // set color map
  Coloring clrng(&geom);
  ColorMap *cmap = colormap_from_name(map_string.c_str());
  clrng.add_cmap(cmap);

  // in case of face coloring C
  Geometry kis;

  // geom is built with face colors from the diagram
  if (face_coloring_method != 'd') {
    if (!face_coloring_method && !face_color.is_set())
      // if no color specified, clear faces
      geom.colors(FACES).clear();
    else if (face_coloring_method == 's') {
      // color faces by symmetry
      Symmetry sym;
      vector<vector<set<int>>> sym_equivs;
      sym.init(geom, &sym_equivs);
      clrng.f_sets(sym_equivs[2], true);
    }
    else if (face_coloring_method == 'c')
      // compound coloring
      clrng.f_parts(true);
    else if (face_coloring_method == 'C') {
      // color by connection
      wythoff_make_tiling(kis, geom, "k", true, false);
      // remove digons
      vector<int> dels;
      for (unsigned int i = 0; i < kis.faces().size(); i++) {
        if (kis.faces(i).size() < 3)
          dels.push_back((int)i);
      }
      kis.del(FACES, dels);
      kis.orient(1); // positive orientation

      // make new verts and edges invisible
      kis.add_missing_impl_edges();
      for (unsigned int i = 0; i < kis.verts().size(); i++) {
        int v_idx = find_vert_by_coords(geom, kis.verts()[i], anti::epsilon);
        if (v_idx == -1) {
          kis.colors(VERTS).set(i, Color::invisible);
          vector<int> edge_idx = find_edges_with_vertex(kis.edges(), i);
          for (unsigned int j = 0; j < edge_idx.size(); j++)
            kis.colors(EDGES).set(edge_idx[j], Color::invisible);
        }
      }
      // the old faces are cleared and kis faces added
      geom.clear(FACES);
      geom.append(kis);
      int blend_type = 1; // first color, invisible edges stay
      merge_coincident_elements(geom, "vef", blend_type, anti::epsilon);

      for (unsigned int i = 0; i < geom.faces().size(); i++) {
        vector<int> face = geom.faces()[i];
        unsigned int fsz = face.size();
        // face to face
        // connections with invisible faces are ignored
        int connections = 0;
        for (unsigned int j = 0; j < fsz; j++) {
          int v1 = face[j];
          int v2 = face[(j + 1) % fsz];
          vector<int> edge = make_edge(v1, v2);
          vector<int> face_idx = find_faces_with_edge(geom.faces(), edge);
          int edge_no = find_edge_in_edge_list(geom.edges(), edge);
          if (!(geom.colors(EDGES).get(edge_no)).is_invisible())
            connections += face_idx.size();
        }
        geom.colors(FACES).set(i, cmap->get_col(connections));
      }
    }
    else
      // use color selected
      clrng.f_one_col(face_color);
  }

  // collect invisible edges
  vector<int> invisible_edges;
  for (unsigned int i = 0; i < geom.edges().size(); i++)
    if ((geom.colors(EDGES).get(i)).is_invisible())
      invisible_edges.push_back(i);

  // color edges
  if (edge_coloring_method == 'f') {
    // if face colors is none clear edges
    if (!face_coloring_method && !face_color.is_set())
      geom.colors(EDGES).clear();
    else
      // edges take colors from faces
      clrng.e_from_adjacent(FACES);
  }
  else if (edge_coloring_method == 'C') {
    // color by connection
    clrng.e_order(true);
  }
  else
    // use color selected
    clrng.e_one_col(edge_color);

  // color vertices
  if (vertex_coloring_method == 'e') {
    // vertices take color from edges
    clrng.v_from_adjacent(EDGES);
  }
  else if (vertex_coloring_method == 'f') {
    // vertices take color from edges
    clrng.v_from_adjacent(FACES);
  }
  else if (vertex_coloring_method == 'n')
    clrng.v_order(true);
  else
    // use color selected
    clrng.v_one_col(vertex_color);

  // reassert invisible edges
  for (unsigned int i = 0; i < invisible_edges.size(); i++)
    geom.colors(EDGES).set(invisible_edges[i], Color::invisible);

  // if using face connection coloring
  // vertices from kis must be made invisible in geom
  if (face_coloring_method == 'C') {
    for (unsigned int i = 0; i < kis.verts().size(); i++) {
      int v_idx = find_vert_by_coords(geom, kis.verts()[i], anti::epsilon);
      if (v_idx != -1) {
        if ((kis.colors(VERTS).get(i)).is_invisible())
          geom.colors(VERTS).set(v_idx, Color::invisible);
      }
    }
  }

  // set transparency
  if (face_opacity > -1) {
    Status stat = Coloring(&geom).apply_transparency(face_opacity);
    if (stat.is_warning())
      fprintf(stderr, "%s: warning: option -T: %s\n", caller.c_str(),
              stat.msg().c_str());
  }
}
