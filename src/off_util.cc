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
   Name: off_util.cc
   Description: utility processing for OFF file, e.g. merge, orient
   Project: Antiprism - http://www.antiprism.com
*/

#include <algorithm>
#include <ctype.h>
#include <map>
#include <math.h>
#include <memory>
#include <set>
#include <stack>
#include <string.h>
#include <string>
#include <vector>

#include "../base/antiprism.h"

#include "help.h"

using std::string;
using std::vector;
using std::set;
using std::map;
using std::pair;
using std::stack;
using std::unique_ptr;

using namespace anti;

// Defined at end of file
int unzip_poly(Geometry &geom, int root, double fract, char centring,
               bool unzip_z_align, char *errmsg);

void triangulate_faces(Geometry &geom, unsigned int winding_rule)
{
  geom.triangulate(Color::invisible, winding_rule);
}

void make_skeleton(Geometry &geom)
{
  geom.add_missing_impl_edges();
  geom.clear(FACES);
}

// Roger Kaufman
// explicit edge type is part of face, or a free edge (not part of face)
void clear_explicit_edges_by_type(Geometry &geom, bool part_of_face)
{
  vector<vector<int>> implicit_edges;
  geom.get_impl_edges(implicit_edges);

  const vector<vector<int>> &edges = geom.edges();

  vector<int> deleted_edges;

  for (unsigned int i = 0; i < edges.size(); i++) {
    int answer = find_edge_in_edge_list(implicit_edges, edges[i]);
    if ((answer > -1 && part_of_face) || (answer == -1 && !part_of_face))
      deleted_edges.push_back(i);
  }

  geom.del(EDGES, deleted_edges);
}

void delete_free_faces(Geometry &geom)
{
  const vector<vector<int>> &faces = geom.faces();
  int fsz = faces.size();
  vector<bool> found(fsz);

  for (int i = 0; i < fsz; i++) {
    if (found[i])
      continue;
    // need to check for faces with lower index than i
    for (int j = 0; j < fsz; j++) {
      if (i == j)
        continue;
      for (int k = 0; k < (int)faces[i].size(); k++) {
        if (vertex_exists_in_face(faces[j], faces[i][k])) {
          found[i] = true;
          found[j] = true;
          break;
        }
      }
      if (found[i])
        break;
    }
  }

  vector<int> face_list;
  for (int i = 0; i < fsz; i++)
    if (!found[i])
      face_list.push_back(i);

  if (face_list.size())
    geom.del(FACES, face_list);
}

void geometry_only(Geometry &geom)
{
  geom.colors(VERTS).clear();
  geom.colors(EDGES).clear();
  geom.colors(FACES).clear();
  clear_explicit_edges_by_type(geom, true);
}

void filter(Geometry &geom, const char *elems)
{
  for (const char *p = elems; *p; p++) {
    switch (*p) {
    case 'V':
      geom.del(VERTS, geom.get_info().get_free_verts());
      break;
    case 'v':
      geom.colors(VERTS).clear();
      break;
    case 'e':
      geom.clear(EDGES);
      break;
    case 'E':
      clear_explicit_edges_by_type(geom, false);
      break;
    case 'D':
      clear_explicit_edges_by_type(geom, true);
      break;
    case 'f':
      geom.clear(FACES);
      break;
    case 'F':
      delete_free_faces(geom);
      break;
    }
  }
}

// get indexes of edge parts
void get_edge_part(vector<int> &edge_part, const Geometry &geom, const int idx,
                   const vector<vector<int>> &vcons, vector<bool> &seen)
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
vector<vector<int>> get_edge_parts(const Geometry &geom)
{
  vector<vector<int>> vcons(geom.verts().size(), vector<int>());

  const vector<vector<int>> &edges = geom.edges();
  for (const auto &edge : edges) {
    vcons[edge[0]].push_back(edge[1]);
    vcons[edge[1]].push_back(edge[0]);
  }

  vector<vector<int>> edge_parts;

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

  return edge_parts;
}

int max_vertex_order(Geometry &geom)
{
  GeometryInfo info(geom);

  vector<int> v_idxs;
  const vector<vector<int>> &vcons = info.get_vert_cons();

  int max_vertex_order = 0;
  for (const auto &vcon : vcons)
    if ((int)vcon.size() > max_vertex_order)
      max_vertex_order = (int)vcon.size();

  return max_vertex_order;
}

void resolve_vert_order_to_verts(Geometry &geom, vector<int> &elem_list)
{
  GeometryInfo info(geom);

  vector<int> v_idxs;
  const vector<vector<int>> &vcons = info.get_vert_cons();

  vector<int> elem_list_verts;
  for (unsigned int i = 0; i < vcons.size(); i++) {
    for (int j : elem_list)
      if ((int)vcons[i].size() == j) {
        elem_list_verts.push_back(i);
        break;
      }
  }

  elem_list.clear();
  elem_list = elem_list_verts;
}

int max_face_sides(Geometry &geom)
{
  const vector<vector<int>> &faces = geom.faces();

  int max_face_sides = 0;
  for (const auto &face : faces)
    if ((int)face.size() > max_face_sides)
      max_face_sides = (int)face.size();

  return max_face_sides;
}

void resolve_face_sides_to_faces(Geometry &geom, vector<int> &elem_list)
{
  const vector<vector<int>> &faces = geom.faces();

  vector<int> elem_list_faces;
  for (unsigned int i = 0; i < faces.size(); i++) {
    for (int j : elem_list)
      if ((int)faces[i].size() == j) {
        elem_list_faces.push_back(i);
        break;
      }
  }

  elem_list.clear();
  elem_list = elem_list_faces;
}

void resolve_color_to_elem(Geometry &geom, const Color &c,
                           const char selection_type_char,
                           vector<vector<int>> &elem_lists)
{
  int elem_type = 0;

  vector<int> elem_list;
  if (selection_type_char == 'x') {
    elem_type = 0;
    for (unsigned int i = 0; i < geom.verts().size(); i++) {
      if (geom.colors(VERTS).get(i) == c)
        elem_list.push_back(i);
    }
  }
  else if (selection_type_char == 'y') {
    elem_type = 1;
    for (unsigned int i = 0; i < geom.edges().size(); i++) {
      if (geom.colors(EDGES).get(i) == c)
        elem_list.push_back(i);
    }
  }
  else if (selection_type_char == 'z') {
    elem_type = 2;
    for (unsigned int i = 0; i < geom.faces().size(); i++) {
      if (geom.colors(FACES).get(i) == c)
        elem_list.push_back(i);
    }
  }

  elem_lists[elem_type] = elem_list;
  elem_list.clear();
}

bool get_del_element_list(Geometry &geom, const string &elem,
                          vector<vector<int>> &elem_lists,
                          vector<vector<int>> &edge_parts,
                          vector<vector<int>> &face_parts, char *errmsg)
{
  *errmsg = '\0';
  if (!elem.size())
    return true;

  unique_ptr<char> elem_str(copy_str(elem.c_str()));
  if (!elem_str.get()) {
    strcpy_msg(errmsg, "could not allocate memory");
    return false;
  }

  const char *elem_type_strs[] = {"vertex",    "edge",      "face",
                                  "edge part", "face part", "vertex orders",
                                  "face sides"};
  char elem_type_char = *elem_str;
  int elem_type;
  int elem_type_name;
  int elems_sz;
  if (elem_type_char == 'v') {
    elem_type_name = 0;
    elem_type = 0;
    elems_sz = geom.verts().size();
  }
  else if (elem_type_char == 'e') {
    elem_type_name = 1;
    elem_type = 1;
    elems_sz = geom.edges().size();
  }
  else if (elem_type_char == 'f') {
    elem_type_name = 2;
    elem_type = 2;
    elems_sz = geom.faces().size();
  }
  else if (elem_type_char == 'E') {
    elem_type_name = 3;
    elem_type = 1;
    elems_sz = edge_parts.size();
  }
  else if (elem_type_char == 'F') {
    elem_type_name = 4;
    elem_type = 2;
    elems_sz = face_parts.size();
  }
  else if (elem_type_char == 'o') {
    elem_type_name = 5;
    elem_type = 0;
    elems_sz = max_vertex_order(geom) + 1;
  }
  else if (elem_type_char == 's') {
    elem_type_name = 6;
    elem_type = 2;
    elems_sz = max_face_sides(geom) + 1;
  }
  else {
    strcpy_msg(errmsg, msg_str("invalid element type '%c', "
                               "should be v, e, f, E, F, o, s",
                               elem_type_char)
                           .c_str());
    return false;
  }

  if (elem_lists[elem_type].size()) {
    strcpy_msg(errmsg, msg_str("list for %s elements already given",
                               elem_type_strs[elem_type_name])
                           .c_str());
    return false;
  }

  Status stat =
      read_idx_list(elem_str.get() + 1, elem_lists[elem_type], elems_sz, false);
  if (stat.is_error()) {
    strcpy_msg(errmsg, msg_str("list for %s elements: %s",
                               elem_type_strs[elem_type], stat.c_msg())
                           .c_str());
    return false;
  }

  return true;
}

// vertex not used. future use?
bool get_del_parts_list(vector<vector<int>> &edge_parts,
                        vector<vector<int>> &face_parts, const string &elem,
                        vector<vector<int>> &elem_lists, char *errmsg)
{
  *errmsg = '\0';
  if (!elem.size())
    return true;

  std::unique_ptr<char> elem_str(copy_str(elem.c_str()));
  if (!elem_str.get()) {
    strcpy_msg(errmsg, "could not allocate memory");
    return false;
  }

  const char *elem_type_strs[] = {"vertex", "edge", "face"};
  char elem_type_char = *elem_str;
  int elem_type;
  int elems_sz = 0;
  if (elem_type_char == 'V') {
    elem_type = 0;
    // elems_sz = parts[0].size();
  }
  else if (elem_type_char == 'E') {
    elem_type = 1;
    elems_sz = edge_parts.size();
  }
  else if (elem_type_char == 'F') {
    elem_type = 2;
    elems_sz = face_parts.size();
  }
  else {
    strcpy_msg(errmsg, msg_str("invalid element type '%c', "
                               "should be E, or F",
                               elem_type_char)
                           .c_str());
    return false;
  }

  if (elem_lists[elem_type].size()) {
    strcpy_msg(errmsg, msg_str("list for %s parts already given",
                               elem_type_strs[elem_type])
                           .c_str());
    return false;
  }

  Status stat =
      read_idx_list(elem_str.get() + 1, elem_lists[elem_type], elems_sz, false);
  if (stat.is_error()) {
    strcpy_msg(errmsg, msg_str("list for %s parts: %s",
                               elem_type_strs[elem_type], stat.c_msg())
                           .c_str());
    return false;
  }

  return true;
}

// if multi is true, old behavior of more than one element type deleted
// collectively is still possible
bool delete_elements(Geometry &geom, vector<string> del_elems, bool keep,
                     string invert_del, bool multi, char *errmsg)
{
  *errmsg = '\0';
  if ((int)del_elems[0].size() <= 1) {
    strcpy_msg(errmsg, "selection string is missing");
    return false;
  }

  // preserve original edge list and complete a new one
  vector<vector<int>> original_edges = geom.edges();
  geom.add_missing_impl_edges();

  // get parts here
  vector<vector<int>> edge_parts;
  edge_parts = get_edge_parts(geom);

  vector<vector<int>> face_parts;
  geom.orient(&face_parts);
  // end get parts here

  vector<vector<int>> elem_lists(3);
  vector<string>::const_iterator vi;
  for (vi = del_elems.begin(); vi != del_elems.end(); ++vi) {
    string vi_str = *vi;
    char elem_type_char = vi_str[0];
    // fprintf(stderr,"elem_type_char = %c\n",elem_type_char);
    string str = "vefEF";
    std::size_t found = str.find(elem_type_char);
    if (found == std::string::npos) {
      strcpy_msg(errmsg, msg_str("invalid element type '%c', "
                                 "should be v, e, f, E , F",
                                 elem_type_char)
                             .c_str());
      return false;
    }

    char selection_type_char = vi_str[1];
    str = "vefEFosxyz";
    found = str.find(selection_type_char);
    if (found == std::string::npos) {
      str = "0123456789-";
      found = str.find(selection_type_char);
      if (found != std::string::npos) {
        selection_type_char = elem_type_char;
      }
      else {
        strcpy_msg(errmsg, msg_str("invalid selection type '%c', "
                                   "should be v, e, f, E, F, o, s, x, y, z",
                                   selection_type_char)
                               .c_str());
        return false;
      }
    }
    else {
      vi_str = vi_str.substr(1);
    }
    // fprintf(stderr,"selection_type_char = %c\n",selection_type_char);
    // fprintf(stderr,"vi_str = %s\n",vi_str.c_str());

    vector<vector<int>> elem_lists_selection(3);

    bool special_selector = false;
    // validate and resolve vertex orders to vertices
    if (selection_type_char == 'o') {
      if (!get_del_element_list(geom, vi_str, elem_lists_selection, edge_parts,
                                face_parts, errmsg))
        return false;
      resolve_vert_order_to_verts(geom, elem_lists_selection[0]);
      // the list is now vertex indexes
      selection_type_char = 'v';
      special_selector = true;
    }
    else
        // validate and resolve face sides to faces
        if (selection_type_char == 's') {
      if (!get_del_element_list(geom, vi_str, elem_lists_selection, edge_parts,
                                face_parts, errmsg))
        return false;
      resolve_face_sides_to_faces(geom, elem_lists_selection[2]);
      // the list is now face indexes
      selection_type_char = 'f';
      special_selector = true;
    }
    else
        // validate and resolve colors to elements
        if (selection_type_char == 'x' || selection_type_char == 'y' ||
            selection_type_char == 'z') {
      // Coloring clrngs[3];
      // if(!read_Colorings(clrngs, vi_str.substr(1), errmsg))
      //   return false;
      Color c;
      char c_name[MSG_SZ];
      if (vi_str.length() > 1)
        strncpy(c_name, vi_str.substr(1).c_str(), MSG_SZ);
      else {
        strcpy_msg(errmsg, "no color specified");
        return false;
      }
      if (!c.read(c_name))
        return false;
      resolve_color_to_elem(geom, c, selection_type_char, elem_lists_selection);
      // the list is now element type
      selection_type_char = (selection_type_char == 'x')
                                ? 'v'
                                : ((selection_type_char == 'y') ? 'e' : 'f');
      special_selector = true;
    }

    if (elem_type_char == 'E' || elem_type_char == 'F') {
      // process parts
      if (elem_type_char == selection_type_char) {
        // resolve part numbers
        if (!get_del_parts_list(edge_parts, face_parts, vi_str,
                                elem_lists_selection, errmsg))
          return false;

        // at this point elem_lists contain part numbers
        vector<int> elem_lists_resolved;
        if (elem_type_char == 'E') {
          for (unsigned int i = 0; i < elem_lists_selection[1].size(); i++) {
            int j = elem_lists_selection[1][i];
            for (int &k : edge_parts[j])
              elem_lists_resolved.push_back(k);
          }
          elem_lists[1].insert(elem_lists[1].end(), elem_lists_resolved.begin(),
                               elem_lists_resolved.end());
        }
        else if (elem_type_char == 'F') {
          for (unsigned int i = 0; i < elem_lists_selection[2].size(); i++) {
            int j = elem_lists_selection[2][i];
            for (int &k : face_parts[j])
              elem_lists_resolved.push_back(k);
          }
          elem_lists[2].insert(elem_lists[2].end(), elem_lists_resolved.begin(),
                               elem_lists_resolved.end());
        }
      }
      // there is a selection type
      else {
        // resolve parts to elements
        if (!special_selector)
          if (!get_del_element_list(geom, vi_str, elem_lists_selection,
                                    edge_parts, face_parts, errmsg))
            return false;

        vector<int> elem_lists_resolved;
        if (elem_type_char == 'E') {
          vector<bool> edge_part_found(edge_parts.size(), false);

          if (selection_type_char == 'v') {
            // add elem_lists for edge parts based on vertex indexes
            for (unsigned int i = 0; i < elem_lists_selection[0].size(); i++) {
              for (unsigned int j = 0; j < edge_parts.size(); j++) {
                if (edge_part_found[j])
                  continue;
                for (unsigned int k = 0; k < edge_parts[j].size(); k++) {
                  if (vertex_exists_in_edge(geom.edges(edge_parts[j][k]),
                                            elem_lists_selection[0][i])) {
                    edge_part_found[j] = true;
                    break;
                  }
                }
              }
            }
          }
          else if (selection_type_char == 'e') {
            // rebuild elem_lists based on edge parts
            for (unsigned int i = 0; i < elem_lists_selection[1].size(); i++) {
              for (unsigned int j = 0; j < edge_parts.size(); j++) {
                if (edge_part_found[j])
                  continue;
                for (unsigned int k = 0; k < edge_parts[j].size(); k++) {
                  if (edge_parts[j][k] == elem_lists_selection[1][i]) {
                    edge_part_found[j] = true;
                    break;
                  }
                }
              }
            }
          }
          else if (selection_type_char == 'f' || selection_type_char == 'F') {
            // if edges selected by face parts
            if (selection_type_char == 'F') {
              vector<int> elem_lists_resolved;
              for (unsigned int i = 0; i < elem_lists_selection[2].size();
                   i++) {
                int j = elem_lists_selection[2][i];
                for (int &k : face_parts[j])
                  elem_lists_resolved.push_back(k);
              }
              elem_lists_selection[2] = elem_lists_resolved;
            }
            // add elem_lists for edge parts based on face indexes
            for (unsigned int i = 0; i < elem_lists_selection[2].size(); i++) {
              for (unsigned int j = 0; j < edge_parts.size(); j++) {
                if (edge_part_found[j])
                  continue;
                for (unsigned int k = 0; k < edge_parts[j].size(); k++) {
                  if (edge_exists_in_face(
                          geom.faces(elem_lists_selection[2][i]),
                          geom.edges(edge_parts[j][k]))) {
                    edge_part_found[j] = true;
                    break;
                  }
                }
              }
            }
          }

          for (unsigned int i = 0; i < edge_part_found.size(); i++)
            if (edge_part_found[i])
              elem_lists[1].insert(elem_lists[1].end(), edge_parts[i].begin(),
                                   edge_parts[i].end());
        }
        else if (elem_type_char == 'F') {
          vector<bool> face_part_found(face_parts.size(), false);

          if (selection_type_char == 'v') {
            // add elem_lists for face parts based on vertex indexes
            for (unsigned int i = 0; i < elem_lists_selection[0].size(); i++) {
              for (unsigned int j = 0; j < face_parts.size(); j++) {
                if (face_part_found[j])
                  continue;
                for (unsigned int k = 0; k < face_parts[j].size(); k++) {
                  if (vertex_exists_in_face(geom.faces(face_parts[j][k]),
                                            elem_lists_selection[0][i])) {
                    face_part_found[j] = true;
                    break;
                  }
                }
              }
            }
          }
          else if (selection_type_char == 'e' || selection_type_char == 'E') {
            // if faces selected by edge parts
            if (selection_type_char == 'E') {
              vector<int> elem_lists_resolved;
              for (unsigned int i = 0; i < elem_lists_selection[1].size();
                   i++) {
                int j = elem_lists_selection[1][i];
                for (int &k : edge_parts[j])
                  elem_lists_resolved.push_back(k);
              }
              elem_lists_selection[1] = elem_lists_resolved;
            }
            // add elem_lists for face parts based on edge indexes
            for (unsigned int i = 0; i < elem_lists_selection[1].size(); i++) {
              for (unsigned int j = 0; j < face_parts.size(); j++) {
                if (face_part_found[j])
                  continue;
                for (unsigned int k = 0; k < face_parts[j].size(); k++) {
                  if (edge_exists_in_face(
                          geom.faces(face_parts[j][k]),
                          geom.edges(elem_lists_selection[1][i]))) {
                    face_part_found[j] = true;
                    break;
                  }
                }
              }
            }
          }
          else if (selection_type_char == 'f') {
            // rebuild elem_lists based on face parts
            for (unsigned int i = 0; i < elem_lists_selection[2].size(); i++) {
              for (unsigned int j = 0; j < face_parts.size(); j++) {
                if (face_part_found[j])
                  continue;
                for (unsigned int k = 0; k < face_parts[j].size(); k++) {
                  if (face_parts[j][k] == elem_lists_selection[2][i]) {
                    face_part_found[j] = true;
                    break;
                  }
                }
              }
            }
          }

          for (unsigned int i = 0; i < face_part_found.size(); i++)
            if (face_part_found[i])
              elem_lists[2].insert(elem_lists[2].end(), face_parts[i].begin(),
                                   face_parts[i].end());
        }
      }
    }
    // process elements
    else {
      if (!special_selector)
        if (!get_del_element_list(geom, vi_str, elem_lists_selection,
                                  edge_parts, face_parts, errmsg))
          return false;

      // previous behaviour
      if (elem_type_char == selection_type_char) {
        for (unsigned int i = 0; i < 3; i++)
          elem_lists[i].insert(elem_lists[i].end(),
                               elem_lists_selection[i].begin(),
                               elem_lists_selection[i].end());
      }
      // when different selection type is specified
      else {
        if (elem_type_char == 'v') {
          if (selection_type_char == 'e' || selection_type_char == 'E') {
            // if vertices selected by edge parts
            if (selection_type_char == 'E') {
              vector<int> elem_lists_resolved;
              for (unsigned int i = 0; i < elem_lists_selection[1].size();
                   i++) {
                int j = elem_lists_selection[1][i];
                for (int &k : edge_parts[j])
                  elem_lists_resolved.push_back(k);
              }
              elem_lists_selection[1] = elem_lists_resolved;
            }
            for (unsigned int i = 0; i < elem_lists_selection[1].size(); i++) {
              vector<int> edge = geom.edges(elem_lists_selection[1][i]);
              elem_lists[0].push_back(edge[0]);
              elem_lists[0].push_back(edge[1]);
            }
          }
          else if (selection_type_char == 'f' || selection_type_char == 'F') {
            // if vertices selected by face parts
            if (selection_type_char == 'F') {
              vector<int> elem_lists_resolved;
              for (unsigned int i = 0; i < elem_lists_selection[2].size();
                   i++) {
                int j = elem_lists_selection[2][i];
                for (int &k : face_parts[j])
                  elem_lists_resolved.push_back(k);
              }
              elem_lists_selection[2] = elem_lists_resolved;
            }
            for (unsigned int i = 0; i < elem_lists_selection[2].size(); i++) {
              vector<int> face = geom.faces(elem_lists_selection[2][i]);
              for (int &j : face)
                elem_lists[0].push_back(j);
            }
          }
        }
        else if (elem_type_char == 'e') {
          if (selection_type_char == 'v') {
            for (unsigned int i = 0; i < elem_lists_selection[0].size(); i++) {
              for (unsigned int j = 0; j < geom.edges().size(); j++) {
                if (vertex_exists_in_edge(geom.edges(j),
                                          elem_lists_selection[0][i]))
                  elem_lists[1].push_back(j);
              }
            }
          }
          else if (selection_type_char == 'E') {
            // if edges selected by edge parts
            vector<int> elem_lists_resolved;
            for (unsigned int i = 0; i < elem_lists_selection[1].size(); i++) {
              int j = elem_lists_selection[1][i];
              for (int &k : edge_parts[j])
                elem_lists_resolved.push_back(k);
            }
            elem_lists_selection[1] = elem_lists_resolved;
            for (unsigned int i = 0; i < elem_lists_selection[1].size(); i++) {
              unsigned int edge_idx1 = elem_lists_selection[1][i];
              for (unsigned int j = 0; j < geom.edges().size(); j++) {
                if (j == edge_idx1)
                  elem_lists[1].push_back(j);
              }
            }
          }
          else if (selection_type_char == 'f' || selection_type_char == 'F') {
            // if edges selected by face parts
            if (selection_type_char == 'F') {
              vector<int> elem_lists_resolved;
              for (unsigned int i = 0; i < elem_lists_selection[2].size();
                   i++) {
                int j = elem_lists_selection[2][i];
                for (int &k : face_parts[j])
                  elem_lists_resolved.push_back(k);
              }
              elem_lists_selection[2] = elem_lists_resolved;
            }
            for (unsigned int i = 0; i < elem_lists_selection[2].size(); i++) {
              vector<int> face = geom.faces(elem_lists_selection[2][i]);
              for (unsigned int j = 0; j < geom.edges().size(); j++) {
                vector<int> edge = geom.edges(j);
                if (edge_exists_in_face(face, edge))
                  elem_lists[1].push_back(j);
              }
            }
          }
        }
        else if (elem_type_char == 'f') {
          if (selection_type_char == 'v') {
            for (unsigned int i = 0; i < elem_lists_selection[0].size(); i++) {
              for (unsigned int j = 0; j < geom.faces().size(); j++) {
                if (vertex_exists_in_edge(geom.faces(j),
                                          elem_lists_selection[0][i]))
                  elem_lists[2].push_back(j);
              }
            }
          }
          else if (selection_type_char == 'e' || selection_type_char == 'E') {
            // if faces selected by edge parts
            if (selection_type_char == 'E') {
              vector<int> elem_lists_resolved;
              for (unsigned int i = 0; i < elem_lists_selection[1].size();
                   i++) {
                int j = elem_lists_selection[1][i];
                for (int &k : edge_parts[j])
                  elem_lists_resolved.push_back(k);
              }
              elem_lists_selection[1] = elem_lists_resolved;
            }
            for (unsigned int i = 0; i < elem_lists_selection[1].size(); i++) {
              vector<int> edge = geom.edges(elem_lists_selection[1][i]);
              for (unsigned int j = 0; j < geom.faces().size(); j++) {
                vector<int> face = geom.faces(j);
                if (edge_exists_in_face(face, edge))
                  elem_lists[2].push_back(j);
              }
            }
          }
          else if (selection_type_char == 'F') {
            // if faces selected by face parts
            vector<int> elem_lists_resolved;
            for (unsigned int i = 0; i < elem_lists_selection[2].size(); i++) {
              int j = elem_lists_selection[2][i];
              for (int &k : face_parts[j])
                elem_lists_resolved.push_back(k);
            }
            elem_lists_selection[2] = elem_lists_resolved;
            for (unsigned int i = 0; i < elem_lists_selection[2].size(); i++) {
              unsigned int face_idx1 = elem_lists_selection[2][i];
              for (unsigned int j = 0; j < geom.faces().size(); j++) {
                if (j == face_idx1)
                  elem_lists[2].push_back(j);
              }
            }
          }
        }
      }
    }
  }

  // keep track if elements were explicitly listed
  vector<bool> explicitly_chosen(3, false);
  for (unsigned int i = 0; i < 3; i++)
    explicitly_chosen[i] = ((int)elem_lists[i].size() ? true : false);

  if (keep) {
    // add in vertices to keep from edges and faces
    for (unsigned int i = 0; i < elem_lists[1].size(); i++) {
      vector<int> edge = geom.edges(elem_lists[1][i]);
      elem_lists[0].push_back(edge[0]);
      elem_lists[0].push_back(edge[1]);
    }

    for (unsigned int i = 0; i < elem_lists[2].size(); i++) {
      vector<int> face = geom.faces(elem_lists[2][i]);
      for (int &j : face)
        elem_lists[0].push_back(j);
    }

    if (!multi) {
      // RK: patch for single element deletion
      // keep edge decorators of faces that are kept
      for (unsigned int i = 0; i < elem_lists[2].size(); i++) {
        vector<int> face = geom.faces(elem_lists[2][i]);
        int sz = face.size();
        for (unsigned int j = 0; j < face.size(); j++) {
          vector<int> edge = make_edge(face[j], face[(j + 1) % sz]);
          elem_lists[1].push_back(find_edge_in_edge_list(geom.edges(), edge));
        }
      }

      // second pass if the face list is inverted
      if (strchr(invert_del.c_str(), 'f')) {
        for (unsigned int i = 0; i < geom.faces().size(); i++) {
          // for any face not found in elem_list[2]
          // remove those edges from elem_list[1]
          if (find(elem_lists[2].begin(), elem_lists[2].end(), i) ==
              elem_lists[2].end()) {
            vector<int> face = geom.faces(i);
            int sz = face.size();
            for (unsigned int j = 0; j < face.size(); j++) {
              vector<int> edge = make_edge(face[j], face[(j + 1) % sz]);
              int answer = find_edge_in_edge_list(geom.edges(), edge);
              elem_lists[1].erase(
                  remove(elem_lists[1].begin(), elem_lists[1].end(), answer),
                  elem_lists[1].end());
            }
          }
        }
      }
    }
  }

  // the lists may contain duplicates and must be unique
  for (unsigned int i = 0; i < 3; i++) {
    sort(elem_lists[i].begin(), elem_lists[i].end());
    auto el = unique(elem_lists[i].begin(), elem_lists[i].end());
    elem_lists[i].resize(el - elem_lists[i].begin());
  }

  vector<bool> invert_deletion(3, false);
  vector<vector<int>> elem_lists_invert(3);

  for (unsigned int i = 0; i < 3; i++) {
    invert_deletion[i] = keep;

    if (multi) {
      if ((i == 0 && strchr(invert_del.c_str(), 'v')) ||
          (i == 1 && strchr(invert_del.c_str(), 'e')) ||
          (i == 2 && strchr(invert_del.c_str(), 'f'))) {
        invert_deletion[i] = !invert_deletion[i];
      }
    }
    else {
      // RK: patch for single element deletion
      // invert deletion for whole model, not just specific element type
      if ((strchr(invert_del.c_str(), 'v')) ||
          (strchr(invert_del.c_str(), 'e')) ||
          (strchr(invert_del.c_str(), 'f'))) {
        invert_deletion[i] = !invert_deletion[i];
      }
    }

    int sz = 0;
    if (i == 0)
      sz = (int)geom.verts().size();
    else if (i == 1)
      sz = (int)geom.edges().size();
    else if (i == 2)
      sz = (int)geom.faces().size();

    // if the inverse is wanted, make a list of elements that were not chosen
    if (invert_deletion[i]) {
      for (int j = 0; j < sz; j++) {
        if (find(elem_lists[i].begin(), elem_lists[i].end(), j) ==
            elem_lists[i].end())
          elem_lists_invert[i].push_back(j);
      }
      elem_lists[i] = elem_lists_invert[i];
    }

    // implicitly selecting all elements?
    if (!explicitly_chosen[i]) {
      // if list is empty and -K then delete all of that element
      // RK: 'multi' is a patch for single element deletion
      if (!elem_lists[i].size() && keep && multi) {
        for (int j = 0; j < sz; j++)
          elem_lists[i].push_back(j);
      }
      else
          // if list is full and -D then don't delete any of that element
          if ((int)elem_lists[i].size() == sz && !keep) {
        elem_lists[i].clear();
      }
    }
  }

  // also delete added edges
  bool explicit_edge_warning = false;
  if (!keep) {
    for (unsigned int i = 0; i < geom.edges().size(); i++) {
      int answer = find_edge_in_edge_list(original_edges, geom.edges(i));
      if (answer == -1) {
        if (find(elem_lists[1].begin(), elem_lists[1].end(), i) !=
            elem_lists[1].end())
          explicit_edge_warning = true;
        elem_lists[1].push_back(i);
      }
    }

    sort(elem_lists[1].begin(), elem_lists[1].end());
    auto el = unique(elem_lists[1].begin(), elem_lists[1].end());
    elem_lists[1].resize(el - elem_lists[1].begin());
  }

  geom.del(FACES, elem_lists[2]);
  geom.del(EDGES, elem_lists[1]);
  if (keep) {
    if (explicitly_chosen[0])
      geom.del(VERTS, elem_lists[0]);
    else
      geom.del(VERTS, geom.get_info().get_free_verts());
  }
  else
    geom.del(VERTS, elem_lists[0]);

  if (explicit_edge_warning)
    strcpy_msg(
        errmsg,
        "some deleted edges were implicit edges. use -e to make them explicit");

  return true;
}

bool add_element(Geometry &geom, const string &elem, char *errmsg)
{
  *errmsg = '\0';
  if (!elem.size())
    return true;

  std::unique_ptr<char> elem_str(copy_str(elem.c_str()));
  if (!elem_str.get()) {
    strcpy_msg(errmsg, "could not allocate memory");
    return false;
  }

  char elem_type_char = *elem_str;
  const char *elem_type_str;
  if (elem_type_char == 'v')
    elem_type_str = "vertex";
  else if (elem_type_char == 'e')
    elem_type_str = "edge";
  else if (elem_type_char == 'f')
    elem_type_str = "face";
  else {
    strcpy_msg(errmsg, msg_str("invalid element type '%c', "
                               "should be v, e, or f",
                               elem_type_char)
                           .c_str());
    return false;
  }

  vector<char *> parts;
  int num_parts = split_line(elem_str.get() + 1, parts, ":");
  if (num_parts > 2) {
    strcpy_msg(errmsg, "more than one ':' used");
    return false;
  }
  if (num_parts == 0 || *parts[0] == '\0') {
    strcpy_msg(errmsg, msg_str("'%s': no element data", elem.c_str()).c_str());
    return false;
  }

  Status stat;
  Color col;
  if (num_parts > 1 && *parts[1] != '\0') {
    if (!(stat = col.read(parts[1]))) {
      strcpy_msg(errmsg, msg_str("%s colour '%s': %s", elem_type_str, parts[1],
                                 stat.c_msg())
                             .c_str());
      return false;
    }
  }

  if (elem_type_char == 'v') {
    Vec3d vec;
    if (!(stat = vec.read(parts[0]))) {
      strcpy_msg(errmsg,
                 msg_str("vertex coordinates '%s': %s", parts[0], stat.c_msg())
                     .c_str());
      return false;
    }
    geom.add_vert(vec, col);
  }
  else if (elem_type_char == 'e' || elem_type_char == 'f') {
    vector<int> idx_list;
    if (!(stat = read_int_list(parts[0], idx_list))) {
      strcpy_msg(errmsg, msg_str("%s index numbers: '%s': %s", elem_type_str,
                                 parts[0], stat.c_msg())
                             .c_str());
      return false;
    }
    int list_sz = idx_list.size();
    for (int i = 0; i < list_sz; i++) {
      int idx = idx_list[i];
      if (idx < 0)
        idx += geom.verts().size();
      if (idx < 0 || idx >= (int)geom.verts().size()) {
        strcpy_msg(errmsg, msg_str("%s index numbers: '%s': %d out of range",
                                   elem_type_str, parts[0], idx_list[i])
                               .c_str());
        return false;
      }
      idx_list[i] = idx;
    }
    if (elem_type_char == 'e') {
      if (list_sz < 2) {
        strcpy_msg(errmsg, msg_str("%s index numbers: '%s': "
                                   "must give at least two numbers",
                                   elem_type_str, parts[0])
                               .c_str());
        return false;
      }
      if (list_sz == 2)
        geom.add_edge(idx_list, col);
      else {
        for (int i = 0; i < list_sz; i++)
          geom.add_edge(make_edge(idx_list[i], idx_list[(i + 1) % list_sz]),
                        col);
      }
    }
    if (elem_type_char == 'f') {
      if (!list_sz) {
        strcpy_msg(errmsg, msg_str("%s index numbers: '%s': "
                                   "no index numbers given",
                                   elem_type_str, parts[0])
                               .c_str());
        return false;
      }
      geom.add_face(idx_list, col);
    }
  }
  return true;
}

bool add_elements(Geometry &geom, vector<string> add_elems, char *errmsg)
{
  vector<string>::const_iterator vi;
  for (vi = add_elems.begin(); vi != add_elems.end(); ++vi) {
    if (!add_element(geom, *vi, errmsg))
      return false;
  }
  return true;
}

void close_poly(Geometry &geom, Color col)
{
  int orig_faces_sz = geom.faces().size();
  close_poly_basic(&geom);
  if (col.is_set()) {
    for (unsigned int i = orig_faces_sz; i < geom.faces().size(); i++)
      geom.colors(FACES).set(i, col);
  }
}

class pr_opts : public ProgramOpts {
public:
  Geometry geom;
  int sig_digits;

  string ofile;

  pr_opts() : ProgramOpts("off_util"), sig_digits(DEF_SIG_DGTS) {}
  void process_command_line(int argc, char **argv);
  void usage();
};

// clang-format off
void pr_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] input_files\n"
"\n"
"Read one or more files in OFF format, combine them into a single file and\n"
"process it. Most operations manipulate elements, such as adding, deleting,\n"
"and combining elements, triangulating and orientimg faces, making an edge\n"
"skeleton, and rounding the precision of coordinates. Other miscellaneous\n"
"operations include projection onto a sphere, creating a net, truncating\n"
"vertices, and converting edges to quadrilaterals. Operations are performeded\n"
"in the order they are given on the command line. input_files is the list of\n"
"files to process, which may include 'null' as an empty geometry, or if not\n"
"given the program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"  -M <args> Sort and merge elements whose coordinates are the same to\n"
"            the number of decimal places given by option -l.  args is 1 or 2\n"
"            comma separated parts. The first part is the elements to merge,\n"
"            which can include: v - vertices, e - edges, f - faces,\n"
"            a - all (vef), b - bond (merge 've' and delete any face\n"
"            coincident with another), s - sort without merging.\n"
"            The second part (default 3) is the merge blend color:\n"
"            first=1, last=2, rgb=3, ryb=4\n"
"  -l <lim>  minimum distance for unique vertex locations as negative\n"
"            exponent (default: %d giving %.0e)\n"
"  -O <opt>  orient the faces first (if possible) then adjust for signed\n"
"            volume, value may be: positive, negative, reverse, or flip\n"
"            (which reverses the orientation of the model as it was input)\n"
"  -T <rat>  truncate vertices by cutting edges at a ratio from each vertex,\n"
"            can also be 'rat,num' to truncate only vertices of order num\n"
"  -E        turn edges into (non-planar) faces\n"
"  -s        skeleton, write the face edges and remove the faces\n"
"  -t <disp> triangulate, include face parts according to winding number\n"
"            from: odd, nonzero, positive, negative, triangulate (synonym\n"
"            for nonzero)\n"
"  -g        geometry only, remove all colours, remove all two-vertex faces\n"
"            (edges) that are also a face edge\n"
"  -x <elms> remove OFF face elements. The element string is processed in\n"
"            order and can include v, e, f to remove OFF faces with one\n"
"            vertex (vertices), two-vertices (edges) and three or more\n"
"            vertices (faces), V to remove vertices that are not part\n"
"            of any face or edge, E to remove two-vertex faces (edges)\n"
"            that are not part of any face, D to remove two-vertex faces (edges)\n"
"            that are also a face edge (decorators), F to remove faces that\n"
"            do not share a vertex with another face\n"
"  -e        Fill in missing explicit edges\n"
"  -D <list> delete a list of elements, list starts with element letter,\n"
"            followed by an index number list, given as index ranges separated\n"
"            by commas. range can be one number or two numbers separated by a\n"
"            hyphen (default range numbers: 0 and largest index). Element\n"
"            letters may also be F or E to delete compound parts by part number.\n"
"            Index number list may be preceded by f, e, v, E or F to find \n"
"            elements based on connectivity to another element type or part.\n"
"            Compound parts may also be found by element number. Only elements\n"
"            specifically specified are deleted. list can have a suffix '%%' to\n"
"            invert results. e or E only act on explicit edges. If some explicit\n"
"            edges are missing, use -e to fill them in\n"
"            special connectivity selectors: o - vertices by vertex order\n"
"               s - faces by number of sides, or to select by a color...\n"
"               x - vertex color, y - edge color, z - face color\n"
"  -K <list> keep a list of elements using the same parameters as -D. Only\n"
"            elements specifically specified are kept along with their vertex\n"
"            and edge decorators if present. Implicit edges which are kept are\n"
"            converted to explicit edges\n"
"  -A <elem> add element, elem is element letter (v, e, f), followed by\n"
"            element data, optionally followed by ':' and a colour. Data is\n"
"               v: three comma separated coordinates\n"
"               e: a comma separated list of index numbers, joined as a ring\n"
"               f: a comma separated list of index numbers\n"
"            negative index numbers are relative to the end of the vertex\n"
"            list, last vertex is -1 (useful to refer to added vertices.)\n"
"  -c <col>  close polyhedron, each hole converted to a face with colour col,\n"
"            holes having a vertex with more than two open edges are not filled\n"
"  -S        project onto unit sphere centred at origin\n"
"  -u <args> unfold a polyhedron into a net, takes up to three comma separated\n"
"            values for base face index, dihedral fraction (normally 1.0 to\n"
"            -1.0, default: 0.0 flat), and final option letters: 'f' centre\n"
"            on centroid of face centres, 'z' align base face normal to z_axis.\n"
"  -d <dgts> number of significant digits (default %d) or if negative\n"
"            then the number of digits after the decimal point\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text, int(-log(::epsilon)/log(10) + 0.5), ::epsilon, DEF_SIG_DGTS);
}
// clang-format on

const char *get_help(const char *name)
{
  map<string, const char *> help;
  help["help"] = help_help;
  help["models"] = help_models;
  help["common_polys"] = help_common_polys;
  help["common_polyhedra"] = help_common_polys;
  help["uniform"] = help_uniform;
  help["archimedean"] = help_uniform;
  help["platonic"] = help_uniform;
  help["uniform_duals"] = help_uniform_duals;
  help["ud"] = help_uniform_duals;
  help["catalans"] = help_uniform_duals;
  help["std_polys"] = help_std_polys;
  help["std_polyhedra"] = help_std_polys;
  help["johnson"] = help_johnson;
  help["polygon"] = help_polygon;
  help["pyramid"] = help_polygon;
  help["prism"] = help_polygon;
  help["antiprism"] = help_polygon;
  help["uc"] = help_uniform_compounds;
  help["geo"] = help_geodesic;
  help["geodesic"] = help_geodesic;
  help["sym_models"] = help_sym_models;
  help["color"] = help_color;
  help["colour"] = help_color;
  help["Color"] = help_color_val;
  help["col_names"] = help_color_names;
  help["col_map"] = help_ColorMap;
  help["expressions"] = help_expreval;
  help["expreval"] = help_expreval;
  help["symmetry"] = help_symmetry;
  help["bowers"] = help_bowers;
  help["schwarz"] = help_schwarz;
  char hname[MSG_SZ];
  to_resource_name(hname, name);
  auto mi = help.find(hname);
  if (mi != help.end())
    return mi->second;
  else
    return nullptr;
}

void pr_opts::process_command_line(int argc, char **argv)
{
  char errmsg[MSG_SZ];
  opterr = 0;
  int c;
  vector<char *> parts;
  Geometry adding;
  Color close_col;
  string arg_id;

  vector<pair<char, char *>> args;

  handle_long_opts(argc, argv);

  bool trailing_option_l = false; // to warn if there is an -l after a -M
  while ((c = getopt(argc, argv, ":hH:st:O:d:x:eD:K:A:c:gT:ESM:l:u:o:")) !=
         -1) {
    if (common_opts(c, optopt))
      continue;

    // Check colour parameter is valid on first pass in case it is omitted
    if (c == 'c')
      print_status_or_exit(close_col.read(optarg), c);
    else if (c == 'l')
      trailing_option_l = true;
    else if (c == 'M')
      trailing_option_l = false;
    else if (c == 'H') {
      const char *help_text = get_help(optarg);
      if (help_text) {
        fprintf(stdout, "\n%s\n", help_text);
        exit(0);
      }
      else
        error(msg_str("no help for '%s' (try 'off_util -H help')", optarg), c);
    }

    args.push_back(pair<char, char *>(c, optarg));
  }
  if (trailing_option_l)
    warning("limit ignored as not followed by a merge operation (-M)", 'l');

  vector<string> ifiles;
  while (argc - optind)
    ifiles.push_back(string(argv[optind++]));

  if (!ifiles.size())
    ifiles.push_back("");

  // Append all input files
  for (auto &ifile : ifiles) {
    Geometry geom_arg;
    if (ifile != "null")
      read_or_error(geom_arg, ifile);
    geom.append(geom_arg);
  }

  int sig_compare = INT_MAX;
  vector<string> add_elems;

  for (auto &arg : args) {
    c = arg.first;
    char *optarg = arg.second;
    switch (c) {
    case 'O':
      print_status_or_exit(get_arg_id(optarg, &arg_id,
                                      "positive=1|negative=2|reverse=3|flip=4",
                                      argmatch_add_id_maps),
                           c);
      print_status_or_exit(geom.orient(atoi(arg_id.c_str())), c);
      break;

    case 'T': {
      vector<char *> parts;
      int parts_sz = split_line(optarg, parts, ",");
      if (parts_sz > 2)
        error("truncation must be 'ratio' or 'ratio,vert_order'", c);
      double trunc_ratio;
      if (!read_double(parts[0], &trunc_ratio))
        error("truncation ratio is not a number", c);

      int trunc_v_ord = 0; // truncate all vertices
      if (parts_sz > 1 && !read_int(parts[1], &trunc_v_ord)) {
        error("truncation vertex order is not an integer", c);
        if (trunc_v_ord < 1)
          error("truncation vertex order is not positive", c);
      }

      truncate_verts(&geom, trunc_ratio, trunc_v_ord);
      break;
    }

    case 'E':
      make_edges_to_faces(&geom);
      break;

    case 's':
      make_skeleton(geom);
      break;

    case 't':
      if (!get_arg_id(optarg, &arg_id,
                      "odd|nonzero|positive|negative|triangulate=1"))
        error(msg_str("invalid winding rule '%s'", optarg).c_str(), c);

      triangulate_faces(geom, TESS_WINDING_ODD + atoi(arg_id.c_str()));
      break;

    case 'd':
      print_status_or_exit(read_int(optarg, &sig_digits), c);
      break;

    case 'x':
      if (strspn(optarg, "vefVEFD") != strlen(optarg))
        error(msg_str("elements to hide are %s, can only include "
                      "vefVEFD\n",
                      optarg),
              c);
      filter(geom, optarg);
      break;

    case 'e':
      geom.add_missing_impl_edges();
      break;

    case 'D': {
      bool invert = (optarg[strlen(optarg) - 1] == '%') ? true : false;
      bool forward_error =
          ((strlen(optarg) > 1) && (optarg[strlen(optarg) - 2] == '%')) ? true
                                                                        : false;
      vector<char *> entries;
      split_line(optarg, entries, "%");
      if (!entries.size())
        error("invalid argument", c);
      if (entries.size() > 1 || forward_error)
        error("extra characters found after %", c);
      vector<string> del_elems;
      del_elems.push_back(entries[0]);
      string invert_del;
      if (invert) {
        if (entries[0][0] == 'f' || entries[0][0] == 'F')
          invert_del = 'f';
        else if (entries[0][0] == 'e' || entries[0][0] == 'E')
          invert_del = 'e';
        else if (entries[0][0] == 'v')
          invert_del = 'v';
      }
      if (!delete_elements(geom, del_elems, false, invert_del, false, errmsg))
        error(errmsg, c);
      if (*errmsg)
        warning(errmsg);
      break;
    }

    case 'K': {
      bool invert = (optarg[strlen(optarg) - 1] == '%') ? true : false;
      bool forward_error =
          ((strlen(optarg) > 1) && (optarg[strlen(optarg) - 2] == '%')) ? true
                                                                        : false;
      vector<char *> entries;
      split_line(optarg, entries, "%");
      if (!entries.size())
        error("invalid argument", c);
      if (entries.size() > 1 || forward_error)
        error("extra characters found after %", c);
      vector<string> del_elems;
      del_elems.push_back(entries[0]);
      string invert_del;
      if (invert) {
        if (entries[0][0] == 'f' || entries[0][0] == 'F')
          invert_del = 'f';
        else if (entries[0][0] == 'e' || entries[0][0] == 'E')
          invert_del = 'e';
        else if (entries[0][0] == 'v')
          invert_del = 'v';
      }
      if (!delete_elements(geom, del_elems, true, invert_del, false, errmsg))
        error(errmsg, c);
      break;
    }

    case 'A':
      add_elems.clear();
      add_elems.push_back(optarg);
      if (!add_elements(geom, add_elems, errmsg))
        error(errmsg, c);
      break;

    case 'g':
      geometry_only(geom);
      break;

    case 'c':
      print_status_or_exit(close_col.read(optarg), c);
      close_poly(geom, close_col);
      break;

    case 'S':
      project_onto_sphere(&geom);
      break;

    case 'M': {
      split_line(optarg, parts, ",");

      // Get merge elements
      char elems[MSG_SZ];
      strncpy(elems, parts[0], MSG_SZ);
      if (strspn(elems, "svefab") != strlen(elems))
        error(msg_str("elements to merge are %s must be v, e, f, a, b "
                      "or s\n",
                      elems),
              c);

      if (strchr(elems, 's') && strlen(elems) > 1)
        error("s is for sorting only, cannot be used with v, e, or f", c);

      if (strchr(elems, 'a') && strlen(elems) > 1)
        error("a includes vef, and must be used alone", c);

      if (strchr(elems, 'b') && strlen(elems) > 1)
        error("b includes vef, and must be used alone", c);

      if (strspn(elems, "ef") && !strchr(elems, 'v'))
        warning("without v, some orphan vertices may result", c);

      if (*elems == 'a')
        strcpy_msg(elems, "vef");

      // Get blend type
      int blend_type = 3;
      if (parts.size() > 1) {
        print_status_or_exit(get_arg_id(parts[1], &arg_id,
                                        "first=1|last=2|rgb=3|ryb=4",
                                        argmatch_add_id_maps),
                             c);
        blend_type = atoi(arg_id.c_str());
      }

      // Process
      double epsilon =
          (sig_compare != INT_MAX) ? pow(10, -sig_compare) : ::epsilon;
      if (*elems == 'b') {
        merge_coincident_elements(&geom, "ve", blend_type, epsilon);
        Geometry tmp = geom;
        vector<map<int, set<int>>> equiv_elems;
        check_congruence(geom, tmp, &equiv_elems, epsilon);
        vector<int> del_faces;
        map<int, set<int>>::iterator mi;
        for (mi = equiv_elems[2].begin(); mi != equiv_elems[2].end(); ++mi) {
          set<int>::iterator si;
          si = mi->second.begin();
          ++si;
          if (si != mi->second.end() && *si < (int)(geom.faces().size())) {
            for (si = mi->second.begin(); si != mi->second.end(); ++si) {
              del_faces.push_back(*si);
            }
          }
        }
        geom.del(FACES, del_faces);
      }
      else
        merge_coincident_elements(&geom, elems, blend_type, epsilon);
      break;
    }

    case 'l':
      print_status_or_exit(read_int(optarg, &sig_compare), c);
      if (sig_compare < 0) {
        warning("limit is negative, and so ignored", c);
      }
      if (sig_compare > DEF_SIG_DGTS) {
        warning("limit is very small, may not be attainable", c);
      }
      break;

    case 'u': {
      int unzip_root;
      split_line(optarg, parts, ",");
      print_status_or_exit(read_int(parts[0], &unzip_root), c);

      double unzip_frac = 0.0;
      if (parts.size() > 1)
        print_status_or_exit(read_double(parts[1], &unzip_frac), c);

      char unzip_centre = 'x';
      char unzip_z_align = false;
      if (parts.size() > 2) {
        if (strspn(parts[2], "zf") != strlen(parts[2]))
          error(msg_str("unzip options are '%s' must include "
                        "only f, u\n",
                        parts[2]),
                c);

        if (strchr(parts[2], 'f'))
          unzip_centre = 'f';
        if (strchr(parts[2], 'z'))
          unzip_z_align = true;
      }

      if (!unzip_poly(geom, unzip_root, unzip_frac, unzip_centre, unzip_z_align,
                      errmsg))
        error(errmsg, 'u');

      break;
    }

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }
}

int main(int argc, char *argv[])
{
  pr_opts opts;
  opts.process_command_line(argc, argv);

  opts.write_or_error(opts.geom, opts.ofile, opts.sig_digits);

  return 0;
}

//------------------------------------------------------------------
// Unzipping and unfolding

struct tree_face {
  int idx;
  int start;
  int cur;
  vector<int> orig_cons;
  vector<int> cons;
  tree_face(int index, vector<int> original_cons, int con_idx = -1);
  int get_next_idx();
};

tree_face::tree_face(int index, vector<int> original_cons, int con_idx)
    : idx(index), orig_cons(original_cons)
{
  if (con_idx == -1) {
    cur = -1;
    start = 0;
  }
  else {
    auto vi = std::find(orig_cons.begin(), orig_cons.end(), con_idx);
    start = cur = vi - orig_cons.begin();
  }
}

int tree_face::get_next_idx()
{
  if (cur - start + 1 >= (int)orig_cons.size()) // all connections seen
    return -1;
  else
    return orig_cons[++cur % orig_cons.size()]; // wrap at end of cons list
}

class unzip_tree {
private:
  int root;
  map<int, vector<int>> tree;

public:
  unzip_tree() : root(-1) {}
  int init_basic(Geometry &geom, int first_face);
  int flatten(const Geometry &geom, Geometry &net_geom, double fract);
};

int unzip_tree::init_basic(Geometry &geom, int first_face)
{
  tree.clear();
  root = first_face;

  Geometry dual;
  get_dual(&dual, geom);
  GeometryInfo info(dual);
  const vector<vector<int>> &f_cons = info.get_vert_cons();
  vector<tree_face> face_list;
  vector<bool> seen(geom.faces().size(), false);

  // start at first vertex of first face
  face_list.push_back(tree_face(root, f_cons[root]));

  int list_sz = -1;
  while (face_list.size() < geom.faces().size()) { // quit when all faces seen
    // while((int)face_list.size()>list_sz) {  // quit when all faces seen
    list_sz = face_list.size();
    for (int i = 0; i < list_sz; i++) { // faces currently in tree
      int idx;
      while ((idx = face_list[i].get_next_idx()) >= 0) { // face conections
        if (!seen[idx]) {                                // unseen face
          seen[idx] = true;
          face_list.push_back(tree_face(idx, f_cons[idx], face_list[i].idx));
          face_list[i].cons.push_back(idx);
          // geom.set_f_col(idx, level);
          continue;
        }
      }
    }
  }

  for (auto &i : face_list)
    if (i.cons.size())
      tree[i.idx] = i.cons;

  return 1;
}

int unzip_tree::flatten(const Geometry &geom, Geometry &net_geom,
                        double fract = 0.0)
{
  net_geom = geom;
  int level = 0;
  int cur_face = root;
  stack<int> face_stack;
  face_stack.push(cur_face);
  unsigned int cur_con_num = 0;
  stack<int> con_num_stack;
  con_num_stack.push(cur_con_num);
  map<int, int> cur_idx_map;
  stack<map<int, int>> idx_map_stack;
  idx_map_stack.push(cur_idx_map);
  stack<Trans3d> trans_stack;
  trans_stack.push(Trans3d());
  stack<Vec3d> norm_stack;

  while (level >= 0) {
    // fprintf(stderr, "new loop: L=%d, cur_con_num=%d, cur_face=%d, prev=%d\n",
    // level, cur_con_num, cur_face, level==0 ? -1 : face_stack.top());

    if (cur_con_num == 0) { // process new face
      // fprintf(stderr, "\tprocess face %d\n", cur_face);
      vector<int> &face = net_geom.raw_faces()[cur_face];
      map<int, int> new_idx_map;
      vector<int> join_edge;
      int first_mapped_i = -1;
      // Trans3d inv_prev_trans = Trans3d::inverse(trans_stack.top());
      for (unsigned int i = 0; i < face.size(); i++) {
        int idx = face[i];
        int new_idx;
        if (cur_face == root) { // new vertices to help deletion later
          net_geom.add_vert(net_geom.verts(idx));
          new_idx = net_geom.verts().size() - 1;
        }
        else {
          // Find vertices common to current and previous faces, and use
          // duplicated vertices for the non-join vertices in the cur face.
          map<int, int>::const_iterator mi = idx_map_stack.top().find(face[i]);
          if (mi == idx_map_stack.top().end()) {
            net_geom.add_vert(trans_stack.top() * net_geom.verts(idx));
            new_idx = net_geom.verts().size() - 1;
          }
          else {
            new_idx = mi->second;
            join_edge.push_back(new_idx);
            // if the join edges don't follow each other they
            // are in the first and last positions, and so are
            // sequential in reverse order
            if (join_edge.size() == 1)
              first_mapped_i = i;
            else if (join_edge.size() == 2) {
              if (i - first_mapped_i > 1)
                std::swap(join_edge[0], join_edge[1]);
            }
          }
        }

        face[i] = new_idx;
        new_idx_map[idx] = new_idx;
      }
      idx_map_stack.push(new_idx_map);

      Trans3d trans;
      Vec3d norm = net_geom.face_norm(cur_face);
      if (cur_face != root) {
        // set join_edge so it is in order for previous face
        const vector<int> &prev_face = net_geom.faces(face_stack.top());
        unsigned int j;
        for (j = 0; j < prev_face.size(); j++)
          if (prev_face[j] == join_edge[0])
            break;
        if (prev_face[(j + 1) % prev_face.size()] == join_edge[1]) {
          std::reverse(face.begin(), face.end());
          norm *= -1.0;
        }
        else
          std::swap(join_edge[0], join_edge[1]);

        Vec3d axis = net_geom.edge_vec(join_edge).unit();
        double ang = angle_around_axis(norm_stack.top(), norm, axis);
        ang *= (fract - 1);
        Trans3d rot = Trans3d::rotate(axis, ang);
        Vec3d offset = net_geom.verts(join_edge[0]);
        trans = Trans3d::translate(offset) * rot * Trans3d::translate(-offset);

        for (int i : face) {
          if (i != join_edge[0] && i != join_edge[1])
            net_geom.raw_verts()[i] = trans * net_geom.verts(i);
        }
        norm = rot * norm; // normal needs rotating because face rotated
      }
      norm_stack.push(norm);
      trans_stack.push(trans * trans_stack.top());
    }

    if (cur_con_num == tree[cur_face].size()) { // finish at this level
      // fprintf(stderr, "\ttree.size()=%d, tree[%d].size()=%d\n", tree.size(),
      // cur_face, tree[cur_face].size());
      // fprintf(stderr, "\tfinish at level\n");
      if (level) {
        cur_face = face_stack.top();
        face_stack.pop();
        cur_con_num = con_num_stack.top() + 1; // reset and increment
        con_num_stack.pop();
        idx_map_stack.pop();
        norm_stack.pop();
        trans_stack.pop();
      }
      level--;
      continue;
    }

    unsigned int next_idx = tree[cur_face][cur_con_num];
    // fprintf(stderr, "\tgo out a level\n");
    level++;
    face_stack.push(cur_face);
    cur_face = next_idx;
    con_num_stack.push(cur_con_num);
    cur_con_num = 0;
  }

  vector<int> del_verts(geom.verts().size());
  for (unsigned int i = 0; i < geom.verts().size(); i++)
    del_verts[i] = i;
  net_geom.del(VERTS, del_verts);
  return 1;
}

int unzip_poly(Geometry &geom, int root, double fract, char centring,
               bool unzip_z_align, char *errmsg)
{
  GeometryInfo info(geom);
  if (!info.is_polyhedron()) {
    strcpy_msg(errmsg, "input not a polyhedron (temporary restriction)");
    return 0;
  }
  if (info.num_parts() > 1) {
    strcpy_msg(errmsg, "input not connected (temporary restriction)");
    return 0;
  }
  if (root < 0 || root >= (int)geom.faces().size()) {
    sprintf(errmsg, "root face '%d' is not a valid face index number", root);
    return 0;
  }
  if (!strchr("fx", centring)) {
    sprintf(errmsg, "invalid centring type '%c', must be f or x", centring);
    return 0;
  }

  if (unzip_z_align)
    geom.transform(Trans3d::rotate(geom.face_norm(root), Vec3d::Z));

  unzip_tree tree;
  tree.init_basic(geom, root);
  Geometry net_geom;
  tree.flatten(geom, net_geom, fract);

  if (centring == 'f') {
    vector<Vec3d> f_cents;
    net_geom.face_cents(f_cents);
    net_geom.transform(Trans3d::translate(-centroid(f_cents)));
  }

  geom = net_geom;
  return 1;
}
