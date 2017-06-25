/*
   Copyright (c) 2017, Roger Kaufman

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
   Name: planar.cc
   Description: shared code originally from src/planar.cc
   Project: Antiprism - http://www.antiprism.com
*/

#include <stdio.h>
#include <stdlib.h>

#include <map>
#include <set>
#include <string>
#include <vector>

#include "planar.h"

using std::string;
using std::vector;
using std::set;
using std::sort;
using std::pair;
using std::make_pair;
using std::map;

namespace anti {

// From planar.cc

// find connections from a vertex v_idx
void find_connections(const Geometry &geom, vector<int> &vcons,
                      const int v_idx) {
  const vector<vector<int>> &edges = geom.edges();

  vcons.clear();
  for (const auto &edge : edges) {
    if (edge[0] == v_idx)
      vcons.push_back(edge[1]);
    else if (edge[1] == v_idx)
      vcons.push_back(edge[0]);
  }
}

// idx will be the dimension not included so that 3D -> 2D projection occurs
void project_using_normal(const Vec3d &normal, int &idx, int &sign) {
  vector<pair<double, int>> normal_info(3);

  for (unsigned int i = 0; i < 3; i++) {
    normal_info[i].second = i;
    normal_info[i].first = fabs(normal[i]); // absolute value for sorting
  }

  stable_sort(normal_info.begin(), normal_info.end());

  idx = normal_info[normal_info.size() - 1].second;
  sign = (normal[idx] >= 0.0) ? 1 : -1;
}

// put faces numbers in face_idxs into fgeom
Geometry faces_to_geom(const Geometry &geom, const vector<int> &face_idxs)
{
  Geometry fgeom;
  fgeom.add_verts(geom.verts());
  fgeom.colors(VERTS) = geom.colors(VERTS);
  for (int j : face_idxs) {
    fgeom.add_face(geom.faces()[j], geom.colors(FACES).get(j));
  }
  fgeom.del(VERTS, fgeom.get_info().get_free_verts());
  return fgeom;
}

// if intersection point does not exist, insert a new one. Otherwise only return
// the index of the existing one
int vertex_into_geom(Geometry &geom, const Vec3d &P, Color vcol,
                     const double eps) {
  int v_idx = find_vert_by_coords(geom, P, eps);
  if (v_idx == -1) {
    geom.add_vert(P, vcol);
    v_idx = geom.verts().size() - 1;
  }

  return v_idx;
}

// if edge already exists, do not create another one and return false. return
// true if new edge created
// check if edge1 and edge2 indexes are equal. If so do not allow an edge length
// of 0 to occur
bool edge_into_geom(Geometry &geom, const int v_idx1, const int v_idx2,
                    Color ecol) {
  // no edge length 0 allowed
  if (v_idx1 == v_idx2)
    return false;

  vector<int> new_edge = make_edge(v_idx1, v_idx2);

  int answer = find_edge_in_edge_list(geom.edges(), new_edge);
  if (answer < 0)
    geom.add_edge(new_edge, ecol);

  return ((answer < 0) ? true : false);
}

// input seperate networks of overlapping edges and merge them into one network
bool mesh_edges(Geometry &geom, const double eps) {
  const vector<Vec3d> &verts = geom.verts();
  const vector<vector<int>> &edges = geom.edges();

  // RK - large diameter meshes are like changing precision
  // standardize on mesh size to radius of 1, restore scale at end
  GeometryInfo info(geom);
  info.set_center(geom.centroid());
  double mesh_radius = info.vert_dist_lims().max;
  geom.transform(Trans3d::scale(1 / mesh_radius));

  // remember original sizes as the geom will be changing size
  int vsz = verts.size();
  int esz = edges.size();

  vector<int> deleted_edges;
  vector<pair<pair<int, int>, int>> new_verts;

  // compare only existing edges
  for (int i = 0; i < esz; i++) {
    vector<pair<double, int>> line_intersections;
    for (int j = 0; j < esz; j++) {
      // don't compare to self
      if (i == j)
        continue;

      // see if the new vertex was already created
      int v_idx = -1;
      for (auto &new_vert : new_verts) {
        if (new_vert.first.first == i && new_vert.first.second == j) {
          v_idx = new_vert.second;
          break;
        }
      }

      // if it doesn't already exist, see if it needs to be created
      if (v_idx == -1) {
        // compare segements P0,P1 with Q0,Q1
        Vec3d intersection_point =
            segments_intersection(verts[edges[i][0]], verts[edges[i][1]],
                                  verts[edges[j][0]], verts[edges[j][1]], eps);
        if (intersection_point.is_set()) {
          // find (or create) index of this vertex
          v_idx =
              vertex_into_geom(geom, intersection_point, Color::invisible, eps);
          // don't include existing vertices
          if (v_idx < vsz)
            v_idx = -1;
          else {
            // store index of vert at i,j. Reverse index i,j so it will be found
            // when encountering edges j,i
            pair<pair<int, int>, int> new_vert;
            new_vert.first.first = j;
            new_vert.first.second = i;
            new_vert.second = v_idx;
            new_verts.push_back(new_vert);
          }
        }
      }

      // if ultimately an intersection was found
      if (v_idx != -1)
        // store distance from P0 to intersection and also the intersection
        // vertex index
        line_intersections.push_back(
            make_pair((verts[edges[i][0]] - verts[v_idx]).len(), v_idx));
    }

    if (line_intersections.size()) {
      // edge i will be replaced. mark it for deletion
      deleted_edges.push_back(i);
      // sort based on distance from P0
      sort(line_intersections.begin(), line_intersections.end());
      // create edgelets from P0 through intersection points to P1 (using
      // indexes)
      edge_into_geom(geom, edges[i][0], line_intersections[0].second,
                     Color::invisible);
      for (unsigned int k = 0; k < line_intersections.size() - 1; k++)
        edge_into_geom(geom, line_intersections[k].second,
                       line_intersections[k + 1].second, Color::invisible);
      edge_into_geom(geom,
                     line_intersections[line_intersections.size() - 1].second,
                     edges[i][1], Color::invisible);
    }
  }

  // delete the replaced edges
  geom.del(EDGES, deleted_edges);

  // restore scale
  geom.transform(Trans3d::scale(mesh_radius));

  return (deleted_edges.size() ? true : false);
}

void build_angle_map(const Geometry &geom,
                     map<pair<int, int>, double> &angle_map,
                     const Vec3d &normal, const double eps) {
  const vector<Vec3d> &verts = geom.verts();

  int idx;
  int sign;

  project_using_normal(normal, idx, sign);

  int idx0 = (idx + 1) % 3;
  int idx1 = (idx + 2) % 3;

  for (unsigned int i = 0; i < verts.size(); i++) {
    vector<int> vcons;
    find_connections(geom, vcons, i);
    for (int k : vcons) {
      double y = verts[k][idx1] - verts[i][idx1];
      double x = verts[k][idx0] - verts[i][idx0];
      double angle = rad2deg(atan2(y, x));
      // eliminate small error
      if (double_eq(angle, 0.0, eps))
        angle = 0.0;
      // make all angles positive
      if (angle < 0.0)
        angle += 360;
      angle_map[make_pair(i, k)] = angle;
    }
  }
}

void build_turn_map(const Geometry &geom, map<pair<int, int>, int> &turn_map,
                    map<pair<int, int>, double> &angle_map) {
  const vector<vector<int>> &edges = geom.edges();

  for (const auto &edge : edges) {
    for (unsigned int j = 0; j < 2; j++) {
      vector<int> vcons;
      int a = edge[!j ? 0 : 1];
      int b = edge[!j ? 1 : 0];

      double base_angle = angle_map[make_pair(b, a)];
      find_connections(geom, vcons, b);
      vector<pair<double, int>> angles;
      for (unsigned int k = 0; k < vcons.size(); k++) {
        int c = vcons[k];
        // can't go backwards unless it is the only possible path
        if (c == a && vcons.size() != 1)
          continue;
        double angle = angle_map[make_pair(b, c)];
        if (angle >= 0 && angle < base_angle)
          angle += 360;
        angles.push_back(make_pair(angle, c));
      }
      // the minimum angle is now the next higher angle than base angle. it will
      // be at the top of the sort
      // map that turn for AB -> C
      sort(angles.begin(), angles.end());
      turn_map[make_pair(a, b)] = angles[0].second;
    }
  }
}

void get_first_unvisited_triad(const Geometry &geom, int &first_unvisited,
                               map<pair<int, int>, bool> &visited,
                               map<pair<int, int>, int> &turn_map,
                               vector<int> &face) {
  const vector<vector<int>> &edges = geom.edges();

  face.clear();
  for (unsigned int i = first_unvisited; i < edges.size(); i++) {
    for (unsigned int j = 0; j < 2; j++) {
      int a = edges[i][!j ? 0 : 1];
      int b = edges[i][!j ? 1 : 0];
      if (visited[make_pair(a, b)])
        continue;
      pair<int, int> edge_pair = make_pair(a, b);
      int c = turn_map[edge_pair];
      if (!visited[make_pair(b, c)]) {
        face.push_back(a);
        face.push_back(b);
        face.push_back(c);
        first_unvisited = i;
        // if j is 1 this edge has been traversed in both directions
        if (j)
          first_unvisited++;
        return;
      }
    }
  }
  return;
}

void construct_faces(Geometry &geom, map<pair<int, int>, int> &turn_map) {
  map<pair<int, int>, bool> visited;

  int first_unvisited = 0;
  vector<int> face;
  get_first_unvisited_triad(geom, first_unvisited, visited, turn_map, face);
  while (face.size()) {
    int fsz = face.size();
    int a = face[fsz - 2];
    int b = face[fsz - 1];
    int c = turn_map[make_pair(a, b)];
    // if c == face[0], still consider face[0] is part of complex polygon if
    // next turn does not lead to face[1]
    if (c != face[0] || turn_map[make_pair(b, c)] != face[1])
      face.push_back(c);
    else {
      // RK: patch. Do not allow faces to have sequential duplicate indexes
      auto fi = unique(face.begin(), face.end());
      face.resize(fi - face.begin());

      geom.add_face(face);
      // mark turns visited and start a new face
      fsz = face.size();
      for (int i = 0; i < fsz; i++)
        visited[make_pair(face[i], face[(i + 1) % fsz])] = true;
      get_first_unvisited_triad(geom, first_unvisited, visited, turn_map, face);
    }
  }
}

void analyze_faces(Geometry &geom, const int planar_merge_type,
                   vector<int> &nonconvex_faces,
                   map<pair<int, int>, double> &angle_map) {
  const vector<vector<int>> &faces = geom.faces();
  vector<int> deleted_faces;
  int deleted_faces_count = 0; // should be only one

  for (unsigned int i = 0; i < faces.size(); i++) {
    double angle_sum = 0;
    bool all_negative_turns = true;
    int fsz = faces[i].size();
    for (int j = 0; j < fsz; j++) {
      int a = faces[i][j];
      int v = faces[i][(j + 1) % fsz];
      int b = faces[i][(j + 2) % fsz];
      double angle =
          angle_map[make_pair(v, b)] - angle_map[make_pair(v, a)] - 180.0;
      if (angle < -180.0)
        angle += 360.0;
      else if (angle > 180.0)
        angle -= 360.0;
      if (angle > 0.0)
        all_negative_turns = false;
      angle_sum += angle;
    }
    int sum = (int)floorf(angle_sum + 0.5);
    // old statement: if ((sum == 360 && planar_merge_type == 1) || (sum != 360
    // && planar_merge_type == 2)) {
    // ideally the interior angle sum of the tiles is -360 degrees. However it
    // is known to be otherwise. Check less than 0
    // the angle sum of the perimeter face is ideally 360. However it is known
    // to sometimes be otherwise. Check greater than or equal to 0
    if ((sum >= 0 && planar_merge_type == 1) ||
        (sum < 0 && planar_merge_type == 2)) {
      deleted_faces.push_back(i);
      deleted_faces_count++;
    } else if (!all_negative_turns)
      nonconvex_faces.push_back(i - deleted_faces_count);
  }
  geom.del(FACES, deleted_faces);
}

void fill_in_faces(Geometry &geom, const int planar_merge_type,
                   vector<int> &nonconvex_faces, const Vec3d &normal,
                   double eps) {
  map<pair<int, int>, double> angle_map;
  map<pair<int, int>, int> turn_map;
  build_angle_map(geom, angle_map, normal, eps);
  build_turn_map(geom, turn_map, angle_map);
  construct_faces(geom, turn_map);
  analyze_faces(geom, planar_merge_type, nonconvex_faces, angle_map);
}

// start of stellate stuff

/* for future use
// http://astronomy.swin.edu.au/~pbourke/geometry/planes/
bool three_plane_intersect(Vec3d Q0, Vec3d n0, Vec3d Q1, Vec3d n1,
                           Vec3d Q2, Vec3d n2, Vec3d &P, double eps)
{
  double tri_prod = vtriple(n0, n1, n2);
  if (fabs(tri_prod) / (n0.len() * n1.len() * n2.len()) < eps)
    return false;
  double d1 = vdot(Q0, n0);
  double d2 = vdot(Q1, n1);
  double d3 = vdot(Q2, n2);
  P = (d1 * vcross(n1, n2) + d2 * vcross(n2, n0) + d3 * vcross(n0, n1)) /
      tri_prod;
  return true;
}

vector<int> neighbor_faces(Geometry &geom, int f_idx)
{
  vector<int> n_faces;

  vector<int> face = geom.faces(f_idx);
  unsigned int sz = face.size();
  for (unsigned int i = 0; i < sz; i++) {
    int v0 = face[i];
    int v1 = face[(i + 1) % sz];
    vector<int> face_idx = find_faces_with_edge(geom.faces(), make_edge(v0,
v1));
    for (unsigned int j = 0; j < face_idx.size(); j++) {
      if (face_idx[j] != f_idx)
        n_faces.push_back(face_idx[j]);
    }
  }

  return n_faces;
}

Vec3d calc_stellation_point(Geometry &geom, int f_idx, const double eps)
{
  //const vector<vector<int>> &faces = geom.faces();
  const vector<Vec3d> &verts = geom.verts();

  vector<int> n_faces = neighbor_faces(geom, f_idx);
//  for (unsigned int i = 0; i < n_faces.size(); i++) {
//    fprintf(stderr,"n_face = %d\n",n_faces[i]);
//  }

  vector<int> face0 = geom.faces(n_faces[0]);
  vector<int> face1 = geom.faces(n_faces[1]);
  vector<int> face2 = geom.faces(n_faces[2]);

  Vec3d P;
  if (!three_plane_intersect(centroid(verts, face0), face_norm(verts, face0),
                             centroid(verts, face1), face_norm(verts, face1),
                             centroid(verts, face2), face_norm(verts, face2),
                             P, eps)) {
//    P.dump("intersection");
    fprintf(stderr,"stellation point not found for face %d\n", f_idx);
  }

  return P;
}
*/

// facelets have unresolved color map numbers
Geometry make_stellation_diagram(Geometry &geom, int f_idx, string sym_string,
                                 int projection_width, const double eps) {
  const vector<vector<int>> &faces = geom.faces();
  const vector<Vec3d> &verts = geom.verts();

  GeometryInfo geom_info(geom);
  geom_info.set_center(geom.centroid());
  double vert_radius = geom_info.vert_dist_lims().max * projection_width;

  Geometry diagram;

  vector<int> face0 = faces[f_idx];
  Vec3d face0_centroid = centroid(verts, face0);
  Vec3d face0_normal = face_norm(verts, face0);

  // Color face_color = geom.colors(FACES).get(f_idx);

  // find where all other faces intersect face(f_idx)
  for (int i = 0; i < (int)faces.size(); i++) {
    // don't intersect self
    if (i == f_idx)
      continue;

    // experimental, didn't work well
    // if using faces of the same color
    // if (geom.colors(FACES).get(i) == face_color)
    //  continue;

    vector<int> face1 = faces[i];
    Vec3d face1_centroid = centroid(verts, face1);
    Vec3d face1_normal = face_norm(verts, face1);
    // digons won't work in plane intersection
    if (face1.size() < 3)
      continue;

    /*
        // if pre-check for parallel
        Vec3d parallel = vcross(face0_normal, face1_normal);
        if (double_eq(parallel[0], 0.0, eps) &&
            double_eq(parallel[1], 0.0, eps) &&
            double_eq(parallel[2], 0.0, eps)) {
          fprintf(stderr,"%d: parallel = %.17lf,%.17lf,%.17lf\n", i,
       parallel[0], parallel[1], parallel[2]);
          continue;
        }
    */

    // will return false if parallel
    Vec3d P, dir;
    if (two_plane_intersect(face0_centroid, face0_normal, face1_centroid,
                            face1_normal, P, dir, eps)) {
      if (!P.is_set())
        continue;

      diagram.add_vert(P - dir * vert_radius);
      diagram.add_vert(P + dir * vert_radius);
      int sz = diagram.verts().size();
      diagram.add_edge(make_edge(sz - 2, sz - 1));
    }
  }

  mesh_edges(diagram, eps);

  // trim off edges with only one connection
  vector<int> del_verts;
  for (unsigned int i = 0; i < diagram.verts().size(); i++) {
    vector<int> edge_idx = find_edges_with_vertex(diagram.edges(), i);
    if (edge_idx.size() == 1)
      del_verts.push_back(i);
  }
  if (del_verts.size())
    diagram.del(VERTS, del_verts);

  // fill in the faces
  int planar_merge_type = 1;   // tile
  vector<int> nonconvex_faces; // dummy
  fill_in_faces(diagram, planar_merge_type, nonconvex_faces, face0_normal, eps);

  // color diagram by symmetry (using map indexes)
  Symmetry diagram_full_sym(diagram);
  Symmetry diagram_sub_sym = diagram_full_sym;

  vector<vector<set<int>>> sym_equivs;
  diagram_sub_sym.init(diagram, &sym_equivs);

  // sub symmetry specified?
  string geom_sym_symbol = sym_string;
  if (!geom_sym_symbol.length()) {
    Symmetry geom_full_sym(geom);
    geom_sym_symbol = geom_full_sym.get_symbol().c_str();
  }

  // check for horizontal symmetry
  if (geom_sym_symbol.find("h") == string::npos) {
    string diagram_sym_symbol = diagram_full_sym.get_symbol().c_str();
    // fprintf(stderr, "diagram_sym_symbol = %s\n", diagram_sym_symbol.c_str());
    // if not vertical symmetry, will change to chiral for diagram
    replace(diagram_sym_symbol.begin(), diagram_sym_symbol.end(), 'D', 'C');
    // if it has vertical symmetry, change to bi-lateral for diagram
    if (geom_sym_symbol.find("v") != string::npos)
      diagram_sym_symbol = "C2v";
    // fprintf(stderr, "diagram_sym_symbol = %s\n", diagram_sym_symbol.c_str());
    diagram_full_sym.get_sub_sym(diagram_sym_symbol, &diagram_sub_sym);
  }

  // color by sub symmetry
  get_equiv_elems(diagram, diagram_sub_sym.get_trans(), &sym_equivs);

  // remap uses indexes which will be resolved per program
  Coloring clrng(&diagram);
  ColorMap *cmap = colormap_from_name("remap");
  clrng.add_cmap(cmap);
  clrng.f_sets(sym_equivs[2], true);

  return diagram;
}

void split_pinched_faces(Geometry &geom, double eps) {
  vector<vector<int>> &faces = geom.raw_faces();

  // have to keep looping in case faces pinched more than once
  bool found;
  do {
    found = false;
    vector<vector<int>> faces_new;
    vector<Color> faces_new_colors;
    vector<int> del_faces;
    for (unsigned int i = 0; i < faces.size(); i++) {
      // record color
      Color c = geom.colors(FACES).get(i);

      // find duplicate vertices in face path
      vector<int> dups;
      for (unsigned int j = 0; j < faces[i].size(); j++) {
        for (unsigned int k = j + 1; k < faces[i].size(); k++) {
          if (faces[i][j] == faces[i][k])
            dups.push_back(faces[i][j]);
        }
      }

      if (dups.size()) {
        del_faces.push_back(i);
        found = true;
      }

      for (unsigned int j = 0; j < dups.size(); j++) {
        vector<int> face = faces[i];
        int sz = face.size();
        for (int k = 0; k < sz; k++) {
          int l = 0;
          if (face[k] == dups[j]) {
            l = k;
            vector<int> facelet;
            do {
              facelet.push_back(face[l % sz]);
              l++;
            } while (face[l % sz] != dups[j]);
            faces_new.push_back(facelet);
            faces_new_colors.push_back(c);
          }
        }
      }
    }

    if (del_faces.size())
      geom.del(FACES, del_faces);
    if (faces_new.size()) {
      for (unsigned int i = 0; i < faces_new.size(); i++)
        geom.add_face(faces_new[i], faces_new_colors[i]);
    }
  } while (found);

  // coincident faces may have occurred
  int blend_type = 1; // first color
  merge_coincident_elements(geom, "f", blend_type, eps);
}

// use diagrams face color indexes to determine inclusion
vector<vector<int>> lists_full(map<int, Geometry> &diagrams,
                               const vector<vector<int>> &idx_lists,
                               bool remove_multiples) {
  vector<vector<int>> idx_lists_full;
  map<int, vector<int>> lists;

  for (unsigned int i = 0; i < idx_lists.size(); i++) {
    vector<int> map_indexes;
    // position 0 is stellation face index
    int stellation_face_idx = idx_lists[i][0];
    // collect color map indexes of chosen faces
    for (unsigned int j = 1; j < idx_lists[i].size(); j++) {
      int face_no = idx_lists[i][j];
      map_indexes.push_back(
          diagrams[stellation_face_idx].colors(FACES).get(face_no).get_index());
    }
    // if list is empty, add stellation face at position 0
    if (!lists[stellation_face_idx].size())
      lists[stellation_face_idx].push_back(stellation_face_idx);
    // find all faces with same color map indexes as those selected
    const vector<vector<int>> &faces = diagrams[stellation_face_idx].faces();
    for (unsigned int j = 0; j < faces.size(); j++) {
      int face_color_index =
          diagrams[stellation_face_idx].colors(FACES).get(j).get_index();
      vector<int>::iterator it;
      it = find(map_indexes.begin(), map_indexes.end(), face_color_index);
      if (it != map_indexes.end())
        lists[stellation_face_idx].push_back(j);
    }
  }

  // loop through map, place contents in vector list, sort
  int i = 0;
  for (auto const &key1 : lists) {
    idx_lists_full.push_back(lists[key1.first]);
    sort(idx_lists_full[i].begin() + 1, idx_lists_full[i].end());
    i++;
  }

  // look for multiple entries
  // faces with multiple entries are not included at all
  if (remove_multiples) {
    for (unsigned int i = 0; i < idx_lists_full.size(); i++) {
      int stellation_face_idx = idx_lists_full[i][0];

      // make hash table of facelet occurrences
      // position 0 is setllation face index
      map<int, int> occurrences;
      for (unsigned int j = 1; j < idx_lists_full[i].size(); j++) {
        occurrences[idx_lists_full[i][j]]++;
      }

      // re-write facelet indexes, elminating multiples
      idx_lists_full[i].clear();
      idx_lists_full[i].push_back(stellation_face_idx);
      for (auto const &key1 : occurrences) {
        if (key1.second == 1)
          idx_lists_full[i].push_back(key1.first);
      }
    }
  }

  return idx_lists_full;
}

// use diagrams and full lists and face color indexes to determine minimal
// resolved faces
vector<vector<int>> lists_resolved(const Geometry &geom,
                                   map<int, Geometry> &diagrams,
                                   const vector<vector<int>> &idx_lists,
                                   const vector<vector<int>> &idx_lists_full,
                                   bool remove_multiples) {
  map<pair<int, int>, int> resolved_faces;
  for (unsigned int i = 0; i < idx_lists_full.size(); i++) {
    // position 0 is stellation face index
    int stellation_face_idx = idx_lists_full[i][0];

    Vec3d face_normal = face_norm(geom.verts(), geom.faces(0));
    if (vdot(face_normal, geom.face_cent(0)) < 0)
      face_normal *= -1.0;

    // get angles of all faces relative to vertex 0 of stellation face
    vector<pair<double, int>> angle_list;
    for (unsigned int j = 1; j < idx_lists_full[i].size(); j++) {
      double angle = angle_around_axis(
          diagrams[stellation_face_idx].face_cent(idx_lists_full[i][j]),
          geom.verts(geom.faces(0)[0]), face_normal);
      angle_list.push_back(make_pair(rad2deg(angle), idx_lists_full[i][j]));
    }

    sort(angle_list.begin(), angle_list.end());
    // reverse(angle_list.begin(), angle_list.end());

    // map resolved faces by color map index, condense by color map indexes
    for (unsigned int j = 0; j < angle_list.size(); j++) {
      int k = angle_list[j].second;
      int face_color_index =
          diagrams[stellation_face_idx].colors(FACES).get(k).get_index();
      pair<int, int> key = make_pair(stellation_face_idx, face_color_index);
      resolved_faces[key] = k;
    }
  }

  vector<vector<int>> idx_lists_resolved;
  // build resolved list using resolved faces with same map index as original
  // faces
  for (unsigned int i = 0; i < idx_lists.size(); i++) {
    vector<int> resolved_list;
    int stellation_face_idx = idx_lists[i][0];
    resolved_list.push_back(stellation_face_idx);
    for (unsigned int j = 1; j < idx_lists[i].size(); j++) {
      int original_color_index = diagrams[stellation_face_idx]
                                     .colors(FACES)
                                     .get(idx_lists[i][j])
                                     .get_index();
      pair<int, int> key = make_pair(stellation_face_idx, original_color_index);
      int resolved_face = resolved_faces[key];
      // zero is a special case, check if really in the map
      key = make_pair(stellation_face_idx, 0);
      if ((resolved_face != 0) ||
          (resolved_faces.find(key) != resolved_faces.end()))
        resolved_list.push_back(resolved_face);
    }

    // don't sort element 0
    sort(resolved_list.begin() + 1, resolved_list.end());
    if (remove_multiples) {
      // a resolved list can itself contain duplicates
      auto last = unique(resolved_list.begin() + 1, resolved_list.end());
      resolved_list.resize(last - resolved_list.begin());
    }

    // a resolved list can be empty, save for the stellation face in position 0
    if (resolved_list.size() > 1)
      idx_lists_resolved.push_back(resolved_list);
  }

  return idx_lists_resolved;
}

// face_index is changed
Geometry faces_to_geom_for_stel(const Geometry &geom, vector<int> face_idxs) {
  // the index list now contains only the diagram's face numbers
  face_idxs.erase(face_idxs.begin());
  return (faces_to_geom(geom, face_idxs));
}

// index lists contains stellation face number in position 0
Geometry make_stellation(const Geometry &geom, map<int, Geometry> &diagrams,
                         const vector<vector<int>> &idx_lists,
                         const string &sym_string, bool merge_faces,
                         bool remove_inline_verts, bool split_pinched,
                         bool resolve_faces, bool remove_multiples,
                         string map_string, double eps) {
  Geometry stellation_full;

  vector<vector<int>> lists = idx_lists;
  if (resolve_faces) {
    vector<vector<int>> idx_lists_full =
        lists_full(diagrams, idx_lists, remove_multiples);
    lists = lists_resolved(geom, diagrams, idx_lists, idx_lists_full,
                           remove_multiples);
  }

  for (unsigned int i = 0; i < lists.size(); i++) {
    int stellation_face_idx = lists[i][0];

    Geometry stellation =
        faces_to_geom_for_stel(diagrams[stellation_face_idx], lists[i]);
    stellation.colors(VERTS).clear();

    if (merge_faces) {
      // save the color map indexes
      vector<Color> cols;
      for (unsigned int j = 0; j < stellation.faces().size(); j++)
        cols.push_back(stellation.colors(FACES).get(j));

      // change to frame, color map indexes are lost
      stellation.add_missing_impl_edges();
      stellation.clear(FACES);

      // merge selected faces
      int planar_merge_type = 2;   // merge
      vector<int> nonconvex_faces; // dummy
      fill_in_faces(stellation, planar_merge_type, nonconvex_faces,
                    face_norm(diagrams[stellation_face_idx].verts(),
                              diagrams[stellation_face_idx].faces(0)),
                    eps);

      // re-apply color map indexes (approximation)
      for (unsigned int j = 0; j < stellation.faces().size(); j++)
        stellation.colors(FACES).set(j, cols[j]);

      // redraw edges
      stellation.raw_edges().clear();
      stellation.add_missing_impl_edges();

      // eliminate in-line vertices
      if (remove_inline_verts) {
        GeometryInfo info(stellation);
        const vector<vector<int>> &vcons = info.get_vert_cons();

        vector<int> del_verts;
        for (unsigned int j = 0; j < stellation.faces().size(); j++) {
          unsigned int fsz = stellation.faces(j).size();
          for (unsigned int k = 0; k < fsz; k++) {
            Vec3d v1 = stellation.verts(stellation.faces(j)[k]);
            int v2_idx = stellation.faces(j)[(k + 1) % fsz];
            // test only when they have 2 connections
            if (vcons[v2_idx].size() != 2)
              continue;
            Vec3d v2 = stellation.verts(v2_idx);
            Vec3d v3 = stellation.verts(stellation.faces(j)[(k + 2) % fsz]);

            // if v2 is in line with v1 and v3 then mark it for deletion
            if (point_in_segment(v2, v1, v3, eps).is_set())
              del_verts.push_back(v2_idx);
          }
        }
        stellation.del(VERTS, del_verts);

        // its possible for free vertices to be formed
        stellation.del(VERTS, stellation.get_info().get_free_verts());
      }

      // split pinched faces into facelets
      if (split_pinched)
        split_pinched_faces(stellation, eps);
    }

    // color faces
    Coloring clrng(&stellation);
    ColorMap *cmap = colormap_from_name(map_string.c_str());
    clrng.add_cmap(cmap);
    clrng.f_apply_cmap();

    // edge and vertices take color from faces
    stellation.add_missing_impl_edges();
    clrng.e_face_color();
    clrng.v_face_color();

    stellation_full.append(stellation);
  }

  // use specified symmetry
  Symmetry sym;
  sym.init(sym_string);

  // mirror face symmetrically
  sym_repeat(stellation_full, stellation_full, sym);

  int blend_type = 3; // rgb
  merge_coincident_elements(stellation_full, "vef", blend_type, eps);

  // orient the result
  stellation_full.orient();

  return stellation_full;
}

// RK - functions for winding number

// Copyright 2001, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.

//    a Point is defined by its coordinates {int x, y;}
//===================================================================

// isLeft(): tests if a point is Left|On|Right of an infinite line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 and P1
//            =0 for P2 on the line
//            <0 for P2 right of the line
//    See: the January 2001 Algorithm "Area of 2D and 3D Triangles and Polygons"

// inline int
// isLeft( Point P0, Point P1, Point P2 )

// RK - this works with wn_PnPoly better if a double is passed back

static double isLeft(const Vec3d &P0, const Vec3d &P1, const Vec3d &P2,
                     const int idx)
{
  int idx1 = (idx + 1) % 3;
  int idx2 = (idx + 2) % 3;

  double P0_x = P0[idx1];
  double P0_y = P0[idx2];

  double P1_x = P1[idx1];
  double P1_y = P1[idx2];

  double P2_x = P2[idx1];
  double P2_y = P2[idx2];

  return ((P1_x - P0_x) * (P2_y - P0_y) - (P2_x - P0_x) * (P1_y - P0_y));
}

//===================================================================

// cn_PnPoly(): crossing number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  0 = outside, 1 = inside
// This code is patterned after [Franklin, 2000]

// bool
// cn_PnPoly( Point P, Point* V, int n )

// RK - try to use same variable names as original pnpoly

/* RK - Not currently used. none of the tests consider epsilon; if used
recommend using functions

bool cn_PnPoly(const Geometry &polygon, const Vec3d &P, const int idx, int
&crossing_number, double eps)
{
   const vector<Vec3d> &verts = polygon.verts();

   int idx1 = (idx+1)%3;
   int idx2 = (idx+2)%3;

   double testx = P[idx1]; // P.x
   double testy = P[idx2]; // P.y

   const vector<int> &face = polygon.faces()[0];
   int n = face.size();

   int cn = 0;    // the crossing number counter

   // loop through all edges of the polygon
   for (int i=0; i<n; i++) {    // edge from face[i] to face[i+1]
      int j = (i+1)%n;

      double vertx_i = verts[face[i]][idx1]; // V[i].x
      double verty_i = verts[face[i]][idx2]; // V[i].y
      double vertx_j = verts[face[j]][idx1]; // V[i+1].x
      double verty_j = verts[face[j]][idx2]; // V[i+1].y

      if (((verty_i <= testy) && (verty_j > testy))          // an upward
crossing
       || ((verty_i > testy) && (verty_j <= testy))) {       // a downward
crossing
         // compute the actual edge-ray intersect x-coordinate
         double vt = (testy - verty_i) / (verty_j - verty_i);
         if (testx < vertx_i + vt * (vertx_j - vertx_i))     // testx <
intersect
            ++cn;   // a valid crossing of y=testy right of testx
       }
   }

   crossing_number = cn;

   //return (cn&1);    // 0 if even (out), and 1 if odd (in)
   return (cn&1 ? true : false); // false if even (out), and true if odd (in)
}
*/
//===================================================================

// wn_PnPoly(): winding number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  wn = the winding number (=0 only if P is outside V[])

// bool
// wn_PnPoly( Point P, Point* V, int n )

// RK - try to use same variable names as original pnpoly
// RK - there is proof that there are mistakes when epsilon is not considered
// RK - wn is n time too large. winding_number = wn/n

static bool wn_PnPoly(const Geometry &polygon, const Vec3d &P, const int idx,
                      int &winding_number, double eps)
{
  winding_number = 0;
  if (!P.is_set() || !polygon.faces().size())
    return false;

  const vector<Vec3d> &verts = polygon.verts();

  int idx2 = (idx + 2) % 3;

  double testy = P[idx2]; // P.y

  const vector<int> &face = polygon.faces()[0];
  int n = face.size();

  int wn = 0; // the winding number counter

  // loop through all edges of the polygon
  for (int i = 0; i < n; i++) { // edge from face[i] to face[i+1]
    int j = (i + 1) % n;

    Vec3d vert_i = verts[face[i]]; // V[i]
    Vec3d vert_j = verts[face[j]]; // V[[i+1]

    double verty_i = vert_i[idx2]; // V[i].y
    double verty_j = vert_j[idx2]; // V[i+1].y

    for (int i = 0; i < n; i++) { // edge from face[i] to face[i+1]
      // if (verty_i <= testy) {   // start y <= P.y
      if (double_le(verty_i, testy, eps)) {
        // if (verty_j > testy)     // an upward crossing
        if (double_gt(verty_j, testy, eps))
          // if (isLeft( vert_i, vert_j, P, idx ) > eps)  // P left of edge
          if (double_gt(isLeft(vert_i, vert_j, P, idx), 0, eps))
            ++wn; // have a valid up intersect
      }
      else { // start y > P.y (no test needed)
        // if (verty_j <= testy)              // a downward crossing
        if (double_le(verty_j, testy, eps))
          // if (isLeft( vert_i, vert_j, P, idx ) < -eps)  // P right of edge
          if (double_lt(isLeft(vert_i, vert_j, P, idx), 0, eps))
            --wn; // have a valid down intersect
      }
    }
  }

  if (n)
    winding_number = wn / n;

  return (wn ? true : false);
}

int get_winding_number(const Geometry &polygon, const Vec3d &point, double eps)
{
  int winding_number;
  wn_PnPoly(polygon, point, 2, winding_number, eps);
  return winding_number;
}

// the winding number is correct when the faces are placed on the xy-plane
// facing z-positive
// geom contains one polygon
// if points is empty, all serial triads of points in a face are tested and the
// highest magnitude winding number is returned
// if complex polygon has both -W and +W winding, +W is returned
int get_winding_number_polygon(const Geometry &polygon,
                               const vector<Vec3d> &points,
                               const Normal &face_normal, bool find_direction,
                               double eps)
{
  // make copies
  Geometry pgon = polygon;
  vector<Vec3d> pts = points;

  int psz = (int)pts.size();
  if (!psz) {
    const vector<Vec3d> &verts = pgon.verts();
    const vector<int> &face = pgon.faces()[0];

    vector<Vec3d> triangle;
    int fsz = (int)face.size();
    for (int i = 0; i < fsz; i++) {
      triangle.push_back(verts[face[i]]);
      triangle.push_back(verts[face[(i + 1) % fsz]]);
      triangle.push_back(verts[face[(i + 2) % fsz]]);
      pts.push_back(centroid(triangle));
      triangle.clear();
    }
    psz = (int)pts.size();
  }

  Vec3d norm = face_normal.outward().unit();
  double z = 1;
  if (find_direction) {
    double test = 0.0;
    for (int i = 2; i >= 0; i--) {
      test = norm[i];
      if (!double_eq(test, 0.0, eps))
        break;
    }
    z = (test > 0.0) ? 1 : -1;
  }

  // rotate polygon to face forward
  Trans3d trans = Trans3d::rotate(norm + pgon.face_cent(0), Vec3d(0, 0, z));
  pgon.transform(trans);

  // also rotate points per matrix
  Geometry vgeom;
  ;
  for (int i = 0; i < psz; i++)
    vgeom.add_vert(pts[i]);
  vgeom.transform(trans);
  pts.clear();
  for (int i = 0; i < psz; i++)
    pts.push_back(vgeom.verts()[i]);
  vgeom.clear_all();

  int winding_number = 0;
  for (int i = 0; i < psz; i++) {
    // projection index = 2;
    int wn = 0;
    wn_PnPoly(pgon, pts[i], 2, wn, eps);
    if ((abs(wn) > abs(winding_number)) ||
        ((wn > 0) && (wn == -winding_number)))
      winding_number = wn;
  }
  return winding_number;
}

int find_polygon_denominator_signed(const Geometry &geom, int face_idx,
                                    double eps)
{
  const vector<int> &face = geom.faces()[face_idx];

  Normal face_normal(geom, face_idx, centroid(geom.verts()), eps);

  vector<int> sface_idxs;
  sface_idxs.push_back(face_idx);
  Geometry polygon = faces_to_geom(geom, sface_idxs);

  vector<Vec3d> points; // empty points to trigger test triangles
  int d = get_winding_number_polygon(polygon, points, face_normal, false, eps);

  int fsz = (int)face.size();
  if (d < 0)
    d += fsz;

  // check for condition we missed the polygon
  if (!d)
    d = 1;

  // hemispherical faces always are treated as positively wound so return them
  // signed
  if (face_normal.is_hemispherical() && d > fsz / 2)
    d = fsz - d;

  return d;
}


} // namespace anti
