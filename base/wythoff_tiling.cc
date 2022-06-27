/*
   Copyright (c) 2012-2020, Adrian Rossiter

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

#include "geometryinfo.h"
#include "programopts.h"
#include "tiling.h"
#include "utils.h"

#include <cstring>
#include <map>
#include <regex>
#include <string>
#include <vector>

using std::map;
using std::pair;
using std::string;
using std::to_string;
using std::vector;

using namespace anti;

//------------------------------------------------------------
// General Wythoff tiling
static ElemProps<Color> get_original_colors(const Geometry &geom, bool is_meta)
{
  // Get the vertex colours first
  ElemProps<Color> orig_colors;
  for (unsigned int i = 0; i < geom.verts().size(); i++) {
    const auto &col = geom.colors(VERTS).get(i);
    if (col.is_set())
      orig_colors.set(i, col);
  }

  // for meta tiling, this is all the colours, for a polyhedron
  // base, add the face colours and then the edge colours.
  if (!is_meta) {
    int f_start = geom.verts().size(); // offset for face index numbers
    for (unsigned int i = 0; i < geom.faces().size(); i++) {
      const auto &col = geom.colors(FACES).get(i);
      if (col.is_set())
        orig_colors.set(i + f_start, col);
    }

    map<vector<int>, Color> e2col;
    for (unsigned int i = 0; i < geom.edges().size(); i++) {
      Color col = geom.colors(EDGES).get(i);
      if (col.is_set())
        e2col[geom.edges(i)] = col;
    }

    // offset for edge index numbers (where index is position in
    // implicicit edge list)
    int e_start = geom.verts().size() + geom.faces().size();
    GeometryInfo info(geom);
    map<vector<int>, int> e2v;
    int e_idx = 0;
    for (const auto &e : info.get_impl_edges()) {
      auto e_it = e2col.find(e);
      if (e_it != e2col.end())
        orig_colors.set(e_idx + e_start, e_it->second);
      e_idx++;
    }
  }
  return orig_colors;
}

static void make_meta(const Geometry &geom, Geometry &meta,
                      double face_ht = 0.0)
{
  meta.clear_all();
  for (const auto &vert : geom.verts())
    meta.add_vert(vert, Color(0));
  int f_start = meta.verts().size();
  for (unsigned int f = 0; f < geom.faces().size(); f++) {
    Vec3d face_pt = geom.face_cent(f);
    if (face_ht)
      face_pt += geom.face_norm(f).with_len(face_ht);
    meta.add_vert(face_pt, Color(2));
  }

  GeometryInfo info(geom);
  map<vector<int>, int> e2v;
  for (const auto &e : info.get_impl_edges())
    e2v[e] = meta.add_vert(geom.edge_cent(e), Color(1));
  for (int f_idx = 0; f_idx < (int)geom.faces().size(); f_idx++) {
    int f_cent_idx = f_start + f_idx;
    for (int v = 0; v < (int)geom.faces(f_idx).size(); v++) {
      int v0 = geom.faces(f_idx, v);
      int v1 = geom.faces_mod(f_idx, v + 1);
      int e_cent_idx = e2v[make_edge(v0, v1)];
      meta.add_face({v0, e_cent_idx, f_cent_idx});
      meta.add_face({v1, e_cent_idx, f_cent_idx});
    }
  }
}

static Status normalize_tri(Geometry &geom, int f_idx, int v0, int v1,
                            Color other_v_col)
{
  vector<int> &face = geom.faces(f_idx);
  if (face.size() != 3)
    return Status::error(msg_str("face %d is not a triangle", f_idx));
  bool found = false;
  for (int i = 0; i < 3; i++) {
    if (face[i] == v0 && face[(i + 1) % face.size()] == v1) {
      found = true;
      break;
    }
  }
  if (!found)
    std::reverse(face.begin(), face.end());

  int other_v_idx;
  for (int i = 0; i < 3; i++) {
    other_v_idx = face[i];
    if (other_v_idx != v0 && other_v_idx != v1)
      break;
  }

  Color this_other_v_col = geom.colors(VERTS).get(other_v_idx);
  if (this_other_v_col.is_set() && this_other_v_col != other_v_col)
    return Status::error("vertices cannot be 3-coloured");
  else
    geom.colors(VERTS).set(other_v_idx, other_v_col);

  for (int i = 0; i < 3; i++)
    if (geom.colors(VERTS).get(face[i]) == Color(0))
      std::rotate(face.begin(), face.begin() + i, face.end());

  return Status::ok();
}

static Status normalize_meta(Geometry &geom)
{
  geom.clear_cols();
  if (geom.faces().size() == 0 || geom.faces().size() % 2)
    return Status::error(
        msg_str("geometry does not have an even number of faces"));

  // int part_num = 0;
  const int done = -1;
  auto edges = geom.get_edge_face_pairs(false);
  vector<int> cur_idx(geom.faces().size(), 0);
  vector<int> prev_face(geom.faces().size(), 0);
  vector<int> e_verts(2);
  vector<int> orig_e_verts(2);
  vector<int> e_faces(2);
  for (unsigned int i = 0; i < geom.faces().size(); i++) {
    if (geom.faces(i).size() != 3)
      return Status::error(msg_str("face %d is not a triangle", i));
    if (cur_idx[i] == done)
      continue;

    // first face in part has colour 0, and original oriented is preserved
    // and vertices are in order VEF. Acts as seed for all other faces
    geom.colors(FACES).set(i, 0);
    geom.colors(VERTS).set(geom.faces_mod(i, 0), Color(0)); // V verts colour 0
    geom.colors(VERTS).set(geom.faces_mod(i, 1), Color(1)); // E verts colour 1
    geom.colors(VERTS).set(geom.faces_mod(i, 2), Color(2)); // F verts colour 2

    int cur_fidx = i;
    // if (parts)
    //  parts->push_back(vector<int>(1, i));
    while (cur_idx[i] != done) {
      int idx = cur_idx[cur_fidx];
      if (idx == done) {
        cur_fidx = prev_face[cur_fidx];
        continue;
      }

      // read off the next edge
      const vector<int> &face = geom.faces(cur_fidx);
      orig_e_verts[0] = face[idx];
      idx = (idx + 1) % face.size();
      orig_e_verts[1] = face[idx];
      cur_idx[cur_fidx] = idx ? idx : done; // set to next idx, or mark done

      e_verts = orig_e_verts;
      if (e_verts[0] > e_verts[1])
        std::swap(e_verts[0], e_verts[1]);
      e_faces = edges.find(e_verts)->second;

      int next_face = (e_faces[0] != cur_fidx) ? e_faces[0] : e_faces[1];
      if (next_face >= 0 && cur_idx[next_face] == 0) { // face not looked at yet
        Color cur_col = geom.colors(FACES).get(cur_fidx);
        // Adjacent faces must be coloured differently
        if (geom.colors(FACES).get(next_face) == cur_col)
          return Status::error("faces cannot be 2-coloured");
        else
          geom.colors(FACES).set(next_face, Color(!cur_col.get_index()));

        int other_v_idx = (idx + 1) % face.size();
        Color other_v_col =
            geom.colors(VERTS).get(geom.faces(cur_fidx, other_v_idx));
        normalize_tri(geom, next_face, orig_e_verts[0], orig_e_verts[1],
                      other_v_col);
        // if (parts)
        //  parts->back().push_back(next_face);
        prev_face[next_face] = cur_fidx;
        cur_fidx = next_face;
      }
    }
    // part_num++;
  }

  // Reorder faces to have colours 0,1,0,1,... will recolour by order later
  vector<int> bad[2];
  for (int i = 0; i < (int)geom.faces().size(); i++)
    if (geom.colors(FACES).get(i).get_index() != i % 2)
      bad[i % 2].push_back(i);
  for (int i = 0; i < (int)bad[0].size(); i++)
    swap(geom.raw_faces()[bad[0][i]], geom.raw_faces()[bad[1][i]]);

  return Status::ok();
}

Color TilingColoring::get_point_color(const TilingPoint &point) const
{
  Color col; // default unset: ct none
  auto type = coloring_types[POINTS];
  if (type == ColoringType::index)
    col.set_index(point.get_index());
  else if (type == ColoringType::component)
    col.set_index(point.get_inclusion());
  else if (type == ColoringType::weight)
    col = Color(point.get_coords().unit()); // unit coordinates as RGB

  return col;
}

ColorMap *TilingColoring::get_default_colormap(int elem) const
{
  string name = "null"; // associated_element and weight
  if (coloring_types[elem] == ColoringType::index)
    name = "spread";
  else if (coloring_types[elem] == ColoringType::component)
    name = "map_red:green:blue:yellow:cyan:magenta:grey80";
  else // none, associated_element, weight
    name = "null";

  return colormap_from_name(name.c_str());
}

Status TilingColoring::set_coloring(int elem, ColoringType clrng_type,
                                    int assoc, Color local)
{
  Status stat;
  bool is_point = (elem == 1);
  if (clrng_type == ColoringType::index || clrng_type == ColoringType::none ||
      (is_point && clrng_type == ColoringType::component) ||
      (is_point && clrng_type == ColoringType::weight) ||
      clrng_type == ColoringType::associated_element)
    coloring_types[elem] = clrng_type;
  else
    return Status::error("invalid colouring type");

  if (clrng_type == ColoringType::associated_element) {
    if (assoc == Tile::V || assoc == Tile::E || assoc == Tile::F ||
        assoc == Tile::P) {
      local_associations[elem] = assoc;
      local_colors[elem] = local;
    }
    else
      return Status::error("invalid associated element type");
  }

  return Status::ok();
}

Status TilingColoring::read_coloring(const string &clrng_str)
{
  Status stat;
  Split parts_elems(clrng_str.c_str(), ":");
  if (parts_elems.size() > ELEMS_SZ)
    return stat.set_error("more than one colon given");

  // 0 - tile, 1 - point
  for (size_t elem = 0; elem < parts_elems.size(); elem++) {
    auto prefix = string((elem == TILES) ? "tile" : "point") + ": ";
    Split parts(parts_elems[elem], ",");
    if (parts.size() > 2)
      return stat.set_error(prefix + "more than one comma given");

    string params = "none=0|index=1|association=2";
    if (elem)
      params += "|component=3|weight=4";
    string arg_id;
    if (!(stat = ProgramOpts::get_arg_id(parts[0], &arg_id, params.c_str())))
      return stat.add_prefix(prefix);

    int id = atoi(arg_id.c_str());
    if (id != 2 && parts.size() > 1)
      return stat.set_error(
          prefix + "comma given, but only valid for associated colouring");

    if (id == 0)
      set_coloring(elem, ColoringType::none);
    else if (id == 1)
      set_coloring(elem, ColoringType::index);
    else if (id == 2) {
      Color local_col;
      int assoc = Tile::F;
      if (parts.size() > 1) {
        if (strcasecmp(parts[1], "V") == 0)
          assoc = Tile::V;
        else if (strcasecmp(parts[1], "E") == 0)
          assoc = Tile::E;
        else if (strcasecmp(parts[1], "F") == 0)
          assoc = Tile::F;
        else {
          if (!(stat = local_col.read(parts[1])))
            return stat;
          assoc = Tile::P;
        }
      }
      set_coloring(elem, ColoringType::associated_element, assoc, local_col);
    }
    else if (elem == POINTS && id == 3)
      set_coloring(elem, ColoringType::component);
    else if (elem == POINTS && id == 4)
      set_coloring(elem, ColoringType::weight);

    if (elem == TILES) { // default point colouring is tile colouring
      coloring_types[POINTS] = coloring_types[TILES];
      local_associations[POINTS] = local_associations[TILES];
      local_colors[POINTS] = local_colors[TILES];
    }
  }

  return stat;
}

string TilingColoring::get_option_help(char op_char)
{
  return msg_str(
      R"(  -%c <mthd> tile colouring method, optionally followed by ':' and
            point colouring method.
               Tile colouring methods
                 none - do not colour tiles
                 index - use the path index, and apply colormap (default)
                 association - colour tiles using corresponding base element
                   colour, optionally followed by a comma and a letter from V,
                   E, F, or a colour, to use for all local tiles (default: F)
               Point colouring methods (default: use tile coloring method)
                 none - do not colour points
                 index - use the point index, and apply colormap
                 association - colour points using corresponding base element
                   colour, optionally followed by a comma and a letter from V,
                   E, F, or a colour, to use for all local points (default: F)
                 component - colour by non-zero components, V=1, E=2, F=4,
                   V&E=3, E&F=6, V&F=5, V&E&F=7, and apply colormap
                 weight - colour using barycentric coordinates as RGB)",
      op_char);
}

bool Tiling::find_nbrs()
{
  auto ef_pairs = meta.get_edge_face_pairs(false);

  // Find the neighbour face opposite each VEF vertex
  nbrs.resize(meta.faces().size(), vector<int>(3));
  for (unsigned int f = 0; f < meta.faces().size(); f++)
    for (int i = 0; i < 3; i++) {
      vector<int> e(2);
      e[0] = meta.faces_mod(f, i + 1);
      e[1] = meta.faces_mod(f, i + 2);
      if (e[0] > e[1])
        std::swap(e[0], e[1]);
      auto ef_i = ef_pairs.find(e);
      if (ef_i == ef_pairs.end())
        return false;
      else if (ef_i->second.size() != 2)
        nbrs[f][i] = -1; // only allow connection for two faces at an edge
      else {
        nbrs[f][i] =
            (ef_i->second[0] != (int)f) ? ef_i->second[0] : ef_i->second[1];
      }
    }
  return true;
}

static Vec3d point_on_face(const Geometry &meta, int f_idx, const Vec3d &crds)
{
  // point coordinates
  Vec3d P = crds[Tile::V] * meta.face_v(f_idx, Tile::V) +
            crds[Tile::E] * meta.face_v(f_idx, Tile::E) +
            crds[Tile::F] * meta.face_v(f_idx, Tile::F);

  // point normal
  // ret[1] =
  //    crds[Tile::V] * vert_norms[meta.faces(f_idx, Tile::V)] +
  //    crds[Tile::E] * vert_norms[meta.faces(f_idx, Tile::E)] +
  //    crds[Tile::F] * vert_norms[meta.faces(f_idx, Tile::F)];
  // ret[1].to_unit();
  return P;
}

Color Tiling::get_associated_element_point_color(int f_idx, int incl) const
{

  int local_incl = coloring.get_local_association(TilingColoring::POINTS);
  int idx = -1;
  if (incl == Tile::V || incl == Tile::E || incl == Tile::F)
    idx = incl;
  else if (local_incl == Tile::V || local_incl == Tile::E ||
           local_incl == Tile::F)
    idx = local_incl;
  else // no direct or indirect association with a base element
    return coloring.get_local_color(TilingColoring::POINTS);

  return orig_colors.get(meta.faces(f_idx, idx));
}

int Tiling::get_associated_element(int start_idx, const string &step,
                                   int assoc_type) const
{
  const map<char, int> elem_idx = {{'V', 0}, {'E', 1}, {'F', 2}};
  int idx;
  if (assoc_type == Tile::VEF)
    idx = -1; // invalid index
  else {
    idx = start_idx;
    for (char op : step)
      idx = nbrs[idx][elem_idx.at(std::toupper(op))]; // move to next tri
  }
  return idx >= 0 ? meta.faces(idx, assoc_type) : idx;
}

// Each pattern point plotted for a meta triangle cooresponds to
// a previously assigned geometry vertex. Get the index of that vertex.
static int
get_index(const vector<int> &face, int f_idx, int pat_pt_idx, int incl,
          const std::vector<std::map<std::vector<int>, std::pair<int, int>>>
              &index_order,
          const std::vector<int> &point_vertex_offsets)
{
  vector<int> incl_key;
  if (incl >= Tile::V && incl <= Tile::F) // meta tile vertex
    incl_key = {face[incl]};
  else if (incl >= Tile::VE && incl <= Tile::FV) // meta tile edge
    incl_key = make_edge(face[incl % 3], face[(incl + 1) % 3]);
  else if (incl == Tile::VEF) // meta tile interior
    incl_key = {f_idx};
  else
    return -1; // invalid inclusion value, shouldn't happen

  const auto it = index_order[incl].find(incl_key);
  if (it == index_order[incl].end())
    return -2; // invalid element key, shouldn't happen

  return point_vertex_offsets[pat_pt_idx] + it->second.first;
}

void Tiling::add_circuit(
    Geometry &geom, int start_idx, const Tile &pat, std::vector<bool> &seen,
    Color col,
    const std::vector<std::map<std::vector<int>, std::pair<int, int>>>
        &index_order,
    const std::vector<int> &point_vertex_offsets) const
{
  // Apply pattern until circuit completes
  vector<int> face;
  int idx = start_idx;
  while (true) {
    seen[idx] = true;
    pat.start_op();
    while (pat.get_op() != Tile::END) {
      if (pat.get_op() == Tile::P) {
        int incl = points[pat.get_idx()].get_inclusion();
        int v_idx = get_index(meta.faces(idx), idx, pat.get_idx(), incl,
                              index_order, point_vertex_offsets);
        face.push_back(v_idx);
      }
      else {
        idx = nbrs[idx][pat.get_op()]; // move to next triangle
        if (idx < 0)
          return; // abandon: circuit tried to cross an open edge
      }
      pat.next_op();
    }
    if (idx == start_idx) // circuit complete
      break;
  }

  geom.add_face(face, col);
}

static void reverse_odd_faces(Geometry &geom)
{
  const int f_sz = geom.faces().size();
  for (int i = 0; i < f_sz; i++)
    if (i % 2)
      reverse(geom.faces(i).begin(), geom.faces(i).end());
}

Status Tiling::set_geom(const Geometry &geom, bool is_meta, double face_ht)
{
  orig_colors = get_original_colors(geom, is_meta);

  if (is_meta) {
    meta = geom;
    Status stat = normalize_meta(meta);
    if (stat.is_error())
      return stat;
  }
  else
    make_meta(geom, meta, face_ht);

  find_nbrs();
  if (is_meta) {
    // Neighbouring faces must have index numbers of opposite parity
    for (int i = 0; i < (int)nbrs.size(); i++)
      for (int j = 0; j < 3; j++)
        if (i % 2 == nbrs[i][j] % 2)
          return Status::error("faces cannot be 2-coloured");
  }

  reverse_odd_faces(meta);
  // vert_norms = meta.get_info().get_vert_norms();
  reverse_odd_faces(meta);
  return Status::ok();
}

Status Tiling::add_tile(const string &pat)
{
  Tile pattern;
  Status stat = pattern.read(pat);
  if (stat)
    pat_paths.push_back(pattern);

  return stat;
}

void Tiling::reverse_pattern()
{
  for (auto &path : pat_paths)
    path.flip_start_faces();
}

void Tiling::start_everywhere()
{
  for (auto &path : pat_paths)
    path.set_start_faces('*');
}

namespace {

void delete_verts(Geometry &geom, const vector<int> &v_nos)
{
  vector<int> dels = v_nos;
  map<int, int> v_map;
  if (!dels.size())
    return;
  sort(dels.begin(), dels.end());
  unsigned int del_verts_cnt = 0;
  int map_to;
  for (unsigned int i = 0; i < geom.verts().size(); i++) {
    if (del_verts_cnt < dels.size() && (int)i == dels[del_verts_cnt]) {
      del_verts_cnt++;
      map_to = -1;
    }
    else {
      map_to = i - del_verts_cnt;
      geom.verts(map_to) = geom.verts(i);
    }
    v_map[i] = map_to;
  }
  geom.raw_verts().resize(geom.verts().size() - del_verts_cnt);
  geom.colors(VERTS).remap(v_map);

  vector<int> del_faces;
  for (unsigned int i = 0; i < geom.faces().size(); i++) {
    int curr_idx = 0;
    for (unsigned int j = 0; j < geom.faces(i).size(); j++) {
      int new_idx = v_map[geom.faces(i, j)];
      if (new_idx >= 0)
        geom.raw_faces()[i][curr_idx++] = new_idx;
    }
    geom.raw_faces()[i].resize(curr_idx);
    if (curr_idx < 2)
      del_faces.push_back(i);
  }
  geom.del(FACES, del_faces);
}

bool valid_start_face(int f, int start_faces)
{
  int pos_tri = f % 2;
  return !((start_faces == '-' && pos_tri) || (start_faces == '+' && !pos_tri));
}

void store_tri(map<vector<int>, pair<int, int>> &elem_to_tri, vector<int> key,
               int tri_idx)
{
  auto pr = elem_to_tri.insert({key, {-1, tri_idx}});
  if (pr.second == false && is_even(pr.first->second.second) &&
      !is_even(tri_idx))
    pr.first->second.second = tri_idx;
}

}; // namespace

Status Tiling::make_tiling(Geometry &geom,
                           vector<Tile::TileReport> *tile_reports) const
{
  geom.clear_all();
  if (tile_reports)
    tile_reports->resize(pat_paths.size());

  // All the possible element inclusion postions V, E, F, VE, EF, FV, VEF.
  // Each entry maps to order (to find index of corresponding point)
  // and example triangle (to generate coordinates of corresponding point)
  auto ef_pairs = meta.get_edge_face_pairs();
  vector<map<vector<int>, pair<int, int>>> index_order(7);
  for (int i = 0; i < (int)meta.faces().size(); i++) {
    const auto &face = meta.faces(i);
    index_order[Tile::VEF][{i}] = {-1, i};
    index_order[Tile::V][{face[Tile::V]}] = {-1, i};
    index_order[Tile::E][{face[Tile::E]}] = {-1, i};
    index_order[Tile::F][{face[Tile::F]}] = {-1, i};
    store_tri(index_order[Tile::VE], make_edge(face[Tile::V], face[Tile::E]),
              i);
    index_order[Tile::EF][make_edge(face[Tile::E], face[Tile::F])] = {-1, i};
    index_order[Tile::FV][make_edge(face[Tile::F], face[Tile::V])] = {-1, i};
  }
  for (int i = (int)Tile::V; i <= (int)Tile::VEF; i++) {
    int pos = 0;
    for (auto &m : index_order[i])
      m.second.first = pos++;
  }

  // Starting offset of vertices corresponding to each pattern point
  vector<int> point_vertex_offsets(points.size());
  for (int i = 0; i < (int)points.size(); i++) {
    point_vertex_offsets[i] = geom.verts().size();
    const auto &pt = points[i];
    int incl = pt.get_inclusion();
    Vec3d crds = pt.get_coords();
    crds /= crds[0] + crds[1] + crds[2];
    for (auto &m : index_order[incl]) {
      const int f_idx = m.second.second;
      Color col = coloring.get_point_color(pt);
      if (coloring.is_associated_element(TilingColoring::POINTS))
        col = get_associated_element_point_color(f_idx, incl);
      geom.add_vert(point_on_face(meta, f_idx, crds), col);
    }
  }

  int faces_sz = meta.faces().size();
  for (unsigned int p_idx = 0; p_idx < pat_paths.size(); p_idx++) {
    const auto &pat = pat_paths[p_idx];
    // Check index range
    auto out_of_range = pat.check_index_range(points.size());
    if (out_of_range.size()) {
      string msg = "Path" + to_string(p_idx) + ": index numbers out of range:";
      for (auto idx : out_of_range)
        msg += " " + to_string(idx) + ",";
      msg.pop_back();
      return Status::error(msg.c_str());
    }

    auto assoc = pat.get_element_association();
    vector<bool> seen(faces_sz, false);
    int start_faces_sz = geom.faces().size();
    unsigned char start_faces = pat.get_start_faces();
    for (int i = 0; i < faces_sz; i++) {
      if (!seen[i] && valid_start_face(i, start_faces)) {
        Color col; // col_type==ColoringType::none
        if (coloring.is_index(TilingColoring::TILES))
          col.set_index(p_idx);
        else if (coloring.is_associated_element(TilingColoring::TILES)) {
          auto type =
              (assoc.assoc_type == Tile::P)
                  ? coloring.get_local_association(TilingColoring::TILES)
                  : assoc.assoc_type;
          if (type == Tile::P)
            col = coloring.get_local_color(TilingColoring::TILES);
          else {
            int col_idx = get_associated_element(i, assoc.step, type);
            if (col_idx >= 0)
              col = orig_colors.get(col_idx);
          }
        }
        add_circuit(geom, i, pat, seen, col, index_order, point_vertex_offsets);
        if (one_of_each_tile)
          break;
      }
    }
    if (tile_reports) {
      assoc.count = geom.faces().size() - start_faces_sz;
      tile_reports->at(p_idx) = assoc;
    }
  }

  delete_verts(geom, geom.get_info().get_free_verts());
  return Status::ok();
}

namespace {

Status read_point(const char *point_str, Vec3d &coords)
{
  coords = Vec3d::zero;
  map<char, int> elem_idx = {{'V', 0}, {'E', 1}, {'F', 2}};
  string pt_string(point_str);
  std::regex re_coord("([-+]?([0-9]*\\.[0-9]+|[0-9]+))?[VEF]");
  std::sregex_token_iterator next(pt_string.begin(), pt_string.end(), re_coord,
                                  {-1, 0});
  std::sregex_token_iterator end;
  if (next == end)
    return Status::error("invalid coordinate string");

  vector<bool> seen(3, false);
  while (next != end) {
    if (next->str() != "")
      return Status::error("invalid characters in coordinates: " + next->str());
    next++;
    int idx = elem_idx[next->str().back()];
    if (seen[idx])
      return Status::error(
          msg_str("coordinates %c given more than once", next->str().back()));
    else
      seen[idx] = true;
    if (next->str().size() < 2)
      coords[idx] = 1;
    else {
      string coord_str = next->str();
      coord_str.pop_back();
      Status stat = read_double(coord_str.c_str(), &coords[idx]);
      if (stat.is_error())
        return stat;
    }
    next++;
  }

  if (coords.len() == 0)
    return Status::error("coordinates cannot all be zero");

  return Status::ok();
}

}; // namespace

Status Tiling::read_pattern(const string &pat)
{

  string pat2 = pat;
  std::smatch m_all;
  std::regex r_all("^\\[(.*)\\](.*)");
  std::regex_match(pat2, m_all, r_all);
  if (m_all.size() < 3)
    return Status::error(
        msg_str("pattern '%s': not in form [Point0,Point1,...]Path0,Path1...",
                pat.c_str()));

  Split parts;
  int num_parts = parts.init(m_all[1].str(), ",");
  points.resize(num_parts);
  for (int i = 0; i < num_parts; i++) {
    Vec3d coords;
    auto stat = read_point(parts[i], coords);
    if (stat.is_error())
      return Status::error(msg_str("Point%d: ", i) + stat.msg());
    points[i] = TilingPoint(coords, i);
  }

  num_parts = parts.init(m_all[2].str(), ",");
  pat_paths.resize(num_parts);
  for (int i = 0; i < num_parts; i++) {
    auto stat = pat_paths[i].read(parts[i]);
    if (stat.is_error())
      return Status::error(msg_str("Path%d: ", i) + stat.msg());
  }
  return Status::ok();
}

Status Tiling::relabel_pattern(string relabel)
{
  for (auto &c : relabel)
    c = toupper(c);

  if (strlen(relabel.c_str()) != 3 || !strchr(relabel.c_str(), 'V') ||
      !strchr(relabel.c_str(), 'E') || !strchr(relabel.c_str(), 'F'))
    return Status::error(
        "relabel string does not contain exactly three letters V, E and F");

  map<char, int> elem_idx = {{'V', 0}, {'E', 1}, {'F', 2}};
  vector<int> relab(3);
  for (int i = 0; i < 3; i++)
    relab[i] = elem_idx[relabel[i]];

  for (auto &pt : points) {
    Vec3d coords;
    for (int i = 0; i < 3; i++)
      coords[relab[i]] = pt.get_coords()[i];
    pt = TilingPoint(coords, pt.get_index());
  }

  for (auto &pat : pat_paths)
    pat.relabel(relab);

  return Status::ok();
}
