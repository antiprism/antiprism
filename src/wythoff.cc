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
   Name: wyt.cc
   Description: make tilings/polyhedra with a generalised Wythoff construction
   Project: Antiprism - http://www.antiprism.com
*/

#include <algorithm>
#include <functional>
#include <memory>
#include <map>
#include <cctype>
#include <math.h>
#include <string.h>
#include <string>
#include <vector>
#include <regex>

#include "../base/antiprism.h"

using std::string;
using std::vector;
using std::map;
using std::swap;
using std::not_equal_to;

using namespace anti;

void make_meta(const Geometry &geom, Geometry &meta, double face_ht=0.0)
{
  meta.clear_all();
  meta.add_verts(geom.verts());
  int f_start = meta.verts().size();
  for (unsigned int f = 0; f < geom.faces().size(); f++) {
    Vec3d face_pt = geom.face_cent(f);
    if(face_ht)
      face_pt += geom.face_norm(f).with_len(face_ht);
    meta.add_vert(face_pt);
  }

  auto ef_pairs = geom.get_edge_face_pairs();
  for (auto &ef : ef_pairs) {
    // The edge and face pair make a quadrilateral
    // Add the centre to this quadrilateral
    int e_idx = meta.add_vert(geom.edge_cent(ef.first));
    // Add four triangles
    if(ef.second[0] >= 0) {
      meta.add_face(ef.first[0], e_idx, ef.second[0] + f_start, -1);
      meta.add_face(ef.first[1], e_idx, ef.second[0] + f_start, -1);
    }
    if(ef.second[1] >= 0) {
      meta.add_face(ef.first[1], e_idx, ef.second[1] + f_start, -1);
      meta.add_face(ef.first[0], e_idx, ef.second[1] + f_start, -1);
    }
  }
}

Status normalize_tri(Geometry &geom, int f_idx, int v0, int v1,
                     Color other_v_col)
{
  vector<int> &face = geom.faces(f_idx);
  if (face.size() != 3)
    return Status::error(
        msg_str("face %d is not a triangle", f_idx));
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
  for(int i=0; i<3; i++) {
    other_v_idx = face[i];
    if(other_v_idx != v0 && other_v_idx != v1)
      break;
  }

  Color this_other_v_col = geom.colors(VERTS).get(other_v_idx);
  if (this_other_v_col.is_set() && this_other_v_col != other_v_col)
      return Status::error("vertices cannot be 3-coloured");
  else
    geom.colors(VERTS).set(other_v_idx, other_v_col);

  for(int i=0; i<3 ; i++)
    if (geom.colors(VERTS).get(face[i]) == Color(0))
      std::rotate(face.begin(), face.begin() + i, face.end());

  return Status::ok();
}

Status normalize_meta(Geometry &geom)
{
  geom.clear_cols();
  if(geom.faces().size() == 0 || geom.faces().size()%2)
    return Status::error(
        msg_str("geometry does not have an even number of faces"));

  //int part_num = 0;
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
    geom.colors(VERTS).set(geom.faces_mod(i, 0), Color(0));  // V verts colour 0
    geom.colors(VERTS).set(geom.faces_mod(i, 1), Color(1));  // E verts colour 1
    geom.colors(VERTS).set(geom.faces_mod(i, 2), Color(2));  // F verts colour 2

    int cur_fidx = i;
    //if (parts)
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
        swap(e_verts[0], e_verts[1]);
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
    //part_num++;
  }

  // Reorder faces to have colours 0,1,0,1,... will recolour by order later
  vector<int> bad[2];
  for(int i=0; i<(int)geom.faces().size(); i++)
    if(geom.colors(FACES).get(i).get_index() != i%2)
      bad[i%2].push_back(i);
  for(int i=0; i<(int)bad[0].size(); i++)
    std::swap(geom.raw_faces()[bad[0][i]], geom.raw_faces()[bad[1][i]]);

  return Status::ok();
}


class Tile {
private:
  vector<int> idxs;     // index numbers of points visited
  vector<int> ops;      // operations between visiting points
  char start_faces;

  mutable vector<int>::const_iterator ops_i;
  mutable vector<int>::const_iterator idxs_i;

public:
  enum { END = -1, V = 0, E, F, P };
  Status read(const string &pat);
  // string get_pattern() const;
  unsigned char get_start_faces() const { return start_faces; }

  void start_op() const
  {
    ops_i = ops.begin();
    idxs_i = idxs.begin();
  }
  void next_op() const
  {
    if (ops_i != ops.end()) {
      ++ops_i;
      if (*ops_i == P)
        ++idxs_i;
    }
  }
  int get_op() const
  {
    if (ops_i == ops.end())
      return END;
    else
      return *ops_i;
  }
  int get_idx() const { return *idxs_i; }

  void relabel(vector<int> relab);
  string tile_string();
};


Status Tile::read(const string &pat)
{
  //initialise
  ops.clear();
  idxs.clear();
  bool has_tris_spec = strchr("+-*", pat[0]);
  start_faces = (has_tris_spec) ? pat[0] : '+';
  unsigned int pos = has_tris_spec;
  if (!std::isdigit(pat[pos]))
      return Status::error("tile specifier: first character (or first "
          "character after +-*)  must be a digit");

  while (pos < pat.size()) {
    int len;
    // point
    if ((len = strspn(pat.substr(pos).c_str(), "0123456789"))) {
      ops.push_back(P);
      int idx;
      read_int(pat.substr(pos, len).c_str(), &idx);
      idxs.push_back(idx);
    }

    // mirrors
    else if ('v' == pat[pos])
      ops.push_back(V);
    else if ('e' == pat[pos])
      ops.push_back(E);
    else if ('f' == pat[pos])
      ops.push_back(F);

    // rotations
    else if ('V' == pat[pos]) {
      ops.push_back(E);
      ops.push_back(F);
    }
    else if ('E' == pat[pos]) {
      ops.push_back(F);
      ops.push_back(V);
    }
    else if ('F' == pat[pos]) {
      ops.push_back(V);
      ops.push_back(E);
    }

    // no op - stay on same triangle
    else if ('_' == pat[pos])
        ;
    else {
      return Status::error(
          msg_str("invalid character '%c' in position %d", pat[pos], pos + 1));
    }

    if (len)
      pos += len;
    else
      pos++;
  }

  return Status::ok();
}

void Tile::relabel(vector<int> relab)
{
  for(auto &op : ops)
    if(0<=op && op<3)
      op = relab[op];
}

string Tile::tile_string()
{
  vector<int> op2;
  string VEF = "VEF";
  string vef = "vef";
  string tile;
  unsigned int p_idx = 0;
  int last_op = -1;
  for(auto op : ops) {
    if(op == P) {
      if(op == last_op)
        tile += "_";
      if(p_idx < idxs.size())
        tile += itostr(idxs[p_idx]);
      else {
        tile += msg_str("ERROR: index %u out of range", p_idx);
        break;
      }
      p_idx++;
    }
    else
      tile += vef[op];
    last_op = op;
  }

  map<char, int> elem_idx = { {'v', 0}, {'e', 1}, {'f', 2} };
  string tile2;
  for(unsigned int i=0; i<tile.size(); i++) {
    // convert pairs of consecutive letters from vef to VEF
    if (i < tile.size() - 1 && strchr("vef", tile[i]) &&
        strchr("vef", tile[i + 1]) &&
        (elem_idx[tile[i]] + 1) % 3 == elem_idx[tile[i + 1]]) {
      tile2.push_back(VEF[(elem_idx[tile[i]]+2) % 3]);
      i++;  // skip second letter of pair
    }
    else
      tile2.push_back(tile[i]);
  }

  return tile2;
}


struct ConwayOperator {
  std::string operator_short;
  std::string operator_name;
  std::string pattern;
};

ConwayOperator conway_operator_list[]
{
  {"d",   "dual",         "[F]0V,0E"},
  {"a",   "ambo",         "[E]0V,0F"},
  {"S",   "seed",         "[V]0E,0F"},
  {"j",   "join",         "[F,V]0_1E"},
  {"k",   "kis",          "[F,V]0_1v1v,1E"},
  {"t",   "truncate",     "[VE]0v0e,0V,0E"},
  {"n",   "needle",       "[V,F]0_1f1f,1E"},
  {"z",   "zip",          "[EF]0e0f,0F,0E"},
  {"e",   "expand",       "[VF]0v0f,0V,0F"},
  {"o",   "ortho",        "[V,E,F]1_0e1_2e"},
  {"g",   "gyro",         "[F,VE,V]1_0F1_2V1E,1E"},
  {"s",   "snub",         "[VEF]0V,0E,0F,0V0E0F"},
  {"b",   "bevel",        "[VEF]0v0e,0e0f,0f0v"},
  {"m",   "meta",         "[V,E,F]*0_1_2"},
  {"c",   "chamfer",      "[V,VF]1F,0_1v1f"},
  {"u",   "subdivide",    "[V,E]1F,0_1e1e"},
  {"l",   "loft",         "[V,VF]1F,0_1v1_0v,0E"},
  {"q",   "quinto",       "[V,E,EF]2F,0_1_2e2_1e"},
  {"L0",  "joined-lace",  "[V,EF2]1F,1e1_0e,1_0E"},
  {"L",   "lace",         "[V,EF2]1F,1e1_0e,1_0v0v,0E"},
  {"K",   "stake",        "[V,EF2,F]0_1_2e1e,1_0v0v,0E"},
  {"M3",  "edge-medial-3","[F,V,VE]0_2_1e2e,2_0v2v"},
  {"M0",  "joined-medial","[F,V,EF]*0_1_2,1_2E"},
  {"m3",  "edge-medial-3","[F,V,VE]*0_1_2,2_0v2v"},
  {"b3",  "bevel3",       "[VEF,E2F]1_0e0v,0e0f,*1_0f0_1f,1E"},
  {"o3",  "ortho3",       "[V,VF,VE]2_0e2_1e,1F,1_2v2_1v,2E"},
  {"X",   "cross",        "[V,E,F,VF]3_1v3_2v,*0_1_3"},
};


class Tiling {
private:
  vector<Vec3d> points;
  vector<Tile> pat_paths;

  Geometry meta;
  vector<vector<int>> nbrs;
  vector<Vec3d> vert_norms;

  bool one_of_each_tile;

  bool find_nbrs();
  vector<Vec3d> point_on_face(int f_idx, const Vec3d &crds) const;
  void add_circuit(Geometry &geom, int start_idx, const Tile &pat,
                   vector<bool> &seen, Color col) const;
  const vector<Tile> &get_pat_paths() const { return pat_paths; }
  void color_meta();

public:
  Tiling(): one_of_each_tile(false) {}
  Status set_geom(const Geometry &geom, bool is_meta=false, double face_ht=0.0);
  void add_tile(const Tile &pattern) { pat_paths.push_back(pattern); }
  Status add_tile(const string &pat);
  void make_tiling(Geometry &geom, vector<int> *tile_counts = nullptr) const;
  Status read_pattern(const string &pat);
  Status relabel_pattern(const string &relabel);
  Status read_conway(const string &op);
  const Geometry &get_meta() const { return meta; }
  void set_one_of_each_tile(bool val=true) { one_of_each_tile = val; };

  string pattern_string();
  void print_conway_list(FILE *ofile=stdout);
};

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
        swap(e[0], e[1]);
      auto ef_i = ef_pairs.find(e);
      if (ef_i == ef_pairs.end())
        return false;
      else if (ef_i->second.size() != 2)
        nbrs[f][i] = -1;  // only allow connection for two faces at an edge
      else {
        nbrs[f][i] =
          (ef_i->second[0] != (int)f) ? ef_i->second[0] : ef_i->second[1];
      }
    }
  return true;
}

inline vector<Vec3d> Tiling::point_on_face(int f_idx, const Vec3d &crds) const
{
  vector<Vec3d> ret(2);
  // point coordinates
  ret[0] = crds[Tile::V] * meta.face_v(f_idx, Tile::V) +
           crds[Tile::E] * meta.face_v(f_idx, Tile::E) +
           crds[Tile::F] * meta.face_v(f_idx, Tile::F);

  // point normal
  ret[1] =
      crds[Tile::V] * vert_norms[meta.faces(f_idx, Tile::V)] +
      crds[Tile::E] * vert_norms[meta.faces(f_idx, Tile::E)] +
      crds[Tile::F] * vert_norms[meta.faces(f_idx, Tile::F)];
  ret[1].to_unit();
  return ret;
}

void Tiling::add_circuit(Geometry &geom, int start_idx, const Tile &pat,
                        vector<bool> &seen, Color col) const
{
  // Apply pattern until circuit completes
  int start_v_sz = geom.verts().size();
  int idx = start_idx;
  while (true) {
    seen[idx] = true;
    pat.start_op();
    while (pat.get_op() != Tile::END) {
      // fprintf(stderr, "op=%d\n", pat.get_op());
      // fprintf(stderr, "nbrs[%d] = %d, %d, %d\n", idx,
      // nbrs[idx][0],nbrs[idx][1],nbrs[idx][2]);
      // fprintf(stderr, "nbrs[%d][%d] = %d, idx=%d, start_idx=%d\n", idx, op,
      // nbrs[idx][op], idx, start_idx);
      if (pat.get_op() == Tile::P) {
        Vec3d crds = points[pat.get_idx()];
        crds /= crds[0] + crds[1] + crds[2];
        vector<Vec3d> pt = point_on_face(idx, crds);
        geom.add_vert(pt[0]);
      }
      else {
        idx = nbrs[idx][pat.get_op()]; // move to next triangle
        if (idx<0)
          return; // abandon: circuit tried to cross an open edge
      }
      pat.next_op();
    }
    if (idx == start_idx) // circuit complete
      break;
  }

  const int f_sz = geom.verts().size() - start_v_sz;
  vector<int> face(f_sz);
  for (int i = 0; i < f_sz; i++)
    face[i] = start_v_sz + i;
  geom.add_face(face, col);
}

static void reverse_odd_faces(Geometry &geom)
{
  const int f_sz = geom.faces().size();
  for (int i = 0; i < f_sz; i++)
    if (i % 2)
      reverse(geom.faces(i).begin(), geom.faces(i).end());
}

void Tiling::color_meta()
{
  meta.clear_cols();
  const Color light(1.0, 0.8, 0.6, 0.5);
  const Color dark(0.1, 0.3, 0.6, 0.5);
  for(int i = 0; i < (int)meta.faces().size(); i++)
    meta.colors(FACES).set(i, (i%2) ? light : dark);
  meta.add_missing_impl_edges();
  Coloring(&meta).e_one_col(Color::invisible);
  Coloring(&meta).v_one_col(Color::invisible);
}

Status Tiling::set_geom(const Geometry &geom, bool is_meta, double face_ht)
{
  if(is_meta) {
    meta = geom;
    Status stat = normalize_meta(meta);
    if (stat.is_error())
      return stat;
  } else
    make_meta(geom, meta, face_ht);
  color_meta();
  find_nbrs();
  reverse_odd_faces(meta);
  vert_norms = meta.get_info().get_vert_norms();
  reverse_odd_faces(meta);
  return Status::ok();
}

Status Tiling::add_tile(const string &pat)
{
  Tile pattern;
  Status stat = pattern.read(pat);
  if (stat)
    add_tile(pattern);

  return stat;
}

bool valid_start_face(int f, int start_faces)
{
  int pos_tri = f%2;
  return !((start_faces == '-' && pos_tri) || (start_faces == '+' && !pos_tri));
}

void Tiling::make_tiling(Geometry &geom, vector<int> *tile_counts) const
{
  // for(unsigned int i=0; i<nbrs.size(); i++)
  //   fprintf(stderr, "nbrs[%d] = %d, %d, %d\n", i,
  //   nbrs[i][0],nbrs[i][1],nbrs[i][2]);
  geom.clear_all();
  if (tile_counts)
    tile_counts->resize(pat_paths.size(), 0);
  int faces_sz = meta.faces().size();
  for (unsigned int p_idx=0; p_idx<pat_paths.size(); p_idx++) {
    const auto &pat = pat_paths[p_idx];
    Color col(p_idx);
    vector<bool> seen(faces_sz, false);
    int start_faces_sz = geom.faces().size();
    unsigned char start_faces = pat.get_start_faces();
    for (int i = 0; i < faces_sz; i++) {
      if (!seen[i] && valid_start_face(i, start_faces)) {
        add_circuit(geom, i, pat, seen, col);
        if(one_of_each_tile)
          break;
      }
    }
    tile_counts->at(p_idx) = geom.faces().size() - start_faces_sz;
  }
}

Status read_point(const char *point_str, Vec3d &coords)
{
  coords = Vec3d::zero;
  map<char, int> elem_idx = { {'V', 0}, {'E', 1}, {'F', 2} };
  string pt_string(point_str);
  std::regex re_coord("[VEF]([-+]?([0-9]*\\.[0-9]+|[0-9]+))?");
  std::sregex_token_iterator next(pt_string.begin(), pt_string.end(), re_coord, {-1, 0});
  std::sregex_token_iterator end;
  if(next == end)
    return Status::error("invalid coordinate string");

  vector<bool> seen(3, false);
  while (next != end) {
    //fprintf(stderr, "'%s'\n", next->str().c_str());
    if(next->str() != "")
      return Status::error("invalid characters in coordinates: " + next->str());
    next++;
    int idx = elem_idx[next->str()[0]];
    //fprintf(stderr, "'%s'\n", next->str().c_str());
    if(seen[idx])
      return Status::error(
          msg_str("coordinates %c given more than once", next->str()[0]));
    else
      seen[idx] = true;
    if(next->str().size() < 2)
      coords[idx] = 1;
    else {
      Status stat = read_double(next->str().c_str()+1, &coords[idx]);
      if(stat.is_error())
        return stat;
    }

    if(coords.len() == 0)
       return Status::error("coordinates cannot all be zero");

    next++;
  }

  return Status::ok();
}


Status Tiling::read_pattern(const string &pat)
{

  string pat2 = pat;
  std::smatch m_all;
  std::regex r_all("^\\[(.*)](.*)");
  std::regex_match(pat2, m_all, r_all);
  if (m_all.size()<3)
    return Status::error(
        msg_str("pattern '%s': not in form [Point0,Point1,...]Path0,Path1...",
                pat.c_str()));

  std::unique_ptr<char> pat_str(copy_str(m_all[1].str().c_str()));
  vector<char *> parts;
  int num_parts = split_line(pat_str.get(), parts, ",");
  points.resize(num_parts);
  for(int i=0; i<num_parts; i++) {
    //auto stat = points[i].read(parts[i]);
    auto stat = read_point(parts[i], points[i]);
    if (stat.is_error())
      return Status::error(msg_str("Point%d: ", i) +
                           stat.msg());
  }

  std::unique_ptr<char> paths(copy_str(m_all[2].str().c_str()));
  num_parts = split_line(paths.get(), parts, ",");
  pat_paths.resize(num_parts);
  for(int i=0; i<num_parts; i++) {
    auto stat = pat_paths[i].read(parts[i]);
    if (stat.is_error())
      return Status::error(msg_str("Path%d: ", i) + stat.msg());
  }
  return Status::ok();
}


Status Tiling::relabel_pattern(const string &relabel)
{
  if (strlen(relabel.c_str()) != 3 || !strchr(relabel.c_str(), 'V') ||
      !strchr(relabel.c_str(), 'E') || !strchr(relabel.c_str(), 'F'))
    return Status::error(
        "relabel string does not contain exactly three letters V, E and F");

  map<char, int> elem_idx = { {'V', 0}, {'E', 1}, {'F', 2} };
  vector<int> relab(3);
  for(int i=0; i<3; i++)
    relab[i] = elem_idx[relabel[i]];

  for(auto &pt : points) {
    Vec3d v = pt;
    for(int i=0; i<3; i++)
      pt[relab[i]] = v[i];
  }

  for(auto &pat : pat_paths)
    pat.relabel(relab);

  return Status::ok();
}

Status Tiling::read_conway(const string &op)
{
  int last_op = sizeof(conway_operator_list)/sizeof(conway_operator_list[0]);
  for(int i=0; i<last_op; i++)
    if(conway_operator_list[i].operator_short == op)
      return read_pattern(conway_operator_list[i].pattern);
  return Status::error("Conway operator not found");
}

static string coord_string(Vec3d v)
{
  string VEF = "VEF";
  string coords;
  for(int i=0; i<3; i++)
    if(v[i]) {
      coords += VEF[i];
      if (v[i] != 1.0)
        coords += msg_str("%g", v[i]);
    }

  return coords;
}

string Tiling::pattern_string()
{
  string pat = "[";
  for(auto v : points)
    pat += coord_string(v) + ",";
  pat.back() = ']';
  for(auto path : pat_paths)
    pat += path.tile_string() + ",";

  pat.pop_back();
  return pat;
}

void Tiling::print_conway_list(FILE *ofile)
{
  fprintf(ofile, "%-5s%-15s%s\n", "Op", "Description", "Tiling Pattern");
  int last_op = sizeof(conway_operator_list)/sizeof(conway_operator_list[0]);
  for (int i = 0; i < last_op; i++) {
    ConwayOperator op = conway_operator_list[i];
    fprintf(ofile, "%-5s%-15s%s\n", op.operator_short.c_str(),
            op.operator_name.c_str(), op.pattern.c_str());
  }
}


class wy_opts : public ProgramOpts {
private:
public:
  bool input_is_meta;
  bool add_meta;
  double face_ht;
  string relabel;
  Tiling tiling;
  bool color_with_value;
  bool quiet;
  string ifile;
  string ofile;

  wy_opts()
      : ProgramOpts("wythoff"), input_is_meta(false), add_meta(false),
      face_ht(0.0), color_with_value(true), quiet(false)
  {
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

// clang-format off
void wy_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Read a file in OFF format and apply a specified pattern to generate polygon\n"
"tiles. The polyhedron faces are divided by a 'meta' operation into triangles\n"
"each having vertices which are a vertex V, edge centre E and face centre F.\n"
"A start point is positioned on one of these triangles, the next point is\n"
"found by using the pattern to step between triangles, leading to a circuit.\n"
"If input_file is not given the program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"  -p <pat>  pattern in form: [Coords0:Coords1:...]Path0,Path1...\n"
"            Coordinates are barycentric, in form VaEbFc:\n"
"              VEF element letters, and a,b,c are barycentric coordinates\n"
"              corresponding to the preceding element letter. Ommiting a\n"
"              an element letter and coordinate sets the coordinate to zero.\n"
"              Ommitting just the coordinate sets the coordinate to 1. E.g\n"
"            V = (1,0,0), VE = (1,1,0), VE2F3 = (1,2,3)\n"
"            Paths are in the form: TrisPidx0Move0Pidx1Move1...\n"
"              Tris: one of +-* (default +) indicating that paths should\n"
"                start for positive, negative or both kinds of triangles.\n"
"              Pidx: an index number of a point from the coordinates list\n"
"              Move: an operation for stepping to the next triangle, given\n"
"                as a series of characters from the following:\n"
"                  _     - no move, stay on the same triangle\n"
"                  v,e,f - step over side opposite V,E,F\n"
"                  V,E,F - step two trianglesi, rotating about V,E,F,\n"
"                          according to: V=ef, E=fv, F=ve\n"
"  -c <op>   Conway polyhedron notation operator, or 'list' to list all\n"
"            available operators with their corresponding patterns\n"
"  -r        relabel pattern, exactly three letters VEF written in any order\n"
"            e.g. EFV relabels the pattern as V->E,v->e,E->F,e->f,F->V,f->v\n"
"  -M        input geometry is a 'meta' tiling, don't apply meta operation\n"
"  -i        tiles are coloured by index number (default, colour by value)\n"
"  -u        output only one example of each type of tile (one per path)\n"
"  -a        add the 'meta'-transformed base\n"
"  -f <ht>   lift the face centres by this height\n"
"  -q        quiet, don't print report\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}
// clang-format on

void wy_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":ho:p:c:f:r:Miuqa")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'p':
      print_status_or_exit(tiling.read_pattern(optarg), 'p');
      break;

    case 'c':
      if(string(optarg) == "list") {
        tiling.print_conway_list();
        exit(0);
      }
      print_status_or_exit(tiling.read_conway(optarg), 'c');
      break;

    case 'r':
      relabel = optarg;
      break;

      case 'a':
      add_meta = true;
      break;

    case 'M':
      input_is_meta = true;
      break;

    case 'u':
      tiling.set_one_of_each_tile(true);
      break;

    case 'f':
      print_status_or_exit(read_double(optarg, &face_ht), 'f');
      break;

    case 'i':
      color_with_value = false;
      break;

    case 'q':
      quiet = true;
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (argc - optind > 1) {
    error("too many arguments");
    exit(1);
  }

  if (argc - optind == 1)
    ifile = argv[optind];
}

int main(int argc, char *argv[])
{
  wy_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);

  GeometryInfo info(geom);
  if (!info.is_orientable())
    opts.error("base polyhedron is not orientable");
  if (!opts.input_is_meta && !info.is_oriented()) {
    opts.warning("base polyhedron is not oriented; it will be oriented.");
    geom.orient();
  }

  Tiling &tiling = opts.tiling;
  if(opts.relabel != "")
    opts.print_status_or_exit(tiling.relabel_pattern(opts.relabel), 'm');

  Status stat = tiling.set_geom(geom, opts.input_is_meta, opts.face_ht);
  if (stat.is_error())
    opts.print_status_or_exit(stat, 'm');
  Geometry ogeom;
  vector<int> tile_counts;
  tiling.make_tiling(ogeom, &tile_counts);
  if (opts.add_meta)
    ogeom.append(tiling.get_meta());

  if(!opts.quiet) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Tiling pattern: %s\n", tiling.pattern_string().c_str());
    fprintf(stderr, "Tile Counts:\n");
    for(unsigned int i=0; i<tile_counts.size(); i++)
      fprintf(stderr, "  %-4u: %d\n", i, tile_counts[i]);
    fprintf(stderr, "\n");
  }


  merge_coincident_elements(ogeom, "vef");
  if (opts.color_with_value) {
    Coloring clrng(&ogeom);
    clrng.add_cmap(colormap_from_name("spread"));
    clrng.e_apply_cmap();
    clrng.f_apply_cmap();
  }
  opts.write_or_error(ogeom, opts.ofile);

  return 0;
}



