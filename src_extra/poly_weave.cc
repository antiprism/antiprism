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
   Name: weave.cc
   Description: make a weave based on edge mid-points
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"
#include <algorithm>
#include <functional>
#include <map>
#include <math.h>
#include <memory>
#include <regex>
#include <string.h>
#include <string>
#include <vector>

using std::map;
using std::not_equal_to;
using std::string;
using std::swap;
using std::vector;

using namespace anti;

void make_meta(const Geometry &geom, Geometry &meta, double face_ht = 0.0)
{
  meta.clear_all();
  meta.add_verts(geom.verts());
  int f_start = meta.verts().size();
  for (unsigned int f = 0; f < geom.faces().size(); f++) {
    Vec3d face_pt = geom.face_cent(f);
    if (face_ht)
      face_pt += geom.face_norm(f).with_len(face_ht);
    meta.add_vert(face_pt);
  }

  Color light(1.0, 0.8, 0.6, 0.5);
  Color dark(0.1, 0.3, 0.6, 0.5);
  auto ef_pairs = geom.get_edge_face_pairs();
  for (auto &ef_pair : ef_pairs) {
    // The edge and face pair make a quadrilateral
    // Add the centre to this quadrilateral
    int e_idx = meta.add_vert(geom.edge_cent(ef_pair.first));
    // Add four triangles
    if (ef_pair.second[0] >= 0) {
      meta.add_face({ef_pair.first[0], e_idx, ef_pair.second[0] + f_start},
                    light);
      meta.add_face({ef_pair.first[1], e_idx, ef_pair.second[0] + f_start},
                    dark);
    }
    if (ef_pair.second[1] >= 0) {
      meta.add_face({ef_pair.first[1], e_idx, ef_pair.second[1] + f_start},
                    light);
      meta.add_face({ef_pair.first[0], e_idx, ef_pair.second[1] + f_start},
                    dark);
    }
  }
  meta.add_missing_impl_edges();
  Coloring(&meta).e_one_col(Color::invisible);
  Coloring(&meta).v_one_col(Color::invisible);
}

class weave_path {
private:
  Vec3d point;
  vector<vector<double>> path_points;
  int path_type;

public:
  enum { LINE, CURVE };
  Status init(const string &path);
  void set_path_type(int typ) { path_type = typ; }
  Vec3d get_point() const { return point; }
  void add_points(vector<Vec3d> *pts, const vector<Vec3d> &start,
                  const vector<Vec3d> &end) const;
};

Status weave_path::init(const string &path)
{
  path_points.clear();
  int buff_sz = path.size() + 1;
  vector<char> path_str(buff_sz);
  strncpy(&path_str[0], path.c_str(), buff_sz - 1);
  path_str[buff_sz - 1] = '\0';
  vector<char *> parts;
  int parts_sz = split_line(&path_str[0], parts, ":", true);
  Status stat = point.read(parts[0]);
  if (stat.is_error())
    return Status::error(msg_str("path '%s': first (join) point: '%s'",
                                 path.c_str(), stat.c_msg()));
  point /= point[0] + point[1] + point[2];
  for (int i = 1; i < parts_sz; i++) {
    vector<double> nums;
    if (!(stat = read_double_list(parts[i], nums)))
      return Status::error(msg_str("path '%s': path point %d: '%s'",
                                   path.c_str(), i + 1, stat.c_msg()));
    if (nums.size() > 3)
      return Status::error(msg_str("path '%s': path point %d: more "
                                   "than three numbers",
                                   path.c_str(), i + 1));
    path_points.push_back(nums);
    // fprintf(stderr, "parts[%d]='%s'\n", i, parts[i]);
  }

  return Status::ok();
}

static Vec3d get_control_point(Vec3d P, Vec3d Q, Vec3d P_up)
{
  Vec3d along = Q - P;
  Vec3d axis = vcross(along, P_up);
  double turn_ang = angle_around_axis(along, P_up, axis) - M_PI / 2;
  Trans3d trans = Trans3d::translate(P) * Trans3d::rotate(axis, turn_ang) *
                  Trans3d::translate(-P);
  Vec3d pt = (2 * P + Q) / 3; // point 1/3 along, chosen experimentally
  return trans * pt;
}

void weave_path::add_points(vector<Vec3d> *pts, const vector<Vec3d> &start,
                            const vector<Vec3d> &end) const
{
  Vec3d along = end[0] - start[0];
  const int path_sz = path_points.size();
  Vec3d Q0 = get_control_point(start[0], end[0], start[1]); // for CURVE
  Vec3d Q1 = get_control_point(end[0], start[0], end[1]);   // for CURVE
  pts->push_back(start[0]);
  for (int i = 0; i < path_sz; i++) {
    double X = (path_points[i].size() > 0) ? path_points[i][0] : 0.0;
    double Y = (path_points[i].size() > 1) ? path_points[i][1] : 0.0;
    double Z = (path_points[i].size() > 2) ? path_points[i][2]
                                           : (i + 1.0) / (path_sz + 1);
    Vec3d up =
        (Z * end[1] + (1 - Z) * start[1]).unit(); // rotation would be better
    Vec3d side = vcross(along, up).unit();
    Vec3d origin;
    if (path_type == LINE)
      origin = (Z * end[0] + (1 - Z) * start[0]);
    else if (path_type == CURVE) {
      const double z = 1 - Z;
      origin = z * z * z * start[0] + 3 * z * z * Z * Q0 + 3 * z * Z * Z * Q1 +
               Z * Z * Z * end[0];
    }

    pts->push_back(origin + up * X + side * Y);
  }
}

class weave_pattern {
private:
  vector<int> ops;
  vector<weave_path> paths;
  unsigned char start_faces;

  mutable vector<int>::const_iterator ops_i;
  mutable vector<weave_path>::const_iterator paths_i;

public:
  enum { END = -1, V = 0, E, F, P };
  Status set_pattern(const string &pat);
  // string get_pattern() const;
  unsigned char get_start_faces() const { return start_faces; }

  void start_op() const
  {
    ops_i = ops.begin();
    paths_i = paths.begin();
  }
  void next_op() const
  {
    if (ops_i != ops.end()) {
      ++ops_i;
      if (*ops_i == P)
        ++paths_i;
    }
  }
  int get_op() const
  {
    if (ops_i == ops.end())
      return END;
    else
      return *ops_i;
  }
  const weave_path &get_path() const { return *paths_i; }
};

Status weave_pattern::set_pattern(const string &pat)
{
  ops.clear();
  paths.clear();
  start_faces = 1;                   // 'left'/'dark' faces
  int path_type = weave_path::CURVE; // path points to follow a curve

  bool reverse = false;
  int pat_sz = pat.size();
  int pos = 0;
  while (pos < pat_sz) {
    int len;
    // path points
    bool add_default_point = !paths.size() && strchr("vefVEF_", pat[pos]);
    if ((len = strspn(pat.substr(pos).c_str(), "0123456789.,-+:")) ||
        add_default_point) {
      ops.push_back(P);
      Status stat;
      weave_path path;
      if (add_default_point)
        path.init("0,1,0");
      else if (!(stat = path.init(pat.substr(pos, len).c_str())))
        return stat;
      path.set_path_type(path_type);
      paths.push_back(path);
      if (add_default_point)
        continue; // reprocess char that triggered adding of default point
    }
    else if (strchr("t", pat[pos])) {
      if ((int)pat.size() == pos + 1)
        return Status::error("'t' is last character, must be "
                             "followed by l, r or b\n");
      char tris = pat[pos + 1];
      if (!strchr("lrb", tris))
        return Status::error(msg_str("'t' followed by '%c', must be "
                                     "followed by l, r or b\n",
                                     pat[pos]));
      start_faces =
          1 * (tris == 'l' || tris == 'b') + 2 * (tris == 'r' || tris == 'b');
      len = 2;
    }

    // path shapes
    else if ('L' == pat[pos]) {
      path_type = weave_path::LINE;
      if (paths.size())
        paths.back().set_path_type(path_type);
    }
    else if ('C' == pat[pos]) {
      path_type = weave_path::CURVE;
      if (paths.size())
        paths.back().set_path_type(path_type);
    }

    // mirrors
    else if ('v' == pat[pos])
      ops.push_back(V);
    else if ('e' == pat[pos])
      ops.push_back(E);
    else if ('f' == pat[pos])
      ops.push_back(F);

    // rotations
    else if ('R' == pat[pos])
      reverse = !reverse;
    else if ('V' == pat[pos]) {
      if (!reverse) {
        ops.push_back(E);
        ops.push_back(F);
      }
      else {
        ops.push_back(F);
        ops.push_back(E);
      }
    }
    else if ('E' == pat[pos]) {
      ops.push_back(F);
      ops.push_back(V);
    }
    else if ('F' == pat[pos]) {
      if (!reverse) {
        ops.push_back(V);
        ops.push_back(E);
      }
      else {
        ops.push_back(E);
        ops.push_back(V);
      }
    }
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

/*
string weave_pattern::get_pattern() const
{
   string ret;
   for(unsigned int i=0; i<ops.size(); i++) {
      if(ops[i] == P)
         ret += "P";
         //ret += vtostr(paths[pt_idx++].get_path().get_point(), ",", 5);
      else if(ops[i] == V)
         ret += "v";
      else if(ops[i] == E)
         ret += "e";
      else if(ops[i] == F)
         ret += "f";
      else {
         ret = "X";
         break;
      }
   }
   ret += "t";
   if(start_faces==1)
      ret += 'l';
   if(start_faces==2)
      ret += 'r';
   if(start_faces==3)
      ret += 'b';

   return ret;
}
*/

class weave {
private:
  vector<weave_pattern> pats;

  Geometry meta;
  vector<vector<int>> nbrs;
  vector<Vec3d> vert_norms;

  bool find_nbrs();
  vector<Vec3d> point_on_face(int f_idx, const Vec3d &crds) const;
  void add_circuit(Geometry &wv, int start_idx, const weave_pattern &pat,
                   vector<bool> &seen) const;
  const vector<weave_pattern> &get_pats() const { return pats; }

public:
  bool set_geom(const Geometry &geom, double face_ht = 0.0);
  void add_pattern(const weave_pattern &pattern) { pats.push_back(pattern); }
  Status add_pattern(const string &pat);
  void make_weave(Geometry &wv) const;

  const Geometry &get_meta() const { return meta; }
};

bool weave::find_nbrs()
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
        nbrs[f][i] = -1; // only allow connection for two faces at an edge
      else {
        nbrs[f][i] =
            (ef_i->second[0] != (int)f) ? ef_i->second[0] : ef_i->second[1];
      }
    }
  return true;
}

inline vector<Vec3d> weave::point_on_face(int f_idx, const Vec3d &crds) const
{
  vector<Vec3d> ret(2);
  // point coordinates
  ret[0] = crds[weave_pattern::V] * meta.face_v(f_idx, weave_pattern::V) +
           crds[weave_pattern::E] * meta.face_v(f_idx, weave_pattern::E) +
           crds[weave_pattern::F] * meta.face_v(f_idx, weave_pattern::F);

  // point normal
  ret[1] =
      crds[weave_pattern::V] * vert_norms[meta.faces(f_idx, weave_pattern::V)] +
      crds[weave_pattern::E] * vert_norms[meta.faces(f_idx, weave_pattern::E)] +
      crds[weave_pattern::F] * vert_norms[meta.faces(f_idx, weave_pattern::F)];
  ret[1].to_unit();
  return ret;
}

void weave::add_circuit(Geometry &wv, int start_idx, const weave_pattern &pat,
                        vector<bool> &seen) const
{
  // Apply pattern until circuit completes
  vector<Vec3d> prev_pt;
  bool finish = false;
  int start_v_sz = wv.verts().size();
  int idx = start_idx;
  while (true) {
    seen[idx] = true;
    pat.start_op();
    while (pat.get_op() != weave_pattern::END) {
      // fprintf(stderr, "op=%d\n", pat.get_op());
      // fprintf(stderr, "nbrs[%d] = %d, %d, %d\n", idx,
      // nbrs[idx][0],nbrs[idx][1],nbrs[idx][2]);
      // fprintf(stderr, "nbrs[%d][%d] = %d, idx=%d, start_idx=%d\n", idx, op,
      // nbrs[idx][op], idx, start_idx);
      if (pat.get_op() == weave_pattern::P) {
        vector<Vec3d> pt = point_on_face(idx, pat.get_path().get_point());
        if (prev_pt.size())
          pat.get_path().add_points(&wv.raw_verts(), prev_pt, pt);
        prev_pt = pt;
      }
      else {
        idx = nbrs[idx][pat.get_op()]; // move to next triangle
        if (idx < 0) { // abandon: circuit tried to cross an open edge
          wv.raw_verts().resize(start_v_sz); // remove any added verts
          return;
        }
      }
      pat.next_op();
    }
    if (finish)
      break;
    if (idx == start_idx && prev_pt.size()) // circuit complete
      finish = true;
  }

  const int f_sz = wv.verts().size() - start_v_sz;
  vector<int> face(f_sz);
  for (int i = 0; i < f_sz; i++)
    face[i] = start_v_sz + i;
  wv.add_face(face);
}

static void reverse_odd_faces(Geometry &geom)
{
  const int f_sz = geom.faces().size();
  for (int i = 0; i < f_sz; i++)
    if (i % 2)
      reverse(geom.faces(i).begin(), geom.faces(i).end());
}

bool weave::set_geom(const Geometry &geom, double face_ht)
{
  make_meta(geom, meta, face_ht);
  find_nbrs();
  reverse_odd_faces(meta);
  vert_norms = meta.get_info().get_vert_norms();
  reverse_odd_faces(meta);
  return true;
}

Status weave::add_pattern(const string &pat)
{
  weave_pattern pattern;
  Status stat = pattern.set_pattern(pat);
  if (stat)
    add_pattern(pattern);

  return stat;
}

bool valid_start_face(int f, int start_faces)
{
  int is_left = f % 2;
  if (start_faces == 3)
    return true;
  else if (start_faces == 2 && is_left)
    return true;
  else if (start_faces == 1 && !is_left)
    return true;
  else
    return false;
}

void weave::make_weave(Geometry &wv) const
{
  // for(unsigned int i=0; i<nbrs.size(); i++)
  //   fprintf(stderr, "nbrs[%d] = %d, %d, %d\n", i,
  //   nbrs[i][0],nbrs[i][1],nbrs[i][2]);

  wv.clear_all();
  int faces_sz = meta.faces().size();
  for (const auto &pat : pats) {
    // fprintf(stderr, "pattern is '%s'\n", pats[p].get_pattern().c_str());
    vector<bool> seen(faces_sz, false);
    unsigned char start_faces = pat.get_start_faces();
    for (int i = 0; i < faces_sz; i++) {
      if (!seen[i] && valid_start_face(i, start_faces))
        add_circuit(wv, i, pat, seen);
    }
  }
}

class wv_opts : public ProgramOpts {
private:
public:
  bool add_meta;
  double face_ht;
  weave wv;
  bool use_default_pattern;
  string ifile;
  string ofile;

  wv_opts()
      : ProgramOpts("poly_weave"), add_meta(false), face_ht(0.0),
        use_default_pattern(true)
  {
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

// clang-format off
void wv_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Read a file in OFF format and make a weave following a specified pattern.\n"
"The polyhedron faces are divided by a 'meta' operation into triangles\n"
"each having vertices which are a vertex V, edge centre E and face centre F.\n"
"A start point is positioned on one of these triangles, the next point is\n"
"found by using the pattern to step between triangles, leading to a circuit.\n"
"Intermediate points may be added along the steps to help with the weave.\n"
"If input_file is not given the program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"  -p <pat>  weave pattern (default 'FEV'), a series of one or more paths.\n"
"            A path is:\n"
"               C,L - path of intermediate points is curve (default),line\n"
"                     (affects following paths, set as first char of pattern)\n"
"               Barycentric coords V,E,F of initial point (default 0,1,0)\n"
"               Optional intermediate points on the path - ':' followed\n"
"                  by 0 to 3 coordinates for 'up' (def: 0), 'side' (def: 0),\n"
"                  'along' (def: equal spacing)\n"
"               Series of characters to step between triangles\n"
"                  -     - stay on the same triangle\n"
"                  V,E,F - step two triangles rotating about V,E,F\n"
"                  R     - R reverse direction of following rotations\n"
"                  v,e,f - step over side opposite V,E,F\n"
"               t followed by l,r,b (def: tl), start circuits from 'left',\n"
"                  'right', both triangles\n"
"  -a        add the 'meta'-transformed base\n"
"  -f <ht>   lift the face centres by this height\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}
// clang-format on

void wv_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":ho:p:f:a")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'p':
      print_status_or_exit(wv.add_pattern(optarg), 'p');
      use_default_pattern = false;
      break;

    case 'a':
      add_meta = true;
      break;

    case 'f':
      print_status_or_exit(read_double(optarg, &face_ht), 'f');
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (use_default_pattern)
    wv.add_pattern("FEV");

  if (argc - optind > 1) {
    error("too many arguments");
    exit(1);
  }

  if (argc - optind == 1)
    ifile = argv[optind];
}

int main(int argc, char *argv[])
{
  wv_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);

  GeometryInfo info(geom);
  if (!info.is_orientable())
    opts.error("base polyhedron is not orientable");
  if (!info.is_oriented()) {
    opts.warning("base polyhedron is not oriented; it will be oriented.");
    geom.orient();
  }

  weave &wv = opts.wv;

  wv.set_geom(geom, opts.face_ht);
  Geometry wv_geom;
  wv.make_weave(wv_geom);
  if (opts.add_meta)
    wv_geom.append(wv.get_meta());

  opts.write_or_error(wv_geom, opts.ofile);

  return 0;
}
