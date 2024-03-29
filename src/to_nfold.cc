/*
   Copyright (c) 2013-2021, Adrian Rossiter

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
   Name: to_nfold.cc
   Description: Generalise an axial model by changing its rotational symmetry.
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <map>
#include <set>
#include <string>
#include <vector>

using std::map;
using std::pair;
using std::set;
using std::string;
using std::swap;
using std::vector;

using namespace anti;

struct PopVertex {
  int index;
  double scale_xy;
  double translate_z;

  Status read(const char *str);
};

Status PopVertex::read(const char *str)
{
  Status stat;
  Split parts(str, ",");
  if (parts.size() < 2 || parts.size() > 3)
    return Status::error(
        "must have either two or three comma separated values");
  if (!(stat = read_int(parts[0], &index)))
    return Status::error("index number: " + stat.msg());
  if (index < 0)
    return Status::error("index number: cannot be negative");
  if (!(stat = read_double(parts[1], &scale_xy)))
    return Status::error("scale-xy: " + stat.msg());
  translate_z = 0.0;
  if (parts.size() > 2 && !(stat = read_double(parts[2], &translate_z)))
    return Status::error("translate-z: " + stat.msg());
  return Status::ok();
}

class nfold_opts : public ProgramOpts {
private:
public:
  int num;
  int denom;
  vector<int> cross_ring_verts;
  vector<PopVertex> pop_vertices;

  string ifile;
  string ofile;

  nfold_opts() : ProgramOpts("to_nfold") {}

  void process_command_line(int argc, char **argv);
  void usage();
};

void nfold_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] fraction [input_file]

Generalise an axial model by changing its rotational symmetry. Read a model,
in OFF format, with an m-fold rotational axis on the z-axis, and create a
new model, generally non-planar, with the same relative connections, but
with an n-fold axis instead. fraction is given as n, or n/d (n and d
integers). Vertices of a face originally separated by x/m of a turn around
the z-axis will be separated by xd/n of a turn in the final model. The new
model will be symmetrically coloured using colours from the base model. If
input_file is not given the program reads from standard input.
%s is based on an idea by Bruce R. Gilson.

Options
%s
  -x <idxs> vertex index numbers, separated by commas, the rings including
            these vertices will be rotated 180 degrees before processing
            and rotated back afterwards
  -p <args> transform ring of vertices of base model, and suppress normal
            to_nfold processing. Arguments are two or three numbers
            separated by commas: vertex index number (specifies the ring),
            the scale-xy factor (ring radius), and an optional translation-z
            (ring height). Can be used multiple times. Index numbers are
            preserved.
  -o <file> write output to file (default: write to standard output)

)",
          prog_name(), prog_name(), help_ver_text);
}

void nfold_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hx:p:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'x':
      print_status_or_exit(read_int_list(optarg, cross_ring_verts, true), c);
      break;

    case 'p': {
      PopVertex pop_vert;
      print_status_or_exit(pop_vert.read(optarg), c);
      pop_vertices.push_back(pop_vert);
      break;
    }

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (argc - optind < 1) {
    error("not enough arguments");
    exit(1);
  }
  if (argc - optind > 2) {
    error("too many arguments");
    exit(1);
  }

  // read fraction
  char *p = strchr(argv[optind], '/');
  denom = 1;
  if (p != nullptr) {
    *p++ = '\0';
    print_status_or_exit(read_int(p, &denom), "fraction n/d, denominator");
  }

  print_status_or_exit(read_int(argv[optind], &num), "fraction n/d, numerator");
  if (num < 2)
    error("must be an integer 2 or greater", "fraction n/d, numerator");
  if (denom < 1)
    error("must be 1 or greater", "fraction n/d, denominator");
  if (denom % num == 0)
    error("numerator cannot divide the denominator", "fraction n/d, numerator");

  if (argc - optind == 2)
    ifile = argv[optind + 1];
}

// ----------------------------------------------------------------------

double to_range_0_2pi(double ang)
{
  return ang - 2 * M_PI * floor(ang / (2 * M_PI));
}

double to_range_plusminus_pi(double ang)
{
  double ret = to_range_0_2pi(ang);
  return (ret > M_PI) ? ret - 2 * M_PI : ret;
}

class cyc_vert {
public:
  double radius;
  double height;
  double cyc_ang;

  cyc_vert(double radius = 0.0, double height = 0.0, double angle = 0.0)
      : radius(radius), height(height), cyc_ang(angle)
  {
  }

  cyc_vert(Vec3d v)
  {
    radius = Vec3d(v[0], v[1], 0).len();
    height = v[2];
    cyc_ang = angle_around_axis(Vec3d::X, v, Vec3d::Z);
    if (fabs(cyc_ang) <
        0.1) // correct for precision problem around v near (0,0,1)
      cyc_ang = angle_around_axis(Vec3d::X - Vec3d::Y, v, Vec3d::Z) - M_PI / 4;
  }

  void to_range_0_2pi() { cyc_ang = ::to_range_0_2pi(cyc_ang); }
  void to_range_plusminus_pi() { cyc_ang = ::to_range_plusminus_pi(cyc_ang); }
  Vec3d to_Vec3d() const
  {
    return Vec3d(radius * cos(cyc_ang), radius * sin(cyc_ang), height);
  }

  // convert absolute angle to relative angle, range (-M_PI,M_PI]
  void angle_abs_to_rel(double from)
  {
    cyc_ang -= from;
    to_range_plusminus_pi();
  }
  bool on_axis(double eps) { return double_eq(radius, 0.0, eps); }

  void dump() const
  {
    Vec3d v = to_Vec3d();
    fprintf(stderr, "r=%g, h=%g, ang=%g (%g) -> (%g, %g, %g)\n", radius, height,
            cyc_ang, rad2deg(cyc_ang), v[0], v[1], v[2]);
  }
};

// ----------------------------------------------------------------------

struct cyc_vert_cmp {
  int n_fold;
  double eps;
  cyc_vert_cmp(int n_fold = 1, double eps = anti::epsilon)
      : n_fold(n_fold), eps(eps)
  {
  }
  bool operator()(const cyc_vert &from, const cyc_vert &to) const
  {
    int ret = double_compare(from.height, to.height, eps);
    if (!ret) {
      ret = double_compare(from.radius, to.radius, eps);
      if (!ret && !double_eq(from.radius, 0.0, eps))
        ret = double_compare(from.cyc_ang, to.cyc_ang, eps);
    }
    return ret < 0;
  }
};

class cyc_verts {
private:
  int n_fold;
  double eps;
  map<cyc_vert, int, cyc_vert_cmp> v2type;

  int get_sector(double ang) const;
  double first_sector(double ang) const;
  double sector_width() const { return 2 * M_PI / n_fold; }

public:
  cyc_verts(int n_fold, double eps) : n_fold(n_fold), eps(eps)
  {
    v2type = map<cyc_vert, int, cyc_vert_cmp>(cyc_vert_cmp(n_fold, eps));
  }

  pair<cyc_vert, int> get_vert_type(cyc_vert cyc_v) const;
  int get_vert_idx(cyc_vert cyc_v) const;
  void add_vert(cyc_vert cyc_v);
  void assign_types_in_order(vector<int> &ring_map);
  void to_vector(vector<cyc_vert> *v_types) const;
};

pair<cyc_vert, int> cyc_verts::get_vert_type(cyc_vert cyc_v) const
{
  cyc_v.to_range_plusminus_pi();
  double ang_orig = cyc_v.cyc_ang;
  cyc_v.cyc_ang = first_sector(ang_orig);
  const map<cyc_vert, int, cyc_vert_cmp>::const_iterator mi =
      v2type.find(cyc_v);
  if (mi == v2type.end())
    return std::make_pair(cyc_v, -1);
  else
    return *mi;
}

int cyc_verts::get_vert_idx(cyc_vert cyc_v) const
{
  pair<cyc_vert, int> v = get_vert_type(cyc_v);
  if (v.second < 0)
    return -1; // vertex not found
  int idx_off = floor((cyc_v.cyc_ang - v.first.cyc_ang) / sector_width() + 0.5);
  idx_off = (2 * n_fold + idx_off) % n_fold;
  if (idx_off < 0)
    idx_off += n_fold;
  return (v.second * n_fold) + idx_off;
}

void cyc_verts::add_vert(cyc_vert cyc_v)
{
  double ang_orig = cyc_v.cyc_ang;
  cyc_v.cyc_ang = first_sector(ang_orig);
  pair<cyc_vert, int> v = get_vert_type(cyc_v);
  if (v.second < 0) {
    int sz = v2type.size();
    v2type[cyc_v] = sz;
  }
}

int cyc_verts::get_sector(double ang) const
{
  ang = to_range_0_2pi(ang);
  if (first_sector(ang) == 0.0)
    ang += 0.5 * sector_width(); // sector boundary! Move to middle of sector
  return floor(ang / sector_width());
}

double cyc_verts::first_sector(double ang) const
{
  ang = to_range_0_2pi(ang);
  double first_ang = fmod(ang, sector_width());
  if (ang < 0 || double_eq(first_ang, sector_width(), eps))
    first_ang = 0.0;
  return first_ang;
}

void cyc_verts::assign_types_in_order(vector<int> &ring_map)
{
  ring_map.resize(v2type.size());
  int i = 0;
  map<cyc_vert, int, cyc_vert_cmp>::iterator mi;
  for (mi = v2type.begin(); mi != v2type.end(); ++mi) {
    ring_map[i] = mi->second;
    mi->second = i++;
  }
}

void cyc_verts::to_vector(vector<cyc_vert> *v_types) const
{
  v_types->resize(v2type.size());
  map<cyc_vert, int, cyc_vert_cmp>::const_iterator mi;
  for (mi = v2type.begin(); mi != v2type.end(); ++mi)
    (*v_types)[mi->second] = mi->first;
}

// ----------------------------------------------------------------------

enum {
  WARN_NONE = 0,
  WARN_EDGE180 = 1,
  WARN_EDGE180_CONSECUTIVE = 2,
  WARN_EDGE180_AMBIGUOUS = 4,
  WARN_NOCHAIN = 8,
  WARN_N_OVER_N = 16,
  WARN_N_OVER_N2 = 32,
  WARN_VERTEX_NOT_FOUND = 64
};

class cyc_chain {
public:
  double start_angle;
  vector<cyc_vert> cyc_vts;
  Color color;
  unsigned char warnings;
  cyc_chain(const Geometry &geom, const vector<int> &face, int len, Color col,
            double eps);
};

cyc_chain::cyc_chain(const Geometry &geom, const vector<int> &face, int len,
                     Color col, double eps)
    : color(col), warnings(WARN_NONE)
{
  // ensure that face starts with a non-axial vertex
  vector<int> face_rot = face;
  int fsz = face_rot.size();
  if (fsz == 0)
    return;
  for (int i = 0; i < fsz; i++) {
    cyc_vert cyc_v = geom.verts(face_rot[i]);
    if (!cyc_v.on_axis(eps)) {
      if (i)
        rotate(face_rot.begin(), face_rot.begin() + i, face_rot.end());
      break;
    }
  }

  double prev_abs_ang = cyc_vert(geom.verts(face_rot.back())).cyc_ang;
  cyc_vts.resize(fsz);
  vector<int> half_turns;
  for (int i = 0; i <= fsz; i++) {
    int idx = i % fsz;
    // fprintf(stderr, "face_rot[%d/%d] = %d", idx, fsz, face_rot[idx]);
    // geom.verts(face_rot[idx]).dump();
    cyc_vert cyc_v(geom.verts(face_rot[idx]));
    bool on_axis = cyc_v.on_axis(eps);
    double cur_abs_ang = on_axis ? prev_abs_ang : cyc_v.cyc_ang;
    if (i == 0) {
      start_angle = cyc_v.cyc_ang;
      // fprintf(stderr, "\nstart_angle=%.17g\n", rad2deg(start_angle));
    }
    else {
      if (on_axis)
        cyc_v.cyc_ang = 0;
      else {
        cyc_v.angle_abs_to_rel(prev_abs_ang);
        if (fsz > 2 && double_eq(fabs(cyc_v.cyc_ang), M_PI, eps))
          half_turns.push_back(idx);
      }
      cyc_vts[idx] = cyc_v;
    }

    prev_abs_ang = cur_abs_ang;
  }

  int hts_sz = half_turns.size();
  if (fsz > 2 && hts_sz)
    warnings |= WARN_EDGE180;

  for (int i = 0; i < hts_sz; i++) {
    int idx = half_turns[i];
    int idx_before = (idx - 1 + fsz) % fsz;
    int idx_after = (idx + 1) % fsz;
    if (hts_sz > 1 && std::abs(half_turns[i] - half_turns[(i + 1) % fsz]) == 1)
      warnings |= WARN_EDGE180_CONSECUTIVE;
    // Find sum of adjacent angles
    double test_ang = cyc_vts[idx_before].cyc_ang + cyc_vts[idx_after].cyc_ang;
    if (double_eq(test_ang, 0.0))
      warnings |= WARN_EDGE180_AMBIGUOUS;

    int sign = 1 - 2 * (test_ang > 0);
    cyc_vts[idx].cyc_ang = sign * M_PI;
  }
  cyc_vts.resize(len);
}

// ----------------------------------------------------------------------

class cyc_geom {
private:
  double eps;
  Geometry geom;
  Geometry chains;
  int from_n;
  vector<vector<set<int>>> sym_equivs;
  unsigned char warnings;
  string error_msg;

  int set_from_n(const Symmetry &sym);
  void add_chain_to_geom(cyc_chain &chain, Geometry *o_geom, int to_n, int to_d,
                         const cyc_verts &cyc_vs);

public:
  cyc_geom(double eps) : eps(eps) {}
  bool set_geom(const Geometry &geo);
  bool make_geom(Geometry *o_geom, int to_n, int to_d = 1,
                 const vector<int> cross_ring_verts = vector<int>());
  const Geometry &get_geom() const { return geom; }
  int get_ring_idx(int v_idx);
  int num_rings() { return sym_equivs[VERTS].size(); }
  Status pop_ring(int ring_idx, double scale_xy, double translate_z);
  const string &get_error_msg() const { return error_msg; }
  vector<string> get_input_warnings();
  vector<string> get_output_warnings();
};

bool cyc_geom::set_geom(const Geometry &geo)
{
  from_n = 0;

  Symmetry sym(geo);
  const set<SymmetryAxis> &axes = sym.get_axes();
  set<SymmetryAxis>::const_iterator ax;
  for (ax = axes.begin(); ax != axes.end(); ++ax) {
    if (double_eq(ax->get_axis()[2], 1.0)) { // aligned with z-axis
      from_n = ax->get_nfold();
      if (ax->get_sym_type() == Symmetry::S)
        from_n /= 2;
    }
  }

  if (!from_n) {
    error_msg = "model does not have cyclic symmetry around the z-axis";
    return false;
  }

  geom = geo;

  Symmetry sym_Cn(Symmetry::C, from_n);
  get_equiv_elems(geom, sym_Cn.get_trans(), &sym_equivs);

  return true;
}

// Is edge is horizontal with centre on the axis
bool is_axial_edge(const vector<int> &edge, const Geometry &geom, double eps)
{
  //  (x, y, z) == (-x, -y, z)
  const Vec3d &v0 = geom.verts(edge[0]);
  const Vec3d &v1 = geom.verts(edge[1]);
  return double_eq(v0[0], -v1[0], eps) && double_eq(v0[1], -v1[1], eps) &&
         double_eq(v0[2], v1[2], eps);
}

// Create a rotated face by stepping through vertex indices
vector<int> get_next_face(vector<int> face, int step, int to_n)
{
  vector<int> f(face.size());
  for (unsigned int i = 0; i < face.size(); i++) {
    f[i] = (face[i] / to_n) * to_n + (face[i] + step) % to_n;
  }
  return f;
}

vector<string> cyc_geom::get_input_warnings()
{
  vector<string> msgs;
  if (warnings & WARN_EDGE180)
    msgs.push_back("includes 180 degree edges");
  if (warnings & WARN_EDGE180_CONSECUTIVE)
    msgs.push_back("includes consecutive 180 degree edges");
  if (warnings & WARN_EDGE180_AMBIGUOUS)
    msgs.push_back("includes a 180 degree edge and its sign cannot be "
                   "unambiguously assigned");
  if (warnings & WARN_NOCHAIN)
    msgs.push_back("includes an edge chain that does not repeat "
                   "symmetrically (invalid result)");
  return msgs;
}

vector<string> cyc_geom::get_output_warnings()
{
  vector<string> msgs;
  if (warnings & WARN_N_OVER_N)
    msgs.push_back("includes a face with a vertex adjacent to itself "
                   "(invalid N/N result?)");
  if (warnings & WARN_N_OVER_N2)
    msgs.push_back("includes a non-axial two-vertex face"
                   "(invalid N/N result?)");
  if (warnings & WARN_VERTEX_NOT_FOUND)
    msgs.push_back("could not determine all vertex positions"
                   "(invalid result)");
  return msgs;
}

void cyc_geom::add_chain_to_geom(cyc_chain &chain, Geometry *o_geom, int to_n,
                                 int to_d, const cyc_verts &cyc_vs)
{
  // fprintf(stderr, "\nADDING CHAIN\n");

  double first_sector_start_angle = fmod(chain.start_angle, 2 * M_PI / from_n);
  double new_start_angle = first_sector_start_angle * from_n * to_d / to_n;
  vector<int> face;
  // fprintf(stderr, "\nstart_angle = %g\n", chain.start_angle);
  // fprintf(stderr, "new_start_angle = %g\n", new_start_angle);
  double prev_angle = new_start_angle;
  bool cycle_completed = false;
  bool all_vertices_valid = true;
  int num_chains_in_face = 0;
  for (int n = 0; n < to_n; n++) {
    // fprintf(stderr, "n = %d\n", n);
    for (unsigned int i = 0; i < chain.cyc_vts.size(); i++) {
      cyc_vert cyc_v = chain.cyc_vts[i];
      if (n == 0 && i == 0)
        cyc_v.cyc_ang = 0.0; // set offset from start angle, not last vert
      // fprintf(stderr, "i=%d, a = %g * pi\n", i, cyc_v.cyc_ang/(M_PI));
      cyc_v.cyc_ang = prev_angle + cyc_v.cyc_ang * from_n * to_d / to_n;
      // cyc_v.dump();
      prev_angle = cyc_v.cyc_ang;
      int idx = cyc_vs.get_vert_idx(cyc_v);
      // fprintf(stderr, "idx = %d\n\n", idx);
      if (i == 0 && face.size() && idx == face[0]) {
        cycle_completed = true;
        break;
      }
      if (face.size() && face.back() == idx)
        warnings |= WARN_N_OVER_N;
      if (idx == -1) {
        warnings |= WARN_VERTEX_NOT_FOUND;
        all_vertices_valid = false;
      }
      face.push_back(idx);
    }
    if (cycle_completed)
      break;
    else
      num_chains_in_face++;
  }

  if (face.size() == 2 && !is_axial_edge(face, *o_geom, eps))
    warnings |= WARN_N_OVER_N2;

  if (all_vertices_valid && num_chains_in_face) {
    for (int i = 0; i < to_n / num_chains_in_face; i++)
      o_geom->add_face(get_next_face(face, i, to_n), chain.color);
  }
}

int cyc_geom::get_ring_idx(int v_idx)
{
  int ring_idx = -1;
  const vector<set<int>> &v_equivs = sym_equivs[VERTS];
  for (int r_idx = 0; r_idx < (int)v_equivs.size(); r_idx++) {
    if (v_equivs[r_idx].count(v_idx) > 0) {
      ring_idx = r_idx;
      break;
    }
  }
  return ring_idx;
}

Status cyc_geom::pop_ring(int ring_idx, double scale_xy, double translate_z)
{
  if (ring_idx < 0 || ring_idx >= (int)sym_equivs[VERTS].size())
    return Status::error(msg_str("invalid ring number %d", ring_idx));

  for (auto v_idx : sym_equivs[VERTS][ring_idx]) {
    geom.verts(v_idx)[0] *= scale_xy;
    geom.verts(v_idx)[1] *= scale_xy;
    geom.verts(v_idx)[2] += translate_z;
  }

  return Status::ok();
}

bool cyc_geom::make_geom(Geometry *o_geom, int to_n, int to_d,
                         const vector<int> cross_ring_verts)
{
  warnings = WARN_NONE;

  // Find the rings that will be crossed
  const vector<set<int>> &v_equivs = sym_equivs[0];
  vector<bool> cross_rings(v_equivs.size(), false);
  for (int v_idx : cross_ring_verts) {
    int ring_idx = get_ring_idx(v_idx);
    if (ring_idx >= 0) {
      if (cross_rings[ring_idx]) {
        error_msg = msg_str("option -x: index number %d is on a ring that has "
                            "already been crossed",
                            v_idx);
        return false;
      }
      else
        cross_rings[ring_idx] = true;
    }
  }

  for (int ring_idx = 0; ring_idx < (int)v_equivs.size(); ring_idx++) {
    if (cross_rings[ring_idx]) {
      // rotate all the vertices in the ring by 180
      for (auto v_idx : v_equivs[ring_idx]) { // cross any crossed ring vertices
        geom.verts(v_idx)[0] *= -1;
        geom.verts(v_idx)[1] *= -1;
      }
    }
  }

  // Set up vertices, as types
  cyc_verts cyc_vs(to_n, sym_eps);        // vertex types in final model.
  cyc_verts cyc_vs_from(from_n, sym_eps); // vertex types in original model
  for (auto &i : sym_equivs[0]) {
    // get equiv vertex in first from_n sector
    cyc_vert cyc_v = cyc_vs_from.get_vert_type(geom.verts(*i.begin())).first;
    // 'scale' to corresponding position in first to_n sector
    cyc_v.cyc_ang *= (double)from_n * to_d / to_n;
    cyc_vs.add_vert(cyc_v);
  }
  vector<int> ring_map;
  cyc_vs.assign_types_in_order(ring_map); // number the vertices nicely

  // Set up final vertices
  vector<cyc_vert> vs;
  cyc_vs.to_vector(&vs);
  for (int i = 0; i < (int)vs.size(); i++) {
    for (int step = 0; step < to_n; step++) {
      cyc_vert v_rot = vs[i];
      v_rot.cyc_ang += step * 2 * M_PI / to_n;
      auto vert = v_rot.to_Vec3d();
      if (cross_rings[ring_map[i]]) { // uncross any crossed ring vertices
        vert[0] *= -1;
        vert[1] *= -1;
      }
      o_geom->add_vert(vert);
    }
  }

  vector<cyc_chain> chains;
  const vector<set<int>> &f_equivs = sym_equivs[FACES];
  for (const auto &f_equiv : f_equivs) {
    const int f_idx = *f_equiv.begin();
    // fprintf(stderr, "\n\nNEW CHAIN face_idx=%d\n", f_idx);
    // fprintf(stderr, "Num faces=%d, from_n=%d\n", (int)f_equivs[i].size(),
    // from_n);
    vector<int> face = geom.faces(f_idx);
    int fsz = face.size();

    // Find the minimum length of chain of edges that make a repeat unit
    int chain_sz = fsz * f_equiv.size() / from_n;
    if (chain_sz < fsz) {
      double steps =
          angle_around_axis(geom.verts(face[0]),
                            geom.verts(face[chain_sz % fsz]), Vec3d::Z) *
          from_n / (2 * M_PI);
      double diff = steps - round(steps);
      if (double_ne(diff, 0.0, eps))
        warnings |= WARN_NOCHAIN;
    }

    cyc_chain cyc_c(geom, face, chain_sz, geom.colors(FACES).get(f_idx), eps);
    chains.push_back(cyc_c);
    warnings |= cyc_c.warnings;
  }

  // infer "horizontal" axial digons as a horizontal edge shared by two
  // equivalent faces
  map<int, int> f2equiv;
  for (unsigned int i = 0; i < f_equivs.size(); i++)
    for (auto si = f_equivs[i].begin(); si != f_equivs[i].end(); ++si)
      f2equiv[*si] = i;

  auto e2f = geom.get_edge_face_pairs(false);
  for (const auto &kp : e2f) {
    if (is_axial_edge(kp.first, geom, eps) && kp.second.size() == 2 &&
        f2equiv[kp.second[0]] == f2equiv[kp.second[1]]) {
      cyc_chain cyc_c(geom, kp.first, 1, Color(), eps);
      chains.push_back(cyc_c);
    }
  }

  for (auto &chain : chains)
    add_chain_to_geom(chain, o_geom, to_n, to_d, cyc_vs);
  return true;
}

// ----------------------------------------------------------------------

int main(int argc, char *argv[])
{
  const double eps = 10e-8; // hardcoded
  nfold_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);

  // merge geometry and check cross rings are valid
  int num_verts_orig = geom.verts().size();
  merge_coincident_elements(geom, "vef", 0, eps);
  if (opts.cross_ring_verts.size() != 0) {
    if ((int)geom.verts().size() != num_verts_orig)
      opts.error("cannot use this option with input model that has coincident "
                 "vertices (merge vertices first)",
                 'x');
    for (int idx : opts.cross_ring_verts) {
      if (idx >= num_verts_orig)
        opts.error(msg_str("index number %d, out of range", idx), 'x');
    }
  }

  for (size_t i = 0; i < geom.edges().size(); i++)
    geom.add_face(geom.edges(i), geom.colors(EDGES).get(i));

  cyc_geom cyc(eps);
  if (!cyc.set_geom(geom))
    opts.error(cyc.get_error_msg(), "input geometry");

  if (opts.pop_vertices.size() > 0) {
    // pop vertices (don't apply nfold processing)
    vector<bool> rings_seen(cyc.num_rings());
    for (auto pop_vert : opts.pop_vertices) {
      int ring_idx = cyc.get_ring_idx(pop_vert.index);
      if (ring_idx < 0)
        opts.error(
            msg_str("vertex index %d is not on a vertex ring", pop_vert.index),
            'p');
      if (rings_seen[ring_idx])
        opts.error(msg_str("vertex index %d is on a vertex ring that has "
                           "already been processed",
                           pop_vert.index),
                   'p');
      cyc.pop_ring(ring_idx, pop_vert.scale_xy, pop_vert.translate_z);
    }
    opts.write_or_error(cyc.get_geom(), opts.ofile);
  }
  else {
    // apply nfold processing
    Geometry o_geom;
    if (!cyc.make_geom(&o_geom, opts.num, opts.denom, opts.cross_ring_verts))
      opts.error(cyc.get_error_msg(), "input geometry");

    o_geom.orient();

    vector<string> msgs = cyc.get_input_warnings();
    for (auto &msg : msgs)
      opts.warning(msg, "input geometry");

    msgs = cyc.get_output_warnings();
    for (auto &msg : msgs)
      opts.warning(msg, "output geometry");

    merge_coincident_elements(o_geom, "vef", 0, eps);
    opts.write_or_error(o_geom, opts.ofile);
  }

  return 0;
}
