/*
   Copyright (c) 2019, Adrian Rossiter

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
   Name: rotegrity.cc
   Description: make rotegrity models
   Project: Antiprism - http://www.antiprism.com
*/

#include <cmath>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>

#include "../base/antiprism.h"

using std::string;
using std::vector;

using namespace anti;

class iter_params {
public:
  int num_iters;
  int num_iters_status;
  int sig_digits;
  FILE *rep_file;
  iter_params()
      : num_iters(10000), num_iters_status(1000), sig_digits(15),
        rep_file(stderr)
  {
  }

  double get_test_val() { return pow(10, -sig_digits); }
  bool quiet() { return rep_file == nullptr; }
  bool check_status(int n)
  {
    return n == num_iters || n == 1 ||
           (num_iters_status > 0 && n % num_iters_status == 0);
  }
  bool checking_status() { return num_iters_status > 0; };
  bool print_progress_dot(int n)
  {
    if (!quiet() && num_iters_status >= 10)
      return (n % num_iters) % (num_iters_status / 10) == 0;
    else
      return false;
  }
};

class rot_opts : public ProgramOpts {
public:
  iter_params it_params;
  double adjust_fact = 98; // generally efficient, chosen from testing
  double end_fraction = 1 / 3.0;
  bool already_twisted = false;
  int method = 0;      // unset
  int output_type = 1; // full
  int col_type = 0;    // unset
  bool twist_program_replacement = false;
  Coloring clrngs[3];
  string ifile;
  string ofile;

  rot_opts() : ProgramOpts("rotegrity") { read_colorings(clrngs, "spread"); }

  void process_command_line(int argc, char **argv);
  void usage();
};

// clang-format off
void rot_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Read a file in OFF format containing a roughly spherical polyhedron, or\n"
"previously twisted model, and try to convert into a rotegrity. Units\n"
"keep original edge colours. If input_file is not given the program reads\n"
"from standard input.\n"
"\n"
"Options\n"
"%s"
"  -f <frac> fraction of length for end sections (default: 1/3)\n"
"  -t        input model is already twisted (produced by this program or\n"
"            'twist' program)\n"
"  -M <mthd> method of conversion from base model - twist, double, join\n"
"            (default: t)\n"
"  -O <type> output type for units: full (face), rotegrity (3 short struts),\n"
"            nexorade (long strut) (default: full). Only 'full' output can \n"
"            be used as input with option -t\n"
"  -m <maps> a comma separated list of colour maps used to transform colour\n"
"            indexes (default: rand), a part consisting of letters from\n"
"            v, e, f, selects the element types to apply the map list to\n"
"            (default 'vef').\n"
"  -c <type> colouring type: edge (base model edges), symmetry, none\n"
"            (default: edge(\n"
"  -n <itrs> number of iterations (default 10000)\n"
"  -s <perc> percentage to adjust corrections on iteration (default: %.0f)\n"
"  -l <lim>  minimum change of vertex distance to terminate, as negative\n"
"            exponent (default: %d giving %.0e)\n"
"  -z <n>    status checking and reporting every n iterations, -1 for no\n"
"            status (default: 1000)\n"
"  -q        quiet, do not print status messages\n"
"  -T        reproduce output of former 'twist' program (see Notes), -f is\n"
"            twist factor, -O is output type, unused options silently ignored\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text,
   adjust_fact, it_params.sig_digits, it_params.get_test_val() );
}
// clang-format on

void rot_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;
  vector<double> nums;
  string arg_id;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hf:tM:O:c:m:n:s:l:z:qTo:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'n':
      print_status_or_exit(read_int(optarg, &it_params.num_iters), c);
      if (it_params.num_iters < 0)
        error("number of iterations must be greater than 0", c);
      break;

    case 'z':
      print_status_or_exit(read_int(optarg, &it_params.num_iters_status), c);
      if (it_params.num_iters_status < -1)
        error("number of iterations must be -1 or greater", c);
      break;

    case 's':
      print_status_or_exit(read_double(optarg, &adjust_fact), c);
      if (adjust_fact < 0 || adjust_fact > 100)
        warning("not in range 0 to 100", c);
      break;

    case 'f':
      print_status_or_exit(read_double(optarg, &end_fraction), c);
      if (end_fraction < 0.0 || end_fraction > 0.5)
        warning("not in range 0 to 0.5", c);
      break;

    case 't':
      already_twisted = true;
      break;

    case 'M':
      print_status_or_exit(get_arg_id(optarg, &arg_id, "twist=1|doube=2|join=3",
                                      argmatch_default),
                           c);
      method = atoi(arg_id.c_str());
      break;

    case 'O':
      print_status_or_exit(get_arg_id(optarg, &arg_id,
                                      "full=1|rotegrity=2|nexorade=3",
                                      argmatch_default),
                           c);
      output_type = atoi(arg_id.c_str());
      break;

    case 'c':
      print_status_or_exit(get_arg_id(optarg, &arg_id,
                                      "edge=1|symmetry=2|none=3",
                                      argmatch_default),
                           c);
      col_type = atoi(arg_id.c_str());
      break;

    case 'm':
      print_status_or_exit(read_colorings(clrngs, optarg), c);
      break;

    case 'l':
      print_status_or_exit(read_int(optarg, &it_params.sig_digits), c);
      if (it_params.sig_digits < 0) {
        warning("termination limit is negative, and so ignored", c);
      }
      if (it_params.sig_digits > DEF_SIG_DGTS) {
        warning("termination limit is very small, may not be attainable", c);
      }
      break;

    case 'q':
      it_params.rep_file = nullptr;
      break;

    case 'T':
      twist_program_replacement = true;
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (argc - optind > 1)
    error("too many arguments");

  if (argc - optind == 1)
    ifile = argv[optind];

  if (already_twisted) {
    if (method)
      error("cannot use option -M with previously twisted input", 't');
    if (col_type)
      error("cannot use option -c to colour previously twisted input", 't');
  }

  if (method == 0)
    method = 1; // default: 'twist'

  if (col_type == 0)
    col_type = 1; // default: 'edge'
}

Status check_polyhedron_model(GeometryInfo &info)
{
  if (!info.is_polyhedron())
    return Status::error("polyhedron model: not connected like polyhedron");
  if (!info.is_orientable())
    return Status::error("polyhedron model: not orientable");
  if (info.genus() != 0)
    return Status::error("polyhedron model: not connected like a sphere");

  return Status::ok();
}

Status check_twist_model(const Geometry &geom)
{
  vector<vector<int>> verts2faces(geom.verts().size());
  for (unsigned int f = 0; f < geom.faces().size(); f++) {
    const auto &face = geom.faces(f);
    if (face.size() != 4)
      return Status::error(
          msg_str("twist model: face %d does not have four vertices", f));
    for (int i = 0; i < 4; i++)
      verts2faces[face[i]].push_back(f);
  }
  for (unsigned int i = 0; i < geom.verts().size(); i++) {
    if (verts2faces[i].size() != 2)
      return Status::error(
          msg_str("twist model: vertex %d does not lie on two faces", i));
  }
  return Status::ok();
}

void make_twist(Geometry &tw_geom, const Geometry &geom, const Symmetry &sym,
                int method, double ratio, int col_type)
{
  double wt = 2 * ratio;
  string pattern;
  if (method == 1)
    pattern = msg_str("[%gV%gE]0V0E0fe", 1 - wt, wt);
  else if (method == 2)
    pattern = msg_str("[%gV%gE,%gE%gF]0V0_1F1", 1 - wt, wt, wt, 1 - wt);
  else if (method == 3)
    pattern = msg_str("[2VE-0.3F,2F3E,3FV]0V0_1F2_1");

  // Model suitable for 'edge' colouring type: col_type == 0
  wythoff_make_tiling(tw_geom, geom, pattern, true, false,
                      Tiling::ColoringType::associated_element);
  if (col_type == 2) { // 'symmetry' colouring
    vector<vector<std::set<int>>> sym_equivs;
    get_equiv_elems(tw_geom, sym.get_trans(), &sym_equivs);
    Coloring coloring(&tw_geom);
    coloring.f_sets(sym_equivs[2], false);
  }
  else if (col_type == 3) // 'none' colouring
    tw_geom.clear_cols();

  int num_faces = tw_geom.faces().size();
  for (int i = 0; i < num_faces; i++) {
    auto face = tw_geom.faces(i);
    int unit_fsz = 4;
    tw_geom.faces(i) = vector<int>(face.begin(), face.begin() + unit_fsz);
    if (method == 2)
      tw_geom.add_face(vector<int>(face.begin() + unit_fsz, face.end()),
                       tw_geom.colors(FACES).get(i));
    else if (method == 3) {
      tw_geom.add_face(
          vector<int>(face.begin() + 5, face.begin() + 5 + unit_fsz),
          tw_geom.colors(FACES).get(i));
      tw_geom.add_face({face[4], face[3], face[8], face[9]},
                       tw_geom.colors(FACES).get(i));
    }
  }

  for (unsigned int i = 0; i < tw_geom.faces().size(); i++) {
    auto &face = tw_geom.faces(i);
    std::rotate(face.begin(), face.begin() + 3, face.end());
    Color col = tw_geom.colors(FACES).get(i);
    tw_geom.add_edge(face[0], face[1], Color::invisible);
    for (int i = 1; i < 4; i++)
      tw_geom.add_edge(face[i], face[(i + 1) % 4], col);
    for (int i = 0; i < 2; i++)
      tw_geom.colors(VERTS).set(face[i], col);
    tw_geom.colors(FACES).set(i, Color::invisible);
  }
}

inline double adjust_vert_to_target(Vec3d &v_new, const Vec3d &vert,
                                    const Vec3d target, double factor)
{
  auto vec_diff = target - vert;
  auto diff = vec_diff.len();
  v_new = (vert + vec_diff * factor).unit();
  return diff;
}

void make_rotegrity(Geometry &rotegrity, const Symmetry &sym,
                    double end_fraction, iter_params it_params, double factor)
{
  const auto test_val = it_params.get_test_val();
  // Read and write to same model (Defered update is slower)
  SymmetricUpdater sym_updater(rotegrity, sym, false);
  const auto &face_orbits = sym_updater.get_equiv_sets(FACES);

  double max_diff = 0.0;
  int cnt;
  for (cnt = 1; cnt <= it_params.num_iters; cnt++) {
    max_diff = 0.0;
    for (const auto &face_orbit : face_orbits) {
      const auto &face =
          sym_updater.get_geom_reading().faces(*face_orbit.begin());
      // Start: make sure vertices are current
      for (int i = 0; i < 4; i++)
        sym_updater.update_from_principal_vertex(face[i]);

      // Face has strut edge (long edge) first
      // Vertex index numbers and coordinates
      const auto &v_idx_end0 = face[0]; // first end
      auto &v_end0 = sym_updater.get_geom_reading().verts(v_idx_end0);
      const auto &v_idx_end1 = face[1]; // second end
      auto &v_end1 = sym_updater.get_geom_reading().verts(v_idx_end1);
      const auto &v_idx_mid1 = face[2]; // adjacent to second end
      auto &v_mid1 = sym_updater.get_geom_reading().verts(v_idx_mid1);
      const auto &v_idx_mid0 = face[3]; // adjacent to first end
      auto &v_mid0 = sym_updater.get_geom_reading().verts(v_idx_mid0);

      // Get normal to strut plane
      auto norm = vcross(v_end0, v_end1);
      auto ang = angle_around_axis(v_end0, v_end1, norm);
      auto rot0 = Trans3d::rotate(norm, ang * end_fraction);
      auto target_v_mid0 = rot0 * v_end0;
      auto rot1 = Trans3d::rotate(norm, -ang * end_fraction);
      auto target_v_mid1 = rot1 * v_end1;

      Vec3d v_new;
      double diff = adjust_vert_to_target(v_new, v_mid0, target_v_mid0, factor);
      // fprintf(stderr, "diff=%10f", diff);
      if (diff > max_diff)
        max_diff = diff;
      sym_updater.update_principal_vertex(v_idx_mid0, v_new);

      diff = adjust_vert_to_target(v_new, v_mid1, target_v_mid1, factor);
      // fprintf(stderr, "\t\t%10f\n", diff);
      if (diff > max_diff)
        max_diff = diff;
      sym_updater.update_principal_vertex(v_idx_mid1, v_new);
    }
    sym_updater.prepare_for_next_iteration();

    if (max_diff <= test_val)
      break;

    if (!it_params.quiet() && it_params.check_status(cnt))
      fprintf(it_params.rep_file, "\niter:%-15d max_diff:%17.15f ", cnt,
              max_diff);
    else if (it_params.print_progress_dot(cnt))
      fprintf(it_params.rep_file, ".");
  }
  if (!it_params.quiet() && it_params.checking_status())
    fprintf(it_params.rep_file, "\nFinal:\niter:%-15d max_diff:%17.15f\n",
            (cnt < it_params.num_iters) ? cnt : it_params.num_iters, max_diff);

  rotegrity = sym_updater.get_geom_final();
}

void to_output_type(Geometry &geom, int type)
{
  if (type == 1)
    return;

  map<vector<int>, Color> edge2col;
  for (unsigned int i = 0; i < geom.edges().size(); i++)
    edge2col[geom.edges(i)] = geom.colors(EDGES).get(i);
  geom.clear(EDGES);

  for (unsigned int i = 0; i < geom.faces().size(); i++) {
    auto &face = geom.faces(i);
    if (type == 2) {
      for (int i = 1; i < 4; i++) {
        auto edge = make_edge(face[i], face[(i + 1) % 4]);
        geom.add_edge_raw(edge, edge2col[edge]);
      }
    }
    else if (type == 3) {
      auto edge = make_edge(face[0], face[1]);
      auto edge_col = make_edge(face[2], face[3]);
      geom.add_edge_raw(edge, edge2col[edge_col]);
    }
  }
  geom.clear(FACES);
}

void make_twist_original(const Geometry &geom, rot_opts &opts);

int main(int argc, char *argv[])
{
  rot_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);

  if (!geom.faces().size())
    opts.error("input file contains no faces");

  if (opts.twist_program_replacement) {
    // Avoid impact of this option on rest of program
    make_twist_original(geom, opts);
    return 0;
  }

  geom.transform(Trans3d::translate(-geom.centroid()));
  Symmetry sym(geom);
  Geometry rotegrity;
  if (opts.already_twisted) {
    opts.print_status_or_exit(check_twist_model(geom));
    rotegrity = geom;
  }
  else {
    GeometryInfo info(geom);
    opts.print_status_or_exit(check_polyhedron_model(info));
    if (!info.is_oriented()) {
      opts.warning("polyhedron model: orienting model as not oriented");
      geom.orient(1);
    }
    make_twist(rotegrity, geom, sym, opts.method, opts.end_fraction,
               opts.col_type);
    project_onto_sphere(rotegrity);
  }

  make_rotegrity(rotegrity, sym.get_max_direct_sub_sym(), opts.end_fraction,
                 opts.it_params, opts.adjust_fact / 100);

  to_output_type(rotegrity, opts.output_type);
  for (int i = 0; i < 3; i++)
    opts.clrngs[i].set_geom(&rotegrity);
  opts.clrngs[VERTS].v_apply_cmap();
  opts.clrngs[EDGES].e_apply_cmap();
  opts.clrngs[FACES].f_apply_cmap();

  opts.write_or_error(rotegrity, opts.ofile);

  return 0;
}

// -----------------------------------------------------------------------
// Original twist program code follows. Likely to be removed in the future.

Geometry twist_original(const Geometry &poly, Geometry &dual, double twist_val,
                        Vec3d centre, bool struts_only)
{
  auto edges = poly.get_edge_face_pairs(true);
  map<vector<int>, vector<int>>::iterator mi, mi_next;
  Geometry twist;

  Vec3d v0, v1, vm; // strut ends and middle
  double ratio;
  Vec3d v0p, v1p, v0d, v1d; // vertices at ends 1 and 2 in poly and dual
  map<vector<int>, vector<int>> attach_edges;
  vector<int> tv_map;
  map<vector<int>, std::pair<Vec3d, Vec3d>> twist_edges;
  map<vector<int>, std::pair<Vec3d, Vec3d>>::iterator mi_tw, mi_tw2;
  for (mi = edges.begin(); mi != edges.end(); mi++) {
    v0p = poly.verts(mi->first[0]);
    v1p = poly.verts(mi->first[1]);
    v0d = dual.verts(mi->second[1]);
    v1d = dual.verts(mi->second[0]);
    // don't need this
    // ratio = (v1p - v0p).len() / (v1d - v0d).len();
    // v0d = centre + (v0d - centre)*ratio;
    // v1d = centre + (v1d - centre)*ratio;
    Vec3d norm = vcross(v1p - v0p, v1d - v0d);
    Trans3d trans = Trans3d::translate(-0.5 * (v0p + v1p));
    trans = Trans3d::rotate(norm, twist_val * M_PI / 2) * trans;
    trans = Trans3d::translate(0.5 * (v0p + v1p)) * trans;
    trans =
        Trans3d::translate((v0d + v1d - v0p - v1p) * 0.5 * twist_val) * trans;

    v0 = trans * v0p;
    v1 = trans * v1p;

    // v0 = v0p + (v0d - v0p)*twist_val;
    // v1 = v1p + (v1d - v1p)*twist_val;

    twist_edges[mi->first] = std::pair<Vec3d, Vec3d>(v0, v1);
    vector<int> ve = mi->first;
    if (ve[1] < ve[0])
      std::swap(ve[1], ve[0]);
    attach_edges[ve] = vector<int>();
  }

  double edge_len = (v1p - v0p).len();
  vector<int>::const_iterator vi;
  for (mi_tw = twist_edges.begin(); mi_tw != twist_edges.end(); mi_tw++) {
    mi = edges.find(mi_tw->first);
    twist.add_face(vector<int>());
    for (int i = 0; i < 2; i++) {
      vector<int> e_next(2);
      e_next[0] = mi_tw->first[(i + 1) % 2];
      const vector<int> &face = poly.faces(mi->second[i]);
      vi = find(face.begin(), face.end(), e_next[0]);
      e_next[1] = ++vi != face.end() ? *vi : face.front();
      // fprintf(stderr, "twist_edge/next (%d, %d)/(%d, %d)\n", mi_tw->first[0],
      // mi_tw->first[1], e_next[0], e_next[1]);
      if (e_next[1] < e_next[0])
        std::swap(e_next[1], e_next[0]);
      mi_tw2 = twist_edges.find(e_next);
      attach_edges[e_next].push_back(twist.verts().size());
      int where;
      Vec3d v = line_plane_intersect(centre, mi_tw2->second.first,
                                     mi_tw2->second.second, mi_tw->second.first,
                                     mi_tw->second.second, &where);
      if (!v.is_set())
        return Geometry(); // Failure in construction

      tv_map.push_back(mi_tw->first[i]);
      twist.add_vert(v);
      twist.faces(twist.faces().size() - 1).push_back(twist.verts().size() - 1);
    }

    Vec3d &v0 = twist.verts(twist.verts().size() - 2);
    Vec3d &v1 = twist.verts(twist.verts().size() - 1);
    ratio = edge_len / (v1 - v0).len();
    v0 = centre + (v0 - centre) * ratio;
    v1 = centre + (v1 - centre) * ratio;
  }

  if (struts_only) {
    for (unsigned int i = 0; i < twist.faces().size(); i++)
      twist.add_edge(make_edge(twist.faces(i, 0), twist.faces(i, 1)), 0);
    twist.clear(FACES);
  }
  else {
    for (unsigned int i = 0; i < twist.faces().size(); i++) {
      vector<int> edge(2);
      edge[0] = tv_map[twist.faces(i, 0)];
      edge[1] = tv_map[twist.faces(i, 1)];
      if (edge[1] < edge[0])
        std::swap(edge[1], edge[0]);
      // fprintf(stderr, "edge = (%d, %d) edge.size = %d\n", edge[0], edge[1],
      // (int)edge.size());
      vector<int> v_idxs = attach_edges.find(edge)->second;
      Vec3d P0 = twist.face_v(i, 0);
      Vec3d edge_vec = twist.face_v(i, 1) - P0;
      if (vdot(edge_vec, twist.verts(v_idxs[0]) - P0) <
          vdot(edge_vec, twist.verts(v_idxs[1]) - P0))
        std::swap(v_idxs[0], v_idxs[1]);
      twist.faces(i).push_back(v_idxs[0]);
      twist.faces(i).push_back(v_idxs[1]);

      for (int j = 0; j < 4; j++)
        twist.add_edge(
            make_edge(twist.faces(i, j), twist.faces(i, (j + 1) % 4)),
            Color(bool(j)));
    }
  }

  return twist;
}

void make_twist_original(const Geometry &geom, rot_opts &opts)
{
  Vec3d origin(0, 0, 0);
  Geometry dual;
  get_dual(dual, geom, 1, origin);

  vector<int> invalid_verts;
  for (unsigned int i = 0; i < dual.verts().size(); i++) {
    if (!dual.verts(i).is_set()) {
      dual.del(VERTS, i);
      int idx = invalid_verts.size() + i;
      invalid_verts.push_back(idx);
      i--;
    }
  }
  if (invalid_verts.size()) {
    string msg(
        "removed invalid vertices (and associated faces) with indices - ");
    char errmsg[MSG_SZ];
    for (unsigned int i = 0; i < invalid_verts.size() - 1; i++) {
      snprintf(errmsg, MSG_SZ, "%d,", invalid_verts[i]);
      msg += string(errmsg);
    }
    snprintf(errmsg, MSG_SZ, "%d", invalid_verts.back());
    msg += string(errmsg);
    opts.warning(msg);
  }

  Geometry twisted =
      twist_original(geom, dual, opts.end_fraction, origin, false);
  if (!twisted.is_set())
    opts.error("failed to construct model", 'T');

  to_output_type(twisted, opts.output_type);
  opts.write_or_error(twisted, opts.ofile);
}
