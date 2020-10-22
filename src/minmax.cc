/*
   Copyright (c) 2003-2020, Adrian Rossiter

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
   Name: minmax.cc
   Description: make equal-edge, face-regular and unscrambled polyhedra
   Project: Antiprism - http://www.antiprism.com
*/

#include <cmath>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>

#include "../base/antiprism.h"

using std::set;
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
      : num_iters(1000), num_iters_status(1000), sig_digits(16),
        rep_file(stderr)
  {
  }

  double get_test_val() { return pow(10, -sig_digits); }
  bool quiet() { return rep_file == nullptr; }
  bool check_status(int n)
  {
    return n == 1 || (num_iters_status > 0 && n % num_iters_status == 0);
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

class mm_opts : public ProgramOpts {
public:
  iter_params it_params;
  char algm;
  char placement;
  double shorten_by;
  double lengthen_by;
  double shorten_rad_by;
  double flatten_by;
  string sym_str;
  Vec4d ellipsoid;

  string ifile;
  string ofile;

  mm_opts()
      : ProgramOpts("minmax"), algm('\0'), placement('\0'), shorten_by(1.0),
        lengthen_by(NAN), shorten_rad_by(NAN), flatten_by(NAN)
  {
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

// clang-format off
void mm_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Read a file in OFF format containing a graph of a polyhedron, with or\n"
"without vertex coordinates, and try to create a spherical or ellipsoidal\n"
"tesselation where the maximum edge is a minimum length, or try to make\n"
"into a regular-faced polyhedron. Option adjustment factors expressed as\n"
"a percentage are approximate. If input_file is not given the program reads\n"
"from standard input.\n"
"\n"
"Options\n"
"%s"
"  -n <itrs> number of iterations (default 1000)\n"
"  -s <perc> percentage to shorten longest edges on iteration (default: 1)\n"
"  -l <perc> percentage to lengthen shortest edges (-a v) on iteration \n"
"            (default: 0.0)\n"
"  -k <perc> -a u - percentage to reduce polygon radius on iteration\n"
"                   (default: value of -s)\n"
"            -a e - percentage to reduce minimum target length on iteration\n"
"                   (default: 1e-6), will oscillate at certain precision,\n"
"                   process output with smaller value to improve solution\n"
"  -f <perc> percentage to reduce distance of vertex from face plane (-a u)\n"
"            on iteration (default: value of -s)\n"
"  -a <alg>  length changing algorithm\n"
"              e - equalise edges on sphere or ellipsoid (default)\n"
"              v - equalise shortest and longest edges attached to a vertex\n"
"                  on sphere or ellipsoid (unscramble algorithm, see -p u)\n"
"              u - make faces into unit-edged regular polygons\n"
"  -p <mthd> method of placement onto a unit sphere:\n"
"              n - project onto the sphere (default for -a v/e)\n"
"              r - random placement\n"
"              u - unscramble: place a small polygon on one side and the\n"
"                  rest of the vertices at a point on the other side\n"
"                  (sets -a v if not specifically set)\n"
"  -y <sub>  maintain symmetry of the base model: sub is symmetry\n"
"            subgroup (Schoenflies notation) or 'full', optionally followed\n"
"            by a ',' and conjugation type (integer)\n"
"  -E <prms> use ellipsoid, three numbers separated by commas are the\n"
"            axis lengths (for a superellipsoid an optional fourth number\n"
"            gives the power)\n"
"  -L <lim>  minimum change of distance/width_of_model to terminate, as \n"
"               negative exponent (default: %d giving %.0e)\n"
"  -z <n>    status checking and reporting every n iterations, -1 for no\n"
"            status (default: 1000)\n"
"  -q        quiet, do not print status messages\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text,
   it_params.sig_digits, it_params.get_test_val() );
}
// clang-format on

void mm_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;
  vector<double> nums;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hn:s:l:k:f:a:p:y:E:L:z:qo:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'o':
      ofile = optarg;
      break;

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
      print_status_or_exit(read_double(optarg, &shorten_by), c);
      if (shorten_by < 0 || shorten_by > 100)
        warning("not in range 0 to 100", c);
      break;

    case 'l':
      print_status_or_exit(read_double(optarg, &lengthen_by), c);
      if (lengthen_by < 0 || lengthen_by > 100) {
        warning("not in range 0 to 100", c);
      }
      break;

    case 'k':
      print_status_or_exit(read_double(optarg, &shorten_rad_by), c);
      if (shorten_rad_by < 0 || shorten_rad_by > 100) {
        warning("not in range 0 to 100", c);
      }
      break;

    case 'f':
      print_status_or_exit(read_double(optarg, &flatten_by), c);
      if (flatten_by < 0 || flatten_by > 100) {
        warning("not in range 0 to 100", c);
      }
      break;

    case 'a':
      if (strlen(optarg) > 1 || !strchr("evu", *optarg))
        error("method is '" + string(optarg) + "' must be e, v or u");
      algm = *optarg;
      break;

    case 'p':
      if (strlen(optarg) > 1 || !strchr("nur", *optarg))
        error("method is '" + string(optarg) + "' must be n, u, or r");
      placement = *optarg;
      break;

    case 'y':
      sym_str = optarg;
      break;

    case 'E':
      print_status_or_exit(read_double_list(optarg, nums), c);
      if (nums.size() < 3 || nums.size() > 4)
        error("must give three or four numbers only", c);
      ellipsoid =
          Vec4d(nums[0], nums[1], nums[2], (nums.size() == 4) ? nums[3] : 2);
      if (ellipsoid[3] <= 0)
        error("superellipsoid power must be greater than zero", c);
      break;

    case 'L':
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

    default:
      error("unknown command line error");
    }
  }

  // Default algorithm is 'e', but if unscrambling it is 'u'
  if (!algm)
    algm = (placement == 'u') ? 'v' : 'e';

  if (placement == 'u' && algm != 'v')
    error(msg_str("algorithm %c is not compatible with unscramble (use -a v)",
                  algm),
          'p');

  if (algm == 'u') {
    if (std::isnan(shorten_rad_by))
      shorten_rad_by = shorten_by;
    if (std::isnan(flatten_by))
      flatten_by = shorten_by;
    if (!std::isnan(lengthen_by))
      warning("set, but not used for this algorithm", 'l');
    if (ellipsoid.is_set())
      warning("set, but not used for this algorithm", 'E');
    if (placement)
      warning(msg_str("placement method used with -a u", placement), 'p');
  }
  else if (algm == 'v') {
    if (std::isnan(lengthen_by))
      lengthen_by = 0.0;
    if (!std::isnan(shorten_rad_by))
      warning("set, but not used for this algorithm", 'k');
    if (!std::isnan(flatten_by))
      warning("set, but not used for this algorithm", 'f');
    if (!placement)
      placement = 'n';
  }
  else if (algm == 'e') {
    if (std::isnan(shorten_rad_by))
      shorten_rad_by = 1e-6;
    if (!std::isnan(lengthen_by))
      warning("set, but not used for this algorithm", 'l');
    if (!std::isnan(flatten_by))
      warning("set, but not used for this algorithm", 'f');
    if (!placement)
      placement = 'n';
  }
  else
    error("invalid algorithm '%c'", algm);

  if (argc - optind > 1)
    error("too many arguments");

  if (argc - optind == 1)
    ifile = argv[optind];
}

Status init_sym(const Geometry &geom, const char *sym_str, Symmetry &sym)
{
  Status stat;
  char sym_cpy[MSG_SZ]; // big enough for normal use
  strcpy_msg(sym_cpy, sym_str);

  Symmetry full_sym(geom);
  vector<char *> parts;
  split_line(sym_cpy, parts, ",");
  if (parts.size() == 0 || parts.size() > 2)
    return Status::error("argument should have 1 or 2 comma separated parts");

  Symmetry sub_sym;
  if (strncmp(parts[0], "full", strlen(parts[0])) == 0)
    sub_sym = full_sym;
  else if (!(stat = sub_sym.init(parts[0], Trans3d())))
    return Status::error(msg_str("sub-symmetry type: %s", stat.c_msg()));

  int sub_sym_conj = 0;
  if (parts.size() > 1 && !(stat = read_int(parts[1], &sub_sym_conj)))
    return Status::error(
        msg_str("sub-symmetry conjugation number: %s", stat.c_msg()));

  if (!(stat = full_sym.get_sub_sym(sub_sym, &sym, sub_sym_conj)))
    return Status::error(msg_str("sub-symmetry: %s", stat.c_msg()));

  return Status::ok();
}

void to_ellipsoid(Vec3d &v, Vec4d ellipsoid)
{
  if (compare(v, Vec3d::zero) == 0)
    v = {0, 0, 1};

  if (ellipsoid.is_set()) {
    double scale = 0;
    for (unsigned int i = 0; i < 3; i++)
      scale += pow(fabs(v[i] / ellipsoid[i]), ellipsoid[3]);
    v *= pow(scale, -1 / ellipsoid[3]);
  }
  else
    v.to_unit();
}

void initial_placement(Geometry &geom, char placement, Vec4d ellipsoid)
{
  switch (placement) {
  case 'n':
    for (unsigned int i = 0; i < geom.verts().size(); i++)
      to_ellipsoid(geom.verts(i), ellipsoid);
    break;

  case 'u':
    for (unsigned int i = 0; i < geom.verts().size(); i++)
      geom.verts(i) = Vec3d(-1, 0, 0);
    if (geom.faces().size()) {
      const vector<int> &face = geom.faces(0);
      for (unsigned int i = 0; i < face.size(); i++) {
        geom.verts(face[i]) = Vec3d(1, 0.01 * cos(i * 2 * M_PI / face.size()),
                                    0.01 * sin(i * 2 * M_PI / face.size()));
        geom.verts(face[i]).to_unit();
      }
    }
    break;

  case 'r': {
    Random rnd;
    rnd.time_seed();
    for (unsigned int i = 0; i < geom.verts().size(); i++) {
      geom.verts(i) = Vec3d::random(rnd);
      to_ellipsoid(geom.verts(i), ellipsoid);
    }
    break;
  }
  }
}

vector<vector<int>> get_edge_connections(const Geometry &geom)
{
  vector<vector<int>> edges(geom.verts().size());
  for (unsigned int i = 0; i < geom.edges().size(); i++) {
    edges[geom.edges(i, 0)].push_back(geom.edges(i, 1));
    edges[geom.edges(i, 1)].push_back(geom.edges(i, 0));
  }
  return edges;
}

//------------------------------------------------------------------
// Unscramble

Status minmax_v(Geometry &geom, iter_params it_params, double shorten_factor,
                double lengthen_factor, Vec4d ellipsoid)
{
  auto eds = GeometryInfo(geom).get_vert_cons();

  int max_edge = 0, min_edge = 0, p0, p1;
  double dist, max_dist = 0, min_dist = 1e100, g_max_dist = 0,
               g_min_dist = 1e100;
  for (int cnt = 1; cnt <= it_params.num_iters; cnt++) {
    g_max_dist = 0;
    g_min_dist = 1e100;
    for (unsigned int v = 0; v < geom.verts().size(); v++) {
      if (eds[v].size() == 0)
        continue;
      max_dist = -1;
      min_dist = 1e100;
      for (unsigned int i = 0; i < eds[v].size(); i++) {
        dist = (geom.verts(v) - geom.verts(eds[v][i])).len2();
        // fprintf(stderr, "dist=%g\n", dist);
        if (dist > max_dist) {
          max_dist = dist;
          max_edge = i;
        }
        if (dist > g_max_dist)
          g_max_dist = dist;
        if (dist < min_dist) {
          min_dist = dist;
          min_edge = i;
        }
        if (dist < g_min_dist)
          g_min_dist = dist;
      }

      p0 = v;
      p1 = eds[v][max_edge];
      Vec3d diff = geom.verts(p1) - geom.verts(p0);
      geom.verts(p0) += diff * shorten_factor;
      to_ellipsoid(geom.verts(p0), ellipsoid);

      p0 = v;
      p1 = eds[v][min_edge];
      diff = geom.verts(p1) - geom.verts(p0);
      geom.verts(p0) -= diff * lengthen_factor;
      to_ellipsoid(geom.verts(p0), ellipsoid);
    }

    if (!it_params.quiet() && it_params.check_status(cnt))
      fprintf(it_params.rep_file, "\niter:%-15d max:%17.15f min:%17.15f ", cnt,
              sqrt(g_max_dist), sqrt(g_min_dist));
    else if (it_params.print_progress_dot(cnt))
      fprintf(it_params.rep_file, ".");
  }
  if (!it_params.quiet() && it_params.checking_status())
    fprintf(it_params.rep_file,
            "\nFinal:\niter:%-15d max:%17.15f min:%17.15f\n",
            it_params.num_iters, sqrt(g_max_dist), sqrt(g_min_dist));

  return Status::ok();
}

//------------------------------------------------------------------
// Equal edge

Status minmax_e_orig(Geometry &geom, iter_params it_params,
                     double shorten_factor, double shrink_factor,
                     Vec4d ellipsoid)
{
  double g_max_dist = 0;
  double g_min_dist = 1e100;
  double g_scale_factor = 1;
  double max_dist = 0;
  double min_dist = 1e100;

  auto eds = GeometryInfo(geom).get_vert_cons();

  for (int cnt = 0; cnt <= it_params.num_iters; cnt++) {
    // Vertex offsets for the iteration.
    const vector<Vec3d> &verts = geom.verts();
    vector<Vec3d> offsets(verts.size(), Vec3d::zero);

    max_dist = 0;
    min_dist = 1e100;
    for (unsigned int v = 0; v < verts.size(); v++) {
      if (eds[v].size() == 0)
        continue;
      for (unsigned int i = 0; i < eds[v].size(); i++) {
        auto vec = geom.verts(v) - geom.verts(eds[v][i]);
        double dist = vec.len();
        if (dist > max_dist)
          max_dist = dist;
        if (dist < min_dist)
          min_dist = dist;
        offsets[v] += (vec / dist) * (g_min_dist - dist) * g_scale_factor *
                      shorten_factor;
      }
    }

    // adjust vertices post-loop (skip first time as global values not set)
    if (cnt > 0) {
      for (unsigned int i = 0; i < offsets.size(); i++) {
        geom.raw_verts()[i] += offsets[i];
        to_ellipsoid(geom.verts(i), ellipsoid);
      }
    }

    g_min_dist = min_dist * (1 - shrink_factor);
    g_max_dist = max_dist;
    g_scale_factor = (g_max_dist - g_min_dist) / g_min_dist;

    if (cnt > 0) {
      if (!it_params.quiet() && it_params.check_status(cnt))
        fprintf(it_params.rep_file, "\niter:%-15d max:%17.15f min:%17.15f ",
                cnt, max_dist, min_dist);
      else if (it_params.print_progress_dot(cnt))
        fprintf(it_params.rep_file, ".");
    }
  }
  if (!it_params.quiet() && it_params.checking_status())
    fprintf(it_params.rep_file,
            "\nFinal:\niter:%-15d max:%17.15f min:%17.15f\n",
            it_params.num_iters, max_dist, min_dist);

  return Status::ok();
}

Status minmax_e_sym(Geometry &base_geom, iter_params it_params,
                    double shorten_factor, double shrink_factor,
                    Vec4d ellipsoid, string sym_str)
{
  Symmetry sym(Symmetry::C1);
  Status stat;
  if (!(stat = init_sym(base_geom, sym_str.c_str(), sym)))
    return stat;

  SymmetricUpdater sym_updater(base_geom, sym, false);
  const Geometry &geom = sym_updater.get_geom_reading();
  auto eds = GeometryInfo(geom).get_vert_cons();

  // Get a list of the vertices that are principal vertices, or are
  // connected to a principal vertex
  vector<int> principal_verts;
  set<int> verts_to_process;
  for (auto &v_orbit : sym_updater.get_equiv_sets(VERTS)) {
    int principal_idx = *v_orbit.begin();
    principal_verts.push_back(principal_idx);
    verts_to_process.insert(principal_idx);
    for (int v_idx : eds[principal_idx])
      verts_to_process.insert(v_idx);
  }

  double g_max_dist = 0;
  double g_min_dist = 1e100;
  double g_scale_factor = 1;
  double max_dist = 0;
  double min_dist = 1e100;

  for (int cnt = 0; cnt <= it_params.num_iters; cnt++) {
    // Vertex offsets for the iteration.
    const vector<Vec3d> &verts = geom.verts();
    vector<Vec3d> offsets(verts.size(), Vec3d::zero);

    // Ensure that the vertices used in adjustment are up to date
    for (int v_idx : verts_to_process)
      sym_updater.update_from_principal_vertex(v_idx);

    max_dist = 0;
    min_dist = 1e100;
    for (int v_idx : principal_verts) {
      if (eds[v_idx].size() == 0)
        continue;
      for (unsigned int i = 0; i < eds[v_idx].size(); i++) {
        auto vec = geom.verts(v_idx) - geom.verts(eds[v_idx][i]);
        double dist = vec.len();
        if (dist > max_dist)
          max_dist = dist;
        if (dist < min_dist)
          min_dist = dist;
        offsets[v_idx] += (vec / dist) * (g_min_dist - dist) * g_scale_factor *
                          shorten_factor;
      }
    }

    // adjust vertices post-loop (skip first time as global values not set)
    if (cnt > 0) {
      for (int v_idx : principal_verts) {
        auto vert = verts[v_idx] + offsets[v_idx];
        to_ellipsoid(vert, ellipsoid);
        sym_updater.update_principal_vertex(v_idx, vert);
      }
    }

    g_min_dist = min_dist * (1 - shrink_factor);
    g_max_dist = max_dist;
    g_scale_factor = (g_max_dist - g_min_dist) / g_min_dist;

    if (cnt > 0) {
      if (!it_params.quiet() && it_params.check_status(cnt))
        fprintf(it_params.rep_file, "\niter:%-15d max:%17.15f min:%17.15f ",
                cnt, max_dist, min_dist);
      else if (it_params.print_progress_dot(cnt))
        fprintf(it_params.rep_file, ".");
    }
  }
  if (!it_params.quiet() && it_params.checking_status())
    fprintf(it_params.rep_file,
            "\nFinal:\niter:%-15d max:%17.15f min:%17.15f\n",
            it_params.num_iters, max_dist, min_dist);

  base_geom = sym_updater.get_geom_final();
  return Status::ok();
}

Status minmax_e(Geometry &base_geom, iter_params it_params,
                double shorten_factor, double shrink_factor,
                Vec4d ellipsoid = Vec4d(), string sym_str = string())
{
  // If no symmetry argument then process with the original algorithm
  if (sym_str.empty())
    return minmax_e_orig(base_geom, it_params, shorten_factor, shrink_factor,
                         ellipsoid);
  else
    return minmax_e_sym(base_geom, it_params, shorten_factor, shrink_factor,
                        ellipsoid, sym_str);
}

//------------------------------------------------------------------
// Unit edge regular polygon

Status minmax_unit_orig(Geometry &geom, iter_params it_params,
                        double shorten_factor, double plane_factor,
                        double radius_factor)
{
  double test_val = it_params.get_test_val();
  const double divergence_test2 = 1e30; // test vertex dist^2 for divergence
  // Scale to get edges close to 1
  GeometryInfo info(geom);
  double scale = info.iedge_length_lims().sum / info.num_iedges();
  if (scale)
    geom.transform(Trans3d::scale(1 / scale));

  const vector<Vec3d> &verts = geom.verts();
  const vector<vector<int>> &faces = geom.faces();

  double max_diff2 = 0;
  Vec3d origin(0, 0, 0);
  vector<double> rads(faces.size());
  for (unsigned int f = 0; f < faces.size(); f++) {
    int N = faces[f].size();
    int D = abs(find_polygon_denominator_signed(geom, f, epsilon));
    if (!D)
      D = 1;
    rads[f] = 0.5 / sin(M_PI * D / N); // circumradius of regular polygon
  }

  bool diverging = false;
  int cnt = 0;
  for (cnt = 1; cnt <= it_params.num_iters; cnt++) {
    // Vertex offsets for the iteration.
    vector<Vec3d> offsets(verts.size(), Vec3d::zero);
    for (unsigned int f_idx = 0; f_idx < faces.size(); f_idx++) {
      const vector<int> &face = faces[f_idx];
      const unsigned int f_sz = face.size();
      Vec3d norm = geom.face_norm(f_idx).unit();
      Vec3d f_cent = geom.face_cent(f_idx);
      if (vdot(norm, f_cent) < 0)
        norm *= -1.0;

      for (unsigned int i = 0; i < f_sz; i++) {
        int v_idx = face[i];
        int v_idx_next = face[(i + 1) % f_sz];

        // offset for unit edges
        vector<int> edge = {v_idx, v_idx_next};
        Vec3d offset =
            (1 - geom.edge_len(edge)) * shorten_factor * geom.edge_vec(edge);
        offsets[v_idx] -= offset;
        offsets[v_idx_next] += offset;

        // offset for planarity
        offsets[v_idx] +=
            vdot(plane_factor * norm, f_cent - verts[v_idx]) * norm;

        // offset for polygon radius
        Vec3d rad_vec = (verts[v_idx] - f_cent);
        offsets[v_idx] +=
            (rads[f_idx] - rad_vec.len()) * radius_factor * rad_vec;
      }
    }

    // adjust vertices post-loop
    for (unsigned int i = 0; i < offsets.size(); i++)
      geom.raw_verts()[i] += offsets[i];

    if (it_params.check_status(cnt)) {
      max_diff2 = 0;
      for (auto &offset : offsets) {
        double diff2 = offset.len2();
        if (diff2 > max_diff2)
          max_diff2 = diff2;
      }

      double width = BoundBox(verts).max_width();
      if (sqrt(max_diff2) / width < test_val)
        break;

      if (!it_params.quiet())
        fprintf(it_params.rep_file, "iter:%-15d max_diff:%17.15e\n", cnt,
                sqrt(max_diff2));

      // see if radius is expanding or contracting unreasonably
      if (divergence_test2 > 0) {
        diverging = false;
        for (auto &offset : offsets)
          if (offset.len2() > divergence_test2) {
            diverging = true;
            break;
          }
      }
    }
  }

  if (!it_params.quiet() && diverging)
    fprintf(it_params.rep_file, "Probably Diverging. Breaking out.\n");
  if (!it_params.quiet() && it_params.checking_status())
    fprintf(it_params.rep_file, "Final:\niter:%-15d max_diff:%17.15e\n",
            cnt - 1 + diverging, sqrt(max_diff2));

  return Status::ok();
}

Status minmax_unit_sym(Geometry &base_geom, iter_params it_params,
                       double shorten_factor, double plane_factor,
                       double radius_factor, const string &sym_str)
{
  Symmetry sym(Symmetry::C1);
  Status stat;
  if (!(stat = init_sym(base_geom, sym_str.c_str(), sym)))
    return stat;

  double test_val = it_params.get_test_val();
  const double divergence_test2 = 1e30; // test vertex dist^2 for divergence
  // Scale to get edges close to 1
  GeometryInfo info(base_geom);
  double scale = info.iedge_length_lims().sum / info.num_iedges();
  if (scale)
    base_geom.transform(Trans3d::scale(1 / scale));
  info.get_vert_cons();
  SymmetricUpdater sym_updater(base_geom, sym, false);
  const Geometry &geom = sym_updater.get_geom_reading();
  const vector<Vec3d> &verts = geom.verts();
  const vector<vector<int>> &faces = geom.faces();

  double max_diff2 = 0;
  Vec3d origin(0, 0, 0);
  vector<double> rads(faces.size());
  for (unsigned int f = 0; f < faces.size(); f++) {
    int N = faces[f].size();
    int D = abs(find_polygon_denominator_signed(geom, f, epsilon));
    if (!D)
      D = 1;
    rads[f] = 0.5 / sin(M_PI * D / N); // circumradius of regular polygon
  }

  // Get a list of the faces that contain a principal vertex of any type
  vector<int> principal_verts;
  vector<int> verts_to_update;
  vector<int> faces_to_process;
  {
    vector<set<int>> vert_faces(geom.verts().size());
    for (int f_idx = 0; f_idx < (int)geom.faces().size(); f_idx++)
      for (int v_idx : faces[f_idx])
        vert_faces[v_idx].insert(f_idx);

    for (auto &v_orbit : sym_updater.get_equiv_sets(VERTS)) {
      int principal_idx = *v_orbit.begin();
      principal_verts.push_back(principal_idx);
      for (int f_idx : vert_faces[principal_idx]) {
        faces_to_process.push_back(f_idx);
        for (int v_idx : faces[f_idx])
          verts_to_update.push_back(v_idx);
      }
    }
  }
  sort(verts_to_update.begin(), verts_to_update.end());
  verts_to_update.erase(unique(verts_to_update.begin(), verts_to_update.end()),
                        verts_to_update.end());
  sort(faces_to_process.begin(), faces_to_process.end());
  faces_to_process.erase(
      unique(faces_to_process.begin(), faces_to_process.end()),
      faces_to_process.end());

  bool diverging = false;
  int cnt = 0;
  for (cnt = 1; cnt <= it_params.num_iters; cnt++) {
    // Vertx offsets for the iteration.
    vector<Vec3d> offsets(verts.size(), Vec3d::zero);

    // Ensure that the vertices used in adjustment are up to date
    for (int v_idx : verts_to_update)
      sym_updater.update_from_principal_vertex(v_idx);

    for (int f_idx : faces_to_process) {
      const vector<int> &face = faces[f_idx];

      Vec3d norm = geom.face_norm(f_idx).unit();
      Vec3d f_cent = geom.face_cent(f_idx);
      if (vdot(norm, f_cent) < 0)
        norm *= -1.0;

      for (int i = 0; i < (int)face.size(); i++) {
        const int v_idx = face[i];
        // offset for planarity
        offsets[v_idx] +=
            vdot(plane_factor * norm, f_cent - verts[v_idx]) * norm;

        // offset for polygon radius
        Vec3d rad_vec = (verts[v_idx] - f_cent);
        offsets[v_idx] +=
            (rads[f_idx] - rad_vec.len()) * radius_factor * rad_vec;
      }
    }

    // offsets for unit edges
    for (int v_idx : principal_verts) {
      for (int v_neigh : info.get_vert_cons()[v_idx]) {
        vector<int> edge = {v_idx, v_neigh};
        Vec3d offset =
            (1 - geom.edge_len(edge)) * shorten_factor * geom.edge_vec(edge);
        offsets[v_idx] -= 2 * offset;
      }
    }

    // adjust principal vertices post-loop
    for (int v_idx : principal_verts)
      sym_updater.update_principal_vertex(v_idx, verts[v_idx] + offsets[v_idx]);

    if (it_params.check_status(cnt)) {
      max_diff2 = 0;
      for (auto &offset : offsets) {
        double diff2 = offset.len2();
        if (diff2 > max_diff2)
          max_diff2 = diff2;
      }

      double width = BoundBox(verts).max_width();
      if (sqrt(max_diff2) / width < test_val)
        break;

      if (!it_params.quiet())
        fprintf(it_params.rep_file, "iter:%-15d max_diff:%17.15e\n", cnt,
                sqrt(max_diff2));

      // see if radius is expanding or contracting unreasonably
      if (divergence_test2 > 0) {
        diverging = false;
        for (auto &offset : offsets)
          if (offset.len2() > divergence_test2) {
            diverging = true;
            break;
          }
      }
    }
  }

  if (!it_params.quiet() && diverging)
    fprintf(it_params.rep_file, "Probably Diverging. Breaking out.\n");
  if (!it_params.quiet() && it_params.checking_status())
    fprintf(it_params.rep_file, "Final:\niter:%-15d max_diff:%17.15e\n",
            cnt - 1 + diverging, sqrt(max_diff2));

  base_geom = sym_updater.get_geom_final();
  return Status::ok();
}

Status minmax_unit(Geometry &base_geom, iter_params it_params,
                   double shorten_factor, double plane_factor,
                   double radius_factor, const string &sym_str)
{
  // If no symmetry argument then process with the original algorithm
  if (sym_str.empty())
    return minmax_unit_orig(base_geom, it_params, shorten_factor, plane_factor,
                            radius_factor);
  else
    return minmax_unit_sym(base_geom, it_params, shorten_factor, plane_factor,
                           radius_factor, sym_str);
}

int main(int argc, char *argv[])
{
  mm_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);

  if (opts.placement)
    initial_placement(geom, opts.placement, opts.ellipsoid);

  if (geom.faces().size() == 0)
    opts.warning("no faces so algorithm not applied");
  else if (opts.algm == 'e')
    minmax_e(geom, opts.it_params, opts.shorten_by / 200,
             opts.shorten_rad_by / 200, opts.ellipsoid, opts.sym_str);
  else if (opts.algm == 'v')
    minmax_v(geom, opts.it_params, opts.shorten_by / 200,
             opts.lengthen_by / 200, opts.ellipsoid);
  else if (opts.algm == 'u')
    opts.print_status_or_exit(minmax_unit(
        geom, opts.it_params, opts.shorten_by / 200, opts.flatten_by / 200,
        opts.shorten_rad_by / 200, opts.sym_str));

  opts.write_or_error(geom, opts.ofile);

  return 0;
}
