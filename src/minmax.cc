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
   Name: minmax.cc
   Description: minimise the maximum edge length on a sphere
   Project: Antiprism - http://www.antiprism.com
*/

#include <cmath>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>

#include "../base/antiprism.h"

using std::isnan;
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

class mm_opts : public ProgramOpts {
public:
  iter_params it_params;
  char algm;
  char placement;
  double shorten_by;
  double lengthen_by;
  double shorten_rad_by;
  double flatten_by;
  Vec4d ellipsoid;

  string ifile;
  string ofile;

  mm_opts()
      : ProgramOpts("minmax"), algm('v'), placement('n'), shorten_by(1.0),
        lengthen_by(NAN), shorten_rad_by(NAN), flatten_by(NAN)
  {
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

void mm_opts::usage()
{
  fprintf(
      stdout,
      "\n"
      "Usage: %s [options] [input_file]\n"
      "\n"
      "Read a file in OFF format containing a graph of a polyhedron, with or\n"
      "without vertex coordinates, and try to create a spherical or "
      "ellipsoidal\n"
      "tesselation where the maximum edge is a minimum length, or try to make\n"
      "into a regular polyhedron. If input_file is not given the program "
      "reads\n"
      "from standard input.\n"
      "\n"
      "Options\n"
      "%s"
      "  -n <itrs> number of iterations (default 1000)\n"
      "  -s <perc> percentage to shorten longest edges on iteration (default: "
      "1)\n"
      "  -l <perc> percentage to lengthen shortest edges (-a a/v) on iteration "
      "\n"
      "            (default: 0.0)\n"
      "  -k <perc> percentage to reduce polygon radius (-a u) on iteration\n"
      "            (default: value of -s)\n"
      "  -f <perc> percentage to reduce distance of vertex from face plane (-a "
      "u)\n"
      "            on iteration (default: value of -s)\n"
      "  -a <alg>  length changing algorithm\n"
      "              v - shortest and longest edges attached to a vertex "
      "(default)\n"
      "              a - shortest and longest of all edges\n"
      "              u - make faces into unit-edged regular polygons (-l "
      "controls\n"
      "                  planarity, ignore -p, -E)\n"
      "  -p <mthd> method of placement onto a unit sphere:\n"
      "              n - project onto the sphere (default)\n"
      "              r - random placement\n"
      "              u - unscramble: place a small polygon on one side and "
      "the\n"
      "                  rest of the vertices at a point on the other side\n"
      "  -E <prms> use ellipsoid, three numbers separated by commas are the\n"
      "            axis lengths (for a superellipsoid an optional fourth "
      "number\n"
      "            gives the power)\n"
      "  -L <lim>  minimum change of distance/width_of_model to terminate, as "
      "\n"
      "               negative exponent (default: %d giving %.0e)\n"
      "  -z <n>    status checking and reporting every n iterations, -1 for "
      "no\n"
      "            status (default: 1000)\n"
      "  -q        quiet, do not print status messages\n"
      "  -o <file> write output to file (default: write to standard output)\n"
      "\n"
      "\n",
      prog_name(), help_ver_text, it_params.sig_digits,
      it_params.get_test_val());
}

void mm_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;
  vector<double> nums;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hn:s:l:k:f:a:p:E:L:z:qo:")) != -1) {
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
      if (strlen(optarg) > 1 || !strchr("avu", *optarg))
        error("method is '" + string(optarg) + "' must be a, v or u");
      algm = *optarg;
      break;

    case 'p':
      if (strlen(optarg) > 1 || !strchr("nur", *optarg))
        error("method is '" + string(optarg) + "' must be n, u, or r");
      placement = *optarg;
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

  if (algm == 'u') {
    if (isnan(shorten_rad_by))
      shorten_rad_by = shorten_by;
    if (isnan(flatten_by))
      flatten_by = shorten_by;
    if (!isnan(lengthen_by))
      warning("set, but not used for this algorithm", 'l');
  }
  else { // algm ia v or a
    if (isnan(shorten_rad_by))
      lengthen_by = 0.0;
    if (!isnan(shorten_rad_by))
      warning("set, but not used for this algorithm", 'k');
    if (!isnan(flatten_by))
      warning("set, but not used for this algorithm", 'f');
  }

  if (argc - optind > 1)
    error("too many arguments");

  if (argc - optind == 1)
    ifile = argv[optind];
}

void to_ellipsoid(Vec3d &v, Vec4d ellipsoid)
{
  if (!ellipsoid.is_set())
    v.to_unit();
  else if (v.len2() > epsilon) {
    double scale = 0;
    for (unsigned int i = 0; i < 3; i++)
      scale += pow(fabs(v[i] / ellipsoid[i]), ellipsoid[3]);
    v *= pow(scale, -1 / ellipsoid[3]);
  }
}

void initial_placement(Geometry &geom, char placement, Vec4d ellipsoid)
{
  const vector<int> &face = geom.faces(0);
  switch (placement) {
  case 'n':
    for (unsigned int i = 0; i < geom.verts().size(); i++)
      to_ellipsoid(geom.verts(i), ellipsoid);
    break;

  case 'u':
    for (unsigned int i = 0; i < geom.verts().size(); i++)
      geom.verts(i) = Vec3d(-1, 0, 0);
    for (unsigned int i = 0; i < face.size(); i++) {
      geom.verts(face[i]) = Vec3d(1, 0.01 * cos(i * 2 * M_PI / face.size()),
                                  0.01 * sin(i * 2 * M_PI / face.size()));
      geom.verts(face[i]).to_unit();
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

void minmax_a(Geometry &geom, iter_params it_params, double shorten_factor,
              double lengthen_factor, Vec4d ellipsoid = Vec4d())
{
  int max_edge = 0, min_edge = 0, p0, p1;
  double dist, max_dist = 0, min_dist = 1e100;
  for (int cnt = 1; cnt <= it_params.num_iters; cnt++) {
    max_dist = 0;
    min_dist = 1e100;
    for (unsigned int i = 0; i < geom.edges().size(); i++) {
      dist = geom.edge_vec(i).len2();
      if (dist > max_dist) {
        max_dist = dist;
        max_edge = i;
      }
      if (dist < min_dist) {
        min_dist = dist;
        min_edge = i;
      }
    }

    p0 = geom.edges(max_edge, 0);
    p1 = geom.edges(max_edge, 1);
    Vec3d diff = geom.edge_vec(max_edge);
    geom.verts(p0) += diff * shorten_factor;
    to_ellipsoid(geom.verts(p0), ellipsoid);
    geom.verts(p1) -= diff * shorten_factor;
    to_ellipsoid(geom.verts(p1), ellipsoid);

    p0 = geom.edges(min_edge, 0);
    p1 = geom.edges(min_edge, 1);
    diff = geom.edge_vec(min_edge);
    geom.verts(p0) -= diff * lengthen_factor;
    to_ellipsoid(geom.verts(p0), ellipsoid);
    geom.verts(p1) += diff * lengthen_factor;
    to_ellipsoid(geom.verts(p1), ellipsoid);

    if (!it_params.quiet() && it_params.check_status(cnt))
      fprintf(it_params.rep_file, "\niter:%-15d max:%17.15f min:%17.15f ", cnt,
              max_dist, min_dist);
    else if (it_params.print_progress_dot(cnt))
      fprintf(it_params.rep_file, ".");
  }
  if (!it_params.quiet() && it_params.checking_status())
    fprintf(it_params.rep_file,
            "\nFinal:\niter:%-15d max:%17.15f min:%17.15f\n",
            it_params.num_iters, max_dist, min_dist);
}

void minmax_v(Geometry &geom, iter_params it_params, vector<vector<int>> &eds,
              double shorten_factor, double lengthen_factor,
              Vec4d ellipsoid = Vec4d())
{
  int max_edge = 0, min_edge = 0, p0, p1;
  double dist, max_dist = 0, min_dist = 1e100, g_max_dist = 0,
               g_min_dist = 1e100;
  for (int cnt = 1; cnt <= it_params.num_iters; cnt++) {
    g_max_dist = 0;
    g_min_dist = 1e100;
    for (unsigned int v = 0; v < geom.verts().size(); v++) {
      if (eds[v].size() == 0)
        continue;
      max_dist = 0;
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
              g_max_dist, g_min_dist);
    else if (it_params.print_progress_dot(cnt))
      fprintf(it_params.rep_file, ".");
  }
  if (!it_params.quiet() && it_params.checking_status())
    fprintf(it_params.rep_file,
            "\nFinal:\niter:%-15d max:%17.15f min:%17.15f\n",
            it_params.num_iters, g_max_dist, g_min_dist);
}

void minmax_unit(Geometry &geom, iter_params it_params, double shorten_factor,
                 double plane_factor, double radius_factor)
{
  double test_val = it_params.get_test_val();
  const double divergence_test2 = 1e30; // test vertex dist^2 for divergence
  // do a scale to get edges close to 1
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
    // fprintf(stderr, "{%d/%d} rad=%g\n", N, D, rads[f]);
  }

  bool diverging = false;
  int cnt = 0;
  for (cnt = 1; cnt <= it_params.num_iters; cnt++) {
    vector<Vec3d> old_verts = verts;

    // Vertx offsets for the iteration.
    vector<Vec3d> offsets(verts.size(), Vec3d::zero);
    for (unsigned int ff = cnt; ff < faces.size() + cnt; ff++) {
      const unsigned int f = ff % faces.size();
      const vector<int> &face = faces[f];
      const unsigned int f_sz = face.size();
      Vec3d norm = geom.face_norm(f).unit();
      Vec3d f_cent = geom.face_cent(f);
      if (vdot(norm, f_cent) < 0)
        norm *= -1.0;

      for (unsigned int vv = cnt; vv < f_sz + cnt; vv++) {
        unsigned int v = vv % f_sz;
        // offset for unit edges
        vector<int> edge = make_edge(face[v], face[(v + 1) % f_sz]);
        Vec3d offset =
            (1 - geom.edge_len(edge)) * shorten_factor * geom.edge_vec(edge);
        offsets[edge[0]] -= offset;
        offsets[edge[1]] += offset;

        // offset for planarity
        offsets[face[v]] +=
            vdot(plane_factor * norm, f_cent - verts[face[v]]) * norm;

        // offset for polygon radius
        Vec3d rad_vec = (verts[face[v]] - f_cent);
        offsets[face[v]] += (rads[f] - rad_vec.len()) * radius_factor * rad_vec;
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
    fprintf(it_params.rep_file, "\nfinal: iter:%-15d max_diff:%17.15e\n", cnt,
            sqrt(max_diff2));
}

int main(int argc, char *argv[])
{
  mm_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);

  if (!geom.edges().size())
    geom.add_missing_impl_edges();

  if (geom.edges().size()) {
    if (opts.algm != 'u')
      initial_placement(geom, opts.placement, opts.ellipsoid);
    if (opts.algm == 'a')
      minmax_a(geom, opts.it_params, opts.shorten_by / 200,
               opts.lengthen_by / 200, opts.ellipsoid);
    else {
      vector<vector<int>> eds(geom.verts().size());
      for (unsigned int i = 0; i < geom.edges().size(); i++) {
        eds[geom.edges(i, 0)].push_back(geom.edges(i, 1));
        eds[geom.edges(i, 1)].push_back(geom.edges(i, 0));
      }
      if (opts.algm == 'v')
        minmax_v(geom, opts.it_params, eds, opts.shorten_by / 200,
                 opts.lengthen_by / 200, opts.ellipsoid);
      else if (opts.algm == 'u')
        minmax_unit(geom, opts.it_params, opts.shorten_by / 200,
                    opts.flatten_by / 200, opts.shorten_rad_by / 200);
    }
  }
  else
    opts.warning("input file contains no edges");

  opts.write_or_error(geom, opts.ofile);

  return 0;
}
