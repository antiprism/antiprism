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
   Name: polar_recip.cc
   Description: make a polar reciprocal
   Project: Antiprism - http://www.antiprism.com
*/

#include <algorithm>
#include <ctype.h>
#include <map>
#include <math.h>
#include <string.h>
#include <string>
#include <vector>

#include "../base/antiprism.h"

using std::string;
using std::vector;
using std::map;
using std::swap;

using namespace anti;

class pr_opts : public ProgramOpts {
private:
public:
  vector<int> space_verts;
  double recip_rad;
  double init_rad;
  char recip_rad_type;
  Vec3d centre;
  Vec3d init_cent;
  char recip_cent_type;
  char init_cent_type;
  bool invert;
  double inf;
  double extra_ideal_elems;
  int num_iters;
  double epsilon;
  bool append;

  string ifile;
  string ofile;

  pr_opts()
      : ProgramOpts("pol_recip"), recip_rad(0), init_rad(0),
        recip_rad_type('x'), recip_cent_type('x'), init_cent_type('x'),
        invert(false), inf(1e15), extra_ideal_elems(true), num_iters(100),
        epsilon(0), append(false)
  {
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

// clang-format off
void pr_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Read a file in OFF format and make a polar reciprocal from the face\n"
"planes. If input_file is not given the program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"  -c <cent> reciprocation centre (default: C)\n"
"              X,Y,Z - centre with these coordinates\n"
"              C - vertex centroid,    M - mid-sphere (approximate)\n"
"              R - circumcentre (avg)\n"
"              e - edge balanced,      E - edge balanced with inversion\n"
"              v - vert/face balanced, V - vert/face balanced with inversion\n"
"              a - v/f/e balanced      A - v/f/e balanced with inversion\n"
"  -C <init> initial value for a centre calculation, in form 'X,Y,Z',\n"
"            or C to use centroid (default), M to calculate approx mid-sphere\n"
"  -r <rad>  reciprocation radius (default: calculated)\n"
"              radius value - use this value as the radius\n"
"              v - nearest vertex distance, V - furthest vertex distance\n"
"              e - nearest edge distance,   E - furthest edge distance\n"
"              f - nearest face distance,   F - furthest face distance\n"
"            or a comma separated list of vertex indices (starting from 0)\n"
"            and the distance is to the space containing those vertices\n"
"  -R <rad>  initial value for a radius calculation (default: calculated)\n"
"  -i        transform dual by reflecting in centre point\n"
"  -I <dist> maximum distance to project any normal or infinite dual vertex\n"
"            (default: %.0e), if 0 then use actual distances and delete\n"
"            infinite points\n"
"  -x        exclude extra elements added to duals with ideal vertices\n"
"  -n <itrs> maximum number of iterations (default: 10000)\n"
"  -l <lim>  minimum distance change to terminate, as negative exponent\n"
"               (default: %d giving %.0e)\n"
"  -a        append dual to original polyhedron\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text, inf, int(-log(::epsilon)/log(10) + 0.5), ::epsilon);
}
// clang-format on

void pr_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  int sig_compare = INT_MAX;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hc:C:r:R:xao:n:l:iI:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'C':
      if (strlen(optarg) == 1 && strchr("CM", *optarg))
        init_cent_type = *optarg;
      else if (init_cent.read(optarg))
        init_cent_type = 'c';
      else
        error("initial centre must be three coordinates, C, or M", c);
      break;

    case 'c':
      if (strlen(optarg) == 1 && strchr("CRMeEvVaA", *optarg))
        recip_cent_type = *optarg;
      else if (centre.read(optarg))
        recip_cent_type = 'c';
      else
        error("centre type must be three coordinates, or letter from CRMeEvVaA",
              c);
      break;

    case 'R':
      print_status_or_exit(read_double(optarg, &init_rad), c);
      if (fabs(init_rad) < epsilon)
        warning("radius is very small (the reciprocal will be very large)", c);
      else if (init_rad < 0)
        warning("radius is negative", c);
      break;

    case 'r':
      if (strlen(optarg) == 1 && strchr("VvEeFf", *optarg)) {
        recip_rad_type = *optarg;
      }
      else if (read_double(optarg, &recip_rad)) {
        if (fabs(recip_rad) < epsilon)
          warning("radius is very small (the reciprocal will be very large)",
                  c);
        else if (recip_rad < 0)
          warning("radius is negative", c);
        recip_rad_type = 'r';
      }
      else if (read_int_list(optarg, space_verts, true))
        recip_rad_type = 'S';
      else
        error("radius must be a radius value, or r,V,v,E,e,F,f or a list of "
              "index numbers",
              c);
      break;

    case 'i':
      invert = true;
      break;

    case 'x':
      extra_ideal_elems = false;
      break;

    case 'a':
      append = true;
      break;

    case 'o':
      ofile = optarg;
      break;

    case 'I':
      print_status_or_exit(read_double(optarg, &inf), c);
      if (inf < 0.0)
        error("distance must be positive, or 0 to disable", c);
      break;

    case 'n':
      print_status_or_exit(read_int(optarg, &num_iters), c);
      if (num_iters < 0)
        error("number of iterations must be greater than 0", c);
      break;

    case 'l':
      print_status_or_exit(read_int(optarg, &sig_compare), c);
      if (sig_compare < 0) {
        warning("limit is negative, and so ignored", c);
      }
      if (sig_compare > DEF_SIG_DGTS) {
        warning("limit is very small, may not be attainable", c);
      }
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

  epsilon = (sig_compare != INT_MAX) ? pow(10, -sig_compare) : ::epsilon;
}

int find_mid_centre(Geometry &geom, double &rad, Vec3d &cent, int n,
                    const double &eps)
{
  vector<vector<int>> g_edges;
  geom.get_impl_edges(g_edges);
  int e_sz = g_edges.size();
  if (!cent.is_set())
    cent = geom.centroid();
  if (fabs(rad) < epsilon) {
    GeometryInfo rep(geom);
    rep.set_center(cent);
    rad = rep.iedge_dist_lims().sum / e_sz;
  }
  Vec3d cur_cent = cent;
  double cur_rad = rad;
  double cent_test = -1, rad_test = -1;
  int cnt;
  for (cnt = 0; cnt < n; cnt++) {
    Vec3d c_diff(0, 0, 0);
    double rad_sum_g = 0;
    for (int e = 0; e < e_sz; ++e) {
      Vec3d Pg = nearest_point(cur_cent, geom.verts(), g_edges[e]);
      double r = (Pg - cur_cent).len();
      rad_sum_g += r;
      c_diff += (1 - r / cur_rad) * (Pg - cur_cent);
    }
    c_diff = c_diff / (double)e_sz;
    cur_cent = cent;
    cent = cur_cent - c_diff / (double)e_sz;
    cur_rad = rad;
    rad = rad_sum_g / e_sz;
    cent_test = (cur_cent - cent).len() / rad;
    rad_test = (cur_rad - rad) / rad;
    // cent.dump("cent");
    // fprintf(stderr, "rad=%g, delt_rad=%g,delta_cent=%g\n", rad, cur_rad-rad,
    //      (cur_cent-cent).len());
    if (fabs(cent_test) < eps && fabs(rad_test) < eps)
      break;
  }
  char str[MSG_SZ];
  fprintf(stderr, "[n=%d, limit=%sachieved, r_test=%g, c_test=%g]\n"
                  "centre=(%s), radius=%.16g\n",
          cnt, cnt == n ? "not " : "", cent_test, rad_test,
          vtostr(str, cent, " "), rad);
  return 1;
}

int find_can_centre(Geometry &geom, char type, double &rad, Vec3d &cent,
                    bool invert, int n, const double &eps)
{
  bool use_e = false;
  bool use_vf = false;
  if (type == 'e')
    use_e = true;
  else if (type == 'v')
    use_vf = true;
  else {
    use_e = true;
    use_vf = true;
  }

  vector<vector<int>> g_edges, d_edges;
  geom.get_impl_edges(g_edges);
  int e_sz = g_edges.size();
  if (!cent.is_set())
    cent = geom.centroid();
  if (fabs(rad) < epsilon) {
    GeometryInfo rep(geom);
    rep.set_center(cent);
    rad = (1 - 2 * (rad < 0)) * rep.iedge_dist_lims().sum / e_sz;
  }
  Geometry dual;
  get_dual(dual, geom, rad, cent);
  if (invert)
    dual.transform(Trans3d::translate(cent) * Trans3d::inversion() *
                   Trans3d::translate(-cent));
  dual.get_impl_edges(d_edges);
  Vec3d cur_cent = cent;
  double cur_rad = rad;
  double cent_test = -1, rad_test = -1;
  int cnt;
  for (cnt = 0; cnt < n; cnt++) {
    Vec3d cent_calc(0, 0, 0);
    double rad_g = 0, rad_d = 0;
    if (use_e) {
      Vec3d e_sum_g(0, 0, 0);
      Vec3d e_sum_d(0, 0, 0);
      double rad_sum_g = 0, rad_sum_d = 0;
      for (int e = 0; e < e_sz; ++e) {
        Vec3d Pg = nearest_point(cur_cent, geom.verts(), g_edges[e]);
        rad_sum_g += (Pg - cur_cent).len();
        e_sum_g += Pg;
        Vec3d Pd = nearest_point(cur_cent, dual.verts(), d_edges[e]);
        rad_sum_d += (Pd - cur_cent).len();
        e_sum_d += Pd;
      }
      cent_calc += (e_sum_g + e_sum_d) / (2.0 * e_sz);
      rad_g += rad_sum_g;
      rad_d += rad_sum_d;
    }
    if (use_vf) {
      const vector<Vec3d> &g_verts = geom.verts();
      const vector<Vec3d> &d_verts = dual.verts();
      int v_sz = g_verts.size();
      int f_sz = geom.faces().size();
      Vec3d f_sum_g(0, 0, 0);
      Vec3d v_sum_d(0, 0, 0);
      Vec3d v_sum_g(0, 0, 0);
      Vec3d f_sum_d(0, 0, 0);
      double v_rad_sum_g = 0, v_rad_sum_d = 0;
      double f_rad_sum_g = 0, f_rad_sum_d = 0;
      for (int f = 0; f < f_sz; ++f) {
        Vec3d Pg = nearest_point(cur_cent, geom.verts(), geom.faces(f));
        f_rad_sum_g += (Pg - cur_cent).len();
        f_sum_g += Pg;
        Vec3d Pd = d_verts[f];
        v_rad_sum_d += (Pd - cur_cent).len();
        v_sum_d += Pd;
      }

      for (int v = 0; v < v_sz; ++v) {
        Vec3d Pg = g_verts[v];
        v_rad_sum_g += (Pg - cur_cent).len();
        v_sum_g += Pg;
        Vec3d Pd = nearest_point(cur_cent, dual.verts(), dual.faces(v));
        f_rad_sum_d += (Pd - cur_cent).len();
        f_sum_d += Pd;
      }

      Vec3d vf_avg_g = (v_sum_g + f_sum_g) / double(v_sz + f_sz);
      Vec3d vf_avg_d = (v_sum_d + f_sum_d) / double(v_sz + f_sz);
      Vec3d vf_cent = (vf_avg_g + vf_avg_d) / 2.0;
      cent_calc += vf_cent;
      rad_g += v_rad_sum_g * f_rad_sum_g;
      rad_d += v_rad_sum_d * f_rad_sum_d;
    }
    if (use_e && use_vf) {
      cent_calc /= 2.0;
    }

    if (is_even(cnt)) {
      cur_rad = rad;
      // rad = sqrt(rad_g/rad_d)*cur_rad;
      rad = 0.5 * cur_rad + 0.5 * sqrt(rad_g / rad_d) * cur_rad;
    }
    else {
      cur_cent = cent;
      cent = 0.5 * cur_cent + 0.5 * cent_calc;
    }

    cent_test = (cur_cent - cent).len() / rad;
    rad_test = (cur_rad - rad) / rad;
    // cent.dump("cent");
    // fprintf(stderr, "rad=%g, delt_rad=%g,delta_cent=%g\n", rad, cur_rad-rad,
    //      (cur_cent-cent).len());
    if (fabs(cent_test) < eps && fabs(rad_test) < eps)
      break;
    get_pol_recip_verts(dual, geom, rad, cur_cent);
    if (invert)
      dual.transform(Trans3d::translate(cur_cent) * Trans3d::inversion() *
                     Trans3d::translate(-cur_cent));
  }
  char str[MSG_SZ];
  fprintf(stderr, "[n=%d, limit=%sachieved, r_test=%g, c_test=%g]\n"
                  "centre=(%s), radius=%.16g%s\n",
          cnt, cnt == n ? "not " : "", cent_test, rad_test,
          vtostr(str, cent, " "), rad, (invert) ? "i" : "");
  //}

  return 1;
}

Vec3d find_circumcenter(const Geometry &geom)
{
  auto vec = Vec4d::zero;
  auto mat = Trans3d::zero();
  mat[0] = geom.verts().size();
  for (const auto &v : geom.verts()) {
    mat[1] += v[0];
    mat[2] += v[1];
    mat[3] += v[2];

    mat[4] += v[0];
    mat[5] += v[0] * v[0];
    mat[6] += v[0] * v[1];
    mat[7] += v[0] * v[2];

    mat[8] += v[1];
    mat[9] += v[1] * v[0];
    mat[10] += v[1] * v[1];
    mat[11] += v[1] * v[2];

    mat[12] += v[2];
    mat[13] += v[2] * v[0];
    mat[14] += v[2] * v[1];
    mat[15] += v[2] * v[2];

    double len = v.len2();
    vec[0] += len;
    vec[1] += len * v[0];
    vec[2] += len * v[1];
    vec[3] += len * v[2];
  }

  auto res = mat.inverse() * vec;
  Vec3d cent;
  cent[0] = res[1] / 2;
  cent[1] = res[2] / 2;
  cent[2] = res[3] / 2;
  double rad = sqrt(cent.len2() + res[0]);
  char str[MSG_SZ];
  fprintf(stderr, "Circumsphere: centre=(%s), radius=%.16g\n",
          vtostr(str, cent, " "), rad);
  return cent;
}

/*
Finds least squares intersection point of edge lines, not midcentre!
Vec3d find_NOTmidcenter(const Geometry &geom)
{
  auto vec = Vec3d::zero;
  auto mat = Trans3d::zero();
  mat[15] = 1; // going to use 3x3 matrix, not row 4 or column 4
  map<vector<int>, vector<int>> e2f;
  geom.get_edge_face_pairs(e2f);
  for (const auto &edge : e2f) {
    Vec3d P = geom.verts(edge.first[0]);
    Vec3d v = geom.edge_vec(edge.first).unit();
    mat[0] += v[0] * v[0] - 1;
    mat[1] += v[0] * v[1];
    mat[2] += v[0] * v[2];

    mat[4] += v[1] * v[0];
    mat[5] += v[1] * v[1] - 1;
    mat[6] += v[1] * v[2];

    mat[8] += v[2] * v[0];
    mat[9] += v[2] * v[1];
    mat[10] += v[2] * v[2] - 1;

    vec[0] += vdot(P, Vec3d(v[0] * v[0] - 1, v[0] * v[1],     v[0] * v[2]));
    vec[1] += vdot(P, Vec3d(v[1] * v[0],     v[1] * v[1] - 1, v[1] * v[2]));
    vec[2] += vdot(P, Vec3d(v[2] * v[0],     v[2] * v[1],     v[2] * v[2] - 1));
  }
  auto cent = mat.inverse() * vec;
  char str[MSG_SZ];
  fprintf(stderr, "Midcentre: centre=(%s)\n", vtostr(str, cent, " "));
  return cent;
}
*/

void find_recip_centre(Geometry &geom, char type, double &rad, Vec3d &cent,
                       int n, const double &eps)
{
  if (fabs(rad) < epsilon)
    rad = epsilon / 2.0;
  bool invert = strchr("EVA", type);
  switch (type) {
  case 'x': // default is C
  case 'C':
    cent = geom.centroid();
    break;
  case 'R':
    cent = find_circumcenter(geom);
    break;
  case 'M':
    find_mid_centre(geom, rad, cent, n, eps);
    break;
  case 'e':
    find_can_centre(geom, 'e', rad, cent, invert, n, eps);
    break;
  case 'E':
    rad = -rad;
    find_can_centre(geom, 'e', rad, cent, invert, n, eps);
    break;
  case 'v':
    find_can_centre(geom, 'v', rad, cent, invert, n, eps);
    break;
  case 'V':
    rad = -rad;
    find_can_centre(geom, 'v', rad, cent, invert, n, eps);
    break;
  case 'a':
    find_can_centre(geom, 'a', rad, cent, invert, n, eps);
    break;
  case 'A':
    rad = -rad;
    find_can_centre(geom, 'a', rad, cent, invert, n, eps);
    break;
  }
}

double find_recip_rad(Geometry &geom, char type, Vec3d cent,
                      vector<int> vidxs = vector<int>())
{
  double rad = 1;
  GeometryInfo rep(geom);
  rep.set_center(cent);
  switch (type) {
  case 'v':
    rad = rep.vert_dist_lims().min;
    break;
  case 'V':
    rad = rep.vert_dist_lims().max;
    break;
  case 'e':
    rad = rep.iedge_dist_lims().min;
    break;
  case 'E':
    rad = rep.iedge_dist_lims().max;
    break;
  case 'f':
    rad = rep.face_dist_lims().min;
    break;
  case 'F':
    rad = rep.face_dist_lims().max;
    break;
  case 'S':
    rad = (nearest_point(cent, geom.verts(), vidxs) - cent).len();
    break;
  }
  return rad;
}

bool strip_free_elements(Geometry &geom)
{
  vector<vector<int>> implicit_edges;
  geom.get_impl_edges(implicit_edges);

  vector<int> free_edges;
  for (unsigned int e_idx = 0; e_idx < geom.edges().size(); e_idx++) {
    if (find_edge_in_edge_list(implicit_edges, geom.edges(e_idx)) == -1)
      free_edges.push_back(e_idx);
  }
  if (free_edges.size())
    geom.del(EDGES, free_edges);

  GeometryInfo info(geom);
  auto free_verts = geom.get_info().get_free_verts();
  if (free_verts.size())
    geom.del(VERTS, free_verts);

  return free_edges.size() || free_verts.size();
}

bool is_polyhedron(const Geometry &geom, char *errmsg = nullptr)
{
  GeometryInfo info(geom);
  if (!info.num_faces()) {
    if (errmsg)
      strcpy_msg(errmsg, "not a polyhedron, has no faces");
    return false;
  }

  if (!info.is_closed()) {
    if (errmsg)
      strcpy_msg(errmsg, "not a polyhedron, is not closed");
    return false;
  }

  if (!info.is_known_connectivity()) {
    if (errmsg)
      strcpy_msg(errmsg, "unknown whether a polyhedron, ambiguous "
                         "connectivity");
    return false;
  }

  return true; // each edge shared by exactly 2 faces
}

int main(int argc, char *argv[])
{
  pr_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);

  char errmsg[MSG_SZ];
  if (!is_polyhedron(geom, errmsg))
    opts.error(errmsg, "input_file");

  if (strip_free_elements(geom))
    opts.warning("stripped vertices or edges which were not part of any face",
                 "input_file");

  // Set up init_centre
  double tmp_rad = 0;
  if (opts.init_cent_type == 'C')
    opts.init_cent = geom.centroid();
  else if (opts.init_cent_type == 'M')
    find_mid_centre(geom, tmp_rad, opts.init_cent, opts.num_iters,
                    opts.epsilon);

  if (opts.recip_rad_type == 'x') {
    if (strchr("CMx", opts.recip_cent_type))
      opts.recip_rad_type = 'e';
    else if (strchr("R", opts.recip_cent_type))
      opts.recip_rad_type = 'v';
  }

  Vec3d centre;
  if (opts.recip_cent_type == 'c') {
    centre = opts.centre;
    opts.init_rad = 1.0;
  }
  else {
    centre = opts.init_cent;
    find_recip_centre(geom, opts.recip_cent_type, opts.init_rad, centre,
                      opts.num_iters, opts.epsilon);
  }

  double radius;
  if (opts.recip_rad_type == 'r')
    radius = opts.recip_rad;
  else if (opts.recip_rad_type == 'x')
    radius = opts.init_rad;
  else
    radius =
        find_recip_rad(geom, opts.recip_rad_type, centre, opts.space_verts);

  Geometry dual;
  get_dual(dual, geom, radius, centre, 1.01 * opts.inf);
  if (opts.invert)
    dual.transform(Trans3d::translate(centre) * Trans3d::inversion() *
                   Trans3d::translate(-centre));

  vector<int> invalid_verts;
  for (unsigned int i = 0; i < dual.verts().size(); i++) {
    if (!dual.verts(i).is_set()) {
      dual.del(VERTS, (int)i);
      int idx = invalid_verts.size() + i;
      invalid_verts.push_back(idx);
      i--;
    }
  }
  if (invalid_verts.size()) {
    string msg(
        "removed invalid vertices (and associated faces) with indices - ");
    for (unsigned int i = 0; i < invalid_verts.size() - 1; i++) {
      snprintf(errmsg, MSG_SZ, "%d,", invalid_verts[i]);
      msg += string(errmsg);
    }
    snprintf(errmsg, MSG_SZ, "%d", invalid_verts.back());
    msg += string(errmsg);
    opts.warning(msg);
  }

  if (opts.extra_ideal_elems)
    add_extra_ideal_elems(dual, centre, 1.005 * opts.inf);

  for (const auto &kv : dual.colors(FACES).get_properties())
    if (kv.second.is_invisible()) {
      opts.warning("dual includes invisible faces (base model included "
                   "invisible vertices)");
      break;
    }

  dual.orient();
  if (opts.append)
    geom.append(dual);

  const Geometry &geom_out = (opts.append) ? geom : dual;
  opts.write_or_error(geom_out, opts.ofile);

  return 0;
}
