/*
   Copyright (c) 2003-2020, Adrian Rossiter, Roger Kaufman
   Includes ideas and algorithms by George W. Hart, http://www.georgehart.com

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
   Name: canonical.cc
   Description: canonicalize a polyhedron
                Uses George Hart's two canonicalization algorithm
                http://library.wolfram.com/infocenter/Articles/2012/
                http://www.georgehart.com/virtual-polyhedra/conway_notation.html
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"

#include <cctype>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using std::string;
using std::vector;

using namespace anti;

class cn_opts : public ProgramOpts {
public:
  IterationControl it_ctrl;

  string ifile;
  string ofile;

  char centering;
  char initial_radius;
  char edge_distribution;
  string shuffle_model_indexes;
  char planarize_method;
  int num_iters_planar;
  char canonical_method;
  int num_iters_canonical;
  double edge_factor;
  double plane_factor;
  bool alternate_algorithm;
  double radius_range_percent;
  string output_parts;
  int face_opacity;
  double offset;
  int roundness;
  char normal_type;

  double epsilon;

  Color ipoints_col;
  Color base_nearpts_col;
  Color dual_nearpts_col;
  Color base_edge_col;
  Color dual_edge_col;
  Color sphere_col;

  cn_opts()
      : ProgramOpts("canonical"), centering('e'), initial_radius('e'),
        edge_distribution('\0'), planarize_method('\0'), num_iters_planar(-1),
        canonical_method('m'), num_iters_canonical(-1), edge_factor(DBL_MAX),
        plane_factor(DBL_MAX), alternate_algorithm(false),
        radius_range_percent(-1), output_parts("b"), face_opacity(-1),
        offset(0), roundness(8), normal_type('n'), epsilon(0),
        ipoints_col(Color(255, 255, 0)), base_nearpts_col(Color(255, 0, 0)),
        dual_nearpts_col(Color(0.0, 0.39216, 0.0)), base_edge_col(Color()),
        dual_edge_col(Color()), sphere_col(Color(255, 255, 255))
  {
    it_ctrl.set_max_iters(-1);
    it_ctrl.set_status_checks("1000,1");
    it_ctrl.set_sig_digits(int(-log(::epsilon) / log(10) + 0.5));
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

void cn_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] [input_file]

Read a polyhedron from a file in OFF format. Canonicalize or planarize it.
Uses algorithms by George W. Hart, http://www.georgehart.com/
http://www.georgehart.com/virtual-polyhedra/conway_notation.html
http://www.georgehart.com/virtual-polyhedra/canonical.html
If input_file is not given the program reads from standard input.

Options
%s
  -e <opt>  edge distribution (default : none)
               s - project vertices onto a sphere
  -s <opt>  shuffle model indexes
               v - vertices, e - edges, f - faces, a - all (default: none)
  -r <opt>  initial radius
               e - average edge near points radius = 1 (default)
               v - average vertex radius = 1
               x - not changed
  -C <opt>  initial centering
               e - edge near points centroid (default)
               v - vertex centroid
               x - not moved
  -p <opt>  planarization (done before canoncalization. default: none)
               q - face centroids magnitude squared
               m - mathematica planarize
               a - sand and fill planarize
               p - fast planarize (poly_form -a p)
  -i <itrs> maximum planarize iterations. -1 for unlimited (default: %d)
            WARNING: unstable models may not finish unless -i is set
  -c <opt>  canonicalization
               m - mathematica version (default)
               b - base/dual version (reciprocate on face normals)
               a - moving edge version
               x - none (default, if -p is set)
  -n <itrs> maximum canonical iterations. -1 for unlimited (default: %d)
            WARNING: unstable models may not finish unless -n is set
  -O <args> output b - base, d - dual, i - intersection points (default: b)
               n - base edge near points, m - dual edge near points
               p - base near points centroid, q - dual near points centroid
               u - minimum tangent sphere, U - maximum; o - origin point
               s - base incircles, S - rings; t - dual incircles, T - rings
  -q <dist> offset for incircles to avoid coplanarity e.g 0.0001 (default: 0)
  -g <opt>  roundness of tangent sphere, positive integer n (default: 8)
  -x <opt>  Normals: n - Newell's, t - triangles, q - quads (default: Newell's)
  -d <perc> radius test. percent difference between minimum and maximum radius
               checks if polyhedron is collapsing. 0 for no test 
               (default: 80 for canonicalizing, not used for planarizing)
  -l <lim>  minimum distance change to terminate, as negative exponent
               (default: %d giving %.0e)
            WARNING: high values can cause non-terminal behaviour. Use -n
  -z <nums> number of iterations between status reports (implies termination
            check) (0 for final report only, -1 for no report), optionally
            followed by a comma and the number of iterations between
            termination checks (0 for report checks only) (default: %d,%d)
  -o <file> write output to file (default: write to standard output)

Extra Options
  for (-c m, -p m, and -p p)
  -E <perc> percentage to scale the edge tangency (default: 50) 
  -P <perc> percentage to scale the face planarity (default: 20) (also -p p)
  -A        alternate algorithm. try if imbalance in result (-c m only)

Coloring Options (run 'off_util -H color' for help on color formats)
  -I <col>  intersection points and/or origin color (default: yellow)
  -N <col>  base near points, centroid, incircles color (default: red)
  -M <col>  dual near points, centroid, incircles color (default: darkgreen)
  -B <col>  base edge color (default: unchanged)
  -D <col>  dual edge color (default: unchanged)
  -U <col>  unit sphere color (default: white)
  -T <tran> base/dual transparency. range from 0 (invisible) to 255 (opaque)

)",
          prog_name(), help_ver_text, it_ctrl.get_max_iters(),
          it_ctrl.get_max_iters(), it_ctrl.get_sig_digits(),
          it_ctrl.get_test_val(), it_ctrl.get_status_check_and_report_iters(),
          it_ctrl.get_status_check_only_iters());
}

void cn_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;
  int num;

  bool p_set = false;
  bool c_set = false;

  handle_long_opts(argc, argv);

  while ((c = getopt(
              argc, argv,
              ":hC:r:e:s:p:i:c:n:O:q:g:E:P:Ad:x:z:I:N:M:B:D:U:T:l:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'e':
      if (strlen(optarg) == 1 && strchr("s", int(*optarg)))
        edge_distribution = *optarg;
      else
        error("edge_distribution method type must be s", c);
      break;

    case 's':
      if (strspn(optarg, "vefa") != strlen(optarg))
        error(msg_str("suffle parts are '%s' must be any or all from "
                      "v, e, f, a",
                      optarg),
              c);
      shuffle_model_indexes = optarg;
      break;

    case 'r':
      if (strlen(optarg) == 1 && strchr("evx", int(*optarg)))
        initial_radius = *optarg;
      else
        error("starting radius type must be e, v, x", c);
      break;

    case 'C':
      if (strlen(optarg) == 1 && strchr("evx", int(*optarg)))
        centering = *optarg;
      else
        error("centering method type must be e, v, x", c);
      break;

    case 'p':
      p_set = true;
      if (strlen(optarg) == 1 && strchr("qmap", int(*optarg)))
        planarize_method = *optarg;
      else
        error("planarize method type must be q, m, a, p", c);
      break;

    case 'i':
      print_status_or_exit(read_int(optarg, &num_iters_planar), c);
      if (num_iters_planar < -1)
        error("number of iterations for planarization must be -1 or greater",
              c);
      break;

    case 'c':
      c_set = true;
      if (strlen(optarg) == 1 && strchr("mbax", int(*optarg)))
        canonical_method = *optarg;
      else
        error("canonical method type must be m, b, a, x", c);
      break;

    case 'n':
      print_status_or_exit(read_int(optarg, &num_iters_canonical), c);
      if (num_iters_canonical < -1)
        error("number of iterations for canonical must be -1 or greater", c);
      break;

    case 'O':
      if (strspn(optarg, "bdinmopqsStTuU") != strlen(optarg))
        error(msg_str("output parts are '%s' must be any or all from "
                      "b, d, i, n, m, o, p, q, s, S, t, T, u, U",
                      optarg),
              c);
      output_parts = optarg;
      break;

    case 'q':
      print_status_or_exit(read_double(optarg, &offset), c);
      break;

    case 'g':
      print_status_or_exit(read_int(optarg, &roundness), c);
      break;

    case 'x':
      if (strlen(optarg) == 1 && strchr("ntq", int(*optarg)))
        normal_type = *optarg;
      else
        error("normal type must be n, t, q", c);
      break;

    case 'E':
      print_status_or_exit(read_double(optarg, &edge_factor), c);
      if (edge_factor <= 0 || edge_factor >= 100)
        warning("edge factor not inside range 0 to 100", c);
      break;

    case 'P':
      print_status_or_exit(read_double(optarg, &plane_factor), c);
      if (plane_factor <= 0 || plane_factor >= 100) {
        warning("plane factor not inside range 0 to 100", c);
      }
      break;

    case 'A':
      alternate_algorithm = true;
      break;

    case 'd':
      print_status_or_exit(read_double(optarg, &radius_range_percent), c);
      if (radius_range_percent < 0)
        error("percentage must not be less than 0", c);
      break;

    case 'z':
      print_status_or_exit(it_ctrl.set_status_checks(optarg), c);
      break;

    case 'I':
      print_status_or_exit(ipoints_col.read(optarg), c);
      break;

    case 'N':
      print_status_or_exit(base_nearpts_col.read(optarg), c);
      break;

    case 'M':
      print_status_or_exit(dual_nearpts_col.read(optarg), c);
      break;

    case 'B':
      print_status_or_exit(base_edge_col.read(optarg), c);
      break;

    case 'D':
      print_status_or_exit(dual_edge_col.read(optarg), c);
      break;

    case 'U':
      print_status_or_exit(sphere_col.read(optarg), c);
      break;

    case 'T':
      print_status_or_exit(read_int(optarg, &face_opacity), c);
      if (face_opacity < 0 || face_opacity > 255) {
        error("face transparency must be between 0 and 255", c);
      }
      break;

    case 'l':
      print_status_or_exit(read_int(optarg, &num), c);
      print_status_or_exit(it_ctrl.set_sig_digits(num), c);
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  // if planarizing only do not canonicalize
  if (p_set && !c_set)
    canonical_method = 'x';

  if (alternate_algorithm && canonical_method != 'm')
    warning("alternate form only has effect in mathematica canonicalization",
            'A');

  // set default variables
  if (canonical_method == 'm' || planarize_method == 'm') {
    if (edge_factor == DBL_MAX)
      edge_factor = 50.0;
  }
  else if (edge_factor != DBL_MAX)
    warning("set, but not used for this algorithm", 'E');

  if (canonical_method == 'm' || planarize_method == 'm' ||
      planarize_method == 'p') {
    if (plane_factor == DBL_MAX)
      plane_factor = 20.0;
  }
  else if (plane_factor != DBL_MAX)
    warning("set, but not used for this algorithm", 'P');

  if ((canonical_method == 'x') && planarize_method &&
      (radius_range_percent > -1))
    warning("set, but not used for planarization", 'd');

  if (argc - optind > 1)
    error("too many arguments");

  if (argc - optind == 1)
    ifile = argv[optind];

  epsilon = it_ctrl.get_test_val();
}

// RK - average radius rather than maximum has more reliability than max
void unitize_vertex_radius(Geometry &geom)
{
  GeometryInfo info(geom);
  info.set_center(geom.centroid());
  // geom.transform(Trans3d::scale(1 / info.vert_dist_lims().max));
  double avg = info.vert_dist_lims().sum / info.num_verts();
  geom.transform(Trans3d::scale(1 / avg));
}

/*
// plane a single face aligned to the z axis
void plane_face(Geometry &polygon)
{
  Vec3d face_normal = polygon.face_norm(0).unit();
  Vec3d face_centroid = polygon.face_cent(0);
  // make sure face_normal points outward
  if (vdot(face_normal, face_centroid) < 0)
    face_normal *= -1.0;

  // this gives the same results (from the mathematica algorithm)
  // for (auto &vert : polygon.raw_verts())
  //  vert += vdot(face_normal, face_centroid - vert) * face_normal;
  // return;

  // rotate face to z axis
  Trans3d trans = Trans3d::rotate(face_normal, Vec3d(0, 0, 1));
  polygon.transform(trans);

  // refresh face centroid
  face_centroid = polygon.face_cent(0);

  // set z of all vertices to height of face centroid
  for (auto &vert : polygon.raw_verts())
    vert[2] = face_centroid[2];

  // rotate face back to original position
  polygon.transform(trans.inverse());
}
*/

void planarity_info(Geometry &geom)
{
  GeometryInfo rep(geom);

  double max_nonplanar = 0;
  double sum_nonplanar = 0;
  unsigned int sz = geom.faces().size();
  for (unsigned int i = 0; i < sz; i++) {
    double nonplanar = rep.get_f_max_nonplanars()[i];
    sum_nonplanar += nonplanar;
    if (nonplanar > max_nonplanar)
      max_nonplanar = nonplanar;
  }
  fprintf(stderr, "\n");
  fprintf(stderr, "maximum_nonplanarity = %.17g\n", max_nonplanar);
  fprintf(stderr, "average_nonplanarity = %.17g\n", sum_nonplanar / sz);
  //  fprintf(stderr, "isoperimetric quotient = %.17g\n",
  //  rep.isoperimetric_quotient());
  fprintf(stderr, "\n");
}

void midradius_info(Geometry &geom, const bool completed)
{
  double min = 0;
  double max = 0;
  Vec3d center;
  double radius = edge_nearpoints_radius(geom, min, max, center);
  fprintf(stderr, "midradius = %.17g (range: %.15g to %.15g)\n", radius, min,
          max);
  fprintf(stderr, "midcenter is the origin\n");
  fprintf(stderr, "near point centroid = (%.17g,%.17g,%.17g)\n", center[0],
          center[1], center[2]);
  double epsilon_local = 1e-8;
  if (!completed)
    fprintf(stderr, "Warning: the calculation did not complete\n");
  else if (double_ne(center[0], 0.0, epsilon_local) ||
           double_ne(center[1], 0.0, epsilon_local) ||
           double_ne(center[2], 0.0, epsilon_local))
    fprintf(stderr,
            "Warning: the result is NOT canonical. edge tangency only\n");
  else
    fprintf(stderr, "the result is canonical\n");
  fprintf(stderr, "\n");
}

void generate_points(const Geometry &base, const Geometry &dual,
                     vector<Vec3d> &ips, vector<Vec3d> &base_nearpts,
                     vector<Vec3d> &dual_nearpts, const cn_opts &opts)
{
  vector<vector<int>> base_edges;
  vector<vector<int>> dual_edges;

  base.get_impl_edges(base_edges);
  dual.get_impl_edges(dual_edges);

  for (auto &base_edge : base_edges)
    base_nearpts.push_back(base.edge_nearpt(base_edge, Vec3d(0, 0, 0)));

  for (auto &dual_edge : dual_edges)
    dual_nearpts.push_back(dual.edge_nearpt(dual_edge, Vec3d(0, 0, 0)));

  if (opts.output_parts.find("i") != string::npos) {
    for (unsigned int i = 0; i < base_edges.size(); i++) {
      Vec3d b0 = base.verts(base_edges[i][0]);
      Vec3d b1 = base.verts(base_edges[i][1]);
      for (unsigned int j = 0; j < dual_edges.size(); j++) {
        Vec3d d0 = dual.verts(dual_edges[j][0]);
        Vec3d d1 = dual.verts(dual_edges[j][1]);
        // does base edge intersect with dual edge
        // use local epsilon
        double epsilon_local = 1e-12;
        Vec3d intersection_point =
            segments_intersection(b0, b1, d0, d1, epsilon_local);
        if (intersection_point.is_set()) {
          ips.push_back(intersection_point);
          break;
        }
      }
    }

    if (ips.size() != base_edges.size()) {
      if (opts.canonical_method != 'x')
        fprintf(stderr,
                "Warning: only %d out of %d intersection points found. try "
                "more precision\n",
                (int)ips.size(), (int)base_edges.size());
      else
        fprintf(stderr,
                "Warning: only canonical models have intersection points\n");
    }
  }
}

void set_edge_colors(Geometry &geom, const Color col)
{
  // it is possible unset faces already exist
  if (col.is_set()) {
    geom.add_missing_impl_edges();
    Coloring(&geom).e_one_col(col);
  }
}

/*
Vec3d face_edge_nearpoints_centroid(const Geometry &geom, double &radius, const
int face_no)
{
  vector<Vec3d> near_pts;

  unsigned int fsz = geom.faces(face_no).size();
  for (unsigned int i = 0; i < fsz; i++) {
    int v1 = geom.faces(face_no)[i];
    int v2 = geom.faces(face_no)[(i + 1) % fsz];

    Vec3d P = geom.edge_nearpt(make_edge(v1,v2), Vec3d(0, 0, 0));
    near_pts.push_back(P);
  }

  Vec3d face_edge_nearpt_centroid = centroid(near_pts);
//  face_edge_nearpt_centroid = geom.face_cent(face_no);

  // get the minimum radius
  double min = DBL_MAX;
  for (unsigned int i = 0; i < fsz; i++) {
    double l = (face_edge_nearpt_centroid - near_pts[i]).len();
    if (l < min)
      min = l;
  }
  radius = min;

  return face_edge_nearpt_centroid;
}

vector<Vec3d> face_edge_nearpoints_centroids(const Geometry &base)
{
  double rad;
  Geometry near_pts;
  unsigned int fsz = base.faces().size();
  for (unsigned int i = 0; i < fsz; i++)
    near_pts.add_vert(face_edge_nearpoints_centroid(base, rad, i));

  return(near_pts.verts());
}
*/

Geometry unit_circle(int polygon_size, const Color &incircle_color, bool filled)
{
  Geometry incircle;
  double arc = deg2rad(360.0 / (double)polygon_size);

  double angle = 0.0;
  for (int i = 0; i < polygon_size; i++) {
    incircle.add_vert(Vec3d(cos(angle), sin(angle), 0.0), Color::invisible);
    incircle.add_edge(make_edge(i, (i + 1) % polygon_size),
                      (filled ? Color::invisible : incircle_color));
    angle += arc;
  }

  if (filled) {
    vector<int> face;
    for (int i = 0; i < polygon_size; i++)
      face.push_back(i);
    incircle.add_face(face, incircle_color);
  }

  return (incircle);
}

// get the minimum incircle radius. for canonicalized they are all equal
double incircle_radius(const Geometry &geom, Vec3d &center, const int face_no)
{
  vector<Vec3d> near_pts;

  // get the near points to measure length from center
  unsigned int fsz = geom.faces(face_no).size();
  for (unsigned int i = 0; i < fsz; i++) {
    int v1 = geom.faces(face_no)[i];
    int v2 = geom.faces(face_no)[(i + 1) % fsz];

    Vec3d P = geom.edge_nearpt(make_edge(v1, v2), Vec3d(0, 0, 0));
    near_pts.push_back(P);
  }

  // get the minimum radius
  double radius = DBL_MAX;
  for (unsigned int i = 0; i < fsz; i++) {
    double l = (center - near_pts[i]).len();
    if (l < radius)
      radius = l;
  }

  return radius;
}

Geometry incircles(const Geometry &geom, const Color &incircle_color,
                   bool filled, double offset)
{
  Geometry incircles;
  Geometry circle = unit_circle(60, incircle_color, filled);

  for (unsigned int i = 0; i < geom.faces().size(); i++) {
    // find incircle rotation place
    Vec3d face_centroid = anti::centroid(geom.verts(), geom.faces(i));
    Vec3d face_normal = face_norm(geom.verts(), geom.faces(i)).unit();
    Vec3d center = face_normal * vdot(face_centroid, face_normal);

    // find radius of incircle, and make incircle of radius
    Geometry incircle = circle;
    double radius = incircle_radius(geom, center, i);
    incircle.transform(Trans3d::scale(radius));

    // set depth of incircle
    double depth = center.len() + offset;
    incircle.transform(Trans3d::translate(Vec3d(0, 0, depth)));

    // rotate incircle into place
    incircle.transform(Trans3d::rotate(Vec3d(0, 0, 1), center));
    incircles.append(incircle);
  }
  return incircles;
}

// RK - models can change size drastically and
// RK - set radius to 1 for get_dual call, if necessary
void reset_model_size(Geometry &base, const cn_opts &opts)
{
  double radius = edge_nearpoints_radius(base);
  if (double_ne(radius, 1.0, opts.epsilon))
    unitize_nearpoints_radius(base);
}

void construct_model(Geometry &base, const cn_opts &opts)
{
  // get statistics before model is changed
  double min = 0;
  double max = 0;
  Vec3d center;
  edge_nearpoints_radius(base, min, max, center);

  Geometry dual;
  get_dual(dual, base, 1, Vec3d(0, 0, 0));

  vector<Vec3d> ips;
  vector<Vec3d> base_nearpts;
  vector<Vec3d> dual_nearpts;
  generate_points(base, dual, ips, base_nearpts, dual_nearpts, opts);

  // base incircles
  Geometry base_incircles;
  if (opts.output_parts.find_first_of("sS") != string::npos) {
    bool filled = (opts.output_parts.find("s") != string::npos) ? true : false;
    base_incircles =
        incircles(base, opts.base_nearpts_col, filled, opts.offset);
  }

  // dual incircles
  Geometry dual_incircles;
  if (opts.output_parts.find_first_of("tT") != string::npos) {
    bool filled = (opts.output_parts.find("t") != string::npos) ? true : false;
    dual_incircles =
        incircles(dual, opts.dual_nearpts_col, filled, opts.offset);
  }

  // clear base if not using
  if (opts.output_parts.find("b") == string::npos)
    base.clear_all();
  // set edge colors here
  else
    set_edge_colors(base, opts.base_edge_col);

  // append dual
  if (opts.output_parts.find("d") != string::npos) {
    // set edge colors here
    set_edge_colors(dual, opts.dual_edge_col);
    base.append(dual);
  }

  // add intersection points
  if (opts.output_parts.find("i") != string::npos) {
    for (auto &ip : ips)
      base.add_vert(ip, opts.ipoints_col);
  }

  // add base near points
  if (opts.output_parts.find("n") != string::npos) {
    for (auto &base_nearpt : base_nearpts)
      base.add_vert(base_nearpt, opts.base_nearpts_col);
  }

  // add dual near points
  if (opts.output_parts.find("m") != string::npos) {
    for (auto &dual_nearpt : dual_nearpts)
      base.add_vert(dual_nearpt, opts.dual_nearpts_col);
  }

  // add base near points centroid
  if (opts.output_parts.find("p") != string::npos) {
    base.add_vert(centroid(base_nearpts), opts.base_nearpts_col);
  }

  // add dual near points centroid
  if (opts.output_parts.find("q") != string::npos) {
    base.add_vert(centroid(dual_nearpts), opts.dual_nearpts_col);
  }

  // add origin point
  if (opts.output_parts.find("o") != string::npos) {
    base.add_vert(Vec3d(0, 0, 0), opts.ipoints_col);
  }

  // apply opacity to faces of base and dual
  if (opts.face_opacity > -1) {
    ColorValuesToRangeHsva valmap(
        msg_str("A%g", (double)opts.face_opacity / 255), Color(255, 255, 255));
    valmap.apply(base, FACES);

    for (const auto &kp : base.colors(FACES).get_properties()) {
      if (kp.second.is_index()) {
        opts.warning("map indexes cannot be made transparent", 'T');
        break;
      }
    }
  }

  // add unit sphere on origin
  if (opts.output_parts.find_first_of("uU") != string::npos) {
    Geometry sgeom;
    string geo_str = "geo_" + std::to_string(opts.roundness) + "_" +
                     std::to_string(opts.roundness);
    sgeom.read_resource(geo_str);
    sgeom.transform(Trans3d::translate(-centroid(sgeom.verts())));
    unitize_vertex_radius(sgeom);

    if (opts.output_parts.find("u") != string::npos)
      sgeom.transform(Trans3d::scale(min));
    else
      sgeom.transform(Trans3d::scale(max));

    Coloring(&sgeom).vef_one_col(Color::invisible, Color::invisible,
                                 opts.sphere_col);
    base.append(sgeom);
  }

  if (base_incircles.verts().size())
    base.append(base_incircles);

  if (dual_incircles.verts().size())
    base.append(dual_incircles);
}

void check_model(const Geometry &geom, string s, const cn_opts &opts)
{
  double epsilon_local1 = 1e-4;         // coincident elements
  double epsilon_local2 = opts.epsilon; // face area (default is 1e-12)

  Geometry geom_merged = geom;
  merge_coincident_elements(geom_merged, "vef", 0, epsilon_local1);

  if (geom.verts().size() != geom_merged.verts().size())
    opts.warning(
        msg_str("possible coincident vertices found in %s", s.c_str()));

  if (geom.faces().size() != geom_merged.faces().size())
    opts.warning(msg_str("possible coincident faces found in %s", s.c_str()));

  GeometryInfo info(geom);
  vector<double> areas = info.get_f_areas();
  for (unsigned int i = 0; i < areas.size(); i++) {
    // fprintf(stderr, "%.17lf\n", areas[i]);
    if (areas[i] < epsilon_local2) {
      opts.warning(msg_str(
          "possible faces of near 0 area found in %s (face %d)", s.c_str(), i));
      break;
    }
  }
}

void check_coincidence(const Geometry &geom, const cn_opts &opts)
{
  string s = "base";
  check_model(geom, s, opts);

  Geometry dual;
  get_dual(dual, geom, 1);
  s = "dual";
  check_model(dual, s, opts);
}

vector<int> geom_deal(Geometry &geom, const int pack_size)
{
  Coloring clrng;
  clrng.set_geom(&geom);
  string map_name = "deal" + std::to_string(pack_size);
  ColorMap *pack = colormap_from_name(map_name.c_str());
  clrng.add_cmap(pack);

  vector<int> deal;
  for (int i = 0; i < pack->effective_size(); i++)
    deal.push_back(pack->get_col(i).get_index());

  return (deal);
}

void shuffle_model_indexes(Geometry &geom, const cn_opts &opts)
{
  bool shuffle_verts =
      (opts.shuffle_model_indexes.find_first_of("va") != string::npos);
  bool shuffle_faces =
      (opts.shuffle_model_indexes.find_first_of("fa") != string::npos);
  bool shuffle_edges =
      (opts.shuffle_model_indexes.find_first_of("ea") != string::npos);

  map<int, int> new_verts;
  if (!shuffle_verts) {
    for (unsigned int i = 0; i < geom.verts().size(); i++)
      new_verts[i] = i;
  }
  else {
    fprintf(stderr, "shuffle model indexes: vertices\n");
    vector<int> deal = geom_deal(geom, geom.verts().size());
    vector<Vec3d> shuffled_verts;
    for (unsigned int i = 0; i < geom.verts().size(); i++) {
      int v_new = deal[i];
      new_verts[v_new] = i;
      shuffled_verts.push_back(geom.verts(v_new));
    }
    geom.raw_verts() = shuffled_verts;
  }

  if (shuffle_faces || shuffle_verts) {
    map<int, int> face_order;
    if (!shuffle_faces) {
      for (unsigned int i = 0; i < geom.faces().size(); i++)
        face_order[i] = i;
    }
    else {
      fprintf(stderr, "shuffle model indexes: faces\n");
      vector<int> deal = geom_deal(geom, geom.faces().size());
      for (unsigned int i = 0; i < geom.faces().size(); i++)
        face_order[deal[i]] = i;
    }
    vector<vector<int>> shuffled_faces;
    for (unsigned int i = 0; i < face_order.size(); i++) {
      vector<int> face;
      for (unsigned int j = 0; j < geom.faces(face_order[i]).size(); j++) {
        int v = geom.faces(face_order[i])[j];
        int v_new = new_verts[v];
        face.push_back(v_new);
      }
      shuffled_faces.push_back(face);
    }
    geom.raw_faces() = shuffled_faces;
  }

  if (!geom.edges().size()) {
    fprintf(stderr, "shuffle model indexes: no explicit edges\n");
  }
  else if (shuffle_edges || shuffle_verts) {
    map<int, int> edge_order;
    if (!shuffle_edges)
      for (unsigned int i = 0; i < geom.edges().size(); i++) {
        edge_order[i] = i;
      }
    else {
      fprintf(stderr, "shuffle model indexes: edges\n");
      vector<int> deal = geom_deal(geom, geom.edges().size());
      for (unsigned int i = 0; i < geom.edges().size(); i++)
        edge_order[deal[i]] = i;
    }
    vector<vector<int>> shuffled_edges;
    for (unsigned int i = 0; i < geom.edges().size(); i++) {
      vector<int> edge;
      for (unsigned int j = 0; j < geom.edges(edge_order[i]).size(); j++) {
        int v = geom.edges(edge_order[i])[j];
        int v_new = new_verts[v];
        edge.push_back(v_new);
      }
      shuffled_edges.push_back(edge);
    }
    geom.raw_edges() = shuffled_edges;
  }
}

int main(int argc, char *argv[])
{
  cn_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);

  if (opts.edge_distribution) {
    fprintf(stderr, "edge distribution: project onto sphere\n");
    if (opts.edge_distribution == 's')
      project_onto_sphere(geom);
  }

  if (opts.shuffle_model_indexes.length())
    shuffle_model_indexes(geom, opts);

  fprintf(stderr, "\n");
  fprintf(stderr, "starting radius: ");
  if (opts.initial_radius == 'e') {
    fprintf(stderr, "average edge near points\n");
    unitize_nearpoints_radius(geom);
  }
  else if (opts.centering == 'v') {
    fprintf(stderr, "average vertex\n");
    unitize_vertex_radius(geom);
  }
  else if (opts.centering == 'x')
    fprintf(stderr, "radius not changed\n");

  fprintf(stderr, "centering: ");
  if (opts.centering == 'e') {
    fprintf(stderr, "edge near points centroid to origin\n");
    geom.transform(
        Trans3d::translate(-edge_nearpoints_centroid(geom, Vec3d(0, 0, 0))));
  }
  else if (opts.centering == 'v') {
    fprintf(stderr, "vertex centroid to origin\n");
    geom.transform(Trans3d::translate(-centroid(geom.verts())));
  }
  else if (opts.centering == 'x')
    fprintf(stderr, "model not moved\n");

  fprintf(stderr, "normals: ");
  if (opts.normal_type == 'n')
    fprintf(stderr, "Newell's method\n");
  else if (opts.normal_type == 't')
    fprintf(stderr, "Triangles method\n");
  else if (opts.normal_type == 'q')
    fprintf(stderr, "Quads method\n");

  bool completed = false;
  if (opts.planarize_method) {
    string planarize_str;
    if (opts.planarize_method == 'q')
      planarize_str = "face centroids magnitude squared";
    else if (opts.planarize_method == 'm')
      planarize_str = "mathematica";
    else if (opts.planarize_method == 'a')
      planarize_str = "sand and fill";
    else if (opts.planarize_method == 'p')
      planarize_str = "poly_form -a p";
    fprintf(stderr, "planarize: %s method\n", planarize_str.c_str());

    bool planarize_only = true;
    opts.it_ctrl.set_max_iters(opts.num_iters_planar);
    double radius_range_pct =
        0; // RK: not used when planarizing
           // (opts.radius_range_percent < 0) ? 0 : opts.radius_range_percent;

    if (opts.planarize_method == 'q') {
      completed = canonicalize_bd(geom, opts.it_ctrl, opts.planarize_method,
                                  radius_range_pct / 100, opts.centering,
                                  opts.normal_type);
    }
    else if (opts.planarize_method == 'm') {
      completed = canonicalize_mm(geom, opts.it_ctrl, opts.edge_factor / 100,
                                  opts.plane_factor / 100,
                                  radius_range_pct / 100, opts.normal_type,
                                  opts.alternate_algorithm, planarize_only);
    }
    else if (opts.planarize_method == 'a') {
      completed =
          canonicalize_unit(geom, opts.it_ctrl, radius_range_pct / 100,
                            opts.centering, opts.normal_type, planarize_only);
    }
    else if (opts.planarize_method == 'p') {
      Symmetry dummy;
      Status stat;
      stat = make_planar(geom, opts.it_ctrl, opts.plane_factor / 100, dummy);
      completed = (stat.is_ok() ? true : false);
    }

    // RK - report planarity
    planarity_info(geom);
  }

  if (opts.canonical_method && opts.canonical_method != 'x') {
    string canonicalize_str;
    if (opts.canonical_method == 'm')
      canonicalize_str = "mathematica";
    else if (opts.canonical_method == 'b')
      canonicalize_str = "base/dual";
    else if (opts.canonical_method == 'a')
      canonicalize_str = "moving edge";
    fprintf(stderr, "canonicalize: %s method\n", canonicalize_str.c_str());

    bool planarize_only = false;
    opts.it_ctrl.set_max_iters(opts.num_iters_canonical);
    double radius_range_pct =
        (opts.radius_range_percent < 0) ? 80 : opts.radius_range_percent;

    if (opts.canonical_method == 'm') {
      completed = canonicalize_mm(geom, opts.it_ctrl, opts.edge_factor / 100,
                                  opts.plane_factor / 100,
                                  radius_range_pct / 100, opts.normal_type,
                                  opts.alternate_algorithm, planarize_only);
    }
    else if (opts.canonical_method == 'b') {
      completed = canonicalize_bd(geom, opts.it_ctrl, opts.canonical_method,
                                  radius_range_pct / 100, opts.centering,
                                  opts.normal_type);
    }
    else if (opts.canonical_method == 'a') {
      completed =
          canonicalize_unit(geom, opts.it_ctrl, radius_range_pct / 100,
                            opts.centering, opts.normal_type, planarize_only);
    }

    // RK - report planarity
    planarity_info(geom);

    // RK - print midradius info
    midradius_info(geom, completed);
  }

  // RK - standardize model radius
  reset_model_size(geom, opts);

  // RK - add coincidence checking the model
  check_coincidence(geom, opts);

  // RK - parts to output
  construct_model(geom, opts);

  opts.write_or_error(geom, opts.ofile);

  return 0;
}
