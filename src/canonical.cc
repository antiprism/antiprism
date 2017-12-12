/*
   Copyright (c) 2003-2017, Adrian Rossiter, Roger Kaufman
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

#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>

#include "../base/antiprism.h"

using std::string;
using std::vector;

using namespace anti;

class cn_opts : public ProgramOpts {
public:
  string ifile;
  string ofile;

  char centering;
  char initial_radius;
  char edge_distribution;
  char planarize_method;
  int num_iters_planar;
  char canonical_method;
  int num_iters_canonical;
  double mm_edge_factor;
  double mm_plane_factor;
  bool alternate_algorithm;
  int rep_count;
  double radius_range_percent;
  string output_parts;
  int face_opacity;
  double offset;
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
        canonical_method('m'), num_iters_canonical(-1), mm_edge_factor(50),
        mm_plane_factor(20), alternate_algorithm(false), rep_count(1000),
        radius_range_percent(80), output_parts("b"), face_opacity(-1),
        offset(0), normal_type('n'), epsilon(0),
        ipoints_col(Color(255, 255, 0)), base_nearpts_col(Color(255, 0, 0)),
        dual_nearpts_col(Color(0.0, 0.39216, 0.0)), base_edge_col(Color()),
        dual_edge_col(Color()), sphere_col(Color(255, 255, 255))
  {
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

// clang-format off
void cn_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Read a polyhedron from a file in OFF format. Canonicalize or planarize it.\n"
"Uses algorithms by George W. Hart, http://www.georgehart.com/\n"
"http://www.georgehart.com/virtual-polyhedra/conway_notation.html\n"
"http://www.georgehart.com/virtual-polyhedra/canonical.html\n"
"If input_file is not given the program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"  -e <opt>  edge distribution (default : none)\n"
"               s - project vertices onto a sphere\n"
"  -r <opt>  initial radius\n"
"               e - average edge near points radius = 1 (default)\n"
"               v - average vertex radius = 1\n"
"               x - not changed\n"
"  -C <opt>  initial centering\n"
"               e - edge near points centroid (default)\n"
"               v - vertex centroid\n"
"               x - not moved\n"
"  -p <opt>  planarization (done before canoncalization. default: none)\n"
"               p - face centroids (magnitude squared)\n"
"               q - face centroids (magnitude)\n"
"               f - face centroids\n"
"               m - mathematica planarize\n"
"               a - sand and fill planarize\n"
"               u - make faces into unit-edged regular polygons (minmax -a u)\n"
"  -i <itrs> maximum number of planarize iterations (default: no limit)\n"
"  -c <opt>  canonicalization\n"
"               m - mathematica version (default)\n"
"               b - base/dual version (reciprocate on face normals)\n"
"               a - moving edge version\n"
"               x - none (default, if -p is set)\n"
"  -n <itrs> maximum number of canonical iterations (default: no limit)\n"
"  -O <args> output b - base, d - dual, i - intersection points (default: b)\n"
"               n - base edge near points, m - dual edge near points\n"
"               p - base near points centeroid, q - dual near points centroid\n"
"               u - minimum tangent sphere, U - maximum, o - origin point\n"
"               s - base incircles, S - rings, t -dual incircles, T -rings\n"
"  -q <dist> offset for incircles to avoid coplanarity e.g 0.0001 (default: 0.0)\n"
"  -x <opt>  Normals: n - Newell's, t - triangles, q - quads (default: Newell's)\n"
"  -d <perc> radius test. precent difference between minumum and maximum radius\n"
"               checks if polyhedron is collapsing. 0 for no test (default: 80)\n"
"  -z <n>    status reporting every n lines. -1 for no status. (default: 1000)\n"
"  -l <lim>  minimum distance change to terminate, as negative exponent\n"
"               (default: %d giving %.0e)\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"Mathematica Canonicalize Options (-c m and -p m)\n"
"  -E <perc> percentage to scale the edge tangency error (default: 50)\n" 
"  -P <perc> percentage to scale the face planarity error (default: 20)\n"
"  -A        alterate algorithm. try if imbalance in result (-c m only)\n" 
"\n"
"Coloring Options (run 'off_util -H color' for help on color formats)\n"
"  -I <col>  intersection points and/or origin color (default: yellow)\n"
"  -N <col>  base near points, centroid, incircles color (default: red)\n"
"  -M <col>  dual near points, centroid, incircles color (default: darkgreen)\n"
"  -B <col>  base edge color (default: unchanged)\n"
"  -D <col>  dual edge color (default: unchanged)\n"
"  -U <col>  unit sphere color (default: white)\n"
"  -T <tran> base/dual transparency. range from 0 (invisible) to 255 (opaque)\n"
"\n"
"\n",prog_name(), help_ver_text, int(-log(::epsilon)/log(10) + 0.5), ::epsilon);
}
// clang-format on 

void cn_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  bool p_set = false;
  bool c_set = false;

  int sig_compare = INT_MAX;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hC:r:e:p:i:c:n:O:q:E:P:Ad:x:z:I:N:M:B:D:U:T:l:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'e':
      if (strlen(optarg) == 1 && strchr("s", int(*optarg)))
        edge_distribution = *optarg;
      else
        error("edge_distribution method type must be s", c);
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
      if (strlen(optarg) == 1 && strchr("pqfmau", int(*optarg)))
        planarize_method = *optarg;
      else
        error("planarize method type must be p, q, f, m, a, u", c);
      break;

    case 'i':
      print_status_or_exit(read_int(optarg, &num_iters_planar), c);
      if (num_iters_planar < 0)
        error(
            "number of iterations for preplanarization must be 0 or greater",
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
      if (num_iters_canonical < 0)
        error("number of iterations for canonical must be 0 or greater", c);
      break;

    case 'O':
      if (strspn(optarg, "bdinmopqsStTuU") != strlen(optarg))
        error(msg_str("output parts are '%s' must be any or all from "
                      "b, d, i, n, m, o, p, q, s, S, t, T, u, U", optarg), c);
      output_parts = optarg;
      break;

    case 'q':
      print_status_or_exit(read_double(optarg, &offset), c);
      break;

    case 'x':
      if (strlen(optarg) == 1 && strchr("ntq", int(*optarg)))
        normal_type = *optarg;
      else
        error("normal type must be n, t, q", c);
      break;

    case 'E':
      print_status_or_exit(read_double(optarg, &mm_edge_factor), c);
      if (mm_edge_factor <= 0 || mm_edge_factor >= 100)
        warning("edge factor not inside range 0 to 100", c);
      break;

    case 'P':
      print_status_or_exit(read_double(optarg, &mm_plane_factor), c);
      if (mm_plane_factor <= 0 || mm_plane_factor >= 100) {
        warning("plane factor not inside range 0 to 100", c);
      }
      break;

    case 'A':
      alternate_algorithm = true;
      break;

    case 'd':
      print_status_or_exit(read_double(optarg, &radius_range_percent), c);
      if (radius_range_percent < 0)
        error("percentage must be greater than 0", c);
      break;

    case 'z':
      print_status_or_exit(read_int(optarg, &rep_count), c);
      if (rep_count < -1)
        error("number of iterations must be -1 or greater", c);
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
      print_status_or_exit(read_int(optarg, &sig_compare), c);
      if (sig_compare < 0) {
        warning("limit is negative, and so ignored", c);
      }
      if (sig_compare > DEF_SIG_DGTS) {
        warning("limit is very small, may not be attainable", c);
      }
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
    warning("alternate form only has effect in mathematica canonicalization", 'A');

  if (argc - optind > 1)
    error("too many arguments");

  if (argc - optind == 1)
    ifile = argv[optind];

  epsilon = (sig_compare != INT_MAX) ? pow(10, -sig_compare) : ::epsilon;
}

// RK - average radius rather than maximum has more reliability than max
void unitize_vertex_radius(Geometry &geom)
{
  GeometryInfo info(geom);
  info.set_center(geom.centroid());
  //geom.transform(Trans3d::scale(1 / info.vert_dist_lims().max));
  double avg = info.vert_dist_lims().sum / info.num_verts();
  geom.transform(Trans3d::scale(1 / avg));
}

// P and Q are modified
void move_line_to_point(Vec3d &P, Vec3d &Q, const Vec3d &X)
{
  Vec3d Y = X + (Q-P);
  Vec3d V = P + X;
  Vec3d P2 = lines_intersection(P, V, X, Y, 0);
  if (P2.is_set()) {
    Q += P2 - P;
    P = P2;
  }
}

// plane a single face aligned to the z axis
void plane_face(Geometry &polygon)
{
  Vec3d face_normal = polygon.face_norm(0).unit();
  Vec3d face_centroid = polygon.face_cent(0);
  // make sure face_normal points outward
  if (vdot(face_normal, face_centroid) < 0)
    face_normal *= -1.0;

  // this gives the same results (from the mathematica algorithm)
  //for (auto &vert : polygon.raw_verts())
  //  vert += vdot(face_normal, face_centroid - vert) * face_normal; 
  //return;

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

// RK - edge near points of base seek 1
bool canonicalize_unit(Geometry &geom, const int num_iters, const double radius_range_percent,
                      const int rep_count, const char centering, 
                      const char normal_type, const bool planar_only, const double eps)
{
  bool completed = false;

  vector<vector<int>> edges;
  geom.get_impl_edges(edges);

  vector<Vec3d> &verts = geom.raw_verts();

  double max_diff2 = 0;
  unsigned int cnt;
  for (cnt = 0; cnt < (unsigned int)num_iters;) {
    vector<Vec3d> verts_last = verts;

    if (!planar_only) {
      for (auto &edge : edges) {
        // unit near point
        Vec3d P = geom.edge_nearpt(edge, Vec3d(0, 0, 0)).unit();
        move_line_to_point(verts[edge[0]], verts[edge[1]], P);
      }

      // re-center for drift
      if (centering == 'e')
        geom.transform(Trans3d::translate(-edge_nearpoints_centroid(geom, Vec3d(0, 0, 0))));
      else
      if (centering == 'v')
        geom.transform(Trans3d::translate(-centroid(geom.verts())));
    }

    //for (unsigned int f = 0; f < geom.faces().size(); f++) {
    for (unsigned int ff = cnt; ff < geom.faces().size() + cnt; ff++) {
      int f = ff % geom.faces().size();
/*
      // give polygon its own geom. face index needs to reside in vector
      vector<int> face_idxs(1);
      face_idxs[0] = f;
      Geometry polygon = faces_to_geom(geom, face_idxs);
      plane_face(polygon);

      // map vertices back into original geom
      // the numerical order of vertex list in polygon geom is preserved
      vector<int> v_idx;
      for (int v : geom.faces(f))
        v_idx.push_back(v);
      sort(v_idx.begin(), v_idx.end());
      int j = 0;
      for (int v : v_idx)
        verts[v] = polygon.verts(j++);
*/
      // RK - this does formulaically what the above does by brute force
      Vec3d face_normal = face_normal_by_type(geom, f, normal_type).unit();
      Vec3d face_centroid = geom.face_cent(f);
      // make sure face_normal points outward
      if (vdot(face_normal, face_centroid) < 0)
        face_normal *= -1.0;
      // place a planar vertex over or under verts[v]
      // adds or subtracts it to get to the planar verts[v]
      for (int v : geom.faces(f))
        verts[v] += vdot(face_normal, face_centroid - verts[v]) *
                    face_normal;
    }

    // len2() for difference value to minimize internal sqrt() calls
    max_diff2 = 0;
    for (unsigned int i = 0; i < verts.size(); i++) {
      double diff2 = (verts[i] - verts_last[i]).len2();
      if (diff2 > max_diff2)
        max_diff2 = diff2;
    }

    // increment count here for reporting
    cnt++;

    if ((rep_count > 0) && (cnt%rep_count == 0))
      fprintf(stderr, "%-15d max_diff=%.17g\n", cnt, sqrt(max_diff2));

    if (sqrt(max_diff2) < eps) {
      completed = true;
      break;
    }

    // if minimum and maximum radius are differing, the polyhedron is crumpling
    if (radius_range_percent && canonical_radius_range_test(geom, radius_range_percent)) {
      fprintf(stderr, "\nbreaking out: radius range detected. try increasing -d\n");
      break;
    }
  }

  if (rep_count > -1) {
    fprintf(stderr, "\n%-15d final max_diff=%.17g\n", cnt, sqrt(max_diff2));
    fprintf(stderr, "\n");
  }

  return completed;
}

void planarity_info(Geometry &geom)
{
  GeometryInfo rep(geom);

  double max_nonplanar = 0;
  double sum_nonplanar = 0;
  int sz = geom.faces().size();
  for (int i = 0; i < sz; i++) {
    double nonplanar = rep.get_f_max_nonplanars()[i];
    sum_nonplanar += nonplanar;
    if (nonplanar > max_nonplanar)
      max_nonplanar = nonplanar;
  }
  fprintf(stderr, "maximum_nonplanarity = %.17g\n", max_nonplanar);
  fprintf(stderr, "average_nonplanarity = %.17g\n", sum_nonplanar/sz);
//  fprintf(stderr, "isoperimetric quotient = %.17g\n", rep.isoperimetric_quotient());
  fprintf(stderr, "\n");
}

void midradius_info(Geometry &geom, const bool completed)
{
  double min = 0;
  double max = 0;
  Vec3d center;
  double radius = edge_nearpoints_radius(geom, min, max, center);
  fprintf(stderr,"midradius = %.17g (range: %.15g to %.15g)\n",radius, min, max);
  fprintf(stderr,"midcenter is the origin\n");
  fprintf(stderr,"near point centroid = (%.17g,%.17g,%.17g)\n",center[0], center[1], center[2]);
  double epsilon_local = 1e-8;
  if (!completed)
    fprintf(stderr,"Warning: the calculation did not complete\n"); 
  else
  if (double_ne(center[0], 0.0, epsilon_local) ||
      double_ne(center[1], 0.0, epsilon_local) ||
      double_ne(center[2], 0.0, epsilon_local))
    fprintf(stderr,"Warning: the result is NOT canonical. edge tangency only\n");
  else
    fprintf(stderr,"the result is canonical\n");
  fprintf(stderr,"\n");
}

void generate_points(const Geometry &base, const Geometry &dual, vector<Vec3d> &ips,
                    vector<Vec3d> &base_nearpts, vector<Vec3d> &dual_nearpts,
                    const cn_opts &opts)
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
        Vec3d intersection_point = segments_intersection(b0, b1, d0, d1, epsilon_local);
        if (intersection_point.is_set()) {
          ips.push_back(intersection_point);
          break;
        }
      }
    }

    if (ips.size() != base_edges.size()) {
      if (opts.canonical_method != 'x')
        fprintf(stderr,"Warning: only %d out of %d intersection points found. try more precision\n",
                (int)ips.size(), (int)base_edges.size());
      else
        fprintf(stderr,"Warning: only canonical models have intersection points\n");
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
Vec3d face_edge_nearpoints_centroid(const Geometry &geom, double &radius, const int face_no)
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
*/

Geometry unit_circle(int polygon_size, const Color &incircle_color, bool filled)
{
  Geometry incircle;
  double arc = deg2rad(360.0 / (double)polygon_size);

  double angle = 0.0;
  for (int i = 0; i < polygon_size; i++) {
    incircle.add_vert(Vec3d(cos(angle), sin(angle), 0.0), Color::invisible);
    incircle.add_edge(make_edge(i, (i + 1) % polygon_size), (filled ? Color::invisible : incircle_color));
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

    Vec3d P = geom.edge_nearpt(make_edge(v1,v2), Vec3d(0, 0, 0));
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

Geometry incircles(const Geometry &geom, const Color &incircle_color, bool filled, double offset)
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

void construct_model(Geometry &base, const cn_opts &opts) {
  // RK - set radius to 1 for get_dual call, if necessary
  double radius = edge_nearpoints_radius(base);
  if (double_ne(radius, 1.0, opts.epsilon))
    unitize_nearpoints_radius(base);

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
    base_incircles = incircles(base, opts.base_nearpts_col, filled, opts.offset);
  }

  // dual incircles
  Geometry dual_incircles;
  if (opts.output_parts.find_first_of("tT") != string::npos) {
    bool filled = (opts.output_parts.find("t") != string::npos) ? true : false;
    dual_incircles = incircles(dual, opts.dual_nearpts_col, filled, opts.offset);
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
    ColorValuesToRangeHsva valmap(msg_str("A%g", (double)opts.face_opacity/255), Color(255, 255, 255));
    valmap.apply(base, FACES);

    for (const auto & kp : base.colors(FACES).get_properties()) {
      if (kp.second.is_index()) {
        opts.warning("map indexes cannot be made transparent",'T');
        break;
      }
    }
  }

  // add unit sphere on origin
  if (opts.output_parts.find_first_of("uU") != string::npos) {
    Geometry sgeom;
    sgeom.read_resource("geo_4_4");
    sgeom.transform(Trans3d::translate(-centroid(sgeom.verts())));
    unitize_vertex_radius(sgeom);

    double min = 0;
    double max = 0;
    Vec3d center;
    edge_nearpoints_radius(base, min, max, center);
    if (opts.output_parts.find("u") != string::npos)
      sgeom.transform(Trans3d::scale(min));
    else
      sgeom.transform(Trans3d::scale(max));

    Coloring(&sgeom).vef_one_col(Color::invisible, Color::invisible, opts.sphere_col);
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
    opts.warning(msg_str("possible coincident vertices found in %s", s.c_str()));

  if (geom.faces().size() != geom_merged.faces().size())
    opts.warning(msg_str("possible coincident faces found in %s", s.c_str()));

  GeometryInfo info(geom);
  vector<double> areas = info.get_f_areas();
  for (unsigned int i = 0; i < areas.size(); i++) {
    //fprintf(stderr, "%.17lf\n", areas[i]);
    if (areas[i] < epsilon_local2) {
      opts.warning(msg_str("possible faces of near 0 area found in %s (face %d)", s.c_str(), i));
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

  fprintf(stderr,"\n");
  fprintf(stderr,"starting radius: ");
  if (opts.initial_radius == 'e') {
    fprintf(stderr, "average edge near points\n");
    unitize_nearpoints_radius(geom);
  }
  else
  if (opts.centering == 'v') {
    fprintf(stderr, "average vertex\n");
    unitize_vertex_radius(geom);
  }
  else
  if (opts.centering == 'x')
    fprintf(stderr, "radius not changed\n");

  fprintf(stderr,"centering: ");
  if (opts.centering == 'e') {
    fprintf(stderr, "edge near points centroid to origin\n");
    geom.transform(Trans3d::translate(-edge_nearpoints_centroid(geom, Vec3d(0, 0, 0))));
  }
  else
  if (opts.centering == 'v') {
    fprintf(stderr, "vertex centroid to origin\n");
    geom.transform(Trans3d::translate(-centroid(geom.verts())));
  }
  else
  if (opts.centering == 'x')
    fprintf(stderr, "model not moved\n");

  fprintf(stderr,"normals: ");
  if (opts.normal_type == 'n')
    fprintf(stderr,"Newell's method\n");
  else
  if (opts.normal_type == 't')
    fprintf(stderr,"Triangles method\n");
  else
  if (opts.normal_type == 'q')
    fprintf(stderr,"Quads method\n");

  bool completed = false;
  if (opts.planarize_method) {
    string planarize_str;
    if (opts.planarize_method == 'p')
      planarize_str = "face centroids magnitude squared";
    else
    if (opts.planarize_method == 'q')
      planarize_str = "face centroids magnitude";
    else
    if (opts.planarize_method == 'f')
      planarize_str = "face centroids";
    else
    if (opts.planarize_method == 'm')
      planarize_str = "mathematica";
    else
    if (opts.planarize_method == 'a')
      planarize_str = "sand and fill";
    else
    if (opts.planarize_method == 'u')
      planarize_str = "minmax -a u";
    fprintf(stderr, "planarize: %s method\n",planarize_str.c_str());

    if (opts.planarize_method == 'm') {
      bool planarize_only = true;
      completed = canonicalize_mm(geom, opts.mm_edge_factor / 100, opts.mm_plane_factor / 100,
                                 opts.num_iters_planar, opts.radius_range_percent / 100, opts.rep_count,
                                 opts.alternate_algorithm, planarize_only, opts.normal_type, opts.epsilon);
    }
    else
    if (opts.planarize_method == 'a') {
      bool planarize_only = true;
      completed = canonicalize_unit(geom, opts.num_iters_planar, opts.radius_range_percent / 100,
                                    opts.rep_count, opts.centering, opts.normal_type, planarize_only, opts.epsilon);
    }
    // case u
    else
    if (opts.planarize_method == 'u')
      completed = minmax_unit_planar(geom, opts.num_iters_planar, opts.radius_range_percent / 100, opts.rep_count, opts.normal_type, opts.epsilon);
    // cases p, q, f
    else
      completed = canonicalize_bd(geom, opts.num_iters_planar, opts.planarize_method,
                                 opts.radius_range_percent / 100, opts.rep_count, opts.centering, opts.normal_type, opts.epsilon);

    // RK - report planarity
    planarity_info(geom);
  }

  if (opts.canonical_method && opts.canonical_method != 'x') {
    string canonicalize_str;
    if (opts.canonical_method == 'm')
      canonicalize_str = "mathematica";
    else
    if (opts.canonical_method == 'b')
      canonicalize_str = "base/dual";
    else
    if (opts.canonical_method == 'a')
      canonicalize_str = "moving edge";
    fprintf(stderr, "canonicalize: %s method\n",canonicalize_str.c_str());
    if (opts.canonical_method == 'm') {
      bool planarize_only = false;
      completed = canonicalize_mm(geom, opts.mm_edge_factor / 100, opts.mm_plane_factor / 100,
                                 opts.num_iters_canonical, opts.radius_range_percent / 100, opts.rep_count,
                                 opts.alternate_algorithm, planarize_only, opts.normal_type, opts.epsilon);
    }
    else
    if (opts.canonical_method == 'b') {
      completed = canonicalize_bd(geom, opts.num_iters_canonical, opts.canonical_method,
                                 opts.radius_range_percent / 100, opts.rep_count, opts.centering, opts.normal_type, opts.epsilon);
    }
    else
    if (opts.canonical_method == 'a') {
      bool planarize_only = false;
      completed = canonicalize_unit(geom, opts.num_iters_canonical, opts.radius_range_percent / 100,
                                    opts.rep_count, opts.centering, opts.normal_type, planarize_only, opts.epsilon);
    }

    // RK - report planarity
    planarity_info(geom);

    // RK - print midradius info
    midradius_info(geom, completed);
  }

  // RK - add coincidence checking the model
  check_coincidence(geom, opts);

  // RK - parts to output
  construct_model(geom, opts);

  opts.write_or_error(geom, opts.ofile);

  return 0;
}
