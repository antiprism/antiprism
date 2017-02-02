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
  string stderr;

  char edge_distribution;
  char planarize_method;
  int num_iters_planar;
  char canonical_method;
  int num_iters_canonical;
  double mm_edge_factor;
  double mm_plane_factor;
  int rep_count;
  double radius_range_percent;

  double epsilon;

  cn_opts()
      : ProgramOpts("canonical"), edge_distribution('\0'),
        planarize_method('\0'), num_iters_planar(10000), canonical_method('m'),
        num_iters_canonical(-1), mm_edge_factor(50), mm_plane_factor(20),
        rep_count(1000), radius_range_percent(80), epsilon(0)
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
"(BETA)\n"
"Read a polyhedron from a file in OFF format. Canonicalize or planarize it.\n"
"Uses algorithms by George W. Hart, http://www.georgehart.com/\n"
"http://www.georgehart.com/virtual-polyhedra/conway_notation.html\n"
"http://www.georgehart.com/virtual-polyhedra/canonical.html\n"
"If input_file is not given the program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"  -e <opt>  edge distribution\n"
"               s - project vertices onto a sphere\n"
"               a - (another method to be implimented?)\n"
"               x - none\n"
"  -p <opt>  planarization (done before canoncalization)\n"
"               p - face centroids (magnitude squared)\n"
"               q - face centroids (magnitude)\n"
"               f - face centroids\n"
"               m - mathematica planarize\n"
"               x - none\n"
"  -i <itrs> maximum number of planarize iterations (default: 10000)\n"
"  -c <opt>  canonicalization\n"
"               m - mathematica version (default)\n"
"               n - conway notation base/dual version\n"
"               x - none\n"
"  -n <itrs> maximum number of iterations (default: no limit)\n"
"  -d <val>  radius test. precent difference between minumum and maximum radius\n"
"               checks if polyhedron is collapsing. 0 for no test (default 10)\n"
"  -z <n>    status reporting every n lines. -1 for no status. (default 1000)\n"
"  -l <lim>  minimum distance change to terminate, as negative exponent\n"
"               (default: %d giving %.0e)\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"Mathematica Canonicalize Options (-M m and -M l)\n"
"  -E <perc> percentage to scale the edge tangency error (default: 50)\n" 
"  -P <perc> percentage to scale the face planarity error (default: 20)\n" 
"\n"
"\n",prog_name(), help_ver_text, int(-log(::epsilon)/log(10) + 0.5), ::epsilon);
}
// clang-format on 

void cn_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  string arg_id;

  int sig_compare = INT_MAX;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":he:p:i:c:E:P:n:d:z:l:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'e':
      if (strlen(optarg) == 1 && strchr("sax", int(*optarg)))
        edge_distribution = *optarg;
      else
        error("edge_distribution method type must be s, a, x", c);
      break;

    case 'p':
      if (strlen(optarg) == 1 && strchr("pqfmx", int(*optarg)))
        planarize_method = *optarg;
      else
        error("planarize method type must be p, q, f, m, x", c);
      break;

    case 'i':
      print_status_or_exit(read_int(optarg, &num_iters_planar), c);
      if (num_iters_planar <= 0)
        error(
            "number of iterations for preplanarization must be greater than 0",
            c);
      break;

    case 'c':
      if (strlen(optarg) == 1 && strchr("mnx", int(*optarg)))
        canonical_method = *optarg;
      else
        error("canonical method type must be m, n, x", c);
      break;

    case 'n':
      print_status_or_exit(read_int(optarg, &num_iters_canonical), c);
      if (num_iters_canonical < 0)
        error("number of iterations must be 0 or greater", c);
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
      stderr = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (argc - optind > 1)
    error("too many arguments");

  if (argc - optind == 1)
    ifile = argv[optind];

  epsilon = (sig_compare != INT_MAX) ? pow(10, -sig_compare) : ::epsilon;
}

// return true if maximum vertex radius is radius_range_percent (0.0 to ...)
// greater than minimum vertex radius
bool radius_range_test(const Geometry &geom, const double &radius_range_percent)
{
  GeometryInfo rep(geom);
  rep.set_center(geom.centroid());

  double min = rep.vert_dist_lims().min;
  double max = rep.vert_dist_lims().max;
    
  // min and max should always be positive, max should always be larger
  return (((max-min)/((max+min)/2.0)) > radius_range_percent) ? true : false;
}

// reciprocalN() is from the Hart's Conway Notation web page
// make array of vertices reciprocal to given planes (face normals)
vector<Vec3d> reciprocalN2(const Geometry &geom)
{
  const vector<vector<int>> &faces = geom.faces();
  const vector<Vec3d> &verts = geom.verts();

  vector<Vec3d> normals;
  for (const auto &face : faces) {
    Vec3d centroid(0, 0, 0);
    Vec3d normal(0, 0, 0);
    double avgEdgeDist = 0;

    int v1 = face.at(face.size() - 2);
    int v2 = face.at(face.size() - 1);
    for (int v3 : face) {
      centroid += verts[v3];
      // orthogonal() was from the Hart's Conway Notation web page. replacement
      // normal += orthogonal(verts[v1], verts[v2], verts[v3]);
      normal += vcross(verts[v3] - verts[v2], verts[v2] - verts[v1]);
      // tangentPoint() was from Hart's Conway Notation web page. replacement
      // avgEdgeDist += tangentPoint(verts[v1], verts[v2]).len();
      Vec3d d = verts[v2] - verts[v1];
      // prevent division by zero
      // avgEdgeDist += (verts[v1] - ((vdot(d,verts[v1])/d.len2()) * d)).len();
      double vdt;
      if (d[0] == 0 && d[1] == 0 && d[2] == 0)
        vdt = 0;
      else
        vdt = vdot(d, verts[v1]) / d.len2();
      avgEdgeDist += (verts[v1] - (vdt * d)).len(); // tangentPoint without call
      v1 = v2;
      v2 = v3;
    }
    centroid *= 1.0 / face.size();
    normal.to_unit();
    avgEdgeDist /= face.size();

    // reciprocal call replace below:
    // prevent division by zero
    // Vec3d ans = reciprocal(normal * vdot(centroid,normal));
    Vec3d v = normal * vdot(centroid, normal);
    Vec3d ans;
    if (v[0] == 0 && v[1] == 0 && v[2] == 0)
      ans = v;
    else {
      ans = v * 1.0 / v.len2();
      ans *= (1 + avgEdgeDist) / 2;
    }
    normals.push_back(ans);
  }

  return normals;
}

// reciprocate on face centers dividing by magnitude squared
vector<Vec3d> reciprocalC2(const Geometry &geom)
{
  vector<Vec3d> centers;
  geom.face_cents(centers);
  for (auto &center : centers)
    center /= center.len2();
  return centers;
}

// reciprocate on face centers dividing by magnitude
vector<Vec3d> reciprocalC_len2(const Geometry &geom)
{
  vector<Vec3d> centers;
  geom.face_cents(centers);
  for (auto &center : centers)
    center /= center.len();
  return centers;
}

// Addition to algorithm by Adrian Rossiter
// Finds the correct centroid for the canonical
Vec3d edge_nearpoints_centroid2(const Geometry &geom, const Vec3d &cent)
{
  vector<vector<int>> edges;
  geom.get_impl_edges(edges);
  int e_sz = edges.size();
  Vec3d e_cent(0, 0, 0);
  for (int e = 0; e < e_sz; ++e)
    e_cent += geom.edge_nearpt(edges[e], cent);
  return e_cent / double(e_sz);
}

// RK - average radius rather than maximum has more reliability than max
void unitize_radius(Geometry &geom)
{
  GeometryInfo info(geom);
  info.set_center(geom.centroid());
  //geom.transform(Trans3d::scale(1 / info.vert_dist_lims().max));
  double avg = info.vert_dist_lims().sum / info.num_verts();
  geom.transform(Trans3d::scale(1 / avg));
}

void centroid_to_origin2(Geometry &geom)
{
  geom.transform(Trans3d::transl(-centroid(geom.verts())));
}

// Implementation of George Hart's planarization and canonicalization algorithms
// http://www.georgehart.com/virtual-polyhedra/conway_notation.html
void canonicalize_cn2(Geometry &base, const int &num_iters, const char &canonical_method, const double &radius_range_percent,
                     const int &rep_count, const double &eps)
{
  Geometry dual;
  get_dual(&dual, base, 0);
  dual.clear_cols();
  const vector<Vec3d> &base_verts = base.verts();
  const vector<Vec3d> &dual_verts = dual.verts();

  double max_diff2 = 0;
  unsigned int cnt;
  for (cnt = 0; cnt < (unsigned int)num_iters;) {
    vector<Vec3d> base_verts_last = base_verts;

    switch (canonical_method) {
    // base/dual canonicalize method
    case 'n': {
      dual.raw_verts() = reciprocalN2(base);
      base.raw_verts() = reciprocalN2(dual);
      Vec3d e_cent = edge_nearpoints_centroid2(base, Vec3d(0, 0, 0));
      base.transform(Trans3d::transl(-0.1 * e_cent));
      break;
    }

    // adjust vertices with side effect of planarization. len2() version
    case 'p':
      // move centroid to origin for balance
      dual.raw_verts() = reciprocalC2(base);
      base.transform(Trans3d::transl(-centroid(dual_verts)));
      base.raw_verts() = reciprocalC2(dual);
      base.transform(Trans3d::transl(-centroid(base_verts)));
      break;

    // adjust vertices with side effect of planarization. len() version
    case 'q':
      // move centroid to origin for balance
      dual.raw_verts() = reciprocalC_len2(base);
      base.transform(Trans3d::transl(-centroid(dual_verts)));
      base.raw_verts() = reciprocalC_len2(dual);
      base.transform(Trans3d::transl(-centroid(base_verts)));
      break;

    // adjust vertices with side effect of planarization. face centroids version
    case 'f':
      base.face_cents(dual.raw_verts());
      dual.face_cents(base.raw_verts());
      break;
    }

    max_diff2 = 0;
    for (unsigned int i = 0; i < base_verts.size(); i++) {
      double diff2 = (base_verts[i] - base_verts_last[i]).len2();
      if (diff2 > max_diff2)
        max_diff2 = diff2;
    }

    // increment count here for reporting
    cnt++;

    if ((rep_count > 0) && (cnt%rep_count == 0))
      fprintf(stderr, "%-15d max_diff=%.17g\n", cnt, sqrt(max_diff2));

    if (sqrt(max_diff2) < eps)
      break;

    // if minimum and maximum radius are differing, the polyhedron is crumpling
    if (radius_range_percent && radius_range_test(base, radius_range_percent)) {
      fprintf(stderr, "\nbreaking out: radius range detected. try increasing -d\n");
      break;
    }
  }

  if (rep_count > -1) {
    fprintf(stderr, "\n%-15d final max_diff=%.17g\n", cnt, sqrt(max_diff2));
    fprintf(stderr, "\n");
  }
}

// Implementation of George Hart's canonicalization algorithm
// http://library.wolfram.com/infocenter/Articles/2012/
void canonicalize_mm2(Geometry &geom, const double &edge_factor, const double &plane_factor,
                     const int &num_iters, const double &radius_range_percent, const int &rep_count,
                     const bool &planar_only, const double &eps)
{
  // RK - functions better when the input polyhedron has a radius near 1
  unitize_radius(geom);

  const vector<Vec3d> &verts = geom.verts();
  const vector<vector<int>> &faces = geom.faces();
  vector<vector<int>> edges;
  geom.get_impl_edges(edges);

  double max_diff2 = 0;
  unsigned int cnt;
  for (cnt = 0; cnt < (unsigned int)num_iters;) {
    vector<Vec3d> verts_last = verts;

    // RK - the model will possibly become non-convex early in the loops.
    // if it contorts too badly, the model will implode. Having the model
    // at a radius of near 1 minimizes this problem
    if (!planar_only) {
      vector<Vec3d> near_pts;
      for (auto &edge : edges) {
        Vec3d P = geom.edge_nearpt(edge, Vec3d(0, 0, 0));
        near_pts.push_back(P);
// RK - these 4 lines cause the near points to be applied in a 2nd loop
// but this causes more problems. The solution is to have the input model
// have an input radius of near 1
//      }
//      int p_cnt = 0;
//      for (auto &edge : edges) {
//        Vec3d P = near_pts[p_cnt++];
        Vec3d offset = edge_factor * (P.len() - 1) * P;
        geom.verts(edge[0]) -= offset;
        geom.verts(edge[1]) -= offset;
      }

      Vec3d cent_near_pts = centroid(near_pts);
      for (unsigned int i = 0; i < verts.size(); i++)
        geom.verts(i) -= cent_near_pts;
    }

    // Make a copy of verts. zero out.
    vector<Vec3d> vs = verts;
    for (auto &v : vs)
      v = Vec3d(0, 0, 0);

    // progressively advances starting face each iteration
    for (unsigned int ff = cnt; ff < faces.size() + cnt; ff++) {
      int f = ff % faces.size();
      if (faces[f].size() == 3)
        continue;
      Vec3d norm = geom.face_norm(f).unit();
      Vec3d f_cent = geom.face_cent(f);
      if (vdot(norm, f_cent) < 0)
        norm *= -1.0;
      for (int v : faces[f])
        vs[v] += vdot(plane_factor * norm, f_cent - verts[v]) * norm;
    }

    // adjust vertices post-loop
    for (unsigned int i = 0; i < vs.size(); i++)
      geom.verts(i) += vs[i];

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

    if (sqrt(max_diff2) < eps)
      break;

    // if minimum and maximum radius are differing, the polyhedron is crumpling
    if (radius_range_percent && radius_range_test(geom, radius_range_percent)) {
      fprintf(stderr, "\nbreaking out: radius range detected. try increasing -d\n");
      break;
    }
  }

  if (rep_count > -1) {
    fprintf(stderr, "\n%-15d final max_diff=%12.10g\n", cnt, sqrt(max_diff2));
    fprintf(stderr, "\n");
  }
}

void planar_info(Geometry &geom)
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

// RK - find nearpoints radius, sets range minimum and maximum
double edge_nearpoints_radius(const Geometry &geom, double &min, double &max)
{
  min = DBL_MAX;
  max = DBL_MIN;

  vector<vector<int>> edges;
  geom.get_impl_edges(edges);
  int e_sz = edges.size();
  double nearpt_radius = 0;
  for (int e = 0; e < e_sz; ++e) {
    double len = geom.edge_nearpt(edges[e], geom.centroid()).len2();
    nearpt_radius += len;
    if (len < min)
      min = len;
    if (len > max)
      max = len;
  }
  return nearpt_radius / double(e_sz);
}

void midradius_info(Geometry &geom)
{
  double min = 0;
  double max = 0;
  double radius = edge_nearpoints_radius(geom, min, max);
  fprintf(stderr,"midradius = %.17g (range: %.15g to %.15g)\n",radius, min, max);
}

/*
void add_tangent_points(const Geometry &geom)
{
  Geometry dual;
  get_dual(&dual, geom, 1, Vec3d(0, 0, 0));

  int sz = geom.faces().size();
  for (int i = 0; i < sz; i++) {
  }
}
*/

int main(int argc, char *argv[])
{
  cn_opts opts;
  opts.process_command_line(argc, argv);

/*
  // test: all johnson solids
  for (int i = 1; i < 93; i++) {
    string str = "j" + std::to_string(i);
    fprintf(stderr,"%s\n",str.c_str());
    Geometry geom;
    opts.read_or_error(geom, str);

    bool planarize_only = false;
    canonicalize_mm2(geom, opts.mm_edge_factor / 100, opts.mm_plane_factor / 100,
                     opts.num_iters_canonical, opts.radius_range_percent / 100, opts.rep_count,
                     planarize_only, opts.epsilon);

    str += ".off";
    opts.write_or_error(geom, str);
  }
  exit(0);
*/

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);

  // RK - the functions expect the model to be centered on the origin
  centroid_to_origin2(geom);

  fprintf(stderr,"\n");
  if (opts.edge_distribution && opts.edge_distribution != 'x') {
    fprintf(stderr, "edge distribution: (project onto sphere)\n");
    if (opts.edge_distribution == 's')
      project_onto_sphere(&geom);
  }

  if (opts.planarize_method && opts.planarize_method != 'x') {
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
    fprintf(stderr, "planarize: (%s method)\n",planarize_str.c_str());

    if (opts.planarize_method == 'm') {
      bool planarize_only = true;
      canonicalize_mm2(geom, opts.mm_edge_factor / 100, opts.mm_plane_factor / 100,
                     opts.num_iters_planar, opts.radius_range_percent / 100, opts.rep_count,
                     planarize_only, opts.epsilon);
    }
    else {
      canonicalize_cn2(geom, opts.num_iters_planar, opts.planarize_method,
                      opts.radius_range_percent / 100, opts.rep_count, opts.epsilon);
    }

    // RK - report planarity
    planar_info(geom);
  }

  if (opts.canonical_method && opts.canonical_method != 'x') {
    fprintf(stderr, "canonicalize: (%s method)\n",((opts.canonical_method == 'm') ? "mathematica" : "base/dual"));
    if (opts.canonical_method == 'm') {
      bool planarize_only = false;
      canonicalize_mm2(geom, opts.mm_edge_factor / 100, opts.mm_plane_factor / 100,
                     opts.num_iters_canonical, opts.radius_range_percent / 100, opts.rep_count,
                     planarize_only, opts.epsilon);
    }
    else
      canonicalize_cn2(geom, opts.num_iters_canonical, opts.canonical_method,
                      opts.radius_range_percent / 100, opts.rep_count, opts.epsilon);

    // RK - report planarity
    planar_info(geom);

    // RK - print midradius before we move it
    midradius_info(geom);

    // RK - questioning to do this?
    //centroid_to_origin2(geom);
  }

  opts.write_or_error(geom, opts.stderr);

  return 0;
}
