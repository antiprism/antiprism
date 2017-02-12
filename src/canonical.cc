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

  char centering;
  char initial_radius;
  char edge_distribution;
  char planarize_method;
  int num_iters_planar;
  char canonical_method;
  int num_iters_canonical;
  double mm_edge_factor;
  double mm_plane_factor;
  bool mm_alternate_loop;
  int rep_count;
  double radius_range_percent;
  string output_parts;
  int face_opacity;

  double epsilon;

  Color ipoints_col;
  Color base_nearpts_col;
  Color dual_nearpts_col;
  Color base_edge_col;
  Color dual_edge_col;
  Color sphere_col;

  cn_opts()
      : ProgramOpts("canonical"), centering('v'), initial_radius('v'), edge_distribution('\0'),
        planarize_method('\0'), num_iters_planar(-1), canonical_method('m'),
        num_iters_canonical(-1), mm_edge_factor(50), mm_plane_factor(20),
        mm_alternate_loop(false), rep_count(1000), radius_range_percent(80),
        output_parts("b"), face_opacity(-1), epsilon(0), ipoints_col(Color(255, 255, 0)),
        base_nearpts_col(Color(255, 0, 0)),
        dual_nearpts_col(Color(0.0, 0.39216, 0.0)), base_edge_col(Color()),
        dual_edge_col(Color()),sphere_col(Color(255, 255, 255))
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
"  -C <opt>  initial centering\n"
"               v - vertex centroid (default)\n"
"               n - edge near points centroid\n"
"               x - not moved\n"
"  -r <opt>  initial radius\n"
"               v - average vertex radius = 1 (default)\n"
"               n - average edge near points radius = 1\n"
"               x - not changed\n"
"  -e <opt>  edge distribution\n"
"               s - project vertices onto a sphere\n"
"               a - (another method to be implimented?)\n"
"  -p <opt>  planarization (done before canoncalization. default: none)\n"
"               p - face centroids (magnitude squared)\n"
"               q - face centroids (magnitude)\n"
"               f - face centroids\n"
"               m - mathematica planarize\n"
"  -i <itrs> maximum number of planarize iterations (default: no limit)\n"
"  -c <opt>  canonicalization\n"
"               m - mathematica version (default)\n"
"               b - base/dual version\n"
"               x - none (default, if -p is set)\n"
"  -n <itrs> maximum number of canonical iterations (default: no limit)\n"
"  -O <args> output b - base, d - dual, i - intersection points (default: b)\n"
"               n - base edge near points, m - dual edge near points\n"
"               p - base near points centeroid, q - dual near points centroid\n"
"               u - unit sphere centered on the origin, o - origin point\n"
"  -d <perc> radius test. precent difference between minumum and maximum radius\n"
"               checks if polyhedron is collapsing. 0 for no test (default: 10)\n"
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
"  -N <col>  base edge near points and/or centroid color (default: red)\n"
"  -M <col>  dual edge near points and/or centroid color (default: darkgreen)\n"
"  -B <col>  base edge color (default: none)\n"
"  -D <col>  dual edge color (default: none)\n"
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

  while ((c = getopt(argc, argv, ":hC:r:e:p:i:c:n:O:E:P:Ad:z:I:N:M:B:D:U:T:l:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'C':
      if (strlen(optarg) == 1 && strchr("vnx", int(*optarg)))
        centering = *optarg;
      else
        error("centering method type must be v, n, x", c);
      break;

    case 'r':
      if (strlen(optarg) == 1 && strchr("vnx", int(*optarg)))
        initial_radius = *optarg;
      else
        error("starting radius type must be v, n, x", c);
      break;

    case 'e':
      if (strlen(optarg) == 1 && strchr("sa", int(*optarg)))
        edge_distribution = *optarg;
      else
        error("edge_distribution method type must be s, a", c);
      break;

    case 'p':
      p_set = true;
      if (strlen(optarg) == 1 && strchr("pqfm", int(*optarg)))
        planarize_method = *optarg;
      else
        error("planarize method type must be p, q, f, m", c);
      break;

    case 'i':
      print_status_or_exit(read_int(optarg, &num_iters_planar), c);
      if (num_iters_planar <= 0)
        error(
            "number of iterations for preplanarization must be greater than 0",
            c);
      break;

    case 'c':
      c_set = true;
      if (strlen(optarg) == 1 && strchr("mbx", int(*optarg)))
        canonical_method = *optarg;
      else
        error("canonical method type must be m, b, x", c);
      break;

    case 'n':
      print_status_or_exit(read_int(optarg, &num_iters_canonical), c);
      if (num_iters_canonical < 0)
        error("number of iterations must be 0 or greater", c);
      break;

    case 'O':
      if (strspn(optarg, "bdinmopqu") != strlen(optarg))
        error(msg_str("output parts are '%s' must be any or all from "
                      "b, d, i, n, m, o, p, q, u", optarg), c);
      output_parts = optarg;
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
      mm_alternate_loop = true;
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
      stderr = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  // if planarizing only do not canonicalize
  if (p_set && !c_set)
    canonical_method = 'x';

  if (mm_alternate_loop && canonical_method != 'm')
    warning("alternate form only has effect in mathematica canonicalization", 'A');  

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
// RK - save of verbatim port code
/*
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
*/

/* RK - tangentPoint() can be replaced by edge_nearpt()
      // tangentPoint() was from Hart's Conway Notation web page
      // point where line v1...v2 tangent to an origin sphere
      // avgEdgeDist += tangentPoint(verts[v1], verts[v2]).len();
      // avgEdgeDist += (verts[v1] - ((vdot(d,verts[v1])/d.len2()) * d)).len();
      Vec3d d = verts[v2] - verts[v1];
      double vdt = 0;
      // prevent division by zero
      if (d[0] != 0 || d[1] != 0 || d[2] != 0)
        vdt = vdot(d, verts[v1]) / d.len2();
      avgEdgeDist += (verts[v1] - (vdt * d)).len();
*/

// reciprocalN() is from the Hart's Conway Notation web page
// make array of vertices reciprocal to given planes (face normals)
vector<Vec3d> reciprocalN2(const Geometry &geom)
{
  const vector<vector<int>> &faces = geom.faces();
  const vector<Vec3d> &verts = geom.verts();

  vector<Vec3d> normals;

  for (unsigned int i = 0; i < faces.size(); i++) {
    // calculate face normal in antiprism
    Vec3d face_normal = face_norm(verts, faces[i]).unit();
    Vec3d face_centroid = anti::centroid(verts, faces[i]);

    vector<int> face = faces[i];
    unsigned int sz = face.size();
    double avgEdgeDist = 0;
    for (unsigned int j = 0; j < sz; j++) {
      int v1 = face[j];
      int v2 = face[(j+1)%sz];

      avgEdgeDist += geom.edge_nearpt(make_edge(v1,v2), Vec3d(0, 0, 0)).len();
    }
    avgEdgeDist /= sz;

    // reciprocal call replace below:
    // Vec3d ans = reciprocal(normal * vdot(centroid,normal));
    Vec3d v = face_normal * vdot(face_centroid, face_normal);
    Vec3d ans = v;
    // prevent division by zero
    if (v[0] != 0 || v[1] != 0 || v[2] != 0)
      ans = v * 1.0 / v.len2();
    // edge correction
    ans *= (1 + avgEdgeDist) / 2;

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
Vec3d edge_nearpoints_centroid(const Geometry &geom, const Vec3d &cent)
{
  vector<vector<int>> edges;
  geom.get_impl_edges(edges);
  int e_sz = edges.size();
  Vec3d e_cent(0, 0, 0);
  for (int e = 0; e < e_sz; e++)
    e_cent += geom.edge_nearpt(edges[e], cent);
  return e_cent / double(e_sz);
}


// Implementation of George Hart's planarization and canonicalization algorithms
// http://www.georgehart.com/virtual-polyhedra/conway_notation.html
bool canonicalize_bd2(Geometry &base, const int &num_iters, const char &canonical_method,
                     const double &radius_range_percent, const int &rep_count, 
                     const char &centering, const double &eps)
{
  bool completed = false;

  Geometry dual;
  // the dual's initial vertex locations are immediately overwritten
  get_dual(&dual, base, 1);
  dual.clear_cols();

  const vector<Vec3d> &base_verts = base.verts();
  const vector<Vec3d> &dual_verts = dual.verts();

  double max_diff2 = 0;
  unsigned int cnt;
  for (cnt = 0; cnt < (unsigned int)num_iters;) {
    vector<Vec3d> base_verts_last = base_verts;

    switch (canonical_method) {
    // base/dual canonicalize method
    case 'b': {
      dual.raw_verts() = reciprocalN2(base);
      base.raw_verts() = reciprocalN2(dual);
      if (centering != 'x') {
        Vec3d e_cent = edge_nearpoints_centroid(base, Vec3d(0, 0, 0));
        base.transform(Trans3d::transl(-0.1 * e_cent));
      }
      break;
    }
/*
    case 'b': {
      dual.raw_verts() = reciprocalN2(base);
      Vec3d d_cent = edge_nearpoints_centroid(dual, Vec3d(0, 0, 0));
      base.raw_verts() = reciprocalN2(dual);
      Vec3d e_cent = edge_nearpoints_centroid(base, Vec3d(0, 0, 0));
double diff = (e_cent - d_cent).len2();
if (diff < eps) {
fprintf(stderr,"%d: converged\n",cnt);
}
//fprintf(stderr,"true\n");
//fprintf(stderr,"%d: e_cent = (%.17g,%.17g,%.17g) len = %.17lf\n",cnt, e_cent[0], e_cent[1], e_cent[2], e_cent.len());
//if (e_cent.len() < eps)
//fprintf(stderr,"true\n");
      base.transform(Trans3d::transl(-0.1 * e_cent));
      break;
    }
*/

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

    if (sqrt(max_diff2) < eps) {
      completed = true;
      break;
    }

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

  return completed;
}

// Implementation of George Hart's canonicalization algorithm
// http://library.wolfram.com/infocenter/Articles/2012/
// RK - the model will possibly become non-convex early in the loops.
// if it contorts too badly, the model will implode. Having the model
// at a radius of near 1 minimizes this problem
bool canonicalize_mm2(Geometry &geom, const double &edge_factor, const double &plane_factor,
                     const int &num_iters, const double &radius_range_percent, const int &rep_count,
                     const bool &planar_only, const bool &alternate_loop, const double &eps)
{
  bool completed = false;

  const vector<Vec3d> &verts = geom.verts();
  const vector<vector<int>> &faces = geom.faces();
  vector<vector<int>> edges;
  geom.get_impl_edges(edges);

  double max_diff2 = 0;
  unsigned int cnt;
  for (cnt = 0; cnt < (unsigned int)num_iters;) {
    vector<Vec3d> verts_last = verts;

    if (!planar_only) {
      vector<Vec3d> near_pts;
      if (!alternate_loop) {
        for (auto &edge : edges) {
          Vec3d P = geom.edge_nearpt(edge, Vec3d(0, 0, 0));
          near_pts.push_back(P);
          Vec3d offset = edge_factor * (P.len() - 1) * P;
          geom.verts(edge[0]) -= offset;
          geom.verts(edge[1]) -= offset;
        }
      }
      // RK - alternate form causes the near points to be applied in a 2nd loop
      // most often not needed unless the model is off balance
      else {
        for (auto &edge : edges) {
          Vec3d P = geom.edge_nearpt(edge, Vec3d(0, 0, 0));
          near_pts.push_back(P);
        // RK - these 4 lines cause the near points to be applied in a 2nd loop
        }
        int p_cnt = 0;
        for (auto &edge : edges) {
          Vec3d P = near_pts[p_cnt++];
          Vec3d offset = edge_factor * (P.len() - 1) * P;
          geom.verts(edge[0]) -= offset;
          geom.verts(edge[1]) -= offset;
        }
      }
/*
      // RK - revolving loop. didn't solve the imbalance problem
      else {
        for (unsigned int ee = cnt; ee < edges.size() + cnt; ee++) {
          int e = ee % edges.size();
          Vec3d P = geom.edge_nearpt(edges[e], Vec3d(0, 0, 0));
          near_pts.push_back(P);
          Vec3d offset = edge_factor * (P.len() - 1) * P;
          geom.verts(edges[e][0]) -= offset;
          geom.verts(edges[e][1]) -= offset;
        }
      }
*/

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

    if (sqrt(max_diff2) < eps) {
      completed = true;
      break;
    }

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

// RK - find nearpoints radius, sets range minimum and maximum
double edge_nearpoints_radius(const Geometry &geom, double &min, double &max, Vec3d &center)
{
  min = DBL_MAX;
  max = DBL_MIN;

  vector<vector<int>> edges;
  geom.get_impl_edges(edges);

  vector<Vec3d> near_pts;

  double nearpt_radius = 0;
  int e_sz = edges.size();
  for (int e = 0; e < e_sz; e++) {
    Vec3d P = geom.edge_nearpt(edges[e], Vec3d(0, 0, 0));
    near_pts.push_back(P);

    double len = P.len2();
    nearpt_radius += len;
    if (len < min)
      min = len;
    if (len > max)
      max = len;
  }

  center = centroid(near_pts);

  return nearpt_radius / double(e_sz);
}

void midradius_info(Geometry &geom, const bool &completed)
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
}

// RK - average of edge near points radius
void unitize_nearpoints_radius(Geometry &geom)
{
  double min = 0;
  double max = 0;
  Vec3d center;
  double avg = edge_nearpoints_radius(geom, min, max, center);
  geom.transform(Trans3d::scale(1 / avg));
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


void generate_points(const Geometry &base, const Geometry &dual, vector<Vec3d> &ip,
                    vector<Vec3d> &base_nearpts, vector<Vec3d> &dual_nearpts,
                    const cn_opts &opts)
{
  vector<vector<int>> base_edges;
  vector<vector<int>> dual_edges;

  if ((opts.output_parts.find("i") != string::npos) ||
      (opts.output_parts.find("n") != string::npos) ||
      (opts.output_parts.find("p") != string::npos))
    base.get_impl_edges(base_edges);

  if ((opts.output_parts.find("i") != string::npos) ||
      (opts.output_parts.find("m") != string::npos) ||
      (opts.output_parts.find("q") != string::npos))
    dual.get_impl_edges(dual_edges);

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
          ip.push_back(intersection_point);
          break;
        }
      }
    }

    if (ip.size() != base_edges.size()) {
      if (opts.canonical_method != 'x')
        fprintf(stderr,"Warning: only %d out of %d intersection points found. try more precision\n",
                (int)ip.size(), (int)base_edges.size());
      else
        fprintf(stderr,"Warning: only canonical models have intersection points\n");
    }
  }

  if ((opts.output_parts.find("n") != string::npos) ||
      (opts.output_parts.find("p") != string::npos)) {
    for (unsigned int e = 0; e < base_edges.size(); e++)
      base_nearpts.push_back(base.edge_nearpt(base_edges[e], Vec3d(0, 0, 0)));
  }

  if ((opts.output_parts.find("m") != string::npos) ||
      (opts.output_parts.find("q") != string::npos)) {
    for (unsigned int e = 0; e < dual_edges.size(); e++)
      dual_nearpts.push_back(dual.edge_nearpt(dual_edges[e], Vec3d(0, 0, 0)));
  }
}

void set_edge_colors(Geometry &geom, const Color &col)
{
  if (col.is_set()) {
    geom.add_missing_impl_edges();
    Coloring clrng(&geom);
    clrng.e_one_col(col);
  }
}

void construct_model(Geometry &base, const cn_opts &opts) {
  Geometry dual;

  // need to generate dual? also needed for tangent points
  if ((opts.output_parts.find("d") != string::npos) ||
      (opts.output_parts.find("i") != string::npos) ||
      (opts.output_parts.find("m") != string::npos) ||
      (opts.output_parts.find("q") != string::npos))
    get_dual(&dual, base, 1, Vec3d(0, 0, 0));

  // need to generate before possible erasure of base
  vector<Vec3d> ip;
  vector<Vec3d> base_nearpts;
  vector<Vec3d> dual_nearpts;
  if ((opts.output_parts.find("i") != string::npos) ||
      (opts.output_parts.find("n") != string::npos) ||
      (opts.output_parts.find("m") != string::npos) ||
      (opts.output_parts.find("p") != string::npos) ||
      (opts.output_parts.find("q") != string::npos))
    generate_points(base, dual, ip, base_nearpts, dual_nearpts, opts);

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
    for (unsigned int i = 0; i < ip.size(); i++)
      base.add_vert(ip[i], opts.ipoints_col);
  }

  // add base near points
  if (opts.output_parts.find("n") != string::npos) {
    for (unsigned int i = 0; i < base_nearpts.size(); i++)
      base.add_vert(base_nearpts[i], opts.base_nearpts_col);
  }

  // add dual near points
  if (opts.output_parts.find("m") != string::npos) {
    for (unsigned int i = 0; i < dual_nearpts.size(); i++)
      base.add_vert(dual_nearpts[i], opts.dual_nearpts_col);
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
    for (int i = 0; i < (int)base.faces().size(); i++) {
      Color col = base.colors(FACES).get(i);
      if (col.is_set()) {
        col = Color(col[0], col[1], col[2], opts.face_opacity);
        base.colors(FACES).set(i, col);
      }
    }
  }

  // add unit sphere on origin
  if (opts.output_parts.find("u") != string::npos) {
    Geometry sgeom;
    sgeom.read_resource("geo_4_4");
    sgeom.transform(Trans3d::transl(-centroid(sgeom.verts())));
    unitize_vertex_radius(sgeom);
    Coloring(&sgeom).vef_one_col(Color::invisible, Color::invisible, opts.sphere_col);
    base.append(sgeom);
  }
}

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

    //geom.transform(Trans3d::transl(-centroid(geom.verts())));
    geom.transform(Trans3d::transl(-edge_nearpoints_centroid(geom, Vec3d(0, 0, 0))));

    unitize_nearpoints_radius(geom);

    //bool planarize_only = false;
    canonicalize_bd2(geom, opts.num_iters_canonical, 'b',
                    opts.radius_range_percent / 100, opts.rep_count, opts.centering, opts.epsilon);

    str += ".off";
    opts.write_or_error(geom, str);
  }
  exit(0);
*/

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);

  fprintf(stderr,"\n");
  fprintf(stderr,"centering: ");
  if (opts.centering == 'v') {
    fprintf(stderr, "(vertex centroid to origin)\n");
    geom.transform(Trans3d::transl(-centroid(geom.verts())));
  }
  else
  if (opts.centering == 'n') {
    fprintf(stderr, "(edge near points centroid to origin)\n");
    geom.transform(Trans3d::transl(-edge_nearpoints_centroid(geom, Vec3d(0, 0, 0))));
  }
  else
  if (opts.centering == 'x')
    fprintf(stderr, "(model not moved)\n");

  fprintf(stderr,"starting radius: ");
  if (opts.centering == 'v') {
    fprintf(stderr, "(average vertex)\n");
    unitize_vertex_radius(geom);
  }
  else
  if (opts.initial_radius == 'n') {
    fprintf(stderr, "(average edge near points)\n");
    unitize_nearpoints_radius(geom);
  }
  else
  if (opts.centering == 'x')
    fprintf(stderr, "(radius not changed)\n");

  if (opts.edge_distribution) {
    fprintf(stderr, "edge distribution: (project onto sphere)\n");
    if (opts.edge_distribution == 's')
      project_onto_sphere(&geom);
  }

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
    fprintf(stderr, "planarize: (%s method)\n",planarize_str.c_str());

    if (opts.planarize_method == 'm') {
      bool planarize_only = true;
      completed = canonicalize_mm2(geom, opts.mm_edge_factor / 100, opts.mm_plane_factor / 100,
                                  opts.num_iters_planar, opts.radius_range_percent / 100, opts.rep_count,
                                  planarize_only, opts.mm_alternate_loop, opts.epsilon);
    }
    else {
      completed = canonicalize_bd2(geom, opts.num_iters_planar, opts.planarize_method,
                                  opts.radius_range_percent / 100, opts.rep_count, opts.centering, opts.epsilon);
    }

    // RK - report planarity
    planarity_info(geom);
  }

  if (opts.canonical_method && opts.canonical_method != 'x') {
    fprintf(stderr, "canonicalize: (%s method)\n",((opts.canonical_method == 'm') ? "mathematica" : "base/dual"));
    if (opts.canonical_method == 'm') {
      bool planarize_only = false;
      completed = canonicalize_mm2(geom, opts.mm_edge_factor / 100, opts.mm_plane_factor / 100,
                                  opts.num_iters_canonical, opts.radius_range_percent / 100, opts.rep_count,
                                  planarize_only, opts.mm_alternate_loop, opts.epsilon);
    }
    else
      completed = canonicalize_bd2(geom, opts.num_iters_canonical, opts.canonical_method,
                                  opts.radius_range_percent / 100, opts.rep_count, opts.centering, opts.epsilon);

    // RK - report planarity
    planarity_info(geom);

    // RK - print midradius info
    midradius_info(geom, completed);
  }

  // RK - parts to output
  construct_model(geom, opts);

  opts.write_or_error(geom, opts.stderr);

  return 0;
}
