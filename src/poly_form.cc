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
   Name: poly_form.cc
   Description: make face-regular, equal-edge and unscrambled polyhedra
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <string>
#include <vector>

using std::set;
using std::string;
using std::vector;

using namespace anti;

class pf_opts : public ProgramOpts {
public:
  IterationControl it_ctrl;
  char algm = '\0';            // algorithm, use default
  double shorten_by = 1.0;     // shorten adjustment, a small amount
  double shorten_rad_by = NAN; // shorten radius adjustment, use default
  double flatten_by = NAN;     // flatten adjustment, use default
  bool use_symmetry = false;   // don't use symmetry
  Vec4d ellipsoid;             // projection ellipsoid, unset, use default

  string ifile;
  string ofile;

  pf_opts() : ProgramOpts("poly_form")
  {
    // The following defaults are suitable for small and medium models
    it_ctrl.set_max_iters(10000); // will finish reasobly quickly
    it_ctrl.set_status_check_and_report_iters(1000);
    // Test for default algorithm -a r not cheap, but enable early termination
    it_ctrl.set_status_check_only_iters(100);
    it_ctrl.set_sig_digits(14); // achievable in reasonable time
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

void pf_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] [input_file]

Read a file in OFF format containing a graph of a polyhedron, with or
without vertex coordinates, and try to create a spherical or ellipsoidal
tesselation where the maximum edge is a minimum length, or try to make
into a regular-faced polyhedron. Option adjustment factors expressed as
a percentage are approximate. If input_file is not given the program reads
from standard input.

Options
%s
  -a <alg>  model forming algorithm
            r - make faces into unit-edged regular polygons (default)
            u - unscramble: place a small polygon on one side of a
                sphere or ellipsoid and the rest of the vertices at a
                point on the other side, then equalise shortest and
                longest edges attached to a vertex
            U - equalise shortest and longest edges attached to a vertex
                  on sphere or ellipsoid (use to continue unscrambling)
            e - equalise edges on sphere or ellipsoid (default)
            p - planarise
  -n <itrs> maximum number of iterations, -1 for unlimited (default: %d)
  -s <perc> percentage to shorten longest edges on iteration (default: 1)
  -k <perc> -a r - percentage to reduce polygon radius on iteration
                   (default: value of -s)
            -a e - percentage to reduce minimum target length on iteration
                   (default: 1e-6), will oscillate at certain precision,
                   process output with smaller value to improve solution
  -f <perc> percentage to reduce distance of vertex from face plane (-a rp)
            on iteration (default: value of -s)
  -y        maintain symmetry of the base model (with -a e/r/p)
  -E <prms> use ellipsoid, three numbers separated by commas are the
            axis lengths (for a superellipsoid an optional fourth number
            gives the power)
  -l <lim>  minimum change of distance/width_of_model to terminate, as 
               negative exponent (default: %d giving %.0e)
  -z <nums> number of iterations between status reports (implies termination
            check) (0 for final report only, -1 for no report), optionally
            followed by a comma and the number of iterations between
            termination checks (0 for report checks only) (default: %d,%d)
  -o <file> write output to file (default: write to standard output)

)",
          prog_name(), help_ver_text, it_ctrl.get_max_iters(),
          it_ctrl.get_sig_digits(), it_ctrl.get_test_val(),
          it_ctrl.get_status_check_and_report_iters(),
          it_ctrl.get_status_check_only_iters());
}

void pf_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;
  int num;
  vector<double> nums;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hn:s:l:k:f:a:yE:z:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'o':
      ofile = optarg;
      break;

    case 'n':
      print_status_or_exit(read_int(optarg, &num), c);
      print_status_or_exit(it_ctrl.set_max_iters(num), c);
      break;

    case 'z':
      print_status_or_exit(it_ctrl.set_status_checks(optarg), c);
      break;

    case 's':
      print_status_or_exit(read_double(optarg, &shorten_by), c);
      if (shorten_by < 0 || shorten_by > 100)
        warning("not in range 0 to 100", c);
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
      if (strlen(optarg) > 1 || !strchr("ruUep", *optarg))
        error("method is '" + string(optarg) + "' must be from ruUep");
      algm = *optarg;
      break;

    case 'y':
      use_symmetry = true;
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

    case 'l':
      print_status_or_exit(read_int(optarg, &num), c);
      print_status_or_exit(it_ctrl.set_sig_digits(num), c);
      break;

    default:
      error("unknown command line error");
    }
  }

  // Default algorithm is 'e', but if unscrambling it is 'u'
  if (!algm)
    algm = 'r';

  if (algm == 'r') {
    if (std::isnan(shorten_rad_by))
      shorten_rad_by = shorten_by;
    if (std::isnan(flatten_by))
      flatten_by = shorten_by;
    if (ellipsoid.is_set())
      warning("set, but not used for this algorithm", 'E');
  }
  else if (algm == 'u' || algm == 'U') {
    if (use_symmetry)
      error(msg_str("cannot maintain symmetry with algorithm '%c'", algm), 'y');
    if (!std::isnan(shorten_rad_by))
      warning("set, but not used for this algorithm", 'k');
    if (!std::isnan(flatten_by))
      warning("set, but not used for this algorithm", 'f');
  }
  else if (algm == 'e') {
    if (std::isnan(shorten_rad_by))
      shorten_rad_by = 1e-6;
    if (!std::isnan(flatten_by))
      warning("set, but not used for this algorithm", 'f');
  }
  else if (algm == 'p') {
    if (std::isnan(flatten_by))
      flatten_by = shorten_by;
    if (!std::isnan(shorten_rad_by))
      warning("set, but not used for this algorithm", 'k');
    if (ellipsoid.is_set())
      warning("set, but not used for this algorithm", 'E');
  }

  if (argc - optind > 1)
    error("too many arguments");

  if (argc - optind == 1)
    ifile = argv[optind];
}

void to_ellipsoid(Vec3d &v, const Vec4d &ellipsoid)
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

void to_ellipsoid(Geometry &geom, const Vec4d &ellipsoid)
{
  for (unsigned int i = 0; i < geom.verts().size(); i++)
    to_ellipsoid(geom.verts(i), ellipsoid);
}

void to_initial_unscramble(Geometry &geom, const Vec4d &ellipsoid)
{
  for (unsigned int i = 0; i < geom.verts().size(); i++)
    geom.verts(i) = Vec3d(-1, 0, 0);
  if (geom.faces().size()) {
    const vector<int> &face = geom.faces(0);
    for (unsigned int i = 0; i < face.size(); i++) {
      geom.verts(face[i]) = Vec3d(1, 0.01 * cos(i * 2 * M_PI / face.size()),
                                  0.01 * sin(i * 2 * M_PI / face.size()));
      to_ellipsoid(geom.verts(face[i]), ellipsoid);
    }
  }
}

//------------------------------------------------------------------
// Unscramble

Status make_unscramble(Geometry &geom, IterationControl it_ctrl,
                       double shorten_factor, Vec4d ellipsoid,
                       bool with_initial_placement)
{
  if (with_initial_placement)
    to_initial_unscramble(geom, ellipsoid);

  // No further processing if no faces, but not an error
  if (geom.faces().size() == 0)
    return Status::ok();

  auto eds = GeometryInfo(geom).get_vert_cons();

  double g_max_dist = 0;
  double g_min_dist = 1e100;
  for (it_ctrl.start_iter(); !it_ctrl.is_done(); it_ctrl.next_iter()) {
    g_max_dist = 0;
    g_min_dist = 1e100;
    for (unsigned int v = 0; v < geom.verts().size(); v++) {
      if (eds[v].size() == 0)
        continue;
      double max_dist = -1;
      double min_dist = 1e100;
      int max_edge = 0;
      for (unsigned int i = 0; i < eds[v].size(); i++) {
        double dist = (geom.verts(v) - geom.verts(eds[v][i])).len2();
        if (dist > max_dist) {
          max_dist = dist;
          max_edge = i;
        }
        if (dist > g_max_dist)
          g_max_dist = dist;
        if (dist < min_dist)
          min_dist = dist;
        if (dist < g_min_dist)
          g_min_dist = dist;
      }

      int p0 = v;
      int p1 = eds[v][max_edge];
      Vec3d diff = geom.verts(p1) - geom.verts(p0);
      geom.verts(p0) += diff * shorten_factor;
      to_ellipsoid(geom.verts(p0), ellipsoid);
    }

    // no status checks
    // if (it_ctrl.is_status_check_iter()) { ; }

    if (it_ctrl.is_status_report_iter()) {
      if (it_ctrl.is_finished())
        it_ctrl.print("Final iteration:\n");

      it_ctrl.print("%-12u max:%17.15f min:%17.15f\n",
                    it_ctrl.get_current_iter(), sqrt(g_max_dist),
                    sqrt(g_min_dist));
    }
  }

  return Status::ok();
}

//------------------------------------------------------------------
// Equal edge

Status make_equal_edges(Geometry &base_geom, IterationControl it_ctrl,
                        double shorten_factor, double shrink_factor,
                        Vec4d ellipsoid, const Symmetry &sym)
{
  to_ellipsoid(base_geom, ellipsoid); // map before finding symmetry

  Status stat;

  bool using_symmetry = (sym.get_sym_type() > Symmetry::C1);

  // No further processing if no faces, but not an error
  if (base_geom.faces().size() == 0)
    return Status::ok();

  SymmetricUpdater sym_updater((using_symmetry) ? base_geom : Geometry(), sym);
  const Geometry &geom =
      (using_symmetry) ? sym_updater.get_geom_working() : base_geom;

  const vector<Vec3d> &verts = geom.verts();
  auto eds = GeometryInfo(geom).get_vert_cons();

  vector<int> principal_verts; // first vertex in each vertex orbit
  vector<int> verts_to_update; // vertices accessed during iteration
  if (using_symmetry) {
    principal_verts = sym_updater.get_principal(VERTS);
    verts_to_update =
        SymmetricUpdater::get_included_verts(principal_verts, eds);
  }
  else {
    principal_verts =
        SymmetricUpdater::sequential_index_list(geom.verts().size());
  }

  double g_max_dist = 0;
  double g_min_dist = 1e100;
  double g_scale_factor = 1;
  double max_dist = 0;
  double min_dist = 1e100;
  double test_val = it_ctrl.get_test_val();

  vector<Vec3d> offsets(verts.size()); // Vertex adjustments

  for (it_ctrl.start_iter_with_setup(); !it_ctrl.is_done();
       it_ctrl.next_iter()) {
    std::fill(offsets.begin(), offsets.end(), Vec3d::zero);

    if (using_symmetry) {
      // Ensure that the vertices used in adjustment are up to date
      for (auto v_idx : verts_to_update)
        sym_updater.update_from_principal_vertex(v_idx);
    }

    max_dist = 0;
    min_dist = 1e100;
    for (int v_idx : principal_verts) {
      if (eds[v_idx].size() == 0)
        continue;
      for (unsigned int i = 0; i < eds[v_idx].size(); i++) {
        auto vec = verts[v_idx] - verts[eds[v_idx][i]];
        double dist = vec.len();
        if (dist > max_dist)
          max_dist = dist;
        if (dist < min_dist)
          min_dist = dist;
        offsets[v_idx] += (vec / dist) * (g_min_dist - dist) * g_scale_factor *
                          shorten_factor;
      }
    }

    // adjust vertices post-loop (skip setup iter as global values not set)
    if (!it_ctrl.is_setup_iter()) {
      for (int v_idx : principal_verts) {
        auto vert = verts[v_idx] + offsets[v_idx];
        to_ellipsoid(vert, ellipsoid);
        if (using_symmetry)
          sym_updater.update_principal_vertex(v_idx, vert);
        else
          base_geom.verts(v_idx) = vert;
      }
    }

    g_min_dist = min_dist * (1 - shrink_factor);
    g_max_dist = max_dist;
    g_scale_factor = (g_max_dist - g_min_dist) / g_min_dist;

    // Do not check status or finish before modifying the model
    if (!it_ctrl.is_setup_iter()) {
      string finish_msg;
      if (it_ctrl.is_status_check_iter()) {

        if ((max_dist - min_dist) < test_val) { // absolute difference
          it_ctrl.set_finished();
          finish_msg = "solved, test value achieved";
        }
        else if (it_ctrl.is_last_iter()) {
          it_ctrl.set_finished();
          finish_msg = "not solved, test value not achieved";
        }
      }

      if (it_ctrl.is_status_report_iter()) {
        if (it_ctrl.is_finished())
          it_ctrl.print("Final iteration (%s):\n", finish_msg.c_str());

        it_ctrl.print("%-12u max:%17.15f min:%17.15f diff:%.11g\n",
                      it_ctrl.get_current_iter(), max_dist, min_dist,
                      max_dist - min_dist);
      }
    }
  }

  if (using_symmetry)
    base_geom = sym_updater.get_geom_final();

  return Status::ok();
}

//------------------------------------------------------------------
// Unit edge regular polygon

Status make_regular_faces(Geometry &base_geom, IterationControl it_ctrl,
                          double shorten_factor, double plane_factor,
                          double radius_factor, const Symmetry &sym)
{
  Status stat;

  bool using_symmetry = (sym.get_sym_type() > Symmetry::C1);

  // No further processing if no faces, but not an error
  if (base_geom.faces().size() == 0)
    return Status::ok();

  double test_val = it_ctrl.get_test_val();
  const double divergence_test2 = 1e30; // test vertex dist^2 for divergence
  // Scale to get edges close to 1
  GeometryInfo info(base_geom);
  double scale = info.iedge_length_lims().sum / info.num_iedges();
  if (scale)
    base_geom.transform(Trans3d::scale(1 / scale));

  SymmetricUpdater sym_updater((using_symmetry) ? base_geom : Geometry(), sym);
  const Geometry &geom =
      (using_symmetry) ? sym_updater.get_geom_working() : base_geom;

  const vector<Vec3d> &verts = geom.verts();
  const vector<vector<int>> &faces = geom.faces();

  double max_diff2 = 0;
  Vec3d origin(0, 0, 0);
  vector<double> rads(faces.size());
  for (unsigned int f = 0; f < faces.size(); f++) {
    int N = faces[f].size();
    int D = std::abs(find_polygon_denominator_signed(geom, f, epsilon));
    if (!D)
      D = 1;
    rads[f] = 0.5 / sin(M_PI * D / N); // circumradius of regular polygon
  }

  vector<int> principal_verts;  // first vertex in each vertex orbit
  vector<int> verts_to_update;  // vertices accessed during iteration
  vector<int> faces_to_process; // faces processed during iteration
  vector<int> edges_to_process; // edges processed during iteration
  if (using_symmetry) {
    principal_verts = sym_updater.get_principal(VERTS);
    faces_to_process = sym_updater.get_associated_elems(info.get_vert_faces());
    edges_to_process =
        sym_updater.get_associated_elems(info.get_vert_impl_edges());
    // All the vertices that need to be updated lie on the faces to be processed
    verts_to_update =
        SymmetricUpdater::get_included_verts(faces_to_process, geom.faces());
  }
  else { // not using_symmetry
    principal_verts =
        SymmetricUpdater::sequential_index_list(geom.verts().size());
    faces_to_process =
        SymmetricUpdater::sequential_index_list(geom.faces().size());
    edges_to_process =
        SymmetricUpdater::sequential_index_list(info.get_impl_edges().size());
  }

  vector<Vec3d> offsets(verts.size()); // Vertex adjustments

  for (it_ctrl.start_iter(); !it_ctrl.is_done(); it_ctrl.next_iter()) {
    std::fill(offsets.begin(), offsets.end(), Vec3d::zero);

    if (using_symmetry) {
      // Ensure that the vertices used in adjustment are up to date
      for (int v_idx : verts_to_update)
        sym_updater.update_from_principal_vertex(v_idx);
    }

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
    for (int e_idx : edges_to_process) {
      const auto &edge = info.get_impl_edges()[e_idx];
      Vec3d offset =
          (1 - geom.edge_len(edge)) * shorten_factor * geom.edge_vec(edge);
      offsets[edge[0]] -= 2 * offset;
      offsets[edge[1]] += 2 * offset;
    }

    // adjust vertices post-loop
    if (using_symmetry) {
      // adjust principal vertices
      for (int v_idx : principal_verts)
        sym_updater.update_principal_vertex(v_idx,
                                            verts[v_idx] + offsets[v_idx]);
    }
    else { // not using_symmetry
      // adjust all vertices
      for (unsigned int i = 0; i < verts.size(); i++)
        base_geom.raw_verts()[i] += offsets[i];
    }

    string finish_msg;
    if (it_ctrl.is_status_check_iter()) {

      max_diff2 = 0;
      for (auto &offset : offsets) {
        double diff2 = offset.len2();
        if (diff2 > max_diff2)
          max_diff2 = diff2;
      }

      double width = BoundBox(verts).max_width();
      if (sqrt(max_diff2) / width < test_val) {
        it_ctrl.set_finished();
        finish_msg = "solved, test value achieved";
      }
      else if (it_ctrl.is_last_iter()) {
        // reached last iteration without solving
        it_ctrl.set_finished();
        finish_msg = "not solved, test value not achieved";
      }

      // check if radius is expanding or contracting unreasonably,
      // but only for the purpose of finishing early
      if (!it_ctrl.is_finished() && divergence_test2 > 0) {
        for (auto &offset : offsets)
          if (offset.len2() > divergence_test2) {
            it_ctrl.set_finished();
            finish_msg = "not solved, quit early as probably diverging";
          }
      }
    }

    if (it_ctrl.is_status_report_iter()) {
      if (it_ctrl.is_finished())
        it_ctrl.print("Final iteration (%s):\n", finish_msg.c_str());

      it_ctrl.print("%-12u max_diff:%17.15e\n", it_ctrl.get_current_iter(),
                    sqrt(max_diff2));
    }
  }

  if (using_symmetry)
    base_geom = sym_updater.get_geom_final();

  return Status::ok();
}

//------------------------------------------------------------------
// Make faces planar

inline Vec3d nearpoint_on_plane(const Vec3d &P, const Vec3d &point_on_plane,
                                const Vec3d &unit_norm)
{
  return P + vdot(point_on_plane - P, unit_norm) * unit_norm;
}

Status make_planar(Geometry &base_geom, IterationControl it_ctrl,
                   double plane_factor, const Symmetry &sym)
{
  // chosen by experiment
  const double intersect_test_val = 1e-5; // test for coplanar faces
  const double diff2_test_val = 10;       // limit for using plane intersection
  const double readjustment = 1.01;       // to adjust adjustment factor
  const double plane_factor_max = 1.1;    // maximum value for adjustment factor
  bool using_symmetry = (sym.get_sym_type() > Symmetry::C1);

  SymmetricUpdater sym_updater((using_symmetry) ? base_geom : Geometry(), sym);
  const Geometry &geom =
      (using_symmetry) ? sym_updater.get_geom_working() : base_geom;

  // No further processing if no faces, but not an error
  if (geom.faces().size() == 0)
    return Status::ok();

  const vector<Vec3d> &verts = geom.verts();
  const vector<vector<int>> &faces = geom.faces();

  // List of faces that a vertex is part of
  GeometryInfo info(geom);
  const auto &vert_faces = info.get_vert_faces();

  // Get a list of the faces that contain a principal vertex of any type
  vector<int> principal_verts;  // first vertex in each vertex orbit
  vector<int> verts_to_update;  // vertices accessed during iteration
  vector<int> faces_to_process; // faces processed during iteration
  if (using_symmetry) {
    principal_verts = sym_updater.get_principal(VERTS);
    faces_to_process = sym_updater.get_associated_elems(vert_faces);
    verts_to_update =
        SymmetricUpdater::get_included_verts(faces_to_process, geom.faces());
  }
  else { // not using_symmetry
    principal_verts =
        SymmetricUpdater::sequential_index_list(geom.verts().size());
    faces_to_process =
        SymmetricUpdater::sequential_index_list(geom.faces().size());
  }

  // Use oversized arrays to avoid mapping
  vector<Vec3d> offsets(verts.size()); // Vertex adjustments
  vector<Vec3d> norms(faces.size());   // Face normals
  vector<Vec3d> cents(faces.size());   // Face centroids

  double test_val = it_ctrl.get_test_val();
  double last_max_diff2 = 0.0;
  for (it_ctrl.start_iter(); !it_ctrl.is_done(); it_ctrl.next_iter()) {
    std::fill(offsets.begin(), offsets.end(), Vec3d::zero);

    if (using_symmetry) {
      // Ensure that the vertices used in adjustment are up to date
      for (auto v_idx : verts_to_update)
        sym_updater.update_from_principal_vertex(v_idx);
    }

    // Initialize face data for just the necessary faces
    for (auto f_idx : faces_to_process) {
      norms[f_idx] = geom.face_norm(f_idx).unit();
      cents[f_idx] = geom.face_cent(f_idx);
    }

    int cnt_proj = 0;
    int cnt_int = 0;
    double max_diff2 = 0.0;
    for (auto v_idx : principal_verts) {
      int intersect_cnt = 0;
      const auto &vfaces = vert_faces[v_idx];
      const int vf_sz = vfaces.size();
      bool good_intersections = (vf_sz >= 3);
      for (int f0 = 0; f0 < vf_sz - 2 && good_intersections; f0++) {
        int f0_idx = vfaces[f0];
        for (int f1 = f0 + 1; f1 < vf_sz - 1 && good_intersections; f1++) {
          int f1_idx = vfaces[f1];
          for (int f2 = f1 + 1; f2 < vf_sz && good_intersections; f2++) {
            int f2_idx = vfaces[f2];
            Vec3d intersection;
            good_intersections = three_plane_intersect(
                cents[f0_idx], norms[f0_idx], cents[f1_idx], norms[f1_idx],
                cents[f2_idx], norms[f2_idx], intersection, intersect_test_val);
            offsets[v_idx] += intersection;
            intersect_cnt++;
          }
        }
      }

      if (good_intersections)
        offsets[v_idx] =
            (offsets[v_idx] / intersect_cnt - verts[v_idx]) * plane_factor;

      // no good 3 plane intersections OR
      // moving to much
      if (!good_intersections ||
          offsets[v_idx].len2() / last_max_diff2 > diff2_test_val) {
        // target vertex is centroid of projection of vertex onto planes
        offsets[v_idx] = Vec3d::zero;
        for (int f0 = 0; f0 < vf_sz; f0++) {
          int f0_idx = vfaces[f0];
          offsets[v_idx] +=
              nearpoint_on_plane(verts[v_idx], cents[f0_idx], norms[f0_idx]);
        }
        offsets[v_idx] = (offsets[v_idx] / vf_sz - verts[v_idx]) * plane_factor;
        cnt_proj++;
      }
      else
        cnt_int++;

      auto diff2 = offsets[v_idx].len2();
      if (diff2 > max_diff2)
        max_diff2 = diff2;
    }

    // adjust vertices post-loop
    if (using_symmetry) {
      // adjust principal vertices
      for (int v_idx : principal_verts)
        sym_updater.update_principal_vertex(v_idx,
                                            verts[v_idx] + offsets[v_idx]);
    }
    else { // not using_symmetry
      // adjust all vertices
      for (unsigned int i = 0; i < verts.size(); i++)
        base_geom.raw_verts()[i] += offsets[i];
    }

    // adjust plane factor
    if (max_diff2 < last_max_diff2)
      plane_factor *= readjustment;
    else
      plane_factor /= readjustment;
    if (plane_factor > plane_factor_max)
      plane_factor = plane_factor_max;
    last_max_diff2 = max_diff2;

    string finish_msg;
    if (it_ctrl.is_status_check_iter()) {
      double width = BoundBox(verts).max_width();
      if (sqrt(max_diff2) / width < test_val) {
        it_ctrl.set_finished();
        finish_msg = "solved, test value achieved";
      }
      else if (it_ctrl.is_last_iter()) {
        // reached last iteration without solving
        it_ctrl.set_finished();
        finish_msg = "not solved, test value not achieved";
      }
    }

    if (it_ctrl.is_status_report_iter()) {
      if (it_ctrl.is_finished())
        it_ctrl.print("Final iteration (%s):\n", finish_msg.c_str());

      it_ctrl.print("%-12u max_diff:%17.15e  -f %-10.5f i/p: %7d/%-7d\n",
                    it_ctrl.get_current_iter(), sqrt(max_diff2),
                    100 * plane_factor, cnt_int, cnt_proj);
    }
  }

  if (using_symmetry)
    base_geom = sym_updater.get_geom_final();

  return Status::ok();
}

//------------------------------------------------------------------
// Main

int main(int argc, char *argv[])
{
  pf_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);

  if (opts.ellipsoid.is_set()) // will only be set for algorithms where valid
    to_ellipsoid(geom, opts.ellipsoid);

  Symmetry sym;
  if (opts.use_symmetry)
    opts.print_status_or_exit(sym.init(geom), 'y');

  if (geom.faces().size() == 0)
    opts.warning("no face date in input, no iterative processing will occur");

  if (opts.algm == 'r') {
    opts.print_status_or_exit(make_regular_faces(
        geom, opts.it_ctrl, opts.shorten_by / 200, opts.flatten_by / 100,
        opts.shorten_rad_by / 100, sym));
  }
  else if (opts.algm == 'u' || opts.algm == 'U') {
    opts.print_status_or_exit(
        make_unscramble(geom, opts.it_ctrl, opts.shorten_by / 200,
                        opts.ellipsoid, (opts.algm == 'u')));
  }
  else if (opts.algm == 'e') {
    opts.print_status_or_exit(
        make_equal_edges(geom, opts.it_ctrl, opts.shorten_by / 200,
                         opts.shorten_rad_by / 100, opts.ellipsoid, sym));
  }
  else if (opts.algm == 'p') {
    opts.print_status_or_exit(
        make_planar(geom, opts.it_ctrl, opts.flatten_by / 100, sym));
  }

  opts.write_or_error(geom, opts.ofile);

  return 0;
}
