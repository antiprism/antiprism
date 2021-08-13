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
   Description: make rotegrity and nexorade models
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using std::string;
using std::vector;

using namespace anti;

class rot_opts : public ProgramOpts {
public:
  IterationControl it_ctrl;
  char algorithm = 'r';    // rotegrity
  double strut_len = 0;    // strut length
  double adjust_fact = 98; // generally efficient, chosen from testing
  double end_fraction = 1 / 3.0;
  bool already_twisted = false;
  int method = 0;         // unset
  char output_type = 'f'; // full
  int col_type = 'x';     // unset
  Coloring clrngs[3];
  string ifile;
  string ofile;

  rot_opts() : ProgramOpts("rotegrity")
  {
    read_colorings(clrngs, "spread");
    it_ctrl.set_max_iters(10000); // will finish reasobly quickly
    it_ctrl.set_status_check_and_report_iters(1000);
    it_ctrl.set_status_check_only_iters(
        1);                     // test cheap, enable early termination
    it_ctrl.set_sig_digits(15); // achievable in reasonable time
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

void rot_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] [input_file]

Read a file in OFF format containing a roughly spherical polyhedron, or
previously twisted model, and try to convert into a rotegrity or nexorade.
Units keep original edge colours. If input_file is not given the program
reads from standard input.

Options
%s
  -a <type> model type: rotegrity, nexorade, for nexorade followed
            by an optional comma and strut length
  -f <frac> fraction of length for end sections (default: 1/3)
  -t        input model is already twisted (produced by -O f)
  -M <mthd> method of conversion from base model - twist, double, joined,
            or X (default: t)
  -O <type> output type for units: full (face), rotegrity (3 short struts),
            nexorade (long strut), Nexorade (long strut, direction vertices)
            (default: full). Only 'full' output can be used as input with
            option -t
  -m <maps> a comma separated list of colour maps used to transform colour
            indexes (default: rand), a part consisting of letters from
            v, e, f, selects the element types to apply the map list to
            (default 'vef')
  -c <type> colouring type: edge (base model edges), symmetry, unit
            (according to shape), none (default: edge)
  -s <perc> percentage to adjust corrections on iteration (default: %.0f)
  -n <itrs> maximum number of iterations, -1 for unlimited (default: %d)
  -l <lim>  minimum change of distance/width_of_model to terminate, as 
               negative exponent (default: %d giving %.0e)
  -z <nums> number of iterations between status reports (implies termination
            check) (0 for final report only, -1 for no report), optionally
            followed by a comma and the number of iterations between
            termination checks (0 for report checks only) (default: %d,%d)
  -o <file> write output to file (default: write to standard output)

)",
          prog_name(), help_ver_text, adjust_fact, it_ctrl.get_max_iters(),
          it_ctrl.get_sig_digits(), it_ctrl.get_test_val(),
          it_ctrl.get_status_check_and_report_iters(),
          it_ctrl.get_status_check_only_iters());
}

void rot_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;
  int num;
  vector<double> nums;
  string arg_id;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":ha:f:tM:O:c:m:n:s:l:z:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'a': {
      Split parts(optarg, ",");
      if (parts.size() < 1 || parts.size() > 2)
        error("argument should have 1-2 comma separated parts", c);
      print_status_or_exit(get_arg_id(parts[0], &arg_id,
                                      "rotegrity=r|nexorade=n|tensegrity=t",
                                      argmatch_add_id_maps),
                           c);
      algorithm = arg_id[0];
      if (parts.size() > 1) {
        if (algorithm == 'n' || algorithm == 't') {
          print_status_or_exit(read_double(parts[1], &strut_len), c);
          if (strut_len < 0)
            error("length of the strut cannot be negative", c);
        }
        else
          error("algorithm does not have a strut length parameter", c);
      }
      break;
    }

    case 'f':
      print_status_or_exit(read_double(optarg, &end_fraction), c);
      if (end_fraction < 0.0 || end_fraction > 0.5)
        warning("not in range 0 to 0.5", c);
      break;

    case 't':
      already_twisted = true;
      break;

    case 'M':
      print_status_or_exit(get_arg_id(optarg, &arg_id,
                                      "twist=1|double=2|joined=3|X=4",
                                      argmatch_add_id_maps),
                           c);
      method = atoi(arg_id.c_str());
      break;

    case 'O':
      print_status_or_exit(
          get_arg_id(optarg, &arg_id,
                     "full=f|rotegrity=r|nexorade=n|Nexorade=N|tensegrity=t",
                     argmatch_case_sensitive),
          c);
      output_type = arg_id[0];
      break;

    case 'c':
      print_status_or_exit(get_arg_id(optarg, &arg_id,
                                      "edge=e|symmetry=s|unit=u|none=n",
                                      argmatch_default),
                           c);
      col_type = arg_id[0];
      break;

    case 'm':
      print_status_or_exit(read_colorings(clrngs, optarg), c);
      break;

    case 's':
      print_status_or_exit(read_double(optarg, &adjust_fact), c);
      if (adjust_fact < 0 || adjust_fact > 100)
        warning("not in range 0 to 100", c);
      break;

    case 'n':
      print_status_or_exit(read_int(optarg, &num), c);
      print_status_or_exit(it_ctrl.set_max_iters(num), c);
      break;

    case 'z':
      print_status_or_exit(it_ctrl.set_status_checks(optarg), c);
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

  if (argc - optind > 1)
    error("too many arguments");

  if (argc - optind == 1)
    ifile = argv[optind];

  if (already_twisted) {
    if (method)
      error("cannot use option -M with previously twisted input", 't');
    if (col_type != 'x')
      error("cannot use option -c to colour previously twisted input", 't');
  }

  if (method == 0)
    method = 1; // default: 'twist'

  if (col_type == 'x') // unset
    col_type = 'e';    // default: 'edge'
}

Status check_polyhedron_model(GeometryInfo &info)
{
  if (!info.is_polyhedron())
    return Status::error("polyhedron model: not connected like polyhedron");
  if (!info.is_orientable())
    return Status::error("polyhedron model: not orientable");
  if (info.genus() != 0)
    return Status::warning("polyhedron model: not connected like a sphere");

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
                int method, double ratio, char col_type)
{
  // set the weight, avoiding models with coincident vertices
  double wt = 2 * ratio;
  int wt_int = round(wt);
  double wt_frac = wt - wt_int;
  if (double_eq(wt_frac, 0.0, 0.001))
    wt += 0.001 * sgn(wt_frac);

  string pattern;
  if (method == 1)
    pattern = msg_str("[%gV%gE]0V0E0fe", 1 - wt, wt);
  else if (method == 2 || method == 4)
    pattern = msg_str("[%gV%gE,%gE%gF]0V0_1F1", 1 - wt, wt, wt, 1 - wt);
  else if (method == 3)
    pattern = msg_str("[2VE-0.3F,2F3E,3FV]0V0_1F2_1");

  // Model suitable for 'edge' colouring type: col_type == 'e'
  wythoff_make_tiling(tw_geom, geom, pattern, true, false, TilingColoring("a"));
  if (col_type == 's') { // 'symmetry' colouring
    vector<vector<std::set<int>>> sym_equivs;
    get_equiv_elems(tw_geom, sym.get_trans(), &sym_equivs);
    Coloring coloring(&tw_geom);
    coloring.f_sets(sym_equivs[2], false);
  }
  else if (col_type == 'n') // 'none' colouring
    tw_geom.clear_cols();

  tw_geom.colors(VERTS).clear();

  int num_faces = tw_geom.faces().size();
  for (int i = 0; i < num_faces; i++) {
    auto face = tw_geom.faces(i);
    auto col = tw_geom.colors(FACES).get(i);
    int unit_fsz = 4;
    if (method == 1) {
      tw_geom.faces(i) = vector<int>(face.begin(), face.begin() + unit_fsz);
    }
    else if (method == 2) {
      tw_geom.faces(i) = vector<int>(face.begin(), face.begin() + unit_fsz);
      tw_geom.add_face(vector<int>(face.begin() + unit_fsz, face.end()), col);
    }
    else if (method == 3) {
      tw_geom.faces(i) = vector<int>(face.begin(), face.begin() + unit_fsz);
      tw_geom.add_face({face[5], face[6], face[7], face[8]}, col);
      tw_geom.add_face({face[4], face[3], face[8], face[9]}, col);
    }
    else if (method == 4) {
      vector<int> squ = {face[1], face[2], face[5], face[6]};
      auto cent = centroid(tw_geom.verts(), squ);
      int v_sz = tw_geom.verts().size();
      for (int j = 0; j < 4; j++)
        tw_geom.add_vert(
            (tw_geom.verts(squ[j]) + tw_geom.verts(squ[(j + 1) % 4]) + cent) /
            3);

      tw_geom.faces(i) = vector<int>({face[0], face[1], v_sz + 0, v_sz + 1});
      tw_geom.add_face({face[3], face[2], v_sz + 1, v_sz + 2}, col);
      tw_geom.add_face({face[4], face[5], v_sz + 2, v_sz + 3}, col);
      tw_geom.add_face({face[7], face[6], v_sz + 3, v_sz + 0}, col);
    }
  }

  for (unsigned int i = 0; i < tw_geom.faces().size(); i++) {
    auto &face = tw_geom.faces(i);
    std::rotate(face.begin(), face.begin() + 3, face.end());
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

//---------------------------------------------------------------------
// Comparison

/// Less, for doubles
struct DoubleLess {
  double limit;

  DoubleLess(double limit = epsilon) : limit(limit) {}
  bool operator()(const double &val1, const double &val2) const
  {
    return double_lt(val1, val2, limit);
  }
};

/// Less, for vectors of doubles
struct VectorDoubleLess {
  double limit;

  VectorDoubleLess(double limit = epsilon) : limit(limit) {}
  bool operator()(const std::vector<double> &vals1,
                  const std::vector<double> &vals2) const;
};

bool VectorDoubleLess::operator()(const vector<double> &vals1,
                                  const vector<double> &vals2) const
{
  int cmp = double_compare(vals1.size(), vals2.size());
  if (!cmp) {
    for (unsigned int i = 0; i < vals1.size(); i++)
      if ((cmp = double_compare(vals1[i], vals2[i], limit)))
        break;
  }
  return cmp < 0;
}
//---------------------------------------------------------------------

string report_rotegrity(Geometry &geom, double end_fraction, bool color)
{
  string report =
      "Two lines per rotegrity unit (angles in degrees)\n"
      "  1: unit_index_number (used for colouring), number_of_units\n"
      "  2: total_unit_angle, outer_segment_angle, inner_segment_angle\n"
      "\n";
  const double eps = 1e-10;
  map<vector<double>, vector<int>, VectorDoubleLess> units(
      (VectorDoubleLess(eps)));
  for (unsigned int i = 0; i < geom.faces().size(); i++) {
    auto ang = rad2deg(
        acos(safe_for_trig(vdot(geom.face_v(i, 1), geom.face_v(i, 0)))));

    units[{ang, ang * end_fraction, ang * (1 - 2 * end_fraction)}].push_back(i);
  }

  int unit_no = 0;
  for (const auto &unit : units) {
    report += msg_str("%-4d, %d\n", unit_no, (int)unit.second.size());
    report += msg_str("    %16.13f,%16.13f,%16.13f\n", unit.first[0],
                      unit.first[1], unit.first[2]);
    if (color) {
      for (auto f : unit.second)
        geom.colors(FACES).set(f, Color(unit_no));
    }
    unit_no++;
  }

  return report;
}

void make_rotegrity(Geometry &base_geom, const Symmetry &sym,
                    double end_fraction, IterationControl it_ctrl,
                    double factor)
{
  // Read and write to same model (Defered update is slower)
  SymmetricUpdater sym_updater(base_geom, sym);
  const Geometry &geom = sym_updater.get_geom_working();
  const vector<Vec3d> &verts = geom.verts();
  const vector<vector<int>> &faces = geom.faces();

  const auto principal_faces = sym_updater.get_principal(FACES);

  double test_val = it_ctrl.get_test_val();
  double max_diff = 0.0;
  for (it_ctrl.start_iter(); !it_ctrl.is_done(); it_ctrl.next_iter()) {
    max_diff = 0.0;
    for (int f_idx : principal_faces) {
      const auto &face = faces[f_idx];
      // Start: make sure vertices are current
      for (int i = 0; i < 4; i++)
        sym_updater.update_from_principal_vertex(face[i]);

      // Face has strut edge (long edge) first
      // Vertex index numbers and coordinates
      auto &v_end0 = verts[face[0]]; // first end
      auto &v_end1 = verts[face[1]]; // second end
      const auto &v_idx_mid1 = face[2];
      auto &v_mid1 = verts[v_idx_mid1]; // adjacent to second end
      const auto &v_idx_mid0 = face[3];
      auto &v_mid0 = verts[v_idx_mid0]; // adjacent to first end

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

    string finish_msg;
    if (it_ctrl.is_status_check_iter()) {

      if (max_diff < test_val) {
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

      it_ctrl.print("%-12u max_diff:%17.15e\n", it_ctrl.get_current_iter(),
                    max_diff);
    }
  }

  base_geom = sym_updater.get_geom_final();

  return;
}

string report_nexorade(Geometry &geom, double strut_len, bool color,
                       const vector<vector<int>> &f2fs)
{
  string report =
      "Radius and range of radius valaues, followed by minimum strut length\n"
      "to ensure contact, then three lines per strut type\n"
      "  1: strut_index_number (used for colouring), number_of_struts, "
      "strut_length\n"
      "  2: length_to_1st_contact, length_to_2nd_contact,\n"
      "     length_to_3rd_contact, length_to_4th_contact\n"
      "  3: angle_of_1st_contact, angle_of_2nd_contact,\n"
      "     angle_of_3rd_contact, angle_of_4th_contact\n"
      "Lengths are measured from either end of the strut\n"
      "Angles (degrees) are measured looking along the strut central line\n"
      "from the white end (display with -c u), anticlockwise, from an outward\n"
      "pointing zero angle.\n"
      "\n";
  const double eps = 1e-8;
  double rad_min = 1e100;
  double rad_max = 0.0;
  double len_max = 0.0;
  map<vector<double>, vector<int>, VectorDoubleLess> units(
      (VectorDoubleLess(eps)));
  vector<bool> f_flipped(geom.faces().size(), false);
  for (unsigned int i = 0; i < geom.faces().size(); i++) {
    // Project onto strut axis
    vector<Vec3d> v(4);     // contact points projected onto axis
    vector<Vec3d> vecs(4);  // contact points vectors perpendicular to axis
    vector<double> angs(4); // contact points vectors perpendicular to axis
    const Vec3d axis = geom.face_v(i, 1) - geom.face_v(i, 0);
    const double axis_len = axis.len();
    if (axis_len > len_max)
      len_max = axis_len;

    for (int j = 0; j < 4; j++) {
      if (j > 1) { // middle contact points
        v[j] = nearest_point(geom.face_v(i, j), geom.face_v(i, 0),
                             geom.face_v(i, 1));
        auto dist = (geom.face_v(i, j) - v[j]).len();
        if (dist < rad_min)
          rad_min = dist;
        else if (dist > rad_max)
          rad_max = dist;
      }
      else // end contact points
        v[j] = geom.face_v(i, j);

      if (j > 1) {
        vecs[j] = geom.face_v(i, j) - v[j];
      }
      else {
        const int f_idx = f2fs[i][j * 2];
        const auto near_pt = nearest_point(
            geom.face_v(i, j), geom.face_v(f_idx, 0), geom.face_v(f_idx, 1));
        vecs[j] = near_pt - geom.face_v(i, j);
      }

      angs[j] = (j == 0) ? 0.0 : angle_around_axis(vecs[0], vecs[j], axis);
    }

    // Choose angle 0 so end angles have same magnitude but different signs
    auto mid = angs[1] / 2;
    if (angs[1] > M_PI) // make the end angles large angles
      mid -= M_PI;
    for (int j = 0; j < 4; j++) {
      double ang = angs[j] - mid;
      if (ang < 0)
        ang += 2 * M_PI;
      else if (ang > 2 * M_PI)
        ang -= 2 * M_PI;
      angs[j] = M_PI - ang;
      // fprintf(stderr, "angs[%d] = %g \n", j, angs[j]);
    }

    // normalise the angles so struts are the same if flipped
    if (fabs(angs[3]) > fabs(angs[2])) {
      f_flipped[i] = true;
      std::swap(angs[0], angs[1]);
      std::swap(angs[2], angs[3]);
      for (int a = 0; a < 4; a++)
        angs[a] *= -1;
    }

    double len1 = (v[1] - v[0]).len();            // first contact length
    double len2 = (v[3] - v[2]).len();            // second contact length
    double len0 = (strut_len) ? strut_len : len1; // strut length
    double dist1 = (len0 - len1) / 2;
    double dist2 = (len0 - len2) / 2;

    units[{len0, dist1, dist2, len0 - dist2, len0 - dist1, angs[0], angs[1],
           angs[2], angs[3]}]
        .push_back(i);
  }

  double radius = (rad_min + rad_max) / 4; // average from max and min diameters
  report += msg_str("radius:            %19.15f (+/-%.15f)\n", radius,
                    radius - rad_min / 2) +
            msg_str("strut length from: %19.15f\n", len_max);

  int unit_no = 0;
  const Color colors[] = {Color(0, 0, 0), Color(255, 255, 255)};
  for (const auto &unit : units) {
    report += msg_str("%-4d, %4d, %16.13f\n", unit_no, (int)unit.second.size(),
                      unit.first[0]);
    report += msg_str("    %16.13f,%16.13f,%16.13f,%16.13f\n", unit.first[1],
                      unit.first[2], unit.first[3], unit.first[4]);
    report += msg_str("    %16.8f,%16.8f,%16.8f,%16.8f\n",
                      rad2deg(unit.first[5]), rad2deg(unit.first[8]),
                      rad2deg(unit.first[7]), rad2deg(unit.first[6]));
    if (color) {
      for (auto f : unit.second) {
        geom.colors(FACES).set(f, Color(unit_no));
        const int v0 = geom.faces(f, 0);
        const int v1 = geom.faces(f, 1);
        bool is_dihedral = (double_eq(unit.first[3], -unit.first[4], eps) &&
                            double_eq(unit.first[5], -unit.first[6], eps));
        geom.colors(VERTS).set(v0, colors[!f_flipped[f] | is_dihedral]);
        geom.colors(VERTS).set(v1, colors[f_flipped[f] | is_dihedral]);
      }
    }
    unit_no++;
  }

  return report;
}

vector<vector<int>> make_nexorade(Geometry &base_geom, const Symmetry &sym,
                                  double end_fraction, IterationControl it_ctrl,
                                  double factor)
{
  // Read and write to same model (Defered update is slower)
  SymmetricUpdater sym_updater(base_geom, sym);
  const Geometry &geom = sym_updater.get_geom_working();
  const vector<Vec3d> &verts = geom.verts();
  const vector<vector<int>> &faces = geom.faces();

  vector<vector<int>> f2fs(faces.size(), vector<int>(4));
  {
    vector<vector<int>> v2fs(verts.size(), vector<int>(4));
    for (unsigned int i = 0; i < faces.size(); i++) {
      // v_idx -> strut_face, is_first_end_con, string_face, is_first_end_con
      for (int j = 0; j < 4; j++) {
        bool is_string = (j > 1); // string or strut
        v2fs[faces[i][j]][2 * is_string] = i;
        v2fs[faces[i][j]][2 * is_string + 1] = (j == 0 || j == 3); // con end
      }
    }
    for (unsigned int i = 0; i < v2fs.size(); i++)
      // fprintf(stderr, "v2fs[%d] = %d,%d,  %d,%d\n", i, v2fs[i][0],
      // v2fs[i][1], v2fs[i][2], v2fs[i][3]);
      for (unsigned int i = 0; i < faces.size(); i++) {
        f2fs[i][0] = v2fs[faces[i][0]][2]; // v0: strut_face -> string_face
        f2fs[i][1] = v2fs[faces[i][0]][3]; // v0: is_first_end_con
        f2fs[i][2] = v2fs[faces[i][1]][2]; // v1: strut_face -> string_face
        f2fs[i][3] = v2fs[faces[i][1]][3]; // v1: is_first_end_con
      }
  }

  const auto principal_faces = sym_updater.get_principal(FACES);

  double rad = -1;
  double max_diff = 0.0;
  double test_val = it_ctrl.get_test_val();

  for (it_ctrl.start_iter_with_setup(); !it_ctrl.is_done();
       it_ctrl.next_iter()) {
    max_diff = 0.0;
    double dist_sum = 0;
    double rad_diff_sum = 0;
    for (int f_idx : principal_faces) {
      const auto &face = faces[f_idx];
      // Start: make sure vertices are current
      for (int i = 0; i < 4; i++)
        sym_updater.update_from_principal_vertex(face[i]);

      for (int i = 0; i < 2; i++) {
        // Edge direction flips with i
        const vector<int> edge = {faces[f_idx][i], faces[f_idx][!i]};

        const int f = f2fs[f_idx][i * 2];
        bool first_end = f2fs[f_idx][i * 2 + 1];
        vector<int> e_con = {faces[f][!first_end], faces[f][first_end]};

        for (auto idx : {edge[0], edge[1], e_con[0], e_con[1]})
          sym_updater.update_from_principal_vertex(idx);

        auto &v_edge0 = verts[edge[0]];
        auto &v_edge1 = verts[edge[1]];
        auto &v_e_con0 = verts[e_con[0]];
        auto &v_e_con1 = verts[e_con[1]];
        Vec3d P; // nearest point on this edge
        Vec3d Q; // nearest point on connected edge
        lines_nearest_points(v_edge0, v_edge1, v_e_con0, v_e_con1, P, Q);
        auto perp = P - Q; // from connected edge to this edge
        dist_sum += perp.len();
        rad_diff_sum += fabs(perp.len() / 2 - rad);

        const auto perp_dir = perp.unit();
        const auto base_point = v_e_con0 + (v_e_con1 - v_e_con0) * end_fraction;
        const auto ideal_point = base_point + perp_dir * (2 * rad);
        const auto diff_vec = ideal_point - v_edge0;
        const auto diff = diff_vec.len();
        if (diff > max_diff)
          max_diff = diff;

        auto v_new = (v_edge0 + diff_vec * factor);
        // adjust vertex (skip setup iter as radius value not set)
        if (!it_ctrl.is_setup_iter())
          sym_updater.update_principal_vertex(edge[0], v_new);
      }
    }
    // Sum of per distances, /2 to make radius, /num_f_orbits for
    // times through loop, /2 for number of ends
    const double target_rad = (dist_sum / 2) / principal_faces.size() / 2;

    if (!it_ctrl.is_setup_iter())
      rad += (target_rad - rad) * factor;
    else // first time through
      rad = target_rad;

    // Do not check status or finish before modifying the model
    if (!it_ctrl.is_setup_iter()) {
      string finish_msg;
      if (it_ctrl.is_status_check_iter()) {

        if (max_diff < test_val) {
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

        it_ctrl.print("%-12u max_diff:%17.15e\n", it_ctrl.get_current_iter(),
                      max_diff);
      }
    }
  }

  base_geom = sym_updater.get_geom_final();

  return f2fs;
}

void to_output_type(Geometry &geom, int out_type, double strut_len)
{
  if (out_type == 'f') // full
    return;

  geom.clear(EDGES);

  for (unsigned int f = 0; f < geom.faces().size(); f++) {
    const auto &face = geom.faces(f);
    const auto &col = geom.colors(FACES).get(f);
    if (out_type == 'r') { // rotegrity
      for (int i = 1; i < 4; i++) {
        geom.add_edge_raw(make_edge(face[i], face[(i + 1) % 4]), col);
        if (i < 3)
          geom.colors(VERTS).set(face[i - 1], col); // face[0] and face[1]
      }
    }
    else if (out_type == 'n' || out_type == 'N') { // nexorade
      auto edge = make_edge(face[0], face[1]);
      geom.add_edge_raw(edge, col);
      if (strut_len) {
        auto mid = geom.edge_cent(edge);
        auto dir = geom.edge_vec(edge).unit();
        geom.verts(edge[0]) = mid - dir * strut_len / 2;
        geom.verts(edge[1]) = mid + dir * strut_len / 2;
      }
      for (int i = 0; i < 2; i++) {
        if (out_type == 'n' || !geom.colors(VERTS).get(edge[0]).is_set()) {
          geom.colors(VERTS).set(edge[0], col);
          geom.colors(VERTS).set(edge[1], col);
        }
      }
    }
    else if (out_type == 't') { // tensegrity
      Color colors[] = {col, Color(1.0, 1.0, 1.0)};
      for (int i = 0; i < 4; i++) {
        auto edge = make_edge(face[i], face[(i + 1) % 4]);
        geom.add_edge_raw(edge, colors[(i > 0)]);
        geom.colors(VERTS).set(face[i], colors[(i < 2)]);
      }
    }
  }
  geom.clear(FACES);
}

int main(int argc, char *argv[])
{
  rot_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);

  if (!geom.faces().size())
    opts.error("input file contains no faces");

  geom.transform(Trans3d::translate(-geom.centroid()));
  Symmetry sym(geom);
  Geometry o_geom;
  if (opts.already_twisted) {
    opts.print_status_or_exit(check_twist_model(geom));
    o_geom = geom;
  }
  else {
    GeometryInfo info(geom);
    opts.print_status_or_exit(check_polyhedron_model(info));
    if (!info.is_oriented()) {
      opts.warning("polyhedron model: orienting model as not oriented");
      geom.orient(1);
    }
    make_twist(o_geom, geom, sym, opts.method, opts.end_fraction,
               opts.col_type);
    project_onto_sphere(o_geom);
  }

  string report;
  if (opts.algorithm == 'r') {
    make_rotegrity(o_geom, sym.get_max_direct_sub_sym(), opts.end_fraction,
                   opts.it_ctrl, opts.adjust_fact / 100);
    report = report_rotegrity(o_geom, opts.end_fraction, opts.col_type == 'u');
  }
  else if (opts.algorithm == 'n') {
    auto f2fs =
        make_nexorade(o_geom, sym.get_max_direct_sub_sym(), opts.end_fraction,
                      opts.it_ctrl, opts.adjust_fact / 100);
    report =
        report_nexorade(o_geom, opts.strut_len, opts.col_type == 'u', f2fs);
  }

  to_output_type(o_geom, opts.output_type, opts.strut_len);
  for (int i = 0; i < 3; i++)
    opts.clrngs[i].set_geom(&o_geom);
  opts.clrngs[VERTS].v_apply_cmap();
  opts.clrngs[EDGES].e_apply_cmap();
  opts.clrngs[FACES].f_apply_cmap();

  if (report.size())
    fprintf(stderr, "\nResult\n------\n%s", report.c_str());
  opts.write_or_error(o_geom, opts.ofile);

  return 0;
}
