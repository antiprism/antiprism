/*
   Copyright (c) 2021, Roger Kaufman

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
   Name: tetra59.cc
   Description: 59 Tetrahedra with Rational Dihedral Angles
                Based on a paper by Kiran S. Kedlaya, Alexander Kolpakov,
                  Bjorn Poonen, AND Mechael Rubinstein
                Space Vectors Forming Rational Angles
                  https://arxiv.org/abs/2011.14232
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <set>
#include <string>
#include <vector>

using std::set;
using std::string;
using std::vector;

using namespace anti;

struct Tetra59Item {
  int N;
  int a12;
  int a34;
  int a13;
  int a24;
  int a14;
  int a23;
  const char *comment;
};

// clang-format off
Tetra59Item tetra59_item_list[] = {
  {12,  3,  4,  3,  4,  6,  8, "H2 (pi/4)"},
  {24,  5,  9,  6,  8, 13, 15, ""},
  {12,  3,  6,  4,  6,  4,  6, "T0"},
  {24,  7, 11,  7, 13,  8, 12, ""},
  {15,  3,  3,  3,  5, 10, 10, "T18"},
  {15,  2,  4,  4,  4, 10, 10, ""},
  {15,  3,  3,  4,  4,  9, 11, ""},
  {15,  3,  3,  5,  5,  9,  9, "T7"},
  {15,  5,  5,  5,  9,  6,  6, "T23"},
  {15,  3,  7,  6,  6,  7,  7, ""},
  {15,  4,  8,  5,  5,  7,  7, ""},
  {21,  3,  9,  7,  7, 12, 12, ""},
  {21,  4, 10,  6,  6, 12, 12, ""},
  {21,  6,  6,  7,  7,  9, 15, ""},
  {30,  6, 12, 10, 15, 10, 20, "T17"},
  {30,  4, 14, 10, 15, 12, 18, ""},
  {60,  8, 28, 19, 31, 25, 35, ""},
  {60, 12, 24, 15, 35, 25, 35, ""},
  {60, 13, 23, 15, 35, 24, 36, ""},
  {60, 13, 23, 19, 31, 20, 40, ""},
  {30,  6, 18, 10, 10, 15, 15, "T13"},
  {30,  4, 16, 12, 12, 15, 15, ""},
  {30,  9, 21, 10, 10, 12, 12, ""},
  {30,  6,  6, 10, 12, 15, 20, "T16"},
  {30,  5,  7, 11, 11, 15, 20, ""},
  {60,  7, 17, 20, 24, 35, 35, ""},
  {60,  7, 17, 22, 22, 33, 37, ""},
  {60, 10, 14, 17, 27, 35, 35, ""},
  {60, 12, 12, 17, 27, 33, 37, ""},
  {30,  6, 10, 10, 15, 12, 18, "T21"},
  {30,  5, 11, 10, 15, 13, 17, ""},
  {60, 10, 22, 21, 29, 25, 35, ""},
  {60, 11, 21, 19, 31, 26, 34, ""},
  {60, 11, 21, 21, 29, 24, 36, ""},
  {60, 12, 20, 19, 31, 25, 35, ""},
  {30,  6, 10,  6, 10, 15, 24, "T6"},
  {60,  7, 25, 12, 20, 35, 43, ""},
  {30,  6, 12,  6, 12, 15, 20, "T2"},
  {60, 12, 24, 13, 23, 29, 41, ""},
  {30,  6, 12, 10, 10, 15, 18, "T3"},
  {30,  7, 13,  9,  9, 15, 18, ""},
  {60, 12, 24, 17, 23, 33, 33, ""},
  {60, 14, 26, 15, 21, 33, 33, ""},
  {60, 15, 21, 20, 20, 27, 39, ""},
  {60, 17, 23, 18, 18, 27, 39, ""},
  {30,  6, 15,  6, 18, 10, 20, "T4"},
  {30,  6, 15,  7, 17,  9, 21, ""},
  {60,  9, 33, 14, 34, 21, 39, ""},
  {60,  9, 33, 15, 33, 20, 40, ""},
  {60, 11, 31, 12, 36, 21, 39, ""},
  {60, 11, 31, 15, 33, 18, 42, ""},
  {30,  6, 15, 10, 15, 12, 15, "T1"},
  {30,  6, 15, 11, 14, 11, 16, ""},
  {30,  8, 13,  8, 17, 12, 15, ""},
  {30,  8, 13,  9, 18, 11, 14, ""},
  {30,  8, 17,  9, 12, 11, 16, ""},
  {30,  9, 12,  9, 18, 10, 15, ""},
  {30, 10, 12, 10, 12, 15, 12, "T5"},
  {60, 19, 25, 20, 24, 29, 25, ""},
};
// clang-format on

class tetra59 {
private:
  int last_tetra59;
  Tetra59Item *tetra59_items;

public:
  tetra59();
  void list_poly(int idx, FILE *fp = stderr);
  void list_polys(FILE *fp = stderr);
  int lookup_sym_no(string sym);
  int get_last_tetra59() { return last_tetra59; }

  int N(int i);
  int a12(int i);
  int a34(int i);
  int a13(int i);
  int a24(int i);
  int a14(int i);
  int a23(int i);
};

tetra59::tetra59()
{
  tetra59_items = tetra59_item_list;
  last_tetra59 = sizeof(tetra59_item_list) / sizeof(tetra59_item_list[0]);
}

void tetra59::list_poly(int idx, FILE *fp)
{
  fprintf(fp, "%2d) N=%2d (%3d,%3d,%3d,%3d,%3d,%3d) %-20s\n", idx + 1,
          tetra59_items[idx].N, tetra59_items[idx].a12, tetra59_items[idx].a34,
          tetra59_items[idx].a13, tetra59_items[idx].a24,
          tetra59_items[idx].a14, tetra59_items[idx].a23,
          tetra59_items[idx].comment);
}

void tetra59::list_polys(FILE *fp)
{
  for (int i = 0; i < last_tetra59; i++)
    list_poly(i, fp);
  fprintf(fp, "\n");
  fprintf(fp, "Entries are given as six dihedral angles of faces 1,2,3,4\n");
  fprintf(fp, "(a12,a34,a13,a24,a14,a23) as multiples of pi/N\n");
  fprintf(fp, "An extra symbol indicates 15 previously known forms from\n");
  fprintf(fp, "V. G. Boltianskii, Hilbert’s third problem, 1978\n");
  fprintf(fp, "No cases of N equal to 21 were previously known\n");
}

int tetra59::lookup_sym_no(string sym)
{
  // remove double spaces and spaces at beginning and end
  string sym_norm;
  bool ignore_if_space = true;
  for (char i : sym) {
    if (i == ' ') {
      if (ignore_if_space)
        continue;
      else
        ignore_if_space = true;
    }
    else
      ignore_if_space = false;
    sym_norm += i;
  }

  if (sym_norm[sym_norm.size() - 1] == ' ')
    sym_norm.resize(sym_norm.size() - 1);

  // remove spaces either side of a punctuation mark
  string sym_norm2;
  for (unsigned int i = 0; i < sym_norm.length(); i++) {
    if (sym_norm[i] == ' ' &&
        ((i > 0 && ispunct(sym_norm[i - 1])) ||
         (i < sym_norm.length() && ispunct(sym_norm[i + 1]))))
      continue;
    sym_norm2 += sym_norm[i];
  }

  // sym_norm2 is now normalised

  // is it blank
  if (sym_norm2 == "")
    return -1;

  // is it the list order number
  char *endptr;
  int idx = strtol(sym_norm2.c_str(), &endptr, 10);
  if (!*endptr) // all of string is an integer
    return idx - 1;

  return idx;
}

int tetra59::N(int i) { return tetra59_items[i].N; }
int tetra59::a12(int i) { return tetra59_items[i].a12; }
int tetra59::a34(int i) { return tetra59_items[i].a34; }
int tetra59::a13(int i) { return tetra59_items[i].a13; }
int tetra59::a24(int i) { return tetra59_items[i].a24; }
int tetra59::a14(int i) { return tetra59_items[i].a14; }
int tetra59::a23(int i) { return tetra59_items[i].a23; }

class tetra59_opts : public ProgramOpts {
public:
  string ifile;
  string ofile;

  string poly;               // polyhedron number 1 to 59
  string case_type;          // multiple alternate cases are possible
  bool verbose = false;      // output math to screen
  int s = 0;                 // for special cases 1 and 2
  double angle = 45;         // angle for special cases (default 30 degrees)
  bool allow_angles = false; // -w switch allows all angles for case s
  bool list_polys = false;   // output the list of models to the screen

  char coloring_method = 'c'; // color method for color by symmetry
  int face_opacity = -1;      // tranparency from 0 to 255

  ColorMapMulti map;

  tetra59_opts() : ProgramOpts("tetra59") {}

  void process_command_line(int argc, char **argv);
  void usage();
};

void tetra59_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] polyhedron (BETA)

Generate 59 Tetrahedra with Rational Dihedral Angles in off format. The 59
Sporadic Tetrahedra is from a paper by Kiran S. Kedlaya, Alexander Kolpakov,
Bjorn Poonen, and Mechael Rubinstein: Space Vectors Forming Rational Angles
The paper can be found at: https://arxiv.org/abs/2011.14232
There are also two special infinite cases included in the findings
The first case was published by M.J.M Hill in 1895. The second case is new.

Options
%s
  -l        display the list of Sporadic Tetrahedra 1 thru 59
  -v        verbose output
  -o <file> write output to file (default: write to standard output)

Special Cases
  -s <int>  special case 1 or 2 (default: none)
              case 1: (pi/2, pi/2, pi - 2x, pi/3, x, x)
                        for pi/6 < x < pi/2 (30 < x < 90 degrees)
              case 2: (5pi/6 - x, pi/6 + x, 2pi/3 - x, 2pi/3 - x, x, x)
                        for pi/6 < x ≤ pi/3 (30 < x <= 60 degrees)
  -a <ang>  angle in degrees (default: 45)
  -w        allow any angle for case s (for testing case 2, 60 < x < 90)

Coloring Options (run 'off_util -H color' for help on color formats)
  -f <mthd> mthd is face coloring method using color in map (default: c)
              keyword: none - sets no color
              c - unique coloring for each compound constituent
              s - symmetric coloring (should always be one color)
  -T <tran> face transparency. valid range from 0 (invisible) to 255 (opaque)
  -m <maps> color maps for all elements to be tried in turn (default: compound)

)",
          prog_name(), help_ver_text);
}

void tetra59_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  string map_file;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hlva:s:wf:T:m:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'l':
      list_polys = true;
      break;

    case 'v':
      verbose = true;
      break;

    case 'a':
      print_status_or_exit(read_double(optarg, &angle), c);
      break;

    case 's':
      print_status_or_exit(read_int(optarg, &s), c);
      if (s < 1 || s > 2)
        error("must be 1 or 2", c);
      break;

    // undocumented switch
    case 'w':
      allow_angles = true;
      warning("using -w allows any angle for case s");
      break;

    case 'f':
      if (!strcasecmp(optarg, "none"))
        coloring_method = '\0';
      else if (strspn(optarg, "cs") != strlen(optarg) || strlen(optarg) > 1)
        error(msg_str("invalid Coloring method '%s'", optarg), c);
      else
        coloring_method = *optarg;
      break;

    case 'T':
      print_status_or_exit(read_int(optarg, &face_opacity), c);
      if (face_opacity < 0 || face_opacity > 255) {
        error("face transparency must be between 0 and 255", c);
      }
      break;

    case 'm':
      map_file = optarg;
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (!allow_angles) {
    if (s == 1) {
      if (angle <= 30 || angle >= 90)
        error("angle for case 1 range is 30 < x < 90 degrees", 'a');
    }
    else if (s == 2) {
      if (angle <= 30 || angle > 60)
        error("angle for case 2 range is 30 < x <= 60 degrees", 'a');
    }
  }

  // convert angle to radians
  angle = deg2rad(angle);

  if (argc == optind && !list_polys && !s)
    error("no polyhedron specified", "polyhedron");

  if (argc - optind > 1)
    error("too many arguments");

  if (argc - optind == 1)
    poly = argv[optind];

  if (!map_file.size())
    map_file = "compound";
  print_status_or_exit(map.init(map_file.c_str()), 'm');
}

// convert cosine to degrees
double cos2deg(double a) { return rad2deg(acos(safe_for_trig(a))); }

// return the sin of an angle given the cosine of an angle
double cos2sin(double a)
{
  // return sqrt(1-a*a);
  return cos((M_PI / 2) - acos(safe_for_trig(a)));
}

// given 3 dihedral angles in radians intersecting, the cosine of the face angle
// opposed to dihedral angle a is returned
double face_cos_a(double a, double b, double c)
{
  double face_angle_a = (cos(a) + cos(b) * cos(c)) / (sin(b) * sin(c));

  // clang-format off
  if (face_angle_a > 1.0) {
    fprintf(stderr, "ERROR: face_cos_a: greater than 1: %.17lf\n", face_angle_a);
    fprintf(stderr, "info: input (radians): a = %g b = %g c = %g\n", a, b, c);
    fprintf(stderr, "info: input (degrees): a = %g b = %g c = %g\n", rad2deg(a), rad2deg(b), rad2deg(c));
    fprintf(stderr, "info: cos(a) = %g cos(b) = %g cos(c) = %g\n", cos(a), cos(b), cos(c));
    fprintf(stderr, "info: numerator: cos(a) + cos(b)*cos(c) = %g\n", cos(a) + cos(b) * cos(c));
    fprintf(stderr, "info: sin(b) = %g sin(c) = %g\n", sin(b), sin(c));
    fprintf(stderr, "info: denominator: sin(b)*sin(c) = %g\n\n", sin(b) * sin(c));
    exit(0);
  }
  // clang-format on

  return face_angle_a;
}

// the dihedral angles are already calculated and are in radians
Geometry make_poly(double a12, double a34, double a13, double a24, double a14,
                   double a23, bool verbose)
{
  // face A angles
  double angle213 = face_cos_a(a14, a12, a13);
  double angle123 = face_cos_a(a24, a12, a23);
  double angle132 = face_cos_a(a34, a13, a23);

  // face B angles
  double angle142 = face_cos_a(a34, a14, a24);
  double angle214 = face_cos_a(a13, a12, a14);
  double angle124 = face_cos_a(a23, a12, a24);

  // face C angles (not needed to build)
  double angle324 = face_cos_a(a12, a23, a24);
  double angle234 = face_cos_a(a13, a23, a34);
  double angle243 = face_cos_a(a23, a24, a34);

  // face D angles (not needed to build)
  double angle134 = face_cos_a(a23, a13, a34);
  double angle143 = face_cos_a(a24, a14, a34);
  double angle314 = face_cos_a(a12, a13, a14);

  // clang-format off
  if (verbose) {
    fprintf(stderr, "P1 angle213 cos = % .17lf (%.17lf deg)\n", angle213, cos2deg(angle213));
    fprintf(stderr, "P2 angle123 cos = % .17lf (%.17lf deg)\n", angle123, cos2deg(angle123));
    fprintf(stderr, "P3 angle132 cos = % .17lf (%.17lf deg)\n", angle132, cos2deg(angle132));
    double sum = cos2deg(angle213) + cos2deg(angle123) + cos2deg(angle132);
    fprintf(stderr, "face A angles sum = %g\n\n", sum);

    fprintf(stderr, "P4 angle142 cos = % .17lf (%.17lf deg)\n", angle142, cos2deg(angle142));
    fprintf(stderr, "P1 angle214 cos = % .17lf (%.17lf deg)\n", angle214, cos2deg(angle214));
    fprintf(stderr, "P2 angle124 cos = % .17lf (%.17lf deg)\n", angle124, cos2deg(angle124));
    sum = cos2deg(angle142) + cos2deg(angle214) + cos2deg(angle124);
    fprintf(stderr, "face B angles sum = %g\n\n", sum);

    fprintf(stderr, "P2 angle324 cos = % .17lf (%.17lf deg)\n", angle324, cos2deg(angle324));
    fprintf(stderr, "P3 angle234 cos = % .17lf (%.17lf deg)\n", angle234, cos2deg(angle234));
    fprintf(stderr, "P4 angle243 cos = % .17lf (%.17lf deg)\n", angle243, cos2deg(angle243));
    sum = cos2deg(angle324) + cos2deg(angle234) + cos2deg(angle243);
    fprintf(stderr, "face C angles sum = %g\n\n", sum);

    fprintf(stderr, "P3 angle134 cos = % .17lf (%.17lf deg)\n", angle134, cos2deg(angle134));
    fprintf(stderr, "P4 angle143 cos = % .17lf (%.17lf deg)\n", angle143, cos2deg(angle143));
    fprintf(stderr, "P1 angle314 cos = % .17lf (%.17lf deg)\n", angle314, cos2deg(angle314));
    sum = cos2deg(angle134) + cos2deg(angle143) + cos2deg(angle314);
    fprintf(stderr, "face D angles sum = %g\n\n", sum);
  }
  // clang-format on

  // edge lengths
  double edge12 = 1;

  double edge13 = (edge12 / cos2sin(angle132)) * cos2sin(angle123);
  double edge23 = (edge12 / cos2sin(angle132)) * cos2sin(angle213); // info only

  double edge14 = (edge12 / cos2sin(angle142)) * cos2sin(angle124);
  double edge24 = (edge12 / cos2sin(angle142)) * cos2sin(angle214); // info only

  double edge34 = (edge14 / cos2sin(angle134)) * cos2sin(angle314); // info only

  if (verbose) {
    fprintf(stderr, "edge12 = %.17lf\n", edge12);
    fprintf(stderr, "edge13 = %.17lf\n", edge13);
    fprintf(stderr, "edge23 = %.17lf\n", edge23);
    fprintf(stderr, "edge14 = %.17lf\n", edge14);
    fprintf(stderr, "edge24 = %.17lf\n", edge24);
    fprintf(stderr, "edge34 = %.17lf\n", edge34);
  }

  Geometry geom;

  // add base unit edge
  geom.add_vert(Vec3d(0, 0, 0)); // P1
  geom.add_vert(Vec3d(0, 1, 0)); // P2

  // rotate P3 on xy plane
  Geometry pgeom;
  Vec3d P = Vec3d(edge13, 0, 0);
  pgeom.add_vert(P);
  pgeom.transform(
      Trans3d::rotate(0, 0, (M_PI / 2) - acos(safe_for_trig(angle213))));
  geom.append(pgeom);

  // rotate P4 on xy plane
  pgeom.clear_all();
  P = Vec3d(edge14, 0, 0);
  pgeom.add_vert(P);
  pgeom.transform(
      Trans3d::rotate(0, 0, (M_PI / 2) - acos(safe_for_trig(angle214))));
  // rotate face C around Y using dihedral angle
  pgeom.transform(Trans3d::rotate(0, a12, 0));
  geom.append(pgeom);

  // add oriented faces, color coded
  geom.add_face({2, 1, 0}, Color(255, 0, 0, 255));   // A red
  geom.add_face({0, 1, 3}, Color(255, 255, 0, 255)); // B yellow
  geom.add_face({1, 2, 3}, Color(0, 255, 0, 255));   // C green
  geom.add_face({3, 2, 0}, Color(0, 0, 255, 255));   // D blue

  return geom;
}

// integer format from the sporadic cases list
Geometry make_poly_sporadic(int N, int A12, int A34, int A13, int A24, int A14,
                            int A23, bool verbose)
{
  // dihedral angles in radians
  double a12 = A12 * M_PI / N;
  double a34 = A34 * M_PI / N;
  double a13 = A13 * M_PI / N;
  double a24 = A24 * M_PI / N;
  double a14 = A14 * M_PI / N;
  double a23 = A23 * M_PI / N;

  // clang-format off
  if (verbose) {
    fprintf(stderr, "\n");
    fprintf(stderr, "dihedral angles at:\n");
    fprintf(stderr, "edge a12 = %2dpi/%d = %.17lf (%.17lf deg)\n", A12, N, a12, rad2deg(a12));
    fprintf(stderr, "edge a34 = %2dpi/%d = %.17lf (%.17lf deg)\n", A34, N, a34, rad2deg(a34));
    fprintf(stderr, "edge a13 = %2dpi/%d = %.17lf (%.17lf deg)\n", A13, N, a13, rad2deg(a13));
    fprintf(stderr, "edge a24 = %2dpi/%d = %.17lf (%.17lf deg)\n", A24, N, a24, rad2deg(a24));
    fprintf(stderr, "edge a14 = %2dpi/%d = %.17lf (%.17lf deg)\n", A14, N, a14, rad2deg(a14));
    fprintf(stderr, "edge a23 = %2dpi/%d = %.17lf (%.17lf deg)\n", A23, N, a23, rad2deg(a23));
    fprintf(stderr, "\n");
  }
  // clang-format on

  return make_poly(a12, a34, a13, a24, a14, a23, verbose);
}

// special cases format
Geometry make_poly_special(double a12, double a34, double a13, double a24,
                           double a14, double a23, bool verbose)
{
  if (verbose) {
    fprintf(stderr, "\n");
    fprintf(stderr, "dihedral angles at:\n");
    fprintf(stderr, "edge a12 = %.17lf (%.17lf deg)\n", a12, rad2deg(a12));
    fprintf(stderr, "edge a34 = %.17lf (%.17lf deg)\n", a34, rad2deg(a34));
    fprintf(stderr, "edge a13 = %.17lf (%.17lf deg)\n", a13, rad2deg(a13));
    fprintf(stderr, "edge a24 = %.17lf (%.17lf deg)\n", a24, rad2deg(a24));
    fprintf(stderr, "edge a14 = %.17lf (%.17lf deg)\n", a14, rad2deg(a14));
    fprintf(stderr, "edge a23 = %.17lf (%.17lf deg)\n", a23, rad2deg(a23));
    fprintf(stderr, "\n");
  }

  return make_poly(a12, a34, a13, a24, a14, a23, verbose);
}

// angle is in radians
Geometry case1(double angle, bool verbose)
{
  // (pi/2, pi/2, pi − 2x, pi/3, x, x)
  // for pi/6 < x < pi/2 (30 < x < 90 degrees)
  return make_poly_special(M_PI / 2, M_PI / 2, M_PI - 2 * angle, M_PI / 3,
                           angle, angle, verbose);
}

// angle is in radians
Geometry case2(double angle, bool verbose)
{
  // (5pi/6 − x, pi/6 + x, 2pi/3 − x, 2pi/3 − x, x, x)
  // for pi/6 < x ≤ pi/3 (30 < x <= 60 degrees)
  return make_poly_special(5 * M_PI / 6 - angle, M_PI / 6 + angle,
                           2 * M_PI / 3 - angle, 2 * M_PI / 3 - angle, angle,
                           angle, verbose);
}

int main(int argc, char *argv[])
{
  // test works for cube, returns 90 degrees
  // fprintf(stderr,"cube dihedral = %.17lf\n", cos2deg(face_cos_a(M_PI/2,
  // M_PI/2, M_PI/2))); exit(0);

  tetra59_opts opts;
  opts.process_command_line(argc, argv);
  tetra59 tetra59s;

  if (opts.list_polys) {
    tetra59s.list_polys();
    exit(0);
  }

  // the list cases
  Geometry geom;
  if (opts.s == 0) {
    int sym_no = tetra59s.lookup_sym_no(opts.poly);
    if (sym_no >= tetra59s.get_last_tetra59())
      opts.error("polyhedron number '" + opts.poly + "' out of range");
    if (sym_no < 0)
      opts.error("unknown polyhedron '" + opts.poly + "'");

    tetra59s.list_poly(sym_no);

    geom = make_poly_sporadic(tetra59s.N(sym_no), tetra59s.a12(sym_no),
                              tetra59s.a34(sym_no), tetra59s.a13(sym_no),
                              tetra59s.a24(sym_no), tetra59s.a14(sym_no),
                              tetra59s.a23(sym_no), opts.verbose);
  }
  // special cases
  else if (opts.s == 1)
    geom = case1(opts.angle, opts.verbose);
  else if (opts.s == 2)
    geom = case2(opts.angle, opts.verbose);

  // compound_coloring(geom, opts);

  opts.write_or_error(geom, opts.ofile);

  return 0;
}
