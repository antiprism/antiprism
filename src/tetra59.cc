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
                  Bjorn Poonen, and Mechael Rubinstein
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

// terms given in paper as (a12,a34,a13,a24,a14,a23)
// values of terms 5 and 6 from the paper are reversed
struct Tetra59Item {
  string comment;
  int N;
  int term1;
  int term2;
  int term3;
  int term4;
  int term6; // reversed
  int term5; // reversed
  int regge;
};

// clang-format off
Tetra59Item tetra59_item_list[] = {
  {"H2",  12,   3,  4,   3,  4,   6,  8,   1},
  {"",    24,   5,  9,   6,  8,  13, 15,   1},
  {"T0",  12,   3,  6,   4,  6,   4,  6,   2},
  {"",    24,   7, 11,   7, 13,   8, 12,   2},
  {"T18", 15,   3,  3,   3,  5,  10, 10,   3},
  {"",    15,   2,  4,   4,  4,  10, 10,   3},
  {"",    15,   3,  3,   4,  4,   9, 11,   3},
  {"T7",  15,   3,  3,   5,  5,   9,  9,   4},
  {"T23", 15,   5,  5,   5,  9,   6,  6,   5},
  {"",    15,   3,  7,   6,  6,   7,  7,   5},
  {"",    15,   4,  8,   5,  5,   7,  7,   5},
  {"",    21,   3,  9,   7,  7,  12, 12,   6},
  {"",    21,   4, 10,   6,  6,  12, 12,   6},
  {"",    21,   6,  6,   7,  7,   9, 15,   6},
  {"T17", 30,   6, 12,  10, 15,  10, 20,   7},
  {"",    30,   4, 14,  10, 15,  12, 18,   7},
  {"",    60,   8, 28,  19, 31,  25, 35,   7},
  {"",    60,  12, 24,  15, 35,  25, 35,   7},
  {"",    60,  13, 23,  15, 35,  24, 36,   7},
  {"",    60,  13, 23,  19, 31,  20, 40,   7},
  {"T13", 30,   6, 18,  10, 10,  15, 15,   8},
  {"",    30,   4, 16,  12, 12,  15, 15,   8},
  {"",    30,   9, 21,  10, 10,  12, 12,   8},
  {"T16", 30,   6,  6,  10, 12,  15, 20,   9},
  {"",    30,   5,  7,  11, 11,  15, 20,   9},
  {"",    60,   7, 17,  20, 24,  35, 35,   9},
  {"",    60,   7, 17,  22, 22,  33, 37,   9},
  {"",    60,  10, 14,  17, 27,  35, 35,   9},
  {"",    60,  12, 12,  17, 27,  33, 37,   9},
  {"T21", 30,   6, 10,  10, 15,  12, 18,  10},
  {"",    30,   5, 11,  10, 15,  13, 17,  10},
  {"",    60,  10, 22,  21, 29,  25, 35,  10},
  {"",    60,  11, 21,  19, 31,  26, 34,  10},
  {"",    60,  11, 21,  21, 29,  24, 36,  10},
  {"",    60,  12, 20,  19, 31,  25, 35,  10},
  {"T6",  30,   6, 10,   6, 10,  15, 24,  11},
  {"",    60,   7, 25,  12, 20,  35, 43,  11},
  {"T2",  30,   6, 12,   6, 12,  15, 20,  12},
  {"",    60,  12, 24,  13, 23,  29, 41,  12},
  {"T3",  30,   6, 12,  10, 10,  15, 18,  13},
  {"",    30,   7, 13,   9,  9,  15, 18,  13},
  {"",    60,  12, 24,  17, 23,  33, 33,  13},
  {"",    60,  14, 26,  15, 21,  33, 33,  13},
  {"",    60,  15, 21,  20, 20,  27, 39,  13},
  {"",    60,  17, 23,  18, 18,  27, 39,  13},
  {"T4",  30,   6, 15,   6, 18,  10, 20,  14},
  {"",    30,   6, 15,   7, 17,   9, 21,  14},
  {"",    60,   9, 33,  14, 34,  21, 39,  14},
  {"",    60,   9, 33,  15, 33,  20, 40,  14},
  {"",    60,  11, 31,  12, 36,  21, 39,  14},
  {"",    60,  11, 31,  15, 33,  18, 42,  14},
  {"T1",  30,   6, 15,  10, 15,  12, 15,  15},
  {"",    30,   6, 15,  11, 14,  11, 16,  15},
  {"",    30,   8, 13,   8, 17,  12, 15,  15},
  {"",    30,   8, 13,   9, 18,  11, 14,  15},
  {"",    30,   8, 17,   9, 12,  11, 16,  15},
  {"",    30,   9, 12,   9, 18,  10, 15,  15},
  {"T5",  30,  10, 12,  10, 12,  15, 12,  16},
  {"",    60,  19, 25,  20, 24,  29, 25,  16},
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
  int reg(int i);
  int term(int term_no, int i);
};

tetra59::tetra59()
{
  tetra59_items = tetra59_item_list;
  last_tetra59 = sizeof(tetra59_item_list) / sizeof(tetra59_item_list[0]);
}

void tetra59::list_poly(int idx, FILE *fp)
{
  fprintf(fp, "%2d) N=%2d (%2d,%3d), (%2d,%3d), (%2d,%3d) %-20s\n", idx + 1,
          N(idx), term(0, idx), term(1, idx), term(2, idx), term(3, idx),
          term(4, idx), term(5, idx), tetra59_items[idx].comment.c_str());
}

void tetra59::list_polys(FILE *fp)
{
  int reg_group = 0;
  for (int i = 0; i < last_tetra59; i++) {
    if (reg_group != reg(i)) {
      reg_group = reg(i);
      fprintf(stderr, "Regge group %d\n", reg_group);
    }
    list_poly(i, fp);
  }
  fprintf(fp, "\n");
  fprintf(fp, "Entries are given as six dihedral angles of faces 1,2,3,4\n");
  fprintf(fp, "(a12,a34,a13,a24,a14,a23) as multiples of pi/N\n");
  fprintf(fp, "An extra symbol indicates 15 previously known forms from\n");
  fprintf(fp, "V. G. Boltianskii, Hilbert’s third problem, 1978\n");
  fprintf(fp, "The cases without a symbol are newly discovered forms\n");
  fprintf(fp, "No cases of N=21 (Regge group 6) were previously known\n");
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
int tetra59::reg(int i) { return tetra59_items[i].regge; }

int tetra59::term(int term_no, int i)
{
  int t = 0;

  if (term_no == 0)
    t = tetra59_items[i].term1;
  else if (term_no == 1)
    t = tetra59_items[i].term2;
  else if (term_no == 2)
    t = tetra59_items[i].term3;
  else if (term_no == 3)
    t = tetra59_items[i].term4;
  else if (term_no == 4)
    t = tetra59_items[i].term5;
  else if (term_no == 5)
    t = tetra59_items[i].term6;

  return t;
}

class tetra59_opts : public ProgramOpts {
public:
  string ifile;
  string ofile;

  string poly;               // polyhedron number 1 to 59
  int s = 0;                 // for special cases 1 and 2
  double angle = 45;         // angle for special cases (default 30 degrees)
  bool allow_angles = false; // -w switch allows all angles for case s
  bool reflect = false;      // make reflection model
  bool scale_volume = false; // scale volume to 1
  char dih_order = 'a';      // order of dihedral pairs
  bool verbose = false;      // output edge math
  bool verbose_face = false; // output face math

  bool list_polys = false; // output the list of models to the screen

  char coloring_method = 'u'; // color method for color by symmetry
  int face_opacity = -1;      // tranparency from 0 to 255

  Color vert_col = Color(255, 215, 0);   // gold
  Color edge_col = Color(211, 211, 211); // lightgray

  ColorMapMulti map;

  tetra59_opts() : ProgramOpts("tetra59") {}

  void process_command_line(int argc, char **argv);
  void usage();
};

void extended_help()
{
  fprintf(stdout, R"(
The project was undertaken in memory of John H. Conway

Abstract. We classify all sets of nonzero vectors in R3 such that the angle
formed by each pair is a rational multiple of Pi .The special case of
four-element subsets lets us classify all tetrahedra whose dihedral angles are
multiples of Pi, solving a 1976 problem of Conway and Jones: there are 2
one-parameter families and 59 sporadic tetrahedra, all but three of which are
related to either the icosidodecahedron or the B3 root lattice. The proof
requires the solution in roots of unity of a W(D6)-symmetric polynomial
equation with 105 monomials (the previous record was 12 monomials).

A brief description of a Regge Symmetry group, a mathematical symmetry.

For all the tetrahedra in a Regge Symmetry group, if given equal volume they
will have the following characteristics

1) For all the tetrahedrons in the group, the sum of the 6 edges will be equal
2) In the group, for a tetrahedron with edges (x,y,a,b,c,d) there may be one
   tetrahedrons with edges x,y,s-a,s-b,s-c,s-d where s = (a+b+c+d)/2
   There may be more than one pairing in a group. Of the pair...
3) The two opposing dihedral angles at edges x and y will be equal
4) The edge lengths of x and y of the two tetrahedra will be equal
5) If i,j,k,l are the dihedral angles at the edges a,b,c,d then the dihedral
   angles at edges s-a,s-b,s-c,s-d are t-i,t-j,t-k,t-l where t = (i+j+k+l)/2
   
For more information see: https://arxiv.org/abs/1903.04929

)");
}

void tetra59_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] polyhedron

Generate 59 Tetrahedra with Rational Dihedral Angles in off format. The 59
Sporadic Tetrahedra is from a paper by Kiran S. Kedlaya, Alexander Kolpakov,
Bjorn Poonen, and Mechael Rubinstein: Space Vectors Forming Rational Angles
The paper can be found at: https://arxiv.org/abs/2011.14232
There are also two special infinite cases included in the findings
The first case was published by M.J.M Hill in 1895. The second case is new.

Options
%s
  -H        abstract from the paper and description of regge symmetry
  -l        display the list of Sporadic Tetrahedra 1 thru 59
  -r        reflect
  -z        scale volume to 1
  -d <mthd> order of dihedral pairs, for matching dihedral angle at position
            term 1 for regge symmetry (also permutes special cases) (default:a)
            pairs are: 1-(a12,a34) 2-(a13,a24) 3-(a14,a23)
            a:1,2,3; b:1,3,2; c:2,1,3; d:2,3,1; e:3,1,2; f:3,2,1
  -v        verbose output of edge math
  -b        verbose output of face math
  -o <file> write output to file (default: write to standard output)

Special Cases
  -s <int>  special case 1 or 2 (default: none)
              case 1: (pi/2, pi/2, pi - 2x, pi/3, x, x)
                        for pi/6 < x < pi/2 (30 < x < 90 degrees)
              case 2: (5pi/6 - x, pi/6 + x, 2pi/3 - x, 2pi/3 - x, x, x)
                        for pi/6 < x <= pi/3 (30 < x <= 60 degrees)
  -a <ang>  angle in degrees (default: 45)
  -w        allow any angle for case s (for testing case 2, 60 < x < 90)

Coloring Options (run 'off_util -H color' for help on color formats)
  -V <col>  vertex color (default: gold)
  -E <col>  edge color   (default: lightgray)
  -f <mthd> mthd is face coloring method using color in map (default: u)
              keyword: none - sets no color
              u - unique coloring
              s - symmetric coloring
              r - regge group (1-16) (not for special cases)
  -T <tran> face transparency. valid range from 0 (invisible) to 255 (opaque)
  -m <maps> color maps for faces to be tried in turn
              (default for -f r: rainbow16, otherwise compound)

)",
          prog_name(), help_ver_text);
}

void tetra59_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  string map_file;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hHlvbrzd:s:a:wf:V:E:T:m:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'H':
      extended_help();
      exit(0);

    case 'l':
      list_polys = true;
      break;

    case 'v':
      verbose = true;
      break;

    case 'b':
      verbose_face = true;
      break;

    case 'r':
      reflect = true;
      break;

    case 'z':
      scale_volume = true;
      break;

    case 'd':
      if (strspn(optarg, "abcdef") != strlen(optarg) || strlen(optarg) > 1)
        error(msg_str("invalid order '%s'", optarg), c);
      else
        dih_order = *optarg;
      break;

    case 's':
      print_status_or_exit(read_int(optarg, &s), c);
      if (s < 1 || s > 2)
        error("must be 1 or 2", c);
      break;

    case 'a':
      print_status_or_exit(read_double(optarg, &angle), c);
      break;

    case 'w':
      allow_angles = true;
      warning("using -w allows any angle for case s");
      break;

    case 'f':
      if (!strcasecmp(optarg, "none"))
        coloring_method = '\0';
      else if (strspn(optarg, "rs") != strlen(optarg) || strlen(optarg) > 1)
        error(msg_str("invalid Coloring method '%s'", optarg), c);
      else
        coloring_method = *optarg;
      break;

    case 'V':
      print_status_or_exit(vert_col.read(optarg), c);
      break;

    case 'E':
      print_status_or_exit(edge_col.read(optarg), c);
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

  if (s && coloring_method == 'r') {
    warning("coloring by regge group has no effect for special cases", 'f');
    coloring_method = 'x';
  }

  if (!map_file.size()) {
    if (coloring_method == 'r')
      map_file = "rainbow16";
    else
      map_file = "compound";
  }
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
  return (cos(a) + cos(b) * cos(c)) / (sin(b) * sin(c));
}

// the dihedral angles are already calculated and are in radians
Geometry make_poly(double a12, double a34, double a13, double a24, double a14,
                   double a23, const tetra59_opts &opts)
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
  double angle243 = face_cos_a(a14, a24, a34);

  // face D angles (not needed to build)
  double angle134 = face_cos_a(a23, a13, a34);
  double angle143 = face_cos_a(a24, a14, a34);
  double angle314 = face_cos_a(a12, a13, a14);

  // clang-format off
  if (opts.verbose_face) {
    fprintf(stderr, "\n");
    fprintf(stderr, "P1 angle213 cos = % .17lf (%g deg)\n", angle213, cos2deg(angle213));
    fprintf(stderr, "P2 angle123 cos = % .17lf (%g deg)\n", angle123, cos2deg(angle123));
    fprintf(stderr, "P3 angle132 cos = % .17lf (%g deg)\n", angle132, cos2deg(angle132));
    double sum = cos2deg(angle213) + cos2deg(angle123) + cos2deg(angle132);
    fprintf(stderr, "face A angles sum = %g\n\n", sum);

    fprintf(stderr, "P4 angle142 cos = % .17lf (%g deg)\n", angle142, cos2deg(angle142));
    fprintf(stderr, "P1 angle214 cos = % .17lf (%g deg)\n", angle214, cos2deg(angle214));
    fprintf(stderr, "P2 angle124 cos = % .17lf (%g deg)\n", angle124, cos2deg(angle124));
    sum = cos2deg(angle142) + cos2deg(angle214) + cos2deg(angle124);
    fprintf(stderr, "face B angles sum = %g\n\n", sum);

    fprintf(stderr, "P2 angle324 cos = % .17lf (%g deg)\n", angle324, cos2deg(angle324));
    fprintf(stderr, "P3 angle234 cos = % .17lf (%g deg)\n", angle234, cos2deg(angle234));
    fprintf(stderr, "P4 angle243 cos = % .17lf (%g deg)\n", angle243, cos2deg(angle243));
    sum = cos2deg(angle324) + cos2deg(angle234) + cos2deg(angle243);
    fprintf(stderr, "face C angles sum = %g\n\n", sum);
      
    fprintf(stderr, "P3 angle134 cos = % .17lf (%g deg)\n", angle134, cos2deg(angle134));
    fprintf(stderr, "P4 angle143 cos = % .17lf (%g deg)\n", angle143, cos2deg(angle143));
    fprintf(stderr, "P1 angle314 cos = % .17lf (%g deg)\n", angle314, cos2deg(angle314));
    sum = cos2deg(angle134) + cos2deg(angle143) + cos2deg(angle314);
    fprintf(stderr, "face D angles sum = %g\n", sum);
  }
  // clang-format on

  if (opts.verbose) {
    fprintf(stderr, "\n");
    fprintf(stderr, "dihedral angles at:\n");
    fprintf(stderr, "edge a12 = %.17lf (%g deg)\n", a12, rad2deg(a12));
    fprintf(stderr, "edge a34 = %.17lf (%g deg)\n", a34, rad2deg(a34));
    fprintf(stderr, "edge a13 = %.17lf (%g deg)\n", a13, rad2deg(a13));
    fprintf(stderr, "edge a24 = %.17lf (%g deg)\n", a24, rad2deg(a24));
    fprintf(stderr, "edge a14 = %.17lf (%g deg)\n", a14, rad2deg(a14));
    fprintf(stderr, "edge a23 = %.17lf (%g deg)\n", a23, rad2deg(a23));
  }

  // edge lengths
  double edge12 = 1;
  double edge13 = (edge12 / cos2sin(angle132)) * cos2sin(angle123);
  double edge14 = (edge12 / cos2sin(angle142)) * cos2sin(angle124);
  double edge23 = (edge12 / cos2sin(angle132)) * cos2sin(angle213); // info only
  double edge24 = (edge12 / cos2sin(angle142)) * cos2sin(angle214); // info only
  double edge34 = (edge14 / cos2sin(angle134)) * cos2sin(angle314); // info only

  if (opts.verbose) {
    fprintf(stderr, "\n");
    fprintf(stderr, "edge12 = %.17lf\n", edge12);
    fprintf(stderr, "edge13 = %.17lf\n", edge13);
    fprintf(stderr, "edge14 = %.17lf\n", edge14);
    fprintf(stderr, "edge23 = %.17lf\n", edge23);
    fprintf(stderr, "edge24 = %.17lf\n", edge24);
    fprintf(stderr, "edge34 = %.17lf\n", edge34);
  }

  Geometry geom;

  double mult = (opts.reflect ? -1 : 1);

  // add base unit edge
  geom.add_vert(Vec3d(0, 0, 0)); // P1
  geom.add_vert(Vec3d(0, 1, 0)); // P2

  // rotate P3 on xy plane
  Geometry pgeom;
  Vec3d P = Vec3d(mult * edge13, 0, 0);
  pgeom.add_vert(P);
  pgeom.transform(
      Trans3d::rotate(0, 0, mult * (M_PI / 2) - acos(safe_for_trig(angle213))));
  geom.append(pgeom);

  // rotate P4 on xy plane
  pgeom.clear_all();
  P = Vec3d(mult * edge14, 0, 0);
  pgeom.add_vert(P);
  pgeom.transform(
      Trans3d::rotate(0, 0, mult * (M_PI / 2) - acos(safe_for_trig(angle214))));
  // rotate face B around Y using dihedral angle
  pgeom.transform(Trans3d::rotate(0, mult * a12, 0));
  geom.append(pgeom);

  // add oriented faces, color coded
  geom.add_face({2, 1, 0}); // A yellow
  geom.add_face({0, 1, 3}); // B red
  geom.add_face({1, 2, 3}); // C green
  geom.add_face({3, 2, 0}); // D blue

  return geom;
}

// integer format from the sporadic cases list
Geometry make_poly_sporadic(int N, int A12, int A34, int A13, int A24, int A14,
                            int A23, int regge, const tetra59_opts &opts)
{
  // dihedral angles in radians
  double a12 = A12 * M_PI / N;
  double a34 = A34 * M_PI / N;
  double a13 = A13 * M_PI / N;
  double a24 = A24 * M_PI / N;
  double a14 = A14 * M_PI / N;
  double a23 = A23 * M_PI / N;

  // output terms as reduced fractions of pi
  if (opts.verbose || opts.verbose_face) {
    fprintf(stderr, "\nregge group = %d\n", regge);

    int a12_gcd = (int)gcd(A12, N);
    int a34_gcd = (int)gcd(A34, N);
    int a13_gcd = (int)gcd(A13, N);
    int a24_gcd = (int)gcd(A24, N);
    int a14_gcd = (int)gcd(A14, N);
    int a23_gcd = (int)gcd(A23, N);

    string a12_num = (A12 / a12_gcd == 1) ? "" : std::to_string(A12 / a12_gcd);
    string a34_num = (A34 / a34_gcd == 1) ? "" : std::to_string(A34 / a34_gcd);
    string a13_num = (A13 / a13_gcd == 1) ? "" : std::to_string(A13 / a13_gcd);
    string a24_num = (A24 / a24_gcd == 1) ? "" : std::to_string(A24 / a24_gcd);
    string a14_num = (A14 / a14_gcd == 1) ? "" : std::to_string(A14 / a14_gcd);
    string a23_num = (A23 / a23_gcd == 1) ? "" : std::to_string(A23 / a23_gcd);

    fprintf(stderr, "\n");
    fprintf(stderr, "%spi/%d, %spi/%d, %spi/%d, %spi/%d, %spi/%d, %spi/%d\n",
            a12_num.c_str(), N / a12_gcd, a34_num.c_str(), N / a34_gcd,
            a13_num.c_str(), N / a13_gcd, a24_num.c_str(), N / a24_gcd,
            a14_num.c_str(), N / a14_gcd, a23_num.c_str(), N / a23_gcd);
  }

  return make_poly(a12, a34, a13, a24, a14, a23, opts);
}

// angle is in radians
Geometry case1(const vector<int> &dihedral_order, const tetra59_opts &opts)
{
  // (pi/2, pi/2, pi − 2x, pi/3, x, x)
  // for pi/6 < x < pi/2 (30 < x < 90 degrees)

  if (opts.verbose || opts.verbose_face)
    fprintf(stderr, "\ncase 1: angle = %g\n", rad2deg(opts.angle));

  vector<double> term(6);
  term[0] = M_PI / 2;
  term[1] = M_PI / 2;
  term[2] = M_PI - 2 * opts.angle;
  term[3] = M_PI / 3;
  term[4] = opts.angle;
  term[5] = opts.angle;

  vector<double> terms(6);
  for (unsigned int i = 0; i < 3; i++) {
    terms[i * 2] = term[dihedral_order[i] * 2];
    terms[i * 2 + 1] = term[dihedral_order[i] * 2 + 1];
  }

  return make_poly(terms[0], terms[1], terms[2], terms[3], terms[4], terms[5],
                   opts);
}

// angle is in radians
Geometry case2(const vector<int> &dihedral_order, const tetra59_opts &opts)
{
  // (5pi/6 − x, pi/6 + x, 2pi/3 − x, 2pi/3 − x, x, x)
  // for pi/6 < x ≤ pi/3 (30 < x <= 60 degrees)

  if (opts.verbose || opts.verbose_face)
    fprintf(stderr, "\ncase 2: angle = %g\n", rad2deg(opts.angle));

  vector<double> term(6);
  term[0] = 5 * M_PI / 6 - opts.angle;
  term[1] = M_PI / 6 + opts.angle;
  term[2] = 2 * M_PI / 3 - opts.angle;
  term[3] = 2 * M_PI / 3 - opts.angle;
  term[4] = opts.angle;
  term[5] = opts.angle;

  vector<double> terms(6);
  for (unsigned int i = 0; i < 3; i++) {
    terms[i * 2] = term[dihedral_order[i] * 2];
    terms[i * 2 + 1] = term[dihedral_order[i] * 2 + 1];
  }

  return make_poly(terms[0], terms[1], terms[2], terms[3], terms[4], terms[5],
                   opts);
}

void face_coloring(Geometry &geom, int regge_grp, const tetra59_opts &opts)
{
  // color by sub-symmetry as map indexes happened by default in sym_repeat()
  if (!opts.coloring_method) {
    // no color, strip colors
    geom.colors(FACES).clear();
    geom.clear(EDGES);
    geom.colors(VERTS).clear();
  }
  else {
    Coloring clrng;
    clrng.add_cmap(opts.map.clone());
    clrng.set_geom(&geom);

    if (opts.coloring_method == 'u') {
      clrng.f_unique(true);
    }
    else if (opts.coloring_method == 's') {
      Symmetry sym;
      vector<vector<set<int>>> sym_equivs;
      sym.init(geom, &sym_equivs);
      clrng.f_sets(sym_equivs[2], true);
    }
    else if (opts.coloring_method == 'r')
      clrng.f_one_col(opts.map.get_col(regge_grp - 1));

    // color vertices
    clrng.v_one_col(opts.vert_col);

    // color edges
    geom.add_missing_impl_edges();
    clrng.e_one_col(opts.edge_col);

    // transparency
    if (opts.face_opacity > -1) {
      ColorValuesToRangeHsva valmap(
          msg_str("A%g", (double)opts.face_opacity / 255));
      valmap.apply(geom, FACES);

      for (const auto &kp : geom.colors(FACES).get_properties()) {
        if (kp.second.is_index()) {
          opts.warning("map indexes cannot be made transparent", 'T');
          break;
        }
      }
    }
  }

  // check if some faces are not set for transparency warning
  if (opts.face_opacity > -1) {
    if (geom.colors(FACES).get_properties().size() < geom.faces().size())
      opts.warning("unset faces cannot be made transparent", 'T');
  }
}

void scale_volume_to_one(Geometry &geom, const tetra59_opts &opts)
{
  GeometryInfo info(geom);
  double scale = pow(fabs(info.volume()), 1.0 / 3.0);
  geom.transform(Trans3d::scale(1 / scale));

  if (opts.verbose) {
    fprintf(stderr, "\nvolume set to 1. edges have been resized\n");
    double sum = 0;
    for (unsigned int i = 0; i < 3; i++) {
      for (unsigned int j = i + 1; j < 4; j++) {
        double len = (geom.verts(i) - geom.verts(j)).len();
        fprintf(stderr, "edge%1d%1d = %.17lf\n", i + 1, j + 1, len);
        sum += len;
      }
    }
    fprintf(stderr, "sum   = %.17lf\n", sum);
  }
}

int main(int argc, char *argv[])
{
  tetra59_opts opts;
  opts.process_command_line(argc, argv);
  tetra59 tetra59s;

  if (opts.list_polys) {
    tetra59s.list_polys();
    exit(0);
  }

  // the list cases
  Geometry geom;
  int regge_grp = 0;

  // a:1,2,3; b:1,3,2; c:2,1,3; d:2,3,1; e:3,1,2; f:3,2,1
  vector<int> dihedral_order(3);
  if (strchr("a", opts.dih_order))
    dihedral_order = {0, 1, 2};
  else if (strchr("b", opts.dih_order))
    dihedral_order = {0, 2, 1};
  else if (strchr("c", opts.dih_order))
    dihedral_order = {1, 0, 2};
  else if (strchr("d", opts.dih_order))
    dihedral_order = {1, 2, 0};
  else if (strchr("e", opts.dih_order))
    dihedral_order = {2, 0, 1};
  else if (strchr("f", opts.dih_order))
    dihedral_order = {2, 1, 0};

  if (opts.s == 0) {
    int sym_no = tetra59s.lookup_sym_no(opts.poly);
    if (sym_no >= tetra59s.get_last_tetra59())
      opts.error("polyhedron number '" + opts.poly + "' out of range");
    if (sym_no < 0)
      opts.error("unknown polyhedron '" + opts.poly + "'");

    // default 'a': 1-(a12,a34) 2-(a13,a24) 3-(a14,a23)
    vector<int> terms(6);
    for (unsigned int i = 0; i < 3; i++) {
      terms[i * 2] = tetra59s.term(dihedral_order[i] * 2, sym_no);
      terms[i * 2 + 1] = tetra59s.term(dihedral_order[i] * 2 + 1, sym_no);
    }

    fprintf(stderr, "%2d) N=%2d (%2d,%3d), (%2d,%3d), (%2d,%3d)\n", sym_no + 1,
            tetra59s.N(sym_no), terms[0], terms[1], terms[2], terms[3],
            terms[4], terms[5]);

    regge_grp = tetra59s.reg(sym_no);
    geom = make_poly_sporadic(tetra59s.N(sym_no), terms[0], terms[1], terms[2],
                              terms[3], terms[4], terms[5], regge_grp, opts);
  }
  // special cases
  else if (opts.s == 1)
    geom = case1(dihedral_order, opts);
  else if (opts.s == 2)
    geom = case2(dihedral_order, opts);

  if (opts.scale_volume)
    scale_volume_to_one(geom, opts);

  if (opts.verbose || opts.verbose_face) {
    Symmetry sym(geom);
    fprintf(stderr, "\n");
    fprintf(stderr, "the symmetry of the tetrahedron is %s\n",
            sym.get_symbol().c_str());
  }

  face_coloring(geom, regge_grp, opts);

  opts.write_or_error(geom, opts.ofile);

  return 0;
}
