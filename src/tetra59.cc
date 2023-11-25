/*
   Copyright (c) 2021-2023, Roger Kaufman

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
#include "color_common.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

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
  void list_poly_int(int idx, FILE *fp = stderr);
  void list_poly_deg(int idx, FILE *fp = stderr);
  void list_polys(int style, FILE *fp = stderr);
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

void tetra59::list_poly_int(int idx, FILE *fp)
{
  fprintf(fp, "%2d) N=%2d (%2d,%3d), (%2d,%3d), (%2d,%3d) %-20s\n", idx + 1,
          N(idx), term(0, idx), term(1, idx), term(2, idx), term(3, idx),
          term(4, idx), term(5, idx), tetra59_items[idx].comment.c_str());
}

void tetra59::list_poly_deg(int idx, FILE *fp)
{
  vector<double> deg_term(6);
  for (int i = 0; i < 6; i++) {
    deg_term[i] = rad2deg(term(i, idx) * M_PI / N(idx));
    deg_term[i] = round(deg_term[i]);
  }

  fprintf(fp, "%2d) (%2g,%4g),  (%2g,%4g),  (%3g,%4g) %-20s\n", idx + 1,
          deg_term[0], deg_term[1], deg_term[2], deg_term[3], deg_term[4],
          deg_term[5], tetra59_items[idx].comment.c_str());
}

void tetra59::list_polys(int style, FILE *fp)
{
  fprintf(stderr, "=====================================================\n");
  if (style == 1) {
    fprintf(stderr, "Dihedral Angles of faces counted as 1,2,3,4    (Pi/N)\n");
    fprintf(stderr, "         (a12,a34) (a13,a24) (a14,a23)\n");
  }
  else {
    fprintf(stderr, "Dihedral Angles of faces counted as 1,2,3,4 (Degrees)\n");
    fprintf(stderr, "    (a12,a34)   (a13,a24)   (a14,a23)\n");
  }
  fprintf(stderr, "=====================================================\n");
  int reg_group = 0;
  for (int i = 0; i < last_tetra59; i++) {
    if (reg_group != reg(i)) {
      reg_group = reg(i);
      fprintf(stderr, "Regge group %d\n", reg_group);
    }
    if (style == 1)
      list_poly_int(i, fp);
    else
      list_poly_deg(i, fp);
  }
  fprintf(fp, "\n");
  fprintf(fp, "Entries are given as six dihedral angles of faces 1,2,3,4\n");
  fprintf(fp, "(a12,a34,a13,a24,a14,a23) %s\n",
          (style == 1) ? "as multiples of pi/N" : "in degrees (rounded)");
  fprintf(fp, "An extra symbol indicates 15 previously known forms from\n");
  fprintf(fp, "V. G. Boltianskii, Hilbert's third problem, 1978\n");
  fprintf(fp, "Values of terms 5 and 6 from the paper are reversed\n");
  fprintf(fp, "The cases without a symbol are newly discovered forms\n");
  fprintf(fp, "No cases of Regge group 6 were previously known\n");
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
  else
    // rk: there are non-digits in the input
    return -1;

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

  string poly;                                    // tetrahedron number 1 to 59
  int s = 0;                                      // for special cases 1 and 2
  double angle = std::numeric_limits<int>::max(); // angle for special cases
                                                  // (default 45 degrees)
  bool allow_angles = false; // -w switch allows all angles for case s

  bool reflect = false;      // make reflection model
  bool scale_volume = false; // scale volume to 1
  char dih_order = 'a';      // order of dihedral pairs
  int show_pair = 0;         // add dihedral pair to model of permutation

  int list_polys = 0; // output the list 1 - integers 2 -degrees

  // output edge and/or face math
  bool verbose_edge_math = false;
  bool verbose_face_math = false;
  bool verbose_pair_math = false;

  // maps are managed
  OffColor off_color = OffColor("");

  int opacity[3] = {-1, -1, -1}; // transparency from 0 to 255, for v,e,f

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
   or more tetrahedrons with edges x,y,s-a,s-b,s-c,s-d where s = (a+b+c+d)/2
   The tetrahedra come in pairs by dihedral angle. There may be more than one
   pairing depending on the three sets of opposing angles. Of the pairings...
3) The two opposing dihedral angles at edges x and y will be equal
4) The edge lengths of x and y of the two tetrahedra will be equal
5) If e,f,g,h are the dihedral angles at the edges a,b,c,d then the dihedral
   angles at edges s-a,s-b,s-c,s-d are t-e,t-f,t-g,t-h where t = (e+f+g+h)/2
   
For more information see: https://arxiv.org/abs/1903.04929

Method of Building Tetrahedrons

The tetrahedrons vertices are numbered 1,2,3 and 4. The dihedral angles given
for the six edges from the paper are (a12,a34,a13,a24,a14,a23). An edge length
of 1 is drawn from the origin to x,y,z = 0,1,0. Face angles for face A and B
are calculated from dihedral angles using trigonometry. Edge a (e13) and
edge b (e14) lengths are calculated from face angles from A and B angles using
trigonometry. Edge a (e13) is drawn rotated into place on the xy plane. Edge
b (e14) is also drawn and rotated on the xy plane and additionally rotated on
the Y axis by dihedral angle a12. At this point 4 vertices exist in 3D space
and the four faces A, B, C and D can be drawn.

Comparing Two Tetrahedra in the Same Regge Group

In a Regge group, two tetrahedra that share the same dihedral angle pair can
be compared and Regge symmetry math demonstrated. Not all angle pairs are in
the same columns pairs. The -d option can choose which order to place the 
dihedral pairs so that the pairs to be compared will be moved to the first
column. The -p option will find the pairing if it exists and move its order
to the first column. -p can take 1 or 2 depending on which order it picks for
the second and third column. The listing -l deg to list the 59 tetrahedra in
terms of degrees will help see which pairing are valid.

)");
}

void tetra59_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] tetrahedron

Generate 59 Tetrahedra with Rational Dihedral Angles in off format. The 59
Sporadic Tetrahedra is from a paper by Kiran S. Kedlaya, Alexander Kolpakov,
Bjorn Poonen, and Mechael Rubinstein: Space Vectors Forming Rational Angles
The paper can be found at: https://arxiv.org/abs/2011.14232
There are also two special infinite cases included in the findings
The first case was published by M.J.M Hill in 1895. The second case is new.

Options
%s
  -H        abstract from the paper and description of regge symmetry
  -v <opt>  verbose output math: e - edges, f - faces, p - pairs, a - all
  -L <int>  display the list of Sporadic Tetrahedra as integer=1, degrees=2
  -o <file> write output to file (default: write to standard output)

Scene Options
  -d <mthd> order of dihedral angles, for matching dihedral angle at position
            term 1 (a12) for regge symmetry (also permutes special cases) 
            angles permuted as pairs: 1-(a12,a34) 2-(a13,a24) 3-(a14,a23)
            for each pair in position one, there are two arrangements
            a:1,2,3  b:1,3,2; c:2,1,3  d:2,3,1; e:3,1,2  f:3,2,1 (default:a)
  -p <int>  build regge pair if found (1 or 2) not for special cases (sets -z)
  -r        reflect
  -z        scale volume to 1

Special Cases
  -s <int>  special case 1 or 2 (default: none)
              case 1: (pi/2, pi/2, pi - 2x, pi/3, x, x)
                        for pi/6 < x < pi/2 (30 < x < 90 degrees)
              case 2: (5pi/6 - x, pi/6 + x, 2pi/3 - x, 2pi/3 - x, x, x)
                        for pi/6 < x <= pi/3 (30 < x <= 60 degrees)
  -a <ang>  angle in degrees (default: 45)
  -w        allow any angle for case s (for testing case 2, 60 < x < 90
              which seems to yield a tetrahedron) 
           
Coloring Options (run 'off_util -H color' for help on color formats)
keyword: none - sets no color
  -F <col>  color the faces according to: (default: u)
              a color value - apply to all faces
              u - unique color
              s - symmetric coloring [,sub_group,conj_type]
              r - regge group (1-16) (not for special cases)
  -E <col>  color the edges according to: (default: lightgray)
              a color value - apply to all edges
              s - symmetric coloring [,sub_group,conj_type]
  -V <col>  color the vertices according to: (default: gold)
              a color value - apply to all vertices
              s - symmetric coloring [,sub_group,conj_type]
  -T <t,e>  transparency. from 0 (invisible) to 255 (opaque). element is any
            or all of, v - vertices, e - edges, f - faces, a - all (default: f)
  -m <maps> a comma separated list of color maps used to transform color
            indexes (default: compound), a part consisting of letters from
            v, e, f, selects the element types to apply the map list to
            (default 'vef'). use map name of 'index' to output index numbers
              compound:   yellow,red,darkgreen,blue,magenta,cyan,darkorange1
              rainbow16:  (special map for -F r)

)",
          prog_name(), help_ver_text);
}

void tetra59_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;
  int num;
  string arg_id;

  Split parts;
  Color col;
  vector<string> map_files;

  off_color.set_f_col_op('u');
  off_color.set_e_col(Color(211, 211, 211)); // lightgray
  off_color.set_v_col(Color(255, 215, 0));   // gold

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hHL:v:rzd:p:s:a:wV:E:F:T:m:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'H':
      extended_help();
      exit(0);

    case 'L':
      print_status_or_exit(
          get_arg_id(optarg, &arg_id, "int=1|deg=2", argmatch_add_id_maps), c);
      list_polys = atoi(arg_id.c_str());
      break;

    case 'v': {
      if (strspn(optarg, "efpa") != strlen(optarg))
        error(msg_str("verbose listings are '%s' must be any or all from "
                      "e, f, p, a",
                      optarg),
              c);
      string verbose = optarg;
      // set booleans once
      verbose_face_math = (verbose.find_first_of("fa") != string::npos);
      verbose_edge_math = (verbose.find_first_of("ea") != string::npos);
      verbose_pair_math = (verbose.find_first_of("pa") != string::npos);
      break;
    }

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

    case 'p':
      print_status_or_exit(read_int(optarg, &show_pair), c);
      if (show_pair < 1 || show_pair > 2)
        error("must be 1 or 2", c);
      scale_volume = true;
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

    case 'V':
      if (col.read(optarg)) {
        off_color.set_v_col(col);
        break;
      }
      parts.init(optarg, ",");
      if (off_color.v_op_check((char *)parts[0], "s"))
        off_color.set_v_col_op(*parts[0]);
      else
        error("invalid coloring", c);

      if (!((strchr("sS", off_color.get_v_col_op()) && parts.size() < 4) ||
            parts.size() < 2))
        error("too many comma separated parts", c);

      if (strchr("sS", off_color.get_v_col_op()))
        off_color.set_v_sub_sym(strlen(optarg) > 2 ? optarg + 2 : "");
      break;

    case 'E':
      if (col.read(optarg)) {
        off_color.set_e_col(col);
        break;
      }
      parts.init(optarg, ",");
      if (off_color.e_op_check((char *)parts[0], "s"))
        off_color.set_e_col_op(*parts[0]);
      else
        error("invalid coloring", c);

      if (!((strchr("sS", off_color.get_e_col_op()) && parts.size() < 4) ||
            parts.size() < 2))
        error("too many comma separated parts", c);

      if (strchr("sS", off_color.get_e_col_op()))
        off_color.set_e_sub_sym(strlen(optarg) > 2 ? optarg + 2 : "");
      break;

    case 'F':
      if (col.read(optarg)) {
        off_color.set_f_col(col);
        break;
      }
      parts.init(optarg, ",");
      if (off_color.f_op_check((char *)parts[0], "usr"))
        off_color.set_f_col_op(*parts[0]);
      else
        error("invalid coloring", c);

      if (!((strchr("sS", off_color.get_f_col_op()) && parts.size() < 4) ||
            parts.size() < 2))
        error("too many comma separated parts", c);

      if (strchr("sS", off_color.get_f_col_op()))
        off_color.set_f_sub_sym(strlen(optarg) > 2 ? optarg + 2 : "");
      break;

    case 'T': {
      int parts_sz = parts.init(optarg, ",");
      if (parts_sz > 2)
        error("the argument has more than 2 parts", c);

      print_status_or_exit(read_int(parts[0], &num), c);
      if (num < 0 || num > 255)
        error("face transparency must be between 0 and 255", c);

      // if only one part, apply to faces as default
      if (parts_sz == 1) {
        opacity[FACES] = num;
      }
      else if (parts_sz > 1) {
        if (strspn(parts[1], "vefa") != strlen(parts[1]))
          error(msg_str("transparency elements are '%s' must be any or all "
                        "from  v, e, f, a",
                        parts[1]),
                c);

        string str = parts[1];
        if (str.find_first_of("va") != string::npos)
          opacity[VERTS] = num;
        if (str.find_first_of("ea") != string::npos)
          opacity[EDGES] = num;
        if (str.find_first_of("fa") != string::npos)
          opacity[FACES] = num;
      }
      break;
    }

    case 'm':
      map_files.push_back(optarg);
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (argc == optind && !list_polys && !s)
    error("no tetrahedron specified", "tetrahedron");

  if (argc - optind > 1)
    error("too many arguments");

  if (argc - optind == 1)
    poly = argv[optind];

  if (s) {
    // angle default if not set
    if (angle == std::numeric_limits<int>::max()) {
      angle = 45;
      warning("angle default of 45 is set", 'a');
    }

    if (show_pair)
      error("make dihedral pair is not valid for special cases", 'p');

    char op = off_color.get_f_col_op();
    if (op && strchr("rR", op)) {
      warning("coloring by regge group has no effect for special cases", 'f');
      off_color.set_f_col_op('u'); // set to default
    }
  }
  else {
    if (angle != std::numeric_limits<int>::max())
      warning("angle has no effect for list cases", 'a');
    if (allow_angles)
      warning("allow angles has no effect for list cases", 'w');
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

  // set all maps in list
  for (unsigned int i = 0; i < map_files.size(); i++)
    print_status_or_exit(read_colorings(off_color.clrngs, map_files[i].c_str()),
                         'm');

  // fill in missing maps
  string default_map_name = "compound";
  for (unsigned int i = 0; i < 3; i++) {
    string map_name = default_map_name;
    // if map is already set, skip
    if (off_color.clrngs[i].get_cmaps().size())
      continue;
    else if (i == FACES) {
      char op = off_color.get_f_col_op();
      if (op && strchr("rR", op))
        map_name = "rainbow16";
    }
    off_color.clrngs[i].add_cmap(colormap_from_name(map_name.c_str()));
  }
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
// opposed to dihedral angle a is returned. Calculation was derived from
// https://en.wikipedia.org/wiki/Dihedral_angle#Geometry
// https://web.archive.org/web/20151125044900/http://www.had2know.com/academics/dihedral-angle-calculator-polyhedron.html
double face_cos_a(double a, double b, double c)
{
  return (cos(a) + cos(b) * cos(c)) / (sin(b) * sin(c));
}

// the dihedral angles are already calculated and are in radians
Geometry make_poly(double a12, double a34, double a13, double a24, double a14,
                   double a23, const tetra59_opts &opts)
{
  // face A angles
  double angle213 = face_cos_a(a14, a12, a13); // used for face A angle
  double angle123 = face_cos_a(a24, a12, a23); // used for edge a (e13)
  double angle132 = face_cos_a(a34, a13, a23); // used for edge a (a13)

  // face B angles
  double angle142 = face_cos_a(a34, a14, a24); // used for edge b (e14)
  double angle214 = face_cos_a(a13, a12, a14); // used for face B angle
  double angle124 = face_cos_a(a23, a12, a24); // used for edge b (e14)

  // face C angles (none needed to build)
  double angle324 = face_cos_a(a12, a23, a24);
  double angle234 = face_cos_a(a13, a23, a34);
  double angle243 = face_cos_a(a14, a24, a34);

  // face D angles (none needed to build)
  double angle134 = face_cos_a(a23, a13, a34);
  double angle143 = face_cos_a(a24, a14, a34);
  double angle314 = face_cos_a(a12, a13, a14);

  // clang-format off
  if (opts.verbose_face_math) {
    fprintf(stderr, "\n====================================================\n");
    fprintf(stderr, "Face Math:\n");
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

  if (opts.verbose_edge_math) {
    fprintf(stderr, "\n====================================================\n");
    fprintf(stderr, "Edge Math:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "dihedral angles at:\n");
    fprintf(stderr, "x (a12) = %.17lf (%g deg)\n", a12, rad2deg(a12));
    fprintf(stderr, "y (a34) = %.17lf (%g deg)\n", a34, rad2deg(a34));
    fprintf(stderr, "a (a13) = %.17lf (%g deg)\n", a13, rad2deg(a13));
    fprintf(stderr, "b (a24) = %.17lf (%g deg)\n", a24, rad2deg(a24));
    fprintf(stderr, "c (a14) = %.17lf (%g deg)\n", a14, rad2deg(a14));
    fprintf(stderr, "d (a23) = %.17lf (%g deg)\n", a23, rad2deg(a23));
  }

  // edge lengths
  double edge12 = 1; // prime edge
  double edge13 = (edge12 / cos2sin(angle132)) * cos2sin(angle123); // face A
  double edge14 = (edge12 / cos2sin(angle142)) * cos2sin(angle124); // face B
  double edge23 = (edge12 / cos2sin(angle132)) * cos2sin(angle213); // info only
  double edge24 = (edge12 / cos2sin(angle142)) * cos2sin(angle214); // info only
  double edge34 = (edge14 / cos2sin(angle134)) * cos2sin(angle314); // info only

  if (opts.verbose_edge_math) {
    fprintf(stderr, "\n");
    fprintf(stderr, "edge lengths:\n");
    fprintf(stderr, "x (e12) = %.17lf\n", edge12);
    fprintf(stderr, "y (e34) = %.17lf\n", edge34);
    fprintf(stderr, "a (e13) = %.17lf\n", edge13);
    fprintf(stderr, "b (e24) = %.17lf\n", edge24);
    fprintf(stderr, "c (e14) = %.17lf\n", edge14);
    fprintf(stderr, "d (e23) = %.17lf\n", edge23);
  }

  Geometry geom;

  // to build reflection, negate
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

  // dihedral angle a12 is used for angle between face A and B
  // rotate face B around Y
  pgeom.transform(Trans3d::rotate(0, mult * a12, 0));
  geom.append(pgeom);

  // add oriented faces, color coded
  geom.add_face({2, 1, 0}); // A yellow
  geom.add_face({0, 1, 3}); // B red
  geom.add_face({1, 2, 3}); // C green
  geom.add_face({3, 2, 0}); // D blue

  return geom;
}

// returns regge_grp, dihedral_angles returned in dih
Geometry case59(int &regge_grp, vector<double> &dih,
                const vector<int> &dihedral_order, const tetra59_opts &opts)
{
  tetra59 tetra59s;

  int sym_no = tetra59s.lookup_sym_no(opts.poly);
  if (sym_no >= tetra59s.get_last_tetra59())
    opts.error("tetrahedron number '" + opts.poly + "' out of range");
  if (sym_no < 0)
    opts.error("unknown tetrahedron '" + opts.poly + "'");

  // integer format from the sporadic cases list
  vector<int> term(6);
  for (unsigned int i = 0; i < 3; i++) {
    term[i * 2] = tetra59s.term(dihedral_order[i] * 2, sym_no);
    term[i * 2 + 1] = tetra59s.term(dihedral_order[i] * 2 + 1, sym_no);
  }

  int N = tetra59s.N(sym_no);

  fprintf(stderr,
          "%2d) N=%2d (%2d,%3d), (%2d,%3d), (%2d,%3d) (permutation %c)\n",
          sym_no + 1, N, term[0], term[1], term[2], term[3], term[4], term[5],
          opts.dih_order);

  // compute dihedral angles in radians
  for (unsigned int i = 0; i < 6; i++)
    dih[i] = term[i] * M_PI / N;

  regge_grp = tetra59s.reg(sym_no);

  // output terms as reduced fractions of pi
  if (opts.verbose_edge_math || opts.verbose_face_math) {
    fprintf(stderr, "\nregge group = %d\n", regge_grp);

    vector<int> gds(6);
    vector<string> numerator(6);
    for (unsigned int i = 0; i < 6; i++) {
      gds[i] = (int)gcd(term[i], N);
      numerator[i] =
          (term[i] / gds[i] == 1) ? "" : std::to_string(term[i] / gds[i]);
    }

    fprintf(stderr, "\n");
    fprintf(stderr, "%spi/%d, %spi/%d, %spi/%d, %spi/%d, %spi/%d, %spi/%d\n",
            numerator[0].c_str(), N / gds[0], numerator[1].c_str(), N / gds[1],
            numerator[2].c_str(), N / gds[2], numerator[3].c_str(), N / gds[3],
            numerator[4].c_str(), N / gds[4], numerator[5].c_str(), N / gds[5]);
  }

  return make_poly(dih[0], dih[1], dih[2], dih[3], dih[4], dih[5], opts);
}

Geometry case1(const vector<int> &dihedral_order, const tetra59_opts &opts)
{
  // (pi/2, pi/2, pi − 2x, pi/3, x, x)
  // for pi/6 < x < pi/2 (30 < x < 90 degrees)

  if (opts.verbose_edge_math || opts.verbose_face_math)
    fprintf(stderr, "\ncase 1: angle = %g\n", rad2deg(opts.angle));

  vector<double> term(6);
  term[0] = M_PI / 2;
  term[1] = M_PI / 2;
  term[2] = M_PI - 2 * opts.angle;
  term[3] = M_PI / 3;
  term[4] = opts.angle;
  term[5] = opts.angle;

  vector<double> dih(6);
  for (unsigned int i = 0; i < 3; i++) {
    dih[i * 2] = term[dihedral_order[i] * 2];
    dih[i * 2 + 1] = term[dihedral_order[i] * 2 + 1];
  }

  return make_poly(dih[0], dih[1], dih[2], dih[3], dih[4], dih[5], opts);
}

Geometry case2(const vector<int> &dihedral_order, const tetra59_opts &opts)
{
  // (5pi/6 − x, pi/6 + x, 2pi/3 − x, 2pi/3 − x, x, x)
  // for pi/6 < x ≤ pi/3 (30 < x <= 60 degrees)

  if (opts.verbose_edge_math || opts.verbose_face_math)
    fprintf(stderr, "\ncase 2: angle = %g\n", rad2deg(opts.angle));

  vector<double> term(6);
  term[0] = 5 * M_PI / 6 - opts.angle;
  term[1] = M_PI / 6 + opts.angle;
  term[2] = 2 * M_PI / 3 - opts.angle;
  term[3] = 2 * M_PI / 3 - opts.angle;
  term[4] = opts.angle;
  term[5] = opts.angle;

  vector<double> dih(6);
  for (unsigned int i = 0; i < 3; i++) {
    dih[i * 2] = term[dihedral_order[i] * 2];
    dih[i * 2 + 1] = term[dihedral_order[i] * 2 + 1];
  }

  return make_poly(dih[0], dih[1], dih[2], dih[3], dih[4], dih[5], opts);
}

void set_dihedral_order(vector<int> &dihedral_order, const tetra59_opts &opts)
{
  // a:1,2,3; b:1,3,2; c:2,1,3; d:2,3,1; e:3,1,2; f:3,2,1
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
}

void face_coloring(Geometry &geom, int regge_grp, tetra59_opts &opts)
{
  char op = opts.off_color.get_f_col_op();
  if (op && strchr("rR", op)) {
    Color col = opts.off_color.clrngs[FACES].get_col(regge_grp - 1);
    opts.off_color.set_f_col(col);
  }

  // any other color options done by class
  Status stat;
  if (!(stat = opts.off_color.off_color_main(geom)))
    opts.error(stat.msg());

  // apply all element transparencies
  apply_transparencies(geom, opts.opacity);
}

void scale_volume_to_one(Geometry &geom, vector<double> &len,
                         const tetra59_opts &opts)
{
  GeometryInfo info(geom);
  double scale = pow(fabs(info.volume()), 1.0 / 3.0);
  geom.transform(Trans3d::scale(1 / scale));

  vector<string> e{"12", "34", "13", "24", "14", "23"};
  vector<char> v{'x', 'y', 'a', 'b', 'c', 'd'};

  double sum = 0;
  for (unsigned int i = 0; i < 6; i++) {
    int a = e[i][0] - '0';
    int b = e[i][1] - '0';
    len[i] = (geom.verts(a - 1) - geom.verts(b - 1)).len();
    sum += len[i];
  }

  if (opts.verbose_edge_math) {
    fprintf(stderr, "\nvolume is set to 1. edges have been resized:\n");

    for (unsigned int i = 0; i < 6; i++)
      fprintf(stderr, "%c (e%s) = %.17lf\n", v[i], e[i].c_str(), len[i]);
    fprintf(stderr, "sum    = %.17lf\n", sum);
  }
}

void find_pair(int &row_no, int &col_no, const vector<int> &dihedral_order,
               const tetra59_opts &opts)
{
  tetra59 tetra59s;

  int sym_no = tetra59s.lookup_sym_no(opts.poly);
  int regge_grp = tetra59s.reg(sym_no);

  vector<int> term(6);
  for (unsigned int i = 0; i < 3; i++) {
    term[i * 2] = tetra59s.term(dihedral_order[i] * 2, sym_no);
    term[i * 2 + 1] = tetra59s.term(dihedral_order[i] * 2 + 1, sym_no);
  }

  // in radians
  vector<double> sym_term(2);
  sym_term[0] = term[0] * M_PI / tetra59s.N(sym_no);
  sym_term[1] = term[1] * M_PI / tetra59s.N(sym_no);

  int last_tetra59 = tetra59s.get_last_tetra59();
  for (int i = 0; i < last_tetra59; i++) {
    if (tetra59s.reg(i) != regge_grp)
      continue;
    if (i == sym_no)
      continue;

    for (int j = 0; j < 3; j++) {
      vector<double> col_term(2);
      col_term[0] = tetra59s.term(j * 2, i) * M_PI / tetra59s.N(i);
      col_term[1] = tetra59s.term(j * 2 + 1, i) * M_PI / tetra59s.N(i);
      if (sym_term[0] == col_term[0] && sym_term[1] == col_term[1]) {
        col_no = j + 1;
        break;
      }
    }

    if (col_no) {
      row_no = i;
      break;
    }
  }
}

// some opts variables are changed
Geometry make_dihedral_pair(int row_no, int col_no, const vector<double> &dih1,
                            const vector<double> &len1, tetra59_opts opts)
{
  string poly_base = opts.poly;
  char dih_order_base = opts.dih_order;

  // reset opts
  opts.poly = std::to_string(row_no + 1);

  // if using the first permutation tuple, the exact match is the second
  // permutation tuple and vice versa
  string first_tuples = "ace";
  if (first_tuples.find(opts.dih_order) != string::npos)
    opts.show_pair = (opts.show_pair == 1 ? 2 : 1);

  if (col_no == 1)
    opts.dih_order = (opts.show_pair == 1 ? 'a' : 'b');
  else if (col_no == 2)
    opts.dih_order = (opts.show_pair == 1 ? 'c' : 'd');
  else if (col_no == 3)
    opts.dih_order = (opts.show_pair == 1 ? 'e' : 'f');

  vector<int> dihedral_order;
  set_dihedral_order(dihedral_order, opts);

  int regge_grp = 0;
  vector<double> dih2(6);
  Geometry geom = case59(regge_grp, dih2, dihedral_order, opts);

  face_coloring(geom, regge_grp, opts);

  // equal volumes required for comparison
  vector<double> len2(6);
  scale_volume_to_one(geom, len2, opts);

  if (opts.verbose_pair_math) {
    fprintf(stderr, "\n====================================================\n");
    fprintf(stderr, "math relations between the two tetrahedra:\n\n");

    vector<string> e{"12", "34", "13", "24", "14", "23"};
    vector<char> v{'x', 'y', 'a', 'b', 'c', 'd'};

    fprintf(stderr, "edge lengths:\n");
    fprintf(stderr, "%c  (e%s) = %.17lf\n", v[0], e[0].c_str(), len1[0]);
    fprintf(stderr, "%c  (e%s) = %.17lf\n", v[1], e[1].c_str(), len1[1]);
    fprintf(stderr, "\n");

    fprintf(stderr, "tetrahedron %s (permutation %c)\n", poly_base.c_str(),
            dih_order_base);
    for (unsigned int i = 2; i < 6; i++)
      fprintf(stderr, "%c  (e%s) = %.17lf\n", v[i], e[i].c_str(), len1[i]);
    fprintf(stderr, "\n");

    fprintf(stderr, "tetrahedron %s (permutation %c)\n", opts.poly.c_str(),
            opts.dih_order);
    for (unsigned int i = 2; i < 6; i++)
      fprintf(stderr, "%c' (e%s) = %.17lf\n", v[i], e[i].c_str(), len2[i]);
    fprintf(stderr, "\n");

    double s = 0;
    for (unsigned int i = 2; i < 6; i++)
      s += len1[i];
    s /= 2;
    fprintf(stderr, "s = (a+b+c+d)/2 = %.17lf\n", s);
    fprintf(stderr, "\n");

    fprintf(stderr, "tetrahedron %s (permutation %c)\n", poly_base.c_str(),
            dih_order_base);
    for (unsigned int i = 2; i < 6; i++)
      fprintf(stderr, "s-%c  = %.17lf\n", v[i], s - len1[i]);
    fprintf(stderr, "\n");

    fprintf(stderr, "tetrahedron %s (permutation %c)\n", opts.poly.c_str(),
            opts.dih_order);
    for (unsigned int i = 2; i < 6; i++)
      fprintf(stderr, "s-%c' = %.17lf\n", v[i], s - len2[i]);
    fprintf(stderr, "\n");

    fprintf(stderr, "difference between a,b,c,d\n");
    for (unsigned int k = 2; k < 6; k++) {
      int j = k;
      for (unsigned int i = 2; i < 6; i++) {
        if (j > 5)
          j = 2;
        else if (j < 2)
          j = 5;
        fprintf(stderr, "%c-%c' = % -.17lf\n", v[i], v[j], len1[i] - len2[j]);
        j += (is_even(k) ? 1 : -1);
      }
      fprintf(stderr, "\n");
    }

    vector<char> d{'x', 'y', 'e', 'f', 'g', 'h'};

    fprintf(stderr, "edge dihedral angles:\n");
    fprintf(stderr, "%c  (a%s) = %.17lf (deg %g)\n", d[0], e[0].c_str(),
            dih1[0], rad2deg(dih1[0]));
    fprintf(stderr, "%c  (a%s) = %.17lf (deg %g)\n", d[1], e[1].c_str(),
            dih1[1], rad2deg(dih1[1]));
    fprintf(stderr, "\n");

    fprintf(stderr, "tetrahedron %s (permutation %c)\n", poly_base.c_str(),
            dih_order_base);
    for (unsigned int i = 2; i < 6; i++)
      fprintf(stderr, "%c  (a%s) = %.17lf (deg %g)\n", d[i], e[i].c_str(),
              dih1[i], rad2deg(dih1[i]));
    fprintf(stderr, "\n");

    fprintf(stderr, "tetrahedron %s (permutation %c)\n", opts.poly.c_str(),
            opts.dih_order);
    for (unsigned int i = 2; i < 6; i++)
      fprintf(stderr, "%c' (a%s) = %.17lf (deg %g)\n", d[i], e[i].c_str(),
              dih2[i], rad2deg(dih2[i]));
    fprintf(stderr, "\n");

    double t = 0;
    for (unsigned int i = 2; i < 6; i++)
      t += dih1[i];
    t /= 2;
    fprintf(stderr, "t = (e+f+g+h)/2 = %.17lf\n", t);
    fprintf(stderr, "\n");

    fprintf(stderr, "tetrahedron %s (permutation %c)\n", poly_base.c_str(),
            dih_order_base);
    for (unsigned int i = 2; i < 6; i++) {
      double rad = t - dih1[i];
      double deg = rad2deg(rad);
      fprintf(stderr, "t-%c  = %.17lf (deg %g)\n", d[i], rad, deg);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "tetrahedron %s (permutation %c)\n", opts.poly.c_str(),
            opts.dih_order);
    for (unsigned int i = 2; i < 6; i++) {
      double rad = t - dih2[i];
      double deg = rad2deg(rad);
      fprintf(stderr, "t-%c' = %.17lf (deg %g)\n", d[i], rad, deg);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "difference between e,f,g,h\n");
    for (unsigned int k = 2; k < 6; k++) {
      int j = k;
      for (unsigned int i = 2; i < 6; i++) {
        if (j > 5)
          j = 2;
        else if (j < 2)
          j = 5;
        double rad = dih1[i] - dih2[j];
        double deg = rad2deg(rad);
        fprintf(stderr, "%c-%c' = % -.17lf (deg %g)\n", d[i], d[j], rad, deg);
        j += (is_even(k) ? 1 : -1);
      }
      fprintf(stderr, "\n");
    }
  }

  return geom;
}

int main(int argc, char *argv[])
{
  tetra59_opts opts;
  opts.process_command_line(argc, argv);

  if (opts.list_polys) {
    tetra59 tetra59s;
    tetra59s.list_polys(opts.list_polys);
    exit(0);
  }

  vector<int> dihedral_order(3);
  set_dihedral_order(dihedral_order, opts);

  Geometry geom;
  int regge_grp = 0;
  vector<double> dih1(6);

  // 59 cases
  if (opts.s == 0)
    geom = case59(regge_grp, dih1, dihedral_order, opts);
  // special cases
  else if (opts.s == 1)
    geom = case1(dihedral_order, opts);
  else if (opts.s == 2)
    geom = case2(dihedral_order, opts);

  vector<double> len1(6);
  if (opts.scale_volume)
    scale_volume_to_one(geom, len1, opts);

  if ((opts.verbose_edge_math || opts.verbose_face_math) & !opts.show_pair) {
    Symmetry sym(geom);
    fprintf(stderr, "\n");
    fprintf(stderr, "the geometric symmetry of the tetrahedron is %s\n",
            sym.get_symbol().c_str());
  }

  face_coloring(geom, regge_grp, opts);

  if (opts.show_pair) {
    int row_no = 0;
    int col_no = 0;
    find_pair(row_no, col_no, dihedral_order, opts);

    if (!col_no)
      fprintf(
          stderr,
          "\nno matching pair for: tetrahedron number: %s permutation: %c\n",
          opts.poly.c_str(), opts.dih_order);
    else {
      fprintf(stderr, "\nbuilding matching pair:\n");
      geom.append(make_dihedral_pair(row_no, col_no, dih1, len1, opts));
    }
  }

  opts.write_or_error(geom, opts.ofile);

  return 0;
}
