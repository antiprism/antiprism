/*
   Copyright (c) 2007-2016, Roger Kaufman
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
   Name: conway.cc
   Description: Conway Notation
                Implementation of George Hart's Conway Notation
                http://www.georgehart.com/virtual-polyhedra/conway_notation.html
   Project: Antiprism - http://www.antiprism.com
*/

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <ctype.h>

#include <map>
#include <set>
#include <string>
#include <vector>

#include "../base/antiprism.h"

using std::string;
using std::vector;
using std::map;
using std::set;

using namespace anti;

#define CN_ONE_THIRD 1 / 3.0
#define CN_ONE_HALF 0.5

class ops {
public:
  int op_pos;
  char op;
  int op_var;
  ops(int n, char o, int v) : op_pos(n), op(o), op_var(v) {}
};

bool cmp_ops(const ops *a, const ops *b) { return a->op_pos > b->op_pos; }

int validate_cn_string(const string &cn_string, vector<ops *> &operations,
                       char &operand, int &poly_size)
{
  char current_op = '\0';
  string number_string;
  int num_val = 0;
  bool delayed_write = false;

  string operators = "abcdegjkmoprstwx";
  string operands = "TCOIDPAY";
  string digits = "0123456789";
  string digits_allowed = "tkPAY";

  operand = '\0';
  poly_size = 0;

  for (unsigned int i = 0; i < cn_string.length(); i++) {
    if (operators.find(cn_string[i]) != string::npos) {
      if (operand != '\0')
        return i;

      if (delayed_write) {
        operations.push_back(new ops(i, current_op, num_val));
        delayed_write = false;
      }

      current_op = cn_string[i];

      // delay if t or k
      if (digits_allowed.find(current_op) == string::npos)
        operations.push_back(new ops(i + 1, current_op, num_val));
      else
        delayed_write = true;
    }
    else if (operands.find(cn_string[i]) != string::npos) {
      if (operand != '\0')
        return i;

      if (delayed_write) {
        operations.push_back(new ops(i, current_op, num_val));
        delayed_write = false;
      }

      operand = cn_string[i];
      current_op = '\0';
    }
    else if (digits.find(digits) != string::npos) {
      if (operand != '\0') {
        if (digits_allowed.find(operand) == string::npos) {
          return i;
        }
      }
      else if (digits_allowed.find(current_op) == string::npos) {
        return i;
      }

      int digits_start = i;
      int digits_end = i;
      while ((digits_end + 1 < (int)cn_string.length()) &&
             (digits.find(cn_string[digits_end + 1]) != string::npos))
        digits_end++;

      number_string =
          cn_string.substr(digits_start, (digits_end - digits_start + 1));

      if (current_op == 'k') {
        num_val = atoi(number_string.c_str());
        if (num_val < 3) {
          fprintf(stderr, "kis k(n), n must be 3 or greater\n");
          return i;
        }
        operations.push_back(new ops(i + 1, current_op, num_val));
        delayed_write = false;
      }
      else if (current_op == 't') {
        num_val = atoi(number_string.c_str());
        if (num_val < 3) {
          fprintf(stderr, "trunc t(n), n must be 3 or greater\n");
          return i;
        }
        operations.push_back(new ops(i + 1, current_op, num_val));
        delayed_write = false;
      }
      else {
        poly_size = atoi(number_string.c_str());
        if (poly_size < 3) {
          fprintf(stderr, "P(n), A(n), or Y(n), n must be 3 or greater\n");
          return i;
        }
      }

      num_val = 0;
      i = digits_end;
    }
  }

  // if t or k was specified alone at end
  if (delayed_write)
    operations.push_back(new ops(cn_string.length(), current_op, num_val));

  // if P, A or Y was specified with no digit n
  if ((digits_allowed.find(operand) != string::npos) && poly_size == 0) {
    fprintf(stderr, "P(n), A(n), or Y(n), n must be 3 or greater\n");
    return cn_string.length();
  }

  return 0;
}

string resolved_cn_string(const string &cn_string,
                          const bool &use_truncate_algorithm)
{
  string resolve_string = cn_string;

  const int num_subst = 25;
  string target[num_subst];
  string resolve[num_subst];
  // G. Hart Commentary
  target[0] = "P4";
  resolve[0] = "C"; // P4 --> C   (C is prism)
  target[1] = "A3";
  resolve[1] = "O"; // A3 --> O   (O is antiprism)
  target[2] = "Y3";
  resolve[2] = "T"; // Y3 --> T   (T is pyramid)
  target[3] = "e";
  resolve[3] = "aa"; // e --> aa   (abbr. for explode)
  target[4] = "b";
  resolve[4] = "ta"; // b --> ta   (abbr. for bevel)
  target[5] = "o";
  resolve[5] = "jj"; // o --> jj   (abbr. for ortho)
  target[6] = "m";
  resolve[6] = "kj"; // m --> kj   (abbr. for meta)
  target[7] = "t";
  resolve[7] = "dk"; // t(n) --> dk(n)d  (dual operations) (see special case)
  target[8] = "j";
  resolve[8] = "dad"; // j --> dad  (dual operations)
  target[9] = "s";
  resolve[9] = "dgd"; // s --> dgd  (dual operations)
  target[10] = "dd";
  resolve[10] = ""; // dd --> null  (order 2)
  target[11] = "ad";
  resolve[11] = "a"; // ad --> a   (a_ = ad_)
  target[12] = "gd";
  resolve[12] = "g"; // gd --> g   (g_ = gd_)
  target[13] = "aY";
  resolve[13] = "A"; // aY --> A   (interesting fact)
  target[14] = "dT";
  resolve[14] = "T"; // dT --> T   (self-dual)
  target[15] = "gT";
  resolve[15] = "D"; // gT --> D   (symm change)
  target[16] = "aT";
  resolve[16] = "O"; // aT --> O   (symm change)
  target[17] = "dC";
  resolve[17] = "O"; // dC --> O   (dual pair)
  target[18] = "dO";
  resolve[18] = "C"; // dO --> C   (dual pair)
  target[19] = "dI";
  resolve[19] = "D"; // dI --> D   (dual pair)
  target[20] = "dD";
  resolve[20] = "I"; // dD --> I   (dual pair)
  target[21] = "aO";
  resolve[21] = "aC"; // aO --> aC  (for uniqueness)
  target[22] = "aI";
  resolve[22] = "aD"; // aI --> aD  (for uniqueness)
  target[23] = "gO";
  resolve[23] = "gC"; // gO --> gC  (for uniqueness)
  target[24] = "gI";
  resolve[24] = "gD"; // gI --> gD  (for uniqueness)

  // first 3 targets are positional
  if (resolve_string.length() > 1) {
    string tmp = resolve_string.substr(resolve_string.length() - 2, 2);
    for (unsigned int i = 0; i < 3; i++) {
      if (tmp == target[i])
        resolve_string.replace(resolve_string.length() - 2, target[i].length(),
                               resolve[i]);
    }
  }

  // 4 to 6
  for (unsigned int i = 3; i < 7; i++) {
    for (size_t pos = resolve_string.find_first_of(target[i], 0);
         pos != string::npos;
         pos = resolve_string.find_first_of(target[i], 0)) {
      resolve_string.replace(pos, target[i].length(), resolve[i]);
    }
  }

  // 7 is special case
  if (!use_truncate_algorithm) {
    for (size_t pos = resolve_string.find_first_of(target[7], 0);
         pos != string::npos;
         pos = resolve_string.find_first_of(target[7], 0)) {
      resolve_string.replace(pos, target[7].length(), resolve[7]);
      size_t pos2 = resolve_string.find_first_not_of("0123456789", pos + 2);
      if (pos2 != string::npos)
        resolve_string.insert(pos2, "d");
      else
        resolve_string.append("d");
    }
  }

  // 8 to end. Because of length greater than 1 of the target, string.find and
  // replace cannot be used
  for (int i = 8; i < num_subst; i++) {
    int j = 0;
    int stop = resolve_string.length() - target[i].length() + 1;
    while (j <= stop) {
      string tmp = resolve_string.substr(j, target[i].length());
      if (tmp == target[i]) {
        resolve_string.replace(j, target[i].length(), resolve[i]);
        j = -1; // we have to begin looking from the begining once again, ok if
                // we modify stop
        stop = resolve_string.length() - target[i].length() + 1;
      }
      j++;
    }
  }

  if (resolve_string != cn_string) {
    fprintf(stderr, "%s resolved to: ", cn_string.c_str());
    if (resolve_string.length() == 0)
      fprintf(stderr, "%s", "NOTHING");
    else
      fprintf(stderr, "%s", resolve_string.c_str());
    fprintf(stderr, "\n");
  }

  return resolve_string;
}

class cn_opts : public ProgramOpts {
public:
  string ifile;
  string ofile;

  string cn_string;
  bool resolve_defeat;
  bool reverse_defeat;
  char operand;
  int poly_size;
  char planarization_method;
  char canonical_method;
  int num_iters_planar;
  int num_iters_canonical;
  int rep_count;
  bool unitize;
  bool verbosity;
  bool use_truncate_algorithm;
  char face_Coloring_method;
  int face_opacity;
  string face_pattern;

  double epsilon;

  Color vert_col;
  Color edge_col;

  ColorMapMulti map;

  vector<ops *> operations;

  cn_opts()
      : ProgramOpts("conway"), cn_string(""), resolve_defeat(false),
        reverse_defeat(false), operand('\0'), poly_size(0),
        planarization_method('p'), canonical_method('\0'),
        num_iters_planar(1000), num_iters_canonical(1000), rep_count(-1),
        unitize(false), verbosity(false), use_truncate_algorithm(false),
        face_Coloring_method('n'), face_opacity(255), face_pattern("1"),
        epsilon(0), vert_col(Color(255, 215, 0)), // gold
        edge_col(Color(211, 211, 211))            // lightgrey
  {
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

void extended_help()
{
  fprintf(
      stdout,
      "\n"
      "The following is a description of Conway Notation edited from the "
      "Conway\n"
      "Notation web page by George W. Hart (http://www.georgehart.com)\n"
      "More detailed information and examples can be found at\n"
      "http://www.georgehart.com/virtual-polyhedra/conway_notation.html\n"
      "and at\n"
      "http://en.wikipedia.org/wiki/Conway_polyhedron_notation\n"
      "\n"
      "Basics: In this notation, one specifies a \"seed\" polyhedron with a "
      "capital\n"
      "letter. Operations to perform on any polyhedron are specified with "
      "lower-case\n"
      "letters preceding it. This program contains a small set of seeds and "
      "operators\n"
      "from which an infinite number of derived polyhedra can be generated.\n"
      "Note: This C++ port of Conway Notation can also operate on OFF files "
      "from\n"
      "standard input if the seed polyhedron is not specified.\n"
      "\n"
      "Seeds: The platonic solids are denoted T, O, C, I, and D, according to "
      "their\n"
      "first letter. Other polyhedra which are implemented here include "
      "prisms, Pn,\n"
      "antiprisms, An, and pyramids, Yn, where n is a number (3 or greater) "
      "which you\n"
      "specify to indicate the size of the base you want, e.g., Y3=T, P4=C, "
      "and A3=O.\n"
      "\n"
      "Operations: Currently, abcdegjkmoprst are defined. They are motivated "
      "by the\n"
      "operations needed to create the Archimedean solids and their duals from "
      "the\n"
      "platonic solids.  Try each on a cube:\n"
      "\n"
      "a = ambo   The ambo operation can be thought of as truncating to the "
      "edge\n"
      "midpoints.  It produces a polyhedron, aX, with one vertex for each edge "
      "of X.\n"
      "There is one face for each face of X and one face for each vertex of "
      "X.\n"
      "Notice that for any X, the vertices of aX are all 4-fold, and that "
      "aX=adX.\n"
      "If two mutually dual polyhedra are in \"dual position,\" with all edges "
      "tangent\n"
      "to a common sphere, the ambo of either is their intersection.  For "
      "example\n"
      "aC=aO is the cuboctahedron.\n"
      "Note: ambo is also known as \"rectifying\" the polyhedron, or "
      "rectification\n"
      "\n"
      "b = bevel  The bevel operation can be defined by bX=taX.  bC is the "
      "truncated\n"
      "cuboctahedron.\n"
      "Note: bevel is also known as \"omnitruncating\" the polyhedron, or "
      "omnitruncation\n"
      "\n"
      "c = chamfer   New hexagonal faces are added in place of edges.\n"
      "\n"
      "d = dual   The dual of a polyhedron has a vertex for each face, and a "
      "face for\n"
      "each vertex, of the original polyhedron, e.g., dC=O.  Duality is an "
      "operation\n"
      "of order two, meaning for any polyhedron X, ddX=X, e.g., ddC=dO=C.\n"
      "\n"
      "e = expand This is Mrs. Stott's expansion operation.  Each face of X "
      "is\n"
      "separated from all its neighbors and reconnected with a new 4-sided "
      "face,\n"
      "corresponding to an edge of X.  An n-gon is then added to connect the "
      "4-sided\n"
      "faces at each n-fold vertex.  For example, eC is the "
      "rhombicuboctahedron.  It\n"
      "turns out that eX=aaX and so eX=edX.\n"
      "Note: expand is also known as \"cantellating\" the polyhedron, or "
      "cantellation\n"
      "\n"
      "g = gyro   The dual operation to s is g. gX=dsdX=dsX, with all 5-sided "
      "faces.\n"
      "The gyrocube, gC=gO=\"pentagonal icositetrahedron,\" is dual to the "
      "snub cube.\n"
      "g is like k but with the new edges connecting the face centers to the "
      "1/3\n"
      "points on the edges rather than the vertices.\n"
      "\n"
      "j = join   The join operator is dual to ambo, so jX=dadX=daX.  jX is "
      "like kX\n"
      "without the original edges of X.  It produces a polyhedron with one "
      "4-sided\n"
      "face for each edge of X.  For example, jC=jO is the rhombic "
      "dodecahedron.\n"
      "\n"
      "k = kis    All faces are processed or kn = just n-sided faces are "
      "processed\n"
      "The kis operation divides each n-sided face into n triangles.  A new "
      "vertex is\n"
      "added in the center of each face, e.g., the kiscube, kC, has 24 "
      "triangular\n"
      "faces.  The k operator is dual to t, meaning kX=dtdX.\n"
      "\n"
      "m = meta   Dual to b, mX=dbX=kjX.  mC has 48 triangular faces.  m is "
      "like k\n"
      "and o combined; new edges connect new vertices at the face centers to "
      "the old\n"
      "vertices and new vertices at the edge midpoints.  mX=mdX.  mC is the\n"
      "\"hexakis octahedron.\"\n"
      "\n"
      "o = ortho  Dual to e, oX=deX=jjX.  oC is the trapezoidal "
      "icositetrahedron, with\n"
      "24 kite-shaped faces.  oX has the effect of putting new vertices in the "
      "middle\n"
      "of each face of X and connecting them, with new edges, to the edge "
      "midpoints of\n"
      "X.\n"
      "\n"
      "p = propellor    Makes each n-gon face into a \"propellor\" of an "
      "n-gon\n"
      "surrounded by n quadrilaterals, e.g., pT is the tetrahedrally "
      "stellated\n"
      "icosahedron. Try pkD and pt6kT. p is a self-dual operation, i.e., "
      "dpdX=pX and\n"
      "dpX=pdX, and p also commutes with a and j, i.e., paX=apX. (This and the "
      "next\n"
      "are extensions were added by George Hart and not specified by Conway)\n"
      "\n"
      "r = reflect   Changes a left-handed solid to right handed, or vice "
      "versa, but\n"
      "has no effect on a reflexible solid. So rC=C, but compare sC and rsC.\n"
      "\n"
      "s = snub   The snub operation produces the snub cube, sC, from C.  It "
      "can be\n"
      "thought of as eC followed by the operation of slicing each of the new "
      "4-fold\n"
      "faces along a diagonal into two triangles.  With a consistent "
      "handedness to\n"
      "these cuts, all the vertices of sX are 5-fold.  Note that sX=sdX.\n"
      "\n"
      "t = truncate  All faces are processed or tn = just n-sided faces are "
      "processed\n"
      "Truncating a polyhedron cuts off each vertex, producing a new n-sided "
      "face for\n"
      "each n-fold vertex.  The faces of the original polyhedron still appear, "
      "but\n"
      "have twice as many sides, e.g., the tC has six octagonal sides "
      "corresponding to\n"
      "the six squares of the C, and eight triangles corresponding to the "
      "cube's eight\n"
      "vertices.\n"
      "\n"
      "w = whirl  Gyro followed by truncation of vertices centered on original "
      "faces.\n"
      "This create 2 new hexagons for every original edge"
      "\n"
      "x = null   Null operation. Nothing is changed. A planarize step is "
      "performed\n"
      "\n");
}

void cn_opts::usage()
{
  fprintf(
      stdout,
      "\n"
      "Usage: %s [options] [Conway Notation string] [input_file]\n"
      "\n"
      "Conway Notation uses algorithms by George W. Hart "
      "(http://www.georgehart.com)\n"
      "http://www.georgehart.com/virtual-polyhedra/conway_notation.html\n"
      "\n"
      "Read a polyhedron from a file in OFF format.\n"
      "If input_file is not given and no seed polyhedron is given in the "
      "notation\n"
      "string then the program reads from standard input.\n"
      "\n"
      "Options\n"
      "%s"
      "  -H        Conway Notation detailed help. seeds and operator "
      "descriptions\n"
      "  -d        don't simplify Conway Notation string\n"
      "  -r        execute operations in reverse order (left to right)\n"
      "  -t        use truncate algorithm instead of simplifying to \"dkd\"\n"
      "              also used in ambo as a truncation of 1/2\n"
      "  -u        make final product be averge unit edge length\n"
      "  -v        verbose output\n"
      "  -o <file> write output to file (default: write to standard output)\n"
      "\n"
      "Canonicalization and Planarization options\n"
      "  -p <mthd> inter-step planarization method\n"
      "            p - conway notation planarize (face centroids reciprocal) "
      "(default)\n"
      "            q - conway notation planarize (face centroids magnitude "
      "reciprocal)\n"
      "            l - mathematica version of planarize\n"
      "  -c <mthd> canonicalize final product using:\n"
      "            n - conway notation version of canonicalization\n"
      "            m - mathematica version of canonicalization\n"
      "            or planarize final product with p, q or l above\n"
      "  -n <itrs> maximum number canonical iterations (default: 1000)\n"
      "  -l <lim>  minimum distance change to terminate canonicalization, as "
      "negative\n"
      "               exponent (default: %d giving %.0e)\n"
      "  -i <itrs> maximum inter-step planarization iterations (default: "
      "1000)\n"
      "               0 = no inter-step planarization\n"
      "  -z <n>    status reporting every n iterations, -1 for no status "
      "(default: -1)\n"
      "\n"
      "\nColoring Options (run 'off_util -H color' for help on color formats)\n"
      "  -V <col>  vertex color (default: gold)\n"
      "  -E <col>  edge color   (default: lightgray)\n"
      "  -f <mthd> mthd is face Coloring method using color in map\n"
      "               key word: none - sets no color (default: n)\n"
      "               n - color by number of sides\n"
      "               s - symmetric Coloring\n"
      "  -T <tran> face transparency. valid range from 0 (invisible) to 255 "
      "(opaque)\n"
      "  -O <strg> face transparency pattern string. valid values\n"
      "               0 - map color alpha value, 1 -T alpha applied (default: "
      "'1')\n"
      "  -m <maps> color maps for faces to be tried in turn (default: m1)\n"
      "               keyword m1: "
      "red,darkorange1,yellow,darkgreen,cyan,blue,magenta,\n"
      "                           white,grey,black\n"
      "               keyword m2: "
      "red,blue,green,yellow,brown,magenta,purple,grue,\n"
      "                           gray,orange (from George Hart\'s original "
      "applet)\n"
      "\n"
      "\n",
      prog_name(), help_ver_text, int(-log(::epsilon) / log(10) + 0.5),
      ::epsilon);
}

void cn_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  int sig_compare = INT_MAX;

  string map_file;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hHdrtuvp:c:n:l:i:z:f:V:E:T:O:m:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'H':
      extended_help();
      exit(0);

    case 'd':
      resolve_defeat = true;
      break;

    case 'r':
      reverse_defeat = true;
      break;

    case 't':
      use_truncate_algorithm = true;
      break;

    case 'u':
      unitize = true;
      break;

    case 'v':
      verbosity = true;
      break;

    case 'p':
      if (strlen(optarg) == 1 && strchr("pql", int(*optarg)))
        planarization_method = *optarg;
      else
        error("planarize method type must be p, q or l", c);
      break;

    case 'c':
      if (strlen(optarg) == 1 && strchr("nmpql", int(*optarg)))
        canonical_method = *optarg;
      else
        error("canonical method type must be n, m, p, q or l", c);
      break;

    case 'n':
      print_status_or_exit(read_int(optarg, &num_iters_canonical), c);
      if (num_iters_canonical <= 0)
        error("number of canonical iterations must be greater than 0", c);
      break;

    case 'l':
      print_status_or_exit(read_int(optarg, &sig_compare), c);
      if (sig_compare < 0) {
        warning("canonicalization limit is negative, and so ignored", c);
      }
      if (sig_compare > DEF_SIG_DGTS) {
        warning("canonicalization limit is very small, may not be attainable",
                c);
      }
      break;

    case 'i':
      print_status_or_exit(read_int(optarg, &num_iters_planar), c);
      if (num_iters_planar < 0)
        error("number of planarization iterations 0 or greater", c);
      break;

    case 'z':
      print_status_or_exit(read_int(optarg, &rep_count), c);
      if (rep_count < -1)
        error("number of iterations must be -1 or greater", c);
      break;

    case 'f':
      if (!strcasecmp(optarg, "none"))
        face_Coloring_method = '\0';
      else if (strspn(optarg, "ns") != strlen(optarg) || strlen(optarg) > 1)
        error(msg_str("invalid face Coloring method '%s'", optarg), c);
      else {
        face_Coloring_method = *optarg;
      }
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

    case 'O':
      if (strspn(optarg, "01") != strlen(optarg))
        error(msg_str("transparency string is '%s' must consist of "
                      "0 and 1's",
                      optarg),
              c);
      face_pattern = optarg;
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

  if (argc - optind > 2)
    error("too many arguments");

  // avoid segfault
  if (argc > 1)
    cn_string = argv[optind];

  if (!strlen(cn_string.c_str()))
    error("no Conway Notation string given");

  if (strspn(cn_string.c_str(), "abcdegjkmoprstwxTCOIDPAY0123456789,") !=
      strlen(cn_string.c_str()))
    error("Conway Notation must consist of abcdegjkmoprstwxTCOIDPAY0123456789");

  if (int pos = validate_cn_string(cn_string, operations, operand, poly_size))
    error(msg_str("Unexpected character in position %d: %s", pos + 1,
                  cn_string.c_str()));

  if (!resolve_defeat) {
    cn_string = resolved_cn_string(cn_string, use_truncate_algorithm);

    for (auto &operation : operations)
      delete operation;
    operations.clear();

    // revalidate (should be valid) to rebuild operations table
    if (int pos = validate_cn_string(cn_string, operations, operand, poly_size))
      error(msg_str("Unexpected character in position %d: %s", pos + 1,
                    cn_string.c_str()));
  }

  optind++;
  if (argc - optind == 1) {
    ifile = argv[optind];
    if (operand)
      error(
          msg_str("operand '%c' was specified so input file '%s' is unexpected",
                  operand, ifile.c_str()));
  }

  // operations are done in reverse order (unless defeated)
  if (!reverse_defeat)
    sort(operations.begin(), operations.end(), cmp_ops);

  if (!map_file.size())
    map_file = "m1";

  if (map_file == "m1" || map_file == "m2") {
    auto *col_map = new ColorMapMap;
    if (map_file == "m1") {
      col_map->set_col(0, Color(1.0, 0.0, 0.0)); // 3-sided faces red
      col_map->set_col(1,
                       Color(1.0, 0.49804, 0.0)); // 4-sided faces darkoranage1
      col_map->set_col(2, Color(1.0, 1.0, 0.0));  // 5-sided faces yellow
      col_map->set_col(3, Color(0.0, 0.39216, 0.0)); // 6-sided faces darkgreen
      col_map->set_col(4, Color(0.0, 1.0, 1.0));     // 7-sided faces cyan
      col_map->set_col(5, Color(0.0, 0.0, 1.0));     // 8-sided faces blue
      col_map->set_col(6, Color(1.0, 0.0, 1.0));     // 9-sided faces magenta
      col_map->set_col(7, Color(1.0, 1.0, 1.0));     // 10-sided faces white
      col_map->set_col(8, Color(0.5, 0.5, 0.5));     // 11-sided faces grey
      col_map->set_col(9, Color(0.0, 0.0, 0.0));     // 12-sided faces black
    }
    else if (map_file == "m2") {
      col_map->set_col(0, Color(0.9, 0.3, 0.3));   // 3-sided faces red
      col_map->set_col(1, Color(0.4, 0.4, 1.0));   // 4-sided faces blue
      col_map->set_col(2, Color(0.4, 0.4, 1.0));   // 5-sided faces green
      col_map->set_col(3, Color(0.9, 0.9, 0.2));   // 6-sided faces yellow
      col_map->set_col(4, Color(0.5, 0.25, 0.25)); // 7-sided faces brown
      col_map->set_col(5, Color(0.8, 0.2, 0.8));   // 8-sided faces magenta
      col_map->set_col(6, Color(0.5, 0.2, 0.8));   // 9-sided faces purple
      col_map->set_col(7, Color(0.1, 0.9, 0.9));   // 10-sided faces grue
      col_map->set_col(8, Color(0.5, 0.5, 0.5));   // 11-sided faces gray
      col_map->set_col(9, Color(1.0, 0.6, 0.1));   // 12-sided faces orange
    }
    col_map->set_wrap();
    map.add_cmap(col_map);
  }
  else
    print_status_or_exit(map.init(map_file.c_str()), 'm');

  epsilon = (sig_compare != INT_MAX) ? pow(10, -sig_compare) : ::epsilon;
}

void verbose(const char &operation, const int &op_var, const bool &verbosity)
{
  if (verbosity) {
    char buf[MSG_SZ];
    buf[0] = '\0';
    switch (operation) {
    case 'a':
      fprintf(stderr, "ambo\n");
      break;
    case '^':
      fprintf(stderr,
              "ambo as truncate to edge midpoints by built in function\n");
      break;
    case 'c':
      fprintf(stderr, "chamfer\n");
      break;
    case 'd':
      fprintf(stderr, "dual\n");
      break;
    case 'g':
      fprintf(stderr, "gyro\n");
      break;
    case 'k':
      if (op_var)
        sprintf(buf, "(%d)", op_var);
      fprintf(stderr, "kis%s\n", buf);
      break;
    case 'p':
      fprintf(stderr, "propellor\n");
      break;
    case 'r':
      fprintf(stderr, "reflect\n");
      break;
    case 'w':
      fprintf(stderr, "whirl as gyro, truncation of vertices centered on "
                      "original faces:\n");
      break;
    case 'b':
      fprintf(stderr, "bevel as truncate, ambo:\n");
      break;
    case 'e':
      fprintf(stderr, "expand as ambo, ambo:\n");
      break;
    case 'j':
      fprintf(stderr, "join as ambo, dual:\n");
      break;
    case 'm':
      fprintf(stderr, "meta as dual, kis, join:\n");
      break;
    case 'o':
      fprintf(stderr, "ortho as join, join:\n");
      break;
    case 's':
      fprintf(stderr, "snub as gyro, dual:\n");
      break;
    case 't':
      if (op_var)
        sprintf(buf, "(%d)", op_var);
      fprintf(stderr, "truncate as dual, kis%s, dual:\n", buf);
      break;
    case '#':
      if (op_var)
        sprintf(buf, "(%d)", op_var);
      fprintf(stderr, "truncate%s by built in function\n", buf);
      break;
    case 'x':
      fprintf(stderr, "null operation\n");
      break;
    case '_':
      fprintf(stderr, "planarizing ...\n");
      break;
    case '@':
      fprintf(stderr, "canonicalizing ...\n");
      break;
    case '$':
      fprintf(stderr, "done.\n");
      break;
    default:
      fprintf(stderr, "unexpected operator to verbose: '%c'\n", operation);
    }
  }
}

void unitize_edges(Geometry &geom)
{
  GeometryInfo info(geom);
  if (info.num_iedges() > 0) {
    double val = info.iedge_length_lims().sum / info.num_iedges();
    geom.transform(Trans3d::scale(1 / val));
  }
}

void centroid_to_origin(Geometry &geom)
{
  geom.transform(Trans3d::transl(-centroid(geom.verts())));
}

void cn_planarize(Geometry &geom, const char &planarization_method,
                  const int &num_iters_planar, const double &eps,
                  const bool &verbosity, const int &rep_count)
{
  int divergence_test = 10;

  if (num_iters_planar != 0) {
    verbose('_', 0, verbosity);
    if (planarization_method == 'p' || planarization_method == 'q')
      canonicalize_cn(&geom, num_iters_planar, planarization_method,
                      divergence_test, rep_count, eps);
    else
      canonicalize_mm(&geom, 50 / 100.0, 20 / 100.0, num_iters_planar,
                      divergence_test, rep_count, true, eps);
  }
}

void canonicalize(Geometry &geom, const char &canonical_method,
                  const int &num_iters_canonical, const double &eps,
                  const bool &verbosity, const int &rep_count)
{
  int divergence_test = 10;

  if (canonical_method == 'n') {
    verbose('@', 0, verbosity);
    centroid_to_origin(geom);
    canonicalize_cn(&geom, num_iters_canonical, canonical_method,
                    divergence_test, rep_count, eps);
  }
  else if (canonical_method == 'm') {
    verbose('@', 0, verbosity);

    // -C s options from canonical.cc found to be helpful with old canonical
    centroid_to_origin(geom);
    project_onto_sphere(&geom);

    canonicalize_mm(&geom, 50 / 100.0, 20 / 100.0, num_iters_canonical,
                    divergence_test, rep_count, false, eps);
  }
  else if (canonical_method)
    // final planarization
    cn_planarize(geom, canonical_method, num_iters_canonical, eps, verbosity,
                 rep_count);
}

void get_operand(Geometry &geom, const char &operand, const int &poly_size)
{
  string uniforms = "TCOID";

  if (uniforms.find(operand) != string::npos) {
    switch (operand) {
    case 'T':
      geom.read_resource("std_tet");
      break;

    case 'C':
      geom.read_resource("std_cube");
      break;

    case 'O':
      geom.read_resource("std_oct");
      break;

    case 'I':
      geom.read_resource("std_ico");
      break;

    case 'D':
      geom.read_resource("std_dod");
      break;
    }
  }
  else {
    Polygon pgon(poly_size, 1);

    switch (operand) {
    case 'P':
      pgon.set_type(Polygon::prism);
      break;

    case 'A':
      pgon.set_type(Polygon::antiprism);
      break;

    case 'Y':
      pgon.set_type(Polygon::pyramid);
      break;
    }

    pgon.set_edge(0, 1.0);

    if (operand == 'Y' && poly_size > 5)
      // Based on circumradius
      pgon.set_height(0, (1 / sin(M_PI / poly_size)) / 2);
    // inradius
    // poly->set_height((1/tan(M_PI/poly_size))/2);
    else
      pgon.set_edge(1, 1.0);

    pgon.make_poly(geom);
  }
}

void build_new_faces(map<string, map<string, string>> &faces_table,
                     map<string, int> verts_table,
                     vector<vector<int>> &faces_new)
{
  map<string, map<string, string>>::iterator ft;
  map<string, string>::iterator ftm;
  string face_name;
  for (ft = faces_table.begin(); ft != faces_table.end(); ft++) {
    for (ftm = ft->second.begin(); ftm != ft->second.end(); ftm++) {
      if (face_name != ft->first) {
        face_name = ft->first;
        string v0 = faces_table[face_name][ftm->first];
        string v = v0;
        vector<int> face;
        do {
          face.push_back(verts_table[v]);
          v = faces_table[face_name][v];
        } while (v != v0);
        if (face.size() > 2) // make sure face is valid
          faces_new.push_back(face);
        face.clear();
      }
    }
  }
}

void cn_ambo(Geometry &geom)
{
  vector<vector<int>> &faces = geom.raw_faces();
  vector<Vec3d> &verts = geom.raw_verts();

  map<string, int> verts_table;
  map<string, map<string, string>> faces_table;
  vector<Vec3d> verts_new;

  char buf1[MSG_SZ], buf2[MSG_SZ], buf3[MSG_SZ];
  unsigned int vert_num = 0;
  for (unsigned int i = 0; i < faces.size(); i++) {
    int v1 = faces[i].at(faces[i].size() - 2);
    int v2 = faces[i].at(faces[i].size() - 1);
    for (unsigned int j = 0; j < faces[i].size(); j++) {
      int v3 = faces[i].at(j);
      if (v1 < v2) {
        sprintf(buf1, "%d_%d", (v1 < v2) ? v1 : v2, (v1 > v2) ? v1 : v2);
        verts_table[buf1] = vert_num++;
        verts_new.push_back((verts[v1] + verts[v2]) * 0.5);
      }
      sprintf(buf1, "f%d", i);
      sprintf(buf2, "%d_%d", (v1 < v2) ? v1 : v2, (v1 > v2) ? v1 : v2);
      sprintf(buf3, "%d_%d", (v2 < v3) ? v2 : v3, (v2 > v3) ? v2 : v3);
      faces_table[buf1][buf2] = buf3;
      sprintf(buf1, "v%d", v2);
      sprintf(buf2, "%d_%d", (v2 < v3) ? v2 : v3, (v2 > v3) ? v2 : v3);
      sprintf(buf3, "%d_%d", (v1 < v2) ? v1 : v2, (v1 > v2) ? v1 : v2);
      faces_table[buf1][buf2] = buf3;
      v1 = v2;
      v2 = v3;
    }
  }

  geom.clear_all();
  verts = verts_new;
  verts_new.clear();

  build_new_faces(faces_table, verts_table, faces);
  faces_table.clear();
  verts_table.clear();

  geom.orient();
}

void cn_chamfer(Geometry &geom)
{
  vector<Vec3d> &verts = geom.raw_verts();
  int sz = verts.size();

  // make all edges a face
  // this is a join operation but retains the original vertex indexes
  make_edges_to_faces(&geom);

  project_onto_sphere(&geom);

  // truncate only on the new vertices
  vector<int> v_idxs;
  for (unsigned int i = sz; i < verts.size(); i++)
    v_idxs.push_back(i);

  truncate_verts(&geom, v_idxs, CN_ONE_HALF);

  geom.orient();
}

void cn_dual(Geometry &geom)
{
  Geometry dual;
  centroid_to_origin(geom);
  get_dual(&dual, geom, 1, Vec3d(0, 0, 0));
  geom = dual;
  geom.orient();
}

void cn_gyro(Geometry &geom)
{
  vector<vector<int>> &faces = geom.raw_faces();
  vector<Vec3d> &verts = geom.raw_verts();

  map<string, int> verts_table;
  map<string, map<string, string>> faces_table;
  vector<Vec3d> verts_new;

  char buf1[MSG_SZ], buf2[MSG_SZ], buf3[MSG_SZ];
  unsigned int vert_num = 0;
  vector<Vec3d> centers;
  geom.face_cents(centers);
  for (unsigned int i = 0; i < faces.size(); i++) {
    sprintf(buf1, "f%d", i);
    verts_table[buf1] = vert_num++;
    verts_new.push_back(centers[i].unit());
  }
  centers.clear();

  for (unsigned int i = 0; i < verts.size(); i++) {
    sprintf(buf1, "v%d", i);
    verts_table[buf1] = vert_num++;
    verts_new.push_back(verts[i]);
  }

  for (unsigned int i = 0; i < faces.size(); i++) {
    int v1 = faces[i].at(faces[i].size() - 2);
    int v2 = faces[i].at(faces[i].size() - 1);
    for (unsigned int j = 0; j < faces[i].size(); j++) {
      int v3 = faces[i].at(j);
      sprintf(buf1, "%d~%d", v1, v2);
      verts_table[buf1] = vert_num++;
      // approx. (2/3)v1 + (1/3)v2
      verts_new.push_back(verts[v1] * 0.7 + verts[v2] * 0.3);

      sprintf(buf1, "%df%d", i, v1);
      sprintf(buf2, "f%d", i);
      sprintf(buf3, "%d~%d", v1, v2);
      faces_table[buf1][buf2] = buf3;
      sprintf(buf2, "%d~%d", v1, v2);
      sprintf(buf3, "%d~%d", v2, v1);
      faces_table[buf1][buf2] = buf3;
      sprintf(buf2, "%d~%d", v2, v1);
      sprintf(buf3, "v%d", v2);
      faces_table[buf1][buf2] = buf3;
      sprintf(buf2, "v%d", v2);
      sprintf(buf3, "%d~%d", v2, v3);
      faces_table[buf1][buf2] = buf3;
      sprintf(buf2, "%d~%d", v2, v3);
      sprintf(buf3, "f%d", i);
      faces_table[buf1][buf2] = buf3;

      v1 = v2;
      v2 = v3;
    }
  }

  geom.clear_all();
  verts = verts_new;
  verts_new.clear();

  build_new_faces(faces_table, verts_table, faces);
  faces_table.clear();
  verts_table.clear();

  geom.orient();
}

void cn_kis(Geometry &geom, const int &n)
{
  vector<vector<int>> &faces = geom.raw_faces();
  vector<Vec3d> &verts = geom.raw_verts();

  vector<Vec3d> centers;
  geom.face_cents(centers);

  vector<vector<int>> faces_new;
  vector<int> face;
  for (unsigned int i = 0; i < faces.size(); i++) {
    if (n > 0 && (int)faces[i].size() != n) {
      faces_new.push_back(faces[i]);
      continue;
    }

    verts.push_back(centers[i]);
    for (unsigned int j = 0; j < faces[i].size() - 1; j++) {
      face.push_back(verts.size() - 1);
      face.push_back(faces[i][j]);
      face.push_back(faces[i][j + 1]);
      faces_new.push_back(face);
      face.clear();
    }

    face.push_back(verts.size() - 1);
    face.push_back(faces[i][faces[i].size() - 1]);
    face.push_back(faces[i][0]);
    faces_new.push_back(face);
    face.clear();
  }
  centers.clear();

  if (faces_new.size() > 0) {
    faces.clear();
    faces = faces_new;
    faces_new.clear();
  }

  geom.orient();
}

void cn_propellor(Geometry &geom)
{
  vector<vector<int>> &faces = geom.raw_faces();
  vector<Vec3d> &verts = geom.raw_verts();

  map<string, int> verts_table;
  map<string, map<string, string>> faces_table;
  vector<Vec3d> verts_new;

  char buf1[MSG_SZ], buf2[MSG_SZ], buf3[MSG_SZ];
  unsigned int vert_num = 0;
  for (unsigned int i = 0; i < verts.size(); i++) {
    sprintf(buf1, "v%d", i);
    verts_table[buf1] = vert_num++;
    verts_new.push_back(verts[i].unit());
  }

  for (unsigned int i = 0; i < faces.size(); i++) {
    int v1 = faces[i].at(faces[i].size() - 2);
    int v2 = faces[i].at(faces[i].size() - 1);
    for (unsigned int j = 0; j < faces[i].size(); j++) {
      int v3 = faces[i].at(j);
      sprintf(buf1, "%d~%d", v1, v2);
      verts_table[buf1] = vert_num++;
      // approx. (2/3)v1 + (1/3)v2
      verts_new.push_back(verts[v1] * 0.7 + verts[v2] * 0.3);

      sprintf(buf1, "v%d", i);
      sprintf(buf2, "%d~%d", v1, v2);
      sprintf(buf3, "%d~%d", v2, v3);
      faces_table[buf1][buf2] = buf3;
      sprintf(buf1, "%df%d", i, v2);
      sprintf(buf2, "%d~%d", v1, v2);
      sprintf(buf3, "%d~%d", v2, v1);
      faces_table[buf1][buf2] = buf3;
      sprintf(buf2, "%d~%d", v2, v1);
      sprintf(buf3, "v%d", v2);
      faces_table[buf1][buf2] = buf3;
      sprintf(buf2, "v%d", v2);
      sprintf(buf3, "%d~%d", v2, v3);
      faces_table[buf1][buf2] = buf3;
      sprintf(buf2, "%d~%d", v2, v3);
      sprintf(buf3, "%d~%d", v1, v2);
      faces_table[buf1][buf2] = buf3;

      v1 = v2;
      v2 = v3;
    }
  }

  geom.clear_all();
  verts = verts_new;
  verts_new.clear();

  build_new_faces(faces_table, verts_table, faces);
  faces_table.clear();
  verts_table.clear();

  geom.orient();
}

void cn_reflect(Geometry &geom) { geom.transform(Trans3d::inversion()); }

void cn_truncate_by_algorithm(Geometry &geom, const double &ratio, const int &n)
{
  truncate_verts(&geom, ratio, n);
  geom.orient();
}

void cn_whirl(Geometry &geom, const char &planarization_method,
              const int &num_iters_planar, const double &eps,
              const bool &verbosity, const int &rep_count)
{
  int num_faces = geom.raw_faces().size();

  verbose('g', 0, verbosity);
  cn_gyro(geom);
  cn_planarize(geom, planarization_method, num_iters_planar, eps, verbosity,
               rep_count);

  verbose('#', 0, verbosity);
  // only truncate on original face centers
  vector<int> v_idxs;
  for (int i = 0; i < num_faces; i++)
    v_idxs.push_back(i);

  truncate_verts(&geom, v_idxs, CN_ONE_HALF, nullptr);

  geom.orient();
}

void cn_expand(Geometry &geom, const bool &use_truncate_algorithm,
               const char &planarization_method, const int &num_iters_planar,
               const double &eps, const bool &verbosity, const int &rep_count)
{
  if (use_truncate_algorithm) {
    verbose('^', 0, verbosity);
    cn_truncate_by_algorithm(geom, CN_ONE_HALF, 0);
  }
  else {
    verbose('a', 0, verbosity);
    cn_ambo(geom);
  }
  cn_planarize(geom, planarization_method, num_iters_planar, eps, verbosity,
               rep_count);
  if (use_truncate_algorithm) {
    verbose('^', 0, verbosity);
    cn_truncate_by_algorithm(geom, CN_ONE_HALF, 0);
  }
  else {
    verbose('a', 0, verbosity);
    cn_ambo(geom);
  }
}

void cn_join(Geometry &geom, const bool &use_truncate_algorithm,
             const char &planarization_method, const int &num_iters_planar,
             const double &eps, const bool &verbosity, const int &rep_count)
{
  //   verbose('d',0,verbosity);
  //   cn_dual(geom);
  //   cn_planarize(geom, planarization_method, num_iters_planar, eps,
  //   verbosity, rep_count);
  if (use_truncate_algorithm) {
    verbose('^', 0, verbosity);
    cn_truncate_by_algorithm(geom, CN_ONE_HALF, 0);
  }
  else {
    verbose('a', 0, verbosity);
    cn_ambo(geom);
  }
  cn_planarize(geom, planarization_method, num_iters_planar, eps, verbosity,
               rep_count);
  verbose('d', 0, verbosity);
  cn_dual(geom);
}

void cn_meta(Geometry &geom, const bool &use_truncate_algorithm,
             const char &planarization_method, const int &num_iters_planar,
             const double &eps, const bool &verbosity, const int &rep_count)
{
  verbose('j', 0, verbosity);
  cn_join(geom, use_truncate_algorithm, planarization_method, num_iters_planar,
          eps, verbosity, rep_count);
  cn_planarize(geom, planarization_method, num_iters_planar, eps, verbosity,
               rep_count);
  verbose('k', 0, verbosity);
  cn_kis(geom, 0);
}

void cn_ortho(Geometry &geom, const bool &use_truncate_algorithm,
              const char &planarization_method, const int &num_iters_planar,
              const double &eps, const bool &verbosity, const int &rep_count)
{
  verbose('j', 0, verbosity);
  cn_join(geom, use_truncate_algorithm, planarization_method, num_iters_planar,
          eps, verbosity, rep_count);
  cn_planarize(geom, planarization_method, num_iters_planar, eps, verbosity,
               rep_count);
  verbose('j', 0, verbosity);
  cn_join(geom, use_truncate_algorithm, planarization_method, num_iters_planar,
          eps, verbosity, rep_count);
}

void cn_truncate(Geometry &geom, const int &n, const char &planarization_method,
                 const int &num_iters_planar, const double &eps,
                 const bool &verbosity, const int &rep_count)
{
  verbose('d', 0, verbosity);
  cn_dual(geom);
  cn_planarize(geom, planarization_method, num_iters_planar, eps, verbosity,
               rep_count);
  verbose('k', n, verbosity);
  cn_kis(geom, n);
  cn_planarize(geom, planarization_method, num_iters_planar, eps, verbosity,
               rep_count);
  verbose('d', 0, verbosity);
  cn_dual(geom);
}

void cn_snub(Geometry &geom, const char &planarization_method,
             const int &num_iters_planar, const double &eps,
             const bool &verbosity, const int &rep_count)
{
  //   verbose('d',0,verbosity);
  //   cn_dual(geom);
  //   cn_planarize(geom, planarization_method, num_iters_planar, eps,
  //   verbosity, rep_count);
  verbose('g', 0, verbosity);
  cn_gyro(geom);
  cn_planarize(geom, planarization_method, num_iters_planar, eps, verbosity,
               rep_count);
  verbose('d', 0, verbosity);
  cn_dual(geom);
}

void cn_bevel(Geometry &geom, const bool &use_truncate_algorithm,
              const char &planarization_method, const int &num_iters_planar,
              const double &eps, const bool &verbosity, const int &rep_count)
{
  if (use_truncate_algorithm) {
    verbose('^', 0, verbosity);
    cn_truncate_by_algorithm(geom, CN_ONE_HALF, 0);
  }
  else {
    verbose('a', 0, verbosity);
    cn_ambo(geom);
  }
  cn_planarize(geom, planarization_method, num_iters_planar, eps, verbosity,
               rep_count);
  if (use_truncate_algorithm) {
    verbose('#', 0, verbosity);
    cn_truncate_by_algorithm(geom, CN_ONE_THIRD, 0);
  }
  else {
    verbose('t', 0, verbosity);
    cn_truncate(geom, 0, planarization_method, num_iters_planar, eps, verbosity,
                rep_count);
  }
}

void do_operations(Geometry &geom, const vector<ops *> &operations,
                   const char &planarization_method,
                   const int &num_iters_planar, const double &eps,
                   const bool &use_truncate_algorithm, const bool &verbosity,
                   const int &rep_count)
{
  for (auto operation : operations) {
    switch (operation->op) {
    // "real" cases
    // ambo
    case 'a':
      if (use_truncate_algorithm) {
        verbose('^', operation->op_var, verbosity);
        cn_truncate_by_algorithm(geom, CN_ONE_HALF, operation->op_var);
      }
      else {
        verbose(operation->op, operation->op_var, verbosity);
        cn_ambo(geom);
      }
      break;

    // chamfer
    case 'c':
      verbose(operation->op, operation->op_var, verbosity);
      cn_chamfer(geom);
      break;

    // dual
    case 'd':
      verbose(operation->op, operation->op_var, verbosity);
      cn_dual(geom);
      break;

    // gyro
    case 'g':
      verbose(operation->op, operation->op_var, verbosity);
      cn_gyro(geom);
      break;

    // kis
    case 'k':
      verbose(operation->op, operation->op_var, verbosity);
      cn_kis(geom, operation->op_var);
      break;

    // propellor
    case 'p':
      verbose(operation->op, operation->op_var, verbosity);
      cn_propellor(geom);
      break;

    // reflect
    case 'r':
      verbose(operation->op, operation->op_var, verbosity);
      cn_reflect(geom);
      break;

    // whirl
    case 'w':
      verbose(operation->op, operation->op_var, verbosity);
      cn_whirl(geom, planarization_method, num_iters_planar, eps, verbosity,
               rep_count);
      break;

    // null
    case 'x':
      verbose(operation->op, operation->op_var, verbosity);
      break;

    // simulated cases for unreduced notation string
    // bevel
    case 'b':
      verbose(operation->op, operation->op_var, verbosity);
      cn_bevel(geom, use_truncate_algorithm, planarization_method,
               num_iters_planar, eps, verbosity, rep_count);
      break;

    // expand
    case 'e':
      verbose(operation->op, operation->op_var, verbosity);
      cn_expand(geom, use_truncate_algorithm, planarization_method,
                num_iters_planar, eps, verbosity, rep_count);
      break;

    // join
    case 'j':
      verbose(operation->op, operation->op_var, verbosity);
      cn_join(geom, use_truncate_algorithm, planarization_method,
              num_iters_planar, eps, verbosity, rep_count);
      break;

    // meta
    case 'm':
      verbose(operation->op, operation->op_var, verbosity);
      cn_meta(geom, use_truncate_algorithm, planarization_method,
              num_iters_planar, eps, verbosity, rep_count);
      break;

    // ortho
    case 'o':
      verbose(operation->op, operation->op_var, verbosity);
      cn_ortho(geom, use_truncate_algorithm, planarization_method,
               num_iters_planar, eps, verbosity, rep_count);
      break;

    // snub
    case 's':
      verbose(operation->op, operation->op_var, verbosity);
      cn_snub(geom, planarization_method, num_iters_planar, eps, verbosity,
              rep_count);
      break;

    // truncate. same as dk(n)d
    case 't':
      if (use_truncate_algorithm) {
        verbose('#', operation->op_var, verbosity);
        cn_truncate_by_algorithm(geom, CN_ONE_THIRD, operation->op_var);
      }
      else {
        verbose(operation->op, operation->op_var, verbosity);
        cn_truncate(geom, operation->op_var, planarization_method,
                    num_iters_planar, eps, verbosity, rep_count);
      }
      break;

    default:
      fprintf(stderr, "unexpected operator: '%c'\n", operation->op);
    }

    cn_planarize(geom, planarization_method, num_iters_planar, eps, verbosity,
                 rep_count);
    //      unitize_edges(geom);
  }
}

void cn_face_Coloring(Geometry &geom, const char &face_Coloring_method,
                      const ColorMapMulti &map, const int &face_opacity,
                      const string &face_pattern)
{
  if (face_Coloring_method == 'n') {
    const vector<vector<int>> &faces = geom.faces();
    for (unsigned int i = 0; i < faces.size(); i++) {
      int fsz = faces[i].size();
      Color col = map.get_col(fsz);
      if (col.is_val()) {
        int opq = (face_pattern[fsz % face_pattern.size()] == '1')
                      ? face_opacity
                      : col[3];
        col = Color(col[0], col[1], col[2], opq);
      }
      geom.colors(FACES).set(i, col);
    }
  }
  else if (face_Coloring_method == 's') {
    Symmetry sym;
    vector<vector<set<int>>> sym_equivs;
    sym.init(geom, &sym_equivs);

    Coloring clrng;
    clrng.add_cmap(map.clone());
    clrng.set_geom(&geom);
    clrng.f_sets(sym_equivs[2], true);

    // transparency
    if (face_opacity != 255) {
      for (unsigned int i = 0; i < geom.faces().size(); i++) {
        Color col = geom.colors(FACES).get(i);
        if (col.is_val())
          col = Color(col[0], col[1], col[2], face_opacity);
        geom.colors(FACES).set(i, col);
      }
    }
  }
}

int main(int argc, char *argv[])
{
  cn_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  if (opts.operand)
    get_operand(geom, opts.operand, opts.poly_size);
  else
    opts.read_or_error(geom, opts.ifile);

  // the program works better with oriented input, centroid at the origin
  geom.orient();
  if (!geom.is_oriented())
    opts.warning("input file contains a non-orientable geometry. output is "
                 "unpredictable");

  centroid_to_origin(geom);

  do_operations(geom, opts.operations, opts.planarization_method,
                opts.num_iters_planar, opts.epsilon,
                opts.use_truncate_algorithm, opts.verbosity, opts.rep_count);
  canonicalize(geom, opts.canonical_method, opts.num_iters_canonical,
               opts.epsilon, opts.verbosity, opts.rep_count);

  if (opts.unitize)
    unitize_edges(geom);

  if (opts.face_Coloring_method)
    cn_face_Coloring(geom, opts.face_Coloring_method, opts.map,
                     opts.face_opacity, opts.face_pattern);

  // color vertices and edges
  Coloring(&geom).vef_one_col(opts.vert_col, opts.edge_col, Color());

  opts.write_or_error(geom, opts.ofile);

  verbose('$', 0, opts.verbosity);

  return 0;
}
