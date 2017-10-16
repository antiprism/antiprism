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

struct ConwayOperator {
  std::string operator_short;
  std::string operator_name;
  bool allow_n;
  bool hart_operator;
};

// clang-format off
ConwayOperator conway_operator_list[]{
    {"a",  "ambo",          false, true  },
    //{"b3", "bevel3",        false, false }, // allows only 3
    {"b",  "bevel",         true,  false },
    {"c",  "chamfer",       false, true  },
    {"d",  "dual",          false, true  },
    {"e",  "expand",        false, false },
    {"g",  "gyro",          false, true  },
    {"j",  "join",          false, false },
    {"k",  "kis",           true,  true  },
    {"K",  "stake",         false, false },
    //{"L0", "joined-lace",   false, false }, // allows only explicit 0
    {"L",  "lace",          true,  false },
    {"l",  "loft",          false, false },
    //{"M0", "joined-medial", false, false }, // allows only explicit 0, 2 or larger
    //{"M3", "edge-medial-3", false, false },
    {"M",  "medial",        true,  false },
    //{"m3", "medial-3",      false, false }, // allows any N
    {"m",  "meta",          true,  false },
    {"n",  "needle",        false, false },
    //{"o3", "ortho3",        false, false }, // allows any N
    {"o",  "ortho",         true,  false },
    {"p",  "propellor",     false, true  },
    {"q",  "quinto",        false, false },
    {"r",  "reflect",       false, false },
    {"S",  "seed",          false, false },
    {"s",  "snub",          false, false },
    {"t",  "truncate",      true,  false },
    {"u",  "subdivide",     false, false },
    {"w",  "whirl",         false, true  },
    {"X",  "cross",         false, false },
    {"z",  "zip",           false, false },

};
// clang-format on

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

  // get allowable operators from table
  string operators;
  int last_op = sizeof(conway_operator_list) / sizeof(conway_operator_list[0]);
  for (int i = 0; i < last_op; i++)
    operators.push_back(conway_operator_list[i].operator_short[0]);

  // get operators which allow N from table
  string digits_allowed = "PAY";
  for (int i = 0; i < last_op; i++)
    if (conway_operator_list[i].allow_n)
      digits_allowed.push_back(conway_operator_list[i].operator_short[0]);

  string operands = "TCOIDPAY";
  string digits = "0123456789";

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

      // patch: for L standing alone to detect set num_val = -1
      if (current_op == 'L')
        num_val = -1;

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
      // will be at end of string as it is an operand
      else if (current_op == 'P' || current_op == 'A' || current_op == 'Y') {
        poly_size = atoi(number_string.c_str());
        if (poly_size < 3) {
          fprintf(stderr, "P(n), A(n), or Y(n), n must be 3 or greater\n");
          return i;
        }
      }
      // b allows only 3. 0 will not be passed to wythoff
      else if (current_op == 'b') {
        num_val = atoi(number_string.c_str());
        if (num_val != 0 && num_val != 3) {
          fprintf(stderr, "b(n), n must be 3\n");
          return i;
        }
        operations.push_back(new ops(i + 1, current_op, num_val));
        delayed_write = false;
      }
      // L allow 0 or 1. 0 will be made explicit in call to wythoff
      else if (current_op == 'L') {
        num_val = atoi(number_string.c_str());
        if (num_val != 0) {
          fprintf(stderr, "L(n), n must 0\n");
          return i;
        }
        operations.push_back(new ops(i + 1, current_op, num_val));
        delayed_write = false;
      }
      // M allow 0, 2 or greater. 0 will be made explicit in call to wythoff
      else if (current_op == 'M') {
        num_val = atoi(number_string.c_str());
        if (num_val == 1) {
          fprintf(stderr, "M(n), n must 0, 2 or greater\n");
          return i;
        }
        // patch: if num_val is 0, use -1 as a place holder
        //operations.push_back(new ops(i + 1, current_op, (num_val == 0) ? -1 : num_val));
        operations.push_back(new ops(i + 1, current_op, num_val));
        delayed_write = false;
      }
      // any other operator with N
      else if (digits_allowed.find(current_op) != string::npos) {
        num_val = atoi(number_string.c_str());
        // for other operators only set num_val if greater than 1
        if (num_val > 1)
          operations.push_back(new ops(i + 1, current_op, num_val));
        delayed_write = false;
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

// G. Hart Commentary
// P4 --> C    (C is prism)
// A3 --> O    (O is antiprism)
// Y3 --> T    (T is pyramid)
// e  --> aa   (abbr. for explode)
// b  --> ta   (abbr. for bevel)
// o  --> jj   (abbr. for ortho)
// m  --> kj   (abbr. for meta)
// t(n) --> dk(n)d  (dual operations) (see special case)
// j  --> dad  (dual operations)
// s  --> dgd  (dual operations)
// dd --> null (order 2)
// ad --> a    (a_ = ad_)
// gd --> g    (g_ = gd_)
// aY --> A    (interesting fact)
// dT --> T    (self-dual)
// gT --> D    (symm change)
// aT --> O    (symm change)
// dC --> O    (dual pair)
// dO --> C    (dual pair)
// dI --> D    (dual pair)
// dD --> I    (dual pair)
// aO --> aC   (for uniqueness)
// aI --> aD   (for uniqueness)
// gO --> gC   (for uniqueness)
// gI --> gD   (for uniqueness)

struct ResolveItem {
  const char *target;
  const char *resolve;
};

// clang-format off
ResolveItem resolve_item_list[] = {
    {  "P4", "C",   },
    {  "A3", "O",   },
    {  "Y3", "T",   },
    {  "e",  "aa",  },
    {  "b",  "ta",  },
    {  "o",  "jj",  },
    {  "m",  "kj",  },
    {  "t",  "dk",  },
    {  "j",  "dad", },
    {  "s",  "dgd", },
    {  "dd", "",    },
    {  "ad", "a",   },
    {  "gd", "g",   },
    {  "aY", "A",   },
    {  "dT", "T",   },
    {  "gT", "D",   },
    {  "aT", "O",   },
    {  "dC", "O",   },
    {  "dO", "C",   },
    {  "dI", "D",   },
    {  "dD", "I",   },
    {  "aO", "aC",  },
    {  "aI", "aD",  },
    {  "gO", "gC",  },
    {  "gI", "gD",  },
};
// clang-format on

string resolved_cn_string(const string &cn_string)
{
  string resolve_string = cn_string;

  int num_subst = sizeof(resolve_item_list) / sizeof(resolve_item_list[0]);

  string target;
  string resolve;

  // first 3 targets are positional
  if (resolve_string.length() > 1) {
    string tmp = resolve_string.substr(resolve_string.length() - 2, 2);
    for (unsigned int i = 0; i < 3; i++) {
      target = resolve_item_list[i].target;
      resolve = resolve_item_list[i].resolve;
      if (tmp == target)
        resolve_string.replace(resolve_string.length() - 2, target.length(),
                               resolve);
    }
  }

  // 4 to 6
  for (unsigned int i = 3; i < 7; i++) {
    target = resolve_item_list[i].target;
    resolve = resolve_item_list[i].resolve;
    for (size_t pos = resolve_string.find_first_of(target, 0);
         pos != string::npos; pos = resolve_string.find_first_of(target, 0)) {
      resolve_string.replace(pos, target.length(), resolve);
    }
  }

  // 7 is special case
  target = resolve_item_list[7].target;
  resolve = resolve_item_list[7].resolve;
  for (size_t pos = resolve_string.find_first_of(target, 0);
       pos != string::npos; pos = resolve_string.find_first_of(target, 0)) {
    resolve_string.replace(pos, target.length(), resolve);
    size_t pos2 = resolve_string.find_first_not_of("0123456789", pos + 2);
    if (pos2 != string::npos)
      resolve_string.insert(pos2, "d");
    else
      resolve_string.append("d");
  }

  // 8 to end. Because of length greater than 1 of the target, string.find and
  // replace cannot be used
  for (int i = 8; i < num_subst; i++) {
    target = resolve_item_list[i].target;
    resolve = resolve_item_list[i].resolve;
    int j = 0;
    int stop = resolve_string.length() - target.length() + 1;
    while (j <= stop) {
      string tmp = resolve_string.substr(j, target.length());
      if (tmp == target) {
        resolve_string.replace(j, target.length(), resolve);
        j = -1; // we have to begin looking from the begining once again, ok if
                // we modify stop
        stop = resolve_string.length() - target.length() + 1;
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
  bool resolve_ops;
  bool hart_mode;
  bool reverse_ops;
  char operand;
  int poly_size;
  char planarization_method;
  int num_iters_planar;
  int rep_count;
  bool unitize;
  bool verbosity;
  char face_coloring_method;
  int face_opacity;
  string face_pattern;

  double epsilon;

  Color vert_col;
  Color edge_col;

  ColorMapMulti map;

  vector<ops *> operations;

  cn_opts()
      : ProgramOpts("conway"), cn_string(""), resolve_ops(false),
        hart_mode(false), reverse_ops(false), operand('\0'), poly_size(0),
        planarization_method('p'), num_iters_planar(1000), rep_count(-1),
        unitize(false), verbosity(false), face_coloring_method('n'),
        face_opacity(-1), face_pattern("1"), epsilon(0),
        vert_col(Color(255, 215, 0)),  // gold
        edge_col(Color(211, 211, 211)) // lightgrey
  {
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

// clang-format off
void extended_help()
{
   fprintf(stdout,
"\n"
"The following is a description of Conway Notation edited from the Conway\n"
"Notation web page by George W. Hart (http://www.georgehart.com)\n"
"More detailed information and examples can be found at\n"
"http://www.georgehart.com/virtual-polyhedra/conway_notation.html\n"
"and at\n"
"http://en.wikipedia.org/wiki/Conway_polyhedron_notation\n"
"\n"
"Basics: In this notation, one specifies a \"seed\" polyhedron with a capital\n"
"letter. Operations to perform on any polyhedron are specified with lower-case\n"
"letters preceding it. This program contains a small set of seeds and operators\n"
"from which an infinite number of derived polyhedra can be generated.\n"
"Note: This C++ port of Conway Notation can also operate on OFF files from\n"
"standard input if the seed polyhedron is not specified.\n"
"\n"
"Seeds: The platonic solids are denoted T, O, C, I, and D, according to their\n"
"first letter. Other polyhedra which are implemented here include prisms, Pn,\n"
"antiprisms, An, and pyramids, Yn, where n is a number (3 or greater) which you\n"
"specify to indicate the size of the base you want, e.g., Y3=T, P4=C, and A3=O.\n"
"\n"
"Operations: Currently, abcdegjkmoprst are defined. They are motivated by the\n"
"operations needed to create the Archimedean solids and their duals from the\n"
"platonic solids.  Try each on a cube:\n"
"\n"
"a = ambo   The ambo operation can be thought of as truncating to the edge\n"
"midpoints.  It produces a polyhedron, aX, with one vertex for each edge of X.\n"
"There is one face for each face of X and one face for each vertex of X.\n"
"Notice that for any X, the vertices of aX are all 4-fold, and that aX=adX.\n"
"If two mutually dual polyhedra are in \"dual position,\" with all edges tangent\n"
"to a common sphere, the ambo of either is their intersection.  For example\n"
"aC=aO is the cuboctahedron.\n"
"Note: ambo is also known as \"rectifying\" the polyhedron, or rectification\n"
"\n"
"b = bevel  The bevel operation can be defined by bX=taX.  bC is the truncated\n"
"cuboctahedron.\n"
"Note: bevel is also known as \"omnitruncating\" the polyhedron, or omnitruncation\n"
"\n"
"c = chamfer   New hexagonal faces are added in place of edges (added)\n"
"\n"
"d = dual   The dual of a polyhedron has a vertex for each face, and a face for\n"
"each vertex, of the original polyhedron, e.g., dC=O.  Duality is an operation\n"
"of order two, meaning for any polyhedron X, ddX=X, e.g., ddC=dO=C.\n" 
"\n"
"e = expand This is Mrs. Stott's expansion operation.  Each face of X is\n"
"separated from all its neighbors and reconnected with a new 4-sided face,\n"
"corresponding to an edge of X.  An n-gon is then added to connect the 4-sided\n"
"faces at each n-fold vertex.  For example, eC is the rhombicuboctahedron.  It\n"
"turns out that eX=aaX and so eX=edX.\n"
"Note: expand is also known as \"cantellating\" the polyhedron, or cantellation\n"
"\n"
"g = gyro   The dual operation to s is g. gX=dsdX=dsX, with all 5-sided faces.\n"
"The gyrocube, gC=gO=\"pentagonal icositetrahedron,\" is dual to the snub cube.\n"
"g is like k but with the new edges connecting the face centers to the 1/3\n"
"points on the edges rather than the vertices.\n"
"\n"
"j = join   The join operator is dual to ambo, so jX=dadX=daX.  jX is like kX\n"
"without the original edges of X.  It produces a polyhedron with one 4-sided\n"
"face for each edge of X.  For example, jC=jO is the rhombic dodecahedron.\n"
"\n"
"k = kis    All faces are processed or kn = just n-sided faces are processed\n"
"The kis operation divides each n-sided face into n triangles.  A new vertex is\n"
"added in the center of each face, e.g., the kiscube, kC, has 24 triangular\n"
"faces.  The k operator is dual to t, meaning kX=dtdX.\n"
"\n"
"m = meta   Dual to b, mX=dbX=kjX.  mC has 48 triangular faces.  m is like k\n"
"and o combined; new edges connect new vertices at the face centers to the old\n"
"vertices and new vertices at the edge midpoints.  mX=mdX.  mC is the\n"
"\"hexakis octahedron.\"\n"
"\n"
"o = ortho  Dual to e, oX=deX=jjX.  oC is the trapezoidal icositetrahedron, with\n"
"24 kite-shaped faces.  oX has the effect of putting new vertices in the middle\n"
"of each face of X and connecting them, with new edges, to the edge midpoints of\n"
"X.\n"
"\n"
"p = propellor    Makes each n-gon face into a \"propellor\" of an n-gon\n"
"surrounded by n quadrilaterals, e.g., pT is the tetrahedrally stellated\n"
"icosahedron. Try pkD and pt6kT. p is a self-dual operation, i.e., dpdX=pX and\n"
"dpX=pdX, and p also commutes with a and j, i.e., paX=apX. (This and the next\n"
"are extensions were added by George Hart and not specified by Conway)\n"
"\n"
"r = reflect   Changes a left-handed solid to right handed, or vice versa, but\n"
"has no effect on a reflexible solid. So rC=C, but compare sC and rsC.\n"
"\n"
"s = snub   The snub operation produces the snub cube, sC, from C.  It can be\n"
"thought of as eC followed by the operation of slicing each of the new 4-fold\n"
"faces along a diagonal into two triangles.  With a consistent handedness to\n"
"these cuts, all the vertices of sX are 5-fold.  Note that sX=sdX.\n"
"\n"
"t = truncate  All faces are processed or tn = just n-sided faces are processed\n"
"Truncating a polyhedron cuts off each vertex, producing a new n-sided face for\n"
"each n-fold vertex.  The faces of the original polyhedron still appear, but\n"
"have twice as many sides, e.g., the tC has six octagonal sides corresponding to\n"
"the six squares of the C, and eight triangles corresponding to the cube's eight\n"
"vertices.\n"
"\n"
"w = whirl  Gyro followed by truncation of vertices centered on original faces.\n"
"This create 2 new hexagons for every original edge (added)\n"
"\n");
}

void cn_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [Conway Notation string] [input_file]\n"
"\n"
"Conway Notation uses algorithms by George W. Hart (http://www.georgehart.com)\n"
"http://www.georgehart.com/virtual-polyhedra/conway_notation.html\n"
"\n"
"Read a polyhedron from a file in OFF format.\n"
"If input_file is not given and no seed polyhedron is given in the notation\n"
"string then the program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"  -H        Conway Notation detailed help. seeds and operator descriptions\n"
"  -s        apply Conway Notation string substitutions\n"
"  -g        use George Hart algorithms (sets -s)\n"
"  -r        execute operations in reverse order (left to right)\n"
"  -u        make final product be averge unit edge length\n"
"  -v        verbose output\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"Planarization options (use canonical program to canonicalize output)\n"
"  -p <mthd> inter-step planarization method\n"
"            p - face centroids (magnitude squared) (default)\n"
"            m - mathematica planarize\n"
"            c - mathematica canonicalize\n"
"  -i <itrs> maximum inter-step planarization iterations (default: 1000)\n"
"  -z <n>    status reporting every n iterations, -1 for no status (default: -1)\n"
"  -l <lim>  minimum distance change to terminate planarization, as negative\n"
"               exponent (default: %d giving %.0e)\n"
"\n"
"\nColoring Options (run 'off_util -H color' for help on color formats)\n"
"  -V <col>  vertex color (default: gold)\n"
"  -E <col>  edge color   (default: lightgray)\n"
"  -f <mthd> mthd is face coloring method using color in map (default: n)\n"
"               key word: none - sets no color\n"
"               n - color by number of sides\n"
"               s - symmetric coloring\n"
"  -T <tran> face transparency. valid range from 0 (invisible) to 255 (opaque)\n"
"  -O <strg> face transparency pattern string (-f n only). valid values\n"
"               0 - map color alpha value, 1 -T alpha applied (default: '1')\n"
"  -m <maps> color maps for faces to be tried in turn (default: m1)\n"
"               keyword m1: red,darkorange1,yellow,darkgreen,cyan,blue,magenta,\n"
"                           white,grey,black\n"
"               keyword m2: red,blue,green,yellow,brown,magenta,purple,grue,\n"
"                           gray,orange (from George Hart\'s original applet)\n"
"\n"
"\n",prog_name(), help_ver_text, int(-log(::epsilon)/log(10) + 0.5), ::epsilon);
}
// clang-format on

void cn_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  int sig_compare = INT_MAX;

  string map_file;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hHsgruvp:l:i:z:f:V:E:T:O:m:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'H':
      extended_help();
      exit(0);

    case 's':
      resolve_ops = true;
      break;

    case 'g':
      hart_mode = true;
      resolve_ops = true;
      break;

    case 'r':
      reverse_ops = true;
      break;

    case 'u':
      unitize = true;
      break;

    case 'v':
      verbosity = true;
      break;

    case 'p':
      if (strlen(optarg) == 1 && strchr("pmc", int(*optarg)))
        planarization_method = *optarg;
      else
        error("planarize method type must be p, m or c", c);
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
        error("number of planar report iterations must be -1 or greater", c);
      break;

    case 'f':
      if (!strcasecmp(optarg, "none"))
        face_coloring_method = '\0';
      else if (strspn(optarg, "ns") != strlen(optarg) || strlen(optarg) > 1)
        error(msg_str("invalid face Coloring method '%s'", optarg), c);
      else {
        face_coloring_method = *optarg;
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

  //  if (strspn(cn_string.c_str(), "abcdegjkmopqrstwTCOIDPAY0123456789,") !=
  //      strlen(cn_string.c_str()))
  //    error("Conway Notation must consist of
  //    abcdegjkmopqrstwTCOIDPAY0123456789");

  if (int pos = validate_cn_string(cn_string, operations, operand, poly_size))
    error(msg_str("Unexpected character in position %d: %s", pos + 1,
                  cn_string.c_str()));

  if (resolve_ops) {
    cn_string = resolved_cn_string(cn_string);

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
  if (!reverse_ops)
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
      auto *col_map0 = new ColorMapMap;
      col_map0->set_col(0, Color(0.9, 0.3, 0.3));   // 3-sided faces red
      col_map0->set_col(1, Color(0.4, 0.4, 1.0));   // 4-sided faces blue
      col_map0->set_col(2, Color(0.2, 0.9, 0.3));   // 5-sided faces green
      col_map0->set_col(3, Color(0.9, 0.9, 0.2));   // 6-sided faces yellow
      col_map0->set_col(4, Color(0.5, 0.25, 0.25)); // 7-sided faces brown
      col_map0->set_col(5, Color(0.8, 0.2, 0.8));   // 8-sided faces magenta
      col_map0->set_col(6, Color(0.5, 0.2, 0.8));   // 9-sided faces purple
      col_map0->set_col(7, Color(0.1, 0.9, 0.9));   // 10-sided faces grue
      col_map0->set_col(8, Color(0.5, 0.5, 0.5));   // 11-sided faces gray
      col_map0->set_col(9, Color(1.0, 0.6, 0.1));   // 12-sided faces orange
      map.add_cmap(col_map0);

      // George Hart had all higher faces at grey
      col_map->set_col(0, Color(0.5, 0.5, 0.5)); // 13-sided faces and higher
    }
    col_map->set_wrap();
    map.add_cmap(col_map);
  }
  else
    print_status_or_exit(map.init(map_file.c_str()), 'm');

  epsilon = (sig_compare != INT_MAX) ? pow(10, -sig_compare) : ::epsilon;
}

void verbose(const char operation, const int op_var, const cn_opts &opts)
{
  if (opts.verbosity) {
    string operator_name;

    int last_op = sizeof(conway_operator_list) / sizeof(conway_operator_list[0]);
    for (int i = 0; i < last_op; i++) {
      if (operation == conway_operator_list[i].operator_short[0]) {
        operator_name = conway_operator_list[i].operator_name;
        break;
      }
    }

    string hart_operators = "acdgkpw";
    string hart_string;
    if (opts.hart_mode && (hart_operators.find(operation) != string::npos))
      hart_string = "(GHart)";

    // special cases
    char buf[MSG_SZ];
    buf[0] = '\0';
    if (operation == '_')
      operator_name = "planarizing ...";
    else if (operation == '$')
      operator_name = "done.";
    // L and op_var is -1 means L stands alone
    else if (operation == 'L' && op_var == -1) {
      buf[0] = '\0';
    }
    // if L or M have a 0
    else if ((operation == 'L' || operation == 'M') && (op_var == 0)) {
      if (operation == 'L')
        operator_name = "joined-lace";
      else if (operation == 'M')
        operator_name = "joined-medial";
    }
    // all other case show op_var when greater than 1
    else if (op_var > 1)
      sprintf(buf, "(%d)", op_var);

    fprintf(stderr, "%s%s %s\n", operator_name.c_str(), buf, hart_string.c_str());
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
  geom.transform(Trans3d::translate(-centroid(geom.verts())));
}

/*
// RK - average radius rather than maximum has more reliability than max
void unitize_vertex_radius(Geometry &geom)
{
  GeometryInfo info(geom);
  info.set_center(geom.centroid());
  //geom.transform(Trans3d::scale(1 / info.vert_dist_lims().max));
  double avg = info.vert_dist_lims().sum / info.num_verts();
  geom.transform(Trans3d::scale(1 / avg));
}
*/

void cn_planarize(Geometry &geom, const cn_opts &opts)
{
  if (opts.num_iters_planar != 0) {
    verbose('_', 0, opts);
    if (opts.planarization_method == 'p')
      planarize_bd(geom, opts.num_iters_planar, opts.rep_count, opts.epsilon);
    else if (opts.planarization_method == 'm')
      planarize_mm(geom, opts.num_iters_planar, opts.rep_count, opts.epsilon);
    else if (opts.planarization_method == 'c') {
      // RK - need?
      // unitize_vertex_radius(geom);
      // geom.transform(Trans3d::translate(-centroid(geom.verts())));
      canonicalize_mm(geom, opts.num_iters_planar, opts.rep_count,
                      opts.epsilon);
    }
  }
}

void get_operand(Geometry &geom, const char operand, const int poly_size)
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

// RK - for hart code
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

// hart_ code ported from George Hart java
void hart_ambo(Geometry &geom)
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

void hart_gyro(Geometry &geom)
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

void hart_kisN(Geometry &geom, const int n)
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

void hart_propellor(Geometry &geom)
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

// chamfer for hart code
void hart_chamfer(Geometry &geom, const cn_opts &opts)
{
  vector<Vec3d> &verts = geom.raw_verts();
  int sz = verts.size();

  // make all edges a face
  // this is a join operation but retains the original vertex indexes
  make_edges_to_faces(geom);

  project_onto_sphere(geom);

  // truncate only on the new vertices
  vector<int> v_idxs;
  for (unsigned int i = sz; i < verts.size(); i++)
    v_idxs.push_back(i);

  verbose('t', 0, opts);
  truncate_verts(geom, v_idxs, CN_ONE_HALF);

  geom.orient();
}

// chamfer for whirl code
void hart_whirl(Geometry &geom, const cn_opts &opts)
{
  int num_faces = geom.raw_faces().size();

  verbose('g', 0, opts);
  hart_gyro(geom);

  cn_planarize(geom, opts);

  // only truncate on original face centers
  vector<int> v_idxs;
  for (int i = 0; i < num_faces; i++)
    v_idxs.push_back(i);

  verbose('t', 0, opts);
  truncate_verts(geom, v_idxs, CN_ONE_HALF, nullptr);

  geom.orient();
}

// functions which can use Antiprism built in features
void hart_dual(Geometry &geom)
{
  Geometry dual;
  centroid_to_origin(geom);
  get_dual(dual, geom, 1, Vec3d(0, 0, 0));
  geom = dual;
  geom.orient();
}

void cn_reflect(Geometry &geom) { geom.transform(Trans3d::inversion()); }

void cn_truncate_by_algorithm(Geometry &geom, const double ratio, const int n)
{
  truncate_verts(geom, ratio, n);
  geom.orient();
}

void do_operations(Geometry &geom, const cn_opts &opts)
{
  for (auto operation : opts.operations) {
    verbose(operation->op, operation->op_var, opts);

    bool hart_operation_done = false;

    // reflection is universal
    if (opts.hart_mode || operation->op == 'r') {
      hart_operation_done = true;

      switch (operation->op) {
      // ambo
      case 'a':
        hart_ambo(geom);
        break;

      // chamfer
      case 'c':
        hart_chamfer(geom, opts);
        break;

      // dual
      case 'd':
        hart_dual(geom);
        break;

      // gyro
      case 'g':
        hart_gyro(geom);
        break;

      // kis
      case 'k':
        hart_kisN(geom, operation->op_var);
        break;

      // propellor
      case 'p':
        hart_propellor(geom);
        break;

      // reflect
      case 'r':
        cn_reflect(geom);
        break;

      // whirl
      case 'w':
        hart_whirl(geom, opts);
        break;

      default:
        hart_operation_done = false;
      }
    }

    // wythoff mode (with exceptions)
    if (!hart_operation_done) {
      // truncate with N>1 uses Hart algorithm
      if (operation->op == 't' && operation->op_var > 1)
        cn_truncate_by_algorithm(geom, CN_ONE_THIRD, operation->op_var);
      // kis with N>1 uses Hart algorithm
      else if (operation->op == 'k' && operation->op_var > 1)
        hart_kisN(geom, operation->op_var);
      // whirl uses Hart algorithm
      else if (operation->op == 'w')
        hart_whirl(geom, opts);
      else {
        // wythoff call requires a string
        string wythoff_op;
        wythoff_op.push_back(operation->op);

        char buf[MSG_SZ];
        if (operation->op_var > 1) {
          sprintf(buf, "%d", operation->op_var);
          wythoff_op += string(buf);
        }
        // L and M are special cases allowing explicit 0
        // if L stands alone op_var with be -1
        else if ((operation->op_var == 0) && operation->op == 'L') {
          sprintf(buf, "0");
          wythoff_op += string(buf);
        }
        else if ((operation->op_var == 0) && operation->op == 'M') {
          sprintf(buf, "0");
          wythoff_op += string(buf);
        }

//fprintf(stderr, "wythoff_op = %s\n", wythoff_op.c_str());
        opts.print_status_or_exit(
            wythoff_make_tiling(geom, geom, wythoff_op, true));

        // remove digon results
        vector<int> dels;
        for(int i=0; i<(int)geom.faces().size(); i++) {
          if(geom.faces(i).size() < 3)
            dels.push_back(i);
        }
        geom.del(FACES, dels); 
      }
      geom.orient();
    }

    // planarize after each step
    cn_planarize(geom, opts);
  }
}

void cn_face_coloring(Geometry &geom, const cn_opts &opts)
{
  if (opts.face_coloring_method == 'n') {
    bool trans_success = true;
    const vector<vector<int>> &faces = geom.faces();
    for (unsigned int i = 0; i < faces.size(); i++) {
      int fsz = faces[i].size();
      Color col = opts.map.get_col(fsz - 3);
      if (opts.face_opacity > -1) {
        int opq = (opts.face_pattern[fsz % opts.face_pattern.size()] == '1')
                      ? opts.face_opacity
                      : col[3];
        // map colors always are set but can be map index or invisible
        if (!col.set_alpha(opq))
          trans_success = false;
      }
      geom.colors(FACES).set(i, col);
    }
    if (!trans_success)
      fprintf(stderr, "warning: some faces could not be made transparent\n");
  }
  else if (opts.face_coloring_method == 's') {
    Symmetry sym;
    vector<vector<set<int>>> sym_equivs;
    sym.init(geom, &sym_equivs);

    Coloring clrng;
    clrng.add_cmap(opts.map.clone());
    clrng.set_geom(&geom);
    clrng.f_sets(sym_equivs[2], true);

    if (opts.face_opacity > -1) {
      ColorValuesToRangeHsva valmap(msg_str("A%g", opts.face_opacity / 255.0));
      valmap.apply(geom, FACES);

      for (const auto &kp : geom.colors(FACES).get_properties()) {
        if (kp.second.is_index()) {
          fprintf(stderr, "warning: map indexes cannot be made transparent\n");
          break;
        }
      }
    }
  }

  // check if some faces are not set for transparency warning
  if (opts.face_opacity > -1) {
    if (geom.colors(FACES).get_properties().size() < geom.faces().size())
      fprintf(stderr, "warning: unset faces cannot be made transparent\n");
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

  do_operations(geom, opts);

  if (opts.unitize)
    unitize_edges(geom);

  cn_face_coloring(geom, opts);

  // color vertices and edges
  Coloring(&geom).vef_one_col(opts.vert_col, opts.edge_col, Color());

  opts.write_or_error(geom, opts.ofile);

  verbose('$', 0, opts);

  return 0;
}
