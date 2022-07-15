/*
   Copyright (c) 2007-2021, Roger Kaufman

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
   Name: n_icons.cc
   Description: Creates Sphericon like Polyhedra. Also known as Streptohedra
   Project: Antiprism - http://www.antiprism.com
*/

#include "n_icons.h"
#include "../base/antiprism.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using std::make_pair;
using std::map;
using std::pair;
using std::set;
using std::sort;
using std::string;
using std::swap;
using std::vector;

using namespace anti;

// Need these here since used in opts
bool half_model(const vector<int> &longitudes)
{
  return (2 * longitudes.back() == longitudes.front());
}

bool full_model(const vector<int> &longitudes)
{
  return (longitudes.front() == longitudes.back());
}

// returns angle in range of 0 to 359.999...
double angle_in_range(double angle, const double eps)
{
  while (angle < 0.0)
    angle += 360.0;
  while (angle >= 360.0)
    angle -= 360.0;

  if (double_eq(angle, 0.0, eps) || double_eq(angle, 360.0, eps))
    angle = 0.0;

  return angle;
}

// angle is either point cut or side cut
bool angle_on_aligned_polygon(const double angle, const double n,
                              const double eps)
{
  double ang = angle_in_range(angle, eps);
  bool ret = double_eq(ang, 180.0, eps);
  if (!ret)
    ret = double_eq(fmod(ang, 180.0 / n), 0.0, eps);
  return ret;
}

ColorMapMap *alloc_default_map()
{
  auto *col_map = new ColorMapMap;

  col_map->set_col(0, Color(255, 0, 0));     // red
  col_map->set_col(1, Color(255, 127, 0));   // darkorange1
  col_map->set_col(2, Color(255, 255, 0));   // yellow
  col_map->set_col(3, Color(0, 100, 0));     // darkgreen
  col_map->set_col(4, Color(0, 255, 255));   // cyan
  col_map->set_col(5, Color(0, 0, 255));     // blue
  col_map->set_col(6, Color(255, 0, 255));   // magenta
  col_map->set_col(7, Color(255, 255, 255)); // white
  col_map->set_col(8, Color(127, 127, 127)); // gray50
  col_map->set_col(9, Color(0, 0, 0));       // black

  col_map->set_wrap();

  return col_map;
}

ColorMapMap *alloc_no_colorMap()
{
  auto *col_map = new ColorMapMap;

  col_map->set_col(0,
                   Color::maximum_index); // index INT_MAX as unset edge color
  // this doesn't work because 'no color' is not allowed in a map
  // col_map->set_col(0, Color());  // unset edge color

  col_map->set_wrap();

  return col_map;
}

class ncon_opts : public ProgramOpts {
public:
  string ofile;

  int ncon_order = 4;                // default n_icon n
  int d = 1;                         // default n_icon d for n/d
  int build_method = 0;              // methods are 1,2 and 3
  bool hide_indent = true;           // for n/d shells, show indented edges
  double inner_radius = NAN;         // for -z 2, overide calculated radius
  double outer_radius = NAN;         // for -z 2, overide calculated radius
  double angle = 0;                  // angle override for -z 3
  bool point_cut = true;             // actually means the dual
  bool hybrid = false;               // hybrid of the base and dual
  bool add_poles = false;            // add a vertex to the point cut model
  int twist = 1;                     // default twist of n_icon is 1
  bool add_symmetry_polygon = false; // add the symmetry polygon for -S coloring
  string closure;                    // if an open model, method of closure
  vector<int> longitudes; // how many longitudes and how many to display
  string hide_elems;      // elements to hide v,e,f or E edges of invisible f
  bool info = false;      // output info of n_icon

  // listing variables
  string ncon_surf;          // what type of surface to list (see help)
  vector<int> ncon_range;    // range of n to list
  bool long_form = false;    // long form of information
  bool filter_case2 = false; // filter out case 2 types

  double eps = anti::epsilon;

  // coloring
  char face_coloring_method = 'S';  // default color by the symmetry algorithm
  int face_opacity = -1;            // tranparency from 0 to 255
  string face_pattern = "1";        // for face transparency patterns
  char edge_coloring_method = '\0'; // edge coloring method is none
  int edge_opacity = -1;            // transparency from 0 to 255
  string edge_pattern = "1";        // for edge transparency patterns
  bool edge_set_no_color = false;   // to set no edge coloring
  Color unused_edge_color = Color::invisible; // edges not part of model
  bool symmetric_coloring = false;            // symmetry color for some modes
  int flood_fill_stop = 0;                    // flood fill early stop
  Color face_default_color = Color(192, 192, 192, 255); // darkgray
  Color edge_default_color = Color(192, 192, 192, 255); // darkgray

  ColorMapMulti face_map;
  ColorMapMulti edge_map;

  // common variables
  bool angle_is_side_cut = false; // if -z 3 angle is side cut
  bool double_sweep = false;      // double sweep is set in build_globe
  bool radius_inversion = false;  // if radius inversion occurs
  int mod_twist = 0;              // mod of abs(twist)

  // former global variable
  bool split = false; // an original variable for splitting caps

  ncon_opts() : ProgramOpts("n_icons") {}

  void process_command_line(int argc, char **argv);
  void usage();
};

void ncon_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options]

Creates Sphericon like Polyhedra. Also known as Streptohedra

Options
%s
  -I        information on current n-icon  
  -l <lim>  minimum distance for unique vertex locations as negative exponent
               (default: %d giving %.0e)
  -o <file> write output to file (default: write to standard output)

Program Options
  -n <n/d>  n-icon of order n. n must be 3 or greater (default: 4)
               use d to make star n-icon. d less than n
  -t <twst> number of twists. Can be negative, positive or 0 (default: 1)
  -s        side-cut of even order n-icon (default is point-cut)
  -H        hybrid of even order n-icon
  -a        angle (-z 3 only)
  -r        override inner radius (-z 2 only)
  -R        override outer radius (-z 2 only)
  -z <mthd> construction method
               1 - n/d must be co-prime. bow-ties can occur (default for d=1)
               2 - n/d compounds allowed. shell model (default for d>1)
               3 - n/d compounds allowed. No bow-ties (default if angle not 0)

Scene Options
  -M <m,m2> longitudes of model of m sides with optional m2 of m sides showing
               m may be odd, 3 or greater if twist is 0 (default: 36,36)
  -A        place a north and south pole in top and bottom if they exist
                only valid if m2<m. Not valid with -c h (-z 1 only)
  -c <clse> close open model if m2<m. Valid values h or v (-z 1,2)
               h = horizontal closure, v = vertical closure   
  -x <elms> v, e and f to remove OFF faces with one vertex (vertices),
               two-vertices (edges) and three or more vertices (faces)
               E - if face is invisible, associated edge is made invisible
  -Y        for n/d shells, when showing edges, show indented edges
  -W        add symmetry polygon (-f S or -e S)

Coloring Options (run 'off_util -H color' for help on color formats)
  -f <mthd> mthd is face coloring method. The coloring is done before twist
               keyword: none - sets no color
               S - color by symmetry polygon (default)
               s - color by circuits algorithm (n/d must be co-prime)
               f - color circuits with flood fill (-z 2,3 any n/d)
               c - color by compound
               a - color by compound, alternate method
               l - color latitudinally
               m - color longitudinally
               b - checkerboard with first two colors in face color list
               n - use each color in succession
               x - first two colors based on sign of x
               y - first two colors based on sign of y
               z - first two colors based on sign of z (z is the twist plane)
               o - use first eight colors per xyz octants
  -S        color circuits symmetrically for coloring method s,f,S (even n)
  -T <tran> face transparency. valid range from 0 (invisible) to 255 (opaque)
  -O <strg> face transparency pattern string. valid values
               0 -T value suppressed, 1 -T value applied  (default: '1')
  -e <mthd> mthd is edge coloring method. The coloring is done before twist
               keyword: none - sets no color
               keyword: Q - defer coloring all edges to option Q  (default)
                  or use the same letter options specified in -f, except c,a
               F - color edges with average adjoining face color
  -U <tran> edge transparency. valid range from 0 (invisible) to 255 (opaque)
  -P <strg> edge transparency pattern string. valid values
               0 -U value suppressed, 1 -U value applied  (default: '1')
  -Q <col>  color given to uncolored edges and vertices of final model
               keyword: none - sets no color (default: invisible)
  -D <c,e>  default color c for uncolored elements e (default: darkgray,ef)
               keyword: none - sets no color. elements e can include e or f
  -X <int>  flood fill stop. used with circuit or compound coloring (-f f,c)
               use 0 (default) to flood fill entire model. if -X is not 0 then
               return 1 from program if entire model has been colored
  -m <maps> color maps to be tried in turn. (default: map_red:darkorange1:
               yellow:darkgreen:cyan:blue:magenta:white:gray50:black%%)
               optionally followed by elements from v, e or f (default: vef)

Surface Count Reporting (options above ignored)
  -L <type> list n-icons with more than one surface. Valid values for type
               n = point cut even order n_icons
               s = side cut even order n-icons (surfaces > 2)
               o = odd order n_icons
               h = hybrids (all)
               i = hybrids (where N/2 is even)
               j = hybrids (where N/2 is odd)
               k = hybrids (where N/4 is even)
               l = hybrids (where N/4 is odd)
  -K <k,k2> range of n-icons to list for multiple surfaces
  -J        long form report
  -Z        filter out case 2 types

)",
          prog_name(), help_ver_text, int(-log(anti::epsilon) / log(10) + 0.5),
          anti::epsilon);
}

void ncon_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  Coloring clrngs[3];

  handle_long_opts(argc, argv);

  while ((c = getopt(
              argc, argv,
              ":hn:t:sHM:x:Ac:z:a:r:R:IL:K:JZm:f:ST:O:e:U:P:Q:YD:X:Wl:o:")) !=
         -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'n': {
      char *p;
      p = strchr(optarg, '/');
      if (p != nullptr) {
        *p++ = '\0';
        print_status_or_exit(read_int(p, &d), "n/d (d part)");
      }

      print_status_or_exit(read_int(optarg, &ncon_order), "n/d (n part)");
      if (ncon_order < 3)
        error("n must be an integer 3 or greater", "n/d (n part)");
      if (d < 1)
        error("d must be 1 or greater", "n/d (d part)");

      if (d >= ncon_order)
        error("d must be less than n", "n/d (d part)");
      break;
    }

    case 't':
      print_status_or_exit(read_int(optarg, &twist), c);
      break;

    case 's':
      point_cut = false;
      break;

    case 'H':
      hybrid = true;
      break;

    case 'M':
      print_status_or_exit(read_int_list(optarg, longitudes, true, 2), c);
      if (longitudes.front() < 3)
        error("m must be 3 or greater", c);
      if (longitudes.size() == 1)
        longitudes.push_back(longitudes.front());
      else if (longitudes.size() == 2) {
        if (longitudes.front() < longitudes.back())
          error("sides shown (m2) must be less than or equal to longitudes (m)",
                c);
        if (longitudes.back() <= 0)
          error("sides shown (m2) must be at least 1", c);
      }
      break;

    case 'x':
      if (strspn(optarg, "vefE") != strlen(optarg))
        error(msg_str("elements to hide are '%s' must be from "
                      "v, e, f and E",
                      optarg),
              c);
      hide_elems = optarg;
      break;

    case 'A':
      add_poles = true;
      break;

    case 'c':
      if (strspn(optarg, "hv") != strlen(optarg) || strlen(optarg) > 1)
        error(msg_str("closure is '%s', must be h or v (not both)", optarg), c);
      closure = optarg;
      break;

    case 'z':
      print_status_or_exit(read_int(optarg, &build_method), c);
      if (build_method < 1 || build_method > 3)
        error("method must be between 1 and 3", c);
      break;

    case 'a':
      print_status_or_exit(read_double(optarg, &angle), c);
      break;

    case 'r':
      print_status_or_exit(read_double(optarg, &inner_radius), c);
      // if (inner_radius <= 0.0)
      //   error("inner radius must be greater than 0", c);
      break;

    case 'R':
      print_status_or_exit(read_double(optarg, &outer_radius), c);
      // if (outer_radius <= 0.0)
      //   error("outer radius must be greater than 0", c);
      break;

    case 'I':
      info = true;
      break;

    case 'L':
      if (strspn(optarg, "nsohijkl") != strlen(optarg) || strlen(optarg) > 1)
        error(msg_str("n-icon type is '%s', must be only one of n, s, "
                      "o, h, i, j, k, or l\n",
                      optarg),
              c);
      ncon_surf = optarg;
      break;

    case 'K':
      print_status_or_exit(read_int_list(optarg, ncon_range, true, 2), c);
      if (ncon_range.front() < 3)
        error("k must be 3 or greater", c);
      if (ncon_range.size() == 1)
        ncon_range.push_back(ncon_range.front());
      else if ((ncon_range.size() == 2) &&
               (ncon_range.front() > ncon_range.back()))
        error("k2 shown must be greater than or equal to k", c);
      break;

    case 'J':
      long_form = true;
      break;

    case 'Z':
      filter_case2 = true;
      break;

    case 'm':
      print_status_or_exit(read_colorings(clrngs, optarg), c);
      break;

    case 'f':
      if (!strcasecmp(optarg, "none"))
        // place holder 'zero'
        face_coloring_method = '0';
      else if (strspn(optarg, "sflmbnxyzocaS") != strlen(optarg) ||
               strlen(optarg) > 1)
        error(msg_str("invalid face coloring method '%c'", *optarg), c);
      else
        face_coloring_method = *optarg;
      break;

    case 'S':
      symmetric_coloring = true;
      break;

    case 'T':
      print_status_or_exit(read_int(optarg, &face_opacity), c);
      if (face_opacity < 0 || face_opacity > 255) {
        error("face transparency must be between 0 and 255", c);
      }
      break;

    case 'O':
      if (strspn(optarg, "01") != strlen(optarg))
        error(msg_str("transparency string is '%s', must consist of "
                      "0 and 1's",
                      optarg),
              c);
      face_pattern = optarg;
      break;

    case 'e':
      if (!strcasecmp(optarg, "none")) {
        edge_set_no_color = true;
        edge_coloring_method = 'n';
      }
      else if (!strcmp(optarg, "Q"))
        edge_coloring_method = '\0';
      else if (strspn(optarg, "sflmbnxyzoFS") != strlen(optarg) ||
               strlen(optarg) > 1)
        error(msg_str("invalid edge coloring method '%s'", optarg), c);
      else
        edge_coloring_method = *optarg;
      break;

    case 'U':
      print_status_or_exit(read_int(optarg, &edge_opacity), c);
      if (edge_opacity < 0 || edge_opacity > 255) {
        error("edge transparency must be between 0 and 255", c);
      }
      break;

    case 'P':
      if (strspn(optarg, "01") != strlen(optarg))
        error(
            msg_str("transparency string %s must consist of 0 and 1's", optarg),
            c);
      edge_pattern = optarg;
      break;

    case 'Q':
      print_status_or_exit(unused_edge_color.read(optarg), c);
      break;

    case 'Y':
      hide_indent = false;
      break;

    case 'D': {
      Color tmp_color;

      bool set_edge = true;
      bool set_face = true;

      Split parts(optarg, ",");
      unsigned int parts_sz = parts.size();

      if (parts_sz > 2)
        error("excess entries for uncolored elements", c);

      for (unsigned int i = 0; i < parts_sz; i++) {
        // first argument is color
        if (i == 0) {
          print_status_or_exit(tmp_color.read(parts[i]), c);
        }
        else if (i == 1) {
          set_edge = false;
          set_face = false;

          if (strspn(parts[i], "ef") != strlen(parts[i]))
            error(msg_str("elements %s must be in e, f", optarg), c);

          if (strchr(parts[i], 'e'))
            set_edge = true;

          if (strchr(parts[i], 'f'))
            set_face = true;
        }
      }

      if (set_edge)
        edge_default_color = tmp_color;
      if (set_face)
        face_default_color = tmp_color;

      break;
    }

    case 'X':
      print_status_or_exit(read_int(optarg, &flood_fill_stop), c);
      if (flood_fill_stop < 0)
        error("flood fill stop must be 0 or greater", c);
      break;

    case 'W':
      add_symmetry_polygon = true;
      break;

    case 'l':
      int sig_compare;
      print_status_or_exit(read_int(optarg, &sig_compare), c);
      if (sig_compare > DEF_SIG_DGTS)
        warning("limit is very small, may not be attainable", c);
      eps = pow(10, -sig_compare);
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (argc - optind > 0)
    error("too many arguments");

  // surfaces subsystem
  if (ncon_surf.length() > 0) {
    if (ncon_range.size() == 0)
      error("for surfaces reporting -K must be specified", "J");
    else if (!is_even(ncon_range.front()) &&
             ncon_range.front() == ncon_range.back()) {
      if (strchr(ncon_surf.c_str(), 'h') || strchr(ncon_surf.c_str(), 'i') ||
          strchr(ncon_surf.c_str(), 'j'))
        error("for hybrid surfaces k must be even", "K");
      else if (strchr(ncon_surf.c_str(), 'n'))
        error("for even order n-icon surfaces k must be even", "K");
      else if (strchr(ncon_surf.c_str(), 's'))
        error("for side cut n-icon surfaces k must be even", "K");
    }
    else {
      if (is_even(ncon_range.front()) &&
          ncon_range.front() == ncon_range.back())
        if (strchr(ncon_surf.c_str(), 'o'))
          error("-K: for odd order n-icons surfaces k must be odd");
    }
  }
  // n_icons option processing
  else {
    // default longitudes to use is 36
    if (longitudes.size() == 0) {
      longitudes.push_back(36);
      longitudes.push_back(36);
    }

    // default build methods (old)
    if (!build_method) {
      // when n/1, build method 1, otherwise build method 2
      if (d == 1 || (ncon_order - d) == 1)
        build_method = 1;
      else
        build_method = 2;

      // radii setting only valid with method 2
      if (!std::isnan(inner_radius) || !std::isnan(outer_radius))
        build_method = 2;

      // angle is only valid with method 3
      if (angle)
        build_method = 3;

      // if flood fill or color by compound, need method 2 (or 3)
      if ((build_method == 1) && strchr("fc", face_coloring_method))
        build_method = 2;
    }

    if (build_method == 1) {
      if (gcd(ncon_order, d) != 1)
        error("when method = 1, n and d must be co-prime", "n");

      if (flood_fill_stop > 0)
        warning("flood fill stop has no effect in construction method 1", "X");
    }

    if (build_method == 2) {
      if (d == 1) {
        if (!std::isnan(inner_radius)) {
          if (hybrid) {
            if (is_even(twist))
              error("for method 2 and d=1, hybrid, inner radius cannot be set "
                    "when twist is even",
                    'r');
          }
          else if (!is_even(ncon_order))
            error("for method 2 and d=1, inner radius cannot be set when N is "
                  "odd",
                  'r');
          else if (!point_cut)
            error("for method 2 and d=1, inner radius cannot be set when side "
                  "cut",
                  'r');
          else if (!is_even(twist))
            error("for method 2 and d=1, inner radius cannot be set when twist "
                  "is odd",
                  'r');
        }
        else if (!std::isnan(outer_radius))
          inner_radius = outer_radius;
      }

      if (strchr(closure.c_str(), 'h')) {
        warning("closure of h not valid in construction method 2", "c");
        closure.clear();
      }
      else if (strchr(closure.c_str(), 'v')) {
        if (2 * longitudes.back() < longitudes.front()) {
          warning("closure of v not working for less than half models", "c");
          closure.clear();
        }
      }

      if (add_poles) {
        warning("poles not valid with construction method 2", "A");
        add_poles = false;
      }
    }
    else {
      if (!std::isnan(inner_radius)) {
        warning("inner radius is only valid in construction method 2", "r");
        inner_radius = NAN;
      }

      if (!std::isnan(outer_radius)) {
        warning("outer radius is only valid in construction method 2", "R");
        outer_radius = NAN;
      }

      if (!hide_indent) {
        warning("show indented edges only valid in construction method 2", "Y");
        hide_indent = true;
      }
    }

    if (build_method == 3) {
      if (angle && !point_cut)
        error("side cut is not valid when angle is not 0", "s");

      if (closure.length()) {
        warning("closure options not valid with construction method 3", "c");
        closure.clear();
      }

      if (add_poles) {
        warning("poles not valid with construction method 3", "A");
        add_poles = false;
      }

      // cannot have odd M with this method
      if (!is_even(longitudes.front()))
        error("-M m must be even with construction method 3", "M");
    }
    else {
      if (angle)
        error("angle is only valid in construction method 3", "a");
    }

    if (!hide_indent && !edge_coloring_method)
      warning("indented edges will not be shown unless an edge coloring is "
              "used (-e)",
              "Y");

    if (ncon_range.size() > 0)
      error("not valid without -J", "K");

    if (long_form)
      error("not valid without -J", "L");

    if (!point_cut) {
      if (!is_even(ncon_order))
        error("side cut is not valid with n-icons of odd n", "s");
    }

    if (hybrid) {
      if (!is_even(ncon_order))
        error("hybrids must be even order", "n");

      if (twist == 0)
        error("hybrids have no twist 0", "t");

      if (!point_cut && build_method != 3)
        error("hybrids and side cut can only be specified for method 3", 's');

      if (2 * longitudes.back() < longitudes.front()) {
        warning("for hybrids m2 cannot be less than half of full model. "
                "setting to half model",
                "M");
        longitudes.back() = longitudes.front() / 2;
      }

      if (closure.length()) {
        warning("closure options not valid with -H", "c");
        closure.clear();
      }

      if (add_poles) {
        warning("poles not valid with -H", "A");
        add_poles = false;
      }

      // This is true only of the face coloring. Edge coloring can be different
      // only in build method 1 it is possible for hybrids to be colored the
      // same using symmetric coloring
      // if (hybrid && symmetric_coloring && build_method == 1)
      //   warning("symmetric coloring is the same an non-symmetric coloring for
      //   hybrids","S");
    }

    // Let us allow globes (Twist = 0) to have uneven number of longitudes
    if ((twist != 0) && (!is_even(longitudes.front())))
      error("if twist -t is not 0 then -M m must be even", "M");

    // poles are only for odd order or even order side cuts
    if (is_even(ncon_order) && point_cut) {
      if (add_poles)
        warning("poles not valid with even point cut n-icons", "a");
      add_poles = false;
    }

    // Don't have poles on half or whole models - Cover longitudinal doesn't
    // need poles
    if ((half_model(longitudes) || full_model(longitudes)) ||
        strchr(closure.c_str(), 'h')) {
      if (add_poles)
        warning("poles not valid with half or full model or when -c h", "a");
      add_poles = false;
    }

    // Twist needs split. Can't have poles. (except for hybrids)
    if (twist != 0 && !hybrid) {
      split = true;
      if (add_poles)
        warning("poles not valid when twist is not 0", "a");
      add_poles = false;
    }

    // Full models have no need for closure
    if (full_model(longitudes)) {
      if (closure.length())
        warning("closure not needed for full model", "c");
      closure.clear();
    }
  }

  // color parameters, defaults...
  if (face_coloring_method == '0')
    face_coloring_method = '\0';
  else if (face_coloring_method == '\0') {
    face_coloring_method = 'S';
  }

  if (add_symmetry_polygon &&
      (face_coloring_method != 'S' && edge_coloring_method != 'S'))
    warning("adding symmetry polygon only has effect for -f S or -e S", 'W');

  // circuit table works with co-prime n/d
  // if n/d is co-prime (and d>1) angle must be 0
  if (face_coloring_method == 's') {
    if (gcd(ncon_order, d) != 1)
      error("circuit coloring only works n and d must be co-prime. Use 'S' or "
            "'f'",
            'f');
    if (d > 1 && (gcd(ncon_order, d) == 1) && double_ne(angle, 0.0, eps))
      error("When n/d is co-prime, angle must be 0. Use 'S' or 'f'", 'f');
  }
  else if (face_coloring_method == 'f') {
    if (build_method == 1)
      error("flood fill face coloring is for build method 2 or 3", 'f');
    if (build_method == 3 && (ncon_order == 2 * d))
      error("flood fill will not work in build method 3 and 2N/N polygons",
            'f');
  }
  else if (face_coloring_method == 'c') {
    if (build_method == 1)
      error("compound coloring is for build method 2 or 3", 'f');
  }

  if (edge_coloring_method == 'f') {
    if (build_method == 1)
      error("flood fill edge coloring is for build method 2 or 3", 'e');
  }

  if (symmetric_coloring && (!strchr("sfS", face_coloring_method) &&
                             !strchr("sfS", edge_coloring_method)))
    error("symmetric coloring is only for coloring methods s,f and S", "S");

  // only for odd order n_icons to be colored the same using symmetric coloring
  if (!is_even(ncon_order) && symmetric_coloring)
    warning("symmetric coloring is the same an non-symmetric coloring for odd "
            "order n-icons",
            "S");

  if ((clrngs[FACES].get_cmaps()).size())
    face_map = clrngs[FACES];
  else
    face_map.add_cmap(alloc_default_map());

  // patch for setting edges with no color
  if (edge_set_no_color)
    edge_map.add_cmap(alloc_no_colorMap());
  else
      // process as was done with face map
      if ((clrngs[1].get_cmaps()).size())
    edge_map = clrngs[EDGES];
  else
    edge_map.add_cmap(alloc_default_map());

  // vertex color map is the same as the edge color map
  // if((clrngs[VERTS].get_cmaps()).size())
  //   warning("vertex map has no effect","m");

  if (build_method == 3) {
    angle_is_side_cut =
        double_eq(fmod(angle_in_range(angle, eps), 360.0 / ncon_order),
                  (180.0 / ncon_order), eps);

    // if it is a point cut 2N/N, method 3 cannot be used
    bool pc =
        !angle_is_side_cut && angle_on_aligned_polygon(angle, ncon_order, eps);
    if (!point_cut)
      pc = false;
    if ((ncon_order == 2 * d) && pc)
      error("when polygon 2N/N and point cut, method 3 cannot be used", 'a');
  }

  // method 2: can't let d > n/2, causes problems
  if ((build_method == 2) && (d > ncon_order / 2))
    d = ncon_order - d;

  // for build method 2, multiply n and twist by 2
  // if d == 1, radii will equal
  if (build_method == 2 && d != 1) {
    ncon_order *= 2;
    twist *= 2;
  }

  mod_twist = std::abs(twist % ncon_order);
}

int longitudinal_faces(const int ncon_order, const bool point_cut)
{
  int lf = (int)floor((double)ncon_order / 2);
  if (is_even(ncon_order) && !point_cut)
    lf--;
  return (lf);
}

int num_lats(const int ncon_order, const bool point_cut)
{
  int lats = 0;
  if (is_even(ncon_order) && point_cut)
    lats = (ncon_order / 2);
  else if (is_even(ncon_order) && !point_cut)
    lats = (ncon_order / 2) + 1;
  else if (!is_even(ncon_order))
    lats = (int)ceil((double)ncon_order / 2);
  return lats;
}

int num_lats_faces(const vector<faceList *> &face_list)
{
  int lats = 0;
  for (auto i : face_list)
    if (i->lat > lats)
      lats = i->lat;
  // account for lat 0
  return lats + 1;
}

int num_lats_edges(const vector<edgeList *> &edge_list)
{
  int lats = 0;
  for (auto i : edge_list)
    if (i->lat > lats)
      lats = i->lat;
  return lats;
}

void add_coord(Geometry &geom, vector<coordList *> &coordinates,
               const Vec3d &vert)
{
  coordinates.push_back(new coordList(geom.add_vert(vert)));
}

void clear_coord(vector<coordList *> &coordinates)
{
  for (auto &coordinate : coordinates)
    delete coordinate;
  coordinates.clear();
}

void add_face(Geometry &geom, vector<faceList *> &face_list,
              const vector<int> &face, const int lat, const int lon)
{
  face_list.push_back(new faceList(geom.add_face(face), lat, lon, 0));
}

void add_face(Geometry &geom, vector<faceList *> &face_list,
              const vector<int> &face, const int lat, const int lon,
              const int polygon_no)
{
  face_list.push_back(new faceList(geom.add_face(face), lat, lon, polygon_no));
}

void clear_faces(vector<faceList *> &face_list)
{
  for (auto &i : face_list)
    delete i;
  face_list.clear();
}

void delete_face_list_items(vector<faceList *> &face_list,
                            const vector<int> &f_nos)
{
  vector<int> dels = f_nos;
  if (!dels.size())
    return;
  sort(dels.begin(), dels.end());
  unsigned int del_faces_cnt = 0;
  int map_to;
  for (unsigned int i = 0; i < face_list.size(); i++) {
    if (del_faces_cnt < dels.size() && (int)i == dels[del_faces_cnt]) {
      del_faces_cnt++;
    }
    else {
      map_to = i - del_faces_cnt;
      face_list[map_to] = face_list[i];
    }
  }
  face_list.resize(face_list.size() - del_faces_cnt);
}

// pass edge by value from make_edge()
// only add_edge_raw can be used else edge count is not correct for n_icons
void add_edge(Geometry &geom, vector<edgeList *> &edge_list,
              const vector<int> &edge, const int lat, const int lon)
{
  edge_list.push_back(new edgeList(geom.add_edge_raw(edge), lat, lon));
}

void clear_edges(vector<edgeList *> &edge_list)
{
  for (auto &i : edge_list)
    delete i;
  edge_list.clear();
}

void delete_edge_list_items(vector<edgeList *> &edge_list,
                            const vector<int> &f_nos)
{
  vector<int> dels = f_nos;
  if (!dels.size())
    return;
  sort(dels.begin(), dels.end());
  unsigned int del_edges_cnt = 0;
  int map_to;
  for (unsigned int i = 0; i < edge_list.size(); i++) {
    if (del_edges_cnt < dels.size() && (int)i == dels[del_edges_cnt]) {
      del_edges_cnt++;
    }
    else {
      map_to = i - del_edges_cnt;
      edge_list[map_to] = edge_list[i];
    }
  }
  edge_list.resize(edge_list.size() - del_edges_cnt);
}

class vertexMap {
public:
  int old_vertex;
  int new_vertex;
  vertexMap(int o, int n) : old_vertex(o), new_vertex(n) {}
};

void remap_elems(vector<vector<int>> &elems,
                 const vector<vertexMap> &coordinate_pairs)
{
  for (auto &elem : elems) {
    for (int &j : elem) {
      for (auto coordinate_pair : coordinate_pairs) {
        if (j == coordinate_pair.old_vertex)
          j = coordinate_pair.new_vertex;
      }
    }
  }
}

void merge_halves(Geometry &geom, vector<polarOrb *> &polar_orbit,
                  const double eps)
{
  const vector<Vec3d> &verts = geom.verts();

  vector<vertexMap> coordinate_pairs;

  for (unsigned int i = 0; i < polar_orbit.size(); i++) {
    int c1 = polar_orbit[i]->coord_no;
    for (unsigned int j = i + 1; j < polar_orbit.size(); j++) {
      int c2 = polar_orbit[j]->coord_no;
      if (!compare(verts[c1], verts[c2], eps)) {
        coordinate_pairs.push_back(vertexMap(c1, c2));
        break;
      }
    }
  }

  remap_elems(geom.raw_faces(), coordinate_pairs);
  remap_elems(geom.raw_edges(), coordinate_pairs);

  geom.del(VERTS, geom.get_info().get_free_verts());
}

// method 3: prime polygon is analog of prime meridian
void build_prime_polygon(Geometry &geom, vector<int> &prime_meridian,
                         vector<coordList *> &coordinates,
                         const vector<poleList *> &pole, const ncon_opts &opts)
{
  // for finding poles, the accuracy must be less than the default
  double epsilon_local = 1e-8;

  int num_polygons = gcd(opts.ncon_order, opts.d);
  int base_polygon =
      (num_polygons == 1) ? opts.ncon_order : opts.ncon_order / num_polygons;

  bool compound = (num_polygons == 1) ? false : true;

  double arc = 360.0 / base_polygon * (opts.d / num_polygons);
  double interior_angle = (180.0 - arc) / 2.0;
  double radius = sin(deg2rad(interior_angle)) / sin(deg2rad(arc));

  // patch for 2N/N
  if (double_eq(arc, 180.0, opts.eps))
    radius = sin(deg2rad(90.0)) / sin(deg2rad(90.0));

  double ang = opts.angle;
  ang -= 90.0;

  // side cut
  if (is_even(opts.ncon_order) && !opts.point_cut)
    ang += 180.0 / opts.ncon_order;

  // ang = angle_in_range(ang,eps);

  for (int i = 0; i < opts.ncon_order; i++) {
    prime_meridian.push_back(i);
    add_coord(
        geom, coordinates,
        Vec3d(cos(deg2rad(ang)) * radius, sin(deg2rad(ang)) * radius, 0.0));
    if (double_eq(fmod(angle_in_range(ang, opts.eps), 360.0), 90.0,
                  epsilon_local)) {
      pole[0]->idx = geom.verts().size() - 1;
      pole[0]->lat = 0;
    }
    else if (double_eq(fmod(angle_in_range(ang, opts.eps), 360.0), 270.0,
                       epsilon_local)) {
      pole[1]->idx = geom.verts().size() - 1;
      // all coloring and circuits considered make it point cut. make it true
      pole[1]->lat = num_lats(opts.ncon_order, true);
    }
    ang += arc;
    if (compound) {
      if (!((i + 1) % base_polygon))
        ang += 360.0 / opts.ncon_order;
    }
  }
}

// reverse polygon indexes of a polygon mirrored on Y
void reverse_poly_indexes_on_y(Geometry &geom, vector<int> &polygon,
                               const double eps)
{
  const vector<Vec3d> &verts = geom.verts();
  vector<bool> swapped(polygon.size());

  for (unsigned int i = 0; i < polygon.size() - 1; i++) {
    if (swapped[i])
      continue;
    for (unsigned int j = i; j < polygon.size(); j++) {
      if (swapped[j])
        continue;
      // if the points have equal Y and have the same fabs(X) then the indexes
      // are mirror/swapped on Y
      if (double_eq(verts[polygon[i]][1], verts[polygon[j]][1], eps) &&
          double_eq(fabs(verts[polygon[i]][0]), fabs(verts[polygon[j]][0]),
                    eps)) {
        swap(polygon[i], polygon[j]);
        swapped[i] = true;
        swapped[j] = true;
      }
    }
  }
}

// bypass is for testing. rotation will not work if true
vector<vector<int>> split_bow_ties(Geometry &geom,
                                   vector<coordList *> &coordinates,
                                   const vector<int> &face, const double eps)
{
  bool bypass = false;

  const vector<Vec3d> &verts = geom.verts();
  vector<vector<int>> faces;

  Vec3d intersection;
  if ((face.size() == 4) && !bypass) {
    for (unsigned int i = 0; i < 2; i++) {
      intersection = segments_intersection(verts[face[i]], verts[face[i + 1]],
                                           verts[face[i + 2]],
                                           verts[face[(i + 3) % 4]], eps);
      if (intersection.is_set()) {
        // make two points, move Z off z-plane plus and minus a little. to be
        // restored to zero later
        intersection[2] = eps * 2.0;
        int v_front = find_vert_by_coords(geom, intersection, eps);
        int v_back = v_front + 1;
        if (v_front == -1) {
          intersection[2] = eps * 2.0;
          add_coord(geom, coordinates, intersection);
          v_front = verts.size() - 1;
          intersection[2] = -eps * 2.0;
          add_coord(geom, coordinates, intersection);
          v_back = verts.size() - 1;
        }

        vector<int> triangle;
        triangle.push_back(face[1]);
        triangle.push_back(face[(i == 0) ? 2 : 0]);
        Vec3d ecent = centroid(verts, triangle);
        triangle.push_back((ecent[2] > 0.0) ? v_front : v_back);
        faces.push_back(triangle);

        triangle.clear();
        triangle.push_back(face[(i == 0) ? 0 : 2]);
        triangle.push_back(face[3]);
        ecent = centroid(verts, triangle);
        triangle.push_back((ecent[2] > 0.0) ? v_front : v_back);
        faces.push_back(triangle);

        break;
      }
    }
  }

  // if no faces were formed just return original face
  if (!faces.size())
    faces.push_back(face);

  return faces;
}

vector<int> find_face_by_lat_lon(const vector<faceList *> &face_list,
                                 const int lat, const int lon)
{
  vector<int> idx;
  for (unsigned int i = 0; i < face_list.size(); i++) {
    if (face_list[i]->lat == lat && face_list[i]->lon == lon)
      idx.push_back(i);
  }
  return idx;
}

vector<int> find_edge_by_lat_lon(const vector<edgeList *> &edge_list,
                                 const int lat, const int lon)
{
  vector<int> idx;
  for (unsigned int i = 0; i < edge_list.size(); i++) {
    if (edge_list[i]->lat == lat && edge_list[i]->lon == lon)
      idx.push_back(i);
  }
  return idx;
}

bool add_edge_wrapper(Geometry &geom, vector<edgeList *> &edge_list,
                      const vector<int> &edge, const int lat,
                      const int lon_front, const int lon_back)
{
  const vector<Vec3d> &verts = geom.verts();

  int edge_no = find_edge_in_edge_list(geom.edges(), edge);
  if (edge_no < 0) {
    double edge_z = centroid(verts, edge)[2];

    int lon = (edge_z > 0.0) ? lon_front : lon_back;
    add_edge(geom, edge_list, make_edge(edge[0], edge[1]), lat, lon);

    if (edge_z > 0.0) // || double_eq(edge_z,0.0,eps))
      edge_list.back()->rotate = true;
  }

  return (edge_no < 0 ? false : true);
}

// split_face_indexes is cleared after use
void apply_latitudes(const Geometry &geom,
                     vector<vector<int>> &split_face_indexes,
                     const vector<faceList *> &face_list,
                     const vector<edgeList *> &edge_list,
                     const vector<poleList *> &pole, const ncon_opts &opts)
{
  const vector<vector<int>> &faces = geom.faces();
  const vector<vector<int>> &edges = geom.edges();
  const vector<Vec3d> &verts = geom.verts();

  int n = opts.ncon_order;
  if (opts.build_method == 2)
    n /= 2;

  // save edge latitudes for restoring them in some cases
  vector<int> edge_lats_save;
  if (opts.build_method == 2) {
    for (unsigned int i = 0; i < edge_list.size(); i++) {
      int j = edge_list[i]->edge_no;
      edge_lats_save.push_back(edge_list[j]->lat);
    }
  }

  // collect Y value of edges and map them
  map<int, vector<int>> levels_edges;
  if (opts.build_method == 2) {
    for (auto i : edge_list) {
      int j = i->edge_no;
      if (i->lat < 0)
        continue;
      int level = i->lat - 1;
      levels_edges[level].push_back(j);
    }

    // if indented edges are used there will be gaps
    map<int, vector<int>> levels_edges_r;
    if (!opts.hide_indent) {
      unsigned int level_l = 0;
      unsigned int level_r = 0;
      while (level_l < levels_edges.size()) {
        if (levels_edges[level_l].size()) {
          levels_edges_r[level_r] = levels_edges[level_l];
          level_r++;
        }
        level_l++;
      }
      // levels_edges.clear();
      levels_edges = levels_edges_r;
    }

    // clear face latitudes
    for (auto i : face_list)
      i->lat = -1;
    // clear edge latitudes
    for (auto i : edge_list) {
      if (i->lat > -1)
        i->lat = -1;
    }
  }
  else if (opts.build_method == 3) {
    vector<pair<double, int>> edge_ys;
    for (auto i : edge_list) {
      int j = i->edge_no;
      double y = verts[edges[j][0]][1];
      edge_ys.push_back(make_pair(y, j));
    }

    sort(edge_ys.begin(), edge_ys.end());
    reverse(edge_ys.begin(), edge_ys.end());

    double last_y = std::numeric_limits<double>::max();
    int level = -1;
    for (auto &edge_y : edge_ys) {
      if (double_ne(edge_y.first, last_y, opts.eps))
        level++;
      last_y = edge_y.first;
      levels_edges[level].push_back(edge_y.second);
    }
  }

  bool do_faces = !(!opts.face_coloring_method ||
                    !strchr("slbfc", opts.face_coloring_method));
  bool do_edges =
      !(!opts.edge_coloring_method || !strchr("l", opts.edge_coloring_method));
  if (do_faces || do_edges) {
    // put split faces into a map
    map<int, int> split_face_map;
    for (auto &split_face_index : split_face_indexes) {
      if (split_face_index.size() == 1) {
        split_face_map[split_face_index[0]] = split_face_index[0];
      }
      else {
        split_face_map[split_face_index[0]] = split_face_index[1];
        split_face_map[split_face_index[1]] = split_face_index[0];
      }
    }
    // split_face_indexes no longer needed
    split_face_indexes.clear();

    // each face is associated with one or two edges
    map<int, vector<int>> faces_edges_map;
    for (unsigned int i = 0; i < faces.size(); i++) {
      vector<int> face = faces[i];
      unsigned int sz = face.size();
      for (unsigned int j = 0; j < sz; j++) {
        vector<int> edge(2);
        edge[0] = face[j];
        edge[1] = face[(j + 1) % sz];
        int ret = find_edge_in_edge_list(edges, edge);
        if (ret > -1)
          faces_edges_map[i].push_back(ret);
      }
    }

    // each edge is associated with two faces
    map<int, vector<int>> edge_faces_map;
    for (unsigned int i = 0; i < faces_edges_map.size(); i++) {
      vector<int> face_idx = faces_edges_map[i];
      for (int k : face_idx) {
        edge_faces_map[k].push_back(i);
      }
    }

    vector<int> last_edges;
    int first_level = -1;
    int lat = 1;
    bool early_ending = false;
    bool early_mode = false;

    int part_number = 1;
    while (part_number) {
      // if south pole is included, colors are backward
      // if ((part_number == 1) && ((pole[0]->lat > -1) || (pole[1]->lat > -1)))
      // {
      //   int v_idx = (pole[0]->lat > -1) ? pole[0]->idx : pole[1]->idx;
      if ((part_number == 1) && (pole[0]->lat > -1)) {
        int v_idx = pole[0]->idx;
        vector<int> faces_with_index = find_faces_with_vertex(faces, v_idx);
        for (int j : faces_with_index) {
          face_list[j]->lat = lat - 1;
          int k = split_face_map[j];
          face_list[k]->lat = lat - 1;
          vector<int> edge_idx = faces_edges_map[k];
          for (int m : edge_idx) {
            if (edge_list[m]->lat == -1) {
              edge_list[m]->lat = lat;
              last_edges.push_back(m);
            }
          }
        }
      }
      // else need to prime side cut for level 0
      else {
        // find first level
        if (part_number == 1) {
          // if n is even, simple math
          if (is_even(n)) {
            first_level = opts.d / 2;
            if (opts.double_sweep)
              first_level *= 2;
          }
          // if n is odd, use vertex number to find level
          // first 'flat' level from the bottom
          else {
            if (opts.d == 1)
              first_level = 0;
            else {
              if (is_even(opts.d))
                first_level = (n / 2) - (opts.d / 2);
              else
                first_level = opts.d / 2;
            }
          }
        }

        // set latitude on first edge level
        for (int j : levels_edges[first_level]) {
          edge_list[j]->lat = lat;
          last_edges.push_back(j);
        }

        // prime system to set latitudes on first set of faces
        for (int &last_edge : last_edges) {
          vector<int> face_idx = edge_faces_map[last_edge];

          // find higher face in terms of Y
          if (centroid(verts, faces[face_idx[1]])[1] >
              centroid(verts, faces[face_idx[0]])[1])
            swap(face_idx[0], face_idx[1]);
          // exception: if odd n, paint upper faces when y level in positive y
          double first_level_y =
              verts[edges[levels_edges[first_level][0]][0]][1];
          if (!is_even(n) && (opts.d > 1) &&
              double_gt(first_level_y, 0, opts.eps))
            swap(face_idx[0], face_idx[1]);

          int face_painted = -1;
          double edge_y = (verts[edges[last_edge][0]])[1];
          for (unsigned int j = 0; j < face_idx.size(); j++) {
            double face_cent_y = centroid(verts, faces[face_idx[j]])[1];
            if (double_eq(edge_y, face_cent_y, opts.eps)) {
              face_painted = j;
              break;
            }
          }

          // if not found: even n, face 0; odd n, face 1
          if (face_painted == -1)
            face_painted = (is_even(n)) ? 0 : 1;

          int para = 0;
          if (part_number == 1)
            para = 0;
          else
            para = lat;

          if (early_mode)
            para--;

          int j = face_idx[face_painted];
          face_list[j]->lat = para;
          int k = split_face_map[j];
          face_list[k]->lat = para;
        }
      }

      // loop to set rest of latitudes
      while (last_edges.size()) {
        vector<int> next_edges;
        for (unsigned int i = 0; i < last_edges.size(); i++) {
          vector<int> face_idx = edge_faces_map[last_edges[i]];

          // test for early ending...
          if (i == 0) {
            early_ending = true;
            for (int k : face_idx) {
              if (face_list[k]->lat == -1)
                early_ending = false;
            }
          }
          if (early_ending)
            break;

          for (int k : face_idx) {
            if (face_list[k]->lat != -1)
              continue;

            int para = 0;
            if (part_number == 1)
              para = lat;
            else
              para = lat + 1;

            if (early_mode)
              para--;

            face_list[k]->lat = para;
            int l = split_face_map[k];
            face_list[l]->lat = para;
            vector<int> edges_idx = faces_edges_map[l];
            for (int n : edges_idx) {
              if (edge_list[n]->lat == -1)
                next_edges.push_back(n);
            }
          }
        }

        lat++;

        for (int next_edge : next_edges)
          edge_list[next_edge]->lat = lat;
        last_edges = next_edges;
      }

      // check for any more unset edges
      bool found = false;
      for (unsigned int i = 0; i < levels_edges.size(); i++) {
        for (unsigned int j = 0; j < levels_edges[i].size(); j++) {
          int k = levels_edges[i][j];
          if (edge_list[k]->lat == -1) {
            found = true;
            first_level = i;
            break;
          }
        }
        // comment out to have prime circuits of parts>1 upside down
        if (found)
          break;
      }

      if (early_ending)
        early_mode = true;
      early_ending = false;

      if (found)
        part_number++;
      else
        part_number = 0;
    }
  }

  // restore edges so circuit table works
  bool reset_done = false;
  if (strchr("sf", opts.edge_coloring_method)) {
    if (opts.build_method == 2) {
      reset_done = true;
      for (unsigned int i = 0; i < edge_lats_save.size(); i++) {
        int lat = edge_lats_save[i];
        if (lat < -1)
          lat = (opts.hide_indent) ? -1 : std::abs(lat + 2);
        edge_list[i]->lat = lat;
      }
    }
    else if (opts.build_method == 3) {
      for (unsigned int i = 0; i < levels_edges.size(); i++) {
        for (unsigned int j = 0; j < levels_edges[i].size(); j++) {
          int k = levels_edges[i][j];
          edge_list[k]->lat = i + 1;
        }
      }
    }
  }

  // marked indented latitudes reset
  if (opts.build_method == 2 && !reset_done) {
    for (unsigned int i = 0; i < edge_lats_save.size(); i++) {
      int lat = edge_lats_save[i];
      if (lat < -1) {
        lat = (opts.hide_indent) ? -1 : std::abs(lat + 2);
        edge_list[i]->lat = lat;
      }
    }
  }
}

// This was the old method of apply_latitudes, it still works for double_sweep
// and d=1
// for method 3: set latitude numbers based on height of Y
// set face latitudes based on edge latitudes
void apply_latitudes(const Geometry &geom,
                     const vector<vector<int>> &original_faces,
                     const vector<vector<int>> &split_face_indexes,
                     const vector<faceList *> &face_list,
                     const vector<edgeList *> &edge_list,
                     const vector<poleList *> &pole, const ncon_opts &opts)
{
  const vector<vector<int>> &faces = geom.faces();
  const vector<vector<int>> &edges = geom.edges();
  const vector<Vec3d> &verts = geom.verts();

  // possible future use when not coloring compounds
  // if false give each latitude unique number, if true pair each layer split by
  // angle
  bool double_sw; // = opts.double_sweep;
  double_sw = false;

  // possible future use
  // hybrid patch: for a hybrid when angle causes a double sweep, one side will
  // be off by 1
  bool hybrid_patch =
      (opts.hybrid && opts.point_cut && double_sw) ? true : false;

  // do edges first and use them as reference latitudes for faces
  // collect Y value of edges and sort them
  vector<pair<double, int>> edge_ys;
  for (unsigned int i = 0; i < edges.size(); i++)
    edge_ys.push_back(make_pair(verts[edges[i][0]][1], i));

  sort(edge_ys.begin(), edge_ys.end());
  reverse(edge_ys.begin(), edge_ys.end());

  int lat = 1;
  int wait = (hybrid_patch) ? 0 : 1;

  double last_y = edge_ys[0].first;
  for (auto &edge_y : edge_ys) {
    if (double_ne(edge_y.first, last_y, opts.eps)) {
      lat++;
      if (double_sw) {
        lat -= wait;
        wait = (wait) ? 0 : 1;
      }
    }
    int j = edge_y.second;
    edge_list[j]->lat = lat;
    last_y = edge_y.first;
  }

  // find faces connected to edges for latitude assignment
  // optimization: if face lats not needed then bypass
  bool do_faces = !(!opts.face_coloring_method ||
                    !strchr("slbfc", opts.face_coloring_method));
  if (!do_faces)
    return;

  // collect minimum Y values of faces
  vector<pair<double, int>> face_ys;
  for (unsigned int i = 0; i < original_faces.size(); i++) {
    vector<int> face = original_faces[i];
    double min_y = std::numeric_limits<double>::max();
    for (int j : face) {
      double y = verts[j][1];
      if (y < min_y)
        min_y = y;
    }
    face_ys.push_back(make_pair(min_y, i));
  }

  lat = 1;
  int max_lat_used = 0;

  // if there is a north pole
  if (pole[0]->idx != -1) {
    Vec3d np = verts[pole[0]->idx];
    for (unsigned int i = 0; i < faces.size(); i++) {
      vector<int> face = faces[i];
      for (int j : face) {
        Vec3d v = verts[j];
        if (!compare(v, np, opts.eps)) {
          face_list[i]->lat = 0;
          break;
        }
      }
    }
    lat++;
  }

  wait = (double_sw) ? 2 : 0;

  last_y = edge_ys[0].first;
  for (auto &edge_y : edge_ys) {
    if (double_ne(edge_y.first, last_y, opts.eps)) {
      lat = max_lat_used + 2;
      if (double_sw) {
        lat -= wait;
        wait = (wait) ? 0 : 2;
      }
      last_y = edge_y.first;
    }

    int j = edge_y.second;
    vector<int> face_idx = find_faces_with_edge(original_faces, edges[j]);
    // edge must have 2 faces
    if (face_idx.size() != 2)
      continue;

    // find higher face in terms of Y
    if (face_ys[face_idx[1]].first > face_ys[face_idx[0]].first)
      swap(face_idx[0], face_idx[1]);

    bool first_lat_used = true;
    for (unsigned int k = 0; k < 2; k++) {
      int l = face_idx[k];
      for (unsigned int m = 0; m < split_face_indexes[l].size(); m++) {
        int n = split_face_indexes[l][m];
        if (face_list[n]->lat == -1) {
          int lat_used = (!k || (k && !first_lat_used)) ? lat - 1 : lat;
          if (lat_used > max_lat_used)
            max_lat_used = lat_used;
          face_list[n]->lat = lat_used;
        }
        else {
          if (!k && (face_list[n]->lat != lat - 1))
            first_lat_used = false;
        }
      }
    }
  }
}

// method 3: fix polygon numbers for compound coloring
void fix_polygon_numbers(const vector<faceList *> &face_list,
                         const ncon_opts &opts)
{
  int lat = 0;
  unsigned int sz = 0;
  do {
    int polygon_min = std::numeric_limits<int>::max();
    for (unsigned int j = 0; j < 2; j++) {
      int lon = (opts.longitudes.front() / 2) - j;
      vector<int> idx = find_face_by_lat_lon(face_list, lat, lon);
      sz = idx.size();
      for (unsigned int k = 0; k < sz; k++) {
        int polygon_no = face_list[idx[k]]->polygon_no;
        if (polygon_no < polygon_min)
          polygon_min = polygon_no;
      }
    }

    if (polygon_min != std::numeric_limits<int>::max()) {
      for (auto j : face_list) {
        if (j->lat == lat)
          j->polygon_no = polygon_min;
      }
    }
    lat++;
  } while (sz);
}

// for method 3: analog to form_globe()
// maximum latitudes is set
void form_angular_model(Geometry &geom, const vector<int> &prime_meridian,
                        vector<coordList *> &coordinates,
                        vector<faceList *> &face_list,
                        vector<edgeList *> &edge_list,
                        const vector<poleList *> &pole,
                        vector<vector<int>> &original_faces,
                        vector<vector<int>> &split_face_indexes,
                        const int polygons_total, const ncon_opts &opts)
{
  const vector<Vec3d> &verts = geom.verts();

  double arc = 360.0 / opts.longitudes.front();

  int num_polygons = gcd(opts.ncon_order, opts.d);
  int base_polygon =
      (num_polygons == 1) ? opts.ncon_order : opts.ncon_order / num_polygons;

  vector<int> meridian_last;
  vector<int> meridian;

  for (int i = 1; i <= polygons_total; i++) {
    // move current meridian one back
    meridian_last = (i == 1) ? prime_meridian : meridian;

    if (i == polygons_total) {
      // if full sweep this works
      meridian = prime_meridian;
      if (!opts.double_sweep)
        reverse_poly_indexes_on_y(geom, meridian, opts.eps);
    }
    else {
      // add one 'meridian' of points
      meridian.clear();
      for (int j = 0; j < opts.ncon_order; j++) {
        int m = prime_meridian[j];
        if (m == pole[0]->idx || m == pole[1]->idx)
          meridian.push_back(m);
        else {
          // Rotate Point Counter-Clockwise about Y Axis (looking down through Y
          // Axis)
          Vec3d v = Trans3d::rotate(0, deg2rad(arc * i), 0) * verts[m];
          add_coord(geom, coordinates, v);
          meridian.push_back(verts.size() - 1);
        }
      }
    }

    // calculate longitude in face list for coloring faces later. latitude
    // calculated after model creation
    int lon_back = (i - 1) % (opts.longitudes.front() / 2);
    int lon_front = lon_back + (opts.longitudes.front() / 2);

    vector<int> top_edge;
    vector<int> bottom_edge;
    vector<int> first_edge;

    int polygon_no = 0;

    for (int p = 0; p < num_polygons; p++) {
      int m_start = base_polygon * p;
      for (int j = 0; j <= base_polygon; j++) {
        // prime bottom edge
        if (j == 0) {
          int m = meridian[j + m_start];
          int ml = meridian_last[j + m_start];
          bottom_edge.clear();
          bottom_edge.push_back(m);
          if (m != ml) // not a pole
            bottom_edge.push_back(ml);
          first_edge = bottom_edge;
          continue;
        }

        // everything else is the top edge
        if (j == base_polygon)
          top_edge = first_edge;
        else {
          int m = meridian[j + m_start];
          int ml = meridian_last[j + m_start];
          top_edge.clear();
          top_edge.push_back(m);
          if (m != ml) // not a pole
            top_edge.push_back(ml);
        }

        vector<int> face;
        for (int &k : bottom_edge)
          face.insert(face.end(), k);
        // top edge is reversed);
        for (int k = top_edge.size() - 1; k >= 0; k--)
          face.insert(face.end(), top_edge[k]);

        vector<vector<int>> face_parts =
            split_bow_ties(geom, coordinates, face, opts.eps);

        vector<int> split_face_idx;

        // if bow tie is split there will be 2 faces (2 triangles)
        for (auto &face_part : face_parts) {
          double face_z = centroid(verts, face_part)[2];

          int lon = (face_z > 0.0) ? lon_front : lon_back;
          add_face(geom, face_list, face_part, -1, lon, polygon_no);

          // keep track of the split face indexes
          split_face_idx.push_back(face_list.size() - 1);

          // the front face is what is rotated
          if (face_z > 0.0)
            face_list.back()->rotate = true;

          // always add edges to discover face latitudes
          if (top_edge.size() == 2)
            add_edge_wrapper(geom, edge_list, top_edge, -1, lon_front,
                             lon_back);
          if (bottom_edge.size() == 2)
            add_edge_wrapper(geom, edge_list, bottom_edge, -1, lon_front,
                             lon_back);
        }

        if (split_face_idx.size()) {
          original_faces.push_back(face);
          split_face_indexes.push_back(split_face_idx);
        }
        split_face_idx.clear();

        bottom_edge = top_edge;
      }

      polygon_no++;
    }
  }

  return;
}

// for method 2: to hide uneeded edges
// note that function used to reverse indented based on manual inner and outer
// radii
void mark_indented_edges_invisible(const vector<edgeList *> &edge_list,
                                   const vector<poleList *> &pole,
                                   const bool radius_reverse,
                                   const ncon_opts &opts)
{
  // for method 2 we used n/2
  int n = opts.ncon_order / 2;

  for (auto i : edge_list) {
    int lat = i->lat;
    bool set_invisible = (is_even(n) && opts.point_cut && !is_even(lat)) ||
                         (is_even(n) && !opts.point_cut && is_even(lat)) ||
                         (!is_even(n) && is_even(lat));
    if (radius_reverse)
      set_invisible = (set_invisible) ? false : true;
    // negative latitudes to temporarily label indented edges
    if (set_invisible) {
      int lat = i->lat;
      lat = -lat - 2;
      i->lat = lat;
    }
  }

  for (unsigned int i = 0; i < 2; i++) {
    if (pole[i]->idx > -1) {
      int lat = pole[i]->lat;
      bool set_invisible = (is_even(n) && opts.point_cut && !is_even(lat)) ||
                           (is_even(n) && !opts.point_cut && is_even(lat)) ||
                           (!is_even(n) && is_even(lat));
      if (radius_reverse)
        set_invisible = (set_invisible) ? false : true;
      if (set_invisible)
        pole[i]->lat = -1;
    }
  }
}

void restore_indented_edges(const vector<edgeList *> &edge_list,
                            const ncon_opts &opts)
{
  for (auto i : edge_list) {
    int lat = i->lat;
    if (lat < -1) {
      lat = (opts.hide_indent) ? -1 : std::abs(lat + 2);
      i->lat = lat;
    }
  }
}

// for method 2
// note: it is called with n and not 2n
// d is a copy
vector<pair<int, int>> get_lat_pairs(const int n, int d, const bool point_cut)
{
  vector<pair<int, int>> lat_pairs;
  pair<int, int> lats;

  // d < n/2
  if (d > n / 2)
    d = n - d;

  int max_lats = n - 1;

  bool toggle = (!is_even(n) || !point_cut) ? true : false;
  for (int i = 0; i < n; i++) {
    int next = (toggle) ? i - (2 * d - 1) : i + (2 * d - 1);
    toggle = (toggle) ? false : true;
    if (next > max_lats)
      next = n - (next % max_lats);
    else if (next < 0)
      next = std::abs(next) - 1;

    lats.first = i;
    lats.second = next;
    if (lats.first > lats.second)
      swap(lats.first, lats.second);

    lat_pairs.push_back(lats);
  }

  sort(lat_pairs.begin(), lat_pairs.end());
  auto li = unique(lat_pairs.begin(), lat_pairs.end());
  lat_pairs.erase(li, lat_pairs.end());

  return lat_pairs;
}

// for method 2: set latitude numbers for pairs of faces of shell models
// then split_face_indexes is filled
void find_split_faces_shell_model(const Geometry &geom,
                                  const vector<faceList *> &face_list,
                                  const vector<edgeList *> &edge_list,
                                  const vector<poleList *> &pole,
                                  vector<vector<int>> &split_face_indexes,
                                  const ncon_opts &opts)
{
  // for method 2 we used n/2
  int n = opts.ncon_order / 2;

  const vector<vector<int>> &faces = geom.faces();
  const vector<Vec3d> &verts = geom.verts();

  int lat = 0;
  if (opts.d != 1) {
    // if radii were inverted, point_cut will have been reversed but not what
    // function is expecting
    bool pc = opts.point_cut;
    if (opts.radius_inversion)
      pc = !pc;
    vector<pair<int, int>> lat_pairs = get_lat_pairs(n, opts.d, pc);

    for (auto &lat_pair : lat_pairs) {
      for (auto j : face_list) {
        if (j->lat == lat_pair.first || j->lat == lat_pair.second)
          j->lat = lat;
      }
      lat++;
    }

    if (opts.hide_indent) {
      for (auto i : edge_list) {
        lat = i->lat;
        if (lat > 1) {
          int adjust = (is_even(lat)) ? 0 : 1;
          i->lat = (int)floor((double)lat / 2) + adjust;
        }
      }

      // fix south pole
      if (pole[1]->idx != -1) {
        lat = pole[1]->lat;
        if (lat > 1) {
          int adjust = (is_even(lat)) ? 0 : 1;
          pole[1]->lat = (int)floor((double)n / 2) + adjust;
        }
      }
    }
  }

  // collect split faces
  lat = 0;
  vector<pair<double, int>> face_zs;
  do {
    face_zs.clear();
    for (unsigned int i = 0; i < face_list.size(); i++) {
      // collect Z centroids of faces
      if (face_list[i]->lat == lat) {
        Vec3d face_cent = centroid(verts, faces[i]);
        double angle = angle_around_axis(face_cent, Vec3d(1, 0, 0), Vec3d::Y);
        face_zs.push_back(make_pair(angle, i));
      }
    }

    sort(face_zs.begin(), face_zs.end());
    reverse(face_zs.begin(), face_zs.end());

    vector<int> split_face_idx;
    double last_z = std::numeric_limits<double>::max();
    for (unsigned int i = 0; i < face_zs.size(); i++) {
      if (double_ne(face_zs[i].first, last_z, opts.eps)) {
        if (i > 0) {
          split_face_indexes.push_back(split_face_idx);
          split_face_idx.clear();
        }
        last_z = face_zs[i].first;
      }
      split_face_idx.push_back(face_zs[i].second);
    }
    if (split_face_idx.size())
      split_face_indexes.push_back(split_face_idx);
    lat++;
  } while (face_zs.size());
}

// if return_calc is true, it returns the calculated values of radius even if
// they are set before hand
// inner_radius, outer_radius, arc, d are changed
void calc_radii(double &inner_radius, double &outer_radius, double &arc,
                const int N, int &d, const ncon_opts &opts,
                const bool return_calc)
{
  // if shell model is created opts.ncon_order needs to be divided by 2 to get
  // the correct outer radius, except when d is 1
  arc = 360.0 / (N / ((opts.build_method == 2 && d != 1) ? 2 : 1)) * d;
  // but this causes a problem for 2N/N so make an exception (always 360/4 = 90)
  if (opts.build_method == 2 && double_eq(arc, 180.0, opts.eps))
    arc = 360.0 / N * d;

  double interior_angle = (180.0 - arc) / 2.0;
  double inner_radius_calc = 0;
  double outer_radius_calc = sin(deg2rad(interior_angle)) / sin(deg2rad(arc));
  if (std::isnan(outer_radius))
    outer_radius = outer_radius_calc;

  if (opts.build_method == 2) {
    // reculate arc
    arc = 360.0 / N;
    // interior angle (not used)
    // interior_angle = (180.0 - arc) / 2.0;

    int n_calc = N / ((d != 1) ? 2 : 1);
    if (2 * d > n_calc)
      d = n_calc - d;
    // formula furnished by Adrian Rossiter
    // r = R * cos(pi*m/n) / cos(pi*(m-1)/n)
    inner_radius_calc = (outer_radius_calc * cos(M_PI * d / n_calc) /
                         cos(M_PI * (d - 1) / n_calc));
    if (std::isnan(inner_radius))
      inner_radius = (d == 1) ? outer_radius_calc : inner_radius_calc;
  }

  // patch: radii cannot be exactly 0
  if (double_eq(inner_radius, 0.0, opts.eps))
    inner_radius += 1e-4;
  if (double_eq(outer_radius, 0.0, opts.eps))
    outer_radius += 1e-4;

  if (return_calc) {
    outer_radius = outer_radius_calc;
    inner_radius = (d == 1) ? outer_radius_calc : inner_radius_calc;
  }
}

// for methods 1 and 2: first longitude of vertices to form globe
// inner_radius, outer_radius set
// point_cut_calc changed from side cut to point cut if build method 2 and d > 1
void build_prime_meridian(Geometry &geom, vector<int> &prime_meridian,
                          vector<coordList *> &coordinates,
                          double &inner_radius, double &outer_radius,
                          bool &point_cut_calc, const ncon_opts &opts)
{
  int n = opts.ncon_order; // pass n and d as const
  int d = opts.d;
  double arc = 0;

  // for build method 2, n will be doubled
  calc_radii(inner_radius, outer_radius, arc, n, d, opts, false);

  bool radii_swapped = false;
  double angle = -90.0;
  // side cut
  if (is_even(opts.ncon_order) && !point_cut_calc) {
    if (opts.build_method == 2 && d != 1) {
      swap(outer_radius, inner_radius);
      radii_swapped = true;
      // now treat it like a point cut
      point_cut_calc = true;
    }
    else {
      // angle += 180.0/opts.ncon_order;
      angle += (arc / 2.0);
    }
  }

  int num_vertices = longitudinal_faces(opts.ncon_order, point_cut_calc) + 1;
  for (int i = 0; i < num_vertices; i++) {
    prime_meridian.push_back(i);
    double radius =
        ((opts.build_method != 2) || is_even(i)) ? outer_radius : inner_radius;
    add_coord(
        geom, coordinates,
        Vec3d(cos(deg2rad(angle)) * radius, sin(deg2rad(angle)) * radius, 0));
    angle += arc;
  }

  // swap radii back for future reference
  if (radii_swapped)
    swap(outer_radius, inner_radius);
}

// methods 1 and 2: calculate polygon numbers
vector<int> calc_polygon_numbers(int n, int d, bool point_cut)
{
  int num_polygons = gcd(n, d);

  vector<int> polygon_numbers(n);
  int polygon_no = 0;
  for (int i = 0; i < n; i++) {
    polygon_numbers[i] = polygon_no;
    polygon_no++;
    if (polygon_no >= num_polygons)
      polygon_no = 0;
  }

  int lats = num_lats(n, point_cut);
  if (!point_cut)
    lats -= 2;

  int start = (point_cut) ? 1 : 0;
  unsigned int end = polygon_numbers.size() - 1;

  for (int i = start; i <= lats; i++) {
    if (polygon_numbers[i] > polygon_numbers[end])
      polygon_numbers[i] = polygon_numbers[end];
    end--;
  }

  polygon_numbers.resize(lats + (is_even(n) ? 1 : 0));
  // will be upside down versus lat in form_globe
  reverse(polygon_numbers.begin(), polygon_numbers.end());

  return (polygon_numbers);
}

// methods 1 and 2
void form_globe(Geometry &geom, const vector<int> &prime_meridian,
                vector<coordList *> &coordinates, vector<faceList *> &face_list,
                vector<edgeList *> &edge_list, const bool point_cut_calc,
                const bool second_half, const ncon_opts &opts)
{
  const vector<vector<int>> &faces = geom.faces();
  const vector<Vec3d> &verts = geom.verts();

  int half = opts.longitudes.front() / 2;

  // patch for polygon number in shell models
  vector<int> polygon_numbers;
  if (opts.build_method == 2)
    polygon_numbers = calc_polygon_numbers(
        opts.ncon_order / (opts.d == 1 ? 1 : 2), opts.d, opts.point_cut);

  // note: the algorithm can make "opts.longitudes.back()" faces but for method
  // 2 needs the whole model
  int longitudes_back = (opts.build_method == 1 || opts.hybrid)
                            ? opts.longitudes.back()
                            : opts.longitudes.front();
  int half_model_marker = 0;

  double arc = 360.0 / opts.longitudes.front();

  int lon_faces = longitudinal_faces(opts.ncon_order, point_cut_calc);
  // used to coordinate latitudes with bands
  int inc1 = (is_even(opts.ncon_order) && point_cut_calc) ? 0 : 1;

  // We close an open model longitudinally only if is desired, not a full model,
  // or only 1 zone is shown
  // only works for build method 1, build method 2 generates the whole model
  int longitudinal_cover = 0;
  if (strchr(opts.closure.c_str(), 'h') && !full_model(opts.longitudes) &&
      longitudes_back != 1)
    longitudinal_cover++;

  for (int i = 1; i <= longitudes_back + longitudinal_cover; i++) {
    for (unsigned int j = 0; j < prime_meridian.size(); j++) {
      if (((is_even(opts.ncon_order) && point_cut_calc) && (j != 0) &&
           ((int)j < prime_meridian.back())) ||
          (is_even(opts.ncon_order) && !point_cut_calc) ||
          (!is_even(opts.ncon_order) && (j > 0))) {
        if ((i != opts.longitudes.front()) && (i != longitudes_back + 1)) {
          // Rotate Point Counter-Clockwise about Y Axis (looking down through Y
          // Axis)
          Vec3d v = Trans3d::rotate(0, deg2rad(arc * i), 0) *
                    verts[prime_meridian[j]];
          add_coord(geom, coordinates, v);
        }
        if ((i == ((opts.longitudes.front() / 2) + 1)) &&
            (half_model_marker == 0))
          half_model_marker = verts.size() - 2;
      }

      bool use_prime;
      int k = 0;
      int l = 0;

      // store in face list for coloring faces later
      int lat = lon_faces - j + inc1;
      int lon = i - 1;

      if ((i == 1) || (i == opts.longitudes.front()) ||
          (i == longitudes_back + 1)) {
        use_prime = true;
        l = 0;

        // Reset coordinate count for last row of faces
        if ((i == opts.longitudes.front()) || (i == longitudes_back + 1)) {
          k = lon_faces - j;
          if (is_even(opts.ncon_order) && point_cut_calc)
            k--;
          // no color for opts.longitudes that span more than one meridian
          if (opts.longitudes.front() - longitudes_back > 1) {
            lat = -1;
            lon = -1;
          }
        }
      }
      else {
        use_prime = false;
        l = faces.size() - lon_faces;
      }

      // for second half of hybrid, adjust longitudes for later mirror image
      lon = (second_half ? std::abs(lon - half) + half - 1 : lon);

      // when j = 0 "prime the system"
      if (j > 0) {

        // patch for polygon number in shell models
        int polygon_no = 0;
        if (opts.build_method == 2) {
          double a = (double)lat / 2.0;
          int idx =
              (is_even(opts.ncon_order / 2)) ? (int)ceil(a) : (int)floor(a);
          if (!opts.point_cut && !is_even(j))
            idx--;
          polygon_no = polygon_numbers[idx];
        }

        if (((is_even(opts.ncon_order) && point_cut_calc) ||
             !is_even(opts.ncon_order)) &&
            (j == 1)) {
          // fprintf(stderr,"South Triangle\n");
          vector<int> face;
          face.push_back(verts.size() - 1 - k);
          face.push_back(prime_meridian.front());

          if (use_prime)
            face.push_back(prime_meridian[1]);
          else
            face.push_back(faces[l][0]);
          add_face(geom, face_list, face, lat, lon, polygon_no);

          // always add edge and can be deleted later. logic for method 1
          // unchanged
          if (opts.edge_coloring_method || opts.build_method > 1)
            add_edge(geom, edge_list, make_edge(face[0], face[2]), lat, lon);
        }
        else if ((is_even(opts.ncon_order) && point_cut_calc) &&
                 (j == prime_meridian.size() - 1)) {
          // fprintf(stderr,"North Triangle\n");
          vector<int> face;
          face.push_back(verts.size() - 1);
          face.push_back(prime_meridian.back());

          if (use_prime)
            face.push_back(prime_meridian.back() - 1);
          else
            face.push_back(faces[l][0]);
          add_face(geom, face_list, face, lat, lon, polygon_no);
        }
        else {
          // fprintf(stderr,"Square\n");
          vector<int> face;
          face.push_back(verts.size() - 1 - k);
          face.push_back(verts.size() - 2 - k);

          if (use_prime) {
            face.push_back(prime_meridian[prime_meridian[j - 1]]);
            face.push_back(prime_meridian[prime_meridian[j]]);
          }
          else {
            face.push_back(faces[l][1]);
            face.push_back(faces[l][0]);
          }
          add_face(geom, face_list, face, lat, lon, polygon_no);

          // always add edge and can be deleted later. logic for method 1
          // unchanged
          if (opts.edge_coloring_method || opts.build_method > 1) {
            add_edge(geom, edge_list, make_edge(face[0], face[3]), lat, lon);

            // patch with the extra edge below, one rotate flag gets missed
            if (half_model_marker != 0)
              edge_list.back()->rotate = true;

            // the last square has two edges if there is no bottom triangle
            if (lat == lon_faces &&
                (is_even(opts.ncon_order) && !point_cut_calc))
              add_edge(geom, edge_list, make_edge(face[1], face[2]), lat + 1,
                       lon);
          }
        }

        if (half_model_marker != 0) {
          face_list.back()->rotate = true;
          if (opts.edge_coloring_method || opts.build_method > 1)
            edge_list.back()->rotate = true;
        }
      }
    }
  }
}

// add caps to method 1 and 2 models
// caps indexes are retained
void add_caps(Geometry &geom, vector<coordList *> &coordinates,
              vector<faceList *> &face_list, const vector<poleList *> &pole,
              vector<int> &caps, const bool point_cut_calc,
              const ncon_opts &opts)
{
  const vector<vector<int>> &faces = geom.faces();
  const vector<Vec3d> &verts = geom.verts();

  // note: faces but for method 2 needs the whole model
  vector<int> lons(2);
  lons.front() = opts.longitudes.front();
  lons.back() = (opts.build_method == 2 && !opts.hybrid)
                    ? opts.longitudes.front()
                    : opts.longitudes.back();

  int lats = num_lats(opts.ncon_order, point_cut_calc);
  // if twist is 0, there is only one cap so lower longitude value by 1 (method
  // 2)
  int lon_front =
      (opts.build_method == 2)
          ? lons.front() / 2 -
                ((opts.mod_twist == 0 || (opts.hybrid && !opts.point_cut)) ? 1
                                                                           : 0)
          : -1;
  int lon_back = (opts.build_method == 2) ? lon_front - 1 : -1;

  pole[0]->idx = -1;
  pole[0]->lat = -1;
  pole[1]->idx = -1;
  pole[1]->lat = lats;

  // Even order and point cut always have poles (and south pole = 0)
  if (is_even(opts.ncon_order) && point_cut_calc) {
    pole[0]->idx = longitudinal_faces(opts.ncon_order, point_cut_calc);
    pole[1]->idx = 0;
  }
  else if (!is_even(opts.ncon_order) || opts.add_poles) {
    if (opts.add_poles)
      pole[0]->idx = 0;
    if (!half_model(lons) || opts.add_poles)
      pole[1]->idx = 0;
  }

  if (pole[0]->idx > -1)
    pole[0]->lat = 0;

  if (opts.add_poles) {
    if (!strchr(opts.hide_elems.c_str(), 't') &&
        ((is_even(opts.ncon_order) && !point_cut_calc) ||
         !is_even(opts.ncon_order))) {
      Vec3d v;
      v[0] = 0;
      v[1] = verts.back()[1];
      v[2] = 0;
      add_coord(geom, coordinates, v);
      pole[0]->idx = verts.size() - 1;

      if (!strchr(opts.hide_elems.c_str(), 'b') &&
          ((is_even(opts.ncon_order) && !point_cut_calc))) {
        Vec3d v;
        v[0] = 0;
        v[1] = verts.front()[1];
        v[2] = 0;
        add_coord(geom, coordinates, v);
        pole[1]->idx = verts.size() - 1;
      }
    }
  }

  // Top
  if (!strchr(opts.hide_elems.c_str(), 't') &&
      ((is_even(opts.ncon_order) && !point_cut_calc) ||
       !is_even(opts.ncon_order))) {
    int j = 1;
    bool split_done = false;
    vector<int> face;
    for (int i = longitudinal_faces(opts.ncon_order, point_cut_calc) - 1;
         j <= lons.back();
         i += longitudinal_faces(opts.ncon_order, point_cut_calc)) {
      if (j == lons.back()) {
        if (!full_model(lons))
          face.push_back(faces[i].back());
        if (opts.split && !split_done && opts.add_poles &&
            (lons.back() > lons.front() / 2))
          face.push_back(pole[0]->idx);
        else if (!opts.split || (opts.split && split_done) ||
                 half_model(lons) ||
                 (!half_model(lons) && (lons.back() != (lons.front() / 2) + 1)))
          face.push_back(faces[i].front());
      }
      else
        face.push_back(faces[i].back());

      if (opts.split && !split_done && (j == (lons.front() / 2) + 1)) {
        add_face(geom, face_list, face, 0, lon_back);
        caps.push_back((int)face_list.size() - 1);
        face.clear();
        i -= longitudinal_faces(opts.ncon_order, point_cut_calc);
        j--;
        split_done = true;
      }
      j++;
    }

    if (split_done && !opts.add_poles)
      face.push_back((faces.back()).front());

    if (opts.add_poles)
      face.push_back(pole[0]->idx);

    add_face(geom, face_list, face, 0, lon_front);
    caps.push_back((int)face_list.size() - 1);
    if (!opts.hybrid)
      face_list.back()->rotate = true;
  }

  // Bottom
  if (!strchr(opts.hide_elems.c_str(), 'b') &&
      ((is_even(opts.ncon_order) && !point_cut_calc))) {
    int j = 1;
    bool split_done = false;
    vector<int> face;

    for (int i = 0; j <= lons.back();
         i += longitudinal_faces(opts.ncon_order, point_cut_calc)) {
      if (j == lons.back()) {
        if (!full_model(lons))
          face.push_back(faces[i][2]);
        if (opts.split && !split_done && opts.add_poles &&
            (lons.back() > lons.front() / 2))
          face.push_back(pole[1]->idx);
        else if (!opts.split || (opts.split && split_done) ||
                 half_model(lons) ||
                 (!half_model(lons) && (lons.back() != (lons.front() / 2) + 1)))
          face.push_back(faces[i][1]);
      }
      else
        face.push_back(faces[i][2]);

      if (opts.split && !split_done && (j == (lons.front() / 2) + 1)) {
        add_face(geom, face_list, face, lats - 1, lon_back);
        caps.push_back((int)face_list.size() - 1);
        face.clear();
        i -= longitudinal_faces(opts.ncon_order, point_cut_calc);
        j--;
        split_done = true;
      }

      j++;
    }

    if (split_done && !opts.add_poles)
      face.push_back((faces.back()).front());

    if (opts.add_poles)
      face.push_back(pole[1]->idx);

    add_face(geom, face_list, face, lats - 1, lon_front);
    caps.push_back((int)face_list.size() - 1);
    if (!opts.hybrid)
      face_list.back()->rotate = true;
  }
}

// for method 1 covering
void close_latitudinal(Geometry &geom, vector<faceList *> &face_list,
                       const vector<poleList *> &pole, const ncon_opts &opts)
{
  bool point_cut_calc = opts.point_cut;
  if (opts.build_method == 2 && opts.d > 1)
    point_cut_calc = true;

  vector<vector<int>> face(3);

  // Cover one side
  if (opts.add_poles && (is_even(opts.ncon_order) && !point_cut_calc))
    face[0].push_back(pole[1]->idx);

  int j = 0;
  for (int i = 1; i <= longitudinal_faces(opts.ncon_order, point_cut_calc) + 1;
       i++) {
    face[0].push_back(j);
    j++;
  }

  // If half model it is all one opts.closure, continue face[0] as face[1]
  if (half_model(opts.longitudes))
    face[1] = face[0];
  else {
    if (opts.add_poles && ((is_even(opts.ncon_order) && !point_cut_calc) ||
                           !is_even(opts.ncon_order)))
      face[0].push_back(pole[0]->idx);

    if (strchr(opts.closure.c_str(), 'v') && (face[0].size() > 2)) {
      add_face(geom, face_list, face[0], -2, -2);
      face_list.back()->rotate = true;
    }
  }

  // Cover other side
  j = geom.verts().size() - 1;

  if (opts.add_poles) {
    if (!is_even(opts.ncon_order))
      j--;
    else if (is_even(opts.ncon_order) && !point_cut_calc)
      j -= 2;
  }

  int k = longitudinal_faces(opts.ncon_order, point_cut_calc) + 1;
  if (opts.add_poles || (is_even(opts.ncon_order) && point_cut_calc)) {
    if (!half_model(opts.longitudes))
      face[1].push_back(pole[0]->idx);
    if (is_even(opts.ncon_order) && point_cut_calc)
      k--;
  }

  for (int i = 1; i < k; i++)
    face[1].push_back(j--);

  if (is_even(opts.ncon_order) && !point_cut_calc)
    face[1].push_back(j);

  if (opts.add_poles || (((is_even(opts.ncon_order) && point_cut_calc) ||
                          !is_even(opts.ncon_order)) &&
                         !half_model(opts.longitudes)))
    face[1].push_back(pole[1]->idx);

  if (strchr(opts.closure.c_str(), 'v') && (face[1].size() > 2)) {
    add_face(geom, face_list, face[1], -2, -2);
    face_list.back()->rotate = true;
  }

  // In cases of odd opts.longitudes sides or even opts.longitudes side cuts
  // which are not half models this generates a third side to cover
  if (!opts.add_poles &&
      ((is_even(opts.ncon_order) && !point_cut_calc) ||
       !is_even(opts.ncon_order)) &&
      (!half_model(opts.longitudes)) && (strchr(opts.closure.c_str(), 'v'))) {
    face[2].push_back(0);
    if (is_even(opts.ncon_order) && !point_cut_calc)
      face[2].push_back(face[1].back());
    face[2].push_back(face[1].front());
    face[2].push_back(face[0].back());
    add_face(geom, face_list, face[2], -2, -2);
    face_list.back()->rotate = true;
  }
}

bool cmp_angle(const pair<pair<double, double>, int> &a,
               const pair<pair<double, double>, int> &b, const double eps)
{
  pair<double, double> ar_a = a.first;
  pair<double, double> ar_b = b.first;
  bool ret = double_eq(ar_a.first, ar_b.first, eps);
  if (ret)
    ret = double_lt(ar_a.second, ar_b.second, eps);
  else
    ret = double_lt(ar_a.first, ar_b.first, eps);
  return ret;
}

class angle_cmp {
public:
  double eps;
  angle_cmp(double ep) : eps(ep) {}
  bool operator()(const pair<pair<double, double>, int> &a,
                  const pair<pair<double, double>, int> &b) const
  {
    return cmp_angle(a, b, eps);
  }
};

// untangle polar orbit
void sort_polar_orbit(Geometry &geom, vector<polarOrb *> &polar_orbit,
                      const double eps)
{
  const vector<Vec3d> &verts = geom.verts();

  Vec3d v0 = verts[polar_orbit[0]->coord_no];
  unsigned int sz = polar_orbit.size();
  vector<pair<pair<double, double>, int>> angles(sz);
  for (unsigned int i = 0; i < sz; i++) {
    int j = polar_orbit[i]->coord_no;
    angles[i].second = j;
    pair<double, double> angle_and_radius;
    angle_and_radius.first = angle_in_range(
        rad2deg(angle_around_axis(v0, verts[j], Vec3d(0, 0, 1))), eps);
    angle_and_radius.second = verts[j].len();
    angles[i].first = angle_and_radius;
  }

  // sort on angles
  sort(angles.begin(), angles.end(), angle_cmp(eps));

  for (unsigned int i = 0; i < sz; i++)
    polar_orbit[i]->coord_no = angles[i].second;
}

void find_polar_orbit(Geometry &geom, vector<polarOrb *> &polar_orbit,
                      const int build_method, const double eps)
{
  vector<Vec3d> &verts = geom.raw_verts();
  for (unsigned int i = 0; i < verts.size(); i++) {
    if (double_eq(verts[i][2], 0.0, eps)) {
      polar_orbit.push_back(new polarOrb(i));
    }
    // in case of build method 3, some points were set aside so not to get
    // into twist plane
    else if ((build_method == 3) && fabs(verts[i][2]) < eps * 3.0) {
      verts[i][2] = 0.0;
    }
  }

  // update for n/m models. works for all, so do it this way now
  sort_polar_orbit(geom, polar_orbit, eps);
}

void ncon_twist(Geometry &geom, const vector<polarOrb *> &polar_orbit,
                const vector<coordList *> &coordinates,
                const vector<faceList *> &face_list,
                const vector<edgeList *> &edge_list, const int ncon_order,
                const int twist, const vector<int> &longitudes)
{
  // this function wasn't designed for twist 0
  if (twist == 0)
    return;

  vector<vector<int>> &faces = geom.raw_faces();
  vector<vector<int>> &edges = geom.raw_edges();
  vector<Vec3d> &verts = geom.raw_verts();

  double arc = (360.0 / ncon_order) * twist * (+1);
  Trans3d rot = Trans3d::rotate(0, 0, deg2rad(-arc));

  for (unsigned int i = 0; i < face_list.size(); i++) {
    if (face_list[i]->rotate) {
      for (int coord_no : faces[i]) {
        if (!coordinates[coord_no]->rotated) {
          verts[coord_no] = rot * verts[coord_no];
          coordinates[coord_no]->rotated = true;
        }
      }
    }
  }

  // Patch - If an open model, some of the polar circle coordinates are not yet
  // rotated
  if (!full_model(longitudes))
    for (auto i : polar_orbit) {
      int j = i->coord_no;
      if (!coordinates[j]->rotated) {
        verts[j] = rot * verts[j];
        coordinates[j]->rotated = true;
      }
    }

  // Create Doubly Circularly Linked List
  for (unsigned int i = 0; i < polar_orbit.size(); i++)
    polar_orbit[i]->forward = i + 1;
  polar_orbit.back()->forward = 0;
  for (unsigned int i = polar_orbit.size() - 1; i > 0; i--)
    polar_orbit[i]->backward = i - 1;
  polar_orbit.front()->backward = polar_orbit.size() - 1;

  for (unsigned int k = 0; k < polar_orbit.size(); k++) {
    int p, q;
    p = polar_orbit[k]->coord_no;
    if (twist > 0)
      q = polar_orbit[k]->forward;
    else
      q = polar_orbit[k]->backward;
    for (int n = 1; n < std::abs(twist); n++) {
      if (twist > 0)
        q = polar_orbit[q]->forward;
      else
        q = polar_orbit[q]->backward;
    }
    q = polar_orbit[q]->coord_no;

    for (unsigned int i = 0; i < faces.size(); i++)
      if (!face_list[i]->rotate)
        for (int &j : faces[i])
          if (j == p)
            j = q * (-1);

    for (unsigned int i = 0; i < edges.size(); i++)
      if (!edge_list[i]->rotate)
        for (int &j : edges[i])
          if (j == p)
            j = q * (-1);
  }

  for (auto &face : faces)
    for (int &j : face)
      if (j < 0)
        j = std::abs(j);

  for (auto &edge : edges)
    for (int &j : edge)
      if (j < 0)
        j = std::abs(j);

  // edges may have been made to be out of numerical order
  for (auto &edge : edges)
    if (edge[0] > edge[1])
      swap(edge[0], edge[1]);
}

// surfaces, case2, case1_twist are changed
void find_surface_count(const vector<surfaceTable *> &surface_table,
                        const int twist, int &surfaces, bool &case2,
                        int &case1_twist)
{
  surfaces = 0;
  case2 = false;
  case1_twist = 0;
  for (auto i : surface_table) {
    if (twist == i->twist) {
      surfaces = i->surfaces;
      case2 = i->case2;
      case1_twist = i->case1_twist;
      break;
    }
  }
}

// ncon_order, point_cut, hybrid are not from opts
void build_surface_table(vector<surfaceTable *> &surface_table,
                         const int max_twist, const int ncon_order,
                         const bool point_cut, const bool hybrid)
{
  // coding idea furnished by Adrian Rossiter
  int axis_edges = 0;
  int n = ncon_order;
  int t_mod = 0;

  if (!is_even(ncon_order) || hybrid) {
    axis_edges = 1;
    if (hybrid) {
      n *= 2;
      t_mod = 1;
    }
  }
  else if (is_even(ncon_order) && !point_cut)
    axis_edges = 2;

  for (int twist = 2; twist <= max_twist; twist++) {
    int total_surfaces = (int)((gcd(n, 2 * twist - t_mod) + axis_edges) / 2);
    // fprintf(stderr,"twist = %d surfaces = %d\n",twist,total_surfaces);

    // subtract out discontinuous surfaces for table
    int continuous_surfaces = total_surfaces - axis_edges;
    bool factor = false;
    if (hybrid)
      factor = (n % (2 * twist - 1) == 0) ? true : false;
    else
      factor = (n % twist == 0) ? true : false;
    int case1_twist = total_surfaces;
    if (is_even(ncon_order) && !point_cut)
      case1_twist--;

    // only store those with more than minimum surface counts
    if ((continuous_surfaces > 1) ||
        (continuous_surfaces > 0 &&
         (!(is_even(ncon_order) && point_cut) || hybrid))) {
      surface_table.push_back(new surfaceTable(twist, continuous_surfaces,
                                               (factor ? twist : case1_twist),
                                               (factor ? false : true)));
    }
  }
}

void model_info(const Geometry &geom, const ncon_opts &opts)
{
  int num_parts = geom.get_info().num_parts();

  if (opts.build_method == 2)
    fprintf(stderr, "shell models are physically non-compounds (but can be "
                    "painted as compounds)\n");
  else if (opts.build_method == 3) {
    if (num_parts > 1)
      fprintf(stderr, "%d physical compound parts were found\n", num_parts);
    else
      fprintf(stderr, "The model is not a compound\n");
  }
  fprintf(stderr, "per actual connectivity, some n_icons may have conflicting "
                  "compound counts\n");
  fprintf(stderr, "-f c and -f a are alternate ways to view compound parts\n");

  if (opts.angle != 0) {
    int actual_order = opts.ncon_order * 2;
    int actual_d = opts.d * 2;
    int actual_twist = opts.twist * 2;
    if (opts.hybrid) {
      actual_twist -= 1;
      fprintf(stderr,
              "With angle, this Hybrid is similar to N = %d/%d  twist "
              "= %d  side cut\n",
              actual_order, actual_d, actual_twist);
    }
    else {
      fprintf(stderr,
              "with angle, this model is similar to N = %d/%d  twist = "
              "%d  side cut\n",
              actual_order, actual_d, actual_twist);
    }
  }

  if (opts.build_method == 2)
    fprintf(stderr, "shell model radii: outer = %.17lf  inner = %.17lf\n",
            opts.outer_radius, opts.inner_radius);

  unsigned long fsz = geom.faces().size();
  unsigned long vsz = geom.verts().size();
  fprintf(stderr,
          "The graphical model shown has %lu faces, %lu vertices, and "
          "%lu implicit edges\n",
          fsz, vsz, (fsz + vsz - 2));

  fprintf(stderr, "========================================\n");
}

// surface_table, sd will be changed
// ncon_order, point_cut, twist, hybrid, info are not from opts
void ncon_info(const int ncon_order, const bool point_cut, const int twist,
               const bool hybrid, const bool info,
               vector<surfaceTable *> &surface_table, surfaceData &sd)
{
  int first, last, forms, chiral, nonchiral, unique;

  if (info) {
    fprintf(stderr, "\n");
    fprintf(stderr, "This is %s n-icon of order %d twisted %d time%s\n",
            (hybrid ? "a hybrid" : (point_cut ? "point cut" : "side cut")),
            ncon_order, twist, ((twist == 1) ? "" : "s"));

    fprintf(stderr, "\n");
    if (hybrid) {
      fprintf(stderr,
              "This is an order %d, hybrid n-icon. By hybrid, it means "
              "it is one half a\n",
              ncon_order);
      fprintf(stderr, "point cut n-icon joined to half of a side cut n-icon. "
                      "Hybrid n-icons are\n");
      fprintf(stderr, "of even order and have at least one discontinuous "
                      "surface and one\n");
      fprintf(stderr, "discontinuous edge. As a non-faceted smooth model it "
                      "would be self dual.\n");
      fprintf(stderr, "Circuit patterns are self dual.\n");
    }
    else if (is_even(ncon_order) && point_cut) {
      fprintf(stderr,
              "The order of this n-icon, %d, is even. It is, what is "
              "termed, a Point Cut.\n",
              ncon_order);
      fprintf(stderr, "Even order, point cut n-icons have at least one "
                      "continuous surface and\n");
      fprintf(stderr, "two discontinuous edges. As a non-faceted smooth model "
                      "it would be dual\n");
      fprintf(stderr, "to the Side Cut n-icon. Circuit patterns are dual to "
                      "Side Cut n-icons.\n");
    }
    else if (is_even(ncon_order) && !point_cut) {
      fprintf(stderr,
              "The order of this n-icon, %d, is even. It is, what is "
              "termed, a Side Cut.\n",
              ncon_order);
      fprintf(stderr, "Even order, side cut n-icons have at least two "
                      "discontinuous surfaces and\n");
      fprintf(stderr, "one continuous edge. As a non-faceted smooth model it "
                      "would be dual to the\n");
      fprintf(stderr, "Point Cut n-icon. Circuit patterns are dual to Point "
                      "Cut n-icons.\n");
    }
    else if (!is_even(ncon_order)) {
      fprintf(stderr,
              "The order of this n-icon, %d, is odd. It could be "
              "termed both a Point Cut\n",
              ncon_order);
      fprintf(stderr, "and a Side Cut depending on the poles. Odd order "
                      "n-icons have at least one\n");
      fprintf(stderr, "discontinuous surface and one discontinuous edge. As a "
                      "non-faceted smooth\n");
      fprintf(stderr,
              "model it would be self dual. Circuit patterns are self dual.\n");
    }
  }

  if (info && twist == 0) {
    fprintf(stderr, "\n");
    fprintf(stderr,
            "Since it is twisted 0 times, it is as a globe with %d latitudes\n",
            num_lats(ncon_order, point_cut));

    if (is_even(ncon_order) && point_cut) {
      fprintf(stderr, "It has a north and south pole\n");
      fprintf(
          stderr,
          "It is the polyhedral dual of the twist 0 side cut order %d n-icon\n",
          ncon_order);
    }
    else if (is_even(ncon_order) && !point_cut) {
      fprintf(stderr, "It has two polar caps\n");
      fprintf(stderr,
              "It is the polyhedral dual of the twist 0 point cut "
              "order %d n-icon\n",
              ncon_order);
    }
    else if (!is_even(ncon_order)) {
      fprintf(stderr, "It has a north polar cap and south pole\n");
      fprintf(stderr, "As a polyhedra it is self dual\n");
    }
  }

  if (info && hybrid) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Note: Hybrids have no twist 0\n");
  }

  int base_form = 1;
  if (hybrid)
    base_form--;

  last = (int)floor((double)ncon_order / 2);
  first = -last;
  if (is_even(ncon_order) && !hybrid)
    first++;

  if (info) {
    fprintf(stderr, "\n");
    fprintf(stderr,
            "Using coloring, there are %d distinct twists ranging from "
            "%d through +%d\n",
            ncon_order, first, last);
  }

  int mod_twist = twist % ncon_order;
  if (hybrid && mod_twist == 0)
    mod_twist = -1;
  int color_twist = mod_twist % last;
  if (is_even(ncon_order) && color_twist == -last)
    color_twist = std::abs(color_twist);
  if (mod_twist != color_twist) {
    if (mod_twist > last)
      color_twist = first + color_twist - 1;
    else if (mod_twist < first)
      color_twist = last + color_twist;
    else
      color_twist = mod_twist;
  }

  if (info) {
    if (color_twist != twist)
      fprintf(stderr, "This twist is the same as color twist %d\n",
              color_twist);
  }

  nonchiral = 0;
  if (is_even(ncon_order)) {
    if (hybrid)
      last = (int)floor((double)(ncon_order + 2) / 4);
    else
      last = (int)floor((double)ncon_order / 4);
    forms = (last * 2) + 1;
    first = -last;
    if ((!hybrid && ncon_order % 4 == 0) ||
        (hybrid && (ncon_order + 2) % 4 == 0)) {
      forms--;
      first++;
    }
    // no twist 0
    if (hybrid)
      forms--;
  }
  else {
    forms = ncon_order;
    last = (int)floor((double)forms / 2);
    first = -last;
  }

  chiral = last;
  if ((!hybrid && ncon_order % 4 == 0) ||
      (hybrid && (ncon_order + 2) % 4 == 0)) {
    nonchiral = 1;
    chiral--;
  }
  unique = forms - last - base_form + nonchiral;

  if (info) {
    fprintf(stderr, "\n");
    fprintf(stderr,
            "Without coloring, there are %d distinct twists ranging "
            "from %d through +%d\n",
            forms, first, last);
    fprintf(stderr, "   %d base model(twist 0)\n", base_form);
    fprintf(stderr, "   %d form%s chiral and %s mirror image%s\n", chiral,
            ((chiral == 1) ? " is" : "s are"),
            ((chiral == 1) ? "has a" : "have"), ((chiral == 1) ? "" : "s"));
    fprintf(stderr, "   %d form%s non-chiral\n", nonchiral,
            ((nonchiral == 1) ? " is" : "s are"));
    fprintf(stderr, "   %d unique form%s of order %d twisted n-icon%s\n\n",
            unique, ((unique == 1) ? "" : "s"), ncon_order,
            ((unique == 1) ? "" : "s"));
  }

  int tmp_forms = forms;
  if (is_even(ncon_order))
    tmp_forms = ncon_order / 2;

  // also correct for no base twist 0 in hybrids
  int base_twist = twist;
  if (base_twist > 0)
    while (base_twist > last) {
      base_twist -= tmp_forms;
      if (hybrid && (base_twist <= 0 && base_twist + tmp_forms > 0))
        base_twist--;
    }
  else if (base_twist < 0)
    while (base_twist < first) {
      base_twist += tmp_forms;
      if (hybrid && (base_twist >= 0 && base_twist - tmp_forms < 0))
        base_twist++;
    }

  // using abs(base_twist) so much it might as well be stored
  int posi_twist = std::abs(base_twist);

  sd.nonchiral =
      (((chiral == 0 || base_twist == 0 ||
         (!hybrid && is_even(ncon_order) && (ncon_order % 4 == 0) &&
          base_twist == last) ||
         ((hybrid && (ncon_order + 2) % 4 == 0) && base_twist == last)))
           ? true
           : false);

  if (info) {
    if (twist < first || twist > last)
      fprintf(stderr, "It is the same as an n-icon of order %d twisted %d %s\n",
              ncon_order, base_twist, ((base_twist == 1) ? "time" : "times"));

    if (sd.nonchiral)
      fprintf(stderr, "It is not chiral\n");
    else
      fprintf(
          stderr,
          "It is the mirror image of an n-icon of order %d twisted %d time%s\n",
          ncon_order, -base_twist, ((-base_twist == 1) ? "" : "s"));
  }

  if (info) {
    fprintf(stderr, "\n");
    double twist_angle = ((double)360 / ncon_order) * twist;
    if (hybrid) {
      if (twist_angle > 0)
        twist_angle -= (double)360 / (ncon_order * 2);
      else
        twist_angle += (double)360 / (ncon_order * 2);
    }
    fprintf(stderr, "It has an absolute twist of %lf degrees, %lf radians\n",
            twist_angle, deg2rad(twist_angle));

    if (twist != base_twist) {
      twist_angle = ((double)360 / ncon_order) * base_twist;
      if (hybrid) {
        if (twist_angle > 0)
          twist_angle -= (double)360 / (ncon_order * 2);
        else
          twist_angle += (double)360 / (ncon_order * 2);
      }
      fprintf(stderr, "which is the same as one at %lf degrees, %lf radians\n",
              twist_angle, deg2rad(twist_angle));
    }
  }

  // if ncon_info() is called multiple times, this keeps surface table from
  // being called only once
  if (!surface_table.size())
    build_surface_table(surface_table, last, ncon_order, point_cut, hybrid);

  sd.c_surfaces = 0;
  sd.c_edges = 0;
  sd.d_surfaces = 0;
  sd.d_edges = 0;
  sd.total_surfaces = 0;
  // int total_edges = 0;

  sd.ncon_case2 = false;
  int case1_twist = 0;

  if (posi_twist == 0) {
    sd.c_surfaces = (int)floor((double)ncon_order + 1) / 2;
    if (is_even(ncon_order) && !point_cut)
      sd.c_surfaces++;
    sd.c_edges = sd.c_surfaces - 1;

    if (!is_even(ncon_order) || !point_cut) {
      sd.c_surfaces--;
      sd.d_surfaces++;
      if (is_even(ncon_order)) {
        sd.c_surfaces--;
        sd.d_surfaces++;
      }
    }
    sd.d_edges = 0;
  }
  else if (hybrid) {
    find_surface_count(surface_table, posi_twist, sd.c_surfaces, sd.ncon_case2,
                       case1_twist);
    sd.d_surfaces = 1;
    sd.c_edges = sd.c_surfaces + sd.d_surfaces - 1;
    sd.d_edges = 1;
  }
  else if (is_even(ncon_order) && point_cut) {
    find_surface_count(surface_table, posi_twist, sd.c_surfaces, sd.ncon_case2,
                       case1_twist);
    if (sd.c_surfaces == 0)
      sd.c_surfaces = 1;
    sd.d_surfaces = 0;
    sd.c_edges = sd.c_surfaces + sd.d_surfaces - 1;
    sd.d_edges = 2;
  }
  else if (is_even(ncon_order) && !point_cut) {
    find_surface_count(surface_table, posi_twist, sd.c_surfaces, sd.ncon_case2,
                       case1_twist);
    sd.d_surfaces = 2;
    sd.c_edges = sd.c_surfaces + sd.d_surfaces - 1;
    sd.d_edges = 0;
  }
  else if (!is_even(ncon_order)) {
    find_surface_count(surface_table, posi_twist, sd.c_surfaces, sd.ncon_case2,
                       case1_twist);
    sd.c_edges = sd.c_surfaces;
    sd.d_surfaces = 1;
    sd.d_edges = 1;
  }
  sd.total_surfaces = sd.c_surfaces + sd.d_surfaces;
  // total_edges = sd.c_edges + sd.d_edges;

  if (info) {
    fprintf(stderr, "\n");
    fprintf(stderr,
            "Treated as conic, it has a total of %d surface%s and a "
            "total of %d edge%s\n",
            sd.total_surfaces, ((sd.total_surfaces == 1) ? "" : "s"),
            sd.c_edges + sd.d_edges,
            ((sd.c_edges + sd.d_edges == 1) ? "" : "s"));
    fprintf(stderr, "   %d surface%s continuous\n", sd.c_surfaces,
            ((sd.c_surfaces == 1) ? " is" : "s are"));
    fprintf(stderr, "   %d surface%s discontinuous\n", sd.d_surfaces,
            ((sd.d_surfaces == 1) ? " is" : "s are"));
    fprintf(stderr, "   %d edge%s continuous\n", sd.c_edges,
            ((sd.c_edges == 1) ? " is" : "s are"));
    fprintf(stderr, "   %d edge%s discontinuous\n", sd.d_edges,
            ((sd.d_edges == 1) ? " is" : "s are"));
  }

  if (info) {
    if ((sd.c_surfaces > 1) ||
        (sd.c_surfaces > 0 &&
         (!(is_even(ncon_order) && point_cut) || hybrid))) {
      fprintf(stderr, "\n");
      if (!sd.ncon_case2) {
        fprintf(stderr, "This is a Case 1 n-icon. Surfaces cannot be colored "
                        "based on an earlier twist\n");
      }
      else {
        char ntype = 'p';
        if (hybrid)
          ntype = 'h';
        else if (is_even(ncon_order) && !point_cut)
          ntype = 's';
        fprintf(stderr, "This is a Case 2 n-icon. Surfaces can be colored "
                        "based on an earlier twist\n");
        char sign = (twist > 0) ? '+' : '-';
        fprintf(stderr,
                "It can be derived by twisting N%d%cT%d%c %c%d increments\n",
                ncon_order, sign, case1_twist, ntype, sign,
                posi_twist - case1_twist);
      }
    }
    fprintf(stderr, "\n");
  }
}

void color_uncolored_faces(Geometry &geom, const ncon_opts &opts)
{
  for (unsigned int i = 0; i < geom.faces().size(); i++) {
    if (!(geom.colors(FACES).get(i)).is_set())
      geom.colors(FACES).set(i, opts.face_default_color);
  }
}

void color_uncolored_edges(Geometry &geom, const ncon_opts &opts)
{
  const vector<vector<int>> &edges = geom.edges();
  const vector<Vec3d> &verts = geom.verts();

  for (unsigned int i = 0; i < edges.size(); i++) {
    if (!(geom.colors(EDGES).get(i)).is_set())
      geom.colors(EDGES).set(i, opts.edge_default_color);
  }

  for (unsigned int i = 0; i < verts.size(); i++) {
    if (!(geom.colors(VERTS).get(i)).is_set())
      geom.colors(VERTS).set(i, opts.edge_default_color);
  }
}

void color_unused_edges(Geometry &geom, const Color &unused_edge_color)
{
  geom.add_missing_impl_edges();
  for (unsigned int i = 0; i < geom.edges().size(); i++) {
    if (!(geom.colors(EDGES).get(i)).is_set())
      geom.colors(EDGES).set(i, unused_edge_color);
  }

  for (unsigned int i = 0; i < geom.verts().size(); i++) {
    if (!(geom.colors(VERTS).get(i)).is_set())
      geom.colors(VERTS).set(i, unused_edge_color);
  }
}

void reassert_colored_verts(Geometry &geom, const Color &default_color,
                            const Color &unused_edge_color)
{
  const vector<vector<int>> &edges = geom.edges();

  for (unsigned int i = 0; i < geom.edges().size(); i++) {
    Color c = geom.colors(EDGES).get(i);
    // if it is an index or it is a value and is completely opaque
    // and not the unused edge color
    if (c.is_set() && (c.is_index() || !c.get_transparency()) &&
        c != unused_edge_color) {
      int v1 = edges[i][0];
      int v2 = edges[i][1];
      Color cv1 = geom.colors(VERTS).get(v1);
      Color cv2 = geom.colors(VERTS).get(v2);
      // let default color be dominant
      if (cv1.is_set() && cv1 != default_color)
        geom.colors(VERTS).set(v1, c);
      if (cv2.is_set() && cv2 != default_color)
        geom.colors(VERTS).set(v2, c);
    }
  }
}

void unset_marked_edges(Geometry &geom)
{
  const vector<vector<int>> &edges = geom.edges();
  const vector<Vec3d> &verts = geom.verts();

  // vector<int> deleted_edges;
  for (unsigned int i = 0; i < edges.size(); i++) {
    Color c = geom.colors(EDGES).get(i);
    if (c.is_maximum_index()) {
      geom.colors(EDGES).set(i, Color());
      // deleted_edges.push_back(i);
    }
  }

  for (unsigned int i = 0; i < verts.size(); i++) {
    Color c = geom.colors(VERTS).get(i);
    if (c.is_maximum_index()) {
      geom.colors(VERTS).set(i, Color());
    }
  }
}

// if map index or invisible, alpha cannot be changed
Color set_opacity(Color &c, const int a)
{
  // in n_icons, if fail just return the input
  Color c_bak = c;
  if (!c.set_alpha(a))
    c = c_bak;
  return c;
}

// safe transparency setting. don't allow alpha set if
// color is index
// color is invisible
// -T or -U have not been set
// c is not changed
void set_vert_color(Geometry &geom, const int i, Color c,
                    const int opacity = -1)
{
  if (opacity > -1)
    c = set_opacity(c, opacity);
  geom.colors(VERTS).set(i, c);
}

void set_edge_col(Geometry &geom, const int i, Color c, const int opacity = -1)
{
  if (opacity > -1)
    c = set_opacity(c, opacity);
  geom.colors(EDGES).set(i, c);
}

void set_face_color(Geometry &geom, const int i, Color c,
                    const int opacity = -1)
{
  if (opacity > -1)
    c = set_opacity(c, opacity);
  geom.colors(FACES).set(i, c);
}

void set_edge_and_verts_col(Geometry &geom, const int i, Color c,
                            const int opacity = -1)
{
  set_edge_col(geom, i, c, opacity);

  set_vert_color(geom, geom.edges()[i][0], c, opacity);
  set_vert_color(geom, geom.edges()[i][1], c, opacity);
}

void set_edge_color(Geometry &geom, const int i, Color c,
                    const int opacity = -1)
{
  set_edge_and_verts_col(geom, i, c, opacity);
}

// point cut is not from opts
void ncon_edge_coloring(Geometry &geom, const vector<edgeList *> &edge_list,
                        const vector<poleList *> &pole,
                        map<int, pair<int, int>> &edge_color_table,
                        const bool point_cut_calc, const ncon_opts &opts)
{
  const vector<vector<int>> &edges = geom.edges();
  const vector<Vec3d> &verts = geom.verts();

  // special situation
  bool pc = opts.point_cut;
  if (opts.build_method == 3 && opts.hybrid && opts.angle_is_side_cut) {
    if (pc != point_cut_calc)
      pc = point_cut_calc;
  }
  if (opts.radius_inversion)
    pc = !pc;

  int opq = 255;

  if (opts.edge_coloring_method == 's') {
    int circuit_count = 0;

    for (auto i : edge_list) {
      int j = i->edge_no;
      int lat = i->lat;

      if (i->lat < 0) {
        set_edge_color(geom, j, Color(), opts.edge_opacity);
      }
      else {
        int col_idx = 0;
        if (i->rotate || (opts.hybrid && pc)) // front side
          col_idx = edge_color_table[lat].second;
        else
          col_idx = edge_color_table[lat].first;

        if (col_idx > circuit_count)
          circuit_count = col_idx;

        opq = opts.edge_pattern[col_idx % opts.edge_pattern.size()] == '1'
                  ? opts.edge_opacity
                  : 255;
        set_edge_color(geom, j, opts.edge_map.get_col(col_idx), opq);
      }
    }

    for (unsigned int i = 0; i < 2; i++) {
      if (pole[i]->idx > -1) {
        int lat = pole[i]->lat;
        if (lat < 0)
          continue;
        int col_idx = edge_color_table[lat].second;
        opq = opts.edge_pattern[col_idx % opts.edge_pattern.size()] == '1'
                  ? opts.edge_opacity
                  : 255;
        set_vert_color(geom, pole[i]->idx, opts.edge_map.get_col(col_idx), opq);
      }
    }

    if (opts.info && !opts.symmetric_coloring)
      fprintf(stderr, "%d edge circuit%s found\n", circuit_count,
              (circuit_count > 1 ? "s were" : " was"));
  }
  else if (opts.edge_coloring_method == 'l') {
    for (auto i : edge_list) {
      int j = i->edge_no;
      if (i->lat < 0) {
        set_edge_color(geom, j, Color(), opts.edge_opacity);
      }
      else {
        int lat = i->lat;
        opq = opts.edge_pattern[lat % opts.edge_pattern.size()] == '1'
                  ? opts.edge_opacity
                  : 255;
        set_edge_color(geom, j, opts.edge_map.get_col(lat), opq);
      }
    }

    for (unsigned int i = 0; i < 2; i++) {
      if (pole[i]->idx > -1) {
        int lat = pole[i]->lat;
        opq = opts.edge_pattern[lat % opts.edge_pattern.size()] == '1'
                  ? opts.edge_opacity
                  : 255;
        set_vert_color(geom, pole[i]->idx, opts.edge_map.get_col(lat), opq);
      }
    }
  }
  else if (opts.edge_coloring_method == 'm') {
    for (auto i : edge_list) {
      int j = i->edge_no;
      if (i->lon < 0) {
        set_edge_color(geom, j, Color(), opts.edge_opacity);
      }
      else {
        int lon = i->lon;
        opq = opts.edge_pattern[lon % opts.edge_pattern.size()] == '1'
                  ? opts.edge_opacity
                  : 255;
        set_edge_color(geom, j, opts.edge_map.get_col(lon), opq);
      }
    }

    // poles don't have any longitude
    for (unsigned int i = 0; i < 2; i++) {
      if (pole[i]->idx > -1) {
        set_vert_color(geom, pole[i]->idx, opts.edge_default_color,
                       opts.edge_opacity);
      }
    }
  }
  else if (opts.edge_coloring_method == 'b') {
    for (auto i : edge_list) {
      int j = i->edge_no;
      if (i->lat < 0)
        set_edge_color(geom, j, Color(), opts.edge_opacity);
      else {
        int n = -1;
        if ((is_even(i->lat) && is_even(i->lon)) ||
            (!is_even(i->lat) && !is_even(i->lon)))
          n = 0;
        else
          n = 1;

        opq = opts.edge_pattern[n % opts.edge_pattern.size()] == '1'
                  ? opts.edge_opacity
                  : 255;
        set_edge_color(geom, j, opts.edge_map.get_col(n), opq);
      }
    }

    // poles will be colored based North/South
    for (unsigned int i = 0; i < 2; i++) {
      if (pole[i]->idx > -1) {
        opq = opts.edge_pattern[i % opts.edge_pattern.size()] == '1'
                  ? opts.edge_opacity
                  : 255;
        set_vert_color(geom, pole[i]->idx, opts.edge_map.get_col(i), opq);
      }
    }
  }
  else if (opts.edge_coloring_method == 'n') {
    // keep track of index beyond loop
    unsigned int k = 0;
    for (unsigned int i = 0; i < edge_list.size(); i++) {
      int j = edge_list[i]->edge_no;
      opq = opts.edge_pattern[i % opts.edge_pattern.size()] == '1'
                ? opts.edge_opacity
                : 255;
      set_edge_color(geom, j, opts.edge_map.get_col(i), opq);
      k++;
    }

    for (unsigned int i = 0; i < 2; i++) {
      int col_idx = k;
      if (pole[i]->idx > -1) {
        opq = opts.edge_pattern[k % opts.edge_pattern.size()] == '1'
                  ? opts.edge_opacity
                  : 255;
        set_vert_color(geom, pole[i]->idx, opts.edge_map.get_col(col_idx), opq);
      }
      k++;
    }
  }
  else if (strchr("xyz", opts.edge_coloring_method)) {
    for (auto i : edge_list) {
      int j = i->edge_no;
      double d = 0.0;
      for (int k : edges[j]) {
        if (opts.edge_coloring_method == 'x')
          d += verts[k][0];
        else if (opts.edge_coloring_method == 'y')
          d += verts[k][1];
        else if (opts.edge_coloring_method == 'z')
          d += verts[k][2];
      }

      // The opts.hybrid base portion will end up in +Z
      if (opts.hybrid && pc)
        d *= -1;

      Color c;
      int n = -1;
      if (d > 0.0)
        n = 0;
      else if (d < 0.0)
        n = 1;
      else {
        n = -1;
        opq = opts.edge_opacity;
        c = opts.edge_default_color;
      }

      if (n > -1) {
        opq = opts.edge_pattern[n % opts.edge_pattern.size()] == '1'
                  ? opts.edge_opacity
                  : 255;
        c = opts.edge_map.get_col(n);
      }
      set_edge_color(geom, j, c, opq);
    }

    // poles are on the axis
    for (unsigned int i = 0; i < 2; i++) {
      if (pole[i]->idx > -1) {
        set_vert_color(geom, pole[i]->idx, opts.edge_default_color,
                       opts.edge_opacity);
      }
    }
  }
  else if (opts.edge_coloring_method == 'o') {
    for (auto i : edge_list) {
      int j = i->edge_no;
      double dx = 0.0;
      double dy = 0.0;
      double dz = 0.0;
      for (int k : edges[j]) {
        dx += verts[k][0];
        dy += verts[k][1];
        dz += verts[k][2];
      }

      // The opts.hybrid base portion will end up in +Z
      if (opts.hybrid && pc)
        dz *= -1;

      Color c;
      int n = -1;
      // by octant number 1 to 8
      if ((dx > 0.0) && (dy > 0.0) && (dz > 0.0))
        n = 0;
      else if ((dx < 0.0) && (dy > 0.0) && (dz > 0.0))
        n = 1;
      else if ((dx < 0.0) && (dy < 0.0) && (dz > 0.0))
        n = 2;
      else if ((dx > 0.0) && (dy < 0.0) && (dz > 0.0))
        n = 3;
      else if ((dx > 0.0) && (dy > 0.0) && (dz < 0.0))
        n = 4;
      else if ((dx < 0.0) && (dy > 0.0) && (dz < 0.0))
        n = 5;
      else if ((dx < 0.0) && (dy < 0.0) && (dz < 0.0))
        n = 6;
      else if ((dx > 0.0) && (dy < 0.0) && (dz < 0.0))
        n = 7;
      else {
        n = -1;
        opq = opts.edge_opacity;
        c = opts.edge_default_color;
      }

      if (n > -1) {
        opq = opts.edge_pattern[n % opts.edge_pattern.size()] == '1'
                  ? opts.edge_opacity
                  : 255;
        c = opts.edge_map.get_col(n);
      }
      set_edge_color(geom, j, c, opq);
    }

    // poles are on the axis
    for (unsigned int i = 0; i < 2; i++) {
      if (pole[i]->idx > -1) {
        set_vert_color(geom, pole[i]->idx, opts.edge_default_color,
                       opts.edge_opacity);
      }
    }
  }
}

void ncon_face_coloring(Geometry &geom, const vector<faceList *> &face_list,
                        map<int, pair<int, int>> &face_color_table,
                        const bool point_cut_calc, const ncon_opts &opts)
{
  const vector<vector<int>> &faces = geom.faces();
  const vector<Vec3d> &verts = geom.verts();

  // special situation
  bool pc = opts.point_cut;
  if (opts.build_method == 3 && opts.hybrid && opts.angle_is_side_cut) {
    if (pc != point_cut_calc)
      pc = point_cut_calc;
  }
  if (opts.radius_inversion)
    pc = !pc;

  int opq = 255;

  if (opts.face_coloring_method == 's') {
    int circuit_count = 0;

    for (auto i : face_list) {
      int j = i->face_no;
      int lat = i->lat;

      if (i->lat < 0) {
        set_face_color(geom, j, Color(), opts.face_opacity);
      }
      else {
        int col_idx = 0;
        if (i->rotate || (opts.hybrid && pc)) // front side
          col_idx = face_color_table[lat].second;
        else
          col_idx = face_color_table[lat].first;

        if (col_idx > circuit_count)
          circuit_count = col_idx;

        opq = opts.face_pattern[col_idx % opts.face_pattern.size()] == '1'
                  ? opts.face_opacity
                  : 255;
        set_face_color(geom, j, opts.face_map.get_col(col_idx), opq);
      }
    }

    if (opts.info && !opts.symmetric_coloring) {
      circuit_count++;
      fprintf(stderr, "%d face circuit%s found\n", circuit_count,
              (circuit_count > 1 ? "s were" : " was"));
    }
  }
  else if (opts.face_coloring_method == 'l') {
    for (auto i : face_list) {
      int j = i->face_no;
      if (i->lat < 0) {
        set_face_color(geom, j, Color(), opts.face_opacity);
      }
      else {
        int lat = i->lat;
        opq = opts.face_pattern[lat % opts.face_pattern.size()] == '1'
                  ? opts.face_opacity
                  : 255;
        set_face_color(geom, j, opts.face_map.get_col(lat), opq);
      }
    }
  }
  else if (opts.face_coloring_method == 'm') {
    for (auto i : face_list) {
      int j = i->face_no;
      if (i->lon < 0) {
        set_face_color(geom, j, Color(), opts.face_opacity);
      }
      else {
        int lon = i->lon;
        opq = opts.face_pattern[lon % opts.face_pattern.size()] == '1'
                  ? opts.face_opacity
                  : 255;
        set_face_color(geom, j, opts.face_map.get_col(lon), opq);
      }
    }
  }
  else if (opts.face_coloring_method == 'b') {
    for (auto i : face_list) {
      int j = i->face_no;
      if (i->lat < 0)
        set_face_color(geom, j, Color(), opts.face_opacity);
      else {
        int n = -1;
        if ((is_even(i->lat) && is_even(i->lon)) ||
            (!is_even(i->lat) && !is_even(i->lon)))
          n = 0;
        else
          n = 1;

        opq = opts.face_pattern[n % opts.face_pattern.size()] == '1'
                  ? opts.face_opacity
                  : 255;
        set_face_color(geom, j, opts.face_map.get_col(n), opq);
      }
    }
  }
  else if (opts.face_coloring_method == 'n') {
    for (unsigned int i = 0; i < face_list.size(); i++) {
      int j = face_list[i]->face_no;
      opq = opts.face_pattern[i % opts.face_pattern.size()] == '1'
                ? opts.face_opacity
                : 255;
      set_face_color(geom, j, opts.face_map.get_col(i), opq);
    }
  }
  else if (strchr("xyz", opts.face_coloring_method)) {
    for (auto i : face_list) {
      int j = i->face_no;
      double d = 0.0;
      for (int k : faces[j]) {
        if (opts.face_coloring_method == 'x')
          d += verts[k][0];
        else if (opts.face_coloring_method == 'y')
          d += verts[k][1];
        else if (opts.face_coloring_method == 'z')
          d += verts[k][2];
      }

      // The opts.hybrid base portion will end up in +Z
      if (opts.hybrid && pc)
        d *= -1;

      Color c;
      int n = -1;
      if (d > 0.0)
        n = 0;
      else if (d < 0.0)
        n = 1;
      else {
        n = -1;
        opq = opts.face_opacity;
        c = opts.face_default_color;
      }

      if (n > -1) {
        opq = opts.face_pattern[n % opts.face_pattern.size()] == '1'
                  ? opts.face_opacity
                  : 255;
        c = opts.face_map.get_col(n);
      }
      set_face_color(geom, j, c, opq);
    }
  }
  else if (opts.face_coloring_method == 'o') {
    for (auto i : face_list) {
      int j = i->face_no;
      double dx = 0.0;
      double dy = 0.0;
      double dz = 0.0;
      for (int k : faces[j]) {
        dx += verts[k][0];
        dy += verts[k][1];
        dz += verts[k][2];
      }

      // The opts.hybrid base portion will end up in +Z
      if (opts.hybrid && pc)
        dz *= -1;

      // by octant number 1 to 8
      Color c;
      int n = -1;
      if ((dx > 0.0) && (dy > 0.0) && (dz > 0.0))
        n = 0;
      else if ((dx < 0.0) && (dy > 0.0) && (dz > 0.0))
        n = 1;
      else if ((dx < 0.0) && (dy < 0.0) && (dz > 0.0))
        n = 2;
      else if ((dx > 0.0) && (dy < 0.0) && (dz > 0.0))
        n = 3;
      else if ((dx > 0.0) && (dy > 0.0) && (dz < 0.0))
        n = 4;
      else if ((dx < 0.0) && (dy > 0.0) && (dz < 0.0))
        n = 5;
      else if ((dx < 0.0) && (dy < 0.0) && (dz < 0.0))
        n = 6;
      else if ((dx > 0.0) && (dy < 0.0) && (dz < 0.0))
        n = 7;
      else {
        n = -1;
        opq = opts.face_opacity;
        c = opts.face_default_color;
      }

      if (n > -1) {
        opq = opts.face_pattern[n % opts.face_pattern.size()] == '1'
                  ? opts.face_opacity
                  : 255;
        c = opts.face_map.get_col(n);
      }
      set_face_color(geom, j, c, opq);
    }
  }
}

vector<int> find_adjacent_face_idx_in_channel(
    const Geometry &geom, const int face_idx,
    const vector<vector<int>> &bare_implicit_edges,
    map<vector<int>, vector<int>> &faces_by_edge, const bool prime)
{
  const vector<vector<int>> &faces = geom.faces();

  vector<int> face_idx_ret;
  vector<vector<int>> adjacent_edges;

  vector<int> face = faces[face_idx];
  unsigned int sz = face.size();
  for (unsigned int i = 0; i < sz; i++) {
    vector<int> edge(2);
    edge[0] = face[i];
    edge[1] = face[(i + 1) % sz];
    if (find_edge_in_edge_list(bare_implicit_edges, edge) > -1) {
      adjacent_edges.push_back(make_edge(edge[0], edge[1]));
    }
  }

  vector<int> adjacent_face_idx;
  for (auto &adjacent_edge : adjacent_edges) {
    vector<int> face_idx_tmp = faces_by_edge[adjacent_edge];
    for (int &j : face_idx_tmp) {
      if (j != face_idx)
        adjacent_face_idx.push_back(j);
    }
  }

  // the first time we "prime" so we return faces in both directions. the second
  // face becomes the "stranded" face
  // there may be faces which would get pinched off (stranded), so there may be
  // more than one
  for (int &i : adjacent_face_idx) {
    if (prime || !(geom.colors(FACES).get(i)).is_set()) {
      face_idx_ret.push_back(i);
    }
  }

  return face_idx_ret;
}

// flood_fill_count is changed
int set_face_colors_by_adjacent_face(
    Geometry &geom, const int start, const Color &c, const int opq,
    const int flood_fill_stop, int &flood_fill_count,
    const vector<vector<int>> &bare_implicit_edges,
    map<vector<int>, vector<int>> &faces_by_edge)
{
  if (flood_fill_stop && (flood_fill_count >= flood_fill_stop))
    return 0;

  vector<int> stranded_faces;

  vector<int> face_idx = find_adjacent_face_idx_in_channel(
      geom, start, bare_implicit_edges, faces_by_edge, true);
  while (face_idx.size()) {
    for (unsigned int i = 0; i < face_idx.size(); i++) {
      if (flood_fill_stop && (flood_fill_count >= flood_fill_stop))
        return 0;
      set_face_color(geom, face_idx[i], c, opq);
      flood_fill_count++;
      if (i > 0)
        stranded_faces.push_back(face_idx[i]);
    }
    face_idx = find_adjacent_face_idx_in_channel(
        geom, face_idx[0], bare_implicit_edges, faces_by_edge, false);
  }

  // check if stranded faces
  for (unsigned int i = 0; i < stranded_faces.size(); i++) {
    face_idx = find_adjacent_face_idx_in_channel(
        geom, stranded_faces[i], bare_implicit_edges, faces_by_edge, false);
    while (face_idx.size()) {
      for (unsigned int i = 0; i < face_idx.size(); i++) {
        if (flood_fill_stop && (flood_fill_count >= flood_fill_stop))
          return 0;
        set_face_color(geom, face_idx[i], c, opq);
        flood_fill_count++;
        if (i > 0)
          stranded_faces.push_back(face_idx[i]);
      }
      face_idx = find_adjacent_face_idx_in_channel(
          geom, face_idx[0], bare_implicit_edges, faces_by_edge, false);
    }
  }

  return (flood_fill_stop ? 1 : 0);
}

void fill_bare_implicit_edges(const Geometry &geom,
                              vector<vector<int>> &bare_implicit_edges)
{
  const vector<vector<int>> &edges = geom.edges();
  vector<vector<int>> implicit_edges;
  geom.get_impl_edges(implicit_edges);

  for (auto &implicit_edge : implicit_edges) {
    if (find_edge_in_edge_list(edges, implicit_edge) < 0)
      bare_implicit_edges.push_back(implicit_edge);
  }
}

void fill_faces_by_edge(const Geometry &geom,
                        map<vector<int>, vector<int>> &faces_by_edge)
{
  const vector<vector<int>> &faces = geom.faces();

  for (unsigned int i = 0; i < faces.size(); i++) {
    unsigned int sz = faces[i].size();
    for (unsigned int j = 0; j < sz; j++)
      faces_by_edge[make_edge(faces[i][j], faces[i][(j + 1) % sz])].push_back(
          i);
  }
}

int ncon_face_coloring_by_adjacent_face(Geometry &geom,
                                        const vector<faceList *> &face_list,
                                        const ncon_opts &opts)
{
  bool debug = false;

  int flood_fill_count = 0;
  int ret = 0;

  // all face colors need to be cleared
  Coloring clrng(&geom);
  clrng.f_one_col(Color());

  vector<vector<int>> bare_implicit_edges;
  fill_bare_implicit_edges(geom, bare_implicit_edges);

  map<vector<int>, vector<int>> faces_by_edge;
  fill_faces_by_edge(geom, faces_by_edge);

  int map_count = 0;
  int lat = 0;
  int lon = opts.longitudes.front() / 2 - 1;
  unsigned int sz = 0;
  do {
    bool painted = false;
    vector<int> idx = find_face_by_lat_lon(face_list, lat, lon);
    sz = idx.size();

    Color c = map_count;

    for (int j : idx) {
      int f_idx = face_list[j]->face_no;
      if ((geom.colors(FACES).get(f_idx)).is_set())
        continue;
      if (opts.flood_fill_stop && (flood_fill_count >= opts.flood_fill_stop))
        break;

      set_face_color(geom, f_idx, c, 255);
      flood_fill_count++;
      ret = set_face_colors_by_adjacent_face(
          geom, f_idx, c, 255, opts.flood_fill_stop, flood_fill_count,
          bare_implicit_edges, faces_by_edge);

      painted = true;
    }
    if (painted)
      map_count++;

    if (opts.flood_fill_stop && (flood_fill_count >= opts.flood_fill_stop))
      break;

    lat++;
  } while (sz);

  if (debug) {
    fprintf(stderr, "flood fill face color map:\n");
    lon = opts.longitudes.front() / 2;
    for (int l = 0; l < 2; l++) {
      lat = 0;
      lon -= l;
      do {
        vector<int> idx = find_face_by_lat_lon(face_list, lat, lon);
        sz = idx.size();
        if (sz) {
          int f_idx = face_list[idx[0]]->face_no;
          int k = geom.colors(FACES).get(f_idx).get_index();
          fprintf(stderr, "%d ", k);
        }
        lat++;
      } while (sz);
      fprintf(stderr, "\n");
    }

    /* useful?
    fprintf(stderr,"flood fill face latitudes:\n");
    lon = opts.longitudes.front()/2;
    for (int l=0;l<2;l++) {
       lat = 0;
       //lon -= l;
       do {
          vector<int> idx = find_face_by_lat_lon(face_list,lat,lon);
          sz = idx.size();
          if (sz) {
             vector<int> face_idx;
             int lat2 = -1;
             if (l==1) {
                face_idx = find_adjacent_face_idx_in_channel(geom, idx[0],
    bare_implicit_edges, faces_by_edge, true);
                sz = face_idx.size();
                //fprintf(stderr,"face idx0 = %d, face idx1 =
    %d\n",face_idx[0],face_idx[1]);
                //fprintf(stderr,"sz = %d lat0 = %d lat1 =
    %d\n",sz,face_list[face_idx[0]]->lat,face_list[face_idx[1]]->lat);
                if (sz) {
                   if (face_idx[1] > (int)face_list.size())
                      face_idx = find_adjacent_face_idx_in_channel(geom,
    face_idx[1], bare_implicit_edges, faces_by_edge, true);
                   lat2 = face_list[face_idx[1]]->lat;
                }
             }
             else
                lat2 = face_list[idx[0]]->lat;
             fprintf(stderr,"%d ",lat2);
          }
          lat++;
       } while( sz );
       fprintf(stderr,"\n");
    }
    */
  }

  vector<int> color_table;
  for (int i = 0; i < map_count; i++)
    color_table.push_back(i);

  // colors for method 2 and 3 are reversed from method 1
  // accept when it is a 1/2 opts.twist
  // don't do this with flood fill stop as the colors change as fill gets larger
  if (!opts.flood_fill_stop) {
    if (opts.mod_twist != 0 && (opts.ncon_order / opts.mod_twist == 2)) {
      reverse(color_table.begin(), color_table.end());
    }
  }

  // equivalent symmetric coloring. odd n_icons are not affected
  if (opts.symmetric_coloring && is_even(opts.ncon_order)) {
    unsigned int sz = color_table.size();
    unsigned int j = sz - 1;
    for (unsigned int i = 0; i < sz / 2; i++) {
      int idx = color_table[i];
      color_table[j] = idx;
      j--;
    }
  }

  // for (unsigned int i=0;i<color_table.size();i++)
  //   fprintf(stderr,"color table = %d\n",color_table[i]);

  // resolve map indexes
  for (unsigned int i = 0; i < geom.faces().size(); i++) {
    Color c_idx = geom.colors(FACES).get(i);
    if (c_idx.is_index()) {
      int idx = c_idx.get_index();
      idx = color_table[idx];
      Color c = opts.face_map.get_col(idx);
      int opq = 255;
      if (opts.face_opacity > -1)
        opq = opts.face_pattern[idx % opts.face_pattern.size()] == '1'
                  ? opts.face_opacity
                  : 255;
      set_face_color(geom, i, c, opq);
    }
  }

  if (!opts.flood_fill_stop && opts.info)
    fprintf(stderr, "%d face circuit%s found\n", map_count,
            (map_count > 1 ? "s were" : " was"));
  return ret;
}

void ncon_edge_coloring_by_adjacent_edge(Geometry &geom,
                                         const vector<edgeList *> &edge_list,
                                         const vector<poleList *> &pole,
                                         const ncon_opts &opts)
{
  bool debug = false;

  bool pc = opts.point_cut;
  if (opts.radius_inversion)
    pc = !pc;

  vector<vector<int>> edges;
  vector<int> edge_no;
  if (!opts.hybrid) {
    for (auto i : edge_list) {
      edges.push_back(geom.edges(i->edge_no));
      edge_no.push_back(i->edge_no);
    }
  }
  else {
    // using edge_list doesn't work for hybrids
    for (unsigned int i = 0; i < geom.edges().size(); i++) {
      Color c = geom.colors(EDGES).get(i);
      // need to include invisible edges or crash
      if (c.is_maximum_index() || c.is_invisible()) {
        edges.push_back(geom.edges(i));
        edge_no.push_back(i);
      }
    }
  }

  // unset colors
  for (unsigned int i = 0; i < geom.edges().size(); i++) {
    Color c = geom.colors(EDGES).get(i);
    if (c.is_maximum_index())
      set_edge_color(geom, i, Color(), 255);
  }

  int map_count = 0;
  int circuit_count = 0;

  // if there is a north pole increment the color map
  // not a circuit, so don't increment circuit count
  // if (is_even(opts.ncon_order) && pc && opts.mod_twist == 0 &&
  // !opts.double_sweep)
  if (pole[0]->idx != -1 && !opts.double_sweep && !opts.hybrid)
    map_count++;

  // edge lats start at 1
  int lat = 1;
  int lon = opts.longitudes.front() / 2 - 1;
  unsigned int sz = 0;
  do {
    bool painted = false;
    vector<int> idx = find_edge_by_lat_lon(edge_list, lat, lon);
    sz = idx.size();

    Color c = map_count;

    for (unsigned int j = 0; j < sz; j++) {
      int e_idx = edge_list[idx[j]]->edge_no;
      vector<int> edge = edges[e_idx];
      vector<int> es = find_edges_with_vertex(edges, edge[0]);
      vector<int> es_tmp = find_edges_with_vertex(edges, edge[1]);
      es.insert(es.end(), es_tmp.begin(), es_tmp.end());

      for (int k = (int)es.size() - 1; k >= 0; k--) {
        if ((geom.colors(EDGES).get(edge_no[es[k]])).is_set())
          es.erase(es.begin() + k);
      }

      while (es.size()) {
        // continue to color in one direction, then the other
        for (int e : es) {
          set_edge_color(geom, edge_no[e], c, 255);
          painted = true;
        }

        vector<int> es_next;
        for (int e : es) {
          edge = edges[e];
          if (!es_next.size())
            es_next = find_edges_with_vertex(edges, edge[0]);
          else {
            es_tmp = find_edges_with_vertex(edges, edge[0]);
            es_next.insert(es_next.end(), es_tmp.begin(), es_tmp.end());
          }

          es_tmp = find_edges_with_vertex(edges, edge[1]);
          es_next.insert(es_next.end(), es_tmp.begin(), es_tmp.end());
        }

        es = es_next;
        for (int k = (int)es.size() - 1; k >= 0; k--) {
          if ((geom.colors(EDGES).get(edge_no[es[k]])).is_set())
            es.erase(es.begin() + k);
        }
      }
    }

    if (painted) {
      map_count++;
      circuit_count++;
    }

    lat++;
  } while (sz);

  // north pole, even opts.point_cut only
  // if there is a pole in build_method 3 it is meant to be there
  int n = (opts.build_method == 2 && opts.d != 1) ? opts.ncon_order / 2
                                                  : opts.ncon_order;
  if ((pole[0]->idx != -1) && ((is_even(n) && pc) || opts.build_method == 3) &&
      !opts.hybrid) {
    int opq = opts.edge_pattern[0 % (int)opts.edge_pattern.size()] == '1'
                  ? opts.edge_opacity
                  : 255;
    Color c = opts.edge_map.get_col(0);
    set_vert_color(geom, pole[0]->idx, c, opq);
  }

  // south pole
  if ((pole[1]->idx != -1) && (!is_even(n) || (is_even(n) && pc)) &&
      !opts.hybrid) {
    int l = (opts.symmetric_coloring && is_even(opts.ncon_order)) ? 0 : lat - 1;
    int opq = opts.edge_pattern[l % (int)opts.edge_pattern.size()] == '1'
                  ? opts.edge_opacity
                  : 255;
    Color c = opts.edge_map.get_col(l);
    set_vert_color(geom, pole[1]->idx, c, opq);
  }

  // some models will have a stranded edge
  bool painted = false;
  for (auto &edge : edges) {
    int edge_no = find_edge_in_edge_list(geom.edges(), edge);
    if (!(geom.colors(EDGES).get(edge_no)).is_set()) {
      set_edge_color(geom, edge_no, Color(0), 255);
      painted = true;
    }
  }
  if (painted) {
    map_count++;
    circuit_count++;
  }

  vector<int> color_table;
  for (int i = 1; i < map_count; i++)
    color_table.push_back(i);

  // equivilent to z=1
  if (is_even(opts.ncon_order) && pc && !opts.hybrid && !opts.double_sweep &&
      ((opts.mod_twist != 0) && (opts.ncon_order % opts.mod_twist == 0)))
    color_table[color_table.size() - 1] = 0;

  // equivalent symmetric coloring. odd n_icons are not affected
  if (opts.symmetric_coloring && is_even(opts.ncon_order)) {
    sz = color_table.size();
    unsigned int j = sz - 1;
    for (unsigned int i = 0; i < sz / 2; i++) {
      int idx = color_table[i];
      color_table[j] = idx;
      j--;
    }
  }

  // for alignment to color map
  color_table.insert(color_table.begin(), 0);

  // colors for method 2 and 3 are reversed from method 1
  // accept when it is a 1/2 opts.twist
  // don't do this with flood fill stop as the colors change as fill gets larger
  if (!opts.flood_fill_stop) {
    if (opts.mod_twist != 0 && (opts.ncon_order / opts.mod_twist == 2)) {
      reverse(color_table.begin(), color_table.end());
    }
  }

  // note: edge_list is a problem for hybrids
  if (debug) {
    fprintf(stderr, "flood fill edge color map:\n");
    lon = opts.longitudes.front() / 2;
    for (int l = 0; l < 2; l++) {
      lon -= l;
      // edge lats start at 1
      lat = 1;
      do {
        vector<int> idx = find_edge_by_lat_lon(edge_list, lat, lon);
        sz = idx.size();
        if (sz) {
          int e_idx = edge_list[idx[0]]->edge_no;
          int k = geom.colors(EDGES).get(e_idx).get_index();
          fprintf(stderr, "%d ", k);
        }
        lat++;
      } while (sz);
      fprintf(stderr, "\n");
    }
  }

  // resolve map indexes
  for (unsigned int i = 0; i < edges.size(); i++) {
    Color c_idx = geom.colors(EDGES).get(edge_no[i]);
    if (c_idx.is_index()) {
      int idx = c_idx.get_index();
      idx = color_table[idx];
      Color c = opts.edge_map.get_col(idx);
      int opq = 255;
      if (opts.edge_opacity > -1)
        opq = opts.edge_pattern[idx % opts.edge_pattern.size()] == '1'
                  ? opts.edge_opacity
                  : 255;
      set_edge_color(geom, edge_no[i], c, opq);
    }
  }

  if (opts.info)
    fprintf(stderr, "%d edge circuit%s found\n", circuit_count,
            (circuit_count > 1 ? "s were" : " was"));
}

int ncon_face_coloring_by_compound(Geometry &geom,
                                   const vector<faceList *> &face_list,
                                   const vector<int> &caps,
                                   const ncon_opts &opts)
{
  int flood_fill_count = 0;
  int ret = 0;

  // all face indexes need to be cleared
  Coloring clrng(&geom);
  clrng.f_one_col(Color());

  int opq = 255;

  /*
     //test polygon_no
     for (unsigned int i=0;i<face_list.size();i++) {
        int face_no = face_list[i]->face_no;
        int polygon_no = face_list[i]->polygon_no;
        if (opts.face_opacity > -1)
           opq = opts.face_pattern[polygon_no%opts.face_pattern.size()] == '1' ?
     opts.face_opacity : 255;
        Color c = set_opacity(opts.face_map.get_col(polygon_no),opq);
        geom.colors(FACES).set(face_no,c);
     }
     return ret;
  */

  vector<pair<int, int>> polygon_table;

  int sz = 0;
  int lat = 0;
  int lon = opts.longitudes.front() / 2;
  if (!opts.hybrid)
    lon--;
  do {
    vector<int> idx = find_face_by_lat_lon(face_list, lat, lon);
    sz = idx.size();
    for (int j = 0; j < sz; j++) {
      int polygon_no = face_list[idx[j]]->polygon_no;
      polygon_table.push_back(make_pair(polygon_no, idx[j]));
    }
    lat++;
  } while (sz);

  sort(polygon_table.begin(), polygon_table.end());
  auto li = unique(polygon_table.begin(), polygon_table.end());
  polygon_table.erase(li, polygon_table.end());

  vector<vector<int>> bare_implicit_edges;
  fill_bare_implicit_edges(geom, bare_implicit_edges);

  map<vector<int>, vector<int>> faces_by_edge;
  fill_faces_by_edge(geom, faces_by_edge);

  int map_count = 0;
  int polygon_no_last = -1;
  for (auto polygons : polygon_table) {
    int polygon_no = polygons.first;
    int face_no = polygons.second;

    if ((geom.colors(FACES).get(face_no)).is_set())
      continue;
    if (opts.flood_fill_stop && (flood_fill_count >= opts.flood_fill_stop))
      return 0;

    Color c = opts.face_map.get_col(polygon_no);
    if (opts.face_opacity > -1)
      opq = opts.face_pattern[polygon_no % opts.face_pattern.size()] == '1'
                ? opts.face_opacity
                : 255;
    set_face_color(geom, face_no, c, opq);

    flood_fill_count++;
    ret = set_face_colors_by_adjacent_face(
        geom, face_no, c, opq, opts.flood_fill_stop, flood_fill_count,
        bare_implicit_edges, faces_by_edge);

    if (polygon_no != polygon_no_last)
      map_count++;
    polygon_no_last = polygon_no;
  }

  // if opts.build_method 1 or 2, twist 0, then end caps will not be colored
  // color them the first map color
  if (opts.build_method < 3 && opts.mod_twist == 0) {
    for (int cap : caps)
      geom.colors(FACES).set(face_list[cap]->face_no, opts.face_map.get_col(0));
  }

  if (!opts.flood_fill_stop && opts.info)
    fprintf(stderr, "%d compound part%s painted\n", map_count,
            (map_count > 1 ? "s were" : " was"));

  return ret;
}

void ncon_alternate_compound_coloring(Geometry &geom, const ncon_opts &opts)
{
  int opq = 255;

  Coloring clrng(&geom);
  clrng.add_cmap(opts.face_map.clone());

  // color by constituents
  clrng.f_parts(false);

  int part_count = -1;
  for (unsigned int i = 0; i < geom.faces().size(); i++) {
    Color c = geom.colors(FACES).get(i);
    if (c.is_index()) {
      int idx = c.get_index();
      if (idx > part_count)
        part_count = idx;

      c = opts.face_map.get_col(idx);
      if (opts.face_opacity > -1)
        opq = opts.face_pattern[idx % (int)opts.face_pattern.size()] == '1'
                  ? opts.face_opacity
                  : 255;
      set_face_color(geom, i, c, opq);
    }
  }

  if (opts.info) {
    part_count++;
    fprintf(stderr, "%d compound part%s painted\n", part_count,
            (part_count > 1 ? "s were" : " was"));
  }
}

/* unused function to line a face with edges of color col
void add_edges_to_face(Geometry &geom, vector<int> face, Color col)
{
   unsigned int sz = face.size();
   for (unsigned int i=0;i<sz;i++)
      geom.add_col_edge(make_edge(face[i],face[(i+1)%sz]),col);
}
*/

// average edge colors with simple RGB
Color ncon_average_edge_color(const vector<Color> &cols)
{
  /* average RGB colors */
  Vec4d col(0.0, 0.0, 0.0, 0.0);

  int j = 0;
  for (auto i : cols) {
    if (i.is_set() && i.is_visible_value()) {
      col += i.get_vec4d();
      j++;
    }
  }

  if (j)
    col /= j;
  return col;
}

// find poles by Y value. Set north and south if found
void find_poles(Geometry &geom, int &north, int &south, const ncon_opts &opts)
{
  // find poles
  north = -1;
  double y_north = std::numeric_limits<double>::min();
  south = -1;
  double y_south = std::numeric_limits<double>::max();
  for (unsigned int i = 0; i < geom.verts().size(); i++) {
    if (double_gt(geom.verts(i)[1], y_north, opts.eps)) {
      y_north = geom.verts(i)[1];
      north = i;
    }
    if (double_lt(geom.verts(i)[1], y_south, opts.eps)) {
      y_south = geom.verts(i)[1];
      south = i;
    }
  }
}

void ncon_edge_coloring_from_faces(Geometry &geom, const ncon_opts &opts)
{
  // save invisible edges
  vector<int> inv_edges;
  for (unsigned int i = 0; i < geom.edges().size(); i++) {
    if ((geom.colors(EDGES).get(i)).is_invisible())
      inv_edges.push_back(i);
  }

  // save invisible verts
  vector<int> inv_verts;
  for (unsigned int i = 0; i < geom.verts().size(); i++) {
    if ((geom.colors(VERTS).get(i)).is_invisible())
      inv_verts.push_back(i);
  }

  Coloring clrng;
  clrng.add_cmap(opts.face_map.clone());
  clrng.set_geom(&geom);
  clrng.e_from_adjacent(FACES);
  clrng.v_from_adjacent(FACES);

  // restore invisible edges
  for (int inv_edge : inv_edges)
    geom.colors(EDGES).set(inv_edge, Color::invisible);

  // restore invisible verts
  for (int inv_vert : inv_verts)
    geom.colors(VERTS).set(inv_vert, Color::invisible);
}

// mark edge circuits for methods 2 and 3 color by adjacent edge
// also for method 1, color by symmetry
void mark_edge_circuits(Geometry &geom, const vector<edgeList *> &edge_list)
{
  for (auto i : edge_list) {
    int j = i->edge_no;
    if (!geom.colors(EDGES).get(j).is_invisible())
      set_edge_color(geom, j, Color::maximum_index, 255);
  }
}

// for method 2, if indented edges are not shown, overwrite them as invisible
void set_indented_edges_invisible(Geometry &geom,
                                  const vector<edgeList *> &edge_list,
                                  const vector<poleList *> &pole)
{
  for (auto i : edge_list) {
    int j = i->edge_no;
    int lat = i->lat;
    if (lat == -1)
      set_edge_color(geom, j, Color::invisible, 255);
  }

  for (unsigned int i = 0; i < 2; i++) {
    if (pole[i]->idx > -1) {
      int lat = pole[i]->lat;
      if (lat == -1)
        set_vert_color(geom, pole[i]->idx, Color::invisible, 255);
    }
  }
}

void make_sequential_map(map<int, pair<int, int>> &color_table)
{
  unsigned int sz = color_table.size();

  // add two to avoid zeroes and make all elements negative.
  for (unsigned int i = 0; i < sz; i++) {
    color_table[i].first = -(color_table[i].first + 2);
    color_table[i].second = -(color_table[i].second + 2);
  }

  int k = 0;
  for (unsigned int i = 0; i < sz; i++) {
    int element = color_table[i].first;
    if (element < 0) {
      for (unsigned int j = i; j < sz; j++) {
        if (color_table[j].first == element)
          color_table[j].first = k;
      }
      for (unsigned int j = 0; j < sz; j++) {
        if (color_table[j].second == element)
          color_table[j].second = k;
      }
      k++;
    }
  }
}

// coding idea for circuit coloring furnished by Adrian Rossiter
void build_circuit_table(const int ncon_order, const int twist,
                         const bool hybrid, const bool symmetric_coloring,
                         map<int, int> &circuit_table)
{
  // use a double size polygon
  int n = 2 * ncon_order;
  int t = 2 * twist;

  if (hybrid)
    t--;

  int t_mult = symmetric_coloring ? 1 : 2;

  int d = gcd(n, t_mult * t);

  for (int i = 0; i <= d / 2; i++) {
    for (int j = i; j < n; j += d)
      circuit_table[j] = i;
    for (int j = d - i; j < n; j += d)
      circuit_table[j] = i;
  }
}

void build_color_table(map<int, pair<int, int>> &color_table,
                       const int ncon_order, const int twist, const bool hybrid,
                       const bool symmetric_coloring, const int increment)
{
  bool debug = false;

  // boolean for testing. this was face coloring options 't'
  bool sequential_colors = true;

  // for the color tables twist is made positive. The positive twist colors work
  // for negative twists.
  unsigned int n = ncon_order;
  unsigned int t = std::abs(twist);

  // for even n, twist 0, symmetric coloring does't happen, so temporarily twist
  // half way
  if (symmetric_coloring && is_even(n) && t == 0) {
    int mod_twist = t % n;
    if (!mod_twist)
      t = n / 2;
  }

  map<int, int> circuit_table;
  build_circuit_table(n, t, hybrid, symmetric_coloring, circuit_table);

  // for (unsigned int i=0;i<circuit_table.size();i++)
  //   fprintf(stderr,"circuit_table[%d] = %d\n",i,circuit_table[i]);

  if (debug)
    fprintf(stderr, "color table:\n");
  for (unsigned int i = 0; i < n; i++) {
    int position = 2 * i + increment;
    // when increment = -1, pos can equal -1, use last position
    if (position < 0)
      position = 2 * n + position;
    unsigned int circuit_idx = position % (int)circuit_table.size();
    color_table[i].first = circuit_table[circuit_idx];
    if (debug)
      fprintf(stderr, "%d ", color_table[i].first);
  }
  if (debug)
    fprintf(stderr, "\n");

  // first is back half.
  // second is front half.
  for (unsigned int i = 0; i < n; i++) {
    color_table[i].second = color_table[(t + i) % n].first;
    if (debug)
      fprintf(stderr, "%d ", color_table[i].second);
  }
  if (debug)
    fprintf(stderr, "\n");

  // circuit table does not start from 1 and is not sequential
  if (sequential_colors)
    make_sequential_map(color_table);
  if (debug) {
    fprintf(stderr, "sequential color table:\n");
    for (unsigned int i = 0; i < n; i++)
      fprintf(stderr, "%d ", color_table[i].first);
    fprintf(stderr, "\n");
    for (unsigned int i = 0; i < n; i++)
      fprintf(stderr, "%d ", color_table[i].second);
    fprintf(stderr, "\n");
  }
}

// point_cut is not that of opts
void ncon_coloring(Geometry &geom, const vector<faceList *> &face_list,
                   const vector<edgeList *> &edge_list,
                   const vector<poleList *> &pole, const bool point_cut_calc,
                   const int lat_mode, const ncon_opts &opts)
{
  map<int, pair<int, int>> edge_color_table;
  map<int, pair<int, int>> face_color_table;

  // for build_color_tables with build_method 2, else use what is passed
  int t = (opts.build_method == 2 && opts.d != 1) ? opts.twist / 2 : opts.twist;
  int n = (opts.build_method == 2 && opts.d != 1) ? opts.ncon_order / 2
                                                  : opts.ncon_order;

  bool hyb = opts.hybrid;

  // works with old method of assigning latitudes with angles in method 3, faces
  // are actually not a hybrid but a normal point cut of n*2 with a twist of
  // t*2-1
  if (lat_mode == 2 && opts.double_sweep) {
    n *= 2;
    t *= 2;
    if (opts.hybrid) {
      hyb = false;
      t--;
    }
  }

  // fprintf(stderr,"opts.point_cut = %s\n",opts.point_cut ? "point" : "side");
  // fprintf(stderr,"point_cut_calc = %s (used for
  // face_increment)\n",point_cut_calc ? "point" : "side");

  bool pc = (!opts.hide_indent || opts.double_sweep || opts.angle_is_side_cut)
                ? point_cut_calc
                : opts.point_cut;

  bool hi = opts.hide_indent;
  if (opts.build_method == 2 && opts.d == 1)
    hi = false;

  if (opts.build_method == 2 && hi) {
    if (opts.radius_inversion)
      pc = !pc;
  }
  // fprintf(stderr,"pc = %s\n",pc ? "point" : "side");

  // increment rules for d = 1
  int face_increment = ((is_even(n) && pc) && !opts.hybrid) ? 1 : 0;
  int edge_increment = face_increment - 1;

  // faces work with old method of apply_latitudes when d=1
  if (lat_mode == 2 && opts.double_sweep)
    face_increment = 1;

  // fprintf(stderr,"face_increment = %d\n", face_increment);
  // fprintf(stderr,"edge_increment = %d\n", edge_increment);

  // fprintf(stderr,"face color table: n = %d t = %d face_increment =
  // %d\n",n,t,face_increment);
  build_color_table(face_color_table, n, t, hyb, opts.symmetric_coloring,
                    face_increment);

  // when not hiding indent of method 2, edges of hybrid
  // are actually not a hybrid but a normal point cut of n*2 with a twist of
  // t*2-1
  if ((opts.build_method == 2 && opts.d != 1) && !hi) {
    n *= 2;
    t *= 2;
    if (opts.hybrid) {
      hyb = false;
      t--;
      // edge increment for point cut
      edge_increment = 0;
    }
  }

  // fprintf(stderr,"edge color table: n = %d t = %d edge_increment =
  // %d\n",n,t,edge_increment);
  build_color_table(edge_color_table, n, t, hyb, opts.symmetric_coloring,
                    edge_increment);

  // front and back twisting is opposite, fix when twisted half way
  if (opts.mod_twist &&
      double_eq((double)opts.ncon_order / opts.mod_twist, 2.0, opts.eps)) {
    for (auto &i : edge_color_table)
      swap(i.second.first, i.second.second);
    for (auto &i : face_color_table)
      swap(i.second.first, i.second.second);
  }

  // point_cut_calc needed below if build_method 3 angle causes side cut for
  // hybrid coloring

  // don't do this when circuits use flood fill in methods 2 or 3 as they will
  // be colored later
  if (opts.face_coloring_method &&
      !(opts.build_method > 1 && strchr("fcaS", opts.face_coloring_method)))
    ncon_face_coloring(geom, face_list, face_color_table, point_cut_calc, opts);

  // don't do this when circuits use flood fill in methods 2 or 3 as they will
  // be colored later
  if (opts.edge_coloring_method &&
      !(opts.build_method > 1 && strchr("fS", opts.edge_coloring_method)))
    ncon_edge_coloring(geom, edge_list, pole, edge_color_table, point_cut_calc,
                       opts);

  // hide indented edges of method 2
  if ((opts.build_method == 2 && opts.d != 1) && hi)
    set_indented_edges_invisible(geom, edge_list, pole);

  // mark edge circuits for methods 2 and 3 flood fill, or color by symmetry
  if ((opts.build_method > 1 && strchr("fS", opts.edge_coloring_method)) ||
      (opts.build_method == 1 && strchr("S", opts.edge_coloring_method)))
    mark_edge_circuits(geom, edge_list);
}

// inner_radius and outer_radius is calculated within
// double sweep is set in build_globe()
// radius_inversion is set
void build_globe(Geometry &geom, vector<coordList *> &coordinates,
                 vector<faceList *> &face_list, vector<edgeList *> &edge_list,
                 vector<poleList *> &pole, vector<int> &caps,
                 double &inner_radius, double &outer_radius,
                 bool &radius_inversion, bool &double_sweep,
                 const bool second_half, const ncon_opts &opts)
{
  // point cut is changed so save a copy
  bool point_cut_calc = opts.point_cut;

  // in methods 2 and 3, most faces are split in two
  vector<vector<int>> split_face_indexes;

  // in method 3, original faces as though they were not split
  vector<vector<int>> original_faces;

  vector<int> prime_meridian;
  if (opts.build_method == 3) {
    build_prime_polygon(geom, prime_meridian, coordinates, pole, opts);

    // if either pole was defined it is a point cut
    point_cut_calc = false;
    for (auto &i : pole)
      if (i->idx != -1)
        point_cut_calc = true;

    // forms with a vertex on Y axis only need 1/2 pass
    // forms with edge parallel with the Y axis only need 1/2 pass
    int polygons_total =
        (point_cut_calc ||
         angle_on_aligned_polygon(opts.angle, opts.ncon_order, opts.eps))
            ? opts.longitudes.front() / 2
            : opts.longitudes.front();
    double_sweep = (polygons_total == opts.longitudes.front()) ? true : false;

    // maximum latitudes is set
    form_angular_model(geom, prime_meridian, coordinates, face_list, edge_list,
                       pole, original_faces, split_face_indexes, polygons_total,
                       opts);
  }
  else {
    // inner_radius and outer_radius is calculated within
    // point_cut_calc changed from side cut to point cut if build method 2 and d
    // > 1
    build_prime_meridian(geom, prime_meridian, coordinates, inner_radius,
                         outer_radius, point_cut_calc, opts);

    if (opts.build_method == 2)
      radius_inversion =
          double_gt(fabs(opts.inner_radius), fabs(opts.outer_radius), opts.eps);

    // need original point cut for polygon numbers
    form_globe(geom, prime_meridian, coordinates, face_list, edge_list,
               point_cut_calc, second_half, opts);

    add_caps(geom, coordinates, face_list, pole, caps, point_cut_calc, opts);
    if (opts.build_method == 2 && opts.d != 1)
      // note that function used to reverse indented based on manual inner and
      // outer radii
      mark_indented_edges_invisible(edge_list, pole, radius_inversion, opts);
    find_split_faces_shell_model(geom, face_list, edge_list, pole,
                                 split_face_indexes, opts);
  }

  // apply latitude numbers at the end
  int lat_mode = 1;

  // old apply_latitudes was specifically for method 3
  if (opts.build_method == 3) {
    // the old apply_latitudes works for double_sweep so completes working with
    // d=1
    if (double_sweep)
      lat_mode = 2;
    // the old apply_latitudes must be used for method 3 compound coloring, new
    // method won't work
    if (opts.face_coloring_method == 'c')
      lat_mode = 2;
    // non-co-prime compounds don't work with new method. Use old method
    // (faster)
    if (gcd(opts.ncon_order, opts.d) != 1)
      lat_mode = 2;
  }

  // if flood fill is going to be used, use old method (faster)
  // but build_method can be 2
  if (opts.face_coloring_method == 'f')
    lat_mode = 2;

  bool do_faces = !(!opts.face_coloring_method ||
                    !strchr("slbfc", opts.face_coloring_method));
  bool do_edges = !(!opts.edge_coloring_method ||
                    !strchr("slbf", opts.edge_coloring_method));

  // latitudes are accessed in methods s,l,b,f,c
  // if not these then there is no need to set
  if (!do_faces && !do_edges)
    lat_mode = 0;

  if (lat_mode != 1)
    if (opts.build_method == 2)
      restore_indented_edges(edge_list, opts);

  if (lat_mode == 1) {
    // new 'sequential' apply_latitudes compatible with both methods 2 and 3
    apply_latitudes(geom, split_face_indexes, face_list, edge_list, pole, opts);
  }
  else if (lat_mode == 2) {
    if (opts.build_method == 3) {
      apply_latitudes(geom, original_faces, split_face_indexes, face_list,
                      edge_list, pole, opts);
      // old apply_latitudes needs to fix polygon numbers for compound coloring
      if ((opts.face_coloring_method == 'c') && !double_sweep)
        fix_polygon_numbers(face_list, opts);
    }
  }

  // sending opts.point_cut
  ncon_coloring(geom, face_list, edge_list, pole, point_cut_calc, lat_mode,
                opts);

  return;
}

double hybrid_twist_angle(const int ncon_order, const int d, const int tw,
                          const int build_method)
{
  int n = ncon_order;
  int t = tw;

  if (build_method == 2 && d != 1) {
    n /= 2;
    t /= 2;
  }

  // there is no twist 0 so positive twists have to be adjusted
  int twist = t;
  if (twist > 0)
    twist -= 1;

  // twist 1, 2, 3 is really twist 0.5, 1.5, 2.5 so add half twist
  double half_twist = 360.0 / (n * 2);
  double angle = half_twist + half_twist * twist * 2;

  return angle;
}

// if partial model, delete appropriate elements
void delete_unused_longitudes(Geometry &geom, vector<faceList *> &face_list,
                              vector<edgeList *> &edge_list,
                              const vector<int> &caps, const bool opposite,
                              const ncon_opts &opts)
{
  if (full_model(opts.longitudes))
    return;

  // don't allow caps to be deleted
  if (opts.build_method == 2) {
    for (int cap : caps) {
      int lon = face_list[cap]->lon;
      if (lon < opts.longitudes.front() / 2)
        face_list[cap]->lon = -1;
    }
  }

  vector<int> delete_list;
  vector<int> delete_elem;
  for (unsigned int i = 0; i < face_list.size(); i++) {
    bool val = (face_list[i]->lon >= opts.longitudes.back());
    if ((!opposite && val) || (opposite && !val)) {
      delete_list.push_back(i);
      // face_no doesn't work for method 3
      // int j = face_list[i]->face_no;
      delete_elem.push_back(i);
    }
  }

  if (delete_list.size()) {
    delete_face_list_items(face_list, delete_list);
    geom.del(FACES, delete_elem);
  }

  delete_list.clear();
  delete_elem.clear();
  for (unsigned int i = 0; i < edge_list.size(); i++) {
    bool val = (edge_list[i]->lon >= opts.longitudes.back());
    if ((!opposite && val) || (opposite && !val)) {
      delete_list.push_back(i);
      // edge_no doesn't work for method 3
      // int j = edge_list[i]->edge_no;
      delete_elem.push_back(i);
    }
  }

  if (delete_list.size()) {
    delete_edge_list_items(edge_list, delete_list);
    geom.del(EDGES, delete_elem);
  }

  geom.del(VERTS, geom.get_info().get_free_verts());

  if (opts.build_method == 2 || (opts.build_method == 1 && opts.hybrid)) {
    vector<int> single_verts;
    for (unsigned int i = 0; i < geom.verts().size(); i++) {
      vector<int> faces_with_index = find_faces_with_vertex(geom.faces(), i);
      if (faces_with_index.size() == 1)
        single_verts.push_back(i);
    }
    geom.del(VERTS, single_verts);
  }
}

// if edges are not specified they need to be cleared
void delete_unused_edges(Geometry &geom, vector<edgeList *> &edge_list,
                         const ncon_opts &opts)
{
  if (!opts.edge_coloring_method) {
    clear_edges(edge_list);
    geom.clear(EDGES);
    geom.colors(VERTS).clear();
  }
}

bool triangle_zero_area(const Geometry &geom, const int idx1, const int idx2,
                        const int idx3, const double eps)
{
  const vector<Vec3d> &verts = geom.verts();
  Vec3d xprod = vcross(verts[idx1] - verts[idx2], verts[idx1] - verts[idx3]);
  return (double_eq(xprod[0], 0.0, eps) && double_eq(xprod[1], 0.0, eps) &&
          double_eq(xprod[2], 0.0, eps));
}

void add_triangles_to_close(Geometry &geom, vector<int> &added_triangles,
                            const double eps)
{
  vector<int> face(3);
  vector<int> face_check(3);

  for (unsigned int i = 0; i < 3; i++)
    face_check[i] = 0;

  vector<vector<int>> unmatched_edges = find_unmatched_edges(geom);
  while (unmatched_edges.size()) {
    vector<bool> used(unmatched_edges.size());
    unsigned int sz = unmatched_edges.size();
    for (unsigned int i = 0; i < sz - 1; i++) {
      if (used[i])
        continue;
      face[0] = unmatched_edges[i][0];
      face[1] = unmatched_edges[i][1];
      bool found = false;
      for (unsigned int j = i + 1; j < sz; j++) {
        if (used[j] || unmatched_edges[j][0] < face[1])
          continue;
        if (unmatched_edges[j][0] > face[1])
          break;
        face[2] = unmatched_edges[j][1];
        found = false;
        if (triangle_zero_area(geom, face[0], face[1], face[2], eps)) {
          geom.add_face(face, Color(0, 0, 0, 0));
          added_triangles.push_back(geom.faces().size() - 1);
          found = true;

          used[i] = true;
          used[j] = true;

          break;
        }
      }
      if (!found) {
        if (face[0] == unmatched_edges[i + 1][0]) {
          face[2] = unmatched_edges[i + 1][1];

          if (triangle_zero_area(geom, face[0], face[1], face[2], eps)) {
            geom.add_face(face, Color(0, 0, 0, 0));
            added_triangles.push_back(geom.faces().size() - 1);
          }

          // check for infinite loop
          if (face[0] == face_check[0] && face[1] == face_check[1] &&
              face[2] == face_check[2]) {
            fprintf(stderr,
                    "warning: face %d %d %d failed. Polygon at limits "
                    "of accuracy (method=3)\n",
                    face[0], face[1], face[2]);
            return;
          }
        }
      }
    }

    face_check[0] = unmatched_edges[0][0];
    face_check[1] = unmatched_edges[0][1];
    face_check[2] = unmatched_edges[1][1];

    unmatched_edges = find_unmatched_edges(geom);
  }
}

// hybrids lose opts.longitudes.front()/2-1 from the edge_list
// be able to back it up and restore it using color elements that are carried
// along with them
void backup_flood_longitude_edges(Geometry &geom, vector<edgeList *> &edge_list,
                                  const ncon_opts &opts)
{
  int lat = 1;
  int lon = opts.longitudes.front() / 2 - 1;
  unsigned int sz = 0;
  do {
    vector<int> idx = find_edge_by_lat_lon(edge_list, lat, lon);
    sz = idx.size();
    for (unsigned int j = 0; j < sz; j++)
      geom.colors(EDGES).set(idx[j], Color(lat));
    lat++;
  } while (sz);
}

// after restore, change color index to INT_MAX for flood fill procedure
void restore_flood_longitude_edges(Geometry &geom,
                                   vector<edgeList *> &edge_list,
                                   const ncon_opts &opts)
{
  const vector<vector<int>> &edges = geom.edges();

  for (unsigned int i = 0; i < edges.size(); i++) {
    Color c = geom.colors(EDGES).get(i);
    if (!c.is_maximum_index() && !c.is_invisible()) {
      edge_list.push_back(
          new edgeList(i, c.get_index(), opts.longitudes.front() / 2 - 1));
      geom.colors(EDGES).set(i, Color::maximum_index);
    }
  }
}

void backup_flood_longitude_faces(Geometry &geom, vector<faceList *> &face_list,
                                  const ncon_opts &opts)
{
  int lat = 0;
  int lon = opts.longitudes.front() / 2 - 1;
  unsigned int sz = 0;
  do {
    vector<int> idx = find_face_by_lat_lon(face_list, lat, lon);
    sz = idx.size();
    for (unsigned int j = 0; j < sz; j++)
      geom.colors(FACES).set(idx[j], Color(lat));
    lat++;
  } while (sz);
}

// after restore, unset face color
void restore_flood_longitude_faces(Geometry &geom,
                                   vector<faceList *> &face_list,
                                   const ncon_opts &opts)
{
  const vector<vector<int>> &faces = geom.faces();

  for (unsigned int i = 0; i < faces.size(); i++) {
    if (geom.colors(FACES).get(i).is_index()) {
      int j = geom.colors(FACES).get(i).get_index();
      face_list.push_back(
          new faceList(i, j, opts.longitudes.front() / 2 - 1, 0));
      geom.colors(FACES).set(i, Color());
    }
  }
}

// opts is not const since some options are calculated
int process_hybrid(Geometry &geom, ncon_opts &opts)
{
  int ret = 0;

  // hybrids which are really point cuts with radii swapped on opposite side
  bool special_hybrids = ((opts.build_method == 2) && (opts.d == 1) &&
                          !std::isnan(opts.inner_radius));

  // retain longitudes settings
  bool full = full_model(opts.longitudes);
  int longitudes_back = opts.longitudes.back();

  // attributes of elements
  vector<coordList *> coordinates;
  vector<faceList *> face_list;
  vector<edgeList *> edge_list;

  // create memory for poles 0 - North Pole 1 - South Pole
  vector<poleList *> pole;
  pole.push_back(new poleList);
  pole.push_back(new poleList);

  // keep track of caps
  vector<int> caps;

  Geometry geom_d;

  // inner and outer radius calculated if not set
  // double sweep is set in build_globe()

  // build side cut half first
  double inner_radius_save = opts.inner_radius;
  double outer_radius_save = opts.outer_radius;
  bool point_cut_save = opts.point_cut;
  opts.point_cut =
      (opts.build_method == 3 && opts.angle_is_side_cut) ? true : false;
  opts.point_cut =
      (opts.build_method == 2 && opts.d == 1 && special_hybrids) ? true : false;
  opts.longitudes.back() = opts.longitudes.front() / 2;
  bool second_half = false;
  build_globe(geom_d, coordinates, face_list, edge_list, pole, caps,
              opts.inner_radius, opts.outer_radius, opts.radius_inversion,
              opts.double_sweep, second_half, opts);
  bool radius_inversion_save = opts.radius_inversion;

  // delete half the model
  if (opts.build_method == 3)
    delete_unused_longitudes(geom_d, face_list, edge_list, caps, false, opts);

  // backup face longitudes that are needed for flood fill
  if (opts.face_coloring_method == 'f')
    backup_flood_longitude_faces(geom_d, face_list, opts);

  // backup edge longitudes that are needed for flood fill
  if (opts.edge_coloring_method == 'f')
    backup_flood_longitude_edges(geom_d, edge_list, opts);

  // start over. build base part second, then rotate it
  clear_coord(coordinates);
  clear_faces(face_list);
  clear_edges(edge_list);
  caps.clear();

  opts.inner_radius = inner_radius_save;
  opts.outer_radius = outer_radius_save;
  opts.point_cut = !opts.point_cut;
  if (opts.build_method == 2 && opts.d == 1 && special_hybrids) {
    opts.point_cut = true;
    swap(opts.inner_radius, opts.outer_radius);
  }
  opts.longitudes.back() = opts.longitudes.front() / 2;
  second_half = true;
  build_globe(geom, coordinates, face_list, edge_list, pole, caps,
              opts.inner_radius, opts.outer_radius, opts.radius_inversion,
              opts.double_sweep, second_half, opts);
  opts.point_cut = point_cut_save;

  // delete half the model
  if (opts.build_method == 3)
    delete_unused_longitudes(geom, face_list, edge_list, caps, true, opts);
  caps.clear();

  // we can do the twist by transforming just one part.
  // negative angle because z reflection
  double twist_angle = hybrid_twist_angle(opts.ncon_order, opts.d, opts.twist,
                                          opts.build_method);
  if (opts.build_method == 2 && opts.d == 1 && special_hybrids) {
    twist_angle = (360.0 / opts.ncon_order) * opts.twist;
    // swap these back for later use
    swap(opts.inner_radius, opts.outer_radius);
    opts.radius_inversion = radius_inversion_save;
  }
  geom.transform(Trans3d::rotate(0, 0, deg2rad(-twist_angle)));
  // methods 1 and 2 need reflection
  if (opts.build_method < 3)
    geom.transform(Trans3d::reflection(Vec3d(0, 0, 1)));

  // merge the two halves and merge the vertices
  geom.append(geom_d);

  // merge by using polar orbit coordinates. this keeps the face_list pointing
  // to the right faces
  vector<polarOrb *> polar_orbit;
  find_polar_orbit(geom, polar_orbit, opts.build_method, opts.eps);
  merge_halves(geom, polar_orbit, opts.eps);
  polar_orbit.clear();

  // for build method 3 there are some unmatched edges
  // make it a valid polyhedron, and be able to flood fill across divide
  vector<int> added_triangles;
  if (opts.build_method == 3)
    add_triangles_to_close(geom, added_triangles, opts.eps);

  opts.longitudes.back() = longitudes_back;

  if (opts.face_coloring_method == 'f') {
    // restore face longitudes needed for flood fill
    restore_flood_longitude_faces(geom, face_list, opts);
    ret = ncon_face_coloring_by_adjacent_face(geom, face_list, opts);
  }
  else if (opts.face_coloring_method == 'c')
    ret = ncon_face_coloring_by_compound(geom, face_list, caps, opts);
  else if (opts.face_coloring_method == 'a')
    ncon_alternate_compound_coloring(geom, opts);

  if (opts.edge_coloring_method == 'f') {
    // restore edge longitudes needed for flood fill
    restore_flood_longitude_edges(geom, edge_list, opts);
    ncon_edge_coloring_by_adjacent_edge(geom, edge_list, pole, opts);
  }

  // if not a full model, added triangles no longer needed
  if (added_triangles.size() && !full_model(opts.longitudes))
    geom.del(FACES, added_triangles);

  // allow for partial open model in hybrids
  if (!full)
    delete_unused_longitudes(geom, face_list, edge_list, caps, false, opts);
  delete_unused_edges(geom, edge_list, opts);

  // in the case of hybrid and side cut
  // model needs to be rotated into position at side_cut angles
  if (opts.build_method == 3 && opts.angle_is_side_cut)
    geom.transform(Trans3d::rotate(0, deg2rad(180.0), 0) *
                   Trans3d::rotate(0, 0, deg2rad(twist_angle)));

  // clean up
  clear_coord(coordinates);
  clear_faces(face_list);
  clear_edges(edge_list);

  return ret;
}

// opts is not const since some options are calculated
int process_normal(Geometry &geom, ncon_opts &opts)
{
  int ret = 0;

  vector<coordList *> coordinates;
  vector<faceList *> face_list;
  vector<edgeList *> edge_list;

  // create memory for poles 0 - North Pole 1 - South Pole
  vector<poleList *> pole;
  pole.push_back(new poleList);
  pole.push_back(new poleList);

  // keep track of caps
  vector<int> caps;

  // inner and outer radius calculated if not set
  // double sweep is set in build_globe()

  bool second_half = false;
  build_globe(geom, coordinates, face_list, edge_list, pole, caps,
              opts.inner_radius, opts.outer_radius, opts.radius_inversion,
              opts.double_sweep, second_half, opts);

  // now we do the twisting
  // twist plane is now determined by points landing on z-plane
  vector<polarOrb *> polar_orbit;
  find_polar_orbit(geom, polar_orbit, opts.build_method, opts.eps);

  // method 1: can't twist when half or less of model is showing
  // method 2 and 3: whole model exists even though some longitudes will later
  // be deleted
  if ((opts.build_method > 1) ||
      (opts.build_method == 1 &&
       (2 * opts.longitudes.back() > opts.longitudes.front()))) {
    // in case of build_method 3, polygon may have been doubled. treat like 2N
    int adjust = (opts.double_sweep) ? 2 : 1;
    ncon_twist(geom, polar_orbit, coordinates, face_list, edge_list,
               opts.ncon_order * adjust, opts.twist * adjust, opts.longitudes);
  }

  // for build method 3 there are some duplicate vertices and unmatched edges
  vector<int> added_triangles;
  if (opts.build_method == 3) {
    polar_orbit.clear();
    find_polar_orbit(geom, polar_orbit, opts.build_method, opts.eps);
    merge_halves(geom, polar_orbit, opts.eps);
    polar_orbit.clear();

    // for build method 3 there are some unmatched edges
    // make it a valid polyhedron, and be able to flood fill across divide
    add_triangles_to_close(geom, added_triangles, opts.eps);
  }

  if (opts.build_method > 1 && opts.face_coloring_method == 'f')
    ret = ncon_face_coloring_by_adjacent_face(geom, face_list, opts);
  else if (opts.face_coloring_method == 'c')
    ret = ncon_face_coloring_by_compound(geom, face_list, caps, opts);
  else if (opts.face_coloring_method == 'a')
    ncon_alternate_compound_coloring(geom, opts);

  if (opts.build_method > 1 && opts.edge_coloring_method == 'f')
    ncon_edge_coloring_by_adjacent_edge(geom, edge_list, pole, opts);

  // if not a full model, added triangles no longer needed
  if (added_triangles.size() && !full_model(opts.longitudes))
    geom.del(FACES, added_triangles);

  delete_unused_longitudes(geom, face_list, edge_list, caps, false, opts);
  delete_unused_edges(geom, edge_list, opts);

  // horizontal closure only works for build method 1 and 2
  if (opts.build_method < 3)
    if (strchr(opts.closure.c_str(), 'v'))
      close_latitudinal(geom, face_list, pole, opts);

  // debug
  if (opts.info && opts.build_method == 3) {
    fprintf(stderr,
            "debug (twist=0): face circuits(d=1) = %d face "
            "circuits(counted) = %d edge circuits(counted) = %d\n",
            num_lats(opts.ncon_order, opts.point_cut),
            num_lats_faces(face_list), num_lats_edges(edge_list));
  }

  // clean up
  clear_coord(coordinates);
  clear_faces(face_list);
  clear_edges(edge_list);

  return ret;
}

void ncon_make_inv_edges_of_inv_faces(Geometry &geom)
{
  const vector<vector<int>> &faces = geom.faces();
  const vector<vector<int>> &edges = geom.edges();

  for (unsigned int i = 0; i < edges.size(); i++) {
    if ((geom.colors(EDGES).get(i)).is_invisible())
      continue;
    vector<int> face_idx = find_faces_with_edge(faces, edges[i]);
    vector<Color> cols;
    for (int j : face_idx)
      cols.push_back(geom.colors(FACES).get(j));
    Color c = ncon_average_edge_color(cols);
    if (c.is_invisible())
      set_edge_color(geom, i, Color::invisible, 255);
  }
}

void filter(Geometry &geom, const char *elems)
{
  for (const char *p = elems; *p; p++) {
    switch (*p) {
    case 'v':
      geom.colors(VERTS).clear();
      break;
    case 'e':
      geom.colors(EDGES).clear();
      break;
    case 'E':
      ncon_make_inv_edges_of_inv_faces(geom);
      break;
    case 'f':
      geom.clear(FACES);
      break;
    }
  }
}

// another coloring method by Adrian Rossiter
struct ht_less {
  static double get_eps() { return 1e-10; }
  bool operator()(const double h0, const double h1) const
  {
    return double_lt(h0, h1, get_eps());
  }
};

Geometry build_gear_polygon(const int N, const int D, const double o_radius,
                            const double i_radius, const double poly_scale,
                            const double eps)
{
  Geometry gear;

  int N2 = N;
  double arc = 180.0 / N2;
  double angle = 0.0;

  if (D > 1)
    N2 *= 2;
  else
    arc *= 2.0;

  vector<int> face;
  for (int i = 0; i < N2; i++) {
    double radius = (is_even(i)) ? o_radius : i_radius;
    gear.add_vert(
        Vec3d(cos(deg2rad(angle)) * radius, sin(deg2rad(angle)) * radius, 0));
    gear.add_edge(make_edge(i, (i + 1) % N2));
    face.push_back(i);
    angle += arc;
  }
  gear.add_face(face);

  if (double_ne(poly_scale, 1.0, eps))
    gear.transform(Trans3d::scale(poly_scale));

  return gear;
}

// transfer colors of the regular polygon to the gear polygon
void transfer_colors(Geometry &gpgon, const Geometry &pgon, const bool digons,
                     ncon_opts &opts)
{
  map<int, int> v_map;
  for (unsigned int i = 0; i < pgon.verts().size(); i++) {
    Vec3d pv = pgon.verts(i);
    for (unsigned int j = 0; j < gpgon.verts().size(); j++) {
      Vec3d gv = gpgon.verts(j);
      if (!compare(pv, gv, opts.eps)) {
        v_map[i] = j;
        break;
      }
    }
  }

  // in case polygons don't line up
  if (!v_map.size()) {
    opts.warning("transfer_colors(): no congruency found");
    return;
  }

  for (unsigned int i = 0; i < v_map.size(); i++) {
    int gp_idx = v_map[i];
    if (opts.hide_indent)
      gpgon.colors(VERTS).set(gp_idx, pgon.colors(VERTS).get(i));
    vector<int> gp_edges = find_edges_with_vertex(gpgon.edges(), gp_idx);
    for (int j : gp_edges) {
      vector<int> gp_edge = gpgon.edges(j);
      for (unsigned int k = 0; k < 2; k++) {
        if (gp_edge[k] != gp_idx) {
          Vec3d gp_vertex =
              (!digons) ? gpgon.verts(gp_edge[k]) : gpgon.edge_cent(gp_edge[k]);
          vector<int> p_edges = find_edges_with_vertex(pgon.edges(), i);
          for (unsigned int m = 0; m < 2; m++) {
            Vec3d v1 = pgon.verts(pgon.edges(p_edges[m])[0]);
            Vec3d v2 = pgon.verts(pgon.edges(p_edges[m])[1]);
            if ((point_in_segment(gp_vertex, v1, v2, opts.eps)).is_set()) {
              gpgon.colors(EDGES).set(j, pgon.colors(EDGES).get(p_edges[m]));
            }
          }
        }
      }
    }
  }

  if (opts.hide_indent && (opts.d != 1)) {
    for (unsigned int i = 1; i < gpgon.verts().size(); i += 2) {
      gpgon.colors(VERTS).set(i, Color(Color::invisible));
    }
  }
}

void pgon_post_process(Geometry &pgon, vector<Vec3d> &axes, const int N,
                       const int twist, const bool hyb, const ncon_opts &opts)
{
  int t_mult = opts.symmetric_coloring ? 1 : 2;

  int Dih_num = N / gcd(t_mult * twist - hyb, N);
  vector<vector<set<int>>> sym_equivs;
  get_equiv_elems(pgon, Symmetry(Symmetry::D, Dih_num).get_trans(),
                  &sym_equivs);

  Coloring e_clrng(&pgon);
  ColorMap *f_map = opts.face_map.clone();
  e_clrng.add_cmap(f_map);
  e_clrng.e_sets(sym_equivs[1]);

  Coloring v_clrng(&pgon);
  ColorMap *e_map = opts.edge_map.clone();
  v_clrng.add_cmap(e_map);
  v_clrng.v_sets(sym_equivs[0]);

  // pgon.transform(Trans3d::rotate(Vec3d::Z, (1-2*(D>1 &&
  // is_even(D%2)))*M_PI/2));
  pgon.transform(Trans3d::rotate(Vec3d::Z, M_PI / 2));
  axes[0] = Vec3d::Y;
  axes[1] =
      Trans3d::rotate(Vec3d::Z, -2 * M_PI * (twist - 0.5 * hyb) / N) * axes[0];
}

void lookup_face_color(Geometry &geom, const int f, const vector<Vec3d> &axes,
                       vector<map<double, Color, ht_less>> &heights,
                       const bool other_axis)
{
  int ax = double_ge(geom.face_cent(f)[2], 0.0,
                     ht_less::get_eps()); // z-coordinate determines axis
  if (other_axis)
    ax = (!ax) ? 1 : 0;
  double hts[3];
  for (int v_idx = 0; v_idx < 3; v_idx++)
    hts[v_idx] = vdot(geom.face_v(f, v_idx), axes[ax]);
  // try to select non-horizontal edge
  int offset = double_eq(hts[0], hts[1], ht_less::get_eps());
  double ht;
  if (offset && double_eq(hts[1], hts[2], ht_less::get_eps()))
    ht = hts[0] / fabs(hts[0]); // horizontal edge (on  horizontal face)
  else {
    // Find nearpoint of swept edge line, make unit, get height on axis
    Vec3d near_pt = nearest_point(Vec3d(0, 0, 0), geom.face_v(f, offset),
                                  geom.face_v(f, offset + 1));
    ht = vdot(near_pt.unit(), axes[ax]);
  }
  map<double, Color, ht_less>::iterator mi;
  if ((mi = heights[ax].find(ht)) != heights[ax].end())
    set_face_color(geom, f, mi->second);
}

void lookup_edge_color(Geometry &geom, const int e, const vector<Vec3d> &axes,
                       vector<map<double, Color, ht_less>> &heights,
                       const bool other_axis)
{
  int ax = geom.edge_cent(e)[2] >= 0.0; // z-coordinate determines axis
  if (other_axis)
    ax = (!ax) ? 1 : 0;
  double hts[2];
  for (int v_idx = 0; v_idx < 2; v_idx++)
    hts[v_idx] = vdot(geom.edge_v(e, v_idx), axes[ax]);
  // select horizontal edges that don't intersect the axis
  if (double_eq(hts[0], hts[1], ht_less::get_eps()) &&
      !lines_intersection(geom.edge_v(e, 0), geom.edge_v(e, 1), Vec3d(0, 0, 0),
                          axes[ax], ht_less::get_eps())
           .is_set()) {
    double ht = vdot(geom.edge_v(e, 0).unit(), axes[ax]);
    map<double, Color, ht_less>::iterator mi;
    if ((mi = heights[ax].find(ht)) != heights[ax].end()) {
      set_edge_color(geom, e, mi->second);
    }
  }
}

void rotate_polygon(Geometry &pgon, const int N, const bool pc, const bool hyb,
                    const ncon_opts &opts)
{
  // rotate polygons
  double rot_angle = 0;
  if ((opts.build_method == 3) && double_ne(opts.angle, 0, opts.eps)) {
    rot_angle = deg2rad(opts.angle);
    if (hyb && pc)
      rot_angle += M_PI / N;
  }
  else if (!is_even(N) || !pc || hyb)
    rot_angle = M_PI / N;

  pgon.transform(Trans3d::rotate(Vec3d::Z, rot_angle));

  // odds are flipped under these rules
  bool is_flipped = false;
  if (!is_even(N)) {
    // method 3 these angles will cause an upward flip
    if ((opts.build_method == 3) && double_ne(opts.angle, 0, opts.eps) &&
        (angle_on_aligned_polygon(opts.angle, N, opts.eps)))
      is_flipped = true;
  }
  else
      // special special hybrids and their kin when N/2 is odd
      if ((opts.build_method == 2) && !is_even(N / 2) && (opts.d == 1) &&
          !std::isnan(opts.inner_radius))
    is_flipped = true;

  // note: the polygon is still on its side
  if (is_flipped)
    pgon.transform(Trans3d::reflection(Vec3d(1, 0, 0)));
}

void color_by_symmetry(Geometry &geom, bool &radius_set, ncon_opts &opts)
{
  // note: N and D are not pre-doubled here for build_method 2
  int N = opts.ncon_order;
  int D = opts.d;
  int twist = opts.twist;
  bool hyb = opts.hybrid;

  bool pc = opts.point_cut;
  if ((opts.build_method == 2) && opts.radius_inversion)
    pc = !pc;

  // when D=1, and radii are not equal, hybrids are really point cuts
  if ((opts.build_method == 2) && (D == 1) && radius_set && hyb) {
    pc = true;
    hyb = false;
  }
  if ((opts.build_method == 3) && opts.angle_is_side_cut)
    pc = false;

  // if build_method 2, and radius set is such a way that it corresponds to
  // another D
  // then change to that D and use normal radius
  if ((opts.build_method == 2 && D > 1) && radius_set) {
    double i_radius = opts.inner_radius;
    double o_radius = opts.outer_radius;
    // inner and outer radius must be both positive or negative
    if ((double_gt(o_radius, 0, opts.eps) &&
         double_gt(i_radius, 0, opts.eps)) ||
        (double_lt(o_radius, 0, opts.eps) &&
         double_lt(i_radius, 0, opts.eps))) {
      i_radius = fabs(i_radius);
      o_radius = fabs(o_radius);
      double ratio = double_gt(o_radius, i_radius, opts.eps)
                         ? o_radius / i_radius
                         : i_radius / o_radius;
      for (int i = 1; i <= N / 2; i++) {
        i_radius = 0;
        o_radius = 0;
        double arc = 0;
        calc_radii(i_radius, o_radius, arc, N * 2, i, opts, true);
        double rat = o_radius / i_radius;
        if (double_eq(rat, ratio, ht_less::get_eps())) {
          if (D != i) {
            D = i;
            opts.warning(msg_str("radii cause D to change to %d", D), 'R');
            radius_set = false;
            if (D == 1)
              N *= 2;
          }
          break;
        }
      }
    }
  }

  Geometry pgon;
  vector<Vec3d> axes(2);
  vector<map<double, Color, ht_less>> heights(2);

  bool gear_polygon_used = false;
  bool digons = (N == 2 * D) ? true : false;

  // create polygon
  Polygon dih(N, D, Polygon::dihedron);
  dih.set_edge(1, 1.00);
  dih.make_poly(pgon);
  pgon.add_missing_impl_edges();

  rotate_polygon(pgon, N, pc, hyb, opts);

  // scale polygon
  double poly_scale = 1.0;
  double i_radius = 0;
  double o_radius = 0;
  double arc = 0;

  // method 3: reflect on Y axis when angled
  if (opts.build_method == 3) {
    // 2N/N polygons need polygon resized
    if (digons)
      pgon.transform(Trans3d::scale(2.0));

    // if it is formed by double sweeping, mirror on Y
    if (opts.double_sweep) {
      Geometry pgon_refl;
      pgon_refl = pgon;
      pgon_refl.transform(Trans3d::reflection(Vec3d(0, 1, 0)));
      pgon.append(pgon_refl);

      // this flip is needed
      if (!is_even(N))
        pgon.transform(Trans3d::reflection(Vec3d(1, 0, 0)));

      // when angle is used
      // a normal becomes a normal side cut of 2N/2D with a twist of 2T
      // a hybrid becomes a normal side cut of 2N/2D with a twist of 2T-1
      N *= 2;
      D *= 2;
      twist *= 2;
      if (hyb) {
        hyb = false;
        twist--;
      }
      pc = false;
    }
  }

  // find inner and outer radius for regular polygon
  if (opts.build_method == 2) {
    int N2 = (D == 1) ? N : N * 2;
    calc_radii(i_radius, o_radius, arc, N2, D, opts, true);

    // patch for 2N/N or D is 1
    if (digons || D == 1)
      radius_set = true;

    // regular polygon that was generated it needs to be rescaled to overlay
    // 2N/N polygons need polygon resized
    if (digons)
      pgon.transform(Trans3d::scale(o_radius / 0.5));

    poly_scale =
        (double_gt(fabs(opts.outer_radius), fabs(opts.inner_radius), opts.eps)
             ? fabs(opts.outer_radius)
             : fabs(opts.inner_radius)) /
        o_radius;
    if (double_ne(poly_scale, 1.0, opts.eps))
      pgon.transform(Trans3d::scale(poly_scale));
  }

  pgon_post_process(pgon, axes, N, twist, hyb, opts);

  // if build method 2, build gear polygon and transfer colors and replace pgon
  // and set radii if they were manually changed
  if ((opts.build_method == 2) && radius_set) {
    gear_polygon_used = true;

    // build gear polygon proper
    Geometry rpgon = build_gear_polygon(N, D, opts.outer_radius,
                                        opts.inner_radius, 1.0, opts.eps);

    rotate_polygon(rpgon, N, pc, hyb, opts);
    pgon_post_process(rpgon, axes, N, twist, hyb, opts);

    // build gear polygon as regular polygon
    Geometry gpgon =
        build_gear_polygon(N, D, o_radius, i_radius, poly_scale, opts.eps);

    rotate_polygon(gpgon, N, pc, hyb, opts);

    if (D > 1) {
      N *= 2;
      twist *= 2;
      // a hybrid becomes a normal of 2N with a twist of 2T-1
      if (hyb) {
        hyb = false;
        twist--;
      }
    }

    pgon_post_process(gpgon, axes, N, twist, hyb, opts);

    // transfer colors to the gear polygon
    transfer_colors(gpgon, pgon, digons, opts);

    // take proper gear polygons vertices
    vector<Vec3d> &verts = gpgon.raw_verts();
    for (unsigned int i = 0; i < rpgon.verts().size(); i++)
      verts[i] = rpgon.verts(i);

    pgon = gpgon;
  }

  // special case, if twisted half way, turn upside down (post coloring)
  // use original N as it may have been doubled
  // catch hybrid case if special hybrids
  if (is_even(opts.ncon_order) && (opts.mod_twist != 0) && (N / twist == 2) &&
      !opts.hybrid)
    pgon.transform(Trans3d::reflection(Vec3d(0, 1, 0)));

  // color faces
  if (opts.face_coloring_method == 'S') {
    // Find nearpoints of polygon edge lines, make unit, get height on each axis
    for (unsigned int i = 0; i < pgon.edges().size(); i++) {
      for (int ax = 0; ax < 2; ax++) {
        Vec3d near_pt =
            nearest_point(Vec3d(0, 0, 0), pgon.edge_v(i, 0), pgon.edge_v(i, 1));
        double ht = vdot(near_pt.unit(), axes[ax]);
        Color c = pgon.colors(EDGES).get(i);
        if (opts.face_opacity > -1) {
          int opq = opts.face_pattern[i % opts.face_pattern.size()] == '1'
                        ? opts.face_opacity
                        : 255;
          c = set_opacity(c, opq);
        }
        heights[ax][ht] = c;
      }
    }

    for (unsigned int f = 0; f < geom.faces().size(); f++) {
      if (geom.faces(f).size() > 2) {
        // if, because negative radii, the face center is shifted onto the wrong
        // axis
        // no color will be found for look up, try the other axis
        lookup_face_color(geom, f, axes, heights, false);
        if (geom.colors(FACES).get(f) == opts.face_default_color)
          lookup_face_color(geom, f, axes, heights, true);
      }
    }
  }

  // color edges
  if (opts.edge_coloring_method == 'S') {
    if (opts.build_method == 2 && !opts.hide_indent && !gear_polygon_used) {
      // build gear polygon proper
      pgon = build_gear_polygon(N, D, opts.outer_radius, opts.inner_radius, 1.0,
                                opts.eps);

      rotate_polygon(pgon, N, pc, hyb, opts);

      if (D > 1) {
        N *= 2;
        twist *= 2;
        // a hybrid becomes a normal of 2N with a twist of 2T-1
        if (hyb) {
          hyb = false;
          twist--;
        }
      }

      pgon_post_process(pgon, axes, N, twist, hyb, opts);

      // special case, if twisted half way, turn upside down (post coloring)
      // use original N as it may have been doubled
      // catch hybrid case if special hybrids
      if (is_even(opts.ncon_order) && (opts.mod_twist != 0) &&
          (N / twist == 2) && !opts.hybrid)
        pgon.transform(Trans3d::reflection(Vec3d(0, 1, 0)));
    }

    // Find nearpoints of polygon vertices, make unit, get height on each axis
    for (int ax = 0; ax < 2; ax++)
      heights[ax].clear();
    for (unsigned int i = 0; i < pgon.verts().size(); i++) {
      for (int ax = 0; ax < 2; ax++) {
        double ht = vdot(pgon.verts(i).unit(), axes[ax]);
        Color c = pgon.colors(VERTS).get(i);
        if (opts.edge_opacity > -1) {
          int opq = opts.edge_pattern[i % opts.edge_pattern.size()] == '1'
                        ? opts.edge_opacity
                        : 255;
          c = set_opacity(c, opq);
        }
        heights[ax][ht] = c;
      }
    }

    for (unsigned int e = 0; e < geom.edges().size(); e++) {
      // marked edges are the ones to be colored
      Color c = geom.colors(EDGES).get(e);
      if (c.is_maximum_index()) {
        lookup_edge_color(geom, e, axes, heights, false);
        c = geom.colors(EDGES).get(e);
        // if, because negative radii, the edge center is shifted onto the wrong
        // axis
        // no color will be found for look up, try the other axis
        if (c.is_maximum_index()) {
          // not found? try again
          lookup_edge_color(geom, e, axes, heights, true);
          c = geom.colors(EDGES).get(e);
          // if still not found, give up and make it invisible
          if (c.is_maximum_index())
            set_edge_color(geom, e, Color::invisible, 255);
        }
      }
    }

    // occasionally some vertices can be missed
    for (unsigned int i = 0; i < geom.verts().size(); i++) {
      Color c = geom.colors(VERTS).get(i);
      if (c.is_maximum_index())
        geom.colors(VERTS).set(i, Color::invisible);
    }

    // try to color poles. in method 2 a model can be forced upside down
    bool pc_p = pc;
    if ((opts.build_method == 2) && (!is_even(opts.ncon_order / 2) && !pc))
      pc_p = true;
    if ((opts.build_method == 3) &&
        (!is_even(opts.ncon_order) &&
         double_eq(angle_in_range(opts.angle, opts.eps), 180.0, opts.eps)))
      pc_p = true;
    if ((!twist || (N / twist == 2)) && pc_p) {
      int north = -1;
      int south = -1;
      find_poles(pgon, north, south, opts);

      // congruent points on model take those colors
      if (north != -1) {
        for (unsigned int i = 0; i < geom.verts().size(); i++) {
          if (!compare(geom.verts(i), pgon.verts(north), opts.eps)) {
            geom.colors(VERTS).set(i, pgon.colors(VERTS).get(north));
            break;
          }
        }
      }
      if (south != -1) {
        for (unsigned int i = 0; i < geom.verts().size(); i++) {
          if (!compare(geom.verts(i), pgon.verts(south), opts.eps)) {
            geom.colors(VERTS).set(i, pgon.colors(VERTS).get(south));
            break;
          }
        }
      }
    }
  }

  // for method 3 and digons, faces have trouble getting colored
  // take colors from attached edge
  if ((opts.build_method == 3) && (opts.face_coloring_method == 'S') &&
      digons) {
    if (!opts.edge_coloring_method)
      opts.warning(
          "digons in method 3 are taken from edge coloring. none was specified",
          'e');
    else {
      for (unsigned int i = 0; i < geom.faces().size(); i++) {
        vector<int> face = geom.faces(i);
        unsigned int sz = face.size();
        for (unsigned int j = 0; j < sz; j++) {
          vector<int> edge = make_edge(face[j], face[(j + 1) % sz]);
          int edge_no = find_edge_in_edge_list(geom.edges(), edge);
          Color c = geom.colors(EDGES).get(edge_no);
          if ((edge_no != -1) && !c.is_invisible()) {
            geom.colors(FACES).set(i, c);
            break;
          }
        }
      }
    }
  }

  // optionally append polygon
  if (opts.add_symmetry_polygon)
    geom.append(pgon);
}

int ncon_subsystem(Geometry &geom, ncon_opts &opts)
{
  int ret = 0;

  bool point_cut_save = opts.point_cut;

  // the best way to do side cut with method 3 is angle
  if (opts.build_method == 3 && opts.hybrid && !opts.point_cut) {
    opts.angle = 180.0 / opts.ncon_order;
    opts.angle_is_side_cut = true;
  }

  bool radius_set =
      ((opts.build_method == 2) &&
       (!std::isnan(opts.inner_radius) || !std::isnan(opts.outer_radius)))
          ? true
          : false;

  if (opts.info)
    fprintf(stderr, "========================================\n");

  if (opts.hybrid)
    ret = process_hybrid(geom, opts);
  else
    ret = process_normal(geom, opts);

  // restore for possible reporting
  if (opts.build_method == 2 && opts.d != 1) {
    opts.ncon_order /= 2;
    opts.twist /= 2;
  }

  opts.point_cut = point_cut_save;

  // if inner radius is greater than outer radius, point cut will change to a
  // side cut and vice-versa
  // note: inner and outer radii may not be completely defined until after
  // construction
  if (opts.radius_inversion) {
    // change point_cut for possible reporting in model_info()
    opts.point_cut = !opts.point_cut;
    opts.warning(msg_str("manual change in radii changed model to %s cut",
                         (opts.point_cut ? "point" : "side")));
  }

  if (opts.info) {
    model_info(geom, opts);

    vector<surfaceTable *> surface_table;
    if ((opts.d > 1 && (opts.ncon_order - opts.d > 1)) ||
        (opts.build_method == 3 && opts.angle != 0.0)) {
      // fprintf(stderr,"face/edge circuit info not available for star n_icons
      // or those with non-zero angles\n");
    }
    else {
      surfaceData sd;
      ncon_info(opts.ncon_order, opts.point_cut, opts.twist, opts.hybrid,
                opts.info, surface_table, sd);
      surface_table.clear();
    }
  }

  if (opts.face_coloring_method == 'S' || opts.edge_coloring_method == 'S')
    color_by_symmetry(geom, radius_set, opts);

  if (opts.edge_coloring_method == 'F')
    ncon_edge_coloring_from_faces(geom, opts);

  // Color post-processing
  // process the uninvolved edges of the n_icon
  if (opts.unused_edge_color.is_set())
    color_unused_edges(geom, opts.unused_edge_color);

  if (opts.edge_set_no_color)
    // now is the point that edges which are supposed to be unset can be done
    unset_marked_edges(geom);
  else
    // when partial model, some vertices can't be rotated, so out of sync with
    // their edges
    // when method 2 or 3, front vertices can get back color because of order of
    // edges in list
    // faster to just recolor all vertices
    reassert_colored_verts(geom, opts.edge_default_color,
                           opts.unused_edge_color);

  // some edges can remain uncolored
  if (opts.edge_coloring_method)
    color_uncolored_edges(geom, opts);

  // some faces can remain uncolored
  if (opts.face_coloring_method)
    color_uncolored_faces(geom, opts);

  geom.orient(1); // positive orientation

  // elements can be chosen to be eliminated completely
  filter(geom, opts.hide_elems.c_str());

  return ret;
}

void surface_subsystem(const ncon_opts &opts)
{
  vector<surfaceTable *> surface_table;
  surfaceData sd;

  char form = opts.ncon_surf[0];

  fprintf(stderr, "\n");

  if (!opts.filter_case2)
    fprintf(stderr,
            "Note: case 2 n-icons are depicted with {curly brackets}\n");

  if ((form != 'o') && (form != 'i'))
    fprintf(stderr,
            "Note: non-chiral n-icons are depicted with [square brackets]\n");

  if (form == 's')
    fprintf(
        stderr,
        "Note: all even order side cut n-icons have at least two surfaces\n");
  else if (form == 'i')
    fprintf(stderr, "Note: only hybrids such that N/2 is even are shown\n");
  else if (form == 'j')
    fprintf(stderr, "Note: only hybrids such that N/2 is odd are shown\n");
  else if (form == 'k')
    fprintf(stderr, "Note: only hybrids such that N/4 is even are shown\n");
  else if (form == 'l')
    fprintf(stderr, "Note: only hybrids such that N/4 is odd are shown\n");

  fprintf(stderr, "\n");

  vector<int> ncon_range = opts.ncon_range;

  int inc = 2;
  if (form == 'o') {
    if (is_even(ncon_range.front()))
      ncon_range.front()++;
  }
  else {
    if (!is_even(ncon_range.front()))
      ncon_range.front()++;
  }

  if ((form == 'i') || (form == 'j')) {
    if (((form == 'i') && !is_even(ncon_range.front() / 2)) ||
        ((form == 'j') && is_even(ncon_range.front() / 2)))
      ncon_range.front() += 2;

    inc += 2;
    form = 'h';
  }
  else if ((form == 'k') || (form == 'l')) {
    if (((form == 'k') && !is_even(ncon_range.front() / 4)) ||
        ((form == 'l') && is_even(ncon_range.front() / 4)))
      ncon_range.front() += 4;

    inc += 6;
    form = 'h';
  }

  if (opts.long_form) {
    fprintf(stderr,
            "                       Surfaces                       Edges\n");
    fprintf(stderr, "                       ------------------------------ "
                    "------------------------\n");
    fprintf(stderr, "%s  %s          %s %s %s %s %s\n\n", "Order", "n-icon",
            "Total", "Continuous", "Discontinuous", "Continuous",
            "Discontinuous");
  }

  int last = 0;
  for (int ncon_order = ncon_range.front(); ncon_order <= ncon_range.back();
       ncon_order += inc) {
    if (is_even(ncon_order)) {
      if (form == 'h')
        last = (int)floor((double)(ncon_order + 2) / 4);
      else
        last = (int)floor((double)ncon_order / 4);
    }
    else
      last = (int)floor((double)ncon_order / 2);

    bool point_cut = false;
    bool hybrid = false;
    bool info = false;

    if (form == 'n' || form == 'o')
      point_cut = true;
    else if (form == 's')
      point_cut = false;
    else if (form == 'h')
      hybrid = true;

    bool none = true;

    if (opts.long_form)
      fprintf(stderr, "%-5d: ", ncon_order);
    else
      fprintf(stderr, "%d: ", ncon_order);

    for (int twist = 2; twist <= last; twist++) {
      ncon_info(ncon_order, point_cut, twist, hybrid, info, surface_table, sd);

      if ((!is_even(ncon_order) && sd.total_surfaces > 1) ||
          (form != 's' && sd.total_surfaces > 1) ||
          (form == 's' && sd.total_surfaces > 2)) {
        if (!sd.ncon_case2 || (sd.ncon_case2 && !opts.filter_case2)) {
          if (!none) {
            if (opts.long_form)
              fprintf(stderr, "%-5d: ", ncon_order);
            else
              fprintf(stderr, ", ");
          }
          string buffer;
          if (sd.nonchiral)
            buffer = "[" + std::to_string(ncon_order) + "+" +
                     std::to_string(twist) + "]";
          else if (sd.ncon_case2)
            buffer = "{" + std::to_string(ncon_order) + "+" +
                     std::to_string(twist) + "}";
          else
            buffer = "(" + std::to_string(ncon_order) + "+" +
                     std::to_string(twist) + ")";
          if (opts.long_form)
            fprintf(stderr, "%-15s %5d %10d %13d %10d %13d\n", buffer.c_str(),
                    sd.total_surfaces, sd.c_surfaces, sd.d_surfaces, sd.c_edges,
                    sd.d_edges);
          else
            fprintf(stderr, "%s", buffer.c_str());
          none = false;
        }
      }
    }
    surface_table.clear();

    if (none) {
      fprintf(stderr, "none");
      if (opts.long_form)
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
  }
}

int main(int argc, char *argv[])
{
  int ret = 0;

  ncon_opts opts;
  opts.process_command_line(argc, argv);

  if (opts.ncon_surf.length())
    surface_subsystem(opts);
  else {
    Geometry geom;
    ret = ncon_subsystem(geom, opts);

    opts.write_or_error(geom, opts.ofile);
  }

  return ret;
}
