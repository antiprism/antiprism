/*
   Copyright (c) 2014-2021, Roger Kaufman, Adrian Rossiter

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
   Name: symmetro.cc
   Description: Make symmetrohedra
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using std::string;
using std::swap;
using std::vector;

using namespace anti;

ColorMapMap *alloc_default_map(const string map_file, const char m)
{
  auto *col_map = new ColorMapMap;
  int n = (m == 'n') ? 3 : 0;

  // color digons gray
  if (n == 3)
    col_map->set_col(2, Color(127, 127, 127)); // digons gray50

  if (map_file == "m1") {
    col_map->set_col(n++, Color(255, 0, 0));   // axis1 red
    col_map->set_col(n++, Color(0, 0, 255));   // axis2 blue
    col_map->set_col(n++, Color(255, 255, 0)); // axis3 yellow
    col_map->set_col(n++, Color(0, 100, 0));   // convex hull darkgreen
  }
  else if (map_file == "m2") {
    // colors from PDF document measured from screen
    col_map->set_col(n++, Color(130, 95, 34));   // 3-sided faces
    col_map->set_col(n++, Color(99, 117, 88));   // 4-sided faces
    col_map->set_col(n++, Color(84, 139, 35));   // 5-sided faces
    col_map->set_col(n++, Color(96, 109, 28));   // 6-sided faces (wrapping)
    col_map->set_col(n++, Color(128, 128, 128)); // 7-sided faces
    col_map->set_col(n++, Color(118, 97, 85));   // 8-sided faces
    // all polygons with higher sides have the same shade
    col_map->set_col(n++, Color(128, 144, 79)); // 9-sided faces and higher
    col_map->set_col(n++, Color(128, 144, 79));
    col_map->set_col(n++, Color(128, 144, 79));
    col_map->set_col(n++, Color(128, 144, 79));
    col_map->set_col(n++, Color(128, 144, 79));
    col_map->set_col(n++, Color(128, 144, 79));
    col_map->set_col(n++, Color(128, 144, 79));
    col_map->set_col(n++, Color(128, 144, 79));
    col_map->set_col(n++, Color(128, 144, 79));
    col_map->set_col(n++, Color(128, 144, 79));
    col_map->set_col(n++, Color(128, 144, 79));
    col_map->set_col(n++, Color(128, 144, 79)); // 20-sides highest needed
  }
  else if (map_file == "m3") {
    col_map->set_col(n++, Color(255, 0, 0));     // 3-sided red
    col_map->set_col(n++, Color(255, 127, 0));   // 4-sided darkoranage1
    col_map->set_col(n++, Color(255, 255, 0));   // 5-sided yellow
    col_map->set_col(n++, Color(0, 100, 0));     // 6-sided darkgreen
    col_map->set_col(n++, Color(0, 255, 255));   // 7-sided cyan
    col_map->set_col(n++, Color(0, 0, 255));     // 8-sided blue
    col_map->set_col(n++, Color(255, 0, 255));   // 9-sided magenta
    col_map->set_col(n++, Color(255, 255, 255)); // 10-sided white
    col_map->set_col(n++, Color(127, 127, 127)); // 11-sided gray50
    col_map->set_col(n++, Color(0, 0, 0));       // 12-sided black
  }

  return col_map;
}

class symmetro_opts : public ProgramOpts {
public:
  string ofile;

  char sym = '\0';          // symmetry entry will be T,O,I or D
  int p = 0;                // rotational order
  int q = 0;                // rotational order
  int dihedral_n = 0;       // dihedral number for D
  int sym_id_no = 1;        // id number for multiple symmetry choices
  vector<int> multipliers;  // l,m,n multipliers for symmetries
  vector<int> d;            // denominators for n/d
  vector<int> d_substitute; // denominators substitutes for n/d
  char sym_mirror = '\0';   // reflection
  int vert_z = std::numeric_limits<int>::max(); // for -c twisters z position
  double rotation = 0;                          // rotation of polygon
  double rotation_as_increment = 0;             // rotation in radians
  bool add_pi = false;                          // add pi to radians rotation
  int rotation_axis = -1;          // axis to apply rotation, defaults to 0
  double angle_between_axes = NAN; // override calculation
  vector<double> scale;            // control scale of polygons 0,1 or 2
  int scale_axis = -1;             // axis for which to set scale
  int convex_hull = 0;             // mode of convex hull, calculated
  string frame_elems;              // axis or sectional tiling or both
  double offset = 0;               // to avoid coplanarity if it occurs
  bool remove_free_faces = false;  // unconnected faces can be removed
  bool verbose = false;            // report math output
  char mode = '\0'; // k - Kaplan, t -twister, s - s symmetry, c - c symmetry

  vector<int> col_axis_idx;               // color method a, colors by axis
  char face_coloring_method = 'a';        // a - by axis, n - by polygon sides
  int face_opacity = -1;                  // transparency from 0 to 255
  bool color_digons = false;              // add an edge to color the digon
  Color vert_col = Color(255, 215, 0);    // gold
  Color edge_col = Color(211, 211, 211);  // lightgray
  Color face_col = Color();               // not set
  Color frame_col = Color(135, 206, 235); // skyblue3
  ColorMapMulti map;

  double eps = anti::epsilon;

  symmetro_opts() : ProgramOpts("symmetro") {}

  void process_command_line(int argc, char **argv);
  void usage();
};

void symmetro_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options]

Symmetrohedra and Twisters are created by placing equilateral polygons centered
on the symmetry axes of Icosahedral, Octahedral, Tetrahedral, or Dihedral
symmetry. The sides of the polygons will be a multiple number of the axis
reflection number. Axes are numbered as 0, 1 and 2

It is possible to generate models such that the polygons intersect. If a
collision is detected, convex hull will be suppressed

options -k, -t, -s and -c cannot be used together but one needs to be specified 

Options
%s
  -v        verbose output
  -l <lim>  minimum distance for unique vertex locations as negative exponent
               (default: %d giving %.0e)
  -o <file> write output to file (default: write to standard output)

Program Options
  -k <s,l,m,n,a> Kaplan-Hart notation. Generate Symmetrohedra based on a study
            by Craig S. Kaplan and George W. Hart (http://www.georgehart.com).
            Url: http://www.cgl.uwaterloo.ca/~csk/projects/symmetrohedra
            s: symmetry type of Symmetrohedra. sets {p,q,2}
               I-icosahedral {5,3,2} O-octahedral {4,3,2} T-tetrahedral {3,3,2}
            l,m,n: multipliers for axis polygons. Separated by commas, two
               multipliers may be * or 0, no more than two positive integers
            a: face rotation type: vertex=1, edge=0  (default: 1)
            example: -k i,2,*,4,e
  -t <s[p,q]i,m1,m2> Twister notation. Generate twister models.
            s: symmetry. I-icosahedral, O-octahedral, T-tetrahedral, D-dihedral
            p,q: rotational order of each of the two axes. may be swapped
            i: (default: 1): integer to select between non-equivalent pairs of
               axes having the same symmetry group and rotational orders
            m1,m2: an integer multiplier for each axis. i.e. m1*p and m2*q
               also can be entered as m1/d, m2/d fractional values
               e.g. T[3,2],1,2  I[5,2]2,1/2,3  D7[7,3],1,2  D11[2,2]5,2,2
                  Axis pairs are from the following
                  T: [3,3], [3,2], [2,2]
                  O: [4,4], [4,3], [4,2]x2, [3,3], [3,2]x2, [2,2]x2
                  I: [5,5], [5,3]x2, [5,2]x3, [3,3]x2, [3,2]x4, [2,2]x4
                  Dn:[n,2], [2,2]x(n/2 rounded down)
  -s <n/d:D>,s S symmetry twisters. denominator d optional (default: 1)
               optional D substitutes a polygon of n/D in place of n/d
               optional s: symmetry override c - C symmetry
  -c <n1/d1:D1,n2/d2:D2,s,v> C symmetry twisters
               optional d denominator (default: 1)
               optional D substitutes a polygon of n/D in place of n/d
               optional n2 can differ from n1. if not specified n2=n1, d2=d1
               optional s: symmetry: c - C (default), h - Ch, v - Cv, d - D
               optional v: vertex index of radial polygon to bring to z plane
                  v of -1 is original position (default: index of largest z)
  -M <opt>  mirroring (may create compound). Can be x, y or z (default: none)
  -a <a,n>  a in degrees of rotation given to polygon applied to optional
               axis n (default: 0)  radians may be entered as 'rad(a)'
               keyword: e - +1 radian, x - +1+pi radians
  -r <r,n>  set the edge length of the polygon on axis n (default: 0)
               to r. Must be non-negative. The default edge length is 1
  -A <a>    a in degrees is angle between axes (default: calculated)
  
Scene Options
  -C <mode> convex hull. polygons=1, suppress=2, force=3, auto=4 (default: 4)
  -q <args> include frame elements in output
               r - rhombic tiling edges, a - rotation axes (default: none)
  -O <dist> amount to offset the first polygon to avoid coplanarity with the
               second polygon, for example 0.0001 (default: 0.0)
  -x        remove any free faces that are produced

Coloring Options (run 'off_util -H color' for help on color formats)
  -V <col>  vertex color (default: gold)
  -E <col>  edge color   (default: lightgray)
  -f <col>  face color. Or use method for using color in map (default: a)
               keyword: none - sets no color
               a - color by axis number
               n - color by number of sides
  -T <tran> face transparency. valid range from 0 (invisible) to 255 (opaque)
  -D        don't cover digons with edge color
  -Q <col>  frame color  (default: skyblue3)
  -m <maps> color maps for faces to be tried in turn (default: m1)
               keyword m1: red,blue,yellow,darkgreen
                  note: position 4 color is for faces added by convex hull
               keyword m2: approximating colors in the symmetrohedra pdf file
               keyword m3: red,darkorange1,yellow,darkgreen,cyan,blue,magenta,
                           white,gray50,black (default when -f n)
               (built in maps start at index 3 when -f n)

)",
          prog_name(), help_ver_text, int(-log(anti::epsilon) / log(10) + 0.5),
          anti::epsilon);
}

// from StackOverflow
// remove blanks from input char*
char *deblank(char *input)
{
  int i, j;
  char *output = input;

  for (i = 0, j = 0; i < (int)strlen(input); i++, j++) {
    if (input[i] != ' ')
      output[j] = input[i];
    else
      j--;
  }
  output[j] = '\0';

  return output;
}

void symmetro_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  string id;
  string map_file;
  vector<int> n;

  handle_long_opts(argc, argv);

  // initialize d
  for (int i = 0; i < 2; i++)
    d.push_back(1);

  // initialize d_substitute
  for (int i = 0; i < 2; i++)
    d_substitute.push_back(0);

  // initialize scale array
  for (int i = 0; i < 3; i++)
    scale.push_back(NAN);

  while ((c = getopt(argc, argv,
                     ":hk:t:m:s:c:M:a:r:A:C:q:O:xvf:Q:V:E:DT:l:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    // Kaplan-Hart notation
    case 'k': {
      if (mode)
        error("-k, -t, -s, -c cannot be used together", c);
      mode = 'k';

      char *opts = deblank(optarg);

      Split pts(opts, ",");
      vector<string> parts;
      for (unsigned int i = 0; i < pts.size(); i++)
        parts.push_back(pts[i]);

      unsigned int parts_sz = parts.size();
      if (parts_sz < 4 || parts_sz > 5)
        error("expecting 4 or 5 parameters for Kaplan-Hart notation", c);
      else if (parts_sz == 4) {
        parts.push_back("v");
        parts_sz++;
      }

      string mult;
      int num_multipliers = 0;
      for (unsigned int i = 0; i < parts_sz; i++) {
        if (i == 0) {
          sym = toupper(parts[i][0]);
          if (!strchr("TOI", sym))
            error(msg_str("invalid symmetry character '%c'", sym), c);
        }
        else if (i == 1 || i == 2) {
          string part = parts[i];
          if (part.length() == 1 && part[0] == '*')
            part = "0";
          mult += part + ",";
        }
        else if (i == 3) {
          string part = parts[i];
          if (part.length() == 1 && part[0] == '*')
            part = "0";
          mult += part;

          print_status_or_exit(
              read_int_list(mult.c_str(), multipliers, true, 3), c);

          // might not be able to happen
          if ((int)multipliers.size() != 3)
            error("3 multipliers must be specified", c);

          for (int multiplier : multipliers) {
            if (multiplier > 0)
              num_multipliers++;
          }
          if (num_multipliers == 0)
            error("at least one axis multiplier must be specified", c);
          else if (num_multipliers == 3)
            error("at least one axis multiplier must be * or zero", c);

          if (multipliers[2] == 1)
            warning("model will contain digons", c);

          int orders[3] = {0, 3, 2};
          orders[0] = (sym == 'T') ? 3 : ((sym == 'O') ? 4 : 5);

          if (num_multipliers == 1) {
            if (multipliers[0]) {
              p = orders[0];
              q = orders[0];
              if (sym == 'T')
                rotation += 120.0;
            }
            else if (multipliers[1]) {
              p = orders[1];
              q = orders[1];
              if (sym == 'T')
                rotation += 120.0;
            }
            else if (multipliers[2]) {
              p = orders[2];
              q = orders[2];
            }
          }
          else if (num_multipliers == 2) {
            if (multipliers[0] && multipliers[1]) {
              p = orders[0];
              q = orders[1];
            }
            else if (multipliers[0] && multipliers[2]) {
              p = orders[0];
              q = orders[2];
            }
            else if (multipliers[1] && multipliers[2]) {
              p = orders[1];
              q = orders[2];
            }
          }
        }
        else if (i == 4) {
          // in the paper, edge connection is shown as 'e', vertex connection is
          // shown a '1'
          print_status_or_exit(get_arg_id(parts[i].c_str(), &id,
                                          "edge=0|vertex=1",
                                          argmatch_add_id_maps),
                               c);
          if (id == "0")
            rotation_as_increment = rad2deg(1.0);
          else if (id == "1")
            rotation_as_increment = rad2deg(0.0);
        }
      }

      // for octahedral and icosahedral, axis2 alone
      if (num_multipliers == 1 && multipliers[2]) {
        if (id == "1") { // vertex connected
          // rotate to coincident faces
          if (sym == 'T') {
            rotation += 45.0; // 45.0 degrees
          }
          else if (sym == 'O') {
            rotation += rad2deg(acos(1.0 / 3.0) /
                                2.0); // 35.26438968275465431577 degrees
          }
          else if (sym == 'I') {
            rotation += rad2deg(acos(2.0 / sqrt(5.0)) /
                                2.0); // 13.28252558853899467604 degrees
            if (!is_even(multipliers[2]))
              rotation += 90.0 / (multipliers[2] * 2.0);
          }
        }
      }

      break;
    }

    // twister notation
    case 't': {
      if (mode)
        error("-k, -t, -s, -c cannot be used together", c);
      mode = 't';

      char *opts = deblank(optarg);
      bool id_num_given = (strstr(opts, "],")) ? false : true;

      Split pts(opts, ",[]");
      vector<string> parts;
      for (unsigned int i = 0; i < pts.size(); i++)
        parts.push_back(pts[i]);

      if (!id_num_given) {
        vector<string>::iterator it;
        it = parts.begin();
        parts.insert(it + 3, "1");
      }

      unsigned int parts_sz = parts.size();
      if (parts_sz != 6)
        error("incorrect format for Twister notation", c);

      for (unsigned int i = 0; i < parts_sz; i++) {
        if (i == 0) {
          sym = toupper(parts[i][0]);
          if (!strchr("TOID", sym))
            error(msg_str("invalid symmetry character '%c'", sym), c);

          // dihedral
          if (sym == 'D') {
            if (parts[i].length() < 2)
              error("No N found after D symmetry specifier", c);
            print_status_or_exit(read_int(parts[i].c_str() + 1, &dihedral_n),
                                 "option t: dihedral symmetry N");
          }
        }
        else if (i == 1) {
          print_status_or_exit(read_int(parts[i].c_str(), &p),
                               "option t: axis 1");
          if (p < 2)
            error("axis 1 rotational order number be greater than 1", c);
        }
        else if (i == 2) {
          print_status_or_exit(read_int(parts[i].c_str(), &q),
                               "option t: axis 2");
          if (q < 2)
            error("axis 2 rotational order number be greater than 1", c);
        }
        else if (i == 3) {
          print_status_or_exit(read_int(parts[i].c_str(), &sym_id_no),
                               "option t: symmetry id number");
          if (sym_id_no <= 0)
            error("symmetry id number must be positive", c);
        }
        else if (i == 4 || i == 5) {
          if (!strchr(parts[i].c_str(), '/')) {
            int mult;
            print_status_or_exit(read_int(parts[i].c_str(), &mult),
                                 "option t: multiplier");
            if (mult <= 0)
              error("multiplier must be positive", c);
            multipliers.push_back(mult);
          }
          else {
            int n_part;
            int d_part = 1;

            Split parts2(parts[i], "/");
            unsigned int parts2_sz = parts2.size();

            for (unsigned int j = 0; j < parts2_sz; j++) {
              if (j == 0) {
                print_status_or_exit(read_int(parts2[j], &n_part),
                                     "option t: n/d (n part)");

                if (n_part <= 0)
                  error("n of n/d must be positive", c);
              }
              else if (j == 1) {
                print_status_or_exit(read_int(parts2[j], &d_part),
                                     "option t: n/d (d part)");

                if (d_part <= 0)
                  error("d of n/d must be positive", c);
              }
            }

            multipliers.push_back(n_part);
            d[(i == 4) ? 0 : 1] = d_part;
          }
        }
      }

      if (sym == 'D') {
        bool reversed = (p < q) ? true : false;
        if ((!reversed && (p != dihedral_n && p != 2)) ||
            (reversed && (q != dihedral_n && q != 2)))
          error(msg_str("when symmetry is D, axis %d rotational order must "
                        "equal 2 or N (%d)",
                        (reversed ? 2 : 1), dihedral_n),
                c);
        if ((!reversed && q != 2) || (reversed && p != 2))
          error(msg_str(
                    "when symmetry is D, axis %d rotational order must equal 2",
                    (reversed ? 2 : 1)),
                c);
      }

      break;
    }

    // S symmetry
    case 's': {
      if (mode)
        error("-k, -t, -s, -c cannot be used together", c);
      mode = 's';

      string sym_override;

      Split parts(optarg, ",");
      unsigned int parts_sz = parts.size();

      for (unsigned int i = 0; i < parts_sz; i++) {
        if (i == 0) {
          int n_part;
          int d_part;
          int d_sub;

          bool found_slash = (strchr(parts[i], '/')) ? true : false;
          bool found_colon = (strchr(parts[i], ':')) ? true : false;

          Split parts2(parts[i], "/:");
          unsigned int parts2_sz = parts2.size();

          for (unsigned int j = 0; j < parts2_sz; j++) {
            if (j == 0) {
              print_status_or_exit(read_int(parts2[j], &n_part),
                                   "option s: n/d (n part)");

              if (n_part < 2)
                error("n must be greater than 1", c);

              n.push_back(n_part);
            }
            else if (j == 1) {
              print_status_or_exit(read_int(parts2[j], &d_part),
                                   "option s: n/d (d part)");

              if (d_part <= 0 && found_slash)
                error("d of n/d must be positive", c);

              d[i] = d_part;
            }
            else if (j == 2) {
              print_status_or_exit(read_int(parts2[j], &d_sub),
                                   "option s: substitute D");

              if (d_sub < 1 || d_sub >= n[i])
                error(
                    msg_str("substitute D must be between 1 and %d", n[i] - 1),
                    c);

              d_substitute[i] = d_sub;
            }
          }

          // if no slash, but colon exists d_substitute will be in the wrong
          // place
          if (!found_slash && found_colon) {
            d_substitute[i] = d[i];
            d[i] = 1;
            if (d_substitute[i] < 1 || d_substitute[i] >= n[i])
              error(msg_str("substitute D must be between 1 and %d", n[i] - 1),
                    c);
          }
        }
        else if (i == 1) {
          if ((strspn(parts[i], "cC") != strlen(parts[i])) ||
              strlen(parts[i]) > 1)
            error(msg_str("symmetry override is '%s' must be c", parts[i]), c);
          sym_override = parts[i];
        }
        else if (i > 1)
          error("too many arguments", c);
      }

      // fill both n/d
      if ((int)n.size() == 1) {
        n.push_back(n[0]);
        d[1] = d[0];
        d_substitute[1] = d_substitute[0];
      }

      if ((double)n[0] / (double)d[0] < 1.5)
        error("polygon: the polygon fraction cannot be less than 3/2 (base "
              "rhombic tiling is not constructible)",
              c);

      if (sym_override.length())
        sym = toupper(sym_override[0]);
      else
        sym = 'S';

      if ((is_even(n[0]) || is_even(d[0])) && !sym_override.length())
        warning("when n or d is even, model will only connect correctly at "
                "certain twist angles. try option -c",
                c);

      p = n[0];
      q = n[1];

      break;
    }

    case 'c': {
      if (mode)
        error("-k, -t, -s, -c cannot be used together", c);
      mode = 'c';

      string sym_override;

      Split parts(optarg, ",");
      unsigned int parts_sz = parts.size();

      for (unsigned int i = 0; i < parts_sz; i++) {
        // if second n/d is not specified and symmetry given move straight to
        // next term
        if ((i == 1) && (strspn(parts[i], "cCvVhHdD") == strlen(parts[i])))
          i++;

        if (i == 0 || i == 1) {
          int n_part;
          int d_part;
          int d_sub;

          bool found_slash = (strchr(parts[i], '/')) ? true : false;
          bool found_colon = (strchr(parts[i], ':')) ? true : false;

          Split parts2(parts[i], "/:");
          unsigned int parts2_sz = parts2.size();

          for (unsigned int j = 0; j < parts2_sz; j++) {
            if (j == 0) {
              print_status_or_exit(
                  read_int(parts2[j], &n_part),
                  msg_str("option c: n/d (n part) (term %d)", i + 1));

              if (n_part < 2)
                error(
                    msg_str("n of n/d must be greater than 1 (term %d)", i + 1),
                    c);
              n.push_back(n_part);
            }
            else if (j == 1) {
              print_status_or_exit(
                  read_int(parts2[j], &d_part),
                  msg_str("option c: n/d (d part) (term %d)", i + 1));

              if (d_part <= 0 && found_slash)
                error(msg_str("d of n/d must be positive (term %d)", i + 1), c);
              d[i] = d_part;
            }
            else if (j == 2) {
              print_status_or_exit(
                  read_int(parts2[j], &d_sub),
                  msg_str("option c: substitute D (term %d)", i + 1));

              if (d_sub < 1 || d_sub >= n[i])
                error(msg_str("substitute D must be between 1 and %d (term %d)",
                              n[i] - 1, i + 1),
                      c);

              d_substitute[i] = d_sub;
            }
          }

          // if no slash, but colon exists d_substitute will be in the wrong
          // place
          if (!found_slash && found_colon) {
            d_substitute[i] = d[i];
            d[i] = 1;
            if (d_substitute[i] < 1 || d_substitute[i] >= n[i])
              error(msg_str("substitute D must be between 1 and %d", n[i] - 1),
                    c);
          }
        }
        else if (i == 2) {
          if (!read_int(parts[i], &vert_z)) {
            if ((strspn(parts[i], "cCvVhHdD") != strlen(parts[i])) ||
                strlen(parts[i]) > 1)
              error(msg_str("symmetry override is '%s' must be c, v, h or d",
                            parts[i]),
                    c);
            sym_override = parts[i];
          }
        }
        else if (i == 3) {
          print_status_or_exit(read_int(parts[i], &vert_z), "option c: vert_z");
        }
        else if (i > 3)
          error("too many arguments", c);
      }

      // fill both n/d
      if ((int)n.size() == 1) {
        n.push_back(n[0]);
        d[1] = d[0];
        d_substitute[1] = d_substitute[0];
      }

      if ((double)n[0] / (double)d[0] < 1.5)
        error("the polygon fraction cannot be less than 3/2 (base rhombic "
              "tiling is not constructible)",
              c);

      if (!sym_override.length())
        sym_override = "C";
      sym = toupper(sym_override[0]);

      if (vert_z != std::numeric_limits<int>::max())
        if (vert_z < -1 || vert_z >= n[1])
          error(msg_str("vert z must be between 0 and %d", n[1] - 1), c);

      p = n[0];
      q = n[1];

      // patch for fill_sym_vec
      if (sym == 'D')
        dihedral_n = p;

      break;
    }

    case 'M':
      if (strspn(optarg, "xyz") != strlen(optarg) || strlen(optarg) > 1)
        error(msg_str("symmetry mirror %s must be x, y or z", optarg), c);
      sym_mirror = optarg[0];
      break;

    // rotation
    case 'a': {
      Split parts(optarg, ",");
      unsigned int parts_sz = parts.size();

      for (unsigned int i = 0; i < parts_sz; i++) {
        if (i == 0) {
          // see if it is built in amount
          char ex = optarg[strlen(parts[i]) - 1];
          if (ex == 'e' || ex == 'x') {
            optarg[strlen(parts[i]) - 1] = '\0';
            double num_part = 0;
            if (strlen(optarg) == 0)
              num_part = 1.0;
            else
              print_status_or_exit(read_double(parts[i], &num_part),
                                   "option a: rotation value");

            rotation_as_increment += rad2deg(num_part);
            if (ex == 'x')
              add_pi = true;
          }
          else {
            // find 'rad' in parts[i], else value is degrees
            bool rotation_as_inc = (strstr(parts[i], "rad")) ? false : true;
            double rot;
            print_status_or_exit(read_double(parts[i], &rot),
                                 "option a: rotation value");
            if (rotation_as_inc)
              rotation_as_increment += rot;
            else
              rotation += rot;
          }
        }
        else if (i == 1) {
          double rotation_axis_tmp;
          print_status_or_exit(read_double(parts[i], &rotation_axis_tmp),
                               "option a: rotation axis");
          // if ( rotation_axis_tmp - a > 0.0 )
          //   error(msg_str("axis numbers must be specified by an integer:
          //   '%g'", rotation_axis_tmp), c);
          rotation_axis = (int)floor(rotation_axis_tmp);

          if (rotation_axis < 0 || rotation_axis > 2)
            error("axis to apply rotation should be 0, 1 or 2", c);
        }
      }

      break;
    }

    // scale
    case 'r': {
      vector<double> scale_tmp;
      print_status_or_exit(read_double_list(optarg, scale_tmp, 2));

      if (scale_tmp.size() > 2)
        error("scale takes 1 or 2 arguments", c);
      else if (scale_tmp.size() == 1)
        scale[0] = scale_tmp[0];
      else if (scale_tmp.size() == 2) {
        scale_axis = (int)floor(scale_tmp[1]);
        // if ( scale_tmp[1] - scale_axis > 0.0 )
        //   error(msg_str("axis numbers must be specified by an integer: '%g'",
        //   scale_tmp[i]), c);

        if (scale_axis < 0 || scale_axis > 2)
          error("axis to apply scale should be 0, 1 or 2", c);
        else
          scale[scale_axis] = scale_tmp[0];
      }
      break;
    }

    case 'A':
      print_status_or_exit(read_double(optarg, &angle_between_axes), c);
      break;

    case 'C':
      print_status_or_exit(get_arg_id(optarg, &id,
                                      "polygons=1|suppress=2|force=3|normal=4",
                                      argmatch_add_id_maps),
                           c);
      convex_hull = atoi(id.c_str());
      break;

    case 'q':
      if (strspn(optarg, "ra") != strlen(optarg))
        error(msg_str("frame elements are '%s' must be from r and a", optarg),
              c);
      frame_elems = optarg;
      break;

    case 'O':
      print_status_or_exit(read_double(optarg, &offset), c);
      break;

    case 'x':
      remove_free_faces = true;
      break;

    case 'v':
      verbose = true;
      break;

    case 'f':
      if (strspn(optarg, "an") && strlen(optarg) == 1)
        face_coloring_method = *optarg;
      else {
        print_status_or_exit(face_col.read(optarg), c);
        face_coloring_method = '\0';
      }
      break;

    case 'Q':
      print_status_or_exit(frame_col.read(optarg));
      break;

    case 'V':
      print_status_or_exit(vert_col.read(optarg));
      break;

    case 'E':
      print_status_or_exit(edge_col.read(optarg));
      break;

    case 'D':
      color_digons = true;
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

  if (!mode)
    error("one of -k, -t, -s, -c must be specified");

  // convert n to multipliers for modes s or c
  if (n.size()) {
    multipliers.clear();
    for (unsigned int i = 0; i < n.size(); i++) {
      multipliers.push_back(n[i]);
      // mode=s or c
      if (mode == 's')
        multipliers[i] /= p;
      else if (mode == 'c') {
        if (i == 0)
          multipliers[i] /= p;
        else if (i == 1)
          multipliers[i] /= q;
      }
    }
  }

  if (!map_file.size())
    map_file = (face_coloring_method == 'n') ? "m3" : "m1";

  if (map_file != "m1" && map_file != "m2" && map_file != "m3")
    print_status_or_exit(map.init(map_file.c_str()), 'm');
  else {
    map.add_cmap(alloc_default_map(map_file, face_coloring_method));
    if (map_file == "m1" || map_file == "m2")
      map.set_wrap();
    else {
      ColorMap *spread_map = colormap_from_name("spread");
      map.add_cmap(spread_map);
    }
  }
}

class symmetro {
public:
  symmetro()
  {
    for (int i = 0; i < 2; i++) {
      mult.push_back(0);
      sym_vec.push_back(Vec3d());
      d.push_back(1);
      d_substitute.push_back(0);
      scale.push_back(1);
    }
  }

  void debug(const char mode);

  void setSym(const char s, const int psym, const int qsym, const int dih_n,
              const int id_no);
  void setMult(const int a, const int m);
  void setScale(const int a, const double s);
  void setD(const int a, const int dee);
  void setD_substitute(const int a, const int dee);

  int getOrder(const int a);
  int getN(const int a);

  double axis_angle(const int n, const int d);
  // double getAngleBetweenAxes( const int axis1, const int axis2 );
  double getAngleBetweenAxesSin(const int axis1, const int axis2);
  void swap_vecs(Vec3d &a, Vec3d &b);
  int fill_sym_vec(const char mode, string *error_msg);

  double angle(const int n, const int d);
  double circumradius(const int n, const int d);

  void substitute_polygon(Geometry &geom, const int axis_no);
  vector<Geometry> calc_polygons(const char mode, const double rotation,
                                 const double rotation_as_increment,
                                 const bool add_pi, const bool swap_axes,
                                 const double offset, const bool verbose,
                                 double &angle_between_axes, string *error_msg);

  ~symmetro() = default;

private:
  char sym;
  int p;
  int q;

  int dihedral_n;
  int sym_id_no;

  vector<int> mult;
  vector<Vec3d> sym_vec;

  vector<int> d;
  vector<int> d_substitute;

  vector<double> scale;
};

void symmetro::debug(const char mode)
{
  fprintf(stderr, "\n");

  if (strchr("kt", mode)) {
    fprintf(stderr, "symmetry = %c%s[%d,%d]%d\n", sym,
            (sym == 'D' ? std::to_string(dihedral_n).c_str() : ""), p, q,
            sym_id_no);
    fprintf(stderr, "\n");
  }

  fprintf(stderr, "vector axis 0: %.17lf %.17lf %.17lf\n", sym_vec[0][0],
          sym_vec[0][1], sym_vec[0][2]);
  fprintf(stderr, "vector axis 1: %.17lf %.17lf %.17lf\n", sym_vec[1][0],
          sym_vec[1][1], sym_vec[1][2]);
  fprintf(stderr, "\n");

  for (int i = 0; i < 2; i++)
    fprintf(stderr, "axis %d: mult = %d  scale = %.17lf\n", i, mult[i],
            scale[i]);
  fprintf(stderr, "\n");

  for (int i = 0; i < 2; i++) {
    string buffer = "";
    if (d_substitute[i] && (d_substitute[i] != d[i]))
      buffer = "(then substituted with " + std::to_string(getN(i)) + "/" +
               std::to_string(d_substitute[i]) + "-gon)";
    if (mult[i])
      fprintf(stderr, "axis %d polygon: %d/%d-gon %s\n", i, getN(i), d[i],
              buffer.c_str());
  }
  fprintf(stderr, "\n");
}

void symmetro::setSym(const char s, const int psym, const int qsym,
                      const int dih_n, const int id_no)
{
  sym = s;

  p = psym;
  q = qsym;

  dihedral_n = dih_n;
  sym_id_no = id_no;
}

void symmetro::setMult(const int a, const int m) { mult[a] = m; }

void symmetro::setScale(const int a, const double s) { scale[a] = s; }

void symmetro::setD(const int a, const int dee) { d[a] = dee; }

void symmetro::setD_substitute(const int a, const int dee)
{
  d_substitute[a] = dee;
}

int symmetro::getOrder(const int a)
{
  switch (a) {
  case 0:
    return p;
  case 1:
    return q;
  default:
    return 0;
  }
}

int symmetro::getN(const int a) { return (getOrder(a) * mult[a]); }

double symmetro::axis_angle(const int n, const int d)
{
  double nn = double(n);
  double dd = double(d);
  return (acos(1.0 / tan(M_PI * dd / nn) / tan(M_PI * (nn - dd) / (2.0 * nn))));
}

// double symmetro::getAngleBetweenAxes( const int axis1, const int axis2 ) {
//   return ( acos(vdot(sym_vec[axis1].unit(), sym_vec[axis2].unit())) );
//}

double symmetro::getAngleBetweenAxesSin(const int axis1, const int axis2)
{
  double sin_angle_between_axes =
      vcross(sym_vec[axis1].unit(), sym_vec[axis2].unit()).len();
  if (fabs(sin_angle_between_axes) > 1.0) {
    sin_angle_between_axes = (sin_angle_between_axes < 0.0) ? -1.0 : 1.0;
  }
  // return ( asin( safe_for_trig(sin_angle_between_axes) ) );
  return (asin(sin_angle_between_axes));
}

void symmetro::swap_vecs(Vec3d &a, Vec3d &b)
{
  a = -a;
  b = -b;
  swap(a, b);
}

int symmetro::fill_sym_vec(const char mode, string *error_msg = nullptr)
{
  int err_no = 0; // 1 - wrong p,q  2 - wrong sym_id_no  3 - wrong sym

  if (sym == 'T') {
    if (p == 3 && q == 3) { // K-H mode +120 degrees
      if (sym_id_no == 1) {
        sym_vec[0] = Vec3d(1, 1, 1);
        sym_vec[1] = Vec3d(-1, -1, 1);
      }
    }
    else if ((p == 3 && q == 2) || (p == 2 && q == 3)) {
      if (sym_id_no == 1) {
        sym_vec[0] = Vec3d(1, 1, 1);
        sym_vec[1] = Vec3d(0, 0, 1);
      }
    }
    else if (p == 2 && q == 2) {
      if (sym_id_no == 1) {
        sym_vec[0] = Vec3d(0, 0, 1);
        sym_vec[1] = Vec3d(1, 0, 0);
      }
    }
    else
      err_no = 1;

    if (p < q)
      swap_vecs(sym_vec[0], sym_vec[1]);
  }
  else if (sym == 'O') {
    if (p == 4 && q == 4) {
      if (sym_id_no == 1) {
        sym_vec[0] = Vec3d(0, 0, 1);
        sym_vec[1] = Vec3d(1, 0, 0);
      }
    }
    else if ((p == 4 && q == 3) || (p == 3 && q == 4)) {
      if (sym_id_no == 1) {
        sym_vec[0] = Vec3d(0, 0, 1);
        sym_vec[1] = Vec3d(1, 1, 1);
      }
    }
    else if ((p == 4 && q == 2) || (p == 2 && q == 4)) {
      if (sym_id_no == 1) {
        sym_vec[0] = Vec3d(0, 0, 1);
        sym_vec[1] = Vec3d(0, 1, 1);
      }
      else if (sym_id_no == 2) {
        sym_vec[0] = Vec3d(0, 0, 1);
        sym_vec[1] = Vec3d(1, 1, 0);
      }
    }
    else if (p == 3 && q == 3) {
      if (sym_id_no == 1) {
        sym_vec[0] = Vec3d(1, 1, 1);
        sym_vec[1] = Vec3d(1, -1, 1);
      }
    }
    else if ((p == 3 && q == 2) || (p == 2 && q == 3)) {
      if (sym_id_no == 1) {
        sym_vec[0] = Vec3d(1, 1, 1);
        sym_vec[1] = Vec3d(0, -1, -1);
      }
      else if (sym_id_no == 2) {
        sym_vec[0] = Vec3d(1, 1, 1);
        sym_vec[1] = Vec3d(1, 0, -1);
      }
    }
    else if (p == 2 && q == 2) {
      if (sym_id_no == 1) {
        sym_vec[0] = Vec3d(0, 1, 1);
        sym_vec[1] = Vec3d(1, 0, 1);
      }
      else if (sym_id_no == 2) {
        sym_vec[0] = Vec3d(0, 1, 1);
        sym_vec[1] = Vec3d(0, 1, -1);
      }
    }
    else
      err_no = 1;

    if (p < q)
      swap_vecs(sym_vec[0], sym_vec[1]);
  }
  else if (sym == 'I') {
    if (p == 5 && q == 5) {
      if (sym_id_no == 1) {
        sym_vec[0] = Vec3d(0, 1, phi);
        sym_vec[1] = Vec3d(0, 1, -phi);
      }
    }
    else if ((p == 5 && q == 3) || (p == 3 && q == 5)) {
      if (sym_id_no == 1) {
        sym_vec[0] = Vec3d(0, 1, phi);
        sym_vec[1] = Vec3d(1, 1, 1);
      }
      else if (sym_id_no == 2) {
        sym_vec[0] = Vec3d(0, 1, phi);
        sym_vec[1] = Vec3d(phi, -1 / phi, 0);
      }
    }
    else if ((p == 5 && q == 2) || (p == 2 && q == 5)) {
      if (sym_id_no == 1) {
        sym_vec[0] = Vec3d(0, 1, phi);
        sym_vec[1] = Vec3d(0, 0, -1);
      }
      else if (sym_id_no == 2) {
        sym_vec[0] = Vec3d(0, 1, phi);
        sym_vec[1] = Vec3d(1, 1 / phi, -phi);
      }
      else if (sym_id_no == 3) {
        sym_vec[0] = Vec3d(0, 1, phi);
        sym_vec[1] = Vec3d(1, 0, 0);
      }
    }
    else if (p == 3 && q == 3) {
      if (sym_id_no == 1) {
        sym_vec[0] = Vec3d(1, 1, 1);
        sym_vec[1] = Vec3d(-1 / phi, 0, -phi);
      }
      else if (sym_id_no == 2) {
        sym_vec[0] = Vec3d(1, 1, 1);
        sym_vec[1] = Vec3d(1, -1, -1);
      }
    }
    else if ((p == 3 && q == 2) || (p == 2 && q == 3)) {
      if (sym_id_no == 1) {
        sym_vec[0] = Vec3d(1, 1, 1);
        sym_vec[1] = Vec3d(-1, -1 / phi, -phi);
      }
      else if (sym_id_no == 2) {
        sym_vec[0] = Vec3d(1, 1, 1);
        sym_vec[1] = Vec3d(-1, 0, 0);
      }
      else if (sym_id_no == 3) {
        sym_vec[0] = Vec3d(1, 1, 1);
        sym_vec[1] = Vec3d(1, -1 / phi, -phi);
      }
      else if (sym_id_no == 4) {
        sym_vec[0] = Vec3d(1, 1, 1);
        sym_vec[1] = Vec3d(1, 1 / phi, -phi);
      }
    }
    else if (p == 2 && q == 2) {
      if (sym_id_no == 1) {
        sym_vec[0] = Vec3d(0, 0, 1);
        sym_vec[1] = Vec3d(1, 1 / phi, phi);
      }
      else if (sym_id_no == 2) {
        sym_vec[0] = Vec3d(0, 0, 1);
        sym_vec[1] = Vec3d(1 / phi, phi, 1);
      }
      else if (sym_id_no == 3) {
        sym_vec[0] = Vec3d(0, 0, 1);
        sym_vec[1] = Vec3d(phi, 1, 1 / phi);
      }
      else if (sym_id_no == 4) {
        sym_vec[0] = Vec3d(0, 0, 1);
        sym_vec[1] = Vec3d(1, 0, 0);
      }
    }
    else
      err_no = 1;

    if (p < q)
      swap_vecs(sym_vec[0], sym_vec[1]);
  }
  else if (sym == 'D' && mode == 't') {
    int p_tmp = (p < q) ? q : p;

    if (p_tmp == 2 && (sym_id_no <= (dihedral_n / 2))) {
      double a = sym_id_no * M_PI / dihedral_n;
      sym_vec[0] = Vec3d(1, 0, 0);
      sym_vec[1] = Vec3d(cos(a), sin(a), 0);
    }
    else if ((p_tmp == dihedral_n) && (sym_id_no == 1)) {
      sym_vec[0] = Vec3d(0, 0, 1);
      sym_vec[1] = Vec3d(1, 0, 0);
    }
    else
      err_no = 2;

    if (p < q)
      swap_vecs(sym_vec[0], sym_vec[1]);
  }
  else
    // D symmetry of option -c included
    if (strchr("SCHVD", sym)) {
      double a = axis_angle(p, d[0]);
      sym_vec[0] = Vec3d(0, 0, 1);
      sym_vec[1] = Vec3d(sin(a), 0, cos(a));
    }
    else
      err_no = 3;

  // sym_vec will only not be set if no id_no was found
  if (!err_no)
    err_no = (sym_vec[0].is_set()) ? 0 : 2;

  if (err_no == 1) {
    if (error_msg)
      *error_msg = msg_str("invalid p,q values: %d,%d", p, q);
  }
  else if (err_no == 2) {
    if (error_msg)
      *error_msg = msg_str("invalid symmetry id no: %d", sym_id_no);
  }
  else if (err_no == 3) {
    if (error_msg)
      *error_msg = msg_str("invalid symmetry: %c", sym);
  }

  return err_no;
}

double symmetro::angle(const int n, const int d)
{
  return ((2.0 * M_PI * double(d) / double(n)));
}

double symmetro::circumradius(const int n, const int d)
{
  double edge_len = 1.0;
  return (edge_len / (2.0 * sin(angle(n, d) / 2.0)));
}

void symmetro::substitute_polygon(Geometry &geom, const int axis_no)
{
  // Make one convex regular polygon
  geom.set_hull("");

  // if star polygon, rethread face
  if (d_substitute[axis_no] > 1) {
    // there will only be one face after convex hull
    vector<int> face = geom.faces(0);
    int n = face.size();
    int d = d_substitute[axis_no];

    // handle compound polygons
    int d_test = (d <= n / 2) ? d : n - d;
    int parts = (!(n % d_test)) ? d_test : 1;

    vector<vector<int>> faces_new;

    int stop = (parts == 1) ? n : n / d_test;
    int k = 0;
    for (int i = 0; i < parts; i++) {
      vector<int> face_new;
      int k2 = k;
      for (int j = 0; j < stop; j++) {
        face_new.push_back(face[k2]);
        k2 = (k2 + d) % n;
      }
      faces_new.push_back(face_new);
      k++;
    }

    // replace with new n/d face(s)
    geom.clear(FACES);
    for (auto &i : faces_new)
      geom.add_face(i);
  }
}

// angle_between_axes in radians, is modified
vector<Geometry> symmetro::calc_polygons(
    const char mode, const double rotation, const double rotation_as_increment,
    const bool add_pi, const bool swap_axes, const double offset,
    const bool verbose, double &angle_between_axes, string *error_msg = nullptr)
{
  // there will be two polygons generated in seperate geoms
  vector<Geometry> pgeom(2);

  // the two axes can be swapped if the angle is applied to the second axis
  vector<int> axis(2);
  axis[0] = 0;
  axis[1] = 1;
  if (swap_axes)
    swap(axis[0], axis[1]);

  double r0 = scale[axis[0]] * circumradius(getN(axis[0]), d[axis[0]]);
  double r1 = scale[axis[1]] * circumradius(getN(axis[1]), d[axis[1]]);

  if (!std::isnan(angle_between_axes))
    angle_between_axes = deg2rad(angle_between_axes);
  else
    angle_between_axes = (mode == 's')
                             ? axis_angle(getN(axis[0]), d[axis[0]])
                             : getAngleBetweenAxesSin(axis[0], axis[1]);
  if (verbose)
    fprintf(stderr, "\nangle between axes: radians = %.17lf degrees = %.17lf\n",
            angle_between_axes, rad2deg(angle_between_axes));

  Trans3d rot = Trans3d::rotate(Vec3d(0, 1, 0), angle_between_axes);
  Trans3d rot_inv = Trans3d::rotate(Vec3d(0, 1, 0), -angle_between_axes);

  double ang = deg2rad(rotation);
  if (rotation_as_increment)
    ang +=
        deg2rad(rotation_as_increment) * angle(getN(axis[0]), d[axis[0]]) / 2.0;
  if (add_pi)
    ang += M_PI;
  if (verbose)
    fprintf(stderr,
            "turn angle: radians = %.17lf degrees = %.17lf on axis %d\n", ang,
            rad2deg(ang), axis[0]);

  Vec3d V = Trans3d::rotate(Vec3d(0, 0, 1), ang) * Vec3d(r0, 0, 0);
  Vec3d q = rot * V;
  Vec3d u = rot * Vec3d(0, 0, 1);

  double a = u[0] * u[0] + u[1] * u[1];
  double b = 2 * (q[0] * u[0] + q[1] * u[1]);
  double c = q[0] * q[0] + q[1] * q[1] - r1 * r1;

  double disc = b * b - 4 * a * c;
  if (disc < -anti::epsilon) {
    if (error_msg)
      *error_msg = "model is not geometrically constructible";

    return pgeom;
  }
  else if (disc < 0)
    disc = 0;

  double sign_flag = -1.0;
  // modes s and c
  if (mode == 's') {
    // AR - The sign flag, which changes for the range 90 to 270 degrees, allows
    // the model to reverse, otherwise the model breaks apart in this range.
    double turn_angle_test_val = fabs(fmod(fabs(ang), 2.0 * M_PI) - M_PI);
    sign_flag = (turn_angle_test_val > M_PI / 2.0) ? -1.0 : 1.0;
  }
  double t = (-b + sign_flag * sqrt(disc)) / (2 * a);

  Vec3d P = V + Vec3d(0, 0, t);
  Vec3d Q = rot * P;

  if (vdot(sym_vec[axis[0]], sym_vec[axis[1]]) > 0.0) {
    sym_vec[axis[1]] *= -1.0;
  }

  for (unsigned int i = 0; i < pgeom.size(); i++) {
    int j = axis[i];
    int n = getN(j);

    // handle compound polygons
    int d_test = (d[j] <= n / 2) ? d[j] : n - d[j];
    int parts = (!(n % d_test)) ? d_test : 1;
    double bump_ang = angle(n, d[j]) / (double)parts;

    if ((n > 0) && scale[j]) {
      double bump_angle = 0.0;
      int vert_idx = 0;

      for (int k = 0; k < parts; k++) {
        for (int idx = 0; idx < n; idx++) {
          if (i == 0) {
            pgeom[j].add_vert(
                Trans3d::rotate(Vec3d(0, 0, 1),
                                (idx * angle(n, d[j])) + bump_angle) *
                    P +
                Vec3d(0, 0, offset));
          }
          else if (i == 1) {
            pgeom[j].add_vert(
                rot_inv *
                Trans3d::rotate(Vec3d(0, 0, 1),
                                (idx * angle(n, d[j])) + bump_angle) *
                Q);
          }
        }

        vector<int> face;
        for (int idx = 0; idx < n; ++idx)
          face.push_back(vert_idx++);
        pgeom[j].add_face(face);
        face.clear();

        bump_angle += bump_ang;
      }

      if (d_substitute[j])
        substitute_polygon(pgeom[j], j);

      pgeom[j].transform(Trans3d::align(Vec3d(0, 0, 1), Vec3d(1, 0, 0),
                                        sym_vec[axis[0]], sym_vec[axis[1]]));
      // this isn't strictly needed. turns model to line up with
      // twister_rhomb.py
      // if ( mode == 's' )
      //   pgeom[j].transform(
      //   Trans3d::rotate(Vec3d(0, 0, deg2rad(180.0/(double)p))) );
    }

    if (!scale[j])
      pgeom[j].clear_all();
  }

  return pgeom;
}

bool is_point_on_polygon_edges(const Geometry &polygon, const Vec3d &P,
                               const double eps)
{
  const vector<int> &face = polygon.faces()[0];
  const vector<Vec3d> &verts = polygon.verts();

  bool answer = false;

  unsigned int fsz = face.size();
  for (unsigned int i = 0; i < fsz; i++) {
    Vec3d v1 = verts[face[i]];
    Vec3d v2 = verts[face[(i + 1) % fsz]];
    if ((point_in_segment(P, v1, v2, eps)).is_set()) {
      answer = true;
      break;
    }
  }

  return answer;
}

bool detect_collision(const Geometry &geom, const symmetro_opts &opts)
{
  const vector<vector<int>> &faces = geom.faces();
  const vector<Vec3d> &verts = geom.verts();

  // don't test digons
  for (unsigned int i = 0; i < faces.size(); i++) {
    vector<int> face0 = faces[i];
    // digons won't work in plane intersection
    if (face0.size() < 3)
      continue;
    for (unsigned int j = i + 1; j < faces.size(); j++) {
      vector<int> face1 = faces[j];
      if (face1.size() < 3)
        continue;

      Vec3d P, dir;
      if (two_plane_intersect(centroid(verts, face0), face_norm(verts, face0),
                              centroid(verts, face1), face_norm(verts, face1),
                              P, dir, opts.eps)) {
        if (!P.is_set())
          continue;
        // if two polygons intersect, see if intersection point is inside
        // polygon
        vector<int> face_idxs;
        face_idxs.push_back(i);
        Geometry polygon = faces_to_geom(geom, face_idxs);
        // get winding number, if not zero, point is on a polygon
        int winding_number = get_winding_number(polygon, P, opts.eps);
        // if point in on an edge set winding number back to zero
        if (winding_number) {
          if (is_point_on_polygon_edges(polygon, P, opts.eps))
            winding_number = 0;
        }
        if (winding_number) {
          return true;
        }
      }
    }
  }

  return false;
}

void delete_free_faces(Geometry &geom)
{
  const vector<vector<int>> &faces = geom.faces();
  unsigned int fsz = faces.size();
  vector<bool> found(fsz);

  for (unsigned int i = 0; i < fsz; i++) {
    if (found[i])
      continue;
    // need to check for faces with lower index than i
    for (unsigned int j = 0; j < fsz; j++) {
      if (i == j)
        continue;
      for (unsigned int k = 0; k < faces[i].size(); k++) {
        if (vertex_exists_in_face(faces[j], faces[i][k])) {
          found[i] = true;
          found[j] = true;
          break;
        }
      }
      if (found[i])
        break;
    }
  }

  vector<int> face_list;
  for (unsigned int i = 0; i < fsz; i++)
    if (!found[i])
      face_list.push_back(i);

  if (face_list.size())
    geom.del(FACES, face_list);
}

void sym_repeat(Geometry &geom, const symmetro_opts &opts)
{
  Symmetry sym;
  if (opts.sym == 'T')
    sym.init(Symmetry::T);
  else if (opts.sym == 'O')
    sym.init(Symmetry::O);
  else if (opts.sym == 'I')
    sym.init(Symmetry::I);
  else if (opts.sym == 'D')
    sym.init(Symmetry::D, (opts.mode == 't') ? opts.dihedral_n : opts.p);
  else if (opts.sym == 'S')
    sym.init(Symmetry::S, opts.p * 2);
  else if (opts.sym == 'C')
    sym.init(Symmetry::C, opts.p);
  else if (opts.sym == 'V')
    sym.init(Symmetry::Cv, opts.p);
  else if (opts.sym == 'H')
    sym.init(Symmetry::Ch, opts.p);

  sym_repeat(geom, geom, sym, ELEM_FACES);

  // reflection
  if (opts.sym_mirror) {
    double x = (opts.sym_mirror == 'x') ? 1 : 0;
    double y = (opts.sym_mirror == 'y') ? 1 : 0;
    double z = (opts.sym_mirror == 'z') ? 1 : 0;

    Geometry geom_refl;
    geom_refl = geom;
    geom_refl.transform(Trans3d::reflection(Vec3d(x, y, z)));
    geom.append(geom_refl);
  }
}

Geometry build_geom(vector<Geometry> &pgeom, const symmetro_opts &opts)
{
  Geometry geom;

  // if option c, align to a z coordinate
  if (opts.mode == 'c') {
    if (!pgeom[1].verts().size()) {
      if (opts.verbose)
        fprintf(
            stderr,
            "option -c: radial polygon not found so model not moved on Z\n");
    }
    else if (opts.vert_z == -1) {
      if (opts.verbose)
        fprintf(stderr,
                "option -c: model not moved on Z (supressed by user)\n");
    }
    else {
      if (opts.verbose)
        fprintf(stderr, "option -c: greatest Z ");

      Vec3d P;
      int vz = opts.vert_z;
      if (vz != std::numeric_limits<int>::max()) {
        P = pgeom[1].verts(vz);

        if (opts.verbose)
          fprintf(stderr, "explicitly set to %d\n", vz);
      }
      else {
        // find greatest Z
        vz = 0;
        P = pgeom[1].verts(vz);
        for (unsigned int i = 1; i < pgeom[1].verts().size(); i++) {
          if (pgeom[1].verts(i)[2] > P[2]) {
            vz = i;
            P = pgeom[1].verts(vz);
          }
        }

        if (opts.verbose)
          fprintf(stderr, "calculated to %d\n", vz);
      }

      for (int i = 0; i < 2; i++) {
        Vec3d P2 = Vec3d(0, 0, -P[2]);
        pgeom[i].transform(Trans3d::translate(P2));
        pgeom[i].transform(Trans3d::rotate(
            0, 0, angle_around_axis(P, Vec3d(1, 0, 0), Vec3d(0, 0, 1))));
      }
    }
  }

  for (int i = 0; i < 2; i++) {
    // if not polygon, repeat for symmetry type
    if (opts.convex_hull > 1)
      sym_repeat(pgeom[i], opts);

    if (opts.face_coloring_method == 'a') {
      Coloring clrng(&pgeom[i]);
      Color col = opts.map.get_col(opts.col_axis_idx[i]);
      clrng.f_one_col(col);
    }

    geom.append(pgeom[i]);
  }

  if (opts.convex_hull > 1)
    merge_coincident_elements(geom, "vf", opts.eps);

  // check for collisions
  bool collision = false;
  if (opts.convex_hull == 4)
    collision = detect_collision(geom, opts);
  if (collision) {
    opts.warning("collision detected. convex hull is suppressed", 'C');
  }

  if ((!collision && opts.convex_hull == 4) || (opts.convex_hull == 3)) {
    Status stat = geom.add_hull("A1");
    // probably never happen
    if (!stat.is_ok())
      if (opts.verbose)
        opts.warning(stat.msg(), 'C');

    // merged faces will retain RGB color
    merge_coincident_elements(geom, "f", opts.eps);

    // after sort merge, only new faces from convex hull will be uncolored
    if (opts.face_coloring_method == 'a') {
      Coloring clrng(&geom);
      for (unsigned int i = 0; i < geom.faces().size(); i++) {
        Color col = geom.colors(FACES).get(i);
        if (!col.is_set()) {
          // convex hull color is map position 3
          col = opts.map.get_col(3);
          geom.colors(FACES).set(i, col);
        }
      }
    }
  }

  if (opts.face_col.is_set())
    Coloring(&geom).f_one_col(opts.face_col);
  else if (opts.face_coloring_method == 'n') {
    Coloring clrng(&geom);
    clrng.add_cmap(opts.map.clone());
    clrng.f_sides(true);
    clrng.f_apply_cmap();
  }

  // apply transparency
  Status stat = Coloring(&geom).apply_transparency(opts.face_opacity);
  if (stat.is_warning())
    opts.warning(stat.msg(), 'T');

  if (opts.remove_free_faces) {
    delete_free_faces(geom);
    geom.del(VERTS, geom.get_info().get_free_verts());
  }

  // if coloring edges, check for digons
  if (opts.color_digons) {
    vector<int> del_faces;
    for (unsigned int i = 0; i < geom.faces().size(); i++) {
      if (geom.faces(i).size() == 2) {
        vector<int> edge = make_edge(geom.faces(i)[0], geom.faces(i)[1]);
        Color col = geom.colors(FACES).get(i);
        geom.add_edge(edge, col);
        del_faces.push_back(i);
      }
    }
    geom.del(FACES, del_faces);
  }

  // color added edges
  geom.add_missing_impl_edges(opts.edge_col);

  // color vertices
  Coloring(&geom).v_one_col(opts.vert_col);

  geom.orient(1); // positive orientation

  return geom;
}

Geometry build_frame(vector<Geometry> &pgeom, const symmetro_opts &opts)
{
  Geometry geom;

  double frame_rad = pgeom[0].verts(0).len();
  Vec3d v0 = centroid(pgeom[0].verts(), pgeom[0].faces(0)).unit() * frame_rad;
  Vec3d v1 = centroid(pgeom[1].verts(), pgeom[1].faces(0)).unit() * frame_rad;

  double num_segs = 10;

  Vec3d ax = vcross(v0, v1).unit();
  double ang = angle_around_axis(v0, v1, ax);
  Trans3d mat = Trans3d::rotate(ax, -ang / num_segs);

  if (strchr(opts.frame_elems.c_str(), 'r')) {
    geom.add_vert(v0);
    for (int i = 0; i < num_segs; i++) {
      geom.add_vert(geom.verts(i) * mat);
      geom.add_edge(make_edge(i, i + 1));
    }

    if (strchr("SCHV", opts.sym)) {
      // Vec3d v2 = Vec3d(v1[0],v1[1],-v1[2]) *
      // Trans3d::rotate(0,0,deg2rad(180.0/opts.p));
      Vec3d v2 = Vec3d(v1[0], v1[1], -v1[2]) *
                 Trans3d::rotate(0, 0,
                                 (M_PI * double(opts.d[0]) / double(opts.p) *
                                  (is_even(opts.p) ? 2.0 : 1.0)));

      ax = vcross(v1, v2).unit();
      ang = angle_around_axis(v1, v2, ax);
      mat = Trans3d::rotate(ax, -ang / num_segs);

      for (int i = num_segs; i < num_segs * 2; i++) {
        geom.add_vert(geom.verts(i) * mat);
        geom.add_edge(make_edge(i, i + 1));
      }
    }
  }

  if (strchr(opts.frame_elems.c_str(), 'a')) {
    unsigned int sz = geom.verts().size();
    geom.add_vert(v0);
    geom.add_vert(-v0);
    geom.add_edge(make_edge(sz, sz + 1));
    geom.add_vert(v1);
    geom.add_vert(-v1);
    geom.add_edge(make_edge(sz + 2, sz + 3));
  }

  // if not polygon, repeat for symmetry type
  if (opts.convex_hull > 1)
    sym_repeat(geom, opts);

  // sort_merge_elems(geom, "ve", opts.eps);

  Coloring(&geom).vef_one_col(opts.frame_col, opts.frame_col, Color());

  return geom;
}

/*
void unitize_edges( vector<Geometry> &pgeom )
{
   vector<double> val(2);
   for( unsigned int i=0; i<pgeom.size(); i++ ) {
      GeometryInfo info(pgeom[i]);
      if (info.num_iedges() > 0)
         val[i] = info.iedge_lengths().sum/info.num_iedges();
   }

   double val_avg = (val[0]+val[1])/2.0;
   pgeom[0].transform(Trans3d::scale(1/val_avg));
   pgeom[1].transform(Trans3d::scale(1/val_avg));
}
*/

int main(int argc, char *argv[])
{
  symmetro_opts opts;
  opts.process_command_line(argc, argv);

  symmetro s;
  s.setSym(opts.sym, opts.p, opts.q, opts.dihedral_n, opts.sym_id_no);

  // indexes will be 0,1 except in the case of mode=k when index 2 can be
  // present
  vector<unsigned int> idx;
  for (unsigned int i = 0; i < opts.multipliers.size(); i++) {
    if (opts.multipliers[i])
      idx.push_back(i);
  }
  // in mode=k there can happen only 1 index
  if ((int)idx.size() == 1)
    idx.push_back(idx[0]);

  // set multipliers, axis index for color, in object
  for (unsigned int i = 0; i < idx.size(); i++) {
    s.setMult(i, opts.multipliers[idx[i]]);
    opts.col_axis_idx.push_back(idx[i]);
  }

  // set d, d_substitute in object
  for (unsigned int i = 0; i < opts.d.size(); i++)
    s.setD(i, opts.d[i]);
  for (unsigned int i = 0; i < opts.d_substitute.size(); i++)
    s.setD_substitute(i, opts.d_substitute[i]);

  // scale will be DBL_MAX when not set
  // if axis not specified and a scale is set, set default scale axis to first
  // axis index
  if (opts.scale_axis == -1 && (!std::isnan(opts.scale[0]))) {
    opts.scale_axis = idx[0];
    swap(opts.scale[0], opts.scale[opts.scale_axis]);
  }
  for (unsigned int i = 0; i < opts.scale.size(); i++) {
    if ((!std::isnan(opts.scale[i])) && (i != idx[0] && i != idx[1]))
      opts.error(
          msg_str("polygon '%d' is not generated so cannot be used for scaling",
                  i),
          'r');
  }
  for (unsigned int i = 0; i < opts.scale.size(); i++) {
    if (opts.scale[i] < 0)
      opts.error("scale cannot be negative", 'r');
  }
  if (opts.mode == 's' || opts.mode == 'c') {
    for (unsigned int i = 0; i < opts.scale.size(); i++) {
      if (!std::isnan(opts.scale[i])) {
        opts.warning(
            "some polygons may not meet when scale is used with -s or -c", 'r');
        break;
      }
    }
  }
  for (unsigned int i = 0; i < idx.size(); i++) {
    double r = (std::isnan(opts.scale[idx[i]])) ? 1 : opts.scale[idx[i]];
    s.setScale(i, r);
  }

  // if not specified, set default rotation axis to first axis index
  if (opts.rotation_axis == -1)
    opts.rotation_axis = idx[0];
  // check rotation axis specifier for zero
  if (opts.mode != 'k' && opts.rotation_axis == 2)
    opts.error("only 0 and 1 are valid for axis when not using option -k", 'a');
  else if (opts.rotation_axis != (int)idx[0] &&
           opts.rotation_axis != (int)idx[1])
    opts.error(
        msg_str("polygon '%d' is not generated so cannot be used for rotation",
                opts.rotation_axis),
        'a');
  // swap axes if we are rotating the axis 1 polygon
  // if -k and p == q then alway use axis 0
  bool swap_axes = false;
  if (opts.mode == 'k' && opts.p == opts.q)
    swap_axes = false;
  else if (opts.rotation_axis == (int)idx[1])
    swap_axes = true;

  // if convex_hull is not set
  if (!opts.convex_hull) {
    for (unsigned int i = 0; i < opts.d.size(); i++) {
      if (opts.d[i] > 1 || opts.d_substitute[i] > 1) {
        // supress convex hull
        opts.convex_hull = 2;
        opts.warning("star polygons detected so convex hull is supressed", 'C');
        break;
      }
    }

    // supress convex hull for models with digons
    if ((opts.multipliers[idx[0]] * opts.p == 2) ||
        (opts.multipliers[idx[1]] * opts.q == 2)) {
      opts.convex_hull = 2;
      opts.warning("model contains digons so convex hull is supressed", 'C');
    }
  }
  // if still not set, convex hull is set to normal
  if (!opts.convex_hull)
    opts.convex_hull = 4;

  // ready to generate

  // fill symmetry axes here
  string error_msg;
  if (s.fill_sym_vec(opts.mode, &error_msg))
    if (!error_msg.empty())
      opts.error(error_msg);

  vector<Geometry> pgeom =
      s.calc_polygons(opts.mode, opts.rotation, opts.rotation_as_increment,
                      opts.add_pi, swap_axes, opts.offset, opts.verbose,
                      opts.angle_between_axes, &error_msg);
  if (!error_msg.empty())
    opts.error(error_msg);

  // if ( opts.d_substitute[0] || opts.d_substitute[1] )
  //   unitize_edges( pgeom );

  if (opts.verbose) {
    s.debug(opts.mode);

    double edge_length[2] = {0};
    for (unsigned int i = 0; i < pgeom.size(); i++) {
      GeometryInfo info(pgeom[i]);
      if (info.num_iedges() > 0) {
        edge_length[i] = info.iedge_length_lims().sum / info.num_iedges();
        fprintf(stderr, "Edge length of polygon %d = %.17lf\n", i,
                edge_length[i]);
      }
    }

    fprintf(stderr, "\n");
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        if (i == j)
          continue;
        // built in epsilon here
        if (edge_length[i] > anti::epsilon && edge_length[j] > anti::epsilon)
          fprintf(stderr, "edge length ratio of polygon %d to %d = %.17lf\n", i,
                  j, edge_length[i] / edge_length[j]);
      }
    }

    fprintf(stderr, "\n");
  }

  Geometry geom;
  geom = build_geom(pgeom, opts);

  if (opts.frame_elems.length()) {
    Geometry geom_frame;
    geom_frame = build_frame(pgeom, opts);
    geom.append(geom_frame);
  }

  opts.write_or_error(geom, opts.ofile);

  return 0;
}
