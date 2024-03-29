/*
   Copyright (c) 2008-2022, Roger Kaufman

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
   Name: bravais.cc
   Description: Generate the 14 Bravais Lattices
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"
#include "color_common.h"
#include "lat_util_common.h"

#include <string>
#include <vector>

using std::string;
using std::swap;
using std::vector;

using namespace anti;

struct bravaisItem {
  const char *crystal_system;
  const char *centering;
};

// clang-format off
bravaisItem bravais_item_list[] = {
   {"Triclinic",  "P"},
   {"Monoclinic", "P"},
   {"Monoclinic", "C"},
   {"Orthorhombic", "P"},
   {"Orthorhombic", "C"},
   {"Orthorhombic", "F"},
   {"Orthorhombic", "I"},
   {"Tetragonal", "P"},
   {"Tetragonal", "I"},
   {"Trigonal", "P"},
   {"Hexagonal", "P"},
   {"Cubic", "P"},
   {"Cubic", "F"},
   {"Cubic", "I"},
};
// clang-format on

class bravais {
private:
  int last_bravais;
  bravaisItem *bravais_items;

public:
  bravais();
  void list_bravais(int idx, FILE *fp = stderr);
  void list_bravais(FILE *fp = stderr);
  int lookup_sym_no(string crystal_system, string centering);
  int get_last_bravais() { return last_bravais; }

  string get_crystal_system(int);
  string get_centering(int);
};

class brav_opts : public ProgramOpts {
public:
  string ofile;
  string cfile;

  vector<double> vecs;       // 3 vectors for edge lengths
  vector<double> angles;     // 3 angles for vectors
  vector<int> grid;          // 3 integers for grid sizes
  vector<int> prim_vec_idxs; // 4 vector vertex indexes for dual

  vector<double> strut_len;            // strut lengths to create if found
  bool cell_struts = false;            // show cell struts
  bool use_centering_for_dual = false; // calculate dual on centering type
  double radius = 0;                   // radius of container
  char radius_default = 's';           // insphere radius is default
  Vec3d radius_by_coord;               // calculate radius by x,y,z coordinate
  Vec3d offset;                        // offset from origin for container
  char container = 'c';                // default container (c all cubic space)
  bool append_container = false;       // append cage of -k container
  bool voronoi_cells = false;          // calculate voronoi cells
  bool voronoi_central_cell = false;   // include voronoi cells only at center
  char auto_grid_type = '\0';          // p,f,i,e or 8 for automatic grid size
  bool grid_for_radius = false;        // for automatic grid size
  bool convex_hull = false;            // convex hull for waterman polyhedra
  bool add_hull = false;               // add lattice to waterman polyhedra
  bool append_lattice = false;         // append lattice to final produc
  bool trans_to_origin = false;        // tranlate lattice centroid to origin
  int r_lattice_type = 0;              // hexagonal cube relation
  bool inversion = false;              // inverts F or I lattices
  bool verbose = false;                // option W gives lattice information
  bool dual_lattice = false;           // calculate dual lattice

  bool list_bravais = false;              // list radii
  Vec3d list_radii_center;                // centroid for calculating radii list
  char list_radii_original_center = '\0'; // c - use centroid o - use q offset
  int list_radii = 0;                     // list radii
  int list_struts = 0;                    // list struts

  double eps = anti::epsilon;

  char color_method = '\0'; // color method for color by symmetry
  int face_opacity = -1;    // transparency from 0 to 255
  Color cent_col;           // color a centroid, add if missing

  // colors for different parts are held in vector
  // 0 - lattice  1 - convex hull  2 - voronoi  3 - hex overlay
  vector<Color> vert_col;
  vector<Color> edge_col;
  vector<Color> face_col;

  string crystal_system; // see help
  string centering;      // see help

  brav_opts() : ProgramOpts("bravais") {}

  void process_command_line(int argc, char **argv);
  void usage();
};

void extended_help()
{
  fprintf(stdout, R"(
Definition: (partly from http://en.wikipedia.org/wiki/Bravais_lattice)
In geometry and crystallography, a Bravais lattice, named after Auguste
Bravais, is an infinite set of points generated by a set of discrete
translation operations. A crystal is made up of one or more atoms (the basis)
which is repeated at each lattice point. The crystal then looks the same when
viewed from any of the lattice points. In all, there are 14 possible Bravais
lattices that fill three-dimensional space.

August Bravais (1811-1863), a French naval officer, adventurer, and physicist
taught a course in applied mathematics for astronomy in the faculty of sciences
in Lyon from 1840. He served as the Chair of Physics, Ecole Polytechnique
between 1845 and 1856. He is best remembered for pointing out in 1845, that
there are 14 unique Bravais lattices in three dimensional crystalline systems,
adjusting the previously result (15 lattices) by Moritz Ludwig Frankenheim
obtained three years before.

A German Crystallographer, Frankenheim (1801-1869) is noted as the first to
enumerate the 32 crystal classes. And he also solved the symmetry systems of
the 7 crystal systems but this work went completely unnoticed at the time.
There is a bit of mystery surrounding what Frankenheim had as the 15th lattice.
Even today, in some texts the Hexagonal lattice with two interior points is
shown in the Trigonal class. But these two lattices use the same set of points
and it is thought that it was this duplication that was eliminated by Bravais.
However, in Bravais' paper, there is no mention of Frankenheim or the
enumeration of lattices he presented.

In this program, the Hexagonal cells and Trigonal cells can be seen together
by using the -R parameter.

Note that End Centered Cubic (would be Cubic C) does not exist but can be
produced by Tetragonal P that has cells of dimensions a,b,c = 1,1,sqrt(2)

Face Centered Cubic (Cubic F or FCC) is duplicated in Body Centered Tetragonal
(Tetragonal I) of dimensions a,b,c = 1,1,sqrt(2). However, the FCC embodied
would be of higher symmetry than the Tetragonal crystal system is allowed.

Similarly, Trigonal at 90 degrees (improper) is SC. Trigonal at 60 degrees is
FCC and Trigonal at acos(-1/3) or 109.47122063449... degrees is BCC.

Also there is no provision for Face Centered Tetragonal (would be Tetragonal F)
or Base Centered Tetragonal (would be Tetragonal C). These would be embodied in
Body Centered Tetragonal (Tetragonal I) and Simple Tetragonal (Tetragonal P)
respectively. This is true at any proportion other than a,b,c = 1,1,sqrt(2)

In Hexagonal, Orthorhombic C can be seen to occur. When Hexagonal of a=b=c is
produced, then Base Centered Tetragonal (would be Tetragonal C) occurs.

Hexagonal is sensitive to which unequal vector corresponds to the non-90 degree
angle. These must be in the same position or a Monoclinic lattice is produced.

Bravais lattices will fall into the following symmetries

Crystal System   Possible Symmetries (32 possible - note no 5 fold symmetries)
Triclinic        C1 Ci
Monoclinic       C2 Cs C2h
Orthorhombic     D2 C2v D2h
Tetragonal       C4 S4 C4h D4 C4v D2d D4h
Trigonal         C3 S6 D3 C3v D3d
Hexagonal        C6 C3h C6h D6 C6v D3h D6h
Cubic            T Th O Td Oh
Of the symbols used for cell centering:
P - stands for Primitive. It is a cube depicted by a vertex at eight corners
C - stands for having a point filled in on the "C" side of the primitive cell
    this is described in some texts as Base Centering. C is most commonly used
    A or B means use the "A" or "B" sides instead. Just a rotation of C
F - stands for Face Centering and fills all three, "A", "B" and "C" sides
I - (from German: innenzentriert, meaning Body Centered) is the primitive cell
    with one point filled in the center of the cell
    
The term Isometric is sometimes used for Cubic. Also allowable in this program.

Monoclinic is defined in this program with angles alpha = gamma = 90 <> beta
as is found in the first volume of International Tables for Crystallography.

Vectors a, b and c correspond to axes x, y and z (before transformations).
)");
}

void brav_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] crystal_system [centering]

Generate Bravais lattices in OFF format. lattice may be specified by its index
number or (start of the) crystal_system name and centering

Crystal System   Centering   Vector Constraints  Angle Constraints
Triclinic        P           no constraints      any not of higher symmetries
Monoclinic       P,C         no constraints      alpha = gamma = 90 <> beta
Orthorhombic     P,C,F,I     a <> b <> c         alpha = beta = gamma = 90
Tetragonal       P,I         a = b <> c          alpha = beta = gamma = 90
Trigonal         P           a = b = c           alpha = beta = gamma <> 90
Hexagonal        P           a = b, c            alpha = beta = 90, gamma = 120
Cubic            P,F,I       a = b = c           alpha = beta = gamma = 90

Centering Types: P - Primitive (Simple) centering on corners (default: P)
                 C - Base (also A or B)  F - Face centering  I - Body centering

Synonyms: Rhombohedral = Trigonal  sc = Cubic P  fcc = Cubic F  bcc = Cubic I

Options
%s
  -H        additional help
  -W        verbose output
  -l <lim>  minimum distance for unique vertex locations as negative exponent
              (default: %d giving %.0e)
  -o <file> write output to file (default: write to standard output)

Lattice Options
  -v <v,n>  vector lengths, non-zero, in form "a,b,c" (default: calculated)
              optional fourth number, vectors taken to root n
  -a <angs> angles in the form "alpha,beta,gamma". Ignored for Orthorhombic,
              Tetragonal, and Cubic. For Hexagonal, any non-90 position may be
              120. Otherwise, if not supplied then random angles are chosen.
              Angles cannot be zero or 180. Angles may be negative values.
              alpha + beta + gamma must be less than 360. Each angle must be
              less than or equal to the sum of the other two angles
  -g <grid> cell grid array. one or three positive integers separated by commas
              a - automatic, of sufficient size for radius (-G required)
              one integer, NxNxN grid. (default: calculated)
              three integers, IxJxK grid. Combinations for grid center:
              even,even,even - on cell corner    odd,odd,odd - in cell body
              even,odd,even - on cell mid-edge   odd,even,odd - on face center
  -G <type> automatic grid center type (type 8 invalid for cell centering = P):
              p - corner, i - body, f - face, e - mid-edge, 8 - eighth cell
  -I        inversion (centering type F or I)
  -d <vrts> output dual of lattice based on primitive vectors
              c - use primitive vectors base on centering type
              four integers - primitive vectors are determined by four vertex
                numbers given by non negative integers. The first vertex
                number is the radial point and the next three vertices are the
                primitive vectors
  -u        add cell struts. Added to cubic grid before transformation
  -s <s,n>  create struts. s is strut length taken to optional root n
              use multiple -s parameters for multiple struts
  -D <opt>  Voronoi (a.k.a Dirichlet) cells (Brillouin zones for duals)
              c - cells only, i - cell(s) touching center only
  -A        append the original lattice to the final product

Container Options
  -c <type> container s - sphere (uses radius) (default: none)
  -k <file> container, convex hull of off file or built in model (uses radius)
  -r <c,n>  radius. c is radius taken to optional root n. n = 2 is sqrt
              or  l - max insphere radius, s - min insphere radius (default)
              or  k - take radius from container specified by -k
  -p <xyz>  radius to lattice, three comma separated coordinates, 0 for origin
  -q <xyz>  center offset, three comma separated coordinates, 0 for origin
  -C <opt>  c - convex hull only, i - keep interior

Coloring Options (run 'off_util -H color' for help on color formats)
  -V <col>  vertex color, (optional) transparency, (optional) elements
              transparency: valid range from 0 (invisible) to 255 (opaque)
              elements to color are l - lattice, c - convex hull, v - voronoi
                                    h - hex relation (default elements: lcvh)
  -E <col>  edge color (same format as for vertices)
  -F <col>  face color (same format as for vertices) or
              special coloring: calculates color from normals, not maps
              lower case outputs map indexes. upper case outputs color values
              y,Y - color by symmetry using face normals
              z,Z - color by symmetry using face normals (chiral)
  -T <tran> face transparency for color by symmetry. valid range from 0 to 255

Scene Options
  -R <opt>  hexagonal/cubic relation (Cubic P or Trigonal only)
              o - hex overlay, f - hex fill, O - cube overlay, F - cube fill
  -O        translate centroid of final product to origin
  -Z <col>  add centroid vertex to final product in color col
  -K        append cage of container of -k to final product

Listing Options
  -B        display the list of Bravais lattices
  -Q <vecs> center for radius calculations in -L (default: centroid)
              c - original center, o - original center + offset in -q
  -L <opt>  list unique radial distances of points (to standard output)
              f - full report, v - values only
  -S <opt>  list every possible strut value (to standard output)
              f - full report, v - values only

)",
          prog_name(), help_ver_text, int(-log(anti::epsilon) / log(10) + 0.5),
          anti::epsilon);
}

void brav_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  bool radius_set = false;
  vector<double> double_parms;

  vert_col.resize(4);
  edge_col.resize(4);
  face_col.resize(4);

  // set some default colors first so they can be optionally unset
  // 0 - lattice  1 - convex hull  2 - voronoi  3 - hex overlay
  // voronoi cells  vcol = gold; ecol = lightgray; fcol = transparent yellow
  vert_col[2] = Color(255, 215, 0);
  edge_col[2] = Color(211, 211, 211);
  face_col[2] = Color(255, 255, 0, 128);

  vert_col[3] = Color(255, 215, 0);
  edge_col[3] = Color(211, 211, 211);

  handle_long_opts(argc, argv);

  while ((c = getopt(
              argc, argv,
              ":hHWc:k:r:p:q:s:uv:a:g:G:Id:l:D:C:AV:E:F:T:Z:KOR:BQ:L:S:o:")) !=
         -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'H':
      extended_help();
      exit(0);

    case 'W':
      verbose = true;
      break;

    case 'c':
      if (strlen(optarg) > 1 || !strchr("s", *optarg))
        error("method is '" + string(optarg) + "' must be s", c);
      container = *optarg;
      break;

    case 'k':
      cfile = optarg;
      break;

    case 'r':
      radius_set = true;
      if (strlen(optarg) == 1 && strchr("lsk", *optarg))
        radius_default = *optarg;
      else {
        print_status_or_exit(read_double_list(optarg, double_parms, 2), c);
        radius = double_parms[0];
        if (double_parms.size() == 2) {
          if (double_parms[1] == 0)
            error("root for radius must be non-zero", c);
          radius = pow(radius, 1 / double_parms[1]);
        }
        if (radius <= 0)
          error("radius cannot be negative or zero", c);
      }
      break;

    case 'p':
      print_status_or_exit(radius_by_coord.read_maths(optarg), c);
      break;

    case 'q':
      print_status_or_exit(offset.read_maths(optarg), c);
      break;

    case 's': {
      double strut_tmp;
      print_status_or_exit(read_double_list(optarg, double_parms, 2), c);
      strut_tmp = double_parms[0];
      if (double_parms.size() == 2) {
        if (double_parms[1] == 0)
          error("root for strut must be non-zero", c);
        strut_tmp = pow(strut_tmp, 1 / double_parms[1]);
      }
      if (strut_tmp <= 0)
        error("strut lengths cannot be negative", "s", c);
      strut_len.push_back(strut_tmp);
      break;
    }

    case 'u':
      cell_struts = true;
      break;

    case 'v':
      print_status_or_exit(read_double_list(optarg, vecs, 4), c);
      if (vecs.size() < 3)
        error(msg_str("three vector lengths needed (%lu were given)",
                      (unsigned long)vecs.size()),
              c);
      if (vecs[0] == 0 || vecs[1] == 0 || vecs[2] == 0)
        error("vector lengths need to be non-zero", c);
      if (vecs.size() == 4) {
        if (vecs[3] == 0)
          error("root for vector must be non-zero", c);
        for (unsigned int i = 0; i < 3; i++)
          vecs[i] = pow(vecs[i], 1 / vecs[3]);
      }
      break;

    case 'a': {
      print_status_or_exit(read_double_list(optarg, angles, 3), c);
      if (angles.size() < 3)
        error(msg_str("three angles needed (%lu were given)",
                      (unsigned long)angles.size()),
              c);

      double alpha = fabs(angles[0]);
      double beta = fabs(angles[1]);
      double gamma = fabs(angles[2]);

      if (fmod(alpha, 180.0) == 0 || fmod(beta, 180.0) == 0 ||
          fmod(gamma, 180.0) == 0)
        error("angles cannot be zero or 180", c);
      if (alpha + beta + gamma >= 360.0)
        error("alpha + beta + gamma must be less than 360.0", c);
      if (alpha >= beta + gamma)
        error("alpha must be less than beta + gamma", c);
      if (beta >= alpha + gamma)
        error("beta must be less than alpha + gamma", c);
      if (gamma >= alpha + beta)
        error("gamma must be less than alpha + beta", c);
      break;
    }

    case 'g':
      if (strchr("a", *optarg))
        grid_for_radius = true;
      else {
        print_status_or_exit(read_int_list(optarg, grid, true, 3), c);
        if (grid.size() == 1) {
          grid.push_back(grid[0]);
          grid.push_back(grid[0]);
        }
        else if (grid.size() != 3)
          error(msg_str("must give one or three numbers (%lu were "
                        "given)",
                        (unsigned long)grid.size()),
                c);
        if (grid[0] == 0 || grid[1] == 0 || grid[2] == 0)
          error("grid requires positive integer(s)", c);
        // if(grid[0] != grid[1] || grid[0] != grid[2] || grid[1] != grid[2])
        //   warning("grid symmetry may be less than cells symmetry", c);
      }
      break;

    case 'G':
      if (strlen(optarg) > 1 || !strchr("pfie8", *optarg))
        error("grid type is '" + string(optarg) + "' must be p, f, i, e, or 8",
              c);
      auto_grid_type = *optarg;
      break;

    case 'I':
      inversion = true;
      break;

    case 'd':
      if (strlen(optarg) == 1) {
        if (strchr("c", *optarg))
          use_centering_for_dual = true;
        else
          error("invalid option", c);
      }
      else {
        print_status_or_exit(read_int_list(optarg, prim_vec_idxs, true, 4), c);
        if (prim_vec_idxs.size() != 4)
          error(msg_str("four lattice vertex indices needed (%lu were"
                        "given)",
                        (unsigned long)prim_vec_idxs.size()),
                c);
        else {
          for (unsigned int i = 0; i < 3; i++) {
            for (unsigned int j = i + 1; j < 4; j++) {
              if (prim_vec_idxs[i] == prim_vec_idxs[j])
                error(msg_str("lattice vertex indices used twice. "
                              "pos: %d & %d",
                              i + 1, j + 1),
                      c);
            }
          }
        }
      }
      break;

    case 'l':
      int sig_compare;
      print_status_or_exit(read_int(optarg, &sig_compare), c);
      if (sig_compare > DEF_SIG_DGTS)
        warning("limit is very small, may not be attainable", c);
      eps = pow(10, -sig_compare);
      break;

    case 'D':
      if (strlen(optarg) > 1 || !strchr("ci", *optarg))
        error("Voronoi cells arg is '" + string(optarg) + "' must be c, i", c);
      voronoi_cells = true;
      if (strchr("i", *optarg))
        voronoi_central_cell = true;
      break;

    case 'C':
      if (strlen(optarg) > 1 || !strchr("ci", *optarg))
        error("convex hull arg is '" + string(optarg) + "' must be c, i", c);
      convex_hull = true;
      if (strchr("i", *optarg))
        add_hull = true;
      break;

    case 'A':
      append_lattice = true;
      break;

    case 'V':
      parse_color_string(this, optarg, c, "lcvh", vert_col);
      break;

    case 'E':
      parse_color_string(this, optarg, c, "lcvh", edge_col);
      break;

    case 'F': {
      Split parts(optarg, ",");
      if (strlen(parts[0]) == 1 && strchr("yYzZ", parts[0][0]))
        color_method = parts[0][0];
      else
        parse_color_string(this, optarg, c, "lcvh", face_col);
      break;
    }

    case 'T':
      print_status_or_exit(read_int(optarg, &face_opacity), c);
      if (face_opacity < 0 || face_opacity > 255) {
        error("face transparency must be between 0 and 255", c);
      }
      break;

    case 'Z':
      print_status_or_exit(cent_col.read(optarg), c);
      break;

    case 'K':
      append_container = true;
      break;

    case 'O':
      trans_to_origin = true;
      break;

    case 'R':
      if (strlen(optarg) > 1 || !strchr("ofOF", *optarg))
        error("convex hull arg is '" + string(optarg) +
                  "' must be o, f, O or F",
              c);
      if (strchr("o", *optarg))
        r_lattice_type = 1;
      else if (strchr("f", *optarg))
        r_lattice_type = 2;
      else if (strchr("O", *optarg))
        r_lattice_type = 3;
      else if (strchr("F", *optarg))
        r_lattice_type = 4;
      break;

    case 'B':
      list_bravais = true;
      break;

    case 'Q':
      if (strlen(optarg) == 1) {
        if (strchr("co", *optarg))
          list_radii_original_center = *optarg;
        else
          error("invalid option", c);
      }
      else
        print_status_or_exit(list_radii_center.read_maths(optarg), c);
      break;

    case 'L':
      if (strlen(optarg) > 1 || !strchr("fv", *optarg))
        error("list radii arg is '" + string(optarg) + "' must be f or v", c);
      if (strchr("f", *optarg))
        list_radii = 1;
      else if (strchr("v", *optarg))
        list_radii = 2;
      break;

    case 'S':
      if (strlen(optarg) > 1 || !strchr("fv", *optarg))
        error("list struts arg is '" + string(optarg) + "' must be f or v", c);
      if (strchr("f", *optarg))
        list_struts = 1;
      else if (strchr("v", *optarg))
        list_struts = 2;
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (argc == optind && !list_bravais)
    error("no lattice specified", "lattice");

  if (argc - optind > 2) {
    error("too many arguments");
    exit(1);
  }

  if (optind < argc)
    crystal_system = argv[optind++];
  if (optind < argc)
    centering = argv[optind++];

  transform(crystal_system.begin(), crystal_system.end(),
            crystal_system.begin(), ::tolower);
  transform(centering.begin(), centering.end(), centering.begin(), ::tolower);

  // don't let centering in when number is entered
  double dummy;
  if (read_double(crystal_system.c_str(), &dummy)) {
    if (centering.length())
      error("too many arguments");
  }

  // have to catch this too
  if (crystal_system == "sc" || crystal_system == "fcc" ||
      crystal_system == "bcc") {
    if (centering.length())
      error("too many arguments");
  }

  if (grid_for_radius && !auto_grid_type)
    error("-G required with -g a", 'g');

  if (radius_set and radius_by_coord.is_set())
    error("-p cannot be used with -r", 'p');

  if (container == 's' && cfile.length())
    error("-c and -k cannot be specified together", 'c');

  if (append_container && !cfile.length())
    error("container can only be appended if one is provided with -k", 'K');

  if (radius_default == 'k' && !cfile.length())
    error("-r k can only be used if -k container is specified", 'r');

  if (container == 'c' && !cfile.length() &&
      (radius != 0 || offset.is_set() || radius_by_coord.is_set()))
    warning("no container in use. Radius parameters ignored", 'c');

  if (use_centering_for_dual || prim_vec_idxs.size())
    dual_lattice = true;

  if (dual_lattice && r_lattice_type > 0)
    error("hex relation struts cannot be done on duals", 'R');

  if (list_radii && list_struts)
    error("cannot list radii and struts at the same time", 'L');

  // for convex hull
  if (face_opacity > -1) {
    // if no face color is set but there is transparency set, use white
    if (!face_col[1].is_set())
      face_col[1] = Color(255, 255, 255);

    // face color can only be made transparent if not index and not invisible
    if (!face_col[1].set_alpha(face_opacity))
      warning("transparency has no effect on map indexes or invisible", "T");

    if (color_method == 's' || color_method == 'c')
      warning("when writing indexes transparency setting ignored", "T");
  }
}

bravais::bravais()
{
  bravais_items = bravais_item_list;
  last_bravais = sizeof(bravais_item_list) / sizeof(bravais_item_list[0]);
}

string bravais::get_crystal_system(int i)
{
  return bravais_items[i].crystal_system;
}

string bravais::get_centering(int i) { return bravais_items[i].centering; }

void bravais::list_bravais(int idx, FILE *fp)
{
  fprintf(fp, "%2d)\t%-s %s\n", idx + 1, bravais_items[idx].crystal_system,
          bravais_items[idx].centering);
}

void bravais::list_bravais(FILE *fp)
{
  for (int i = 0; i < last_bravais; i++)
    list_bravais(i, fp);
}

int bravais::lookup_sym_no(string crystal_system, string centering)
{
  // everything expects lower case
  transform(crystal_system.begin(), crystal_system.end(),
            crystal_system.begin(), ::tolower);
  transform(centering.begin(), centering.end(), centering.begin(), ::tolower);

  // is it blank
  if (!crystal_system.length())
    return -1;

  // is it the list order number
  char *endptr;
  int idx = strtol(crystal_system.c_str(), &endptr, 10);
  if (!*endptr) // all of string is an integer
    return idx - 1;

  idx = -1;

  // is it a partial name
  for (int i = 0; i < last_bravais; i++) {
    string name = get_crystal_system(i);
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    if (strncmp(crystal_system.c_str(), name.c_str(), crystal_system.size()) ==
        0) {
      crystal_system = get_crystal_system(i);
      break;
    }
  }

  // is it rhombohedral or isometric
  if (strncmp(crystal_system.c_str(), "rhombohedral", crystal_system.size()) ==
      0)
    crystal_system = "trigonal";
  else if (strncmp(crystal_system.c_str(), "isometric",
                   crystal_system.size()) == 0)
    crystal_system = "cubic";

  // well known synonyms
  if (crystal_system == "sc") {
    crystal_system = "cubic";
    centering = "p";
  }
  else if (crystal_system == "fcc") {
    crystal_system = "cubic";
    centering = "f";
  }
  else if (crystal_system == "bcc") {
    crystal_system = "cubic";
    centering = "i";
  }

  // default centering is simple
  if (!centering.length())
    centering = "p";

  string full_name = crystal_system + centering;

  // is it a lattice name
  transform(full_name.begin(), full_name.end(), full_name.begin(), ::tolower);
  for (int i = 0; i < last_bravais; i++) {
    string name = get_crystal_system(i) + get_centering(i);
    transform(name.begin(), name.end(), name.begin(), ::tolower);
    if (full_name == name)
      return i;
  }

  return idx;
}

// formulas from http://en.wikipedia.org/wiki/Bravais_lattice
// crystal_system is not changed
double bravais_volume(string crystal_system, const vector<double> &vecs,
                      const vector<double> &angles)
{
  transform(crystal_system.begin(), crystal_system.end(),
            crystal_system.begin(), ::tolower);

  double volume = 0;
  // find crystal system
  if (crystal_system == "triclinic") {
    // a * b * c *
    // sqrt(1-cos(alpha)^2-cos(beta)^2-cos(gamma)^2+2*cos(alpha)*cos(beta)*cos(gamma))
    double cos_a = cos(angles[0] * M_PI / 180);
    double cos_b = cos(angles[0] * M_PI / 180);
    double cos_g = cos(angles[0] * M_PI / 180);
    double tmp = 1 - cos_a * cos_a - cos_b * cos_b - cos_g * cos_g +
                 2 * cos_a * cos_b * cos_g;
    volume = vecs[0] * vecs[1] * vecs[2] * sqrt(tmp);
  }
  else if (crystal_system == "monoclinic" || crystal_system == "hexagonal") {
    double a = 0;
    for (unsigned int i = 0; i < 3; i++) {
      if (angles[i] != 90) {
        a = angles[i];
        break;
      }
    }
    volume = vecs[0] * vecs[1] * vecs[2] *
             sin(a * M_PI / 180); // a * b * c * sin(alpha)
  }
  else if (crystal_system == "orthorhombic")
    volume = vecs[0] * vecs[1] * vecs[2]; // a * b * c
  else if (crystal_system == "tetragonal")
    volume = vecs[0] * vecs[0] * vecs[2]; // a^2 * c
  else if (crystal_system == "trigonal") {
    // a^3 * sqrt(1-3*cos(alpha)^2+2*cos(alpha)^3)
    double cos_a = cos(angles[0] * M_PI / 180);
    double tmp = 1 - 3 * cos_a * cos_a + 2 * cos_a * cos_a * cos_a;
    volume = vecs[0] * vecs[0] * vecs[0] * sqrt(tmp);
  }
  else
    // incorrect
    // if ( crystal_system == "hexagonal" )
    //   volume = (3.0*pow(3.0,(1.0/3.0))*vecs[0]*vecs[0]*vecs[2])/2; //
    //   (3*cube_rt(3)*a^2*c)/2
    // else
    if (crystal_system == "cubic")
      volume = vecs[0] * vecs[0] * vecs[0]; // a^3
    else {
      fprintf(stderr, "error in bravais_volume: unknown crystal system: %s\n",
              crystal_system.c_str());
    }

  return volume;
}

// return random angle
double bravais_random_angle(Random &ran, const double max_angle)
{
  return ran.ran_in_range_exclude_end(0.001,
                                      max_angle); // avoid 0 and max_angle
}

// sort the three elements without altering original order
void sort_three(double &x, double &y, double &z, vector<double> v)
{
  sort(v.begin(), v.end());

  x = v[0];
  y = v[1];
  z = v[2];
}

int bravais_check(string &crystal_system, string &centering,
                  vector<double> &vecs, vector<double> &angles,
                  const int strictness, bool verbose)
{
  // strictness = 0 any change in crystal system allowed
  // strictness = 1 only upgrades in crystal system allowed
  // strictness = 2 no change in crystal system allowed

  // Crystal System       Vector Constraints      Angle Constraints
  // Triclinic        (1) doesn't matter      (1) alpha <> beta <> gamma
  //                                          (1) alpha = beta <> 90 <> gamma
  //                  (2) a <> b <> c         (4) alpha = beta = gamma <> 90
  //                  (special case)
  //                  (3) a = b <> c          (4) alpha = beta = gamma <> 90
  //                  (special case)
  // Monoclinic       (1) doesn't matter      (2) alpha = gamma = 90 <> beta
  //                  (2) a <> b <> c         (5) alpha = beta = 90, gamma = 120
  //                  (special case)
  // Orthorhombic     (2) a <> b <> c         (3) alpha = beta = gamma = 90
  // Tetragonal       (3) a = b <> c          (3) alpha = beta = gamma = 90
  // Trigonal         (4) a = b = c           (4) alpha = beta = gamma <> 90
  // Note on Hexagonal: odd vector position must line up with odd angle
  // position, else monoclinic
  // Hexagonal        (3) a = b <> c          (5) alpha = beta = 90, gamma = 120
  // (or 60)
  //                  (4) a = b = c           (5) alpha = beta = 90, gamma = 120
  //                  (or 60)
  // Cubic            (4) a = b = c           (3) alpha = beta = gamma = 90

  // load makeshift "database"
  vector<string> crystal_systems;
  crystal_systems.push_back("unknown");
  crystal_systems.push_back("triclinic");
  crystal_systems.push_back("monoclinic");
  crystal_systems.push_back("orthorhombic");
  crystal_systems.push_back("tetragonal");
  crystal_systems.push_back("trigonal");
  crystal_systems.push_back("hexagonal");
  crystal_systems.push_back("cubic");

  vector<string> vec_cases;
  vec_cases.push_back("unknown vector case");
  vec_cases.push_back("no constraints");
  vec_cases.push_back("a <> b <> c");
  vec_cases.push_back("a = b <> c");
  vec_cases.push_back("a = b = c");

  vector<string> angle_cases;
  angle_cases.push_back("unknown angle case");
  angle_cases.push_back("anything not in higher symmetries");
  angle_cases.push_back("alpha = gamma = 90 <> beta");
  angle_cases.push_back("alpha = beta = gamma = 90");
  angle_cases.push_back("alpha = beta = gamma <> 90");
  angle_cases.push_back("alpha = beta = 90, gamma = 120");

  vector<int> vec_rule;
  vec_rule.push_back(0);
  vec_rule.push_back(1);
  vec_rule.push_back(1);
  vec_rule.push_back(2);
  vec_rule.push_back(3);
  vec_rule.push_back(4);
  vec_rule.push_back(3);
  vec_rule.push_back(4);

  vector<int> angle_rule;
  angle_rule.push_back(0);
  angle_rule.push_back(1);
  angle_rule.push_back(2);
  angle_rule.push_back(3);
  angle_rule.push_back(3);
  angle_rule.push_back(4);
  angle_rule.push_back(5);
  angle_rule.push_back(3);

  // in case function is called with upper case lettering
  transform(crystal_system.begin(), crystal_system.end(),
            crystal_system.begin(), ::tolower);
  transform(centering.begin(), centering.end(), centering.begin(), ::tolower);

  // find chosen crystal system
  int crystal_system_index = 0;
  for (int i = 1; i <= 7; i++) {
    if (crystal_system == crystal_systems[i]) {
      crystal_system_index = i;
      break;
    }
  }

  // in case of misspellings
  if (!crystal_system_index) {
    fprintf(stderr, "error in bravais_check: unknown crystal system: %s\n",
            crystal_system.c_str());
    return 0;
  }

  double a = 0;
  double b = 0;
  double c = 0;

  // if vecs not set, give valid values
  if (!vecs.size()) {
    // hexagonal case gets bumped up to all equal vectors
    if (vec_rule[crystal_system_index] == 4 || crystal_system_index == 6) {
      a = 2;
      b = 2;
      c = 2;
    }
    else if (vec_rule[crystal_system_index] == 3) {
      a = 2;
      b = 2;
      c = 4;
    }
    else if (vec_rule[crystal_system_index] == 2 ||
             vec_rule[crystal_system_index] == 1) {
      a = 2;
      b = 3;
      c = 4;
    }

    // if simple cubic then space them 1 apart instead
    if (centering == "p") {
      a /= 2.0;
      b /= 2.0;
      c /= 2.0;
    }

    vecs.push_back(a);
    vecs.push_back(b);
    vecs.push_back(c);
  }

  if (verbose)
    fprintf(stderr, "info: vector a = %g, b = %g, c = %g\n", vecs[0], vecs[1],
            vecs[2]);

  double alpha = 0;
  double beta = 0;
  double gamma = 0;

  if (angles.size())
    sort_three(alpha, beta, gamma, angles);

  // fixed angle cases, otherwise if not given, random angles assigned
  if (angle_rule[crystal_system_index] == 3) {
    if ((alpha && alpha != 90.0) || (beta && beta != 90.0) ||
        (gamma && gamma != 90.0))
      fprintf(stderr,
              "warning: for crystal system %s, alpha, beta and gamma "
              "must be 90 degrees. -a ignored\n",
              crystal_system.c_str());
    alpha = 90;
    beta = 90;
    gamma = 90;
    angles.clear();
    angles.push_back(alpha);
    angles.push_back(beta);
    angles.push_back(gamma);
  }
  else
    // hexagonal, 120 can be in any position (alpha, beta, gamma are sorted
    // order)
    if (angle_rule[crystal_system_index] == 5 &&
        !(alpha == 90 && beta == 90 && gamma == 120)) {
      fprintf(stderr,
              "warning: for crystal system %s, of alpha, beta and gamma, "
              "two values must be 90 degrees and one value 120 degrees. "
              "-a ignored\n",
              crystal_system.c_str());
      alpha = 90;
      beta = 90;
      gamma = 120;
      angles.clear();
      angles.push_back(alpha);
      angles.push_back(beta);
      angles.push_back(gamma);
    }
    else
      // angles can be set to random values for other crystal systems
      if (!angles.size()) {
        Random ran;
        ran.time_seed();

        if (crystal_system == "monoclinic") {
          alpha = 90.0;
          beta = 90.0;
          gamma = 90.0;
          // don't let beta equal 90
          while (beta == 90.0)
            beta = bravais_random_angle(ran, 180.0);
        }
        else if (crystal_system == "trigonal") {
          alpha = bravais_random_angle(ran, 60.0);
          beta = alpha;
          gamma = alpha;
        }
        else if (crystal_system == "triclinic") {
          alpha = 360.0;
          beta = 360.0;
          gamma = 360.0;
          while (alpha + beta + gamma >= 360 || gamma >= alpha + beta ||
                 beta >= alpha + gamma || alpha >= beta + gamma) {
            alpha = bravais_random_angle(ran, 180.0);
            beta = bravais_random_angle(ran, 180.0);
            gamma = bravais_random_angle(ran, 180.0);
          }
        }

        angles.push_back(alpha);
        angles.push_back(beta);
        angles.push_back(gamma);

        fprintf(stderr, "warning: RANDOM angles generated\n");
      }

  if (verbose)
    fprintf(stderr, "info: angles alpha = %5.5g, beta = %5.5g, gamma = %5.5g\n",
            angles[0], angles[1], angles[2]);

  // for evaulation, sort vector lengths and angles, low to high (if parameters
  // are set)
  sort_three(a, b, c, vecs);
  sort_three(alpha, beta, gamma, angles);

  int vec_case = 0;
  if (a != b && b != c)
    vec_case = 2;
  else if ((a == b && b != c) || (a != b && b == c)) // e.g. 1,1,2 or 1,2,2
    vec_case = 3;
  else if (a == b && b == c)
    vec_case = 4;

  int angle_case = 0;
  if ((alpha == 90 && beta == 90 && gamma == 120) ||
      (alpha == 60 && beta == 90 && gamma == 90))
    angle_case = 5;
  else if (alpha == 90 && beta == 90 && gamma == 90)
    angle_case = 3;
  else if ((alpha == 90 && beta == 90 && gamma != 90) ||
           (alpha != 90 && beta == 90 &&
            gamma == 90)) // e.g. 90,90,100 or 45,90,90
    angle_case = 2;
  else if (alpha == beta && beta == gamma && alpha != 90 && beta != 90 &&
           gamma != 90)
    angle_case = 4;
  else
    // anything else will base case 1
    angle_case = 1;

  // odd vector position must line up with odd angle position, else monoclinic
  if ((vec_case == 3 || vec_case == 4) && angle_case == 5) {
    int odd_angle = 0;
    for (; odd_angle < 3; odd_angle++) {
      if (angles[odd_angle] != 90.0)
        break;
    }
    if (vecs[(odd_angle + 1) % 3] != vecs[(odd_angle + 2) % 3]) {
      angle_case = 2;
      vec_case = 1;
    }
  }

  // find actual crystal system by vec and angle input
  int csystem = 0;
  // triclinic
  if ((angle_case == 1) || // vec_case can be any
      (vec_case == 2 && angle_case == 4) ||
      (vec_case == 3 &&
       angle_case == 4)) { // have to allow for triclinic with unequal vectors
    csystem = 1;
    vec_case = 1;
  }
  else
    // monoclinic
    if ((angle_case == 2) ||                  // vec_case can be any
        (vec_case == 2 && angle_case == 5)) { // have to allow for hexagonal
                                              // angles but three unequal
                                              // vectors
      csystem = 2;
      vec_case = 1;
      angle_case = 2; // in case hex angles were present
    }
    else
      // orthorhombic
      if (vec_case == 2 && angle_case == 3)
        csystem = 3;
      else
        // tetragonal
        if (vec_case == 3 && angle_case == 3)
          csystem = 4;
        else
          // trigonal
          if (vec_case == 4 && angle_case == 4)
            csystem = 5;
          else
            // hexagonal
            if ((vec_case == 3 || vec_case == 4) && angle_case == 5)
              csystem = 6;
            else
              // cubic
              if (vec_case == 4 && angle_case == 3)
                csystem = 7;

  if (crystal_system_index != csystem) {
    fprintf(stderr, "\n");
    fprintf(stderr, "warning in bravais_check: MISMATCH ...\n\n");
    fprintf(stderr, "asked for '%s' which requires:\n",
            crystal_systems[crystal_system_index].c_str());
    fprintf(stderr, "vector lengths: '%s'\n",
            vec_cases[vec_rule[crystal_system_index]].c_str());
    fprintf(stderr, "angles: '%s'\n",
            angle_cases[angle_rule[crystal_system_index]].c_str());
    fprintf(stderr, "\nbut instead found '%s' which has:\n",
            crystal_systems[csystem].c_str());
    fprintf(stderr, "vector lengths: '%s'\n", vec_cases[vec_case].c_str());
    fprintf(stderr, "angles: '%s'\n\n", angle_cases[angle_case].c_str());
  }

  // this should not be able to happen
  if (csystem == 0) {
    fprintf(stderr, "error in bravais_check: unknown crystal system case\n");
    return 0;
  }

  // check if system change occured. if strict is other than zero, may return
  // error state
  if (csystem != crystal_system_index) {
    // have to update crystal system name for other functioning
    crystal_system = crystal_systems[csystem];

    if (strictness == 1 && csystem < crystal_system_index) {
      fprintf(stderr, "error in bravais_check: crystal system downgrade "
                      "disallowed. cannot continue\n");
      return 0;
    }
    else if (strictness == 2) {
      fprintf(stderr, "error in bravais_check: strict type matching in force. "
                      "cannot continue\n");
      return 0;
    }
  }

  // only "monoclinic" or "orthorhombic" can have Base Centering (C)
  if (csystem != 2 && csystem != 3 && centering == "c") {
    fprintf(stderr,
            "bravais_check: warning: %s cannot have Base Centering (C)\n",
            crystal_systems[csystem].c_str());
    if (strictness > 0) {
      fprintf(stderr, "bravais_check: error: crystal system downgrade "
                      "disallowed. cannot continue\n");
      return 0;
    }
    else {
      fprintf(stderr,
              "bravais_check: warning: downgrade to Simple Centering (P)\n");
      centering = "p";
    }
  }

  return 1;
}

void bravais_centering_p(Geometry &geom)
{
  geom.add_vert(Vec3d(1, 1, 1));    // 0
  geom.add_vert(Vec3d(1, 1, -1));   // 1
  geom.add_vert(Vec3d(1, -1, 1));   // 2
  geom.add_vert(Vec3d(1, -1, -1));  // 3
  geom.add_vert(Vec3d(-1, 1, 1));   // 4
  geom.add_vert(Vec3d(-1, 1, -1));  // 5
  geom.add_vert(Vec3d(-1, -1, 1));  // 6
  geom.add_vert(Vec3d(-1, -1, -1)); // 7
}

void bravais_centering_a(Geometry &geom)
{
  bravais_centering_p(geom);

  geom.add_vert(Vec3d(1, 0, 0));  // 8
  geom.add_vert(Vec3d(-1, 0, 0)); // 9
}

void bravais_centering_b(Geometry &geom)
{
  bravais_centering_p(geom);

  geom.add_vert(Vec3d(0, 1, 0));  // 8
  geom.add_vert(Vec3d(0, -1, 0)); // 9
}

void bravais_centering_c(Geometry &geom)
{
  bravais_centering_p(geom);

  geom.add_vert(Vec3d(0, 0, 1));  // 8
  geom.add_vert(Vec3d(0, 0, -1)); // 9
}

void bravais_centering_f(Geometry &geom)
{
  bravais_centering_p(geom);

  geom.add_vert(Vec3d(1, 0, 0));  // 8
  geom.add_vert(Vec3d(-1, 0, 0)); // 9
  geom.add_vert(Vec3d(0, 1, 0));  // 10
  geom.add_vert(Vec3d(0, -1, 0)); // 11
  geom.add_vert(Vec3d(0, 0, 1));  // 12
  geom.add_vert(Vec3d(0, 0, -1)); // 13
}

void bravais_centering_i(Geometry &geom)
{
  bravais_centering_p(geom);
  geom.add_vert(Vec3d(0, 0, 0)); // 8
}

void bravais_centering_f_inverted(Geometry &geom)
{
  geom.add_vert(Vec3d(-1, -1, 0)); // 1
  geom.add_vert(Vec3d(-1, 0, -1)); // 2
  geom.add_vert(Vec3d(-1, 0, 1));  // 3
  geom.add_vert(Vec3d(-1, 1, 0));  // 4
  geom.add_vert(Vec3d(0, -1, -1)); // 5
  geom.add_vert(Vec3d(0, -1, 1));  // 6
  geom.add_vert(Vec3d(0, 0, 0));   // 7
  geom.add_vert(Vec3d(0, 1, -1));  // 8
  geom.add_vert(Vec3d(0, 1, 1));   // 9
  geom.add_vert(Vec3d(1, -1, 0));  // 10
  geom.add_vert(Vec3d(1, 0, -1));  // 11
  geom.add_vert(Vec3d(1, 0, 1));   // 12
  geom.add_vert(Vec3d(1, 1, 0));   // 13
}

void bravais_centering_i_inverted(Geometry &geom)
{
  geom.add_vert(Vec3d(-1, 0, 0));  // 1
  geom.add_vert(Vec3d(0, -1, -1)); // 2
  geom.add_vert(Vec3d(0, -1, 1));  // 3
  geom.add_vert(Vec3d(0, 1, -1));  // 4
  geom.add_vert(Vec3d(0, 1, 1));   // 5
  geom.add_vert(Vec3d(1, 0, 0));   // 6
}

void bravais_cell_struts(Geometry &geom, const Color &edge_col)
{
  int f[] = {0, 1, 0, 2, 0, 4, 3, 1, 3, 2, 3, 7,
             5, 1, 5, 4, 5, 7, 6, 2, 6, 4, 6, 7};
  vector<int> edge(2);
  for (int i = 0; i < 12; i++) {
    for (int j = 0; j < 2; j++)
      edge[j] = f[i * 2 + j];
    geom.add_edge(edge, edge_col);
  }
}

void bravais_cell(Geometry &geom, const string &centering,
                  const bool cell_struts, const Color &vert_col,
                  const Color &edge_col, const bool inversion)
{
  if (inversion) {
    if (centering == "f")
      bravais_centering_f_inverted(geom);
    else if (centering == "i")
      bravais_centering_i_inverted(geom);
  }
  else if (centering == "p")
    bravais_centering_p(geom);
  else if (centering == "a")
    bravais_centering_a(geom);
  else if (centering == "b")
    bravais_centering_b(geom);
  else if (centering == "c")
    bravais_centering_c(geom);
  else if (centering == "f")
    bravais_centering_f(geom);
  else if (centering == "i")
    bravais_centering_i(geom);

  if (vert_col.is_set()) {
    Coloring vc(&geom);
    vc.v_one_col(vert_col);
  }

  if (cell_struts & !inversion)
    bravais_cell_struts(geom, edge_col);
}

void bravais_scale(Geometry &geom, const vector<double> &vecs,
                   const bool inverse)
{
  // divide by 2 since original scale is 2
  Trans3d m = Trans3d::scale(vecs[0] / 2, vecs[1] / 2, vecs[2] / 2);
  geom.transform((inverse) ? m.inverse() : m);
}

// applies angles alpha,beta and gamma without changing scale
void bravais_warp(Geometry &geom, const vector<double> &angles,
                  const bool inverse)
{
  double yz = deg2rad(angles[0]); // alpha
  double zx = deg2rad(angles[1]); // beta
  double xy = deg2rad(angles[2]); // gamma

  Trans3d m = Trans3d::angles_between_axes(yz, zx, xy);
  geom.transform((inverse) ? m.inverse() : m);
}

// vecs is not changed
double bravais_radius(const vector<int> &grid, vector<double> vecs,
                      const vector<double> &angles, const char radius_default)
{
  Geometry tgeom;
  tgeom.read_resource("std_cube");

  // compensate for grid size
  for (unsigned int i = 0; i < 3; i++)
    vecs[i] *= grid[i];

  bravais_scale(tgeom, vecs, false);
  bravais_warp(tgeom, angles, false);

  return lattice_radius(tgeom, radius_default);
}

double bravais_auto_grid_size(const double radius, const vector<double> &vecs,
                              const vector<double> &angles)
{
  vector<int> grid(3);
  for (unsigned int i = 0; i < 3; i++)
    grid[i] = 1;
  return (radius / bravais_radius(grid, vecs, angles, 's'));
}

// radius_by_coord is not changed
double bravais_radius_by_coord(Vec3d radius_by_coord, const Vec3d &offset,
                               const vector<double> &vecs,
                               const vector<double> &angles)
{
  Geometry tgeom;
  tgeom.add_vert(radius_by_coord);

  // do same scale and warp as happened to lattice
  bravais_scale(tgeom, vecs, false);
  bravais_warp(tgeom, angles, false);

  const vector<Vec3d> &tverts = tgeom.verts();
  radius_by_coord = tverts[0];
  // radius_by_coord.dump("radius_by_coord");

  Vec3d cent = Vec3d(0, 0, 0);
  if (offset.is_set())
    cent += offset;

  return ((cent - radius_by_coord).len());
}

// makes local copy of tgeom
void geom_to_grid_translate(Geometry &geom, Geometry tgeom,
                            const Trans3d &transl_matrix)
{
  tgeom.transform(transl_matrix);
  geom.append(tgeom);
}

void geom_to_grid(Geometry &geom, const vector<int> &grid,
                  const vector<double> &cell_size, const double eps)
{
  Geometry tgeom = geom;
  geom.clear_all();

  for (int x = 0; x < grid[0]; x++) {
    for (int y = 0; y < grid[1]; y++) {
      for (int z = 0; z < grid[2]; z++) {
        geom_to_grid_translate(
            geom, tgeom,
            Trans3d::translate(
                Vec3d(cell_size[0] * x, cell_size[1] * y, cell_size[2] * z)));
      }
    }
  }

  merge_coincident_elements(geom, "vef", eps);
}

void bravais_eighth_cell_grid(Geometry &geom)
{
  const vector<Vec3d> &verts = geom.verts();

  vector<int> del_verts;
  for (unsigned int i = 0; i < verts.size(); i++) {
    for (unsigned int j = 0; j < 3; j++) {
      if (verts[i][j] < 0) {
        del_verts.push_back(i);
        break;
      }
    }
  }

  if (del_verts.size())
    geom.del(VERTS, del_verts);
}

Trans3d r_lattice_trans_mat(const bool inverse)
{
  Trans3d trans_m;
  // transformation matrix by Adrian Rossiter
  // off_trans -X 1,1,1,0,-1,1,-1,0,1
  trans_m[0] = 1;
  trans_m[1] = 1;
  trans_m[2] = 1;
  trans_m[4] = 0;
  trans_m[5] = -1;
  trans_m[6] = 1;
  trans_m[8] = -1;
  trans_m[9] = 0;
  trans_m[10] = 1;

  /*
  // inverse of above matrix will be like this below
  // off_trans -X 1/3,1/3,-2/3,1/3,-2/3,1/3,1/3,1/3,1/3
  double one_third = 1.0/3.0;
  double two_thirds = 2.0/3.0;
  trans_m[0]=one_third; trans_m[1]=one_third;  trans_m[2]=-two_thirds;
  trans_m[4]=one_third; trans_m[5]=-two_thirds; trans_m[6]=one_third;
  trans_m[8]=one_third; trans_m[9]=one_third;  trans_m[10]=one_third;
  */

  return (inverse) ? trans_m.inverse() : trans_m;
}

void r_lattice_overlay(Geometry &geom, const vector<int> &grid,
                       const vector<double> &cell_size, const Color &vert_col,
                       const Color &edge_col, const double eps)
{
  Geometry hgeom;

  bravais_cell(hgeom, "p", true, vert_col, edge_col,
               false); // hard code struts and colors for the "R" lattice
  geom_to_grid(hgeom, grid, cell_size, eps);

  hgeom.transform(r_lattice_trans_mat(false));

  // translate hexagonal onto cube
  Vec3d adjust = Vec3d(0, 0, 0);
  if (!is_even(grid[1]))
    adjust += Vec3d(0, 1, 0);
  if (!is_even(grid[0]))
    adjust += Vec3d(0, 0, 1);
  hgeom.transform(Trans3d::translate(-centroid(hgeom.verts()) +
                                     centroid(geom.verts()) - adjust));

  geom.append(hgeom);
}

void bravais_grid_type(vector<int> &grid, const char auto_grid_type)
{
  string pattern; // 0 - even 1 - odd
  if (auto_grid_type == 'p')
    pattern = "000";
  else if (auto_grid_type == 'i')
    pattern = "111";
  else if (auto_grid_type == 'e')
    pattern = "010";
  else if (auto_grid_type == 'f')
    pattern = "101";

  for (unsigned int i = 0; i < 3; i++) {
    if ((pattern[i] == '0' && !is_even(grid[i])) ||
        (pattern[i] == '1' && is_even(grid[i])))
      grid[i]++;
  }
}

void bravais_primitive_vectors(vector<Vec3d> &primitive_vectors,
                               const string &centering)
{
  double len = sqrt(2);
  primitive_vectors.push_back(Vec3d(0, 0, 0));
  if (centering == "p") {
    primitive_vectors.push_back(Vec3d(2, 0, 0));
    primitive_vectors.push_back(Vec3d(0, 2, 0));
    primitive_vectors.push_back(Vec3d(0, 0, 2));
  }
  else if (strchr("fcba", *(centering.c_str()))) {
    primitive_vectors.push_back(Vec3d(len, len, 0));
    primitive_vectors.push_back(Vec3d(len, 0, len));
    primitive_vectors.push_back(Vec3d(0, len, len));
  }
  else if (centering == "i") {
    primitive_vectors.push_back(Vec3d(len, len, -len));
    primitive_vectors.push_back(Vec3d(len, -len, len));
    primitive_vectors.push_back(Vec3d(-len, len, len));
  }
}

double tetrahedral_volume(const Vec3d &a0, const Vec3d &a1, const Vec3d &a2,
                          const Vec3d &a3)
{
  return -vtriple(a1 - a0, a2 - a0, a3 - a0) / factorial(3);
}

void bravais_dual(Geometry &geom, const vector<Vec3d> &primitive_vectors,
                  const vector<int> &prim_vec_idxs, const vector<double> &vecs,
                  const vector<double> &angles, const double eps)
{
  Geometry tgeom;

  if (primitive_vectors.size()) {
    for (unsigned int i = 0; i < 4; i++)
      tgeom.add_vert(primitive_vectors[i]);
  }
  else {
    const vector<Vec3d> &verts = geom.verts();
    for (unsigned int i = 0; i < 4; i++) {
      if (prim_vec_idxs[i] >= (int)verts.size()) {
        fprintf(stderr,
                "error in bravais_dual: vertex %d does not exist. "
                "cannot perform dual\n",
                prim_vec_idxs[i]);
        return;
      }
      tgeom.add_vert(verts[prim_vec_idxs[i]]);
    }

    // for(unsigned int i=0; i<4; i++) {
    //   tgeom.add_vert(verts[prim_vec_idxs[i]]);
    //}
  }

  // changing orientation
  // if ( vecs[0] != vecs[1] )
  //   swap( vecs[0], vecs[1] );
  // vecs[2] = sqrt( vecs[0] * vecs[1] * sin(angles[2]*M_PI/180) );

  // do same scale and warp as happened to lattice
  bravais_scale(tgeom, vecs, false);
  bravais_warp(tgeom, angles, false);

  const vector<Vec3d> &tverts = tgeom.verts();
  Vec3d a0 = tverts[0]; // radial vector
  Vec3d a1 = tverts[1];
  Vec3d a2 = tverts[2];
  Vec3d a3 = tverts[3];

  double vecs_volume = tetrahedral_volume(a0, a1, a2, a3);
  if (double_eq(vecs_volume, 0, eps)) {
    fprintf(stderr, "error in bravais_dual: vectors are on one plane. cannot "
                    "perform dual\n");
    return;
  }

  double det = vtriple(a1, a2, a3);

  Vec3d b1 = vcross(a2, a3) / det;
  Vec3d b2 = vcross(a3, a1) / det;
  Vec3d b3 = vcross(a1, a2) / det;

  Trans3d m = Trans3d(b1, b2, b3).transpose();

  // if the dual lattice is desired to always have the spacing of the original
  // uncomment the code below and and at restore scale
  // store size to restore scale after dual
  // GeometryInfo info(geom);
  // double val = 1;
  // if (info.num_iedges() > 0)
  // double val = info.iedge_lengths().sum/info.num_iedges();

  // undo angles and vectors of the input geom
  bravais_warp(geom, angles, true);
  bravais_scale(geom, vecs, true);

  // apply dual transform
  geom.transform(m);

  // restore scale
  // geom.transform(Trans3d::scale(1/val));
}

void do_bravais(Geometry &geom, Geometry &container, brav_opts &opts)
{
  int strictness = 0;
  if (!bravais_check(opts.crystal_system, opts.centering, opts.vecs,
                     opts.angles, strictness, opts.verbose))
    return;

  // if hexagonal/cubic fill option, force struts and colors for cage
  bool struts = (opts.r_lattice_type == 2 || opts.r_lattice_type == 4)
                    ? true
                    : opts.cell_struts;
  Color vertex_col = (opts.r_lattice_type == 2 || opts.r_lattice_type == 4)
                         ? opts.vert_col[3]
                         : opts.vert_col[0];
  Color edge_col = (opts.r_lattice_type == 2 || opts.r_lattice_type == 4)
                       ? opts.edge_col[3]
                       : opts.edge_col[0];

  // for duals, force to primitive centering
  string centering = opts.dual_lattice ? "p" : opts.centering;
  // in case centering changed, check inversion
  if (centering != "f" && centering != "i")
    opts.inversion = false;
  bravais_cell(geom, centering, struts, vertex_col, edge_col, opts.inversion);

  // if hexagonal/cubic fill option, fill cell
  if (opts.r_lattice_type == 2 || opts.r_lattice_type == 4) {
    double one_third = 1.0 / 3.0;
    geom.add_vert(Vec3d(one_third, one_third, -one_third), opts.vert_col[0]);
    geom.add_vert(Vec3d(-one_third, -one_third, one_third), opts.vert_col[0]);
  }

  // correct for cell size for simple lattices
  if (opts.radius_by_coord.is_set()) {
    // if (opts.auto_grid_type == '8')
    //   opts.radius_by_coord *= sqrt(0.5);
    // else
    if (opts.centering == "p")
      opts.radius_by_coord *= 2.0;
  }

  // find automatic grid size
  if (opts.grid_for_radius) {
    double radius = opts.radius;
    if (opts.radius_default == 'k')
      radius = lattice_radius(container, opts.radius_default);
    else if (opts.radius_by_coord.is_set()) {
      vector<double> tangles(3, 90);
      radius = bravais_radius_by_coord(opts.radius_by_coord, opts.offset,
                                       opts.vecs, tangles);
    }

    radius = bravais_auto_grid_size(radius, opts.vecs, opts.angles);

    // grid big enough for radius
    if (radius != 0) {
      for (unsigned int i = 0; i < 3; i++)
        opts.grid.push_back((int)ceil(radius));
    }
  }

  // if grid size is still not set
  if (!opts.grid.size()) {
    // make grid big enough for most default spherical grids, when radius has
    // not been set.
    if (opts.container == 's') {
      int rep = (opts.centering == "p") ? 6 : 3;
      for (unsigned int i = 0; i < 3; i++)
        opts.grid.push_back(rep);
    }
    // default grid size
    else {
      for (unsigned int i = 0; i < 3; i++)
        opts.grid.push_back(2);
    }
  }

  if (opts.auto_grid_type)
    bravais_grid_type(opts.grid, opts.auto_grid_type);

  // original cell size is 2x2x2
  vector<double> cell_size(3, 2.0);

  geom_to_grid(geom, opts.grid, cell_size, opts.eps);

  if (opts.r_lattice_type == 1 || opts.r_lattice_type == 3) {
    r_lattice_overlay(geom, opts.grid, cell_size, opts.vert_col[3],
                      opts.edge_col[3], opts.eps);
    if (opts.r_lattice_type == 3)
      geom.transform(r_lattice_trans_mat(true));
  }

  if (opts.auto_grid_type == '8')
    bravais_eighth_cell_grid(geom);
  else
    // translate grid into the positive realm
    geom.transform(Trans3d::translate(Vec3d(1, 1, 1)));

  if (opts.verbose) {
    fprintf(stderr, "info: grid size is %d x %d x %d", opts.grid[0],
            opts.grid[1], opts.grid[2]);
    if (opts.auto_grid_type == '8')
      fprintf(stderr, " (clipped to: %d.5 x %d.5 x %d.5)", opts.grid[0] - 1,
              opts.grid[1] - 1, opts.grid[2] - 1);
    fprintf(stderr, "\n");
  }

  // if hexagonal/cubic fill option, position for true scale and warp
  if (opts.r_lattice_type == 2 || opts.r_lattice_type == 4) {
    vector<double> vecs(3);
    vecs[0] = sqrt(2);
    vecs[1] = sqrt(3);
    vecs[2] = sqrt(2);
    bravais_scale(geom, vecs, false);

    vector<double> angles(3);
    angles[0] = 90.0;
    angles[1] = 120.0;
    angles[2] = 90.0;
    bravais_warp(geom, angles, false);

    // rotate into place
    geom.transform(Trans3d::rotate(acos(1.0 / 3.0) / 2.0, 0, 0));
    geom.transform(Trans3d::rotate(0, 0, -45.0 * (M_PI / 180.0)));

    if (opts.cell_struts)
      add_color_struts(geom, 1.0, opts.edge_col[0], opts.eps);

    if (opts.r_lattice_type == 4)
      geom.transform(r_lattice_trans_mat(true));
  }

  bravais_scale(geom, opts.vecs, false);
  bravais_warp(geom, opts.angles, false);

  // save original center
  Vec3d original_center = centroid(geom.verts());

  // save lattice in case if adding back in end
  Geometry tgeom;
  if (opts.append_lattice)
    tgeom = geom;

  // if dual graph do that first
  if (opts.dual_lattice) {
    vector<Vec3d> primitive_vectors;
    if (opts.use_centering_for_dual)
      bravais_primitive_vectors(primitive_vectors, opts.centering);
    bravais_dual(geom, primitive_vectors, opts.prim_vec_idxs, opts.vecs,
                 opts.angles, opts.eps);
  }

  for (unsigned int i = 0; i < opts.strut_len.size(); i++)
    add_color_struts(geom, opts.strut_len[i] * opts.strut_len[i],
                     opts.edge_col[0], opts.eps);

  // radius calculation if needed
  if (opts.radius_by_coord.is_set())
    opts.radius = bravais_radius_by_coord(opts.radius_by_coord, opts.offset,
                                          opts.vecs, opts.angles);
  else if (!opts.radius)
    opts.radius =
        bravais_radius(opts.grid, opts.vecs, opts.angles, opts.radius_default);

  // scoop
  if (opts.cfile.length() > 0) {
    geom_container_clip(geom, container,
                        (opts.radius_default == 'k')
                            ? lattice_radius(container, opts.radius_default)
                            : opts.radius,
                        opts.offset, opts.verbose, opts.eps);
  }
  else if (opts.container == 's') {
    geom_spherical_clip(geom, opts.radius, opts.offset, opts.verbose, opts.eps);
  }

  if (opts.voronoi_cells) {
    Geometry vgeom;
    if (get_voronoi_geom(geom, vgeom, opts.voronoi_central_cell, false,
                         opts.eps)) {
      Coloring(&vgeom).vef_one_col(opts.vert_col[2], opts.edge_col[2],
                                   opts.face_col[2]);
      geom = vgeom;
    }
  }

  if (opts.convex_hull) {
    Status stat = (opts.add_hull) ? geom.add_hull() : geom.set_hull();
    if (stat.is_error())
      fprintf(stderr, "%s\n", stat.c_msg());
    else {
      geom.orient(1); // positive orientation
      if (opts.verbose)
        convex_hull_report(geom, opts.add_hull);
      Coloring(&geom).vef_one_col(opts.vert_col[1], opts.edge_col[1],
                                  opts.face_col[1]);
    }
  }

  if (opts.append_lattice) {
    geom.append(tgeom);
    tgeom.clear_all();
  }

  if (opts.color_method)
    color_by_symmetry_normals(geom, opts.color_method, opts.face_opacity,
                              opts.eps);

  if (opts.trans_to_origin)
    geom.transform(Trans3d::translate(-centroid(geom.verts())));

  if (opts.list_radii) {
    if (opts.list_radii_original_center) {
      opts.list_radii_center = original_center;
      if (opts.list_radii_original_center == 'o')
        opts.list_radii_center += opts.offset;
    }
    list_grid_radii(opts.ofile, geom, opts.list_radii_center, opts.list_radii,
                    opts.eps);
  }
  else if (opts.list_struts)
    list_grid_struts(opts.ofile, geom, opts.list_struts, opts.eps);

  // last so not to alter listing outcomes
  if (opts.cent_col.is_set())
    color_centroid(geom, opts.cent_col, opts.eps);

  if (opts.append_container) {
    container.add_missing_impl_edges();
    container.clear(FACES);
    geom.append(container);
  }

  // fprintf(stderr,"The volume of a single cell is %lf\n",
  // bravais_volume(opts.crystal_system, opts.vecs, opts.angles));
}

int main(int argc, char **argv)
{
  brav_opts opts;
  opts.process_command_line(argc, argv);
  bravais ubravais;

  if (opts.list_bravais) {
    ubravais.list_bravais();
    exit(0);
  }

  // read in container file if using. Check existance
  Geometry container;
  if (opts.cfile.length() > 0)
    opts.read_or_error(container, opts.cfile);

  // list lookups are based on centering type "C"
  string centering_temp;
  if (opts.centering == "a" || opts.centering == "b") {
    centering_temp = opts.centering;
    opts.centering = "c";
  }

  int sym_no = ubravais.lookup_sym_no(opts.crystal_system, opts.centering);

  if (sym_no < 0) {
    // restore "a" or "b" if it exists
    if (centering_temp.length())
      opts.centering = centering_temp;
    opts.error("unknown lattice '" + opts.crystal_system + " " +
               opts.centering + "'");
  }
  else if (sym_no >= ubravais.get_last_bravais())
    opts.error("lattice number '" + opts.crystal_system + "' out of range");
  else {
    opts.crystal_system = ubravais.get_crystal_system(sym_no);
    // restore "a" or "b" if it exists
    if (centering_temp.length())
      opts.centering = centering_temp;
    else
      opts.centering = ubravais.get_centering(sym_no);
  }

  // from here on out everything expects lower case
  transform(opts.crystal_system.begin(), opts.crystal_system.end(),
            opts.crystal_system.begin(), ::tolower);
  transform(opts.centering.begin(), opts.centering.end(),
            opts.centering.begin(), ::tolower);

  // inversion only allowed for f or i centering
  if (opts.inversion) {
    if (opts.centering != "f" && opts.centering != "i")
      opts.error("inversion can only be used with centering F or I", 'I');
    if (opts.cell_struts)
      opts.warning("-u has no effect with inversion", 'u');
  }

  // in place of ubravais.list_bravais(sym_no);
  string crystal_system_print_str = opts.crystal_system;
  if (opts.crystal_system == "trigonal")
    crystal_system_print_str = "trigonal(rhombohedral)";
  string centering_print_str = opts.centering;
  if (opts.centering == "a" || opts.centering == "b")
    centering_print_str = "c(" + opts.centering + ")";

  if (opts.verbose)
    fprintf(stderr, "%2d)\t%-s %s\n", sym_no + 1,
            crystal_system_print_str.c_str(), centering_print_str.c_str());

  if (opts.auto_grid_type == '8' && opts.centering == "p")
    opts.error(
        "grid type of 8 cannot be used with primitive (simple) bravais cells",
        'G');

  if (opts.r_lattice_type > 0 &&
      !(opts.crystal_system == "cubic" && opts.centering == "p") &&
      opts.crystal_system != "trigonal")
    opts.warning("hex relation is only proper on Cubic P or Trigonal", 'R');

  // if (opts.dual_lattice && (!strchr("pfi", *(opts.centering.c_str()))))
  //   opts.error("dual lattices are only only allowed for centering types p, f
  //   or i", 'd');

  Geometry geom;
  do_bravais(geom, container, opts);

  if (!opts.list_radii && !opts.list_struts)
    opts.write_or_error(geom, opts.ofile);

  return 0;
}
