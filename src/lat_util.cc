/*
   Copyright (c) 2008-2023, Roger Kaufman

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
   Name: lat_util.cc
   Description: Do various things to lattices
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"
#include "color_common.h"
#include "lat_util_common.h"

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using std::string;
using std::vector;

using namespace anti;

class lutil_opts : public ProgramOpts {
public:
  vector<string> ifiles;
  string ofile;
  string cfile;
  string rfile; // geom to be repeated

  vector<double> strut_len;          // strut lengths to create if found
  double radius = 0;                 // radius of container
  char radius_default = 's';         // insphere radius is default
  Vec3d offset;                      // offset from origin for container
  char container = 'c';              // default container (c all cubic space)
  bool append_container = false;     // append cage of -k container
  bool voronoi_cells = false;        // calculate voronoi cells
  bool voronoi_central_cell = false; // include voronoi cells only at center
  bool convex_hull = false;          // convex hull for waterman polyhedra
  bool add_hull = false;             // add lattice to waterman polyhedra
  bool append_lattice = false;       // append lattice to final produc
  bool trans_to_origin = false;      // tranlate lattice centroid to origin
  bool verbose = false;              // option W gives lattice information

  bool strip_faces = true;      // strip faces off input model
  string R_merge_elems = "vef"; // merge for repeat elements
  int blend_type = 3;           // RGB blend for repeat elements

  Vec3d list_radii_center;                // centroid for calculating radii list
  char list_radii_original_center = '\0'; // c - use centroid o - use q offset
  int list_radii = 0;                     // list radii
  int list_struts = 0;                    // list struts

  double eps = anti::epsilon;

  char color_method = '\0';        // color method for color by symmetry
  int face_opacity = -1;           // transparency from 0 to 255
  Color cent_col;                  // color a centroid, add if missing
  char color_edges_by_sqrt = '\0'; // will be r or R, color edges base on sqrt

  // colors for different parts are held in vector
  // 0 - lattice  1 - convex hull  2 - voronoi
  vector<Color> vert_col;
  vector<Color> edge_col;
  vector<Color> face_col;

  lutil_opts() : ProgramOpts("lat_util") {}

  void process_command_line(int argc, char **argv);
  void usage();
};

void lutil_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] [input_files]

Read one or more files in OFF format, combine them into a single file and
process it. Operations take place in the order listed below. input_files is the
list of files to process. If they don't exist, explicit edges are created.
If the input possesses faces they are stripped by default.

Options
%s
  -I        verbose output
  -l <lim>  minimum distance for unique vertex locations as negative exponent
              (default: %d giving %.0e)
  -o <file> write output to file (default: write to standard output)

Scene Options
  -z        suppress stripping of faces
  -c <type> container, c - cube (default), s - sphere (uses radius)
  -k <file> container, convex hull of off file or built in model (uses radius)
  -r <c,n>  radius. c is radius taken to optional root n. n = 2 is sqrt
              or  l - max insphere radius, s - min insphere radius (default)
              or  k - take radius from container specified by -k
  -q <xyz>  center offset, three comma separated coordinates, 0 for origin
  -s <s,n>  create struts. s is strut length taken to optional root n
              use multiple -s parameters for multiple struts
  -D <opt>  Voronoi (a.k.a Dirichlet) cells (Brillouin zones for duals)
              c - cells only, i - cell(s) touching center only
  -C <opt>  c - convex hull only, i - keep interior
  -A        append the original lattice to the final product
  -R <fi,s> repeat off file fi at every vertex in lattice. If optional s is
            set, sort and merge elements whose coordinates are the same to
            the number of decimal places given by option -l.  elements can
            include: v - vertices, e - edges, f - faces,  a - all (vef)
            n - no merging  (default 'a'. Colors blended as RGB)
  -K        append cage of container of -k to final product
  -Z <col>  add center vertex to final product in color col
  -O        translate center of final product to origin

Listing Options
  -Q <vecs> center for radius calculations in -L (default: centroid)
              c - original center, o - original center + offset in -q
  -L <opt>  list unique radial distances of points (to standard output)
              f - full report, v - values only
  -S <opt>  list every possible strut value (to standard output)
              f - full report, v - values only

Coloring Options (run 'off_util -H color' for help on color formats)
  -V <col>  vertex color, (optional) transparency, (optional) elements
              transparency: valid range from 0 (invisible) to 255 (opaque)
              elements to color are l - lattice, c - convex hull, v - voronoi
                 (default elements: lcv)
              Note: input element colors from -R input are not changed
  -E <col>  edge color (same format as for vertices)
              lower case outputs map indexes. upper case outputs color values
              t,T - for color edges by root value of final product
  -F <col>  face color (same format as for vertices)
              special coloring: calculates color from normals, not maps
              lower case outputs map indexes. upper case outputs color values
              y,Y - color by symmetry using face normals
              z,Z - color by symmetry using face normals (chiral)
  -T <tran> face transparency for color by symmetry. valid range from 0 to 255

)",
          prog_name(), help_ver_text, int(-log(anti::epsilon) / log(10) + 0.5),
          anti::epsilon);
}

void lutil_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  vector<double> double_parms;
  Split parts;
  Color col_tmp;

  vert_col.resize(4);
  edge_col.resize(4);
  face_col.resize(4);

  // set some default colors
  // 0 - lattice  1 - convex hull  2 - voronoi
  // voronoi cells  vcol = gold; ecol = lightgray; fcol = transparent yellow
  vert_col[2] = Color(255, 215, 0);
  edge_col[2] = Color(211, 211, 211);
  face_col[2] = Color(255, 255, 0, 128);

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv,
                     ":hIzc:k:r:q:s:D:C:AV:E:F:T:Z:KOR:Q:L:S:l:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'I':
      verbose = true;
      break;

    case 'z':
      strip_faces = false;
      break;

    case 'c':
      if (strlen(optarg) > 1 || !strchr("cs", *optarg))
        error("method is '" + string(optarg) + "' must be c or s", c);
      container = *optarg;
      break;

    case 'k':
      cfile = optarg;
      break;

    case 'r':
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
          error("radius cannot be negative or zero", "s", c);
      }
      break;

    case 'q':
      print_status_or_exit(offset.read_maths(optarg), c);
      break;

    case 's': {
      double strut_tmp;
      print_status_or_exit(read_double_list(optarg, double_parms, 2), c);
      strut_tmp = double_parms[0];
      if (double_parms.size() == 2)
        strut_tmp = pow(strut_tmp, 1 / double_parms[1]);
      if (strut_tmp <= 0)
        error("strut lengths cannot be negative", "s");
      strut_len.push_back(strut_tmp);
      break;
    }

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
      parse_color_string(this, optarg, c, "lcv", vert_col);
      break;

    case 'E':
      parts.init(optarg, ",");
      if (strlen(parts[0]) == 1 && strchr("tT", parts[0][0]))
        color_edges_by_sqrt = parts[0][0];
      else
        parse_color_string(this, optarg, c, "lcv", edge_col);
      break;

    case 'F':
      parts.init(optarg, ",");
      if (strlen(parts[0]) == 1 && strchr("yYzZ", parts[0][0]))
        color_method = parts[0][0];
      else
        parse_color_string(this, optarg, c, "lcv", face_col);
      break;

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

    case 'R': {
      int parts_sz = parts.init(optarg, ",");
      if (parts_sz > 2)
        error("the argument has more than 2 parts", c);

      rfile = parts[0];
      if (parts_sz > 1) {
        if (strspn(parts[1], "vefan") != strlen(parts[1]))
          error(msg_str("elements to merge are %s must be v, e, f, a or n",
                        parts[1])
                    .c_str(),
                c);
        if (strchr(parts[1], 'a') && strlen(parts[1]) > 1) {
          error("a includes vef, and must be used alone", c);
        }
        if (strchr(parts[1], 'n') && strlen(parts[1]) > 1) {
          error("n must be used alone", c);
        }
        if (strspn(parts[1], "ef") && !strchr(parts[1], 'v')) {
          warning("without v, some orphan vertices may result", c);
        }
        R_merge_elems = parts[1];
        if (R_merge_elems == "a")
          R_merge_elems = "vef";
        else if (R_merge_elems == "n")
          R_merge_elems = "";
      }
      break;
    }

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

  while (argc - optind)
    ifiles.push_back(string(argv[optind++]));

  if (container == 's' && cfile.length())
    error("-c and -k cannot be specified together", 'c');

  if (append_container && !cfile.length())
    error("container can only be appended if one is provided with -k", 'K');

  if (radius_default == 'k' && !cfile.length())
    error("-r k can only be used if -k container is specified", 'r');

  if (container == 'c' && !cfile.length() && (radius != 0 || offset.is_set()))
    warning("cubic container in use. Radius and offset ignored");

  // if ((container == 's' || (cfile.length() && !radius_default)) && !radius)
  //   error("radius not set");

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

void make_skeleton(Geometry &geom, const bool strip_faces)
{
  geom.add_missing_impl_edges();
  if (strip_faces)
    geom.clear(FACES);
}

void process_lattices(Geometry &geom, Geometry &container,
                      const Geometry &repeater, lutil_opts &opts)
{
  // add explicit edges and remove faces if necessary
  make_skeleton(geom, opts.strip_faces);

  // save original center
  Vec3d original_center = centroid(geom.verts());

  // save lattice in case if adding back in end
  Geometry tgeom;
  if (opts.append_lattice)
    tgeom = geom;

  Coloring(&geom).vef_one_col(opts.vert_col[0], opts.edge_col[0],
                              opts.face_col[0]);

  for (unsigned int i = 0; i < opts.strut_len.size(); i++)
    add_color_struts(geom, opts.strut_len[i] * opts.strut_len[i],
                     opts.edge_col[0], opts.eps);

  if (!opts.radius && opts.radius_default != 'k')
    opts.radius = lattice_radius(geom, opts.radius_default);

  // scoop
  if (opts.cfile.length())
    geom_container_clip(geom, container,
                        (opts.radius_default == 'k')
                            ? lattice_radius(container, opts.radius_default)
                            : opts.radius,
                        opts.offset, opts.verbose, opts.eps);
  else if (opts.container == 's')
    geom_spherical_clip(geom, opts.radius, opts.offset, opts.verbose, opts.eps);

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
      geom.orient(1);
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

  // face color by symmetry normals
  if (opts.color_method)
    color_by_symmetry_normals(geom, opts.color_method, opts.face_opacity,
                              opts.eps);

  // if color by sqrt was used, override all edges of all structure
  if (opts.color_edges_by_sqrt)
    color_edges_by_sqrt(geom, opts.color_edges_by_sqrt);

  // place geom at every vertex in lattice
  if (opts.rfile.length()) {
    Geometry geom2;
    for (const auto &i : geom.verts()) {
      Geometry rep = repeater;
      rep.transform(Trans3d::translate(i));
      geom2.append(rep);
    }
    geom = geom2;
    if (opts.R_merge_elems != "")
      merge_coincident_elements(geom, opts.R_merge_elems, opts.blend_type,
                                opts.eps);
  }

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

  if (opts.append_container) {
    container.add_missing_impl_edges();
    container.clear(FACES);
    geom.append(container);
  }

  if (opts.cent_col.is_set())
    color_centroid(geom, opts.cent_col, opts.eps);

  if (opts.trans_to_origin)
    geom.transform(Trans3d::translate(-centroid(geom.verts())));
}

int main(int argc, char *argv[])
{
  lutil_opts opts;
  opts.process_command_line(argc, argv);
  if (!opts.ifiles.size())
    opts.ifiles.push_back("");

  // read in container file if using. Check existance
  Geometry container;
  if (opts.cfile.length())
    opts.read_or_error(container, opts.cfile);

  Geometry repeater;
  if (opts.rfile.length())
    opts.read_or_error(repeater, opts.rfile);

  Geometry geoms;
  for (unsigned int i = 0; i < opts.ifiles.size(); i++) {
    Geometry geom;
    opts.read_or_error(geom, opts.ifiles[i]);
    geoms.append(geom);
  }

  process_lattices(geoms, container, repeater, opts);

  if (!opts.list_radii && !opts.list_struts)
    opts.write_or_error(geoms, opts.ofile);

  return 0;
}
