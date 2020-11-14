/*
   Copyright (c) 2008-2016, Roger Kaufman

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
#include "lattice_grid.h"

#include <cctype>
#include <climits>
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
  string rfile;

  vector<double> strut_len;
  bool strip_faces;
  char color_edges_by_sqrt;
  char container;
  bool append_container;
  double radius;
  char radius_default;
  Vec3d offset;
  bool voronoi_cells;
  bool voronoi_central_cell;
  bool convex_hull;
  bool add_hull;
  bool append_lattice;
  char color_method;
  int face_opacity;
  Vec3d list_radii_center;
  char list_radii_original_center;
  int list_radii;
  int list_struts;
  Color cent_col;
  bool trans_to_origin;
  string R_merge_elems;
  int blend_type;

  bool verbose;
  double epsilon;

  // 0 - lattice  1 - convex hull  2 - voronoi
  vector<Color> vert_col;
  vector<Color> edge_col;
  vector<Color> face_col;

  lutil_opts()
      : ProgramOpts("lat_util"), strip_faces(true), color_edges_by_sqrt('\0'),
        container('c'), append_container(false), radius(0), radius_default('s'),
        voronoi_cells(false), voronoi_central_cell(false), convex_hull(false),
        add_hull(false), append_lattice(false), color_method('\0'),
        face_opacity(-1), list_radii_original_center('\0'), list_radii(0),
        list_struts(0), trans_to_origin(false), R_merge_elems("vef"),
        blend_type(3), verbose(false), epsilon(0)
  {
  }

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
  -h        this help message
  -I        verbose output
  -z        suppress stripping of faces
  -c <type> container, c - cube (default), s - sphere (uses radius)
  -k <file> container, convex hull of off file or built in model (uses radius)
  -r <c,n>  radius. c is radius taken to optional root n. n = 2 is sqrt
               or  l - max insphere radius, s - min insphere radius (default)
               or  k - take radius from container specified by -k
  -q <vecs> center offset, in form \a_val,b_val,c_val\ (default: none)
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
  -l <lim>  minimum distance for unique vertex locations as negative exponent
               (default: %d giving %.0e)
  -o <file> write output to file (default: write to standard output)

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
               Note: input elements from -R input are not changed
  -E <col>  edge color (same format as for vertices) or
               key word: r,R for color edges by root value of final product
               lower case outputs map indexes. upper case outputs color values
  -F <col>  face color (same format as for vertices) or
               key word: s,S color by symmetry using face normals
               key word: c,C color by symmetry using face normals (chiral)
  -T <tran> face transparency for color by symmetry. valid range from 0 to 255

)",
          prog_name(), help_ver_text, int(-log(::epsilon) / log(10) + 0.5),
          ::epsilon);
}

void lutil_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;
  char errmsg[MSG_SZ] = {0};

  int sig_compare = INT_MAX;
  vector<double> double_parms;
  Split parts;
  Color col_tmp;

  vert_col.resize(4);
  edge_col.resize(4);
  face_col.resize(4);

  // set some default colors
  // 0 - lattice  1 - convex hull  2 - voronoi
  // voronoi cells  vcol = gold; ecol = lightgrey; fcol = transparent yellow
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
      print_status_or_exit(offset.read(optarg), c);
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
      if (!strcasecmp(parts[0], "r"))
        color_edges_by_sqrt = parts[0][0];
      else
        parse_color_string(this, optarg, c, "lcv", edge_col);
      break;

    case 'F':
      parts.init(optarg, ",");
      if ((!strcasecmp(parts[0], "s")) || (!strcasecmp(parts[0], "c")))
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
        if (strspn(parts[1], "vefan") != strlen(parts[1])) {
          strcpy_msg(errmsg,
                     msg_str("elements to merge are %s must be v, e, f, a or n",
                             parts[1])
                         .c_str());
          error(errmsg, c);
        }
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
        print_status_or_exit(list_radii_center.read(optarg), c);
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
      print_status_or_exit(read_int(optarg, &sig_compare), c);
      if (sig_compare < 0) {
        warning("limit is negative, and so ignored", c);
      }
      if (sig_compare > DEF_SIG_DGTS) {
        warning("limit is very small, may not be attainable", c);
      }
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
    error("-c and -k cannot be specified together");

  if (append_container && !cfile.length())
    error("container can only be appended if one is provided with -k", 'K');

  if (radius_default == 'k' && !cfile.length())
    error("-r k can only be used if -k container is specified", 'r');

  if (container == 'c' && !cfile.length() && (radius != 0 || offset.is_set()))
    warning("cubic container in use. Radius and offset ignored");

  // if ((container == 's' || (cfile.length() && !radius_default)) && !radius)
  //   error("radius not set");

  if (list_radii && list_struts)
    error("cannot list radii and struts at the same time");

  epsilon = (sig_compare != INT_MAX) ? pow(10, -sig_compare) : ::epsilon;
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
                     opts.edge_col[0]);

  if (!opts.radius && opts.radius_default != 'k')
    opts.radius = lattice_radius(geom, opts.radius_default);

  // scoop
  if (opts.cfile.length())
    geom_container_clip(geom, container,
                        (opts.radius_default == 'k')
                            ? lattice_radius(container, opts.radius_default)
                            : opts.radius,
                        opts.offset, opts.verbose, opts.epsilon);
  else if (opts.container == 's')
    geom_spherical_clip(geom, opts.radius, opts.offset, opts.verbose,
                        opts.epsilon);

  if (opts.voronoi_cells) {
    Geometry vgeom;
    if (get_voronoi_geom(geom, vgeom, opts.voronoi_central_cell, false,
                         opts.epsilon)) {
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
      geom.orient();
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
                              opts.epsilon);

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
                                opts.epsilon);
  }

  if (opts.list_radii) {
    if (opts.list_radii_original_center) {
      opts.list_radii_center = original_center;
      if (opts.list_radii_original_center == 'o')
        opts.list_radii_center += opts.offset;
    }
    list_grid_radii(opts.ofile, geom, opts.list_radii_center, opts.list_radii,
                    opts.epsilon);
  }
  else if (opts.list_struts)
    list_grid_struts(opts.ofile, geom, opts.list_struts, opts.epsilon);

  if (opts.append_container) {
    container.add_missing_impl_edges();
    container.clear(FACES);
    geom.append(container);
  }

  if (opts.cent_col.is_set())
    color_centroid(geom, opts.cent_col, opts.epsilon);

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
