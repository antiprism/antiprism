/*
   Copyright (c) 2011-2016, Roger Kaufman

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
   Name: off_normals.cc
   Description: Display normals of faces, implicit edges, or vertices
   Project: Antiprism - http://www.antiprism.com
*/

#include <cstdio>
#include <cstdlib>

#include <cctype>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "../base/antiprism.h"

using std::pair;
using std::string;
using std::swap;
using std::vector;

using namespace anti;

class off_normals_opts : public ProgramOpts {
public:
  string ifile;
  string ofile;

  bool unit_normals;
  char exclude_normals_elems;
  char force_normals_polarity;
  bool elem_normal_vecs;
  char edge_normal_method;
  char base_normal_method;
  string show_pointing;
  string show_elems;
  string average_pattern;
  bool alternate_calculation;

  Color outward_normal_col;
  Color inward_normal_col;
  Color hemispherical_normal_col;
  Color edge_normal_col;
  Color base_normal_col;

  Vec3d center;

  int sig_compare;
  double epsilon;

  off_normals_opts()
      : ProgramOpts("off_normals"), unit_normals(false),
        exclude_normals_elems('\0'), force_normals_polarity('\0'),
        elem_normal_vecs(false), edge_normal_method('\0'),
        base_normal_method('b'), show_pointing("oih"), show_elems("f"),
        average_pattern("r"), alternate_calculation(false),
        hemispherical_normal_col(Color(127, 127, 127)), sig_compare(INT_MAX),
        epsilon(0)
  {
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

void off_normals_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] [input_file]

Display normals of faces, implicit edges, and vertices
If input_file is not given the program reads from standard input.

Options
%s
  -u        unit normals  (positional normals otherwise)
  -e        connect to element centroid
  -p <opt>  force polarity. o - set all outward,  i - set all inward
               r - reverse both inward and outward
  -i <elms> include normals. The element string can include o, i and h
               to show, respectively, outward, inward and hemispherical
               note: exclusion occurs before -p  (default: oih)
  -s <elms> include elements. The element string can include v, e and f
               to show, respectively, vertices, edges and faces  (default: f)
  -d <opt>  delete elements.  f - delete faces of excluded normals
               a - delete all of original model
  -c <opts> average pattern string for edge and vertex normals. Done before -p
               r - raw,  o - outward,  i - inward,  u - unit  (default: r)
  -a        alternate calculation for vertex normals
  -C <xyz>  center of model, in form 'X,Y,Z'  (default: centroid)
  -l <lim>  minimum distance for unique vertex locations as negative exponent
               (default: %d giving %.0e)
  -o <file> write output to file  (default: write to standard output)

Coloring Options (run 'off_util -H color' for help on color formats)
  -O <col>  outward normal vertex color
  -I <col>  inward normal vertex color
               default: vertex color is negative of outward col
  -H <col>  hemispherical normal vertex color  (default: gray50)
  -E <col>  normal vector color. connected to element centroid
               default: color of normal vertex
               key word: r take random color
  -B <col>  normal vector base color. color at element centroid
               key word: b take color of element (default)
               key word: n take color of normal vertex
)",
prog_name(), help_ver_text, int(-log(::epsilon)/log(10) + 0.5),::epsilon);
}

void off_normals_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":huep:i:s:d:c:aO:I:H:E:B:C:l:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'u':
      unit_normals = true;
      break;

    case 'e':
      elem_normal_vecs = true;
      break;

    case 'p':
      if (strlen(optarg) > 1 || !strchr("oir", *optarg))
        error(msg_str("force polarity is '%s', must be o, i, or r", optarg), c);
      force_normals_polarity = *optarg;
      break;

    case 'i':
      if (strspn(optarg, "oih") != strlen(optarg))
        error(msg_str("pointing to include are '%s', must be from o, i, and h",
                      optarg),
              c);
      show_pointing = optarg;
      break;

    case 's':
      if (strspn(optarg, "vef") != strlen(optarg))
        error(msg_str("elements to hide are '%s', must be from v, e, and f",
                      optarg),
              c);
      show_elems = optarg;
      break;

    case 'd':
      if (strlen(optarg) > 1 || !strchr("fa", *optarg))
        error(msg_str("delete elements is '%s', must be f or a", optarg), c);
      exclude_normals_elems = *optarg;
      break;

    case 'c':
      if (strspn(optarg, "roiu") != strlen(optarg))
        error(msg_str("average string is '%s', must be from r, o, i, and u",
                      optarg),
              c);
      average_pattern = optarg;
      break;

    case 'a':
      alternate_calculation = true;
      break;

    case 'O':
      print_status_or_exit(outward_normal_col.read(optarg), c);
      break;

    case 'I':
      print_status_or_exit(inward_normal_col.read(optarg), c);
      break;

    case 'H':
      print_status_or_exit(hemispherical_normal_col.read(optarg), c);
      break;

    case 'E':
      if (strchr("r", *optarg))
        edge_normal_method = *optarg;
      else
        print_status_or_exit(edge_normal_col.read(optarg), c);
      break;

    case 'B':
      if (strchr("bn", *optarg))
        base_normal_method = *optarg;
      else
        print_status_or_exit(base_normal_col.read(optarg), c);
      break;

    case 'C':
      print_status_or_exit(center.read(optarg), c);
      break;

    case 'l':
      print_status_or_exit(read_int(optarg, &sig_compare), c);
      if (sig_compare < 0) {
        warning("limit is negative, and so ignored", c);
      }
      if (sig_compare > 16) {
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

  if (argc - optind > 1)
    error("too many arguments");

  if (argc - optind == 1)
    ifile = argv[optind];

  epsilon = (sig_compare != INT_MAX) ? pow(10, -sig_compare) : ::epsilon;
}

void add_normals(Geometry &geom, const off_normals_opts &opts)
{
  Color outward_col = opts.outward_normal_col;
  Color inward_col = opts.inward_normal_col;
  if (!inward_col.is_set() && opts.outward_normal_col.is_set()) {
    inward_col = opts.outward_normal_col;
    inward_col.set_complement();
  }

  Color col;

  Geometry ngeom;
  vector<int> deleted_faces;

  FaceNormals x_normals(geom, opts.center, opts.epsilon);

  if (strchr(opts.show_elems.c_str(), 'f')) {
    std::unique_ptr<ColorMap> cmap(colormap_from_name("rnd")); // for -E r
    for (unsigned int i = 0; i < x_normals.size(); i++) {
      Normal x_normal = x_normals[i];

      bool plotted = false;

      Vec3d normal;
      if (x_normal.is_hemispherical()) {
        if (strchr(opts.show_pointing.c_str(), 'h')) {
          plotted = true;

          normal = x_normal.raw();
          col = opts.hemispherical_normal_col;
        }
      }
      else {
        if (x_normal.is_inward()) { // normal points inward
          if (strchr(opts.show_pointing.c_str(), 'i')) {
            plotted = true;

            if ((opts.force_normals_polarity == 'o') ||
                (opts.force_normals_polarity == 'r')) { // force it outwards
              normal = x_normal.outward();
              col = outward_col;
            }
            else {
              normal = x_normal.inward();
              col = inward_col;
            }
          }
        }

        if (x_normal.is_outward()) { // normal points outward
          if (strchr(opts.show_pointing.c_str(), 'o')) {
            plotted = true;

            if ((opts.force_normals_polarity == 'i') ||
                (opts.force_normals_polarity == 'r')) { // force it inwards
              normal = x_normal.inward();
              col = inward_col;
            }
            else {
              normal = x_normal.outward();
              col = outward_col;
            }
          }
        }
      }

      if (!plotted)
        deleted_faces.push_back(i); // if deleting faces of unplotted normals
      else {
        if (!normal.is_set())
          normal = x_normal.raw();

        if (opts.unit_normals)
          normal = normal.unit();

        ngeom.add_vert(normal, col);

        if (opts.elem_normal_vecs) {
          // get base color
          Color bcol = (opts.base_normal_col.is_set())
                           ? opts.base_normal_col
                           : ((opts.base_normal_method == 'b')
                                  ? geom.colors(FACES).get((int)i)
                                  : col);
          // add point at centroid
          ngeom.add_vert(geom.face_cent(i), bcol);

          // get edge color
          Color ecol;
          int sz = ngeom.verts().size();
          if (opts.edge_normal_method == 'r')
            ecol = cmap->get_col(sz - 1);
          else
            ecol = (opts.edge_normal_col.is_set()) ? opts.edge_normal_col : col;
          // edge from face centroid to normal
          ngeom.add_edge(make_edge(sz - 1, sz - 2), ecol);
        }
      }
    }
  }

  if (strchr(opts.show_elems.c_str(), 'e')) {
    vector<vector<int>> implicit_edges;
    geom.get_impl_edges(implicit_edges);

    for (auto edge : implicit_edges) {
      Normal x_normal =
          x_normals.edge_normal(edge[0], edge[1], opts.average_pattern);

      if (x_normal.is_hemispherical()) {
        if (!strchr(opts.show_pointing.c_str(), 'h')) {
          continue;
        }
      }

      Vec3d normal;
      if (x_normal.is_inward()) { // normal points inward
        if (!strchr(opts.show_pointing.c_str(), 'i')) {
          continue;
        }
        if ((opts.force_normals_polarity == 'o') ||
            (opts.force_normals_polarity == 'r')) { // force it outwards
          normal = x_normal.outward();
          col = outward_col;
        }
        else {
          normal = x_normal.inward();
          col = inward_col;
        }
      }

      if (x_normal.is_outward()) { // normal points outward
        if (!strchr(opts.show_pointing.c_str(), 'o')) {
          continue;
        }
        if ((opts.force_normals_polarity == 'i') ||
            (opts.force_normals_polarity == 'r')) { // force it inwards
          normal = x_normal.inward();
          col = inward_col;
        }
        else {
          normal = x_normal.outward();
          col = outward_col;
        }
      }

      if (!normal.is_set())
        normal = x_normal.raw();

      if (opts.unit_normals)
        normal = normal.unit();

      ngeom.add_vert(normal, col);

      if (opts.elem_normal_vecs) {
        // get base color for edge. might not explicitly be colored
        int e_idx = find_edge_in_edge_list(geom.edges(), edge);
        Color expl_col = (e_idx > -1) ? geom.colors(EDGES).get(e_idx) : Color();

        Color bcol = (opts.base_normal_col.is_set())
                         ? opts.base_normal_col
                         : ((opts.base_normal_method == 'b') ? expl_col : col);
        // add point at centroid
        ngeom.add_vert(centroid(geom.verts(), edge), bcol);
        // get edge color
        Color ecol =
            (opts.edge_normal_col.is_set()) ? opts.edge_normal_col : col;
        // edge from face centroid to normal
        ngeom.add_edge(
            make_edge(ngeom.verts().size() - 1, ngeom.verts().size() - 2),
            ecol);
      }
    }
  }

  if (strchr(opts.show_elems.c_str(), 'v')) {
    const vector<Vec3d> &verts = geom.verts();

    // in case this is needed for alternate vertex normal calculation
    GeometryInfo info(geom);
    const vector<Vec3d> &v_norms = info.get_vert_norms();

    for (unsigned int i = 0; i < verts.size(); i++) {
      Normal x_normal;
      if (opts.alternate_calculation)
        x_normal = Normal(geom, v_norms[i], i, opts.center, opts.epsilon);
      else
        x_normal = x_normals.vert_normal(i, opts.average_pattern);

      if (x_normal.is_hemispherical()) {
        if (!strchr(opts.show_pointing.c_str(), 'h')) {
          continue;
        }
      }

      Vec3d normal;
      if (x_normal.is_inward()) { // normal points inward
        if (!strchr(opts.show_pointing.c_str(), 'i')) {
          continue;
        }
        if ((opts.force_normals_polarity == 'o') ||
            (opts.force_normals_polarity == 'r')) { // force it outwards
          normal = x_normal.outward();
          col = outward_col;
        }
        else {
          normal = x_normal.inward();
          col = inward_col;
        }
      }

      if (x_normal.is_outward()) { // normal points outward
        if (!strchr(opts.show_pointing.c_str(), 'o')) {
          continue;
        }
        if ((opts.force_normals_polarity == 'i') ||
            (opts.force_normals_polarity == 'r')) { // force it inwards
          normal = x_normal.inward();
          col = inward_col;
        }
        else {
          normal = x_normal.outward();
          col = outward_col;
        }
      }

      if (!normal.is_set())
        normal = x_normal.raw();

      if (opts.unit_normals)
        normal = normal.unit();

      ngeom.add_vert(normal, col);

      if (opts.elem_normal_vecs) {
        // get base color
        Color bcol = (opts.base_normal_col.is_set())
                         ? opts.base_normal_col
                         : ((opts.base_normal_method == 'b')
                                ? geom.colors(VERTS).get((int)i)
                                : col);
        // add point at centroid
        ngeom.add_vert(verts[i], bcol);
        // get edge color
        Color ecol =
            (opts.edge_normal_col.is_set()) ? opts.edge_normal_col : col;
        // edge from face centroid to normal
        ngeom.add_edge(
            make_edge(ngeom.verts().size() - 1, ngeom.verts().size() - 2),
            ecol);
      }
    }
  }

  if (opts.exclude_normals_elems == 'f')
    geom.del(FACES, deleted_faces);
  else if (opts.exclude_normals_elems == 'a')
    geom.clear_all();
  geom.append(ngeom);
}

/*
Vec3d line_nearest_point(Vec3d P, Vec3d A, Vec3d B)
{
   Vec3d v1 = P-A;
   Vec3d v2 = B-A;

   double v2m = v2.mag2();

   double D = vdot(v1,v2);

   double dist = D/v2m;

   return (A+(v2*dist));
}
*/

int main(int argc, char *argv[])
{
  off_normals_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);

  add_normals(geom, opts);

  opts.write_or_error(geom, opts.ofile);

  return 0;
}
