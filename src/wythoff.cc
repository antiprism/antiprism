/*
   Copyright (c) 2003-2017, Adrian Rossiter

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
   Name: wythoff.cc
   Description: make tilings/polyhedra with a generalised Wythoff construction
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <functional>
#include <memory>
#include <string>
#include <vector>

using std::string;
using std::vector;

using namespace anti;

class wy_opts : public ProgramOpts {
private:
public:
  bool input_is_meta = false;
  bool add_meta = false;
  double face_ht = 0.0;
  string relabel;
  bool reverse = false;
  Tiling tiling;
  TilingColoring col_type;
  Coloring clrngs[3];
  bool quiet = false;
  string ifile;
  string ofile;

  wy_opts() : ProgramOpts("wythoff") {}

  void process_command_line(int argc, char **argv);
  void usage();
};

void wy_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] [input_file]

Read a file in OFF format and apply a specified pattern to generate polygon
tiles. The polyhedron faces are divided by a 'meta' operation into triangles
each having vertices which are a vertex V, edge centre E and face centre F.
A start point is positioned on one of these triangles, the next point is
found by using the pattern to step between triangles, leading to a circuit.
If input_file is not given the program reads from standard input.

Options
%s
  -p <pat>  pattern in form: [Coords0:Coords1:...]Path0,Path1...
            Coordinates are barycentric, in form aVbEcF:
              VEF element letters, and a,b,c are barycentric coordinates
              corresponding to the following element letter. Ommiting a
              an element letter and coordinate sets the coordinate to zero.
              Ommitting just the coordinate sets the coordinate to 1. E.g
            V = (1,0,0), VE = (1,1,0), V2E3F = (1,2,3)
            Paths are in the form: TrisPidx0Move0Pidx1Move1...
              Tris: one of +-* (default +) indicating that paths should
                start for positive, negative or both kinds of triangles.
              Pidx: an index number of a point from the coordinates list
              Move: an operation for stepping to the next triangle, given
                as a series of characters from the following:
                  _     - no move, stay on the same triangle
                  v,e,f - step over side opposite V,E,F
                  V,E,F - step two trianglesi, rotating about V,E,F,
                          according to: V=ef, E=fv, F=ve
              Paths can start or end with either a move or a point
  -c <op>   Conway polyhedron notation operator, or 'list' to list all
            available operators with their corresponding patterns
  -R        reverse pattern, exchanges the signs of the start triangles
  -r <elms> relabel pattern, exactly three letters VEF written in any order
            e.g. EFV relabels the pattern as V->E,v->e,E->F,e->f,F->V,f->v
  -M        input geometry is a 'meta' tiling, don't apply meta operation
%s
  -m <maps> a comma separated list of colour maps used to transform colour
            indexes (default: spread), a part consisting of letters from
            v, e, f, selects the element types to apply the map list to
  -u        output only one example of each type of tile (one per path)
  -a        add the 'meta'-transformed base
  -f <ht>   lift the face centres by this height
  -q        quiet, don't print report
  -o <file> write output to file (default: write to standard output)

)",
          prog_name(), help_ver_text,
          TilingColoring::get_option_help('C').c_str());
}

void wy_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":ho:p:c:f:Rr:MC:m:uqa")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'p':
      print_status_or_exit(tiling.read_pattern(optarg), c);
      break;

    case 'c':
      if (string(optarg) == "list") {
        tiling.print_conway_list();
        exit(0);
      }
      print_status_or_exit(tiling.read_conway(optarg), c);
      break;

    case 'R':
      reverse = true;
      break;

    case 'r':
      relabel = optarg;
      break;

    case 'a':
      add_meta = true;
      break;

    case 'M':
      input_is_meta = true;
      break;

    case 'u':
      tiling.set_one_of_each_tile(true);
      break;

    case 'f':
      print_status_or_exit(read_double(optarg, &face_ht), c);
      break;

    case 'C':
      print_status_or_exit(col_type.read_coloring(optarg), c);
      break;

    case 'm':
      print_status_or_exit(read_colorings(clrngs, optarg), c);
      break;

    case 'q':
      quiet = true;
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (argc - optind > 1) {
    error("too many arguments");
    exit(1);
  }

  if (argc - optind == 1)
    ifile = argv[optind];
}

void color_meta(Geometry &meta)
{
  const Color light(1.0, 0.8, 0.6, 0.5);
  const Color dark(0.1, 0.3, 0.6, 0.5);
  for (int i = 0; i < (int)meta.faces().size(); i++)
    meta.colors(FACES).set(i, (i % 2) ? light : dark);
  meta.add_missing_impl_edges();
  Coloring clrng(&meta);
  clrng.e_one_col(Color::invisible);
  clrng.v_one_col(Color::invisible);
}

int main(int argc, char *argv[])
{
  wy_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);

  GeometryInfo info(geom);
  bool orientable = true;
  if (!info.is_orientable()) {
    orientable = false;
    opts.warning("base polyhedron is not orientable: tiles will start from "
                 "every meta triangle");
  }
  if (!opts.input_is_meta && orientable && !info.is_oriented()) {
    opts.warning("base polyhedron is not oriented: it will be oriented.");
    geom.orient(1); // positive orientation
  }

  Tiling &tiling = opts.tiling;
  if (opts.relabel != "")
    opts.print_status_or_exit(tiling.relabel_pattern(opts.relabel), 'r');
  if (opts.reverse) {
    tiling.reverse_pattern();
    opts.warning("base polyhedron is not oriented: reverse has no effect", 'R');
  }

  if (!orientable)
    tiling.start_everywhere();

  Status stat = tiling.set_geom(geom, opts.input_is_meta, opts.face_ht);
  if (stat.is_error())
    opts.print_status_or_exit(stat, 'm');
  tiling.set_coloring(opts.col_type);
  Geometry ogeom;
  vector<Tile::TileReport> tile_reports;
  opts.print_status_or_exit(tiling.make_tiling(ogeom, &tile_reports));
  if (!orientable)
    merge_coincident_elements(ogeom, "f");

  if (!opts.quiet) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Tiling pattern: %s\n", tiling.pattern_string().c_str());
    fprintf(stderr, "Tile Counts:\n");
    fprintf(stderr, "No.:  Tiles:  Type:  Full Association:\n");
    string assoc_elem_str = "VEF345X1"; // VEF=6 -> X
    for (unsigned int i = 0; i < tile_reports.size(); i++) {
      const auto &rep = tile_reports[i];
      string full_assoc = rep.step.size() ? "(" + rep.step + ")" : "";
      full_assoc += rep.assoc.size() ? rep.assoc : "1";
      full_assoc += rep.step_back.size() ? "(" + rep.step_back + ")" : "";
      fprintf(stderr, "%3u, %6d, %5c,  %s\n", i, rep.count,
              assoc_elem_str[rep.assoc_type], full_assoc.c_str());
    }
    fprintf(stderr, "\n");
  }

  // colour the final model
  auto &clrngs = opts.clrngs;
  for (int i = 0; i < 3; i++) {
    auto &clrng = clrngs[i];
    clrng.set_geom(&ogeom);
    if (!clrng.get_cmaps().size()) { // no colouring set for this element
      if (i == VERTS)
        clrng.add_cmap(tiling.get_default_point_colormap());
      else // FACES and EDGES
        clrng.add_cmap(tiling.get_default_tile_colormap());
    }
  }

  clrngs[FACES].f_apply_cmap();
  clrngs[EDGES].e_apply_cmap();
  clrngs[VERTS].v_apply_cmap();

  if (opts.add_meta) {
    Geometry meta = tiling.get_meta();
    color_meta(meta);
    ogeom.append(meta);
  }

  opts.write_or_error(ogeom, opts.ofile);

  return 0;
}
