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

#include <algorithm>
#include <functional>
#include <memory>
#include <cctype>
#include <math.h>
#include <string.h>
#include <string>
#include <vector>

#include "../base/antiprism.h"

using std::string;
using std::vector;

using namespace anti;

class wy_opts : public ProgramOpts {
private:
public:
  bool input_is_meta;
  bool add_meta;
  double face_ht;
  string relabel;
  Tiling tiling;
  bool color_with_value;
  bool quiet;
  string ifile;
  string ofile;

  wy_opts()
      : ProgramOpts("wythoff"), input_is_meta(false), add_meta(false),
      face_ht(0.0), color_with_value(true), quiet(false)
  {
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

// clang-format off
void wy_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Read a file in OFF format and apply a specified pattern to generate polygon\n"
"tiles. The polyhedron faces are divided by a 'meta' operation into triangles\n"
"each having vertices which are a vertex V, edge centre E and face centre F.\n"
"A start point is positioned on one of these triangles, the next point is\n"
"found by using the pattern to step between triangles, leading to a circuit.\n"
"If input_file is not given the program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"  -p <pat>  pattern in form: [Coords0:Coords1:...]Path0,Path1...\n"
"            Coordinates are barycentric, in form VaEbFc:\n"
"              VEF element letters, and a,b,c are barycentric coordinates\n"
"              corresponding to the preceding element letter. Ommiting a\n"
"              an element letter and coordinate sets the coordinate to zero.\n"
"              Ommitting just the coordinate sets the coordinate to 1. E.g\n"
"            V = (1,0,0), VE = (1,1,0), VE2F3 = (1,2,3)\n"
"            Paths are in the form: TrisPidx0Move0Pidx1Move1...\n"
"              Tris: one of +-* (default +) indicating that paths should\n"
"                start for positive, negative or both kinds of triangles.\n"
"              Pidx: an index number of a point from the coordinates list\n"
"              Move: an operation for stepping to the next triangle, given\n"
"                as a series of characters from the following:\n"
"                  _     - no move, stay on the same triangle\n"
"                  v,e,f - step over side opposite V,E,F\n"
"                  V,E,F - step two trianglesi, rotating about V,E,F,\n"
"                          according to: V=ef, E=fv, F=ve\n"
"  -c <op>   Conway polyhedron notation operator, or 'list' to list all\n"
"            available operators with their corresponding patterns\n"
"  -r        relabel pattern, exactly three letters VEF written in any order\n"
"            e.g. EFV relabels the pattern as V->E,v->e,E->F,e->f,F->V,f->v\n"
"  -M        input geometry is a 'meta' tiling, don't apply meta operation\n"
"  -i        tiles are coloured by index number (default, colour by value)\n"
"  -u        output only one example of each type of tile (one per path)\n"
"  -a        add the 'meta'-transformed base\n"
"  -f <ht>   lift the face centres by this height\n"
"  -q        quiet, don't print report\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}
// clang-format on

void wy_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":ho:p:c:f:r:Miuqa")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'p':
      print_status_or_exit(tiling.read_pattern(optarg), 'p');
      break;

    case 'c':
      if(string(optarg) == "list") {
        tiling.print_conway_list();
        exit(0);
      }
      print_status_or_exit(tiling.read_conway(optarg), 'c');
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
      print_status_or_exit(read_double(optarg, &face_ht), 'f');
      break;

    case 'i':
      color_with_value = false;
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
  for(int i = 0; i < (int)meta.faces().size(); i++)
    meta.colors(FACES).set(i, (i%2) ? light : dark);
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
  if (!info.is_orientable())
    opts.error("base polyhedron is not orientable");
  if (!opts.input_is_meta && !info.is_oriented()) {
    opts.warning("base polyhedron is not oriented; it will be oriented.");
    geom.orient();
  }

  Tiling &tiling = opts.tiling;
  if(opts.relabel != "")
    opts.print_status_or_exit(tiling.relabel_pattern(opts.relabel), 'm');

  Status stat = tiling.set_geom(geom, opts.input_is_meta, opts.face_ht);
  if (stat.is_error())
    opts.print_status_or_exit(stat, 'm');
  Geometry ogeom;
  vector<int> tile_counts;
  opts.print_status_or_exit(tiling.make_tiling(ogeom, &tile_counts));

  if(!opts.quiet) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Tiling pattern: %s\n", tiling.pattern_string().c_str());
    fprintf(stderr, "Tile Counts:\n");
    for(unsigned int i=0; i<tile_counts.size(); i++)
      fprintf(stderr, "  %-4u: %d\n", i, tile_counts[i]);
    fprintf(stderr, "\n");
  }

  if (opts.color_with_value) {
    Coloring clrng(&ogeom);
    clrng.add_cmap(colormap_from_name("spread"));
    clrng.e_apply_cmap();
    clrng.f_apply_cmap();
    Coloring v_clrng(&ogeom);
    v_clrng.add_cmap(
        colormap_from_name("map_red:green:blue:yellow:cyan:magenta:grey80"));
    v_clrng.v_apply_cmap();
  }

  if (opts.add_meta) {
    Geometry meta = tiling.get_meta();
    if(opts.color_with_value)
      color_meta(meta);
    ogeom.append(meta);
  }

  opts.write_or_error(ogeom, opts.ofile);

  return 0;
}

