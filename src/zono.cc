/*
   Copyright (c) 2003-2021, Adrian Rossiter

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
   Name: zono.cc
   Description: make zonohedra
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using std::string;
using std::vector;

using namespace ::anti;

class zo_opts : public ProgramOpts {
public:
  char method = 'v';
  Vec3d center = Vec3d(0, 0, 0);
  bool centroid = false;
  bool out_star = false;
  bool unit_len = false;
  Color zone_col;
  int pol_num = 0;
  int pol_spiral_step = 0;
  Polygon pgon = 2;
  bool non_polar_opt = false;
  Coloring clrngs[3];
  string ifile;
  Geometry seed_geom;
  Geometry star_translation;
  string open_flgs;
  string ofile;

  zo_opts() : ProgramOpts("zono") { read_colorings(clrngs, "spread"); }
  void process_command_line(int argc, char **argv);
  void usage();
};

void zo_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] [star_file] 

Make a zonohedron or add zones to a convex seed polyhdron. The zones are
created from a star of vectors, which can be based on a polyhdron (input
model and option -m) or initialised (option -P) to make a polar zonohedron.
If input_file is not given the program reads from standard input

Options
%s
  -m <mthd> method to create star from input, can be
               v - centre to vertices are vectors (default)
               a - all vertex to vertex are vectors
               i - implicit edges (face sides) are vectors
               e - explicit edges are vectors
  -c <cent> centre of points for method v, in form 'X,Y,Z', or C to use
            centroid (default '0,0,0')
  -s        output the star (instead of the zonohedron)
  -S <poly> seed model to add zones to, must be convex
  -u        make vectors unit length
  -C <col>  colour for new zone faces
  -P <arg>  polar zonohedron from ordered star, arg can be an offset polygon
            given as an integer or fraction (e.g. 5, 7/2) or 's' to use
            star_file. Optionally, follow by a comma and an integer to
            to make a spirallohedron with that spiral width (default:0, a
            polar zonohedron). Any further comma separated parts are colour
            maps to colour the faces.
  -T <star> translation surface made from two ordered stars (argument star
            and star_file) optionally followed by a comma and two characters
            (or one to use for both) to indicate if the stars of vectors
            are open, or close to form a loop: d-detect, c-force close,
             o-leave open.
  -o <file> write output to file (default: write to standard output)

)",
          prog_name(), help_ver_text);
}

void zo_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;
  Split parts;
  string map_str;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hm:c:S:suC:P:T:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'm':
      if (!(strlen(optarg) == 1 && strchr("vaei", *optarg)))
        error("unknown method '" + string(optarg) + "'", c);
      method = *optarg;
      non_polar_opt = true;
      break;

    case 'c':
      if ((strlen(optarg) == 1 && strchr("Cc", *optarg)))
        centroid = true;
      else if (!center.read_maths(optarg))
        error("invalid centre '" + string(optarg) + "'", c);
      non_polar_opt = true;
      break;

    case 's':
      out_star = true;
      break;

    case 'S': {
      read_or_error(seed_geom, optarg);
      if (get_convexity(seed_geom) != Convexity::convex_strict)
        warning("seed geometry is not strictly convex", c);
      break;
    }

    case 'u':
      unit_len = true;
      break;

    case 'C':
      print_status_or_exit(zone_col.read(optarg));
      break;

    case 'P':
      parts.init(optarg, ",");
      if (parts.size() > 1)
        print_status_or_exit(read_int(parts[1], &pol_spiral_step), c);
      else
        pol_spiral_step = 0; // default

      if (strcmp(parts[0], "s") == 0)
        pol_num = -1; // use star argument
      else {
        int pol_denom = 0;
        print_status_or_exit(read_fraction(parts[0], &pol_num, &pol_denom), c);
        if (pol_num < 2)
          error("number of sides must be 2 or greater", c);
        if (pol_denom < 1)
          error("denominator must be 1 or greater", c);
        if (pol_denom % pol_num == 0)
          error("denominator cannot be a multiple of the number of "
                "sides",
                c);
        pgon = Polygon(pol_num, pol_denom);

        if (pol_spiral_step % pol_denom &&
            pol_spiral_step % (pol_num - pol_denom))
          error("spiral step must be divisible by the polygon denominator, or "
                "by the numerator minus the denominator",
                c);
        if (pol_spiral_step && pol_spiral_step % pol_num == 0)
          error("spiral step must not be divisible by the polygon numerator",
                c);
      }

      map_str = "";
      if (parts.size() > 2) {
        for (int i = 2; i < (int)parts.size(); i++)
          map_str += parts[i] + string(",");
      }
      else
        map_str = "spread";

      print_status_or_exit(read_colorings(clrngs, map_str.c_str()), c);
      break;

    case 'T':
      parts.init(optarg, ",");
      if (parts.size() > 1) {
        if (strspn(parts[1], "doc") != strlen(parts[1]))
          error("close setting can contain only the characters d, o, c", c);
        if (strlen(parts[1]) > 2)
          error("close setting cannot be longer than two characters", c);
        open_flgs = parts[1];
      }

      read_or_error(star_translation, parts[0]);
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (pol_num > 0 && star_translation.is_set())
    error("options -P and -T cannot be used together");

  if (pol_num > 0 && (non_polar_opt || argc - optind))
    error("option -P polygon parameter is not compatible with options m, c or "
          "an input file");

  if (argc - optind > 1)
    error("too many arguments");

  if (argc - optind == 1)
    ifile = argv[optind];
}

int main(int argc, char **argv)
{
  zo_opts opts;
  opts.process_command_line(argc, argv);

  vector<Vec3d> star;
  if (opts.pol_num > 0) {
    Geometry gstar;
    Polygon pgon2(opts.pol_num);
    pgon2.add_polygon(gstar, sqrt(0.5));
    star = gstar.verts();
  }
  else {
    Geometry geom;
    opts.read_or_error(geom, opts.ifile);

    if (opts.centroid)
      opts.center = geom.centroid();

    star = get_star(geom, opts.method, opts.center);
  }

  // Set star vectors to unit, and remove any parallel vectors
  if (opts.unit_len) {
    for (auto &i : star)
      i.to_unit();
    int star_final_sz = star.size();
    for (int i = 0; i < star_final_sz; i++)
      for (int j = i + 1; j < star_final_sz; j++) {
        if (vcross(star[i], star[j]).len() < anti::epsilon)
          std::swap(star[j], star[--star_final_sz]);
      }
    star.resize(star_final_sz);
  }

  Geometry zono;
  if (opts.out_star)
    zono.add_verts(star);
  else if (opts.seed_geom.is_set())
    opts.print_status_or_exit(make_zonohedrified_polyhedron(
        zono, opts.seed_geom, star, opts.zone_col));
  else if (opts.pol_num) {
    make_polar_zonohedron(zono, star,
                          opts.pgon.get_step() * opts.pgon.get_parts(),
                          opts.pol_spiral_step);
    opts.clrngs[FACES].set_geom(&zono);
    opts.clrngs[FACES].f_apply_cmap();
  }
  else if (opts.star_translation.is_set()) {
    Geometry g_star;
    g_star.add_verts(star);
    make_translation_surface(zono, g_star, opts.star_translation,
                             opts.open_flgs);
  }
  else {
    opts.print_status_or_exit(make_zonohedron(zono, star));
    if (opts.zone_col.is_set())
      Coloring(&zono).f_one_col(opts.zone_col);
  }

  opts.write_or_error(zono, opts.ofile);

  return 0;
}
