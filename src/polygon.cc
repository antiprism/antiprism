/*
   Copyright (c) 2003-2016, Adrian Rossiter

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
   Name: polygon.cc
   Description: program to generate polyhedra based on polygons
   Project: Antiprism - http://www.antiprism.com
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <map>
#include <string>
#include <vector>

#include "../base/antiprism.h"

using std::string;
using std::vector;
using std::map;

using namespace anti;

class pg_opts : public ProgramOpts {
public:
  Polygon pgon;
  string subtype_str;
  int num_sides;
  int step;
  string ofile;

  pg_opts() : ProgramOpts("polygon") {}
  void process_command_line(int argc, char **argv);
  void usage();
};

// clang-format off
void pg_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] type num_sides\n"
"\n"
"Make polyhedra based on polygons. The polyhedron type and subtype may be\n"
"given by (partial) name or number\n"
"   1. prism (angle: polygons twist, forming an antiprism)\n"
"        subtypes: 1. antiprism, from triangulating side faces\n"
"                  2. trapezohedron, with this antiprism-ised belt\n"
"                  3. crown (-A integer to specify a second polygon step)\n"
"   2. antiprism (angle: polygons twist)\n"
"        subtypes: 1. trapezohedron, with this antiprism belt\n"
"                  2. antihermaphrodite, with this antiprism base\n"
"                  3. scalenohedron (-L for apex height)\n"
"                  4. subdivided_scalenohedron (-L for apex height)\n"
"                  5. crown (-A integer to specify a second polygon step)\n"
"   3. pyramid (angle: base separates, polygons twist to antihermaphrodite)\n"
"        subtypes: 1. antihermaphrodite, with this antiprism base\n"
"                  2. elongated (-L for prism height)\n"
"                  3. gyroelongated (-L for antiprism height)\n"
"   4. dipyramid (angle: base separates, polygons twist to trapezohedron)\n"
"        subtypes: 1. trapezohedron, with this pyramid apex\n"
"                  2. elongated (-L for prism height)\n"
"                  3. gyroelongated (-L for antiprism height)\n"
"                  4. dipyramid_scalenohedron (-R for alternate vertex radius)\n"
"   5. cupola\n"
"        subtypes: 1. elongated (-L for prism height)\n"
"                  2. gyroelongated (-L for antiprism height)\n"
"                  3. cupoloid\n"
"   6. orthobicupola\n"
"        subtypes: 1. elongated (-L for prism height)\n"
"                  2. gyroelongated (-L for antiprism height)\n"
"   7. gyrobicupola\n"
"        subtypes: 1. elongated (-L for prism height)\n"
"                  2. gyroelongated (-L for antiprism height)\n"
"   8. snub-antiprism\n"
"        subtypes: 1. inverted, triangle band inverted\n"
"   9. dihedron\n"
"        subtypes: 1. polygon\n"
"  10. crown polyhedron (either -A integer to specify a second polygon step,\n"
"                            or -a angle to specify a twist\n"
"\n"
"num_sides is a number (N) optionally followed by / and a second\n"
"number (N/D). N is the number of vertices spaced equally on a\n"
"circle, and each vertex is joined to the Dth (default: 1) vertex\n"
"moving around the circle.\n"
"\n"
"Options\n"
"%s"
"  -s <subt> a number or name (see type list above) indicting a subtype\n"
"            or modification of a polyhedron\n"
"  -e <len>  length of polygon edges (default: 1)\n"
"  -E <len>  length of non-polygon edges (default: calculate from -l)\n"
"  -r <rad>  circumradius of polygon (default: calculate from -e)\n"
"  -R <rad>  a second radius value\n"
"  -l <hgt>  height of upper vertices above base polygon\n"
"            (default: circumradius of polygon)\n"
"  -L <hgt>  a second height or length value\n"
"  -a <ang>  twist angle (degrees)\n"
"  -A <ang>  a second twist angle (degrees)\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}
// clang-format on

void pg_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;
  double val;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":he:E:r:R:l:L::a:A:s:o")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 's':
      subtype_str = optarg;
      break;

    case 'r':
      print_status_or_exit(read_double(optarg, &val), c);
      pgon.set_radius(0, val);
      break;

    case 'R':
      print_status_or_exit(read_double(optarg, &val), c);
      pgon.set_radius(1, val);
      break;

    case 'l':
      print_status_or_exit(read_double(optarg, &val), c);
      pgon.set_height(0, val);
      break;

    case 'L':
      print_status_or_exit(read_double(optarg, &val), c);
      pgon.set_height(1, val);
      break;

    case 'e':
      print_status_or_exit(read_double(optarg, &val), c);
      pgon.set_edge(0, val);
      break;

    case 'E':
      print_status_or_exit(read_double(optarg, &val), c);
      pgon.set_edge(1, val);
      break;

    case 'a':
      print_status_or_exit(read_double(optarg, &val), c);
      pgon.set_twist_angle(0, deg2rad(val));
      break;

    case 'A':
      print_status_or_exit(read_double(optarg, &val), c);
      pgon.set_twist_angle(1, deg2rad(val));
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (argc - optind != 2)
    error("must give two arguments");

  // Determine type
  string params = "|prism|antiprism|pyramid|dipyramid|cupola|orthobicupola"
                  "|gyrobicupola|snub-antiprism|dihedron|crown";

  string arg_id;
  print_status_or_exit(get_arg_id(argv[optind], &arg_id, params.c_str(),
                                  argmatch_default | argmatch_add_id_maps),
                       "polyhedron type");

  int type = atoi(arg_id.c_str());
  pgon.set_type(type);

  // Determine subtype
  if (subtype_str != "") {
    params = "0";
    if (type == Polygon::prism)
      params += "|antiprism|trapezohedron|crown";
    else if (type == Polygon::antiprism)
      params += "|trapezohedron|antihermaphrodite|scalenohedron"
                "|subdivided_scalenohedron|crown";
    //"|subdivided_scalenohedron|crown";
    else if (type == Polygon::pyramid)
      params += "|antihermaphrodite|elongated|gyroelongated";
    else if (type == Polygon::dipyramid)
      params += "|trapezohedron|elongated|gyroelongated"
                "|dipyramid_scalenohedron";
    else if (type == Polygon::orthobicupola || type == Polygon::gyrobicupola)
      params += "|elongated|gyroelongated";
    else if (type == Polygon::cupola)
      params += "|elongated|gyroelongated|cupoloid";
    else if (type == Polygon::snub_antiprism)
      params += "|inverted";
    else if (type == Polygon::dihedron)
      params += "|polygon";

    print_status_or_exit(get_arg_id(subtype_str.c_str(), &arg_id,
                                    params.c_str(),
                                    argmatch_default | argmatch_add_id_maps),
                         "polyhedron subtype");

    pgon.set_subtype(atoi(arg_id.c_str()));
  }

  // read polygon
  int step = 1;
  char *p = strchr(argv[optind + 1], '/');
  if (p != nullptr) {
    *p++ = '\0';
    print_status_or_exit(read_int(p, &step),
                         "number of sides (denominator of polygon fraction)");
  }

  int num_sides;
  print_status_or_exit(read_int(argv[optind + 1], &num_sides),
                       "number of sides");
  print_status_or_exit(pgon.set_fraction(num_sides, step), "number of sides");

  unsigned int param_flags;
  vector<unsigned int> conflict_flags;
  pgon.get_acceptable_params(param_flags, conflict_flags);

  for (int i = 0; i < 2; i++)
    if (pgon.value_is_set(pgon.get_edge(i)) && pgon.get_edge(i) < 0.0)
      error("edge length cannot be negative", (i) ? 'E' : 'e');

  map<int, char> flag_char;
  flag_char[Polygon::R0] = 'r';
  flag_char[Polygon::R1] = 'R';
  flag_char[Polygon::H0] = 'l';
  flag_char[Polygon::H1] = 'L';
  flag_char[Polygon::E0] = 'e';
  flag_char[Polygon::E1] = 'E';
  flag_char[Polygon::A0] = 'a';
  flag_char[Polygon::A1] = 'A';

  unsigned int params_set = pgon.get_params_set();
  if (!(params_set & (Polygon::R0 | Polygon::E0)))
    pgon.set_edge(0, 1.0); // default for polygon size

  for (unsigned int conflict_flag : conflict_flags) {
    if ((params_set & conflict_flag) == conflict_flag) {
      fprintf(stderr, "flags=%04x\n", (params_set & conflict_flag));
      string conflict_opts = "options ";
      for (auto &mi : flag_char) {
        if (mi.first & conflict_flag)
          conflict_opts += msg_str("-%c, ", mi.second);
      }
      if (conflict_opts.size())
        conflict_opts.resize(conflict_opts.size() - 2); // clear last ", "
      error("can only set one of these options", conflict_opts);
    }
  }

  unsigned int unused_params = params_set & (~param_flags);
  if (unused_params) {
    string unused_opts = "option(s) ";
    for (auto &mi : flag_char) {
      if (mi.first & unused_params)
        unused_opts += msg_str("-%c, ", mi.second);
    }
    if (unused_opts.size())
      unused_opts.resize(unused_opts.size() - 2); // clear last ", "

    warning("set but not used", unused_opts);
  }
}

int main(int argc, char *argv[])
{
  pg_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  opts.print_status_or_exit(opts.pgon.make_poly(geom),
                            "could not construct polyhedron");

  opts.write_or_error(geom, opts.ofile);

  return 0;
}
