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
   Name: lat_grid.cc
   Description: program to make lattices and grids with intger coordinates
   Project: Antiprism - http://www.antiprism.com
*/

#include <stdio.h>
#include <stdlib.h>

#include "../base/antiprism.h"
#include "lattice_grid.h"
#include <string>
#include <vector>

using std::string;
using std::vector;

using namespace anti;

class lg_opts : public ProgramOpts {
public:
  double o_width;
  double i_width;
  Vec3d centre;
  int strut_len2;
  COORD_TEST_F coord_test;
  char container;

  string ofile;

  lg_opts()
      : ProgramOpts("lat_grid"), o_width(6), i_width(-1), strut_len2(0),
        coord_test(sc_test), container('c')
  {
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

// clang-format off
void lg_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [lat_type [outer_width [inner_width]]]\n"
"\n"
"Make a lattice or grid with integer coordinates. Lattice types (default: sc)\n"
"are followed by name and then and valid strut arguments (squares of length,\n"
"first is square of the radius of balls packed at the vertex positions)\n"
"  sc                 - Simple Cubic                         (1, 2, 3)\n"
"  fcc                - Face Centred Cubic                   (2, 4, 6, 12)\n"
"  bcc                - Body Centred Cubic                   (3, 4, 8)\n"
"  hcp                - Hexagonal Close Packing              (18)\n"
"  rh_dodec           - Rhombic Dodecahedra                  (3, 8)\n" 
"  cubo_oct           - Cuboctahedron / Octahedron           (2)\n"
"  tr_oct             - Truncated Octahedron                 (2)\n"
"  tr_tet_tet         - Truncated Tetrahedron / Tetrahedron  (2)\n"
"  tr_tet_tr_oct_cubo - Truncated Tetrahedron / \n"
"                       Truncated Octahedron / Cuboctahedron (4)\n"
"  diamond            - Diamond                              (3)\n"
"  hcp_diamond        - HCP Diamond                          (27)\n"
"  k_4                - K_4 Crystal                          (2)\n"
""
"Inner and outer widths are the sizes of the inner and outer containers.\n"
"For cubes, these are the length of a side, for spheres these are the\n"
"squares of the radii.\n"
"\n"
"Options\n"
"%s"
"  -C <cent> centre of lattice, in form \"x_val,y_val,z_val\"\n"
"  -c <type> container, c - cube (default), s - sphere\n"
"  -s <len2> create struts, the value is the square of the strut length\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}
// clang-format on

void lg_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hc:s:C:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 's':
      print_status_or_exit(read_int(optarg, &strut_len2), c);
      if (strut_len2 < 0)
        error("square of strut length cannot be negative", "s");
      break;

    case 'C':
      print_status_or_exit(centre.read(optarg), c);
      break;

    case 'c':
      if (strlen(optarg) > 1 || !strchr("cs", *optarg))
        error("method is '" + string(optarg) + "' must be c or s");
      container = *optarg;
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (argc - optind > 3) {
    error("too many arguments");
    exit(1);
  }

  if (optind < argc) {
    string lt = argv[optind];
    if (lt == "sc")
      coord_test = sc_test;
    else if (lt == "fcc")
      coord_test = fcc_test;
    else if (lt == "bcc")
      coord_test = bcc_test;
    else if (lt == "hcp")
      coord_test = hcp_test;
    else if (lt == "rh_dodec")
      coord_test = rh_dodec_test;
    else if (lt == "cubo_oct")
      coord_test = cubo_oct_test;
    else if (lt == "tr_oct")
      coord_test = tr_oct_test;
    else if (lt == "tr_tet_tet")
      coord_test = tr_tet_tet_test;
    else if (lt == "tr_tet_tr_oct_cubo")
      coord_test = tr_tet_tr_oct_cubo_test;
    else if (lt == "diamond")
      coord_test = diamond_test;
    else if (lt == "hcp_diamond")
      coord_test = hcp_diamond_test;
    else if (lt == "k_4")
      coord_test = k_4_test;
    else
      error("unknown lattice type '" + lt + "' (try " + prog_name() + " -h)",
            "lat_type");
  }

  optind++;
  if (optind < argc) {
    print_status_or_exit(read_double(argv[optind], &o_width), "outer_width");
    if (o_width < 0)
      error("cannot be negative", "outer_width");
  }

  optind++;
  if (optind < argc) {
    print_status_or_exit(read_double(argv[optind], &i_width), "inner_width");
    if (i_width > o_width)
      error("cannot be greater than outer container width", "inner_width");
  }
}

int main(int argc, char **argv)
{
  lg_opts opts;
  opts.process_command_line(argc, argv);

  int_lat_grid *lat;
  if (opts.container == 'c')
    lat = new int_lat_grid();
  else
    lat = new sph_lat_grid();

  lat->set_o_width(opts.o_width);
  lat->set_i_width(opts.i_width);
  lat->set_centre(opts.centre);
  lat->set_coord_test(opts.coord_test);

  Geometry geom;
  lat->make_lattice(geom);
  delete lat;

  if (opts.strut_len2 > 0)
    add_struts(geom, opts.strut_len2);

  opts.write_or_error(geom, opts.ofile);

  return 0;
}
