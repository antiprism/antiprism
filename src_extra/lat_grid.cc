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

#include "../base/antiprism.h"

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using std::string;
using std::vector;

using namespace anti;

bool sc_test(int /*x*/, int /*y*/, int /*z*/) // dist2 = 1, 2, 3
{
  return 1;
}

bool fcc_test(int x, int y, int z) // dist2 = 2, 4, 6, 12
{
  return (x + y + z) % 2 == 0;
}

bool bcc_test(int x, int y, int z) // dist2 = 3, 4, 8
{
  return ((x % 2 && y % 2 && z % 2) || !(x % 2 || y % 2 || z % 2));
}

bool hcp_test(int x, int y, int z) // dist2 = 18
{
  int m = x + y + z;
  int n = (m < 0); // for integer division
  return (m % 6 == 0) && ((x - m / 12 + n) % 3 == 0) &&
         ((y - m / 12 + n) % 3 == 0);
}

bool rh_dodec_test(int x, int y, int z) // dist2 = 3 (8)
{
  return ((x % 2 && y % 2 && z % 2) ||
          !(((x + y + z) / 2) % 2 == 0 || (x % 2 || y % 2 || z % 2)));
}

bool cubo_oct_test(int x, int y, int z) // dist2 = 2
{
  return (std::abs(x) % 2 + std::abs(y) % 2 + std::abs(z) % 2) == 2;
}

bool tr_oct_test(int x, int y, int z) // dist2 = 2
{
  return ((z % 4 == 0 && x % 4 && y % 4 && (x - y) % 2) ||
          (y % 4 == 0 && z % 4 && x % 4 && (z - x) % 2) ||
          (x % 4 == 0 && y % 4 && z % 4 && (y - z) % 2));
}

bool tr_tet_tet_test(int x, int y, int z) // dist2 = 2
{
  return ((x % 2 == 0 && (y % 4 + 4) % 4 == (((x + z)) % 4 + 4) % 4) ||
          (y % 2 == 0 && (z % 4 + 4) % 4 == (((y + x)) % 4 + 4) % 4) ||
          (z % 2 == 0 && (x % 4 + 4) % 4 == (((z + y)) % 4 + 4) % 4));
}

bool tr_tet_tr_oct_cubo_test(int x, int y, int z) // dist2 = 4
{
  x = std::abs(x) % 6;
  if (x > 3)
    x = 6 - x;
  y = std::abs(y) % 6;
  if (y > 3)
    y = 6 - y;
  z = std::abs(z) % 6;
  if (z > 3)
    z = 6 - z;
  int dist2 = x * x + y * y;
  return ((z % 6 == 0 && (dist2 == 2 || dist2 == 8)) ||
          (z % 6 == 1 && (dist2 == 1 || dist2 == 13)) ||
          (z % 6 == 2 && (dist2 == 4 || dist2 == 10)) ||
          (z % 6 == 3 && dist2 == 5));
}

bool diamond_test(int x, int y, int z) //  dist2 = 3
{
  return (((x % 2 + 2) % 2 + (y % 2 + 2) % 2 + (z % 2 + 2) % 2) % 3 == 0 &&
          (x / 2 + (x % 2 < 0) + y / 2 + (y % 2 < 0) + z / 2 + (z % 2 < 0)) %
                  2 ==
              0);
}

// Coordinates from Wendy Krieger
bool hcp_diamond_test(int x, int y, int z) //  dist2 = 27
{
  int pt[][3] = {{0, 0, 0}, {3, 3, 3}, {6, 0, 6}, {9, 3, 9}};
  for (auto &i : pt) {
    int tx = x - i[0];
    int ty = y - i[1];
    int tz = z - i[2];
    int sum = tx + ty + tz;
    if (sum % 24 == 0) {
      int n8 = sum / 3;
      if ((tx - n8) % 6 == 0 && (ty - n8) % 6 == 0 && (tz - n8) % 6 == 0)
        return true;
    }
  }
  return false;
}

// Coordinates from Vladimir Bulatov
bool k_4_test(int x, int y, int z) //  dist2 = 2
{
  if ((x + y + z) % 2)
    return false;
  x = (x % 4) + ((x < 0) ? 4 : 0);
  y = (y % 4) + ((y < 0) ? 4 : 0);
  z = (z % 4) + ((z < 0) ? 4 : 0);
  if ((x == 0 && y == 0 && z == 0) || (x == 0 && y == 1 && z == 3) ||
      (x == 1 && y == 0 && z == 1) || (x == 1 && y == 1 && z == 2) ||
      (x == 2 && y == 2 && z == 2) || (x == 2 && y == 3 && z == 1) ||
      (x == 3 && y == 2 && z == 3) || (x == 3 && y == 3 && z == 0))
    return true;
  return false;
}

void add_struts(Geometry &geom, int len2)
{
  const vector<Vec3d> &verts = geom.verts();
  for (unsigned int i = 0; i < verts.size(); i++)
    for (unsigned int j = i; j < verts.size(); j++) {
      if (fabs((verts[i] - verts[j]).len2() - len2) < anti::epsilon)
        geom.add_edge(make_edge(i, j));
    }
}

typedef bool (*COORD_TEST_F)(int, int, int);

class int_lat_grid {
protected:
  double o_width;
  double i_width;
  anti::Vec3d centre;
  COORD_TEST_F coord_test;

public:
  // enum { l_sc, l_fcc, l_bcc, l_rh_dodec, l_cubo_oct,
  //   l_tr_oct, l_tr_tet_tet, l_tr_oct_tr_tet_cubo, l_diamond }
  int_lat_grid() {}
  virtual ~int_lat_grid() = default;
  void set_coord_test(COORD_TEST_F func) { coord_test = func; }
  virtual void set_o_width(double w) { o_width = w; }
  virtual void set_i_width(double w) { i_width = w; }
  virtual void set_centre(anti::Vec3d cent) { centre = cent; }
  virtual void make_lattice(anti::Geometry &geom);
  // void add_struts(anti::Geometry &geom, int len2);
};

class sph_lat_grid : public int_lat_grid {
public:
  sph_lat_grid() {}
  virtual void set_o_width(double w) { o_width = w; }
  virtual void set_i_width(double w) { i_width = w; }
  virtual void make_lattice(anti::Geometry &geom);
};

void int_lat_grid::make_lattice(Geometry &geom)
{
  if (!centre.is_set())
    centre = Vec3d(1, 1, 1) * (o_width / 2.0);
  double o_off = o_width / 2.0 + anti::epsilon;
  double i_off = i_width / 2.0 - anti::epsilon;
  int i, j, k;
  for (k = int(ceil(centre[2] - o_off)); k <= centre[2] + o_off; k++)
    for (j = int(ceil(centre[1] - o_off)); j <= centre[1] + o_off; j++)
      for (i = int(ceil(centre[0] - o_off)); i <= centre[0] + o_off; i++) {
        if (i > centre[0] - i_off && i < centre[0] + i_off &&
            j > centre[1] - i_off && j < centre[1] + i_off &&
            k > centre[2] - i_off && k < centre[2] + i_off)
          continue;
        if (coord_test(i, j, k))
          geom.add_vert(Vec3d(i, j, k));
      }
}

void sph_lat_grid::make_lattice(Geometry &geom)
{
  if (!centre.is_set())
    centre = Vec3d(0, 0, 0);
  double o_off = o_width + anti::epsilon;
  double i_off = i_width - anti::epsilon;
  int i, j, k;
  for (k = int(ceil(centre[2] - o_off)); k <= centre[2] + o_off; k++)
    for (j = int(ceil(centre[1] - o_off)); j <= centre[1] + o_off; j++)
      for (i = int(ceil(centre[0] - o_off)); i <= centre[0] + o_off; i++) {
        double dist2 = (Vec3d(i, j, k) - centre).len2();
        if (o_off < dist2 || i_off > dist2)
          continue;
        if (coord_test(i, j, k))
          geom.add_vert(Vec3d(i, j, k));
      }
}

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

void lg_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] [lat_type [outer_width [inner_width]]]

Make a lattice or grid with integer coordinates. Lattice types (default: sc)
are followed by name and then and valid strut arguments (squares of length,
first is square of the radius of balls packed at the vertex positions)
  sc                 - Simple Cubic                         (1, 2, 3)
  fcc                - Face Centred Cubic                   (2, 4, 6, 12)
  bcc                - Body Centred Cubic                   (3, 4, 8)
  hcp                - Hexagonal Close Packing              (18)
  rh_dodec           - Rhombic Dodecahedra                  (3, 8)\
  cubo_oct           - Cuboctahedron / Octahedron           (2)
  tr_oct             - Truncated Octahedron                 (2)
  tr_tet_tet         - Truncated Tetrahedron / Tetrahedron  (2)
  tr_tet_tr_oct_cubo - Truncated Tetrahedron / 
                       Truncated Octahedron / Cuboctahedron (4)
  diamond            - Diamond                              (3)
  hcp_diamond        - HCP Diamond                          (27)
  k_4                - K_4 Crystal                          (2)

Inner and outer widths are the sizes of the inner and outer containers.
For cubes, these are the length of a side, for spheres these are the
squares of the radii.

Options
%s
  -C <cent> centre of lattice, in form \"x_val,y_val,z_val\"
  -c <type> container, c - cube (default), s - sphere
  -s <len2> create struts, the value is the square of the strut length
  -o <file> write output to file (default: write to standard output)

)",
          prog_name(), help_ver_text);
}

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
