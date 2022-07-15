/*
   Copyright (c) 2003-2022, Adrian Rossiter, Roger Kaufman

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
   Name: common.cc
   Description: code shared by source code in /src
   Project: Antiprism - http://www.antiprism.com
*/

#include "common.h"
#include "../base/antiprism.h"

#include <cstdio>
#include <memory>
#include <string>
#include <vector>

using std::set;
using std::string;
using std::vector;

using namespace anti;

// for lat_grid.cc

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

// for lat_util.cc and bravais.cc

// parse special color entry strings for bravais and lat_util
//
// 1 arg  - color name/index
// 2 args - color name/index, transparency
// 2 args - color name/index, assignment
// 3 args - color name/index, transparency, assignment
// 3 args - color r,g,b
// 4 args - color r,g,b, transparency
// 4 args - color r,g,b, assignment
// 5 args - color r,g,b, transparency, assignment
// vcol is a vector size 4 of pre-allocated colors
void parse_color_string(const ProgramOpts *opts, const char *optarg,
                        const char c, const string &allowed_chars,
                        vector<Color> &vcol)
{
  Split parts(optarg, ",");
  unsigned int parts_sz = parts.size();

  int dummy;
  Status is_int = read_int(parts[0], &dummy);
  if (is_int.is_ok() && parts_sz > 5)
    opts->error("the argument is a color value and has more than 5 parts", c);
  else if (!is_int.is_ok() && parts_sz > 3)
    opts->error("the argument is a color name and has more than 3 parts", c);

  Color col;
  bool valid_color = false;
  unsigned int next_parms_idx = 1;

  // see if entry is a valid r,g,b
  if (parts_sz >= 3) {
    string color_str_tmp = string(parts[0]) + "," + parts[1] + "," + parts[2];
    if (col.read(color_str_tmp.c_str())) {
      if (col.is_set())
        valid_color = true;
      next_parms_idx = 3;
    }
  }

  // check to see if it is a valid color name
  if (!valid_color) {
    if (col.read(parts[0])) {
      if (col.is_set())
        valid_color = true;
      next_parms_idx = 1;
    }
  }

  // if the color is not valid, output the color error
  if (!valid_color)
    opts->print_status_or_exit(col.read(parts[0]), c);

  // check for transparency
  if (parts_sz > next_parms_idx) {
    // if the next part is an integer
    int opacity = 255;
    if (read_int(parts[next_parms_idx], &opacity)) {
      if (opacity < 0 || opacity > 255)
        opts->error(
            msg_str("transparency is '%d' and must be between 0 and 255",
                    opacity),
            c);
      if (!col.set_alpha(opacity))
        opts->warning("transparency has no effect on map indexes or invisible",
                      c);
      next_parms_idx++;
    }
  }

  unsigned int conv_elems = 15;
  if (parts_sz > next_parms_idx) {
    if (strspn(parts[next_parms_idx], allowed_chars.c_str()) !=
        strlen(parts[next_parms_idx]))
      opts->error(msg_str("elements to map are '%s' must be from '%s'",
                          parts[next_parms_idx], allowed_chars.c_str()),
                  c);
    conv_elems = 8 * (strchr(parts[next_parms_idx], 'h') != nullptr) +
                 4 * (strchr(parts[next_parms_idx], 'v') != nullptr) +
                 2 * (strchr(parts[next_parms_idx], 'c') != nullptr) +
                 1 * (strchr(parts[next_parms_idx], 'l') != nullptr);
    next_parms_idx++;
  }

  for (int i = 0; i < 4; i++) {
    if (conv_elems & (1 << i))
      vcol[i] = col;
  }

  if (parts_sz > next_parms_idx)
    opts->error(msg_str("extraneous input: '%s'", parts[next_parms_idx]), c);
}

double lattice_radius(const Geometry &geom, const char radius_type)
{
  Geometry tgeom = geom;

  // allow radius to be calculated for some polygons
  int dimensions;
  Status stat = tgeom.set_hull("", &dimensions);
  if (stat.is_error()) {
    fprintf(stderr, "%s\n", stat.c_msg());
    fprintf(stderr,
            "lattice_radius: warning: convex hull could not be created\n");
    return 0;
  }
  tgeom.orient(1); // positive orientation

  GeometryInfo rep(tgeom);
  rep.set_center(centroid(tgeom.verts()));

  double radius = 0;
  if (radius_type == 'k')
    radius = rep.vert_dist_lims().max;
  else if (radius_type == 's' && dimensions > 1)
    radius = rep.face_dist_lims().min;
  else if (radius_type == 'l' && dimensions > 1)
    radius = rep.face_dist_lims().max;

  return radius;
}

// utility function to convert a single Vec3d to a vector of Vec3ds
static inline vector<Vec3d> as_vector(const Vec3d &v)
{
  vector<Vec3d> ret;
  ret.push_back(v);
  return ret;
}

void geom_container_clip(Geometry &geom, Geometry &container,
                         const double radius, const Vec3d &offset,
                         const bool verbose, const double eps)
{
  // container has to be convex and 3 dimensional
  Status stat = container.set_hull();
  if (stat.is_error()) {
    fprintf(stderr, "%s\n", stat.c_msg());
    fprintf(stderr,
            "geom_container_clip: warning: convex hull could not be created\n");
    return;
  }
  container.orient(1); // positive orientation

  // standardize radius of 1 on maximum vertex. Then set radius
  Trans3d trans_m =
      Trans3d::scale((1.0 / lattice_radius(container, 'k')) * radius);
  container.transform(trans_m);

  const vector<Vec3d> &verts = geom.verts();
  Vec3d grid_cent = centroid(verts);
  if (offset.is_set())
    grid_cent += offset;

  if (verbose)
    fprintf(stderr, "info: radius = %g (square root of %g)\n", radius,
            radius * radius);

  // translate container to center of grid
  Vec3d container_cent = centroid(container.verts());
  trans_m = Trans3d::translate(-container_cent + grid_cent);
  container.transform(trans_m);

  vector<int> del_verts;
  for (unsigned int i = 0; i < verts.size(); i++) {
    if (!are_points_in_hull(as_vector(verts[i]), container,
                            INCLUSION_IN | INCLUSION_ON, eps))
      del_verts.push_back(i);
  }

  if (del_verts.size())
    geom.del(VERTS, del_verts);

  if (!verts.size())
    fprintf(
        stderr,
        "bravais_container_clip: warning: all vertices were clipped out!\n");
}

void geom_spherical_clip(Geometry &geom, const double radius,
                         const Vec3d &offset, const bool verbose,
                         const double eps)
{
  const vector<Vec3d> &verts = geom.verts();
  Vec3d cent = centroid(verts);
  if (offset.is_set())
    cent += offset;

  if (verbose)
    fprintf(stderr, "info: radius = %g (square root of %g)\n", radius,
            radius * radius);

  vector<int> del_verts;
  for (unsigned int i = 0; i < verts.size(); i++) {
    double len = (cent - verts[i]).len();
    if (double_gt(len, radius, eps))
      del_verts.push_back(i);
  }

  if (del_verts.size())
    geom.del(VERTS, del_verts);

  if (!verts.size())
    fprintf(
        stderr,
        "bravais_spherical_clip: warning: all vertices were clipped out!\n");
}

void list_grid_radii(const string &file_name, const Geometry &geom,
                     const Vec3d &list_radii_center, int report_type,
                     const double eps)
{
  FILE *ofile = stdout; // write to stdout by default
  if (file_name.length())
    ofile = fopen(file_name.c_str(), "w");
  if (!ofile) {
    fprintf(stderr, "could not output file \'%s\'\n", file_name.c_str());
    return;
  }

  const vector<Vec3d> &verts = geom.verts();
  Vec3d cent = list_radii_center;
  if (!cent.is_set())
    cent = centroid(verts);

  vector<double> radii;
  for (const auto &vert : verts)
    radii.push_back((cent - vert).len());

  sort(radii.begin(), radii.end());

  // eliminate 0 radius
  int start = 0;
  if (double_eq(radii[0], 0, eps))
    start++;

  // prime things
  double comp = radii[start];
  int occur_total = 0;
  int occur = 1;
  int rank = 1;

  if (report_type == 1) {
    string buffer;
    if (!list_radii_center.is_set())
      buffer = "centroid";
    else
      buffer = msg_str("%g,%g,%g", list_radii_center[0], list_radii_center[1],
                       list_radii_center[2]);

    fprintf(ofile,
            "\nList of unique radial distances in grid using center: %s\n\n",
            buffer.c_str());

    fprintf(ofile, "Rank\tDistance\tD Squared\tOccurrence\n");
    fprintf(ofile, "----\t--------\t---------\t----------\n");
  }
  for (unsigned int i = start + 1; i < radii.size(); i++) {
    if (double_eq(radii[i], comp, eps))
      occur++;
    else {
      occur_total += occur;
      if (report_type == 1)
        fprintf(ofile, "%d\t%-8g\t%-8g\t%d\n", rank, comp, comp * comp, occur);
      else if (report_type == 2)
        fprintf(ofile, "%.17g\n", comp);
      comp = radii[i];
      occur = 1;
      rank++;
    }
  }
  occur_total += occur;
  if (report_type == 1)
    fprintf(ofile, "%d\t%-8g\t%-8g\t%d\n\n", rank, comp, comp * comp, occur);
  else if (report_type == 2)
    fprintf(ofile, "%.17g\n", comp);

  if (report_type == 1)
    fprintf(ofile, "Total occurrences = %d\n\n", occur_total);
}

void list_grid_struts(const string &file_name, const Geometry &geom,
                      int report_type, const double eps)
{
  FILE *ofile = stdout; // write to stdout by default
  if (file_name.length())
    ofile = fopen(file_name.c_str(), "w");
  if (!ofile) {
    fprintf(stderr, "could not output file \'%s\'\n", file_name.c_str());
    return;
  }

  const vector<Vec3d> &verts = geom.verts();

  vector<double> struts;
  for (unsigned int i = 0; i < verts.size(); i++)
    for (unsigned int j = i + 1; j < verts.size(); j++)
      struts.push_back((verts[i] - verts[j]).len());

  sort(struts.begin(), struts.end());

  // eliminate 0 radius
  int start = 0;
  if (double_eq(struts[0], 0, eps))
    start++;

  // prime things
  double comp = struts[start];
  int occur_total = 0;
  int occur = 1;
  int rank = 1;

  if (report_type == 1) {
    fprintf(ofile, "\nList of unique strut lengths in grid\n\n");

    fprintf(ofile, "Rank\tDistance\tD Squared\tOccurrence\n");
    fprintf(ofile, "----\t--------\t---------\t----------\n");
  }
  for (unsigned int i = start + 1; i < struts.size(); i++) {
    if (double_eq(struts[i], comp, eps))
      occur++;
    else {
      occur_total += occur;
      if (report_type == 1)
        fprintf(ofile, "%d\t%-8g\t%-8g\t%d\n", rank, comp, comp * comp, occur);
      else if (report_type == 2)
        fprintf(ofile, "%.17g\n", comp);
      comp = struts[i];
      occur = 1;
      rank++;
    }
  }
  occur_total += occur;
  if (report_type == 1)
    fprintf(ofile, "%d\t%-8g\t%-8g\t%d\n\n", rank, comp, comp * comp, occur);
  else if (report_type == 2)
    fprintf(ofile, "%.17g\n", comp);

  if (report_type == 1)
    fprintf(ofile, "Total occurrences = %d\n\n", occur_total);
}

void add_color_struts(Geometry &geom, const double len2, Color &edge_col,
                      const double eps)
{
  const vector<Vec3d> &verts = geom.verts();

  for (unsigned int i = 0; i < verts.size(); i++)
    for (unsigned int j = i; j < verts.size(); j++) {
      if (fabs((verts[i] - verts[j]).len2() - len2) < eps)
        geom.add_edge(make_edge(i, j), edge_col);
    }
}

void color_centroid(Geometry &geom, Color &cent_col, const double eps)
{
  const vector<Vec3d> &verts = geom.verts();
  Vec3d cent = centroid(verts);
  int cent_idx = find_vert_by_coords(geom, cent, eps);

  if (cent_idx == -1)
    geom.add_vert(cent, cent_col);
  else
    geom.colors(VERTS).set(cent_idx, cent_col);
}

int get_voronoi_geom(Geometry &geom, Geometry &vgeom, const bool central_cells,
                     const bool one_cell_only, const double eps)
{
  // do this in case compound lattice was sent. Simultaneous points cause
  // problems for Voronoi Cells
  merge_coincident_elements(geom, "vef", eps);

  // store convex hull of lattice as speed up
  // has to be 3 dimensional
  Geometry hgeom = geom;
  Status stat = hgeom.set_hull();
  if (stat.is_error()) {
    fprintf(stderr, "%s\n", stat.c_msg());
    fprintf(stderr,
            "get_voronoi_geom: warning: convex hull could not be created\n");
    return 0;
  }
  hgeom.orient(1); // positive orientation

  // Add centroid to a vector. Needed in this form for are_points_in_hull()
  vector<Vec3d> cent = as_vector(centroid(hgeom.verts()));

  vector<Geometry> cells;
  get_voronoi_cells(geom.verts(), &cells);

  for (auto &cell : cells) {
    if (central_cells &&
        !are_points_in_hull(cent, hgeom, INCLUSION_IN | INCLUSION_ON, eps)) {
      continue;
    }
    else if (!are_points_in_hull(cell.verts(), hgeom,
                                 INCLUSION_IN | INCLUSION_ON, eps)) {
      continue;
    }
    vgeom.append(cell);
    if (one_cell_only)
      break;
  }

  if (!(vgeom.verts()).size()) {
    fprintf(stderr,
            "get_voronoi_geom: warning: after Voronoi cells, geom is empty\n");
    return 0;
  }

  return 1;
}

// for lat_util.cc, bravais.cc and waterman.cc

// Rotational octahedral by Adrian Rossiter
Vec3d sort_Vec3d_chiral(const Vec3d &v, const double eps)
{
  Vec3d c = v;
  // Rotate into positive octant
  if (c[0] < 0) {
    c[0] = -c[0];
    c[2] = -c[2];
  }
  if (c[1] < 0) {
    c[1] = -c[1];
    c[2] = -c[2];
  }
  if (c[2] < 0) {
    std::swap(c[0], c[1]);
    c[2] = -c[2];
  }

  // if c[1] is maximum rotate to first place: 1,2,0
  if (c[1] > c[0] && c[1] > c[2] - eps)
    c = Vec3d(c[1], c[2], c[0]);
  else
      // if c[2] is maximum rotate to first place: 2,0,1
      if (c[2] > c[0] && c[2] > c[1] + eps)
    c = Vec3d(c[2], c[0], c[1]);
  else
    // if c[0] is maximum do nothing
    c = Vec3d(c[0], c[1], c[2]);

  // Check whether c is near negative triangle external boundary, and
  // rotate to corresponding positive triangle boundary if so.
  if (double_eq(c[0], c[1], eps))
    c = Vec3d(c[1], c[0], c[2]);
  if (double_eq(c[2], 0, eps))
    c = Vec3d(c[0], -c[2], c[1]);

  return c;
}

// sort an absolute value of Vec3d without altering original Vec3d
Vec3d sort_Vec3d(Vec3d &v)
{
  vector<double> c;
  c.push_back(fabs(v[0]));
  c.push_back(fabs(v[1]));
  c.push_back(fabs(v[2]));

  sort(c.begin(), c.end());

  return (Vec3d(c[0], c[1], c[2]));
}

void color_by_symmetry_normals(Geometry &geom, const char color_method,
                               const int face_opacity, const double eps)
{
  const vector<vector<int>> &faces = geom.faces();
  const vector<Vec3d> &verts = geom.verts();

  string map_name = "rnd";
  if (face_opacity > -1)
    map_name += msg_str("_A%g", (double)face_opacity / 255);
  std::unique_ptr<ColorMap> cmap(colormap_from_name(map_name.c_str()));

  for (unsigned int i = 0; i < faces.size(); i++) {
    Vec3d norm = face_norm(verts, faces[i]).unit();
    if (color_method == 's' || color_method == 'S')
      norm = sort_Vec3d(norm);
    else if (color_method == 'c' || color_method == 'C')
      norm = sort_Vec3d_chiral(norm, eps);

    long idx = (long)(norm[0] * 1000000) + (long)(norm[1] * 10000) +
               (long)norm[2] * 100;
    if (color_method == 'S' || color_method == 'C')
      geom.colors(FACES).set(i, cmap->get_col(idx));
    else
      geom.colors(FACES).set(i, idx);
  }
}

void color_edges_by_sqrt(Geometry &geom, const char color_method)
{
  geom.add_missing_impl_edges();

  std::unique_ptr<ColorMap> cmap(colormap_from_name("rnd"));
  // e_coloring clrg(&geom);
  for (unsigned int i = 0; i < geom.edges().size(); i++) {
    // geom.colors(EDGES).set(i, int(floor(pow(geom.edge_len(i),2)+0.5)));
    int idx = int(floor(pow(geom.edge_len(i), 2) + 0.5));
    if (color_method == 'R')
      geom.colors(EDGES).set(i, cmap->get_col(idx));
    // geom.colors(EDGES).set(i,clrg.idx_to_rand_val(idx));
    else
      geom.colors(EDGES).set(i, idx);
  }
}

void convex_hull_report(const Geometry &geom, const bool add_hull)
{
  GeometryInfo rep(geom);
  fprintf(stderr, "\n");
  fprintf(stderr, "convex hull information:\n");
  if (!add_hull)
    fprintf(stderr, "num_verts = %d\n", rep.num_verts());
  fprintf(stderr, "num_faces = %d\n", rep.num_faces());
  if (!add_hull && rep.num_verts() > 1)
    fprintf(stderr, "num_edges = %d\n", rep.num_verts() + rep.num_faces() - 2);
  if (rep.num_verts() > 2) {
    double area = rep.face_areas().sum;
    fprintf(stderr, "area      = %.17g\n", area);
    fprintf(stderr, "volume    = %.17g\n", fabs(rep.volume()));
    if (area) {
      fprintf(stderr, "isoperimetric quotient (spherical nature: v^2/a^3 x 36 "
                      "x PI = 1 is sphere)\n");
      fprintf(stderr, "  iq      = %.17g\n", rep.isoperimetric_quotient());
    }
  }
  fprintf(stderr, "end convex hull information\n");
  fprintf(stderr, "\n");
}

// for stellate.cc and miller.cc

// use anti::epsilon so it is not effected by opt epsilon
void color_stellation(Geometry &geom, const char face_coloring_method,
                      const char edge_coloring_method,
                      const char vertex_coloring_method,
                      const Color &face_color, const Color &edge_color,
                      const Color &vertex_color, const int face_opacity,
                      const string &map_string, const string caller)
{
  // set color map
  Coloring clrng(&geom);
  ColorMap *cmap = colormap_from_name(map_string.c_str());
  clrng.add_cmap(cmap);

  // in case of face coloring C
  Geometry kis;

  // geom is built with face colors from the diagram
  if (face_coloring_method != 'd') {
    if (!face_coloring_method && !face_color.is_set())
      // if no color specified, clear faces
      geom.colors(FACES).clear();
    else if (face_coloring_method == 's') {
      // color faces by symmetry
      Symmetry sym;
      vector<vector<set<int>>> sym_equivs;
      sym.init(geom, &sym_equivs);
      clrng.f_sets(sym_equivs[2], true);
    }
    else if (face_coloring_method == 'c')
      // compound coloring
      clrng.f_parts(true);
    else if (face_coloring_method == 'C') {
      // color by connection
      wythoff_make_tiling(kis, geom, "k", true, false);
      // remove digons
      vector<int> dels;
      for (unsigned int i = 0; i < kis.faces().size(); i++) {
        if (kis.faces(i).size() < 3)
          dels.push_back((int)i);
      }
      kis.del(FACES, dels);
      kis.orient(1); // positive orientation

      // make new verts and edges invisible
      kis.add_missing_impl_edges();
      for (unsigned int i = 0; i < kis.verts().size(); i++) {
        int v_idx = find_vert_by_coords(geom, kis.verts()[i], anti::epsilon);
        if (v_idx == -1) {
          kis.colors(VERTS).set(i, Color::invisible);
          vector<int> edge_idx = find_edges_with_vertex(kis.edges(), i);
          for (unsigned int j = 0; j < edge_idx.size(); j++)
            kis.colors(EDGES).set(edge_idx[j], Color::invisible);
        }
      }
      // the old faces are cleared and kis faces added
      geom.clear(FACES);
      geom.append(kis);
      int blend_type = 1; // first color, invisible edges stay
      merge_coincident_elements(geom, "vef", blend_type, anti::epsilon);

      for (unsigned int i = 0; i < geom.faces().size(); i++) {
        vector<int> face = geom.faces()[i];
        unsigned int fsz = face.size();
        // face to face
        // connections with invisible faces are ignored
        int connections = 0;
        for (unsigned int j = 0; j < fsz; j++) {
          int v1 = face[j];
          int v2 = face[(j + 1) % fsz];
          vector<int> edge = make_edge(v1, v2);
          vector<int> face_idx = find_faces_with_edge(geom.faces(), edge);
          int edge_no = find_edge_in_edge_list(geom.edges(), edge);
          if (!(geom.colors(EDGES).get(edge_no)).is_invisible())
            connections += face_idx.size();
        }
        geom.colors(FACES).set(i, cmap->get_col(connections));
      }
    }
    else
      // use color selected
      clrng.f_one_col(face_color);
  }

  // collect invisible edges
  vector<int> invisible_edges;
  for (unsigned int i = 0; i < geom.edges().size(); i++)
    if ((geom.colors(EDGES).get(i)).is_invisible())
      invisible_edges.push_back(i);

  // color edges
  if (edge_coloring_method == 'f') {
    // if face colors is none clear edges
    if (!face_coloring_method && !face_color.is_set())
      geom.colors(EDGES).clear();
    else
      // edges take colors from faces
      clrng.e_from_adjacent(FACES);
  }
  else if (edge_coloring_method == 'C') {
    // color by connection
    clrng.e_order(true);
  }
  else
    // use color selected
    clrng.e_one_col(edge_color);

  // color vertices
  if (vertex_coloring_method == 'e') {
    // vertices take color from edges
    clrng.v_from_adjacent(EDGES);
  }
  else if (vertex_coloring_method == 'f') {
    // vertices take color from edges
    clrng.v_from_adjacent(FACES);
  }
  else if (vertex_coloring_method == 'n')
    clrng.v_order(true);
  else
    // use color selected
    clrng.v_one_col(vertex_color);

  // reassert invisible edges
  for (unsigned int i = 0; i < invisible_edges.size(); i++)
    geom.colors(EDGES).set(invisible_edges[i], Color::invisible);

  // if using face connection coloring
  // vertices from kis must be made invisible in geom
  if (face_coloring_method == 'C') {
    for (unsigned int i = 0; i < kis.verts().size(); i++) {
      int v_idx = find_vert_by_coords(geom, kis.verts()[i], anti::epsilon);
      if (v_idx != -1) {
        if ((kis.colors(VERTS).get(i)).is_invisible())
          geom.colors(VERTS).set(v_idx, Color::invisible);
      }
    }
  }

  // set transparency
  if (face_opacity > -1) {
    Status stat = Coloring(&geom).apply_transparency(face_opacity);
    if (stat.is_warning())
      fprintf(stderr, "%s: warning: option -T: %s\n", caller.c_str(),
              stat.msg().c_str());
  }
}

// for canonical.cc and conway.cc

// RK - find nearpoints radius, sets range minimum and maximum
double edge_nearpoints_radius(const Geometry &geom, double &min, double &max,
                              Vec3d &center)
{
  min = std::numeric_limits<double>::max();
  max = std::numeric_limits<double>::min();

  vector<vector<int>> edges;
  geom.get_impl_edges(edges);

  vector<Vec3d> near_pts;

  double nearpt_radius = 0;
  for (auto &edge : edges) {
    Vec3d P = geom.edge_nearpt(edge, Vec3d(0, 0, 0));
    near_pts.push_back(P);

    double l = P.len();
    nearpt_radius += l;
    if (l < min)
      min = l;
    if (l > max)
      max = l;
  }

  center = centroid(near_pts);

  return nearpt_radius / double(edges.size());
}

// RK - wrapper
double edge_nearpoints_radius(const Geometry &geom)
{
  double min = 0;
  double max = 0;
  Vec3d center;
  return edge_nearpoints_radius(geom, min, max, center);
}

// sets radius of geom to average of edge near points radius
void unitize_nearpoints_radius(Geometry &geom)
{
  double avg = edge_nearpoints_radius(geom);
  geom.transform(Trans3d::scale(1 / avg));
}

// return true if maximum vertex radius is radius_range_percent (0.0 to ...)
// greater than minimum vertex radius (visible for canonical.cc)
bool canonical_radius_range_test(const Geometry &geom,
                                 const double radius_range_percent)
{
  GeometryInfo rep(geom);
  rep.set_center(geom.centroid());

  double min = rep.vert_dist_lims().min;
  double max = rep.vert_dist_lims().max;

  // min and max should always be positive, max should always be larger
  return (((max - min) / ((max + min) / 2.0)) > radius_range_percent) ? true
                                                                      : false;
}

// Addition to algorithm by Adrian Rossiter
// Finds the edge near points centroid
Vec3d edge_nearpoints_centroid(Geometry &geom, const Vec3d cent)
{
  vector<vector<int>> edges;
  geom.get_impl_edges(edges);
  Vec3d e_cent(0, 0, 0);
  for (auto &edge : edges)
    e_cent += geom.edge_nearpt(edge, cent);
  return e_cent / double(edges.size());
}

// reciprocalN() is from the Hart's Conway Notation web page
// make array of vertices reciprocal to given planes (face normals)
// RK - save of verbatim port code
/*
vector<Vec3d> reciprocalN_old(const Geometry &geom)
{
  const vector<vector<int>> &faces = geom.faces();
  const vector<Vec3d> &verts = geom.verts();

  vector<Vec3d> normals;
  for (const auto &face : faces) {
    Vec3d centroid(0, 0, 0);
    Vec3d normal(0, 0, 0);
    double avgEdgeDist = 0;

    int v1 = face.at(face.size() - 2);
    int v2 = face.at(face.size() - 1);
    for (int v3 : face) {
      centroid += verts[v3];
      // orthogonal() was from the Hart's Conway Notation web page. replacement
      // normal += orthogonal(verts[v1], verts[v2], verts[v3]);
      normal += vcross(verts[v3] - verts[v2], verts[v2] - verts[v1]);
      // tangentPoint() was from Hart's Conway Notation web page. replacement
      // avgEdgeDist += tangentPoint(verts[v1], verts[v2]).len();
      Vec3d d = verts[v2] - verts[v1];
      // prevent division by zero
      // avgEdgeDist += (verts[v1] - ((vdot(d,verts[v1])/d.len2()) * d)).len();
      double vdt;
      if (d[0] == 0 && d[1] == 0 && d[2] == 0)
        vdt = 0;
      else
        vdt = vdot(d, verts[v1]) / d.len2();
      avgEdgeDist += (verts[v1] - (vdt * d)).len(); // tangentPoint without call
      v1 = v2;
      v2 = v3;
    }
    centroid *= 1.0 / face.size();
    normal.to_unit();
    avgEdgeDist /= face.size();

    // reciprocal call replace below:
    // prevent division by zero
    // Vec3d ans = reciprocal(normal * vdot(centroid,normal));
    Vec3d v = normal * vdot(centroid, normal);
    Vec3d ans;
    if (v[0] == 0 && v[1] == 0 && v[2] == 0)
      ans = v;
    else {
      ans = v * 1.0 / v.len2();
      ans *= (1 + avgEdgeDist) / 2;
    }
    normals.push_back(ans);
  }

  return normals;
}
*/

/* RK - save of simplified tangent code
      Vec3d d = verts[v2] - verts[v1];
      double vdt = 0;
      // prevent division by zero
      if (d[0] != 0 || d[1] != 0 || d[2] != 0)
        vdt = vdot(d, verts[v1]) / d.len2();
      avgEdgeDist += (verts[v1] - (vdt * d)).len();
*/

/*
// RK - Gives the same answer as the built in face_norm().unit()
// even when nonplanar and measuring all edges
Vec3d face_norm_newell(const Geometry &geom, vector<int> &face)
{
  const vector<Vec3d> &v = geom.verts();

  Vec3d face_normal(0, 0, 0);

  unsigned int sz = face.size();
  for (unsigned int i = 0; i < sz; i++) {
    int v1 = face[i];
    int v2 = face[(i + 1) % sz];

    face_normal[0] += (v[v1][1] - v[v2][1]) * (v[v1][2] + v[v2][2]);
    face_normal[1] += (v[v1][2] - v[v2][2]) * (v[v1][0] + v[v2][0]);
    face_normal[2] += (v[v1][0] - v[v2][0]) * (v[v1][1] + v[v2][1]);
  }

  return face_normal.to_unit();
}

// return the unit normal of all perimeter triangles
Vec3d face_norm_newell(const Geometry &geom, const int f_idx)
{
  vector<int> face = geom.faces(f_idx);
  return face_norm_newell(geom, face);
}
*/

// reciprocalN() is from the Hart's Conway Notation web page
// make array of vertices reciprocal to given planes (face normals)
// RK - has accuracy issues and will have trouble with -l 16
vector<Vec3d> reciprocalN(const Geometry &geom)
{
  vector<Vec3d> normals;

  for (const auto &face : geom.faces()) {
    // RK - the algoritm was written to use triangles for measuring
    // non-planar faces. Now method can be chosen
    Vec3d face_normal = face_norm(geom.verts(), face).unit();
    Vec3d face_centroid = anti::centroid(geom.verts(), face);
    // make sure face_normal points outward
    if (vdot(face_normal, face_centroid) < 0)
      face_normal *= -1.0;

    // RK - find the average lenth of the edge near points
    unsigned int sz = face.size();
    double avgEdgeDist = 0;
    for (unsigned int j = 0; j < sz; j++) {
      int v1 = face[j];
      int v2 = face[(j + 1) % sz];

      avgEdgeDist += geom.edge_nearpt(make_edge(v1, v2), Vec3d(0, 0, 0)).len2();
    }

    // RK - sqrt of length squared here
    avgEdgeDist = sqrt(avgEdgeDist / sz);

    // the face normal height set to intersect face at v
    Vec3d v = face_normal * vdot(face_centroid, face_normal);

    // adjust v to the reciprocal value
    Vec3d ans = v;
    // prevent division by zero
    if (v[0] != 0 || v[1] != 0 || v[2] != 0)
      ans = v * 1.0 / v.len2();

    // edge correction (of v based on all edges of the face)
    ans *= (1 + avgEdgeDist) / 2;

    normals.push_back(ans);
  }

  return normals;
}

// reciprocate on face centers dividing by magnitude squared
vector<Vec3d> reciprocalC_len2(const Geometry &geom)
{
  vector<Vec3d> centers;
  geom.face_cents(centers);
  for (auto &center : centers)
    center /= center.len2();
  return centers;
}

// Implementation of George Hart's planarization and canonicalization algorithms
// http://www.georgehart.com/virtual-polyhedra/conway_notation.html
bool canonicalize_bd(Geometry &base, IterationControl it_ctrl,
                     double radius_range_percent, const bool planarize_only)
{
  bool completed = false;
  it_ctrl.set_finished(false);

  Geometry dual;
  // the dual's initial vertex locations are immediately overwritten
  get_dual(dual, base, 1);
  dual.clear_cols();

  double test_val = it_ctrl.get_test_val();
  double max_diff2 = 0;

  for (it_ctrl.start_iter(); !it_ctrl.is_done(); it_ctrl.next_iter()) {
    vector<Vec3d> base_verts_last = base.verts();

    if (planarize_only) {
      // adjust vertices with side effect of planarization. len2() version
      dual.raw_verts() = reciprocalC_len2(base);
      base.raw_verts() = reciprocalC_len2(dual);
      // move centroid to origin for balance
      base.transform(Trans3d::translate(-centroid(base.verts())));
    }
    else {
      // base/dual canonicalize method
      dual.raw_verts() = reciprocalN(base);
      base.raw_verts() = reciprocalN(dual);
      // re-center for drift
      base.transform(
          Trans3d::translate(-edge_nearpoints_centroid(base, Vec3d(0, 0, 0))));
    }

    // reduces size imbalance problem with this algorithm
    unitize_nearpoints_radius(base);

    string finish_msg;
    if (it_ctrl.is_status_check_iter()) {
      // len2() for difference value to minimize internal sqrt() calls
      max_diff2 = 0;
      for (unsigned int i = 0; i < base.verts().size(); i++) {
        double diff2 = (base.verts(i) - base_verts_last[i]).len2();
        if (diff2 > max_diff2)
          max_diff2 = diff2;
      }

      // break out if volume goes to zero, restore last good vertices
      if (std::isnan(base.verts(0)[0])) {
        it_ctrl.set_finished();
        finish_msg = "error, volume went to zero";
        base.raw_verts() = base_verts_last;
        radius_range_percent = 0; // bypass range test
      }
      else if (sqrt(max_diff2) < test_val) {
        completed = true;
        it_ctrl.set_finished();
        finish_msg = "solved, test value achieved";
      }
      else if (it_ctrl.is_last_iter()) {
        // reached last iteration without solving
        it_ctrl.set_finished();
        finish_msg = "not solved, test value not achieved";
      }

      // check if radius is expanding or contracting unreasonably,
      // but only for the purpose of finishing early
      // if minimum and maximum radius are differing, the polyhedron is
      // crumpling
      if (radius_range_percent &&
          canonical_radius_range_test(base, radius_range_percent)) {
        if (!it_ctrl.is_finished())
          it_ctrl.set_finished();
        finish_msg = "breaking out: radius range detected. try increasing -d";
      }
    }

    if (it_ctrl.is_status_report_iter()) {
      if (it_ctrl.is_finished())
        it_ctrl.print("Final iteration (%s):\n", finish_msg.c_str());

      it_ctrl.print("%-12u max_diff:%17.15e\n", it_ctrl.get_current_iter(),
                    sqrt(max_diff2));
    }
  }

  return completed;
}

// RK - wrapper for basic planarization with base/dual algorithm
// meant to be called with finite num_iters (not -1)
bool planarize_bd(Geometry &geom, IterationControl it_ctrl)
{
  double radius_range_percent = 0;
  bool planarize_only = true;
  return canonicalize_bd(geom, it_ctrl, radius_range_percent, planarize_only);
}

/*
// plane a single face aligned to the z axis
void plane_face(Geometry &polygon)
{
  Vec3d face_normal = polygon.face_norm(0).unit();
  Vec3d face_centroid = polygon.face_cent(0);
  // make sure face_normal points outward
  if (vdot(face_normal, face_centroid) < 0)
    face_normal *= -1.0;

  // this gives the same results (from the mathematica algorithm)
  // for (auto &vert : polygon.raw_verts())
  //  vert += vdot(face_normal, face_centroid - vert) * face_normal;
  // return;

  // rotate face to z axis
  Trans3d trans = Trans3d::rotate(face_normal, Vec3d(0, 0, 1));
  polygon.transform(trans);

  // refresh face centroid
  face_centroid = polygon.face_cent(0);

  // set z of all vertices to height of face centroid
  for (auto &vert : polygon.raw_verts())
    vert[2] = face_centroid[2];

  // rotate face back to original position
  polygon.transform(trans.inverse());
}
*/

// P and Q are modified
void move_line_to_point(Vec3d &P, Vec3d &Q, const Vec3d &X)
{
  Vec3d Y = X + (Q - P);
  Vec3d V = P + X;
  Vec3d P2 = lines_intersection(P, V, X, Y, 0);
  if (P2.is_set()) {
    Q += P2 - P;
    P = P2;
  }
}

// RK - edge near points of base seek 1
bool canonicalize_unit(Geometry &geom, IterationControl it_ctrl,
                       const double radius_range_percent,
                       const bool planarize_only)
{
  bool completed = false;
  it_ctrl.set_finished(false);

  vector<vector<int>> edges;
  geom.get_impl_edges(edges);

  vector<Vec3d> &verts = geom.raw_verts();

  double test_val = it_ctrl.get_test_val();
  double max_diff2 = 0;

  for (it_ctrl.start_iter(); !it_ctrl.is_done(); it_ctrl.next_iter()) {
    vector<Vec3d> verts_last = verts;

    if (!planarize_only) {
      vector<Vec3d> near_pts;
      for (auto &edge : edges) {
        // unit near point
        Vec3d P = geom.edge_nearpt(edge, Vec3d(0, 0, 0)).unit();
        near_pts.push_back(P);
        move_line_to_point(verts[edge[0]], verts[edge[1]], P);
      }

      // re-center for drift
      Vec3d cent_near_pts = centroid(near_pts);
      for (unsigned int i = 0; i < verts.size(); i++)
        verts[i] -= cent_near_pts;
    }

    for (unsigned int f = 0; f < geom.faces().size(); f++) {
      /*
            // give polygon its own geom. face index needs to reside in vector
            vector<int> face_idxs(1);
            face_idxs[0] = f;
            Geometry polygon = faces_to_geom(geom, face_idxs);
            plane_face(polygon);

            // map vertices back into original geom
            // the numerical order of vertex list in polygon geom is preserved
            vector<int> v_idx;
            for (int v : geom.faces(f))
              v_idx.push_back(v);
            sort(v_idx.begin(), v_idx.end());
            int j = 0;
            for (int v : v_idx)
              verts[v] = polygon.verts(j++);
      */
      // RK - this does formulaically what the above does by brute force
      Vec3d face_normal = face_norm(geom.verts(), geom.faces(f)).unit();
      Vec3d face_centroid = geom.face_cent(f);
      // make sure face_normal points outward
      if (vdot(face_normal, face_centroid) < 0)
        face_normal *= -1.0;
      // place a planar vertex over or under verts[v]
      // adds or subtracts it to get to the planar verts[v]
      for (int v : geom.faces(f))
        verts[v] += vdot(face_normal, face_centroid - verts[v]) * face_normal;
    }

    string finish_msg;
    if (it_ctrl.is_status_check_iter()) {
      // len2() for difference value to minimize internal sqrt() calls
      max_diff2 = 0;
      for (unsigned int i = 0; i < verts.size(); i++) {
        double diff2 = (verts[i] - verts_last[i]).len2();
        if (diff2 > max_diff2)
          max_diff2 = diff2;
      }

      if (sqrt(max_diff2) < test_val) {
        completed = true;
        it_ctrl.set_finished();
        finish_msg = "solved, test value achieved";
      }
      else if (it_ctrl.is_last_iter()) {
        // reached last iteration without solving
        it_ctrl.set_finished();
        finish_msg = "not solved, test value not achieved";
      }

      // check if radius is expanding or contracting unreasonably,
      // but only for the purpose of finishing early
      // if minimum and maximum radius are differing, the polyhedron is
      // crumpling
      if (radius_range_percent &&
          canonical_radius_range_test(geom, radius_range_percent)) {
        if (!it_ctrl.is_finished())
          it_ctrl.set_finished();
        finish_msg = "breaking out: radius range detected. try increasing -d";
      }
    }

    if (it_ctrl.is_status_report_iter()) {
      if (it_ctrl.is_finished())
        it_ctrl.print("Final iteration (%s):\n", finish_msg.c_str());

      it_ctrl.print("%-12u max_diff:%17.15e\n", it_ctrl.get_current_iter(),
                    sqrt(max_diff2));
    }
  }

  return completed;
}

// RK - wrapper for basic planarization with unit algorithm
// meant to be called with finite num_iters (not -1)
bool planarize_unit(Geometry &geom, IterationControl it_ctrl)
{
  double radius_range_percent = 0;
  bool planarize_only = true;
  return canonicalize_unit(geom, it_ctrl, radius_range_percent, planarize_only);
}
