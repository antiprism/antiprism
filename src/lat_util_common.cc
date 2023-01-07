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
   Name: lat_util_common.cc
   Description: lat_util code shared by source code in /src
   Project: Antiprism - http://www.antiprism.com
*/

#include "lat_util_common.h"
#include "../base/antiprism.h"

#include <cstdio>
#include <memory>
#include <string>
#include <vector>

using std::set;
using std::string;
using std::vector;

using namespace anti;

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
