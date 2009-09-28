/*
   Copyright (c) 2003-2008, Adrian Rossiter

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

/*!\file geom_utils.h
 * \brief Utilities associated with geometries
*/

#ifndef GEOM_UTILS_H
#define GEOM_UTILS_H

void orient_face(vector<int> &face, int v0, int v1);

void triangulate_basic(geom_if &geom, bool sq_diag=true, col_val inv=col_val(),
      vector<int> *fmap=0);

double minimum_distance(const geom_if &geom, double sig_dist=0, char *errmsg=0);

void make_polar_zono(geom_if &zono, int star_n, bool out_star);

bool make_resource_geom(geom_if &geom, string name, char *errmsg=0);

// test points versus hull functions
bool test_points_vs_hull(const vector<vec3d> &, col_geom_v &, bool, bool, bool, double);
bool is_geom_fully_outside_hull(col_geom_v &, col_geom_v &, double);
bool is_geom_outside_hull(col_geom_v &, col_geom_v &, double);
bool is_geom_on_surface_hull(col_geom_v &, col_geom_v &, double);
bool is_geom_inside_hull(col_geom_v &, col_geom_v &, double);
bool is_geom_fully_inside_hull(col_geom_v &, col_geom_v &, double);
bool is_point_fully_outside_hull(const vec3d &, col_geom_v &, double);
bool is_point_outside_hull(const vec3d &, col_geom_v &, double);
bool is_point_on_surface_hull(const vec3d &, col_geom_v &, double);
bool is_point_inside_hull(const vec3d &, col_geom_v &, double);
bool is_point_fully_inside_hull(const vec3d &, col_geom_v &, double);

#endif // GEOM_UTILS_H

