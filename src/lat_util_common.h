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
   Name: lat_util_common.h
   Description: lat_util_code shared by source code in /src
   Project: Antiprism - http://www.antiprism.com
*/

#ifndef LAT_UTIL_COMMON_H
#define LAT_UTIL_COMMON_H

#include <cstdio>
#include <memory>
#include <string>
#include <vector>

using std::string;
using std::vector;

#include "../base/antiprism.h"

// for lat_util.cc and bravais.cc

void parse_color_string(const anti::ProgramOpts *, const char *, const char,
                        const string &, vector<anti::Color> &);

double lattice_radius(const anti::Geometry &, const char);

void geom_container_clip(anti::Geometry &, anti::Geometry &, const double,
                         const anti::Vec3d &, const bool,
                         double eps = anti::epsilon);

void geom_spherical_clip(anti::Geometry &, const double, const anti::Vec3d &,
                         const bool, double eps = anti::epsilon);

void list_grid_radii(const string &, const anti::Geometry &,
                     const anti::Vec3d &, int report_type = 1,
                     double eps = anti::epsilon);

void list_grid_struts(const string &, const anti::Geometry &,
                      int report_type = 1, double eps = anti::epsilon);

void add_color_struts(anti::Geometry &, const double, anti::Color &,
                      double eps = anti::epsilon);

void color_centroid(anti::Geometry &, anti::Color &,
                    double eps = anti::epsilon);

int get_voronoi_geom(anti::Geometry &, anti::Geometry &, const bool, const bool,
                     double eps = anti::epsilon);

// for lat_util.cc, bravais.cc and waterman.cc

void convex_hull_report(const anti::Geometry &, const bool);

#endif // LAT_UTIL_COMMON_H
