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
   Name: lattice_grid.cc
   Description: grids and lattices with integer coordinates
   Project: Antiprism - http://www.antiprism.com
*/

#ifndef LATTICE_GRID_H
#define LATTICE_GRID_H

#include <string>
#include <vector>

#include "../base/antiprism.h"

using std::vector;

typedef bool (*COORD_TEST_F)(int, int, int);

bool sc_test(int x, int y, int z);                 // dist2 = 1, 2, 3
bool fcc_test(int x, int y, int z);                // dist2 = 2, 4, 6, 12
bool bcc_test(int x, int y, int z);                // dist2 = 3, 4, 8
bool hcp_test(int x, int y, int z);                // dist2 = 18
bool rh_dodec_test(int x, int y, int z);           // dist2 = 3 (8) ,
bool cubo_oct_test(int x, int y, int z);           // dist2 = 2
bool tr_oct_test(int x, int y, int z);             // dist2 = 2
bool tr_tet_tet_test(int x, int y, int z);         // dist2 = 2
bool tr_tet_tr_oct_cubo_test(int x, int y, int z); // dist2 = 4
bool diamond_test(int x, int y, int z);            // dist2 = 3
bool k_4_test(int x, int y, int z);                // dist2 = 2
bool hcp_diamond_test(int x, int y, int z);        // dist2 = 27

void add_struts(anti::Geometry &geom, int len2);

// for lattice code only
void parse_color_string(const anti::ProgramOpts *, const char *, const char,
                        const std::string &, vector<anti::Color> &);

double lattice_radius(const anti::Geometry &, const char);
void geom_container_clip(anti::Geometry &, anti::Geometry &, const double,
                         const anti::Vec3d &, const bool,
                         double eps = anti::epsilon);
void geom_spherical_clip(anti::Geometry &, const double, const anti::Vec3d &,
                         const bool, double eps = anti::epsilon);
void list_grid_radii(const std::string &, const anti::Geometry &,
                     const anti::Vec3d &, int report_type = 1,
                     double eps = anti::epsilon);
void list_grid_struts(const std::string &, const anti::Geometry &,
                      int report_type = 1, double eps = anti::epsilon);
void add_color_struts(anti::Geometry &, const double, anti::Color &,
                      double eps = anti::epsilon);
void color_centroid(anti::Geometry &, anti::Color &,
                    double eps = anti::epsilon);

// color functions for lattice programs
void color_by_symmetry_normals(anti::Geometry &, const char, const int,
                               double eps = anti::epsilon);
void color_edges_by_sqrt(anti::Geometry &, const char);

// convex hull and voronoi wrappers
void convex_hull_report(const anti::Geometry &, const bool);
int get_voronoi_geom(anti::Geometry &, anti::Geometry &, const bool, const bool,
                     double eps = anti::epsilon);

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
  // void add_struts(Geometry &geom, int len2);
};

class sph_lat_grid : public int_lat_grid {
public:
  sph_lat_grid() {}
  virtual void set_o_width(double w) { o_width = w; }
  virtual void set_i_width(double w) { i_width = w; }
  virtual void make_lattice(anti::Geometry &geom);
};

#endif // LATTICE_GRID_H
