/*
   Copyright (c) 2003, Adrian Rossiter

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

#include <vector>
#include "../base/geom.h"
using std::vector;

typedef bool (* COORD_TEST_F)(int, int, int);

bool sc_test(int x, int y, int z);                  // dist2 = 1, 2, 3
bool fcc_test(int x, int y, int z);                 // dist2 = 2, 4, 6, 12
bool bcc_test(int x, int y, int z);                 // dist2 = 3, 4, 8
bool hcp_test(int x, int y, int z);                 // dist2 = 18
bool rh_dodec_test(int x, int y, int z);            // dist2 = 3 (8) , 
bool cubo_oct_test(int x, int y, int z);            // dist2 = 2
bool tr_oct_test(int x, int y, int z);              // dist2 = 2
bool tr_tet_tet_test(int x, int y, int z);          // dist2 = 2
bool tr_tet_tr_oct_cubo_test(int x, int y, int z);  // dist2 = 4
bool diamond_test(int x, int y, int z);             // dist2 = 3
bool k_4_test(int x, int y, int z);                 // dist2 = 2
bool hcp_diamond_test(int x, int y, int z);         // dist2 = 27

void add_struts(geom_if &geom, int len2);

// for lattice code only
double lattice_radius(const geom_if &, char);
void geom_container_clip(col_geom_v &, col_geom_v &, double, vec3d, double);
void geom_spherical_clip(col_geom_v &, double, vec3d, double);
void list_grid_radii(const col_geom_v &, vec3d, double);
void list_grid_struts(const col_geom_v &, double);
void add_color_struts(col_geom_v &, double, col_val);
void color_centroid(col_geom_v &, col_val, double);

// color functions for lattice programs
void color_by_symmetry_normals(col_geom_v &, char, int);
void color_edges_by_sqrt(col_geom_v &, char);

// convex hull and voronoi wrappers
void convex_hull_report(const geom_v &, bool);
int get_voronoi_geom(col_geom_v &, col_geom_v &, bool, bool, double);


class int_lat_grid
{
   protected:
      double o_width;
      double i_width;
      vec3d centre;
      COORD_TEST_F coord_test;
   public:
      //enum { l_sc, l_fcc, l_bcc, l_rh_dodec, l_cubo_oct,
      //   l_tr_oct, l_tr_tet_tet, l_tr_oct_tr_tet_cubo, l_diamond }
      int_lat_grid() {}
      virtual ~int_lat_grid() {}
      void set_coord_test(COORD_TEST_F func) {coord_test=func;}
      virtual void set_o_width(double w)   {o_width = w;}
      virtual void set_i_width(double w)   {i_width = w;}
      virtual void set_centre(vec3d cent)   {centre = cent;}
      virtual void make_lattice(geom_if &geom);
      //void add_struts(geom_if &geom, int len2);
};

class sph_lat_grid: public int_lat_grid
{
   public:
      sph_lat_grid() {}
      virtual void set_o_width(double w)   {o_width = w;}
      virtual void set_i_width(double w)   {i_width = w;}
      virtual void make_lattice(geom_if &geom);
};
   
#endif // LATTICE_GRID_H

