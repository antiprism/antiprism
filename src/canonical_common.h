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
   Name: canonical_common.h
   Description: canonical code shared in /src
   Project: Antiprism - http://www.antiprism.com
*/

#ifndef CANONICAL_COMMON_H
#define CANONICAL_COMMON_H

#include <cstdio>
#include <string>
#include <vector>

using std::string;
using std::vector;

#include "../base/antiprism.h"

// for lat_grid.cc

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

void color_by_symmetry_normals(anti::Geometry &, const char, const int,
                               double eps = anti::epsilon);

void color_edges_by_sqrt(anti::Geometry &, const char);

void convex_hull_report(const anti::Geometry &, const bool);

// for stellate.cc and miller.cc

void color_stellation(anti::Geometry &geom, const char, const char, const char,
                      const anti::Color &, const anti::Color &,
                      const anti::Color &, const int, const string &,
                      const string);

// for canonical.cc and conway.cc

/// find nearpoints radius, sets range minimum and maximum
/**\param geom geometry.
 * \param min returns the minimum nearpoints radius.
 * \param max returns the maximum nearpoints radius.
 * \param center returns the cente rof the nearpoints.
 * \returns the average radius of the nearpoints. */
double edge_nearpoints_radius(const anti::Geometry &geom, double &min,
                              double &max, anti::Vec3d &center);

/// wrapper for above.
/**\param geom geometry. */
double edge_nearpoints_radius(const anti::Geometry &geom);

/// sets radius of geom to average of edge near points radius
/**\param geom geometry. */
void unitize_nearpoints_radius(anti::Geometry &geom);

/// return true if maximum vertex radius is radius_range_percent (0.0 to ...)
/**greater than minimum vertex radius (visible for canonical.cc)
 * \param geom geometry to measure.
 * \param radius_range_percent limit to maximum radius over minimum radius */
bool canonical_radius_range_test(const anti::Geometry &geom,
                                 const double radius_range_percent);

/// returns the edge near points centroid
/**\param geom geometry to measure
 * \param cent centre from which to calculate nearpoints on edges
 * \return the centroid of the nearpoints. */
anti::Vec3d edge_nearpoints_centroid(anti::Geometry &geom,
                                     const anti::Vec3d cent = anti::Vec3d(0, 0,
                                                                          0));

/// Canonicalize (George Hart "Conway Notation" algorithm)
/**See http://www.georgehart.com/virtual-polyhedra/conway_notation.html
 * \param base geometry to canonicalise.
 * \param it_ctrl interation control.
 * \param radius_range_percent if the model outer radius increases this
 *  much over the inner radius then it is growing too much, terminate.
 * \param planarize_only planarise only.
 * \return \c true if success, otherwise \c false */
bool canonicalize_bd(anti::Geometry &base, anti::IterationControl it_ctrl,
                     double radius_range_percent, const bool planarize_only);

/// an abbreviated wrapper for planarization with the base/dual method
/**\param geom geometry to planarize.
 * \param it_ctrl interation control.
 * \return \c true if success, otherwise \c false */
bool planarize_bd(anti::Geometry &geom, anti::IterationControl it_ctrl);

/// RK - edge near points of base seek 1
/**\param geom geometry to canonicalise.
 * \param it_ctrl interation control.
 * \param radius_range_percent if the model outer radius increases this
 *  much over the inner radius then it is growing too much, terminate.
 * \param planarize_only planarise only.
 * \return \c true if success, otherwise \c false */
bool canonicalize_unit(anti::Geometry &geom, anti::IterationControl it_ctrl,
                       const double radius_range_percent,
                       const bool planarize_only);

/// an abbreviated wrapper for planarization with canonicalize_unit
/**\param geom geometry to planarize.
 * \param it_ctrl interation control
 * \return \c true if success, otherwise \c false */
bool planarize_unit(anti::Geometry &geom, anti::IterationControl it_ctrl);

#endif // CANONICAL_COMMON_H
