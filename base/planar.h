/*
   Copyright (c) 2017, Roger Kaufman

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
   Name: planar.h
   Description: shared code originally from src/planar.cc
   Project: Antiprism - http://www.antiprism.com
*/

#ifndef PLANAR_H
#define PLANAR_H

#include "geometry.h"
#include "geometryinfo.h"
#include "geometryutils.h"
#include "vec3d.h"

using std::string;
using std::vector;
using std::map;

namespace anti {

// from planar.cc

/// Make a geometry from a specified set of faces
/**\param geom the geometry.
 * \param face_idxs the index numbers of the faces to use
 * \return A geometry made of those faces. */
Geometry faces_to_geom(const Geometry &geom, const std::vector<int> &face_idxs);

/// add a vector P into the geom unless a point already occupies that point
/**\param geom the geometry.
 * \param P a point.
 * \param vcol color of the new point.
 * \param eps value for contolling the limit of precision.
 * \return the index of the new point, or the occupying point. */
int vertex_into_geom(Geometry &geom, const Vec3d &P, Color vcol,
                     const double eps);

/// add an edge v1, v2 into the geom unless an edge of v1, v2 already exists
/**\param geom the geometry.
 * \param v_idx1 is the first index.
 * \param v_idx2 is the second index.
 * \param ecol color of the new edge.
 * \return the true if new edge was created. */
bool edge_into_geom(Geometry &geom, const int v_idx1, const int v_idx2,
                    Color ecol);

/// add intersection points to all crossing edges. The input contains on faces
/**\param geom the geometry.
 * \param eps value for contolling the limit of precision.
 * \return the true if original edges were replaced. */
bool mesh_edges(Geometry &geom, const double eps);

/// Project a normal from 3D -> 2D
/**\param normal input normal.
 * \param idx will be the dimension not included (x:0 y:1 z:2).
 * \param sign will be the sign of the projection. */
void project_using_normal(const Vec3d &normal, int &idx, int &sign);

/// fill in faces into a faceless mesh
/**\param geom the geometry.
 * \param planar_merge_type (1: tile, 2: merge).
 * \param nonconvex_faces if any faces are nonconvex, and index list,
 * \param normal of mesh.
 * \param eps value for contolling the limit of precision. */
void fill_in_faces(Geometry &geom, const int planar_merge_type,
                   vector<int> &nonconvex_faces, const Vec3d &normal,
                   double eps = epsilon);

// start of stellate stuff

/* for future use
bool three_plane_intersect(Vec3d Q0, Vec3d n0, Vec3d Q1, Vec3d n1,
                           Vec3d Q2, Vec3d n2, Vec3d &P, double eps=epsilon);

vector<int> neighbor_faces(Geometry &geom, int f_idx);

Vec3d calc_stellation_point(Geometry &geom, int f_idx, double eps=epsilon);
*/

/// make stellation diagram for a face of a geom
/**\param geom the geometry.
 * \param f_idx is face to make diagram.
 * \param sym_string is sub-symmetry of stellation.
 * \param projection_width is length of line extents of diagram.
 * \param eps value for contolling the limit of precision. */
Geometry make_stellation_diagram(Geometry &geom, int f_idx,
                                 string sym_string = "",
                                 int projection_width = 500,
                                 double eps = epsilon);

/// if faces are pinched (revisited vertices) in a geom, split them
/**\param geom the geometry.
 * \param eps value for contolling the limit of precision. */
void split_pinched_faces(Geometry &geom, double eps = epsilon);

/// return lists of index for full stellation diagram
/// Note: index lists first position is reference to the map of stellation
/// diagrams
/**\param diagrams a map of stellation diagrams.
 * \param idx_lists the partial lists.
 * \param remove_multiples if true any duplicates in lists are removed */
vector<vector<int>> lists_full(map<int, Geometry> &diagrams,
                               const vector<vector<int>> &idx_lists,
                               bool remove_multiples = true);

/// return lists which have been standardized
/// Note: index lists first position is reference to the map of stellation
/// diagrams
/**\param geom the geometry.
 * \param sym_string the symmetry symbol.
 * \param diagrams a map of stellation diagrams.
 * \param idx_lists the partial lists.
 * \param idx_lists_full the full lists (from lists_full)
 * \param remove_multiples if true any duplicates in lists are removed */
vector<vector<int>> lists_resolved(const Geometry &geom,
                                   const string &sym_string,
                                   map<int, Geometry> &diagrams,
                                   const vector<vector<int>> &idx_lists,
                                   const vector<vector<int>> &idx_lists_full,
                                   bool remove_multiples = true);

/// take fully connected compound and produce unconnected model
/**\param geom the geometry. */
void rebuild_compound(Geometry &geom);

/// make a stellation
/// Note: index lists first position is reference to the map of stellation
/// diagrams
/**\param geom the geometry.
 * \param diagrams a map of stellation diagrams.
 * \param idx_lists the index lists used
 * \param sym_string is sub-symmetry of stellation.
 * \param merge_faces try to make full faces from adjacent facelets
 * \param remove_inline_verts remove vertices which are on a line between
 *        two other vertices (if merge_faces is true)
 * \param split_pinched split faces with revisited vertices (if merge_faces is
 * true)
 * \param resolve_faces to standard
 * \param remove_multiples if true any duplicates in lists are removed (if
 * resolved_faces is true)
 * \param map_string color map for coloring from diagrams
 * \param eps value for contolling the limit of precision. */
Geometry make_stellation(const Geometry &geom, map<int, Geometry> &diagrams,
                         const vector<vector<int>> &idx_lists,
                         const string &sym_string, bool merge_faces = true,
                         bool remove_inline_verts = true,
                         bool split_pinched = true, bool resolve_faces = false,
                         bool remove_multiples = false,
                         string map_string = "compound", double eps = epsilon);

// RK - functions for winding number

/// Get winding number of a point in a polygon
/**\param polygon geometry containing the polygon
 * \param point point to test
 * \param eps a small number, coordinates differing by less than eps are
 *  the same.
 * \return The winding number. */
int get_winding_number(const Geometry &polygon, const Vec3d &point,
                       double eps = epsilon);

/// Get largest winding number of any point in a polygon
/**\param polygon geometry containing the polygon
 * \param points if not empty, test just these points
 * \param face_normal the face normal
 * \param find_direction \c true find suitable z-orientation for normal,
 *  \c false don't reorient
 * \param eps a small number, coordinates differing by less than eps are
 *  the same.
 * \return The winding number. */
int get_winding_number_polygon(const Geometry &polygon,
                               const std::vector<Vec3d> &points,
                               const Normal &face_normal,
                               bool find_direction = false,
                               const double eps = epsilon);

/// Find the (signed) denominator of a wound polygon
/**\param geom the geometry.
 * \param face_idx face index.
 * \param eps a small number, coordinates differing by less than eps are
 *  the same.
 * \return The signed olygon denominator. */
int find_polygon_denominator_signed(const Geometry &geom, int face_idx,
                                    double eps = epsilon);

} // namespace anti

#endif // PLANAR_H
