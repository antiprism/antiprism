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

/*!\file vec_utils.h
   \brief Vector utilities
*/


#ifndef VEC_UTILS_H
#define VEC_UTILS_H

#include <vector>
#include "vec3d.h"

using std::vector;

///Get the centroid of a set of points
/**\param pts the points
 * \param idxs the index numbers of the points to use, or if none (the default)
 * then use all the points.
 * \return The centroid. */
vec3d centroid(const vector<vec3d> &pts, const vector<int> &idxs = vector<int>());


///Get the point of intersection of a line and a plane.
/**\param Q a point on the plane.
 * \param n the normal to the plane
 * \param P0 a point on the line
 * \param P1 a second point on the line
 * \param where to return the the intersection point was relative to the
 * two points on the line, the meaning of the values is
 * <ul>
 * <li>0 - outside P0
 * <li>1 - on P0
 * <li>2 - between P0 and P1
 * <li>3 - on P1
 * <li>4 - outside P1
 * </ul>
 * \return the point of intersection (it will be unset if the line doesn't
 * intersect the plane. */
vec3d line_plane_intersect(vec3d Q, vec3d n, vec3d P0, vec3d P1, int *where=0);

///Get the point of intersection of a line and a plane.
/**\param Q0 a point on the plane.
 * \param Q1 a second point on the plane.
 * \param Q2 a third points on the plane, not on the line Q0Q1.
 * \param P0 a point on the line
 * \param P1 a second point on the line
 * \param where to return the the intersection point was relative to the
 * two points on the line, the meaning of the values is
 * <ul>
 * <li>0 - outside P0
 * <li>1 - on P0
 * <li>2 - between P0 and P1
 * <li>3 - on P1
 * <li>4 - outside P1
 * </ul>
 * \return The point of intersection (it will be unset if the line doesn't
 * intersect the plane*/
vec3d line_plane_intersect(vec3d Q0, vec3d Q1, vec3d Q2, vec3d P0, vec3d P1, int *where=0);

///Get the line where two planes intersect.
/**\param Q0 a point on the first plane.
 * \param n0 the normal to the first plane.
 * \param Q1 a point on the second plane.
 * \param n1 the normal to the second plane.
 * \param P to return a point on the line of intersection.
 * \param dir to return the direction of the line of intersection.
 * \return \c true if the planes intersect, otherwise \c false. */
bool two_plane_intersect(vec3d Q0, vec3d n0, vec3d Q1, vec3d n1,
      vec3d &P, vec3d &dir);

//unused
bool three_plane_intersect(vec3d Q0, vec3d n0, vec3d Q1, vec3d n1,
      vec3d Q2, vec3d &n2, vec3d &P);

///Get the nearest point on a line to another (skew) line
/**\param P0 a point on the first line.
 * \param P1 another point on the first line.
 * \param Q0 a point on the second line.
 * \param Q1 another point on the second line.
 * \param P to return the nearest point on the first line.
 * \param Q to return the nearest point on the second line.
 * \return \c true if there is one nesarpoint per line, otherwise \c false (the
 * lines are parallel.) */
bool lines_nearest_points(vec3d P0, vec3d P1, vec3d Q0, vec3d Q1,
      vec3d &P, vec3d &Q);

///Get the intersection point of two lines
/**If the lines do not exactly intersect then the point returned
 * is half way along the shortest segment that joins the lines.
 * \param P a point on the first line.
 * \param P_dir the direction of the first line.
 * \param Q a point on the second line.
 * \param Q_dir the direction of the second line.
 * \return The point of intersection. */
vec3d lines_intersection(vec3d P, vec3d P_dir, vec3d Q, vec3d Q_dir);


///Get the nearest point on a line to a particular point.
/**\param P a point.
 * \param Q0 a point on the line.
 * \param Q1 another point on the line.
 * \return The point on the line through Q0 and Q1 that is closest to P. */
inline vec3d nearest_point(vec3d P, vec3d Q0, vec3d Q1)
{
   vec3d u = (Q1 - Q0).unit();
   return Q0 + vdot(u, P - Q0)*u;
}

///Get the nearest point on a plane to a particular point.
/**\param P a point.
 * \param Q0 a point on the plane.
 * \param Q1 a second point on the plane.
 * \param Q2 a third points on the plane, not on the line Q0Q1.
 * \return The point on the plane through Q0, Q1 and Q2 that is closest to P. */
vec3d nearest_point(vec3d P, vec3d Q0, vec3d Q1, vec3d Q2);

///Get the nearest point on a space to a particular point.
/**\param P a point.
 * \param points independant points determining the space (one, two
 * or three points as the space is a point, line or plane.)
 * \return The point on the space that is closest to P. */
vec3d nearest_point(vec3d P, const vector <vec3d> &points);

///Get the nearest point on a space to a particular point.
/**\param P a point.
 * \param points a set of points.
 * \param idxs the index numbers of independant points from \a points that
 * determine the space (one, two or three index numbers as the space
 * is a point, line or plane.)
 * \return The point on the space that is closest to P. */
vec3d nearest_point(vec3d P, const vector <vec3d> &points,
      const vector<int> &idxs);

///Get a face normal and face area
/**\param verts a set of vertices
 * \param face the index numbers of the vertices in \a verts that make the face.
 * \param allow_zero if \c true then the length of the returned normal
 * is the area of the face, if \c false then this will not be true for
 * faces with a signed area close to zero.
 * \return A normal to the face. */
vec3d face_norm(const vector<vec3d> &verts, const vector<int> &face, bool allow_zero=false);

#endif // VEC_UTILS_H

