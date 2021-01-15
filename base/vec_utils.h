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

/*!\file vec_utils.h
   \brief Vector utilities
*/

#ifndef VEC_UTILS_H
#define VEC_UTILS_H

#include "vec3d.h"

#include <vector>

namespace anti {

/// Get the centroid of a set of points
/**\param pts the points
 * \param idxs the index numbers of the points to use, or if none (the default)
 *  then use all the points.
 * \return The centroid. */
Vec3d centroid(const std::vector<Vec3d> &pts,
               const std::vector<int> &idxs = std::vector<int>());

/// Get the point of intersection of a line and a plane.
/**\param Q a point on the plane.
 * \param n the normal to the plane
 * \param P0 a point on the line
 * \param P1 a second point on the line
 * \param where to return the the intersection point was relative to the
 *  two points on the line, the meaning of the values is
 *  <ul>
 *  <li>0 - outside P0
 *  <li>1 - on P0
 *  <li>2 - between P0 and P1
 *  <li>3 - on P1
 *  <li>4 - outside P1
 *  </ul>
 * \param eps value for contolling the limit of precision
 * \return the point of intersection (it will be unset if the line doesn't
 *  intersect the plane. */
Vec3d line_plane_intersect(Vec3d Q, Vec3d n, Vec3d P0, Vec3d P1,
                           int *where = nullptr, double eps = epsilon);

/// Get the point of intersection of a line and a plane.
/**\param Q0 a point on the plane.
 * \param Q1 a second point on the plane.
 * \param Q2 a third points on the plane, not on the line Q0Q1.
 * \param P0 a point on the line
 * \param P1 a second point on the line
 * \param where to return the the intersection point was relative to the
 *  two points on the line, the meaning of the values is
 *  <ul>
 *  <li>0 - outside P0
 *  <li>1 - on P0
 *  <li>2 - between P0 and P1
 *  <li>3 - on P1
 *  <li>4 - outside P1
 *  </ul>
 * \param eps value for contolling the limit of precision
 * \return The point of intersection (it will be unset if the line doesn't
 *  intersect the plane*/
Vec3d line_plane_intersect(Vec3d Q0, Vec3d Q1, Vec3d Q2, Vec3d P0, Vec3d P1,
                           int *where = nullptr, double eps = epsilon);

/// Get the line where two planes intersect.
/**\param Q0 a point on the first plane.
 * \param n0 the normal to the first plane.
 * \param Q1 a point on the second plane.
 * \param n1 the normal to the second plane.
 * \param P to return a point on the line of intersection.
 * \param dir to return the direction of the line of intersection.
 * \param eps value for contolling the limit of precision
 * \return \c true if the planes intersect, otherwise \c false. */
bool two_plane_intersect(Vec3d Q0, Vec3d n0, Vec3d Q1, Vec3d n1, Vec3d &P,
                         Vec3d &dir, double eps = epsilon);

/// Get the point where three planes intersect.
/**\param Q0 a point on the first plane.
 * \param n0 the normal to the first plane.
 * \param Q1 a point on the second plane.
 * \param n1 the normal to the second plane.
 * \param Q2 a point on the second plane.
 * \param n2 the normal to the second plane.
 * \param P to return the point of intersection.
 * \param eps value for contolling the limit of precision
 * \return \c true if the planes intersect, otherwise \c false. */
bool three_plane_intersect(Vec3d Q0, Vec3d n0, Vec3d Q1, Vec3d n1, Vec3d Q2,
                           Vec3d n2, Vec3d &P, double eps = epsilon);

/// Get the nearest point on a line to another (skew) line
/**\param P0 a point on the first line.
 * \param P1 another point on the first line.
 * \param Q0 a point on the second line.
 * \param Q1 another point on the second line.
 * \param P to return the nearest point on the first line.
 * \param Q to return the nearest point on the second line.
 * \param eps value for contolling the limit of precision
 * \return \c true if there is one nearpoint per line, otherwise \c false (the
 *  lines are parallel.) */
bool lines_nearest_points(Vec3d P0, Vec3d P1, Vec3d Q0, Vec3d Q1, Vec3d &P,
                          Vec3d &Q, double eps = epsilon);

/// Get the point of intersection between two lines
/**\param P0 a point on the first line.
 * \param P1 another point on the first line.
 * \param Q0 a point on the second line.
 * \param Q1 another point on the second line.
 * \param eps value for controlling the limit of precision. If eps is 0 then
 *  always return an intersection point, even if lines are skew or parallel.
 * \return The point of intersection (it will be unset if eps>0 and the point
 *  is not within distance eps from both lines*/
Vec3d lines_intersection(const Vec3d &P0, const Vec3d &P1, const Vec3d &Q0,
                         const Vec3d &Q1, double eps = epsilon);

/// Get the point of intersection of two segments (or best approximation if
/// segments are skew).
/**\param P0 one end of first segment.
 * \param P1 other end of first segment.
 * \param Q0 one end of second segment.
 * \param Q1 other end of second segment.
 * \param eps value for contolling the limit of precision
 * \return The intersection point, or unset if no intersection point.*/
Vec3d segments_intersection(const Vec3d &P0, const Vec3d &P1, const Vec3d &Q0,
                            const Vec3d &Q1, double eps = epsilon);

/// Get the nearest point on a line to a particular point.
/**\param P a point.
 * \param Q0 a point on the line.
 * \param Q1 another point on the line.
 * \return The point on the line through Q0 and Q1 that is closest to P. */
inline Vec3d nearest_point(Vec3d P, Vec3d Q0, Vec3d Q1)
{
  Vec3d u = (Q1 - Q0).unit();
  return Q0 + vdot(u, P - Q0) * u;
}

/// Get the nearpoint on a plane
/**\param P a point on the plane, usually a vertex
 * \param point_on_plane a point on the plane, usually the centroid
 * \param unit_norm unit normal of plane
 * \return nearpoint vector */
inline Vec3d nearpoint_on_plane(const Vec3d &P, const Vec3d &point_on_plane,
                                const Vec3d &unit_norm)
{
  return P + vdot(point_on_plane - P, unit_norm) * unit_norm;
}

/// Get the nearest point on a plane to a particular point.
/**\param P a point.
 * \param Q0 a point on the plane.
 * \param Q1 a second point on the plane.
 * \param Q2 a third points on the plane, not on the line Q0Q1.
 * \param eps value for contolling the limit of precision
 * \return The point on the plane through Q0, Q1 and Q2 that is closest to P. */
Vec3d nearest_point(Vec3d P, Vec3d Q0, Vec3d Q1, Vec3d Q2,
                    double eps = epsilon);

/// Get the nearest point on a space to a particular point.
/**\param P a point.
 * \param points independant points determining the space (one, two
 *  or three points as the space is a point, line or plane.)
 * \param eps value for contolling the limit of precision
 * \return The point on the space that is closest to P. */
Vec3d nearest_point(Vec3d P, const std::vector<Vec3d> &points,
                    double eps = epsilon);

/// Get the nearest point on a space to a particular point.
/**\param P a point.
 * \param points a set of points.
 * \param idxs the index numbers of independant points from \a points that
 *  determine the space (one, two or three index numbers as the space
 *  is a point, line or plane.)
 * \param eps value for contolling the limit of precision
 * \return The point on the space that is closest to P. */
Vec3d nearest_point(Vec3d P, const std::vector<Vec3d> &points,
                    const std::vector<int> &idxs, double eps = epsilon);

/// Check if a point lies on a segment.
/**\param P a point.
 * \param Q0 one end of the segment.
 * \param Q1 the other end of the segment.
 * \param eps value for contolling the limit of precision
 * \return \c true if the point is on the segment, otherwise \c false. */
bool in_segment(const Vec3d &P, const Vec3d &Q0, const Vec3d &Q1,
                double eps = epsilon);

/// Adjust a point on a segment to actually lie on the segment.
/**\param P a point.
 * \param Q0 one end of the segment.
 * \param Q1 the other end of the segment.
 * \param eps value for contolling the limit of precision
 * \return The point adjusted to lie on the segment, unset if the point was
 *  not on the segment */
Vec3d point_in_segment(const Vec3d &P, const Vec3d &Q0, const Vec3d &Q1,
                       double eps = epsilon);

/// Get a face normal and face area
/**\param verts a set of vertices
 * \param face the index numbers of the vertices in \a verts that make the face.
 * \param allow_zero if \c true then the length of the returned normal
 *  is the area of the face, if \c false then this will not be true for
 *  faces with a signed area close to zero.
 * \return A normal to the face. */
Vec3d face_norm(const std::vector<Vec3d> &verts, const std::vector<int> &face,
                bool allow_zero = false);

/// Get the angle required to rotate one vector onto another around an axis
/**\param v0 vector to rotate (perpendicular to axis)
 * \param v1 vector to rotate onto (perpendicular to axis)
 * \param axis axis to rotate around (perpendicular to v0 and v1)
 * \return angle, in range 0 <= ang < 2PI */
double angle_around_axis(const Vec3d &v0, const Vec3d &v1, const Vec3d &axis);

} // namespace anti

#endif // VEC_UTILS_H
