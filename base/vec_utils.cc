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
   Name: vec_utils.cc
   Description: vector utilities
   Project: Antiprism - http://www.antiprism.com
*/

#include <vector>

#include "boundbox.h"
#include "mathutils.h"
#include "vec_utils.h"

using std::vector;

namespace anti {

Vec3d line_plane_intersect(Vec3d Q, Vec3d n, Vec3d P0, Vec3d P1, int *where,
                           double eps)
{
  double seg_proj = vdot(n, P1 - P0);
  if (fabs(seg_proj) < eps)
    return Vec3d();
  double t = vdot(n, Q - P0) / seg_proj;
  if (where) {
    *where = 0; // X < P0
    if (t > 0)
      *where = 2; // P0 < X < P1
    if (t > 1)
      *where = 4; // P1 < X
    if (t > 0 - eps && t < 0 + eps)
      *where = 1; // X = P0
    if (t > 1 - eps && t < 1 + eps)
      *where = 3; // X = P1
  }
  return P0 + t * (P1 - P0);
}

Vec3d line_plane_intersect(Vec3d Q0, Vec3d Q1, Vec3d Q2, Vec3d P0, Vec3d P1,
                           int *where, double eps)
{
  Vec3d n = vcross(Q1 - Q0, Q2 - Q0);
  // n.unit();
  return line_plane_intersect(Q0, n, P0, P1, where, eps);
}

// http://astronomy.swin.edu.au/~pbourke/geometry/planeplane/
bool two_plane_intersect(Vec3d Q0, Vec3d n0, Vec3d Q1, Vec3d n1, Vec3d &P,
                         Vec3d &dir, double eps)
{
  dir = vcross(n0, n1) / (n0.len() * n1.len());
  if (fabs(dir.len()) < eps)
    return false;
  dir.to_unit();
  double d1 = vdot(Q0, n0);
  double d2 = vdot(Q1, n1);
  double det = vdot(n0, n0) * vdot(n1, n1) - pow(vdot(n0, n1), 2);
  double c1 = (d1 * vdot(n1, n1) - d2 * vdot(n0, n1)) / det;
  double c2 = (d2 * vdot(n0, n0) - d1 * vdot(n0, n1)) / det;
  P = c1 * n0 + c2 * n1;
  return true;
}

// http://astronomy.swin.edu.au/~pbourke/geometry/planes/
bool three_plane_intersect(Vec3d Q0, Vec3d n0, Vec3d Q1, Vec3d n1, Vec3d Q2,
                           Vec3d n2, Vec3d &P, double eps)
{
  double tri_prod = vtriple(n0, n1, n2);
  if (fabs(tri_prod) / (n0.len2() * n1.len2() * n2.len2()) < eps * eps)
    return false;
  double d1 = vdot(Q0, n0);
  double d2 = vdot(Q1, n1);
  double d3 = vdot(Q2, n2);
  P = (d1 * vcross(n1, n2) + d2 * vcross(n2, n0) + d3 * vcross(n0, n1)) /
      tri_prod;
  return true;
}

// http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm
bool lines_nearest_points(Vec3d P0, Vec3d P1, Vec3d Q0, Vec3d Q1, Vec3d &P,
                          Vec3d &Q, double eps)
{
  Vec3d u = (P1 - P0);
  Vec3d v = (Q1 - Q0);
  Vec3d w = P0 - Q0;
  double a = vdot(u, u);
  double b = vdot(u, v);
  double c = vdot(v, v);
  double d = vdot(u, w);
  double e = vdot(v, w);
  double D = a * c - b * b;
  double s = (b * e - c * d) / D;
  double t = (a * e - b * d) / D;
  P = P0 + s * u;
  Q = Q0 + t * v;
  return (D >= eps);
}

Vec3d nearest_point(Vec3d P, Vec3d Q0, Vec3d Q1, Vec3d Q2, double eps)
{
  Vec3d norm = vcross(Q0 - Q1, Q2 - Q1);
  if (norm.len2() <= (eps * eps))
    return Vec3d();

  norm.to_unit();
  return P + vdot(Q0 - P, norm) * norm;
}

Vec3d nearest_point(Vec3d P, const vector<Vec3d> &points, double eps)
{
  switch (points.size()) {
  case (0):
    return P;
  case (1):
    return points[0];
  case (2):
    if ((points[1] - points[0]).len2() <= eps * eps)
      return points[0];
    else
      return nearest_point(P, points[0], points[1]);
  default:
    Vec3d nearest;
    unsigned int p2;
    for (p2 = 1; p2 < points.size(); p2++)
      if ((points[p2] - points[0]).len2() > eps * eps)
        break;
    if (p2 == points.size())
      return points[0];
    if (p2 == points.size() - 1)
      return nearest_point(P, points[0], points[1]);

    for (unsigned int i = p2; i < points.size(); i++) {
      nearest = nearest_point(P, points[0], points[p2], points[i], eps);
      if (nearest.is_set())
        return nearest;
    }

    return nearest_point(P, points[0], points[1]);
  }
}

Vec3d nearest_point(Vec3d P, const vector<Vec3d> &points,
                    const vector<int> &idxs, double eps)
{
  vector<Vec3d> points2;
  for (int idx : idxs)
    points2.push_back(points[idx]);

  return nearest_point(P, points2, eps);
}

double angle_around_axis(const Vec3d &v0, const Vec3d &v1, const Vec3d &axis)
{
  Vec3d n0 = vcross(v0, axis).unit();
  Vec3d n1 = vcross(v1, axis).unit();
  double ang = acos(safe_for_trig(vdot(n0, n1)));
  if (vdot(axis, vcross(n0, n1)) < 0)
    ang = 2 * M_PI - ang;
  return ang;
}

Vec3d lines_intersection(const Vec3d &P0, const Vec3d &P1, const Vec3d &Q0,
                         const Vec3d &Q1, double eps)
{
  Vec3d P, Q;
  if (eps == 0.0) { // always return an intersection point
    if (!lines_nearest_points(P0, P1, Q0, Q1, P, Q, 0)) {
      P = (P0 + P1) / 2.0;
      Q = (Q0 + Q1) / 2.0;
    }
    return (P + Q) / 2.0;
  }

  // lines might not be parallel and still miss so check if nearest points is
  // not zero
  if (lines_nearest_points(P0, P1, Q0, Q1, P, Q, eps) && ((P - Q).len() < eps))
    return (P + Q) / 2.0;
  else
    return Vec3d();
}

bool in_segment(const Vec3d &P, const Vec3d &Q0, const Vec3d &Q1, double eps)
{
  BoundBox bb;
  vector<Vec3d> points(2);

  points[0] = Q0;
  points[1] = Q1;

  bb.add_points(points);
  Vec3d min = bb.get_min();
  Vec3d max = bb.get_max();

  return (compare(min, P, eps) <= 0 && compare(max, P, eps) >= 0);
}

Vec3d segments_intersection(const Vec3d &P0, const Vec3d &P1, const Vec3d &Q0,
                            const Vec3d &Q1, double eps)
{
  Vec3d intersection = lines_intersection(P0, P1, Q0, Q1, eps);

  if (intersection.is_set() && (!in_segment(intersection, P0, P1, eps) ||
                                !in_segment(intersection, Q0, Q1, eps)))
    intersection = Vec3d();

  return intersection;
}

Vec3d point_in_segment(const Vec3d &P, const Vec3d &Q0, const Vec3d &Q1,
                       double eps)
{
  Vec3d intersection = nearest_point(P, Q0, Q1);

  if (intersection.is_set() &&
      (compare(P, intersection, eps) || !in_segment(intersection, Q0, Q1, eps)))
    intersection = Vec3d();

  return intersection;
}

} // namespace anti
