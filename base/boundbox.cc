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

/*\file bbox.h
   \brief Bounding box and sphere
*/

#include "boundbox.h"

#include <vector>

using std::vector;

namespace anti {

const double big_double = 1e100;

BoundBox::BoundBox(const vector<Vec3d> &points, double cutoff)
    : min_coords(Vec3d(big_double, big_double, big_double)),
      max_coords(Vec3d(-big_double, -big_double, -big_double))
{
  if (points.size())
    add_points(points, cutoff);
}

void BoundBox::add_points(const vector<Vec3d> &points, double cutoff)
{
  double cut_off2 = cutoff * cutoff;
  for (const auto &point : points) {
    for (int j = 0; j < 3; j++) {
      if (cutoff < 0 || point.len2() < cut_off2) {
        if (point[j] < min_coords[j])
          min_coords[j] = point[j];
        if (point[j] > max_coords[j])
          max_coords[j] = point[j];
      }
    }
  }
}

void BoundBox::add_b_box(const BoundBox &b_box)
{
  for (int j = 0; j < 3; j++) {
    if (b_box.min_coords[j] < min_coords[j])
      min_coords[j] = b_box.min_coords[j];
    if (b_box.max_coords[j] > max_coords[j])
      max_coords[j] = b_box.max_coords[j];
  }
}

BoundSphere::BoundSphere(const vector<Vec3d> &points, double cutoff)
    : radius(-1), cut_off(cutoff), cut_off_cnt(0)
{
  if (points.size())
    add_points(points, cutoff);
}

void BoundSphere::add_points(const vector<Vec3d> &points, double cutoff)
{
  set_cut_off(cutoff);
  find_radius_centre(points);
}

void BoundSphere::add_b_sphere(const BoundSphere &b_sphere)
{
  // Set larger cut_off
  if (cut_off < 0 || b_sphere.cut_off < 0)
    set_cut_off(-1);
  else if (b_sphere.cut_off > cut_off)
    set_cut_off(b_sphere.cut_off);

  cut_off_cnt += b_sphere.cut_off_cnt;

  add_b_sphere(b_sphere.centre, b_sphere.radius);
}

void BoundSphere::find_radius_centre(const vector<Vec3d> &pts)
{
  double cut_off2 = cut_off * cut_off;
  double rad2 = 0;
  BoundBox bbox(pts, cut_off);
  Vec3d cent = bbox.get_centre();
  Vec3d centrd(0, 0, 0); // centroid
  int n = 0;             // number of points
  for (const auto &pt : pts) {
    double dist2 = (pt - cent).len2();
    if (cut_off < 0 || dist2 < cut_off2) {
      centrd += pt;
      n++;
      if (dist2 > rad2)
        rad2 = dist2;
    }
    else
      cut_off_cnt++;
  }

  centrd /= (double)n;
  double rad3 = 0;
  for (const auto &pt : pts) {
    double dist2 = (pt - centrd).len2();
    if (cut_off < 0 || dist2 < cut_off2) {
      if (dist2 > rad3)
        rad3 = dist2;
    }
  }

  double new_radius;
  Vec3d new_centre;
  if (rad3 < rad2 + epsilon) {
    new_radius = sqrt(rad3);
    new_centre = centrd;
  }
  else {
    new_radius = sqrt(rad2);
    new_centre = cent;
  }

  add_b_sphere(new_centre, new_radius);
}

void BoundSphere::add_b_sphere(Vec3d cent, double rad)
{
  if (centre.is_set()) {
    Vec3d offset = cent - centre;
    if (offset.len() <
        (radius + rad) * epsilon) // relatively very close centres
      radius = (radius > rad) ? radius : rad;
    else { // use container of two spheres
      Vec3d u = offset.unit();
      centre = 0.5 * ((centre - u * radius) + (cent + u * rad));
      radius = ((cent + u * rad) - centre).len();
    }
  }
  else {
    centre = cent;
    radius = rad;
  }
}

void BoundSphere::set_cut_off(double cutoff)
{
  if (cutoff > 0 && (cut_off < 0 || cutoff < cut_off))
    cut_off = cutoff;
}

} // namespace anti
