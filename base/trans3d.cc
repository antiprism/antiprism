/*
   Copyright (c) 2003-2021, Adrian Rossiter

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

/* \file Trans3d.cc
 *\brief Matrix transformations for 3D geometry.
 */

#include "trans3d.h"
#include "mathutils.h"

#include <cstdio>
#include <vector>

using std::vector;

namespace anti {

/*
bool less_than(const Trans3d &m1, const Trans3d &m2, double eps)
{
   for(int i=0; i<12; i++) {
      if(fabs(m1[i] - m2[i]) > eps) {
         if(m1[i] < m2[i])
            return 1;
         else
            return 0;
      }
   }
   return 0;
}
*/

int compare(const Trans3d &trans1, const Trans3d &trans2, double eps)
{
  for (int i = 0; i < 12; i++) {
    if (fabs(trans1[i] - trans2[i]) > eps) {
      if (trans1[i] < trans2[i])
        return 1;
      else
        return -1;
    }
  }
  return 0;
}

double Trans3d::det() const { return determinant(m, 4); }

Trans3d Trans3d::rotate(Vec3d axis, double angle)
{
  Trans3d trans; // unit matrix.
  axis.to_unit();
  double c = cos(angle);
  double s = sin(angle);
  double t = 1.0 - c;
  trans[0] = c + axis[0] * axis[0] * t;
  trans[5] = c + axis[1] * axis[1] * t;
  trans[10] = c + axis[2] * axis[2] * t;

  double tmp1 = axis[0] * axis[1] * t;
  double tmp2 = axis[2] * s;
  trans[4] = tmp1 + tmp2;
  trans[1] = tmp1 - tmp2;

  tmp1 = axis[0] * axis[2] * t;
  tmp2 = axis[1] * s;
  trans[8] = tmp1 - tmp2;
  trans[2] = tmp1 + tmp2;

  tmp1 = axis[1] * axis[2] * t;
  tmp2 = axis[0] * s;
  trans[9] = tmp1 + tmp2;
  trans[6] = tmp1 - tmp2;
  return trans;
}

Trans3d Trans3d::rotate(Vec3d v_from, Vec3d v_to)
{
  v_from.to_unit();
  v_to.to_unit();
  Vec3d axis = vcross(v_from, v_to);
  double cos_a = vdot(v_from, v_to);
  if (fabs(cos_a) >= 1 - anti::epsilon) {
    cos_a = cos_a > 0 ? 1 : -1;
    axis = vcross(v_from, Vec3d(1.2135, 2.09865, 3.23784)); // fix this
  }

  return rotate(axis, acos(cos_a));
}

Trans3d Trans3d::reflection(Vec3d norm)
{
  Trans3d trans;
  norm.to_unit();
  const double x = norm[0];
  const double y = norm[1];
  const double z = norm[2];
  trans[0] = -x * x + y * y + z * z;
  trans[1] = -2 * x * y;
  trans[2] = -2 * x * z;
  trans[4] = -2 * y * x;
  trans[5] = x * x - y * y + z * z;
  trans[6] = -2 * y * z;
  trans[8] = -2 * z * x;
  trans[9] = -2 * z * y;
  trans[10] = x * x + y * y - z * z;
  return trans;
}

Trans3d Trans3d::angles_between_axes(double yz_ang, double zx_ang,
                                     double xy_ang, bool *valid)
{
  double sq = 1 - cos(xy_ang) * cos(xy_ang) - cos(yz_ang) * cos(yz_ang) -
              cos(zx_ang) * cos(zx_ang) +
              2 * cos(xy_ang) * cos(yz_ang) * cos(zx_ang);
  if (sq < -anti::epsilon) {
    if (valid)
      *valid = false;
    return zero();
  }
  else if (sq < 0)
    sq = 0;

  if (valid)
    *valid = true;
  Vec3d new_x(1, 0, 0);
  Vec3d new_y(cos(xy_ang), sin(xy_ang), 0);
  Vec3d new_z(cos(zx_ang),
              -cos(zx_ang) / tan(xy_ang) + cos(yz_ang) / sin(xy_ang),
              sqrt(sq) / sin(xy_ang));
  return Trans3d(new_x, new_y, new_z).transpose();
}

Trans3d Trans3d::inverse() const
{
  Trans3d inv;
  double determ = determinant(m, 4);
  if (fabs(determ) > anti::epsilon) {
    // http://www.euclideanspace.com/maths/algebra/matrix/functions
    // /inverse/fourD/index.htm
    inv[0] = m[6] * m[11] * m[13] - m[7] * m[10] * m[13] + m[7] * m[9] * m[14] -
             m[5] * m[11] * m[14] - m[6] * m[9] * m[15] + m[5] * m[10] * m[15];
    inv[1] = m[3] * m[10] * m[13] - m[2] * m[11] * m[13] - m[3] * m[9] * m[14] +
             m[1] * m[11] * m[14] + m[2] * m[9] * m[15] - m[1] * m[10] * m[15];
    inv[2] = m[2] * m[7] * m[13] - m[3] * m[6] * m[13] + m[3] * m[5] * m[14] -
             m[1] * m[7] * m[14] - m[2] * m[5] * m[15] + m[1] * m[6] * m[15];
    inv[3] = m[3] * m[6] * m[9] - m[2] * m[7] * m[9] - m[3] * m[5] * m[10] +
             m[1] * m[7] * m[10] + m[2] * m[5] * m[11] - m[1] * m[6] * m[11];
    inv[4] = m[7] * m[10] * m[12] - m[6] * m[11] * m[12] - m[7] * m[8] * m[14] +
             m[4] * m[11] * m[14] + m[6] * m[8] * m[15] - m[4] * m[10] * m[15];
    inv[5] = m[2] * m[11] * m[12] - m[3] * m[10] * m[12] + m[3] * m[8] * m[14] -
             m[0] * m[11] * m[14] - m[2] * m[8] * m[15] + m[0] * m[10] * m[15];
    inv[6] = m[3] * m[6] * m[12] - m[2] * m[7] * m[12] - m[3] * m[4] * m[14] +
             m[0] * m[7] * m[14] + m[2] * m[4] * m[15] - m[0] * m[6] * m[15];
    inv[7] = m[2] * m[7] * m[8] - m[3] * m[6] * m[8] + m[3] * m[4] * m[10] -
             m[0] * m[7] * m[10] - m[2] * m[4] * m[11] + m[0] * m[6] * m[11];
    inv[8] = m[5] * m[11] * m[12] - m[7] * m[9] * m[12] + m[7] * m[8] * m[13] -
             m[4] * m[11] * m[13] - m[5] * m[8] * m[15] + m[4] * m[9] * m[15];
    inv[9] = m[3] * m[9] * m[12] - m[1] * m[11] * m[12] - m[3] * m[8] * m[13] +
             m[0] * m[11] * m[13] + m[1] * m[8] * m[15] - m[0] * m[9] * m[15];
    inv[10] = m[1] * m[7] * m[12] - m[3] * m[5] * m[12] + m[3] * m[4] * m[13] -
              m[0] * m[7] * m[13] - m[1] * m[4] * m[15] + m[0] * m[5] * m[15];
    inv[11] = m[3] * m[5] * m[8] - m[1] * m[7] * m[8] - m[3] * m[4] * m[9] +
              m[0] * m[7] * m[9] + m[1] * m[4] * m[11] - m[0] * m[5] * m[11];
    inv[12] = m[6] * m[9] * m[12] - m[5] * m[10] * m[12] - m[6] * m[8] * m[13] +
              m[4] * m[10] * m[13] + m[5] * m[8] * m[14] - m[4] * m[9] * m[14];
    inv[13] = m[1] * m[10] * m[12] - m[2] * m[9] * m[12] + m[2] * m[8] * m[13] -
              m[0] * m[10] * m[13] - m[1] * m[8] * m[14] + m[0] * m[9] * m[14];
    inv[14] = m[2] * m[5] * m[12] - m[1] * m[6] * m[12] - m[2] * m[4] * m[13] +
              m[0] * m[6] * m[13] + m[1] * m[4] * m[14] - m[0] * m[5] * m[14];
    inv[15] = m[1] * m[6] * m[8] - m[2] * m[5] * m[8] + m[2] * m[4] * m[9] -
              m[0] * m[6] * m[9] - m[1] * m[4] * m[10] + m[0] * m[5] * m[10];

    for (double &val : inv.m)
      val /= determ;
  }
  else // det very small
    inv = zero();

  return inv;
}

Trans3d &Trans3d::operator*=(const Trans3d &trans)
{
  Trans3d new_m = zero();
  for (int i = 0; i < 16; i++)
    for (int j = 0; j < 4; j++)
      new_m[i] += m[(i / 4) * 4 + j] * trans[(j)*4 + (i % 4)];

  *this = new_m;
  return *this;
}

Vec3d operator*(const Trans3d &trans, const Vec3d &vec)
{
  auto new_v = Vec3d::zero;
  for (int i = 0; i < 12; i++) {
    if (i % 4 != 3)
      new_v[i / 4] += trans[i] * vec[i % 4];
    else
      new_v[i / 4] += trans[i];
  }
  return new_v;
}

Vec4d operator*(const Trans3d &trans, const Vec4d &vec)
{
  auto new_v = Vec4d::zero;
  for (int i = 0; i < 16; i++)
    new_v[i / 4] += trans[i] * vec[i % 4];
  return new_v;
}

// http://www.j3d.org/matrix_faq/matrfaq_latest.html  Q55.
Vec4d Trans3d::get_quaternion() const
{
  Vec4d quat; // X, Y, Z, W
  double T = 1 + m[0] + m[5] + m[10];
  double S;

  // If the trace of the matrix is equal to zero then identify
  // which major diagonal element has the greatest value.
  // Depending on this, calculate the following:
  if (T > anti::epsilon) {
    S = sqrt(T) * 2;
    quat[0] = (m[9] - m[6]) / S;
    quat[1] = (m[2] - m[8]) / S;
    quat[2] = (m[4] - m[1]) / S;
    quat[3] = 0.25 * S;
  }
  else if (m[0] > m[5] && m[0] > m[10]) { // Column 0:
    S = sqrt(1.0 + m[0] - m[5] - m[10]) * 2;
    quat[0] = 0.25 * S;
    quat[1] = (m[4] + m[1]) / S;
    quat[2] = (m[2] + m[8]) / S;
    quat[3] = (m[9] - m[6]) / S;
  }
  else if (m[5] > m[10]) { // Column 1:
    S = sqrt(1.0 + m[5] - m[0] - m[10]) * 2;
    quat[0] = (m[4] + m[1]) / S;
    quat[1] = 0.25 * S;
    quat[2] = (m[9] + m[6]) / S;
    quat[3] = (m[2] - m[8]) / S;
  }
  else { // Column 2:
    S = sqrt(1.0 + m[10] - m[0] - m[5]) * 2;
    quat[0] = (m[2] + m[8]) / S;
    quat[1] = (m[9] + m[6]) / S;
    quat[2] = 0.25 * S;
    quat[3] = (m[4] - m[1]) / S;
  }
  return quat;
}

// http://www.gamedev.net/community/forums/topic.asp?topic_id=166424
Vec3d quat2euler(const Vec4d &quat)
{
  double sqx = quat[0] * quat[0];
  double sqy = quat[1] * quat[1];
  double sqz = quat[2] * quat[2];
  double sqw = quat[3] * quat[3];

  return -Vec3d(
      atan2(2 * (quat[1] * quat[2] + quat[0] * quat[3]),
            -sqx - sqy + sqz + sqw),
      asin(safe_for_trig(-2 * (quat[2] * quat[0] - quat[1] * quat[3]))),
      atan2(2 * (quat[0] * quat[1] + quat[2] * quat[3]),
            sqx - sqy - sqz + sqw));
}

Vec3d Trans3d::get_euler() const { return quat2euler(get_quaternion()); }

void Trans3d::dump(const char *var, FILE *file) const
{
  fprintf(file, "%s", var);
  for (int i = 0; i < 16; i++) {
    if (!(i % 4))
      fprintf(file, "\n\t");
    fprintf(file, "\t%g", m[i]);
  }
  fprintf(file, "\n");
}

Trans3d Trans3d::align(vector<Vec3d> from, vector<Vec3d> to)
{
  const Trans3d r = Trans3d::rotate(0.1, 0, 0) * Trans3d::rotate(0, 0.2, 0) *
                    Trans3d::rotate(0, 0, 0.3);
  const Trans3d inv_r = Trans3d::rotate(0, 0, -0.3) *
                        Trans3d::rotate(0, -0.2, 0) *
                        Trans3d::rotate(-0.1, 0, 0);

  if (to.size() > 1)
    transform(to, r);

  Trans3d trans = Trans3d::translate(to[0] - from[0]);

  if (to.size() > 1) {
    trans = Trans3d::translate(to[0]) *
            Trans3d::rotate(from[1] - from[0], to[1] - to[0]) *
            Trans3d::translate(-to[0]) * trans;
  }

  if (to.size() > 2) {
    for (unsigned int j = 0; j < 3; j++)
      from[j] = trans * from[j];

    Vec3d norm1 = vcross(from[2] - from[0], from[1] - from[0]);
    Vec3d norm2 = vcross(to[2] - to[0], to[1] - to[0]);
    // Maybe test form norm size here
    norm1.to_unit();
    norm2.to_unit();
    trans = Trans3d::translate(to[0]) * Trans3d::rotate(norm1, norm2) *
            Trans3d::translate(-to[0]) * trans;
  }

  return (to.size() > 1) ? inv_r * trans : trans;
}

Trans3d Trans3d::align(Vec3d from1, Vec3d from2, Vec3d from3, Vec3d to1,
                       Vec3d to2, Vec3d to3)
{
  vector<Vec3d> from(3);
  from[0] = from1;
  from[1] = from2;
  from[2] = from3;
  vector<Vec3d> to(3);
  to[0] = to1;
  to[1] = to2;
  to[2] = to3;
  return align(from, to);
}

Trans3d Trans3d::align(Vec3d from1, Vec3d from2, Vec3d to1, Vec3d to2)
{

  const Trans3d r = Trans3d::rotate(0.1, 0, 0) * Trans3d::rotate(0, 0.2, 0) *
                    Trans3d::rotate(0, 0, 0.3);
  const Trans3d inv_r = Trans3d::rotate(0, 0, -0.3) *
                        Trans3d::rotate(0, -0.2, 0) *
                        Trans3d::rotate(-0.1, 0, 0);
  to1 = r * to1;
  to2 = r * to2;

  Trans3d trans = Trans3d::rotate(from1, to1);
  from1 = trans * from1;
  from2 = trans * from2;

  Vec3d norm1 = vcross(from2, from1);
  Vec3d norm2 = vcross(to2, to1);
  // Maybe test form norm size here
  norm1.to_unit();
  norm2.to_unit();

  // find the angle to rotate abot to1, in the range -180<ang<=180
  double ang = acos(safe_for_trig(vdot(norm1, norm2)));
  if (vtriple(to1, norm1, norm2) < 0)
    ang *= -1;

  // trans = Trans3d::rotate(to1, acos(vdot(norm1, norm2))) * trans;
  //*this = inv_r*trans;

  trans = Trans3d::rotate(to1, ang) * trans;
  return inv_r * trans;
}

} // namespace anti
