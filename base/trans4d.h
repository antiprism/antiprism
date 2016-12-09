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

/*!\file trans4d.h
  \brief Matrix transformations for 4D geometry.
*/

#ifndef TRANS4D_H
#define TRANS4D_H

#include <math.h>

#include "mathutils.h"
#include "vec4d.h"

namespace anti {

/// Matrix for transformations in 4D
/** A 4x4 matrix generally multiplying left to right on a column vector. */
class Trans4d {
private:
  double m[25];

public:
  /// Constructor
  /**Initialises to the identity */
  Trans4d() { to_unit(); }

  /// Constructor
  /**Initialise the rows of the rotation part.
   * \param r1 first row.
   * \param r2 second row.
   * \param r3 third row.
   * \param r4 fourth row. */
  Trans4d(Vec4d r1, Vec4d r2, Vec4d r3, Vec4d r4);

  /// Read access to the element array.
  /**\param idx index into the 16 element array.
   * \return The element. */
  double operator[](int idx) const { return m[idx]; }

  /// Write access to the element array.
  /**\param idx index into the 16 element array.
   * \return The element. */
  double &operator[](int idx) { return m[idx]; }

  /// Compose this transformation with another transformation
  /**Multiply this matrix by another matrix.
   * \param trans the transformation to compose with.
   * \return A reference to the result, held in this transformation. */
  Trans4d &operator*=(const Trans4d &trans);

  /// Multiply the matrix of this transformation by a scalar.
  /**\param n the scalar.
   * \return A reference to the result, held in this transformation. */
  Trans4d &operator*=(double n);

  /// Add another transformation matrix to this transformation matrix.
  /**\param trans the transformation whose matrix will be added.
   * \return A reference to the result, held in this transformation. */
  Trans4d &operator+=(const Trans4d &trans);

  /// Transpose the matrix
  /**\return A reference this transformation with its matrix transposed. */
  Trans4d &transpose();

  /// Set the transformation to the identity.
  /**\return A reference to this transformation. */
  Trans4d &to_unit();

  /// Get an identity transformation.
  /**\return An identity transformation. */
  static Trans4d unit();

  /// Set the transformation matrix to zero
  /**\return A reference to this transformation. */
  Trans4d &to_zero();

  /// Get a transformation with zero matrix.
  /**\return A transformation with zero matrix. */
  static Trans4d zero();

  /// Set a rotation by plane axis and angle.
  /**Rotate around an arbitrary plane axis through the origin.
   * \param u the first basis vector of the axis plane to rotate around.
   * \param v the second basis vector of the axis plane to rotate around.
   * \param angle the angle to rotate, in radians.
   * \return  A reference to this matrix with the resulting rotation set. */
  Trans4d &set_rot(Vec4d u, Vec4d v, double angle);

  /// Get a rotation by plane axis and angle.
  /**Rotate around an arbitrary plane axis through the origin.
   * \param u the first basis vector of the axis plane to rotate around.
   * \param v the second basis vector of the axis plane to rotate around.
   * \param angle the angle to rotate, in radians.
   * \return  A matrix with the resulting rotation set. */
  static Trans4d rot(Vec4d u, Vec4d v, double angle);

  /// Set a rotation by rotating one direction vector onto another.
  /**\param v_from the direction vector to rotate.
   * \param v_to the direction \a v_from should point after the rotation.
   * \return  A reference to this matrix with the resulting rotation set. */
  Trans4d &set_rot(Vec4d v_from, Vec4d v_to);

  /// Get a rotation by rotating one direction vector onto another.
  /**\param v_from the direction vector to rotate.
   * \param v_to the direction \a v_from should point after the rotation.
   * \return  A matrix with the resulting rotation set. */
  static Trans4d rot(Vec4d v_from, Vec4d v_to);

  /// Set a rotation by rotating around the coordinate-plane axes.
  /**Rotate around the x-axis, y-axis and then z-axis.
   * \param xy_ang angle in radians to rotate around xy-axis.
   * \param yz_ang angle in radians to rotate around yz-axis.
   * \param zw_ang angle in radians to rotate around zw-axis.
   * \param wx_ang angle in radians to rotate around wx-axis.
   * \param xz_ang angle in radians to rotate around xz-axis.
   * \param yw_ang angle in radians to rotate around yw-axis.
   * \return  A reference to this matrix with the resulting rotation set. */
  Trans4d &set_rot(double xy_ang, double yz_ang, double zw_ang, double wx_ang,
                   double xz_ang, double yw_ang);

  /// Get a rotation by rotating around the coordinate-plane axes.
  /**Rotate around the x-axis, y-axis and then z-axis.
   * \param xy_ang angle in radians to rotate around xy-axis.
   * \param yz_ang angle in radians to rotate around yz-axis.
   * \param zw_ang angle in radians to rotate around zw-axis.
   * \param wx_ang angle in radians to rotate around wx-axis.
   * \param xz_ang angle in radians to rotate around xz-axis.
   * \param yw_ang angle in radians to rotate around yw-axis.
   * \return  A matrix with the resulting rotation set. */
  static Trans4d rot(double xy_ang, double yz_ang, double zw_ang, double wx_ang,
                     double xz_ang, double yw_ang);

  /// Set a translation.
  /**\param transl the vector to translate by.
   * \return  A reference to this matrix with the translation set. */
  Trans4d &set_transl(Vec4d transl);

  /// Get a translation.
  /**\param transl the vector to translate by.
   * \return  A matrix with the translation set. */
  static Trans4d transl(Vec4d transl);

  /// Set a uniform scaling.
  /** Scale uniformly in all directions, away from the origin.
   * \param scal the amount to scale.
   * \return  A reference to this matrix with the scaling set. */
  Trans4d &set_scale(double scal);

  /// Get a uniform scaling.
  /**Scale uniformly in all directions, away from the origin.
   * \param scal the amount to scale.
   * \return  A matrix with the scaling set. */
  static Trans4d scale(double scal);

  /// Set a scaling by components.
  /**Scale in the direction of the coordinate axes, away from the origin.
   * \param x_scal the amount to scale in the direction of the x-axis.
   * \param y_scal the amount to scale in the direction of the y-axis.
   * \param z_scal the amount to scale in the direction of the z-axis.
   * \param w_scal the amount to scale in the direction of the z-axis.
   * \return  A reference to this matrix with the scaling set. */
  Trans4d &set_scale(double x_scal, double y_scal, double z_scal,
                     double w_scal);

  /// Get a scaling by components.
  /**Scale in the direction of the coordinate axes, away from the origin.
   * \param x_scal the amount to scale in the direction of the x-axis.
   * \param y_scal the amount to scale in the direction of the y-axis.
   * \param z_scal the amount to scale in the direction of the z-axis.
   * \param w_scal the amount to scale in the direction of the z-axis.
   * \return  A reference to this matrix with the scaling set. */
  static Trans4d scale(double x_scal, double y_scal, double z_scal,
                       double w_scal);

  /// Set an inversion.
  /** Inversion through the origin.
   * \return  A reference to this transformation set to an inversion. */
  Trans4d &set_inversion();

  /// Get an inversion.
  /** Inversion through the origin.
   * \return  A transformation set to an inversion. */
  static Trans4d inversion();

  /// Debugging print of the matrix of a transformation.
  /**\param var a string to identify the matrix variable.
   * \param file file stream to print the variable. */
  void dump(const char *var = "", FILE *file = stderr); // print matrix
};

/// Transform a column vector.
/**\param trans the transformation.
 * \param vec the column vector.
 * \return The transformed vector (left-multiplied by the
 *  transformation matrix). */
Vec4d operator*(const Trans4d &trans, const Vec4d &vec);

/// Transform a row vector.
/**\param vec the column vector.
 * \param trans the transformation matrix.
 * \return The transformed vector (right-multiplied by the
 *  transformation matrix). */
Vec4d operator*(const Vec4d &vec, const Trans4d &trans);

/// Compose two transformations
/**\param trans1 the first transformation.
 * \param trans2 the second transformation.
 * \return The composition of the two transformations (matrix of the first
 *  multiplying the matrix of the second). */
Trans4d operator*(const Trans4d &trans1, const Trans4d &trans2);

/// Multiply the matrix of this transformation by a scalar.
/**\param n the scalar.
 * \param trans the transformation whose matris is to be multiplied
 * \return The resulting transformation. */
Trans4d operator*(double n, const Trans4d &trans);

/// Multiply the matrix of this transformation by a scalar.
/**\param trans the transformation whose matris is to be multiplied
 * \param n the scalar.
 * \return The resulting transformation. */
Trans4d operator*(const Trans4d &trans, double n);

/// Add the matrices of two transformations.
/**\param trans1 the first transformation.
 * \param trans2 the second transformation.
 * \return The resulting transformation. */
Trans4d operator+(const Trans4d &trans1, const Trans4d &trans2);

// inline functions
inline Trans4d::Trans4d(Vec4d r1, Vec4d r2, Vec4d r3, Vec4d r4)
{
  to_unit();
  Vec4d *pvec[] = {&r1, &r2, &r3, &r4};
  for (int i = 0; i < 20; i++)
    if ((i % 5) != 4)
      m[i] = (*pvec[i / 5])[i % 5];
}

inline Trans4d &Trans4d::to_unit()
{
  for (int i = 0; i < 25; i++)
    m[i] = ((i % 5) == (i / 5)) ? 1 : 0;
  return *this;
}

inline Trans4d Trans4d::unit()
{
  Trans4d trans;
  return trans.to_unit();
}

inline Trans4d &Trans4d::to_zero()
{
  for (auto &entry : m)
    entry = 0;
  return *this;
}

inline Trans4d Trans4d::zero()
{
  Trans4d trans;
  return trans.to_zero();
}

inline Trans4d &Trans4d::set_rot(Vec4d u, Vec4d v, double angle)
{
  u.to_unit();
  v.to_unit();
  Trans4d A(Vec4d(0, +(u[2] * v[3] - u[3] * v[2]), -(u[1] * v[3] - u[3] * v[1]),
                  +(u[1] * v[2] - u[2] * v[1])),
            Vec4d(-(u[2] * v[3] - u[3] * v[2]), 0, +(u[0] * v[3] - u[3] * v[0]),
                  -(u[0] * v[2] - u[2] * v[0])),
            Vec4d(+(u[1] * v[3] - u[3] * v[1]), -(u[0] * v[3] - u[3] * v[0]), 0,
                  +(u[0] * v[1] - u[1] * v[0])),
            Vec4d(-(u[1] * v[2] - u[2] * v[1]), +(u[0] * v[2] - u[2] * v[0]),
                  -(u[0] * v[1] - u[1] * v[0]), 0));
  *this = (unit() + (sin(angle) * A) + ((1 - cos(angle)) * (A * A)));
  m[24] = 1;
  return *this;
}

inline Trans4d Trans4d::rot(Vec4d u, Vec4d v, double angle)
{
  Trans4d trans;
  return trans.set_rot(u, v, angle);
}

inline Trans4d &Trans4d::set_rot(Vec4d v_from, Vec4d v_to)
{
  v_from.unit();
  v_to.unit();
  Vec4d n2, n3, orth;
  Random rnd;
  rnd.time_seed();
  n2 = vcross(Vec4d::random(rnd).unit(), v_from, v_to);
  n3 = vcross(n2, v_from, v_to);
  orth = vcross(v_from, n2, n3);
  set_rot(n2, n3, acos(safe_for_trig(vdot(v_from, v_to))));
  return *this;
}

inline Trans4d Trans4d::rot(Vec4d v_from, Vec4d v_to)
{
  Trans4d trans;
  return trans.set_rot(v_from, v_to);
}

inline Trans4d &Trans4d::set_rot(double xy_ang, double yz_ang, double zw_ang,
                                 double wx_ang, double xz_ang, double yw_ang)
{
  int planes[] = {2, 3, 3, 1, 0, 1, 1, 2, 1, 3, 2, 0};
  double angs[] = {xy_ang, yz_ang, zw_ang, wx_ang, xz_ang, yw_ang};
  to_unit();
  for (int p = 0; p < 6; p++) {
    Vec4d p1(0, 0, 0, 0);
    p1[planes[2 * p]] = 1;
    Vec4d p2(0, 0, 0, 0);
    p2[planes[2 * p + 1]] = 1;
    *this = rot(p1, p2, angs[p]) * (*this);
  }
  return *this;
}

inline Trans4d Trans4d::rot(double xy_ang, double yz_ang, double zw_ang,
                            double wx_ang, double xz_ang, double yw_ang)
{
  Trans4d trans;
  return trans.set_rot(xy_ang, yz_ang, zw_ang, wx_ang, xz_ang, yw_ang);
}

inline Trans4d &Trans4d::set_transl(Vec4d transl)
{
  to_unit();
  m[4] = transl[0];
  m[9] = transl[1];
  m[14] = transl[2];
  m[19] = transl[3];
  return *this;
}

inline Trans4d Trans4d::transl(Vec4d transl)
{
  Trans4d trans;
  return trans.set_transl(transl);
}

inline Trans4d &Trans4d::set_scale(double scal)
{
  return set_scale(scal, scal, scal, scal);
}

inline Trans4d Trans4d::scale(double scal)
{
  Trans4d trans;
  return trans.set_scale(scal);
}

inline Trans4d &Trans4d::set_scale(double x_scal, double y_scal, double z_scal,
                                   double w_scal)
{
  m[0] = x_scal;
  m[6] = y_scal;
  m[12] = z_scal;
  m[18] = w_scal;
  return *this;
}

inline Trans4d Trans4d::scale(double x_scal, double y_scal, double z_scal,
                              double w_scal)
{
  Trans4d trans;
  return trans.set_scale(x_scal, y_scal, z_scal, w_scal);
}

inline Trans4d &Trans4d::transpose()
{
  for (int i = 0; i < 25; i++)
    if (i / 5 < i % 5) {
      double tmp = m[i];
      m[i] = m[(i % 5) * 5 + (i / 5)];
      m[(i % 5) * 5 + (i / 5)] = tmp;
    }

  return *this;
}

inline Trans4d &Trans4d::set_inversion() { return set_scale(-1); }

inline Trans4d Trans4d::inversion()
{
  Trans4d trans;
  return trans.set_inversion();
}

inline Trans4d &Trans4d::operator*=(const Trans4d &trans)
{
  Trans4d new_m;
  new_m.to_zero();
  for (int i = 0; i < 25; i++)
    for (int j = 0; j < 5; j++)
      new_m[i] += m[(i / 5) * 5 + j] * trans[(j)*5 + (i % 5)];

  *this = new_m;
  return *this;
}

inline Trans4d &Trans4d::operator*=(double n)
{
  for (auto &entry : m)
    entry *= n;
  return *this;
}

inline Trans4d &Trans4d::operator+=(const Trans4d &trans)
{
  for (int i = 0; i < 25; i++)
    m[i] += trans[i];
  return *this;
}

inline Vec4d operator*(const Trans4d &trans, const Vec4d &vec)
{
  Vec4d new_v(0, 0, 0, 0);
  for (int i = 0; i < 20; i++) {
    if (i % 5 != 4)
      new_v[i / 5] += trans[i] * vec[i % 5];
    else
      new_v[i / 5] += trans[i];
  }
  return new_v;
}

inline Vec4d operator*(const Vec4d &vec, const Trans4d &trans)
{
  Trans4d m_ret = trans;
  m_ret.transpose();
  return m_ret * vec;
}

inline Trans4d operator*(const Trans4d &trans1, const Trans4d &trans2)
{
  Trans4d m_ret = trans1;
  return m_ret *= trans2;
}

inline Trans4d operator*(double n, const Trans4d &trans)
{
  Trans4d m_ret = trans;
  return m_ret *= n;
}

inline Trans4d operator*(const Trans4d &trans, double n)
{
  Trans4d m_ret = trans;
  return m_ret *= n;
}

inline Trans4d operator+(const Trans4d &trans1, const Trans4d &trans2)
{
  Trans4d m_ret = trans1;
  return m_ret += trans2;
}

} // namespace anti

#endif // TRANS4D_H
