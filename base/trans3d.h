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

/*!\file trans3d.h
 *\brief Matrix transformations for 3D geometry.
*/

#ifndef TRANS3D_H
#define TRANS3D_H

#include "vec3d.h"
#include "vec4d.h"
#include <math.h>
#include <vector>

namespace anti {

/// Transformations in 3D
/** A 4x4 matrix generally multiplying left to right on a column vector. */
class Trans3d {
private:
  static inline double det2(const double &a11, const double &a12,
                            const double &a21, const double &a22)
  {
    return a11 * a22 - a12 * a21;
  }

  double m[16];

public:
  /// Constructor
  /** Initialises to the identity. */
  Trans3d() { to_unit(); }

  /// Constructor
  /** Initialise the rows of the rotation part.
   * \param r1 first row.
   * \param r2 second row.
   * \param r3 third row. */
  Trans3d(const Vec3d &r1, const Vec3d &r2, const Vec3d &r3);

  /// Constructor
  /** Initialise the rows.
   * \param r1 first row.
   * \param r2 second row.
   * \param r3 third row.
   * \param r4 fourth row. */
  Trans3d(const Vec4d &r1, const Vec4d &r2, const Vec4d &r3,
          Vec4d r4 = Vec4d(0, 0, 0, 1));

  /// Set the columns.
  /** Initialise the top three rows by columns
   * \param c1 first column.
   * \param c2 second column.
   * \param c3 third column.
   * \param c4 fourth column.
   * \return A reference to the result, held in this matrix. */
  Trans3d &set_cols(const Vec3d &c1, const Vec3d &c2, const Vec3d &c3,
                    const Vec3d &c4 = Vec3d(0, 0, 0));

  /// Read access to the element array.
  /**\param idx index into the 16 element array.
   * \return The element. */
  double operator[](int idx) const { return m[idx]; }

  /// Write access to the element array.
  /**\param idx index into the 16 element array.
   * \return The element. */
  double &operator[](int idx) { return m[idx]; }

  /// Read access to the start of the element array.
  /**\return The start of the array. */
  const double *get_m() { return m; }

  /// Compose this transformation with another transformation
  /**Multiply this matrix by another matrix.
   * \param trans the transformation to compose with.
   * \return A reference to the result, held in this transformation. */
  Trans3d &operator*=(const Trans3d &trans);

  /// Transpose the matrix
  /**\return A reference this transformation with its matrix transposed. */
  Trans3d &transpose();

  /// Set to the inverse.
  /**\return  A reference to this transformation set to its inverse. */
  Trans3d &set_inverse();

  /// Get the inverse.
  /**\param trans the transformation to get the inverse of.
   * \return  The inverse of \p trans. */
  static Trans3d inverse(const Trans3d &trans);

  /// Determinant
  /**\return The determinant of the transformation matrix. */
  double det() const;

  /// Set the transformation to the identity.
  /**\return A reference to this transformation. */
  Trans3d &to_unit();

  /// Get an identity transformation.
  /**\return An identity transformation. */
  static Trans3d unit();

  /// Set the transformation matrix to zero
  /**\return A reference to this transformation. */
  Trans3d &to_zero();

  /// Get a transformation with zero matrix.
  /**\return A transformation with zero matrix. */
  static Trans3d zero();

  /// Set a rotation by rotating around the coordinate axes.
  /** Rotate around the x-axis, y-axis and then z-axis.
   * \param x_ang angle in radians to rotate around x-axis.
   * \param y_ang angle in radians to rotate around y-axis.
   * \param z_ang angle in radians to rotate around z-axis.
   * \return  A reference to this matrix with the resulting rotation set. */
  Trans3d &set_rot(double x_ang, double y_ang, double z_ang);

  /// Get a rotation by rotating around the coordinate axes.
  /** Rotate around the x-axis, y-axis and then z-axis.
   * \param x_ang angle in radians to rotate around x-axis.
   * \param y_ang angle in radians to rotate around y-axis.
   * \param z_ang angle in radians to rotate around z-axis.
   * \return  A matrix with the resulting rotation set. */
  static Trans3d rot(double x_ang, double y_ang, double z_ang);

  /// Set a rotation by rotating around the coordinate axes.
  /** Rotate around the x-axis, y-axis then z-axis.
   * \param angles a vector whose components are the the angles in radians
   *  to rotate around he x-axis, y-axis and z-axis.
   * \return  A reference to this matrix with the resulting rotation set. */
  Trans3d &set_rot(Vec3d angles);

  /// Get a rotation by rotating around the coordinate axes.
  /** Rotate around the x-axis, y-axis then z-axis.
   * \param angles a vector whose components are the the angles in radians
   *  to rotate around he x-axis, y-axis and z-axis.
   * \return  A matrix with the resulting rotation set. */
  static Trans3d rot(Vec3d angles);

  /// Set a rotation by axis and angle.
  /** Rotate around an arbitrary axis through the origin.
   * \param axis the axis o rotate around.
   * \param angle the angle to rotate, in radians.
   * \return  A reference to this matrix with the resulting rotation set. */
  Trans3d &set_rot(Vec3d axis, double angle);

  /// Get a rotation by axis and angle.
  /** Rotate around an arbitrary axis through the origin.
   * \param axis the axis o rotate around.
   * \param angle the angle to rotate, in radians.
   * \return  A matrix with the resulting rotation set. */
  static Trans3d rot(Vec3d axis, double angle);

  /// Set a rotation by rotating one direction vector onto another.
  /**\param v_from the direction vector to rotate.
   * \param v_to the direction \a v_from should point after the rotation.
   * \return  A reference to this matrix with the resulting rotation set. */
  Trans3d &set_rot(Vec3d v_from, Vec3d v_to);

  /// Get a rotation by rotating one direction vector onto another.
  /**\param v_from the direction vector to rotate.
   * \param v_to the direction \a v_from should point after the rotation.
   * \return  A reference to this matrix with the resulting rotation set. */
  static Trans3d rot(Vec3d v_from, Vec3d v_to);

  /// Set a translation.
  /**\param transl the vector to translate by.
   * \return  A reference to this matrix with the translation set. */
  Trans3d &set_transl(Vec3d transl);

  /// Get a translation.
  /**\param transl the vector to translate by.
   * \return  A matrix with the translation set. */
  static Trans3d transl(Vec3d transl);

  /// Set a reflection.
  /** Reflect in a mirror plane passing through the origin.
   * \param norm a normal to the mirror plane.
   * \return  A reference to this matrix with the reflection set. */
  Trans3d &set_refl(Vec3d norm);

  /// Get a reflection.
  /** Reflect in a mirror plane passing through the origin.
   * \param norm a normal to the mirror plane.
   * \return  A matrix with the reflection set. */
  static Trans3d refl(Vec3d norm);

  /// Set a uniform scaling.
  /** Scale uniformly in all directions, away from the origin.
   * \param scal the amount to scale.
   * \return  A reference to this matrix with the scaling set. */
  Trans3d &set_scale(double scal);

  /// Get a uniform scaling.
  /** Scale uniformly in all directions, away from the origin.
   * \param scal the amount to scale.
   * \return  A matrix with the scaling set. */
  static Trans3d scale(double scal);

  /// Set a scaling by components.
  /** Scale in the direction of the coordinate axes, away from the origin.
   * \param x_scal the amount to scale in the direction of the x-axis.
   * \param y_scal the amount to scale in the direction of the y-axis.
   * \param z_scal the amount to scale in the direction of the z-axis.
   * \return  A reference to this matrix with the scaling set. */
  Trans3d &set_scale(double x_scal, double y_scal, double z_scal);

  /// Get a scaling by components.
  /** Scale in the direction of the coordinate axes, away from the origin.
   * \param x_scal the amount to scale in the direction of the x-axis.
   * \param y_scal the amount to scale in the direction of the y-axis.
   * \param z_scal the amount to scale in the direction of the z-axis.
   * \return  A reference to this matrix with the scaling set. */
  static Trans3d scale(double x_scal, double y_scal, double z_scal);

  /// Set a scaling in a single direction.
  /** Scale in a single arbitrary direction, away from the origin.
   * \param dir the direction to scale.
   * \param scal the amount to scale.
   * \return  A reference to this matrix with the scaling set. */
  Trans3d &set_scale(Vec3d dir, double scal);

  /// Get a scaling in a single direction.
  /** Scale in a single arbitrary direction, away from the origin.
   * \param dir the direction to scale.
   * \param scal the amount to scale.
   * \return  A matrix with the scaling set. */
  static Trans3d scale(Vec3d dir, double scal);

  /// Set an inversion.
  /** Inversion through the origin.
   * \return  A reference to this transformation set to an inversion. */
  Trans3d &set_inversion();

  /// Get an inversion.
  /** Inversion through the origin.
   * \return  A transformation set to an inversion. */
  static Trans3d inversion();

  /// Set a transformation that aligns two pairs of vectors
  /**\param from1 first vector to align.
   * \param from2 second vector to align.
   * \param to1 align \p from1 with this vector.
   * \param to2 align \p from2 with this vector.
   * \return reference to this transformation set to the alignment. */
  Trans3d &set_alignment(Vec3d from1, Vec3d from2, Vec3d to1, Vec3d to2);

  /// Get a transformation that aligns two pairs of vectors
  /**\param from1 first vector to align.
   * \param from2 second vector to align.
   * \param to1 align \p from1 with this vector.
   * \param to2 align \p from2 with this vector.
   * \return the alignment transformation. */
  static Trans3d alignment(Vec3d from1, Vec3d from2, Vec3d to1, Vec3d to2);

  /// Set a transformation that aligns two sets of three points
  /**\param from1 first point to align.
   * \param from2 second point to align.
   * \param from3 third point to align.
   * \param to1 align \p from1 with this point.
   * \param to2 align \p from2 with this point.
   * \param to3 align \p from3 with this point.
   * \return reference to this transformation set to the alignment. */
  Trans3d &set_alignment(Vec3d from1, Vec3d from2, Vec3d from3, Vec3d to1,
                         Vec3d to2, Vec3d to3);

  /// Get a transformation that aligns two sets of three points
  /**\param from1 first point to align.
   * \param from2 second point to align.
   * \param from3 third point to align.
   * \param to1 align \p from1 with this point.
   * \param to2 align \p from2 with this point.
   * \param to3 align \p from3 with this point.
   * \return the alignment transformation. */
  static Trans3d alignment(Vec3d from1, Vec3d from2, Vec3d from3, Vec3d to1,
                           Vec3d to2, Vec3d to3);

  /// Set a transformation that aligns two sets of one, two or three points
  /**\param from one, two or three points to align.
   * \param to corresponding points for the \p from points to align with.
   * \return reference to this transformation set to the alignment. */
  Trans3d &set_alignment(std::vector<Vec3d> from, std::vector<Vec3d> to);

  /// Get a transformation that aligns two sets of one, two or three points
  /**\param from one, two or three points to align.
   * \param to corresponding points for the \p from points to align with.
   * \return the alignment transformation. */
  static Trans3d alignment(std::vector<Vec3d> from, std::vector<Vec3d> to);

  /// Set a transformation that makes particular angles between mapped axes.
  /** The x-axis is left unchanged. The y-axis is mapped to the y'-axis
   *  in the xy-plane. The z-axis is mapped to the z'-axis.
   * \param yz_ang the angle between the y'-axis and the z'axis.
   * \param zx_ang the angle between the z'-axis and the x-axis.
   * \param xy_ang the angle between the x-axis and the y'-axis.
   * \param valid to return whether a transformation based on
   *  the angles is valid.
   * \return  A reference to this matrix with the transformation set. */
  Trans3d &set_trans_by_angles(double yz_ang, double zx_ang, double xy_ang,
                               bool *valid = nullptr);

  /// Get a transformation that makes particular angles between mapped axes.
  /** The x-axis is left unchanged. The y-axis is mapped to the y'-axis
   *  in the xy-plane. The z-axis is mapped to the z'-axis.
   * \param zx_ang the angle between the z'-axis and the x-axis.
   * \param xy_ang the angle between the x-axis and the y'-axis.
   * \param yz_ang the angle between the y'-axis and the z'axis.
   * \param valid to return whether a transformation based on
   *  the angles is valid.
   * \return  A matrix with the transformation set. */
  static Trans3d trans_by_angles(double yz_ang, double zx_ang, double xy_ang,
                                 bool *valid = nullptr);

  /// Get the translation vector.
  /* \return the translation component of the transformation. */
  Vec3d get_transl() const;

  /// Get quaternion
  /** Convert the rotation part of the transformation into a quaterinon.
   * \return A vector representing a quaternion
   *  (components X, Y, Z, and W.) */
  Vec4d get_quaternion() const;

  /// Get Euler angles.
  /**\return A vector representing rotations about the X, Y and Z
   *  axis, in that order. */
  Vec3d get_euler() const;

  /// Debugging print of the matrix of a transformation.
  /**\param var a string to identify the matrix variable.
   * \param file file stream to print the variable. */
  void dump(const char *var = "", FILE *file = stderr) const;
};

/// Transform a column vector.
/**\param trans the transformation.
 * \param vec the column vector.
 * \return The transformed vector (left-multiplied by the
 *  transformation matrix). */
Vec3d operator*(const Trans3d &trans, const Vec3d &vec);

/// Transform a row vector.
/**\param vec the column vector.
 * \param trans the transformation matrix.
 * \return The transformed vector (right-multiplied by the
 *  transformation matrix). */
Vec3d operator*(const Vec3d &vec, const Trans3d &trans);

/// Compose two transformations
/**\param trans1 the first transformation.
 * \param trans2 the second transformation.
 * \return The composition of the two transformations (matrix of the first
 *  multiplying the matrix of the second). */
Trans3d operator*(const Trans3d &trans1, const Trans3d &trans2);

/// Compare two transformations (so they may be ordered)
/**\param trans1 the first transformation.
 * \param trans2 the second transformation.
 * \param eps a small number, numbers differing by less than eps are the same
 * \return  -1, 0, or 1 to indicate trans1 is less, equal or greater
 *  than trans2 */
int compare(const Trans3d &trans1, const Trans3d &trans2, double eps = epsilon);

/// Less for two transformations (so they may be ordered)
/**\param trans1 the first transformation.
 * \param trans2 the second transformation.
 * \return \c true if trans1 < trans 2, otherwise \c false.*/
bool operator<(const Trans3d &trans1, const Trans3d &trans2);

/// Transform a set of vectors
/**\param vecs the (column) vectors to transform.
 * \param trans the transformation to apply. */
void transform(std::vector<Vec3d> &vecs, const Trans3d &trans);

// inline functions
inline Trans3d::Trans3d(const Vec3d &r1, const Vec3d &r2, const Vec3d &r3)
{
  to_unit();
  const Vec3d *pvec[] = {&r1, &r2, &r3};
  for (int i = 0; i < 12; i++)
    if ((i % 4) != 3)
      m[i] = (*pvec[i / 4])[i % 4];
}

inline Trans3d::Trans3d(const Vec4d &r1, const Vec4d &r2, const Vec4d &r3,
                        Vec4d r4)
{
  const Vec4d *pvec[] = {&r1, &r2, &r3, &r4};
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      m[i] = (*pvec[i])[j];
}

inline Trans3d &Trans3d::to_unit()
{
  for (int i = 0; i < 16; i++)
    m[i] = ((i % 4) == (i / 4)) ? 1 : 0;
  return *this;
}

inline Trans3d Trans3d::unit()
{
  Trans3d trans;
  return trans.to_unit();
}

inline Trans3d &Trans3d::to_zero()
{
  for (auto &entry : m)
    entry = 0;
  return *this;
}

inline Trans3d Trans3d::zero()
{
  Trans3d trans;
  return trans.to_zero();
}

inline Trans3d &Trans3d::set_transl(Vec3d transl)
{
  to_unit();
  m[3] = transl[0];
  m[7] = transl[1];
  m[11] = transl[2];
  return *this;
}

inline Trans3d Trans3d::transl(Vec3d transl)
{
  Trans3d trans;
  return trans.set_transl(transl);
}

inline Trans3d Trans3d::rot(Vec3d axis, double angle)
{
  Trans3d trans;
  return trans.set_rot(axis, angle);
}

inline Trans3d &Trans3d::set_rot(double x_ang, double y_ang, double z_ang)
{
  return set_rot(Vec3d(0, 0, 1), z_ang) *=
         rot(Vec3d(0, 1, 0), y_ang) * rot(Vec3d(1, 0, 0), x_ang);
}

inline Trans3d Trans3d::rot(double x_ang, double y_ang, double z_ang)
{
  Trans3d trans;
  return trans.set_rot(x_ang, y_ang, z_ang);
}

inline Trans3d &Trans3d::set_rot(Vec3d angles)
{
  return set_rot(angles[0], angles[1], angles[2]);
}

inline Trans3d Trans3d::rot(Vec3d angles)
{
  Trans3d trans;
  return trans.set_rot(angles);
}

inline Trans3d &Trans3d::set_cols(const Vec3d &c1, const Vec3d &c2,
                                  const Vec3d &c3, const Vec3d &c4)
{
  const Vec3d *cols[] = {&c1, &c2, &c3, &c4};
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 3; j++)
      m[j * 4 + i] = (*cols[i])[j];
  return *this;
}

inline Trans3d Trans3d::rot(Vec3d v_from, Vec3d v_to)
{
  Trans3d trans;
  return trans.set_rot(v_from, v_to);
}

inline Trans3d Trans3d::refl(Vec3d norm)
{
  Trans3d trans;
  return trans.set_refl(norm);
}

inline Trans3d &Trans3d::set_scale(double scal)
{
  return set_scale(scal, scal, scal);
}

inline Trans3d Trans3d::scale(double scal)
{
  Trans3d trans;
  return trans.set_scale(scal);
}

inline Trans3d &Trans3d::set_scale(double x_scal, double y_scal, double z_scal)
{
  m[0] = x_scal;
  m[5] = y_scal;
  m[10] = z_scal;
  return *this;
}

inline Trans3d Trans3d::scale(double x_scal, double y_scal, double z_scal)
{
  Trans3d trans;
  return trans.set_scale(x_scal, y_scal, z_scal);
}

inline Trans3d &Trans3d::set_scale(Vec3d dir, double scal)
{
  Vec3d x(1, 0, 0);
  return set_rot(x, dir) *= scale(scal, 1, 1) * rot(dir, x);
}

inline Trans3d Trans3d::scale(Vec3d dir, double scal)
{
  Trans3d trans;
  return trans.set_scale(dir, scal);
}

inline Trans3d &Trans3d::set_inversion() { return set_scale(-1); }

inline Trans3d Trans3d::inversion()
{
  Trans3d trans;
  return trans.set_inversion();
}

inline Trans3d Trans3d::alignment(Vec3d from1, Vec3d from2, Vec3d to1,
                                  Vec3d to2)
{
  Trans3d trans;
  return trans.set_alignment(from1, from2, to1, to2);
}

inline Trans3d Trans3d::alignment(Vec3d from1, Vec3d from2, Vec3d from3,
                                  Vec3d to1, Vec3d to2, Vec3d to3)
{
  Trans3d trans;
  return trans.set_alignment(from1, from2, from3, to1, to2, to3);
}

inline Trans3d Trans3d::alignment(std::vector<Vec3d> from,
                                  std::vector<Vec3d> to)
{
  Trans3d trans;
  return trans.set_alignment(from, to);
}

inline Trans3d Trans3d::trans_by_angles(double yz_ang, double zx_ang,
                                        double xy_ang, bool *valid)
{
  Trans3d trans;
  return trans.set_trans_by_angles(yz_ang, zx_ang, xy_ang, valid);
}

inline Trans3d &Trans3d::transpose()
{
  for (int i = 0; i < 16; i++)
    if (i / 4 < i % 4) {
      double tmp = m[i];
      m[i] = m[(i % 4) * 4 + (i / 4)];
      m[(i % 4) * 4 + (i / 4)] = tmp;
    }
  return *this;
}

inline Trans3d Trans3d::inverse(const Trans3d &trans)
{
  Trans3d inv = trans;
  return inv.set_inverse();
}

inline double Trans3d::det() const
{
  return +m[0] * det2(m[5], m[6], m[9], m[10]) -
         m[1] * det2(m[4], m[6], m[8], m[10]) +
         m[2] * det2(m[4], m[5], m[8], m[9]);
}

inline Vec3d Trans3d::get_transl() const { return Vec3d(m[3], m[7], m[11]); }

inline Vec3d operator*(const Vec3d &vec, const Trans3d &trans)
{
  Trans3d m_ret = trans;
  m_ret.transpose();
  return m_ret * vec;
}

inline Trans3d operator*(const Trans3d &trans1, const Trans3d &trans2)
{
  Trans3d m_ret = trans1;
  return m_ret *= trans2;
}

inline bool operator<(const Trans3d &trans1, const Trans3d &trans2)
{
  return compare(trans1, trans2, epsilon) == -1;
}

inline bool operator==(const Trans3d &trans1, const Trans3d &trans2)
{
  for (int i = 0; i < 16; i++)
    if (!double_eq(trans1[i], trans2[i], epsilon))
      return false;
  return true;
}

inline void transform(std::vector<Vec3d> &vecs, const Trans3d &trans)
{
  for (auto &v : vecs)
    v = trans * v;
}

} // namespace anti

#endif // TRANS3D_H
