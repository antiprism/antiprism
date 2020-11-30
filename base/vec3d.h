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

/*! \file vec3d.h
 *  \brief Vector for 3D geometry
 *
 *  A vector class with common vector operations.
 */

#ifndef VEC3D_H
#define VEC3D_H

#include "const.h"
#include "mathutils.h"
#include "random.h"
#include "status.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace anti {

/// Vector with 3 components
class Vec3d {
private:
  double v[3];

public:
  /// Constructor
  /** The vector is initialised to the unset state */
  Vec3d() { unset(); }

  /// Constructor
  /**\param x_val first component
   * \param y_val second component
   * \param z_val third component */
  Vec3d(double x_val, double y_val, double z_val);

  /// Constructor
  /** Initialise from an array
   * \param vals pointer to an array of at least three values
   * to use as the x, y and z components */
  Vec3d(double *vals);

  /// Constructor
  /** Initialise from an array
   * \param vals pointer to an array of at least three values
   * to use as the x, y and z components */
  Vec3d(float *vals);

  /// Convert into a unit vector with the same direction.
  /**\return A reference to the vector. */
  Vec3d &to_unit();

  /// Get a unit vector with the same direction.
  /**\return The unit vector. */
  Vec3d unit() const;

  /// Get a vector in the same direction with a given length.
  /**\param length length for the vector
   * \return A vector in the same direction with the given length. */
  Vec3d with_len(double length);

  /// Get the length of the vector
  /**\return The length. */
  double len() const;

  /// Get the square of the length of the vector
  /**\return The square of the length. */
  double len2() const;

  /// Multiply this vector by a scalar
  /**\param n number to multiply by
   * \return A reference to this vector. */
  Vec3d &operator*=(double n);

  /// Divide this vector by a scalar
  /**\param n number to divide by
   * \return A reference to this vector. */
  Vec3d &operator/=(double n);

  /// Add a vector to this vector
  /**\param vec the vector to add
   * \return A reference to this vector. */
  Vec3d &operator+=(Vec3d vec);

  /// Subtract a vector from this vector
  /**\param vec the vector to subtract
   * \return A reference to this vector. */
  Vec3d &operator-=(Vec3d vec);

  /// Get the vector components
  /**\param idx the component index.
   * \return The value of the component. */
  double operator[](int idx) const { return v[idx]; }

  /// Get the vector components
  /**\param idx the component index.
   * \return A reference to the component. */
  double &operator[](int idx) { return v[idx]; }

  /// Get the x component
  /**\return The value of the x component. */
  double x() const { return v[0]; }

  /// Set the x component
  /**\param val the value of the x component. */
  void x(double val) { v[0] = val; }

  /// Get the y component
  /**\return The value of the y component. */
  double y() const { return v[1]; }

  /// Set the y component
  /**\param val the value of the y component. */
  void y(double val) { v[1] = val; }

  /// Get the z component
  /**\return The value of the z component. */
  double z() const { return v[2]; }

  /// Set the z component
  /**\param val the value of the z component. */
  void z(double val) { v[2] = val; }

  /// Get a pointer to the vector component array.
  /** \return A pointer to the underlying component array. */
  const double *get_v() const { return v; }

  /// Get the component along another vector.
  /**\param along component direction
   * \return The component vector. */
  Vec3d component(Vec3d along);

  /// Unset the vector.
  /**Put the vector into the initial unset state. The vector will return
   * \c false if tested */
  void unset();

  /// Get a random vector.
  // Uses random numbers provided by the Random argument.
  /**\param rnd to provide the random numbers
   * \return A random vector with length less then or equal to one. */
  static Vec3d random(Random &rnd);

  /// Check whether a vector has been set
  /**\return \c true if set, otherwise \c false */
  bool is_set() const { return !std::isnan(v[0]); }

  /// Read a vector from a string
  /**\param str a string containing three decimals separated by
   * commas or spaces.
   * \return status, which evaluates to \c true if a valid vector was read,
   * otherwise \c false to indictate an error. */
  Status read(const char *str);

  /// Print a vector variable
  /**\return a string representation of the variable */
  // std::string str() const;

  /// Convert to a coordinate string
  /**\param sep the separator between the numbers.
   * \param sig_dgts the number of significant digits in the conversion,
   *  or if negative then the number of digits after the decimal point.
   * \return The string. */
  std::string to_str(const char *sep = ", ", int sig_dgts = 17) const;

  /// Debugging print of a vector variable
  /**\param var a string to identify the vector variable.
   * \param file file stream to print the variable. */
  void dump(const char *var = "", FILE *file = stderr) const;

  static Vec3d X;    //<Unit vector in the direction of the x-axis
  static Vec3d Y;    //<Unit vector in the direction of the y-axis
  static Vec3d Z;    //<Unit vector in the direction of the z-axis
  static Vec3d zero; //<Zero vector
};

/// Add two vectors
/**\param vec1 a vector
 * \param vec2 a vector to add
 * \return The resulting vector (\c vec1 + \c vec2). */
Vec3d operator+(Vec3d vec1, Vec3d vec2);

/// Subtract one vector from another
/**\param vec1 a vector
 * \param vec2 a vector to subtract
 * \return The resulting vector (\c vec1 - \c vec2). */
Vec3d operator-(Vec3d vec1, Vec3d vec2);

/// The negative of a vector
/**\param vec a vector
 * \return The negative of the vector (-\c v). */
Vec3d operator-(Vec3d vec);

/// Multiply a vector by a scalar
/**\param vec the vector
 * \param n the scalar
 * \return The resulting vector (\c n * \c v). */
Vec3d operator*(Vec3d vec, double n);

/// Multiply a vector by a scalar
/**\param n the scalar
 * \param vec the vector
 * \return The resulting vector (\c n * \c v). */
Vec3d operator*(double n, Vec3d vec);

/// Divide a vector by a scalar
/**\param vec the vector
 * \param n the scalar
 * \return The resulting vector (1/\c n * \c v). */
Vec3d operator/(Vec3d vec, double n);

/// The cross product (vector product)
/**\param vec1 the first vector
 * \param vec2 the second vector
 * \return The cross product (\c vec1 x \c vec2). */
inline Vec3d vcross(const Vec3d &vec1, const Vec3d &vec2)
{
  Vec3d vprod;
  vprod[0] = vec2[2] * vec1[1] - vec2[1] * vec1[2];
  vprod[1] = vec2[0] * vec1[2] - vec2[2] * vec1[0];
  vprod[2] = vec2[1] * vec1[0] - vec2[0] * vec1[1];
  return vprod;
}

/// The dot product (scalar product)
/**\param vec1 the first vector
 * \param vec2 the second vector
 * \return The dot product (\c vec1 . \c vec2). */
inline double vdot(const Vec3d &vec1, const Vec3d &vec2)
{
  return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

/// The triple product
/**\param vec1 the first vector
 * \param vec2 the second vector
 * \param vec3 the third vector
 * \return The dot product (\c vec1 . (\c vec2 x \c vec3). */
inline double vtriple(Vec3d vec1, Vec3d vec2, Vec3d vec3)
{
  return vdot(vec1, vcross(vec2, vec3));
}

inline Vec3d::Vec3d(double x_val, double y_val, double z_val)
{
  v[0] = x_val;
  v[1] = y_val;
  v[2] = z_val;
}

inline Vec3d::Vec3d(double *vals)
{
  v[0] = vals[0];
  v[1] = vals[1];
  v[2] = vals[2];
}

inline Vec3d::Vec3d(float *vals)
{
  v[0] = vals[0];
  v[1] = vals[1];
  v[2] = vals[2];
}

inline Vec3d Vec3d::random(Random &rnd)
{
  Vec3d u;
  do {
    u[0] = 1.0 - 2.0 * rnd.ranf();
    u[1] = 1.0 - 2.0 * rnd.ranf();
    u[2] = 1.0 - 2.0 * rnd.ranf();
  } while (u.len2() > 1);
  return u;
}

inline Vec3d Vec3d::unit() const
{
  Vec3d ret = *this;
  return ret.to_unit();
}

inline Vec3d &Vec3d::to_unit()
{
  double ln = len();
  if (ln > 1e-20)
    operator*=(1 / ln);
  else {
    v[0] = 0;
    v[1] = 0;
    v[2] = 1;
  }
  return *this;
}

inline Vec3d Vec3d::with_len(double length) { return unit() * length; }

inline double Vec3d::len2() const
{
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

inline double Vec3d::len() const { return sqrt(len2()); }

inline Vec3d Vec3d::component(Vec3d along)
{
  along.to_unit();
  return vdot(*this, along) * along;
}

inline Vec3d &Vec3d::operator*=(double n)
{
  v[0] *= n;
  v[1] *= n;
  v[2] *= n;
  return *this;
}

inline Vec3d &Vec3d::operator/=(double n) { return *this *= (1.0 / n); }

inline Vec3d &Vec3d::operator+=(Vec3d vec)
{
  v[0] += vec.v[0];
  v[1] += vec.v[1];
  v[2] += vec.v[2];
  return *this;
}

inline Vec3d &Vec3d::operator-=(Vec3d vec)
{
  v[0] -= vec.v[0];
  v[1] -= vec.v[1];
  v[2] -= vec.v[2];
  return *this;
}

inline Vec3d operator+(Vec3d vec1, Vec3d vec2)
{
  Vec3d ret = vec1;
  return ret += vec2;
}

inline Vec3d operator-(Vec3d vec1, Vec3d vec2)
{
  Vec3d ret = vec1;
  return ret -= vec2;
}

inline Vec3d operator-(Vec3d vec)
{
  Vec3d ret = vec;
  return ret *= -1;
}

inline Vec3d operator*(Vec3d vec, double n)
{
  Vec3d ret = vec;
  return ret *= n;
}

inline Vec3d operator*(double n, Vec3d vec)
{
  Vec3d ret = vec;
  return ret *= n;
}

inline Vec3d operator/(Vec3d vec, double n)
{
  Vec3d ret = vec;
  return ret /= n;
}

inline int compare(const Vec3d &vec1, const Vec3d &vec2, double eps = epsilon)
{
  if (!vec1.is_set() && !vec2.is_set())
    return 0;
  if (!vec1.is_set())
    return -1;
  if (!vec2.is_set())
    return 1;
  for (int i = 0; i < 3; i++) {
    /*
          if(fabs(vec1[i]-vec2[i])>eps) {
             if(vec1[i]<vec2[i])
                return -1;
             if(vec1[i]>vec2[i])
                return 1;
          }
    */
    int ret = double_compare(vec1[i], vec2[i], eps);
    if (ret)
      return ret;
  }
  return 0;
}

} // namespace anti

#endif // VEC3D_H
