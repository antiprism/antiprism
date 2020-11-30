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

/*! \file vec4d.h
 *  \brief Vector for 4D geometry
 *
 *  A vector class with common vector operations.
 */

#ifndef VEC4D_H
#define VEC4D_H

#include "random.h"
#include "status.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

namespace anti {

/// Vector with 4 components
class Vec4d {
private:
  double v[4];

public:
  /// Constructor
  /** The vector is initialised to the unset state */
  Vec4d() { unset(); }

  /// Constructor
  /**\param x first component
   * \param y second component
   * \param z third component
   * \param w fourth component */
  Vec4d(double x, double y, double z, double w);

  /// Convert into a unit vector with the same direction.
  /**\return A reference to the vector. */
  Vec4d &to_unit();

  /// Get a unit vector with the same direction.
  /**\return The unit vector. */
  Vec4d unit() const;

  /// Get a vector in the same direction with a given length.
  /**\param length length for the vector
   * \return A vector in the same direction with the given length. */
  Vec4d with_len(double length);

  /// Get the length of the vector
  /**\return The length. */
  double len() const;

  /// Get the square of the length of the vector
  /**\return The square of the length. */
  double len2() const;

  /// Multiply this vector by a scalar
  /**\param n number to multiply by
   * \return A reference to this vector. */
  Vec4d &operator*=(double n);

  /// Divide this vector by a scalar
  /**\param n number to divide by
   * \return A reference to this vector. */
  Vec4d &operator/=(double n);

  /// Add a vector to this vector
  /**\param vec the vector to add
   * \return A reference to this vector. */
  Vec4d &operator+=(Vec4d vec);

  /// Subtract a vector from this vector
  /**\param vec the vector to subtract
   * \return A reference to this vector. */
  Vec4d &operator-=(Vec4d vec);

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

  /// Get the w component
  /**\return The value of the z component. */
  double w() const { return v[3]; }

  /// Set the w component
  /**\param val the value of the w component. */
  void w(double val) { v[3] = val; }

  /// Get a pointer to the vector component array.
  /** \return A pointer to the underlying component array. */
  const double *get_v() const { return v; }

  /// Unset the vector.
  /**Put the vector into the initial unset state. The vector will return
   * \c false if tested */
  void unset();

  /// Get a random vector.
  // Uses random numbers provided by the Random argument.
  /**\param rnd to provide the random numbers
   * \return A random vector with length less then or equal to one. */
  static Vec4d random(Random &rnd);

  /// Check whether a vector has been set
  /**\return \c true if set, otherwise \c false */
  bool is_set() const { return !std::isnan(v[0]); }

  /// Read a vector from a string
  /**\param str a string containing three decimals separated by
   * commas or spaces.
   * \return status, which evaluates to \c true if a valid vector was read,
   * otherwise \c false to indictate an error. */
  Status read(const char *str);

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

  static Vec4d X;
  static Vec4d Y;
  static Vec4d Z;
  static Vec4d W;
  static Vec4d zero;
};

Vec4d operator+(Vec4d vec1, Vec4d vec2);
Vec4d operator-(Vec4d vec1, Vec4d vec2);
Vec4d operator-(Vec4d vec);
Vec4d operator*(Vec4d vec, double n);
Vec4d operator*(double n, Vec4d vec);
Vec4d operator/(Vec4d vec, double n);

/// The dot product (scalar product)
/**\param vec1 the first vector
 * \param vec2 the second vector
 * \return The dot product (\c vec1 . \c vec2). */
inline double vdot(const Vec4d &vec1, const Vec4d &vec2)
{
  return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2] +
         vec1[3] * vec2[3];
}

/// The "cross" product (vector product)
/**\param vec1 the first vector
 * \param vec2 the second vector
 * \param vec3 the third vector
 * \return A vector normal to the three vectors. */
Vec4d vcross(const Vec4d &vec1, const Vec4d &vec2, const Vec4d &vec3);

// inline member functions
inline Vec4d::Vec4d(double x, double y, double z, double w)
{
  v[0] = x;
  v[1] = y;
  v[2] = z;
  v[3] = w;
}

inline Vec4d Vec4d::random(Random &rnd)
{
  Vec4d u;
  do {
    u[0] = 1.0 - 2.0 * rnd.ranf();
    u[1] = 1.0 - 2.0 * rnd.ranf();
    u[2] = 1.0 - 2.0 * rnd.ranf();
    u[3] = 1.0 - 2.0 * rnd.ranf();
  } while (u.len2() > 1);
  return u;
}

inline Vec4d Vec4d::unit() const
{
  Vec4d ret = *this;
  return ret.to_unit();
}

inline Vec4d &Vec4d::to_unit()
{
  double ln = len();
  if (ln > 1e-20)
    operator*=(1 / ln);
  else {
    v[0] = 0;
    v[1] = 0;
    v[2] = 0;
    v[3] = 1;
  }
  return *this;
}

inline Vec4d Vec4d::with_len(double length) { return unit() * length; }

inline double Vec4d::len2() const
{
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3];
}

inline double Vec4d::len() const { return sqrt(len2()); }

inline Vec4d &Vec4d::operator*=(double n)
{
  v[0] *= n;
  v[1] *= n;
  v[2] *= n;
  v[3] *= n;
  return *this;
}

inline Vec4d &Vec4d::operator/=(double n) { return *this *= (1.0 / n); }

inline Vec4d &Vec4d::operator+=(Vec4d vec)
{
  v[0] += vec.v[0];
  v[1] += vec.v[1];
  v[2] += vec.v[2];
  v[3] += vec.v[3];
  return *this;
}

inline Vec4d &Vec4d::operator-=(Vec4d vec)
{
  v[0] -= vec.v[0];
  v[1] -= vec.v[1];
  v[2] -= vec.v[2];
  v[3] -= vec.v[3];
  return *this;
}

inline Vec4d operator+(Vec4d vec1, Vec4d vec2)
{
  Vec4d ret = vec1;
  return ret += vec2;
}

inline Vec4d operator-(Vec4d vec1, Vec4d vec2)
{
  Vec4d ret = vec1;
  return ret -= vec2;
}

inline Vec4d operator-(Vec4d vec)
{
  Vec4d ret = vec;
  return ret *= -1;
}

inline Vec4d operator*(Vec4d vec, double n)
{
  Vec4d ret = vec;
  return ret *= n;
}

inline Vec4d operator*(double n, Vec4d vec)
{
  Vec4d ret = vec;
  return ret *= n;
}

inline Vec4d operator/(Vec4d vec, double n)
{
  Vec4d ret = vec;
  return ret /= n;
}

} // namespace anti

#endif // VEC4D_H
