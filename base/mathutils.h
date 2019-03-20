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

/*!\file mathutils.h
   \brief utility routines for maths operations.
*/

#ifndef MATHUTILS_H
#define MATHUTILS_H

#include <algorithm>
#include <math.h>

#include "const.h"

#ifndef NAN
#define NAN (0.0F / 0.0F)
#endif

namespace anti {

/// Get the Greatest Common Divisor.
/**\param m the first number.
 * \param n the second number.
 * \return the greatest common divisor. */
long gcd(long m, long n);

/// Get a factorial.
/**\param n an integer.
 * \return the factorial of the integer. */
int factorial(int n);

/// Convert a floating point number to a close rational number.
/**Can easily fail if <code> MAX_LONG*eps^2 > 1 </code>.
 * \param f the floating point number.
 * \param num to return the numerator of the rational number.
 * \param denom to return the denominator of the rational number.
 * \param eps the maximum difference between the floating point
 * number and the rational approximation.
 * \param max_steps the maximum number of iterations, in case the process
 * does not converge (default 100 should exceed any reasonable number of steps).
 */
void double2rational(double f, long &num, long &denom, double eps = 0.00001,
                     int max_steps = 100);

/// Convert degrees to radians.
/**\param ang angle in degres, to convert.
 * \return The angle in radians. */
double deg2rad(double ang = 1.0);

/// Convert radians to degrees.
/**\param ang angle in radians, to convert.
 * \return The angle in degrees. */
double rad2deg(double ang = 1.0);

/// Check whether an integer is even.
/**\param n the number to check.
 * \return true if the number is even, otherwise false. \n */
bool is_even(int n);

/// Solve a quartic equation
/**\param coeffs the coefficients of the quartic, indexed by power of x.
 * \param sol the real parts of the solutions (in increasing order).
 * \param sol_i the imaginary parts of the solutions (if sol_i is not zero).
 * \return The number of real solutions. */
int quartic(double coeffs[5], double sol[4], double sol_i[4] = nullptr);

/// Solve a cubic equation
/**\param coeffs the coefficients of the cubic, indexed by power of x.
 * \param sol the real solutions (in ascending order).
 * \return The number of real solutions. */
int cubic(double coeffs[4], double sol[3]);

/// Determinant of a matrix
/**\param mat the (n*n) matrix coefficients (left to right, top to bottom)
 * \param n the number of rows/columns
 * \return the determinant, or NAN if the matrix is singular. */
double determinant(const double *mat, int n);

/// Compare two doubles for <, =, >
/**\param d1 first number
 * \param d2 second number
 * \param eps a small number, numbers differing by less than eps are the same
 * \return -1: d1<d2, 0:d1=d2, 1:d1>d2. */
int double_compare(double d1, double d2, double eps = epsilon);

/// Compare two doubles for equality
/**\param d1 first number
 * \param d2 second number
 * \param eps a small number, numbers differing by less than eps are the same
 * \return true if d1=d2, otherwise false */
bool double_eq(double d1, double d2, double eps = epsilon);

/// Compare two doubles for inequality
/**\param d1 first number
 * \param d2 second number
 * \param eps a small number, numbers differing by less than eps are the same
 * \return true if d1=/=d2, otherwise false */
bool double_ne(double d1, double d2, double eps = epsilon);

/// Compare two doubles for greater than
/**\param d1 first number
 * \param d2 second number
 * \param eps a small number, numbers differing by less than eps are the same
 * \return true if d1>d2, otherwise false */
bool double_gt(double d1, double d2, double eps = epsilon);

/// Compare two doubles for greater than or equal
/**\param d1 first number
 * \param d2 second number
 * \param eps a small number, numbers differing by less than eps are the same
 * \return true if d1>=d2, otherwise false */
bool double_ge(double d1, double d2, double eps = epsilon);

/// Compare two doubles for less than
/**\param d1 first number
 * \param d2 second number
 * \param eps a small number, numbers differing by less than eps are the same
 * \return true if d1<d2, otherwise false */
bool double_lt(double d1, double d2, double eps = epsilon);

/// Compare two doubles for less than or equal
/**\param d1 first number
 * \param d2 second number
 * \param eps a small number, numbers differing by less than eps are the same
 * \return true if d1<=d2, otherwise false */
bool double_le(double d1, double d2, double eps = epsilon);

/// Make a value safe to use as an argument with acos and asin
/**Map the value to the nearest value in the range
 * -1.0<=val<=1.0 to ensure that it is safe to use with
 *  \c acos() and \c asin().
 * \param val value to make safe.
 * \return A safe value. */
double safe_for_trig(double val);

/// Get sign, as -1, 1
/**\param val the value to get the sign of
 * \return 1 for x>=0 else -1 */
template <typename T> int sgn(T val) { return (val < 0) ? -1 : 1; }

// inline function definitions

inline double deg2rad(double ang) { return ang * M_PI / 180; }

inline double rad2deg(double ang) { return ang * 180 / M_PI; }

inline bool is_even(int n) { return (n % 2 == 0); }

inline double safe_for_trig(double val)
{
  return (std::min(1.0, std::max(-1.0, val)));
}

// considering epsilon, return 0 if equal, -1 if d1 < d2, 1 if d1 > d2
inline int double_compare(double d1, double d2, double eps)
{
  const double diff = d1 - d2;
  return (diff < eps ? (diff > -eps ? 0 : -1) : 1); // excludes epsilon as zero
}

// true if d1 == d2 are considering epsilon
inline bool double_eq(double d1, double d2, double eps)
{
  // const double diff = d1 - d2;
  // return diff < eps && diff > -eps; // excludes epsilon as zero
  return (!double_compare(d1, d2, eps));
}

// true if d1 != d2 considering epsilon
inline bool double_ne(double d1, double d2, double eps)
{
  return (double_compare(d1, d2, eps));
}

// true if d1 > d2 considering epsilon
inline bool double_gt(double d1, double d2, double eps)
{
  return (double_compare(d1, d2, eps) > 0);
}

// true if d1 >= d2 considering epsilon
inline bool double_ge(double d1, double d2, double eps)
{
  return (double_compare(d1, d2, eps) > -1);
}

// true if d1 < d2 considering epsilon
inline bool double_lt(double d1, double d2, double eps)
{
  return (double_compare(d1, d2, eps) < 0);
}

// true if d1 >= d2 considering epsilon
inline bool double_le(double d1, double d2, double eps)
{
  return (double_compare(d1, d2, eps) < 1);
}

} // namespace anti

#endif // MATHUTILS_H
