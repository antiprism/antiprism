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

/*!\file utils.h
   \brief utility routines for maths operations.
*/


#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <algorithm>
#include <math.h>

#include "const.h"

///Get the Greatest Common Divisor.
/**\param m the first number.
 * \param n the second number.
 * \return the greatest common divisor. */
long gcd(long m, long n);

///Get a factorial.
/**\param n an integer.
 * \return the factorial of the integer. */
int factorial(int n);

///Convert a floating point number to a close rational number.
/**Can easily fail if <code> MAX_LONG*eps^2 > 1 </code>.
 * \param f the floating point number.
 * \param num to return the numerator of the rational number.
 * \param denom to return the denominator of the rational number.
 * \param eps the maximum difference between the floating point
 * number and the rational approximation. */
void double2rational(double f, long &num, long &denom, double eps = 0.00001);

///Convert degrees to radians.
/**\param ang angle in degres, to convert.
 * \return The angle in radians. */
double deg2rad(double ang=1.0);

///Convert radians to degrees.
/**\param ang angle in radians, to convert.
 * \return The angle in degrees. */
double rad2deg(double ang=1.0);

///Check whether an integer is even.
/**\param n the number to check.
 * \return true if the number is even, otherwise false. \n */
bool is_even(int n);

///Solve a quartic equation
/**\param coeffs the coefficients of the quartic, indexed by power of x.
 * \param sol the real parts of the solutions (in increasing order).
 * \param sol_i the imaginary parts of the solutions (if sol_i is not zero).
 * \return The number of real solutions. */
int quartic(double coeffs[5], double sol[4], double sol_i[4]=0);

///Solve a cubic equation
/**\param coeffs the coefficients of the cubic, indexed by power of x.
 * \param sol the real solutions (in ascending order).
 * \return The number of real solutions. */
int cubic(double coeffs[4], double sol[3]);

bool double_equality(const double &d1, const double &d2,
      const double &eps=epsilon);

///Make a value safe to use as an argument with acos and asin
/**Map the value to the nearest value in the range
 * -1.0<=val<=1.0 to ensure that it is safe to use with
 *  \c acos() and \c asin().
 * \param val value to make safe.
 * \return A safe value. */
double safe_for_trig(double val);

// inline function definitions

inline double deg2rad(double ang)
{ 
   return ang * M_PI/180;
}

inline double rad2deg(double ang)
{ 
   return ang * 180/M_PI;
}

inline bool is_even(int n)
{
   return (n%2==0);
}

inline double safe_for_trig(double val) {
   return (std::min(1.0, std::max(-1.0, val)));
}

inline bool double_equality(const double &d1, const double &d2,
      const double &eps)
{
   const double diff = d1 - d2;
   return diff < eps && diff > -eps;
}



#endif // MATH_UTILS_H

