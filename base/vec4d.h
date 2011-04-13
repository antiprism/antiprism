/*
   Copyright (c) 2003, Adrian Rossiter

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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rand_gen.h"

///Vector with 4 components
class vec4d
{
   public:
     double v[4];

     ///Constructor
     /** The vector is initialised to the unset state */
     vec4d() { unset(); }
     
     ///Constructor
     /**\param x first component 
      * \param y second component 
      * \param z third component 
      * \param w fourth component */
     vec4d(double x, double y, double z, double w);

     ///Convert into a unit vector with the same direction.
     /**\return A reference to the vector. */
     vec4d &to_unit();

     ///Get a unit vector with the same direction.
     /**\return The unit vector. */
     vec4d unit() const;

     ///Get the magnitude of the vector
     /**\return The magnitude. */
     double mag() const;

     ///Get the square of the magnitude of the vector
     /**\return The square of the magnitude. */
     double mag2() const;
     
     ///Multiply this vector by a scalar
     /**\return A reference to this vector. */
     vec4d &operator *=(double n);
     
     ///Divide this vector by a scalar
     /**\return A reference to this vector. */
     vec4d &operator /=(double n);           // divide by scalar n
     
     ///Add a vector to this vector
     /**\return A reference to this vector. */
     vec4d &operator +=(vec4d v);				// add vector v
     
     ///Subtract a vector from this vector
     /**\return A reference to this vector. */
     vec4d &operator -=(vec4d v);
 

     ///Get the vector components
     /**\param idx the component index.
      * \return The value of the component. */
     inline double operator [](int idx) const { return v[idx]; }
 
     ///Get the vector components
     /**\param idx the component index.
      * \return A reference to the component. */
     inline double &operator [](int idx) { return v[idx]; }
           
     ///Unset the vector.
     /**Put the vector into the initial unset state. The vector will return
      * \c false if tested */
     void unset();
     
     ///Get a random vector.
     /**\return A random vector with magnitude less then or equal to one. */
     static vec4d random();
     
     ///Get a random vector.
     //Uses random numbers provided by the rand_gen argument.
     /**\return A random vector with magnitude less then or equal to one. */
     static vec4d random(rand_gen &rnd);
     
     ///Check whether a vector has been set
     /**\return \c true if set, otherwise \c false */
     bool is_set() const { return !isnan(v[0]); }
     
     ///Read a vector from a string
     /**\param str a string containing three decimals separated by
      * commas or spaces.
      * \param errmsg an array at least \c MSG_SZ chars long to
      * return any error message.
      * \return true if a valid vector was read, otherwise false
      * and the error is detailed in \a errmsg. */
     bool read(const char *str, char *errmsg=0);
     

     ///Debugging print of a vector variable
     /**\param var a string to identify the vector variable.
      * \param file file stream to print the variable. */
     void dump(const char *var="", FILE *file=stderr) const;

};

vec4d operator +(vec4d v1, vec4d v2);
vec4d operator -(vec4d v1, vec4d v2);
vec4d operator -(vec4d v1);
vec4d operator *(vec4d v1, double n);
vec4d operator *(double n, vec4d v1);
vec4d operator /(vec4d v1, double n);


///The dot product (scalar product)
/**\param v1 the first vector
 * \param v2 the second vector
 * \return The dot product (\c v1 . \c v2). */
inline double vdot(const vec4d &v1, const vec4d &v2)
{
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2] + v1[3]*v2[3];
}


///The "cross" product (vector product)
/**\param v1 the first vector
 * \param v2 the second vector
 * \param v3 the third vector
 * \return A vector normal to the three vectors. */
vec4d vcross(const vec4d &v1, const vec4d &v2, const vec4d &v3);

// inline member functions
inline vec4d::vec4d(double i, double j, double k, double l)
{
	v[0] = i;
	v[1] = j;
	v[2] = k;
	v[3] = l;
}

inline vec4d vec4d::random()
{
   vec4d u;
   do {
      u[0] = 1.0 - 2.0*rand()/(RAND_MAX-1);
      u[1] = 1.0 - 2.0*rand()/(RAND_MAX-1);
      u[2] = 1.0 - 2.0*rand()/(RAND_MAX-1);
      u[3] = 1.0 - 2.0*rand()/(RAND_MAX-1);
   } while(u.mag2()>1);
	return u;
}

inline vec4d vec4d::random(rand_gen &rnd)
{
   vec4d u;
   do {
      u[0] = 1.0 - 2.0*rnd.ranf();
      u[1] = 1.0 - 2.0*rnd.ranf();
      u[2] = 1.0 - 2.0*rnd.ranf();
      u[3] = 1.0 - 2.0*rnd.ranf();
   } while(u.mag2()>1);
	return u;
}

inline vec4d vec4d::unit() const
{
   vec4d ret = *this;
   return ret.to_unit();
}

inline vec4d& vec4d::to_unit()
{
	double magn = mag();
	if(magn > 1e-20)
		operator *=(1/magn);
	else {
		v[0] = 0;
		v[1] = 0;
		v[2] = 0;
		v[3] = 1;
	}
	return *this;
}


inline double vec4d::mag2() const
{
	return v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3];
}


inline double vec4d::mag() const
{
	return sqrt(mag2());
}


inline vec4d &vec4d::operator *=(double n)
{
	v[0] *= n;
	v[1] *= n;
	v[2] *= n;
	v[3] *= n;
	return *this;
}

inline vec4d &vec4d::operator /=(double n)
{
   return *this *=(1.0/n);
}

inline vec4d &vec4d::operator +=(vec4d v2)
{
	v[0] += v2.v[0];
	v[1] += v2.v[1];
	v[2] += v2.v[2];
	v[3] += v2.v[3];
	return *this;
}


inline vec4d &vec4d::operator -=(vec4d v2)
{
	v[0] -= v2.v[0];
	v[1] -= v2.v[1];
	v[2] -= v2.v[2];
	v[3] -= v2.v[3];
	return *this;
}



inline vec4d operator +(vec4d v1, vec4d v2)
{
   vec4d ret = v1;
   return ret += v2;
}

inline vec4d operator -(vec4d v1, vec4d v2)
{ 
   vec4d ret = v1;
   return ret -= v2;
}

inline vec4d operator -(vec4d v1)
{ 
   vec4d ret = v1;
   return ret *= -1;
}

inline vec4d operator *(vec4d v1, double n)
{ 
   vec4d ret = v1;
   return ret *= n; 
}

inline vec4d operator *(double n, vec4d v1)
{
   vec4d ret = v1;
   return ret *= n;
}

inline vec4d operator /(vec4d v1, double n)
{ 
   vec4d ret = v1;
   return ret /= n; 
}


#endif // VEC4D_H

