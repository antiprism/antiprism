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

/*! \file vec3d.h
 *  \brief Vector for 3D geometry
 *
 *  A vector class with common vector operations.
 */

#ifndef VEC3D_H
#define VEC3D_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "const.h"

///Vector with 3 components
class vec3d
{
   private:
     double v[3];

   public:
     ///Constructor
     /** The vector is initialised to the unset state */
     vec3d() { unset(); }

     ///Constructor
     /**\param x_val first component 
      * \param y_val second component 
      * \param z_val third component */
     vec3d(double x_val, double y_val, double z_val);
     
     ///Constructor
     /** Initialise from an array
      * \param vals pointer to an array of at least three values
      * to use as the x, y and z components */
     vec3d(double *vals);
     
     ///Constructor
     /** Initialise from an array
      * \param vals pointer to an array of at least three values
      * to use as the x, y and z components */
     vec3d(float *vals);

     ///Convert into a unit vector with the same direction.
     /**\return A reference to the vector. */
     vec3d &to_unit();

     ///Get a unit vector with the same direction.
     /**\return The unit vector. */
     vec3d unit() const;

     ///Get the magnitude of the vector
     /**\return The magnitude. */
     double mag() const;

     ///Get the square of the magnitude of the vector
     /**\return The square of the magnitude. */
     double mag2() const;
     
     ///Multiply this vector by a scalar
     /**\return A reference to this vector. */
     vec3d &operator *=(double n);
     
     ///Divide this vector by a scalar
     /**\return A reference to this vector. */
     vec3d &operator /=(double n);           // divide by scalar n
     
     ///Add a vector to this vector
     /**\return A reference to this vector. */
     vec3d &operator +=(vec3d v);				// add vector v
     
     ///Subtract a vector from this vector
     /**\return A reference to this vector. */
     vec3d &operator -=(vec3d v);
 
     ///Get the vector components
     /**\param idx the component index.
      * \return The value of the component. */
     inline double operator [](int idx) const { return v[idx]; }
 
     ///Get the vector components
     /**\param idx the component index.
      * \return A reference to the component. */
     inline double &operator [](int idx) { return v[idx]; }
     
     ///Get a pointer to the vector component array.
     /** \return A pointer to the underlying component array. */
     const double *get_v() const { return v; }
     
     ///Get the component along another vector.
     /**\return The component vector. */
     vec3d component(vec3d along);
     
     ///Unset the vector.
     /**Put the vector into the initial unset state. The vector will return
      * \c false if tested */
     void unset();
     
     ///Get a random vector.
     /**\return A random vector with magnitude less then or equal to one. */
     static vec3d random();
     
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
     bool read(char *str, char *errmsg=0);
     

     ///Debugging print of a vector variable
     /**\param var a string to identify the vector variable.
      * \param file file stream to print the variable. */
     void dump(const char *var="", FILE *file=stderr) const;

     static vec3d x;
     static vec3d y;
     static vec3d z;
     static vec3d zero;

};



///Add two vectors
/**\param v1 a vector
 * \param v2 a vector to add
 * \return The resulting vector (\c v1 + \c v2). */
vec3d operator +(vec3d v1, vec3d v2);

///Subtract one vector from another
/**\param v1 a vector
 * \param v2 a vector to subtract
 * \return The resulting vector (\c v1 - \c v2). */
vec3d operator -(vec3d v1, vec3d v2);

///The negative of a vector
/**\param v a vector
 * \return The negative of the vector (-\c v). */
vec3d operator -(vec3d v);

///Multiply a vector by a scalar
/**\param v the vector
 * \param n the scalar
 * \return The resulting vector (\c n * \c v). */
vec3d operator *(vec3d v, double n);

///Multiply a vector by a scalar
/**\param n the scalar
 * \param v the vector
 * \return The resulting vector (\c n * \c v). */
vec3d operator *(double n, vec3d v);

///Divide a vector by a scalar
/**\param v the vector
 * \param n the scalar
 * \return The resulting vector (1/\c n * \c v). */
vec3d operator /(vec3d v, double n);

///The cross product (vector product)
/**\param v1 the first vector
 * \param v2 the second vector
 * \return The cross product (\c v1 x \c v2). */
inline vec3d vcross(const vec3d &v1, const vec3d &v2)
{  
	vec3d vprod;
	vprod[0] = v2[2]*v1[1] - v2[1]*v1[2];
	vprod[1] = v2[0]*v1[2] - v2[2]*v1[0];
	vprod[2] = v2[1]*v1[0] - v2[0]*v1[1];
	return vprod;
}

///The dot product (scalar product)
/**\param v1 the first vector
 * \param v2 the second vector
 * \return The dot product (\c v1 . \c v2). */
inline double vdot(const vec3d &v1, const vec3d &v2)
{
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

///The triple product
/**\param v1 the first vector
 * \param v2 the second vector
 * \param v3 the third vector
 * \return The dot product (\c v1 . (\c v2 x \c v3). */
inline double vtriple(vec3d v1, vec3d v2, vec3d v3)
{
   return vdot(v1, vcross(v2, v3));
}


inline vec3d::vec3d(double x_val, double y_val, double z_val)
{
	v[0] = x_val;
	v[1] = y_val;
	v[2] = z_val;
}

inline vec3d::vec3d(double *vals)
{
	v[0] = vals[0];
	v[1] = vals[1];
	v[2] = vals[2];
}

inline vec3d::vec3d(float *vals)
{
	v[0] = vals[0];
	v[1] = vals[1];
	v[2] = vals[2];
}

inline vec3d vec3d::random()
{
   vec3d u;
   do {
      u[0] = 1.0 - 2.0*rand()/(RAND_MAX-1);
      u[1] = 1.0 - 2.0*rand()/(RAND_MAX-1);
      u[2] = 1.0 - 2.0*rand()/(RAND_MAX-1);
   } while(u.mag2()>1);
	return u;
}


inline vec3d vec3d::unit() const
{
   vec3d ret = *this;
   return ret.to_unit();
}


inline vec3d& vec3d::to_unit()
{
   /*
	double magn = mag();
	operator *=(1/magn);
   for(int i=0; i<3; i++)
	   if(isinf(v[i]))
         return x;   // return some unit vector
	return *this;
   */

   double magn = mag();
   if(magn > 1e-20)
      operator *=(1/magn);
   else {
      v[0] = 0;
      v[1] = 0;
      v[2] = 1;
   }
   return *this;
}


inline double vec3d::mag2() const
{
	return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}


inline double vec3d::mag() const
{
	return sqrt(mag2());
}

inline vec3d vec3d::component(vec3d along)
{
   along.to_unit();
   return vdot(*this, along)*along;
}

inline vec3d &vec3d::operator *=(double n)
{
	v[0] *= n;
	v[1] *= n;
	v[2] *= n;
	return *this;
}


inline vec3d &vec3d::operator /=(double n)
{
   return *this *=(1.0/n);
}


inline vec3d &vec3d::operator +=(vec3d v2)
{
	v[0] += v2.v[0];
	v[1] += v2.v[1];
	v[2] += v2.v[2];
	return *this;
}


inline vec3d &vec3d::operator -=(vec3d v2)
{
	v[0] -= v2.v[0];
	v[1] -= v2.v[1];
	v[2] -= v2.v[2];
	return *this;
}



inline vec3d operator +(vec3d v1, vec3d v2)
{
   vec3d ret = v1;
   return ret += v2;
}

inline vec3d operator -(vec3d v1, vec3d v2)
{ 
   vec3d ret = v1;
   return ret -= v2;
}

inline vec3d operator -(vec3d v1)
{ 
   vec3d ret = v1;
   return ret *= -1;
}

inline vec3d operator *(vec3d v1, double n)
{ 
   vec3d ret = v1;
   return ret *= n; 
}

inline vec3d operator *(double n, vec3d v1)
{
   vec3d ret = v1;
   return ret *= n;
}

inline vec3d operator /(vec3d v1, double n)
{ 
   vec3d ret = v1;
   return ret /= n; 
}


inline int compare(const vec3d &v1, const vec3d &v2, double eps=epsilon)
{
   if(!v1.is_set() && !v2.is_set())
      return 0;
   if(!v1.is_set())
      return -1;
   if(!v2.is_set())
      return 1;
   for(int i=0; i<3; i++)
      if(fabs(v1[i]-v2[i])>eps) {
         if(v1[i]<v2[i])
            return -1;
         if(v1[i]>v2[i])
            return 1;
      }
   return 0;
}


#endif // VEC3D_H

