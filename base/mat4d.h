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


/*!\file mat4d.h
  \brief Matrix transformations for 4D geometry.
*/


#ifndef MAT4D_H
#define MAT4D_H

#include <math.h>

#include "math_utils.h"
#include "vec4d.h"


///Matrix for transformations in 4D
class mat4d
{
   private:
     double m[25];

   public:

     ///Constructor
     /**Initialises to the identity */
     mat4d() { to_unit(); }

     ///Constructor
     /**Initialise the rows of the rotation part.
      * \param r1 first row.
      * \param r2 second row.
      * \param r3 third row.
      * \param r4 fourth row. */
     mat4d(vec4d r1, vec4d r2, vec4d r3, vec4d r4);


     ///Read access to the element array.
     /**\param idx index into the 16 element array.
      * \return The element. */
     double operator [](int idx) const { return m[idx]; }
     
     ///Write access to the element array.
     /**\param idx index into the 16 element array.
      * \return The element. */
     double &operator [](int idx) { return m[idx]; }
     
     ///Multiply this matrix by another matrix.
     /**\param mat the matrix to multiply.
      * \return A reference to the result, held in this matrix. */
     mat4d &operator *=(const mat4d &mat);
     
     ///Multiply this matrix by a scalar.
     /**\param n the scalar.
      * \return A reference to the result, held in this matrix. */
     mat4d &operator *=(double n);

     ///Add another matrix to this matrix.
     /**\param mat the matrix to add.
      * \return A reference to the result, held in this matrix. */
     mat4d &operator +=(const mat4d &mat);
     
        
     ///Transpose the matrix
     //* \return A reference this matrix, transposed. */
     mat4d &transpose();
     
     ///Set the matrix to the unit matrix.
     //* \return A reference to this matrix. */
     mat4d &to_unit();
     
     ///Get a unit matrix.
     //* \return A unit matrix. */
     static mat4d unit();
     
     ///Set the matrix to zero
     //* \return A reference to this matrix. */
     mat4d &to_zero();
     
     ///Get a zero matrix.
     //* \return A zero matrix. */
     static mat4d zero();
     

     ///Set a rotation by plane axis and angle.
     /**Rotate around an arbitrary plane axis through the origin.
      * \param u the first basis vector of the axis plane to rotate around.
      * \param v the second basis vector of the axis plane to rotate around.
      * \param angle the angle to rotate, in radians.
      * \return  A reference to this matrix with the resulting rotation set. */
     mat4d &set_rot(vec4d u, vec4d v, double angle);

     ///Get a rotation by plane axis and angle.
     /**Rotate around an arbitrary plane axis through the origin.
      * \param u the first basis vector of the axis plane to rotate around.
      * \param v the second basis vector of the axis plane to rotate around.
      * \param angle the angle to rotate, in radians.
      * \return  A matrix with the resulting rotation set. */
     static mat4d rot(vec4d u, vec4d v, double angle);
     
     ///Set a rotation by rotating one direction vector onto another.
     /**\param v_from the direction vector to rotate.
      * \param v_to the direction \a v_from should point after the rotation.
      * \return  A reference to this matrix with the resulting rotation set. */
     mat4d &set_rot(vec4d v_from, vec4d v_to);
     
     ///Get a rotation by rotating one direction vector onto another.
     /**\param v_from the direction vector to rotate.
      * \param v_to the direction \a v_from should point after the rotation.
      * \return  A matrix with the resulting rotation set. */
     static mat4d rot(vec4d v_from, vec4d v_to);
     
     ///Set a rotation by rotating around the coordinate-plane axes.
     /**Rotate around the x-axis, y-axis and then z-axis.
      * \param xy_ang angle in radians to rotate around xy-axis.
      * \param yz_ang angle in radians to rotate around yz-axis.
      * \param zw_ang angle in radians to rotate around zw-axis.
      * \param wx_ang angle in radians to rotate around wx-axis.
      * \param xz_ang angle in radians to rotate around xz-axis.
      * \param yw_ang angle in radians to rotate around yw-axis.
      * \return  A reference to this matrix with the resulting rotation set. */
     mat4d &set_rot(double xy_ang, double yz_ang, double zw_ang,
                      double wx_ang, double xz_ang, double yw_ang);

     ///Get a rotation by rotating around the coordinate-plane axes.
     /**Rotate around the x-axis, y-axis and then z-axis.
      * \param xy_ang angle in radians to rotate around xy-axis.
      * \param yz_ang angle in radians to rotate around yz-axis.
      * \param zw_ang angle in radians to rotate around zw-axis.
      * \param wx_ang angle in radians to rotate around wx-axis.
      * \param xz_ang angle in radians to rotate around xz-axis.
      * \param yw_ang angle in radians to rotate around yw-axis.
      * \return  A matrix with the resulting rotation set. */
     static mat4d rot(double xy_ang, double yz_ang, double zw_ang,
                      double wx_ang, double xz_ang, double yw_ang);
     
     ///Set a translation.
     /**\param trans the vector to translate by.
      * \return  A reference to this matrix with the translation set. */
     mat4d &set_transl(vec4d trans);
     
     ///Get a translation.
     /**\param trans the vector to translate by.
      * \return  A matrix with the translation set. */
     static mat4d transl(vec4d trans);

     
     ///Set a uniform scaling.
     /**Scale uniformly in all directions, away from the origin.
      * \param scal the amount to scale.
      * \return  A reference to this matrix with the scaling set. */
     mat4d &set_scale(double scal);

     ///Get a uniform scaling.
     /**Scale uniformly in all directions, away from the origin.
      * \param scal the amount to scale.
      * \return  A matrix with the scaling set. */
     static mat4d scale(double scal);
     
     ///Set a scaling by components.
     /**Scale in the direction of the coordinate axes, away from the origin.
      * \param x_scal the amount to scale in the direction of the x-axis.
      * \param y_scal the amount to scale in the direction of the y-axis.
      * \param z_scal the amount to scale in the direction of the z-axis.
      * \param w_scal the amount to scale in the direction of the z-axis.
      * \return  A reference to this matrix with the scaling set. */
     mat4d &set_scale(double x_scal, double y_scal, double z_scal,
           double w_scal);
     
     ///Get a scaling by components.
     /**Scale in the direction of the coordinate axes, away from the origin.
      * \param x_scal the amount to scale in the direction of the x-axis.
      * \param y_scal the amount to scale in the direction of the y-axis.
      * \param z_scal the amount to scale in the direction of the z-axis.
      * \param w_scal the amount to scale in the direction of the z-axis.
      * \return  A reference to this matrix with the scaling set. */
     static mat4d scale(double x_scal, double y_scal, double z_scal,
           double w_scal);
 
     ///Set an inversion.
     /**Inversion through the origin.
      * \return  A reference to this matrix with the inversion set. */
     mat4d &set_inversion();

     ///Get an inversion.
     /**Inversion through the origin.
      * \return  A matrix with the inversion set. */
     static mat4d inversion();
     
      ///Debugging print of a matrix
     /**\param var a string to identify the matrix variable.
      * \param file file stream to print the variable. */
     void dump(const char *var="", FILE *file=stderr);         // print matrix
};


///Transform a column vector.
/**\param mat the transformation matrix.
 * \param v the column vector.
 * \return The transformed vector (left-multiplied by the matrix). */
vec4d operator *(const mat4d &mat, const vec4d &v);

///Transform a row vector.
/**\param v the column vector.
 * \param mat the transformation matrix.
 * \return The transformed vector (right-multiplied by the matrix). */
vec4d operator *(const vec4d &v, const mat4d &mat);

///Multiply two matrices.
/**\param m1 the first matrix.
 * \param m2 the second matrix.
 * \return The result of the first matrix multiplying the second. */
mat4d operator *(const mat4d &m1, const mat4d &m2);

///Multiply a matrix by a scalar.
/**\param n the scalar.
 * \param mat the matrix.
 * \return the result of multiplying the matrix by the scalar. */
mat4d operator *(double n, const mat4d &mat);

///Multiply a matrix by a scalar.
/**\param mat the matrix.
 * \param n the scalar.
 * \return the result of multiplying the matrix by the scalar. */
mat4d operator *(const mat4d &mat, double n);

///Add two matrices.
/**\param m1 the first matrix.
 * \param m2 the second matrix.
 * \return The result of adding the two matrices. */
mat4d operator +(const mat4d &m1, const mat4d &m2);


// inline functions
inline mat4d::mat4d(vec4d r1, vec4d r2, vec4d r3, vec4d r4)
{
   to_unit();
   vec4d *pvec[] = {&r1, &r2, &r3, &r4};
   for(int i=0; i<20; i++)
      if((i%5) != 4)
         m[i] = (*pvec[i/5])[i%5];
}

inline mat4d &mat4d::to_unit()
{
   for(int i=0; i<25; i++)
      m[i]=((i%5)==(i/5))?1:0;
   return *this;
}

inline mat4d mat4d::unit()
{ 
   mat4d mat;
   return mat.to_unit();
}

inline mat4d &mat4d::to_zero()
{
   for(int i=0; i<25; i++)
      m[i]=0;
   return *this;
}

inline mat4d mat4d::zero()
{ 
   mat4d mat;
   return mat.to_zero();
}


inline mat4d &mat4d::set_rot(vec4d u, vec4d v, double angle)
{
   u.unit();
   v.unit();
   mat4d A(vec4d(0, +(u[2]*v[3]-u[3]*v[2]),
              -(u[1]*v[3]-u[3]*v[1]), +(u[1]*v[2]-u[2]*v[1])),
           vec4d(-(u[2]*v[3]-u[3]*v[2]), 0,
              +(u[0]*v[3]-u[3]*v[0]), -(u[0]*v[2]-u[2]*v[0])),
           vec4d(+(u[1]*v[3]-u[3]*v[1]), -(u[0]*v[3]-u[3]*v[0]),
              0, +(u[0]*v[1]-u[1]*v[0])),
           vec4d(-(u[1]*v[2]-u[2]*v[1]), +(u[0]*v[2]-u[2]*v[0]),
              -(u[0]*v[1]-u[1]*v[0]), 0)  );
   *this = (unit() + (sin(angle)*A) + ((1-cos(angle))*(A*A))); 
   m[24] = 1;
   return *this;
}

inline mat4d mat4d::rot(vec4d u, vec4d v, double angle)
{ 
   mat4d mat;
   return mat.set_rot(u, v, angle);
}


inline mat4d &mat4d::set_rot(vec4d v_from, vec4d v_to)
{
   v_from.unit();
   v_to.unit();
   vec4d n2, n3, orth;
   rand_gen rnd;
   rnd.time_seed();
   n2 = vcross(vec4d::random(rnd).unit(), v_from, v_to);
   n3 = vcross(n2, v_from, v_to);
   orth = vcross(v_from, n2, n3);
   set_rot(n2, n3, acos(safe_for_trig(vdot(v_from, v_to))) );
   return *this;
}

inline mat4d mat4d::rot(vec4d v_from, vec4d v_to)
{
   mat4d mat;
   return mat.set_rot(v_from, v_to);
}



inline mat4d &mat4d::set_rot(double xy_ang, double yz_ang, double zw_ang,
                      double wx_ang, double xz_ang, double yw_ang)
{
   int planes[]  = {2,3,  3,1,  0,1,  1,2,  1,3,  2,0 };
   double angs[] = {xy_ang, yz_ang, zw_ang, wx_ang, xz_ang, yw_ang};
   to_unit();
   for(int p=0; p<6; p++) {
      vec4d p1(0,0,0,0);
      p1[planes[2*p]] = 1;
      vec4d p2(0,0,0,0);
      p2[planes[2*p+1]] = 1;
      *this = rot(p1, p2, angs[p]) * (*this);
   }
   return *this;
}

inline mat4d mat4d::rot(double xy_ang, double yz_ang, double zw_ang,
                      double wx_ang, double xz_ang, double yw_ang)
{
   mat4d mat;
   return mat.set_rot(xy_ang, yz_ang, zw_ang, wx_ang, xz_ang, yw_ang);
}
 
inline mat4d &mat4d::set_transl(vec4d trans)
{
   to_unit();
   m[4] = trans[0];
   m[9] = trans[1];
   m[14]= trans[2];
   m[19]= trans[3];
   return *this;
}
     
inline mat4d mat4d::transl(vec4d trans)
{
   mat4d mat;
   return mat.set_transl(trans);
}



inline mat4d &mat4d::set_scale(double scal)
{
   return set_scale(scal, scal, scal, scal);
}
     
inline mat4d mat4d::scale(double scal)
{
   mat4d mat;
   return mat.set_scale(scal);
}

inline mat4d &mat4d::set_scale(double x_scal, double y_scal, double z_scal,
      double w_scal)
{
   m[0] = x_scal;
   m[6] = y_scal;
   m[12] = z_scal;
   m[18] = w_scal;
   return *this;
}
     
inline mat4d mat4d::scale(double x_scal, double y_scal, double z_scal, 
      double w_scal)
{
   mat4d mat;
   return mat.set_scale(x_scal, y_scal, z_scal, w_scal);
}


inline mat4d &mat4d::transpose()
{
   for(int i=0; i<25; i++)
      if(i/5 < i%5) {
         double tmp = m[i];
         m[i] = m[(i%5)*5 + (i/5)];
         m[(i%5)*5 + (i/5)] = tmp;
      }

   return *this;
}


inline mat4d &mat4d::set_inversion()
{ 
   return set_scale(-1);
}   

inline mat4d mat4d::inversion()
{ 
   mat4d mat;
   return mat.set_inversion();
}   




inline mat4d &mat4d::operator *=(const mat4d &mat)
{
   mat4d new_m;
   new_m.to_zero();
   for(int i=0; i<25; i++)
      for(int j=0; j<5; j++)
         new_m[i] += m[(i/5)*5+j]*mat[(j)*5 + (i%5)];

   *this = new_m;
   return *this;
}

inline mat4d &mat4d::operator *=(double n)
{
	for(int i=0; i<25; i++)
      m[i] *= n;
	return *this;
}

inline mat4d &mat4d::operator +=(const mat4d &mat)
{
   for(int i=0; i<25; i++)
      m[i] += mat[i];
   return *this;
}


inline vec4d operator *(const mat4d &mat, const vec4d &v)
{
   vec4d new_v(0, 0, 0, 0);
   for(int i=0; i<20; i++) {
      if(i%5 != 4)
         new_v[i/5] += mat[i]*v[i%5];
      else
         new_v[i/5] += mat[i];
   }
   return new_v;
}

inline vec4d operator *(const vec4d &v, const mat4d &m)
{
   mat4d m_ret = m;
   m_ret.transpose();
   return m_ret*v;
}


inline mat4d operator *(const mat4d &m1, const mat4d &m2)
{
   mat4d m_ret = m1;
   return m_ret*=m2;
}

inline mat4d operator *(double n, const mat4d &mat)
{
   mat4d m_ret = mat;
   return m_ret *= n;
}

inline mat4d operator *(const mat4d &mat, double n)
{
   mat4d m_ret = mat;
   return m_ret *= n;
}

inline mat4d operator +(const mat4d &m1, const mat4d &m2)
{
   mat4d m_ret = m1;
   return m_ret+=m2;
}


#endif // MAT4D_H


