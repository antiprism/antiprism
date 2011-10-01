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

/*!\file mat3d.h
 *\brief Matrix transformations for 3D geometry.
*/


#ifndef MAT3D_H
#define MAT3D_H

#include <math.h>
#include <vector>
#include "vec3d.h"
#include "vec4d.h"

using std::vector;


///Matrix for transformations in 3D
class mat3d
{
   private:
      static double det2(const double &a11, const double &a12,
                         const double &a21, const double &a22)
         { return a11*a22 - a12*a21; }
     
      double m[16];

   public:

     ///Constructor
     /**Initialises to the identity. */
     mat3d() { to_unit(); }

     ///Constructor
     /**Initialise the rows of the rotation part.
      * \param r1 first row.
      * \param r2 second row.
      * \param r3 third row. */
     mat3d(const vec3d &r1, const vec3d &r2, const vec3d &r3);

     ///Constructor
     /**Initialise the rows.
      * \param r1 first row.
      * \param r2 second row.
      * \param r3 third row.
      * \param r4 fourth row. */
     mat3d(const vec4d &r1, const vec4d &r2, const vec4d &r3,
           vec4d r4=vec4d(0,0,0,1));

     ///Set the columns.
     /**Initialise the top three rows by columns
      * \param c1 first column.
      * \param c2 second column.
      * \param c3 third column.
      * \param c4 fourth column. */
     mat3d &set_cols(const vec3d &c1, const vec3d &c2, const vec3d &c3,
           const vec3d &c4=vec3d(0,0,0));

     ///Read access to the element array.
     /**\param idx index into the 16 element array.
      * \return The element. */
     double operator [](int idx) const { return m[idx]; }
     
     ///Write access to the element array.
     /**\param idx index into the 16 element array.
      * \return The element. */
     double &operator [](int idx) { return m[idx]; }
    
     ///Read access to the start of the element array.
     /**\return The start of the array. */
     const double *get_m() {return m;}
    
     ///Multiply this matrix by another matrix.
     /**\param mat the matrix to multiply.
      * \return A reference to the result, held in this matrix. */
     mat3d &operator *=(const mat3d &mat);

     ///Transpose the matrix
     //* \return A reference this matrix, transposed. */
     mat3d &transpose();
     
     ///Set to the inverse of the matrix.
     /**\return  A reference to this matrix set to its inverse. */
     mat3d &set_inverse();

     ///Get the inverse of a matrix.
     /**\param mat a matrix to get the inverse of.
      * \return  The inverse of \p mat. */
     static mat3d inverse(const mat3d &mat);

     ///Determinant
     /**\return The determinant of the matrix. */
     double det() const;
     
     ///Set the matrix to the unit matrix.
     //* \return A reference to this matrix. */
     mat3d &to_unit();
     
     ///Get a unit matrix.
     //* \return A unit matrix. */
     static mat3d unit();
     
     ///Set the matrix to zero
     //* \return A reference to this matrix. */
     mat3d &to_zero();
     
     ///Get a zero matrix.
     //* \return A zero matrix. */
     static mat3d zero();
     
     ///Set a rotation by rotating around the coordinate axes.
     /**Rotate around the x-axis, y-axis and then z-axis.
      * \param x_ang angle in radians to rotate around x-axis.
      * \param y_ang angle in radians to rotate around y-axis.
      * \param z_ang angle in radians to rotate around z-axis.
      * \return  A reference to this matrix with the resulting rotation set. */
     mat3d &set_rot(double x_ang, double y_ang, double z_ang);
     
     ///Get a rotation by rotating around the coordinate axes.
     /**Rotate around the x-axis, y-axis and then z-axis.
      * \param x_ang angle in radians to rotate around x-axis.
      * \param y_ang angle in radians to rotate around y-axis.
      * \param z_ang angle in radians to rotate around z-axis.
      * \return  A matrix with the resulting rotation set. */
     static mat3d rot(double x_ang, double y_ang, double z_ang);

     ///Set a rotation by rotating around the coordinate axes.
     /**Rotate around the x-axis, y-axis then z-axis.
      * \param angles a vector whose components are the the angles in radians
      * to rotate around he x-axis, y-axis and z-axis.
      * \return  A reference to this matrix with the resulting rotation set. */
     mat3d &set_rot(vec3d angles);

     ///Get a rotation by rotating around the coordinate axes.
     /**Rotate around the x-axis, y-axis then z-axis.
      * \param angles a vector whose components are the the angles in radians
      * to rotate around he x-axis, y-axis and z-axis.
      * \return  A matrix with the resulting rotation set. */
     static mat3d rot(vec3d angles);
     
     ///Set a rotation by axis and angle.
     /**Rotate around an arbitrary axis through the origin.
      * \param axis the axis o rotate around.
      * \param angle the angle to rotate, in radians.
      * \return  A reference to this matrix with the resulting rotation set. */
     mat3d &set_rot(vec3d axis, double angle);
     
     ///Get a rotation by axis and angle.
     /**Rotate around an arbitrary axis through the origin.
      * \param axis the axis o rotate around.
      * \param angle the angle to rotate, in radians.
      * \return  A matrix with the resulting rotation set. */
     static mat3d rot(vec3d axis, double angle);
     
     ///Set a rotation by rotating one direction vector onto another.
     /**\param v_from the direction vector to rotate.
      * \param v_to the direction \a v_from should point after the rotation.
      * \return  A reference to this matrix with the resulting rotation set. */
     mat3d &set_rot(vec3d v_from, vec3d v_to);
     
     ///Get a rotation by rotating one direction vector onto another.
     /**\param v_from the direction vector to rotate.
      * \param v_to the direction \a v_from should point after the rotation.
      * \return  A reference to this matrix with the resulting rotation set. */
     static mat3d rot(vec3d v_from, vec3d v_to);
     
     ///Set a translation.
     /**\param trans the vector to translate by.
      * \return  A reference to this matrix with the translation set. */
     mat3d &set_transl(vec3d trans);
     
     ///Get a translation.
     /**\param trans the vector to translate by.
      * \return  A matrix with the translation set. */
     static mat3d transl(vec3d trans);
     
     ///Set a reflection.
     /**Reflect in a mirror plane passing through the origin.
      * \param norm a normal to the mirror plane.
      * \return  A reference to this matrix with the reflection set. */
     mat3d &set_refl(vec3d norm);
     
     ///Get a reflection.
     /**Reflect in a mirror plane passing through the origin.
      * \param norm a normal to the mirror plane.
      * \return  A matrix with the reflection set. */
     static mat3d refl(vec3d norm);

     ///Set a uniform scaling.
     /**Scale uniformly in all directions, away from the origin.
      * \param scal the amount to scale.
      * \return  A reference to this matrix with the scaling set. */
     mat3d &set_scale(double scal);

     ///Get a uniform scaling.
     /**Scale uniformly in all directions, away from the origin.
      * \param scal the amount to scale.
      * \return  A matrix with the scaling set. */
     static mat3d scale(double scal);
     
     ///Set a scaling by components.
     /**Scale in the direction of the coordinate axes, away from the origin.
      * \param x_scal the amount to scale in the direction of the x-axis.
      * \param y_scal the amount to scale in the direction of the y-axis.
      * \param z_scal the amount to scale in the direction of the z-axis.
      * \return  A reference to this matrix with the scaling set. */
     mat3d &set_scale(double x_scal, double y_scal, double z_scal);
     
     ///Get a scaling by components.
     /**Scale in the direction of the coordinate axes, away from the origin.
      * \param x_scal the amount to scale in the direction of the x-axis.
      * \param y_scal the amount to scale in the direction of the y-axis.
      * \param z_scal the amount to scale in the direction of the z-axis.
      * \return  A reference to this matrix with the scaling set. */
     static mat3d scale(double x_scal, double y_scal, double z_scal);
     
     ///Set a scaling in a single direction.
     /**Scale in a single arbitrary direction, away from the origin.
      * \param dir the direction to scale.
      * \param scal the amount to scale.
      * \return  A reference to this matrix with the scaling set. */
     mat3d &set_scale(vec3d dir, double scal);
     
     ///Get a scaling in a single direction.
     /**Scale in a single arbitrary direction, away from the origin.
      * \param dir the direction to scale.
      * \param scal the amount to scale.
      * \return  A matrix with the scaling set. */
     static mat3d scale(vec3d dir, double scal);

     ///Set an inversion.
     /**Inversion through the origin.
      * \return  A reference to this matrix with the inversion set. */
     mat3d &set_inversion();

     ///Get an inversion.
     /**Inversion through the origin.
      * \return  A matrix with the inversion set. */
     static mat3d inversion();

     ///Set a transformation that aligns two pairs of vectors
     mat3d &set_alignment(vec3d from1, vec3d from2, vec3d to1, vec3d to2);

     ///Get a transformation that aligns two pairs of vectors
     static mat3d alignment(vec3d from1, vec3d from2, vec3d to1, vec3d to2);

     ///Set a transformation that aligns two sets of three points
     mat3d &set_alignment(vec3d from1, vec3d from2, vec3d from3,
           vec3d to1, vec3d to2, vec3d to3);

     ///Get a transformation that aligns two sets of three points
     static mat3d alignment(vec3d from1, vec3d from2, vec3d from3,
           vec3d to1, vec3d to2, vec3d to3);

     ///Set a transformation that aligns two sets of one, two or three points
     mat3d &set_alignment(vector<vec3d> from, vector<vec3d> to);

     ///Get a transformation that aligns two sets of one, two or three points
     static mat3d alignment(vector<vec3d> from, vector<vec3d> to);


     ///Set a transformation that makes particular angles between mapped axes.
     /**The x-axis is left unchanged. The y-axis is mapped to the y'-axis
      * in the xy-plane. The z-axis is mapped to the z'-axis.
      * \param yz_ang the angle between the y'-axis and the z'axis.
      * \param zx_ang the angle between the z'-axis and the x-axis.
      * \param xy_ang the angle between the x-axis and the y'-axis.
      * \param valid to return whether a transformation based on
      * the angles is valid.
      * \return  A reference to this matrix with the transformation set. */
     mat3d &set_trans_by_angles(double yz_ang, double zx_ang, double xy_ang,
           bool *valid=0);
     
     ///Get a transformation that makes particular angles between mapped axes.
     /**The x-axis is left unchanged. The y-axis is mapped to the y'-axis
      * in the xy-plane. The z-axis is mapped to the z'-axis.
      * \param zx_ang the angle between the z'-axis and the x-axis.
      * \param xy_ang the angle between the x-axis and the y'-axis.
      * \param yz_ang the angle between the y'-axis and the z'axis.
      * \param valid to return whether a transformation based on
      * the angles is valid.
      * \return  A matrix with the transformation set. */
     static mat3d trans_by_angles(double yz_ang, double zx_ang, double xy_ang,
           bool *valid=0);
     
     ///Get the translation vector.
     /* \return the translation component of the matrix. */
     vec3d get_transl() const;

     ///Get quaternion
     /**Convert the rotation part of the matrix into a quaterinon.
      * \return A vector representing a quaternion
      * (components X, Y, Z, and W.) */
     vec4d get_quaternion() const;

     ///Get Euler angles.
     /**\return A vector representing rotations about the X, Y and Z
      * axis, in that order. */
     vec3d get_euler() const;


     ///Debugging print of a matrix
     /**\param var a string to identify the matrix variable.
      * \param file file stream to print the variable. */
     void dump(const char *var="", FILE *file=stderr) const;
};


///Transform a column vector.
/**\param mat the transformation matrix.
 * \param v the column vector.
 * \return The transformed vector (left-multiplied by the matrix). */
vec3d operator *(const mat3d &mat, const vec3d &v);
     
///Transform a row vector.
/**\param v the column vector.
 * \param mat the transformation matrix.
 * \return The transformed vector (right-multiplied by the matrix). */
vec3d operator *(const vec3d &v, const mat3d &mat);

///Multiply two matrices.
/**\param m1 the first matrix.
 * \param m2 the second matrix.
 * \return The result of the first matrix multiplying the second. */
mat3d operator *(const mat3d &m1, const mat3d &m2);

int compare(const mat3d &m1, const mat3d &m2, double eps=epsilon);
bool operator <(const mat3d &m1, const mat3d &m2);

///Transform a set of vectors
/**\param vecs the (column) vectors to transform.
 * \param mat the matrix transformation. */
void transform(vector<vec3d> &vecs, const mat3d &mat);

// inline functions
inline mat3d::mat3d(const vec3d &r1, const vec3d &r2, const vec3d &r3)
{
   to_unit();
   const vec3d *pvec[] = {&r1, &r2, &r3};
   for(int i=0; i<12; i++)
      if((i%4) != 3)
         m[i] = (*pvec[i/4])[i%4];
}
     
inline mat3d::mat3d(const vec4d &r1, const vec4d &r2, const vec4d &r3, vec4d r4)
{
   const vec4d *pvec[] = {&r1, &r2, &r3, &r4};
   for(int i=0; i<4; i++)
      for(int j=0; j<4; j++)
         m[i] = (*pvec[i])[j];
}
 
inline mat3d &mat3d::to_unit()
{
   for(int i=0; i<16; i++)
      m[i]=((i%4)==(i/4))?1:0;
   return *this;
}

inline mat3d mat3d::unit()
{ 
   mat3d mat;
   return mat.to_unit();
}

inline mat3d &mat3d::to_zero()
{ 
   for(int i=0; i<16; i++)
      m[i]=0;
   return *this;
}

inline mat3d mat3d::zero()
{ 
   mat3d mat;
   return mat.to_zero();
}

inline mat3d &mat3d::set_transl(vec3d trans)
{
   to_unit();
   m[3] = trans[0];
   m[7] = trans[1];
   m[11] = trans[2];
   return *this;
}
     
inline mat3d mat3d::transl(vec3d trans)
{
   mat3d mat;
   return mat.set_transl(trans);
}


inline mat3d &mat3d::set_rot(vec3d axis, double a)
{
   to_zero();
   m[15] = 1;
   axis.to_unit();
   double c = cos(a);
   double s = sin(a);
   double t = 1.0 - c;
   m[0] = c + axis[0]*axis[0]*t;
   m[5] = c + axis[1]*axis[1]*t;
   m[10] = c + axis[2]*axis[2]*t;


   double tmp1 = axis[0]*axis[1]*t;
   double tmp2 = axis[2]*s;
   m[4] = tmp1 + tmp2;
   m[1] = tmp1 - tmp2;

   tmp1 = axis[0]*axis[2]*t;
   tmp2 = axis[1]*s;
   m[8] = tmp1 - tmp2;
   m[2] = tmp1 + tmp2;

   tmp1 = axis[1]*axis[2]*t;
   tmp2 = axis[0]*s;
   m[9] = tmp1 + tmp2;
   m[6] = tmp1 - tmp2;
   return *this;
}

inline mat3d mat3d::rot(vec3d axis, double a)
{ 
   mat3d mat;
   return mat.set_rot(axis, a);
}


inline mat3d &mat3d::set_rot(double x_ang, double y_ang, double z_ang)
{
	return set_rot(vec3d(0,0,1), z_ang) *=
      rot(vec3d(0,1,0), y_ang) * rot(vec3d(1,0,0), x_ang);
}	
     
inline mat3d mat3d::rot(double x_ang, double y_ang, double z_ang)
{
   mat3d mat;
   return mat.set_rot(x_ang, y_ang, z_ang);
}
     
inline mat3d &mat3d::set_rot(vec3d angles)
{
   return set_rot(angles[0], angles[1], angles[2]);
}

inline mat3d mat3d::rot(vec3d angles)
{ 
   mat3d mat;
   return mat.set_rot(angles);
}

inline mat3d &mat3d::set_cols(const vec3d &c1, const vec3d &c2, const vec3d &c3,
           const vec3d &c4)
{
   const vec3d *cols[] = {&c1, &c2, &c3, &c4};
   for(int i=0; i<4; i++)
      for(int j=0; j<3; j++)
         m[j*4+i] = (*cols[i])[j];
   return *this;
}

inline mat3d &mat3d::set_rot(vec3d v_from, vec3d v_to)
{
   v_from.to_unit();
   v_to.to_unit();
   vec3d axis = vcross(v_from, v_to);
   double cos_a = vdot(v_from, v_to);
   if(fabs(cos_a) >= 1-epsilon) {
      cos_a = cos_a>0 ? 1 : -1;
      axis = vcross(v_from, vec3d(1.2135,2.09865,3.23784));  // fix this
   }
      
   return set_rot(axis, acos(cos_a));
}

inline mat3d mat3d::rot(vec3d v_from, vec3d v_to)
{
   mat3d mat;
   return mat.set_rot(v_from, v_to);
}
         
   
inline mat3d &mat3d::set_refl(vec3d norm)
{
   norm.to_unit();
   double x = norm[0];
   double y = norm[1];
   double z = norm[2];
   m[0] = -x*x + y*y + z*z;
   m[1] = -2*x*y;
   m[2] = -2*x*z;
   m[4] = -2*y*x;
   m[5] = x*x - y*y + z*z;
   m[6] = -2*y*z;
   m[8] = -2*z*x;
   m[9] = -2*z*y;
   m[10]= x*x + y*y - z*z;
   return *this;
}

inline mat3d mat3d::refl(vec3d norm)
{
   mat3d mat;
   return mat.set_refl(norm);
}

inline mat3d &mat3d::set_scale(double scal)
{
   return set_scale(scal, scal, scal);
}
     
inline mat3d mat3d::scale(double scal)
{
   mat3d mat;
   return mat.set_scale(scal);
}

inline mat3d &mat3d::set_scale(double x_scal, double y_scal, double z_scal)
{
   m[0] = x_scal;
   m[5] = y_scal;
   m[10] = z_scal;
   return *this;
}
     
inline mat3d mat3d::scale(double x_scal, double y_scal, double z_scal)
{
   mat3d mat;
   return mat.set_scale(x_scal, y_scal, z_scal);
}

inline mat3d &mat3d::set_scale(vec3d v, double scal)
{
   vec3d x(1,0,0);
   return set_rot(x, v) *= scale(scal, 1, 1) * rot(v, x);
}
     
inline mat3d mat3d::scale(vec3d v, double scal)
{
   mat3d mat;
   return mat.set_scale(v, scal);
}

inline mat3d &mat3d::set_inversion()
{ 
   return set_scale(-1);
}   

inline mat3d mat3d::inversion()
{ 
   mat3d mat;
   return mat.set_inversion();
}   


inline mat3d mat3d::alignment(vec3d from1, vec3d from2, vec3d to1, vec3d to2)
{ 
   mat3d mat;
   return mat.set_alignment(from1, from2, to1, to2);
}


inline mat3d mat3d::alignment(vec3d from1, vec3d from2, vec3d from3,
           vec3d to1, vec3d to2, vec3d to3)
{ 
   mat3d mat;
   return mat.set_alignment(from1, from2, from3, to1, to2, to3);
}


inline mat3d mat3d::alignment(vector<vec3d> from, vector<vec3d> to)
{ 
   mat3d mat;
   return mat.set_alignment(from, to);
}


inline mat3d &mat3d::set_trans_by_angles(double yz_ang, double zx_ang,
      double xy_ang, bool *valid)
{
   double sq = 1-cos(xy_ang)*cos(xy_ang)-cos(yz_ang)*cos(yz_ang) -
      cos(zx_ang)*cos(zx_ang) + 2*cos(xy_ang)*cos(yz_ang)*cos(zx_ang);
   if(sq<-epsilon) {
      if(valid)
         *valid = false;
      to_zero();
      return *this;
   }
   else if(sq<0)
      sq = 0;

   if(valid)
      *valid = true;
   vec3d new_x(1, 0, 0);
   vec3d new_y(cos(xy_ang), sin(xy_ang), 0);
   vec3d new_z(cos(zx_ang), -cos(zx_ang)/tan(xy_ang) + cos(yz_ang)/sin(xy_ang),
               sqrt(sq)/sin(xy_ang));
   set_cols(new_x, new_y, new_z);
   return *this;
}
   



inline mat3d mat3d::trans_by_angles(double yz_ang, double zx_ang, double xy_ang,
      bool *valid)
{
   mat3d mat;
   return mat.set_trans_by_angles(yz_ang, zx_ang, xy_ang, valid);
}


inline mat3d &mat3d::transpose()
{
   for(int i=0; i<16; i++)
      if(i/4 < i%4) {
         double tmp = m[i];
         m[i] =  m[(i%4)*4 + (i/4)];
         m[(i%4)*4 + (i/4)] = tmp;
      }
   return *this;
}


inline mat3d mat3d::inverse(const mat3d &mat)
{ 
   mat3d inv = mat;
   return inv.set_inverse();
}   


inline double mat3d::det() const
{
   return +m[0]*det2(m[5], m[6], m[9], m[10])
          -m[1]*det2(m[4], m[6], m[8], m[10])
          +m[2]*det2(m[4], m[5], m[8], m[9]);
}

inline vec3d mat3d::get_transl() const
{
   return vec3d(m[3], m[7], m[11]);
}


inline mat3d &mat3d::operator *=(const mat3d &mat)
{
   mat3d new_m;
   new_m.to_zero();
   for(int i=0; i<16; i++)
      for(int j=0; j<4; j++)
         new_m[i] += m[(i/4)*4+j]*mat[(j)*4 + (i%4)];

   *this = new_m;
   return *this;
}


inline vec3d operator *(const mat3d &mat, const vec3d &v)
{
   vec3d new_v(0, 0, 0);
   for(int i=0; i<12; i++) {
      if(i%4 != 3)
         new_v[i/4] += mat[i]*v[i%4];
      else
         new_v[i/4] += mat[i];
   }
   return new_v;
}

inline vec3d operator *(const vec3d &v, const mat3d &m)
{
   mat3d m_ret = m;
   m_ret.transpose();
   return m_ret*v;
}

inline mat3d operator *(const mat3d &m1, const mat3d &m2)
{
   mat3d m_ret = m1;
   return m_ret*=m2;
}

inline bool operator <(const mat3d &m1, const mat3d &m2)
{
   return compare(m1, m2, epsilon)==-1;
}

inline bool operator ==(const mat3d &m1, const mat3d &m2)
{
   for(int i=0; i<16; i++)
      if(!double_eq(m1[i], m2[i], epsilon))
         return false;
   return true;
}

inline void transform(vector<vec3d> &verts, const mat3d &trans)
{
   for(unsigned int i=0; i<verts.size(); i++)
      verts[i] = trans * verts[i];
}

#endif // MAT3D_H


