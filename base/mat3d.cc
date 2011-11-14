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

/* \file mat3d.cc
 *\brief Matrix transformations for 3D geometry.
*/




#include <stdio.h>

#include "math_utils.h"
#include "mat3d.h"

/*     
bool less_than(const mat3d &m1, const mat3d &m2, double eps)
{
   for(int i=0; i<12; i++) {
      if(fabs(m1[i] - m2[i]) > eps) {
         if(m1[i] < m2[i])
            return 1;
         else
            return 0;
      }
   }
   return 0;
}
*/

int compare(const mat3d &m1, const mat3d &m2, double eps)
{
   for(int i=0; i<12; i++) {
      if(fabs(m1[i] - m2[i]) > eps) {
         if(m1[i] < m2[i])
            return 1;
         else
            return -1;
      }
   }
   return 0;
}

mat3d &mat3d::set_rot(vec3d axis, double a)
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

mat3d &mat3d::set_rot(vec3d v_from, vec3d v_to)
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

mat3d &mat3d::set_refl(vec3d norm)
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


mat3d &mat3d::set_trans_by_angles(double yz_ang, double zx_ang,
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
   

mat3d &mat3d::set_inverse()
{
   mat3d inv;
   double determ = det();
   if(fabs(determ)>epsilon) {
      // http://www.euclideanspace.com/maths/algebra/matrix/functions
      // /inverse/fourD/index.htm
      inv.m[0]  = m[6]*m[11]*m[13] - m[7]*m[10]*m[13] + m[7]*m[9]*m[14]
                - m[5]*m[11]*m[14] - m[6]*m[9]*m[15]  + m[5]*m[10]*m[15];
      inv.m[1]  = m[3]*m[10]*m[13] - m[2]*m[11]*m[13] - m[3]*m[9]*m[14]
                + m[1]*m[11]*m[14] + m[2]*m[9]*m[15]  - m[1]*m[10]*m[15];
      inv.m[2]  = m[2]*m[7]*m[13]  - m[3]*m[6]*m[13]  + m[3]*m[5]*m[14]
                - m[1]*m[7]*m[14]  - m[2]*m[5]*m[15]  + m[1]*m[6]*m[15];
      inv.m[3]  = m[3]*m[6]*m[9]   - m[2]*m[7]*m[9]   - m[3]*m[5]*m[10]
                + m[1]*m[7]*m[10]  + m[2]*m[5]*m[11]  - m[1]*m[6]*m[11];
      inv.m[4]  = m[7]*m[10]*m[12] - m[6]*m[11]*m[12] - m[7]*m[8]*m[14]
                + m[4]*m[11]*m[14] + m[6]*m[8]*m[15]  - m[4]*m[10]*m[15];
      inv.m[5]  = m[2]*m[11]*m[12] - m[3]*m[10]*m[12] + m[3]*m[8]*m[14]
                - m[0]*m[11]*m[14] - m[2]*m[8]*m[15]  + m[0]*m[10]*m[15];
      inv.m[6]  = m[3]*m[6]*m[12]  - m[2]*m[7]*m[12]  - m[3]*m[4]*m[14]
                + m[0]*m[7]*m[14]  + m[2]*m[4]*m[15]  - m[0]*m[6]*m[15];
      inv.m[7]  = m[2]*m[7]*m[8]   - m[3]*m[6]*m[8]   + m[3]*m[4]*m[10]
                - m[0]*m[7]*m[10]  - m[2]*m[4]*m[11]  + m[0]*m[6]*m[11];
      inv.m[8]  = m[5]*m[11]*m[12] - m[7]*m[9]*m[12]  + m[7]*m[8]*m[13]
                - m[4]*m[11]*m[13] - m[5]*m[8]*m[15]  + m[4]*m[9]*m[15];
      inv.m[9]  = m[3]*m[9]*m[12]  - m[1]*m[11]*m[12] - m[3]*m[8]*m[13]
                + m[0]*m[11]*m[13] + m[1]*m[8]*m[15]  - m[0]*m[9]*m[15];
      inv.m[10] = m[1]*m[7]*m[12]  - m[3]*m[5]*m[12]  + m[3]*m[4]*m[13]
                - m[0]*m[7]*m[13]  - m[1]*m[4]*m[15]  + m[0]*m[5]*m[15];
      inv.m[11] = m[3]*m[5]*m[8]   - m[1]*m[7]*m[8]   - m[3]*m[4]*m[9]
                + m[0]*m[7]*m[9]   + m[1]*m[4]*m[11]  - m[0]*m[5]*m[11];
      inv.m[12] = m[6]*m[9]*m[12]  - m[5]*m[10]*m[12] - m[6]*m[8]*m[13]
                + m[4]*m[10]*m[13] + m[5]*m[8]*m[14]  - m[4]*m[9]*m[14];
      inv.m[13] = m[1]*m[10]*m[12] - m[2]*m[9]*m[12]  + m[2]*m[8]*m[13]
                - m[0]*m[10]*m[13] - m[1]*m[8]*m[14]  + m[0]*m[9]*m[14];
      inv.m[14] = m[2]*m[5]*m[12]  - m[1]*m[6]*m[12]  - m[2]*m[4]*m[13]
                + m[0]*m[6]*m[13]  + m[1]*m[4]*m[14]  - m[0]*m[5]*m[14];
      inv.m[15] = m[1]*m[6]*m[8]   - m[2]*m[5]*m[8]   + m[2]*m[4]*m[9]
                - m[0]*m[6]*m[9]   - m[1]*m[4]*m[10]  + m[0]*m[5]*m[10];
      
      for(int i=0; i<16; i++)
         inv.m[i] /= determ;
      *this = inv;
   }
   else // det very small 
      to_zero();

   return *this;
}

mat3d &mat3d::operator *=(const mat3d &mat)
{
   mat3d new_m;
   new_m.to_zero();
   for(int i=0; i<16; i++)
      for(int j=0; j<4; j++)
         new_m[i] += m[(i/4)*4+j]*mat[(j)*4 + (i%4)];

   *this = new_m;
   return *this;
}


vec3d operator *(const mat3d &mat, const vec3d &v)
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



// http://www.j3d.org/matrix_faq/matrfaq_latest.html  Q55.
vec4d mat3d::get_quaternion() const
{
   vec4d quat; // X, Y, Z, W
   double T = 1 + m[0] + m[5] + m[10];
   double S;

   // If the trace of the matrix is equal to zero then identify
   // which major diagonal element has the greatest value.
   // Depending on this, calculate the following:
   if(T > epsilon) {
      S = sqrt(T) * 2;
      quat[0] = (m[9] - m[6])/S;
      quat[1] = (m[2] - m[8])/S;
      quat[2] = (m[4] - m[1])/S;
      quat[3] = 0.25*S;
   }
   else if(m[0]>m[5] && m[0]>m[10]) {    // Column 0: 
      S = sqrt(1.0 + m[0] - m[5] - m[10])*2;
      quat[0] = 0.25 * S;
      quat[1] = (m[4] + m[1])/S;
      quat[2] = (m[2] + m[8])/S;
      quat[3] = (m[9] - m[6])/S;
   }
   else if (m[5] > m[10]) {              // Column 1: 
      S  = sqrt(1.0 + m[5] - m[0] - m[10])*2;
      quat[0] = (m[4] + m[1])/S;
      quat[1] = 0.25*S;
      quat[2] = (m[9] + m[6])/S;
      quat[3] = (m[2] - m[8])/S;
   }
   else {                                // Column 2:
      S  = sqrt(1.0 + m[10] - m[0] - m[5])*2;
      quat[0] = (m[2] + m[8])/S;
      quat[1] = (m[9] + m[6])/S;
      quat[2] = 0.25*S;
      quat[3] = (m[4] - m[1])/S;
   }
   return quat;
}

// http://www.gamedev.net/community/forums/topic.asp?topic_id=166424
vec3d quat2euler(const vec4d &quat)
{
   double sqx = quat[0]*quat[0];
   double sqy = quat[1]*quat[1];
   double sqz = quat[2]*quat[2];
   double sqw = quat[3]*quat[3];
   
   return -vec3d(
         atan2(2*(quat[1]*quat[2] + quat[0]*quat[3]), -sqx-sqy+sqz+sqw),
         asin(safe_for_trig(-2*(quat[2]*quat[0] - quat[1]*quat[3]))),
         atan2(2*(quat[0]*quat[1] + quat[2]*quat[3]),  sqx-sqy-sqz+sqw)  );
}
 
vec3d mat3d::get_euler() const
{
   return quat2euler(get_quaternion());
}



void mat3d::dump(const char *var, FILE *strm) const
{
   fprintf(strm, "%s", var);
   for(int i=0; i<16; i++) {
      if(!(i%4))
         fprintf(strm, "\n\t");
      fprintf(strm, "\t%g", m[i]);
   }
   fprintf(strm, "\n");
}
      
mat3d &mat3d::set_alignment(vector<vec3d> from, vector<vec3d> to)
{
   const mat3d r = mat3d::rot(0.1,0,0)
                 * mat3d::rot(0,0.2,0)
                 * mat3d::rot(0,0,0.3);
   const mat3d inv_r = mat3d::rot(0,0,-0.3)
                 * mat3d::rot(0,-0.2,0)
                 * mat3d::rot(-0.1,0,0);
       
   if(to.size()>1)
      transform(to, r);

   mat3d trans = mat3d::transl(to[0] - from[0]);

   if(to.size()>1) {
      trans = mat3d::transl(to[0])
            * mat3d::rot(from[1]-from[0], to[1]-to[0]) 
            * mat3d::transl(-to[0])
            * trans;
   }
   
   if(to.size()>2) {
      for(unsigned int j=0; j<3; j++)
         from[j] = trans * from[j];
         
      vec3d norm1 = vcross(from[2]-from[0], from[1]-from[0]);
      vec3d norm2 = vcross(to[2]-to[0], to[1]-to[0]);
      // Maybe test form norm size here
      norm1.to_unit();
      norm2.to_unit();
      trans = mat3d::transl(to[0])
            * mat3d::rot(norm1, norm2)
            * mat3d::transl(-to[0])
            * trans;
   }
  
   *this = (to.size()>1) ? inv_r * trans : trans;
   return *this;
}

mat3d &mat3d::set_alignment(vec3d from1, vec3d from2, vec3d from3,
      vec3d to1, vec3d to2, vec3d to3)
{
   vector<vec3d> from(3);
   from[0] = from1;
   from[1] = from2;
   from[2] = from3;
   vector<vec3d> to(3);
   to[0] = to1;
   to[1] = to2;
   to[2] = to3;
   return set_alignment(from, to);
}


mat3d &mat3d::set_alignment(vec3d from1, vec3d from2, vec3d to1, vec3d to2)
{
   
   const mat3d r = mat3d::rot(0.1,0,0)
                 * mat3d::rot(0,0.2,0)
                 * mat3d::rot(0,0,0.3);
   const mat3d inv_r = mat3d::rot(0,0,-0.3)
                 * mat3d::rot(0,-0.2,0)
                 * mat3d::rot(-0.1,0,0);
   to1 = r * to1;
   to2 = r * to2;
   

   mat3d trans = mat3d::rot(from1, to1);
   from1 = trans * from1;
   from2 = trans * from2;

   vec3d norm1 = vcross(from2, from1);
   vec3d norm2 = vcross(to2, to1);
   // Maybe test form norm size here
   norm1.to_unit();
   norm2.to_unit();

   // find the angle to rotate abot to1, in the range -180<ang<=180
   double ang = acos(safe_for_trig(vdot(norm1, norm2)));
   if(vtriple(to1, norm1, norm2)<0)
      ang *= -1;

   //trans = mat3d::rot(to1, acos(vdot(norm1, norm2))) * trans;
   //*this = inv_r*trans;
   
   trans = mat3d::rot(to1, ang) * trans;
   *this = inv_r*trans;
   //*this = mat3d::rot(to1, ang) * trans;
   return *this;
}
 
