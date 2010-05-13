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

/* \file vec_utils_norm.cc
   \brief Vector utilities, face normals
*/


#include "vec_utils.h"

// Get a normal to the plane containing three points.
/* \param Q0 a point on the plane.
 * \param Q1 a second point on the plane.
 * \param Q2 a third points on the plane, not on the line Q0Q1.
 * \return the normal. */
inline vec3d find_norm(vec3d Q0, vec3d Q1, vec3d Q2)
{ return vcross((Q0-Q1).unit(), (-Q1+Q2).unit()); }



vec3d face_norm_largest(const vector<vec3d> &verts, const vector<int> &face)
{
   unsigned int sz = face.size();
   vec3d norm = vec3d(0,0,0);
   vec3d n = vec3d(0,0,0);
   double max = 0;
   for(unsigned int i=0; i<sz; i++) {
      n = find_norm(verts[face[i]], verts[face[(i+1)%sz]],
                    verts[face[(i+2)%sz]]);
      double mag2 = n.mag2();
      if(mag2 > max) {
         max = mag2;
         norm = n;
      }
   }
   return norm;
}

// adapted from http://jgt.akpeters.com/papers/Sunday02/
double findArea(const vector<vec3d> &verts, const vector<int> &face, int idx0, int idx1)
{
   int sz = face.size();
   double sum = 0.0;
   for (int i=1; i <= sz; i++)
      sum += verts[face[i%sz]][idx0]*
         (verts[face[(i+1)%sz]][idx1]-verts[face[(i-1+sz)%sz]][idx1]);
   return (sum / 2.0);
}

vec3d face_norm(const vector<vec3d> &verts, const vector<int> &face, bool allow_zero)
{
   // Newell normal
   vec3d norm(findArea(verts, face, 1, 2),
              findArea(verts, face, 2, 0),
              findArea(verts, face, 0, 1));
   return (allow_zero || norm.mag()>1e-8)? norm : face_norm_largest(verts, face);
}


