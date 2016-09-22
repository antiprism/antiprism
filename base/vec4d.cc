/*
   Copyright (c) 2003-2011, Adrian Rossiter
   
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

/* \file vec4d.cc
 *  \brief Vectors for 4D geometry
 *
 *  A vector class with common vector operations.
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "vec4d.h"

vec4d vec4d::X(1, 0, 0, 0);
vec4d vec4d::Y(0, 1, 0, 0);
vec4d vec4d::Z(0, 0, 1, 0);
vec4d vec4d::W(0, 0, 0, 1);
vec4d vec4d::zero(0, 0, 0, 0);


bool vec4d::read(const char *str, char *errmsg)
{
   int i;
   double f[4];
   char buff;
   char buff2;
   if( sscanf(str, " %lf , %lf , %lf , %lf %c", f, f+1, f+2, f+3, &buff) != 4 &&
         sscanf(str, " %lf %lf %lf %lf %c", f, f+1, f+2, f+3, &buff2) != 4 ) {
      if(errmsg)
         strcpy(errmsg, "didn't find four numbers");
      return false;
   }

   for(i=0; i<4; i++) {
      if(isinf(f[i])) {
         if(errmsg) {
            const char *pos[] = {"first", "second", "third", "fourth"};
            sprintf(errmsg, "%s number too large\n", pos[i]);
         }
         return false;
      }
   }

   v[0] = f[0];
   v[1] = f[1];
   v[2] = f[2];
   v[3] = f[3];

   return true;
}


void vec4d::dump(const char *var, FILE *file) const
{
   if(var)
      fprintf(file, "%s=", var);
   if(is_set())  // ihas been set
      fprintf(file, "(%f,%f,%f,%f)\n", v[0], v[1], v[2], v[3]);
   else
      fprintf(file, "(not set)\n");
}

// make vector unusable
void vec4d::unset()
{
   v[0]=NAN;
}


inline double det2(double a11, double a12, double a21, double a22)
{
   return a11*a22 - a12*a21;
}


inline double det3(double a11, double a12, double a13,
                   double a21, double a22, double a23,
                   double a31, double a32, double a33)
{
   return +a11*det2(a22, a23, a32, a33)
          -a12*det2(a21, a23, a31, a33)
          +a13*det2(a21, a22, a31, a32);
}


vec4d vcross(const vec4d &v1, const vec4d &v2, const vec4d &v3)
{
	vec4d vprod;
   vprod[0] =  det3(v1[1],v1[2],v1[3], v2[1],v2[2],v2[3], v3[1],v3[2],v3[3]);
	vprod[1] = -det3(v1[0],v1[2],v1[3], v2[0],v2[2],v2[3], v3[0],v3[2],v3[3]);
	vprod[2] =  det3(v1[0],v1[1],v1[3], v2[0],v2[1],v2[3], v3[0],v3[1],v3[3]);
	vprod[3] = -det3(v1[0],v1[1],v1[2], v2[0],v2[1],v2[2], v3[0],v3[1],v3[2]);
	return vprod;
}


