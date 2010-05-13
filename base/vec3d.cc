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

/* \file vec3d.cc
 *  \brief Vectors for 3D geometry
 *
 *  A vector class with common vector operations.
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "vec3d.h"

#ifndef NAN
#define NAN (0.0F/0.0F)
#endif

vec3d vec3d::x(1, 0, 0);
vec3d vec3d::y(0, 1, 0);
vec3d vec3d::z(0, 0, 1);
vec3d vec3d::zero(0, 0, 0);

bool vec3d::read(char *str, char *errmsg)
{
   int i;
   double f[3];
   char buff;
   char buff2;
   if( sscanf(str, " %lf , %lf , %lf %c", f, f+1, f+2, &buff) != 3 &&
       sscanf(str, " %lf %lf %lf %c", f, f+1, f+2, &buff2) != 3 ) {
      if(errmsg)
         strcpy(errmsg, "didn't find three numbers");
      return false;
   }

   for(i=0; i<3; i++) {
      if(isinf(f[i])) {
         if(errmsg) {
            const char *pos[] = {"first", "second", "third"};
            sprintf(errmsg, "%s number too large\n", pos[i]);
         }
         return false;
      }
   }

   v[0] = f[0];
   v[1] = f[1];
   v[2] = f[2];

   return true;
}

void vec3d::dump(const char *var, FILE *file) const
{
   if(var)
      fprintf(file, "%s=", var);
   if(is_set())  // ihas been set
      fprintf(file, "(%f,%f,%f)\n", v[0], v[1], v[2]);
   else
      fprintf(file, "(not set)\n");
}

// make vector unusable
void vec3d::unset()
{
   v[0]=NAN;
}


