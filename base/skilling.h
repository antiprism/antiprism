/*
   Copyright (c) 2009, Roger Kaufman

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

/*!\file skilling.h
 * \Uniform Compounds catalogued by John Skilling 
 */

#ifndef SKILLING_H
#define SKILLING_H

#include <stdio.h>
#include <string>

#include "geom.h"

using std::string;


double v_ang_at_ax(const vec3d &v0, const vec3d &v1, const vec3d &ax);
int make_resource_uniform_compound(geom_if &geom, string pname, char *errmsg=0);


struct UCItem {
   int uc_case;
   const char *constituent;
   const char *sym_from;
   const char *sym_to;
   double angle;
   const char *name;
   const char *description;
};

class uc_poly
{
   private:
      int last_uc;
      UCItem* uc_items;

   public:     
      uc_poly();
      int get_poly(col_geom_v &geom, int sym, double angle, int n, int d, int k);
      void list_poly(int idx, FILE *fp=stderr);
      void list_polys(FILE *fp=stderr);
      int lookup_sym_no(string sym);
      int get_last_uc() { return last_uc; }

};

#endif // SKILLING_H  

