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

/*!\file std_polys.h
 * \brief Platonic and Archimedean polyhedra
 */

#ifndef STD_POLYS_H
#define STD_POLYS_H


#include <stdio.h>
#include <string.h>

#include "geom.h"
#include "transforms.h"

using std::string;


///A basic function type that makes a model.
typedef void (*model_func) (geom_if &);

///Put face list into a normalised form.
/**Sometimes a model is created and the face list may not be consistent
 * between machines, perhaps for numeric reasons. This function puts
 * the list in a consistent form.
 * \param geom the geometry whose faces list will be put into the
 * normalized form */
void normalised_face_list(geom_if &geom);

bool make_resource_geom(geom_if &geom, string name, char *errmsg=0);

void set_resource_polygon_color(geom_if &geom);
string expand_abbrevs(const string &name,
      const char *abbrevs[][2], size_t last);

struct Uniform {
	string Wythoff;
	short Kaleido, Coxeter, Wenninger;
   string name;
   string dual;
};

class uni_poly
{
   private:
      int last_uniform;
      Uniform* uniform;

   public:      
      uni_poly() { last_uniform = get_uniform_table(&uniform); }
      int get_poly(geom_if &geom, int sym_no);
      int lookup_sym_no(string sym, int is_dual);
      int get_last_uniform() { return last_uniform; }
      int get_uniform_table(Uniform **uniform);
};



struct JohnsonItem {
   model_func pfunc;
   const char *name;
};

class j_poly
{
   private:
      int last_J;
      JohnsonItem* J_items;

   public:
     
      j_poly();
      int get_poly(geom_if &geom, int sym);
      int lookup_sym_no(string sym);
      int get_last_J() { return last_J; }

};



#endif // STD_POLYS_H  

