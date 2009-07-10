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
#include <string>

#include "geom.h"
#include "transforms.h"
#include "kaleido.h"

using std::string;


///A basic function type that makes a model.
typedef void (*model_func) (geom_if &);

///Get a model.
/**A utility funcion that returns a model.
 * \param pfunc a pointer to a function that uses a geometry argument
 * to return a model.
 * \return The model. */
geom_v model(model_func pfunc);

///Make a tetrahedron.
/**\param geom used to return the model */
void tetrahedron(geom_if &geom);

///Make a cube.
/**\param geom used to return the model */
void cube(geom_if &geom);

///Make an octahedron.
/**\param geom used to return the model */
void octahedron(geom_if &geom);

///Make a dodecahedron.
/**\param geom used to return the model */
void dodecahedron(geom_if &geom);

///Make an icosahedron.
/**\param geom used to return the model */
void icosahedron(geom_if &geom);


///Make a truncated tetrahedron.
/**\param geom used to return the model */
void tr_tetrahedron(geom_if &geom);

///Make a cuboctahedron.
/**\param geom used to return the model */
void cuboctahedron(geom_if &geom);

///Make a truncated cube.
/**\param geom used to return the model */
void tr_cube(geom_if &geom);

///Make a truncated octahedron.
/**\param geom used to return the model */
void tr_octahedron(geom_if &geom);

///Make a rhombicuboctahedron.
/**\param geom used to return the model */
void rhombicuboctahedron(geom_if &geom);

///Make a truncated cuboctahedron.
/**\param geom used to return the model */
void tr_cuboctahedron(geom_if &geom);

///Make a snub cube.
/**\param geom used to return the model */
void snub_cube(geom_if &geom);

///Make an icosidodecahedron.
/**\param geom used to return the model */
void icosidodecahedron(geom_if &geom);

///Make a truncated dodecahedron.
/**\param geom used to return the model */
void tr_dodecahedron(geom_if &geom);

///Make a truncated icosahedron.
/**\param geom used to return the model */
void tr_icosahedron(geom_if &geom);

///Make a rhombicosidodecahedron.
/**\param geom used to return the model */
void rhombicosidodecahedron(geom_if &geom);

///Make a truncated icosidodecahedron.
/**\param geom used to return the model */
void tr_icosidodecahedron(geom_if &geom);

///Make a snub dodecahedron.
/**\param geom used to return the model */
void snub_dodecahedron(geom_if &geom);

///Put face list into a normalised form.
/**Sometimes a model is created and the face list may not be consistent
 * between machines, perhaps for numeric reasons. This function puts
 * the list in a consistent form.
 * \param geom the geometry whose faces list will be put into the
 * normalized form */
void normalised_face_list(geom_if &geom);

// Inline function definitions

inline geom_v model(model_func pfunc)
{ 
   geom_v g;
   pfunc(g);
   return g;
}


class uni_poly
{
   private:
      int last_uniform;
      Uniform *uniform; 
      Polyhedron *poly;
      
   public:
      uni_poly() { last_uniform = get_uniform_list(&uniform); }
      int get_poly(geom_if &geom, int sym_no);
      bool ok() { return poly != 0; }
      void list_poly(int idx, FILE *fp=stderr);
      void list_polys(FILE *fp=stderr);
      int lookup_sym_no(string sym);
      int get_last_uniform() { return last_uniform; }

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
      void list_poly(int idx, FILE *fp=stderr);
      void list_polys(FILE *fp=stderr);
      int lookup_sym_no(string sym);
      int get_last_J() { return last_J; }

};



#endif // STD_POLYS_H  

