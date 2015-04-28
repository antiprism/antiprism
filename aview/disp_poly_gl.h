/*
   Copyright (c) 2008, Adrian Rossiter

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

/**\file disp_poly.h
   \brief support for displaying a model as a polyhedron
*/


#ifndef DISP_POLY_GL_H
#define DISP_POLY_GL_H

#include "../base/antiprism.h"

class disp_poly_gl : public virtual disp_poly
{
   private:
      bool show_orientation;

   protected:
      void gl_verts(const scene &scen);
      void gl_edges(const scene &scen);
      void gl_faces(const scene &scen);

   public:
      disp_poly_gl();

      geom_disp *clone() const { return new disp_poly_gl(*this); };
      void gl_geom(const scene &scen);
      void set_show_orientation(bool show=true) { show_orientation = show; }
      bool get_show_orientation() { return show_orientation; }
};


class disp_num_labels_gl : public virtual disp_num_labels
{
   private:
      void gl_faces(const scene &scen);
      void gl_verts(const scene &scen);
      void gl_edges(const scene &scen);

   public:
      geom_disp *clone() const { return new disp_num_labels_gl(*this); };
      void gl_geom(const scene &scen);
 };


class disp_sym_gl: public virtual disp_sym, public virtual disp_poly_gl
{
   public:
      geom_disp *clone() const { return new disp_sym_gl(*this); }
      void gl_geom(const scene &scen);

};


#endif // DISP_POLY_GL_H


