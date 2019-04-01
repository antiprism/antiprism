/*
   Copyright (c) 2008-2016, Adrian Rossiter

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

/**\file DisplayPoly.h
   \brief support for displaying a model as a polyhedron
*/

#ifndef DISP_POLY_GL_H
#define DISP_POLY_GL_H

#include "../base/antiprism.h"

using namespace anti;

class DisplayPoly_gl : public virtual DisplayPoly {
private:
  bool show_orientation;
  int transparency_type;

protected:
  void gl_verts(const Scene &scen);
  void gl_edges(const Scene &scen);
  void gl_faces(const Scene &scen);

public:
  enum { trans_model = 0, trans_50pc, trans_0pc };
  DisplayPoly_gl();

  GeometryDisplay *clone() const { return new DisplayPoly_gl(*this); };
  void gl_geom(const Scene &scen);
  void set_show_orientation(bool show = true) { show_orientation = show; }
  bool get_show_orientation() { return show_orientation; }
  void set_transparency_type(int typ) { transparency_type = typ; }
  int get_transparency_type() { return transparency_type; }
};

class DisplayNumLabels_gl : public virtual DisplayNumLabels {
private:
  bool use_alt_labels = false;
  float text_scale = 1.0;
  void gl_faces(const Scene &scen);
  void gl_verts(const Scene &scen);
  void gl_edges(const Scene &scen);

public:
  GeometryDisplay *clone() const { return new DisplayNumLabels_gl(*this); };
  void gl_geom(const Scene &scen);
  void write_label(const Scene &scen, char *label, Vec3d pos,
                   Vec3d norm = Vec3d());
  void set_use_alt_labels(bool use_alt = true) { use_alt_labels = use_alt; }
  bool get_use_alt_labels() { return use_alt_labels; }
  void set_text_scale(float scale) { text_scale = scale; }
  float get_text_scale() { return text_scale; }
};

class DisplaySymmetry_gl : public virtual DisplaySymmetry,
                           public virtual DisplayPoly_gl {
public:
  GeometryDisplay *clone() const { return new DisplaySymmetry_gl(*this); }
  void gl_geom(const Scene &scen);
};

#endif // DISP_POLY_GL_H
