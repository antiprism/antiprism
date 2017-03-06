/*
   Copyright (c) 2003-2016, Adrian Rossiter

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

/*
   Name: gl_writer.h
   Description: write a OpenGL output
   Project: Antiprism - http://www.antiprism.com
*/

#ifndef GL_WRITER_H
#define GL_WRITER_H

#include "../base/antiprism.h"
#include "gl_stroke_font.h"

using namespace anti;

// conversion functions
inline void to_gl_matrix(double *gl_m, const Trans3d &trans)
{
  Trans3d m = trans.transpose();
  for (int i = 0; i < 16; i++)
    gl_m[i] = m[i];
}

class gl_writer {
private:
  double text_size;
  Color text_colour;

  void cur_camera(const Camera &cam);
  void geometry_objects(const Scene &scen);

public:
  gl_writer() : text_size(1.0), text_colour(Color(10, 10, 10)) {}

  void set_text_size(double size) { text_size = size; }
  void set_text_colour(Vec3d col) { text_colour = col; }

  void init(const Scene &scen);
  void write(const Scene &scen);
};

#endif // GL_WRITER_H
