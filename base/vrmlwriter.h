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

/**\file vrmlwriter.h
   \brief write a VRML file
*/

#ifndef VRMLWRITER_H
#define VRMLWRITER_H

#include "scene.h"

#include <cstdio>
#include <string>

namespace anti {

// conversion functions
std::string vrml_vec(double x, double y, double z, int sig_digits = 10);
inline std::string vrml_vec(const Vec3d &v, int sig_digits = 10)
{
  return vrml_vec(v[0], v[1], v[2], sig_digits);
}
std::string vrml_col(const Color &col);

class VrmlWriter {
private:
  void header(FILE *ofile);
  void scene_header(FILE *ofile, Scene &scen);
  void cameras(FILE *ofile, Scene &scen);
  void geometry_headers(FILE *ofile, Scene &scen);
  void geometry_objects(FILE *ofile, const Scene &scen, int sig_digits);

public:
  void write(FILE *ofile, Scene &scen, int sig_digits);
};

} // namespace anti

#endif // VRMLWRITER_H
