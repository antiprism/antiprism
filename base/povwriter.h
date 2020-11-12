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
   Name: povwriter.h
   Description: write a POV file
   Project: Antiprism - http://www.antiprism.com
*/

#ifndef POVWRITER_H
#define POVWRITER_H

#include "scene.h"
#include "utils.h"

#include <cstdio>
#include <string>
#include <vector>

namespace anti {

// conversion functions
std::string pov_vec(double x, double y, double z, int sig_digits = 10);
inline std::string pov_vec(const Vec3d &v, int sig_digits = 10)
{
  return pov_vec(v[0], v[1], v[2], sig_digits);
}
std::string pov_col(const Color &col);

class PovWriter : public SceneItem {
private:
  char o_type;
  int stereo_type;
  int shadow;
  double text_size;
  Color text_colour;

  std::vector<std::string> includes;
  std::vector<std::string> obj_includes;

  void scene_header(FILE *ofile, const Scene &scen);
  void common_macros(FILE *ofile);
  void include_files(FILE *ofile, const Scene &scen);
  void cameras(FILE *ofile, const Scene &scen);
  void scene_setup(FILE *ofile);
  void stereo_setup(FILE *ofile);
  void geometry_objects(FILE *ofile, const Scene &scen, int sig_dgts);
  void stereo_clipping(FILE *ofile);
  void camera_lights(FILE *ofile);

public:
  PovWriter()
      : o_type('a'), stereo_type(0), shadow(1), text_size(1.0),
        text_colour(Color(10, 10, 10))
  {
  }

  void set_file_name(std::string fname) { set_name(dots2underscores(fname)); }
  void set_o_type(char otype);
  int get_o_type() { return o_type; }

  void set_stereo_type(int s_type) { stereo_type = s_type; }
  int get_stereo_type() { return stereo_type; }
  void set_shadow(int shad) { shadow = shad; }
  int get_shadow() { return shadow; }

  void set_text_size(double size) { text_size = size; }
  void set_text_colour(Vec3d col) { text_colour = col; }

  void set_includes(std::vector<std::string> incs) { includes = incs; }
  void set_obj_includes(std::vector<std::string> incs) { obj_includes = incs; }

  void write(FILE *ofile, const Scene &scen, int sig_dgts = 10);
};

} // namespace anti

#endif // POVWRITER_H
