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

/* \file Vec3d.cc
 *  \brief Vectors for 3D geometry
 *
 *  A vector class with common vector operations.
 */

#include "vec3d.h"
#include "utils.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace anti {

Vec3d Vec3d::X(1, 0, 0);
Vec3d Vec3d::Y(0, 1, 0);
Vec3d Vec3d::Z(0, 0, 1);
Vec3d Vec3d::zero(0, 0, 0);

Status Vec3d::read(const char *str)
{
  int i;
  double f[3];
  char buff;
  char buff2;
  if (sscanf(str, " %lf , %lf , %lf %c", f, f + 1, f + 2, &buff) != 3 &&
      sscanf(str, " %lf %lf %lf %c", f, f + 1, f + 2, &buff2) != 3)
    return Status::error("didn't find three numbers");

  for (i = 0; i < 3; i++) {
    if (std::isinf(f[i])) {
      const char *pos[] = {"first", "second", "third"};
      return Status::error(msg_str("%s number too large\n", pos[i]));
    }
  }

  v[0] = f[0];
  v[1] = f[1];
  v[2] = f[2];

  return Status::ok();
}

Status Vec3d::read_maths(const char *str)
{
  std::vector<double> coords;
  Status stat;
  if ((stat = read_double_list(str, coords))) {
    if (coords.size() == 1 && coords[0] == 0)
      *this = {0, 0, 0}; // set to zero/null/origin
    else if (coords.size() == 3)
      *this = coords.data(); // set to values calculated from maths expressions
    else
      stat.set_error(
          msg_str("must give 3 coordinates (%d given)", (int)coords.size()));
  }

  return stat;
}

std::string Vec3d::to_str(const char *sep, int sig_dgts) const
{
  if (!is_set())
    return "not set";

  const auto *v = get_v();
  if (sig_dgts > 0)
    return msg_str("%.*g%s%.*g%s%.*g", sig_dgts, v[0], sep, sig_dgts, v[1], sep,
                   sig_dgts, v[2]);
  else
    return msg_str("%.*f%s%.*f%s%.*f", -sig_dgts, v[0], sep, -sig_dgts, v[1],
                   sep, -sig_dgts, v[2]);
}

void Vec3d::dump(const char *var, FILE *file) const
{
  if (var)
    fprintf(file, "%s=", var);
  fprintf(file, "(%s)\n", to_str().c_str());
}

// make vector unusable
void Vec3d::unset() { v[0] = NAN; }

} // namespace anti
