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

/* \file crds_read.cc
   \brief Read and write OFF and coordinate files
*/

#include "private_off_file.h"
#include "utils.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

using std::string;
using std::vector;

void crds_file_read(FILE *ifile, Geometry &geom, char *first_line)
{
  int read_ret = 0;
  char *line = first_line;
  if (!line)
    read_ret = read_off_line(ifile, &line);

  while (read_ret == 0) {
    char *l = line;
    while (*l) {
      if (*l == ',' || isspace(*l))
        *l = ' ';
      l++;
    }

    double cs[4];
    int num = sscanf(line, "%lf %lf %lf %lf", &cs[0], &cs[1], &cs[2], &cs[3]);
    if (num == 3)
      geom.add_vert(Vec3d(cs[0], cs[1], cs[2]));

    free(line);
    read_ret = read_off_line(ifile, &line);
  }
  free(line);
}
