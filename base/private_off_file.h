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

/* !\file private_off_file.h
   \brief Read and write OFF and coordinate files
*/

#ifndef PRIVATE_OFF_FILE_H
#define PRIVATE_OFF_FILE_H

#include "geometry.h"
#include <stdio.h>

using namespace anti;

int read_off_line(FILE *fp, char **line);

bool crds_file_read(std::string file_name, anti::Geometry &geom,
                    char *errmsg = nullptr);
void crds_file_read(FILE *ifile, anti::Geometry &geom,
                    char *first_line = nullptr);

bool crds_write(std::string file_name, const anti::Geometry &geom,
                char *errmsg = nullptr, const char *sep = " ",
                int sig_dgts = DEF_SIG_DGTS);
void crds_write(FILE *ofile, const anti::Geometry &geom, const char *sep = " ",
                int sig_dgts = DEF_SIG_DGTS);

bool obj_write(std::string file_name, const Geometry &geom,
               char *errmsg = nullptr, const char *sep = " ",
               int sig_dgts = DEF_SIG_DGTS);
void obj_write(FILE *ofile, const Geometry &geom, const char *sep = " ",
               int sig_dgts = DEF_SIG_DGTS);

bool off_file_read(std::string file_name, anti::Geometry &geom,
                   char *errmsg = nullptr);
bool off_file_read(FILE *ifile, anti::Geometry &geom, char *errmsg = nullptr);

bool off_file_write(std::string file_name, const anti::Geometry &geom,
                    char *errmsg = nullptr, int sig_dgts = DEF_SIG_DGTS);
void off_file_write(FILE *ofile, const anti::Geometry &geom,
                    int sig_dgts = DEF_SIG_DGTS);

bool off_file_write(std::string file_name,
                    const std::vector<const anti::Geometry *> &geoms,
                    char *errmsg = nullptr, int sig_dgts = DEF_SIG_DGTS);
void off_file_write(FILE *ofile,
                    const std::vector<const anti::Geometry *> &geoms,
                    int sig_dgts = DEF_SIG_DGTS);

#endif // PRIVATE_OFF_FILE_H
