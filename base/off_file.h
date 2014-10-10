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

/*!\file off_file.h
   \brief Read and write OFF and coordinate files
*/


#ifndef OFF_FILE_H
#define OFF_FILE_H

#include <stdio.h>
#include "geom.h"


int read_off_line(FILE *fp, char **line);


bool crds_file_read(string file_name, geom_if &geom, char *errmsg=0);
void crds_file_read(FILE *ifile, geom_if &geom, char *first_line=0);

bool crds_write(string file_name,  const geom_if &geom, char *errmsg=0,
      const char *sep=" ", int sig_dgts=DEF_SIG_DGTS);
void crds_write(FILE *ofile,  const geom_if &geom,
       const char *sep=" ", int sig_dgts=DEF_SIG_DGTS);
       
bool obj_write(string file_name,  const geom_if &geom, char *errmsg=0,
      const char *sep=" ", int sig_dgts=DEF_SIG_DGTS);
void obj_write(FILE *ofile,  const geom_if &geom,
       const char *sep=" ", int sig_dgts=DEF_SIG_DGTS);

bool off_file_read(string file_name, geom_if &geom, char *errmsg=0);
bool off_file_read(FILE *ifile, geom_if &geom, char *errmsg=0);

bool off_file_write(string file_name, const geom_if &geom, char *errmsg=0,
      int sig_dgts=DEF_SIG_DGTS);
void off_file_write(FILE *ofile, const geom_if &geom,
      int sig_dgts=DEF_SIG_DGTS);

bool off_file_write(string file_name, const vector<const geom_if *> &geoms,
      char *errmsg=0, int sig_dgts=DEF_SIG_DGTS);
void off_file_write(FILE *ofile, const vector<const geom_if *> &geoms,
      int sig_dgts=DEF_SIG_DGTS);

#endif // OFF_FILE_H
   
