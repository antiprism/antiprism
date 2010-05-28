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

/*!\file const.h
   \brief Global constants
*/

#ifndef CONST_H
#define CONST_H

#include <stddef.h>

///Less than this magnitude may be taken as zero
const double epsilon = 1e-12;

/**Less than this magnitude may be taken as zero when determining
 * symmetry elements */
const double sym_eps = 1e-5;

///The size of the character array used for returning messages from functions
const size_t MSG_SZ = 256;

///The default number of significant digits when writing numbers
const int DEF_SIG_DGTS = 16;

///The default "infinity" distance to ignore far points when positioning camera
const int DEF_CAMERA_INF_DIST = 1000.0;

///Characters that separate parts of a resource name
const char RES_SEPARATOR[] = "_ ";

///For selecting which elements a function acts upon
const char ELEM_NONE = 0;
const char ELEM_VERTS = 1;
const char ELEM_EDGES = 2;
const char ELEM_FACES = 4;
const char ELEM_ALL = ELEM_VERTS | ELEM_EDGES | ELEM_FACES;

#endif // CONST_H

