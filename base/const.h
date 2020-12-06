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

/*!\file const.h
   \brief Global constants
*/

#ifndef CONST_H
#define CONST_H

#include <cmath>
#include <cstddef>

namespace anti {
// frequently used constants
static const double phi = (sqrt(5) + 1) / 2;         // 1.61803
static const double sqrt_phi_plus_2 = sqrt(phi + 2); // 1.90211
static const double sqrt_2 = sqrt(2);                // 1.41421
static const double sqrt_3 = sqrt(3);                // 1.73205

/// Less than this magnitude may be taken as zero
const double epsilon = 1e-12;

/**Less than this magnitude may be taken as zero when determining
 * symmetry elements */
const double sym_eps = 1e-5;

/// The default number of significant digits when writing numbers
const int DEF_SIG_DGTS = 16;

/// The default "infinity" distance to ignore far points when positioning camera
const double DEF_CAMERA_INF_DIST = 1000.0;

/// Flags for selecting which elements a function acts upon
enum {
  ELEM_NONE = 0,                                   ///< None
  ELEM_VERTS = 1,                                  ///< Vertices
  ELEM_EDGES = 2,                                  ///< Edges
  ELEM_FACES = 4,                                  ///< Faces
  ELEM_ALL = ELEM_VERTS | ELEM_EDGES | ELEM_FACES, ///< All elements
};

/// Indicate element type when an array holds data on all three element types
enum {
  VERTS = 0, ///< Vertices
  EDGES = 1, ///< Edges
  FACES = 2, ///< Faces
};

/// Flags to select inclusion
enum {
  INCLUSION_IN = 1, ///< In
  INCLUSION_ON = 2, ///< On
  INCLUSION_OUT = 4 ///< Out
};

/// Flags to select tesselator winding rules
enum {
  TESS_WINDING_ODD = 100130,         ///< Odd winding number
  TESS_WINDING_NONZERO = 100131,     ///< Nonzero winding number
  TESS_WINDING_POSITIVE = 100132,    ///< Positive winding number
  TESS_WINDING_NEGATIVE = 100133,    ///< Negative winding number
  TESS_WINDING_ABS_GEQ_TWO = 100134, ///< Absolute value of winding number >= 2
};

} // namespace anti

#endif // CONST_H
