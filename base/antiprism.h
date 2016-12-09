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

/*!\file antiprism.h
 * \brief Includes all the headers need to use the Antiprism library
 */

/**\mainpage Antiprism Library Documentation
 * \section features Features
 *
 * The Antiprism library provides a framework for working
 * with polyhedra. It includes:
 * - import of OFF, LG3D and simple text coordinate files
 * - export of OFF, VRML, POV, LG3D and text coordinate files
 * - vector and matrix operations
 * - geometric utilities
 * - analysis of polyhedra
 * - polyhedron operations
 * - colouring operations
 * - symmetry operations
 * - models of commonly used polyhedra
 *
 */

#ifndef ANTIPRISM_H
#define ANTIPRISM_H

#include "boundbox.h"
#include "color.h"
#include "coloring.h"
#include "colormap.h"
#include "const.h"
#include "displaypoly.h"
#include "elemprops.h"
#include "geometry.h"
#include "geometryinfo.h"
#include "geometryutils.h"
#include "getopt.h"
#include "mathutils.h"
#include "normal.h"
#include "polygon.h"
#include "povwriter.h"
#include "random.h"
#include "scene.h"
#include "status.h"
#include "symmetry.h"
#include "timer.h"
#include "trans3d.h"
#include "trans4d.h"
#include "utils.h"
#include "vec3d.h"
#include "vec4d.h"
#include "vec_utils.h"
#include "vrmlwriter.h"

#endif // ANTIPRISM_H
