/*
   Copyright (c) 2003-2022, Roger Kaufman

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
   Name: stellate_common.h
   Description: stellate code shared in /src
   Project: Antiprism - http://www.antiprism.com
*/

#ifndef STELLATE_COMMON_H
#define STELLATE_COMMON_H

#include <cstdio>
#include <set>
#include <string>
#include <vector>

using std::set;
using std::string;
using std::vector;

#include "../base/antiprism.h"

void color_stellation(anti::Geometry &geom, const char, const char, const char,
                      const anti::Color &, const anti::Color &,
                      const anti::Color &, const int, const string &,
                      const string);

#endif // STELLATE_COMMON_H
