/*
   Copyright (c) 2003, Adrian Rossiter

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
   Name: col_geom.cc
   Description: representation of colours in geometries
   Project: Antiprism - http://www.antiprism.com
*/

#include <stdlib.h>
#include <algorithm>

#include "col_geom.h"

void col_geom::remap_cols(map<int, col_val> &cols, const map<int, int> &chg_map)
{
   if(!chg_map.size())
      return;
   map<int, col_val> new_cols;
   map<int, int>::const_iterator mi;
   for(mi=chg_map.begin(); mi!=chg_map.end(); mi++) {
      if(mi->second != -1) {
         map<int, col_val>::iterator cmi = cols.find(mi->first);
         if(cmi!=cols.end())
            new_cols[mi->second] = cmi->second;
      }
   }

   cols = new_cols;
}


void col_geom::append(const col_geom &geom, int v_offset, int e_offset, int f_offset)
{
   int offs[] = {v_offset, e_offset, f_offset};
   for(int i=0; i<3; i++) {
      map<int, col_val>::const_iterator mi;
      for(mi=geom.elem_cols[i].begin(); mi!=geom.elem_cols[i].end(); mi++)
         elem_cols[i][mi->first + offs[i]] = mi->second;
   }
}
 
