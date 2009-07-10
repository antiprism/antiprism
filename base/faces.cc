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
   Name: faces.cc
   Description: face manipulations
   Project: Antiprism - http://www.antiprism.com
*/

#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <functional>
#include <vector>
#include <map>
#include <string>

#include "geom.h"


using std::vector;
using std::map;
using std::pair;
using std::swap;
using std::string;



// close a polyhedron, no more than two open edges per vertex
bool close_poly_basic(geom_if &geom)
{
   map<vector<int>, int> edges;
   map<vector<int>, int>::iterator mi, mi2;
   vector<int> edge(2);
   vector<int> r_edge(2);
   for(unsigned int i=0; i<geom.faces().size(); ++i)
      for(unsigned int j=0; j<geom.faces(i).size(); ++j) {
         edge[0] = geom.faces(i, j);
         edge[1] = geom.faces(i, (j+1)%geom.faces(i).size());
         mi = edges.find(edge);
         if(mi!=edges.end())
            mi->second++;
         else {
            r_edge[0] = edge[1];
            r_edge[1] = edge[0];
            mi2 = edges.find(r_edge);
            if(mi2!=edges.end())
               mi2->second++;
            else
               edges[r_edge] = 1;  // edge is oriented for missing face
         }
      }
  
   map<int, vector<int> > neighbours;
   map<int, vector<int> >::iterator ni;
   for(mi=edges.begin(); mi!=edges.end(); mi++) {
      if(mi->second == 1) {
         for(int i=0; i<2; i++) {
            const int from = mi->first[i];
            const int to = mi->first[(i+1)%2];
            ni = neighbours.find(from);
            if(ni == neighbours.end()) {
               vector<int> idxs(2, -1);
               idxs[i] = to;
               neighbours[from] = idxs;
            }
            else {
               if(ni->second[0]<0)
                  ni->second[0] = to;
               else if(ni->second[1]<0)
                  ni->second[1] = to;
               else                      // three open edges at a vertex
                  return false;
            }
         }
      }
   }
   
   while(neighbours.size()) {
      vector<int> face;
      ni = neighbours.begin();
      int prev = -1;
      int curr = ni->first;
      face.push_back(curr);
      while(true) {
         int new_curr = (ni->second[0]!=prev) ? ni->second[0] : ni->second[1];
         prev = curr;
         curr = new_curr;
         neighbours.erase(ni);
         ni = neighbours.find(curr);
         if(ni==neighbours.end())  // completed face
            break;
         face.push_back(curr);
      }
      geom.add_face(face);
   }
   return true;
}


bool face_bond(geom_if &geom, geom_if &bgeom, int f, int b_f, int off, bool merge)
{
   //transform(bgeom, mat3d::rot(1.0001, 0.0002, 0.0003));

   vector<vec3d> pts;
   pts.push_back(geom.verts(geom.faces(f, 0)));
   pts.push_back(geom.verts(geom.faces(f, 1)));
   pts.push_back(geom.verts(geom.faces(f, 2)));
   int fsz = bgeom.faces(b_f).size();
   vector<vec3d> bpts;
   bpts.push_back(bgeom.verts(bgeom.faces(b_f, off%fsz)));
   bpts.push_back(bgeom.verts(bgeom.faces(b_f, (off-1+fsz)%fsz)));
   bpts.push_back(bgeom.verts(bgeom.faces(b_f, (off-2+fsz)%fsz)));
   
   bgeom.transform(mat3d::alignment(bpts, pts));
   
   map<int, int> vmap;
   for(unsigned int i=0; i<bgeom.faces(b_f).size(); i++)
      vmap[bgeom.faces(b_f, (off-i+fsz)%fsz)+geom.verts().size()] =
            geom.faces(f, i);
   
   if(merge) {
      vector<int> del_faces(1);
      del_faces[0] = f;
      geom.delete_faces(del_faces);
      del_faces[0] = b_f;
      bgeom.delete_faces(del_faces);
      geom.append(bgeom);

      geom.verts_merge(vmap);
   }

   return true;
}


