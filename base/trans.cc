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
   Name: trans.cc
   Description: polyhedron transformations
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
#include "info.h"
#include "transforms.h"


using std::vector;
using std::map;
using std::pair;
using std::swap;
using std::string;


void edges_to_faces(const geom_if &geom, geom_if &egeom, bool face_verts)
{
   egeom.add_verts(geom.verts());
   for(unsigned int f=0; f<geom.faces().size(); f++)
      egeom.add_vert(face_verts ? geom.face_cent(f) : vec3d(0,0,0));
   
   map<vector<int>, vector<int> > edges;
   geom.get_edge_face_pairs(edges, false);
   map<vector<int>, vector<int> >::iterator mi;
   
   for(mi=edges.begin(); mi!=edges.end(); mi++) {
      vector<int> face(4);
      face[0] = mi->first[0];
      face[1] = mi->second[0]+geom.verts().size();
      face[2] = mi->first[1];
      face[3] = mi->second[1]+geom.verts().size();
      egeom.add_face(face);
   }
}


// truncate vertices specified by index number
void truncate_verts(geom_if &geom, vector<int> &v_idxs, double ratio, geom_info *info)
{
   geom_info tmp_info(geom);
   if(!info)
      info = &tmp_info;
  
   bool merge = (ratio==0.5);
   map<vector<int>, vector<int> > e_to_vs;  // vertex pairs to merge
   const vector<vec3d> &verts = geom.verts();
   map<vector<int>, vector<int> > v3maps;
   for(unsigned int i=0; i<v_idxs.size(); i++) {
      int idx = v_idxs[i];
      const vector<int> &cons = info->get_vert_cons()[idx];
      int orig_vsz = verts.size();
      vector<int> tface;
      for(unsigned int j=0; j<cons.size(); j++) {
         geom.add_vert(verts[idx] + ratio*(verts[cons[j]] - verts[idx]));
         tface.push_back(orig_vsz+j);
         // set up maps for faces conversion of old index to new index pair
         vector<int> v3(3);
         v3[0] = cons[j];
         v3[1] = idx;
         v3[2] = cons[(j+1)%cons.size()];
         vector<int> v2(2);
         v2[0] = orig_vsz + j;
         v2[1] = orig_vsz + (j+1)%cons.size();
         v3maps[v3] = v2;
         if(merge)
            e_to_vs[make_edge(v3[0], v3[1])].push_back(orig_vsz+j);
      }
      geom.add_face(tface);
   }

   for(unsigned int i=0; i<geom.faces().size(); i++) {
      vector<int> tface;
      int fsz = geom.faces(i).size();
      vector<vector<int> > v3s(fsz, vector<int>(3));
      for(int j=0; j<fsz; j++) {
         v3s[j][0] = geom.faces(i, j);
         v3s[j][1] = geom.faces(i, (j+1)%fsz);
         v3s[j][2] = geom.faces(i, (j+2)%fsz);
      }
      for(unsigned int j=0; j<v3s.size(); j++) {
         map<vector<int>, vector<int> >::iterator mi;
         if((mi=v3maps.find(v3s[j])) != v3maps.end()) {    // v3 in right order
            tface.push_back(mi->second[0]);
            tface.push_back(mi->second[1]);
         }
         else {
            swap(v3s[j][0], v3s[j][2]);
            if((mi=v3maps.find(v3s[j])) != v3maps.end()) { // v3 in rev order
               tface.push_back(mi->second[1]);
               tface.push_back(mi->second[0]);
            }
            else                                           // v3 not found
               tface.push_back(v3s[j][1]);
         }
      }
      geom.raw_faces()[i] = tface;
   }

   if(merge) {
      map<int, int> idx_map;
      map<vector<int>, vector<int> >::iterator mi;
      for(mi=e_to_vs.begin(); mi!=e_to_vs.end(); mi++) {
         if(mi->second.size()==2) // two vertices on an edge
            idx_map[mi->second[1]] = mi->second[0];
      }
      geom.verts_merge(idx_map);
   }
   geom.delete_verts(v_idxs);
}           

// truncate vertices with vertex order (default order = 0, truncate all)
void truncate_verts(geom_if &geom, double ratio, int order, geom_info *info)
{      
   geom_info tmp_info(geom);
   if(!info)
      info = &tmp_info;

   vector<int> v_idxs;
   const vector<vector<int> > &vcons = info->get_vert_cons();
   for(unsigned int i=0; i<vcons.size(); i++)
      if(order==0 || (int)vcons[i].size()==order)
         v_idxs.push_back(i);
   truncate_verts(geom, v_idxs, ratio, info);
}

