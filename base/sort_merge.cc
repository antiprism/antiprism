/*
   Copyright (c) 2007,2011 Roger Kaufman

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
   Name: sort_merge.cc
   Description: Sort vertices and faces indexes of an OFF file.
                Merge coincident vertices and/or faces.
   Project: Antiprism - http://www.antiprism.com
*/

#include <string.h>
#include <string>
#include <vector>
#include <algorithm>

#include "geom.h"
#include "transforms.h"
#include "math_utils.h"
#include "coloring.h"

using std::string;
using std::vector;


bool polygon_sort(vector<int> &polygon)
{
   bool reversed = false;
   vector<int> polygon_sorted;

   // Find lowest index make it the first polygon vertex index
   vector<int>::iterator iter = min_element(polygon.begin(), polygon.end()); 
   polygon_sorted.push_back(*iter);
   
   // Find rest of polygon indexes
   if ( iter < polygon.end() )
      polygon_sorted.insert(polygon_sorted.end(), (iter+1), polygon.end());
   if ( iter > polygon.begin() )
      polygon_sorted.insert(polygon_sorted.end(), polygon.begin(), iter);

   // reverse them if necessary
   iter = polygon_sorted.begin()+1;
   if ( *iter > polygon_sorted.back() ) {
      reverse(iter, polygon_sorted.end());
      reversed = true;
   }

   polygon.clear();
   polygon.assign(polygon_sorted.begin(), polygon_sorted.end());
   return reversed;
}

class vertexMap
{
public:
   int old_vertex;
   int new_vertex;
   vertexMap(int o, int n) : old_vertex(o), new_vertex(n) {}
};

bool cmp_vertexMap(const vertexMap &a, const vertexMap &b)
{
      return a.old_vertex < b.old_vertex;
}

void remap_faces(vector<vector<int> > &faces, const vector<vertexMap> &vm)
{
   for(unsigned int i=0;i<faces.size();i++)
      for(unsigned int j=0;j<faces[i].size();j++)
         faces[i][j] = vm[faces[i][j]].new_vertex;
}

class facesSort
{
public:
   vector<int> face;
   vector<int> face_with_cv;
   vector<int> equiv_faces;
   col_val color;
   int face_no;
   bool deleted;
   bool reversed;
   facesSort(vector<int> f, int i, bool rev) :
      face(f), face_no(i), reversed(rev) { deleted = false; }
   facesSort(vector<int> f, vector<int> fcv, int i, bool rev) :
      face(f), face_with_cv(fcv), face_no(i), reversed(rev) { deleted = false; }
};

bool cmp_faces(const facesSort &a, const facesSort &b)
{
   if ( a.face.size() != b.face.size() )
      return a.face.size() < b.face.size();
   else
      return a.face < b.face;
}

col_val average_color(const vector<col_val> &cols, const int &blend_type)
{
   if (blend_type == 1)
      return cols[0]; // first color
   else
   if (blend_type == 2)
      return cols[cols.size()-1]; // last color
   else
   if (blend_type == 3)
      return blend_RGB_centroid(cols); // simple RGB
   else
   if (blend_type == 4)
      return blend_RGB_centroid(cols,3,true); // RYB mode

   // should not get here
   return col_val();
}

col_val average_face_color(const col_geom &cg, const vector<facesSort> &fs, const char &elem, const int &begin, const int &end, const int &blend_type)
{
   // if only one instance, return its own color
   // if blend_type is 0 then do not continue with averaging. return the first color if more than one
   if (!(end-begin) || !blend_type) {
      if (elem=='f')
         return cg.get_f_col(fs[begin].face_no);
      else
         return cg.get_e_col(fs[begin].face_no);
   }
   
   vector<col_val> cols;  
   bool unset_found = false;
   for(int i=begin;i<=end;i++) {
      col_val col;
      if(elem=='f')
         col = cg.get_f_col(fs[i].face_no);
      else
         col = cg.get_e_col(fs[i].face_no);

      if (col.is_inv())
         continue;
      else
      if (!col.is_set() || col.is_idx())
         unset_found = true;
      else
         cols.push_back(col);
   }
   
   return (cols.size() ? average_color(cols, blend_type) : ((unset_found) ? col_val() : col_val(col_val::invisible)));
}

void sort_faces(geom_if &geom, const vector<vertexMap> &vm_no_cv, const vector<vertexMap> &vm_with_cv,
                const string &delete_elems, const char &elem,
                map<int, set<int> > *equiv_elems = 0,
                int blend_type = 1)
{
   vector<vector<int> > &faces = (elem=='f') ? geom.raw_faces() : geom.raw_edges();
   col_geom *cg = dynamic_cast<col_geom *>(&geom);

   vector<facesSort> fs;
   vector<col_val> cols;

   // make a copy of unmapped faces with coincident vertices if necessary
   vector<vector<int> > faces_with_cv;
   if ( vm_with_cv.size() ) {
      faces_with_cv = faces;
      remap_faces(faces_with_cv, vm_with_cv);
   }

   // use geom's face's and not a copy for no_cv
   if ( vm_no_cv.size() )
      remap_faces(faces, vm_no_cv);

   // load face sort vector
   for(unsigned int i=0;i<faces.size();i++) {
      bool reversed;
      reversed = polygon_sort(faces[i]);
      if ( !faces_with_cv.size() )
         fs.push_back(facesSort(faces[i],i, reversed));
      else
      {
         polygon_sort(faces_with_cv[i]);
         fs.push_back(facesSort(faces[i],faces_with_cv[i],i, reversed));
      }
   }
   
   // check to see we are actually deleting elements
   // erase coincident elements if any
   bool deleting_faces = ( (elem=='e' && strchr(delete_elems.c_str(), 'e')) || (elem=='f' && strchr(delete_elems.c_str(), 'f')) );

   // don't need to find coinicident faces/edges if no map provided
   if ( vm_no_cv.size() ) {
      stable_sort(fs.begin(), fs.end(), cmp_faces);

      if ( deleting_faces ) {
         if(equiv_elems && fs.size())
            fs[0].equiv_faces.push_back(fs[0].face_no);
         int cur_undeleted = 0;  // the first face in a set of equivalent faces
         unsigned int j = 0;
         for(unsigned int i=0;i<fs.size()-1;i++) {
            j=i+1;
            if ( fs[i].face == fs[j].face )
               fs[j].deleted = true;
            else {
               if(cg && !equiv_elems)
                  fs[cur_undeleted].color = average_face_color(*cg, fs, elem, cur_undeleted, i, blend_type);
               cur_undeleted = j;
            }

            if(equiv_elems)
               fs[cur_undeleted].equiv_faces.push_back(fs[j].face_no);
         }
         if(cg && !equiv_elems)
            fs[cur_undeleted].color = average_face_color(*cg, fs, elem, cur_undeleted, j, blend_type);
      }
   }

   // if coincient vertices are included we need to use/re-sort by those faces
   // but deleted flag will be preserved
   if ( !strchr(delete_elems.c_str(), 'v') ) {
      for(unsigned int i=0;i<fs.size();i++)
         fs[i].face = fs[i].face_with_cv;
      stable_sort(fs.begin(), fs.end(), cmp_faces);
   }

   // rewrite sorted face list, record colors for rebuilding color table
   bool restore_orient = (strchr(delete_elems.c_str(), 's')==0);
   faces.clear();
   for(unsigned int i=0;i<fs.size();i++) {
      if ( !fs[i].deleted ) {
         if(cg && !equiv_elems) {
            if (deleting_faces)
               cols.push_back(fs[i].color);
            else {
               if(elem=='f')
                  cols.push_back(cg->get_f_col(fs[i].face_no));
               else
                  cols.push_back(cg->get_e_col(fs[i].face_no));
            }
         }

         if(restore_orient && fs[i].reversed)
            reverse(fs[i].face.begin(), fs[i].face.end());
         faces.push_back(fs[i].face);
         if(equiv_elems) {
            (*equiv_elems)[faces.size()-1].insert(fs[i].equiv_faces.begin(), fs[i].equiv_faces.end());
         }
      }
   }

   // rebuild color table
   if(cg && !equiv_elems) {
      if(elem=='f') {
         cg->clear_f_cols();
         for(unsigned int i=0;i<cols.size();i++)
            cg->set_f_col(i,cols[i]);
      }
      else {
         cg->clear_e_cols();
         for(unsigned int i=0;i<cols.size();i++)
            cg->set_e_col(i,cols[i]);
      }
   }
}

class vertSort
{
public:
   vec3d vert;
   col_val color;
   int vert_no;
   bool deleted;
   vertSort(vec3d v, int i) : vert(v), vert_no(i) { deleted = false; }
};

bool cmp_verts(const vertSort &a, const vertSort &b, const double &eps)
{
   return (compare(a.vert,b.vert,eps)<0);
}

class vert_cmp
{
public:
   double eps;
   vert_cmp(double ep): eps(ep) {}
   bool operator() (const vertSort &a, const vertSort &b) { return cmp_verts(a, b, eps); }
};

col_val average_vert_color(const col_geom &cg, const vector<vertSort> &vs, const int &begin, const int &end, const int &blend_type)
{
   // if only one instance, return its own color
   // if blend_type is 0 then do not continue with averaging. return the first color if more than one
   if (!(end-begin) || !blend_type)
      return cg.get_v_col(vs[begin].vert_no);

   vector<col_val> cols;  
   bool unset_found = false;
   for(int i=begin;i<=end;i++) {
      col_val col = cg.get_v_col(vs[i].vert_no);

      if (col.is_inv())
         continue;
      else
      if (!col.is_set() || col.is_idx())
         unset_found = true;
      else
         cols.push_back(col);
   }
   
   return (cols.size() ? average_color(cols, blend_type) : ((unset_found) ? col_val() : col_val(col_val::invisible)));
}

void sort_vertices(geom_if &geom, vector<vertexMap> &vm_no_cv, vector<vertexMap > &vm_with_cv, 
                   const string &delete_elems,
                   map<int, set<int> > *equiv_elems = 0,
                   int blend_type = 1,  double eps = epsilon)
{
   vector<vec3d> &verts = geom.raw_verts();
   col_geom *cg = dynamic_cast<col_geom *>(&geom);

   vector<vertSort> vs;
   vector<col_val> cols;

   for(unsigned int i=0;i<verts.size();i++)
      vs.push_back(vertSort(verts[i],i));
   stable_sort(vs.begin(), vs.end(), vert_cmp(eps));

   // build map from old to new vertices with coincident vertices included
   // if v is specified we will not be using it
   bool merge_verts = strchr(delete_elems.c_str(), 'v');
   if ( !merge_verts ) {
      for(unsigned int i=0;i<vs.size();i++)
         vm_with_cv.push_back(vertexMap(vs[i].vert_no,i));
      stable_sort(vm_with_cv.begin(), vm_with_cv.end(), cmp_vertexMap);
   }

   // build map for coincident vertices if any deletion is to be done in vef.
   if ( strlen(delete_elems.c_str()) ) {
      int v = 0;
      int cur_undeleted = 0;  // the first face in a set of equivalent verts
      unsigned int j = 0;
      for(unsigned int i=0;i<vs.size()-1;i++) {
         j=i+1;
         vm_no_cv.push_back(vertexMap(vs[i].vert_no,v));
         if(!compare(vs[i].vert, vs[j].vert, eps)) {
            // don't set delete flag if we will not be using it
            if ( merge_verts )
               vs[j].deleted = true;
         }
         else {
            v++;
            if(cg && !equiv_elems)
               vs[cur_undeleted].color = average_vert_color(*cg, vs, cur_undeleted, i, blend_type);
            cur_undeleted = j;
         }
      }
      if(cg && !equiv_elems)
         vs[cur_undeleted].color = average_vert_color(*cg, vs, cur_undeleted, j, blend_type);

      vm_no_cv.push_back(vertexMap(vs.back().vert_no,v));
      stable_sort(vm_no_cv.begin(), vm_no_cv.end(), cmp_vertexMap);
   }

   // rewrite sorted vertex list, record colors for rebuilding color table
   verts.clear();
   for(unsigned int i=0;i<vs.size();i++) {
      if ( !vs[i].deleted ) {
         if(cg && !equiv_elems) {
            if (merge_verts)
               cols.push_back(vs[i].color);
            else
               cols.push_back(cg->get_v_col(vs[i].vert_no));
         }
         
         verts.push_back(vs[i].vert);
      }
      if(equiv_elems)
         (*equiv_elems)[verts.size()-1].insert(vs[i].vert_no);
   }

   // rebuild color table
   if(cg && !equiv_elems) {
      cg->clear_v_cols();
      for(unsigned int i=0;i<cols.size();i++)
         cg->set_v_col(i,cols[i]);
   }
}

bool sort_merge_elems(geom_if &geom, const string &merge_elems, vector<map<int, set<int> > > *equiv_elems, bool chk_congruence, int blend_type, double eps)
{
   // an empty geom cannot be processed
   if (!geom.verts().size())
      return false;
      
   // Congruence check expects a polyhedron with elements occurring
   // in coincident pairs, and exits early if this is not the case
   vector<map<int, set<int> > > equiv_elems_tmp;
   if(chk_congruence && !equiv_elems)
      equiv_elems = &equiv_elems_tmp;

   if(equiv_elems) {
      equiv_elems->clear();
      equiv_elems->resize(3);
   }

   unsigned int num_verts = geom.verts().size();
   unsigned int num_edges = geom.edges().size();
   unsigned int num_faces = geom.faces().size();

   vector<vertexMap> vm_no_cv, vm_with_cv;
   sort_vertices(geom, vm_no_cv, vm_with_cv, merge_elems, (equiv_elems ? &(*equiv_elems)[0] : 0), blend_type, eps);
   if(chk_congruence && (*equiv_elems)[0].size()*2 != num_verts)
      return false;

   if(geom.edges().size())
      sort_faces(geom, vm_no_cv, vm_with_cv, merge_elems, 'e', (equiv_elems ? &(*equiv_elems)[1] : 0), blend_type);
   if(chk_congruence && (*equiv_elems)[1].size()*2 != num_edges)
      return false;
   
   if(geom.faces().size())
      sort_faces(geom, vm_no_cv, vm_with_cv, merge_elems, 'f', (equiv_elems ? &(*equiv_elems)[2] : 0), blend_type);
   if(chk_congruence && (*equiv_elems)[2].size()*2 != num_faces)
      return false;
   
   return true;
}

void sort_merge_elems(geom_if &geom, const string &merge_elems, vector<map<int, set<int> > > *equiv_elems, double eps)
{
   sort_merge_elems(geom, merge_elems, equiv_elems, false, 0, eps);
}

void sort_merge_elems(geom_if &geom, const string &merge_elems, const int &blend_type, double eps)
{
   vector<map<int, set<int> > > *equiv_elems=0;
   sort_merge_elems(geom, merge_elems, equiv_elems, false, blend_type, eps);
}

void sort_merge_elems(geom_if &geom, const string &merge_elems, double eps)
{
   vector<map<int, set<int> > > *equiv_elems=0;
   sort_merge_elems(geom, merge_elems, equiv_elems, false, 0, eps);
}

bool check_congruence(const geom_if &geom1, const geom_if &geom2, vector<map<int, set<int> > > *equiv_elems, double eps)
{
   col_geom_v geom = geom1;
   geom.append(geom2);
   int ret = sort_merge_elems(geom, "vef", equiv_elems, true, 0, eps);
   /*
   for(int i=0; i<3; i++) {
      const char *elems[] = { "verts", "edges", "faces" };
      fprintf(stderr, "\n%s\n", elems[i]);
      map<int, set<int> >::iterator mi;
      for(mi=(*equiv_elems)[i].begin(); mi!=(*equiv_elems)[i].end(); ++mi) {
         fprintf(stderr, "%d <- ", mi->first);
         for(set<int>::iterator j=mi->second.begin(); j!=mi->second.end(); j++)
            fprintf(stderr, "%d  ", *j);
         fprintf(stderr, "\n");
      }
   }
   */
   return ret;
}

void get_congruence_maps(const geom_if &geom, mat3d trans,
      vector<vector<int > > &elem_maps, double eps)
{
   elem_maps.resize(3);
   col_geom_v tmp = geom;
   tmp.transform(trans);
   vector<map<int, set<int> > > equiv_elems;
   check_congruence(geom, tmp, &equiv_elems, eps);
   int cnts[3];
   cnts[0]=geom.verts().size();
   cnts[1]=geom.edges().size();
   cnts[2]=geom.faces().size();
   for(int i=0; i<3; i++) {
      elem_maps[i].resize(cnts[i]);
      map<int, set<int> >::iterator mi;
      for(mi=equiv_elems[i].begin(); mi!=equiv_elems[i].end(); ++mi) {
         int to = *mi->second.begin();
         if(to>=cnts[i])
            to -= cnts[i];
         elem_maps[i][to] = to;
         set<int>::iterator si;
         for(si=mi->second.begin(); si!=mi->second.end(); ++si) {
            int from = (*si<cnts[i]) ? *si : *si-cnts[i];
            if(to != from) {
               elem_maps[i][to] = from;
               break;
            }
         }
      }
   }
   /*
   for(int i=0; i<3; i++) {
      fprintf(stderr, "%d:\n", i);
      for(unsigned int v=0; v<elem_maps[i].size(); ++v) {
         fprintf(stderr, "   %d -> %d", v, elem_maps[i][v]);
         fprintf(stderr, "\n");
      }
      fprintf(stderr, "\n");
   }
   */
}


