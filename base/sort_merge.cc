/*
   Copyright (c) 2007-2012 Roger Kaufman

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
#include "info.h"

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
   if (iter < polygon.end())
      polygon_sorted.insert(polygon_sorted.end(), (iter+1), polygon.end());
   if (iter > polygon.begin())
      polygon_sorted.insert(polygon_sorted.end(), polygon.begin(), iter);

   // reverse them if necessary
   iter = polygon_sorted.begin()+1;
   if (*iter > polygon_sorted.back()) {
      reverse(iter, polygon_sorted.end());
      reversed = true;
   }

   polygon = polygon_sorted;
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
   int face_no;
   vector<int> face;
   vector<int> face_all_verts;
   col_val col;
   col_val average_col;
   bool deleted;
   bool reversed;

   // for congruence test
   vector<int> equiv_faces;

   facesSort(int i, vector<int> f, col_val c, bool rev) :
      face_no(i), face(f), col(c), reversed(rev) { deleted = false; }
   facesSort(int i, vector<int> f, vector<int> fcv, col_val c, bool rev) :
      face_no(i), face(f), face_all_verts(fcv), col(c), reversed(rev) { deleted = false; }
};

bool cmp_faces(const facesSort &a, const facesSort &b)
{
   if (a.face.size() != b.face.size())
      return a.face.size() < b.face.size();
   else
      return a.face < b.face;
}

bool cmp_face_no(const facesSort &a, const facesSort &b)
{
   return a.face_no < b.face_no;
}

col_val average_color(const vector<col_val> &cols, const int &blend_type)
{
   if (blend_type <= 1)
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

col_val average_face_color(const vector<facesSort> &fs, const char &elem, const int &begin, const int &end, const int &blend_type)
{
   // quick decision, if only one instance, return its own color
   if (!(end-begin)) {
      if (elem=='f')
         return fs[begin].col;
      else
         return fs[begin].col;
   }
   
   // collect colors
   vector<col_val> cols;  
   for(int i=begin;i<=end;i++) {
      col_val col;
      if (elem=='f')
         col = fs[i].col;
      else
         col = fs[i].col;
      cols.push_back(col);
   }
   
   return (cols.size() ? average_color(cols, blend_type) : col_val());
}

void sort_faces(geom_if &geom, const vector<vertexMap> &vm_all_verts, const vector<vertexMap> &vm_merged_verts,
                const string &delete_elems, const char &elem,
                map<int, set<int> > *equiv_elems = 0,
                int blend_type = 1)
{
   vector<vector<int> > &faces = (elem=='f') ? geom.raw_faces() : geom.raw_edges();
   col_geom *cg = dynamic_cast<col_geom *>(&geom);

   // sort only when 's' is set
   bool sort_only = strchr(delete_elems.c_str(), 's');
   bool merge_verts = strchr(delete_elems.c_str(), 'v');

   vector<facesSort> fs;

   // make a copy of faces with all vertices if necessary
   vector<vector<int> > faces_all_verts;
   if (!merge_verts) {
      faces_all_verts = faces;
      remap_faces(faces_all_verts, vm_all_verts);
   }

   // use geom's face's and not a copy for merged vertices
   if (vm_merged_verts.size())
      remap_faces(faces, vm_merged_verts);

   // load face sort vector
   for(unsigned int i=0;i<faces.size();i++) {
      col_val col;
      if (cg) {
         if (elem=='f')
            col = cg->get_f_col(i);
         else
            col = cg->get_e_col(i);
      }

      // use the reverse flag from the faces with merged vertices
      bool reversed = polygon_sort(faces[i]);
      // if using faces with all verts mapped, use different declaration
      if (!merge_verts) {
         polygon_sort(faces_all_verts[i]);
         fs.push_back(facesSort(i, faces[i], faces_all_verts[i], col, reversed));
      }
      else
         fs.push_back(facesSort(i, faces[i], col, reversed));
   }
   faces_all_verts.clear();

   // clear some memory
   faces.clear();
   if (cg) {
      if (elem=='f')
         cg->clear_f_cols();
      else
         cg->clear_e_cols();
   }

   // sort on faces with merged vertices
   stable_sort(fs.begin(), fs.end(), cmp_faces);

   // check to see we are actually deleting elements
   // mark coincident faces for skipping if any
   bool deleting_faces = ( (elem=='e' && strchr(delete_elems.c_str(), 'e')) || (elem=='f' && strchr(delete_elems.c_str(), 'f')) );
   if (deleting_faces) {
      if (equiv_elems)
         fs[0].equiv_faces.push_back(fs[0].face_no);

      int cur_undeleted = 0;  // the first face in a set of equivalent faces
      unsigned int j = 0;
      for(unsigned int i=0;i<fs.size()-1;i++) {
         j=i+1;
         if (fs[i].face == fs[j].face)
            fs[j].deleted = true;
         else {
            if (!equiv_elems)
               fs[cur_undeleted].average_col = average_face_color(fs, elem, cur_undeleted, i, blend_type);
            cur_undeleted = j;
         }

         if (equiv_elems)
            fs[cur_undeleted].equiv_faces.push_back(fs[j].face_no);
      }
      if (!equiv_elems)
         fs[cur_undeleted].average_col = average_face_color(fs, elem, cur_undeleted, j, blend_type);
   }

   // if all vertices use/re-sort by those faces
   // but deleted flag will be preserved
   if (!merge_verts) {
      for(unsigned int i=0;i<fs.size();i++)
         fs[i].face = fs[i].face_all_verts;
      stable_sort(fs.begin(), fs.end(), cmp_faces);
   }

   // restore original sort of faces unless sort only or congruency check
   if (!sort_only && !equiv_elems)
      sort(fs.begin(), fs.end(), cmp_face_no);

   // write out sorted faces and colors
   for(unsigned i=0; i<fs.size(); i++) {
      if (!fs[i].deleted) {
        // restore orientation unless sort only
        if (!sort_only && fs[i].reversed)
            reverse(fs[i].face.begin(), fs[i].face.end());
         faces.push_back(fs[i].face);

         if (equiv_elems)
            (*equiv_elems)[faces.size()-1].insert(fs[i].equiv_faces.begin(), fs[i].equiv_faces.end());
         else
         if (cg) {
            if (deleting_faces) {
               if (elem=='f')
                  cg->set_f_col(i,fs[i].average_col);
               else
                  cg->set_e_col(i,fs[i].average_col);
            }
            else {
               if (elem=='f')
                  cg->set_f_col(i,fs[i].col);
               else
                  cg->set_e_col(i,fs[i].col);
            }
         }
      }
   }
}

class vertSort
{
public:
   int vert_no;
   vec3d vert;
   col_val col;
   col_val average_col;
   bool deleted;

   int vert_no_new;

   vertSort(int i, vec3d v, col_val c) : vert_no(i), vert(v), col(c) { deleted = false; }
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

class vertSortPostMerge
{
public:
   int vert_no;
   vec3d vert;
   col_val col;
   int vert_new;
   vertSortPostMerge(int i, vec3d v, col_val c, int n) :
      vert_no(i), vert(v), col(c), vert_new(n) {}
};

bool cmp_vert_no(const vertSortPostMerge &a, const vertSortPostMerge &b)
{
   return a.vert_no < b.vert_no;
}

col_val average_vert_color(const vector<vertSort> &vs, const int &begin, const int &end, const int &blend_type)
{
   // quick decision, if only one instance, return its own color
   if (!(end-begin))
      return vs[begin].col;

   // collect colors
   vector<col_val> cols;  
   for(int i=begin;i<=end;i++)
      cols.push_back(vs[i].col);
   
   return (cols.size() ? average_color(cols, blend_type) : col_val());
}

void sort_vertices(geom_if &geom, vector<vertexMap> &vm_all_verts, vector<vertexMap > &vm_merged_verts, 
                   const string &delete_elems,
                   map<int, set<int> > *equiv_elems = 0,
                   int blend_type = 1, double eps = epsilon)
{
   vector<vec3d> &verts = geom.raw_verts();
   col_geom *cg = dynamic_cast<col_geom *>(&geom);

   // sort only when 's' is set
   bool sort_only = strchr(delete_elems.c_str(), 's');
   bool merge_verts = strchr(delete_elems.c_str(), 'v');

   vector<vertSort> vs;

   // load vertex sort vector
   for(unsigned int i=0;i<verts.size();i++) {
      col_val col;
      if (cg)
         col = cg->get_v_col(i);
      vs.push_back(vertSort(i,verts[i],col));
   }

   // clear some memory
   verts.clear();
   if (cg)
      cg->clear_v_cols();

   stable_sort(vs.begin(), vs.end(), vert_cmp(eps));

   // if not merging vertices, build map from old to new vertices for all vertices
   if (!merge_verts) {
      for(unsigned int i=0;i<vs.size();i++)
         vm_all_verts.push_back(vertexMap(vs[i].vert_no,i));
      stable_sort(vm_all_verts.begin(), vm_all_verts.end(), cmp_vertexMap);
   }

   // if merging vertices, build map for merged vertices
   // mark coincident vertices for skipping if any
   // this is always done
   if (true) {
      int v = 0;
      int cur_undeleted = 0;  // the first face in a set of equivalent verts
      unsigned int j = 0;
      for(unsigned int i=0;i<vs.size()-1;i++) {
         j=i+1;

         vm_merged_verts.push_back(vertexMap(vs[i].vert_no,v));

         if (!compare(vs[i].vert, vs[j].vert, eps)) {
            // don't set delete flag if we will not be using it
            if (merge_verts)
               vs[j].deleted = true;
         }
         else {
            v++;
            if (cg && !equiv_elems)
               vs[cur_undeleted].average_col = average_vert_color(vs, cur_undeleted, i, blend_type);
            cur_undeleted = j;
         }
      }
      if (cg && !equiv_elems)
         vs[cur_undeleted].average_col = average_vert_color(vs, cur_undeleted, j, blend_type);

      vm_merged_verts.push_back(vertexMap(vs.back().vert_no,v));

      stable_sort(vm_merged_verts.begin(), vm_merged_verts.end(), cmp_vertexMap);
   }

   // the vertices to be written out have to be put into a second structure
   // so that no deleted vertices (if any) exist in the list
   vector<vertSortPostMerge> vspm;

   for(unsigned int i=0;i<vs.size();i++) {
      if ( !vs[i].deleted ) {
         col_val col;
         if(cg && !equiv_elems) {
            if (merge_verts)
               col = vs[i].average_col;
            else
               col = cg->get_v_col(vs[i].vert_no);
         }

         vspm.push_back(vertSortPostMerge(vs[i].vert_no, vs[i].vert, col, (int)vspm.size()));
      }

      if (equiv_elems)
         (*equiv_elems)[vspm.size()-1].insert(vs[i].vert_no);
   }

   // restore original sort of vertices unless sort only or congruency check
   if (!sort_only && !equiv_elems) {
      sort(vspm.begin(), vspm.end(), cmp_vert_no);

      // adjust the vertex maps
      for(unsigned i=0; i<vm_merged_verts.size(); i++) {
         for(unsigned j=0; j<vspm.size(); j++) {
            if (vspm[j].vert_new == vm_merged_verts[i].new_vertex) {
               vm_merged_verts[i].new_vertex = j;
               break;
            }
         }
      }

      for(unsigned i=0; i<vm_all_verts.size(); i++)
         vm_all_verts[i].new_vertex = vm_all_verts[i].old_vertex;
   }

   // write out sorted vertices and colors
   for(unsigned i=0; i<vspm.size(); i++) {
      verts.push_back(vspm[i].vert);
      if(cg && !equiv_elems)
         cg->set_v_col(i,vspm[i].col);
   }
}

bool sort_merge_elems(geom_if &geom, const string &merge_elems, vector<map<int, set<int> > > *equiv_elems, bool chk_congruence, int blend_type, double eps)
{
   // an empty geom cannot be processed
   if (!geom.verts().size())
      return false;

   // if there is nothing to merge return
   if (!merge_elems.size())
      return false;
      
   // Congruence check expects a polyhedron with elements occurring
   // in coincident pairs, and exits early if this is not the case
   vector<map<int, set<int> > > equiv_elems_tmp;
   if (chk_congruence && !equiv_elems)
      equiv_elems = &equiv_elems_tmp;

   if (equiv_elems) {
      equiv_elems->clear();
      equiv_elems->resize(3);
   }

   unsigned int num_verts = geom.verts().size();
   unsigned int num_edges = geom.edges().size();
   unsigned int num_faces = geom.faces().size();

   vector<vertexMap> vm_all_verts, vm_merged_verts;
   sort_vertices(geom, vm_all_verts, vm_merged_verts, merge_elems, (equiv_elems ? &(*equiv_elems)[0] : 0), blend_type, eps);
   if (chk_congruence && (*equiv_elems)[0].size()*2 != num_verts)
      return false;

   if (geom.edges().size())
      sort_faces(geom, vm_all_verts, vm_merged_verts, merge_elems, 'e', (equiv_elems ? &(*equiv_elems)[1] : 0), blend_type);
   if (chk_congruence && (*equiv_elems)[1].size()*2 != num_edges)
      return false;
   
   if (geom.faces().size())
      sort_faces(geom, vm_all_verts, vm_merged_verts, merge_elems, 'f', (equiv_elems ? &(*equiv_elems)[2] : 0), blend_type);
   if (chk_congruence && (*equiv_elems)[2].size()*2 != num_faces)
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
         if (to>=cnts[i])
            to -= cnts[i];
         elem_maps[i][to] = to;
         set<int>::iterator si;
         for(si=mi->second.begin(); si!=mi->second.end(); ++si) {
            int from = (*si<cnts[i]) ? *si : *si-cnts[i];
            if (to != from) {
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


