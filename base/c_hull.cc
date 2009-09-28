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

/* \file c_hull.cc
   \brief wrapper around qhull for convex hulls and delaunay edges.
*/

#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <functional>
#include <vector>
#include <map>
#include <string>

#include "geom.h"
#include "geom_utils.h"
#include "utils.h"

#include "qh_qhull_a.h"

#define NUM_TO_USE_DELAUNAY 20

using std::vector;
using std::map;
using std::pair;
using std::swap;
using std::string;



bool make_hull(geom_if &geom, bool append, string qh_args, char *errmsg)
{
   vector<vec3d> pts = geom.verts();
   
   map<int, col_val> vcols;
   col_geom *cg = dynamic_cast<col_geom *>(&geom);
   if(cg)
      vcols = cg->vert_cols();

   if(!append)
      geom.clear_all();
   
   const int dim=3;
   coordT *points = new coordT[pts.size()*dim];
   
   for(unsigned i=0; i<pts.size(); i++)
      for(int j=0; j<dim; j++)
         points[i*dim+j] = pts[i][j];

   qh_args.insert(0, "qhull o ");

   boolT ismalloc= False;    // don't free points in qh_freeqhull() or realloc
   FILE *outfile= NULL;      // suppress output from qh_produce_output()
   FILE *errfile = fopen("/dev/null", "w");    // suppress qhull error messages
   if(!errfile)              // must be a valid pointer
      errfile = stderr;

   int ret=qh_new_qhull(dim, pts.size(), points, ismalloc, (char *)qh_args.c_str(), outfile, errfile); 

   if(ret) {
      if(errmsg)
         snprintf(errmsg, MSG_SZ, "error calculating convex hull");
      return false;
   }
  
   map <size_t, int> vert_order;
   vertexT *vertex;
   
   if(!append) {
      int i=0;
      FORALLvertices {
         size_t idx = (vertex->point-points)/dim;
         vert_order[idx] = i++;
         int v_idx = geom.add_vert(pts[idx]);
         if(cg)
            cg->set_v_col(v_idx, vcols[idx]);
      }
   }

   facetT *facet;
   FORALLfacets {
      vector<int> face, ordered_face;
      vertexT *vid; //, **vidp;
      int vid_i, vid_n;
      FOREACHsetelement_i_(vertexT, facet->vertices, vid) {
         if(!append)
            face.push_back(vert_order[(vid->point-points)/dim]);
         else
            face.push_back((vid->point-points)/dim);
      }
      if(face.size()>3) {
         // order the face vertices so joining them sequentially will
         // form the polygon
         unsigned int r_cnt=0;
         vector<int> lns;
         ridgeT *ridge; //, **ridgep;
         int ridge_i, ridge_n;
         FOREACHsetelement_i_(ridgeT, facet->ridges, ridge) {
            r_cnt++; 
            FOREACHsetelement_i_(vertexT, ridge->vertices, vid) {
               if(!append)
                  lns.push_back(vert_order[(vid->point-points)/dim]);
               else
                  lns.push_back((vid->point-points)/dim);
            }
         }
         if(lns.size()!=2*r_cnt) {
            snprintf(errmsg, MSG_SZ, "ridge included more than two points");
         }

         ordered_face.push_back(lns[0]);
         int pt = lns[1];
         for(unsigned int j=1; j<r_cnt; j++) {
            ordered_face.push_back(pt);
            for(unsigned int k=j; k<r_cnt; k++) {
               if(lns[k*2]==pt) {
                  pt=lns[k*2+1];
                  swap(lns[j*2], lns[k*2+1]);
                  swap(lns[j*2+1], lns[k*2]);
                  break;
               }
               else if(lns[k*2+1]==pt) {
                  pt=lns[k*2];
                  swap(lns[j*2], lns[k*2]);
                  swap(lns[j*2+1], lns[k*2+1]);
                  break;
               }
            }
         }
      }
      else {
         ordered_face = face;
      }
      geom.add_face(ordered_face);
   }
   
   // clean up
   qh_freeqhull(!qh_ALL);                 // free long memory
   int curlong, totlong;
   qh_memfreeshort (&curlong, &totlong);  // free short mem and mem allocator
  
   delete[] points;

   return true;
}


bool add_hull(geom_if &geom, string qh_args, char *errmsg)
{
   return make_hull(geom, true, qh_args, errmsg);
}

int set_hull(geom_if &geom, string qh_args, char *errmsg)
{
   int ret = make_hull(geom, false, qh_args, errmsg);
   if(!ret)
      geom.clear_all();
   return ret;
}

      
        
         

bool get_delaunay_edges(const geom_if &geom, map<pair<int, int>, int> &edges,
      string qh_args, char *errmsg)
{
   vector<vec3d> pts = geom.verts();
   
   const int dim=3;
   coordT *points = new coordT[pts.size()*dim];
   
   for(unsigned i=0; i<pts.size(); i++)
      for(int j=0; j<dim; j++)
         points[i*dim+j] = pts[i][j];

   qh_args.insert(0, "qhull d Qbb QJ o ");

   boolT ismalloc= False;    // don't free points in qh_freeqhull() or realloc
   FILE *outfile= NULL;      // suppress output from qh_produce_output()
   FILE *errfile = fopen("/dev/null", "w");    // suppress qhull error messages
   if(!errfile)              // must be a valid pointer
      errfile = stderr;

   int ret=qh_new_qhull(dim, pts.size(), points, ismalloc, (char *)qh_args.c_str(), outfile, errfile); 

   if(ret) {
      if(errmsg)
         snprintf(errmsg, MSG_SZ, "error calculating delaunay triangulation");
      return false;
   }
  
   
   edges.clear();
   facetT *facet;
   FORALLfacets {
      if(!facet->good || facet->upperdelaunay)
         continue;
      vector<int> tet;
      vertexT *vid;
      int vid_i, vid_n;
      FOREACHsetelement_i_(vertexT, facet->vertices, vid) {
         tet.push_back((vid->point - qh first_point )/(dim+1));
      }
      
      if(tet.size()<4)
         continue;
      
      bool coplanar=true;
      vec3d v1=pts[tet[1]]-pts[tet[0]];
      vec3d v2=pts[tet[2]]-pts[tet[0]];
      for(unsigned int i=3; i<tet.size(); i++)
         if(fabs(vtriple(v1,v2,pts[tet[i]]-pts[tet[0]]))>epsilon) {
            coplanar=false;
            break;
         }
      if(coplanar)
         continue;

      sort(tet.begin(), tet.end());
      for(unsigned int i=0; i<tet.size()-1; i++)
         for(unsigned int j=i+1; j<tet.size(); j++) {
            edges[pair<int, int>(tet[i], tet[j])]++;
         }
   }
   
   // clean up
   qh_freeqhull(!qh_ALL);                 // free long memory
   int curlong, totlong;
   qh_memfreeshort (&curlong, &totlong);  // free short mem and mem allocator
  
   delete[] points;

   return true;
}

double minimum_distance(const geom_if &geom, double sig_dist, char *errmsg)
{
   double sig_dist2 = sig_dist*sig_dist;
   const vector<vec3d> &verts = geom.verts();
   //double min_dist = 1e100;
   double min_sig_dist = 1e100;
   double mag2;
   
   if(!verts.size())
      min_sig_dist = 0;
   else if(verts.size() < NUM_TO_USE_DELAUNAY) {
      vec3d dvec;
      for(unsigned int i=0; i<verts.size()-1; i++)
         for(unsigned int j=i+1; j<verts.size(); j++) {
            mag2 = (verts[j] - verts[i]).mag2();
            if(mag2>=sig_dist2 && mag2<min_sig_dist)
               min_sig_dist = mag2;
         }
   }
   else {
      map<pair<int, int>, int> edges;
      if(!get_delaunay_edges(geom, edges, "", errmsg))
         return -1;

      map<pair<int, int>, int>::iterator mi;
      for(mi=edges.begin(); mi!=edges.end(); ++mi) {
         mag2 = (verts[mi->first.second] -  verts[mi->first.first]).mag2();
         if(mag2>=sig_dist2 && mag2<min_sig_dist)
            min_sig_dist = mag2;
      }
   }
   return min_sig_dist<1e99 ? sqrt(min_sig_dist) : sig_dist;
}




int get_voronoi_cells(geom_if &geom, vector<col_geom_v> &cells,
      string qh_args, char *errmsg)
{
   vector<vec3d> pts = *geom.get_verts();

   const int dim=3;
   coordT *points = new coordT[pts.size()*dim];

   for(unsigned i=0; i<pts.size(); i++)
      for(int j=0; j<dim; j++)
         points[i*dim+j] = pts[i][j];

   qh_args.insert(0, "qhull v o ");

   boolT ismalloc= False;    /* don't free points in qh_freeqhull() or realloc*/
   FILE *outfile= NULL;    /* output from qh_produce_output()   NULL ?*/
   FILE *errfile= stderr;    /* error messages from qhull code */
   //FILE *errfile= NULL;    /* error messages from qhull code */

   int ret=qh_new_qhull(dim, pts.size(), points, ismalloc, (char *)qh_args.c_str(), outfile, errfile); 


   if(ret) {
      if(errmsg)
         snprintf(errmsg, MSG_SZ, "error calculating voronoi cell");
      return false;
   }

   //qh_setvoronoi_all();
   qh_clearcenters (qh_ASvoronoi);
   qh_vertexneighbors();
   qh_findgood_all (qh facet_list);
   qh RANDOMdist= False;

   facetT *facetlist = qh facet_list;
   setT *facets = NULL;
   boolT printall = !qh_ALL;

   //qh_printfacets(fp, qh_PRINToff, facetlist, facets, printall);
   int numcenters, numvertices= 0, numneighbors, numinf, vertex_i, vertex_n;
   facetT *facet, *neighbor, **neighborp;
   setT *vertices;
   vertexT *vertex;
   boolT islower;
   unsigned int numfacets= (unsigned int) qh num_facets;

   vertices= qh_markvoronoi (facetlist, facets, printall, &islower, &numcenters);
   FOREACHvertex_i_(vertices) {
      if (vertex) {
         numvertices++;
         numneighbors = numinf = 0;
         FOREACHneighbor_(vertex) {
            if (neighbor->visitid == 0)
               numinf= 1;
            else if (neighbor->visitid < numfacets)
               numneighbors++;
         }
         if (numinf && !numneighbors) {
            SETelem_(vertices, vertex_i)= NULL;
            numvertices--;
         }
      }
   }

   col_geom_v vcells;
   FORALLfacet_(facetlist) {
      if (facet->visitid && facet->visitid < numfacets) {
         if (!facet->normal || !facet->upperdelaunay || !qh ATinfinity) {
            if (!facet->center)
               facet->center = qh_facetcenter (facet->vertices);
            vcells.add_vert(vec3d(facet->center[0],
                     facet->center[1], facet->center[2]));
         }
         else
            vcells.add_vert(vec3d(1000,1000,1000));
      }
   }
   FOREACHvertex_i_(vertices) {
      vector<int> face;
      if (vertex) {
         //qh_order_vertexneighbors(vertex);
         qsort (SETaddr_(vertex->neighbors, vertexT), 
               qh_setsize (vertex->neighbors),
               sizeof (facetT *), qh_compare_facetvisit);
      }
      if (vertex) {
         FOREACHneighbor_(vertex) {
            if (neighbor->visitid < numfacets)
               face.push_back(neighbor->visitid-1);
         }
      }
      if(face[0] >= 0)
         vcells.add_face(face);
   }
   qh_settempfree (&vertices);


   // clean up
   qh_freeqhull(!qh_ALL);                 // free long memory
   int curlong, totlong;
   qh_memfreeshort (&curlong, &totlong);  // free short mem and mem allocator

   delete[] points;
   
   const vector<vector<int> > &faces = *vcells.get_faces();
   for(unsigned int i=0; i<faces.size(); i++) {
      col_geom_v cell;
      for(unsigned int j=0; j<faces[i].size(); j++)
         cell.add_vert((*vcells.get_verts())[faces[i][j]]);
      cell.add_hull();
      cells.push_back(cell);
   }
   return true;
}

 
// test points versus hull functions

bool test_points_vs_hull(const vector<vec3d> &P, col_geom_v &hull, bool inside, bool surface, bool outside, double epsilon)
{
   const vector<vec3d> &verts = hull.verts();
   const vector<vector<int> > &faces = hull.faces();

   vec3d C = centroid(verts);

   bool answer = true;
   for(unsigned int i=0;i<faces.size();i++) {
      vec3d n = face_norm(verts, faces[i]).unit();
      double D = vdot(verts[faces[i][0]]-C, n);
      if(D < 0) { // Make sure the normal points outwards
         D = -D;
         n = -n;
      }

      for(unsigned int j=0;j<P.size();j++) {
         double t = vdot(P[j]-C, n);
         if (t < D-epsilon && !inside)
            answer = false;
         else
         if (t > D+epsilon && !outside)
            answer = false;
         else
         if (!surface)
            answer = false;
 
      if(!answer)
         break;
      }

   if(!answer)
      break;
   }

   return answer;
}

bool is_geom_fully_outside_hull(col_geom_v &geom, col_geom_v &hull, double epsilon)
{
   return test_points_vs_hull(geom.verts(),hull,false,false,true,epsilon);
}

bool is_geom_outside_hull(col_geom_v &geom, col_geom_v &hull, double epsilon)
{
   return test_points_vs_hull(geom.verts(),hull,false,true,true,epsilon);
}

bool is_geom_on_surface_hull(col_geom_v &geom, col_geom_v &hull, double epsilon)
{
   return test_points_vs_hull(geom.verts(),hull,false,true,false,epsilon);
}

bool is_geom_inside_hull(col_geom_v &geom, col_geom_v &hull, double epsilon)
{
   return test_points_vs_hull(geom.verts(),hull,true,true,false,epsilon);
}

bool is_geom_fully_inside_hull(col_geom_v &geom, col_geom_v &hull, double epsilon)
{
   return test_points_vs_hull(geom.verts(),hull,true,false,false,epsilon);
}

bool is_point_fully_outside_hull(const vec3d &P, col_geom_v &hull, double epsilon)
{
   col_geom_v tgeom;
   tgeom.add_vert(P);
   return is_geom_fully_outside_hull(tgeom, hull, epsilon);
}

bool is_point_outside_hull(const vec3d &P, col_geom_v &hull, double epsilon)
{
   col_geom_v tgeom;
   tgeom.add_vert(P);
   return is_geom_outside_hull(tgeom, hull, epsilon);
}

bool is_point_on_surface_hull(const vec3d &P, col_geom_v &hull, double epsilon)
{
   col_geom_v tgeom;
   tgeom.add_vert(P);
   return is_geom_on_surface_hull(tgeom, hull, epsilon);
}

bool is_point_inside_hull(const vec3d &P, col_geom_v &hull, double epsilon)
{
   col_geom_v tgeom;
   tgeom.add_vert(P);
   return is_geom_inside_hull(tgeom, hull, epsilon);
}

bool is_point_fully_inside_hull(const vec3d &P, col_geom_v &hull, double epsilon)
{
   col_geom_v tgeom;
   tgeom.add_vert(P);
   return is_geom_fully_inside_hull(tgeom, hull, epsilon);
}

