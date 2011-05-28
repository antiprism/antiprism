/*
   Copyright (c) 2003-2007, Adrian Rossiter, Roger Kaufman
   Includes ideas and algorithms by George W. Hart, http://www.georgehart.com

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
   Name: canonic.cc
   Description: canonicalize a polyhedron
                Implementation of George Hart's canonicalization algorithm
                http://library.wolfram.com/infocenter/Articles/2012/
   Project: Antiprism - http://www.antiprism.com
*/

#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "geom.h"
#include "info.h"

using std::string;
using std::vector;


// Implementation of George Hart's canonicalization algorithm
// http://library.wolfram.com/infocenter/Articles/2012/
void canonicalize_mm(geom_if &geom, double edge_factor, double plane_factor, int n,
                     int divergence_test, int rep_count, bool planar_only, double eps)
{
   // do a scale to get edges close to 1
   geom_info info(geom);
   double scale = info.iedge_lengths().sum/info.num_iedges();
   if (scale)
      geom.transform(mat3d::scale(1/scale));

   const vector<vec3d> &verts = geom.verts();
   const vector<vector<int> > &faces = geom.faces();
   vector<vector<int> > edges;
   geom.get_impl_edges(edges);

   geom_info rep(geom);
   rep.set_center(geom.centroid());
   double starting_radius = rep.vert_dists().max;

   double max_diff2=0;
   vec3d origin(0,0,0);
   unsigned int cnt;
   for(cnt=0; cnt<(unsigned int)n;) {
      vector<vec3d> old_verts = verts;

      if (!planar_only) {
         vector<vec3d> near_pts;
         for(unsigned int e=0; e<edges.size(); e++) {
            vec3d P = geom.edge_nearpt(edges[e], origin);
            near_pts.push_back(P);
            vec3d offset = edge_factor*(P.mag()-1)*P;
            geom.raw_verts()[edges[e][0]] -= offset;
            geom.raw_verts()[edges[e][1]] -= offset;
         }
         
         vec3d cent_near_pts = centroid(near_pts);
         for(unsigned int i=0; i<verts.size(); i++)
            geom.raw_verts()[i] -= cent_near_pts;
      }

      // Make a copy of verts. zero out.
      vector<vec3d> vs = verts;
      for(unsigned int i=0; i<vs.size(); i++)
         vs[i] = vec3d(0,0,0);

      for(unsigned int ff=cnt; ff<faces.size()+cnt; ff++) {
         int f = ff%faces.size();
         if(faces[f].size()==3)
            continue;
         vec3d norm = geom.face_norm(f).unit();
         vec3d f_cent = geom.face_cent(f);
         if(vdot(norm, f_cent)<0)
            norm *= -1.0;
         for(unsigned int v=0; v<faces[f].size(); v++)
            vs[faces[f][v]] +=
               vdot(plane_factor*norm, f_cent - verts[faces[f][v]])*norm;
      }

      // adjust vertices post-loop
      for(unsigned int i=0; i<vs.size(); i++)
         geom.raw_verts()[i] += vs[i];
      
      max_diff2 = 0;
      for(unsigned int i=0; i<verts.size(); i++) {
         double diff2 = (verts[i] - old_verts[i]).mag2();
         if(diff2>max_diff2)
            max_diff2 = diff2;
      }

      // increment count here for reporting
      cnt++;

      if(rep_count > 0 && (cnt)%rep_count == 0)
         fprintf(stderr, "%-15d max_diff=%12.10g\n", cnt, sqrt(max_diff2));

      if(sqrt(max_diff2)<eps)
         break;

      // see if radius is expanding or contracting unreasonably
      if (divergence_test > 0) {
         rep.set_center(geom.centroid());
         if((rep.vert_dists().max > starting_radius * divergence_test) ||
            (rep.vert_dists().max < starting_radius / divergence_test)) {
            fprintf(stderr,"Probably Diverging. Breaking out.\n");
            break;
         }
      }
   }
   if(rep_count > -1) {
      fprintf(stderr, "\n%-15d final max_diff=%12.10g\n", cnt, sqrt(max_diff2));
      fprintf(stderr, "\n");
   }
}
            
vector<vec3d> reciprocalN(geom_if &geom)
{
   const vector<vector<int> > &faces = geom.faces();
   const vector<vec3d> &verts = geom.verts();

   vector<vec3d> normals;
   for(unsigned int i=0;i<faces.size();i++) {
      vec3d centroid(0,0,0);
      vec3d normal(0,0,0);
      double avgEdgeDist = 0;

      int v1 = faces[i].at(faces[i].size()-2);
      int v2 = faces[i].at(faces[i].size()-1);
      for(unsigned int j=0;j<faces[i].size();j++) {
         int v3 = faces[i].at(j);
         centroid += verts[v3];
         //normal += orthogonal(verts[v1], verts[v2], verts[v3]);
         normal += vcross(verts[v3]-verts[v2], verts[v2]-verts[v1]); // orthogonal without call
         //avgEdgeDist += tangentPoint(verts[v1], verts[v2]).mag();
         vec3d d = verts[v2] - verts[v1];
         // prevent division by zero
         //avgEdgeDist += (verts[v1] - ((vdot(d,verts[v1])/d.mag2()) * d)).mag(); // tangentPoint without call
         double vdt;
         if (d[0] == 0 && d[1] == 0 && d[2] == 0)
            vdt = 0;
         else
            vdt = vdot(d,verts[v1])/d.mag2();
         avgEdgeDist += (verts[v1] - (vdt * d)).mag(); // tangentPoint without call
         v1 = v2;
         v2 = v3;
      }
      centroid *= 1.0/faces[i].size();
      normal.to_unit();
      avgEdgeDist /= faces[i].size();

      // reciprocal call replace below:
      //vec3d ans = reciprocal(normal * vdot(centroid,normal));
      vec3d v = normal * vdot(centroid,normal);
      // prevent division by zero
      vec3d ans;
      if (v[0] == 0 && v[1] == 0 && v[2] == 0)
         ans = v;
      else {
         ans = v * 1.0/v.mag2();
         ans *= (1+avgEdgeDist)/2;
      }
      normals.push_back(ans);
   }

   return normals;
}

vector<vec3d> reciprocalC(geom_if &geom)
{
   vector<vec3d> centers;
   geom.face_cents(centers);
   for(unsigned int i=0; i<centers.size(); i++)
      centers[i] /= centers[i].mag2();
   return centers;
}

vector<vec3d> reciprocalC_mag(geom_if &geom)
{
   vector<vec3d> centers;
   geom.face_cents(centers);
   for(unsigned int i=0; i<centers.size(); i++)
      centers[i] /= centers[i].mag();
   return centers;
}

// Addition to algorithm by Adrian Rossiter. Finds the correct centroid for the canonical
vec3d edge_nearpoints_centroid(geom_if &geom, vec3d cent)
{
   vector<vector<int> > edges;
   geom.get_impl_edges(edges);
   int e_sz = edges.size();
   vec3d e_cent(0,0,0);
   for(int e=0; e<e_sz; ++e)
      e_cent += geom.edge_nearpt(edges[e], cent);
   return e_cent / double(e_sz);
}

// Implementation of George Hart's planarization and canonicalization algorithms
// http://www.georgehart.com/virtual-polyhedra/conway_notation.html
void canonicalize_cn(geom_if &geom, int n, char method, int divergence_test, int rep_count, double eps)
{
   // do a scale to get edges close to 1
   geom_info info(geom);
   double scale = info.iedge_lengths().sum/info.num_iedges();
   if (scale)
      geom.transform(mat3d::scale(1/scale));

   geom_v dual;
   get_dual(geom, dual, 0);
   const vector<vec3d> &verts = geom.verts();
   const vector<vec3d> &d_verts = dual.verts();

   geom_info rep(geom);
   rep.set_center(geom.centroid());
   double starting_radius = rep.vert_dists().max;

   double max_diff2=0;
   unsigned int cnt;
   for(cnt=0; cnt<(unsigned int)n;) {
      vector<vec3d> old_verts = verts;
      vec3d e_cent;

      switch(method) {
         // base/dual canonicalize method
         case 'n':
            dual.raw_verts() = reciprocalN(geom);
            geom.raw_verts() = reciprocalN(dual);
            e_cent = edge_nearpoints_centroid(geom, vec3d(0,0,0));
            geom.transform(mat3d::transl(-0.1*e_cent));
            break;

         // adjust vertices with side effect of planarization. mag2() version
         case 'p':
            // move centroid to origin for balance
            dual.raw_verts() = reciprocalC(geom);
            geom.transform(mat3d::transl(-centroid(d_verts)));
            geom.raw_verts() = reciprocalC(dual);
            geom.transform(mat3d::transl(-centroid(verts)));
            break;

         // adjust vertices with side effect of planarization. mag() version
         case 'q':
            // move centroid to origin for balance
            dual.raw_verts() = reciprocalC_mag(geom);
            geom.transform(mat3d::transl(-centroid(d_verts)));
            geom.raw_verts() = reciprocalC_mag(dual);
            geom.transform(mat3d::transl(-centroid(verts)));
            break;

         case 'x':
            geom.face_cents(dual.raw_verts());
            dual.face_cents(geom.raw_verts());
            break;
         }
      
      max_diff2 = 0;
      for(unsigned int i=0; i<verts.size(); i++) {
         double diff2 = (verts[i] - old_verts[i]).mag2();
         if(diff2>max_diff2)
            max_diff2 = diff2;
      }

      // increment count here for reporting
      cnt++;

      if(rep_count > 0 && (cnt)%rep_count == 0)
         fprintf(stderr, "%-15d max_diff=%12.10g\n", cnt, sqrt(max_diff2));

      if(sqrt(max_diff2)<eps)
         break;

      // see if radius is expanding or contracting unreasonably
      if (divergence_test > 0) {
         rep.set_center(geom.centroid());
         if((rep.vert_dists().max > starting_radius * divergence_test) ||
            (rep.vert_dists().min < starting_radius / divergence_test)) {
            fprintf(stderr,"Probably Diverging. Breaking out.\n");
            break;
         }
      }
   }

   if(rep_count > -1) {
      fprintf(stderr, "\n%-15d final max_diff=%12.10g\n", cnt, sqrt(max_diff2));
      fprintf(stderr, "\n");
   }
}
