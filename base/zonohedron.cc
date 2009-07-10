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
   Name: zono.cc
   Description: creating zonohedra
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
#include "utils.h"
#include "polygons.h"

using std::vector;
using std::map;
using std::pair;
using std::swap;
using std::logical_not;
using std::string;


struct vec_less {
   bool operator ()(const vec3d &v1, const vec3d &v2)
   { return compare(v1, v2, 1e-10)<0; }
};


// Normalised direction of the line of a vector
vec3d normalised_dir(const vec3d v) {
   if(v[0]<-epsilon)
      return -v;
   else if(v[0]<epsilon) {
      if(v[1]<-epsilon)
         return -v;
      else if(v[1]<epsilon) {
         if(v[2]<-epsilon)
            return -v;
      }
   }
   return v;
}


void get_star(const geom_if &geom, vector<vec3d> &star, char type,
      vec3d centre)
{
   star.clear();
   const vector<vec3d> &verts = geom.verts();
   switch(type) {
      case 'v':      // vertices
         for(unsigned int i=0; i<verts.size(); i++) {
            vec3d v = verts[i] - centre;
            if(v.mag2()>epsilon*epsilon)
               star.push_back(v);
         }
         break;
      case 'a':      // any vector between two vertices
         for(unsigned int i=0; i<verts.size(); i++)
            for(unsigned int j=i+1; j<verts.size(); j++)
               star.push_back(geom.edge_vec(i,j));
         break;
      case 'e':      // explicit edge vectors
         for(unsigned int i=0; i<geom.edges().size(); i++)
            star.push_back(geom.edge_vec(i));
         break;
      case 'i': {    // implict edge vectors
         geom_info info(geom);
         const vector<vector<int> > edges = info.get_impl_edges();
         for(unsigned int i=0; i<edges.size(); i++)
            star.push_back(geom.edge_vec(edges[i]));
         break;
      }
   }
}


void make_zono_1d(geom_if &zono, const vector<vec3d> &star)
{
   zono.clear_all();
   if(!star.size())
      return;
   vector<vec3d> &verts = *zono.get_verts();
   vec3d pos = vec3d(0, 0, 0);
   vec3d neg = vec3d(0, 0, 0);
   for(unsigned int i=0; i<star.size(); i++) {
      double cos_a = vdot(star[0], star[i]);
      if(cos_a > 0)
         pos += star[i];
      else
         neg += star[i];
   }
   verts.push_back(pos);
   verts.push_back(neg);
}


void make_zono_2d(geom_if &zono, const vector<vec3d> &star, vec3d fnorm)
{
   zono.clear_all();
   map<vec3d, set<int>, vec_less> stars;
   for(unsigned int i=0; i<star.size(); i++)
      stars[normalised_dir(star[i]).unit()].insert(i);

   vector<vec3d> &verts = *zono.get_verts();
   map<vec3d, set<int> >::iterator mi;
   for(mi=stars.begin(); mi!=stars.end(); ++mi) {
      vector<vec3d> z_face_star;
      vec3d pos = vec3d(0, 0, 0);
      vec3d neg = vec3d(0, 0, 0);
      const vec3d &norm = vcross(mi->first, fnorm).unit();
      for(int i=0; i<(int)star.size(); i++) {
         if(mi->second.find(i) != mi->second.end())
            z_face_star.push_back(star[i]);
         else {
            double cos_a = vdot(norm, star[i]);
            if(cos_a > 0)
               pos += star[i];
            else
               neg += star[i];
         }
      }
      geom_v zono_face;
      make_zono_1d(zono_face, z_face_star);
      vector<vec3d> &zf_verts = *zono_face.get_verts();
      for(unsigned int k=0; k<zf_verts.size(); k++) {
         verts.push_back(pos + zf_verts[k]);
         verts.push_back(neg + zf_verts[k]);
      }
   }

}


bool make_zono(geom_if &zono, const vector<vec3d> &star, char *errmsg)
{
   zono.clear_all();
   
   vector<double> star_mags(star.size());
   for(unsigned int i=0; i< star.size(); i++)
      star_mags[i] = star[i].mag();

   map<vec3d, set<int>, vec_less > stars2d;
   for(unsigned int i=0; i<star.size()-1; i++)
      for(unsigned int j=i+1; j<star.size(); j++) {
         vec3d n = vcross(star[i], star[j])/(star_mags[i]*star_mags[j]);
         if(n.mag2()>epsilon*epsilon) {
            vec3d norm = normalised_dir(n).unit();
            stars2d[norm].insert(i);
            stars2d[norm].insert(j);
         }
      }

   vector<vec3d> &verts = *zono.get_verts();
   map<vec3d, set<int> >::iterator mi;
   for(mi=stars2d.begin(); mi!=stars2d.end(); ++mi) {
      vector<vec3d> z_face_star;
      vec3d pos = vec3d(0, 0, 0);
      vec3d neg = vec3d(0, 0, 0);
      const vec3d &norm = mi->first;
      for(int i=0; i<(int)star.size(); i++) {
         if(mi->second.find(i) != mi->second.end())
            z_face_star.push_back(star[i]);
         else {
            double cos_a = vdot(norm, star[i]);
            if(cos_a > 0)
               pos += star[i];
            else
               neg += star[i];
         }
      }
      
      geom_v zono_face;
      make_zono_2d(zono_face, z_face_star, norm);
      vector<vec3d> &zf_verts = *zono_face.get_verts();
      for(unsigned int k=0; k<zf_verts.size(); k++) {
         verts.push_back(pos + zf_verts[k]);
         verts.push_back(neg + zf_verts[k]);
      }
   }
   
   return zono.set_hull("A0.9999999", errmsg);
}


/*

void make_zono_1d(geom_if &zono, const vector<vec3d> &star)
{
   zono.clear_all();
   vector<double> star_mags(star.size());
   for(unsigned int i=0; i< star.size(); i++)
      star_mags[i] = star[i].mag();
   
   vector<vec3d> &verts = *zono.get_verts();
   for(unsigned int i=0; i<star.size(); i++) {
      vec3d pos = vec3d(0, 0, 0);
      vec3d neg = vec3d(0, 0, 0);
      for(unsigned int k=0; k<star.size(); k++) {
         double cos_a = vdot(star[0], star[k])/star_mags[k];
         if(cos_a > 0)
            pos += star[k];
         else
            neg += star[k];
      }
      
      //if(pos==vec3d(0, 0, 0))
      //   continue;
      verts.push_back(pos);
      verts.push_back(neg);
   }
}

void make_zono_2d(geom_if &zono, const vector<vec3d> &star, vec3d fnorm)
{
   zono.clear_all();
   vector<double> star_mags(star.size());
   for(unsigned int i=0; i< star.size(); i++)
      star_mags[i] = star[i].mag();
  
   vector<bool> seen(star.size(), false);
   vector<vec3d> &verts = *zono.get_verts();
   for(unsigned int i=0; i<star.size(); i++) {
      if(seen[i])
         continue;
      vec3d pos = vec3d(0, 0, 0);
      vec3d neg = vec3d(0, 0, 0);
      vec3d norm = vcross(fnorm, star[i]).unit();
      vector<vec3d> perp;
      for(unsigned int k=i; k<star.size(); k++) {
         double cos_a = vdot(norm, star[k])/star_mags[k];
         //if(fabs(cos_a) < epsilon && star[k] != vec3d(0, 0, 0))
         if(fabs(cos_a) < epsilon) {
            seen[k] = true;
            perp.push_back(star[k]);
         }
         else if(cos_a > 0)
            pos += star[k];
         else
            neg += star[k];
      }
      
      //verts.push_back(pos);
      //verts.push_back(neg);
      if(perp.size()) {
         geom_v zono_face;
         make_zono_1d(zono_face, perp);
         vector<vec3d> &zf_verts = *zono_face.get_verts();
         for(unsigned int k=0; k<zf_verts.size(); k++) {
            verts.push_back(pos + zf_verts[k]);
            verts.push_back(neg + zf_verts[k]);
         }
      }
   }
}


bool make_zono(geom_if &zono, const vector<vec3d> &star, char *errmsg)
{
   zono.clear_all();
   vector<double> star_mags(star.size());
   for(unsigned int i=0; i< star.size(); i++)
      star_mags[i] = star[i].mag();
   
   vector<vec3d> &verts = *zono.get_verts();
   for(unsigned int i=0; i<star.size()-1; i++) {
      for(unsigned int j=i+1; j<star.size(); j++) {
         //fprintf(stderr, "face no. %d\n", i*(star.size()-1)+j - (j>i) );
         vec3d pos = vec3d(0, 0, 0);
         vec3d neg = vec3d(0, 0, 0);
         vec3d norm = vcross(star[i], star[j]).unit();
         vector<vec3d> perp;
         for(unsigned int k=0; k<star.size(); k++) {
            double cos_a = vdot(norm, star[k])/star_mags[k];
            //if(fabs(cos_a) < epsilon && star[k] != vec3d(0, 0, 0))
            if(fabs(cos_a) < epsilon)
               perp.push_back(star[k]);
            else if(cos_a > 0)
               pos += star[k];
            else
               neg += star[k];
         }
      
         //verts.push_back(pos);
         //verts.push_back(neg);
         if(perp.size()) {
            geom_v zono_face;
            make_zono_2d(zono_face, perp, norm);
            vector<vec3d> &zf_verts = *zono_face.get_verts();
            for(unsigned int k=0; k<zf_verts.size(); k++) {
               verts.push_back(pos + zf_verts[k]);
               verts.push_back(neg + zf_verts[k]);
            }
         }
      }
   }
   
   return 1; //zono.set_hull("A0.9999999", errmsg);
}
*/

void make_polar_zono(geom_if &zono, int star_n, bool out_star)
{
   zono.clear_all();
   geom_v star;
   if(is_even(star_n))
      uni_pgon(star, prism(star_n));
   else
      uni_pgon(star, antiprism(star_n));
      
   if(out_star)
      zono = star;
   else
      zono.set_zono(star.get_star());
}


