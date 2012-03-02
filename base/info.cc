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
   Name: info.cc
   Description: information from OFF file
   Project: Antiprism - http://www.antiprism.com
*/

#include <string.h>
#include <set>
#include <vector>
#include <utility>
#include <algorithm>
#include "info.h"
#include "transforms.h"
#include "math_utils.h"
#include "geom_utils.h"

using std::set;
using std::vector;
using std::pair;



int cmp_angles(const double &a, const double &b, double min_diff)
{
   if(fabs(a-b) > min_diff) {
      if(a < b)
         return -1;
      else
         return 1;
   }
   return 0;
}


int cmp_face_angles(const vector<double> &f1, const vector<double> &f2, double min_diff)
{
   if(f1.size() != f2.size())
      return 0;

   for(unsigned int i=0; i<f1.size(); i++) {
      if(int cmp = cmp_angles(f1[i], f2[i], min_diff))
         return cmp;
   }
   return 0;
}

void geom_info::reset()
{
   oriented = -1;
   orientable = -1;
   found_connectivity = false;
   genus_val = INT_MAX;
   dual.clear_all();
   sym = sch_sym();
   efpairs.clear();
   face_angles.clear();
   vert_dihed.clear();
   dihedral_angles.clear();
   e_lengths.clear();
   sol_angles.clear();
   vertex_angles.clear();
   f_areas.clear();
   f_perimeters.clear();
   vert_cons.clear();
   vert_norms.clear();
   found_free_verts = false;
   free_verts.clear();
}


bool geom_info::is_oriented()
{
   if(oriented<0)
      oriented = geom.is_oriented();
   return oriented;
}

bool geom_info::is_orientable()
{
   if(orientable<0) {
      geom_v geom2 = geom;
      number_parts = geom2.orient();
      orientable = geom2.is_oriented();
   }
   return orientable;
}

void geom_info::find_edge_face_pairs()
{
   if(is_oriented())
      geom.get_edge_face_pairs(efpairs, true);
   else
      geom.get_edge_face_pairs(efpairs, false);
}

void geom_info::find_connectivity()
{
   // get a copy of edge face pairs with all faces around an edge
   map<vector<int>, vector<int> > tmp_efpairs;
   if(is_oriented())
      geom.get_edge_face_pairs(tmp_efpairs, false);
   else if(efpairs.size()==0)
      find_edge_face_pairs();

   const map<vector<int>, vector<int> > &pairs =
      is_oriented() ? tmp_efpairs : efpairs;

   known_connectivity = true;
   even_connectivity = true;
   polyhedron = true;
   closed = true;
   map<vector<int>, vector<int> >::const_iterator ei;
   for(ei=pairs.begin(); ei!=pairs.end(); ++ei) {
      if(ei->second.size()==1)     // One faces at an edge
         closed = false;
      if(ei->second.size()!=2)    // Edge not met be exactly 2 faces
         polyhedron = false;
      if(ei->second.size()%2)     // Odd number of faces at an edge
         even_connectivity = false;
      if(ei->second.size()>2)     // More than two faces at an edge
         known_connectivity = false;
   }

   found_connectivity = true;
}

void geom_info::find_f_areas()
{
   int fsz = geom.faces().size();
   f_areas.resize(fsz);
   area.init();
   vol = 0;
   for(int i=0; i<fsz; i++) {
      f_areas[i] = face_area(i);
      if(f_areas[i]<area.min) {
         area.min = f_areas[i];
         area.idx[elem_lims::IDX_MIN] = i;
      }
      if(f_areas[i]>area.max) {
         area.max = f_areas[i];
         area.idx[elem_lims::IDX_MAX] = i;
      }
      area.sum += f_areas[i];
      vol += face_vol(i);
   }
}

void geom_info::find_f_perimeters()
{
   f_perimeters.resize(num_faces());
   const vector<vector<int> > &faces = geom.faces();
   for(unsigned int i=0; i<faces.size(); i++) {
      double perim = 0.0;
      unsigned int fsz = faces[i].size();
      for(unsigned int j=0; j<fsz; j++)
         perim += geom.edge_vec(faces[i][j], faces[i][(j+1)%fsz]).mag();
      f_perimeters[i] = perim;
   }
}


double geom_info::face_vol(int f_no)
{
   double f_vol = 0;
   const vector<vec3d> &verts = geom.verts();
   const vector<int> &face = geom.faces(f_no);
   vec3d V = verts[0];
   vec3d v0 = verts[face[0]];
   for(unsigned int i=1; i<face.size()-1; i++)
      f_vol += vtriple(v0-V, verts[face[i]]-V, verts[face[i+1]]-V);
   return f_vol/6;
}

double geom_info::face_area(int f_no)
{
   return geom.face_norm(f_no, true).mag();
}

void geom_info::face_angles_lengths(int f_no, vector<double> &angs,
      vector<double> *lens)
{
   vector<double> l;
   if(!lens)
      lens = &l;
   const vector<vec3d> &verts = geom.verts();
   const vector<vector<int> > &faces = geom.faces();
   unsigned int fsz = faces[f_no].size();
   angs.resize(fsz);
   lens->resize(fsz);
   vec3d norm;
   double ang_sum=0;
   for(unsigned int j=0; j<fsz; j++) {
      vec3d v0 = verts[faces[f_no][j]]-verts[faces[f_no][(j+fsz-1)%fsz]];
      vec3d v1 = verts[faces[f_no][j]]-verts[faces[f_no][(j+1)%fsz]];
      if(j==0)
         (*lens)[fsz-1] = v0.mag();
      (*lens)[j] = v1.mag();
      if(!norm.is_set())
         norm = vcross(v0,v1);
      angs[j] = acos(safe_for_trig(
                       vdot(v0, v1)/((*lens)[(j+fsz-1)%fsz]*(*lens)[j])  ));
      double sign = vdot(norm, vcross(v0, v1));
      if(sign<0)
         angs[j] = 2*M_PI - angs[j];
      ang_sum += angs[j];
   }
   if(ang_sum/fsz>M_PI)
      for(unsigned int j=0; j<fsz; j++)
         angs[j] = 2*M_PI - angs[j];
}
 
void geom_info::find_face_angles()
{
   ang.init();
   num_angs = 0;
   map<vector<double>, int, ang_vect_less>::iterator fi;
   map<double, int, ang_less>::iterator ai;
   
   for(unsigned int i=0; i<geom.faces().size(); i++) {
      unsigned int fsz = geom.faces(i).size();
      vector<double> f1(fsz), f2(fsz), f_min(fsz);
      face_angles_lengths(i, f1);
      f_min = f1;
      for(unsigned int offset=0; offset<fsz; offset++) {
         pair<int, int> vf_pr(geom.faces(i, offset), i);
         vertex_plane_angs[vf_pr] = f1[offset];
         ai = plane_angles.find(f1[offset]);
         if(ai==plane_angles.end())
            plane_angles[f1[offset]]=1;
         else
            ai->second += 1;
         if(ang.max<f1[offset])
            ang.max = f1[offset];
         if(ang.min>f1[offset])
            ang.min = f1[offset];
         ang.sum += f1[offset];
         num_angs++;
         
         for(unsigned int k=0; k<fsz; k++)
            f2[(k+offset)%fsz] = f1[k];
         if(cmp_face_angles(f2, f_min) < 0)
            f_min = f2;
         reverse(f2.begin(), f2.end());
         if(cmp_face_angles(f2, f_min) < 0)
            f_min = f2;
      }
      fi = face_angles.find(f_min);
      if(fi==face_angles.end())
         face_angles[f_min]=1;
      else
         fi->second += 1;
   }
}

void geom_info::find_dihedral_angles()
{
   if(efpairs.size()==0)
      find_edge_face_pairs();
   edge_dihedrals.resize(efpairs.size());

   dih_angles.init();
   map<vector<int>, vector<int> >::iterator ei; 
   map<double, int, ang_less>::iterator di;
   double cos_a=1, sign=1;
   int e_idx=-1;
   for(ei=efpairs.begin(); ei!=efpairs.end(); ++ei) {
      e_idx++;
      if(ei->second[0]>=0 && ei->second[1]>=0) { // pair of faces
         vec3d n0;
         vec3d n1;
         if(is_oriented()) {
            n0 = geom.face_norm(ei->second[0]).unit();
            n1 = geom.face_norm(ei->second[1]).unit();
            vec3d e_dir = geom.verts(ei->first[1]) - geom.verts(ei->first[0]);
            sign = vdot(e_dir, vcross(n0, n1));
         }
         else {
            vector<int> f0 = geom.faces(ei->second[0]);
            vector<int> f1 = geom.faces(ei->second[1]);
            orient_face(f0, ei->first[0], ei->first[1]);
            orient_face(f1, ei->first[1], ei->first[0]);
            n0 = face_norm(geom.verts(), f0).unit();
            n1 = face_norm(geom.verts(), f1).unit();
            sign = 1;
         }
         cos_a =-vdot(n0, n1);
      }
      else
         cos_a = 1;  // one face on an edge
      double ang = acos(safe_for_trig(cos_a));
      if(sign<0) // in oriented polyhedron
          ang = 2*M_PI - ang;

      edge_dihedrals[e_idx]=ang;
  
      if(ang>dih_angles.max) {
         dih_angles.max = ang;
         dih_angles.idx[elem_lims::IDX_MAX] = ei->first[0];
         dih_angles.idx[elem_lims::IDX_MAX2] = ei->first[1];
      }
      if(ang<dih_angles.min) {
         dih_angles.min = ang;
         dih_angles.idx[elem_lims::IDX_MIN] = ei->first[0];
         dih_angles.idx[elem_lims::IDX_MIN2] = ei->first[1];
      }
      if(fabs(ang - M_PI)<fabs(dih_angles.mid)) {
         dih_angles.mid = ang;
         dih_angles.idx[elem_lims::IDX_MID] = ei->first[0];
         dih_angles.idx[elem_lims::IDX_MID2] = ei->first[1];
      }
      
      di = dihedral_angles.find(ang);
      if(di==dihedral_angles.end())
         dihedral_angles[ang]=1;
      else
         di->second += 1;
   }
}


void geom_info::find_vert_cons()
{
   vert_cons.resize(num_verts(), vector<int>());
   if(!num_faces())
      return;
   geom_v &dual = get_dual();

   vector<int> del_faces;
   const vector<vector<int> > &dfaces = dual.faces(); 
   
   // this needs checking
   for(unsigned int i=0; i<dfaces.size(); i++)
      if(dfaces[i].size()<2)
         del_faces.push_back(i);
   //dual.delete_faces(del_faces);

   if(del_faces.size()) {   // set up an unordered list
      vector<vector<int> > es = geom.edges();
      geom.get_impl_edges(es);
      for(unsigned int i=0; i<es.size(); i++) {
         vert_cons[es[i][0]].push_back(es[i][1]);
         vert_cons[es[i][1]].push_back(es[i][0]);
      }
      return;
   }
   
   /* if(!dual.get_faces()->size()) {   // set up an unordered list
      for(unsigned int i=0; i<geom.edges().size(); i++) {
         vert_cons[geom.edges(i, 0)].push_back(geom.edges(i, 1));
         vert_cons[geom.edges(i, 1)].push_back(geom.edges(i, 0));
      }
      return;
   } */

   dual.orient();
   const vector<vector<int> > &faces = geom.faces();
   for(unsigned int i=0; i<dfaces.size(); i++) {
      for(unsigned int j=0; j<dfaces[i].size(); j++) {
         int f = dfaces[i][j];
         int sz = faces[f].size();
         for(int v=0; v<sz; v++) {
            if(faces[f][v]==(int)i) { // found vertex in surrounding face
               if(j==0) {  // first vertex, must be in following face
                  int idx = faces[f][(v+1)%sz];
                  const vector<int> &nface =
                              faces[dfaces[i][(j+1)%dfaces[i].size()]];
                  if(find(nface.begin(), nface.end(), idx)!=nface.end())
                     vert_cons[i].push_back(faces[f][(v+1)%sz]);
                  else
                     vert_cons[i].push_back(faces[f][(v+sz-1)%sz]);
               }
               else {
                  if(faces[f][(v+1)%sz]==vert_cons[i][j-1])
                     vert_cons[i].push_back(faces[f][(v+sz-1)%sz]);
                  else 
                     vert_cons[i].push_back(faces[f][(v+1)%sz]);
               }
               break;
            }
         }     
      }
   }
   
   // first dual vertex on dual face 0 is a face containing vertex 0
   int f0 = dfaces[0][0];
   // find vertices before and after 0 
   vector<int> edge(2);
   int sz = faces[f0].size();
   for(int i=0; i<sz; i++)
      if(faces[f0][i]==0) {
         edge[0] = faces[f0][(i-1+sz)%sz];
         edge[1] = faces[f0][(i+1)%sz];
         break;
      }
 
   for(unsigned int j=0; j<vert_cons[0].size(); j++) {
      if(vert_cons[0][j] == edge[0]) {
         if(vert_cons[0][(j+1)%vert_cons[0].size()]!=edge[1]) {
            for(unsigned int v=0; v<vert_cons.size(); v++)
               reverse(vert_cons[v].begin(), vert_cons[v].end());
         }
         break;
      }
   }
}

void geom_info::find_vert_norms()
{
   const geom_if &geom = get_geom();
   const unsigned int verts_sz = geom.verts().size();
   vert_norms.assign(verts_sz, vec3d::zero);
   vector<int> vert_face_cnt(geom.verts().size(), 0);
   for(unsigned int i=0; i<geom.faces().size(); i++) {
      vec3d norm = geom.face_norm(i).unit();
      for(unsigned int j=0; j<geom.faces(i).size(); j++) {
         const int v_idx = geom.faces(i, j);
         vert_norms[v_idx] += norm;
         vert_face_cnt[v_idx]++;
      }
   }

   for(unsigned int i=0; i<verts_sz; i++) {
      if(int cnt = vert_face_cnt[i])
         vert_norms[i] /= cnt;
      else
         vert_norms[i].unset();
   }
}
   

void geom_info::find_free_verts()
{
   free_verts.clear();
   const geom_if &geom = get_geom();
   vector<int> cnt(geom.verts().size());
   for(unsigned int i=0; i<geom.faces().size(); i++)
      for(unsigned int j=0; j<geom.faces(i).size(); j++)
         cnt[geom.faces(i, j)]++;
   for(unsigned int i=0; i<geom.edges().size(); i++)
      for(unsigned int j=0; j<geom.edges(i).size(); j++)
         cnt[geom.edges(i, j)]++;
  
   for(int i=0; i<(int)cnt.size(); i++)
      if(cnt[i]==0)
         free_verts.push_back(i);

   found_free_verts = true;
}


static double sph_tri_area(vec3d u0, vec3d u1, vec3d u2)
{
   double sign = 1 - 2*(vtriple(u0, u1, u2)>0);
   vec3d u[3] = {u0, u1, u2};
   double ang, area=0;
   for(int i=0; i<3; i++) {
      ang = acos(safe_for_trig(vdot(vcross(u[(i-1+3)%3], u[i]).unit(),
                       vcross(u[i], u[(i+1)%3]).unit()) ));
      if(sign<0)
         ang = 2*M_PI - ang;
      area += ang;
   }
   return area - M_PI;
}

void geom_info::find_solid_angles()
{
   if(!vert_cons.size())
      find_vert_cons();
   
   vertex_angles = vector<double> (num_verts(), 0);
  
   so_angles.init();
   map<double, int, ang_less>::iterator si;
   for(unsigned int i=0; i<vert_cons.size(); i++) {
      vector<vec3d> dirs(vert_cons[i].size());
      for(unsigned int j=0; j<vert_cons[i].size(); j++)
         dirs[j] = geom.verts(i) - geom.verts(vert_cons[i][j]);

      for(unsigned int j=1; j<vert_cons[i].size()-1; j++)
         vertex_angles[i] += sph_tri_area(dirs[0], dirs[j], dirs[j+1]);

      //if(!is_oriented()) {
      //   fmod(vertex_angles[i], 4*M_PI);
      //   if(vertex_angles[i]>2*M_PI)
      //      vertex_angles[i] = 4*M_PI - vertex_angles[i];
     // }

      
      if(vertex_angles[i]>so_angles.max) {
         so_angles.max = vertex_angles[i];
         so_angles.idx[elem_lims::IDX_MAX] = i;
      }
      if(vertex_angles[i]<so_angles.min) {
         so_angles.min = vertex_angles[i];
         so_angles.idx[elem_lims::IDX_MIN] = i;
      }
      if(fabs(vertex_angles[i])<fabs(so_angles.mid)) {
         so_angles.mid = vertex_angles[i];
         so_angles.idx[elem_lims::IDX_MID] = i;
      }
      
      si = sol_angles.find(vertex_angles[i]);
      if(si==sol_angles.end())
         sol_angles[vertex_angles[i]]=1;
      else
         si->second += 1;
   }
}
   
void geom_info::find_e_lengths(map<double, int, ang_less> &e_lens,
      const vector<vector<int> > &edges, elem_lims &lens)
{
   lens.init();
   map<double, int, ang_less>::iterator ei;
   for(unsigned int i=0; i<edges.size(); i++) {
      double dist = geom.edge_len(edges[i]);
      if(dist>lens.max) {
         lens.max = dist;
         lens.idx[elem_lims::IDX_MAX] = edges[i][0];
         lens.idx[elem_lims::IDX_MAX2] = edges[i][1];
      }
      if(dist<lens.min) {
         lens.min = dist;
         lens.idx[elem_lims::IDX_MIN] = edges[i][0];
         lens.idx[elem_lims::IDX_MIN2] = edges[i][1];
      }
      lens.sum += dist; 
      ei = e_lens.find(dist);
      if(ei==e_lens.end())
         e_lens[dist]=1;
      else
         ei->second += 1;
   }
}


void geom_info::find_f_max_nonplanars()
{
   f_max_nonplanars.resize(geom.faces().size());
   for(unsigned int f=0; f<geom.faces().size(); f++) {
      if(geom.faces(f).size()==3) {
         f_max_nonplanars[f] = 0;
         continue;
      }
      vec3d norm = geom.face_norm(f).unit();
      vec3d f_cent = geom.face_cent(f);
      double max = 0;
      for(unsigned int v=0; v<geom.faces(f).size(); v++) {
         double dist = fabs(vdot(norm, f_cent - geom.verts(geom.faces(f, v))));
         if(dist>max)
            max = dist;
      }
      f_max_nonplanars[f] = max;
   }
}

void geom_info::find_symmetry()
{
   sym.init(geom);
}


bool geom_info::is_closed()
{
   if(!found_connectivity)
      find_connectivity();
   return closed;
}


bool geom_info::is_polyhedron()
{
   if(!found_connectivity)
      find_connectivity();
   return polyhedron;
}


bool geom_info::is_even_connectivity()
{
   if(!found_connectivity)
      find_connectivity();
   return even_connectivity;
}


bool geom_info::is_known_connectivity()
{
   if(!found_connectivity)
      find_connectivity();
   return known_connectivity;
}


void v_distances::set_values()
{
   init();
   for(unsigned int i=0; i<geom.verts().size(); i++) {
      double dist = (geom.verts(i) - center).mag();
      if(dist < min) {
         min = dist;
         idx[elem_lims::IDX_MIN] = i;
      }
      if(dist > max) {
         max = dist;
         idx[elem_lims::IDX_MAX] = i;
      }
      sum += dist;
   }
}

bool geom_info::is_known_genus()
{
   return genus()<INT_MAX-1;
}

int geom_info::genus()
{
   if(genus_val==INT_MAX) {
      genus_val = INT_MAX-1;                    // 'not known' value
      if(num_parts()==1  && is_known_connectivity()) {
         int euler_char = num_verts() -  num_edges() + num_faces();
         geom_v geom2 = geom;
         if(close_poly_basic(geom2)) {
            int num_boundaries  = geom2.faces().size() - num_faces();
            genus_val = 2 - num_boundaries - euler_char;
            if(is_orientable())
               genus_val /= 2;
            else
               genus_val *= -1; // to indicate non-orientable genus (demigenus)
         }
      }
   }
   return genus_val;
}

void e_distances::set_values()
{
   init();
   for(unsigned int i=0; i< geom.edges().size(); i++) {
      double dist = geom.edge_nearpt(i, center).mag();
      if(dist < min) {
         min = dist;
         idx[elem_lims::IDX_MIN] =  geom.edges(i, 0);
         idx[elem_lims::IDX_MIN2] = geom.edges(i, 1);
      }
      if(dist > max) {
         max = dist;
         idx[elem_lims::IDX_MAX] =  geom.edges(i, 0);
         idx[elem_lims::IDX_MAX2] = geom.edges(i, 1);
      }
      sum += dist;
   }
}


void ie_distances::set_values()
{
   init();
   vector<vector<int> > edges;
   geom.get_impl_edges(edges);
   for(unsigned int i=0; i< edges.size(); i++) {
      double dist = (geom.edge_nearpt(edges[i], center) - center).mag();
      if(dist < min) {
         min = dist;
         idx[elem_lims::IDX_MIN] =  edges[i][0];
         idx[elem_lims::IDX_MIN2] = edges[i][1];
      }
      if(dist > max) {
         max = dist;
         idx[elem_lims::IDX_MAX] =  edges[i][0];
         idx[elem_lims::IDX_MAX2] = edges[i][1];
      }
      sum += dist;
   }
}


void f_distances::set_values()
{
   init();
   for(unsigned int i=0; i<geom.faces().size(); i++) {
      double dist = geom.face_nearpt(i, center).mag();
      if(dist < min) {
         min = dist;
         idx[elem_lims::IDX_MIN] = i;
      }
      if(dist > max) {
         max = dist;
         idx[elem_lims::IDX_MAX] = i;
      }
      sum += dist;
   }
}



