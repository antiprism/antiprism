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

/*\file geom.cc
 * \brief Classes to represent a geometry
*/


#include <stdlib.h>
#include <stdarg.h>
#include <algorithm>

#include "geom.h"
#include "info.h"
#include "off_file.h"
#include "geodesic.h"
#include "std_polys.h"
#include "symmetry.h"
#include "coloring.h"

int add_hull(geom_if &geom, string qh_args="", char *errmsg=0);
int set_hull(geom_if &geom, string qh_args="", char *errmsg=0);
int orient_geom(geom_if &geom, vector<vector<int> > *parts=0);
bool orient_geom(geom_if &geom, int type, char *errmsg=0);
void orient_reverse(geom_if &geom);
int triangulate(geom_if &geom, col_val inv=col_val(),
      unsigned int winding_rule=TESS_WINDING_NONZERO, vector<int> *fmap=0);
void get_star(const geom_if &geom, vector<vec3d> &star, char type='v',
      vec3d centre=vec3d(0,0,0));
bool make_zono(geom_if &zono, const vector<vec3d> &star, char *errmsg=0);

using std::swap;
using std::pair;

// geom_if implementation

int geom_if::add_vert(vec3d vert)
{
   raw_verts().push_back(vert);
   return verts().size()-1;
}
      
int geom_if::add_verts(const vector<vec3d> &vs)
{
   for(vector<vec3d>::const_iterator vi=vs.begin(); vi!=vs.end(); vi++)
      add_vert(*vi);
   return verts().size()-1;
}
         
int geom_if::add_face(const vector<int> &face)
{
   raw_faces().push_back(face);
   return faces().size()-1;
}

int geom_if::add_face(int v1, ...)
{
   vector<int> face;
   va_list ap;
   va_start(ap, v1); 
   for(int i=v1; i!=-1; i=va_arg(ap, int))
      face.push_back(i);
   va_end(ap);
   return add_face(face);
}



int geom_if::add_faces(const vector<vector<int> > &fs)
{
   for(vector<vector<int> >::const_iterator vi=fs.begin(); vi!=fs.end(); vi++)
      add_face(*vi);
   return faces().size()-1;
}
         
int geom_if::add_edge_raw(const vector<int> &edge)
{ 
   raw_edges().push_back(edge);
   return edges().size()-1;
}
      
int geom_if::add_edges_raw(const vector<vector<int> > &es)
{
   for(vector<vector<int> >::const_iterator vi=es.begin(); vi!=es.end(); vi++)
      add_edge_raw(*vi);
   return edges().size()-1;
}
         
      

int geom_if::add_edge(vector<int> edge)
{
   if(edge[0]>edge[1])
      swap(edge[0], edge[1]);
   vector<vector<int> >::const_iterator ei = find(edges().begin(),
         edges().end(), edge);
   if(ei != edges().end())
      return ei - edges().begin();
   raw_edges().push_back(edge);
   return edges().size()-1;
}
      

int geom_if::add_edge(int v_idx1, int v_idx2)
{
   return add_edge(make_edge(v_idx1, v_idx2));
}


int geom_if::add_edges(const vector<vector<int> > &es)
{
   for(vector<vector<int> >::const_iterator vi=es.begin(); vi!=es.end(); vi++)
      add_edge(*vi);
   return edges().size()-1;
}
         

void geom_if::delete_vert(int vert, map<int, int> *vert_map)
{
   vector<int> del_verts;
   del_verts.push_back(vert);
   delete_verts(del_verts, vert_map);
}

 
void geom_if::delete_verts(const vector<int> &v_nos, map<int, int> *vert_map)
{
   vector<int> dels = v_nos;
   map<int, int> tmp;
   map<int, int> *v_map = (vert_map) ? vert_map : &tmp;
   v_map->clear();
   if(!dels.size())
      return;
   sort(dels.begin(), dels.end());
   unsigned int del_verts_cnt=0;
   int map_to;
   for(unsigned int i=0; i<verts().size(); i++) {
      if(del_verts_cnt<dels.size() && (int)i==dels[del_verts_cnt]) {
         del_verts_cnt++;
         map_to = -1;
      }
      else {
         map_to = i - del_verts_cnt;
         raw_verts()[map_to] = verts(i);
      }
      (*v_map)[i] = map_to;
   }
   raw_verts().resize(verts().size()-del_verts_cnt);

   vector<int> del_faces;
   for(unsigned int i=0; i<faces().size(); i++) {
      int curr_idx = 0;
      for(unsigned int j=0; j<faces(i).size(); j++) {
         int new_idx = (*v_map)[faces(i, j)];
         if(new_idx >= 0)
            raw_faces()[i][curr_idx++] = new_idx;
      }
      raw_faces()[i].resize(curr_idx);
      if(curr_idx<3)
         del_faces.push_back(i);
   }
   delete_faces(del_faces);
   
   vector<int> del_edges;
   for(unsigned int i=0; i<edges().size(); i++)
      for(unsigned int j=0; j<2; j++) {
         int new_idx = (*v_map)[edges(i, j)];
         if(new_idx >= 0)
            raw_edges()[i][j] = new_idx;
         else {
            del_edges.push_back(i);
            break;
         }
      }
   delete_edges(del_edges);
}
   


void geom_if::delete_faces(const vector<int> &f_nos, map<int, int> *face_map)
{
   vector<int> dels = f_nos;
   if(face_map)
      face_map->clear();
   if(!dels.size())
      return;
   sort(dels.begin(), dels.end());
   unsigned int del_faces_cnt=0;
   int map_to;
   for(unsigned int i=0; i<faces().size() ; i++) {
      if(del_faces_cnt<dels.size() && (int)i==dels[del_faces_cnt]) {
         del_faces_cnt++;
         map_to = -1;
      }
      else {
         map_to = i - del_faces_cnt;
         raw_faces()[map_to] = faces(i);
      }
      if(face_map)
         (*face_map)[i] = map_to;
   }
   raw_faces().resize(faces().size()-del_faces_cnt);
}
   

void geom_if::delete_edges(const vector<int> &e_nos, map<int, int> *edge_map)
{
   vector<int> dels = e_nos;
   if(edge_map)
      edge_map->clear();
   if(!dels.size())
      return;
   sort(dels.begin(), dels.end());
   unsigned int del_edges_cnt=0;
   int map_to;
   for(unsigned int i=0; i<edges().size() ; i++) {
      if(del_edges_cnt<dels.size() && (int)i==dels[del_edges_cnt]) {
         del_edges_cnt++;
         map_to = -1;
      }
      else {
         map_to = i - del_edges_cnt;
         raw_edges()[map_to] = edges(i);
      }
      if(edge_map)
         (*edge_map)[i] = map_to;
   }
   raw_edges().resize(edges().size()-del_edges_cnt);
}
   

void remap_shift(vector<vector<int> > &elems, int offset)
{
   for(unsigned int i=0; i<elems.size(); i++)
      for(unsigned int j=0; j<elems[i].size(); j++)
         elems[i][j]+=offset;
}



void geom_if::append(const geom_if &geom)
{
   vector<vec3d> g_verts = geom.verts();
   vector<vector<int> > g_faces = geom.faces();
   vector<vector<int> > g_edges = geom.edges();
   
   remap_shift(g_edges, verts().size());
   remap_shift(g_faces, verts().size());
   
   raw_verts().insert(raw_verts().end(), g_verts.begin(), g_verts.end());
   raw_edges().insert(raw_edges().end(), g_edges.begin(), g_edges.end());
   raw_faces().insert(raw_faces().end(), g_faces.begin(), g_faces.end());
}

void geom_if::clear_verts() { raw_verts().clear(); }
void geom_if::clear_edges() { raw_edges().clear(); }
void geom_if::clear_faces() { raw_faces().clear(); }
void geom_if::clear_all() { clear_verts(); clear_edges(); clear_faces(); }


void geom_if::get_impl_edges(vector<vector<int> > &edgs) const
{
   // edgs hasn't been cleared
   for(unsigned int i=0; i<faces().size(); ++i)
      for(unsigned int j=0; j<faces(i).size(); ++j)
         edgs.push_back(
               make_edge(faces(i, j), faces(i, (j+1)%faces(i).size())) );

   // Clear duplicate edges (each edge appears in two faces)
   sort(edgs.begin(), edgs.end());
   vector<vector<int> >::iterator vit = unique(edgs.begin(), edgs.end());
   edgs.erase(vit, edgs.end());
}

int geom_if::add_hull(string qh_args, char *errmsg)
{
   return ::add_hull(*this, qh_args, errmsg);
}

int geom_if::set_hull(string qh_args, char *errmsg)
{
   return ::set_hull(*this, qh_args, errmsg);
}

vector<vec3d> geom_if::get_star(char type, vec3d centre)
{
   vector<vec3d> star;
   ::get_star(*this, star, type, centre);
   return star;
}

bool geom_if::set_zono(const vector<vec3d> &star, char *errmsg)
{
   return make_zono(*this, star, errmsg);
}
      
int geom_if::orient(vector<vector<int> > *parts)
{
   return orient_geom(*this, parts);
}

bool geom_if::orient(int type, char *errmsg)
{
   return orient_geom(*this, type, errmsg);
}

bool geom_if::set_geodesic_planar(const geom_if &base, int m, int n)
{
   if(m<0 || n<0 || (m==0&&n==0))
      return false;                         // invalid pattern
   geodesic geod(base, m, n, 'p');
   geod.make_geo(*this);
   return true;                             // valid pattern
}
      
bool geom_if::set_geodesic_sphere(const geom_if &base, int m, int n, vec3d cent)
{
   if(m<0 || n<0 || (m==0&&n==0))
      return false;                         // invalid pattern
   geodesic geod(base, m, n, 's', cent);
   geod.make_geo(*this);
   return true;                             // valid pattern
}


void geom_if::sphere_projection(vec3d centre, double radius)
{
   vector<vec3d> &verts = *get_verts();
   for(unsigned int i=0; i<verts.size(); i++)
      verts[i] = (verts[i]-centre).unit()*radius;
}


void geom_if::orient_reverse()
{
   return ::orient_reverse(*this);
}


void geom_if::sym_align()
{
   transform(sch_sym(*this).get_to_std());
}


void geom_if::triangulate(col_val col, unsigned int winding, vector<int> *fmap)
{
   ::triangulate(*this, col, winding, fmap);
}

bool geom_if::read(string file_name, char *errmsg)
{
   return off_file_read(file_name, *this, errmsg);
}

bool geom_if::read(FILE *file, char *errmsg)
{
   return off_file_read(file, *this, errmsg);
}

bool geom_if::read_resource(string res_name, char *errmsg)
{
   return make_resource_geom(*this, res_name, errmsg);
}

bool geom_if::write(string file_name, char *errmsg, int sig_dgts) const
{
   return off_file_write(file_name, *this, errmsg, sig_dgts);
}

void geom_if::write(FILE *file, int sig_dgts) const
{
   off_file_write(file, *this, sig_dgts);
}

bool geom_if::write_crds(string file_name, char *errmsg, const char *sep,
      int sig_dgts) const
{
   return crds_write(file_name, *this, errmsg, sep, sig_dgts);
}

void geom_if::write_crds(FILE *file, const char *sep, int sig_dgts) const
{
   return crds_write(file, *this, sep, sig_dgts);
}

geom_info geom_if::get_info() const
{
   return geom_info(*this);
}


vector<int> make_edge(int v_idx1, int v_idx2)
{
   vector<int> edge(2);
   edge[0] = v_idx1;
   edge[1] = v_idx2;
   if(edge[0]>edge[1])
      swap(edge[0], edge[1]);
   return edge;
}

bool geom_if::is_oriented() const
{
   set<pair<int, int> > edgs;
   pair<int, int> edge;
   for(unsigned int f=0; f<faces().size(); f++)
      for(unsigned int i=0; i<faces(f).size(); i++) {
         edge = pair<int, int>(faces(f, i), faces(f, (i+1)%faces(f).size()));
         if(edgs.find(edge)!=edgs.end())
            return false;
         else
            edgs.insert(edge);
      }
   return true;
}



void geom_if::add_missing_impl_edges()
{
   // save original edges
   vector<vector<int> > e_edges = edges();
   
   clear_edges();
   get_impl_edges(raw_edges());

   // restore original edges and colours
   for(unsigned int e=0; e<e_edges.size(); ++e)
      add_edge(e_edges[e]);
}


void geom_if::get_edge_face_pairs(map<vector<int>, vector<int> > &edge2facepr,
            bool oriented) const
{
   map<vector<int>, vector<int> >::iterator mi;
   vector<int> vrts(2);
   for(unsigned int i=0; i<faces().size(); ++i) {
      for(unsigned int j=0; j<faces(i).size(); ++j) {
         vrts[0]=faces(i, j);
         vrts[1]=faces(i, (j+1)%faces(i).size());
         int face_pos=0;
         if(vrts[0] > vrts[1]) {
            swap(vrts[0], vrts[1]);
            face_pos=1;
         }
         if(oriented) {
            mi = edge2facepr.find(vrts);
            if(mi==edge2facepr.end()) {
               edge2facepr[vrts].resize(2);
               edge2facepr[vrts][(face_pos+1)%2]=-1;
            }
            edge2facepr[vrts][face_pos]=i;
         }
         else
            edge2facepr[vrts].push_back(i);
      }
   }
}


void geom_if::verts_merge(map<int, int> &vmap)
{
   map<int, int>::iterator vmi;
   vector<int> del_faces;
   for(unsigned int i=0; i<faces().size(); i++) {
      for(unsigned int j=0; j<faces(i).size(); j++) {
         int idx = faces(i, j);
         vmi = vmap.find(idx);
         if(vmi != vmap.end())
            raw_faces()[i][j] = vmap[idx];
      }
      vector<int>::iterator vi = unique(raw_faces()[i].begin(), raw_faces()[i].end());
      if(faces(i, 0)==*(vi-1))
         vi--;
      raw_faces()[i].resize(vi-faces(i).begin());
      if(faces(i).size()<3)
         del_faces.push_back(i);
   }
   delete_faces(del_faces);

   vector<int> del_edges;
   for(unsigned int i=0; i<edges().size(); i++) {
      for(unsigned int j=0; j<edges(i).size(); j++) {
         int idx = edges(i, j);
         vmi = vmap.find(idx);
         if(vmi != vmap.end())
            raw_edges()[i][j] = vmap[idx];
      }
      if(edges(i, 0)==edges(i, 1))
         del_edges.push_back(i);
   }
   delete_edges(del_edges);

   vector<int> del_verts;
   for(vmi=vmap.begin(); vmi!=vmap.end(); vmi++)
      del_verts.push_back(vmi->first);
   delete_verts(del_verts);
}


col_geom_v &col_geom_v::operator =(const geom_if &geom)
{
   clear_all();
   append(geom);
   return *this;
}


void col_geom_v::add_missing_impl_edges()
{
   // save original edges and colours
   vector<vector<int> > e_edges = edges();
   map<int, col_val> cols = edge_cols();
   
   clear_edges();
   get_impl_edges(raw_edges());

   // restore original edges and colours
   for(unsigned int e=0; e<e_edges.size(); ++e)
      add_col_edge(e_edges[e], get_col(cols, e));
}



void col_geom_v::color_vef(col_val vert_col, col_val edge_col, col_val face_col)
{
   coloring clrng(this);
   
   if (vert_col.is_set())
      clrng.v_one_col(vert_col);

   if (edge_col.is_set()) {
      add_missing_impl_edges();
      clrng.e_one_col(edge_col);
   }

   if (face_col.is_set())
      clrng.f_one_col(face_col);
}



void col_geom_v::delete_verts(const vector<int> &v_nos, map<int, int> *vert_map)
{
   map<int, int> tmp;
   map<int, int> *v_map = (vert_map) ? vert_map : &tmp;
   geom_if::delete_verts(v_nos, v_map);
   remap_vert_cols(*v_map);
}

void col_geom_v::delete_faces(const vector<int> &f_nos, map<int, int> *face_map)
{
   map<int, int> tmp;
   map<int, int> *f_map = (face_map) ? face_map : &tmp;
   geom_if::delete_faces(f_nos, f_map);
   remap_face_cols(*f_map);
}

void col_geom_v::delete_edges(const vector<int> &e_nos, map<int, int> *edge_map)
{
   map<int, int> tmp;
   map<int, int> *e_map = (edge_map) ? edge_map : &tmp;
   geom_if::delete_edges(e_nos, e_map);
   remap_edge_cols(*e_map);
}

void col_geom_v::append(const geom_if &geom)
{
   const col_geom *cg = dynamic_cast<const col_geom *>(&geom);
   if(cg)
      col_geom::append(*cg, verts().size(), edges().size(), faces().size());
   geom_if::append(geom);
}




