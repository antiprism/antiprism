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
   Name: geodesics.h
   Description: Generate geodesic spheres and polyhedra 
   Project: Antiprism - http://www.antiprism.com
*/

#ifndef GEODESICS_H
#define GEODESICS_H


#include <string>
#include <vector>
#include <map>

#include "geom.h"

using std::string;
using std::vector;
using std::map;
using std::pair;

typedef pair<int, int> int_pr;
inline int_pr mk_int_pr(int x, int y)
   { int_pr pr; pr.first=x; pr.second=y; return pr; }


class ij_pos {
   private:
      int pos;
   public:
      enum {face=0, v0=3, v1=6, v2=5, e0=2, e1=4, e2=1, out=7};
      ij_pos(int p): pos(p) {}
      bool operator == (int pos2) { return pos==pos2; }
      bool operator == (ij_pos p2) { return pos==p2.pos; }
      bool is_vert() {return (pos==v0||pos==v1||pos==v2);}
      bool is_edge() {return (pos==e0||pos==e1||pos==e2);}
      bool is_face() {return (pos==face);}
      bool is_out() {return (pos==out);}
      int to_int() { return pos; }
      string dump() { const char *txt[] = {"face", "e2", "e0", "v0", "e1", "v2", "v1", "out"}; return txt[pos];}
};


class geodesic
{
   private:
      enum { noindex=-1 };
      
      col_geom_v base;
      int freq;
      int F;
      int m;
      int n;
      char method;
      vec3d centre;
      
      map<vector<int>, int> edge_idx;
      map<vector<int>, vector<int> > edge_faces;
      //map<vector<int>, int> face_idx;
      map<int_pr, int> grid_idxs;

      void init();
      void sphere_projection(geom_if &geom);
      void make_grid_idxs();
      int grid_x(int i, int j) { return i*(-n) + j*(m+n); }
      int grid_y(int i, int j) { return i*(m+n) + j*(-m); }
      int_pr rot_e0(int_pr crds) // half-rot about centre e0
         { return mk_int_pr(freq-crds.first, -crds.second); }
      int_pr rot_f(int_pr crds)   // third-rot about centre f
         { return mk_int_pr(freq-crds.first-crds.second, crds.first); }
      int_pr normal_crds(int_pr crds);
      
      int coord_i(int_pr crds)
         { return (m*crds.first + (m+n)*crds.second)/(m*m+m*n+n*n); }
      int coord_j(int_pr crds)
         { return ((m+n)*crds.first + n*crds.second)/(m*m+m*n+n*n); }
      
      void grid_to_points(vector<int> indx, vector<vec3d> &gverts);
      bool tri_test(int i, int j, int di, int dj);
      void grid_to_tris(vector<int> indx, vector<vector<int> > &new_tris, vector<vector<int> > &orig_edges);
      vector<int> make_face_indexes(int i, const vector<int> &face);
      int index_map(int i, int j, const vector<int> &indx, int p_idx=noindex);

      int get_edge_index(int v0, int v1);
      vector<int> get_face_indexes(vector<int> face);

      ij_pos get_pos_xy(int x, int y) {
         if(x<0 || y<0 || x+y>freq) return ij_pos::out;
         else return ij_pos((x==0) + 2*(y==0) + 4*(x+y==freq)); }
      ij_pos get_pos(int i, int j)
         { return get_pos_xy(grid_x(i, j), grid_y(i, j)); }

   public:
      enum { err_not_tri=1 };
      geodesic(const geom_if &base_poly, int mm, int nn=0, char mthd='s',
            vec3d cen=vec3d(0,0,0));
      void make_geo(geom_if &geo);
};


#endif // GEODESICS_H
  

