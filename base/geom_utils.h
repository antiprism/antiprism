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

/*!\file geom_utils.h
 * \brief Utilities associated with geometries
*/

#ifndef GEOM_UTILS_H
#define GEOM_UTILS_H

void orient_face(vector<int> &face, int v0, int v1);

void triangulate_basic(geom_if &geom, bool sq_diag=true, col_val inv=col_val(),
      vector<int> *fmap=0);

double minimum_distance(const geom_if &geom, double sig_dist=0, char *errmsg=0);

void make_polar_zono(geom_if &zono, int star_n, bool out_star);


// RK - test points versus hull functions
bool is_geom_fully_outside_hull(const geom_if &, const geom_if &, double);
bool is_geom_outside_hull(const geom_if &, const geom_if &, double);
bool is_geom_on_surface_hull(const geom_if &, const geom_if &, double);
bool is_geom_inside_hull(const geom_if &, const geom_if &, double);
bool is_geom_fully_inside_hull(const geom_if &, const geom_if &, double);
bool is_point_fully_outside_hull(const vec3d &, const geom_if &, double);
bool is_point_outside_hull(const vec3d &, const geom_if &, double);
bool is_point_on_surface_hull(const vec3d &, const geom_if &, double);
bool is_point_inside_hull(const vec3d &, const geom_if &, double);
bool is_point_fully_inside_hull(const vec3d &, const geom_if &, double);


// RK - Various find functions for geom
int find_vertex_by_coordinate(geom_if &geom, vec3d v, double eps=epsilon);
bool edge_exists_in_face(const vector<int> &face, const vector<int> &edge);
bool are_edges_equal(const vector<int> &edge2, const vector<int> &edge1);
vector<int> find_faces_with_edge(const vector<vector<int> > &faces, const vector<int> &edge);
bool vertex_exists_in_face(const vector<int> &face, int v_idx);
bool vertex_exists_in_edge(const vector<int> &edge, int v_idx);
vector<int> find_faces_with_vertex(const vector<vector<int> > &faces, int v_idx);
vector<int> find_edges_with_vertex(const vector<vector<int> > &edges, int v_idx);
int find_edge_in_edge_list(const vector<vector<int> > &edges, const vector<int> &edge);
vector<vector<int> > find_unmatched_edges(col_geom_v &geom);

// RK - the normals classes. xnormals for one normal. fnormals for the list of face normals.
class xnormal
{
   private:
      vec3d normal;
      int direction; // 1 = outward  0 = hemispherical  -1 = inward
      
   public:
      xnormal() { direction = 0; } // unset normal
      xnormal(const geom_if &geom, int face_idx, vec3d C=vec3d(), double eps=epsilon); // face normal
      xnormal(const geom_if &geom, vec3d norm, int v_idx, vec3d C=vec3d(), double eps=epsilon); // for precalculated normals

      virtual ~xnormal() {}

      bool is_set() const { return normal.is_set(); }

      vec3d raw() const { return normal; }  
      vec3d unit() const { return normal.unit(); }  // unit() of raw normal  
      vec3d outward() const { return (direction>0 ? normal : -normal); }
      vec3d inward() const { return (direction<0 ? normal : -normal); }

      bool is_outward() const { return (direction==1); }
      bool is_inward() const { return (direction==-1); }      
      bool is_hemispherical() const { return (direction==0); }
};

class fnormals
{
   private:
      const geom_if *ngeom;
      vector<xnormal> normals;

      vec3d average_normals(vector<int> &face_idx, string average_pattern);
      
   public:
      fnormals() {};
      fnormals(const geom_if &geom, vec3d C=vec3d(), double eps=epsilon) { refresh(geom, C, eps); }
      void refresh(const geom_if &, vec3d C=vec3d(), double eps=epsilon);
      virtual ~fnormals() {}

      xnormal edge_normal(int idx1, int idx2, string average_pattern, vec3d C=vec3d(), double eps=epsilon);
      xnormal vertex_normal(int idx, string average_pattern, vec3d C=vec3d(), double eps=epsilon);
      
      unsigned int size() const { return (normals.size()); }
      bool in_range(unsigned int idx) const { return (idx < (unsigned int)size()); }
      bool is_set() const { return (size() > 0); }

      xnormal operator [](unsigned int idx) const { return ((is_set() && in_range(idx)) ? normals[idx] : xnormal()); }
};

#endif // GEOM_UTILS_H

