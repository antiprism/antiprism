/*
   Copyright (c) 2010-2011, Adrian Rossiter
   
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
   \brief utilities for geometries
*/

#include <stdio.h>
#include <stdlib.h>

#include <vector>

#include "geom.h"
#include "math_utils.h"
#include "geom_utils.h"


// RK - test points versus hull functions

bool test_points_vs_hull(const vector<vec3d> &P, const geom_if &hull, const bool &inside, const bool &surface, const bool &outside, const double &eps)
{
   const vector<vec3d> &verts = hull.verts();
   const vector<vector<int> > &faces = hull.faces();

   vec3d C = centroid(verts);

   bool answer = true;
   for(unsigned int i=0;i<faces.size();i++) {
      vec3d n = face_norm(verts, faces[i]).unit();
      double D = vdot(verts[faces[i][0]]-C, n);
      //if(D < 0)
      if(double_compare(D,0,eps) < 0) { // Make sure the normal points outwards
         D = -D;
         n = -n;
      }

      for(unsigned int j=0;j<P.size();j++) {
         double t = vdot(P[j]-C, n);
         //if (t < D-eps && !inside)
         if ((double_compare(t,D,eps) < 0) && !inside)
            answer = false;
         else
         //if (t > D+eps && !outside)
         if ((double_compare(t,D,eps) > 0) && !outside)
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

bool is_geom_fully_outside_hull(const geom_if &geom, const geom_if &hull, double eps)
{
   return test_points_vs_hull(geom.verts(),hull,false,false,true,eps);
}

bool is_geom_outside_hull(const geom_if &geom, const geom_if &hull, double eps)
{
   return test_points_vs_hull(geom.verts(),hull,false,true,true,eps);
}

bool is_geom_on_surface_hull(const geom_if &geom, const geom_if &hull, double eps)
{
   return test_points_vs_hull(geom.verts(),hull,false,true,false,eps);
}

bool is_geom_inside_hull(const geom_if &geom, const geom_if &hull, double eps)
{
   return test_points_vs_hull(geom.verts(),hull,true,true,false,eps);
}

bool is_geom_fully_inside_hull(const geom_if &geom, const geom_if &hull, double eps)
{
   return test_points_vs_hull(geom.verts(),hull,true,false,false,eps);
}

bool is_point_fully_outside_hull(const vec3d &P, const geom_if &hull, double eps)
{
   geom_v tgeom;
   tgeom.add_vert(P);
   return is_geom_fully_outside_hull(tgeom, hull, eps);
}

bool is_point_outside_hull(const vec3d &P, const geom_if &hull, double eps)
{
   geom_v tgeom;
   tgeom.add_vert(P);
   return is_geom_outside_hull(tgeom, hull, eps);
}

bool is_point_on_surface_hull(const vec3d &P, const geom_if &hull, double eps)
{
   geom_v tgeom;
   tgeom.add_vert(P);
   return is_geom_on_surface_hull(tgeom, hull, eps);
}

bool is_point_inside_hull(const vec3d &P, const geom_if &hull, double eps)
{
   geom_v tgeom;
   tgeom.add_vert(P);
   return is_geom_inside_hull(tgeom, hull, eps);
}

bool is_point_fully_inside_hull(const vec3d &P, const geom_if &hull, double eps)
{
   geom_v tgeom;
   tgeom.add_vert(P);
   return is_geom_fully_inside_hull(tgeom, hull, eps);
}


// RK - Various find functions for geom

int find_vertex_by_coordinate(geom_if &geom, const vec3d &v, double eps)
{
   const vector<vec3d> &verts = geom.verts();      
   int v_idx = -1;
   for(unsigned int i=0; i<verts.size(); i++) {
      if (!compare(verts[i], v, eps)) {
         v_idx = i;
         break;
      }
   }
   return v_idx;
}

// elem could be face or another edge
bool edge_exists_in_elem(const vector<int> &elem, const vector<int> &edge)
{
   vector<int> edge1 = make_edge(edge[0], edge[1]);

   bool found = false;

   int sz = elem.size();
   for (int i=0;i<sz;i++) {
      vector<int> edge2 = make_edge(elem[i], elem[(i+1)%sz]);

      if (edge1 == edge2) {
         found = true;
         break;
      }
   }

   return found;
}

bool edge_exists_in_face(const vector<int> &face, const vector<int> &edge)
{
   return edge_exists_in_elem(face, edge);
}

bool are_edges_equal(const vector<int> &edge2, const vector<int> &edge1)
{
   return edge_exists_in_elem(edge2, edge1);
}

vector<int> find_faces_with_edge(const vector<vector<int> > &faces, const vector<int> &edge)
{
   vector<int> face_idxs;
   for (unsigned int i=0;i<faces.size();i++) {
      if (edge_exists_in_elem(faces[i], edge))
         face_idxs.push_back(i);
   }

   return face_idxs;
}

bool vertex_exists_in_elem(const vector<int> &elem, const int &v_idx)
{
   bool found = false;

   for (unsigned int i=0;i<elem.size();i++) {
      if (v_idx == elem[i]) {
         found = true;
         break;
      }
   }

   return found;
}

bool vertex_exists_in_face(const vector<int> &face, const int &v_idx)
{
   return vertex_exists_in_elem(face, v_idx);
}

bool vertex_exists_in_edge(const vector<int> &edge, const int &v_idx)
{
   return vertex_exists_in_elem(edge, v_idx);
}

vector<int> find_elems_with_vertex(const vector<vector<int> > &elems, const int &v_idx)
{
   vector<int> elem_idxs;
   for (unsigned int i=0;i<elems.size();i++) {
      if (vertex_exists_in_elem(elems[i], v_idx))
         elem_idxs.push_back(i);
   }

   return elem_idxs;
}

vector<int> find_faces_with_vertex(const vector<vector<int> > &faces, const int &v_idx)
{
   return find_elems_with_vertex(faces, v_idx);
}

vector<int> find_edges_with_vertex(const vector<vector<int> > &edges, const int &v_idx)
{
   return find_elems_with_vertex(edges, v_idx);
}

// find index of an edge in an edge list. -1 if not found
int find_edge_in_edge_list(const vector<vector<int> > &edges, const vector<int> &edge)
{
   int found = -1;
   for (unsigned int i=0;i<edges.size();i++) {
      if (edge_exists_in_elem(edges[i], edge)) {
         found = i;
         break;
      }
   }

   return found;
}

vector<vector<int> > find_unmatched_edges(col_geom_v &geom)
{
   const vector<vector<int> > &faces = geom.faces();
   vector<vector<int> > edges;

   // can't use get_impl_edges() here because we need to know the duplicates
   for (unsigned int i=0;i<faces.size();i++) {
      vector<int> face = faces[i];
      int sz = face.size();
      for (int j=0;j<sz;j++)
         edges.push_back(make_edge(face[j],face[(j+1)%sz]));
   }

   sort(edges.begin(), edges.end());

   vector<vector<int> > unmatched_edges;

   int sz = edges.size();
   for (int i=1;i<sz;) {
      if (edges[i-1] == edges[i])
         i+=2;
      else {
         unmatched_edges.push_back(edges[i-1]);
         if (i==sz-1)
            unmatched_edges.push_back(edges[i]);
         i++;
      }
   }

   return unmatched_edges;
}


// RK - the normals classes. xnormals for one normal. fnormals for the list of face normals.

// face normal
xnormal::xnormal(const geom_if &geom, const int &face_idx, vec3d C, double eps)
{
   const vector<int> &face = geom.faces()[face_idx];
   const vector<vec3d> &verts = geom.verts();

   normal = face_norm(verts, face);
   
   if (!C.is_set())
      C = centroid(verts);

   double D = vdot(verts[face[0]]-C, normal);

   // test case centroid is on the hemi's face
   if (double_eq(vdot(centroid(verts, face)-C, normal), 0.0, eps))
      direction = 0;
   else
      direction = double_compare(D, 0.0, eps);

}

// edge and vertex normals which have already been calculated from fnormals
xnormal::xnormal(const geom_if &geom, const vec3d &norm, const int &v_idx, vec3d C, double eps)
{
   const vector<vec3d> &verts = geom.verts();

   normal = norm;

   if (!C.is_set())
      C = centroid(verts);

   double D = vdot(verts[v_idx]-C, normal);

   direction = double_compare(D, 0.0, eps);
}

void fnormals::refresh(const geom_if &geom, vec3d C, double eps)
{
   ngeom = &geom;

   if (!C.is_set())
      vec3d C = centroid((*ngeom).verts());
   
   normals.clear();
   for(unsigned int i=0;i<geom.faces().size();i++)
      normals.push_back(xnormal(*ngeom, i, C, eps));
}

// private common code for averaging edge and vertex normals
vec3d fnormals::average_normals(vector<int> &face_idx, const string &average_pattern)
{
   vec3d norm;

   if (face_idx.size()) {
      vector<vec3d> f_normals;
      for (unsigned int i=0;i<face_idx.size();i++) {
         vec3d normal = normals[face_idx[i]].raw(); // default
         for (unsigned int j=0;j<average_pattern.size();j++) {
            if (average_pattern[j] == 'r') 
               normal = normals[face_idx[i]].raw();
            else
            if (average_pattern[j] == 'o') 
               normal = normals[face_idx[i]].outward();
            else
            if (average_pattern[j] == 'i') 
               normal = normals[face_idx[i]].inward();
            else
            if (average_pattern[j] == 'u') 
               normal = normal.unit();
         }
         f_normals.push_back(normal);
      }

      norm = centroid(f_normals);
   }

   return norm;
}

// the edge normal is centroid of all the face normals of which edge is a part of those faces
xnormal fnormals::edge_normal(const int &idx1, const int &idx2, const string &average_pattern, vec3d C, double eps)
{
   vector<int> edge = make_edge(idx1,idx2);

   vector<int> face_idx = find_faces_with_edge((*ngeom).faces(), edge);

   vec3d norm = average_normals(face_idx, average_pattern);

   return (norm.is_set() ? xnormal(*ngeom, norm, idx1, C, eps) : xnormal());
}

// the vertex normal is centroid of all the face normals of which vertex is a part of those faces
xnormal fnormals::vertex_normal(const int &idx, const string &average_pattern, vec3d C, double eps)
{
   vector<int> face_idx = find_faces_with_vertex((*ngeom).faces(), idx);

   vec3d norm = average_normals(face_idx, average_pattern);

   return (norm.is_set() ? xnormal(*ngeom, norm, idx, C, eps) : xnormal());
}


