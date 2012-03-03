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
   Name: transforms.h
   Description: transforming one polyhedron into another
   Project: Antiprism - http://www.antiprism.com
*/



#ifndef TRANSFORMS_H
#define TRANSFORMS_H

#include <map>
#include "geom.h"
#include "coloring.h"
#include "symmetry.h"

using std::pair;

class geom_info;


void sym_repeat(geom_if &geom, const geom_if &part, const t_set &ts,
      char col_part_elems=ELEM_NONE, coloring *clrngs=0);
bool sym_repeat(geom_if &geom, const geom_if &part, const sch_sym &sym,
      char col_part_elems=ELEM_NONE, coloring *clrngs=0);

int get_voronoi_cells(geom_if &geom, vector<col_geom_v> &cells,
      string qh_args="", char *errmsg=0);



void get_pol_recip_verts(const geom_if &geom, geom_if &dual, double recip_rad, vec3d centre, double inf = 1e15);

void get_dual(const geom_if &geom, geom_if &dual, double recip_rad=0, vec3d centre=vec3d(0, 0, 0), double inf = 1e20);
void add_extra_ideal_elems(geom_if &geom, vec3d centre, double inf);
void limit_distance(geom_if &geom, vec3d centre, double inf);

void edges_to_faces(const geom_if &geom, geom_if &egeom, bool face_verts=false);

void truncate_verts(geom_if &geom, vector<int> &v_idxs, double ratio, geom_info *info=0);
void truncate_verts(geom_if &geom, double ratio, int order=0,geom_info *info=0);


//delete_elem contains v, e, or f or "" to just sort.
void sort_merge_elems(geom_if &geom, const string &merge_elems, vector<map<int, set<int> > > *equiv_elems=0, double eps=epsilon);
void sort_merge_elems(geom_if &geom, const string &merge_elems, const int &blend_type, double eps=epsilon);
void sort_merge_elems(geom_if &geom, const string &merge_elems, double eps=epsilon);
bool check_congruence(const geom_if &geom1, const geom_if &geom2, vector<map<int, set<int> > > *equiv_elems=0, double eps=epsilon);
void get_congruence_maps(const geom_if &geom, mat3d trans,
      vector<vector<int > > &elem_maps, double eps=epsilon);


void canonicalize_mm(geom_if &geom, double edge_factor, double plane_factor, int n, int divergence_test, int rep_count, bool planar_only, double eps=epsilon);
void canonicalize_cn(geom_if &geom, int n, char method, int divergence_test, int rep_count, double eps=epsilon);

bool close_poly_basic(geom_if &geom);
bool face_bond(geom_if &geom, geom_if &bgeom, int f=0, int b_f=0, int off=0, bool merge=true);
bool face_bond_direct(geom_if &geom, geom_if &bgeom, int f=0, int b_f=0, bool merge=true);

void transform_and_repeat(geom_if &geom, string sym_to, string sym_from, mat3d pos=mat3d());

#endif // TRANSFORMS_H
