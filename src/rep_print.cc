/*
   Copyright (c) 2003-2009, Adrian Rossiter

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
   Name: info_print.cc
   Description: information from OFF file - print functions 
   Project: Antiprism - http://www.antiprism.com
*/

#include <stdio.h>
#include <string.h>
#include <set>
#include <vector>
#include <utility>
#include "rep_print.h"

using std::set;
using std::vector;
using std::pair;

char *rep_printer::idx2s(char *buf, int idx, int elems_sz)
{ 
   buf[MSG_SZ-1] = 0;
   if(idx<elems_sz)
      snprintf(buf, MSG_SZ-1, "%d", idx);
   else 
      snprintf(buf, MSG_SZ-1, "E%d", idx-elems_sz);
   return buf;
}


void rep_printer::general_sec()
{
   char s1[MSG_SZ];
   fprintf(ofile, "[general]\n");
   fprintf(ofile, "num_verts = %d\n", num_verts());
   fprintf(ofile, "num_faces = %d\n", num_faces());
   fprintf(ofile, "num_edges = %d\n", num_edges());
   fprintf(ofile, "centroid = (%s)\n",v2s(s1, geom.centroid()));
   fprintf(ofile, "oriented = %s\n", is_oriented()?"yes":"no");
   fprintf(ofile, "orientable = %s\n", is_orientable()?"yes":"no");
   fprintf(ofile, "closed = %s\n", is_closed()?"yes":"no");
   fprintf(ofile, "num_parts = %d\n", num_parts());
   fprintf(ofile, "area = %s\n", d2s(s1, face_areas().sum));
   fprintf(ofile, "volume = ");
   if(is_closed())
      fprintf(ofile, "%s\n", d2s(s1, volume()));
   else
      fprintf(ofile, "n/a (not closed)\n");
fprintf(ofile, "\n");
}

void rep_printer::faces_sec()
{
   char s1[MSG_SZ];
   fprintf(ofile, "[faces]\n");
   fprintf(ofile, "num_faces = %d\n", num_faces());
   fprintf(ofile, "area = %s\n", d2s(s1, face_areas().sum));
   fprintf(ofile, "face_area_max = %s (%d)\n",
         d2s(s1, face_areas().max), face_areas().idx[elem_lims::IDX_MAX]);
   fprintf(ofile, "face_area_min = %s (%d)\n",
         d2s(s1, face_areas().min), face_areas().idx[elem_lims::IDX_MIN]);
   fprintf(ofile, "face_area_avg = %s\n",
         d2s(s1, face_areas().sum/num_faces()));
   fprintf(ofile, "volume = ");
   if(is_closed())
      fprintf(ofile, "%s\n", d2s(s1, volume()));
   else
      fprintf(ofile, "n/a (not closed)\n");
   fprintf(ofile, "vol2_by_area3 = ");
   if(is_closed())
      fprintf(ofile, "%s\n", d2s(s1, vol2_by_area3()));
   else
      fprintf(ofile, "n/a (not closed)\n");
   fprintf(ofile, "\n");
}


void rep_printer::angles_sec()
{
   char s1[MSG_SZ];
   fprintf(ofile, "[angles]\n");
   fprintf(ofile, "angle_max = %s\n",
         d2s(s1, rad2deg(angles().max)));
   fprintf(ofile, "angle_min = %s\n",
         d2s(s1, rad2deg(angles().min)));
   fprintf(ofile, "angle_avg = %s\n",
         d2s(s1, rad2deg(angles().sum/num_angles())));
   fprintf(ofile, "angle_defect = %s\n",
         d2s(s1, rad2deg(angle_defect())));
   fprintf(ofile, "\n");
}



void rep_printer::solid_angles_sec()
{
   char s1[MSG_SZ], s2[MSG_SZ];
   fprintf(ofile, "[solid_angles]\n");
   fprintf(ofile, "dihed_angle_max = %s (%d, %d)\n",
         d2s(s1, rad2deg(dihed_angles().max)),
         dihed_angles().idx[elem_lims::IDX_MAX], dihed_angles().idx[elem_lims::IDX_MAX2]);
   fprintf(ofile, "dihed_angle_min = %s (%d, %d)\n",
         d2s(s1, rad2deg(dihed_angles().min)),
         dihed_angles().idx[elem_lims::IDX_MIN],
         dihed_angles().idx[elem_lims::IDX_MIN2]);
   fprintf(ofile, "dihed_angle_flattest = %s (%d, %d)\n",
         d2s(s1, rad2deg(dihed_angles().mid)),
         dihed_angles().idx[elem_lims::IDX_MID],
         dihed_angles().idx[elem_lims::IDX_MID2]);
   fprintf(ofile, "solid_angle_max = %s [%sx4PI] (%d)\n",
         d2s(s1, solid_angles().max),
         d2s(s2, solid_angles().max/(4*M_PI)),
         solid_angles().idx[elem_lims::IDX_MAX]);
   fprintf(ofile, "solid_angle_min = %s [%sx4PI] (%d)\n",
         d2s(s1, solid_angles().min),
         d2s(s2, solid_angles().min/(4*M_PI)),
         solid_angles().idx[elem_lims::IDX_MIN]);
   fprintf(ofile, "\n");
}

void rep_printer::edges_sec()
{
   char s1[MSG_SZ];
   fprintf(ofile, "[edges]\n");
   fprintf(ofile, "num_edges = %d\n", num_edges());
   fprintf(ofile, "perimeter = %s\n",
         d2s(s1, edge_lengths().sum));
   fprintf(ofile, "edge_length_max = %s (%d,%d)\n",
         d2s(s1, edge_lengths().max),
         edge_lengths().idx[elem_lims::IDX_MAX],
         edge_lengths().idx[elem_lims::IDX_MAX2]);
   fprintf(ofile, "edge_length_min = %s (%d,%d)\n",
         d2s(s1, edge_lengths().min),
         edge_lengths().idx[elem_lims::IDX_MIN],
         edge_lengths().idx[elem_lims::IDX_MIN2]);
   fprintf(ofile, "edge_length_avg = %s\n",
         d2s(s1, edge_lengths().sum/num_edges()));
   fprintf(ofile, "\n");
}



void rep_printer::distances_sec()
{
   char s1[MSG_SZ];
   fprintf(ofile, "[distances]\n");
   fprintf(ofile, "given_center = (%s)\n",
         v2s(s1, center()));
   fprintf(ofile, "vert_min = %s (%d)\n",
         d2s(s1, vert_dists().min), vert_dists().idx[elem_lims::IDX_MIN]);
   fprintf(ofile, "vert_max = %s (%d)\n",
         d2s(s1, vert_dists().max), vert_dists().idx[elem_lims::IDX_MAX]);
   fprintf(ofile, "vert_avg = %s\n",
         d2s(s1, vert_dists().sum/num_verts()));
   fprintf(ofile, "face_min = %s (%d)\n", d2s(s1, face_dists().min),
          face_dists().idx[elem_lims::IDX_MIN]);
   fprintf(ofile, "face_max = %s (%d)\n", d2s(s1, face_dists().max),
         face_dists().idx[elem_lims::IDX_MAX]);
   fprintf(ofile, "face_avg = %s\n",
         d2s(s1, face_dists().sum/num_faces()));
   fprintf(ofile, "edge_min = %s (%d,%d)\n",
         d2s(s1, edge_dists().min),
         edge_dists().idx[elem_lims::IDX_MIN],
         edge_dists().idx[elem_lims::IDX_MIN2]);
   fprintf(ofile, "edge_max = %s (%d,%d)\n",
         d2s(s1, edge_dists().max),
         edge_dists().idx[elem_lims::IDX_MAX],
         edge_dists().idx[elem_lims::IDX_MAX2]);
   fprintf(ofile, "edge_avg = %s\n",
         d2s(s1, edge_dists().sum/num_edges()));
   fprintf(ofile, "\n");
}


void rep_printer::symmetry()
{
   char s1[MSG_SZ];
   fprintf(ofile, "[symmetry]\n");
   fprintf(ofile, "type = %s\n", get_symmetry_type_name().c_str());
   fprintf(ofile, "alignment_to_std =");
   mat3d m =  get_symmetry_alignment_to_std();
   for(int i=0; i<12; i++)
      fprintf(ofile, "%c%s", (i==0)?' ':',', d2s(s1, m[i]));
   fprintf(ofile, "\n\n");
   fprintf(ofile, "[symmetry_axes]\n");
   const set<sch_axis> &axes = get_symmetry_axes();
   set<sch_axis>::const_iterator ai;
   for(ai=axes.begin(); ai!=axes.end(); ++ai) {
      fprintf(ofile, "%s,",
         sch_sym(ai->get_sym_type(), ai->get_nfold()).get_symbol().c_str());
      if(ai->get_axis().is_set())
         fprintf(ofile, "%s,", v2s(s1, ai->get_axis()));
      else
         fprintf(ofile, "n/a,");
      if(ai->get_perp().is_set())
         fprintf(ofile, "%s,", v2s(s1, ai->get_perp()));
      else
         fprintf(ofile, "n/a");
      fprintf(ofile, "\n");
   }
   fprintf(ofile, "\n");
}


void rep_printer::edge_lengths_cnts()
{
   char s1[MSG_SZ];
   fprintf(ofile, "[edge_lengths_cnts]\n");
   map<double, int, ang_less> edge_lengths = get_e_lengths();
   map<double, int, ang_less>::iterator ei;
   for(ei=edge_lengths.begin(); ei!=edge_lengths.end(); ++ei)
      fprintf(ofile, "%s = %d\n",
            d2s(s1, ei->first), ei->second);
   fprintf(ofile, "\n");
}


void rep_printer::dihedral_angles_cnts()
{
   char s1[MSG_SZ];
   fprintf(ofile, "[dihedral_angles_cnts]\n");
   map<double, int, ang_less> &dihedrals = get_dihedral_angles();
   map<double, int, ang_less>::iterator di;
   for(di=dihedrals.begin(); di!=dihedrals.end(); ++di)
      fprintf(ofile, "%s = %d\n",
            d2s(s1, rad2deg(di->first)),
            di->second);
   fprintf(ofile, "\n");
}

void rep_printer::solid_angles_cnts()
{
   char s1[MSG_SZ];
   fprintf(ofile, "[solid_angles_cnts]\n");
   map<double, int, ang_less> &solid_angs = get_solid_angles();
   map<double, int, ang_less>::iterator si;
   for(si=solid_angs.begin(); si!=solid_angs.end(); ++si)
      fprintf(ofile, "%s = %d\n",
            d2s(s1, si->first), si->second);
   fprintf(ofile, "\n");
}



void rep_printer::vert_order_cnts()
{
   fprintf(ofile, "[vert_order_cnts]\n");
   map<int, int>::iterator mi;
   map<int, int> cnts;
   const vector<vector<int> > &v_cons = get_vert_cons();
   for(unsigned int i=0; i<v_cons.size(); i++) {
      mi = cnts.find(v_cons[i].size());
      if(mi == cnts.end())
         cnts[v_cons[i].size()]=1;
      else
         mi->second += 1;
   }
   for(mi=cnts.begin(); mi!=cnts.end(); ++mi)
      fprintf(ofile, "%d = %d\n", mi->first, mi->second);
   fprintf(ofile, "\n");
}

void rep_printer::face_sides_cnts()
{
   fprintf(ofile, "[face_sides_cnts]\n");
   map<int, int>::iterator mi;
   map<int, int> cnts;
   for(unsigned int i=0; i<geom.faces().size(); i++) {
      mi = cnts.find(geom.faces(i).size());
      if(mi == cnts.end())
         cnts[geom.faces(i).size()]=1;
      else
         mi->second += 1;
   }
   for(mi=cnts.begin(); mi!=cnts.end(); ++mi)
      fprintf(ofile, "%d = %d\n", mi->first, mi->second);
   fprintf(ofile, "\n");
}

void rep_printer::face_angles_cnts()
{
   char s1[MSG_SZ];
   fprintf(ofile, "[face_angles_cnts]\n");
   map<vector<double>, int, ang_vect_less> &face_angs = get_face_angles();
   map<vector<double>, int, ang_vect_less>::iterator fi;
   for(fi=face_angs.begin(); fi!=face_angs.end(); ++fi) {
      for(unsigned int i=0; i<fi->first.size(); i++)
         fprintf(ofile, "%s%s",
               d2s(s1, rad2deg(fi->first[i])),
               (i<fi->first.size()-1)?",":"");
      fprintf(ofile, " = %d\n", fi->second);
   }
   fprintf(ofile, "\n");
}


void rep_printer::v_index(int v_idx)
{
   char str[MSG_SZ];
   fprintf(ofile, "%s", vidx2s(str, v_idx));
}


void rep_printer::v_coords(int v_idx)
{
   char str[MSG_SZ];
   fprintf(ofile, "%s", v2s(str, geom.verts(v_idx)));
}
   
void rep_printer::v_neighbours(int v_idx)
{
   char str[MSG_SZ];
   const vector<int> &vcons = get_vert_cons()[v_idx];
   for(unsigned int i=0; i<vcons.size(); i++)
      fprintf(ofile, "%s%s", vidx2s(str, vcons[i]), (i<vcons.size()-1)?" ":"");
}
   
void rep_printer::v_face_idxs(int v_idx)
{
   char str[MSG_SZ];
   geom_v &dual = get_dual();
   const vector<int> &fcons = (*dual.get_faces())[v_idx];
   for(unsigned int i=0; i<fcons.size(); i++)
      fprintf(ofile, "%s%s", fidx2s(str, fcons[i]), (i<fcons.size()-1)?" ":"");
}
   
void rep_printer::v_solid_angle(int v_idx)
{
   char str[MSG_SZ];
   fprintf(ofile, "%s", d2s(str, get_vertex_angles()[v_idx]));
}
   
void rep_printer::v_order(int v_idx)
{
   fprintf(ofile, "%lu", (unsigned long)get_vert_cons()[v_idx].size());
}
   
void rep_printer::v_distance(int v_idx)
{
   char str[MSG_SZ];
   fprintf(ofile, "%s", d2s(str, (geom.verts(v_idx)-center()).mag()) );
}

void rep_printer::v_angles(int v_idx)
{
   char str[MSG_SZ];
   geom_v &dual = get_dual();
   const vector<int> &fcons = (*dual.get_faces())[v_idx];
   for(unsigned int i=0; i<fcons.size(); i++) {
      pair<int, int> vf_pr(v_idx, fcons[i]);
      fprintf(ofile, "%s%s", d2s(str, rad2deg(get_vertex_plane_angs()[vf_pr])),
            (i<fcons.size()-1)?" ":"");
   }
}


void rep_printer::e_index(int e_idx)
{
   char str[MSG_SZ];
   fprintf(ofile, "%s", eidx2s(str, e_idx));
}

void rep_printer::e_vert_idxs(int e_idx)
{
   char str[MSG_SZ], str2[MSG_SZ];
   vector<int> edge = geom.edges(e_idx);
   fprintf(ofile, "%s %s", vidx2s(str, edge[0]), vidx2s(str2, edge[1]));   
}
   
void rep_printer::e_face_idxs(int e_idx)
{
   char str[MSG_SZ], str2[MSG_SZ];
   vector<int> edge = geom.edges(e_idx);
   map<vector<int>, vector<int> >::iterator ei =
      get_edge_face_pairs().find(edge);
   if(ei != get_edge_face_pairs().end())
      fprintf(ofile, "%s %s", fidx2s(str, ei->second[0]),
                              fidx2s(str2,ei->second[1]));
   else
      fprintf(ofile, "-1 -1");
}
   
void rep_printer::e_dihedral_angle(int e_idx)
{
   char str[MSG_SZ];
   fprintf(ofile, "%s", d2s(str, rad2deg(get_edge_dihedrals()[e_idx])));
}
   
void rep_printer::e_central_angle(int e_idx)
{
   vec3d v0 = geom.verts(geom.edges(e_idx, 0)) - center();
   vec3d v1 = geom.verts(geom.edges(e_idx, 1)) - center();
   char str[MSG_SZ];
   fprintf(ofile, "%s", d2s(str, rad2deg(acos(vdot(v0.unit(), v1.unit())))));
}
   
void rep_printer::e_distance(int e_idx)
{
   vec3d v0 = geom.verts(geom.edges(e_idx, 0)) - center();
   vec3d v1 = geom.verts(geom.edges(e_idx, 1)) - center();
   double dist = (nearest_point(center(), v0, v1)-center()).mag();
   char str[MSG_SZ];
   fprintf(ofile, "%s", d2s(str, dist));
}
   
void rep_printer::e_direction(int e_idx)
{
   vec3d v0 = geom.verts(geom.edges(e_idx, 0));
   vec3d v1 = geom.verts(geom.edges(e_idx, 1));
   char str[MSG_SZ];
   fprintf(ofile, "%s", v2s(str, (v1-v0).unit()));
}

void rep_printer::e_length(int e_idx)
{
   vec3d v0 = geom.verts(geom.edges(e_idx, 0));
   vec3d v1 = geom.verts(geom.edges(e_idx, 1));
   char str[MSG_SZ];
   fprintf(ofile, "%s", d2s(str, (v1-v0).mag()));
}


void rep_printer::f_index(int f_idx)
{
   char str[MSG_SZ];
   fprintf(ofile, "%s", fidx2s(str, f_idx));
}

void rep_printer::f_vert_idxs(int f_idx)
{
   char str[MSG_SZ];
   const vector<int> &face = geom.faces(f_idx);
   for(unsigned int i=0; i<face.size(); i++)
      fprintf(ofile, "%s%s", vidx2s(str, face[i]), (i<face.size()-1)?" ":"");
}
   
void rep_printer::f_neighbours(int f_idx)
{
   char str[MSG_SZ];
   map<vector<int>, vector<int> >::iterator ei;
   vector<int> edge(2);
   const vector<int> &face = geom.faces(f_idx);
   unsigned int sz = face.size();
   for(unsigned int i=0; i<sz; i++) {
      ei = get_edge_face_pairs().find(make_edge(face[i], face[(i+1)%sz]));
      int neigh = (ei->second[0]!=f_idx)? ei->second[0]: ei->second[1];
      fprintf(ofile, "%s%s", fidx2s(str, neigh), (i<sz-1)?" ":"");
   }
}
   
void rep_printer::f_normal(int f_idx)
{
   char str[MSG_SZ];
   fprintf(ofile, "%s", v2s(str, geom.face_norm(f_idx).unit()));
}
   
void rep_printer::f_angles(int f_idx)
{
   vector<double> angs;
   face_angles_lengths(f_idx, angs);
   char str[MSG_SZ];
   for(unsigned int i=0; i<angs.size(); i++)
      fprintf(ofile, "%s%s", d2s(str, rad2deg(angs[i])),
            (i<angs.size()-1)?" ":"");
}
   
void rep_printer::f_sides(int f_idx)
{
   fprintf(ofile, "%lu", (unsigned long)geom.faces(f_idx).size());
}
   
void rep_printer::f_distance(int f_idx)
{
   double dist = (nearest_point(center(), geom.verts(), geom.faces(f_idx))
         -center()).mag();
   char str[MSG_SZ];
   fprintf(ofile, "%s", d2s(str, dist));
}
   
void rep_printer::f_area(int f_idx)
{
   char str[MSG_SZ];
   fprintf(ofile, "%s", d2s(str, get_f_areas()[f_idx]));
}
   
void rep_printer::f_centroid(int f_idx)
{
   char str[MSG_SZ];
   fprintf(ofile, "%s", v2s(str, geom.face_cent(f_idx)));
}

void rep_printer::f_lengths(int f_idx)
{
   char str[MSG_SZ];
   const vector<vec3d> &verts = geom.verts();
   const vector<int> &face = geom.faces(f_idx);
   unsigned int sz = face.size();
   for(unsigned int i=0; i<sz; i++)
      fprintf(ofile, "%s%s",
            d2s(str,(verts[face[i]]-verts[face[(i+1)%sz]]).mag()),
            (i<face.size()-1)?" ":"");
}

void rep_printer::f_max_nonplanar(int f_idx)
{
   char str[MSG_SZ];
   fprintf(ofile, "%s", d2s(str, get_f_max_nonplanars()[f_idx]));
}


