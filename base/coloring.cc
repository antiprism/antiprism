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

/* \file coloring.cc
 * \brief classes to color all elements of a type.
 */

#include <limits.h>
#include <string.h>
#include <algorithm>
#include "rand_gen.h"
#include "prop_col.h"
#include "coloring.h"
#include "transforms.h"
#include "bbox.h"
#include "utils.h"
#include "math_utils.h"
#include "info.h"

coloring::coloring(col_geom_v *geo): geom(geo), cycle_msecs(0)
{
}

coloring::~coloring()
{
}

coloring::coloring(const coloring &clrng): color_map_multi(clrng),
   geom(clrng.geom), cycle_msecs(clrng.cycle_msecs)
{
}

coloring &coloring::operator=(const coloring &clrng)
{
   if(this!=&clrng) {
      color_map_multi::operator=(clrng);
      geom = clrng.geom;
      cycle_msecs = clrng.cycle_msecs;
   }

   return *this;
}

 
void coloring::cycle_map_cols()
{
   set_shift(get_shift()+1);
}
       
 
void coloring::set_all_idx_to_val(map<int, col_val> &cols)
{
   map<int, col_val>::iterator mi;
   for(mi=cols.begin(); mi!=cols.end(); mi++)
      if(mi->second.is_idx())
         mi->second = get_col(mi->second.get_idx());
}


inline double fract(double rng[], double frac)
{
   return fmod(rng[0] + (rng[1]-rng[0])*frac, 1+epsilon);
}


int coloring::y_gradient(vec3d vec, vec3d cent, double height, int def_sz)
{
   int sz = def_sz;
   const vector<color_map *> &cmaps = get_cmaps();
   if(cmaps.size()>0) {
      if(cmaps[0]->effective_size()>0)
         sz = cmaps[0]->effective_size();
      else
         sz = INT_MAX;
   }

   return (int)floor(sz * (0.5*height+(vec-cent)[1])/(height+epsilon));
}


void coloring::setup_lights(col_geom_v &lts)
{
   if(lts.verts().size()==0) {
      lts.add_col_vert(vec3d(1,0,0), vec3d(1,0,0));
      lts.add_col_vert(vec3d(0,1,0), vec3d(0,1,0));
      lts.add_col_vert(vec3d(0,0,1), vec3d(0,0,1));
      lts.add_col_vert(vec3d(-1,0,0), vec3d(0,1,1));
      lts.add_col_vert(vec3d(0,-1,0), vec3d(1,0,1));
      lts.add_col_vert(vec3d(0,0,-1), vec3d(1,1,0));
   }
   else
      for(unsigned int l=0; l<lts.verts().size(); l++)
         lts.raw_verts()[l].to_unit();
}



col_val coloring::light(vec3d vec, col_geom_v &lts)
{
   vec3d col_sum(0,0,0);
   vec.to_unit();
   double dot;
   for(unsigned int l=0; l<lts.verts().size(); l++)
      if((dot = vdot(vec, lts.verts(l)))>0)
         col_sum += dot * lts.get_v_col(l).get_vec3d();
   
   for(int j=0; j<3; j++)
      if(col_sum[j]>1)
         col_sum[j] = 1;

   return col_val(col_sum);
}





void coloring::v_apply_cmap()
{
   set_all_idx_to_val(get_geom()->raw_vert_cols());
}


void coloring::v_one_col(col_val col)
{
   for(unsigned int i=0; i<get_geom()->verts().size(); i++)
      get_geom()->set_v_col(i, col);
}


void coloring::v_unique(bool apply_map)
{
   for(unsigned int i=0; i<get_geom()->verts().size(); i++) {
      if(apply_map)
         get_geom()->set_v_col(i, get_col(i));
      else
         get_geom()->set_v_col(i, i);
   }
}


void coloring::v_sets(const vector<set<int> > &equivs, bool apply_map)
{
   for(unsigned int i=0; i<equivs.size(); i++) {
      for(set<int>::iterator si=equivs[i].begin(); si!=equivs[i].end(); ++si) {
         if(apply_map)
            get_geom()->set_v_col(*si, get_col(i));
         else
            get_geom()->set_v_col(*si, i);
      }
   }
}


void coloring::v_proper(bool apply_map)
{
   prop_color prop(geom->verts().size());

   const vector<vector<int> > &faces = geom->faces();
   for(unsigned int i=0; i<faces.size(); ++i)
      for(unsigned int j=0; j<faces[i].size(); ++j)
         prop.set_adj(faces[i][j], faces[i][(j+1)%faces[i].size()]);
   
   prop.find_colors();
   for(unsigned int i=0; i<get_geom()->verts().size(); i++)
      if(apply_map)
         get_geom()->set_v_col(i, get_col(prop.get_color(i)));
      else
         get_geom()->set_v_col(i,prop.get_color(i));
}


void coloring::v_order(bool apply_map)
{
   vector<int> f_cnt(get_geom()->verts().size());
   for(unsigned int i=0; i<get_geom()->faces().size(); i++)
      for(unsigned int j=0; j<get_geom()->faces(i).size(); j++)
         f_cnt[get_geom()->faces(i, j)]++;
   for(unsigned int i=0; i<get_geom()->verts().size(); i++)
      if(apply_map)
         get_geom()->set_v_col(i, get_col(f_cnt[i]));
      else
         get_geom()->set_v_col(i, f_cnt[i]);
}


void coloring::v_position(bool apply_map)
{
   bound_box bb(get_geom()->verts());
   vec3d cent = bb.get_centre();
   double height = bb.get_max()[1] - bb.get_min()[1];
   for(unsigned int i=0; i<get_geom()->verts().size(); i++) {
      int idx = y_gradient(get_geom()->verts(i), cent, height);
      if(apply_map)
         get_geom()->set_v_col(i, get_col(idx));
      else
         get_geom()->set_v_col(i, idx);
   }
}

void coloring::v_avg_angle(bool apply_map)
{
   geom_info info(*get_geom());
   col_geom_v tmp;
   tmp.raw_verts() = get_geom()->verts();
   for(int i=0; i<info.num_verts(); i++)
      tmp.add_face(info.get_vert_figs()[i].size()? info.get_vert_figs()[i][0]
                                                 : vector<int>() );
   col_geom_v *orig_geom = get_geom();
   set_geom(&tmp);
   f_avg_angle(apply_map);
   orig_geom->raw_vert_cols() = tmp.face_cols();
   set_geom(orig_geom);
}



void coloring::v_lights(col_geom_v lts)
{
   setup_lights(lts);
   vec3d cent = get_geom()->centroid();
   for(unsigned int i=0; i<get_geom()->verts().size(); i++)
      get_geom()->set_v_col(i, light(get_geom()->verts(i) - cent, lts));
}

void coloring::face_edge_color(const vector<vector<int> > &elems,
      const map<int, col_val> &cmap)
{
   vector<vector<int> > v_elems(get_geom()->verts().size());
   for(unsigned int i=0; i<elems.size(); ++i)
      for(unsigned int j=0; j<elems[i].size(); ++j)
         v_elems[elems[i][j]].push_back(i);
   
   for(unsigned int i=0; i<v_elems.size(); ++i) {
      vec4d col(0,0,0,0);
      int val_cnt = 0;
      int first_idx = -1;
      for(unsigned int j=0; j<v_elems[i].size(); ++j) {
         col_val ecol = col_geom::get_col(cmap, v_elems[i][j]);
         if(ecol.is_val()) {
            col += ecol.get_vec4d();
            val_cnt++;
         }
         else if(first_idx==-1 && ecol.is_idx())
            first_idx = ecol.get_idx();
      }
      if(val_cnt)
         get_geom()->set_v_col(i, col_val(col/double(val_cnt)));
      else if(first_idx!=-1)
         get_geom()->set_v_col(i, col_val(first_idx));
   }
}


void coloring::v_face_color()
{
   face_edge_color(get_geom()->faces(), get_geom()->face_cols());
}


void coloring::v_edge_color()
{
   face_edge_color(get_geom()->edges(), get_geom()->edge_cols());
}





void coloring::f_apply_cmap()
{
   set_all_idx_to_val(get_geom()->raw_face_cols());
}


void coloring::f_one_col(col_val col)
{
   for(unsigned int i=0; i<get_geom()->faces().size(); i++)
      get_geom()->set_f_col(i, col);
}


void coloring::f_sets(const vector<set<int> > &equivs, bool apply_map)
{
   for(unsigned int i=0; i<equivs.size(); i++) {
      col_val col(i);
      for(set<int>::iterator si=equivs[i].begin(); si!=equivs[i].end(); ++si) {
         if(apply_map)
            get_geom()->set_f_col(*si, get_col(i));
         else
            get_geom()->set_f_col(*si, i);
      }
   }
}


void coloring::f_unique(bool apply_map)
{
   for(unsigned int i=0; i<get_geom()->faces().size(); i++) {
      if(apply_map)
         get_geom()->set_f_col(i, get_col(i));
      else
         get_geom()->set_f_col(i, i);
   }
}


void coloring::f_proper(bool apply_map)
{
   prop_color prop(geom->faces().size());

   map<vector<int>, vector<int> > ef_prs;
   geom->get_edge_face_pairs(ef_prs, false);
   map<vector<int>, vector<int> >::const_iterator mi;
   for(mi=ef_prs.begin(); mi!=ef_prs.end(); ++mi) {
      const vector<int> &f_idxs = mi->second; 
      for(unsigned int i=0; i<f_idxs.size(); ++i)
         for(unsigned int j=i+1; j<f_idxs.size(); ++j)
            prop.set_adj(f_idxs[i], f_idxs[j]);
   }
   
   prop.find_colors();
   for(unsigned int i=0; i<get_geom()->faces().size(); i++)
      if(apply_map)
         get_geom()->set_f_col(i, get_col(prop.get_color(i)));
      else
         get_geom()->set_f_col(i,prop.get_color(i));
}


void coloring::f_sides(bool apply_map)
{
   for(unsigned int i=0; i<get_geom()->faces().size(); i++)
      if(apply_map)
         get_geom()->set_f_col(i, get_col(get_geom()->faces(i).size()));
      else
         get_geom()->set_f_col(i, get_geom()->faces(i).size());
}


void coloring::f_avg_angle(bool apply_map)
{
   geom_info info(*get_geom());
   int faces_sz = get_geom()->faces().size();
   for(int i=0; i<faces_sz; i++) {
      vector<double> f_angs;
      info.face_angles_lengths(i, f_angs);
      double ang_sum = 0;
      for(unsigned int j=0; j<f_angs.size(); j++)
         ang_sum += f_angs[j];
      // set invalid faces to have impossible angle 181
      int idx = f_angs.size() ? (int)(rad2deg(ang_sum/f_angs.size())+0.5) : 181;
      if(apply_map)
         get_geom()->set_f_col(i, get_col(idx));
      else
         get_geom()->set_f_col(i, idx);
   }
}


void coloring::f_parts(bool apply_map)
{
   vector<vector<int> > parts;
   geom_v gtmp = *get_geom();
   gtmp.orient(&parts);
   for(unsigned int i=0; i<parts.size(); i++)
      for(unsigned int j=0; j<parts[i].size(); j++)
         if(apply_map)
            get_geom()->set_f_col(parts[i][j], get_col(i));
         else
            get_geom()->set_f_col(parts[i][j], i);
}

void coloring::f_normal(bool apply_map)
{
   for(unsigned int i=0; i<get_geom()->faces().size(); i++) {
      int idx = y_gradient(get_geom()->face_norm(i).unit());
      if(apply_map)
         get_geom()->set_f_col(i, get_col(idx));
      else
         get_geom()->set_f_col(i, idx);
   }
}


void coloring::f_centroid(bool apply_map)
{
   bound_box bb(get_geom()->verts());
   vec3d cent = bb.get_centre();
   double height = bb.get_max()[1] - bb.get_min()[1];
   for(unsigned int i=0; i<get_geom()->faces().size(); i++) {
      int idx = y_gradient(get_geom()->face_cent(i), cent, height);
      if(apply_map)
         get_geom()->set_f_col(i, get_col(idx));
      else
         get_geom()->set_f_col(i, idx);
   }
}

void coloring::f_lights(col_geom_v lts)
{
   setup_lights(lts);
   for(unsigned int i=0; i<get_geom()->faces().size(); i++)
      get_geom()->set_f_col(i, light(get_geom()->face_norm(i), lts));
}

void coloring::f_lights2(col_geom_v lts)
{
   setup_lights(lts);
   for(unsigned int i=0; i<get_geom()->faces().size(); i++)
      get_geom()->set_f_col(i, light(get_geom()->face_cent(i), lts));
   
}




void coloring::e_apply_cmap()
{ 
   set_all_idx_to_val(get_geom()->raw_edge_cols());
}

void coloring::e_one_col(col_val col)
{
   for(unsigned int i=0; i<get_geom()->edges().size(); i++)
      get_geom()->set_e_col(i, col);
}


void coloring::e_sets(const vector<set<int> > &equivs, bool apply_map)
{
   for(unsigned int i=0; i<equivs.size(); i++) {
      for(set<int>::iterator si=equivs[i].begin(); si!=equivs[i].end(); ++si) {
         if(apply_map)
            get_geom()->set_e_col(*si, get_col(i));
         else
            get_geom()->set_e_col(*si, i);
      }
   }
}


void coloring::e_unique(bool apply_map)
{
   for(unsigned int i=0; i<get_geom()->edges().size(); i++) {
      if(apply_map)
         get_geom()->set_e_col(i, get_col(i));
      else
         get_geom()->set_e_col(i, i);
   }
}


void coloring::e_proper(bool apply_map)
{
   map<vector<int>, vector<int> > edges;
   geom->get_edge_face_pairs(edges, false);
   map<vector<int>, vector<int> >::iterator mi;
   int idx=0;
   for(mi=edges.begin(); mi!=edges.end(); ++mi)
      mi->second[0] = idx++;   // assign index numbers to the impl edges
   prop_color prop(idx);       // idx is total number of impl edges


   const vector<vector<int> > &faces = geom->faces();
   vector<int> e_next(2);
   for(unsigned int i=0; i<faces.size(); ++i) {
      const vector<int> &face = faces[i];
      int f_sz = face.size();
      for(int j=0; j<f_sz; ++j) {
         vector<int> e = make_edge(face[j], face[(j+1)%f_sz]);
         vector<int> e_next = make_edge(face[(j+1)%f_sz], face[(j+2)%f_sz]);
         // An edge is adjacent to the edge that follows it on a face
         prop.set_adj( edges[e][0], edges[e_next][0]);
      }
   }
   
   prop.find_colors();
   for(mi=edges.begin(); mi!=edges.end(); ++mi) {
      int e_idx = edges[mi->first][0];
      int col_idx = prop.get_color(e_idx);
      col_val col = (apply_map) ? get_col(col_idx) : col_val(col_idx);
      get_geom()->add_col_edge(mi->first[0], mi->first[1], col);
   }
}



void coloring::e_face_color()
{
   const vector<vector<int> > &faces = get_geom()->faces();
   vector<vector<int> > efaces(get_geom()->edges().size());
   for(unsigned int i=0; i<faces.size(); ++i) {
      for(unsigned int j=0; j<faces[i].size(); ++j) {
         vector<int> edge =
            make_edge(faces[i][j],faces[i][(j+1)%faces[i].size()]);
         vector<vector<int> >::const_iterator ei = get_geom()->edges().begin();
         while((ei=find(ei, get_geom()->edges().end(), edge)) !=
               get_geom()->edges().end()) {
            efaces[ei-get_geom()->edges().begin()].push_back(i);
            ++ei;
         }
      }
   }
   
   for(unsigned int i=0; i<efaces.size(); ++i) {
      vec4d col(0,0,0,0);
      int val_cnt = 0;
      int first_idx = -1;
      for(unsigned int j=0; j<efaces[i].size(); ++j) {
         col_val fcol = get_geom()->get_f_col(efaces[i][j]);
         if(fcol.is_val()) {
            col += fcol.get_vec4d();
            val_cnt++;
         }
         else if(first_idx==-1 && fcol.is_idx())
            first_idx = fcol.get_idx();
      }
      if(val_cnt)
         get_geom()->set_e_col(i, col_val(col/double(val_cnt)));
      else if(first_idx!=-1)
         get_geom()->set_e_col(i, col_val(first_idx));
   }
}


void coloring::edge_color_and_branch(int idx, int part, bool apply_map,
      vector<vector<int> > &vcons, vector<bool> &seen)
{
   if(seen[idx])
      return;
   else
      seen[idx] = true;

   for(unsigned int i=0; i<vcons[idx].size(); i++) {
      int next_idx = vcons[idx][i];
      if(idx == next_idx)
         continue;
      vector<int> edge(2);
      edge[0] = idx;
      edge[1] = next_idx;
      int e_idx = get_geom()->add_edge(edge);
      if(apply_map)
         get_geom()->set_e_col(e_idx, get_col(part));
      else
         get_geom()->set_e_col(e_idx, part);
      edge_color_and_branch(next_idx, part, apply_map, vcons, seen);
   }
}



void coloring::e_parts(bool apply_map)
{
   vector<vector<int> > vcons(get_geom()->get_verts()->size(), vector<int>());
   const vector<vector<int> > &edges = get_geom()->edges();
   for(unsigned int i=0; i<edges.size(); i++) {
      vcons[edges[i][0]].push_back(edges[i][1]);
      vcons[edges[i][1]].push_back(edges[i][0]);
   }

   int part=0;
   vector<bool> seen(vcons.size(), false);
   for(unsigned int i=0; i<vcons.size(); i++)
      if(!seen[i])
         edge_color_and_branch(i, part++, apply_map, vcons, seen);
}


void coloring::e_direction(bool apply_map)
{
   for(unsigned int i=0; i<get_geom()->edges().size(); i++) {
      vec3d v = 2.0*(get_geom()->edge_vec(i)).unit();
      if(v[1]<0)
         v = -v;
      v -= vec3d::Y;  // put v[1] in the range -1.0 to 1.0;
      int idx = y_gradient(v);
      if(apply_map)
         get_geom()->set_e_col(i, get_col(idx));
      else
         get_geom()->set_e_col(i, idx);
   }
}

void coloring::e_mid_point(bool apply_map)
{
   bound_box bb(get_geom()->verts());
   vec3d cent = bb.get_centre();
   double height = bb.get_max()[1] - bb.get_min()[1];
   for(unsigned int i=0; i<get_geom()->edges().size(); i++) {
      int idx = y_gradient(get_geom()->edge_cent(i), cent, height);
      if(apply_map)
         get_geom()->set_e_col(i, get_col(idx));
      else
         get_geom()->set_e_col(i, idx);
   }
}



void coloring::e_lights(col_geom_v lts)
{
   setup_lights(lts);
   vec3d cent = get_geom()->centroid();
   for(unsigned int i=0; i<get_geom()->edges().size(); i++)
      get_geom()->set_e_col(i, light(get_geom()->edge_nearpt(i,cent)-cent,lts));
}


static bool get_cycle_rate(const char *str, double *cps)
{
   size_t len = strlen(str);
   if(len>=3 && str[len-2]=='h' && str[len-1]=='z') {
      char str_copy[MSG_SZ];
      strncpy(str_copy, str, len-2);
      str_copy[len-2] = '\0';
      double cycs;
      if(read_double(str_copy, &cycs) && cycs>=0.0) {
         *cps = cycs;
         return true;
      }
   }

   return false;
}




bool read_colorings(coloring clrngs[], const char *line, char *errmsg)
{
   if(errmsg)
      *errmsg = '\0';

   char line_copy[MSG_SZ];
   strncpy(line_copy, line, MSG_SZ);
   line_copy[MSG_SZ-1] = '\0';

   vector<char *> parts;
   int parts_sz = split_line(line_copy, parts, ",");

   char errmsg2[MSG_SZ];
      *errmsg2 = '\0';
   vector<char *> map_names;
   coloring clrng;
   unsigned int conv_elems = 0;
   
   for(int i=0; i<parts_sz; i++) {
      color_map *col_map = init_color_map(parts[i], errmsg2);
      double cps;
      if(get_cycle_rate(parts[i], &cps)) {
         clrng.set_cycle_msecs((int)(1000/cps));
         if(col_map && errmsg)
            snprintf(errmsg, MSG_SZ,
                  "cycle_rate '%s' is also a valid colour map name", parts[i]);
      }
      else if(strspn(parts[i], "vef") == strlen(parts[i])) {
         conv_elems |= 4*(strchr(parts[i], 'f')!=0) +
                       2*(strchr(parts[i], 'e')!=0) +
                       1*(strchr(parts[i], 'v')!=0);
         if(col_map && errmsg)
            snprintf(errmsg, MSG_SZ,
                  "conversion elements '%s' is also a valid colour map name",
                  parts[i]);
      }
      else if(col_map) {
         clrng.add_cmap(col_map);
         if(errmsg && *errmsg2)
            snprintf(errmsg, MSG_SZ, "map %d: %s", i+1, errmsg2);
      }
      else {
         if(errmsg)
            strcpy(errmsg, errmsg2);
         return 0;
      }
   }

   if(!conv_elems)
      conv_elems = 7;  // no elements specified, use "fev"
   for(int i=0; i<3; i++) {
      if((conv_elems & (1<<i)))
         clrngs[i] = clrng;
   }
   
   return 1;
}


/*
RGB->RYB
0 >= 60 degrees, multiply by 2
60 >= 120 degrees, add 60
120 >= 180 maps to 180 to 210
180 >= 240 maps to 210 to 240

So, compression happens
at 120 it is adding 60
at 180 it is adding 30
at 240 it is adding 0

Then
120 >= 240 degrees, N+((240-N)/2)


RYB->RGB
0 >= 120 degrees, divide by 2
120 >= 180 degrees, subtract 60
180 maps to 120
210 maps to 180
240 maps to itself

So, decompression happens
at 180 it is subtracting 60
at 210 it is subtracting 30
at 240 it is subtracting 0

Then
180 >= 240 degrees, N-(240-N)
*/

// angle represented by 0 to 360 degrees
// input: HSV/HSL angle
// output: angle adjusted for RYB mode
double hsx_to_ryb(double angle)
{
   if (angle > 0.0 && angle <= 60.0)
      angle *= 2.0;
   else
   if (angle > 60.0 && angle <= 120.0)
      angle += 60.0;
   else
   if (angle > 120.0 && angle <= 240.0)
      angle += (240.0-angle)/2;
   return angle;
}

// angle represented by 0 to 360 degrees
// input: angle adjusted for RYB mode
// output: HSV/HSL angle
double ryb_to_hsx(double angle)
{
   if (angle > 0.0 && angle <= 120.0)
      angle /= 2.0;
   else
   if (angle > 120.0 && angle <= 180.0)
      angle -= 60.0;
   else
   if (angle > 180.0 && angle <= 240.0)
      angle -= (240.0-angle);
   return angle;
}

col_val rgb_complement(const col_val &col, bool ryb_mode)
{
   if (!col.is_val() || col.is_inv())
      return col;

   // only need hue so algorithm doesn't matter
   vec4d hsxa = col.get_hsva();
   double angle = rad2deg(2*M_PI*hsxa[0]);
 
   if (ryb_mode)
      angle = hsx_to_ryb(angle);

   angle += 180.0;
   if (angle >= 360.0)
      angle -= 360.0;
      
   if (ryb_mode)
      angle = ryb_to_hsx(angle);
   angle /= 360.0;
   
   col_val rcol;
   rcol.set_hsva(angle,hsxa[1],hsxa[2],hsxa[3]);
   return rcol;
}

// wrapper for multiple HSV/HSL algorithms
vec4d get_hsxa(const col_val &col, int color_system_mode)
{
   return (color_system_mode == 1) ? col.get_hsva() : col.get_hsla();
}

col_val set_hsxa(double hue, double sat, double val, double alpha, int color_system_mode)
{
   col_val col;
   if (color_system_mode <= 1)
      col.set_hsva(hue,sat,val,alpha);
   else
   if (color_system_mode == 2)
      col.set_hsla(hue,sat,val,alpha);
   return col;
}

// core code furnished by Adrian Rossiter
col_val blend_HSX_centroid(const vector<col_val> &cols, int color_system_mode, double sat_power, double sat_threshold, double value_power, double value_advance, int alpha_mode, bool ryb_mode)
{
   // no colors, return unset color. one color, return that color
   if (!cols.size())
      return col_val();
   else
   if (cols.size() == 1)
      return cols[0];
      
   // saturation power can't be 0 or less
   if (sat_power <= 0.0)
      sat_power = 1.0;
   
   // can't blend two or less colors to black
   if (cols.size() < 3)
      value_power = 0.0;
   
   double saturation_sum = 0.0;

   double alpha_min = 1.1;
   double alpha_max = -0.1;

   // check for unset, map index, or invisible
   col_val map_found;
   bool invisible_found = false;
   bool unset_found = false;

   int cols_sz = cols.size();
 
   vec4d sum(0.0, 0.0, 0.0, 0.0);
   for(unsigned int i=0; i<cols.size(); i++) {
      // indexes, invisible or unset are not averaged in
      if (cols[i].is_idx()) {
         if (!map_found.is_set())
            map_found = cols[i];
         cols_sz--;
         continue;
      }
      else
      if (cols[i].is_inv()) {
         invisible_found = true;
         cols_sz--;
         continue;
      }
      else
      if (!cols[i].is_set()) {
         unset_found = true;
         cols_sz--;
         continue;
      }

      vec4d hsxa = get_hsxa(cols[i], color_system_mode);
         
      double S = pow(hsxa[1], sat_power); // map onto distorted disc
      double angle = 2*M_PI*hsxa[0];
      
      // RYB mode
      if (ryb_mode)
         angle = deg2rad(hsx_to_ryb(rad2deg(angle)));

      // if value_power is set, simulate subtractive coloring for 3 or more colors
      double V = (value_power <= 0.0) ? hsxa[2] : pow(fabs(60.0-fmod(rad2deg(angle)+(ryb_mode ? 60.0 : 0.0)+value_advance,120.0))/60, value_power);

      alpha_min = (hsxa[3] < alpha_min) ? hsxa[3] : alpha_min;
      alpha_max = (hsxa[3] > alpha_max) ? hsxa[3] : alpha_max;

      sum += vec4d(S*cos(angle), S*sin(angle), V, hsxa[3]); // point in cylinder

      if (sat_threshold < 1.0)
         saturation_sum+=hsxa[1];
   }

   // if no colors are being averaged, all are a mix of unset, indexes, and/or invisible
   // heirarchy: first map index, invisible, unset
   if (!cols_sz) {
      if (map_found.is_set())
         return(map_found);
      else
      if (invisible_found)
         return(col_val(col_val::invisible));
      else
      if (unset_found)
         return(col_val());
   }

   // average
   sum /= cols_sz;

   double H = 0.0;
   double S = 0.0;
   
   // saturation
   S = pow(sum[0]*sum[0]+sum[1]*sum[1], 0.5/sat_power); // map back from distorted disc
   
   // saturations less than 1/255 will happen due to inaccuracy of HSx->RGB conversion
   // note that 255,255,254 has a saturation of 1/255 = 0.00392157..., which is the smallest valid saturation
   if (S<1/255.0)
      S = 0.0;

   // saturation of color centroid is higher than sat_threshold, use average saturation         
   if (S > sat_threshold)
      S = saturation_sum/cols_sz;
   
   // hue
   // if saturation is 0, no need to calculate hue
   //if (double_ne(sum[0],0.0,epsilon) && double_ne(sum[1],0.0,epsilon)) { // old error, made nice pattern, but wrong
   if (S != 0.0) {
      H = atan2(sum[1], sum[0])/(2*M_PI);
      if(H<0)
         H += 1.0;

      // RYB mode
      if (ryb_mode)
         H = deg2rad(ryb_to_hsx(rad2deg(H*2*M_PI)))/(2*M_PI);
   }

   double A = (alpha_mode <= 1) ? sum[3] : ((alpha_mode == 2) ? alpha_min : alpha_max);
      
   col_val col = set_hsxa(H, S, sum[2], A, color_system_mode);
 
   return col;
}

col_val blend_RGB_centroid(const vector<col_val> &cols, int alpha_mode, bool ryb_mode)
{
   // no colors, return unset color. one color, return that color
   if (!cols.size())
      return col_val();
   else
   if (cols.size() == 1)
      return cols[0];

   double alpha_min = 1.1;
   double alpha_max = -0.1;

   // check for unset, map index, or invisible
   col_val map_found;
   bool invisible_found = false;
   bool unset_found = false;

   int cols_sz = cols.size();

   vec4d col(0.0, 0.0, 0.0, 0.0);
   for(unsigned int i=0; i<cols.size(); i++) {
     // indexes, invisible or unset are not averaged in
      if (cols[i].is_idx()) {
         if (!map_found.is_set())
            map_found = cols[i];
         cols_sz--;
         continue;
      }
      else
      if (cols[i].is_inv()) {
         invisible_found = true;
         cols_sz--;
         continue;
      }
      else
      if (!cols[i].is_set()) {
         unset_found = true;
         cols_sz--;
         continue;
      }

      col_val rcol = cols[i];
      if (ryb_mode) {
         // only need hue so algorithm doesn't matter
         vec4d hsxa = rcol.get_hsva();
         hsxa[0] = hsx_to_ryb(rad2deg(2*M_PI*hsxa[0]))/360.0;
         rcol.set_hsva(hsxa);
      }

      double A = 1.0 - rcol.get_transd();
      alpha_min = (A < alpha_min) ? A : alpha_min;
      alpha_max = (A > alpha_max) ? A : alpha_max;

      col += rcol.get_vec4d();
   }

   // if no colors are being averaged, all are a mix of unset, indexes, and/or invisible
   // heirarchy: first map index, invisible, unset
   if (!cols_sz) {
      if (map_found.is_set())
         return(map_found);
      else
      if (invisible_found)
         return(col_val(col_val::invisible));
      else
      if (unset_found)
         return(col_val());
   }

   col /= cols_sz;

   if (ryb_mode) {
      col_val rcol(col);
      // only need hue so algorithm doesn't matter
      vec4d hsxa = rcol.get_hsva();
      hsxa[0] = ryb_to_hsx(rad2deg(2*M_PI*hsxa[0]))/360.0;
      rcol.set_hsva(hsxa);
      col = rcol.get_vec4d();
   }

   if (alpha_mode > 1)
      col[3] = (alpha_mode == 2) ? alpha_min : alpha_max;

   return col;
}

