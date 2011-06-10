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
   Name: lattice_grid.cc
   Description: grids and lattices with integer coordinates
   Project: Antiprism - http://www.antiprism.com
*/

#include <stdio.h>
#include <stdlib.h>

#include <ctype.h>
#include <math.h>

#include <string>
#include <vector>
#include <algorithm>

#include "../base/antiprism.h"
#include "lattice_grid.h"

bool sc_test(int /*x*/, int /*y*/, int /*z*/)        // dist2 = 1, 2, 3
{ return 1; }

bool fcc_test(int x, int y, int z)                   // dist2 = 2, 4, 6, 12
{ return (x+y+z)%2==0; }

bool bcc_test(int x, int y, int z)                   //dist2 = 3, 4, 8
{   return ( (x%2&&y%2&&z%2) || !(x%2||y%2||z%2) ); }

bool hcp_test(int x, int y, int z)                   // dist2 = 18
{ 
   int m = x+y+z;
   int n = (m < 0); // for integer division  
   return (m%6==0) && ((x-m/12+n)%3==0) && ((y-m/12+n)%3==0);
}

bool rh_dodec_test(int x, int y, int z)              // dist2 = 3 (8)
{   return ( (x%2&&y%2&&z%2) || !(((x+y+z)/2)%2==0||(x%2||y%2||z%2)) ); }

bool cubo_oct_test(int x, int y, int z)              // dist2 = 2
{ return (abs(x)%2 + abs(y)%2 + abs(z)%2)==2; }

bool tr_oct_test(int x, int y, int z)                // dist2 = 2
{ return ( (z%4==0 && x%4 && y%4 && (x-y)%2) ||
           (y%4==0 && z%4 && x%4 && (z-x)%2) ||
           (x%4==0 && y%4 && z%4 && (y-z)%2) ); }

bool tr_tet_tet_test(int x, int y, int z)            // dist2 = 2
{ return ( ( x%2==0 && (y%4+4)%4 == (((x + z))%4+4)%4) ||
           ( y%2==0 && (z%4+4)%4 == (((y + x))%4+4)%4) ||
           ( z%2==0 && (x%4+4)%4 == (((z + y))%4+4)%4) ); }

bool tr_tet_tr_oct_cubo_test(int x, int y, int z)    // dist2 = 4
{ x=abs(x)%6; if(x>3) x=6-x;
   y=abs(y)%6; if(y>3) y=6-y;
   z=abs(z)%6; if(z>3) z=6-z;
   int dist2 = x*x + y*y;
   return ( (z%6==0 && (dist2==2 || dist2==8 )) ||
            (z%6==1 && (dist2==1 || dist2==13)) ||
            (z%6==2 && (dist2==4 || dist2==10)) ||
            (z%6==3 && dist2==5) ); }

bool diamond_test(int x, int y, int z)               //  dist2 = 3
{ return ( ((x%2+2)%2+(y%2+2)%2+(z%2+2)%2)%3==0 &&
            (x/2+(x%2<0)+y/2+(y%2<0)+z/2+(z%2<0))%2==0); }

// Coordinates from Wendy Krieger
bool hcp_diamond_test(int x, int y, int z)            //  dist2 = 27
{
   int pt[][3] = {{0,0,0}, {3,3,3}, {6,0,6}, {9,3,9}};
   for(int i=0; i<4; i++) {
      int tx = x-pt[i][0];
      int ty = y-pt[i][1];
      int tz = z-pt[i][2];
      int sum = tx+ty+tz;
      if(sum%24==0) {
         int n8 = sum/3;
         if( (tx-n8)%6==0 && (ty-n8)%6==0 && (tz-n8)%6==0)
            return true;
      }
   }     
   return false;
}        
 
// Coordinates from Vladimir Bulatov
bool k_4_test(int x, int y, int z)                   //  dist2 = 2
{
   if((x+y+z)%2)
      return false;
   x = (x%4) + ((x<0)?4:0);
   y = (y%4) + ((y<0)?4:0);
   z = (z%4) + ((z<0)?4:0);
   if( (x==0&&y==0&&z==0) || 
       (x==0&&y==1&&z==3) || 
       (x==1&&y==0&&z==1) || 
       (x==1&&y==1&&z==2) || 
       (x==2&&y==2&&z==2) || 
       (x==2&&y==3&&z==1) || 
       (x==3&&y==2&&z==3) || 
       (x==3&&y==3&&z==0))
      return true;
   return false;
}

void add_struts(geom_if &geom, int len2)
{
   const vector<vec3d> &verts = geom.verts();
   for(unsigned int i=0; i<verts.size(); i++)
      for(unsigned int j=i; j<verts.size(); j++) {
         if(fabs((verts[i]-verts[j]).mag2() - len2) < epsilon)
            geom.add_edge(make_edge(i, j));
      }
}

void int_lat_grid::make_lattice(geom_if &geom)
{  
   if(!centre.is_set())
      centre = vec3d(1,1,1)*(o_width/2.0);
   double o_off = o_width/2.0 + epsilon;
   double i_off = i_width/2.0 - epsilon;
   int i, j, k;
   for (k=int(ceil(centre[2]-o_off)); k<=centre[2]+o_off; k++)
      for (j=int(ceil(centre[1]-o_off)); j<=centre[1]+o_off; j++)
         for (i=int(ceil(centre[0]-o_off)); i<=centre[0]+o_off; i++) {
            if(i>centre[0]-i_off && i<centre[0]+i_off &&
               j>centre[1]-i_off && j<centre[1]+i_off &&
               k>centre[2]-i_off && k<centre[2]+i_off)
               continue;
            if(coord_test(i, j, k))
               geom.add_vert(vec3d(i, j, k));
         }
}

void sph_lat_grid::make_lattice(geom_if &geom)
{
   if(!centre.is_set())
      centre = vec3d(0,0,0);
   double o_off = o_width + epsilon;
   double i_off = i_width - epsilon;
   int i, j, k;
   for (k=int(ceil(centre[2]-o_off)); k<=centre[2]+o_off; k++)
      for (j=int(ceil(centre[1]-o_off)); j<=centre[1]+o_off; j++)
         for (i=int(ceil(centre[0]-o_off)); i<=centre[0]+o_off; i++) {
            double dist2 = (vec3d(i, j, k)-centre).mag2();
            if(o_off<dist2 || i_off>dist2)
               continue;
            if(coord_test(i, j, k))
               geom.add_vert(vec3d(i, j, k));
         }
}

// for lattice code only

double lattice_radius(const geom_if &geom, const char &radius_type)
{
   geom_v tgeom = geom;
   
   char errmsg[MSG_SZ]="";
   // allow radius to be calculated for some polygons
   int dimensions = tgeom.set_hull("",errmsg);
   if(dimensions < 0) {
      fprintf(stderr,"%s\n",errmsg);
      fprintf(stderr,"lattice_radius: warning: convex hull could not be created\n");
      return 0;
   }
   tgeom.orient();
   
   geom_info rep(tgeom);
   rep.set_center(centroid(tgeom.verts()));
   
   double radius = 0;
   if (radius_type == 'k')
      radius = rep.vert_dists().max;
   else
   if (radius_type == 's' && dimensions > 1)
      radius = rep.face_dists().min;
   else
   if (radius_type == 'l' && dimensions > 1)
      radius = rep.face_dists().max;
      
   return radius;
}

void geom_container_clip(col_geom_v &geom, col_geom_v &container, const double &radius, const vec3d &offset, double eps)
{
   // container has to be convex and 3 dimensional
   char errmsg[MSG_SZ]=""; 
   if(container.set_hull("",errmsg) < 0) {
      fprintf(stderr,"%s\n",errmsg);
      fprintf(stderr,"geom_container_clip: warning: convex hull could not be created\n");
      return;
   }
   container.orient();

   // standardize radius of 1 on maximum vertex. Then set radius
   mat3d trans_m = mat3d::scale((1.0/lattice_radius(container,'k'))*radius);
   container.transform(trans_m);

   const vector<vec3d> &verts = geom.verts();
   vec3d grid_cent = centroid(verts);
   if (offset.is_set())
      grid_cent += offset;

   fprintf(stderr,"info: radius = %g (square root of %g)\n",radius,radius*radius);

   // translate container to center of grid
   vec3d container_cent = centroid(container.verts());
   trans_m = mat3d::transl(-container_cent+grid_cent);
   container.transform(trans_m);

   vector<int> del_verts;
   for(unsigned int i=0; i<verts.size(); i++) {
      if (!is_point_inside_hull(verts[i], container, eps))
         del_verts.push_back(i);
   }

   if (del_verts.size())
      geom.delete_verts(del_verts);

   if (!verts.size())
      fprintf(stderr,"bravais_container_clip: warning: all vertices were clipped out!\n");
}

void geom_spherical_clip(col_geom_v &geom, const double &radius, const vec3d &offset, double eps)
{
   const vector<vec3d> &verts = geom.verts();
   vec3d cent = centroid(verts);
   if (offset.is_set())
      cent += offset;

   fprintf(stderr,"info: radius = %g (square root of %g)\n",radius,radius*radius);

   vector<int> del_verts;
   for(unsigned int i=0; i<verts.size(); i++) {
      double len = (cent-verts[i]).mag();
      if (double_gt(len, radius, eps))
         del_verts.push_back(i);
   }

   if (del_verts.size())
      geom.delete_verts(del_verts);

   if (!verts.size())
      fprintf(stderr,"bravais_spherical_clip: warning: all vertices were clipped out!\n");
}

void list_grid_radii(const col_geom_v &geom, const vec3d &offset, double eps)
{
   const vector<vec3d> &verts = geom.verts();
   vec3d cent = centroid(verts);
   if (offset.is_set())
      cent += offset;

   vector<double> radii;
   for(unsigned int i=0;i<verts.size();i++)
      radii.push_back((cent-verts[i]).mag());

   sort( radii.begin(), radii.end() );

   // eliminate 0 radius
   int start = 0;
   if (double_eq(radii[0], 0, eps))
      start++;

   // prime things
   double comp = radii[start];
   int occur_total = 0;
   int occur = 1;
   int rank = 1;

   fprintf(stderr,"\nList of unique radial distances from center (and offset) in grid\n\n");

   fprintf(stderr,"Rank\tDistance\tD Squared\tOccurrence\n");
   fprintf(stderr,"----\t--------\t---------\t----------\n");
   for(unsigned int i=start+1;i<radii.size();i++) {
      if (double_eq(radii[i], comp, eps))
         occur++;
      else {
         occur_total += occur;
         fprintf(stderr,"%d\t%-8g\t%-8g\t%d\n",rank,comp,comp*comp,occur);
         comp = radii[i];
         occur = 1;
         rank++;
      }
   }
   occur_total += occur;
   fprintf(stderr,"%d\t%-8g\t%-8g\t%d\n\n",rank,comp,comp*comp,occur);
   fprintf(stderr,"Total occurrences = %d\n\n",occur_total);
}

void list_grid_struts(const col_geom_v &geom, double eps)
{
   const vector<vec3d> &verts = geom.verts();

   vector<double> struts;
   for(unsigned int i=0;i<verts.size();i++)
      for(unsigned int j=i+1;j<verts.size();j++)
         struts.push_back((verts[i]-verts[j]).mag());

   sort( struts.begin(), struts.end() );

   // eliminate 0 radius
   int start = 0;
   if (double_eq(struts[0], 0, eps))
      start++;

   // prime things
   double comp = struts[start];
   int occur_total = 0;
   int occur = 1;
   int rank = 1;

   fprintf(stderr,"\nList of unique strut lengths in grid\n\n");

   fprintf(stderr,"Rank\tDistance\tD Squared\tOccurrence\n");
   fprintf(stderr,"----\t--------\t---------\t----------\n");
   for(unsigned int i=start+1;i<struts.size();i++) {
      if (double_eq(struts[i], comp, eps))
         occur++;
      else {
         occur_total += occur;
         fprintf(stderr,"%d\t%-8g\t%-8g\t%d\n",rank,comp,comp*comp,occur);
         comp = struts[i];
         occur = 1;
         rank++;
      }
   }
   occur_total += occur;
   fprintf(stderr,"%d\t%-8g\t%-8g\t%d\n\n",rank,comp,comp*comp,occur);
   fprintf(stderr,"Total occurrences = %d\n\n",occur_total);
}

void add_color_struts(col_geom_v &geom, const double &len2, col_val &edge_col, double eps)
{
   const vector<vec3d> &verts = geom.verts();

   for(unsigned int i=0; i<verts.size(); i++)
      for(unsigned int j=i; j<verts.size(); j++) {
         if(fabs((verts[i]-verts[j]).mag2() - len2) < eps)
            geom.add_col_edge(make_edge(i, j), edge_col);
      }
}

void color_centroid(col_geom_v &geom, col_val &cent_col, double eps)
{
   const vector<vec3d> &verts = geom.verts();
   vec3d cent = centroid(verts);
   int cent_idx = find_vertex_by_coordinate(geom, cent, eps);
   
   if (cent_idx == -1)
      geom.add_col_vert(cent, cent_col);
   else 
      geom.set_v_col(cent_idx, cent_col);
}

// color functions
         
// Rotational octahedral by Adrian Rossiter
vec3d sort_vec3d_chiral(const vec3d &v, double eps)
{
   vec3d c = v;
   // Rotate into positive octant
   if(c[0]<0) {
      c[0] = -c[0];
      c[2] = -c[2];
   }
   if(c[1]<0) {
      c[1] = -c[1];
      c[2] = -c[2];
   }
   if(c[2]<0) {
      std::swap(c[0], c[1]);
      c[2] = -c[2];
   }

   // if c[1] is maximum rotate to first place: 1,2,0
   if (c[1] > c[0] && c[1] > c[2]-eps)
      c = vec3d(c[1],c[2],c[0]);
   else
   // if c[2] is maximum rotate to first place: 2,0,1
   if (c[2] > c[0] && c[2] > c[1]+eps)
      c = vec3d(c[2],c[0],c[1]);
   else
   // if c[0] is maximum do nothing
      c = vec3d(c[0],c[1],c[2]);

   // Check whether c is near negative triangle external boundary, and
   // rotate to corresponding positive triangle boundary if so.
   if(double_eq(c[0], c[1], eps))
      c = vec3d(c[1], c[0], c[2]);
   if(double_eq(c[2], 0, eps))
      c = vec3d(c[0], -c[2], c[1]);

   return c;
}

// sort an absolute value of vec3d without altering original vec3d
vec3d sort_vec3d(vec3d &v)
{
   vector<double> c;
   c.push_back(fabs(v[0]));
   c.push_back(fabs(v[1]));
   c.push_back(fabs(v[2]));

   sort( c.begin(), c.end() );

   return(vec3d(c[0],c[1],c[2]));
}

void color_by_symmetry_normals(col_geom_v &geom, const char &color_method, const int &face_opacity, double eps)
{
   const vector<vector<int> > &faces = geom.faces();
   const vector<vec3d> &verts = geom.verts();

   color_map_range_rand_hsv cmap;
   //f_coloring clrg(&geom);

   // transparency
   if (face_opacity != -1)
      cmap.set_range(3, vector<double>(1, (double)face_opacity/255));

   for(unsigned int i=0; i<faces.size(); i++) {
      vec3d norm = face_norm(verts, faces[i]).unit();
      if(color_method == 's' || color_method == 'S')
         norm = sort_vec3d(norm);
      else if(color_method == 'c' || color_method == 'C')
         norm = sort_vec3d_chiral(norm, eps);

      long idx = (long)(norm[0]*1000000) + (long)(norm[1]*10000) + (long)norm[2]*100;
      if(color_method == 'S' || color_method == 'C')
         geom.set_f_col(i,cmap.get_col(idx));
      else
         geom.set_f_col(i,idx);
   }
}

void color_edges_by_sqrt(col_geom_v &geom, const char &color_method)
{
   geom.add_missing_impl_edges();

   color_map_range_rand_hsv cmap;
   //e_coloring clrg(&geom);
   for(unsigned int i=0; i<geom.edges().size(); i++) {
      // geom.set_e_col(i, int(floor(pow(geom.edge_len(i),2)+0.5)));
      int idx = int(floor(pow(geom.edge_len(i),2)+0.5));
      if (color_method == 'R')
         geom.set_e_col(i,cmap.get_col(idx));
         //geom.set_e_col(i,clrg.idx_to_rand_val(idx));
      else
         geom.set_e_col(i,idx);
   }
}

// convex hull and voronoi wrappers

void convex_hull_report(const geom_v &geom, const bool &add_hull)
{
   geom_info rep(geom);
   fprintf(stderr,"\n");
   fprintf(stderr, "convex hull information:\n");
   if (!add_hull)
      fprintf(stderr, "num_verts = %d\n", rep.num_verts());
   fprintf(stderr, "num_faces = %d\n", rep.num_faces());
   if (!add_hull && rep.num_verts() > 1)
      fprintf(stderr, "num_edges = %d\n", rep.num_verts()+rep.num_faces()-2);
   if (rep.num_verts() > 2) {
      double area = rep.face_areas().sum;
      fprintf(stderr, "area      = %.17g\n", area);
      fprintf(stderr, "volume    = %.17g\n", fabs(rep.volume()));
      if (area) {
         fprintf(stderr, "isoperimetric quotient (spherical nature: v^2/a^3 x 36 x PI = 1 is sphere)\n"); 
         fprintf(stderr, "  iq      = %.17g\n", rep.isoperimetric_quotient());
      }
   }
   fprintf(stderr, "end convex hull information\n");
   fprintf(stderr,"\n");
}

int get_voronoi_geom(col_geom_v &geom, col_geom_v &vgeom, const bool &central_cells, const bool &one_cell_only, double eps)
{
   // do this in case compound lattice was sent. Simultaneous points cause problems for Voronoi Cells
   sort_merge_elems(geom, "vef", eps);
   
   // store convex hull of lattice as speed up to is_geom_inside_hull()
   // has to be 3 dimensional
   col_geom_v hgeom = geom;
   char errmsg[MSG_SZ]=""; 
   if(hgeom.set_hull("",errmsg) < 0) {
      fprintf(stderr,"%s\n",errmsg);
      fprintf(stderr,"get_voronoi_geom: warning: convex hull could not be created\n");
      return 0;
   }
   hgeom.orient();
   
   vec3d cent = centroid(hgeom.verts());

   vector<col_geom_v> cells;
   get_voronoi_cells(geom, cells);
   
   for(unsigned int i=0; i<cells.size(); i++) {      
      if (central_cells && !is_point_inside_hull(cent, cells[i], eps)) {
         continue;
      }
      else
      if (!is_geom_inside_hull(cells[i], hgeom, eps)) {
         continue;
      }
      vgeom.append(cells[i]);
      if (one_cell_only)
         break;
   }
 
   if (!(vgeom.verts()).size()) {
      fprintf(stderr,"get_voronoi_geom: warning: after Voronoi cells, geom is empty\n");
      return 0;
   }
   
   return 1;
}
