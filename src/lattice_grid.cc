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


#include <math.h>
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
{ int m=x+y+z; return (m%6==0) && ((x-m/12)%3==0) && ((y-m/12)%3==0); }

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
   vector<vec3d> &verts = *geom.get_verts();
   for(unsigned int i=0; i<verts.size(); i++)
      for(unsigned int j=i; j<verts.size(); j++) {
         if(fabs((verts[i]-verts[j]).mag2() - len2) < epsilon)
            geom.add_edge(make_edge(i, j));
      }
}

void int_lat_grid::make_lattice(geom_if &geom)
{
   if(!centre.is_set())
      centre = vec3d(1,1,1)*(o_width/2);
   double o_off = o_width/2 + epsilon;
   double i_off = i_width/2 - epsilon;
   int i, j, k;
   for (k=int(ceil(centre[2]-o_off)); k<=centre[2]+o_off; k++)
      for (j=int(ceil(centre[1]-o_off)); j<=centre[1]+o_off; j++)
         for (i=int(ceil(centre[0]-o_off)); i<=centre[0]+o_off; i++) {
            if(i>centre[0]-i_off && i<centre[0]+i_off &&
               j>centre[1]-i_off && j<centre[1]+i_off &&
               k>centre[2]-i_off && k<centre[2]+i_off)
               continue;
            if(coord_test(i, j, k))
               geom.get_verts()->push_back(vec3d(i, j, k));
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
               geom.get_verts()->push_back(vec3d(i, j, k));
         }
}


// new stuff starts here

// some actual stuff for lattice code only

void geom_container_clip(col_geom_v &geom, col_geom_v &container, double radius, vec3d offset, double epsilon)
{
   if(!container.set_hull()) {
      fprintf(stderr,"warning from geom_container_clip: Convex hull could not be created\n");
      return;
   }
   container.orient();

   // standardize radius of 1 on maximum vertex. Then set radius
   geom_info rep(container);
   rep.set_center(centroid(container.verts()));
   mat3d trans_m = mat3d::scale((1.0/rep.vert_dists().max)*radius);
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
      if (!is_point_inside_hull(verts[i], container, epsilon))
         del_verts.push_back(i);
   }

   if (del_verts.size())
      geom.delete_verts(del_verts);

   if (!verts.size())
      fprintf(stderr,"warning from bravais_container_clip: all vertices were clipped out!\n");
}

void geom_spherical_clip(col_geom_v &geom, double radius, vec3d offset, double epsilon)
{
   const vector<vec3d> &verts = geom.verts();
   vec3d cent = centroid(verts);
   if (offset.is_set())
      cent += offset;

   fprintf(stderr,"info: radius = %g (square root of %g)\n",radius,radius*radius);

   vector<int> del_verts;
   for(unsigned int i=0; i<verts.size(); i++) {
      double len = (cent-verts[i]).mag();
      if (!double_equality(len, radius, epsilon) && len > radius)
         del_verts.push_back(i);
   }

   if (del_verts.size())
      geom.delete_verts(del_verts);

   if (!verts.size())
      fprintf(stderr,"warning from bravais_spherical_clip: all vertices were clipped out!\n");
}

void list_grid_radii(col_geom_v &geom, vec3d offset, double epsilon)
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
   if (double_equality(radii[0], 0, epsilon))
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
      if (double_equality(radii[i], comp, epsilon))
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

void list_grid_struts(col_geom_v &geom, double epsilon)
{
   const vector<vec3d> &verts = geom.verts();

   vector<double> struts;
   for(unsigned int i=0;i<verts.size();i++)
      for(unsigned int j=i+1;j<verts.size();j++)
         struts.push_back((verts[i]-verts[j]).mag());

   sort( struts.begin(), struts.end() );

   // eliminate 0 radius
   int start = 0;
   if (double_equality(struts[0], 0, epsilon))
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
      if (double_equality(struts[i], comp, epsilon))
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

void add_color_struts(col_geom_v &geom, double len2, col_val edge_col)
{
   const vector<vec3d> &verts = geom.verts();

   for(unsigned int i=0; i<verts.size(); i++)
      for(unsigned int j=i; j<verts.size(); j++) {
         if(fabs((verts[i]-verts[j]).mag2() - len2) < epsilon) {
            vector<int> edge(2);
            edge[0] = i;
            edge[1] = j;
            geom.add_col_edge(edge, edge_col);
         }
      }
}

void color_centroid(col_geom_v &geom, col_val cent_col)
{
   const vector<vec3d> &verts = geom.verts();
   vec3d cent = centroid(verts);
   int cent_idx = -1;
   for(unsigned int i=0; i<verts.size(); i++) {
      if (double_equality(verts[i][0], cent[0], epsilon) &&
          double_equality(verts[i][1], cent[1], epsilon) &&
          double_equality(verts[i][2], cent[2], epsilon)) {
         cent_idx = i;
         break;
      }
   }
   
   if (cent_idx == -1)
      geom.add_col_vert(cent, cent_col);
   else 
      geom.set_v_col(cent_idx, cent_col);
}

// Rotational octahedral by Adrian Rossiter
vec3d sort_vec3d_chiral(const vec3d &v)
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
   if (c[1] > c[0] && c[1] > c[2]-epsilon)
      c = vec3d(c[1],c[2],c[0]);
   else
   // if c[2] is maximum rotate to first place: 2,0,1
   if (c[2] > c[0] && c[2] > c[1]+epsilon)
      c = vec3d(c[2],c[0],c[1]);
   else
   // if c[0] is maximum do nothing
      c = vec3d(c[0],c[1],c[2]);

   // Check whether c is near negative triangle external boundary, and
   // rotate to corresponding positive triangle boundary if so.
   if(double_equality(c[0], c[1], epsilon))
      c = vec3d(c[1], c[0], c[2]);
   if(double_equality(c[2], 0, epsilon))
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

void color_by_symmetry_normals(col_geom_v &geom, char color_method, int face_opacity)
{
   const vector<vector<int> > &faces = geom.faces();
   const vector<vec3d> &verts = geom.verts();

   color_map_range_rand_hsv cmap;
   //f_coloring clrg(&geom);

   // transparency
   if (face_opacity != -1)
      cmap.set_range(3, vector<double>(1, (double)face_opacity/255));
      //clrg.set_alp_range((double)face_opacity/255);

   for(unsigned int i=0; i<faces.size(); i++) {
      vec3d norm = face_norm(verts, faces[i]).unit();
      if(color_method == 's' || color_method == 'S')
         norm = sort_vec3d(norm);
      else if(color_method == 'c' || color_method == 'C')
         norm = sort_vec3d_chiral(norm);

      long idx = (long)(norm[0]*1000000) + (long)(norm[1]*10000) + (long)norm[2]*100;
      if(color_method == 'S' || color_method == 'C')
         geom.set_f_col(i,cmap.get_col(idx));
         //geom.set_f_col(i,clrg.idx_to_rand_val(idx));
      else
         geom.set_f_col(i,idx);
   }
}

void color_edges_by_sqrt(col_geom_v &geom, char color_method)
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

int do_convex_hull(col_geom_v &geom, bool add_hull, bool verbose)
{
   const vector<vec3d> &verts = geom.verts();

   if (verts.size() < 4) {
      fprintf(stderr,"warning from do_convex_hull: Convex hull not possible with only %lu points\n", (unsigned long)verts.size());
      return 0;
   }

   if(!(add_hull ? geom.add_hull() : geom.set_hull())) {
      fprintf(stderr,"warning from do_convex_hull: Convex hull could not be created\n");
      return 0;
   }
   
   geom.orient();
      
   if (verbose && !add_hull) {
      geom_info rep(geom);
      fprintf(stderr, "convex hull information:\n");
      fprintf(stderr, "num_verts = %d\n", rep.num_verts());
      fprintf(stderr, "num_faces = %d\n", rep.num_faces());
      fprintf(stderr, "num_edges = %d\n", rep.num_verts()+rep.num_faces()-2);
      fprintf(stderr, "area      = %.17g\n", rep.face_areas().sum);
      fprintf(stderr, "volume    = %.17g\n", fabs(rep.volume()));
      fprintf(stderr, "isoperimetric quotient (spherical nature: v^2/a^3 x 36 x PI = 1 is sphere)\n"); 
      fprintf(stderr, "  iq      = %.17g\n", rep.vol2_by_area3()*M_PI*36.0);
      fprintf(stderr, "end convex hull information\n");
   }
   
   return 1;
}

int get_voronoi_geom(col_geom_v &geom, col_geom_v &vgeom, bool central_cells, bool one_cell_only, double epsilon)
{
   // do this in case compound lattice was sent. Simultaneous points cause problems for Voronoi Cells
   sort_merge_elems(geom, "vef", epsilon);
   
   // store convex hull as speed up to is_geom_inside_hull()
   col_geom_v hgeom;
   hgeom = geom;
   hgeom.set_hull();
   hgeom.orient();
   
   vec3d cent = centroid(hgeom.verts());

   vector<col_geom_v> cells;
   get_voronoi_cells(geom, cells);
   
   for(unsigned int i=0; i<cells.size(); i++) {      
      if (central_cells && !is_point_inside_hull(cent, cells[i], epsilon)) {
         continue;
      }
      else
      if (!is_geom_inside_hull(cells[i], hgeom, epsilon)) {
         continue;
      }
      vgeom.append(cells[i]);
      if (one_cell_only)
         break;
   }
 
   if (!(vgeom.verts()).size()) {
      fprintf(stderr,"warning: after Voronoi cells, geom is empty.\n");
      return 0;
   }
   
   return 1;
}


