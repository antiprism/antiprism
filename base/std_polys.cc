/*
   Copyright (c) 2003, Adrian Rossiter

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

/* \file std_polys.cc
 * \brief Platonic and Archimedean polyhedra
 */


#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#include "polygons.h"
#include "std_polys.h"
#include "math_utils.h"
#include "utils.h"
#include "coloring.h"

using std::swap;

void normalised_face_list(geom_if &geom)
{
   geom.orient();
   if(vdot(geom.face_norm(0), geom.face_cent(0))<0)
      geom.orient_reverse();
   sort(geom.raw_faces().begin(), geom.raw_faces().end());
}


void tetrahedron(geom_if &geom)
{
   geom.clear_all();
   geom.add_vert(vec3d(-1,  1, -1)); // 0
   geom.add_vert(vec3d(-1, -1,  1)); // 1
   geom.add_vert(vec3d( 1,  1,  1)); // 2
   geom.add_vert(vec3d( 1, -1, -1)); // 3

   int f[] = { 0,1,2,   0,3,1,   0,2,3,   2,1,3};
   for(int i=0; i<4; i++) {
      vector<int> face;
      for(int j=0; j<3; j++)
         face.push_back(f[i*3+j]);
      geom.add_face(face);
   }
}



void cube(geom_if &geom)
{
   geom.clear_all();
   geom.add_vert(vec3d( 1,  1,  1)); // 0
   geom.add_vert(vec3d( 1,  1, -1)); // 1
   geom.add_vert(vec3d( 1, -1,  1)); // 2
   geom.add_vert(vec3d( 1, -1, -1)); // 3
   geom.add_vert(vec3d(-1,  1,  1)); // 4
   geom.add_vert(vec3d(-1,  1, -1)); // 5
   geom.add_vert(vec3d(-1, -1,  1)); // 6
   geom.add_vert(vec3d(-1, -1, -1)); // 7

   int f[] = { 0,4,6,2,  1,3,7,5,  0,2,3,1,  4,5,7,6,  0,1,5,4,  2,6,7,3 };
   for(int i=0; i<6; i++) {
      vector<int> face;
      for(int j=0; j<4; j++)
         face.push_back(f[i*4+j]);
      geom.add_face(face);
   }
}

void octahedron(geom_if &geom)
{
   geom.clear_all();
   geom.add_vert(vec3d( 1,  0,  0)); // 0
   geom.add_vert(vec3d(-1,  0,  0)); // 1
   geom.add_vert(vec3d( 0,  1,  0)); // 2
   geom.add_vert(vec3d( 0, -1,  0)); // 3
   geom.add_vert(vec3d( 0,  0,  1)); // 4
   geom.add_vert(vec3d( 0,  0, -1)); // 5

   int f[] = { 0,2,4,  0,4,3,  0,3,5,  0,5,2,  1,4,2,  1,3,4,  1,5,3,  1,2,5 };
   for(int i=0; i<8; i++) {
      vector<int> face;
      for(int j=0; j<3; j++)
         face.push_back(f[i*3+j]);
      geom.add_face(face);
   }
}


void dodecahedron(geom_if &geom)
{
   geom.clear_all();
   double phi = (sqrt(5)+1)/2;
   double iphi = 1/phi;
   geom.add_vert(vec3d( 1,  1,  1));     //  0
   geom.add_vert(vec3d( 1,  1, -1));     //  1
   geom.add_vert(vec3d( 1, -1,  1));     //  2
   geom.add_vert(vec3d( 1, -1, -1));     //  3
   geom.add_vert(vec3d(-1,  1,  1));     //  4
   geom.add_vert(vec3d(-1,  1, -1));     //  5
   geom.add_vert(vec3d(-1, -1,  1));     //  6
   geom.add_vert(vec3d(-1, -1, -1));     //  7
   geom.add_vert(vec3d(0,  iphi,  phi)); //  8
   geom.add_vert(vec3d(0,  iphi, -phi)); //  9
   geom.add_vert(vec3d(0, -iphi, -phi)); // 10
   geom.add_vert(vec3d(0, -iphi,  phi)); // 11
   geom.add_vert(vec3d( phi, 0,  iphi)); // 12
   geom.add_vert(vec3d(-phi, 0,  iphi)); // 13
   geom.add_vert(vec3d(-phi, 0, -iphi)); // 14
   geom.add_vert(vec3d( phi, 0, -iphi)); // 15
   geom.add_vert(vec3d( iphi,  phi, 0)); // 16
   geom.add_vert(vec3d( iphi, -phi, 0)); // 17
   geom.add_vert(vec3d(-iphi, -phi, 0)); // 18
   geom.add_vert(vec3d(-iphi,  phi, 0)); // 19


   int f[] = { 12,15, 1,16, 0,  8,11, 2,12, 0, 16,19, 4, 8, 0,
               15,12, 2,17, 3, 10, 9, 1,15, 3, 17,18, 7,10, 3,
               19,16, 1, 9, 5, 14,13, 4,19, 5,  9,10, 7,14, 5,
               13,14, 7,18, 6, 11, 8, 4,13, 6, 18,17, 2,11, 6 };
   for(int i=0; i<12; i++) {
      vector<int> face;
      for(int j=0; j<5; j++)
         face.push_back(f[i*5+j]);
      geom.add_face(face);
   }
}

void icosahedron(geom_if &geom)
{
   geom.clear_all();
   double phi = (sqrt(5)+1)/2;
   geom.add_vert(vec3d(0,  phi,  1)); //  0
   geom.add_vert(vec3d(0, -phi,  1)); //  1
   geom.add_vert(vec3d(0, -phi, -1)); //  2
   geom.add_vert(vec3d(0,  phi, -1)); //  3
   geom.add_vert(vec3d( 1, 0,  phi)); //  4
   geom.add_vert(vec3d( 1, 0, -phi)); //  5
   geom.add_vert(vec3d(-1, 0, -phi)); //  6
   geom.add_vert(vec3d(-1, 0,  phi)); //  7
   geom.add_vert(vec3d( phi,  1, 0)); //  8
   geom.add_vert(vec3d(-phi,  1, 0)); //  9
   geom.add_vert(vec3d(-phi, -1, 0)); // 10
   geom.add_vert(vec3d( phi, -1, 0)); // 11

   int f[] = { 0, 4, 7,   1, 7, 4,   2, 5, 6,   3, 6, 5,
               4, 8,11,   5,11, 8,   6, 9,10,   7,10, 9,
               8, 0, 3,   9, 3, 0,  10, 1, 2,  11, 2, 1,
               0, 8, 4,   0, 7, 9,   1, 4,11,   1,10, 7,
               2,11, 5,   2, 6,10,   3, 5, 8,   3, 9, 6 };
   for(int i=0; i<20; i++) {
      vector<int> face;
      for(int j=0; j<3; j++)
         face.push_back(f[i*3+j]);
      geom.add_face(face);
   }
}

void tr_tetrahedron(geom_if &geom)
{
   tetrahedron(geom);
   truncate_verts(geom, 1/3.0);
   geom.orient();
}
  

void cuboctahedron(geom_if &geom)
{
   cube(geom);
   truncate_verts(geom, 1/2.0);
   geom.orient();
}


void tr_cube(geom_if &geom)
{
   cube(geom);
   truncate_verts(geom, 1-sqrt(1/2.0));
   geom.orient();
}


void tr_octahedron(geom_if &geom)
{
   octahedron(geom);
   truncate_verts(geom, 1/3.0);
   geom.orient();
}


void rhombicuboctahedron(geom_if &geom)
{
   geom.clear_all();
   for(int i=0; i<2; i++)
      for(int j=0; j<2; j++)
         for(int k=0; k<2; k++) {
            vec3d v((1-2*i), (1-2*j), (1-2*k)*(1+sqrt(2)));
            for(int l=0; l<3; l++)
               geom.add_vert(vec3d(v[l], v[(l+1)%3], v[(l+2)%3]));
         }
   geom.add_hull();
   normalised_face_list(geom);
}


void tr_cuboctahedron(geom_if &geom)
{
   geom.clear_all();
   for(int i=0; i<2; i++)
      for(int j=0; j<2; j++)
         for(int k=0; k<2; k++) {
            vec3d v((1-2*i), (1-2*j)*(1+sqrt(2)), (1-2*k)*(1+sqrt(8)));
            for(int l=0; l<3; l++) {
               geom.add_vert(vec3d(v[l], v[(l+1)%3], v[(l+2)%3]));
               geom.add_vert(vec3d(v[(l+1)%3], v[l], v[(l+2)%3]));
            }
         }
   geom.add_hull();
   normalised_face_list(geom);
}


void snub_cube(geom_if &geom)
{
   geom.clear_all();
   double third = 1/3.0;
   double K = third*(pow(17+sqrt(297), third) - pow(-17+sqrt(297), third) -1);
   for(int i=0; i<2; i++)
      for(int j=0; j<2; j++)
         for(int k=0; k<2; k++) {
            vec3d v((1-2*i), (1-2*j)*K, (1-2*k)/K);
            if(!is_even(i+j+k))
               swap(v[0], v[1]);
            for(int l=0; l<3; l++)
               geom.add_vert(vec3d(v[l], v[(l+1)%3], v[(l+2)%3]));
         }
   geom.add_hull();
   normalised_face_list(geom);
}


void icosidodecahedron(geom_if &geom)
{
   dodecahedron(geom);
   truncate_verts(geom, 1/2.0);
   geom.orient();
}


void tr_dodecahedron(geom_if &geom)
{
   double phi = (sqrt(5)+1)/2;
   dodecahedron(geom);
   truncate_verts(geom, 1/(phi+2));
   geom.orient();
}


void tr_icosahedron(geom_if &geom)
{
   icosahedron(geom);
   truncate_verts(geom, 1/3.0);
   geom.orient();
}


void rhombicosidodecahedron(geom_if &geom)
{
   geom.clear_all();
   double phi = (sqrt(5)+1)/2;
   vec3d v[3];
   for(int i=0; i<2; i++)
      for(int j=0; j<2; j++)
         for(int k=0; k<2; k++) {
            v[0] = vec3d ((1-2*i), (1-2*j), (1-2*k)*(phi*phi*phi));
            v[1] = vec3d((1-2*i)*(phi*phi), (1-2*j)*phi, (1-2*k)*(2*phi));
            v[2] = vec3d((1-2*i)*(phi*phi), (1-2*j)*(2+phi), 0);
            for(int l=0; l<3; l++)
               for(int m=0; m<3; m++)
                  if(m!=2 || k)
                     geom.add_vert(vec3d(v[m][l], v[m][(l+1)%3],
                              v[m][(l+2)%3]));
         }
   geom.add_hull();
   normalised_face_list(geom);
}


void tr_icosidodecahedron(geom_if &geom)
{
   geom.clear_all();
   double phi = (sqrt(5)+1)/2;
   vec3d v[5];
   for(int i=0; i<2; i++)
      for(int j=0; j<2; j++)
         for(int k=0; k<2; k++) {
            v[0] = vec3d((1-2*i)/phi, (1-2*j)/phi, (1-2*k)*(3+phi));
            v[1] = vec3d((1-2*i)*2/phi, (1-2*j)*phi, (1-2*k)*(1+2*phi));
            v[2] = vec3d((1-2*i)/phi, (1-2*j)*(phi*phi), (1-2*k)*(3*phi-1));
            v[3] = vec3d((1-2*i)*(2*phi-1), (1-2*j)*2, (1-2*k)*(2+phi));
            v[4] = vec3d((1-2*i)*phi, (1-2*j)*3, (1-2*k)*2*phi);
            for(int l=0; l<5; l++)
               for(int m=0; m<3; m++)
                  geom.add_vert(vec3d(v[l][m], v[l][(m+1)%3], v[l][(m+2)%3]));
         }
   geom.add_hull();
   normalised_face_list(geom);
}


void snub_dodecahedron(geom_if &geom)
{
   geom.clear_all();
   double phi = (sqrt(5)+1)/2;
   double third = 1/3.0;
   double K = pow(phi/2+sqrt(phi-5/27.0)/2, third) +
      pow(phi/2-sqrt(phi-5/27.0)/2, third);
   double A = K - 1/K;
   double B = K*phi + phi*phi + phi/K;
            
   vec3d v[5];
   for(int i=0; i<2; i++)
      for(int j=0; j<2; j++)
         for(int k=0; k<2; k++) {
            if(!is_even(i+j+k))
               continue;
            v[0] = vec3d((1-2*i)*2*A, (1-2*j)*2, (1-2*k)*2*B);
            v[1] = vec3d((1-2*i)*(A+B/phi+phi), (1-2*j)*(-A*phi+B+1/phi),
                     (1-2*k)*(A/phi+B*phi-1));
            v[2] = vec3d((1-2*i)*(-A/phi+B*phi+1), (1-2*j)*(-A+B/phi-phi),
                     (1-2*k)*(A*phi+B-1/phi));
            v[3] = vec3d((1-2*i)*(-A/phi+B*phi-1), (1-2*j)*(A-B/phi-phi),
                     (1-2*k)*(A*phi+B+1/phi));
            v[4] = vec3d((1-2*i)*(A+B/phi-phi), (1-2*j)*(A*phi-B+1/phi),
                     (1-2*k)*(A/phi+B*phi+1));
            for(int l=0; l<5; l++)
               for(int m=0; m<3; m++) {
                  geom.add_vert(vec3d(v[l][m], v[l][(m+1)%3], v[l][(m+2)%3]));
               }
         }
   geom.add_hull();
   normalised_face_list(geom);
}

string expand_abbrevs(const string &name, const char *abbrevs[][2], size_t last)
{
   string expanded;
   char name_cpy[MSG_SZ];
   strncpy(name_cpy, name.c_str(), MSG_SZ);
   vector<char *> parts;
   split_line(name_cpy, parts, "_");
   for(unsigned int i=1; i<parts.size(); i++) {
      size_t j;
      for(j=0; j<last; j++) {
         if(strcmp(parts[i], abbrevs[j][0])==0)
            break;
      }
      if(i>1)
         expanded += ' ';
      expanded += (j<last) ? abbrevs[j][1] : parts[i];
   }
   return expanded;
}


const char *alt_names[][2] = {
   {"tet", "u1" }, 
   {"tetrahedron", "u1" }, 
   {"tr_tetrahedron", "u2"}, 
   {"tr_tet", "u2"}, 
   {"octahedron", "u5" }, 
   {"oct", "u5" }, 
   {"cube", "u6" }, 
   {"cuboctahedron", "u7" }, 
   {"cubo", "u7" }, 
   {"tr_octahedron", "u8" }, 
   {"tr_oct", "u8" }, 
   {"tr_cube", "u9" }, 
   {"rhombicuboctahedron", "u10" }, 
   {"rhombicubo", "u10" }, 
   {"tr_cuboctahedron", "u11" }, 
   {"tr_cubo", "u11" }, 
   {"snub_cube", "u12" }, 
   {"icosahedron", "u22" }, 
   {"ico", "u22" }, 
   {"icosa", "u22" }, 
   {"dodecahedron", "u23" }, 
   {"dod", "u23" }, 
   {"icosidodecahedron", "u24" }, 
   {"icosid", "u24" }, 
   {"tr_icosahedron", "u25" }, 
   {"tr_ico", "u25" }, 
   {"tr_icosa", "u25" }, 
   {"tr_dodecahedron", "u26" }, 
   {"tr_dod", "u26" }, 
   {"rhombicosidodecahedron", "u27" }, 
   {"rhombicosid", "u27" }, 
   {"tr_icosidodecahedron", "u28" }, 
   {"tr_icosid", "u28" }, 
   {"snub_dodecahedron", "u29" },
   {"snub_dod", "u29" }
};

const char *u_abbrevs[][2] = {
   {"tr", "truncated"},
   {"sm", "small"},
   {"gr", "great"},
   {"st", "stellated"},
   {"sn", "snub"},
   {"tet", "tetrahedron"},
   {"ico", "icosahedron"},
   {"icosa", "icosahedron"},
   {"dod", "dodecahedron"},
   {"oct", "octahedron"},
   {"cubo", "cuboctahedron"},
   {"icosid", "icosidodecahedron"}
};



int make_resource_uniform(geom_if &geom, string name, char *errmsg)
{
   if(name.size()<2 || !strchr("uU", name[0]) || name.find('.')!=string::npos)
      return -1; // not uniform name

   uni_poly uni;
   int sym_no;
   if(read_int(name.c_str()+1, &sym_no)) {
      sym_no--;
      if(sym_no<0 || sym_no >= uni.get_last_uniform()) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "uniform polyhedron number out of range");
         return 1; // fail
      }
   }
   else if(name[1] == '_') {
      string expanded = expand_abbrevs(name, u_abbrevs,
            sizeof(u_abbrevs)/sizeof(u_abbrevs[0]));
      sym_no = uni.lookup_sym_no(expanded.c_str());
      if(sym_no == -1) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "invalid uniform polyhedron name");
         return 1; // fail
      }
   }
   else
      return -1; // not uniform name

   
   uni.get_poly(geom, sym_no);

   double e_len = geom.edge_vec(geom.faces(0,0), geom.faces(0,1)).mag();
   geom.transform(mat3d::scale(1/e_len));
   return 0; // name found
}

const char *j_abbrevs[][2] = {
   {"tri", "triangular"},
   {"sq", "square"},
   {"squ", "square"},
   {"pe", "pentagonal"},
   {"pen", "pentagonal"},
   {"el", "elongated"},
   {"ge", "gyroelongated"},
   {"tr", "truncated"},
   {"au", "augmented"},
   {"ba", "biaugmenbted"},
   {"ta", "triaugmented"},
};



int make_resource_johnson(geom_if &geom, string name, char *errmsg)
{
   if(name.size()<2 || !strchr("jJ", name[0]) || name.find('.')!=string::npos)
      return -1; // not johnson name

   j_poly json;
   int sym_no;
   if(read_int(name.c_str()+1, &sym_no)) {
      sym_no--;
      if(sym_no<0 || sym_no >= json.get_last_J()) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "Johnson polyhedron number out of range");
         return 1; // fail
      }
   }
   else if(name[1] == '_') {
      string expanded = expand_abbrevs(name, j_abbrevs,
            sizeof(j_abbrevs)/sizeof(j_abbrevs[0]));
      sym_no = json.lookup_sym_no(expanded.c_str());
      if(sym_no == -1) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "invalid Johnson polyhedron name");
         return 1; // fail
      }
   }
   else
      return -1; // not johnson name
   
   json.get_poly(geom, sym_no);
   return 0; // name found
}


bool make_resource_geom(geom_if &geom, string name, char *errmsg)
{
   if(errmsg)
      *errmsg = '\0';
   geom.clear_all();
   
   if(!name.size() || name.find('.')!=string::npos)
      return false;

   // Look for an internal alternative name
   char alt_name[MSG_SZ];
   to_resource_name(alt_name, name.c_str());

   for(unsigned int i=0; i<sizeof(alt_names)/sizeof(alt_names[0]); i++) {
      if(strcmp(alt_name, alt_names[i][0])==0) {
         name = alt_names[i][1];  // set name to the usual name for the model
         break;
      }
   }

   char errmsg2[MSG_SZ];
   bool geom_ok = false;
   
   if(!geom_ok) {
      int ret = make_resource_pgon(geom, name, errmsg2);
      if(ret==0)
         geom_ok = true;
      else if(ret > 0) {
         if(errmsg)
            strcpy(errmsg, errmsg2);
         return false;
      }
   }

   if(!geom_ok) {
      int ret = make_resource_uniform(geom, name, errmsg2);
      if(ret==0)
         geom_ok = true;
      else if(ret > 0) {
         if(errmsg)
            strcpy(errmsg, errmsg2);
         return false;
      }
   }

   if(!geom_ok) {
      int ret = make_resource_johnson(geom, name, errmsg2);
      if(ret==0)
         geom_ok = true;
      else if(ret > 0) {
         if(errmsg)
            strcpy(errmsg, errmsg2);
         return false;
      }
   }


   if(geom_ok)
      if(col_geom_v *cg = dynamic_cast<col_geom_v *>(&geom)) {
         coloring clrng(cg);
         color_map_map *overrides = new color_map_map;
         overrides->set_col(60, col_val(0.9,0.45,0.0)); // triangle
         overrides->set_col(36, col_val(0.7,0.1,0.2));  // pentagram
         clrng.add_cmap(overrides);
         color_map *cmap = init_color_map("spread+53*12");
         clrng.add_cmap(cmap);

         clrng.f_avg_angle(true);
      }

   return geom_ok;
}

 

