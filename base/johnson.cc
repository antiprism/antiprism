/*
   Copyright (c) 2006-2009, Adrian Rossiter

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
   Name: johnson.cc
   Description: johnson polyhedra 
   Project: Antiprism - http://www.antiprism.com
*/

#include "std_polys.h"
#include "polygons.h"

void unit_edge(geom_if &geom)
{
   vec3d v0 = (*geom.get_verts())[(*geom.get_faces())[0][0]];
   vec3d v1 = (*geom.get_verts())[(*geom.get_faces())[0][1]];
   geom.transform(mat3d::scale(1.0/(v1-v0).mag()));
}


// face bond two J_polyhedra
void bond(geom_if &geom, model_func base, model_func brick,
      int f=0, int b_f=0, int off=0)
{
   base(geom);
   geom_v geom2;
   brick(geom2);
   face_bond(geom, geom2, f, b_f, off);
}

// face bond two J_polyhedra, with a fixed twist
void tbond(geom_if &geom, model_func base, model_func brick, int f=0, int b_f=0)
{ bond(geom, base, brick, f, b_f, 1); }

// ------------------------------------------------------------------------
// Non-Johnson elementary parts

void J_tetrahedron(geom_if &geom)
                  { uni_pgon(geom, pyramid(3)); }
void J_icosahedron(geom_if &geom)
                  { icosahedron(geom); unit_edge(geom); }
void J_dodecahedron(geom_if &geom)
                  { dodecahedron(geom); unit_edge(geom); }
void J_icosidodecahedron(geom_if &geom)
                  { icosidodecahedron(geom); unit_edge(geom); }
void J_tr_tetrahedron(geom_if &geom)
                  { tr_tetrahedron(geom); unit_edge(geom); }
void J_tr_cube(geom_if &geom)
                  { tr_cube(geom); unit_edge(geom); }
void J_tr_icosahedron(geom_if &geom)
                  { tr_icosahedron(geom); unit_edge(geom); }
void J_tr_dodecahedron(geom_if &geom)
                  { tr_dodecahedron(geom); unit_edge(geom); }
void J_rhombicosidodecahedron(geom_if &geom)
                  { rhombicosidodecahedron(geom); unit_edge(geom); }
void J_prism3(geom_if &geom) { uni_pgon(geom, prism(3)); }
void J_prism4(geom_if &geom) { uni_pgon(geom, prism(4)); }
void J_prism5(geom_if &geom) { uni_pgon(geom, prism(5)); }
void J_prism6(geom_if &geom) { uni_pgon(geom, prism(6)); }
void J_prism8(geom_if &geom) { uni_pgon(geom, prism(8)); }
void J_prism10(geom_if &geom) { uni_pgon(geom, prism(10)); }
void J_antiprism4(geom_if &geom) { uni_pgon(geom, antiprism(4)); }
void J_antiprism5(geom_if &geom) { uni_pgon(geom, antiprism(5)); }
void J_antiprism6(geom_if &geom) { uni_pgon(geom, antiprism(6)); }
void J_antiprism8(geom_if &geom) { uni_pgon(geom, antiprism(8)); }
void J_antiprism10(geom_if &geom) { uni_pgon(geom, antiprism(10)); }
// ------------------------------------------------------------------------


// square pyramid
void J1(geom_if &geom) { uni_pgon(geom, pyramid(4)); }

// pentagonal pyramid
void J2(geom_if &geom) { uni_pgon(geom, pyramid(5)); }

// triangular cupola
void J3(geom_if &geom) { uni_pgon(geom, cupola(3)); }

// square cupola
void J4(geom_if &geom) { uni_pgon(geom, cupola(4)); }

// pentagonal cupola
void J5(geom_if &geom) { uni_pgon(geom, cupola(5)); }

// pentagonal rotunda
void J6(geom_if &geom)
{
   J_icosidodecahedron(geom);
   vector<vec3d> &verts = *geom.get_verts();
   vector<vector<int> > &faces = *geom.get_faces();
   geom.transform(mat3d::rot(centroid(verts, faces[0]), vec3d(0,1,0)));
   vector<int> del_verts;
   for(unsigned int i=0; i<verts.size(); i++)
      if(verts[i][1]<-epsilon)
         del_verts.push_back(i);
   geom.delete_verts(del_verts);
   close_poly_basic(geom);
   swap(faces[0], faces[16]); // large face to 0, the default bonding face
}

// elongated triangular pyramid
void J7(geom_if &geom) { bond(geom, J_prism3, J_tetrahedron); }

// elongated square pyramid
void J8(geom_if &geom) { bond(geom, J_prism4, J1); }

// elongated pentagonal pyramid
void J9(geom_if &geom) { bond(geom, J_prism5, J2); }

// gyroelongated square pyramid
void J10(geom_if &geom) { bond(geom, J_antiprism4, J1); }

// gyroelongated pentagonal pyramid
void J11(geom_if &geom) { bond(geom, J_antiprism5, J2); }

// triangular dipyramid
void J12(geom_if &geom) { bond(geom, J_tetrahedron, J_tetrahedron); }

// pentagonal dipyramid
void J13(geom_if &geom) { bond(geom, J2, J2); }

// elongated triangular dipyramid
void J14(geom_if &geom) { bond(geom, J7, J_tetrahedron); }

// elongated square dipyramid
void J15(geom_if &geom) { bond(geom, J8, J1); }

// elongated pentagonal dipyramid
void J16(geom_if &geom) { bond(geom, J9, J2); }

// gyroelongated square dipyramid
void J17(geom_if &geom) { bond(geom, J10, J1); }

// elongated triangular cupola
void J18(geom_if &geom) { bond(geom, J_prism6, J3); }

// elongated square cupola
void J19(geom_if &geom) { bond(geom, J_prism8, J4); }

// elongated pentagonal cupola
void J20(geom_if &geom) { bond(geom, J_prism10, J5); }

// elongated pentagonal rotunda
void J21(geom_if &geom) { bond(geom, J_prism10, J6); }

// gyroelongated triangular cupola
void J22(geom_if &geom) { bond(geom, J_antiprism6, J3); }

// gyroelongated square cupola
void J23(geom_if &geom) { bond(geom, J_antiprism8, J4); }

// gyroelongated pentagonal cupola
void J24(geom_if &geom) { bond(geom, J_antiprism10, J5); }

// gyroelongated pentagonal rotunda
void J25(geom_if &geom) { bond(geom, J_antiprism10, J6); }

// gyrobifastigium
void J26(geom_if &geom) { bond(geom, J_prism3, J_prism3, 4, 4); }

// triangular orthobicupola
void J27(geom_if &geom) { tbond(geom, J3, J3); }

// square orthobicupola
void J28(geom_if &geom) { tbond(geom, J4, J4); }

// square gyrobicupola
void J29(geom_if &geom) { bond(geom, J4, J4); }

// pentagonal orthobicupola
void J30(geom_if &geom) { tbond(geom, J5, J5); }

// pentagonal gyrobicupola
void J31(geom_if &geom) { bond(geom, J5, J5); }

// pentagonal orthocupolarotunda
void J32(geom_if &geom) { tbond(geom, J6, J5); }

// pentagonal gyrocupolarotunda
void J33(geom_if &geom) { bond(geom, J6, J5); }

// pentagonal orthobirotunda
void J34(geom_if &geom) { tbond(geom, J6, J6); }

// elongated triangular orthobicupola
void J35(geom_if &geom) { bond(geom, J18, J3); }

// elongated triangular gyrobicupola
void J36(geom_if &geom) { tbond(geom, J18, J3); }

// elongated square gyrobicupola
void J37(geom_if &geom) { tbond(geom, J19, J4); }

// elongated pentagonal orthobicupola
void J38(geom_if &geom) { bond(geom, J20, J5); }

// elongated pentagonal gyrobicupola
void J39(geom_if &geom) { tbond(geom, J20, J5); }

// elongated pentagonal orthocupolarotunda
void J40(geom_if &geom) { bond(geom, J21, J5); }

// elongated pentagonal gyrocupolarotunda
void J41(geom_if &geom) { tbond(geom, J21, J5); }

// elongated pentagonal orthobirotunda
void J42(geom_if &geom) { bond(geom, J21, J6); }

// elongated pentagonal gyrobirotunda
void J43(geom_if &geom) { tbond(geom, J21, J6); }

// gyroelongated triangular bicupola
void J44(geom_if &geom) { bond(geom, J22, J3); }

// gyroelongated square bicupola
void J45(geom_if &geom) { bond(geom, J23, J4); }

// gyroelongated pentagonal bicupola
void J46(geom_if &geom) { bond(geom, J24, J5); }

// gyroelongated pentagonal cupolarotunda
void J47(geom_if &geom) { bond(geom, J25, J5); }

// gyroelongated pentagonal birotunda
void J48(geom_if &geom) { bond(geom, J25, J6); }

// augmented triangular prism
void J49(geom_if &geom) { bond(geom, J_prism3, J1, 4); }

// biaugmented triangular prism
void J50(geom_if &geom) { bond(geom, J49, J1, 3); }

// triaugmented triangular prism
void J51(geom_if &geom) { bond(geom, J50, J1, 2); }

// augmented pentagonal prism
void J52(geom_if &geom) { bond(geom, J_prism5, J1, 4); }

// biaugmented pentagonal prism
void J53(geom_if &geom) { bond(geom, J52, J1, 2); }

// augmented hexagonal prism
void J54(geom_if &geom) { bond(geom, J_prism6, J1, 7); }

// parabiaugmented hexagonal prism
void J55(geom_if &geom) { bond(geom, J54, J1, 4); }

// metabiaugmented hexagonal prism
void J56(geom_if &geom) { bond(geom, J54, J1, 5); }

// triaugmented hexagonal prism
void J57(geom_if &geom) { bond(geom, J56, J1, 3); }

// augmented dodecahedron
void J58(geom_if &geom) { bond(geom, J_dodecahedron, J2, 9); }

// parabiaugmented dodecahedron
void J59(geom_if &geom) { bond(geom, J58, J2, 0); }

// metabiaugmented dodecahedron
void J60(geom_if &geom) { bond(geom, J58, J2, 3); }

// triaugmented dodecahedron
void J61(geom_if &geom) { bond(geom, J60, J2, 2); }

// tridiminished icosahedron
void J63(geom_if &geom)
{
   J_icosahedron(geom);
   // largest first to avoid remapping
   geom.delete_vert(9);
   close_poly_basic(geom);
   geom.delete_vert(8);
   close_poly_basic(geom);
   geom.delete_vert(1);
   close_poly_basic(geom);
}

// metabidiminished icosahedron
void J62(geom_if &geom) { bond(geom, J63, J2, 5); }

// augmented tridiminished icosahedron
void J64(geom_if &geom) { bond(geom, J63, J_tetrahedron); }

// augmented truncated tetrahedron
void J65(geom_if &geom) { bond(geom, J_tr_tetrahedron, J3); }

// augmented truncated cube
void J66(geom_if &geom) { bond(geom, J_tr_cube, J4, 4); }

// biaugmented truncated cube
void J67(geom_if &geom) { bond(geom, J66, J4, 4); }

// augmented truncated dodecahedron
void J68(geom_if &geom) { bond(geom, J_tr_dodecahedron, J5); }

// parabiaugmented truncated dodecahedron
void J69(geom_if &geom) { bond(geom, J68, J5, 8); }

// metabiaugmented truncated dodecahedron
void J70(geom_if &geom) { bond(geom, J68, J5, 6); }

// triaugmented truncated dodecahedron
void J71(geom_if &geom) { bond(geom, J70, J5, 9); }

// parabidiminished rhombicosidodecahedron
void J80(geom_if &geom)
{
   J_rhombicosidodecahedron(geom);
   int dels[] = {4,5,8,20,19,  58,57,42,43,53};
   vector<int> del_verts;
   for(int i=0; i<10; i++)
      del_verts.push_back(dels[i]);
   geom.delete_verts(del_verts);
   close_poly_basic(geom);
}

// paragyrate diminished rhombicosidodecahedron
void J77(geom_if &geom) { bond(geom, J80, J5, 41); }

// gyrate rhombicosidodecahedron
void J72(geom_if &geom) { tbond(geom, J77, J5, 40); }

// parabigyrate rhombicosidodecahedron
void J73(geom_if &geom) { bond(geom, J77, J5, 40); }

// tridiminished rhombicosidodecahedron
void J83(geom_if &geom)
{
   J_rhombicosidodecahedron(geom);
   int dels[] = {4,5,8,20,19,  13,12,27,28,38,  9,10,44,25,24};
   vector<int> del_verts;
   for(int i=0; i<15; i++)
      del_verts.push_back(dels[i]);
   geom.delete_verts(del_verts);
   close_poly_basic(geom);
}

// gyrate bidiminished rhombicosidodecahedron
void J82(geom_if &geom) { tbond(geom, J83, J5, 30); }

// bigyrate diminished rhombicosidodecahedron
void J79(geom_if &geom) { tbond(geom, J82, J5, 30); }

// metabigyrate rhombicosidodecahedron
void J74(geom_if &geom) { tbond(geom, J79, J5, 29); }

// trigyrate rhombicosidodecahedron
void J75(geom_if &geom) { bond(geom, J79, J5, 29); }

// diminished rhombicosidodecahedron
void J76(geom_if &geom) { tbond(geom, J80, J5, 41); }

// metagyrate diminished rhombicosidodecahedron
void J78(geom_if &geom) { bond(geom, J82, J5, 30); }

// metabidiminished rhombicosidodecahedron
void J81(geom_if &geom) { bond(geom, J83, J5, 30); }

// snub disphenoid
void J84(geom_if &geom)
{
   make_resource_pgon(geom, "snu2");
  
   /*
   geom.clear_all();
   double sol[3];
   
   double coeffs_x[] = { 8, -4, -3, 1 };
   cubic(coeffs_x, sol);
   double x = sol[1];         // 1.28916854644831
   
   double h = sqrt(4-2*x*x)/2;
   double H = sqrt(3-x*x) + h;
 
   vector<vec3d> &verts = *geom.get_verts();
   verts.push_back(vec3d(1, -H, 0));
   verts.push_back(vec3d(0, -h, x));
   verts.push_back(vec3d(x, h, 0));
   verts.push_back(vec3d(0, H, 1));
   for(int i=0; i<4; i++) {
      verts[i] *= 1/2.0;
      verts.push_back(vec3d(-verts[i][0], verts[i][1], -verts[i][2]));
   }

   geom.add_hull();
   normalised_face_list(geom);
   */
}

// snub square antiprism
void J85(geom_if &geom)
{
   make_resource_pgon(geom, "snu4");
  
   // vertices arranged in four layers, middle layers have heights
   // h1 and h2 and vertices lie at radius r from the main axis
   // 2r^4 - 2sqrt(2)r^3 + (2sqrt(2)-7)r^2 + 2(1+sqrt(2))r + 1 = 0
   // h1 = sqrt(1 - (r-sqrt(1/2)^2)
   // h2 = sqrt(3/4 - (r-1/2)^2)
   
   /*
   geom.clear_all();
   double coeffs[] = { 1, 2*(1+sqrt(2)), 2*sqrt(2)-7, -2*sqrt(2), 2 };
   double sol[4];
   quartic(coeffs, sol);
   double r = sol[1]; // 1.213205545866313 calculation from quartic formula
   
   double h1 = sqrt(1 - pow(r-sqrt(1/2.0), 2));
   double h2 = sqrt(3/4.0 - pow(r-1/2.0, 2));
   double H = (h1+h2)/2;  // new top layer height
   double h = H-h1;       // new upper inner layer height

   vector<vec3d> &verts = *geom.get_verts();
   vector<vec3d> vs(4);
   vs[0] = vec3d(sqrt(2)/2, H, 0);
   vs[1] = vec3d(r, h, 0);
   vs[2] = vec3d(0.5, -H, 0.5);
   vs[3] = vec3d(r/sqrt(2), -h, r/sqrt(2));
   for(int i=0; i<4; i++) {
         verts.push_back(vec3d(vs[i][0], vs[i][1], vs[i][2]));
         verts.push_back(vec3d(-vs[i][0], vs[i][1], -vs[i][2]));
         verts.push_back(vec3d(-vs[i][2], vs[i][1], vs[i][0]));
         verts.push_back(vec3d(vs[i][2], vs[i][1], -vs[i][0]));
   }

   geom.add_hull();
   normalised_face_list(geom);
   */
}


// sphenocorona
void J86(geom_if &geom)
{
   // solution of the coordinate distance equations by Wolfram|Alpha
   // http://www.wolframalpha.com/
   double A = (9.0-sqrt(6.0)+2.0*sqrt(3.0*(71.0-19.0*sqrt(6.0))));  
   double a = 1.0/15.0*A;
   double b = sqrt(3.0 +2.0/15.0*A - 1.0/225.0*A*A);
   double c = 1.0/900.0*(
                 -65.0*b
                 -124.0*A*b
                 +1444.0/225.0*A*A*b
                 +548.0/3375.0*A*A*A*b
                 +b*(2907.0+204.0*A -1394.0/125.0*A*A -204.0/625.0*A*A*A)
                 +b*(-912.0 -268.0/3.0*A + 1936.0/1125.0*A*A
                    + 3368.0/16875.0*A*A*A));

   double d = 49.0/162.0 + 181.0/2430.0*A
              - 47.0/60750.0*A*A + 1.0/182250.0*A*A*A;

   double e = 1.0/3240.0*(
                 8855.0*b
                 -77.0/3.0*A*b
                 -611.0/75.0*A*A*b
                 +13.0/225.0*A*A*A*b);

   double crds[] = {
      0, 1, e,
      0, -1, e,
      a, 0, b,
      -a, 0, b,
      1, d, c,
      -1, d, c,
      -1, -d, c,
      1, -d, c,
      -1, 0, 0,
      1, 0, 0 };

   for(int i=0; i<10; i++)
         geom.add_vert(0.5*vec3d(crds[i*3], crds[i*3+1], crds[i*3+2]));

   geom.add_hull();
   normalised_face_list(geom);
}

// augmented sphenocorona
void J87(geom_if &geom) { bond(geom, J86, J1, 7); }

// sphenomegacorona
void J88(geom_if &geom)
{
   // use coordinates for now
   double crds[] = {
-0.49999999036257647, 0.80401400055670513, -0.59462870024862657,
-0.85474303991789613, 1.525545959813478, -2.8834319214551709e-05,
3.2400751663637069e-08, 1.664857813358968, -0.50003532206212287,
-1.2831023174837659, 0.62193745633186093, 1.3105631615658241e-05,
0.50000000963742264, 0.80401397706789002, -0.59462872713345483,
-0.49999999325924788, 4.4580538388935278e-05, 4.1941094472124707e-05,
5.9286669848764452e-08, 1.6649042144649679, 0.49996467686134521,
-0.49999995838805011, 0.8040691838455809, 0.59463796973636562,
0.85474312506003369, 1.525545919659673, -2.8880278456533051e-05,
0.50000000674075196, 4.4557049573881892e-05, 4.1914209643821378e-05,
0.500000041611949, 0.80406916035676579, 0.59463794285153726,
1.283102360178771, 0.62193739605475362, 1.303663964354256e-05 };
      
   vector<vec3d> &verts = *geom.get_verts();
   for(int i=0; i<12; i++)
         verts.push_back(vec3d(crds[i*3], crds[i*3+1], crds[i*3+2]));

   // align it
   vec3d mid1 = (verts[5]+verts[9])/2.0;
   vec3d mid2 = (verts[6]+verts[2])/2.0;
   geom.transform(mat3d::rot(mid2-mid1, vec3d(0,1,0)) * mat3d::transl(-mid1));
   vec3d v(verts[11][0], 0, verts[11][2]);
   geom.transform(mat3d::rot(v, vec3d(1,0,0)));
                 
   geom.add_hull();
   unit_edge(geom);
   normalised_face_list(geom);
}

// hebesphenomegacorona
void J89(geom_if &geom)
{
   // use coordinates for now
   double crds[] = {
0.50000064613933071, 0.97620515606108926, 0.71684542256978023,
1.1012961263957419, 0.62325160696545301, -1.239770378138855e-07,
0.83565925281672204, 1.587324763626395, 8.6130658756724206e-07,
0.50000038894916266, -7.5401639915758517e-07, 0.49999980668207922,
5.7105086951525469e-07, 1.8146434114432239, 0.50000166125003365,
-0.49999935386039468, 0.9762052564116076, 0.71684615685900477,
-0.49999961105056279, -6.5366588096466728e-07, 0.50000054097130375,
0.49999965466002039, 6.5661389633121093e-08, -0.50000019331731527,
0.49999959339660088, 0.97620633122463629, -0.71684420885626532,
-1.6323827283651201e-07, 1.8146442311210129, -0.49999833874936078,
-0.83565889062755272, 1.587324931344037, 2.0885374910873009e-06,
-1.1012959576982091, 0.62325182799671008, 1.493362595557812e-06,
-0.50000034533970505, 1.660119078365313e-07, -0.49999945902809068,
-0.5000004066031245, 0.97620643157515452, -0.71684347456704078 };
   
   vector<vec3d> &verts = *geom.get_verts();
   for(int i=0; i<14; i++)
         verts.push_back(vec3d(crds[i*3], crds[i*3+1], crds[i*3+2]));

   // align it
   vec3d mid1 = (verts[3]+verts[6]+verts[12]+verts[7])/4.0;
   vec3d mid2 = (verts[4]+verts[9])/2.0;
   geom.transform(mat3d::rot(mid2-mid1, vec3d(0,1,0)) * mat3d::transl(-mid1));
   vec3d v(verts[1][0], 0, verts[1][2]);
   geom.transform(mat3d::rot(v, vec3d(1,0,0)));
   
   geom.add_hull();
   unit_edge(geom);
   normalised_face_list(geom);
}

// disphenocingulum
void J90(geom_if &geom)
{
   // use coordinates for now
   double crds[] = {
-0.089880753431450133, 0.1646625164023697, -1.156160363750983,
0.61473692232257793, 0.60765411531969293, -0.60183881960565155,
-0.34978464587621588, 0.70710659443430013, -0.35728340895747512,
-0.74458688834140663, -0.54244461894771501, -0.88902860590803612,
0.65470819525853052, -0.37778008769805038, -0.76713167747514521,
-1.004490780786172, -5.4091578465311604e-07, -0.090151651114527998,
0.34978558799105458, 0.70710784735046706, 0.35728043555561162,
-0.61473580163880148, 0.6076540641487097, 0.60183602014396054,
2.0603485740081029e-06, -1.084887223048135, -0.49999991963219809,
1.0044929925172139, 1.420328567327977e-06, 0.090149914418578184,
-0.86029393220798744, -0.98543442698263761, 2.6913674317440649e-08,
0.86029791329427485, -0.98543288231409509, -4.0118625950398793e-08,
-0.65470530521950809, -0.37777992154306239, 0.76713060150445223,
0.089882669508380697, 0.16466469986945739, 1.156158339064705,
2.099306650979977e-06, -1.084886348564962, 0.50000008036741883,
0.74459007403453981, -0.54244172715244232, 0.8890278179276716 };

   vector<vec3d> &verts = *geom.get_verts();
   for(int i=0; i<16; i++)
         verts.push_back(vec3d(crds[i*3], crds[i*3+1], crds[i*3+2]));

   // align it
   vec3d mid1 = (verts[3]+verts[0])/2.0;
   vec3d mid2 = (verts[15]+verts[13])/2.0;
   geom.transform(mat3d::rot(mid2-mid1, vec3d(0,1,0)) * mat3d::transl(-mid1));
   vec3d v(verts[10][0], 0, verts[10][2]);
   geom.transform(mat3d::rot(v, vec3d(1,0,0)));
                 
   geom.add_hull();
   unit_edge(geom);
   normalised_face_list(geom);
}

// bilunabirotunda
void J91(geom_if &geom)
{
   geom.clear_all();
   double phi = (sqrt(5)+1)/2;
   vector<vec3d> &verts = *geom.get_verts();
   verts.push_back(vec3d( phi, 1, 1));
   verts.push_back(vec3d(-phi, 1, 1));
   verts.push_back(vec3d( phi, 1, -1));
   verts.push_back(vec3d(-phi, 1, -1));
   verts.push_back(vec3d( 1, 1+phi, 0));
   verts.push_back(vec3d(-1, 1+phi, 0));
   verts.push_back(vec3d(0, 0, phi));
   for(int i=0; i<7; i++) {
      verts[i] /= 2.0;
      verts.push_back(-verts[i]);
   }
   geom.add_hull();
   normalised_face_list(geom);
}

// triangular hebesphenorotunda
void J92(geom_if &geom)
{
   // hexagon is base at height h2.
   // irregular hexagon section though edges and pentagon diagonals at
   // height 0 has sides alternating 1 and phi. The vertices lie
   // at distance r from the main axis.
   // top triangle is at height h1.
   geom.clear_all();
   double phi = (sqrt(5)+1)/2;
   double r_2 = 2*phi*phi/3;           // r^2
   double x = sqrt(r_2-1/4.0);         // axis to mid-edge on irregular hexagon
   double x1 = 1/sqrt(3);              // axis to mid-edge on top triangle
   double x2 = sqrt(3)/2;              // axis to mid-edge on base hexagon

   double h1 =  sqrt(3/4.0 - (x-x1)*(x-x1));       // height of top tri
   double h2 = -sqrt(1 - (x-x2)*(x-x2));           // height of base hex
   
   vec3d tri_mid(-x1/2, h1, 0);                    // top tri mid edge
   vec3d phi_mid(-phi*sqrt(5/12.0), 0, 0);         // pentagon mid-diagonal 
   vec3d p_apex = tri_mid + phi*(phi_mid-tri_mid); // pentagon apex
   
   vector<vec3d> &verts = *geom.get_verts();
   vector<vec3d> vs(8);
   vs[0] = vec3d(x1, h1, 0);       // top triangle
   vs[1] = vec3d(x, 0, 1/2.0);     // irregular hexagon
   vs[2] = vec3d(x, 0, -1/2.0);    // irregular hexagon
   vs[3] = p_apex;                 // pentagon apex
   vs[4] = vec3d(x2, h2, 1/2.0);   // irregular hexagon
   vs[5] = vec3d(x2, h2, -1/2.0);  // irregular hexagon
   for(int i=0; i<6; i++) {
      double rot120_x = - vs[i][0]/2         - vs[i][2]*sqrt(3)/2;
      double rot120_z =   vs[i][0]*sqrt(3)/2 - vs[i][2]/2;
      verts.push_back(vs[i]);
      verts.push_back(vec3d(rot120_x, vs[i][1],  rot120_z));
      verts.push_back(vec3d(rot120_x, vs[i][1], -rot120_z));
   }

   geom.add_hull();
   normalised_face_list(geom);
}



JohnsonItem j_item_list[] = {
   {J1, "square pyramid"},
   {J2, "pentagonal pyramid"},
   {J3, "triangular cupola"},
   {J4, "square cupola"},
   {J5, "pentagonal cupola"},
   {J6, "pentagonal rotunda"},
   {J7, "elongated triangular pyramid"},
   {J8, "elongated square pyramid"},
   {J9, "elongated pentagonal pyramid"},
   {J10, "gyroelongated square pyramid"},
   {J11, "gyroelongated pentagonal pyramid"},
   {J12, "triangular dipyramid"},
   {J13, "pentagonal dipyramid"},
   {J14, "elongated triangular dipyramid"},
   {J15, "elongated square dipyramid"},
   {J16, "elongated pentagonal dipyramid"},
   {J17, "gyroelongated square dipyramid"},
   {J18, "elongated triangular cupola"},
   {J19, "elongated square cupola"},
   {J20, "elongated pentagonal cupola"},
   {J21, "elongated pentagonal rotunda"},
   {J22, "gyroelongated triangular cupola"},
   {J23, "gyroelongated square cupola"},
   {J24, "gyroelongated pentagonal cupola"},
   {J25, "gyroelongated pentagonal rotunda"},
   {J26, "gyrobifastigium"},
   {J27, "triangular orthobicupola"},
   {J28, "square orthobicupola"},
   {J29, "square gyrobicupola"},
   {J30, "pentagonal orthobicupola"},
   {J31, "pentagonal gyrobicupola"},
   {J32, "pentagonal orthocupolarotunda"},
   {J33, "pentagonal gyrocupolarotunda"},
   {J34, "pentagonal orthobirotunda"},
   {J35, "elongated triangular orthobicupola"},
   {J36, "elongated triangular gyrobicupola"},
   {J37, "elongated square gyrobicupola"},
   {J38, "elongated pentagonal orthobicupola"},
   {J39, "elongated pentagonal gyrobicupola"},
   {J40, "elongated pentagonal orthocupolarotunda"},
   {J41, "elongated pentagonal gyrocupolarotunda"},
   {J42, "elongated pentagonal orthobirotunda"},
   {J43, "elongated pentagonal gyrobirotunda"},
   {J44, "gyroelongated triangular bicupola"},
   {J45, "gyroelongated square bicupola"},
   {J46, "gyroelongated pentagonal bicupola"},
   {J47, "gyroelongated pentagonal cupolarotunda"},
   {J48, "gyroelongated pentagonal birotunda"},
   {J49, "augmented triangular prism"},
   {J50, "biaugmented triangular prism"},
   {J51, "triaugmented triangular prism"},
   {J52, "augmented pentagonal prism"},
   {J53, "biaugmented pentagonal prism"},
   {J54, "augmented hexagonal prism"},
   {J55, "parabiaugmented hexagonal prism"},
   {J56, "metabiaugmented hexagonal prism"},
   {J57, "triaugmented hexagonal prism"},
   {J58, "augmented dodecahedron"},
   {J59, "parabiaugmented dodecahedron"},
   {J60, "metabiaugmented dodecahedron"},
   {J61, "triaugmented dodecahedron"},
   {J62, "metabidiminished icosahedron"},
   {J63, "tridiminished icosahedron"},
   {J64, "augmented tridiminished icosahedron"},
   {J65, "augmented truncated tetrahedron"},
   {J66, "augmented truncated cube"},
   {J67, "biaugmented truncated cube"},
   {J68, "augmented truncated dodecahedron"},
   {J69, "parabiaugmented truncated dodecahedron"},
   {J70, "metabiaugmented truncated dodecahedron"},
   {J71, "triaugmented truncated dodecahedron"},
   {J72, "gyrate rhombicosidodecahedron"},
   {J73, "parabigyrate rhombicosidodecahedron"},
   {J74, "metabigyrate rhombicosidodecahedron"},
   {J75, "trigyrate rhombicosidodecahedron"},
   {J76, "diminished rhombicosidodecahedron"},
   {J77, "paragyrate diminished rhombicosidodecahedron"},
   {J78, "metagyrate diminished rhombicosidodecahedron"},
   {J79, "bigyrate diminished rhombicosidodecahedron"},
   {J80, "parabidiminished rhombicosidodecahedron"},
   {J81, "metabidiminished rhombicosidodecahedron"},
   {J82, "gyrate bidiminished rhombicosidodecahedron"},
   {J83, "tridiminished rhombicosidodecahedron"},
   {J84, "snub disphenoid"},
   {J85, "snub square antiprism"},
   {J86, "sphenocorona"},
   {J87, "augmented sphenocorona"},
   {J88, "sphenomegacorona"},
   {J89, "hebesphenomegacorona"},
   {J90, "disphenocingulum"},
   {J91, "bilunabirotunda"},
   {J92, "triangular hebesphenorotunda"}
};



j_poly::j_poly()
{ 
   J_items = j_item_list;
   last_J = sizeof (j_item_list) / sizeof (j_item_list[0]);
}


void j_poly::list_poly(int idx, FILE *fp)
{
   fprintf(fp, "%2d) %s\n", idx+1, J_items[idx].name);
}

void j_poly::list_polys(FILE *fp)
{
   for(int i=0; i<last_J; i++)
      list_poly(i, fp);
}

int j_poly::get_poly(geom_if &geom, int sym)
{
   J_items[sym].pfunc(geom);
   return 1;
}


int j_poly::lookup_sym_no(string sym)
{
   // remove double spaces and spaces at beginning and end
   string sym_norm;
   bool ignore_if_space = true;
   for(unsigned int i=0; i<sym.length(); i++) {
      if(sym[i]==' ') {
         if(ignore_if_space)
            continue;
         else
            ignore_if_space = true;
      }
      else
         ignore_if_space = false;
      sym_norm+=sym[i];
   }

   if(sym_norm[sym_norm.size()-1]==' ')
      sym_norm.resize(sym_norm.size()-1);
         
   // is it blank
   if(sym_norm=="")
      return -1;

   // is it a J number
   int offset = (sym_norm[0]=='j'||sym_norm[0]=='J');
   char *endptr;
   int idx = strtol(sym_norm.c_str()+offset, &endptr, 10);
   if(!*endptr)     // all of string is an integer
      return idx-1;

   // is it a poly name
   idx= -1;
   for(unsigned int i=0; i<sym_norm.size(); i++)
      if(isalpha(sym_norm[i]))
         sym_norm[i] = tolower(sym_norm[i]);
   for(int i=0; i<last_J; i++) {
      if(sym_norm==J_items[i].name)
         return i;
      
      if(idx<0 && strncmp(sym_norm.c_str(), J_items[i].name,sym_norm.size())==0)
         idx = i;
   }

   return idx;
}

