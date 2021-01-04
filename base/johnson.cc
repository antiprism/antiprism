/*
   Copyright (c) 2006-2021, Adrian Rossiter

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

/*
   Name: johnson.cc
   Description: johnson polyhedra
   Project: Antiprism - http://www.antiprism.com
*/

/* These notes are from this thread on the internet: http://snipurl.com/hh13

   Norman Johnson says there are 17 simple Johnson solids (although he uses
   the word elementary rather than simple).

   (1-6, 63, 80, 83-86, 88-92).

   However, as Paul Gailunas has pointed out, there are 11 more elementary
   convex solids whose faces are regular polygons:
   tetrahedron
   dodecahedron
   truncated tetrahedron
   truncated cube
   truncated dodecahedron
   truncated icosahedron
   rhombitruncated cuboctahedron
   rhombitruncated icosidodecahedron
   snub cube
   snub dodecahedron
   truncated octahedron

   The cube's not in the above list because it's a prism.
*/

#include "polygon.h"
#include "private_std_polys.h"

#include <cctype>
#include <string>
#include <vector>

using std::string;
using std::vector;

using namespace anti;

void unit_edge(Geometry &geom)
{
  Vec3d v0 = geom.face_v(0, 0);
  Vec3d v1 = geom.face_v(0, 1);
  geom.transform(Trans3d::scale(1.0 / (v1 - v0).len()));
}

// face bond two J_polyhedra
void bond(Geometry &geom, model_func base, model_func brick, int f = 0,
          int b_f = 0, int off = 0)
{
  base(geom);
  Geometry geom2;
  brick(geom2);
  geom2.clear_cols();
  face_bond(geom, geom2, f, b_f, off);
}

// face bond two J_polyhedra, with a fixed twist
void tbond(Geometry &geom, model_func base, model_func brick, int f = 0,
           int b_f = 0)
{
  bond(geom, base, brick, f, b_f, 1);
}

// ------------------------------------------------------------------------
// Non-Johnson elementary parts

void J_tetrahedron(Geometry &geom)
{
  uni_pgon(geom, Polygon::as_type(Polygon::pyramid, 3));
}
void J_icosahedron(Geometry &geom)
{
  geom.read_resource("std_ico");
  unit_edge(geom);
}
void J_dodecahedron(Geometry &geom)
{
  geom.read_resource("std_dod");
  unit_edge(geom);
}
void J_icosidodecahedron(Geometry &geom)
{
  geom.read_resource("std_icosid");
  unit_edge(geom);
}
void J_tr_tetrahedron(Geometry &geom)
{
  geom.read_resource("std_tr_tet");
  unit_edge(geom);
}
void J_tr_cube(Geometry &geom)
{
  geom.read_resource("std_tr_cube");
  unit_edge(geom);
}
void J_tr_icosahedron(Geometry &geom)
{
  geom.read_resource("std_tr_icosa");
  unit_edge(geom);
}
void J_tr_dodecahedron(Geometry &geom)
{
  geom.read_resource("std_tr_dod");
  unit_edge(geom);
}
void J_rhombicosidodecahedron(Geometry &geom)
{
  geom.read_resource("std_rhombicosid");
  unit_edge(geom);
}
void J_prism3(Geometry &geom)
{
  uni_pgon(geom, Polygon::as_type(Polygon::prism, 3));
}
void J_prism4(Geometry &geom)
{
  uni_pgon(geom, Polygon::as_type(Polygon::prism, 4));
}
void J_prism5(Geometry &geom)
{
  uni_pgon(geom, Polygon::as_type(Polygon::prism, 5));
}
void J_prism6(Geometry &geom)
{
  uni_pgon(geom, Polygon::as_type(Polygon::prism, 6));
}
void J_prism8(Geometry &geom)
{
  uni_pgon(geom, Polygon::as_type(Polygon::prism, 8));
}
void J_prism10(Geometry &geom)
{
  uni_pgon(geom, Polygon::as_type(Polygon::prism, 10));
}
void J_antiprism4(Geometry &geom)
{
  uni_pgon(geom, Polygon::as_type(Polygon::antiprism, 4));
}
void J_antiprism5(Geometry &geom)
{
  uni_pgon(geom, Polygon::as_type(Polygon::antiprism, 5));
}
void J_antiprism6(Geometry &geom)
{
  uni_pgon(geom, Polygon::as_type(Polygon::antiprism, 6));
}
void J_antiprism8(Geometry &geom)
{
  uni_pgon(geom, Polygon::as_type(Polygon::antiprism, 8));
}
void J_antiprism10(Geometry &geom)
{
  uni_pgon(geom, Polygon::as_type(Polygon::antiprism, 10));
}
// ------------------------------------------------------------------------

// square pyramid
void J1(Geometry &geom)
{
  uni_pgon(geom, Polygon::as_type(Polygon::pyramid, 4));
}

// pentagonal pyramid
void J2(Geometry &geom)
{
  uni_pgon(geom, Polygon::as_type(Polygon::pyramid, 5));
}

// triangular cupola
void J3(Geometry &geom)
{
  uni_pgon(geom, Polygon::as_type(Polygon::cupola, 3));
}

// square cupola
void J4(Geometry &geom)
{
  uni_pgon(geom, Polygon::as_type(Polygon::cupola, 4));
}

// pentagonal cupola
void J5(Geometry &geom)
{
  uni_pgon(geom, Polygon::as_type(Polygon::cupola, 5));
}

// pentagonal rotunda (elementary)
void J6(Geometry &geom)
{
  J_icosidodecahedron(geom);
  geom.transform(Trans3d::rotate(geom.face_cent(29), Vec3d(0, 1, 0)));
  vector<int> del_verts;
  for (unsigned int i = 0; i < geom.verts().size(); i++)
    if (geom.verts(i)[1] < -anti::epsilon)
      del_verts.push_back(i);
  geom.del(VERTS, del_verts);
  close_poly_basic(geom);
  swap(geom.faces(0), geom.faces(16)); // large face to 0, the def bonding face
}

// elongated triangular pyramid
void J7(Geometry &geom) { bond(geom, J_prism3, J_tetrahedron); }

// elongated square pyramid
void J8(Geometry &geom) { bond(geom, J_prism4, J1); }

// elongated pentagonal pyramid
void J9(Geometry &geom) { bond(geom, J_prism5, J2); }

// gyroelongated square pyramid
void J10(Geometry &geom) { bond(geom, J_antiprism4, J1); }

// gyroelongated pentagonal pyramid
void J11(Geometry &geom) { bond(geom, J_antiprism5, J2); }

// triangular dipyramid
void J12(Geometry &geom) { bond(geom, J_tetrahedron, J_tetrahedron); }

// pentagonal dipyramid
void J13(Geometry &geom) { bond(geom, J2, J2); }

// elongated triangular dipyramid
void J14(Geometry &geom) { bond(geom, J7, J_tetrahedron); }

// elongated square dipyramid
void J15(Geometry &geom) { bond(geom, J8, J1); }

// elongated pentagonal dipyramid
void J16(Geometry &geom) { bond(geom, J9, J2); }

// gyroelongated square dipyramid
void J17(Geometry &geom) { bond(geom, J10, J1); }

// elongated triangular cupola
void J18(Geometry &geom) { bond(geom, J_prism6, J3); }

// elongated square cupola
void J19(Geometry &geom) { bond(geom, J_prism8, J4); }

// elongated pentagonal cupola
void J20(Geometry &geom) { bond(geom, J_prism10, J5); }

// elongated pentagonal rotunda
void J21(Geometry &geom) { bond(geom, J_prism10, J6); }

// gyroelongated triangular cupola
void J22(Geometry &geom) { bond(geom, J_antiprism6, J3); }

// gyroelongated square cupola
void J23(Geometry &geom) { bond(geom, J_antiprism8, J4); }

// gyroelongated pentagonal cupola
void J24(Geometry &geom) { bond(geom, J_antiprism10, J5); }

// gyroelongated pentagonal rotunda
void J25(Geometry &geom) { bond(geom, J_antiprism10, J6); }

// gyrobifastigium
void J26(Geometry &geom) { bond(geom, J_prism3, J_prism3, 4, 4); }

// triangular orthobicupola
void J27(Geometry &geom) { tbond(geom, J3, J3); }

// square orthobicupola
void J28(Geometry &geom) { tbond(geom, J4, J4); }

// square gyrobicupola
void J29(Geometry &geom) { bond(geom, J4, J4); }

// pentagonal orthobicupola
void J30(Geometry &geom) { tbond(geom, J5, J5); }

// pentagonal gyrobicupola
void J31(Geometry &geom) { bond(geom, J5, J5); }

// pentagonal orthocupolarotunda
void J32(Geometry &geom) { tbond(geom, J6, J5); }

// pentagonal gyrocupolarotunda
void J33(Geometry &geom) { bond(geom, J6, J5); }

// pentagonal orthobirotunda
void J34(Geometry &geom) { tbond(geom, J6, J6); }

// elongated triangular orthobicupola
void J35(Geometry &geom) { bond(geom, J18, J3); }

// elongated triangular gyrobicupola
void J36(Geometry &geom) { tbond(geom, J18, J3); }

// elongated square gyrobicupola
void J37(Geometry &geom) { tbond(geom, J19, J4); }

// elongated pentagonal orthobicupola
void J38(Geometry &geom) { bond(geom, J20, J5); }

// elongated pentagonal gyrobicupola
void J39(Geometry &geom) { tbond(geom, J20, J5); }

// elongated pentagonal orthocupolarotunda
void J40(Geometry &geom) { bond(geom, J21, J5); }

// elongated pentagonal gyrocupolarotunda
void J41(Geometry &geom) { tbond(geom, J21, J5); }

// elongated pentagonal orthobirotunda
void J42(Geometry &geom) { bond(geom, J21, J6); }

// elongated pentagonal gyrobirotunda
void J43(Geometry &geom) { tbond(geom, J21, J6); }

// gyroelongated triangular bicupola
void J44(Geometry &geom) { bond(geom, J22, J3); }

// gyroelongated square bicupola
void J45(Geometry &geom) { bond(geom, J23, J4); }

// gyroelongated pentagonal bicupola
void J46(Geometry &geom) { bond(geom, J24, J5); }

// gyroelongated pentagonal cupolarotunda
void J47(Geometry &geom) { bond(geom, J25, J5); }

// gyroelongated pentagonal birotunda
void J48(Geometry &geom) { bond(geom, J25, J6); }

// augmented triangular prism
void J49(Geometry &geom) { bond(geom, J_prism3, J1, 4); }

// biaugmented triangular prism
void J50(Geometry &geom) { bond(geom, J49, J1, 3); }

// triaugmented triangular prism
void J51(Geometry &geom) { bond(geom, J50, J1, 2); }

// augmented pentagonal prism
void J52(Geometry &geom) { bond(geom, J_prism5, J1, 4); }

// biaugmented pentagonal prism
void J53(Geometry &geom) { bond(geom, J52, J1, 2); }

// augmented hexagonal prism
void J54(Geometry &geom) { bond(geom, J_prism6, J1, 7); }

// parabiaugmented hexagonal prism
void J55(Geometry &geom) { bond(geom, J54, J1, 4); }

// metabiaugmented hexagonal prism
void J56(Geometry &geom) { bond(geom, J54, J1, 5); }

// triaugmented hexagonal prism
void J57(Geometry &geom) { bond(geom, J56, J1, 3); }

// augmented dodecahedron
void J58(Geometry &geom) { bond(geom, J_dodecahedron, J2, 6); }

// parabiaugmented dodecahedron
void J59(Geometry &geom) { bond(geom, J58, J2, 2); }

// metabiaugmented dodecahedron
void J60(Geometry &geom) { bond(geom, J58, J2, 4); }

// triaugmented dodecahedron
void J61(Geometry &geom) { bond(geom, J60, J2, 1); }

// tridiminished icosahedron (elementary)
void J63(Geometry &geom)
{
  J_icosahedron(geom);
  // largest first to avoid remapping
  geom.del(VERTS, 8);
  close_poly_basic(geom);
  geom.del(VERTS, 7);
  close_poly_basic(geom);
  geom.del(VERTS, 6);
  close_poly_basic(geom);
}

// metabidiminished icosahedron
void J62(Geometry &geom) { bond(geom, J63, J2, 6); }

// augmented tridiminished icosahedron
void J64(Geometry &geom) { bond(geom, J63, J_tetrahedron, 1); }

// augmented truncated tetrahedron
void J65(Geometry &geom) { bond(geom, J_tr_tetrahedron, J3, 7); }

// augmented truncated cube
void J66(Geometry &geom) { bond(geom, J_tr_cube, J4, 10); }

// biaugmented truncated cube
void J67(Geometry &geom) { bond(geom, J66, J4, 9); }

// augmented truncated dodecahedron
void J68(Geometry &geom) { bond(geom, J_tr_dodecahedron, J5, 25, 0, 1); }

// parabiaugmented truncated dodecahedron
void J69(Geometry &geom) { bond(geom, J68, J5, 23, 0, 1); }

// metabiaugmented truncated dodecahedron
void J70(Geometry &geom) { bond(geom, J68, J5, 25, 0, 1); }

// triaugmented truncated dodecahedron
void J71(Geometry &geom) { bond(geom, J70, J5, 24); }

// parabidiminished rhombicosidodecahedron (elementary)
void J80(Geometry &geom)
{
  J_rhombicosidodecahedron(geom);
  int dels[] = {4, 5, 44, 45, 58, 2, 3, 42, 43, 57};
  vector<int> del_verts;
  for (int &del : dels)
    del_verts.push_back(del);
  geom.del(VERTS, del_verts);
  close_poly_basic(geom);
}

// paragyrate diminished rhombicosidodecahedron
void J77(Geometry &geom) { bond(geom, J80, J5, 41); }

// gyrate rhombicosidodecahedron
void J72(Geometry &geom) { tbond(geom, J77, J5, 40); }

// parabigyrate rhombicosidodecahedron
void J73(Geometry &geom) { bond(geom, J77, J5, 40); }

// tridiminished rhombicosidodecahedron (elementary)
void J83(Geometry &geom)
{
  J_rhombicosidodecahedron(geom);
  int dels[] = {4, 5, 44, 45, 58, 16, 18, 32, 34, 52, 6, 7, 46, 47, 59};
  vector<int> del_verts;
  for (int &del : dels)
    del_verts.push_back(del);
  geom.del(VERTS, del_verts);
  close_poly_basic(geom);
}

// gyrate bidiminished rhombicosidodecahedron
void J82(Geometry &geom) { tbond(geom, J83, J5, 31); }

// bigyrate diminished rhombicosidodecahedron
void J79(Geometry &geom) { tbond(geom, J82, J5, 30); }

// metabigyrate rhombicosidodecahedron
void J74(Geometry &geom) { tbond(geom, J79, J5, 29); }

// trigyrate rhombicosidodecahedron
void J75(Geometry &geom) { bond(geom, J79, J5, 29); }

// diminished rhombicosidodecahedron
void J76(Geometry &geom) { tbond(geom, J80, J5, 41); }

// metagyrate diminished rhombicosidodecahedron
void J78(Geometry &geom) { bond(geom, J82, J5, 30); }

// metabidiminished rhombicosidodecahedron
void J81(Geometry &geom) { bond(geom, J83, J5, 31); }

// snub disphenoid (elementary)
void J84(Geometry &geom)
{
  geom.read_resource("std_snu2");

  /*
  geom.clear_all();
  double sol[3];

  double coeffs_x[] = { 8, -4, -3, 1 };
  cubic(coeffs_x, sol);
  double x = sol[1];         // 1.28916854644831

  double h = sqrt(4-2*x*x)/2;
  double H = sqrt(3-x*x) + h;

  vector<Vec3d> &verts = *geom.get_verts();
  verts.push_back(Vec3d(1, -H, 0));
  verts.push_back(Vec3d(0, -h, x));
  verts.push_back(Vec3d(x, h, 0));
  verts.push_back(Vec3d(0, H, 1));
  for(int i=0; i<4; i++) {
     verts[i] *= 1/2.0;
     verts.push_back(Vec3d(-verts[i][0], verts[i][1], -verts[i][2]));
  }

  geom.add_hull();
  normalised_face_list(geom);
  */
}

// snub square antiprism (elementary)
void J85(Geometry &geom)
{
  geom.read_resource("std_snu4");

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

  vector<Vec3d> &verts = *geom.get_verts();
  vector<Vec3d> vs(4);
  vs[0] = Vec3d(sqrt(2)/2, H, 0);
  vs[1] = Vec3d(r, h, 0);
  vs[2] = Vec3d(0.5, -H, 0.5);
  vs[3] = Vec3d(r/sqrt(2), -h, r/sqrt(2));
  for(int i=0; i<4; i++) {
        verts.push_back(Vec3d(vs[i][0], vs[i][1], vs[i][2]));
        verts.push_back(Vec3d(-vs[i][0], vs[i][1], -vs[i][2]));
        verts.push_back(Vec3d(-vs[i][2], vs[i][1], vs[i][0]));
        verts.push_back(Vec3d(vs[i][2], vs[i][1], -vs[i][0]));
  }

  geom.add_hull();
  normalised_face_list(geom);
  */
}

// sphenocorona (elementary)
void J86(Geometry &geom)
{
  // solution of the coordinate distance equations by Wolfram|Alpha
  // http://www.wolframalpha.com/
  double A = (9.0 - sqrt(6.0) + 2.0 * sqrt(3.0 * (71.0 - 19.0 * sqrt(6.0))));
  double a = 1.0 / 15.0 * A;
  double b = sqrt(3.0 + 2.0 / 15.0 * A - 1.0 / 225.0 * A * A);
  double c = 1.0 / 900.0 *
             (-65.0 * b - 124.0 * A * b + 1444.0 / 225.0 * A * A * b +
              548.0 / 3375.0 * A * A * A * b +
              b * (2907.0 + 204.0 * A - 1394.0 / 125.0 * A * A -
                   204.0 / 625.0 * A * A * A) +
              b * (-912.0 - 268.0 / 3.0 * A + 1936.0 / 1125.0 * A * A +
                   3368.0 / 16875.0 * A * A * A));

  double d = 49.0 / 162.0 + 181.0 / 2430.0 * A - 47.0 / 60750.0 * A * A +
             1.0 / 182250.0 * A * A * A;

  double e = 1.0 / 3240.0 *
             (8855.0 * b - 77.0 / 3.0 * A * b - 611.0 / 75.0 * A * A * b +
              13.0 / 225.0 * A * A * A * b);

  double crds[] = {0,  1, e, 0,  -1, e, a, 0,  b, -a, 0, b, 1, d, c,
                   -1, d, c, -1, -d, c, 1, -d, c, -1, 0, 0, 1, 0, 0};

  for (int i = 0; i < 10; i++)
    geom.add_vert(0.5 * Vec3d(crds[i * 3], crds[i * 3 + 1], crds[i * 3 + 2]));

  geom.add_hull();
  normalised_face_list(geom);
}

// augmented sphenocorona
void J87(Geometry &geom) { bond(geom, J86, J1, 7); }

// clang-format off
// sphenomegacorona (elementary)
void J88(Geometry &geom)
{
  // there is no known rational formula for sphenomegacorona (j88) coordinates

  // coordinates are from njohnson.py (see antiprism_python in github)
  double a = 0.5;
  double b = 0.5946333356326384;
  double c = 0.8039970125282817;
  double d = 1.2831023388312692;
  double e = 0.6218928580688124;
  double f = 1.6648364469102344;
  double g = 0.8547430824889655;
  double h = 1.5255013725848436;

  double crds[] = {
      0.0,a,0.0,
      0.0,-a,0.0,
      -b,a,c,
      b,-a,c,
      b,a,c,
      -b,-a,c,
      0.0,d,e,
      0.0,-d,e,
      -a,0.0,f,
      a,0.0,f,
      0.0,g,h,
      0.0,-g,h};

  for (int i = 0; i < 12; i++)
    geom.add_vert(Vec3d(crds[i * 3], crds[i * 3 + 1], crds[i * 3 + 2]));

  geom.add_hull();
  normalised_face_list(geom);
}

// hebesphenomegacorona (elementary)
void J89(Geometry &geom)
{
  // there is no known rational formula for hebesphenomegacorona (j89) coordinates

  // coordinates are from njohnson.py (see antiprism_python in github)
  double a = 0.5;
  double b = 0.7168448157134573;
  double c = 0.9762060878207003;
  double d = 1.1012960420472777;
  double e = 0.6232520114835398;
  double f = 1.8146441152849828;
  double g = 0.8356590717223666;
  double h = 1.5873251414880025;

  double crds[] = {
      a,a,0.0,
      -a,a,0.0,
      -a,-a,0.0,
      a,-a,0.0,
      a,b,c,
      -a,-b,c,
      -a,b,c,
      a,-b,c,
      d,0.0,e,
      -d,0.0,e,
      0.0,a,f,
      0.0,-a,f,
      g,0.0,h,
      -g,0.0,h};

  for (int i = 0; i < 14; i++)
    geom.add_vert(Vec3d(crds[i * 3], crds[i * 3 + 1], crds[i * 3 + 2]));

  geom.add_hull();
  normalised_face_list(geom);
}

// disphenocingulum (elementary)
void J90(Geometry &geom)
{
  // there is no known rational formula for disphenocingulum (j90) coordinates

  // coordinates are from njohnson.py (see antiprism_python in github)
  double a = 0.3535533905932738;
  double b = 1.104437942079934;
  double c = 0.1888902221636221;
  double d = 0.8959970033501695;
  double e = 0.4629476039153646;
  double f = 0.7965438721919097;
  double g = 0.3250029759499704;

  double crds[] = {
      a,a,-b,
      -a,-a,-b,
      -c,d,-e,
      c,-d,-e,
      d,-c,-e,
      -d,c,-e,
      f,f,-g,
      -f,-f,-g,
      a,-a,b,
      -a,a,b,
      d,c,e,
      -d,-c,e,
      -c,-d,e,
      c,d,e,
      f,-f,g,
      -f,f,g};

  for (int i = 0; i < 16; i++)
    geom.add_vert(Vec3d(crds[i * 3], crds[i * 3 + 1], crds[i * 3 + 2]));

  geom.add_hull();
  normalised_face_list(geom);
}
// clang-format on

// bilunabirotunda (elementary)
void J91(Geometry &geom)
{
  geom.clear_all();
  double phi = (sqrt(5) + 1) / 2;
  geom.add_vert(Vec3d(phi, 1, 1));
  geom.add_vert(Vec3d(-phi, 1, 1));
  geom.add_vert(Vec3d(phi, 1, -1));
  geom.add_vert(Vec3d(-phi, 1, -1));
  geom.add_vert(Vec3d(1, 1 + phi, 0));
  geom.add_vert(Vec3d(-1, 1 + phi, 0));
  geom.add_vert(Vec3d(0, 0, phi));
  for (int i = 0; i < 7; i++) {
    geom.verts(i) /= 2.0;
    geom.add_vert(-geom.verts(i));
  }
  geom.add_hull();
  normalised_face_list(geom);
}

// triangular hebesphenorotunda (elementary)
void J92(Geometry &geom)
{
  // hexagon is base at height h2.
  // irregular hexagon section though edges and pentagon diagonals at
  // height 0 has sides alternating 1 and phi. The vertices lie
  // at distance r from the main axis.
  // top triangle is at height h1.
  geom.clear_all();
  double phi = (sqrt(5) + 1) / 2;
  double r_2 = 2 * phi * phi / 3; // r^2
  double x = sqrt(r_2 - 1 / 4.0); // axis to mid-edge on irregular hexagon
  double x1 = 1 / sqrt(3);        // axis to mid-edge on top triangle
  double x2 = sqrt(3) / 2;        // axis to mid-edge on base hexagon

  double h1 = sqrt(3 / 4.0 - (x - x1) * (x - x1)); // height of top tri
  double h2 = -sqrt(1 - (x - x2) * (x - x2));      // height of base hex

  Vec3d tri_mid(-x1 / 2, h1, 0);                      // top tri mid edge
  Vec3d phi_mid(-phi * sqrt(5 / 12.0), 0, 0);         // pentagon mid-diagonal
  Vec3d p_apex = tri_mid + phi * (phi_mid - tri_mid); // pentagon apex

  vector<Vec3d> vs(8);
  vs[0] = Vec3d(x1, h1, 0);        // top triangle
  vs[1] = Vec3d(x, 0, 1 / 2.0);    // irregular hexagon
  vs[2] = Vec3d(x, 0, -1 / 2.0);   // irregular hexagon
  vs[3] = p_apex;                  // pentagon apex
  vs[4] = Vec3d(x2, h2, 1 / 2.0);  // irregular hexagon
  vs[5] = Vec3d(x2, h2, -1 / 2.0); // irregular hexagon
  for (int i = 0; i < 6; i++) {
    double rot120_x = -vs[i][0] / 2 - vs[i][2] * sqrt(3) / 2;
    double rot120_z = vs[i][0] * sqrt(3) / 2 - vs[i][2] / 2;
    geom.add_vert(vs[i]);
    geom.add_vert(Vec3d(rot120_x, vs[i][1], rot120_z));
    geom.add_vert(Vec3d(rot120_x, vs[i][1], -rot120_z));
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
    {J92, "triangular hebesphenorotunda"}};

Johnson::Johnson()
{
  J_items = j_item_list;
  last_J = sizeof(j_item_list) / sizeof(j_item_list[0]);
}

int Johnson::get_poly(Geometry &geom, int sym)
{
  J_items[sym].pfunc(geom);
  return 1;
}

int Johnson::lookup_sym_no(string sym)
{
  // remove double spaces and spaces at beginning and end
  string sym_norm;
  bool ignore_if_space = true;
  for (char i : sym) {
    if (i == ' ') {
      if (ignore_if_space)
        continue;
      else
        ignore_if_space = true;
    }
    else
      ignore_if_space = false;
    sym_norm += i;
  }

  if (sym_norm[sym_norm.size() - 1] == ' ')
    sym_norm.resize(sym_norm.size() - 1);

  // is it blank
  if (sym_norm == "")
    return -1;

  // is it a J number
  int offset = (sym_norm[0] == 'j' || sym_norm[0] == 'J');
  char *endptr;
  int idx = strtol(sym_norm.c_str() + offset, &endptr, 10);
  if (!*endptr) { // all of string is an integer
    idx -= 1;
    if (idx < 0 || idx >= last_J) // out of range
      return -1;
    else
      return idx - 1;
  }

  // is it a poly name
  idx = -1;
  for (char &i : sym_norm)
    if (isalpha(i))
      i = tolower(i);
  for (int i = 0; i < last_J; i++) {
    if (sym_norm == J_items[i].name)
      return i;

    if (idx < 0 &&
        strncmp(sym_norm.c_str(), J_items[i].name, sym_norm.size()) == 0)
      idx = i;
  }

  return idx;
}
