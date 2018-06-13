/*
   Copyright (c) 2003-2016, Adrian Rossiter

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
   Name: unitile2d.cc
   Description: uniform tilings in 2d
   Project: Antiprism - http://www.antiprism.com
*/

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>

#include "../base/antiprism.h"

using std::string;
using std::vector;

using namespace anti;

void to_nearest(double *num, double step) { *num = ceil(*num / step) * step; }

class unitile : public Geometry {
private:
  int pat; // pattern id number
  double x_end;
  double x_inc;
  double y_end;
  double y_inc;
  bool to_tile; // adjust size to make pattern tilable
  Trans3d trans_m;
  Vec3d shear;
  vector<Vec3d> poly; // current polygon

  void ut_4444();   //  1
  void ut_333333(); //  2
  void ut_666();    //  3
  void ut_3636();   //  4
  void ut_33344();  //  5
  void ut_33434();  //  6
  void ut_33336();  //  7
  void ut_31212();  //  8
  void ut_488();    //  9
  void ut_3464();   // 10
  void ut_4612();   // 11

  void set_polygon(int num_sides, double rot_ang = 0);
  void add_polygon(double x_start, double y_start);
  void set_incs(double xinc, double yinc);

public:
  enum { ut_open, ut_join, ut_join2, ut_twist, ut_twist2, ut_twist3 };

  unitile(int patt, double width, double height, bool totile = true,
          Trans3d trans = Trans3d(), Vec3d shr = Vec3d(0, 0, 0))
      : pat(patt), x_end(width), y_end(height), to_tile(totile), trans_m(trans),
        shear(shr)
  {
  }

  void plane(int lr_join, int tb_join);
  void conic_frust(double top_rad, double bot_rad, double ht);
  void mobius(double sect_rad, double ring_rad);
  void torus(double sect_rad, double ring_rad);
  void torus_trefoil(double sect_rad, double ring_rad, vector<double> pq);
  void klein(double sect_rad, double ring_rad);
  void klein2(double sect_rad, double ring_rad);
  void roman_boy(double t);
  void roman();
  void cross_cap();
  void cross_cap2();
  void torus4d(Trans4d rot4d_m);
  void klein4d(Trans4d rot4d_m);
  void proj4d(Trans4d rot4d_m);
};

void unitile::set_incs(double xinc, double yinc)
{
  x_inc = xinc;
  y_inc = yinc;
  if (to_tile) { // assumes starts are 0
    x_end = ceil(x_end / x_inc) * x_inc;
    y_end = ceil(y_end / y_inc) * y_inc;
  }
}

void unitile::set_polygon(int num_sides, double rot_ang)
{
  poly.clear();
  double ang = 2 * M_PI / num_sides;
  double rad = 0.5 / sin(ang / 2); // radius for unit edge polygon

  for (int i = 0; i < num_sides; i++)
    poly.push_back(
        Vec3d(rad * cos(i * ang + rot_ang), rad * sin(i * ang + rot_ang), 0));
}

void unitile::add_polygon(double x_start, double y_start)
{
  vector<Vec3d> &verts = raw_verts();
  vector<vector<int>> &faces = raw_faces();
  int x_steps = (x_end - epsilon - x_start) / x_inc;
  int y_steps = (y_end - epsilon - y_start) / y_inc;
  for (int i = 0; i <= x_steps; i++) {
    double x = x_start + i * x_inc;
    if (x < -epsilon)
      continue;
    for (int j = 0; j <= y_steps; j++) {
      double y = y_start + j * y_inc;
      if (y < -epsilon)
        continue;
      Vec3d cent(x, y, 0);
      vector<Vec3d> tile;
      for (auto &i : poly)
        tile.push_back(i + cent);
      vector<int> face;
      for (unsigned int j = 0; j < tile.size(); j++)
        face.push_back(verts.size() + j);
      faces.push_back(face);
      verts.insert(verts.end(), tile.begin(), tile.end());
    }
  }
}

void unitile::ut_4444()
{
  set_incs(1, 1);
  set_polygon(4, M_PI / 4);
  add_polygon(0.5, 0.5);
}

void unitile::ut_333333()
{
  set_incs(sqrt(3), 1);
  // set_incs(1, 1);
  set_polygon(3, M_PI);
  add_polygon(0, 0);
  add_polygon(sqrt(3) / 2, 0.5);
  set_polygon(3);
  add_polygon(sqrt(3) / 3, 0);
  add_polygon(sqrt(3) * 5 / 6, 0.5);
}

void unitile::ut_666()
{
  set_incs(sqrt(3), 3);
  set_polygon(6, M_PI / 6);
  add_polygon(0, 0);
  add_polygon(sqrt(3) / 2, 1.5);
}

void unitile::ut_3636()
{
  set_incs(sqrt(3) * 2, 2);
  set_polygon(3);
  add_polygon(sqrt(3) * 2 / 3, 0);
  add_polygon(sqrt(3) * 5 / 3, 1);
  set_polygon(3, M_PI);
  add_polygon(sqrt(3) * 1 / 3, 1);
  add_polygon(sqrt(3) * 4 / 3, 0);
  set_polygon(6, M_PI / 6);
  add_polygon(0, 0);
  add_polygon(sqrt(3), 1);
}

void unitile::ut_33344()
{
  set_incs(1, 2 + sqrt(3));
  set_polygon(3, M_PI / 2);
  add_polygon(0, 0.5 + sqrt(3) / 6);
  add_polygon(0.5, 1.5 + sqrt(3) * 2 / 3);
  set_polygon(3, -M_PI / 2);
  add_polygon(0, 1.5 + sqrt(3) * 5 / 6);
  add_polygon(0.5, 0.5 + sqrt(3) * 1 / 3);
  set_polygon(4, M_PI / 4);
  add_polygon(0, 0);
  add_polygon(0.5, 1 + sqrt(3) / 2);
}

void unitile::ut_33434()
{
  double l = 1 + sqrt(3);
  set_incs(l, l);
  set_polygon(4, -M_PI / 12);
  add_polygon(l / 4, l / 4);
  add_polygon(3 * l / 4, 3 * l / 4);
  set_polygon(4, M_PI / 12);
  add_polygon(l / 4, 3 * l / 4);
  add_polygon(3 * l / 4, l / 4);
  double d = 1 / sqrt(12);
  set_polygon(3, 0);
  add_polygon(d, l / 2);
  add_polygon(l / 2 + d, 0);
  set_polygon(3, M_PI);
  add_polygon(l - d, l / 2);
  add_polygon(l / 2 - d, 0);
  set_polygon(3, M_PI / 2);
  add_polygon(0, d);
  add_polygon(l / 2, l / 2 + d);
  set_polygon(3, -M_PI / 2);
  add_polygon(0, l - d);
  add_polygon(l / 2, l / 2 - d);
}

void unitile::ut_33336()
{
  double rot = acos(5 / (sqrt(7) * 2));
  set_incs(sqrt(7), sqrt(21));
  set_polygon(6, rot);
  add_polygon(0, 0);
  add_polygon(sqrt(7) / 2, sqrt(21) / 2);
  for (int i = 0; i < 6; i++) {
    double rot2 = rot + M_PI / 6 + i * M_PI / 3;
    Trans3d trans = Trans3d::rotate(0, 0, rot2) *
                    Trans3d::translate(Vec3d(sqrt(3) - 1 / sqrt(3), 0, 0));
    set_polygon(3, rot2);
    Vec3d cent = trans * Vec3d(0, 0, 0);
    add_polygon(cent[0], cent[1]);
    add_polygon(cent[0] + sqrt(7) / 2, cent[1] + sqrt(21) / 2);
    if (i < 3) {
      trans = trans * Trans3d::translate(Vec3d(0, 1, 0));
      cent = trans * Vec3d(0, 0, 0);
      if (i == 1)
        add_polygon(cent[0], cent[1]);
      add_polygon(cent[0] + sqrt(7) / 2, cent[1] + sqrt(21) / 2);
    }
  }
}

void unitile::ut_31212()
{
  set_incs(2 + sqrt(3), 3 + sqrt(3) * 2);
  set_polygon(12, M_PI / 12);
  add_polygon(0, 0);
  add_polygon(1 + sqrt(3) / 2, 1.5 + sqrt(3));
  set_polygon(3, -M_PI / 6);
  add_polygon(1 + sqrt(3) / 2, 2.5 + sqrt(3) * 5 / 3);
  add_polygon(0, 1 + sqrt(3) * 2 / 3);
  set_polygon(3, M_PI / 6);
  add_polygon(1 + sqrt(3) / 2, 0.5 + sqrt(3) / 3);
  add_polygon(0, 2 + sqrt(3) * 4 / 3);
}

void unitile::ut_488()
{
  double l = 1 + sqrt(2);
  set_incs(l, l);
  set_polygon(4);
  add_polygon(l / 2, l / 2);
  set_polygon(8, M_PI / 8);
  add_polygon(0, 0);
}

void unitile::ut_3464()
{
  set_incs(1 + sqrt(3), 3 + sqrt(3));
  set_polygon(3, M_PI / 6);
  add_polygon(0, 1 + sqrt(3) / 3);
  add_polygon(0.5 + sqrt(3) / 2, 2.5 + sqrt(3) * 5 / 6);
  set_polygon(3, -M_PI / 6);
  add_polygon(0, 2 + sqrt(3) * 2 / 3);
  add_polygon(0.5 + sqrt(3) / 2, 0.5 + sqrt(3) / 6);
  set_polygon(4, M_PI / 12);
  add_polygon(0.25 + sqrt(3) / 4, 0.75 + sqrt(3) / 4);
  add_polygon(-0.25 - sqrt(3) / 4, 2.25 + sqrt(3) * 3 / 4);
  set_polygon(4, -M_PI / 12);
  add_polygon(-0.25 - sqrt(3) / 4, 0.75 + sqrt(3) / 4);
  add_polygon(0.25 + sqrt(3) / 4, 2.25 + sqrt(3) * 3 / 4);
  set_polygon(4, M_PI / 4);
  add_polygon(0.5 + sqrt(3) / 2, 0);
  add_polygon(0, 1.5 + sqrt(3) / 2);
  set_polygon(6, M_PI / 6);
  add_polygon(0, 0);
  add_polygon(0.5 + sqrt(3) / 2, 1.5 + sqrt(3) / 2);
}

void unitile::ut_4612()
{
  set_incs(3 + sqrt(3), 3 + sqrt(3) * 3);
  set_polygon(4, M_PI / 12);
  add_polygon(0.75 + sqrt(3) / 4, 0.75 + sqrt(3) * 3 / 4);
  add_polygon(2.25 + sqrt(3) * 3 / 4, 2.25 + sqrt(3) * 9 / 4);
  set_polygon(4, -M_PI / 12);
  add_polygon(2.25 + sqrt(3) * 3 / 4, 0.75 + sqrt(3) * 3 / 4);
  add_polygon(0.75 + sqrt(3) / 4, 2.25 + sqrt(3) * 9 / 4);
  set_polygon(4, M_PI / 4);
  add_polygon(1.5 + sqrt(3) / 2, 0);
  add_polygon(0, 1.5 + sqrt(3) * 3 / 2);
  set_polygon(6);
  add_polygon(0, 1 + sqrt(3));
  add_polygon(0, 2 + sqrt(3) * 2);
  add_polygon(1.5 + sqrt(3) / 2, 2.5 + sqrt(3) * 5 / 2);
  add_polygon(1.5 + sqrt(3) / 2, 0.5 + sqrt(3) / 2);
  set_polygon(12, M_PI / 12);
  add_polygon(0, 0);
  add_polygon(1.5 + sqrt(3) / 2, 1.5 + sqrt(3) * 3 / 2);
}

void unitile::plane(int lr_join, int tb_join)
{
  double epsilon = 1e-6;
  clear_all();
  void (unitile::*pat_funcs[])() = {nullptr,
                                    &unitile::ut_4444,
                                    &unitile::ut_333333,
                                    &unitile::ut_666,
                                    &unitile::ut_3636,
                                    &unitile::ut_33344,
                                    &unitile::ut_33434,
                                    &unitile::ut_33336,
                                    &unitile::ut_31212,
                                    &unitile::ut_488,
                                    &unitile::ut_3464,
                                    &unitile::ut_4612};
  (this->*(pat_funcs[pat]))();

  double x_sh_inc = shear[0] * x_inc;
  double y_sh_inc = shear[1] * y_inc;
  double x_sh_inc2 = x_sh_inc * (1 + y_sh_inc / y_end);
  double y_sh_inc2 = y_sh_inc * (1 + x_sh_inc / x_end);
  double x_cross = x_end + x_sh_inc - x_sh_inc2;
  double y_cross = y_end + y_sh_inc - y_sh_inc2;
  vector<Vec3d> &verts = raw_verts();
  for (auto &vert : verts) {
    double x = vert[0];
    vert[0] += (shear[0] * x_inc) * (vert[1] / y_end);
    vert[0] *= x_end / x_cross;
    vert[1] += (shear[1] * y_inc) * (x / x_end);
    vert[1] *= y_end / y_cross;
    ;
  }
  for (auto &vert : verts) {
    if (lr_join == ut_twist &&
        (vert[0] < -epsilon || vert[0] > x_end - epsilon))
      vert[1] = -vert[1];
    if (lr_join == ut_twist2 &&
        (vert[0] < -epsilon || vert[0] > x_end - epsilon))
      vert[1] = y_end - vert[1];
    if (lr_join == ut_twist3 &&
        (vert[0] < -epsilon || vert[0] > x_end - epsilon))
      vert[1] = fmod(1.5 * y_end - vert[1], y_end);
    if (tb_join == ut_twist &&
        (vert[1] < -epsilon || vert[1] > y_end - epsilon))
      vert[0] = -vert[0];
    if (tb_join == ut_twist2 &&
        (vert[1] < -epsilon || vert[1] > y_end - epsilon))
      vert[0] = x_end - vert[0];
    if (tb_join == ut_twist3 &&
        (vert[1] < -epsilon || vert[1] > y_end - epsilon))
      vert[0] = fmod(1.5 * x_end - vert[0], x_end);
    if (tb_join == ut_join2 &&
        (vert[1] < -epsilon || vert[1] > y_end - epsilon))
      vert[0] = fmod(0.5 * x_end - vert[0], x_end);
    if (lr_join != ut_open) {
      vert[0] = fmod(vert[0] + x_end, x_end - epsilon);
    }
    if (tb_join != ut_open) {
      vert[1] = fmod(vert[1] + y_end, y_end - epsilon);
    }
  }
  merge_coincident_elements(*this, "v", 1e-5);

  transform(trans_m);
}

void unitile::conic_frust(double top_rad, double bot_rad, double ht)
{
  plane(ut_join, ut_open);
  double a0, h, rad;
  vector<Vec3d> &verts = raw_verts();
  for (auto &vert : verts) {
    a0 = 2 * M_PI * vert[0] / x_end;
    h = ht * vert[1] / y_end;
    rad = (top_rad - bot_rad) * vert[1] / y_end + bot_rad;
    vert = Vec3d(rad * cos(a0), rad * sin(a0), h);
  }
}

void unitile::mobius(double sect_rad, double ring_rad)
{
  plane(ut_twist2, ut_open);
  double a0;
  vector<Vec3d> &verts = raw_verts();
  for (auto &vert : verts) {
    a0 = 2 * M_PI * vert[0] / x_end;
    vert[1] -= y_end / 2;
    vert =
        Vec3d(sin(a0) * (sect_rad * vert[1] / y_end * cos(a0 / 2) + ring_rad),
              sect_rad * vert[1] / y_end * sin(a0 / 2),
              cos(a0) * (sect_rad * vert[1] / y_end * cos(a0 / 2) + ring_rad));
  }
}

void unitile::torus(double sect_rad, double ring_rad)
{
  plane(ut_join, ut_join);
  double a0, a1;
  vector<Vec3d> &verts = raw_verts();
  for (auto &vert : verts) {
    a0 = 2 * M_PI * vert[0] / x_end;
    a1 = 2 * M_PI * vert[1] / y_end;
    vert = Vec3d(sin(a1) * (sect_rad * cos(a0) + ring_rad), sect_rad * sin(a0),
                 cos(a1) * (sect_rad * cos(a0) + ring_rad));
  }
}

void unitile::torus_trefoil(double sect_rad, double ring_rad,
    vector<double> pq)
{
  plane(ut_join, ut_join);
  vector<Vec3d> &verts = raw_verts();
  int i=-1;
  for (auto &vert : verts) {
    i++;
    double a0 = 2 * M_PI * vert[0] / x_end;
    double a1 = 2 * M_PI * vert[1] / y_end;

    int q = 3;
    int p = 2;

    if(pq.size()>0) {
      if(floor(pq[0]) > 0)
         q = floor(pq[0]);
      p = 1;
      if(pq.size()>1 && floor(pq[1]) > 0)
         p = floor(pq[1]);
    }

    double a = ring_rad;
    double b = 1.0;
    double x = (a + b*cos(q*a1))*cos(p*a1);
    double y = (a + b*cos(q*a1))*sin(p*a1);
    double z = b*sin(q*a1);
    Vec3d P(x, y, z);

    double a11 = a1 + 0.0001;
    double x0 = (a + b*cos(q*a11))*cos(p*a11);
    double y0 = (a + b*cos(q*a11))*sin(p*a11);
    double z0 = b*sin(q*a11);
    Vec3d P0(x0, y0, z0);
    Vec3d Q(a*cos(p*a1), a*sin(p*a1), 0);

    Vec3d dir = (P0 - P).unit();
    Vec3d norm = P - Q;

    double ang = a0+(p+1)*a1; // stops the tube twist
    Vec3d v(sect_rad * cos(ang), sect_rad * sin(ang), 0);
    vert = Trans3d::translate(P) *
           Trans3d::align(Vec3d::Z, Vec3d::Y, dir, norm) * v;
  }
}

// http://en.wikipedia.org/wiki/Image:KleinBottle-01.png
void unitile::klein(double sect_rad, double /*ring_rad*/)
{
  plane(ut_twist3, ut_join);
  double a0, a1;
  vector<Vec3d> &verts = raw_verts();
  for (auto &vert : verts) {
    // a0 = 2*M_PI*fmod(1+verts[i][0]/x_end, 1);
    // a1 = 2*M_PI*fmod(1+verts[i][1]/y_end, 1);
    a0 = 2 * M_PI * vert[0] / x_end;
    a1 = 2 * M_PI * vert[1] / y_end;
    if (a0 < M_PI)
      vert = Vec3d(6 * cos(a0) * (1 + sin(a0)) +
                       4 * sect_rad * (1 - 0.5 * cos(a0)) * cos(a0) * cos(a1),
                   16 * sin(a0) +
                       4 * sect_rad * (1 - 0.5 * cos(a0)) * sin(a0) * cos(a1),
                   4 * sect_rad * (1 - 0.5 * cos(a0)) * sin(a1));
    else
      vert = Vec3d(6 * cos(a0) * (1 + sin(a0)) -
                       4 * sect_rad * (1 - 0.5 * cos(a0)) * cos(a1),
                   16 * sin(a0), 4 * sect_rad * (1 - 0.5 * cos(a0)) * sin(a1));
  }
}

void unitile::klein2(double sect_rad, double ring_rad)
{
  plane(ut_twist, ut_join);
  double a0, a1;
  vector<Vec3d> &verts = raw_verts();
  for (auto &vert : verts) {
    a0 = 2 * M_PI * vert[0] / x_end;
    a1 = 2 * M_PI * vert[1] / y_end;
    vert = Vec3d(
        (ring_rad +
         (cos(0.5 * a0) * sin(a1) - sin(0.5 * a0) * sin(2 * a1)) * sect_rad) *
            cos(a0),
        (ring_rad +
         (cos(0.5 * a0) * sin(a1) - sin(0.5 * a0) * sin(2 * a1)) * sect_rad) *
            sin(a0),
        (sin(0.5 * a0) * sin(a1) + cos(0.5 * a0) * sin(2 * a1)) * sect_rad);
  }
}

void unitile::roman()
{
  plane(ut_open, ut_open);
  double a0, a1;
  vector<Vec3d> &verts = raw_verts();
  for (auto &vert : verts) {
    a0 = M_PI * vert[0] / x_end;
    a1 = M_PI * (vert[1] / y_end - 0.5);
    vert = Vec3d(0.5 * sin(2 * a0) * sin(a1) * sin(a1),
                 0.5 * sin(a0) * cos(2 * a1), 0.5 * cos(a0) * sin(2 * a1));
    // verts[i] = Vec3d(
    //   sin(2*a0)*cos(a1)*cos(a1),
    //   sin(a0)*sin(2*a1),
    //   cos(a0)*cos(2*a1) );
    vert = Vec3d(0.5 * cos(a0) * sin(2 * a1), 0.5 * sin(a0) * sin(2 * a1),
                 0.5 * sin(2 * a0) * cos(a1) * cos(a1));
  }
}

void unitile::roman_boy(double t)
{
  plane(ut_open, ut_open);
  double a0, a1;
  vector<Vec3d> &verts = raw_verts();
  for (auto &vert : verts) {
    a0 = M_PI * (0.5 - vert[0] / x_end);
    a1 = M_PI * vert[1] / y_end;
    double x =
        (sqrt(2) * cos(2 * a0) * cos(a1) * cos(a1) + cos(a0) * sin(2 * a1)) /
        (2 - t * sqrt(2) * sin(3 * a0) * sin(2 * a1));
    double y =
        (sqrt(2) * sin(2 * a0) * cos(a1) * cos(a1) - sin(a0) * sin(2 * a1)) /
        (2 - t * sqrt(2) * sin(3 * a0) * sin(2 * a1));
    double z =
        (3 * cos(a1) * cos(a1)) / (2 - t * sqrt(2) * sin(3 * a0) * sin(2 * a1));
    vert = Vec3d(x, y, z);
  }
}

void unitile::cross_cap()
{
  plane(ut_twist2, ut_twist2);
  vector<Vec3d> &verts = raw_verts();
  for (auto &vert : verts) {
    double x = 1 - 2.0 * vert[0] / x_end;
    double y = 1 - 2.0 * vert[1] / y_end;
    double a0 = atan2(y, x);
    // x *= 0.9;
    // y *= 0.9;
    double dist_to_edge;
    if (fabs(x) < epsilon && fabs(y) < epsilon)
      dist_to_edge = 1;
    else if (fabs(x) > fabs(y))
      dist_to_edge = sqrt(1 + fabs(y * y / (x * x)));
    else
      dist_to_edge = sqrt(1 + fabs(x * x / (y * y)));

    x /= dist_to_edge;
    y /= dist_to_edge;
    double r = 2.0 / 3 + (1 / (1 + cos(2 * (a0 - M_PI / 4)) / 2));
    double a1 = M_PI * sqrt(x * x + y * y);

    x = r * cos(a0) * sin(a1);
    y = r * sin(a0) * sin(a1);
    double z = r * (cos(a1) - 1);
    vert = Vec3d(x, y, z);
  }
}

void unitile::cross_cap2()
{
  plane(ut_open, ut_open);
  double a0, a1;
  vector<Vec3d> &verts = raw_verts();
  for (auto &vert : verts) {
    a0 = 2 * M_PI * vert[0] / x_end;
    a1 = 0.5 * M_PI * vert[1] / y_end;
    // verts[i] = Vec3d(cos(a0)*sin(2*a1), sin(a0)*sin(2*a1),
    //                 cos(a1)*cos(a1)-cos(a0)*cos(a0)*sin(a1)*sin(a1));

    vert = Vec3d(sin(a0) * sin(2 * a1), sin(2 * a0) * sin(a1) * sin(a1),
                 cos(2 * a0) * sin(a1) * sin(a1));
  }
}

void unitile::torus4d(Trans4d rot4d_m)
{
  plane(ut_join, ut_join);
  double a0, a1;
  vector<Vec3d> &verts = raw_verts();
  for (auto &vert : verts) {
    a0 = 2 * M_PI * vert[0] / x_end;
    a1 = 2 * M_PI * vert[1] / y_end;
    Vec4d v4d = rot4d_m * Vec4d(cos(a0), sin(a0), cos(a1), sin(a1));
    // Vec4d v4d = rot * Vec4d(cos(a0+a1), sin(a0+a1), cos(a0-a1), sin(a0-a1));
    vert = Vec3d(v4d[0], v4d[1], v4d[2]);
  }
}

void unitile::klein4d(Trans4d rot4d_m)
{
  plane(ut_join, ut_join2);
  double a0, a1;
  vector<Vec3d> &verts = raw_verts();
  for (auto &vert : verts) {
    a0 = 2 * M_PI * vert[0] / x_end;
    a1 = 2 * M_PI * vert[1] / y_end;

    double x = cos(a0) / (1 + pow(sin(a0), 2));
    double y = sin(a0) * x;
    double a2 = a1 / 2;
    double x2 = x * cos(a2) - y * sin(a2);
    double y2 = x * sin(a2) + y * cos(a2);
    Vec4d v4d = rot4d_m * Vec4d(x2, y2, cos(a1), sin(a1));
    // Vec4d v4d = rot * Vec4d(cos(a0+a1), sin(a0+a1), cos(a0-a1), sin(a0-a1));
    vert = Vec3d(v4d[0], v4d[1], v4d[2]);
  }
}

void unitile::proj4d(Trans4d rot4d_m)
{
  plane(ut_twist, ut_twist);
  double a0, a1;
  vector<Vec3d> &verts = raw_verts();
  for (auto &vert : verts) {
    a0 = -M_PI * (2 * vert[0] / x_end - 0.5);
    a1 = M_PI * (2 * vert[1] / y_end - 0.5);

    double x = cos(a0) / (1 + pow(sin(a0), 2));
    double y = sin(a0) * x;
    double a2 = a1 / 2;
    double x2 = x * cos(a2) - y * sin(a2);
    double y2 = x * sin(a2) + y * cos(a2);
    double z = cos(a1) / (1 + pow(sin(a1), 2));
    double w = sin(a1) * z;
    double a3 = a0 / 2;
    double z2 = z * cos(a3) - w * sin(a3);
    double w2 = z * sin(a3) + w * cos(a3);
    Vec4d v4d = rot4d_m * Vec4d(x2, y2, z2, w2);
    // Vec4d v4d = rot * Vec4d(cos(a0+a1), sin(a0+a1), cos(a0-a1), sin(a0-a1));
    vert = Vec3d(v4d[0], v4d[1], v4d[2]);
  }
}

class ut_opts : public ProgramOpts {
public:
  char surface;
  int pattern;
  double width;
  double height;
  double to_tile;
  double r;
  double R;
  vector <double> d;
  Trans3d trans_m;
  Trans4d rot4d_m;
  Vec3d shear;

  string ofile;

  ut_opts()
      : ProgramOpts("unitile2d"), surface('p'), pattern(1), width(20),
        height(-1), to_tile(true), r(1), R(3), shear(Vec3d(0, 0, 0))
  {
  }
  void process_command_line(int argc, char **argv);
  void usage();
};

// clang-format off
void ut_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [pattern]\n"
"\n"
"Makes uniform tilings on a surface. Surfaces include plane or cylinder\n"
"section, torus, Klein bottle, and many others.\n"
"Patterns are (default: 1)\n"
"   1 - 4,4,4,4      {4,4}    5 - 3,3,3,4,4      9  - 4,8,8\n"
"   2 - 3,3,3,3,3,3  {3,6}    6 - 3,3,4,3,4      10 - 3,4,6,4\n"
"   3 - 6,6,6        {6,3}    7 - 3,3,3,3,6      11 - 4,6,12\n"
"   4 - 3,6,3,6               8 - 3,12,12\n"
"\n"
"Options\n"
"%s"
"  -s <type> surface type:\n"
"        p - plane (default)\n"
"        c - conical frustum (-R bottom radius, -r top radius, -d height,\n"
"            -T X,Y,0)\n"
"        m - mobius strip (-R ring radius, -r strip width, -T X,0,0)\n"
"        t - torus (-R ring radius, -r tube radius. -T X,Y,0)\n"
"        T - trefoil (-R ring radius, -r tube radius. -T X,Y,0, -d N,D)\n"
"        k - Klein bottle (-r tube radius, -T X,0,0)\n"
"        K - figure-8 Klein bottle (-R ring radius, -r tube radius, -T X,0,0)\n"
"        C - cross-cap (needs tile vertices at tiling rectangle corners)\n"
"        w - cross-cap (does not preserve tiling, only works with pattern 1!)\n"
"        r - Roman (does not preserve tiling)\n"
"        b - Boy's (does not presere tiling)\n"
"        R - Roman to Boy's (-d stage of tansformation, range 0.0-1.0\n"
"            does not preserve tiling\n"
"        x - 4D torus (-W rotatation in 4D before projection onto xyz, -T X,Y,0\n"
"        y - 4D Klein bottle (-W rotatation in 4D before projection onto xyz)\n"
"  -w <wdth> width of tiling (default: 20)\n"
"  -l <ht>   height of tiling (default: width)\n"
"  -g        use given height and width, don't increase to make\n"
"            rectangular tile\n"
"  -r <rad>  'minor' radius of surface\n"
"  -R <rad>  'major' radius of surface\n"
"  -d <dist> height of conic frustrum (-s C, default: 5.0) or parameter for\n"
"            Roman to Boy's surface (-s R, default 0.0)\n"
"  -T <tran> translate pattern, three numbers separated by commas which are\n"
"            used as the x, y and z displacements\n"
"  -S <X,Y>  \"shear\" the base rectangular tiling pattern by X units in the\n"
"            w direction and Y units in the l direction (used with -s t)\n"
"  -W <rot>  rotation of 4D surface before projection, six angles\n"
"            separated by commas to rotate in planes xy,yz,zw,wx,xz,yw\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}
// clang-format on

void ut_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;
  vector<double> nums;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hp:s:l:w:r:R:d:gT:S:W:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 's':
      if (strlen(optarg) != 1 || !strchr("pcmtTkKCrbRwxy", *optarg))
        error("surface type must be one of pcmtTkKCwrbRxy", c);
      surface = *optarg;
      break;

    case 'w':
      print_status_or_exit(read_double(optarg, &width), c);
      if (width < 0)
        error("width cannot be negative", c);
      break;

    case 'l':
      print_status_or_exit(read_double(optarg, &height), c);
      if (height < 0)
        error("length cannot be negative", c);
      break;

    case 'g':
      to_tile = false;
      break;

    case 'r':
      print_status_or_exit(read_double(optarg, &r), c);
      if (r < 0)
        warning("radius is negative", c);
      break;

    case 'R':
      print_status_or_exit(read_double(optarg, &R), c);
      if (R < 0)
        warning("radius is negative", c);
      break;

    case 'd':
      print_status_or_exit(read_double_list(optarg, d), c);
      break;

    case 'T':
      print_status_or_exit(read_double_list(optarg, nums), c);
      if (nums.size() != 3)
        error(msg_str("must give exactly three numbers (%lu were given)",
                      (unsigned long)nums.size()),
              c);
      trans_m = Trans3d::translate(Vec3d(nums[0], nums[1], nums[2]));
      break;

    case 'S':
      print_status_or_exit(read_double_list(optarg, nums), c);
      if (nums.size() != 2)
        error(msg_str("must give exactly thwo numbers (%lu were given)",
                      (unsigned long)nums.size()),
              c);
      // trans_m = Trans3d::transl(Vec3d(nums[0], nums[1], nums[2]));
      shear[0] = nums[0];
      shear[1] = nums[1];
      break;

    case 'W':
      print_status_or_exit(read_double_list(optarg, nums), c);
      if (nums.size() != 6)
        error(msg_str("4d rotation must be exactly six angles "
                      "(%lu were given)",
                      (unsigned long)nums.size()),
              c);
      for (auto &num : nums)
        num = deg2rad(num);
      rot4d_m =
          Trans4d::rotate(nums[0], nums[1], nums[2], nums[3], nums[4], nums[5]);
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (argc - optind > 1)
    error("too many arguments");

  pattern = 1;
  if (argc - optind > 0) {
    print_status_or_exit(read_int(argv[optind], &pattern), "pattern");
    if (pattern < 1 || pattern > 11)
      error("must be a number from 1 to 11", "pattern");
  }

  if (height < 0)
    height = width;
}

int main(int argc, char *argv[])
{
  ut_opts opts;
  opts.process_command_line(argc, argv);

  unitile ut(opts.pattern, opts.width, opts.height, opts.to_tile, opts.trans_m,
             opts.shear);

  switch (opts.surface) {
  case 'p':
    ut.plane(unitile::ut_open, unitile::ut_open);
    break;
  case 'c':
    ut.conic_frust(opts.r, opts.R, opts.d.size()==0 ? 5.0 : opts.d[0]);
    break;
  case 'm':
    ut.mobius(opts.r, opts.R);
    break;
  case 't':
    ut.torus(opts.r, opts.R);
    break;
  case 'T':
    ut.torus_trefoil(opts.r, opts.R, opts.d);
    break;
  case 'k':
    ut.klein(opts.r, opts.R);
    break;
  case 'K':
    ut.klein2(opts.r, opts.R);
    break;
  case 'C':
    ut.cross_cap();
    break;
  case 'r':
    // ut.roman_boy(0);
    ut.roman();
    break;
  case 'b':
    ut.roman_boy(1);
    break;
  case 'R':
    ut.roman_boy(opts.d.size()==0 ? 0.0 : opts.d[0]);
    break;
  case 'w':
    ut.cross_cap2();
    break;
  case 'x':
    ut.torus4d(opts.rot4d_m);
    break;
  case 'y':
    ut.klein4d(opts.rot4d_m);
    break;
  case 'z':
    ut.proj4d(opts.rot4d_m);
    break;
  }

  opts.write_or_error(ut, opts.ofile);

  return 0;
}
