/*
   Copyright (c) 2003-2022, Adrian Rossiter

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
   Name: lat_grid_common.cc
   Description: lat_grid code shared in /src
   Project: Antiprism - http://www.antiprism.com
*/

#include "lat_grid_common.h"
#include "../base/antiprism.h"

#include <cstdio>
#include <string>
#include <vector>

using std::string;
using std::vector;

using namespace anti;

// for lat_grid.cc

bool sc_test(int /*x*/, int /*y*/, int /*z*/) // dist2 = 1, 2, 3
{
  return 1;
}

bool fcc_test(int x, int y, int z) // dist2 = 2, 4, 6, 12
{
  return (x + y + z) % 2 == 0;
}

bool bcc_test(int x, int y, int z) // dist2 = 3, 4, 8
{
  return ((x % 2 && y % 2 && z % 2) || !(x % 2 || y % 2 || z % 2));
}

bool hcp_test(int x, int y, int z) // dist2 = 18
{
  int m = x + y + z;
  int n = (m < 0); // for integer division
  return (m % 6 == 0) && ((x - m / 12 + n) % 3 == 0) &&
         ((y - m / 12 + n) % 3 == 0);
}

bool rh_dodec_test(int x, int y, int z) // dist2 = 3 (8)
{
  return ((x % 2 && y % 2 && z % 2) ||
          !(((x + y + z) / 2) % 2 == 0 || (x % 2 || y % 2 || z % 2)));
}

bool cubo_oct_test(int x, int y, int z) // dist2 = 2
{
  return (std::abs(x) % 2 + std::abs(y) % 2 + std::abs(z) % 2) == 2;
}

bool tr_oct_test(int x, int y, int z) // dist2 = 2
{
  return ((z % 4 == 0 && x % 4 && y % 4 && (x - y) % 2) ||
          (y % 4 == 0 && z % 4 && x % 4 && (z - x) % 2) ||
          (x % 4 == 0 && y % 4 && z % 4 && (y - z) % 2));
}

bool tr_tet_tet_test(int x, int y, int z) // dist2 = 2
{
  return ((x % 2 == 0 && (y % 4 + 4) % 4 == (((x + z)) % 4 + 4) % 4) ||
          (y % 2 == 0 && (z % 4 + 4) % 4 == (((y + x)) % 4 + 4) % 4) ||
          (z % 2 == 0 && (x % 4 + 4) % 4 == (((z + y)) % 4 + 4) % 4));
}

bool tr_tet_tr_oct_cubo_test(int x, int y, int z) // dist2 = 4
{
  x = std::abs(x) % 6;
  if (x > 3)
    x = 6 - x;
  y = std::abs(y) % 6;
  if (y > 3)
    y = 6 - y;
  z = std::abs(z) % 6;
  if (z > 3)
    z = 6 - z;
  int dist2 = x * x + y * y;
  return ((z % 6 == 0 && (dist2 == 2 || dist2 == 8)) ||
          (z % 6 == 1 && (dist2 == 1 || dist2 == 13)) ||
          (z % 6 == 2 && (dist2 == 4 || dist2 == 10)) ||
          (z % 6 == 3 && dist2 == 5));
}

bool diamond_test(int x, int y, int z) //  dist2 = 3
{
  return (((x % 2 + 2) % 2 + (y % 2 + 2) % 2 + (z % 2 + 2) % 2) % 3 == 0 &&
          (x / 2 + (x % 2 < 0) + y / 2 + (y % 2 < 0) + z / 2 + (z % 2 < 0)) %
                  2 ==
              0);
}

// Coordinates from Wendy Krieger
bool hcp_diamond_test(int x, int y, int z) //  dist2 = 27
{
  int pt[][3] = {{0, 0, 0}, {3, 3, 3}, {6, 0, 6}, {9, 3, 9}};
  for (auto &i : pt) {
    int tx = x - i[0];
    int ty = y - i[1];
    int tz = z - i[2];
    int sum = tx + ty + tz;
    if (sum % 24 == 0) {
      int n8 = sum / 3;
      if ((tx - n8) % 6 == 0 && (ty - n8) % 6 == 0 && (tz - n8) % 6 == 0)
        return true;
    }
  }
  return false;
}

// Coordinates from Vladimir Bulatov
bool k_4_test(int x, int y, int z) //  dist2 = 2
{
  if ((x + y + z) % 2)
    return false;
  x = (x % 4) + ((x < 0) ? 4 : 0);
  y = (y % 4) + ((y < 0) ? 4 : 0);
  z = (z % 4) + ((z < 0) ? 4 : 0);
  if ((x == 0 && y == 0 && z == 0) || (x == 0 && y == 1 && z == 3) ||
      (x == 1 && y == 0 && z == 1) || (x == 1 && y == 1 && z == 2) ||
      (x == 2 && y == 2 && z == 2) || (x == 2 && y == 3 && z == 1) ||
      (x == 3 && y == 2 && z == 3) || (x == 3 && y == 3 && z == 0))
    return true;
  return false;
}

void add_struts(Geometry &geom, int len2)
{
  const vector<Vec3d> &verts = geom.verts();
  for (unsigned int i = 0; i < verts.size(); i++)
    for (unsigned int j = i; j < verts.size(); j++) {
      if (fabs((verts[i] - verts[j]).len2() - len2) < anti::epsilon)
        geom.add_edge(make_edge(i, j));
    }
}

void int_lat_grid::make_lattice(Geometry &geom)
{
  if (!centre.is_set())
    centre = Vec3d(1, 1, 1) * (o_width / 2.0);
  double o_off = o_width / 2.0 + anti::epsilon;
  double i_off = i_width / 2.0 - anti::epsilon;
  int i, j, k;
  for (k = int(ceil(centre[2] - o_off)); k <= centre[2] + o_off; k++)
    for (j = int(ceil(centre[1] - o_off)); j <= centre[1] + o_off; j++)
      for (i = int(ceil(centre[0] - o_off)); i <= centre[0] + o_off; i++) {
        if (i > centre[0] - i_off && i < centre[0] + i_off &&
            j > centre[1] - i_off && j < centre[1] + i_off &&
            k > centre[2] - i_off && k < centre[2] + i_off)
          continue;
        if (coord_test(i, j, k))
          geom.add_vert(Vec3d(i, j, k));
      }
}

void sph_lat_grid::make_lattice(Geometry &geom)
{
  if (!centre.is_set())
    centre = Vec3d(0, 0, 0);
  double o_off = o_width + anti::epsilon;
  double i_off = i_width - anti::epsilon;
  int i, j, k;
  for (k = int(ceil(centre[2] - o_off)); k <= centre[2] + o_off; k++)
    for (j = int(ceil(centre[1] - o_off)); j <= centre[1] + o_off; j++)
      for (i = int(ceil(centre[0] - o_off)); i <= centre[0] + o_off; i++) {
        double dist2 = (Vec3d(i, j, k) - centre).len2();
        if (o_off < dist2 || i_off > dist2)
          continue;
        if (coord_test(i, j, k))
          geom.add_vert(Vec3d(i, j, k));
      }
}
