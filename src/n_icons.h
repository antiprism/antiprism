/*
   Copyright (c) 2007-2020, Roger Kaufman

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
   Name: n_icons.h
   Description: Creates Sphericon like Polyhedra
   Project: Antiprism - http://www.antiprism.com
*/

#ifndef NCONS_H
#define NCONS_H

class coordList {
public:
  int coord_no;
  bool rotated;
  coordList(int c) : coord_no(c) { rotated = false; }
};

// pole 0 - North Pole
// pole 1 - South Pole
class poleList {
public:
  int idx;
  int lat;
  poleList()
  {
    idx = -1;
    lat = -1;
  }
};

class edgeList {
public:
  int edge_no;
  int lat;
  int lon;
  bool rotate;
  edgeList(int f, int lat, int lon) : edge_no(f), lat(lat), lon(lon)
  {
    rotate = false;
  }
};

class faceList {
public:
  int face_no;
  int lat;
  int lon;
  int polygon_no;
  bool rotate;
  faceList(int f, int lat, int lon, int polygon_no)
      : face_no(f), lat(lat), lon(lon), polygon_no(polygon_no)
  {
    rotate = false;
  }
};

class polarOrb {
public:
  int coord_no;
  int forward;
  int backward;
  polarOrb(int c) : coord_no(c) {}
};

struct surfaceData {
  int c_surfaces;
  int c_edges;
  int d_surfaces;
  int d_edges;
  int total_surfaces;
  int total_edges;
  bool nonchiral;
  int case1_twist;
  bool case2;
  int compound_parts; // borrow for compounds
};

#endif
