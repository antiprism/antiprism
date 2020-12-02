/*
   Copyright (c) 2009-2020, Roger Kaufman

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
   Name: skilling.cc
   Description: Uniform Compounds catalogued by John Skilling
   Project: Antiprism - http://www.antiprism.com
*/

#include "coloring.h"
#include "geometryutils.h"
#include "mathutils.h"
#include "private_std_polys.h"
#include "random.h"
#include "utils.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using std::string;
using std::vector;

using namespace anti;

// clang-format off
// oct_case = 0 // UC19 20 Tetrahemihexahedra 45.6696674755810220... degrees (using 0 degrees as starting point)
// oct_case = 1 // UC14 20 Octahedra -23.430911382616046... degrees (using UC10 which starts at -45 degrees) (negate the returned angle)
//                 case 1 not currently used. Instead use oct_case = 2 and add 45.6696674755810220... = 67.908423568545970 degrees
// oct_case = 2 // UC15 10 Octahedra  22.238756092964941... degrees (using UC10 white starts at -45 degrees)
//              // UC16 10 Octahedra  82.238756092964941... degrees (use oct_case = 2 and add 60 degrees)
//                 note: information only. UC15 and UC16 are generated with 0 and 60 degrees using the same tranformation used by 20 Tetrahemihexahedra
// returns radians
// clang-format on
static double octahedra_angle(int oct_case)
{
  double a = 1 / (phi * phi) - 1 / sqrt(phi);
  double b = -1 + 1 / sqrt(phi * phi * phi);
  double c = 1 / phi + sqrt(phi);

  Vec3d v_ax = Vec3d(1, 1, 1); // unit axis vector

  // start out with 45.6696674755810220... = 23.430911382616046... +
  // 22.238756092964941...
  Vec3d v0 = Vec3d(a, -b, c);
  Vec3d v1 = Vec3d(-a, b, c);

  if (oct_case == 1)
    v1 = Vec3d(0, 0, 1); // this angle is 23.430911382616046...
  else if (oct_case == 2)
    v0 = Vec3d(0, 0, 1); // this angle is 22.238756092964941...

  return angle_around_axis(v0, v1, v_ax);
}

static void build_uniform_compound(Geometry &geom, int uc_case, int uc_num,
                                   string sym_from, string sym_to, double angle)
{
  // UC04 2 tetrahedra
  // UC05 5 tetrahedra
  // UC06 10 tetrahedra
  // UC09 5 cubes
  // UC17 5 octahedra
  // UC18 5 tetrahemihexahedra
  // UC42 3 square antiprisms
  // UC43 6 square antiprisms
  // UC46 2 icosahedra
  // UC48 2 great dodecahedra
  // UC50 2 small stellated dodecahedra
  // UC52 2 great icosahedra
  // UC54 2 truncated tetrahedra
  // UC55 5 truncated tetrahedra
  // UC56 10 truncated tetrahedra
  // UC57 5 truncated cubes
  // UC58 5 stellated truncated hexahedra
  // UC59 5 cuboctahedra
  // UC60 5 cubohemioctahedra
  // UC61 5 octahemioctahedra
  // UC62 5 rhombicuboctahedra
  // UC63 5 small rhombihexahedra
  // UC64 5 small cubicuboctahedra
  // UC65 5 great cubicuboctahedra
  // UC66 5 great rhombihexahedra
  // UC67 5 great rhombicuboctahedra
  // UC68 2 snub cubes
  // UC69 2 snub dodecahedra
  // UC70 2 great snub icosidodecahedra
  // UC71 2 great inverted snub icosidodecahedra
  // UC72 2 great retrosnub icosidodecahedra
  // UC73 2 snub dodecadodecahedra
  // UC74 2 inverted snub dodecadodecahedra
  // UC75 2 snub icosidodecadodecahedra
  if (uc_case == 0) {
    if (uc_num == 6 || uc_num == 56) {
      transform_and_repeat(geom, "Oh", sym_from);
      sym_from = "Oh";
    }
    else if (uc_num == 43) {
      transform_and_repeat(geom, "D8", sym_from);
      sym_from = "D8";
    }
    transform_and_repeat(geom, sym_to, sym_from);
  }
  else if (uc_case == 1) {
    // UC01 6 tetrahedra rotational
    // UC02 12 tetrahedra rotational
    // UC03 6 tetrahedra
    // RK - for UC02 & UC03, generate this model actually as 6 stella
    // octangula
    if (uc_num == 2 || uc_num == 3) {
      transform_and_repeat(geom, "S2", sym_from);
      if (uc_num == 3)
        sym_from = "Oh";
    }
    transform_and_repeat(geom, sym_to, sym_from, Trans3d::rotate(0, 0, angle));
  }
  else if (uc_case == 2) {
    // UC07 6 cubes rotational
    geom.transform(Trans3d::rotate(0, 0, angle));
    transform_and_repeat(geom, "D2", sym_from);
    transform_and_repeat(geom, sym_to, "D2");
  }
  else if (uc_case == 3) {
    // UC10 4 octahedra rotational
    // UC11 8 octahedra rotational
    // UC12 4 octahedra
    // UC13 20 octahedra rotational
    // UC14 20 octahedra
    // UC30 4 triangular prisms
    // UC31 8 triangular prisms
    // UC38 4 hexagonal prisms
    string sym_to_local = sym_to;
    if (uc_num == 10 || uc_num == 11 || uc_num == 13 || uc_num == 14) {
      angle -= M_PI / 4; // variable (subtract 45 degrees so 0 causes
                         // simultaneous octahedra)
      // make UC10 first, then use that to make UC11, UC13 and UC14
      if (uc_num == 11 || uc_num == 13 || uc_num == 14)
        sym_to_local = "T";
    }
    else if (uc_num == 31)
      transform_and_repeat(geom, "D3v", "D3h");
    transform_and_repeat(geom, sym_to_local, sym_from,
                         Trans3d::rotate(Vec3d(0, 0, 1), Vec3d(1, 1, 1)) *
                             Trans3d::rotate(0, 0, angle));
    if (uc_num == 11 || uc_num == 13 || uc_num == 14)
      transform_and_repeat(geom, sym_to, "T");
  }
  else if (uc_case == 4) {
    // UC15 10 octahedra 1 (0 degrees)
    // UC16 10 octahedra 2 (60 degrees)
    // UC19 20 tetrahemihexahedra
    if (uc_num == 19)
      geom.transform(
          Trans3d::rotate(Vec3d(0, 0, 1),
                          Vec3d(0.5, 1 / sqrt(12), 1 / sqrt(6))) *
          Trans3d::rotate(0, 0, -M_PI / 12)); // -15 degrees, placement of
                                              // tetrahemihexahedra like a
                                              // 3-antiprism
    transform_and_repeat(geom, sym_to, sym_from,
                         Trans3d::align(Vec3d(0, 0, 1), Vec3d(1, 0, 0),
                                        Vec3d(1, 1, 1),
                                        Vec3d(-phi, -1, phi + 1)) *
                             Trans3d::rotate(0, 0, angle));
  }
  else if (uc_case == 5) {
    // UC20 2k n d gonal prisms rotational
    // UC21 k n d gonal prisms
    // UC22 2k n odd d gonal antiprisms rotational
    // UC23 k n odd d gonal antiprisms
    // UC24 2k n even d gonal antiprisms rotational
    // UC25 k n even d gonal antiprisms
    if (uc_num == 20 || uc_num == 22 || uc_num == 24)
      transform_and_repeat(geom, sym_from, sym_from,
                           Trans3d::rotate(0, 0, angle)); // variable
    transform_and_repeat(geom, sym_to, sym_from);
  }
  else if (uc_case == 6) {
    // UC26 12 pentagonal antiprisms rotational
    // UC28 12 pentagrammic crossed antiprisms rotational
    angle += (uc_num == 26) ? M_PI / 5 : 0; // UC26 needs 36 degrees to
                                            // correspond to the 0 degree point
                                            // of UC28
    geom.transform(Trans3d::rotate(0, 0, angle));
    transform_and_repeat(geom, "D10h", sym_from);
    transform_and_repeat(geom, sym_to, "D10h",
                         Trans3d::rotate(Vec3d(0, 0, 1), Vec3d(0, 1, phi)));
  }
  else if (uc_case == 7) {
    // UC27 6 pentagonal antiprisms (angle = 0)
    // UC29 6 pentagrammic crossed antiprisms
    // UC34 6 pentagonal prisms
    // UC35 12 pentagonal prisms
    // UC36 6 pentagrammic prisms
    // UC37 12 pentagrammic prisms
    // UC40 6 decagonal prisms
    // UC41 6 decagrammic prisms
    // UC44 6 pentagrammic antiprisms
    // UC45 12 pentagrammic antiprisms
    if (uc_num == 35 || uc_num == 37 || uc_num == 45) {
      transform_and_repeat(geom, "D10h", sym_from);
      sym_from = "D10h";
    }
    transform_and_repeat(geom, sym_to, sym_from,
                         Trans3d::rotate(Vec3d(0, 0, 1), Vec3d(0, 1, phi)) *
                             Trans3d::rotate(0, 0, angle));
  }
  else if (uc_case == 8) {
    // UC32 10 triangular prisms
    // UC33 20 triangular prisms
    // UC39 10 hexagonal prisms
    if (uc_num == 33) {
      transform_and_repeat(geom, "D6h", sym_from);
      sym_from = "D6h";
    }
    transform_and_repeat(
        geom, sym_to, sym_from,
        Trans3d::rotate(Vec3d(0, 0, 1), Vec3d(1 / phi, 0, phi)) *
            Trans3d::rotate(0, 0, angle));
  }
  else if (uc_case == 9) {
    // UC08 3 cubes
    // UC47 5 icosahedra
    // UC49 5 great dodecahedra
    // UC51 5 small stellated dodecahedra
    // UC53 5 great icosahedra
    transform_and_repeat(geom, sym_to, sym_from, Trans3d::rotate(angle, 0, 0));
  }
}

// M_PI/12 = 15 degrees
// M_PI/6  = 30 degrees
// M_PI/5  = 36 degrees
// M_PI/4  = 45 degrees
// M_PI/3  = 60 degrees
// M_PI/2  = 90 degrees

// clang-format off
UCItem uc_item_list[] = {
   { 1,   "u1",       "S4",    "T",     -1,      "UC1",   "6 tetrahedra rotational"},
   { 1,   "u1",       "S4",    "T",     -1,      "UC2",   "12 tetrahedra rotational"},
   { 1,   "u1",       "S4",    "T",     M_PI/4,  "UC3",   "6 tetrahedra"},
   { 0,   "u1",       "Td",    "Oh",    -1,      "UC4",   "2 tetrahedra"},
   { 0,   "u1",       "Td",    "I",     -1,      "UC5",   "5 tetrahedra"},
   { 0,   "u1",       "Td",    "I",     -1,      "UC6",   "10 tetrahedra"},
   { 2,   "u6",       "S4",    "T",     -1,      "UC7",   "6 cubes rotational"},
   { 9,   "u6",       "Oh",    "T",     M_PI/4,  "UC8",   "3 cubes"},
   { 0,   "u6",       "Oh",    "I",     -1,      "UC9",   "5 cubes"},
   { 3,   "ant3",     "D3v",   "T",     -1,      "UC10",  "4 octahedra rotational"},
   { 3,   "ant3",     "D3v",   "O",     -1,      "UC11",  "8 octahedra rotational"},
   { 3,   "ant3",     "D3v",   "T",     M_PI/12, "UC12",  "4 octahedra"},
   { 3,   "ant3",     "D3v",   "I",     -1,      "UC13",  "20 octahedra rotational"},
   { 3,   "ant3",     "D3v",   "I",     -1,      "UC14",  "20 octahedra"},
   { 4,   "ant3",     "D3v",   "I",      0,      "UC15",  "10 octahedra 1"},
   { 4,   "ant3",     "D3v",   "I",     M_PI/3,  "UC16",  "10 octahedra 2"},
   { 0,   "u5",       "Oh",    "I",     -1,      "UC17",  "5 octahedra"},
   { 0,   "u4",       "Oh",    "I",     -1,      "UC18",  "5 tetrahemihexahedra"},
   { 4,   "u4",       "D3v",   "I",     -1,      "UC19",  "20 tetrahemihexahedra"},
   { 5,   "pri",      "D",     "D",     -1,      "UC20",  "2k n d gonal prisms rotational"},
   { 5,   "pri",      "D",     "D",     -1,      "UC21",  "k n d gonal prisms"},
   { 5,   "ant",      "D",     "D",     -1,      "UC22",  "2k n odd d gonal antiprisms rotational"},
   { 5,   "ant",      "D",     "D",     -1,      "UC23",  "k n odd d gonal antiprisms"},
   { 5,   "ant",      "D",     "D",     -1,      "UC24",  "2k n even d gonal antiprisms rotational"},
   { 5,   "ant",      "D",     "D",     -1,      "UC25",  "k n even d gonal antiprisms"},
   { 6,   "ant5",     "D5h",   "I",     -1,      "UC26",  "12 pentagonal antiprisms rotational"},
   { 7,   "ant5",     "D5h",   "I",      0,      "UC27",  "6 pentagonal antiprisms"},
   { 6,   "ant5/3",   "D5h",   "I",     -1,      "UC28",  "12 pentagrammic crossed antiprisms rotational"},
   { 7,   "ant5/3",   "D5h",   "I",     M_PI/5,  "UC29",  "6 pentagrammic crossed antiprisms"},
   { 3,   "pri3",     "D3h",   "O",     M_PI/12, "UC30",  "4 triangular prisms"},
   { 3,   "pri3",     "D3h",   "O",     M_PI/12, "UC31",  "8 triangular prisms"},
   { 8,   "pri3",     "D3h",   "I",     M_PI/6,  "UC32",  "10 triangular prisms"},
   { 8,   "pri3",     "D3h",   "I",     M_PI/6,  "UC33",  "20 triangular prisms"},
   { 7,   "pri5",     "D5h",   "I",     M_PI/5,  "UC34",  "6 pentagonal prisms"},
   { 7,   "pri5",     "D5h",   "I",     M_PI/5,  "UC35",  "12 pentagonal prisms"},
   { 7,   "pri5/2",   "D5h",   "I",     M_PI/5,  "UC36",  "6 pentagrammic prisms"},
   { 7,   "pri5/2",   "D5h",   "I",     M_PI/5,  "UC37",  "12 pentagrammic prisms"},
   { 3,   "pri6",     "D6h",   "O",     M_PI/12, "UC38",  "4 hexagonal prisms"},
   { 8,   "pri6",     "D6h",   "I",     M_PI/6,  "UC39",  "10 hexagonal prisms"},
   { 7,   "pri10",    "D10h",  "I",     M_PI/5,  "UC40",  "6 decagonal prisms"},
   { 7,   "pri10/3",  "D10h",  "I",     M_PI/5,  "UC41",  "6 decagrammic prisms"},
   { 0,   "ant4",     "D4h",   "T",     -1,      "UC42",  "3 square antiprisms"},
   { 0,   "ant4",     "D4h",   "T",     -1,      "UC43",  "6 square antiprisms"},
   { 7,   "ant5/2",   "D5h",   "I",     M_PI/5,  "UC44",  "6 pentagrammic antiprisms"},
   { 7,   "ant5/2",   "D5h",   "I",     M_PI/5,  "UC45",  "12 pentagrammic antiprisms"},
   { 0,   "u22",      "Ih",    "O",     -1,      "UC46",  "2 icosahedra"},
   { 9,   "u22",      "Th",    "I",     M_PI/2,  "UC47",  "5 icosahedra"},
   { 0,   "u35",      "Ih",    "O",     -1,      "UC48",  "2 great dodecahedra"},
   { 9,   "u35",      "Th",    "I",     M_PI/2,  "UC49",  "5 great dodecahedra"},
   { 0,   "u34",      "Ih",    "O",     -1,      "UC50",  "2 small stellated dodecahedra"},
   { 9,   "u34",      "Th",    "I",     M_PI/2,  "UC51",  "5 small stellated dodecahedra"},
   { 0,   "u53",      "Ih",    "O",     -1,      "UC52",  "2 great icosahedra"},
   { 9,   "u53",      "Th",    "I",     M_PI/2,  "UC53",  "5 great icosahedra"},
   { 0,   "u2",       "Td",    "Oh",    -1,      "UC54",  "2 truncated tetrahedra"},
   { 0,   "u2",       "Td",    "I",     -1,      "UC55",  "5 truncated tetrahedra"},
   { 0,   "u2",       "Td",    "I",     -1,      "UC56",  "10 truncated tetrahedra"},
   { 0,   "u9",       "Td",    "I",     -1,      "UC57",  "5 truncated cubes"},
   { 0,   "u19",      "Td",    "I",     -1,      "UC58",  "5 stellated truncated hexahedra"},
   { 0,   "u7",       "Td",    "I",     -1,      "UC59",  "5 cuboctahedra"},
   { 0,   "u15",      "Td",    "I",     -1,      "UC60",  "5 cubohemioctahedra"},
   { 0,   "u3",       "Td",    "I",     -1,      "UC61",  "5 octahemioctahedra"},
   { 0,   "u10",      "Td",    "I",     -1,      "UC62",  "5 rhombicuboctahedra"},
   { 0,   "u18",      "Td",    "I",     -1,      "UC63",  "5 small rhombihexahedra"},
   { 0,   "u13",      "Td",    "I",     -1,      "UC64",  "5 small cubicuboctahedra"},
   { 0,   "u14",      "Td",    "I",     -1,      "UC65",  "5 great cubicuboctahedra"},
   { 0,   "u21",      "Td",    "I",     -1,      "UC66",  "5 great rhombihexahedra"},
   { 0,   "u17",      "Td",    "I",     -1,      "UC67",  "5 great rhombicuboctahedra"},
   { 0,   "u12",      "I",     "Ih",    -1,      "UC68",  "2 snub cubes"},
   { 0,   "u29",      "I",     "Ih",    -1,      "UC69",  "2 snub dodecahedra"},
   { 0,   "u57",      "I",     "Ih",    -1,      "UC70",  "2 great snub icosidodecahedra"},
   { 0,   "u69",      "I",     "Ih",    -1,      "UC71",  "2 great inverted snub icosidodecahedra"},
   { 0,   "u74",      "I",     "Ih",    -1,      "UC72",  "2 great retrosnub icosidodecahedra"},
   { 0,   "u40",      "I",     "Ih",    -1,      "UC73",  "2 snub dodecadodecahedra"},
   { 0,   "u60",      "I",     "Ih",    -1,      "UC74",  "2 inverted snub dodecadodecahedra"},
   { 0,   "u46",      "I",     "Ih",    -1,      "UC75",  "2 snub icosidodecadodecahedra"},
};
// clang-format on

UniformCompound::UniformCompound()
{
  uc_items = uc_item_list;
  last_uc = sizeof(uc_item_list) / sizeof(uc_item_list[0]);
}

int UniformCompound::get_poly(Geometry &geom, int sym, double angle, int n,
                              int d, int k, bool is_std)
{
  string constituent_str = uc_items[sym].constituent;
  if (is_std)
    constituent_str = "std_" + constituent_str;

  string sym_from = uc_items[sym].sym_from;
  string sym_to = uc_items[sym].sym_to;

  int uc_num = sym + 1;
  if (uc_num >= 20 && uc_num <= 25) {
    string tmp_str = std::to_string(n) + "/" + std::to_string(d);
    constituent_str += tmp_str;

    string sym_from_local =
        "D" + std::to_string(n) + ((uc_num == 22 || uc_num == 23) ? "v" : "h");
    string sym_to_local =
        "D" + std::to_string(k * n) +
        (((uc_num == 22 || uc_num == 23) && !is_even(k)) ? "v" : "h");

    sym_from = sym_from_local;
    sym_to = sym_to_local;
  }
  else if (uc_num == 14)
    angle = octahedra_angle(2) +
            octahedra_angle(0); // add 45.6696674755810220... degrees
  else if (uc_num == 19)
    angle = octahedra_angle(0);
  else if (uc_items[sym].angle != -1)
    angle = uc_items[sym].angle;

  geom.read_resource(constituent_str);

  build_uniform_compound(geom, uc_items[sym].uc_case, uc_num, sym_from, sym_to,
                         angle);

  return 1;
}

int UniformCompound::lookup_sym_no(string sym)
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

  // remove spaces either side of a punctuation mark
  string sym_norm2;
  for (unsigned int i = 0; i < sym_norm.length(); i++) {
    if (sym_norm[i] == ' ' &&
        ((i > 0 && ispunct(sym_norm[i - 1])) ||
         (i < sym_norm.length() && ispunct(sym_norm[i + 1]))))
      continue;
    sym_norm2 += sym_norm[i];
  }

  // sym_norm2 is now normalised

  // is it blank
  if (sym_norm2 == "")
    return -1;

  int idx = -1;

  // is it a poly name or description
  for (int i = 0; i < get_last_uc(); i++) {
    if (!strncasecmp(sym_norm2.c_str(), uc_item_list[i].name,
                     sym_norm2.size()) ||
        !strncasecmp(sym_norm2.c_str(), uc_item_list[i].description,
                     sym_norm2.size())) {
      idx = i;
      break;
    }
  }

  return idx;
}

void UniformCompound::assign_uc_value(char operand, const char *digits_str,
                                      double &angle, int &n, int &d, int &k)
{
  if (operand == 'a')
    angle = atof(digits_str);
  else if (operand == 'n')
    n = atoi(digits_str);
  else if (operand == '/')
    d = atoi(digits_str);
  else if (operand == 'k')
    k = atoi(digits_str);
}

int UniformCompound::parse_uc_args(string &name, double &angle, int &n, int &d,
                                   int &k, string *error_msg = nullptr)
{
  int ret = 0;

  int loc = 0;
  int loc1 = name.rfind("_");
  int loc2 = name.rfind(" ");
  if (loc1 == (int)string::npos || loc1 < loc2)
    loc = loc2;
  else if (loc2 == (int)string::npos || loc1 > loc2)
    loc = loc1;

  string uc_name = name.substr(0, loc);

  if (loc + 1 >= (int)name.length()) {
    if (error_msg)
      *error_msg = "argument string not found";
  }
  // process uc args
  else {
    string and_the_rest = name.substr((loc + 1));

    string operators = "an/k";
    string digits = "0123456789.-";
    char operand = '\0';
    string digits_str;

    for (char i : and_the_rest) {
      if (operators.find(i) != string::npos) {
        if (operand) {
          // don't accept empty digit string or "." at end of string
          if (!digits_str.length() ||
              digits_str[digits_str.length() - 1] == '.') {
            if (error_msg)
              *error_msg = "no digits found, or decimal point at end";
            break;
          }
          else
            assign_uc_value(operand, digits_str.c_str(), angle, n, d, k);
        }
        digits_str.clear();
        operand = i;
      }
      else if (digits.find(i) != string::npos) {
        if (!operand) {
          if (error_msg)
            *error_msg = "operator expected";
          break;
        }
        else if (operand == '/' && n < 1) {
          if (error_msg)
            *error_msg = "d of n/d supplied but n is zero or not set";
          break;
        }
        else if (operand != 'a' && (i == '.' || i == '-')) {
          if (error_msg)
            *error_msg =
                msg_str("operator %c should have a positive integer", operand);
          break;
        }
        if ((digits_str.find('.') != string::npos) && i == '.') {
          if (error_msg)
            *error_msg = "decimal point encountered more than once";
          break;
        }
        else
          digits_str += i;
      }
      else {
        if (error_msg)
          *error_msg = msg_str("unexpected character: %c", i);
        break;
      }
    }

    if (operand) {
      // don't accept empty digit string or "." at end of string
      if (!digits_str.length() || digits_str[digits_str.length() - 1] == '.') {
        if (error_msg)
          *error_msg = "no digits found, or decimal point at end";
      }
      else
        assign_uc_value(operand, digits_str.c_str(), angle, n, d, k);
    }
  }

  if (error_msg) {
    if (error_msg->empty()) {
      if (n == 0)
        *error_msg = "operator n must not be 0";
      else if (d == 0)
        *error_msg = "operator / must not be 0";
      else if (k == 0)
        *error_msg = "operator k must not be 0";
      else if (n > 0 && d > 0 && gcd(n, d) != 1)
        *error_msg = "n and d must be co-prime";
    }
  }

  if (error_msg) {
    if (!error_msg->empty())
      ret = 1; // fail
  }

  if (!ret) {
    name = uc_name;
    if (n != -1 &&
        d == -1) // if n is set, set d or it will be randomly selected
      d = 1;
    if (angle != INFINITY)
      angle = deg2rad(angle);
  }

  return ret;
}

int UniformCompound::set_uc_args(int sym, double &angle, int &n, int &d, int &k,
                                 string *error_msg = nullptr)
{
  int need_angle[] = {1,  2,  7,  10, 11, 13,
                      20, 22, 24, 26, 28}; // 11 occurrances
  bool needs_angle = false;
  for (int i : need_angle) {
    if (sym == i) {
      needs_angle = true;
      break;
    }
  }

  Random ran;
  ran.time_seed();

  if (needs_angle && angle == INFINITY) {
    angle = ran.ran_in_range_exclude_end(0, 360); // 0 to 360
    angle = deg2rad(angle);
  }
  else if (!needs_angle && angle != INFINITY) {
    if (error_msg)
      *error_msg = msg_str("for UC%d, angle is not needed", sym);
    return 1;
  }

  // n,d,k are used only in 20 to 25
  bool uc20_25 = (sym >= 20 && sym <= 25);
  if (uc20_25 && n == -1)
    n = ran.ran_int_in_range(2, 20); // range 2 to 20
  else if (!uc20_25 && n != -1) {
    if (error_msg)
      *error_msg = msg_str("for UC%d, n is not needed", sym);
    return 1;
  }

  if (uc20_25 && d == -1) {
    if (sym == 20 || sym == 21) { // n and d must be co-prime
      if (n < 3)                  // n = 2 would be a square
        n = 3;
      while (!(gcd(n, d) == 1 && (double)n / d > 2.0))
        d = ran.ran_int_in_range(1, 10); // range 1 to 10
    }
    else if (sym == 22 || sym == 23) { // n and d must be co-prime, d must be
                                       // odd, so n can be 2 or greater
      if (n < 5) // n = 2 is a tetrahedron and is allowed (2/1, 3/1, 4/1)
        d = 1;
      else {
        while (!(!is_even(d) && gcd(n, d) == 1 && (double)n / d > 3.0 / 2))
          d = ran.ran_int_in_range(1, 10); // range 1 to 10
      }
    }
    else if (sym == 24 || sym == 25) { // n and d must be co-prime, d must be
                                       // even, so n must be odd (recalculate),
                                       // 5 or greater, ant3/2 would be a
                                       // triangle
      while (is_even(n) || n < 5)
        n = ran.ran_int_in_range(5, 20); // range 5 to 20
      while (!(is_even(d) && gcd(n, d) == 1 && (double)n / d > 3.0 / 2))
        d = ran.ran_int_in_range(2, 10); // range 2 to 10
    }
  }
  else if (!uc20_25 && d != -1) {
    if (error_msg)
      *error_msg = msg_str("for UC%d, d is not needed", sym);
    return 1;
  }

  int k_min = 0;
  if (uc20_25)
    k_min = (is_even(sym)) ? 1 : 2;

  if (uc20_25 && k == -1)
    k = ran.ran_int_in_range(k_min, 4); // range k_min to 4
  else if (!uc20_25 && k != -1) {
    if (error_msg)
      *error_msg = msg_str("for UC%d, k is not needed", sym);
    return 1;
  }

  // a few more checks for prismatic, now that sym is known, n/d and k are
  // checked
  if (uc20_25) {
    if (k < k_min) {
      if (error_msg)
        *error_msg = msg_str("for UC%d, k must be greater than %d", sym, k_min);
      return 1;
    }

    if (sym == 20 || sym == 21) {
      // RK - this constraint is not necessary for prisms
      // if ((double)n/d <= 2.0) {
      //   if (error_msg)
      //     *error_msg = msg_str("for UC%d, n/d (%d/%d) must be greater
      //   than 2",sym,n,d).c_str()); return 1;
      //}
    }
    else if (sym >= 22) {
      if ((sym == 22 || sym == 23) && is_even(d)) {
        if (error_msg)
          *error_msg = msg_str("for UC%d, d must be odd", sym);
        return 1;
      }
      else if ((sym == 24 || sym == 25) && !is_even(d)) {
        if (error_msg)
          *error_msg = msg_str("for UC%d, d must be even", sym);
        return 1;
      }
      if ((double)n / d <= 3.0 / 2) {
        if (error_msg)
          *error_msg = msg_str("for UC%d, n/d (%d/%d) must be greater than 3/2",
                               sym, n, d);
        return 1;
      }
    }
  }

  return 0;
}
