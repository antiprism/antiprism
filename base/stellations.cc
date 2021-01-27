/*
   Copyright (c) 2017-2021, Roger Kaufman

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
   Name: stellations.cc
   Description: stellations added for Wenninger numbers
   Project: Antiprism - http://www.antiprism.com
*/

#include "planar.h"
#include "private_std_polys.h"
#include "utils.h"

#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>

using std::map;
using std::string;
using std::vector;

using namespace anti;

// clang-format off
// Table of Wenninger Models which are not Uniforms
// Number, base polyhedron, diagram numbers 1, diagram numbers 2, symmetry, remove inline vertices
WenningerItem wenninger_item_list[] = {
   { 19, "oct",     "0,1",                  "",                        "",  "",  "Oh", true  }, // compound of two tetrahedra, also uc4
   { 23, "ico",     "0,30",                 "",                        "",  "",  "Ih", true  }, // compound of five octahedra, also uc17
   { 24, "ico",     "0,16,31,43",           "0,32,57",                 "",  "",  "I",  true  }, // compound of five tetrahedra, also uc5
   { 25, "ico",     "0,16",                 "0,32,57",                 "",  "",  "Ih", true  }, // compound of ten tetrahedra, also uc6
   { 26, "ico",     "0,46",                 "",                        "",  "",  "Ih", true  }, // small triambic icosahedron, also u30_d
   { 27, "ico",     "0,32,34,57",           "",                        "",  "",  "Ih", true  }, // second stellation of icosahedron
   { 28, "ico",     "0,16",                 "0,57",                "0,11",  "",  "Ih", true  }, // third stellation of icosahedron, excavated dodecahedron
   { 29, "ico",     "0,11,32",              "",                        "",  "",  "Ih", true  }, // fourth stellation of icosahedron
   { 30, "ico",     "0,16",                 "0,33,57",                 "",  "",  "Ih", true  }, // fifth stellation of icosahedron
   { 31, "ico",     "0,34,57",              "0,11",                    "",  "",  "Ih", true  }, // sixth stellation of icosahedron
   { 32, "ico",     "0,29",                 "0,43",                    "",  "",  "Ih", true  }, // seventh stellation of icosahedron
   { 33, "ico",     "0,32",                 "0,33",                    "",  "",  "Ih", true  }, // eighth stellafion of icosahedron, great triambic icosahedron, also u47_d
   { 34, "ico",     "0,44",                 "0,31",                "0,34",  "",  "Ih", true  }, // ninth stellation of icosahedron
   { 35, "ico",     "0,14,32,57,62",        "",                        "",  "",  "I",  true  }, // tenth stellation of icosahedron
   { 36, "ico",     "0,43,44",              "0,32,57,62",              "",  "",  "I",  true  }, // eleventh stellation of icosahedron
   { 37, "ico",     "0,43,44",              "",                        "",  "",  "Ih", true  }, // twelfth stellation of icosahedron
   { 38, "ico",     "0,17,29,31",           "0,14,32,57",          "0,34",  "",  "I",  true  }, // thirteenth stellation of icosahedron
   { 39, "ico",     "0,43,44",              "0,57,62",          "0,11,61",  "",  "I",  true  }, // fourteenth stellation of icosahedron
   { 40, "ico",     "0,16,17,29,31",        "0,14,32,57",              "",  "",  "I",  false }, // fifteenth stellation of icosahedron
   { 42, "ico",     "0,8,12",               "",                        "",  "",  "Ih", true  }, // final stellation of the icosahedron
   { 43, "cubo",    "0,19",                 "8,25",                    "",  "",  "Oh", true  }, // compound of cube and octahedron
   { 44, "cubo",    "0,2",                  "8,26",                    "",  "",  "Oh", true  }, // second stellation of cuboctahedron
   { 45, "cubo",    "0,1,24",               "8,9,10",                  "",  "",  "Oh", true  }, // third stellation of cuboctahedron
   { 46, "cubo",    "0,5,20,24",            "8,7,8,9",                 "",  "",  "Oh", true  }, // fourth stellation of cuboctahedron
   { 47, "icosid",  "0,142",                "20,125",                  "",  "",  "Ih", true  }, // compound of dodecahedron and icosahedron
   { 48, "icosid",  "0,143",                "20,132",                  "",  "",  "Ih", true  }, // second stellation of icosidodecahedron
   { 49, "icosid",  "0,94,223",             "20,127,180",              "",  "",  "Ih", true  }, // third stellation of icosidodecahedron
   { 50, "icosid",  "0,28,223",             "20,127,182,211",    "20,181",  "",  "Ih", true  }, // compound of small stellated dodecahedron and triakis icosahedron
   { 51, "icosid",  "0,93,160",             "20,182,211",        "20,181",  "",  "Ih", true  }, // compound of small stellated dodecahedron and five octahedra
   { 52, "icosid",  "0,93",                 "20,130,180,182,211","20,181",  "",  "Ih", true  }, // sixth stellation of icosidodecahedron
   { 53, "icosid",  "0,195,196",            "20,210,228",              "",  "",  "Ih", true  }, // seventh stellation of icosidodecahedron
   { 54, "icosid",  "0,23,24,112,113,114,115,117,164,165,197",
                    "20,23,228",                                       "",  "",  "I",  true  }, // Compound of five tetrahedra and great dodecahedron
   { 55, "icosid",  "0,61,62,109,110,197",  "0,113",
                    "20,130,180",           "20,131,210",                        "Ih", true  }, // ninth stellation of icosidodecahedron
   { 56, "icosid",  "0,108,109,110,197",    "0,111,113,222",
                    "20,132,181",           "20,127,182",                        "Ih", true  }, // tenth stellation of icosidodecahedron
   { 57, "icosid",  "0,26,118",             "0,57",            "20,24,25",  "",  "Ih", false }, // eleventh stellation of icosidodecahedron
   { 58, "icosid",  "0,63,108,110",         "0,111,113",       "20,24,26",  "",  "Ih", false }, // twelfth stellation of icosidodecahedron
   { 59, "icosid",  "0,58,59",              "0,168,227",
                    "20,28,30,86",          "20,26,27,29,87,215",                "Ih", true  }, // thirteenth stellation of icosidodecahedron
   { 60, "icosid",  "0,60,61,168,220",      "",
                    "20,28,30,86",          "20,29,87,151,215",                  "Ih", false }, // fourteenth stellation of icosidodecahedron
   { 61, "icosid",  "0,60,90,111,220",      "",
                    "20,30,86",             "20,29,87,215",                      "Ih", true  }, // compound of great stellated dodecahedron and great icosahedron
   { 62, "icosid",  "0,111",                "0,60,63,108,227", "20,209,214", "", "Ih", false }, // fifteenth stellation of icosidodecahedron
   { 63, "icosid",  "0,92,95",              "20,210,228",              "",  "",  "Ih", true  }, // sixteenth stellation of icosidodecahedron
   { 64, "icosid",  "0,62,63",              "",
                    "20,28,30,86",          "20,29,87,215",                      "Ih", false }, // seventeenth stellation of icosidodecahedron
   { 65, "icosid",  "0,60,61,90,220",       "20,24,25,26,27",          "",  "",  "Ih", true  }, // eighteenth stellation of icosidodecahedron
   { 66, "icosid",  "0,67,70",              "0,66,68,69,103",
                    "20,82,87",             "20,30,86",                          "Ih", false }, // nineteenth stellation of icosidodecahedron

};
// clang-format on

Wenninger::Wenninger() { Wenninger_items = wenninger_item_list; }

// idx changed locally
bool Wenninger::test_valid(int idx)
{
  // Wenninger number range has gaps
  idx++;
  return (idx == 19 || (idx >= 23 && idx <= 40) || (idx >= 42 && idx <= 66))
             ? true
             : false;
}

// sym changed locally
int Wenninger::get_poly(Geometry &geom, int sym)
{
  // Wenninger number range has gaps, reset to table indexes
  sym++;
  if (sym == 19)
    sym = 0;
  else if (sym >= 23 && sym <= 40)
    sym -= 22;
  else if (sym >= 42 && sym <= 66)
    sym -= 23;

  geom.read_resource(Wenninger_items[sym].poly);

  vector<string> diagram_list_strings;
  if (strlen(Wenninger_items[sym].diagram_str1))
    diagram_list_strings.push_back(Wenninger_items[sym].diagram_str1);
  if (strlen(Wenninger_items[sym].diagram_str2))
    diagram_list_strings.push_back(Wenninger_items[sym].diagram_str2);
  if (strlen(Wenninger_items[sym].diagram_str3))
    diagram_list_strings.push_back(Wenninger_items[sym].diagram_str3);
  if (strlen(Wenninger_items[sym].diagram_str4))
    diagram_list_strings.push_back(Wenninger_items[sym].diagram_str4);

  string sym_str = Wenninger_items[sym].sub_sym;

  int sz = diagram_list_strings.size();
  vector<vector<int>> idx_lists(sz);

  map<int, Geometry> diagrams;

  // data sizes are verified
  for (unsigned int i = 0; i < diagram_list_strings.size(); i++) {
    if (!diagram_list_strings[i].length())
      continue;

    read_idx_list(diagram_list_strings[i].c_str(), idx_lists[i],
                  std::numeric_limits<int>::max(), false);

    // stellation face index is in the first position
    int stellation_face_idx = idx_lists[i][0];

    // construct the diagrams
    if (!diagrams[stellation_face_idx].verts().size())
      diagrams[stellation_face_idx] =
          make_stellation_diagram(geom, stellation_face_idx, sym_str);
  }

  bool merge_faces = true;
  bool remove_inline_verts = Wenninger_items[sym].remove_inline_verts;
  bool split_pinched = true;
  bool resolve_faces = true;
  bool remove_multiples = true;
  string map_string = "compound";

  // these models need lower precision
  int number = Wenninger_items[sym].number;
  double local_epsilon = epsilon;
  if (number == 59 || number == 60 || number == 66)
    local_epsilon = 1e-11;

  geom = make_stellation(geom, diagrams, idx_lists, sym_str, merge_faces,
                         remove_inline_verts, split_pinched, resolve_faces,
                         remove_multiples, map_string, local_epsilon);

  return 1;
}
