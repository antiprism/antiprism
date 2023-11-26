/*
   Copyright (c) 2003-2023, Adrian Rossiter, Roger Kaufman

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
   Name: color_common.h
   Description: color code shared in /src
   Project: Antiprism - http://www.antiprism.com
*/

#ifndef COLOR_COMMON_H
#define COLOR_COMMON_H

#include <cstdio>
#include <set>
#include <string>
#include <vector>

using std::set;
using std::string;
using std::vector;

#include "../base/antiprism.h"

class OffColor {
public:
  OffColor() { read_colorings(clrngs, "spread"); }
  OffColor(string s)
  {
    if (s.length())
      read_colorings(clrngs, s.c_str());
  }

  void set_v_col_op(char v);
  void set_v_col(const anti::Color v);
  void set_v_sub_sym(const string v);

  void set_e_col_op(char e);
  void set_e_col(const anti::Color e);
  void set_e_sub_sym(const string e);
  void set_e_min_len_diff(const double e);

  void set_f_col_op(char f);
  void set_f_col(const anti::Color f);
  void set_f_sub_sym(const string f);

  char get_v_col_op();
  anti::Color get_v_col();
  string get_v_sub_sym();

  char get_e_col_op();
  anti::Color get_e_col();
  string get_e_sub_sym();
  double get_e_min_len_diff();

  char get_f_col_op();
  anti::Color get_f_col();
  string get_f_sub_sym();

  // array of private colorings is problematical. public for now
  anti::Coloring clrngs[3];

  /* not sure about this code. direct access for now
  void set_clrng(anti::Coloring &clrng, int n);
  anti::Coloring get_clrng(int n);
  */

  // if op_str is not set, all options would be available
  bool v_op_check(char *v_col_op, const char *op_str = "uUpPsSnNaAFEcCLM");
  bool e_op_check(char *e_col_op,
                  const char *op_str = "uUpPsSnNkKFVdDjJgGcCLlM");
  bool f_op_check(char *f_col_op, const char *op_str = "uUpPsSnNaAkKEVgGcCLlM");

  anti::Status off_color_main(anti::Geometry &geom);

  ~OffColor() = default;

private:
  char v_col_op = '\0';
  string v_sub_sym;
  vector<set<int>> v_equivs;
  anti::Color v_col;

  char e_col_op = '\0';
  string e_sub_sym;
  vector<set<int>> e_equivs;
  anti::Color e_col;
  double e_min_len_diff;

  char f_col_op = '\0';
  string f_sub_sym;
  vector<set<int>> f_equivs;
  anti::Color f_col;

  // operator when value is present
  char val_op = '$';

  // if true, when lower case input then upper case operations are issued
  bool upper = true;
};

// global transparency call
anti::Status apply_transparency(anti::Geometry &geom, const int opacity,
                                const int elem = anti::FACES);

// global transparency call
void apply_transparencies(anti::Geometry &geom, const int (&opacity)[3]);

anti::ColorMapMap *alloc_no_colorMap();

// color edges by dihedral angle
void color_edges_by_dihedral(anti::Geometry &geom, anti::Coloring &clrng,
                             bool apply_map = true, double eps = anti::epsilon);

// color faces by convexity compared to other faces
void color_faces_by_convexity(anti::Geometry &geom, anti::Coloring &clrng,
                              bool apply_map = true,
                              double eps = anti::epsilon);

// color faces by connection to other faces
void color_faces_by_connection(anti::Geometry &geom, anti::Coloring &clrng,
                               bool apply_map = true);

// color faces by connection to other faces
// when vertices and edges are also handled. wraps off_color_main()
anti::Status color_faces_by_connection_vef(anti::Geometry &geom,
                                           OffColor &off_color);

// for lat_util.cc, bravais.cc and waterman.cc
void color_by_symmetry_normals(anti::Geometry &, const char, const int,
                               double eps = anti::epsilon);

void color_edges_by_sqrt(anti::Geometry &, const char);

#endif // COLOR_COMMON_H
