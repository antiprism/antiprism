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

/* !\file private_std_polys.h
 * \brief Platonic and Archimedean polyhedra
 */

#ifndef PRIVATE_STD_POLYS_H
#define PRIVATE_STD_POLYS_H

#include "geometry.h"
#include "geometryutils.h"

#include <cstdio>
#include <cstring>
#include <string>

/// Characters that separate parts of a resource name
const char RES_SEPARATOR[] = "_ ";

/// A basic function type that makes a model.
typedef void (*model_func)(anti::Geometry &);

/// Put face list into a normalized form.
/**Sometimes a model is created and the face list may not be consistent
 * between machines, perhaps for numeric reasons. This function puts
 * the list in a consistent form.
 * \param geom the geometry whose faces list will be put into the
 * normalized form */
void normalized_face_list(anti::Geometry &geom);

anti::Status make_resource_geom(anti::Geometry &geom, std::string name);

struct UniformItem {
  model_func pfunc;
  std::string Wythoff;
  short Kaleido, Coxeter, Wenninger;
  std::string name;
  std::string dual;
};

class Uniform {
private:
  int last_uniform;
  UniformItem *uniform_items;

public:
  Uniform();
  int get_poly(anti::Geometry &geom, int sym);
  int lookup_sym_no(std::string sym, int is_dual);
  int get_last_uniform() { return last_uniform; }
};

struct UCItem {
  int uc_case;
  const char *constituent;
  const char *sym_from;
  const char *sym_to;
  double angle;
  const char *name;
  const char *description;
};

class UniformCompound {
private:
  int last_uc;
  UCItem *uc_items;

public:
  UniformCompound();
  int get_poly(anti::Geometry &geom, int sym, double angle, int n, int d, int k,
               bool is_std);
  int lookup_sym_no(std::string sym);
  int get_last_uc() { return last_uc; }

  void assign_uc_value(char operand, const char *digits_str, double &angle,
                       int &n, int &d, int &k);
  int parse_uc_args(std::string &name, double &angle, int &n, int &d, int &k,
                    std::string *error_msg);
  int set_uc_args(int sym, double &angle, int &n, int &d, int &k,
                  std::string *error_msg);
};

struct JohnsonItem {
  model_func pfunc;
  const char *name;
};

class Johnson {
private:
  int last_J;
  JohnsonItem *J_items;

public:
  Johnson();
  int get_poly(anti::Geometry &geom, int sym);
  int lookup_sym_no(std::string sym);
  int get_last_J() { return last_J; }
};

class Wythoff {
private:
  std::vector<int> fracs;
  std::vector<anti::Vec3d> verts;
  int bar_pos;

  anti::Status read_symbol(const char *sym = nullptr);
  std::vector<int> map_to_min_triangle();
  bool assign_verts();

public:
  Wythoff(const char *sym, anti::Status *status);
  anti::Status make_poly(anti::Geometry &geom);
  bool make_tri_poly(anti::Geometry &geom);
  bool make_tri(anti::Geometry &geom);
  bool is_set() { return bar_pos != -1; }

  std::string to_str();
  std::string get_tri_sym();
};

struct WenningerItem {
  int number;
  const char *poly;
  const char *diagram_str1;
  const char *diagram_str2;
  const char *diagram_str3;
  const char *diagram_str4;
  const char *sub_sym;
  bool remove_inline_verts;
};

class Wenninger {
private:
  WenningerItem *Wenninger_items;

public:
  Wenninger();
  int get_poly(anti::Geometry &geom, int sym);
  bool test_valid(int idx);
};

#endif // PRIVATE_STD_POLYS_H
