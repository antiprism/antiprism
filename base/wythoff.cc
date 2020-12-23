/*
   Copyright (c) 2012-2019, Adrian Rossiter

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

#include "geometry.h"
#include "geometryinfo.h"
#include "private_std_polys.h"
#include "programopts.h"
#include "symmetry.h"
#include "tiling.h"
#include "utils.h"

#include <algorithm>
#include <cstring>
#include <functional>
#include <regex>
#include <set>
#include <string>
#include <vector>

using std::map;
using std::pair;
using std::string;
using std::swap;
using std::to_string;
using std::vector;

using namespace anti;

// Schwarz triangles
const int num_schwarz_tris = 44;

// clang-format off
int schwarz_triangles[num_schwarz_tris][6] = {
    {2, 1, 3, 1, 3, 1}, //  0
    {2, 1, 3, 1, 3, 2}, //  1
    {2, 1, 3, 1, 4, 1}, //  2
    {2, 1, 3, 1, 4, 3}, //  3
    {2, 1, 3, 1, 5, 1}, //  4
    {2, 1, 3, 1, 5, 2}, //  5
    {2, 1, 3, 1, 5, 3}, //  6
    {2, 1, 3, 1, 5, 4}, //  7
    {2, 1, 3, 2, 3, 2}, //  8
    {2, 1, 3, 2, 4, 1}, //  9
    {2, 1, 3, 2, 4, 3}, // 10
    {2, 1, 3, 2, 5, 1}, // 11
    {2, 1, 3, 2, 5, 2}, // 12
    {2, 1, 3, 2, 5, 3}, // 13
    {2, 1, 3, 2, 5, 4}, // 14
    {2, 1, 5, 1, 5, 2}, // 15
    {2, 1, 5, 1, 5, 3}, // 16
    {2, 1, 5, 2, 5, 4}, // 17
    {2, 1, 5, 3, 5, 4}, // 18
    {3, 1, 3, 1, 3, 2}, // 19
    {3, 1, 3, 1, 5, 2}, // 20
    {3, 1, 3, 1, 5, 4}, // 21
    {3, 1, 3, 2, 5, 1}, // 22
    {3, 1, 3, 2, 5, 3}, // 23
    {3, 1, 4, 1, 4, 3}, // 24
    {3, 1, 5, 1, 5, 3}, // 25
    {3, 1, 5, 1, 5, 4}, // 26
    {3, 1, 5, 2, 5, 3}, // 27
    {3, 1, 5, 2, 5, 4}, // 28
    {3, 2, 3, 2, 3, 2}, // 29
    {3, 2, 3, 2, 5, 2}, // 30
    {3, 2, 3, 2, 5, 4}, // 31
    {3, 2, 4, 1, 4, 1}, // 32
    {3, 2, 4, 3, 4, 3}, // 33
    {3, 2, 5, 1, 5, 1}, // 34
    {3, 2, 5, 1, 5, 2}, // 35
    {3, 2, 5, 2, 5, 2}, // 36
    {3, 2, 5, 3, 5, 3}, // 37
    {3, 2, 5, 3, 5, 4}, // 38
    {3, 2, 5, 4, 5, 4}, // 39
    {5, 1, 5, 1, 5, 4}, // 40
    {5, 2, 5, 2, 5, 2}, // 41
    {5, 2, 5, 3, 5, 3}, // 42
    {5, 4, 5, 4, 5, 4}, // 43
};

static double schwarz_triangles_verts[num_schwarz_tris][9] = {
    { //  0
      0, 1, 0,
      1/sqrt_3, 1/sqrt_3, -1/sqrt_3,
      1/sqrt_3, 1/sqrt_3, 1/sqrt_3, },
    { //  1
      0, 1, 0,
      1/sqrt_3, -1/sqrt_3, 1/sqrt_3,
      1/sqrt_3, 1/sqrt_3, -1/sqrt_3, },
    { //  2
      1/sqrt_2, 1/sqrt_2, 0,
      1/sqrt_3, 1/sqrt_3, 1/sqrt_3,
      1, 0, 0, },
    { //  3
      1/sqrt_2, 1/sqrt_2, 0,
      -1/sqrt_3, -1/sqrt_3, -1/sqrt_3,
      1, 0, 0, },
    { //  4
      0.5/phi, phi/2, 0.5,
      1/sqrt_3, 1/sqrt_3, 1/sqrt_3,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
    { //  5
      0.5, -0.5/phi, phi/2,
      phi/sqrt_3, (phi-1)/sqrt_3, 0,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
    { //  6
      0, 1, 0,
      phi/sqrt_3, -(phi-1)/sqrt_3, 0,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
    { //  7
      0.5/phi, phi/2, 0.5,
      0, -phi/sqrt_3, -(phi-1)/sqrt_3,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
    { //  8
      0, -1, 0,
      1/sqrt_3, 1/sqrt_3, -1/sqrt_3,
      1/sqrt_3, 1/sqrt_3, 1/sqrt_3, },
    { //  9
      -1/sqrt_2, 1/sqrt_2, 0,
      -1/sqrt_3, 1/sqrt_3, -1/sqrt_3,
      1, 0, 0, },
    { // 10
      -1/sqrt_2, 0, -1/sqrt_2,
      1/sqrt_3, 1/sqrt_3, 1/sqrt_3,
      1, 0, 0, },
    { // 11
      0.5/phi, -phi/2, -0.5,
      1/sqrt_3, -1/sqrt_3, -1/sqrt_3,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
    { // 12
      0.5, 0.5/phi, -phi/2,
      phi/sqrt_3, -(phi-1)/sqrt_3, 0,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
    { // 13
      0, -1, 0,
      phi/sqrt_3, (phi-1)/sqrt_3, 0,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
    { // 14
      -0.5, -0.5/phi, -phi/2,
      1/sqrt_3, 1/sqrt_3, 1/sqrt_3,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
    { // 15
      0.5/phi, phi/2, 0.5,
      -1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, 0,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
    { // 16
      0.5/phi, phi/2, 0.5,
      1/sqrt_phi_plus_2, -phi/sqrt_phi_plus_2, 0,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
    { // 17
      0.5/phi, -phi/2, -0.5,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
      -1/sqrt_phi_plus_2, -phi/sqrt_phi_plus_2, 0, },
    { // 18
      0.5/phi, -phi/2, -0.5,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
      1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, 0, },
    { // 19
      1/sqrt_3, 1/sqrt_3, -1/sqrt_3,
      1/sqrt_3, -1/sqrt_3, 1/sqrt_3,
      1/sqrt_3, 1/sqrt_3, 1/sqrt_3, },
    { // 20
      1/sqrt_3, 1/sqrt_3, 1/sqrt_3,
      (phi-1)/sqrt_3, 0, phi/sqrt_3,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
    { // 21
      phi/sqrt_3, (phi-1)/sqrt_3, 0,
      -1/sqrt_3, -1/sqrt_3, 1/sqrt_3,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
    { // 22
      phi/sqrt_3, -(phi-1)/sqrt_3, 0,
      phi/sqrt_3, (phi-1)/sqrt_3, 0,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
    { // 23
      0, -phi/sqrt_3, -(phi-1)/sqrt_3,
      1/sqrt_3, 1/sqrt_3, 1/sqrt_3,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
    { // 24
      1/sqrt_3, 1/sqrt_3, -1/sqrt_3,
      0, 0, 1,
      0, 1, 0, },
    { // 25
      1/sqrt_3, 1/sqrt_3, 1/sqrt_3,
      0, -1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
    { // 26
      1/sqrt_3, 1/sqrt_3, 1/sqrt_3,
      -1/sqrt_phi_plus_2, -phi/sqrt_phi_plus_2, 0,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
    { // 27
      phi/sqrt_3, -(phi-1)/sqrt_3, 0,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
      1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, 0, },
    { // 28
      1/sqrt_3, -1/sqrt_3, -1/sqrt_3,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
      1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, 0, },
    { // 29
      1/sqrt_3, 1/sqrt_3, -1/sqrt_3,
      1/sqrt_3, -1/sqrt_3, 1/sqrt_3,
      -1/sqrt_3, 1/sqrt_3, 1/sqrt_3, },
    { // 30
      1/sqrt_3, -1/sqrt_3, -1/sqrt_3,
      (phi-1)/sqrt_3, 0, -phi/sqrt_3,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
    { // 31
      phi/sqrt_3, -(phi-1)/sqrt_3, 0,
      -1/sqrt_3, 1/sqrt_3, -1/sqrt_3,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
    { // 32
      1/sqrt_3, 1/sqrt_3, -1/sqrt_3,
      0, 1, 0,
      1, 0, 0, },
    { // 33
      1/sqrt_3, 1/sqrt_3, -1/sqrt_3,
      0, -1, 0,
      0, 0, 1, },
    { // 34
      1/sqrt_3, 1/sqrt_3, 1/sqrt_3,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
      1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, 0, },
    { // 35
      1/sqrt_3, 1/sqrt_3, 1/sqrt_3,
      1/sqrt_phi_plus_2, -phi/sqrt_phi_plus_2, 0,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
    { // 36
      phi/sqrt_3, (phi-1)/sqrt_3, 0,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
      1/sqrt_phi_plus_2, -phi/sqrt_phi_plus_2, 0, },
    { // 37
      1/sqrt_3, 1/sqrt_3, -1/sqrt_3,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
      1/sqrt_phi_plus_2, -phi/sqrt_phi_plus_2, 0, },
    { // 38
      1/sqrt_3, -1/sqrt_3, -1/sqrt_3,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
      0, 1/sqrt_phi_plus_2, -phi/sqrt_phi_plus_2, },
    { // 39
      1/sqrt_3, -1/sqrt_3, -1/sqrt_3,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
      -1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, 0, },
    { // 40
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
      1/sqrt_phi_plus_2, -phi/sqrt_phi_plus_2, 0,
      0, -1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, },
    { // 41
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
      1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, 0,
      -1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, 0, },
    { // 42
      1/sqrt_phi_plus_2, -phi/sqrt_phi_plus_2, 0,
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
      1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2, 0, },
    { // 43
      0, 1/sqrt_phi_plus_2, phi/sqrt_phi_plus_2,
      1/sqrt_phi_plus_2, -phi/sqrt_phi_plus_2, 0,
      0, 1/sqrt_phi_plus_2, -phi/sqrt_phi_plus_2, },
};
// clang-format on

static void frac_swap(vector<int> &fracs, vector<Vec3d> &vecs, int frac0,
                      int frac1)
{
  swap(fracs[2 * frac0], fracs[2 * frac1]);
  swap(fracs[2 * frac0 + 1], fracs[2 * frac1 + 1]);
  swap(vecs[frac0], vecs[frac1]);
}

static int frac_cmp(const int &n0, const int &d0, const int &n1, const int &d1)
{
  if (n0 > n1)
    return 1;
  else if (n0 < n1)
    return -1;
  else { // n0==n1
    if (d0 > d1)
      return 1;
    else if (d0 < d1)
      return -1;
    else
      return 0; // d0==d1
  }
}

static bool frac_less(const vector<int> &fracs, int frac0, int frac1)
{
  const int f0 = 2 * frac0;
  const int f1 = 2 * frac1;
  return frac_cmp(fracs[f0], fracs[f0 + 1], fracs[f1], fracs[f1 + 1]) == -1;
}

class frac_vect_less {
public:
  bool operator()(const vector<int> &f0, const vector<int> &f1) const
  {
    if (f0.size() != f1.size())
      return (f0.size() < f1.size());
    for (unsigned int i = 0; i < f0.size() / 2; i++) {
      int cmp = frac_cmp(f0[2 * i], f0[2 * i + 1], f1[2 * i], f1[2 * i + 1]);
      if (cmp == -1)
        return true; // less
      else if (cmp == 1)
        return false; // greater
    }
    return false; // equal
  }
};

static void tri_normalise(vector<int> &fracs, vector<Vec3d> &vecs)
{
  // bubble sort
  if (frac_less(fracs, 2, 1))
    frac_swap(fracs, vecs, 2, 1);
  if (frac_less(fracs, 1, 0))
    frac_swap(fracs, vecs, 1, 0);
  // smallest fraction is now in first place
  if (frac_less(fracs, 2, 1))
    frac_swap(fracs, vecs, 2, 1);
  // sorting complete
}

Wythoff::Wythoff(const char *sym, Status *status)
{
  if (!(*status = read_symbol(sym)))
    bar_pos = -1;
  else if (!assign_verts()) {
    status->set_error("symbol for non-finite construction (unsupported)");
    bar_pos = -1;
  }
}

string Wythoff::to_str()
{
  string sym;
  if (is_set()) {
    for (int i = 0; i < 3; i++) {
      if (bar_pos == i)
        sym += "|";
      else if (i)
        sym += " ";
      sym += to_string(fracs[2 * i]);
      if (fracs[2 * i + 1] > 1)
        sym += "/" + to_string(fracs[2 * i + 1]);
    }
    if (bar_pos == 3)
      sym += "|";
  }
  else {
    sym = "no symbol set";
  }
  return sym;
}

static string get_tri_symmetry(const vector<int> &fracs)
{
  vector<int> fs = fracs;
  vector<Vec3d> tmp(6);
  tri_normalise(fs, tmp);
  string sym;
  if (fs[2] == 2)
    sym = "D" + to_string(fs[4]);
  else if (fs[4] == 5)
    sym = "I";
  else if (fs[4] == 4)
    sym = "O";
  else if (fs[4] == 3)
    sym = "T";
  return sym;
}

string Wythoff::get_tri_sym() { return get_tri_symmetry(fracs); }

Status Wythoff::read_symbol(const char *sym)
{
  fracs = vector<int>(6);
  bar_pos = -1; // initialise to failure value

  // remove double spaces and spaces at beginning and end
  string sym_norm;
  int sym_len = strlen(sym);
  bool ignore_if_space = true;
  for (int i = 0; i < sym_len; i++) {
    if (sym[i] == ' ') {
      if (ignore_if_space)
        continue;
      else
        ignore_if_space = true;
    }
    else
      ignore_if_space = false;
    sym_norm += sym[i];
  }

  if (sym_norm[sym_norm.size() - 1] == ' ')
    sym_norm.resize(sym_norm.size() - 1);

  // remove spaces either side of a punctuation mark
  // record space counts and position of bar
  int bar_cnt = 0;
  int bar_off = -1;
  int space_before_bar_cnt = 0;
  int space_after_bar_cnt = 0;
  bool last_char_was_bar = false;
  string sym_norm2;
  for (unsigned int i = 0; i < sym_norm.length(); i++) {
    if (sym_norm[i] == ' ' &&
        ((i > 0 && ispunct(sym_norm[i - 1])) ||
         (i < sym_norm.length() && ispunct(sym_norm[i + 1]))))
      continue;

    if (sym_norm[i] == ' ') {
      if (bar_cnt)
        space_after_bar_cnt++;
      else
        space_before_bar_cnt++;
    }

    if (sym_norm[i] == '|') {
      last_char_was_bar = true;
      bar_cnt++;
      bar_off = sym_norm2.size();
      if (bar_off)
        sym_norm2 += ' ';
    }
    else {
      last_char_was_bar = false;
      sym_norm2 += sym_norm[i];
    }
  }
  if (last_char_was_bar)
    sym_norm2.resize(sym_norm2.size() - 1);

  size_t off = strspn(sym_norm2.c_str(), "0123456789/| ");
  if (off != sym_norm2.size())
    return Status::error(
        msg_str("unrecognised character '%c' in symbol", sym_norm2[off]));

  if (bar_cnt == 0)
    return Status::error("no bar in symbol");
  else if (bar_cnt > 1)
    return Status::error("more than one bar in symbol");

  int bar_pstn;
  if (bar_off == 0)
    bar_pstn = 0;
  else if (bar_off == (int)sym_norm2.size())
    bar_pstn = 3;
  else if (space_after_bar_cnt)
    bar_pstn = 1;
  else
    bar_pstn = 2;

  int total_spaces = space_before_bar_cnt + space_after_bar_cnt;
  // was | converted to space between two fractions
  if (bar_pstn == 1 || bar_pstn == 2)
    total_spaces++;

  if (total_spaces < 2)
    return Status::error("less than three fractions in symbol");
  else if (total_spaces > 2)
    return Status::error("more than three fractions in symbol");

  Status stat;
  char *fracs_str = &sym_norm2[0]; // destructive, not using a copy
  char *frac_p = strtok(fracs_str, " ");
  for (int f = 0; f < 3; f++) {
    if (!frac_p)
      return Status::error("internal symbol parsing error");

    char *p = strchr(frac_p, '/');
    int denominator = 1;
    int numerator;
    if (p != nullptr) {
      *p++ = '\0';
      if (!(stat = read_int(p, &denominator)))
        return Status::error(
            msg_str("denominator of fraction %d: %s", f + 1, stat.c_msg()));
    }

    if (!(stat = read_int(frac_p, &numerator)))
      return Status::error(
          msg_str("numerator of fraction %d: %s", f + 1, stat.c_msg()));

    if (numerator < 2)
      return Status::error(msg_str("numerator of fraction %d: must be an "
                                   "integer 2 or greater",
                                   f + 1));

    if (denominator % numerator == 0)
      return Status::error(msg_str("denominator of fraction %d: cannot be "
                                   "a multiple of the numerator",
                                   f + 1));

    fracs[2 * f] = numerator;
    fracs[2 * f + 1] = denominator % numerator;

    frac_p = strtok(nullptr, " ");
  }

  bar_pos = bar_pstn; // clears failure value
  return Status::ok();
}

static bool get_tri_verts(const vector<int> &norm_fracs,
                          vector<Vec3d> &norm_verts)
{
  norm_verts.resize(6);
  bool found = false;
  if (norm_fracs[2] == 2) { // Dihedral
    found = true;
    norm_verts[0] = Vec3d::X;
    norm_verts[1] = Trans3d::rotate(Vec3d::Z, M_PI * (double)norm_fracs[5] /
                                                  (double)norm_fracs[4]) *
                    Vec3d::X;
    norm_verts[2] = Vec3d::Z;
  }
  else { // Check other triangles in Schwarz list
    for (int tri = 0; tri < num_schwarz_tris; tri++) {
      bool are_equal = true;
      for (int i = 0; i < 6; i++) {
        if (schwarz_triangles[tri][i] != norm_fracs[i]) {
          are_equal = false;
          break;
        }
      }
      if (are_equal) {
        found = true;
        for (int i = 0; i < 3; i++)
          for (int j = 0; j < 3; j++)
            norm_verts[i][j] = schwarz_triangles_verts[tri][i * 3 + j];

        break;
      }
    }
  }

  return found;
}

static bool assign_vertices(const vector<int> &fracs, vector<Vec3d> &verts)
{
  verts.resize(3);
  vector<Vec3d> v_map(3);
  for (int i = 0; i < 3; i++)
    v_map[i] = Vec3d((double)i, 0.0, 0.0);

  vector<int> norm_fracs = fracs;
  tri_normalise(norm_fracs, v_map);
  vector<Vec3d> norm_verts;
  bool ret = get_tri_verts(norm_fracs, norm_verts);
  if (ret) {
    for (int i = 0; i < 3; i++)
      verts[(int)floor(v_map[i][0] + 0.5)] = norm_verts[i];
  }

  return ret;
}

bool Wythoff::assign_verts() { return assign_vertices(fracs, verts); }

static Vec3d get_angle_bisector_norm(Vec3d v0, Vec3d v1, Vec3d v2)
{
  double ang = angle_around_axis(v1, v2, v0);
  return Trans3d::rotate(v0, ang / 2) * vcross(v0, v1);
}

static Vec3d get_fermat_point(Vec3d v0, Vec3d v1, Vec3d v2, bool degenerate,
                              string &msg)
{
  Vec3d v[3];
  v[0] = v0;
  v[1] = v1;
  v[2] = v2;
  Vec3d pt = (v0 + v1 + v2).unit(); // initialise to approx centroid
  // Use a fixed large number of iterations with small change. Degenerates are
  // sensitive and may produce different results with different params
  int iters = degenerate ? 50000 : 1000;
  double off_factor = degenerate ? 0.01 : 0.1;
  for (int n = 0; n < iters; n++) {
    Vec3d offset = Vec3d::zero;
    for (auto &i : v)
      offset += (i.component(pt) - i).unit();

    pt = (pt + off_factor * offset).unit();
  }

  msg.clear();
  double max_ang = 0.0;
  for (int i = 0; i < 3; i++) {
    double ang = angle_around_axis(v[i], v[(i + 1) % 3], pt);
    if (ang > M_PI)
      ang = 2 * M_PI - ang;
    double ang_diff = fabs(2 * M_PI / 3 - ang);
    max_ang = ang > max_ang ? ang_diff : max_ang;
  }
  if (max_ang > epsilon)
    msg =
        msg_str("inaccurate calculation of fermat point (angle difference %g)",
                max_ang);

  return pt;
}

static void add_faces(Geometry &geom, Vec3d pt, int num, int denom,
                      const Vec3d &axis, Color col, const Symmetry &sym)
{
  // avoid extra windings
  int gr_fact = gcd(num, denom);
  num /= gr_fact;
  denom /= gr_fact;

  double ang = 2 * M_PI * (double)denom / (double)num;

  int sides = num;
  Geometry face_geom, sym_face_geom;
  vector<int> face(sides);
  for (int i = 0; i < sides; i++) {
    face_geom.add_vert(Trans3d::rotate(axis, ang * i) * pt);
    face[i] = i;
  }

  if (sides > 2)
    face_geom.add_face(face, col);
  else
    face_geom.add_edge(face, col);
  sym_repeat(sym_face_geom, face_geom, sym);
  merge_coincident_elements(sym_face_geom, "vf", epsilon);
  geom.append(sym_face_geom);
}

static void add_faces(Geometry &geom, Vec3d pt, int num, int denom,
                      const vector<Vec3d> &axes, int idx, const Symmetry &sym)
{
  add_faces(geom, pt, num, denom, axes[idx], Color(idx), sym);
}

Status Wythoff::make_poly(Geometry &geom)
{
  geom.clear_all();
  string message;
  Symmetry sym(get_tri_sym());
  if (bar_pos == 0) {
    int max_fract = 0;
    for (int i = 0; i < 3; i++) { // find smallest fraction (largest angle)
      if (double(fracs[2 * i]) / fracs[2 * i + 1] <=
          double(fracs[2 * max_fract]) / fracs[2 * max_fract + 1])
        max_fract = i;
    }
    // Check for invalid antiprism
    if (2 * fracs[2 * max_fract] < 3 * fracs[2 * max_fract + 1] &&
        fracs[(2 * max_fract + 2) % 6] == 2 &&
        fracs[(2 * max_fract + 4) % 6] == 2) {
      return Status::error("symbol leads to nonconstructible antiprism");
    }

    Vec3d f_pt; // Fermat point
    Vec3d pt;   // Final construction point

    bool degenerate = false;

    // triangles with a 3/2 vertex are generally a problem, if there
    // is only one of these then the solution (generally) lies at this
    // vertex and requires special processing the case of
    int cnt_3_2 = 0;
    int pos_3_2 = 0;
    for (int i = 0; i < 3; i++)
      if (fracs[2 * i] == 3 && fracs[2 * i + 1] == 2) {
        cnt_3_2++;
        pos_3_2 = i;
      }
    if (cnt_3_2 == 1)
      degenerate = true;

    // first check for non-dihedral isoscelese triangle with 3/2 apex
    if (cnt_3_2 == 1 &&
        fracs[(2 * pos_3_2 + 2) % 6] == fracs[(2 * pos_3_2 + 4) % 6] &&
        fracs[(2 * pos_3_2 + 3) % 6] == fracs[(2 * pos_3_2 + 5) % 6] &&
        fracs[(2 * pos_3_2 + 2) % 6] != 2) {
      // |3/2 5/3 5/3 or |3/2 5/4 5/4 have a different construction point
      if (fracs[(2 * pos_3_2 + 2) % 6] == 5 &&
          (fracs[(2 * pos_3_2 + 3) % 6] == 3 ||
           fracs[(2 * pos_3_2 + 3) % 6] == 4)) {
        pt = Trans3d::reflection(
                 vcross(verts[(pos_3_2 + 2) % 3], verts[pos_3_2])) *
             verts[(pos_3_2 + 1) % 3];
      }
      else { // take apex as Fermat point and use smallest circumcentre
        pt = verts[(pos_3_2 + 1) % 3] + verts[(pos_3_2 + 2) % 3];
      }
    }
    else { // general case
      f_pt =
          get_fermat_point(verts[0], verts[1], verts[2], degenerate, message);

      // Reflect in sides of triangle
      Vec3d u0 = Trans3d::reflection(vcross(verts[1], verts[2])) * f_pt;
      Vec3d u1 = Trans3d::reflection(vcross(verts[2], verts[0])) * f_pt;
      Vec3d u2 = Trans3d::reflection(vcross(verts[0], verts[1])) * f_pt;
      pt = vcross(u0 - u1, u1 - u2); // circumcentre
    }

    pt.to_unit();

    add_faces(geom, pt, fracs[0], fracs[1], verts, 0, sym);
    add_faces(geom, pt, fracs[2], fracs[3], verts, 1, sym);
    add_faces(geom, pt, fracs[4], fracs[5], verts, 2, sym);

    // Add snub triangle faces
    int dir = 1 - 2 * (vtriple(verts[0], verts[1], verts[2]) > 0);
    Vec3d tri_cent = pt;
    tri_cent += Trans3d::rotate(verts[0], dir * 2 * M_PI * (double)fracs[1] /
                                              (double)fracs[0]) *
                pt;
    tri_cent += Trans3d::rotate(verts[1], -dir * 2 * M_PI * (double)fracs[3] /
                                              (double)fracs[2]) *
                pt;
    add_faces(geom, pt, 3, 2, tri_cent, Color(3), sym);
  }
  else if (bar_pos == 1) {
    Vec3d pt = verts[0];
    if (fracs[2] == 2 && fracs[4] == 2) { // P|2 2 is degenerate
      geom.add_vert(pt);
      geom.add_vert(-pt);
      geom.add_face({0, 1}); // add as face, for sizing by edge length
    }
    else { // usual construction
      add_faces(geom, pt, fracs[2], fracs[3], verts, 1, sym);
      add_faces(geom, pt, fracs[4], fracs[5], verts, 2, sym);
    }
  }
  else if (bar_pos == 2) {
    Vec3d n0 = get_angle_bisector_norm(verts[2], verts[0], verts[1]);
    Vec3d n1 = vcross(verts[0], verts[1]);
    Vec3d pt = vcross(n0, n1).unit();

    add_faces(geom, pt, fracs[0], fracs[1], verts, 0, sym);
    add_faces(geom, pt, fracs[2], fracs[3], verts, 1, sym);
    // All hemis apart from 3/2 3 | 3 have duplicated faces
    merge_coincident_elements(geom, "vf", epsilon);
    add_faces(geom, pt, 2 * fracs[4], fracs[5], verts, 2, sym);
  }
  else if (bar_pos == 3) {
    Vec3d n0 = get_angle_bisector_norm(verts[1], verts[2], verts[0]);
    Vec3d n1 = get_angle_bisector_norm(verts[2], verts[0], verts[1]);
    Vec3d pt = vcross(n0, n1).unit();

    add_faces(geom, pt, 2 * fracs[0], fracs[1], verts, 0, sym);
    add_faces(geom, pt, 2 * fracs[2], fracs[3], verts, 1, sym);
    add_faces(geom, pt, 2 * fracs[4], fracs[5], verts, 2, sym);
  }
  else
    return Status::error("invalid symbol");

  merge_coincident_elements(geom, "v", epsilon);
  return (message.empty()) ? Status::ok() : Status::warning(message);
}

bool Wythoff::make_tri(Geometry &geom)
{
  geom.clear_all();
  if (is_set()) {
    geom.add_verts(verts);
    geom.add_face({0, 1, 2});
  }
  return is_set();
}

bool Wythoff::make_tri_poly(Geometry &geom)
{
  geom.clear_all();
  if (is_set()) {
    Symmetry sym(get_tri_sym());

    if (sym.get_sym_type() == Symmetry::D) {
      // N/D with D even is double wrapped surface and cannot be merged. Use
      // specific construction rather than symmetry repeat with merge
      int N, D;
      for (int i = 0; i < 3; i++) {
        N = fracs[2 * i];
        D = fracs[2 * i + 1];
        if (N != 2 || D != 1)
          break;
      };
      geom.clear_all();
      geom.add_vert(Vec3d::Z);
      geom.add_vert(-Vec3d::Z);
      for (int i = 0; i < 2 * N; i++) {
        double ang = i * M_PI * D / N;
        geom.add_vert(Vec3d(cos(ang), sin(ang), 0));
        geom.add_face({2 + i, 2 + (i + 1) % (2 * N), 0}, Color(i % 2));
        geom.add_face({1, 2 + (i + 1) % (2 * N), 2 + i}, Color((i + 1) % 2));
      }
    }
    else {
      Geometry tri;
      make_tri(tri);
      Geometry geom_tri;
      sym_repeat(geom_tri, tri, sym);
      Coloring clrng(&geom_tri);
      clrng.f_one_col(0);
      geom.append(geom_tri);
      Vec3d norm =
          (sym.get_sym_type() == Symmetry::T) ? Vec3d(1, 1, 0) : Vec3d::Z;
      geom_tri.transform(Trans3d::reflection(norm));
      clrng.f_one_col(1);
      geom.append(geom_tri);
      merge_coincident_elements(geom, "v", epsilon);
    }
  }
  return is_set();
}

//------------------------------------------------------------
// General Wythoff tiling
static ElemProps<Color> get_original_colors(const Geometry &geom, bool is_meta)
{
  // Get the vertex colours first
  ElemProps<Color> orig_colors;
  for (unsigned int i = 0; i < geom.verts().size(); i++) {
    const auto &col = geom.colors(VERTS).get(i);
    if (col.is_set())
      orig_colors.set(i, col);
  }

  // for meta tiling, this is all the colours, for a polyhedron
  // base, add the face colours and then the edge colours.
  if (!is_meta) {
    int f_start = geom.verts().size(); // offset for face index numbers
    for (unsigned int i = 0; i < geom.faces().size(); i++) {
      const auto &col = geom.colors(FACES).get(i);
      if (col.is_set())
        orig_colors.set(i + f_start, col);
    }

    map<vector<int>, Color> e2col;
    for (unsigned int i = 0; i < geom.edges().size(); i++) {
      Color col = geom.colors(EDGES).get(i);
      if (col.is_set())
        e2col[geom.edges(i)] = col;
    }

    // offset for edge index numbers (where index is position in
    // implicicit edge list)
    int e_start = geom.verts().size() + geom.faces().size();
    GeometryInfo info(geom);
    map<vector<int>, int> e2v;
    int e_idx = 0;
    for (const auto &e : info.get_impl_edges()) {
      auto e_it = e2col.find(e);
      if (e_it != e2col.end())
        orig_colors.set(e_idx + e_start, e_it->second);
      e_idx++;
    }
  }
  return orig_colors;
}

static void make_meta(const Geometry &geom, Geometry &meta,
                      double face_ht = 0.0)
{
  meta.clear_all();
  for (const auto &vert : geom.verts())
    meta.add_vert(vert, Color(0));
  int f_start = meta.verts().size();
  for (unsigned int f = 0; f < geom.faces().size(); f++) {
    Vec3d face_pt = geom.face_cent(f);
    if (face_ht)
      face_pt += geom.face_norm(f).with_len(face_ht);
    meta.add_vert(face_pt, Color(2));
  }

  GeometryInfo info(geom);
  map<vector<int>, int> e2v;
  for (const auto &e : info.get_impl_edges())
    e2v[e] = meta.add_vert(geom.edge_cent(e), Color(1));
  for (int f_idx = 0; f_idx < (int)geom.faces().size(); f_idx++) {
    int f_cent_idx = f_start + f_idx;
    for (int v = 0; v < (int)geom.faces(f_idx).size(); v++) {
      int v0 = geom.faces(f_idx, v);
      int v1 = geom.faces_mod(f_idx, v + 1);
      int e_cent_idx = e2v[make_edge(v0, v1)];
      meta.add_face({v0, e_cent_idx, f_cent_idx});
      meta.add_face({v1, e_cent_idx, f_cent_idx});
    }
  }
}

static Status normalize_tri(Geometry &geom, int f_idx, int v0, int v1,
                            Color other_v_col)
{
  vector<int> &face = geom.faces(f_idx);
  if (face.size() != 3)
    return Status::error(msg_str("face %d is not a triangle", f_idx));
  bool found = false;
  for (int i = 0; i < 3; i++) {
    if (face[i] == v0 && face[(i + 1) % face.size()] == v1) {
      found = true;
      break;
    }
  }
  if (!found)
    std::reverse(face.begin(), face.end());

  int other_v_idx;
  for (int i = 0; i < 3; i++) {
    other_v_idx = face[i];
    if (other_v_idx != v0 && other_v_idx != v1)
      break;
  }

  Color this_other_v_col = geom.colors(VERTS).get(other_v_idx);
  if (this_other_v_col.is_set() && this_other_v_col != other_v_col)
    return Status::error("vertices cannot be 3-coloured");
  else
    geom.colors(VERTS).set(other_v_idx, other_v_col);

  for (int i = 0; i < 3; i++)
    if (geom.colors(VERTS).get(face[i]) == Color(0))
      std::rotate(face.begin(), face.begin() + i, face.end());

  return Status::ok();
}

static Status normalize_meta(Geometry &geom)
{
  geom.clear_cols();
  if (geom.faces().size() == 0 || geom.faces().size() % 2)
    return Status::error(
        msg_str("geometry does not have an even number of faces"));

  // int part_num = 0;
  const int done = -1;
  auto edges = geom.get_edge_face_pairs(false);
  vector<int> cur_idx(geom.faces().size(), 0);
  vector<int> prev_face(geom.faces().size(), 0);
  vector<int> e_verts(2);
  vector<int> orig_e_verts(2);
  vector<int> e_faces(2);
  for (unsigned int i = 0; i < geom.faces().size(); i++) {
    if (geom.faces(i).size() != 3)
      return Status::error(msg_str("face %d is not a triangle", i));
    if (cur_idx[i] == done)
      continue;

    // first face in part has colour 0, and original oriented is preserved
    // and vertices are in order VEF. Acts as seed for all other faces
    geom.colors(FACES).set(i, 0);
    geom.colors(VERTS).set(geom.faces_mod(i, 0), Color(0)); // V verts colour 0
    geom.colors(VERTS).set(geom.faces_mod(i, 1), Color(1)); // E verts colour 1
    geom.colors(VERTS).set(geom.faces_mod(i, 2), Color(2)); // F verts colour 2

    int cur_fidx = i;
    // if (parts)
    //  parts->push_back(vector<int>(1, i));
    while (cur_idx[i] != done) {
      int idx = cur_idx[cur_fidx];
      if (idx == done) {
        cur_fidx = prev_face[cur_fidx];
        continue;
      }

      // read off the next edge
      const vector<int> &face = geom.faces(cur_fidx);
      orig_e_verts[0] = face[idx];
      idx = (idx + 1) % face.size();
      orig_e_verts[1] = face[idx];
      cur_idx[cur_fidx] = idx ? idx : done; // set to next idx, or mark done

      e_verts = orig_e_verts;
      if (e_verts[0] > e_verts[1])
        swap(e_verts[0], e_verts[1]);
      e_faces = edges.find(e_verts)->second;

      int next_face = (e_faces[0] != cur_fidx) ? e_faces[0] : e_faces[1];
      if (next_face >= 0 && cur_idx[next_face] == 0) { // face not looked at yet
        Color cur_col = geom.colors(FACES).get(cur_fidx);
        // Adjacent faces must be coloured differently
        if (geom.colors(FACES).get(next_face) == cur_col)
          return Status::error("faces cannot be 2-coloured");
        else
          geom.colors(FACES).set(next_face, Color(!cur_col.get_index()));

        int other_v_idx = (idx + 1) % face.size();
        Color other_v_col =
            geom.colors(VERTS).get(geom.faces(cur_fidx, other_v_idx));
        normalize_tri(geom, next_face, orig_e_verts[0], orig_e_verts[1],
                      other_v_col);
        // if (parts)
        //  parts->back().push_back(next_face);
        prev_face[next_face] = cur_fidx;
        cur_fidx = next_face;
      }
    }
    // part_num++;
  }

  // Reorder faces to have colours 0,1,0,1,... will recolour by order later
  vector<int> bad[2];
  for (int i = 0; i < (int)geom.faces().size(); i++)
    if (geom.colors(FACES).get(i).get_index() != i % 2)
      bad[i % 2].push_back(i);
  for (int i = 0; i < (int)bad[0].size(); i++)
    swap(geom.raw_faces()[bad[0][i]], geom.raw_faces()[bad[1][i]]);

  return Status::ok();
}

void Tile::start_op() const
{
  ops_i = ops.begin();
  idxs_i = idxs.begin();
}

void Tile::next_op() const
{
  if (ops_i != ops.end()) {
    ++ops_i;
    if (*ops_i == P)
      ++idxs_i;
  }
}

int Tile::get_op() const
{
  if (ops_i == ops.end())
    return END;
  else
    return *ops_i;
}

int Tile::get_idx() const { return *idxs_i; }

// https://codereview.stackexchange.com/questions/187212/remove-all-adjacent-duplicates-in-a-string-using-a-stack
template <typename Iter> Iter remove_duplicates(Iter begin, Iter end)
{
  using namespace std::placeholders;

  auto dest = begin;

  do {
    Iter pair = std::adjacent_find(begin, end);
    if (dest != begin) {
      std::copy(begin, pair, dest);
      dest += std::distance(begin, pair);
    }
    else {
      dest = pair;
    }
    begin = pair;
    if (pair != end) {
      begin += 2;
    }
  } while (begin != end);

  return dest;
}

template <typename Iter> Iter repeatedly_remove_duplicates(Iter begin, Iter end)
{
  Iter new_end;
  while ((new_end = remove_duplicates(begin, end)) != end) {
    end = new_end;
  }
  return end;
}

std::string repeatedly_remove_duplicates(std::string s)
{
  s.erase(repeatedly_remove_duplicates(s.begin(), s.end()), s.end());
  return s;
}

Tile::TileReport Tile::get_element_association() const
{
  TileReport rep;
  string elems = "vef";
  string association = "VEFX";
  string ops_str;
  for (auto op : ops)
    if (op != P)
      ops_str.push_back(elems[op]);

  auto reduced = repeatedly_remove_duplicates(ops_str);
  int sz = reduced.size();
  int mismatch_idx;
  for (mismatch_idx = 0; mismatch_idx < sz; mismatch_idx++)
    if (reduced[mismatch_idx] != reduced[sz - 1 - mismatch_idx])
      break;

  rep.step = reduced.substr(0, mismatch_idx);
  rep.assoc = reduced.substr(mismatch_idx, sz - 2 * mismatch_idx);
  rep.step_back = reduced.substr(sz - mismatch_idx, sz);

  vector<bool> contains(3, false); // contains which of v, e, f
  for (int i = 0; i < 3; i++)
    if (rep.assoc.find(elems[i]) != string::npos)
      contains[i] = true;

  int elem_tri_idx;
  if (contains[0] && contains[1] && contains[2]) // v, e and f
    elem_tri_idx = VEF;
  else if (contains[0] && contains[1]) // v and e
    elem_tri_idx = F;
  else if (contains[1] && contains[2]) // e and f
    elem_tri_idx = V;
  else if (contains[2] && contains[0]) // f and v
    elem_tri_idx = E;
  else if (contains[0]) // v
    elem_tri_idx = P;   // triangle-like
  else if (contains[1]) // e
    elem_tri_idx = P;   // triangle-like
  else if (contains[2]) // f
    elem_tri_idx = P;   // triangle-like
  else                  // empty
    elem_tri_idx = P;   // assign to triangle

  rep.assoc_type = elem_tri_idx;

  // string assoc_elem_str = "VEF345X";   // VEF=6 -> X
  // fprintf(stderr, "final=(%s)%s(%s) association=%c\n\n", rep.step.c_str(),
  //        rep.assoc.c_str(), rep.step_back.c_str(),
  //        assoc_elem_str[elem_tri_idx]);
  return rep;
}

Status Tile::read(const string &pat)
{
  // initialise
  ops.clear();
  idxs.clear();
  bool has_tris_spec = strchr("+-*", pat[0]);
  start_faces = (has_tris_spec) ? pat[0] : '+';
  unsigned int pos = has_tris_spec;
  if (!std::isdigit(pat[pos]) && !std::isdigit(pat.back()))
    return Status::error(
        "tile specifier: first character (or first character after +-*), "
        "or last character must be a digit");

  while (pos < pat.size()) {
    int len;
    // point
    if ((len = strspn(pat.substr(pos).c_str(), "0123456789"))) {
      ops.push_back(P);
      int idx;
      read_int(pat.substr(pos, len).c_str(), &idx);
      idxs.push_back(idx);
    }

    // mirrors
    else if ('v' == pat[pos])
      ops.push_back(V);
    else if ('e' == pat[pos])
      ops.push_back(E);
    else if ('f' == pat[pos])
      ops.push_back(F);

    // rotations
    else if ('V' == pat[pos]) {
      ops.push_back(E);
      ops.push_back(F);
    }
    else if ('E' == pat[pos]) {
      ops.push_back(F);
      ops.push_back(V);
    }
    else if ('F' == pat[pos]) {
      ops.push_back(V);
      ops.push_back(E);
    }

    // no op - stay on same triangle
    else if ('_' == pat[pos])
      ;
    else {
      return Status::error(
          msg_str("invalid character '%c' in position %d", pat[pos], pos + 1));
    }

    if (len)
      pos += len;
    else
      pos++;
  }

  get_element_association();
  return Status::ok();
}

void Tile::relabel(vector<int> relab)
{
  for (auto &op : ops)
    if (0 <= op && op < 3)
      op = relab[op];
}

void Tile::flip_start_faces()
{
  // Flip +/-, * is left unchanged
  if (start_faces == '+')
    start_faces = '-';
  else if (start_faces == '-')
    start_faces = '+';
}

vector<int> Tile::check_index_range(int num_points) const
{
  vector<int> out_of_range;
  for (auto &idx : idxs)
    if (idx < 0 || idx >= num_points)
      out_of_range.push_back(idx);
  return out_of_range;
}

string Tile::tile_string()
{
  vector<int> op2;
  string VEF = "VEF";
  string vef = "vef";
  string tile;
  if (start_faces != '+')
    tile += start_faces;
  unsigned int p_idx = 0;
  int last_op = -1;
  for (auto op : ops) {
    if (op == P) {
      if (op == last_op)
        tile += "_";
      if (p_idx < idxs.size())
        tile += to_string(idxs[p_idx]);
      else {
        tile += msg_str("ERROR: index %u out of range", p_idx);
        break;
      }
      p_idx++;
    }
    else
      tile += vef[op];
    last_op = op;
  }

  map<char, int> elem_idx = {{'v', 0}, {'e', 1}, {'f', 2}};
  string tile2;
  for (unsigned int i = 0; i < tile.size(); i++) {
    // convert pairs of consecutive letters from vef to VEF
    if (i < tile.size() - 1 && strchr("vef", tile[i]) &&
        strchr("vef", tile[i + 1]) &&
        (elem_idx[tile[i]] + 1) % 3 == elem_idx[tile[i + 1]]) {
      tile2.push_back(VEF[(elem_idx[tile[i]] + 2) % 3]);
      i++; // skip second letter of pair
    }
    else
      tile2.push_back(tile[i]);
  }

  return tile2;
}

struct ConwayOperator {
  std::string operator_short;
  std::string operator_name;
  std::string pattern;
  int num_params;
};

// clang-format off
ConwayOperator conway_operator_list[]{
    // Equivalent: d, a, s
    {"d",   "dual",           "[F]0V,0E", 0},
    {"a",   "ambo",           "[E]0F,0V", 0},
    {"S",   "seed",           "[V]0E,0F", 0},

    {"j",   "join",           "[F,V]0_1E", 0},

    // Equivalent: k, n, u
    {"k",   "kis",            "[F,V]0_1v1v,1E", 0},
    {"n",   "needle",         "[V,F]0_1f1f,1E", 0},
    {"u",   "subdivide",      "[V,E]1F,0_1e1e", 2},

    // Equivalent: t, z, e (tile order to match e0=z and e1=e
    {"t",   "truncate",       "[VE]0V0E,0V,0E", 0},
    {"z",   "zip",            "[EF]0E0F,0F,0E", 0},
    {"e",   "expand",         "[FV]0V,0F,0F0V", 2},

    // Symmetric: s, m, b
    {"s",   "snub",           "[VEF]0V,0E,0F,0V0E0F", 1},
    {"m",   "meta",           "[V,E,F]*0_1_2", 1},
    {"b",   "bevel",          "[VEF]0e0f,0v0e,0f0v", 1},

    {"o",   "ortho",          "[V,E,F]1_0e1_2e", 2},
    {"g",   "gyro",           "[F,VE,V]1_0F1_2V1E,1E", 1},
    {"c",   "chamfer",        "[V,VF]1F,0_1v1f", 0},
    {"l",   "loft",           "[V,VF]1F,0_1v1_0v,0E", 1},
    {"p",   "propellor",      "[V,VEF]1F,1_0V1E1F,1E", 0},
    {"q",   "quinto",         "[V,E,EF]2F,0_1_2e2_1e", 0},
    {"L_0", "joined-lace",    "[V,EF]1F,1e1_0e,1_0E", 0},
    {"G",   "opposite-lace",  "[V,EF]1F,1e1_0e,0_1f1f,1E", 0},
    {"L",   "lace",           "[V,EF]1F,1e1_0e,1_0v0v,0E", 1},
    {"K",   "stake",          "[V,EF,F]0_1_2e1e,1_0v0v,0E", 0},
    {"M",   "edge-medial",    "[F,3V,V2E]0_2_1e2e,2_0v2v,2E", 1},
    {"J",   "joined-medial",  "[F,V,EF]*0_1_2,1_2E", 0},
    {"X",   "cross",          "[V,E,F,VF]3_1v3_2v,*0_1_3", 0},
    {"w",   "whirl",          "[VF,VE,V]0F,0_1V2_1E1_0F,1E", 0},
    {"E",   "ethel",          "[V,VE,VF]0_1_2e1e,2F,1_2v2f", 0},
    {"W",   "waffle",         "[V,E,F,V2E,VF]0_4_3f4f,2_4_3v3_4v,3E", 0},
    {"B",   "bowtie",         "[V,E,F,VE,EF]1_3_4,0_3_4_2e4_1_3e", 0},
};
// clang-format on

Status TilingColoring::set_type_associated(int elem, Color col)
{
  col_type = ColoringType::associated_element;
  if (elem == Tile::V || elem == Tile::E || elem == Tile::F ||
      elem == Tile::P) {
    local_tile_assoc = elem;
    color = col;
  }
  else
    return Status::error("invalid associated element type");

  return Status::ok();
}

Status TilingColoring::read(const string &str)
{
  Status stat;
  Split parts(str.c_str(), ",");
  if (parts.size() > 2)
    return stat.set_error("more than one comma given");

  string arg_id;
  if (!(stat = ProgramOpts::get_arg_id(parts[0], &arg_id,
                                       "none=0|index=1|value=2|association=3")))
    return stat;

  int id = atoi(arg_id.c_str());
  if (id != 3 && parts.size() > 1)
    return stat.set_error(
        "comma given, but only valid for associated colouring");

  if (id == 0)
    set_type_none();
  else if (id == 1 || id == 2)
    set_type_path_index((id == 2));
  else {
    if (parts.size() > 1) {
      if (strcasecmp(parts[1], "V") == 0)
        set_type_associated(Tile::V);
      else if (strcasecmp(parts[1], "E") == 0)
        set_type_associated(Tile::E);
      else if (strcasecmp(parts[1], "F") == 0)
        set_type_associated(Tile::F);
      else {
        Color col;
        if (!(stat = col.read(parts[1])))
          return stat;
        set_type_associated(Tile::P, col);
      }
    }
    else
      set_type_associated(Tile::F); // Default association for local tiles
  }

  return stat;
}

string TilingColoring::get_option_help(char op_char)
{
  return msg_str(R"(  -%c <mthd> colouring method for tiles, method can be
               none - do not colour tiles
               index - use the path index, colour with index numbers
               value - use the path index, colour with colour values (default)
               association - colour tiles using corresponding base element
                 colour, optionally followed by a comma and a letter from V,
                 E, F, or a colour, to colour local tiles by that element type
                 (default: F), or colour all local tiles with a single colour)",
                 op_char);
}

bool Tiling::find_nbrs()
{
  auto ef_pairs = meta.get_edge_face_pairs(false);

  // Find the neighbour face opposite each VEF vertex
  nbrs.resize(meta.faces().size(), vector<int>(3));
  for (unsigned int f = 0; f < meta.faces().size(); f++)
    for (int i = 0; i < 3; i++) {
      vector<int> e(2);
      e[0] = meta.faces_mod(f, i + 1);
      e[1] = meta.faces_mod(f, i + 2);
      if (e[0] > e[1])
        swap(e[0], e[1]);
      auto ef_i = ef_pairs.find(e);
      if (ef_i == ef_pairs.end())
        return false;
      else if (ef_i->second.size() != 2)
        nbrs[f][i] = -1; // only allow connection for two faces at an edge
      else {
        nbrs[f][i] =
            (ef_i->second[0] != (int)f) ? ef_i->second[0] : ef_i->second[1];
      }
    }
  return true;
}

static Vec3d point_on_face(const Geometry &meta, int f_idx, const Vec3d &crds)
{
  // point coordinates
  Vec3d P = crds[Tile::V] * meta.face_v(f_idx, Tile::V) +
            crds[Tile::E] * meta.face_v(f_idx, Tile::E) +
            crds[Tile::F] * meta.face_v(f_idx, Tile::F);

  // point normal
  // ret[1] =
  //    crds[Tile::V] * vert_norms[meta.faces(f_idx, Tile::V)] +
  //    crds[Tile::E] * vert_norms[meta.faces(f_idx, Tile::E)] +
  //    crds[Tile::F] * vert_norms[meta.faces(f_idx, Tile::F)];
  // ret[1].to_unit();
  return P;
}

Color Tiling::get_associated_element_point_color(int f_idx, int incl) const
{
  int idx = -1;
  if (incl == Tile::V) // on a vertex
    idx = Tile::V;
  else if (incl == Tile::E) // on an edge
    idx = Tile::E;
  else // on a face
    idx = Tile::F;

  return orig_colors.get(meta.faces(f_idx, idx));
}

int Tiling::get_associated_element(int start_idx, const string &step,
                                   int assoc_type) const
{
  const map<char, int> elem_idx = {{'V', 0}, {'E', 1}, {'F', 2}};
  int idx;
  if (assoc_type == Tile::VEF)
    idx = -1; // invalid index
  else {
    idx = start_idx;
    for (char op : step)
      idx = nbrs[idx][elem_idx.at(std::toupper(op))]; // move to next tri
  }
  return idx >= 0 ? meta.faces(idx, assoc_type) : idx;
}

// Each pattern point plotted for a meta triangle cooresponds to
// a previously assigned geometry vertex. Get the index of that vertex.
static int
get_index(const vector<int> &face, int f_idx, int pat_pt_idx, int incl,
          const std::vector<std::map<std::vector<int>, std::pair<int, int>>>
              &index_order,
          const std::vector<int> &point_vertex_offsets)
{
  vector<int> incl_key;
  if (incl >= Tile::V && incl <= Tile::F) // meta tile vertex
    incl_key = {face[incl]};
  else if (incl >= Tile::VE && incl <= Tile::FV) // meta tile edge
    incl_key = make_edge(face[incl % 3], face[(incl + 1) % 3]);
  else if (incl == Tile::VEF) // meta tile interior
    incl_key = {f_idx};
  else
    return -1; // invalid inclusion value, shouldn't happen

  const auto it = index_order[incl].find(incl_key);
  if (it == index_order[incl].end())
    return -2; // invalid element key, shouldn't happen

  return point_vertex_offsets[pat_pt_idx] + it->second.first;
}

void Tiling::add_circuit(
    Geometry &geom, int start_idx, const Tile &pat, std::vector<bool> &seen,
    Color col,
    const std::vector<std::map<std::vector<int>, std::pair<int, int>>>
        &index_order,
    const std::vector<int> &point_vertex_offsets) const
{
  // Apply pattern until circuit completes
  vector<int> face;
  int idx = start_idx;
  while (true) {
    seen[idx] = true;
    pat.start_op();
    while (pat.get_op() != Tile::END) {
      if (pat.get_op() == Tile::P) {
        int incl = points[pat.get_idx()].second.get_index();
        int v_idx = get_index(meta.faces(idx), idx, pat.get_idx(), incl,
                              index_order, point_vertex_offsets);
        face.push_back(v_idx);
      }
      else {
        idx = nbrs[idx][pat.get_op()]; // move to next triangle
        if (idx < 0)
          return; // abandon: circuit tried to cross an open edge
      }
      pat.next_op();
    }
    if (idx == start_idx) // circuit complete
      break;
  }

  geom.add_face(face, col);
}

static void reverse_odd_faces(Geometry &geom)
{
  const int f_sz = geom.faces().size();
  for (int i = 0; i < f_sz; i++)
    if (i % 2)
      reverse(geom.faces(i).begin(), geom.faces(i).end());
}

Status Tiling::set_geom(const Geometry &geom, bool is_meta, double face_ht)
{
  orig_colors = get_original_colors(geom, is_meta);

  if (is_meta) {
    meta = geom;
    Status stat = normalize_meta(meta);
    if (stat.is_error())
      return stat;
  }
  else
    make_meta(geom, meta, face_ht);

  find_nbrs();
  if (is_meta) {
    // Neighbouring faces must have index numbers of opposite parity
    for (int i = 0; i < (int)nbrs.size(); i++)
      for (int j = 0; j < 3; j++)
        if (i % 2 == nbrs[i][j] % 2)
          return Status::error("faces cannot be 2-coloured");
  }

  reverse_odd_faces(meta);
  // vert_norms = meta.get_info().get_vert_norms();
  reverse_odd_faces(meta);
  return Status::ok();
}

Status Tiling::add_tile(const string &pat)
{
  Tile pattern;
  Status stat = pattern.read(pat);
  if (stat)
    pat_paths.push_back(pattern);

  return stat;
}

void Tiling::reverse_pattern()
{
  for (auto &path : pat_paths)
    path.flip_start_faces();
}

void Tiling::start_everywhere()
{
  for (auto &path : pat_paths)
    path.set_start_faces('*');
}

static void delete_verts(Geometry &geom, const vector<int> &v_nos)
{
  vector<int> dels = v_nos;
  map<int, int> v_map;
  if (!dels.size())
    return;
  sort(dels.begin(), dels.end());
  unsigned int del_verts_cnt = 0;
  int map_to;
  for (unsigned int i = 0; i < geom.verts().size(); i++) {
    if (del_verts_cnt < dels.size() && (int)i == dels[del_verts_cnt]) {
      del_verts_cnt++;
      map_to = -1;
    }
    else {
      map_to = i - del_verts_cnt;
      geom.verts(map_to) = geom.verts(i);
    }
    v_map[i] = map_to;
  }
  geom.raw_verts().resize(geom.verts().size() - del_verts_cnt);
  geom.colors(VERTS).remap(v_map);

  vector<int> del_faces;
  for (unsigned int i = 0; i < geom.faces().size(); i++) {
    int curr_idx = 0;
    for (unsigned int j = 0; j < geom.faces(i).size(); j++) {
      int new_idx = v_map[geom.faces(i, j)];
      if (new_idx >= 0)
        geom.raw_faces()[i][curr_idx++] = new_idx;
    }
    geom.raw_faces()[i].resize(curr_idx);
    if (curr_idx < 2)
      del_faces.push_back(i);
  }
  geom.del(FACES, del_faces);
}

static bool valid_start_face(int f, int start_faces)
{
  int pos_tri = f % 2;
  return !((start_faces == '-' && pos_tri) || (start_faces == '+' && !pos_tri));
}

static void store_tri(map<vector<int>, pair<int, int>> &elem_to_tri,
                      vector<int> key, int tri_idx)
{
  auto pr = elem_to_tri.insert({key, {-1, tri_idx}});
  if (pr.second == false && is_even(pr.first->second.second) &&
      !is_even(tri_idx))
    pr.first->second.second = tri_idx;
}

Status Tiling::make_tiling(Geometry &geom, TilingColoring col_type,
                           vector<Tile::TileReport> *tile_reports) const
{
  geom.clear_all();
  if (tile_reports)
    tile_reports->resize(pat_paths.size());

  // All the possible element inclusion postions V, E, F, VE, EF, FV, VEF.
  // Each entry maps to order (to find index of corresponding point)
  // and example triangle (to generate coordinates of corresponding point)
  auto ef_pairs = meta.get_edge_face_pairs();
  vector<map<vector<int>, pair<int, int>>> index_order(7);
  for (int i = 0; i < (int)meta.faces().size(); i++) {
    const auto &face = meta.faces(i);
    index_order[Tile::VEF][{i}] = {-1, i};
    index_order[Tile::V][{face[Tile::V]}] = {-1, i};
    index_order[Tile::E][{face[Tile::E]}] = {-1, i};
    index_order[Tile::F][{face[Tile::F]}] = {-1, i};
    store_tri(index_order[Tile::VE], make_edge(face[Tile::V], face[Tile::E]),
              i);
    index_order[Tile::EF][make_edge(face[Tile::E], face[Tile::F])] = {-1, i};
    index_order[Tile::FV][make_edge(face[Tile::F], face[Tile::V])] = {-1, i};
  }
  for (int i = (int)Tile::V; i <= (int)Tile::VEF; i++) {
    int pos = 0;
    for (auto &m : index_order[i])
      m.second.first = pos++;
  }

  // Starting offset of vertices corresponding to each pattern point
  vector<int> point_vertex_offsets(points.size());
  for (int i = 0; i < (int)points.size(); i++) {
    point_vertex_offsets[i] = geom.verts().size();
    const auto &pt = points[i];
    int incl = pt.second.get_index();
    Vec3d crds = pt.first;
    crds /= crds[0] + crds[1] + crds[2];
    for (auto &m : index_order[incl]) {
      const int f_idx = m.second.second;
      Color col; // col_type==ColoringType::none
      if (col_type.is_path_index())
        col = pt.second; // Colour by element type inclusion in coords
      else if (col_type.is_associated_element())
        col = get_associated_element_point_color(f_idx, incl);
      geom.add_vert(point_on_face(meta, f_idx, crds), col);
    }
  }

  int faces_sz = meta.faces().size();
  for (unsigned int p_idx = 0; p_idx < pat_paths.size(); p_idx++) {
    const auto &pat = pat_paths[p_idx];
    // Check index range
    auto out_of_range = pat.check_index_range(points.size());
    if (out_of_range.size()) {
      string msg = "Path" + to_string(p_idx) + ": index numbers out of range:";
      for (auto idx : out_of_range)
        msg += " " + to_string(idx) + ",";
      msg.pop_back();
      return Status::error(msg.c_str());
    }

    auto assoc = pat.get_element_association();
    vector<bool> seen(faces_sz, false);
    int start_faces_sz = geom.faces().size();
    unsigned char start_faces = pat.get_start_faces();
    for (int i = 0; i < faces_sz; i++) {
      // for(int f=0; f<(int)seen.size(); f++)
      //  fprintf(stderr, "i=%d, seen[%d]=%d\n", i, f, (int)seen[f]);
      // fprintf(stderr, "\n");
      if (!seen[i] && valid_start_face(i, start_faces)) {
        Color col; // col_type==ColoringType::none
        if (col_type.is_path_index())
          col.set_index(p_idx);
        else if (col_type.is_associated_element()) {
          auto type = (assoc.assoc_type == Tile::P)
                          ? col_type.get_local_tile_assoc()
                          : assoc.assoc_type;
          if (type == Tile::P)
            col = col_type.get_color();
          else {
            int col_idx = get_associated_element(i, assoc.step, type);
            if (col_idx >= 0)
              col = orig_colors.get(col_idx);
          }
        }
        add_circuit(geom, i, pat, seen, col, index_order, point_vertex_offsets);
        if (one_of_each_tile)
          break;
      }
    }
    if (tile_reports) {
      assoc.count = geom.faces().size() - start_faces_sz;
      tile_reports->at(p_idx) = assoc;
    }
  }

  delete_verts(geom, geom.get_info().get_free_verts());
  return Status::ok();
}

static void color_point(std::pair<Vec3d, Color> &point)
{
  // val {1, 2, 4, 3, 6, 5, 7, 0} = v, e, f, ve, ef, fv, vef, 0
  int col_idx_map[] = {7, 0, 1, 3, 2, 5, 4, 6}; // v, e, f, ve, ef, fv, vef, 0
  const Vec3d &pt = point.first;
  int val = (pt[0] != 0) + (pt[1] != 0) * 2 + (pt[2] != 0) * 4;
  point.second.set_index(col_idx_map[val]);
}

Status read_point(const char *point_str, std::pair<Vec3d, Color> &point)
{
  Vec3d &coords = point.first;
  coords = Vec3d::zero;
  map<char, int> elem_idx = {{'V', 0}, {'E', 1}, {'F', 2}};
  string pt_string(point_str);
  std::regex re_coord("([-+]?([0-9]*\\.[0-9]+|[0-9]+))?[VEF]");
  std::sregex_token_iterator next(pt_string.begin(), pt_string.end(), re_coord,
                                  {-1, 0});
  std::sregex_token_iterator end;
  if (next == end)
    return Status::error("invalid coordinate string");

  vector<bool> seen(3, false);
  while (next != end) {
    if (next->str() != "")
      return Status::error("invalid characters in coordinates: " + next->str());
    next++;
    int idx = elem_idx[next->str().back()];
    if (seen[idx])
      return Status::error(
          msg_str("coordinates %c given more than once", next->str().back()));
    else
      seen[idx] = true;
    if (next->str().size() < 2)
      coords[idx] = 1;
    else {
      string coord_str = next->str();
      coord_str.pop_back();
      Status stat = read_double(coord_str.c_str(), &coords[idx]);
      if (stat.is_error())
        return stat;
    }
    next++;
  }

  if (coords.len() == 0)
    return Status::error("coordinates cannot all be zero");

  color_point(point);

  return Status::ok();
}

Status Tiling::read_pattern(const string &pat)
{

  string pat2 = pat;
  std::smatch m_all;
  std::regex r_all("^\\[(.*)\\](.*)");
  std::regex_match(pat2, m_all, r_all);
  if (m_all.size() < 3)
    return Status::error(
        msg_str("pattern '%s': not in form [Point0,Point1,...]Path0,Path1...",
                pat.c_str()));

  Split parts;
  int num_parts = parts.init(m_all[1].str(), ",");
  points.resize(num_parts);
  for (int i = 0; i < num_parts; i++) {
    auto stat = read_point(parts[i], points[i]);
    if (stat.is_error())
      return Status::error(msg_str("Point%d: ", i) + stat.msg());
  }

  num_parts = parts.init(m_all[2].str(), ",");
  pat_paths.resize(num_parts);
  for (int i = 0; i < num_parts; i++) {
    auto stat = pat_paths[i].read(parts[i]);
    if (stat.is_error())
      return Status::error(msg_str("Path%d: ", i) + stat.msg());
  }
  return Status::ok();
}

Status Tiling::relabel_pattern(string relabel)
{
  for (auto &c : relabel)
    c = toupper(c);

  if (strlen(relabel.c_str()) != 3 || !strchr(relabel.c_str(), 'V') ||
      !strchr(relabel.c_str(), 'E') || !strchr(relabel.c_str(), 'F'))
    return Status::error(
        "relabel string does not contain exactly three letters V, E and F");

  map<char, int> elem_idx = {{'V', 0}, {'E', 1}, {'F', 2}};
  vector<int> relab(3);
  for (int i = 0; i < 3; i++)
    relab[i] = elem_idx[relabel[i]];

  for (auto &pt : points) {
    Vec3d v = pt.first;
    for (int i = 0; i < 3; i++)
      pt.first[relab[i]] = v[i];
    color_point(pt);
  }

  for (auto &pat : pat_paths)
    pat.relabel(relab);

  return Status::ok();
}

// nomalize point, if integer coordinates remove any common factor
static Vec3d normalize_point(const Vec3d &Q)
{
  vector<long> coords(3);
  for (int i = 0; i < 3; i++) {
    coords[i] = std::lround(Q[i]);
    if (!double_eq(coords[i], Q[i]))
      return Vec3d(Q); // not integer coordinates
  }

  long div = gcd(coords[0], gcd(coords[1], coords[2]));
  if (div == 0) // will not happen for barycentric coordinates
    div = 1;

  return Vec3d(coords[0] / div, coords[1] / div, coords[2] / div);
}

static string coord_string(Vec3d v)
{
  string VEF = "VEF";
  string coords;
  for (int i = 0; i < 3; i++)
    if (v[i]) {
      if (v[i] != 1.0) {
        auto str = msg_str("%g", v[i]);
        if (str.find('e') != string::npos)
          str = msg_str("%f", v[i]);
        coords += str;
      }
      coords += VEF[i];
    }

  return coords;
}

static string M_pattern(int N)
{
  N += 2; // FIX CODE BELOW TO AVOID THIS

  string pat = "[F";
  for (int i = 0; i < N + 1; i += 2) {
    double e_num = i;
    double v_num = N - e_num;
    pat += "," + coord_string(normalize_point(Vec3d(v_num, e_num, 0)));
  }
  int last_idx = N / 2 + 1;

  pat += "]0_2_1e2e";

  for (int i = 2; i < last_idx; i++)
    pat += msg_str(",*0_%d_%d", i, i + 1);

  if (!is_even(N)) {
    pat += msg_str(",%d_0v%dv", last_idx, last_idx);
    pat += msg_str(",%dE", last_idx);
  }

  return pat;
}

static string m_pattern(int N)
{
  N += 1; // FIX CODE BELOW TO AVOID THIS

  string pat = "[F";
  for (int i = 0; i < N + 1; i += 2) {
    double e_num = i;
    double v_num = N - e_num;
    pat += "," + coord_string(normalize_point(Vec3d(v_num, e_num, 0)));
  }
  int last_idx = N / 2 + 1;

  pat += "]";

  for (int i = 1; i < last_idx; i++)
    pat += msg_str(",*0_%d_%d", i, i + 1);

  if (!is_even(N)) {
    pat += msg_str(",%d_0v%dv", last_idx, last_idx);
    pat += msg_str(",%dE", last_idx);
  }

  return pat;
}

static string b_pattern(int N)
{
  N += 1; // FIX CODE BELOW TO AVOID THIS

  string pat = "[";
  for (int b = 1; b <= N + N % 2; b += 2)
    pat += coord_string(Vec3d((N - b), b, 1)) + ",";
  pat.back() = ']';

  pat += "0e0f,";

  for (int b = 0; b < N + N % 2; b += 2)
    pat += msg_str("%d_", b / 2);
  pat.back() = 'v';
  for (int b = 0; b < N - 1; b += 2)
    pat += msg_str("%d_", N / 2 - b / 2 - 1);
  if (pat.back() == '_')
    pat.pop_back();
  pat += 'e';

  for (int b = 0; b < N - 2 + N % 2; b += 2)
    pat += msg_str(",%d_%df%d_%df", b / 2, b / 2 + 1, b / 2 + 1, b / 2);
  if (N % 2)
    pat += msg_str(",%dE", N / 2);
  else
    pat += msg_str(",%dv%df", N / 2 - 1, N / 2 - 1);

  return pat;
}

static string g_pattern(int N)
{
  if (N < 1)
    return ""; // number out of range

  string pat = "[V";
  int divs = 2 * N + 1;
  for (int b = 0; b < N; b++) {
    const int e_coord = 2 * (b + 1);
    pat += "," + coord_string(Vec3d(divs - e_coord, e_coord, 0));
  }
  pat += ",F]";

  int F_idx = N + 1;

  if (N == 1)
    pat += msg_str("1_2F1_0V1E");
  else
    pat += msg_str("%d_1_0e1_2e", F_idx);

  for (int b = 0; b < N - 1; b++) {
    pat += msg_str(",%d", F_idx);
    int div_start = 2 * b + 1;
    bool past_center = false;
    for (int i = 0; i < 3; i++) {
      int div = div_start + i;
      char op = '_';
      if (div > N && !past_center) {
        op = 'v';
        past_center = true;
      }
      pat += msg_str("%c%d", op, (div <= N) ? div : 2 * N + 1 - div);
    }
    if (past_center)
      pat += 'v';
  }

  pat += msg_str(",%dE", N);

  return pat;
}

static string s_pattern(int N)
{
  if (N < 1)
    return ""; // number out of range

  string pat = "[";
  int divs = N;
  for (int b = 0; b < N / 2 + 1; b++) {
    const int e_coord = 2 * b;
    pat += coord_string(Vec3d(divs - e_coord, e_coord, 1)) + ',';
  }
  pat.back() = ']';

  auto div2idx = [](int idx, int N) { return (idx <= N / 2) ? idx : N - idx; };

  pat += msg_str("0V,%dE,", N / 2);

  bool past_center = false;
  for (int b = 0; b < N; b++) {
    string op = (b) ? "_" : "";
    if (2 * b > N && !past_center) {
      op = "v";
      past_center = true;
    }
    pat += msg_str("%s%d", op.c_str(), div2idx(b, N));
  }
  if (past_center)
    pat += 'v';

  pat += 'F';

  for (int b = 0; b < N / 2; b++) {
    pat += msg_str(",%d_%df%df", div2idx(b, N), div2idx(b + 1, N),
                   div2idx(N - (b + 1), N));
    pat += msg_str(",%d_f%d_%df", div2idx(b, N), div2idx(N - (b + 1), N),
                   div2idx(N - b, N));
  }
  if (!is_even(N))
    pat += msg_str(",%df%dv%dvf", div2idx(N / 2, N), div2idx(N - N / 2, N),
                   div2idx(N - (N / 2 + 1), N));

  return pat;
}

static string l_pattern(int N)
{
  if (N < 0)
    return ""; // number out of range

  string pat = "[";
  int divs = N + 1;
  for (int i = 0; i < divs; i++) {
    pat += coord_string(normalize_point(Vec3d(divs - i, 0, i))) + ',';
  }
  pat.back() = ']';

  pat += msg_str("%dF,", divs - 1);

  for (int i = 0; i < divs - 1; i++)
    pat += msg_str("%d_%dv%d_%dv,", i, i + 1, i + 1, i);

  pat += msg_str("0E");

  return pat;
}

static string L_pattern(int N)
{
  if (N < 0)
    return ""; // number out of range

  if (N == 0)
    return "[V,EF]1F,1e1_0e,1_0E"; // L0 is its own operator

  string pat = "[";
  int divs = N + 1;
  for (int i = 0; i < divs; i++) {
    if (i == 0)
      pat += coord_string(Vec3d(1, 0, 0)) + ',';
    else if (i == 1)
      pat += coord_string(Vec3d(0, 5, 1)) + ',';
    else if (i % 2 == 0)
      pat += coord_string(Vec3d(1, 0, pow(5, i - 1))) + ',';
    else if (i % 2 == 1)
      pat += coord_string(Vec3d(0, 1, pow(2, i / 2) * pow(5, i - 2))) + ',';
  }
  pat.back() = ']';

  pat += msg_str("%dF,", divs - 1);

  for (int i = 0; i < divs - 1; i++) {
    const char *c = (i % 2) ? "ve" : "ev";
    pat += msg_str("%d%c%d_%d%c,", i + 1, c[0], i + 1, i, c[0]);
    pat += msg_str("%d%c%d_%d%c,", i, c[1], i, i + 1, c[1]);
  }

  pat += msg_str("0E");

  return pat;
}

class Tri {
public:
  vector<Vec3d> verts;
  string op;
  Tri(const vector<Vec3d> &verts, string op) : verts(verts), op(op) {}
  bool cart2bary(const Vec3d &P, Vec3d &Q);
};

bool Tri::cart2bary(const Vec3d &P, Vec3d &Q)
{
  auto &A = verts[0];
  auto &B = verts[1];
  auto &C = verts[2];
  double det = (B[1] - C[1]) * (A[0] - C[0]) + (C[0] - B[0]) * (A[1] - C[1]);
  Q[0] = ((B[1] - C[1]) * (P[0] - C[0]) + (C[0] - B[0]) * (P[1] - C[1])) / det;
  Q[1] = ((C[1] - A[1]) * (P[0] - C[0]) + (A[0] - C[0]) * (P[1] - C[1])) / det;
  Q[2] = 1 - Q[0] - Q[1];

  for (int i = 0; i < 3; i++) {
    if (double_eq(Q[i], 0))
      Q[i] = 0;
    if (double_lt(Q[i], 0))
      return false;
  }

  return true;
}

static string invert_op(char op)
{
  const char *c = strchr("VEF", op);
  if (c == nullptr) // mirror or no_op
    return string(1, op);
  else if (*c == 'V')
    return "fe";
  else if (*c == 'E')
    return "vf";
  else if (*c == 'F')
    return "ev";
  else
    return string(); // error
}

static string invert_op(const string &op)
{
  string inv;
  for (auto it = op.rbegin(); it < op.rend(); it++)
    inv += invert_op(*it);
  return inv;
}

// remove repeated mirrors and trailing _
static string normalize_pattern(const string &pattern)
{
  const char all_ops[] = "_VEFvef";
  string new_pattern;

  const char *p = pattern.c_str();
  bool points_sec = true;
  while (*p) {
    size_t len = 1; // number of characters consumed
    if (points_sec) {
      if (*p == ']')
        points_sec = false;
      new_pattern.push_back(*p);
    }
    else { // on path
      if (strchr(all_ops, *p)) {
        len = strspn(p, all_ops);
        string new_ops; // expand rotations
        for (size_t i = 0; i < len; i++) {
          if (p[i] == 'V')
            new_ops += "ef";
          else if (p[i] == 'E')
            new_ops += "fv";
          else if (p[i] == 'F')
            new_ops += "ve";
          else
            new_ops.push_back(p[i]);
        }

        new_ops = repeatedly_remove_duplicates(new_ops);
        if (new_ops.empty() && *(p + len) != '\0' && *(p + len) != ',')
          new_ops = "_"; // need this to seperate point index numbers
        new_pattern += new_ops;
      }
      else
        new_pattern.push_back(*p);
    }
    p += len;
  }

  return new_pattern;
}

static string u_pattern(int N, int M)
{
  // One of the parameters can be 0, but not both
  if (N < 1 && M < 1)
    return "";

  // Colour indexes for labelling tiles
  enum {
    tile_invalid = -1,   // invalid
    tile_unused = 0,     // outside the meta triangles of interest
    tile_center = 1,     // centred on the face centre
    tile_inner_edge = 2, // on an internal edge between meta triangles
    tile_normal = 3      // any other tile in the meta triangles
  };

  // Make a geodesic division of an equilateral triangle
  // Base model includes a triangle on each edge to capture the points
  // of the geodesic division outside the central triangle
  Geometry base;
  base.add_vert({-2, 0, 0});       // central triangle has vertex to the left
  base.add_vert({1, sqrt(3), 0});  // label V0
  base.add_vert({1, -sqrt(3), 0}); // label V1
  for (int i = 0; i < 3; i++)
    base.add_vert(-2 * base.verts(i));

  base.add_face({0, 1, 2}, Color(tile_unused)); // only central face has index
  base.add_face({3, 2, 1});
  base.add_face({4, 0, 2});
  base.add_face({5, 1, 0});
  Geometry geo;
  make_geodesic_planar(geo, base, M, N);

  // Find tiles that lie on a positive and negative meta triangle
  for (size_t i = 0; i < geo.faces().size(); i++) {
    auto cent = geo.face_cent(i); // tile centre, use for testing position
    bool on_central_face = (geo.colors(FACES).get(i).is_index());
    if (on_central_face && double_ge(cent[0], 0)) { // in sector to the right
      int tile_type = tile_invalid;
      if (double_eq(cent.len(), 0.0))
        tile_type = tile_center;
      else if (double_eq(sqrt(3) * cent[0], cent[1]))
        tile_type = tile_inner_edge; // top edge only (avoid duplication)
      else if (double_gt(sqrt(3) * fabs(cent[0]), fabs(cent[1])))
        tile_type = tile_normal; // inside sector

      // The right edge of tiles has no tiles shared with right original face
      // and so there is no possibility of duplication

      if (tile_type > tile_unused)
        geo.colors(FACES).set(i, tile_type);
    }
  }

  // Label some points of the meta triangle tiling of the base triangles
  const Vec3d V0 = base.verts(1);  // 60
  const Vec3d V1 = base.verts(2);  // -60
  const Vec3d V2 = base.verts(0);  // 180
  const Vec3d E01 = (V0 + V1) / 2; // 0
  const Vec3d E02 = (V0 + V2) / 2; // 120
  const Vec3d E12 = (V1 + V2) / 2; // -120
  const Vec3d F(0, 0, 0);          // centre (origin)
  const Vec3d F01 = 2 * E01;       // 0
  auto F_r = F - V2;               // Face centre to right

  // These are all the meta triangles which might include points
  vector<Tri> meta_triangles({
      Tri({V0, E01, F}, "_"),               // base
      Tri({V1, E01, F}, "v"),               // base negative
      Tri({V0, E02, F}, "e"),               // fan anticlockwise from base
      Tri({V0, E01, F01}, "f"),             // outside base
      Tri({V0, F_r + 0.5 * V0, F_r}, "ef"), // outside outside base
      Tri({V1, E01, F01}, "vf"),            // outside base negative
      Tri({V1, E12, F}, "ve"),              // fan clockwise from base negative
  });

  // base vertex index -> meta triangle index and barycentric coords
  map<int, pair<int, Vec3d>> vert_idx2tri_pt;

  // barycentric coords -> pattern point index
  auto vec_cmp = [](const Vec3d &v1, const Vec3d &v2) {
    return compare(v1, v2) == -1;
  };
  map<Vec3d, size_t, decltype(vec_cmp)> b_coords2pat_pt(vec_cmp);

  const int to_int_factor = M * M + M * N + N * N; // for integer coordinates
  string pattern_paths;
  string pattern_edge_paths;

  for (size_t i = 0; i < geo.faces().size(); i++) {
    auto col = geo.colors(FACES).get(i);
    int tile_type = col.is_index() ? col.get_index() : tile_invalid;
    int edge_start_idx = -1;       // face offset where an edge tile starts
    if (tile_type > tile_unused) { // a valid tile
      vector<int> tile_pts;
      vector<string> tile_ops;
      for (int v_idx = 0; v_idx < 3; v_idx++) {
        const auto P = geo.face_v(i, v_idx);

        Vec3d Q;    // point in barycentric coordinates for containing meta
        size_t tri; // will hold index of meta triangle where point found
        for (tri = 0; tri < meta_triangles.size(); tri++)
          if (meta_triangles[tri].cart2bary(P, Q)) // P is in meta triangle
            break;

        if (tri < meta_triangles.size()) { // point found in a meta triangle
          Q = normalize_point(Q * to_int_factor); // to low integer coordinates
          // If Q has already been seen, get the pattern index number,
          // otherwise set it
          auto it_bool = b_coords2pat_pt.insert({Q, b_coords2pat_pt.size()});
          auto pat_pt_idx = it_bool.first->second; // idx is iterator value

          string op = meta_triangles[tri].op;
          if (tile_type == tile_center) {
            tile_ops.push_back((op == "_" || is_even(op.size())) ? "_" : "v");
            tile_pts.push_back(pat_pt_idx);
            break; // only need one point for the central triangle
          }
          else { // tile_inner_edge or tile_normal
            tile_ops.push_back(op);
            tile_pts.push_back(pat_pt_idx);
          }
        }
        else {                          // point not found in a meta triangle
          if (tile_type != tile_center) // not an error for centre triangle
            return "";                  // INTERNAL ERROR
        }

        // test for edge on E, to convert to tile that winds E
        if (compare(geo.edge_cent(
                        {geo.faces(i, v_idx), geo.faces_mod(i, v_idx + 1)}),
                    E01) == 0)
          edge_start_idx = v_idx;
      }

      // Have found tiles as paths of points on meta triangles, need to
      // convert to paths of points flipping between meta triangles
      string tile_str;
      if (tile_type == tile_center) {
        // No extra handling required (handled when point was found)
        if (tile_ops[0] == "v")
          tile_str += '-';
        tile_str += to_string(tile_pts[0]) + "F";
      }
      else { // tile_inner_edge or tile_normal
        // handling will depend on whether path is positive or negative
        // get the index of the first _, otherwise the first v
        int first = -1; // index of first _ if found, or first v
        for (int i = 0; i < 3; i++) {
          if (tile_ops[i] == "_") { // if found then done
            first = i;
            break;
          }
          else if (tile_ops[i] == "v" && first == -1) // only record first v
            first = i;
        }

        if (first == -1) { // shouldn't happen
          return "";       // INTERNAL ERROR
        }
        else if (tile_ops[first] == "_") {
          for (int i = 0; i < 3; i++) {
            int idx = (i + first) % 3; // index starts at operator '_'
            string start_op = (i) ? tile_ops[idx] : ""; // skip leading op
            tile_str +=
                start_op + to_string(tile_pts[idx]) + invert_op(start_op);
          }
        }
        else if (tile_ops[first] == "v") {
          tile_str = "-";
          for (int i = 0; i < 3; i++) {
            int idx = (i + first) % 3; // index starts at operator 'v'
            // conjugate the operation with 'v' to reuse the base operators
            string start_op = (i) ? "v" + tile_ops[idx] : ""; // skip leading op
            tile_str +=
                start_op + to_string(tile_pts[idx]) + invert_op(start_op);
          }
        }
        else {
          return ""; // INTERNAL ERROR, should not occur
        }

        if (edge_start_idx != -1) { // edge tile found
          int e0 = edge_start_idx;
          int e1 = (e0 + 1) % 3;
          // choose e0 to prefer '_' over 'v'
          if ((tile_ops[e1] == "_") ||
              ((tile_ops[e1] == "v") && (tile_ops[e0] != "_")))
            swap(e0, e1);
          pattern_edge_paths = ((tile_ops[e0] == "v") ? "-" : "") +
                               to_string(tile_pts[e0]) + "E,";
        }
      }
      if (!tile_str.empty())
        pattern_paths += tile_str + ",";
    }
  }

  if (M + N == 1)
    pattern_edge_paths = "0E,"; // isn't set by the algorithm for 1,0

  // make an array of the pattern points in index order
  vector<Vec3d> pat_points(b_coords2pat_pt.size());
  for (const auto &kv : b_coords2pat_pt)
    pat_points[kv.second] = kv.first;

  // write the points in the pattern format
  string pattern = "[";
  for (const auto &coords : pat_points)
    pattern += coord_string(coords) + ",";
  pattern.back() = ']';

  pattern += pattern_paths + pattern_edge_paths;
  pattern.pop_back(); // remove the final ','
  pattern = normalize_pattern(pattern);

  return pattern;
}

static string oe_pattern(char op, int N, int M)
{
  // One of the parameters can be 0, but not both
  if (N < 1 && M < 1)
    return "";

  // Colour indexes for labelling tiles
  enum {
    tile_invalid = -1, // invalid
    tile_face = 0,     // centred on a face centre
    tile_edge = 1,     // centred on an edge centre
    tile_vert = 2,     // centred on a vertex centre
    tile_normal = 3,   // any other tile in a meta triangle
  };

  // Label some points of the meta triangle tiling of the base triangles
  const Vec3d V0(1, 1, 0);         // 45
  const Vec3d V1(1, -1, 0);        // -45
  const Vec3d V2(-1, -1, 0);       // -135
  const Vec3d V3(-1, 1, 0);        // 135
  const Vec3d E01 = (V0 + V1) / 2; // 0
  const Vec3d E30 = (V3 + V0) / 2; // 90
  const Vec3d E12 = (V1 + V2) / 2; // -90
  const Vec3d F(0, 0, 0);          // centre (origin)
  const Vec3d F01 = 2 * E01;       // 0

  // These are all the meta triangles which might include points
  vector<Tri> meta_triangles({
      Tri({V0, E01, F}, "_"),    // base
      Tri({V1, E01, F}, "v"),    // base negative
      Tri({V0, E30, F}, "e"),    // fan anticlockwise from base
      Tri({V0, E01, F01}, "f"),  // outside base
      Tri({V1, E01, F01}, "vf"), // outside base negative
      Tri({V1, E12, F}, "ve"),   // fan clockwise from base negative
  });

  // barycentric coords -> pattern point index
  auto vec_cmp = [](const Vec3d &v1, const Vec3d &v2) {
    return compare(v1, v2) == -1;
  };
  map<Vec3d, size_t, decltype(vec_cmp)> b_coords2pat_pt(vec_cmp);

  const int to_int_factor = 2 * (M * M + N * N); // XXXfor integer coordinates
  string pattern_paths;
  string pattern_edge_paths;

  Geometry geom;
  geom.add_verts({V0, V1, V2, V3, F});
  geom.add_face({0, 1, 4});
  double tile_edge_len = 2 / sqrt(N * N + M * M);
  Vec3d vx = Vec3d(N, M, 0).with_len(tile_edge_len);
  Vec3d vy = Vec3d(-M, N, 0).with_len(tile_edge_len);
  for (int i = -(N + M + 1); i < N + M + 1; i++) {   // enough
    for (int j = -(N + M + 1); j < N + M + 1; j++) { // enough
      Vec3d off = -0.5 * (vx + vy);
      auto C = is_even(M + N + (op == 'e')) * off + i * vx + j * vy;
      const auto x = C[0];
      const auto y = C[1];
      int tile_type = tile_invalid;
      if (double_ge(x, 0) && double_le(x, 1) && // 0 <= x <= 1
          double_le(std::abs(y), 1)) {          // -1 <= y <= 1
        if (compare(C, Vec3d(0, 0, 0)) == 0)
          tile_type = tile_face;
        else if (compare(C, Vec3d(1, 0, 0)) == 0)
          tile_type = tile_edge;
        else if (compare(C, Vec3d(1, 1, 0)) == 0)
          tile_type = tile_vert;
        else if ((double_lt(x, 1) && double_gt(x, std::abs(y))) || // interior
                 double_eq(x, y) ||                                // e mirror
                 (double_eq(x, 1) && double_ge(y, 0)))             // f mirror
          tile_type = tile_normal;
      }
      if (tile_type == tile_invalid)
        continue;

      // The tile is valid
      vector<Vec3d> points = {C + (vx + vy) / 2, C + (vx - vy) / 2,
                              C + (-vx - vy) / 2, C + (-vx + vy) / 2};

      int vsz = geom.verts().size();
      geom.add_verts(points);
      geom.add_face({vsz, vsz + 1, vsz + 2, vsz + 3});

      int edge_start_idx = -1; // face offset where an edge tile starts
      vector<int> tile_pts;
      vector<string> tile_ops;
      for (int v_idx = 0; v_idx < 4; v_idx++) {
        const auto P = points[v_idx];

        Vec3d Q;    // point in barycentric coordinates for containing meta
        size_t tri; // will hold index of meta triangle where point found
        for (tri = 0; tri < meta_triangles.size(); tri++)
          if (meta_triangles[tri].cart2bary(P, Q)) // P is in meta triangle
            break;

        if (tri < meta_triangles.size()) { // point found in a meta triangle
          Q = normalize_point(Q * to_int_factor); // to low integer coords
          // If Q has already been seen, get the pattern index number,
          // otherwise set it
          auto it_bool = b_coords2pat_pt.insert({Q, b_coords2pat_pt.size()});
          auto pat_pt_idx = it_bool.first->second; // idx is iterator value

          string op = meta_triangles[tri].op;
          if (tile_type == tile_normal) {
            tile_ops.push_back(op);
            tile_pts.push_back(pat_pt_idx);
          }
          else { // tile cycles an element
            tile_ops.push_back((op == "_" || is_even(op.size())) ? "_" : "v");
            tile_pts.push_back(pat_pt_idx);
            if (tile_type != tile_edge || tile_ops.size() > 1)
              break; // need one point for F and V, and two for E
          }
        }
        else {                            // point not found in a meta triangle
          if (tile_type == tile_normal) { // not an error for centre triangle
            fprintf(stderr, "INTERNAL ERROR: wythoff o_%d_%d (point: %s)\n", N,
                    M, P.to_str().c_str());
            return ""; // INTERNAL ERROR
          }
        }

        // test for edge on E, to convert to tile that winds E
        if (compare((points[v_idx] + points[(v_idx + 1) % 4]) / 2, E01) == 0)
          edge_start_idx = v_idx;
      }

      // Have found tiles as paths of points on meta triangles, need to
      // convert to paths of points flipping between meta triangles
      string tile_str;
      if (tile_type != tile_normal) {
        // No extra handling required (handled when point was found)
        if (tile_ops[0] == "v")
          tile_str += '-';
        string elems("FEV");
        tile_str += to_string(tile_pts[0]);
        if (tile_type == tile_edge)
          tile_str += "v" + to_string(tile_pts[1]) + "v";
        tile_str += elems[tile_type];
      }
      else { // tile_normal
        // handling will depend on whether path is positive or negative
        // get the index of the first _, otherwise the first v
        int first = -1; // index of first _ if found, or first v
        for (int i = 0; i < 4; i++) {
          if (tile_ops[i] == "_") { // if found then done
            first = i;
            break;
          }
          else if (tile_ops[i] == "v" && first == -1) // only record first v
            first = i;
        }

        if (first == -1) { // shouldn't happen
          string t;
          for (int i = 0; i < 4; i++)
            t += tile_ops[i] + to_string(tile_pts[i]);
          fprintf(stderr, "INTERNAL ERROR: wythoff o_%d_%d (path: %s)\n", N, M,
                  t.c_str());
          return ""; // INTERNAL ERROR
        }
        else if (tile_ops[first] == "_") {
          for (int i = 0; i < 4; i++) {
            int idx = (i + first) % 4; // index starts at operator '_'
            string start_op = (i) ? tile_ops[idx] : ""; // skip leading op
            tile_str +=
                start_op + to_string(tile_pts[idx]) + invert_op(start_op);
          }
        }
        else if (tile_ops[first] == "v") {
          tile_str = "-";
          for (int i = 0; i < 4; i++) {
            int idx = (i + first) % 4; // index starts at operator 'v'
            // conjugate the operation with 'v' to reuse the base
            // operators
            string start_op = (i) ? "v" + tile_ops[idx] : ""; // skip leading op
            tile_str +=
                start_op + to_string(tile_pts[idx]) + invert_op(start_op);
          }
        }
        else {
          string t;
          for (int i = 0; i < 4; i++)
            t += tile_ops[i] + to_string(tile_pts[i]);
          fprintf(stderr, "INTERNAL ERROR: wythoff o_%d_%d (path: %s)\n", N, M,
                  t.c_str());
          return ""; // INTERNAL ERROR, should not occur
        }

        if (edge_start_idx != -1) { // edge tile found
          int e0 = edge_start_idx;
          int e1 = (e0 + 1) % 4;
          // choose e0 to prefer '_' over 'v'
          if ((tile_ops[e1] == "_") ||
              ((tile_ops[e1] == "v") && (tile_ops[e0] != "_")))
            swap(e0, e1);
          pattern_edge_paths = ((tile_ops[e0] == "v") ? "-" : "") +
                               to_string(tile_pts[e0]) + "E,";
        }
      }
      if (!tile_str.empty())
        pattern_paths += tile_str + ",";
    }
  }

  if (M + N == 1)
    pattern_edge_paths = "0E,"; // isn't set by the algorithm for 1,0

  // make an array of the pattern points in index order
  vector<Vec3d> pat_points(b_coords2pat_pt.size());
  for (const auto &kv : b_coords2pat_pt)
    pat_points[kv.second] = kv.first;

  // write the points in the pattern format
  string pattern = "[";
  for (const auto &coords : pat_points)
    pattern += coord_string(coords) + ",";
  pattern.back() = ']';

  pattern += pattern_paths + pattern_edge_paths;
  pattern.pop_back(); // remove the final ','
  pattern = normalize_pattern(pattern);

  return pattern;

  // geom.write("tmp.off");
}

Status Tiling::read_conway(const string &op)
{
  const auto unset = std::numeric_limits<int>::max(); // unset parameter value
  int param1 = unset;
  int param2 = unset;
  char op_char; // Conway operation letter
  char buff;
  bool valid_format =
      (sscanf(op.c_str(), "%c%c", &op_char, &buff) == 1 ||             // X
       sscanf(op.c_str(), "%c%d%c", &op_char, &param1, &buff) == 2 ||  // X1
       sscanf(op.c_str(), "%c_%d%c", &op_char, &param1, &buff) == 2 || // X_1
       sscanf(op.c_str(), "%c_%d_%d%c", &op_char, &param1, &param2,
              &buff) == 3 // X_1_2
      );
  if (!valid_format)
    return Status::error("Conway operator '" + op + "': invalid format");

  // Look up operator
  int op_idx = -1;             // index number of operator in list
  vector<string> param_ops(3); // op letters for each number of paramieters
  int last_op = sizeof(conway_operator_list) / sizeof(conway_operator_list[0]);
  for (int i = 0; i < last_op; i++) {
    const char c = conway_operator_list[i].operator_short[0];
    param_ops[conway_operator_list[i].num_params] += c;
    if (c == op_char)
      op_idx = i;
  }

  auto prefix = msg_str("Conway operator %c: ", op_char);
  if (op_idx < 0) // not found
    return Status::error(prefix + "operator not recognised");

  int number_of_params = (param1 != unset) + (param2 != unset);
  if (number_of_params > conway_operator_list[op_idx].num_params) {
    if (conway_operator_list[op_idx].num_params == 0)
      return Status::error(prefix + "does not take a parameter");
    else
      return Status::error(prefix + "too many parameters for operator");
  }
  if (param1 < 0 || param2 < 0)
    return Status::error(prefix + "parameter cannot be negative");

  string pat;            // final wythoff pattern
  if (param1 == unset) { // no parameter given, use base pattern
    pat = conway_operator_list[op_idx].pattern;
  }
  else { // parameter given, use sequence processing
    if (strchr(param_ops[1].c_str(), op_char)) {
      if (op_char == 'M')
        pat = M_pattern(param1);
      else if (op_char == 'm')
        pat = m_pattern(param1);
      else if (op_char == 'b')
        pat = b_pattern(param1);
      else if (op_char == 'g')
        pat = g_pattern(param1);
      else if (op_char == 's')
        pat = s_pattern(param1);
      else if (op_char == 'l')
        pat = l_pattern(param1);
      else if (op_char == 'L')
        pat = L_pattern(param1);
      else
        return Status::error(
            prefix + "one-parameter operator not handled (internal error)");

      if (pat == "")
        return Status::error(prefix + msg_str("invalid parameter %d", param1));
    }
    else if (strchr(param_ops[2].c_str(), op_char)) {
      if (op_char == 'u')
        pat = u_pattern((param1 == unset) ? 3 : param1,
                        (param2 == unset) ? 0 : param2);
      if (op_char == 'o')
        pat = oe_pattern('o', (param1 == unset) ? 2 : param1,
                         (param2 == unset) ? 0 : param2);
      if (op_char == 'e')
        pat = oe_pattern('e', (param1 == unset) ? 2 : param1,
                         (param2 == unset) ? 0 : param2);
      if (pat == "") {
        if (param2 == unset)
          return Status::error(prefix +
                               msg_str("invalid parameter %d", param1));
        else
          return Status::error(
              prefix + msg_str("invalid parameters %d, %d", param1, param2));
      }
    }
  }

  return read_pattern(pat);
}

string Tiling::pattern_string()
{
  string pat = "[";
  for (auto v : points)
    pat += coord_string(v.first) + ",";
  pat.back() = ']';
  for (auto path : pat_paths)
    pat += path.tile_string() + ",";

  pat.pop_back();
  return pat;
}

void Tiling::print_conway_list(FILE *ofile)
{
  fprintf(ofile, "%-4s%-5s%-15s%s\n", "Op", "Args", "Description",
          "Tiling Pattern");
  fprintf(ofile, "%-4s%-5s%-15s%s\n", "--", "----", "-----------",
          "--------------");
  int last_op = sizeof(conway_operator_list) / sizeof(conway_operator_list[0]);
  for (int i = 0; i < last_op; i++) {
    ConwayOperator op = conway_operator_list[i];
    string params = op.num_params ? msg_str("[%d]", op.num_params) : "   ";
    fprintf(ofile, "%-4s%-5s%-15s%s\n", op.operator_short.c_str(),
            params.c_str(), op.operator_name.c_str(), op.pattern.c_str());
  }
  fprintf(ofile, R"(
Some operators are part of a sequence having an integer parameter.
  m, g, b, M, s, l, L:  1 is the base, 0 is a lower level operator
                        (except L_0, which is a standalone operator)
  u, o, e:              2 is the base, 1 is a lower level operator

Operators g, o, e take two integer parameters.

Examples e, M, m_3, g_3, s_1, u_3, u_2_3, o_3_3, e_0_4

Some Conway operators, like t and k, take a number to filter the elements
that the pattern will be applied to, but this is not supported.
)");
}

namespace anti {
// export these functions

Status wythoff_make_tiling(Geometry &tiled_geom, const Geometry &base_geom,
                           const std::string &pat, bool oriented, bool reverse,
                           TilingColoring col_type)
{
  Tiling tiling;
  Status stat =
      (pat[0] == '[') ? tiling.read_pattern(pat) : tiling.read_conway(pat);
  if (!stat.is_error()) {
    tiling.set_geom(base_geom); // not meta, so will not fail
    if (!oriented)
      tiling.start_everywhere();
    if (reverse)
      tiling.reverse_pattern();
    tiling.make_tiling(tiled_geom, col_type);
    if (!oriented) // some tiles may be doubled
      merge_coincident_elements(tiled_geom, "f");
  }
  return stat;
}

/// Get vertex points of a Schwarz triangle, and its symmetry group
bool get_schwarz_tri_verts(const vector<int> &fracs, vector<Vec3d> &verts,
                           Symmetry *sym)
{
  int ret = assign_vertices(fracs, verts);
  if (ret && sym) {
    *sym = Symmetry(get_tri_symmetry(fracs));
  }
  return ret;
};

bool get_schwarz_tri_fracs(int tri_idx, vector<int> &fracs)
{
  if (tri_idx < 0 || tri_idx >= num_schwarz_tris)
    return false;
  fracs.resize(6);
  for (int i = 0; i < 6; i++)
    fracs[i] = schwarz_triangles[tri_idx][i];
  return true;
};

} // namespace anti
