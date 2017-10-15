/*
   Copyright (c) 2012-2016, Adrian Rossiter

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

#include <regex>
#include <string.h>
#include <string>
#include <vector>

#include "geometry.h"
#include "geometryinfo.h"
#include "symmetry.h"
#include "tiling.h"
#include "utils.h"

#include "private_std_polys.h"

using std::string;
using std::vector;
using std::map;
using std::pair;
using std::swap;

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
      sym += itostr(fracs[2 * i]);
      if (fracs[2 * i + 1] > 1)
        sym += "/" + itostr(fracs[2 * i + 1]);
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
    sym = "D" + itostr(fs[4]);
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
  char fracs_str[MSG_SZ];
  strncpy(fracs_str, sym_norm2.c_str(), MSG_SZ);
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
                              char *msg = nullptr)
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

  if (msg) {
    *msg = '\0';
    double max_ang = 0.0;
    for (int i = 0; i < 3; i++) {
      double ang = angle_around_axis(v[i], v[(i + 1) % 3], pt);
      if (ang > M_PI)
        ang = 2 * M_PI - ang;
      double ang_diff = fabs(2 * M_PI / 3 - ang);
      max_ang = ang > max_ang ? ang_diff : max_ang;
    }
    if (max_ang > epsilon)
      sprintf(msg, "innacurate calculation of fermat point "
                   "(angle difference %g)",
              max_ang);
  }

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

bool Wythoff::make_poly(Geometry &geom, char *errmsg)
{
  if (errmsg)
    *errmsg = '\0';
  geom.clear_all();
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
      if (errmsg)
        sprintf(errmsg, "symbol leads to nonconstructible antiprism");
      return false;
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
      f_pt = get_fermat_point(verts[0], verts[1], verts[2], degenerate, errmsg);

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
      geom.add_face(0, 1, -1); // add as face, for sizing by edge length
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
    return false;

  merge_coincident_elements(geom, "v", epsilon);
  return true;
}

bool Wythoff::make_tri(Geometry &geom)
{
  geom.clear_all();
  if (is_set()) {
    geom.add_verts(verts);
    geom.add_face(0, 1, 2, -1);
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

  auto ef_pairs = geom.get_edge_face_pairs();
  for (auto &ef : ef_pairs) {
    // The edge and face pair make a quadrilateral
    // Add the centre to this quadrilateral
    int e_idx = meta.add_vert(geom.edge_cent(ef.first), Color(1));
    // Add four triangles
    if (ef.second[0] >= 0) {
      meta.add_face(ef.first[0], e_idx, ef.second[0] + f_start, -1);
      meta.add_face(ef.first[1], e_idx, ef.second[0] + f_start, -1);
    }
    if (ef.second[1] >= 0) {
      meta.add_face(ef.first[1], e_idx, ef.second[1] + f_start, -1);
      meta.add_face(ef.first[0], e_idx, ef.second[1] + f_start, -1);
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

Status Tile::read(const string &pat)
{
  // initialise
  ops.clear();
  idxs.clear();
  bool has_tris_spec = strchr("+-*", pat[0]);
  start_faces = (has_tris_spec) ? pat[0] : '+';
  unsigned int pos = has_tris_spec;
  if (!std::isdigit(pat[pos]))
    return Status::error("tile specifier: first character (or first "
                         "character after +-*)  must be a digit");

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

  return Status::ok();
}

void Tile::relabel(vector<int> relab)
{
  for (auto &op : ops)
    if (0 <= op && op < 3)
      op = relab[op];
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
        tile += itostr(idxs[p_idx]);
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
};

ConwayOperator conway_operator_list[]{
    // Equivalent: d, a, s
    {"d", "dual", "[F]0V,0E"},
    {"a", "ambo", "[E]0F,0V"},
    {"S", "seed", "[V]0E,0F"},

    {"j", "join", "[F,V]0_1E"},

    // Equivalent: k, n, u
    {"k", "kis", "[F,V]0_1v1v,1E"},
    {"n", "needle", "[V,F]0_1f1f,1E"},
    {"u", "subdivide", "[V,E]0_1e1e,1F"},

    // Equivalent: t, z, e
    {"t", "truncate", "[VE]0v0e,0V,0E"},
    {"z", "zip", "[EF]0e0f,0E,0F"},
    {"e", "expand", "[FV]0f0v,0F,0V"},

    // Symmetric: s, m, b
    {"s", "snub", "[VEF]0V,0E,0F,0V0E0F"},
    {"m", "meta", "[V,E,F]*0_1_2"},
    {"b", "bevel", "[VEF]0v0e,0e0f,0f0v"},

    {"o", "ortho", "[V,E,F]1_0e1_2e"},
    {"g", "gyro", "[F,VE,V]1_0F1_2V1E,1E"},
    {"c", "chamfer", "[V,VF]1F,0_1v1f"},
    {"l", "loft", "[V,VF]1F,0_1v1_0v,0E"},
    {"p", "propeller", "[V,VEF]1F,1_0V1E1F,1E"},
    {"q", "quinto", "[V,E,EF]2F,0_1_2e2_1e"},
    {"L0", "joined-lace", "[V,EF2]1F,1e1_0e,1_0E"},
    {"L", "lace", "[V,EF2]1F,1e1_0e,1_0v0v,0E"},
    {"K", "stake", "[V,EF2,F]0_1_2e1e,1_0v0v,0E"},
    {"M0", "joined-medial", "[F,V,EF]*0_1_2,1_2E"},
    {"M3", "edge-medial-3", "[F,V,VE2]0_2_1e2e,2_0v2v"},
    {"m3", "medial-3", "[F,V,VE]*0_1_2,2_0v2v"},
    {"b3", "bevel3", "[VEF,E2F]1_0e0v,0e0f,*1_0f0_1f,1E"},
    {"o3", "ortho3", "[VF2,V,VE2]0_2e1_2e,0_2v2_0v,0F,2E"},
    {"X", "cross", "[V,E,F,VF]3_1v3_2v,*0_1_3"},
};

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

static bool valid_start_face(int f, int start_faces)
{
  int pos_tri = f % 2;
  return !((start_faces == '-' && pos_tri) || (start_faces == '+' && !pos_tri));
}

Status Tiling::make_tiling(Geometry &geom, vector<int> *tile_counts) const
{
  geom.clear_all();
  if (tile_counts)
    tile_counts->resize(pat_paths.size(), 0);

  // All the possible element inclusion postions V, E, F, VE, EF, FV, VEF.
  // Each entry maps to order (to find index of corresponding point)
  // and example triangle (to generate coordinates of corresponding point)
  auto ef_pairs = meta.get_edge_face_pairs();
  vector<map<vector<int>, pair<int, int>>> index_order(7);
  for (int i = 0; i < (int)meta.faces().size(); i++) {
    const auto &face = meta.faces(i);
    index_order[Tile::V][{face[Tile::V]}] = {-1, i};
    index_order[Tile::E][{face[Tile::E]}] = {-1, i};
    index_order[Tile::F][{face[Tile::F]}] = {-1, i};
    index_order[Tile::VE][make_edge(face[Tile::V], face[Tile::E])] = {-1, i};
    index_order[Tile::EF][make_edge(face[Tile::E], face[Tile::F])] = {-1, i};
    index_order[Tile::FV][make_edge(face[Tile::F], face[Tile::V])] = {-1, i};
    index_order[Tile::VEF][{i}] = {-1, i};
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
    for (auto &m : index_order[incl])
      geom.add_vert(point_on_face(meta, m.second.second, crds), pt.second);
  }

  int faces_sz = meta.faces().size();
  for (unsigned int p_idx = 0; p_idx < pat_paths.size(); p_idx++) {
    const auto &pat = pat_paths[p_idx];
    // Check index range
    auto out_of_range = pat.check_index_range(points.size());
    if (out_of_range.size()) {
      string msg = "Path" + itostr(p_idx) + ": index numbers out of range:";
      for (auto idx : out_of_range)
        msg += " " + itostr(idx) + ",";
      msg.pop_back();
      return Status::error(msg.c_str());
    }

    Color col(p_idx);
    vector<bool> seen(faces_sz, false);
    int start_faces_sz = geom.faces().size();
    unsigned char start_faces = pat.get_start_faces();
    for (int i = 0; i < faces_sz; i++) {
      if (!seen[i] && valid_start_face(i, start_faces)) {
        add_circuit(geom, i, pat, seen, col, index_order, point_vertex_offsets);
        if (one_of_each_tile)
          break;
      }
    }
    if (tile_counts)
      tile_counts->at(p_idx) = geom.faces().size() - start_faces_sz;
  }

  geom.del(VERTS, geom.get_info().get_free_verts());
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
  std::regex re_coord("[VEF]([-+]?([0-9]*\\.[0-9]+|[0-9]+))?");
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
    int idx = elem_idx[next->str()[0]];
    if (seen[idx])
      return Status::error(
          msg_str("coordinates %c given more than once", next->str()[0]));
    else
      seen[idx] = true;
    if (next->str().size() < 2)
      coords[idx] = 1;
    else {
      Status stat = read_double(next->str().c_str() + 1, &coords[idx]);
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
  std::regex r_all("^\\[(.*)](.*)");
  std::regex_match(pat2, m_all, r_all);
  if (m_all.size() < 3)
    return Status::error(
        msg_str("pattern '%s': not in form [Point0,Point1,...]Path0,Path1...",
                pat.c_str()));

  std::unique_ptr<char> pat_str(copy_str(m_all[1].str().c_str()));
  vector<char *> parts;
  int num_parts = split_line(pat_str.get(), parts, ",");
  points.resize(num_parts);
  for (int i = 0; i < num_parts; i++) {
    auto stat = read_point(parts[i], points[i]);
    if (stat.is_error())
      return Status::error(msg_str("Point%d: ", i) + stat.msg());
  }

  std::unique_ptr<char> paths(copy_str(m_all[2].str().c_str()));
  num_parts = split_line(paths.get(), parts, ",");
  pat_paths.resize(num_parts);
  for (int i = 0; i < num_parts; i++) {
    auto stat = pat_paths[i].read(parts[i]);
    if (stat.is_error())
      return Status::error(msg_str("Path%d: ", i) + stat.msg());
  }
  return Status::ok();
}

Status Tiling::relabel_pattern(const string &relabel)
{
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

static string coord_string(Vec3d v)
{
  string VEF = "VEF";
  string coords;
  for (int i = 0; i < 3; i++)
    if (v[i]) {
      coords += VEF[i];
      if (v[i] != 1.0)
        coords += msg_str("%g", v[i]);
    }

  return coords;
}

static Status medial_pattern(char op, int N, string &pat)
{
  if (N < 2)
    return Status::error(msg_str("invalid number %d", N));

  pat = "[F";
  for (int i = 0; i < N + 1; i += 2) {
    double e_num = i;
    double v_num = N - e_num;
    pat += "," + coord_string(Vec3d(v_num, e_num, 0));
  }
  int last_idx = N / 2 + 1;

  if (op == 'M')
    pat += "]0_2_1e2e";
  else if (op == 'm')
    pat += "]*0_1_2";

  for (int i = 2; i < last_idx; i++)
    pat += msg_str(",*0_%d_%d", i, i + 1);

  if (!is_even(N)) {
    pat += msg_str(",%d_0v%dv", last_idx, last_idx);
    pat += msg_str(",%dE", last_idx);
  }

  return Status::ok();
}

static Status ortho_pattern(char /*op*/, int N, string &pat)
{
  if (N < 2)
    return Status::error(msg_str("invalid number %d", N));

  pat = "[";
  for (int a = 0; a <= N; a += 2)
    for (int b = 0; b <= a; b += 2)
      pat += coord_string(Vec3d(((a + N % 2) - b), b, (N - (a + N % 2)))) + ",";
  pat.back() = ']';

  auto crds2idx = [](int a, int b) { return (a / 2 + 1) * a / 4 + b / 2; };
  for (int a = 0; a < N - N % 2; a += 2)
    for (int b = 0; b < a; b += 2)
      pat += msg_str("*%d_%d_%d_%d,", crds2idx(a, b), crds2idx(a, b + 2),
                     crds2idx(a + 2, b + 4), crds2idx(a + 2, b + 2));

  for (int a = 0; a < N - N % 2; a += 2)
    pat += msg_str("%d_%de%d_%de,", crds2idx(a, 0), crds2idx(a + 2, 2),
                   crds2idx(a + 2, 0), crds2idx(a + 2, 2));

  if (N % 2) {
    for (int a = 0; a < N - N % 2; a += 2)
      pat += msg_str("%d_%dv%d_%dv,", crds2idx(a, a), crds2idx(a + 2, a + 2),
                     crds2idx(a + 2, a + 2), crds2idx(a, a));
    pat += msg_str("0F,%dE,", crds2idx(N - 1, N - 1));
  }
  pat.pop_back();

  return Status::ok();
}

Status Tiling::read_conway(const string &op)
{
  int last_op = sizeof(conway_operator_list) / sizeof(conway_operator_list[0]);
  for (int i = 0; i < last_op; i++)
    if (conway_operator_list[i].operator_short == op)
      return read_pattern(conway_operator_list[i].pattern);

  Status stat;
  string pat;
  char op_char;
  int op_int;
  char buff;
  if (sscanf(op.c_str(), "%c%d%c", &op_char, &op_int, &buff) == 2) {
    if (op_int < 1)
      return Status::error("Conway operater number must be a positive integer");
    if (op_char == 'M' || op_char == 'm')
      stat = medial_pattern(op_char, op_int, pat);
    else if (op_char == 'o')
      stat = ortho_pattern(op_char, op_int, pat);
    else
      stat.set_error("Conway operator '" + op + "' not known");
  }
  else
    stat.set_error("Conway operator '" + op + "' not known");

  if (stat.is_error())
    return stat;

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
  fprintf(ofile, "%-5s%-15s%s\n", "Op", "Description", "Tiling Pattern");
  fprintf(ofile, "%-5s%-15s%s\n", "--", "-----------", "--------------");
  int last_op = sizeof(conway_operator_list) / sizeof(conway_operator_list[0]);
  for (int i = 0; i < last_op; i++) {
    ConwayOperator op = conway_operator_list[i];
    fprintf(ofile, "%-5s%-15s%s\n", op.operator_short.c_str(),
            op.operator_name.c_str(), op.pattern.c_str());
  }
  fprintf(
      ofile,
      "\nOperators M, m and o accept a general positive integer, e.g. M5\n");
}

namespace anti {
// export these functions

Status wythoff_make_tiling(Geometry &tiled_geom, const Geometry &base_geom,
                           const std::string &pat, bool pat_is_conway_op)
{
  Tiling tiling;
  Status stat =
      (pat_is_conway_op) ? tiling.read_conway(pat) : tiling.read_pattern(pat);
  if (!stat.is_error()) {
    tiling.set_geom(base_geom); // not meta, so will not fail
    tiling.make_tiling(tiled_geom);
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
