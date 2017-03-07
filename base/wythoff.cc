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

#include <string.h>
#include <string>
#include <vector>

#include "geometry.h"
#include "symmetry.h"
#include "utils.h"

#include "private_std_polys.h"

using std::string;
using std::vector;
using std::swap;

using namespace anti;

// Schwarz triangles
const int num_schwarz_tris = 44;

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
    {
        //  0
        0, 1, 0, 1 / sqrt_3, 1 / sqrt_3, -1 / sqrt_3, 1 / sqrt_3, 1 / sqrt_3,
        1 / sqrt_3,
    },
    {
        //  1
        0, 1, 0, 1 / sqrt_3, -1 / sqrt_3, 1 / sqrt_3, 1 / sqrt_3, 1 / sqrt_3,
        -1 / sqrt_3,
    },
    {
        //  2
        1 / sqrt_2, 1 / sqrt_2, 0, 1 / sqrt_3, 1 / sqrt_3, 1 / sqrt_3, 1, 0, 0,
    },
    {
        //  3
        1 / sqrt_2, 1 / sqrt_2, 0, -1 / sqrt_3, -1 / sqrt_3, -1 / sqrt_3, 1, 0,
        0,
    },
    {
        //  4
        0.5 / phi, phi / 2, 0.5, 1 / sqrt_3, 1 / sqrt_3, 1 / sqrt_3, 0,
        1 / sqrt_phi_plus_2, phi / sqrt_phi_plus_2,
    },
    {
        //  5
        0.5, -0.5 / phi, phi / 2, phi / sqrt_3, (phi - 1) / sqrt_3, 0, 0,
        1 / sqrt_phi_plus_2, phi / sqrt_phi_plus_2,
    },
    {
        //  6
        0, 1, 0, phi / sqrt_3, -(phi - 1) / sqrt_3, 0, 0, 1 / sqrt_phi_plus_2,
        phi / sqrt_phi_plus_2,
    },
    {
        //  7
        0.5 / phi, phi / 2, 0.5, 0, -phi / sqrt_3, -(phi - 1) / sqrt_3, 0,
        1 / sqrt_phi_plus_2, phi / sqrt_phi_plus_2,
    },
    {
        //  8
        0, -1, 0, 1 / sqrt_3, 1 / sqrt_3, -1 / sqrt_3, 1 / sqrt_3, 1 / sqrt_3,
        1 / sqrt_3,
    },
    {
        //  9
        -1 / sqrt_2, 1 / sqrt_2, 0, -1 / sqrt_3, 1 / sqrt_3, -1 / sqrt_3, 1, 0,
        0,
    },
    {
        // 10
        -1 / sqrt_2, 0, -1 / sqrt_2, 1 / sqrt_3, 1 / sqrt_3, 1 / sqrt_3, 1, 0,
        0,
    },
    {
        // 11
        0.5 / phi, -phi / 2, -0.5, 1 / sqrt_3, -1 / sqrt_3, -1 / sqrt_3, 0,
        1 / sqrt_phi_plus_2, phi / sqrt_phi_plus_2,
    },
    {
        // 12
        0.5, 0.5 / phi, -phi / 2, phi / sqrt_3, -(phi - 1) / sqrt_3, 0, 0,
        1 / sqrt_phi_plus_2, phi / sqrt_phi_plus_2,
    },
    {
        // 13
        0, -1, 0, phi / sqrt_3, (phi - 1) / sqrt_3, 0, 0, 1 / sqrt_phi_plus_2,
        phi / sqrt_phi_plus_2,
    },
    {
        // 14
        -0.5, -0.5 / phi, -phi / 2, 1 / sqrt_3, 1 / sqrt_3, 1 / sqrt_3, 0,
        1 / sqrt_phi_plus_2, phi / sqrt_phi_plus_2,
    },
    {
        // 15
        0.5 / phi, phi / 2, 0.5, -1 / sqrt_phi_plus_2, phi / sqrt_phi_plus_2, 0,
        0, 1 / sqrt_phi_plus_2, phi / sqrt_phi_plus_2,
    },
    {
        // 16
        0.5 / phi, phi / 2, 0.5, 1 / sqrt_phi_plus_2, -phi / sqrt_phi_plus_2, 0,
        0, 1 / sqrt_phi_plus_2, phi / sqrt_phi_plus_2,
    },
    {
        // 17
        0.5 / phi, -phi / 2, -0.5, 0, 1 / sqrt_phi_plus_2,
        phi / sqrt_phi_plus_2, -1 / sqrt_phi_plus_2, -phi / sqrt_phi_plus_2, 0,
    },
    {
        // 18
        0.5 / phi, -phi / 2, -0.5, 0, 1 / sqrt_phi_plus_2,
        phi / sqrt_phi_plus_2, 1 / sqrt_phi_plus_2, phi / sqrt_phi_plus_2, 0,
    },
    {
        // 19
        1 / sqrt_3, 1 / sqrt_3, -1 / sqrt_3, 1 / sqrt_3, -1 / sqrt_3,
        1 / sqrt_3, 1 / sqrt_3, 1 / sqrt_3, 1 / sqrt_3,
    },
    {
        // 20
        1 / sqrt_3, 1 / sqrt_3, 1 / sqrt_3, (phi - 1) / sqrt_3, 0, phi / sqrt_3,
        0, 1 / sqrt_phi_plus_2, phi / sqrt_phi_plus_2,
    },
    {
        // 21
        phi / sqrt_3, (phi - 1) / sqrt_3, 0, -1 / sqrt_3, -1 / sqrt_3,
        1 / sqrt_3, 0, 1 / sqrt_phi_plus_2, phi / sqrt_phi_plus_2,
    },
    {
        // 22
        phi / sqrt_3, -(phi - 1) / sqrt_3, 0, phi / sqrt_3, (phi - 1) / sqrt_3,
        0, 0, 1 / sqrt_phi_plus_2, phi / sqrt_phi_plus_2,
    },
    {
        // 23
        0, -phi / sqrt_3, -(phi - 1) / sqrt_3, 1 / sqrt_3, 1 / sqrt_3,
        1 / sqrt_3, 0, 1 / sqrt_phi_plus_2, phi / sqrt_phi_plus_2,
    },
    {
        // 24
        1 / sqrt_3, 1 / sqrt_3, -1 / sqrt_3, 0, 0, 1, 0, 1, 0,
    },
    {
        // 25
        1 / sqrt_3, 1 / sqrt_3, 1 / sqrt_3, 0, -1 / sqrt_phi_plus_2,
        phi / sqrt_phi_plus_2, 0, 1 / sqrt_phi_plus_2, phi / sqrt_phi_plus_2,
    },
    {
        // 26
        1 / sqrt_3, 1 / sqrt_3, 1 / sqrt_3, -1 / sqrt_phi_plus_2,
        -phi / sqrt_phi_plus_2, 0, 0, 1 / sqrt_phi_plus_2,
        phi / sqrt_phi_plus_2,
    },
    {
        // 27
        phi / sqrt_3, -(phi - 1) / sqrt_3, 0, 0, 1 / sqrt_phi_plus_2,
        phi / sqrt_phi_plus_2, 1 / sqrt_phi_plus_2, phi / sqrt_phi_plus_2, 0,
    },
    {
        // 28
        1 / sqrt_3, -1 / sqrt_3, -1 / sqrt_3, 0, 1 / sqrt_phi_plus_2,
        phi / sqrt_phi_plus_2, 1 / sqrt_phi_plus_2, phi / sqrt_phi_plus_2, 0,
    },
    {
        // 29
        1 / sqrt_3, 1 / sqrt_3, -1 / sqrt_3, 1 / sqrt_3, -1 / sqrt_3,
        1 / sqrt_3, -1 / sqrt_3, 1 / sqrt_3, 1 / sqrt_3,
    },
    {
        // 30
        1 / sqrt_3, -1 / sqrt_3, -1 / sqrt_3, (phi - 1) / sqrt_3, 0,
        -phi / sqrt_3, 0, 1 / sqrt_phi_plus_2, phi / sqrt_phi_plus_2,
    },
    {
        // 31
        phi / sqrt_3, -(phi - 1) / sqrt_3, 0, -1 / sqrt_3, 1 / sqrt_3,
        -1 / sqrt_3, 0, 1 / sqrt_phi_plus_2, phi / sqrt_phi_plus_2,
    },
    {
        // 32
        1 / sqrt_3, 1 / sqrt_3, -1 / sqrt_3, 0, 1, 0, 1, 0, 0,
    },
    {
        // 33
        1 / sqrt_3, 1 / sqrt_3, -1 / sqrt_3, 0, -1, 0, 0, 0, 1,
    },
    {
        // 34
        1 / sqrt_3, 1 / sqrt_3, 1 / sqrt_3, 0, 1 / sqrt_phi_plus_2,
        phi / sqrt_phi_plus_2, 1 / sqrt_phi_plus_2, phi / sqrt_phi_plus_2, 0,
    },
    {
        // 35
        1 / sqrt_3, 1 / sqrt_3, 1 / sqrt_3, 1 / sqrt_phi_plus_2,
        -phi / sqrt_phi_plus_2, 0, 0, 1 / sqrt_phi_plus_2,
        phi / sqrt_phi_plus_2,
    },
    {
        // 36
        phi / sqrt_3, (phi - 1) / sqrt_3, 0, 0, 1 / sqrt_phi_plus_2,
        phi / sqrt_phi_plus_2, 1 / sqrt_phi_plus_2, -phi / sqrt_phi_plus_2, 0,
    },
    {
        // 37
        1 / sqrt_3, 1 / sqrt_3, -1 / sqrt_3, 0, 1 / sqrt_phi_plus_2,
        phi / sqrt_phi_plus_2, 1 / sqrt_phi_plus_2, -phi / sqrt_phi_plus_2, 0,
    },
    {
        // 38
        1 / sqrt_3, -1 / sqrt_3, -1 / sqrt_3, 0, 1 / sqrt_phi_plus_2,
        phi / sqrt_phi_plus_2, 0, 1 / sqrt_phi_plus_2, -phi / sqrt_phi_plus_2,
    },
    {
        // 39
        1 / sqrt_3, -1 / sqrt_3, -1 / sqrt_3, 0, 1 / sqrt_phi_plus_2,
        phi / sqrt_phi_plus_2, -1 / sqrt_phi_plus_2, phi / sqrt_phi_plus_2, 0,
    },
    {
        // 40
        0, 1 / sqrt_phi_plus_2, phi / sqrt_phi_plus_2, 1 / sqrt_phi_plus_2,
        -phi / sqrt_phi_plus_2, 0, 0, -1 / sqrt_phi_plus_2,
        phi / sqrt_phi_plus_2,
    },
    {
        // 41
        0, 1 / sqrt_phi_plus_2, phi / sqrt_phi_plus_2, 1 / sqrt_phi_plus_2,
        phi / sqrt_phi_plus_2, 0, -1 / sqrt_phi_plus_2, phi / sqrt_phi_plus_2,
        0,
    },
    {
        // 42
        1 / sqrt_phi_plus_2, -phi / sqrt_phi_plus_2, 0, 0, 1 / sqrt_phi_plus_2,
        phi / sqrt_phi_plus_2, 1 / sqrt_phi_plus_2, phi / sqrt_phi_plus_2, 0,
    },
    {
        // 43
        0, 1 / sqrt_phi_plus_2, phi / sqrt_phi_plus_2, 1 / sqrt_phi_plus_2,
        -phi / sqrt_phi_plus_2, 0, 0, 1 / sqrt_phi_plus_2,
        -phi / sqrt_phi_plus_2,
    },
};

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
/*
   verts.resize(3);
   vector<Vec3d> v_map(3);
   for(int i=0; i<3; i++)
      v_map[i] = Vec3d((double)i, 0.0, 0.0);

   vector<int> norm_fracs = fracs;
   tri_normalise(norm_fracs, v_map);
   vector<Vec3d> norm_verts;
   bool ret = get_tri_verts(norm_fracs, norm_verts);
   if(ret) {
      for(int i=0; i<3; i++)
         verts[(int)floor(v_map[i][0]+0.5)] = norm_verts[i];
   }

   return ret;
}

*/

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
  if (gr_fact % 2 == 0) // even number of windings cancel out
    return;
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

  if(sides>2)
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
  return is_set();
}

namespace anti {

// export this function
bool get_schwarz_tri_verts(const vector<int> &fracs, vector<Vec3d> &verts,
                           Symmetry *sym)
{
  int ret = assign_vertices(fracs, verts);
  if (ret && sym) {
    *sym = Symmetry(get_tri_symmetry(fracs));
  }
  return ret;
};

// export this function
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
