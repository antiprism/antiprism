/*
   Copyright (c) 2008-2016, Roger Kaufman

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
   Name: iso_delta.cc
   Description: Isohedral Deltahedra
                Based on a paper by G. C. Shephard
                Periodica Mathematica Hungarica, Volume 39, Numbers 1-3, 2000 ,
   pp. 83-106(24).
                with enhancements by Jim McNeill
   (http://www.orchidpalms.com/polyhedra)
                and Adrian Rossiter (http://www.antiprism.com)
   Project: Antiprism - http://www.antiprism.com
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <ctype.h>

#include <set>
#include <string>
#include <vector>

#include "../base/antiprism.h"

using std::string;
using std::vector;
using std::set;
using std::swap;

using namespace anti;

struct IsoDeltaItem {
  const char *name;
  const char *sym_type;
  const char *symbol;
  const char *comment;
};

struct IsoDeltaVector {
  const char *name;
  Vec3d A;
  Vec3d B;
  Vec3d C;
  double alpha;
  double beta;
  double gamma;
};

// clang-format off
IsoDeltaItem iso_delta_item_list[] = {
   {"T1(1)",   "Oh", "[1/3,2/3,1/3]", "Augmented Cube (=O5(1))"},
   {"T1(2)",   "Oh", "[2/3,2/3,2/3]", "Excavated Cube (=O5(2))"},
   {"T1(3)",   "Td", "[1/3,2/3,1/3]", "Additonal Isomer of T1(1) (Mobius)"},
   {"T2(1)",   "Td", "[1/3,1/3,2/3]", "Augmented Tetrahedron"},
   {"T2(2)",   "Td", "[2/3,1/3,1/3]", "Tetrahedron"},
   {"O1",      "Oh", "[1/2,1/2,2/3]", "COMPOUND of 3 Octahedra (UD08)"},
   {"O2(1)",   "Oh", "[1/2,1/3,1/4]", "Relaxed Excavated Rhombic Dodecahedron (Mobius)"},
   {"O2(2)",   "Oh", "[1/2,1/3,1/4]", "Relaxed Augmented Rhombic Dodecahedron (Mobius)"},
   {"O3",      "Oh", "[1/2,1/4,2/4]", "COMPOUND of 3 8/3 Star Bipyramids"},
   {"O4(1)",   "Oh", "[1/3,1/3,2/3]", "COMPOUND of 2 Augmented Tetrahedra T2(1) (UC54_d)"},
   {"O4(2)",   "Oh", "[2/3,1/3,1/3]", "COMPOUND of 2 Tetrahedra (Stella Octangula UC04)"},
   // next 2 redundent
   {"O5(1)",   "Oh", "[1/3,2/3,1/3]", "Augmented Cube (=T1(1)) (Repeat)"},
   {"O5(2)",   "Oh", "[2/3,2/3,2/3]", "Excavated Cube (=T1(2)) (Repeat)"},
   {"O5(3)",   "Oh", "[1/3,2/3,1/3]", "COMPOUND of 2 T1(3)"},
   {"O6(1)",   "Oh", "[2/3,1/4,1/4]", "Augmented Octahedron"},
   {"O6(2)",   "Oh", "[2/3,3/4,3/4]", "Excavated Octahedron"},
   {"O7",      "Oh", "[2/4,2/4,2/4]", "Octahedron"},
   {"I1(1)",   "Ih", "[1/2,1/3,1/5]", "Relaxed Excavated Rhombic Triacontahedron (Mobius)"},
   {"I1(2)",   "Ih", "[1/2,1/3,1/5]", "Relaxed Augmented Rhombic Triacontahedron (Mobius)"},
   {"I2(1)",   "Ih", "[1/2,2/5,4/5]", "Relaxed Augmented Did (U36) *"},
   {"I2(2)",   "Ih", "[1/2,3/5,4/5]", "Relaxed Excavated Did (U36) *"},
   {"I3(1)",   "Ih", "[1/2,1/3,2/5]", "Relaxed Augmented Gid (U54) *"},
   {"I3(2)",   "Ih", "[1/2,1/3,2/5]", "Relaxed Excavated Gid (U54) *"},
   {"I3(3)",   "Ih", "[1/2,2/3,3/5]", "Relaxed Dual of Gaquatid (UD73) *"},
   {"I3(4)",   "Ih", "[1/2,2/3,2/5]", "Additional Isomer of Augmented RT I1(2) **"},
   {"I4",      "Ih", "[1/2,1/2,1/2]", "COMPOUND of 5 Octahedra (UC17)"},
   {"I5(1)",   "Ih", "[2/3,3/5,4/5]", "Relaxed Augmented Dual of Ditdid (UD41)"},
   {"I5(2)",   "Ih", "[1/3,2/5,4/5]", "Relaxed Excavated Dual of Ditdid (UD41)"},
   {"I6(1)",   "Ih", "[2/3,1/5,1/5]", "Augmented Icosahedron"},
   {"I6(2)",   "Ih", "[2/3,4/5,4/5]", "Excavated Icosahedron"},
   {"I7(1)",   "Ih", "[1/3,2/5,3/5]", "Excavated Gike (U53)"},
   {"I7(2)",   "Ih", "[2/3,2/5,3/5]", "Augmented Gike (U53)"},
   {"I7(3)",   "Ih", "[1/3,2/5,3/5]", "Additional Isomer of Augmented Icosahedron I6(1)"},
   {"I8(1)",   "Ih", "[1/3,1/3,4/5]", "Augmented Gissid (U52)"},
   {"I8(2)",   "Ih", "[1/3,1/3,4/5]", "Excavated Gissid (U52)"},
   {"I9(1)",   "Ih", "[1/3,1/3,2/5]", "Augmented Dodecahedron"},
   {"I9(2)",   "Ih", "[1/3,1/3,2/5]", "Excavated Dodecahedron"},
   {"I9(3)",   "Ih", "[1/3,1/3,2/5]", "Relaxed Excavated Dual of Sidtid (UD30)"},
   // redundant 
   {"I9(4)",   "Ih", "[1/3,1/3,2/5]", "Relaxed Excavated Dual of Sidtid (UD30) (Repeat)"},
   {"I10(1)",  "Ih", "[4/5,4/5,4/5]", "Great Icosahedron or Gike (U53)"},
   {"I10(2)",  "Ih", "[1/5,1/5,4/5]", "Augmented Sissid (U34)"},
   // redundant 
   {"I10(3)",  "Ih", "[1/5,1/5,4/5]", "Augmented Sissid (U34) (Repeat)"},
   {"I11(1)",  "Ih", "[2/5,2/5,2/5]", "Icosahedron"},
   {"I11(2)",  "Ih", "[2/5,3/5,3/5]", "Excavated Gad (U35)"},
};


// notes added as to what changed from what was listed in the paper

// values found to be formulas
// 0.20710678118654752 = -(1-sqrt(2))/2
// 0.237248933904040   = (1.0/6)*sqrt(21-6*sqrt(10))
// 0.33079226912480375 = sqrt(3)/(sqrt(5)+3)
// 0.5352331346596349  = sqrt(3)/(sqrt(5)+1)
// 0.541196100146197   = sqrt((2-sqrt(2))/2)
// 0.58778525229247312 = 0.5*sqrt(0.5*(5-sqrt(5)))
// 0.84089641525371454 = sqrt(sqrt(2))
// 0.95105651629515357 = 0.5*sqrt(0.5*(5+sqrt(5)))
// 1.02062072615965754 = sqrt(25/24)
// 1.05374551483176583 = (1.0/6)*sqrt(21+6*sqrt(10))
// 1.11803398874989485 = sqrt(5)/2
// 1.20710678118654752 = (1+sqrt(2))/2
// 1.40125853844407354 = -sqrt(3)/(-sqrt(5)+1)
// 1.90211303259030714 = sqrt(phi)*5^(1/4) = sqrt(0.5*(5+sqrt(5))) = sqrt(phi+2) = sqrt_phi_plus_2 (built in constant)

IsoDeltaVector iso_delta_vector_list[] = {
   {"T1(1)",   Vec3d(1,-1,1)/sqrt_3, Vec3d(1,0,0), Vec3d(1,-1,-1)/sqrt_3,
               sqrt_3/2, (1+sqrt_2)/2, sqrt_3/2},
   {"T1(2)",   Vec3d(1,-1,1)/sqrt_3, Vec3d(-1,0,0), Vec3d(1,-1,-1)/sqrt_3,
               sqrt_3/2, -(1-sqrt_2)/2, sqrt_3/2},
   {"T1(3)",   Vec3d(1,-1,1)/sqrt_3, Vec3d(1,0,0), Vec3d(1,-1,-1)/sqrt_3,
               (1.0/6)*sqrt(21-6*sqrt(10)), sqrt(5)/2, (1.0/6)*sqrt(21+6*sqrt(10))},
   // B changed
   {"T2(1)",   Vec3d(1,-1,1)/sqrt_3, Vec3d(1,-1,-1)/sqrt_3, Vec3d(-1,-1,-1)/sqrt_3,
               sqrt(6)/4, sqrt(25.0/24), sqrt(6)/4},
   {"T2(2)",   Vec3d(1,-1,1)/sqrt_3, Vec3d(-1,1,1)/sqrt_3, Vec3d(-1,-1,-1)/sqrt_3,
               sqrt(6)/4, sqrt(6)/4, sqrt(6)/4},
   {"O1",      Vec3d(1,0,1)/sqrt_2, Vec3d(1,0,-1)/sqrt_2, Vec3d(0,1,0),
               sqrt_2/2, sqrt_2/2, sqrt_2/2},
   // alpha changed
   {"O2(1)",   Vec3d(1,0,0), Vec3d(1,0,1)/sqrt_2, Vec3d(1,1,1)/sqrt_3,
               1.080488239084420, 0.118828664898468, 1.094667047615130},
   {"O2(2)",   Vec3d(1,0,0), Vec3d(1,0,1)/sqrt_2, Vec3d(1,1,1)/sqrt_3,
               1.204738778767500, 1.375617671088820, 0.515547971500924},
   {"O3",      Vec3d(1,0,1)/sqrt_2, Vec3d(-1,0,0), Vec3d(0,1,0),
               sqrt((2-sqrt_2)/2), sqrt((2-sqrt_2)/2), sqrt(sqrt_2)},
   {"O4(1)",   Vec3d(1,-1,1)/sqrt_3, Vec3d(1,-1,-1)/sqrt_3, Vec3d(-1,-1,-1)/sqrt_3,
               sqrt(6)/4, sqrt(25.0/24), sqrt(6)/4},
   // B changed
   {"O4(2)",   Vec3d(-1,1,-1)/sqrt_3, Vec3d(1,-1,-1)/sqrt_3, Vec3d(1,1,1)/sqrt_3,
               sqrt(6)/4, sqrt(6)/4, sqrt(6)/4},
   {"O5(1)",   Vec3d(1,-1,1)/sqrt_3, Vec3d(1,0,0), Vec3d(1,-1,-1)/sqrt_3,
               sqrt_3/2, (1+sqrt_2)/2, sqrt_3/2},
   {"O5(2)",   Vec3d(1,-1,1)/sqrt_3, Vec3d(-1,0,0), Vec3d(1,-1,-1)/sqrt_3,
               sqrt_3/2, -(1-sqrt_2)/2, sqrt_3/2},
   {"O5(3)",   Vec3d(1,-1,1)/sqrt_3, Vec3d(1,0,0), Vec3d(1,-1,-1)/sqrt_3,
               (1.0/6)*sqrt(21-6*sqrt(10)), sqrt(5)/2, (1.0/6)*sqrt(21+6*sqrt(10))},
   // triangle flipped for algorithm: A<-->C, C<-->A, alpha<-->gamma, gamma<-->alpha
   {"O6(1)",   Vec3d(1,0,0), Vec3d(0,-1,0), Vec3d(1,-1,1)/sqrt(3),
               sqrt(2)/2, sqrt(2)/2, sqrt(6)/2},
   {"O6(2)",   Vec3d(-1,1,-1)/sqrt_3, Vec3d(1,0,0), Vec3d(0,-1,0),
               sqrt(6)/6, sqrt_2/2, sqrt_2/2},
   {"O7",      Vec3d(1,0,0), Vec3d(0,1,0), Vec3d(0,0,1),
               sqrt_2/2, sqrt_2/2, sqrt_2/2},
   // alpha changed
   {"I1(1)",   Vec3d(phi,phi*phi,1)/(2*phi), Vec3d(1,1,1)/sqrt_3, Vec3d(0,phi,1)/sqrt_phi_plus_2,
               0.674181480291942, 1.600435268943570, 1.508572470741690},
   {"I1(2)",   Vec3d(phi,phi*phi,1)/(2*phi), Vec3d(1,1,1)/sqrt_3, Vec3d(0,phi,1)/sqrt_phi_plus_2,
               1.901892201462340, 1.042221422390510, 1.602608615469960},
   // the following six added to the table (not in the paper)
   // I2(1) thru I3(3) added by Jim McNeill
   {"I2(1)",   Vec3d(-phi,phi*phi,-1)/(2*phi), Vec3d(-1,0,-phi)/sqrt_phi_plus_2, Vec3d(0,phi,-1)/sqrt_phi_plus_2,
               1.035310785747180, 1.017991958789870, 0.041794129091506},
   {"I2(2)",   Vec3d(phi,phi*phi,1)/(2*phi), Vec3d(-1,0,-phi)/sqrt_phi_plus_2, Vec3d(0,-phi,-1)/sqrt_phi_plus_2,
               0.069590240312534, 0.961660565811817, 0.940133523007424},
   {"I3(1)",   Vec3d(-phi,phi*phi,1)/(2*phi), Vec3d(1,1,1)/sqrt_3, Vec3d(0,phi,-1)/sqrt_phi_plus_2,
               1.067540004877190, 0.454812301059250, 0.979984198587005},
   {"I3(2)",   Vec3d(phi,-phi*phi,-1)/(2*phi), Vec3d(1,1,1)/sqrt_3, Vec3d(0,phi,-1)/sqrt_phi_plus_2,
               0.365481361740863, 0.809498075533857, 0.758298681854245},
   {"I3(3)",   Vec3d(-phi,phi*phi,1)/(2*phi), Vec3d(1,1,1)/sqrt_3, Vec3d(0,phi,-1)/sqrt_phi_plus_2,
               1.037389492284540, 0.123497623343534, 1.015782486189710},
   // I3(4) result added from programmatic search by Adrian Rossiter
   {"I3(4)",   Vec3d(-phi,phi*phi,1)/(2*phi), Vec3d(1,1,1)/sqrt_3, Vec3d(0,-phi,1)/sqrt_phi_plus_2,
               0.840436103166278, 0.919239756651211, 0.257365254781315},
   // A, B changed
   {"I4",      Vec3d(phi,-phi*phi,-1)/(2*phi), Vec3d(-1,-phi,phi*phi)/(2*phi), Vec3d(phi*phi,1,phi)/(2*phi),
               1/sqrt_2, 1/sqrt_2, 1/sqrt_2},
   // C changed
   {"I5(1)",   Vec3d(1,1,1)/sqrt_3, Vec3d(-1,0,-phi)/sqrt_phi_plus_2, Vec3d(phi,-1,0)/sqrt_phi_plus_2,
               0.740118744863774, 0.305243296610152, 0.825499999909858},
   // alpha changed
   {"I5(2)",   Vec3d(1,1,1)/sqrt_3, Vec3d(-1,0,-phi)/sqrt_phi_plus_2, Vec3d(-phi,1,0)/sqrt_phi_plus_2,
               0.095586833624572, 0.922356501413017, 0.977651218274709},
   {"I6(1)",   Vec3d(1,1,1)/sqrt_3, Vec3d(1,0,phi)/sqrt_phi_plus_2, Vec3d(phi,1,0)/sqrt_phi_plus_2,
               1.572257895003900, 0.5*sqrt_phi_plus_2, 0.5*sqrt_phi_plus_2},
   {"I6(2)",   Vec3d(1,1,1)/sqrt_3, Vec3d(-1,0,-phi)/sqrt_phi_plus_2, Vec3d(-phi,-1,0)/sqrt_phi_plus_2,
               0.060735266851555, 0.5*sqrt_phi_plus_2, 0.5*sqrt_phi_plus_2},
   // C changed
   {"I7(1)",   Vec3d(1,1,1)/sqrt_3, Vec3d(0,-phi,1)/sqrt_phi_plus_2, Vec3d(-phi,1,0)/sqrt_phi_plus_2,
               0.706232491219458, 0.5*sqrt(0.5*(5-sqrt(5))), 0.5*sqrt(0.5*(5-sqrt(5)))},
   {"I7(2)",   Vec3d(1,1,1)/sqrt_3, Vec3d(0,phi,-1)/sqrt_phi_plus_2, Vec3d(phi,-1,0)/sqrt_phi_plus_2,
               0.926760670635994, 0.5*sqrt(0.5*(5-sqrt(5))), 0.5*sqrt(0.5*(5-sqrt(5)))},
   {"I7(3)",   Vec3d(1,1,1)/sqrt_3, Vec3d(0,-phi,1)/sqrt_phi_plus_2, Vec3d(phi,-1,0)/sqrt_phi_plus_2,
               sqrt_3/(sqrt(5)+3), 0.883687468829836, 1.007795749176520},
   // C, beta changed
   {"I8(1)",   Vec3d(1,1,1)/sqrt_3, Vec3d(0,-1/phi,-phi)/sqrt_3, Vec3d(phi,-1,0)/sqrt_phi_plus_2,
               sqrt_3/(sqrt(5)+1), sqrt_3/(sqrt(5)+1), 0.5*sqrt_phi_plus_2},
   // B changed
   {"I8(2)",   Vec3d(1,1,1)/sqrt_3, Vec3d(0,-1/phi,-phi)/sqrt_3, Vec3d(-phi,1,0)/sqrt_phi_plus_2,
               sqrt_3/(sqrt(5)+1), sqrt_3/(sqrt(5)+1), 0.750245100408926},
   {"I9(1)",   Vec3d(1,1,1)/sqrt_3, Vec3d(1,0,phi)/sqrt_phi_plus_2, Vec3d(phi,0,1/phi)/sqrt_3,
               -sqrt_3/(-sqrt(5)+1), 1.639247476530740, -sqrt_3/(-sqrt(5)+1)},
   {"I9(2)",   Vec3d(1,1,1)/sqrt_3, Vec3d(1,0,phi)/sqrt_phi_plus_2, Vec3d(phi,0,1/phi)/sqrt_3,
               -sqrt_3/(-sqrt(5)+1), 0.5*sqrt(0.5*(5-sqrt(5))), -sqrt_3/(-sqrt(5)+1)},
   {"I9(3)",   Vec3d(1,1,1)/sqrt_3, Vec3d(1,0,phi)/sqrt_phi_plus_2, Vec3d(phi,0,1/phi)/sqrt_3,
               1.056261606816530, 1.606723212103540, 1.497317965649610},
   // Redundant
   {"I9(4)",   Vec3d(1,1,1)/sqrt_3, Vec3d(1,0,phi)/sqrt_phi_plus_2, Vec3d(phi,0,1/phi)/sqrt_3,
               1.497317965649610, 1.606723212103540, 1.056261606816530},
   {"I10(1)",  Vec3d(1,0,phi)/sqrt_phi_plus_2, Vec3d(0,-phi,-1)/sqrt_phi_plus_2, Vec3d(-phi,1,0)/sqrt_phi_plus_2,
               0.5*sqrt(0.5*(5-sqrt(5))), 0.5*sqrt(0.5*(5-sqrt(5))), 0.5*sqrt(0.5*(5-sqrt(5)))},
   // C changed
   {"I10(2)",  Vec3d(1,0,phi)/sqrt_phi_plus_2, Vec3d(0,-phi,-1)/sqrt_phi_plus_2, Vec3d(phi,-1,0)/sqrt_phi_plus_2,
               0.5*sqrt(0.5*(5-sqrt(5))), 0.5*sqrt(0.5*(5-sqrt(5))), 1.113516364411610},
   // Redundant
   {"I10(3)",  Vec3d(1,0,phi)/sqrt_phi_plus_2, Vec3d(0,phi,1)/sqrt_phi_plus_2, Vec3d(phi,-1,0)/sqrt_phi_plus_2,
               1.113516364411610, 0.5*sqrt(0.5*(5-sqrt(5))), 0.5*sqrt(0.5*(5-sqrt(5)))},
   {"I11(1)",  Vec3d(1,0,phi)/sqrt_phi_plus_2, Vec3d(0,phi,1)/sqrt_phi_plus_2, Vec3d(phi,1,0)/sqrt_phi_plus_2,
               0.5*sqrt_phi_plus_2, 0.5*sqrt_phi_plus_2, 0.5*sqrt_phi_plus_2},
   {"I11(2)",  Vec3d(1,0,phi)/sqrt_phi_plus_2, Vec3d(0,-phi,-1)/sqrt_phi_plus_2, Vec3d(phi,1,0)/sqrt_phi_plus_2,
               0.5*sqrt_phi_plus_2, 0.100405707943114, 0.5*sqrt_phi_plus_2},
};
// clang-format on

class id_poly {
private:
  int last_iso_delta;
  IsoDeltaItem *iso_delta_items;
  IsoDeltaVector *iso_delta_vectors;

public:
  id_poly();
  void list_poly(int idx, FILE *fp = stderr);
  void list_polys(FILE *fp = stderr);
  int lookup_sym_no(string sym);
  int get_last_iso_delta() { return last_iso_delta; }

  const char *get_sym_type(int i);
  Vec3d A(int i);
  Vec3d B(int i);
  Vec3d C(int i);
  double alpha(int i);
  double beta(int i);
  double gamma(int i);
};

id_poly::id_poly()
{
  iso_delta_items = iso_delta_item_list;
  iso_delta_vectors = iso_delta_vector_list;
  last_iso_delta = sizeof(iso_delta_item_list) / sizeof(iso_delta_item_list[0]);
}

void id_poly::list_poly(int idx, FILE *fp)
{
  fprintf(fp, "%2d) %-7s %-2s %13s %-45s\n", idx + 1, iso_delta_items[idx].name,
          iso_delta_items[idx].sym_type, iso_delta_items[idx].symbol,
          iso_delta_items[idx].comment);
}

void id_poly::list_polys(FILE *fp)
{
  for (int i = 0; i < last_iso_delta; i++)
    list_poly(i, fp);
  fprintf(fp, "\n");
  fprintf(fp,
          "34 Isohedral Deltahedra + 4 Repeats + 6 Compounds = 44 Total\n\n");
  fprintf(fp, "*  added to the table by Jim McNeill "
              "(http://www.orchidpalms.com/polyhedra)\n");
  fprintf(
      fp,
      "** added to the table by Adrian Rossiter (http://www.antiprism.com)\n");
}

int id_poly::lookup_sym_no(string sym)
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

  // is it the list order number
  char *endptr;
  int idx = strtol(sym_norm2.c_str(), &endptr, 10);
  if (!*endptr) // all of string is an integer
    return idx - 1;

  idx = -1;

  // is it a poly name
  for (char &i : sym_norm2)
    if (isalpha(i))
      i = toupper(i);
  for (int i = 0; i < last_iso_delta; i++) {
    if (sym_norm2 == iso_delta_item_list[i].name)
      return i;

    if (idx < 0 &&
        strncmp(sym_norm2.c_str(), iso_delta_item_list[i].name,
                sym_norm2.size()) == 0)
      idx = i;
  }

  return idx;
}

const char *id_poly::get_sym_type(int i)
{
  return iso_delta_item_list[i].sym_type;
}

Vec3d id_poly::A(int i) { return iso_delta_vectors[i].A; }

Vec3d id_poly::B(int i) { return iso_delta_vectors[i].B; }

Vec3d id_poly::C(int i) { return iso_delta_vectors[i].C; }

double id_poly::alpha(int i) { return iso_delta_vectors[i].alpha; }

double id_poly::beta(int i) { return iso_delta_vectors[i].beta; }

double id_poly::gamma(int i) { return iso_delta_vectors[i].gamma; }

class id_opts : public ProgramOpts {
public:
  string ifile;
  string ofile;

  bool list_polys;
  string poly;

  bool make_dipyramid;
  bool triangle_only;
  bool verbose;
  bool allow_angles;
  double angle;
  int n;
  int d;
  int k;
  int s;
  char Coloring_method;
  int face_opacity;

  ColorMapMulti map;
  string case_type;

  id_opts()
      : ProgramOpts("iso_delta"), list_polys(false), make_dipyramid(false),
        triangle_only(false), verbose(false), allow_angles(false),
        angle(INFINITY), n(0), d(1), k(0), s(1), Coloring_method('c'),
        face_opacity(255)
  {
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

// clang-format off
void id_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] polyhedron\n"
"\n"
"Make Isohedral Deltahedra in OFF format. The polyhedron may be specified\n"
"by its list number, or the (start of the) name of the polyhedron.\n"
"Additional examples can be generated using -d or -c paramters.\n"
"Based on a paper by G. C. Shephard\n"
"Periodica Mathematica Hungarica, Volume 39, Numbers 1-3, 2000 , pp. 83-106(24).\n"
"with enhancements by Jim McNeill (http://www.orchidpalms.com/polyhedra)\n"
"and Adrian Rossiter (http://www.antiprism.com)\n"
"\n"
"Options\n"
"%s"
"  -l        display the list of Isohedral Deltahedra 1 thru 42\n"
"  -t        generate triangle only (Isohedral Deltahedra 1 to 42 and option -d)\n"
"  -v        verbose output (Isohedral Deltahedra 1 thru 42 and option -d)\n"
"  -o <file> write output to file (default: write to standard output)\n"
// undocumented switch
//"  -w        allow angle on b,e,f,m,n,o,p,q creates bi-hedral forms\n"
"\nIsohedral Deltahedra Options\n"
"  -a <ang>  angle\n"
"  -n <n/d>  n/d (d is optional)\n"
"              note: for option -d and compound cases c, g, h:\n"
"              n and d must be such that (2 < n/d < 6)\n"
"  -k <int>  in special cases, for specifying number of constituents\n"
"  -s <int>  in special cases, for subtypes (default: 1)\n"
"\nIsohedral Deltahedra Special Cases\n"
"  -d        dipyramid of n/d using -n n/d (infinite set)\n"
"  -c <type> compound cases (a thru f from Shephard's paper)\n"
"              a - tetrahedron repeated k times, evenly spaced\n"
"                     when k=1 tetrahedron, when k=2 Stella Octangula\n"
"                     Uniform Compound Set UC23 when n/d is 2/1\n"
"              b - 5 or 10 tetrahedra\n"
"                     s=1 icosahedral, s=2 with horizontal reflection\n"
"                     Uniform Compounds UC05 and UC06\n"
"              c - 2 dipyramids of n/d using -a angle (default: calculated)\n"
"                     relaxed dual of Uniform Compound Set UC20 and k=1\n"
"              d - 6 octahedra using -a angle (default: 22.5)\n"
"                     At 45.0 degrees is 3 Octahedra\n"
"                     dual of Uniform Compound Set UC07\n"
"              e - 4 or 8 triangular dipyramids\n"
"                     s=1 octahedral, s=2 with horizontal reflection\n"
"                     relaxed duals of Uniform Compounds UC30 & UC31\n"
"              f - 6 or 12 5/1 pentagonal or 5/2 star dipyramids\n"
"                     5/1: s=1 icosahedral, s=2 with horizontal reflection\n"
"                          relaxed duals of Uniform Compounds UC34 & UC35\n"
"                     5/2: s=3 icosahedral, s=4 with horizontal reflection\n"
"                          relaxed duals of Uniform Compounds UC36 & UC37\n"
"\n         additional cases:\n"
"              g - 2 tetrahedra using -a angle (default: 45.0)\n"
"                     At 45.0 degrees is Uniform Compound UC04\n"
"                     Uniform Compound Set UC23 when n/d is 2/1 and k=1\n"
"              h - 2 tetrahedra repeated k times, evenly spaced\n"
"                     using -a angle (default: 1.0)\n"
"                     Uniform Compound Set UC23 when n/d is 2/1 for any k\n"
"              i - 6 tetrahedra using -a angle (default: 45.0)\n"
"                     Uniform Compound UC01. At 45.0 degrees is UC03\n"
"              j - 12 tetrahedra using -a angle (default: 30.0)\n"
"                     Uniform Compound UC02. At 45.0 degrees is UC03\n"
"              k - 2 dipyramids of n/d repeated k times, evenly spaced\n"
"                     using -a angle (default: 1.0)\n"
"                     relaxed dual of Uniform Compound Set UC20\n"
"              l - k dipyramids of n/d using -n n/d, evenly spaced\n"
"                     relaxed dual of Uniform Compound Set UC21\n"
"              m - 10 or 20 triangular dipyramids\n"
"                     s=1 icosahedral, s=2 with horizontal reflection\n"
"                     relaxed duals of Uniform Compounds UC32 & UC33\n"
"              n - 6 10/3 star dipyramids\n"
"                     relaxed dual of Uniform Compounds UC41\n"
"              o - 5 or 10 Augmented Tetrahedra T2(1)\n"
"                     s=1 icosahedral, s=2 with horizontal reflection\n"
"                     relaxed duals of Uniform Compounds UC55 & UC56\n"
"              p - 5 Augmented Octahedra O6(1)\n"
"                     relaxed dual of Uniform Compounds UC57\n"
"              q - 5 Excavated Octahedra O6(2)\n"
"                     relaxed dual of Uniform Compounds UC58\n"
"\nColoring Options (run 'off_util -H color' for help on color formats)\n"
"  -f <opt> compound coloring\n"
"              key word: none - sets no color (default: c)\n"
"              c - unique coloring for each constituent\n"
"              s - symmetric coloring (should always be one color)\n"
"  -T <tran> face transparency. valid range from 0 (invisible) to 255 (opaque)\n"
"  -m <maps> color maps for all elements to be tried in turn (default: compound)\n"
"\n"
"\n", prog_name(), help_ver_text);
}
// clang-format on

void id_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  bool s_is_set = false;
  string map_file;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hldtvwc:n:a:k:s:f:T:m:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'l':
      list_polys = true;
      break;

    case 'd':
      make_dipyramid = true;
      break;

    case 't':
      triangle_only = true;
      break;

    case 'v':
      verbose = true;
      break;

    // undocumented switch
    case 'w':
      allow_angles = true;
      warning("using -w to allow angles creates bihedral forms");
      break;

    case 'c':
      if (strspn(optarg, "abcdefghijklmnopq") != strlen(optarg) ||
          strlen(optarg) > 1)
        error(msg_str("case type is '%s' must be only one of a thru q", optarg),
              c);
      case_type = optarg;
      break;

    case 'a':
      print_status_or_exit(read_double(optarg, &angle), c);
      angle = deg2rad(angle);
      break;

    case 'n': {
      char *p;
      p = strchr(optarg, '/');
      if (p != nullptr) {
        *p++ = '\0';
        print_status_or_exit(read_int(p, &d), "n/d (d part)");
      }

      print_status_or_exit(read_int(optarg, &n), "n/d (n part)");
      if (n < 1)
        error("must be an integer 1 or greater", "n/d (n part)");
      if (d < 1)
        error("d must be 1 or greater", "n/d (d part)");
      break;
    }

    case 'k':
      print_status_or_exit(read_int(optarg, &k), c);
      if (k < 1)
        error("must be an integer 1 or greater", c);
      break;

    case 's':
      print_status_or_exit(read_int(optarg, &s), c);
      if (s < 1)
        error("must be an integer 1 or greater", c);
      s_is_set = true;
      break;

    case 'f':
      if (!strcasecmp(optarg, "none"))
        Coloring_method = '\0';
      else if (strspn(optarg, "cs") != strlen(optarg) || strlen(optarg) > 1)
        error(msg_str("invalid Coloring method '%s'", optarg), c);
      else
        Coloring_method = *optarg;
      break;

    case 'T':
      print_status_or_exit(read_int(optarg, &face_opacity), c);
      if (face_opacity < 0 || face_opacity > 255) {
        error("face transparency must be between 0 and 255", c);
      }
      break;

    case 'm':
      map_file = optarg;
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (case_type == "c" || case_type == "k" || case_type == "l" ||
      make_dipyramid) {
    if (n <= 0)
      error("n must be greater than 0");

    if (d >= n)
      error("d must be less than number of sides", "n/d (d part)");
    if (((double)n / d <= 2.0) || ((double)n / d >= 6.0))
      error("2 < n/d < 6 is enforced");
  }
  else if (n != 0)
    warning("n/d is ignored");

  if ((case_type == "a" || case_type == "l" || make_dipyramid) ||
      (!allow_angles &&
       (case_type == "b" || case_type == "e" || case_type == "f" ||
        case_type == "m" || case_type == "n" || case_type == "o" ||
        case_type == "p" || case_type == "q"))) {
    if (angle != INFINITY) {
      warning("-a angle is ignored");
      angle = INFINITY;
    }
  }

  // check k
  if (case_type == "a" || case_type == "h" || case_type == "k" ||
      case_type == "l") {
    if (k <= 0)
      error("k must be greater then 0");
  }
  else if (k != 0)
    warning("k is ignored");

  // check s
  if (case_type == "b" || case_type == "e" || case_type == "m" ||
      case_type == "o") {
    if (s != 1 && s != 2)
      error("s must 1 or 2");
  }
  else if (case_type == "f") {
    if (s < 1 || s > 4)
      error("s must 1, 2, 3 or 4");
  }
  else if (s_is_set)
    warning("s is ignored");

  // if one from the table or -d then force color by constituent
  if (!case_type.length() &&
      (Coloring_method == 'T' || Coloring_method == 't')) {
    if (Coloring_method == 'T')
      Coloring_method = 'C';
    else if (Coloring_method == 't')
      Coloring_method = 'c';
  }

  if (argc == optind && !list_polys && !case_type.length() && !make_dipyramid)
    error("no polyhedron specified", "polyhedron");

  while (optind < argc) {
    poly += argv[optind];
    poly += " ";
    optind++;
  }

  if (poly.length() != 0 && (case_type.length() || make_dipyramid))
    warning("polyhedron specifier ignored if using -c or -d");

  if (!map_file.size())
    map_file = "compound";
  print_status_or_exit(map.init(map_file.c_str()), 'm');
}

void verbose_output(const Vec3d &A, const Vec3d &B, const Vec3d &C,
                    const double alpha, const double beta, const double gamma)
{
  fprintf(stderr, "Schwarz triangle:\n");
  fprintf(stderr, "A = (% .15lf,% .15lf,% .15lf )\n", A[0], A[1], A[2]);
  fprintf(stderr, "B = (% .15lf,% .15lf,% .15lf )\n", B[0], B[1], B[2]);
  fprintf(stderr, "C = (% .15lf,% .15lf,% .15lf )\n\n", C[0], C[1], C[2]);

  fprintf(stderr, "alpha = % .15lf\n", alpha);
  fprintf(stderr, "beta  = % .15lf\n", beta);
  fprintf(stderr, "gamma = % .15lf\n\n", gamma);

  fprintf(stderr, "Deltahedra triangle:\n");
  fprintf(stderr, "A = (% .15lf,% .15lf,% .15lf )\n", A[0] * alpha,
          A[1] * alpha, A[2] * alpha);
  fprintf(stderr, "B = (% .15lf,% .15lf,% .15lf )\n", B[0] * beta, B[1] * beta,
          B[2] * beta);
  fprintf(stderr, "C = (% .15lf,% .15lf,% .15lf )\n\n", C[0] * gamma,
          C[1] * gamma, C[2] * gamma);
}

void make_triangle(Geometry &geom, const Vec3d &A, const Vec3d &B,
                   const Vec3d &C, Color c)
{
  int sz = geom.verts().size();

  geom.add_vert(A);
  geom.add_vert(B);
  geom.add_vert(C);

  vector<int> face(3);
  for (unsigned int i = 0; i < 3; i++)
    face[i] = sz++;

  geom.add_face(face, c);
}

// quartic formula by Adrian Rossiter
int find_equ_tris(Vec3d A, Vec3d B, Vec3d C, vector<double> &factors,
                  const bool verbose)
{
  /*
     A.dump("\nA");
     B.dump("B");
     C.dump("C");
     // check for unit vectors
     fprintf(stderr, "lens = %g, %g, %g\n\n", A.len(), B.len(), C.len());
  */

  double AB = vdot(A, B);
  double BC = vdot(B, C);
  double CA = vdot(C, A);
  double A2 = vdot(A, A);
  double B2 = vdot(B, B);
  double C2 = vdot(C, C);

  /*
     fprintf(stderr,"AB = %.15lf\n",AB);
     fprintf(stderr,"BC = %.15lf\n",BC);
     fprintf(stderr,"CA = %.15lf\n",CA);
     fprintf(stderr,"A2 = %.15lf\n",A2);
     fprintf(stderr,"B2 = %.15lf\n",B2);
     fprintf(stderr,"C2 = %.15lf\n\n",C2);
  */

  // C.C*(4*B.C*B.C - B.B*C.C) c^4
  // -4*B.C*(2*B.C*C.A + A.B*C.C) c^3
  // +2* (8*A.B*B.C*C.A + A.A*B.B*C.C) c^2
  // -4*A.B*(2*A.B*C.A + B.C*A.A) c
  // + A.A*(4*A.B*A.B - A.A*B.B) = 0

  double coeffs[5];
  coeffs[0] = A2 * (4 * AB * AB - A2 * B2);
  coeffs[1] = -4 * AB * (2 * AB * CA + BC * A2);
  coeffs[2] = 2 * (8 * AB * BC * CA + A2 * B2 * C2);
  coeffs[3] = -4 * BC * (2 * BC * CA + AB * C2);
  coeffs[4] = C2 * (4 * BC * BC - B2 * C2);

  double sol[4] = {0, 0, 0, 0};
  int num_roots = quartic(coeffs, sol);
  if (verbose)
    fprintf(stderr, "\nquartic formula found %d roots\n\n", num_roots);

  factors.clear();

  double c = 1e100;
  for (int i = 0; i < num_roots; i++) {
    /* RK: since matching is used, don't reject any quartic solution
          if(fabs(c-sol[i])<1e-10) {
             if (verbose)
                fprintf(stderr,"rejecting quartic solution: %.15lf\n\n",sol[i]);
             continue;
          }
    */

    c = sol[i];
    if (verbose)
      fprintf(stderr, "c = %.15lf\n", c);

    // b = (A.B +/- sqrt(A.B*A.B + B.B*(C.C*c^2 - 2*C.A*c)) / B.B

    double rt = AB * AB + B2 * (C2 * c * c - 2 * CA * c);
    if (verbose)
      fprintf(stderr, "   rt = %.15lf\n", rt);

    // if root is 0 or negative continue
    if (rt < -sqrt(epsilon)) {
      if (verbose)
        fprintf(stderr, "   rejecting negative root: %.15lf\n\n", rt);
      continue;
    }

    rt = (rt > sqrt(epsilon)) ? sqrt(rt) : 0;

    double b = (AB + rt) / B2;

    int num_factors = factors.size();

    // AB and AC should be almost the same length, if not reject
    // RK: was 1e-3 (.001), but with reversing B and C method, we can reject
    // almost all cases outside epsilon
    if (fabs((b * B - A).len() - (b * B - c * C).len()) < epsilon) {
      double AB_len = (b * B - A).len();
      if (verbose)
        fprintf(stderr, "   |AbB| = %.15lf, |AcC| = %.15lf, |bBcC| = %.15lf\n",
                1.0, (c * C - A).len() / AB_len,
                (c * C - b * B).len() / AB_len);
      factors.push_back(c);
      factors.push_back(b);
      if (verbose)
        fprintf(stderr, "   b = %.15lf\n\n", factors.back());
    }

    if (rt != 0) {
      b = (AB - rt) / B2;
      // AB and AC should be almost the same length, if not reject
      // RK: was 1e-3 (.001), but with reversing B and C method, we can reject
      // almost all cases outside epsilon
      if (fabs((b * B - A).len() - (b * B - c * C).len()) < epsilon) {
        double AB_len = (b * B - A).len();
        if (verbose)
          fprintf(stderr,
                  "   |AbB| = %.15lf, |AcC| = %.15lf, |bBcC| = %.15lf\n", 1.0,
                  (c * C - A).len() / AB_len, (c * C - b * B).len() / AB_len);
        factors.push_back(c);
        factors.push_back(b);
        if (verbose)
          fprintf(stderr, "   b = %.15lf\n\n", b);
      }
    }
    else if (verbose)
      fprintf(stderr, "   Double root\n\n");

    if (factors.size() - num_factors == 4) {
      if (verbose)
        fprintf(stderr,
                "   ***** DOUBLE VALUE for non-double root!!!! *****\n\n");
      // factors.resize(num_factors);
    }
  }

  if (verbose)
    fprintf(stderr, "total of %d possible triangles\n\n",
            (int)factors.size() / 2);

  return num_roots;
}

// find closest match to alpha, beta and gamma values in the table
// if found, use the formula derived values. If not found, table values
// unchanged
bool refine_abg(/* const string &sym_type, */ bool reverse, const bool verbose,
                const Vec3d &A, const Vec3d &B, const Vec3d &C, double &alpha,
                double &beta, double &gamma)
{
  bool found = false;

  /* geom2 is for debug OFF files
     // schwarz triangle and Shephards transformation
     Geometry geom2;
     make_triangle(geom2, A, B, C, Color(0.0,0.0,1.0));
     make_triangle(geom2, A*alpha, B*beta, C*gamma, Color(1.0,1.0,0.0));
     if (sym_type[0]=='D' || sym_type[0]=='I')
        geom2.transform(Trans3d::rotate(0, M_PI/2, 0));
     geom2.write("tri1.off");
     geom2.clear_all();
  */

  vector<double> factors;
  find_equ_tris(A, (reverse ? C : B), (reverse ? B : C), factors, verbose);

  /*
     // schwarz triangle and calculated values
     make_triangle(geom2, A, B, C, Color(0.0,0.0,1.0));
  */

  for (unsigned int k = 0; k < factors.size() / 2; k++) {
    double a = 1;
    double b = factors[2 * k + (reverse ? 0 : 1)];
    double g = factors[2 * k + (reverse ? 1 : 0)];

    if (verbose) {
      fprintf(stderr, "calculated from formula (%d):\n", k + 1);
      fprintf(stderr, "alpha: % .15lf\n", a);
      fprintf(stderr, "beta:  % .15lf\n", b);
      fprintf(stderr, "gamma: % .15lf\n\n", g);
    }

    Geometry triangle;
    make_triangle(triangle, A * a, B * b, C * g, Color(1.0, 0.0, 0.0));

    // unit edges
    GeometryInfo info(triangle);
    if (info.num_iedges() > 0) {
      double val = info.iedge_length_lims().sum / info.num_iedges();
      triangle.transform(Trans3d::scale(1 / val));
    }
    /*
          geom2.append(triangle);
    */

    // unitized a,b,g
    int i = 0;
    int j = 0;

    Vec3d t = triangle.verts()[0];
    for (i = 0; i < 3; i++)
      if (double_ne(t[i], 0, epsilon))
        break;
    for (j = 0; j < 3; j++)
      if (double_ne(A[i], 0, epsilon))
        break;
    a = t[i] / A[j];

    t = triangle.verts()[1];
    for (i = 0; i < 3; i++)
      if (double_ne(t[i], 0, epsilon))
        break;
    for (j = 0; j < 3; j++)
      if (double_ne(B[j], 0, epsilon))
        break;
    b = t[i] / B[j];

    t = triangle.verts()[2];
    for (i = 0; i < 3; i++)
      if (double_ne(t[i], 0, epsilon))
        break;
    for (j = 0; j < 3; j++)
      if (double_ne(C[j], 0, epsilon))
        break;
    g = t[i] / C[j];

    if (verbose) {
      fprintf(stderr, "unit edge values:\n");
      fprintf(stderr, "alpha: % .15lf\n", a);
      fprintf(stderr, "beta:  % .15lf\n", b);
      fprintf(stderr, "gamma: % .15lf\n\n", g);
    }

    // hard code epsilon value since hard coded table values are only 15 places
    // compare to table
    double eps = 1e-12;
    if (double_eq(alpha, a, eps) && double_eq(beta, b, eps) &&
        double_eq(gamma, g, eps)) {
      if (verbose)
        fprintf(stderr, "***** found match!! *****\n\n");
      found = true;

      alpha = a;
      beta = b;
      gamma = g;
    }

    /*
          // make polyhedron
          Geometry poly;
          if (sym_type[0]=='D' || sym_type[0]=='I')
             triangle.transform(Trans3d::rotate(0, M_PI/2, 0));
          sym_repeat(poly, triangle, sym_type);
          sort_merge_elems(poly, "vef", epsilon);
          char filename[80];
          sprintf(filename,"triangle%d_poly.off",k+1);
          poly.write(filename);
    */
  }

  /*
     if (sym_type[0]=='D' || sym_type[0]=='I')
        geom2.transform(Trans3d::rotate(0, M_PI/2, 0));
     geom2.write("tri2.off");
  */

  return found;
}

void make_poly(Geometry &geom, const string &sym_type, const bool triangle_only,
               const bool verbose, const Vec3d &A, const Vec3d &B,
               const Vec3d &C, const double alpha, const double beta,
               const double gamma)
{
  double a = alpha;
  double b = beta;
  double g = gamma;

  bool found = refine_abg(/* sym_type, */ false, verbose, A, B, C, a, b, g);
  if (verbose)
    fprintf(stderr, "1st try: match %sfound\n\n", (found ? "" : "not "));
  if (!found) {
    found = refine_abg(/* sym_type, */ true, verbose, A, B, C, a, b, g);
    if (verbose)
      fprintf(stderr, "2nd try: match %sfound%s\n\n", (found ? "" : "not "),
              (found ? "" : " (using table values)"));
  }

  if (verbose)
    verbose_output(A, B, C, a, b, g);

  // make poly called twice from dipyramid code, don't remake triangle
  if (geom.verts().size() == 0)
    make_triangle(geom, A * a, B * b, C * g, Color());

  if (sym_type[0] == 'D' || sym_type[0] == 'I')
    geom.transform(Trans3d::rotate(0, M_PI / 2, 0));

  if (!triangle_only) {
    sym_repeat(&geom, geom, sym_type);
    merge_coincident_elements(&geom, "vef", epsilon);
  }
}

void make_delta_dipyramid(Geometry &geom, const int n, const int d,
                          bool triangle_only = false, bool verbose = false)
{
  char buf1[MSG_SZ];
  char buf2[MSG_SZ];

  sprintf(buf1, "D%dh", n);
  sprintf(buf2, "[%d/%d,1/2,1/2]", d, n);
  fprintf(stderr, "Dihedral Group:  %s  %s  %d/%d %sdipyramid\n", buf1, buf2, n,
          d, ((d == 1) ? "" : "star "));

  string sym_type = buf1;

  Vec3d A = Vec3d(1, 0, 0);
  Vec3d B = Vec3d(cos((2.0 * M_PI * d) / n), sin((2.0 * M_PI * d) / n), 0);
  Vec3d C = Vec3d(0, 0, 1);
  double alpha = 1 / (2.0 * sin((M_PI * d) / n));
  double beta = alpha;
  double gamma = sqrt(1.0 - alpha * alpha);

  // force triangle_only, verbose to false
  make_poly(geom, sym_type, true, (triangle_only ? verbose : false), A, B, C,
            alpha, beta, gamma);

  if (!triangle_only)
    make_poly(geom, sym_type, false, verbose, A, B, C, alpha, beta, gamma);
  else if (verbose)
    verbose_output(A, B, C, alpha, beta, gamma);
}

void tet_to_dihedral(Geometry &geom, const string &sym_from, const int k,
                     Trans3d pos = Trans3d())
{
  char sym_to[MSG_SZ];
  sprintf(sym_to, "D%d%s", k, ((k % 4 == 0) ? "h" : "v"));
  transform_and_repeat(&geom, sym_to, sym_from, pos);
}

void case_a_star_tetrahedron(Geometry &geom, const int k)
{
  geom.read_resource("u1");
  tet_to_dihedral(geom, "Td", 2 * k);
}

void case_b_5_or_10_tetrahedra(Geometry &geom, double angle, const int k)
{
  if (angle == INFINITY)
    angle = 0;

  geom.read_resource("u1");

  // to construct in one statement for Ih
  // transform_and_repeat(&geom, (k == 1 ? "I" : "Ih"), "Td",
  // Trans3d::rotate(0,angle,0));

  if (k == 1)
    transform_and_repeat(&geom, "I", "Td", Trans3d::rotate(0, angle, 0));
  else if (k == 2) {
    transform_and_repeat(&geom, "Oh", "Td", Trans3d::rotate(0, angle, 0));
    transform_and_repeat(&geom, "I", "Oh");
  }
}

void case_c_2_dipyramids(Geometry &geom, double angle, const int n, const int d)
{
  if (angle == INFINITY) {
    angle = (M_PI / 2) / n; // 90/n degrees
    fprintf(stderr, "angle calculated is %g\n", rad2deg(angle));
  }

  fprintf(stderr, "Using: ");
  make_delta_dipyramid(geom, n, d);

  char sym_from[MSG_SZ];
  sprintf(sym_from, "D%dh", n);
  char sym_to[MSG_SZ];
  sprintf(sym_to, "D%dh", 2 * n);

  // advance angle so that angle = 0 is coincident constituents
  geom.transform(Trans3d::rotate(0, 0, angle + M_PI / (2 * n)));
  transform_and_repeat(&geom, sym_to, sym_from);
}

void case_d_6_octahedra(Geometry &geom, double angle)
{
  if (angle == INFINITY)
    angle = (M_PI / 8); // 22.5 degrees

  geom.read_resource("u5");

  // at 0 degrees, produced 3 coincident octahedra
  transform_and_repeat(&geom, "D2h", "Oh", Trans3d::rotate(0, 0, angle));
  transform_and_repeat(&geom, "T", "D2h");
}

void case_e_4_or_8_triangular_dipyramids(Geometry &geom, double angle,
                                         const int k)
{
  if (angle == INFINITY)
    angle = 0;

  fprintf(stderr, "Using: ");
  make_delta_dipyramid(geom, 3, 1);

  // to construct in one statement for Oh
  // transform_and_repeat(&geom, (k == 1 ? "O" : "Oh"), "D3h",
  //    Trans3d::rotate(Vec3d(0,0,1),Vec3d(1,1,1)) *
  //    Trans3d::rotate(0,0,angle+M_PI/12));

  if (k == 1)
    transform_and_repeat(&geom, "O", "D3h",
                         Trans3d::rotate(Vec3d(0, 0, 1), Vec3d(1, 1, 1)) *
                             Trans3d::rotate(0, 0, angle + M_PI / 12));
  else if (k == 2) {
    transform_and_repeat(&geom, "D6h", "D3h");
    transform_and_repeat(&geom, "O", "D6h",
                         Trans3d::rotate(Vec3d(0, 0, 1), Vec3d(1, 1, 1)) *
                             Trans3d::rotate(0, 0, angle + M_PI / 12));
  }
}

void case_f_6_or_12_pentagonal_dipyramids(Geometry &geom, double angle,
                                          const int k)
{
  if (angle == INFINITY)
    angle = 0;

  fprintf(stderr, "Using: ");
  if (k == 1 || k == 2)
    make_delta_dipyramid(geom, 5, 1);
  else if (k == 3 || k == 4)
    make_delta_dipyramid(geom, 5, 2);

  // to construct in one statement for Ih
  // transform_and_repeat(&geom, ((k == 1 || k == 3) ? "I" : "Ih"), "D5h",
  //    Trans3d::rotate(Vec3d(0,0,1),Vec3d(0,1,phi)) *
  //    Trans3d::rotate(0,0,angle+M_PI/5));

  if (k == 1 || k == 3)
    transform_and_repeat(&geom, "I", "D5h",
                         Trans3d::rotate(Vec3d(0, 0, 1), Vec3d(0, 1, phi)) *
                             Trans3d::rotate(0, 0, angle + M_PI / 5));
  else if (k == 2 || k == 4) {
    transform_and_repeat(&geom, "D10h", "D5h");
    transform_and_repeat(&geom, "I", "D10h",
                         Trans3d::rotate(Vec3d(0, 0, 1), Vec3d(0, 1, phi)) *
                             Trans3d::rotate(0, 0, angle + M_PI / 5));
  }
}

void case_g_2_tetrahedra(Geometry &geom, double angle)
{
  if (angle == INFINITY)
    angle = (M_PI / 4); // 45 degrees

  geom.read_resource("u1");

  // advance angle so that angle = 0 is coincident constituents
  geom.transform(Trans3d::rotate(0, 0, angle));
  tet_to_dihedral(geom, "S4", 2); // 2*k=4
}

void case_h_2k_tetrahedra(Geometry &geom, double angle, const int k)
{
  if (angle == INFINITY)
    angle = deg2rad(1.0);

  case_g_2_tetrahedra(geom, angle);
  tet_to_dihedral(geom, "D2v", 2 * k);
}

void case_i_6_tetrahedra(Geometry &geom, double angle)
{
  if (angle == INFINITY)
    angle = (M_PI / 4); // 45 degrees

  geom.read_resource("u1");

  transform_and_repeat(&geom, "T", "T", Trans3d::rotate(0, 0, angle));
}

void case_j_12_tetrahedra(Geometry &geom, double angle)
{
  if (angle == INFINITY)
    angle = (M_PI / 6); // 30 degrees;

  geom.read_resource("u1");

  transform_and_repeat(&geom, "S2", "T");
  transform_and_repeat(&geom, "T", "Oh", Trans3d::rotate(0, 0, angle));
}

void case_k_2k_dipyramids(Geometry &geom, double angle, const int k,
                          const int n, const int d)
{
  if (angle == INFINITY)
    angle = deg2rad(1.0);

  fprintf(stderr, "Using: ");
  make_delta_dipyramid(geom, n, d);

  char sym_from[MSG_SZ];
  sprintf(sym_from, "D%dh", n);
  char sym_to[MSG_SZ];
  sprintf(sym_to, "D%dh", k * n);

  transform_and_repeat(&geom, sym_from, sym_from,
                       Trans3d::rotate(0, 0, angle + M_PI));
  transform_and_repeat(&geom, sym_to, sym_from);
}

void case_l_k_dipyramids(Geometry &geom, const int k, const int n, const int d)
{
  fprintf(stderr, "Using: ");
  make_delta_dipyramid(geom, n, d);

  char sym_from[MSG_SZ];
  sprintf(sym_from, "D%dh", n);
  char sym_to[MSG_SZ];
  sprintf(sym_to, "D%dh", k * n);

  transform_and_repeat(&geom, sym_to, sym_from);
}

void case_m_10_or_20_triangular_dipyramids(Geometry &geom, double angle,
                                           const int k)
{
  if (angle == INFINITY)
    angle = 0;

  fprintf(stderr, "Using: ");
  make_delta_dipyramid(geom, 3, 1);

  // to construct in one statement for Ih
  // transform_and_repeat(&geom, (k == 1 ? "I" : "Ih"), "D3h",
  //   Trans3d::rotate(Vec3d(0,0,1), Vec3d(1/phi,0,phi)) *
  //   Trans3d::rotate(0,0,angle+M_PI/6));

  if (k == 1)
    transform_and_repeat(
        &geom, "I", "D3h",
        Trans3d::rotate(Vec3d(0, 0, 1), Vec3d(1 / phi, 0, phi)) *
            Trans3d::rotate(0, 0, angle + M_PI / 6));
  else if (k == 2) {
    transform_and_repeat(&geom, "D6h", "D3h");
    transform_and_repeat(
        &geom, "I", "D6h",
        Trans3d::rotate(Vec3d(0, 0, 1), Vec3d(1 / phi, 0, phi)) *
            Trans3d::rotate(0, 0, angle + M_PI / 6));
  }
}

void case_n_6_10_3_star_dipyramids(Geometry &geom, double angle)
{
  if (angle == INFINITY)
    angle = 0;

  fprintf(stderr, "Using: ");
  make_delta_dipyramid(geom, 10, 3);

  transform_and_repeat(&geom, "I", "D10h",
                       Trans3d::rotate(Vec3d(0, 0, 1), Vec3d(0, 1, phi)) *
                           Trans3d::rotate(0, 0, angle + M_PI / 10));
}

void case_o_5_or_10_augmented_tetrahedra(Geometry &geom, double angle,
                                         const int k)
{
  if (angle == INFINITY)
    angle = 0;

  // to construct in one statement for Ih
  // transform_and_repeat(&geom, (k == 1 ? "I" : "Ih"), "Td",
  // Trans3d::rotate(0,angle,0));

  if (k == 1)
    transform_and_repeat(&geom, "I", "Td", Trans3d::rotate(0, angle, 0));
  else if (k == 2) {
    transform_and_repeat(&geom, "Oh", "Td", Trans3d::rotate(0, angle, 0));
    transform_and_repeat(&geom, "I", "Oh");
  }
}

void case_p_5_augmented_octahedra(Geometry &geom, double angle)
{
  if (angle == INFINITY)
    angle = 0;

  transform_and_repeat(&geom, "I", "Oh", Trans3d::rotate(0, angle, 0));
}

void case_q_5_excavated_octahedra(Geometry &geom, double angle)
{
  if (angle == INFINITY)
    angle = 0;

  transform_and_repeat(&geom, "I", "Oh", Trans3d::rotate(0, angle, 0));
}

void compound_Coloring(Geometry &geom, const char Coloring_method,
                       const ColorMapMulti &map, const int face_opacity)
{
  // color by sub-symmetry  as map indexes happened by default in sym_repeat()
  if (!Coloring_method) {
    // no color, strip colors
    geom.colors(FACES).clear();
    geom.clear(EDGES);
    geom.colors(VERTS).clear();
  }
  else {
    Coloring clrng;
    clrng.add_cmap(map.clone());
    clrng.set_geom(&geom);

    if (Coloring_method == 'c') {
      // color by constituents
      clrng.f_parts(true);
    }
    else if (Coloring_method == 's') {
      Symmetry sym;
      vector<vector<set<int>>> sym_equivs;
      sym.init(geom, &sym_equivs);
      clrng.f_sets(sym_equivs[2], true);
    }

    // blend edges
    geom.add_missing_impl_edges();
    clrng.e_face_color();
    clrng.v_face_color();

    // transparency
    if (face_opacity != 255) {
      for (unsigned int i = 0; i < geom.faces().size(); i++) {
        Color col = geom.colors(FACES).get(i);
        if (col.is_value())
          col = Color(col[0], col[1], col[2], face_opacity);
        geom.colors(FACES).set(i, col);
      }
    }
  }
}

int main(int argc, char *argv[])
{
  id_opts opts;
  opts.process_command_line(argc, argv);
  id_poly id_polys;

  if (opts.list_polys) {
    id_polys.list_polys();
    exit(0);
  }

  Geometry geom;
  if (opts.make_dipyramid)
    make_delta_dipyramid(geom, opts.n, opts.d, opts.triangle_only,
                         opts.verbose);
  else if (opts.case_type.length()) {
    if (opts.case_type == "a")
      case_a_star_tetrahedron(geom, opts.k);
    else if (opts.case_type == "b")
      case_b_5_or_10_tetrahedra(geom, opts.angle, opts.s);
    else if (opts.case_type == "c")
      case_c_2_dipyramids(geom, opts.angle, opts.n, opts.d);
    else if (opts.case_type == "d")
      case_d_6_octahedra(geom, opts.angle);
    else if (opts.case_type == "e")
      case_e_4_or_8_triangular_dipyramids(geom, opts.angle, opts.s);
    else if (opts.case_type == "f")
      case_f_6_or_12_pentagonal_dipyramids(geom, opts.angle, opts.s);
    else if (opts.case_type == "g")
      case_g_2_tetrahedra(geom, opts.angle);
    else if (opts.case_type == "h")
      case_h_2k_tetrahedra(geom, opts.angle, opts.k);
    else if (opts.case_type == "i")
      case_i_6_tetrahedra(geom, opts.angle);
    else if (opts.case_type == "j")
      case_j_12_tetrahedra(geom, opts.angle);
    else if (opts.case_type == "k")
      case_k_2k_dipyramids(geom, opts.angle, opts.k, opts.n, opts.d);
    else if (opts.case_type == "l")
      case_l_k_dipyramids(geom, opts.k, opts.n, opts.d);
    else if (opts.case_type == "m")
      case_m_10_or_20_triangular_dipyramids(geom, opts.angle, opts.s);
    else if (opts.case_type == "n")
      case_n_6_10_3_star_dipyramids(geom, opts.angle);
    else if (opts.case_type == "o") {
      int sym_no = 3;
      string sym_type = id_polys.get_sym_type(sym_no);
      make_poly(geom, sym_type, false, false, id_polys.A(sym_no),
                id_polys.B(sym_no), id_polys.C(sym_no), id_polys.alpha(sym_no),
                id_polys.beta(sym_no), id_polys.gamma(sym_no));
      case_o_5_or_10_augmented_tetrahedra(geom, opts.angle, opts.s);
    }
    else if (opts.case_type == "p") {
      int sym_no = 14;
      string sym_type = id_polys.get_sym_type(sym_no);
      make_poly(geom, sym_type, false, false, id_polys.A(sym_no),
                id_polys.B(sym_no), id_polys.C(sym_no), id_polys.alpha(sym_no),
                id_polys.beta(sym_no), id_polys.gamma(sym_no));
      case_p_5_augmented_octahedra(geom, opts.angle);
    }
    else if (opts.case_type == "q") {
      int sym_no = 15;
      string sym_type = id_polys.get_sym_type(sym_no);
      make_poly(geom, sym_type, false, false, id_polys.A(sym_no),
                id_polys.B(sym_no), id_polys.C(sym_no), id_polys.alpha(sym_no),
                id_polys.beta(sym_no), id_polys.gamma(sym_no));
      case_q_5_excavated_octahedra(geom, opts.angle);
    }
  }
  else {
    int sym_no = id_polys.lookup_sym_no(opts.poly);
    if (sym_no >= id_polys.get_last_iso_delta())
      opts.error("polyhedron number '" + opts.poly + "' out of range");
    if (sym_no < 0)
      opts.error("unknown polyhedron '" + opts.poly + "'");

    id_polys.list_poly(sym_no);

    // patch for 9 and 14. Made with make_poly they will have merged vertices
    // between constituents
    // make parts and use transform_and_repeat on results
    if (sym_no + 1 == 9 || sym_no + 1 == 14) {
      geom.clear_all();
      if (sym_no + 1 == 9) {
        make_delta_dipyramid(geom, 8, 3, opts.triangle_only, opts.verbose);
        if (!opts.triangle_only)
          transform_and_repeat(&geom, "Oh", "D8h");
      }
      else if (sym_no + 1 == 14) {
        sym_no = 2; // consituent of case 14, sent to make_poly
        string sym_type = id_polys.get_sym_type(sym_no);
        make_poly(geom, sym_type, opts.triangle_only, opts.verbose,
                  id_polys.A(sym_no), id_polys.B(sym_no), id_polys.C(sym_no),
                  id_polys.alpha(sym_no), id_polys.beta(sym_no),
                  id_polys.gamma(sym_no));
        if (!opts.triangle_only)
          transform_and_repeat(&geom, "Oh", "Td");
      }
      geom.transform(Trans3d::rotate(
          0, M_PI / 2,
          0)); // as it did in make_poly for same color order as before
    }
    // process all the rest as normal
    else {
      string sym_type = id_polys.get_sym_type(sym_no);
      make_poly(geom, sym_type, opts.triangle_only, opts.verbose,
                id_polys.A(sym_no), id_polys.B(sym_no), id_polys.C(sym_no),
                id_polys.alpha(sym_no), id_polys.beta(sym_no),
                id_polys.gamma(sym_no));
    }
  }

  // orient for positive volume
  geom.orient();

  compound_Coloring(geom, opts.Coloring_method, opts.map, opts.face_opacity);

  opts.write_or_error(geom, opts.ofile);

  return 0;
}
