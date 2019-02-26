/*
   Copyright (c) 2016, Adrian Rossiter

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
   Name: iso_kite.cc
   Description: make isohedral kite-faced polyhedra
   Project: Antiprism - http://www.antiprism.com
*/

#include <algorithm>
#include <map>
#include <string>
#include <vector>

#include "../base/antiprism.h"

using std::map;
using std::string;
using std::vector;

using namespace anti;

// clang-format off
const char *base_models[] = { "T1", "T2", "O1", "O2", "O2B", "I1", "I2", "I3",
   "I4", "I4B", "I5", "I5B", "I6", "I6B", "I6C", "I7", "I7B", "I8", "I8B",
   "I9", "I9B", "I10", "n/d", "end" };

const char *compounds[][4] = {
   { "T1"   , "Oh", "Td",    ""},      //excluding A=B uc4
   { "T1"   , "I" , "T" ,    ""},      //excluding A=B, which are Th uc5
   { "T1"   , "Ih", "T" ,    ""},      //excluding A=B, which are Th uc6
   { "T2"   , "Oh", "Td",    ""},      //excluding A=B uc4
   { "T2"   , "I" , "T" ,    ""},      //excluding A=B, which are Th uc5
   { "T2"   , "Ih", "T" ,    ""},      //excluding A=B, which are Th uc6

   { "O1"   , "Ih", "Th",    ""},      // uc17
   { "O2B"  , "Ih", "Th",    ""},      // uc17

   { "3"    , "Th", "S6",    ""},      //uc10
   { "3"    , "Oh", "S6",    ""},      // uc11
   { "3"    , "Oh", "D3v",   ""},      // uc12
   { "3"    , "Oh", "D3v,1", ""},      // uc12
   { "6/2"  , "Oh", "D3v",   "1"},     // uc38

   { "3"    , "Ih", "S6",    ""},      // uc13
   { "3"    , "Ih", "D3v,1", ""},      // uc15
   { "3"    , "Ih", "D3v",   ""},      // uc16
   { "6/2"  , "Ih", "D3v",   "1"},     // uc39

   { "3/2"  , "O" , "D3",    ""},      // uc30
   { "3/2"  , "Oh", "D3",    ""},      // uc31
   { "3/2"  , "I" , "D3",    ""},      // uc32
   { "3/2"  , "Ih", "D3",    ""},      // uc33

   { "4"    , "O" , "D4",    ""},      // uc42
   { "4"    , "Oh", "D4",    ""},      // uc43

   { "4/3"  , "O" , "D4",    ""},      // uc42
   { "4/3"  , "Oh", "D4",    ""},      // uc43


   { "5"    , "Ih", "S10",   ""},      // uc26
   { "5"    , "Ih", "D5v",   ""},      // uc27
   { "5"    , "Ih", "D5v,1", ""},      // uc27
   { "10/2" , "Ih", "D5v",   "1"},     // uc40

   { "5/2"  , "I" , "D5",    ""},      // uc44
   { "5/2"  , "Ih", "D5",    ""},      // uc45

   { "5/3"  , "Ih", "S10",   ""},      // uc28
   { "5/3"  , "Ih", "D5v",   ""},      // uc29
   { "5/3"  , "Ih", "D5v,1", ""},      // uc29
   { "10/6" , "Ih", "D5v",   "1"},     // uc41

   { "5/4"  , "I" , "D5",    ""},      // uc44
   { "5/4"  , "Ih", "D5",    ""},      // uc45


   { "n/d_odd" , "DN_odd*nv" , "S2n", ""},    // uc22
   { "n/d_odd" , "DN_even*nh", "S2n", ""},    // uc22
   { "n/d_odd" , "DN_odd*nv" , "Dnv", ""},    // uc23
   { "n/d_odd" , "DN_even*nh", "Dnv", ""},    // uc23
   { "n/d_even", "DN*nh"     , "Cnh", ""},    // uc24
   { "n/d_even", "DN*nh"     , "Dnh", ""},    // uc25
   { "end"     ,  "end"      , "end", "end"}  // end
};
// clang-format on

Status read_fraction(const char *frac_str, int &num, int &denom)
{
  char frac_str_cpy[MSG_SZ];
  strncpy(frac_str_cpy, frac_str, MSG_SZ);

  Status stat;
  denom = 1;
  char *p = strchr(frac_str_cpy, '/');
  if (p != nullptr) {
    *p++ = '\0';
    if (!(stat = read_int(p, &denom)))
      return Status::error(msg_str("denominator, %s", stat.c_msg()));
  }
  if (!(stat = read_int(frac_str_cpy, &num)))
    return Status::error(msg_str("numerator, %s", stat.c_msg()));

  return Status::ok();
}

Status read_triangle_fraction(char *frac_str, int &num, int &denom)
{
  Status stat = read_fraction(frac_str, num, denom);
  if (stat.is_error())
    return stat;

  if (num < 2)
    return Status::error("numerator must be 2 or greater");

  if (denom == 0)
    return Status::error("denominator cannot be 0");

  if (denom >= num)
    return Status::error("denominator must be less than the numerator");

  return Status::ok();
}

class kt_opts : public ProgramOpts {
private:
  void read_vert_height(int idx, char c);

public:
  vector<int> triangle;
  string model_name;
  int num_fracs;
  unsigned int heights_set;
  double heights[3];
  double angle;
  int num_parts;
  int color_type;
  int list_idx; // -1: ignore, 0: print list, >0: make model with index
  bool kite_only;
  int verb; // verbosity - 0:no report, 1:print report

  string ofile;

  kt_opts()
      : ProgramOpts("iso_kite"), heights_set(0), angle(NAN), num_parts(0),
        color_type(-1), list_idx(-1), kite_only(false), verb(1)
  {
    for (double &height : heights)
      height = 1.0;
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

// clang-format off
void kt_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] model_args...\n"
"\n"
"Create all isohedral polyhedra and polyhedron compounds whose face type\n"
"is a kite (or dart)\n"
"\n"
"Model_args can be a Schwarz triangle name (which can also be given as two\n"
"or three fractions A, B, C (default 2) which specify the Schwarz triangle\n"
"vertices, where n/d corresponds to an angle PI*d/n), optionally followed by\n"
"a permutation letter from AaBbCc, to rotate the fraction A, B, or C to the\n"
"front (and if lowercase then follow by swapping the last two fractions)\n"
"  T1:  (3 3 2)         T2:  (3 3 3/2)\n"
"  O1:  (4 3 2)         O2:  (4/3 4 3)\n"
"  I1:  (5 3 2)         I2:  (5/2 3 2)\n"
"  I3:  (5/2 5 2)       I4:  (5/2 3 3)\n"
"  I5:  (5/4 3 3)       I6:  (5/3 5 3)\n"
"  I7:  (5/4 5 3)       I8:  (5/3 5/2 3)\n"
"  I9:  (5/4 5 5)       I10: (5/2 5/2 5/2)\n"
"If model_args is a single fraction then make a trapezohedron (with B and C\n"
"as antiprism vertices).\n"
"The kite has vertices along A, C, B, and C reflected in AB. The height of up\n"
"to two vertices may be given, and the third is calculated (B and C have\n"
"the same height in a trapezohedron). If less than two heights are given\n"
"the calculated vertex will be B, if not specified, otherwise A. The kite\n"
"is repeated with the symmetry corresponding to the base Schwarz triangle,\n"
"or trapezohedron model.\n"
"Notes: the tetrahedral triangles T1 and T2 both produce a model with\n"
"octahedral symmetry when height A=-B (cube and rhombic dodecahedron,\n"
"respectively). These models do not form an octahedral compound, and have\n"
"common symmetry subgroup of Th in the icosahedral compounds. The 3\n"
"trapezohedron also has octahedral symmetry at height A=+/-B (cube), and\n"
"likewise has a common symmetry subgroup of Th in icosahedral compounds.\n"
"\n"
"Options\n"
"  -h        this help message\n"
"  -l num    'list' will list all compounds of the specified base model, or\n"
"            all models and compounds if no base model specified. The list\n"
"            includes a model number and base model, and if a compound, the\n"
"            alignment symmetry (possibly with realignment number) and final\n"
"            symmetry (and possibly a final pre-realignment number).\n"
"            Specifying a number will make the the model with that number\n"
"            in the list.\n"
"  -A ht     height of kite apex on OA (default 1.0 or calculated)\n"
"  -B ht     height of kite apex on OB (only magnitude used if trapezohedron,\n"
"            default 1.0 or calculated)\n"
"  -C ht     height of kite side vertex on OC (Schwarz triangle only,\n"
"            default 1.0 or calculated)\n"
"  -a ang    angle, for compounds where this is a parameter\n"
"  -N parts  number of parts, for compounds where this is a parameter\n"
"  -c type   colour the faces around each vertex of a type, from AaBbCc\n"
"            (colouring by value/index for upper/lower case) using a\n"
"            different colour for each set (Schwarz models only), or Kk\n"
"            to colour by base component (which may be compound)\n"
"  -k        output a single kite (colours not applied)\n"
"  -q        quiet, don't print final report\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name());
}
// clang-format on

void kt_opts::read_vert_height(int idx, char c)
{
  heights_set |= (1 << idx);
  if (heights_set == 7)
    error("can only set two of options A, B, C", c);
  print_status_or_exit(read_double(optarg, &heights[idx]), c);
}

int tri_from_str(const char *str, vector<int> &fracs, char *errmsg = nullptr)
{
  if (errmsg)
    *errmsg = '\0';
  char str_cpy[MSG_SZ];
  strncpy(str_cpy, str, MSG_SZ);

  char sym_char;
  int tri_idx;
  char perm;
  char buf;
  int num = sscanf(str_cpy, "%c%d%c%c", &sym_char, &tri_idx, &perm, &buf);
  if (num < 2 || num > 3) {
    if (errmsg)
      sprintf(errmsg, "invalid format for triangle description");
    return -1; // -1 for skipping if multiple formats considered
  }

  if (!strchr("TOI", sym_char)) {
    if (errmsg)
      sprintf(errmsg, "first letter is not T, O or I");
    return 1; // fail
  }
  if (tri_idx < 1 || (sym_char == 'T' && tri_idx > 2) ||
      (sym_char == 'O' && tri_idx > 2) || (sym_char == 'I' && tri_idx > 10)) {
    if (errmsg)
      sprintf(errmsg, "number out of range for triangle type");
    return 1; // fail
  }

  map<string, int> axes;
  axes["T1"] = 0;
  axes["T2"] = 19;
  axes["O1"] = 2;
  axes["O2"] = 24;
  axes["I1"] = 4;
  axes["I2"] = 5;
  axes["I3"] = 15;
  axes["I4"] = 20;
  axes["I5"] = 21;
  axes["I6"] = 25;
  axes["I7"] = 26;
  axes["I8"] = 27;
  axes["I9"] = 40;
  axes["I10"] = 41;

  string base_str = msg_str("%c%d", sym_char, tri_idx);
  auto mi = axes.find(base_str);
  if (mi == axes.end()) {
    if (errmsg)
      sprintf(errmsg, "triangle not found (internal error)");
    return 1; // fail
  }

  get_schwarz_tri_fracs(mi->second, fracs);
  // reverse triangle so that any dihedral axis is last
  std::swap(fracs[0], fracs[4]);
  std::swap(fracs[1], fracs[5]);

  if (num == 3) { // permutation specified
    if (!strchr("ABCabc", perm)) {
      if (errmsg)
        sprintf(errmsg, "permutation letter is not from ABCabc");
      return 1; // fail
    }
    if (strchr("Bb", perm))
      std::rotate(fracs.begin(), fracs.begin() + 2, fracs.end());
    if (strchr("Cc", perm))
      std::rotate(fracs.begin(), fracs.begin() + 4, fracs.end());
    if (strchr("abc", perm)) {
      std::swap(fracs[2], fracs[4]);
      std::swap(fracs[3], fracs[5]);
    }
  }

  return 0; // success
}

// convert model name to fractions, return number of fractions or 0
int args2fracs(const char *model, vector<int> &fracs, char *errmsg = nullptr)
{
  if (errmsg)
    *errmsg = '\0';
  fracs.resize(6);
  if (*model == '\0') {
    if (errmsg)
      strcpy_msg(errmsg, "model name is blank");
    return 0;
  }
  char errmsg2[MSG_SZ];
  char model_cpy[MSG_SZ];
  strncpy(model_cpy, model, MSG_SZ);
  vector<char *> parts;
  int num_parts = split_line(model_cpy, parts);

  if (num_parts == 1) {
    int ret = tri_from_str(model, fracs, errmsg2);
    if (errmsg)
      strcpy_msg(errmsg, errmsg2);
    if (ret == 0) // success
      return 3;
    else if (ret == 1) // failure
      return 0;

    // else ret==-1, invalid format so ignore
  }

  Status stat;
  if (num_parts > 0) {
    if (!(stat = read_triangle_fraction(parts[0], fracs[0], fracs[1]))) {
      if (errmsg)
        snprintf(errmsg, MSG_SZ, "fraction A is '%s': %s", parts[0],
                 stat.c_msg());
      return 0;
    }
  }
  if (num_parts > 1) {
    if (!(stat = read_triangle_fraction(parts[1], fracs[2], fracs[3]))) {
      if (errmsg)
        snprintf(errmsg, MSG_SZ, "fraction B is '%s': %s", parts[1],
                 stat.c_msg());
      return 0;
    }
    if (num_parts > 2) {
      if (!(stat = read_triangle_fraction(parts[2], fracs[4], fracs[5]))) {
        if (errmsg)
          snprintf(errmsg, MSG_SZ, "fraction C is '%s': %s", parts[2],
                   stat.c_msg());
        return 0;
      }
    }
    else {
      fracs[4] = 2;
      fracs[5] = 1;
    }
  }
  return (num_parts == 1) ? 1 : 3;
}

void kt_opts::process_command_line(int argc, char **argv)
{
  char errmsg[MSG_SZ];
  opterr = 0;
  int c;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hA:B:C:a:N:c:l:kqo:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'A':
      read_vert_height(0, c);
      break;

    case 'B':
      read_vert_height(1, c);
      break;

    case 'C':
      read_vert_height(2, c);
      break;

    case 'a':
      print_status_or_exit(read_double(optarg, &angle), c);
      // angle = deg2rad(angle);
      break;

    case 'N':
      print_status_or_exit(read_int(optarg, &num_parts), c);
      if (num_parts < 1)
        error("number of parts must be 1 or greater", c);
      break;

    case 'c':
      if (strlen(optarg) > 1)
        error(msg_str("colour letter has extra letters"), c);
      else {
        const char *v_letters = "AaBbCcKk";
        const char *p = strchr(v_letters, *optarg);
        if (p)
          color_type = p - v_letters;
        else
          error(
              msg_str("colour letter is '%c', must be from AaBbCcKk", *optarg),
              c);
      }
      break;

    case 'l':
      if (strncmp(optarg, "list", strlen(optarg)) == 0)
        list_idx = 0;
      else {
        if (!read_int(optarg, &list_idx))
          error(msg_str("list number is '%s', must be an integer "
                        "or 'list'",
                        optarg),
                c);
        if (list_idx < 1)
          error("list number cannot be less than 1", c);
      }
      break;

    case 'k':
      kite_only = true;
      break;

    case 'q':
      verb = 0;
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (kite_only && color_type >= 0) {
    color_type = -1;
    warning("colouring ignored for single kite model", 'c');
  }

  triangle.resize(6);
  num_fracs = argc - optind;
  if (num_fracs > 3)
    error(msg_str("cannot give more than 3 fractions (%d given)", num_fracs));

  if (num_fracs > 0)
    model_name = argv[optind];
  if (num_fracs > 1)
    model_name += string(" ") + argv[optind + 1];
  if (num_fracs > 2)
    model_name += string(" ") + argv[optind + 2];

  // args2fracs may change num_fracs
  num_fracs = args2fracs(model_name.c_str(), triangle, errmsg);
  if (model_name.size() && !num_fracs)
    error(errmsg);

  if (num_fracs == 1) { // trapezohedron
    if (heights_set & 4)
      warning("C height cannot be set for a trapezohedron", 'C');
    if (color_type > 1 && color_type < 6) // BbCc (Schwarz only)
      warning("BbCc colouring not valid for trapezohedron, ignored", 'c');
  }

  // print model list if no model is specified
  if (model_name == "" && list_idx < 0)
    list_idx = 0;

  if (list_idx < 0) { // base model only
    if (!std::isnan(angle))
      warning("given, but not used for this model", 'a');
    if (num_parts)
      warning("given, but not used for this model", 'N');
  }
}

bool get_B(const Vec3d &A, Vec3d &B, const Vec3d &C, const Vec3d &norm_AB)
{
  norm_AB.dump();
  Vec3d pivot = (C + Trans3d::reflection(norm_AB) * C) / 2;
  if (vcross((A - pivot).unit(), B.unit()).len() >
      epsilon) // lines not parallel
    B = lines_intersection(A, pivot, Vec3d(0, 0, 0), B, 0);
  else // lines parallel
    B *= (A.len() + C.len()) *
         10000; // large relative dist to represent infinity

  for (int i = 0; i < 3; i++)
    if (std::isnan(B[i]))
      return false;
  return true;
}

bool get_C(const Vec3d &A, const Vec3d &B, Vec3d &C)
{
  // AOB not colinear in Schwarz triangle
  Vec3d n = nearest_point(Vec3d(0, 0, 0), A, B);
  C = line_plane_intersect(A, n, Vec3d(0, 0, 0), C);
  for (int i = 0; i < 3; i++)
    if (std::isnan(C[i]))
      return false;
  return true;
}

bool is_degenerate_model(const vector<Vec3d> &kite, char *errmsg)
{

  // kite vertices - A:0, C:1, B:2, C':3
  Vec3d k012_offset = kite[1] - nearest_point(kite[1], kite[0], kite[2]);
  if (k012_offset.len() < epsilon) {
    strcpy_msg(errmsg, "kite shape is a line segment");
    return true;
  }

  Vec3d k123_offset = kite[2] - nearest_point(kite[2], kite[1], kite[3]);
  Vec3d k301_offset = kite[0] - nearest_point(kite[0], kite[3], kite[1]);
  if (k123_offset.len() < epsilon || k301_offset.len() < epsilon) {
    strcpy_msg(errmsg, "kite shape is a triangle");
    return true;
  }

  if (fabs(vtriple(kite[0], kite[1], kite[2])) < epsilon) {
    strcpy_msg(errmsg, "model does not enclose any space");
    return true;
  }

  return false;
}

bool fraction_to_trap_kite(vector<int> tri, Symmetry &sym, vector<Vec3d> &kite,
                           unsigned int hts_set, double hts[],
                           vector<double> *hts_used)
{
  int n = tri[0] / gcd(tri[0], tri[1]); // base polygon number of sides
  // Index numbers used depend on construction of trapezohedron model,
  // none of the vertices will coincide with the origin.
  int A_idx = 2 * n + 1;
  int B_idx = n + 1;
  int C_idx = 1;
  Polygon trap(tri[0], tri[1], Polygon::antiprism,
               Polygon::sub_antiprism_trapezohedron);
  Geometry geom;
  const vector<Vec3d> &verts = geom.verts();
  trap.make_poly(geom);
  // make the default a bit squashed
  geom.transform(
      Trans3d::scale(1, 1, 0.5 * verts[B_idx].len() / verts[A_idx][2]));
  sym = Symmetry(msg_str("D%d", tri[0]));
  if (hts_set & 1 || !(hts_set & 2)) { // All cases except B only
    if (hts_set & 1)                   // A was set
      geom.transform(Trans3d::scale(1, 1, hts[0] / verts[A_idx][2]));
    if (hts_set & 2) { // A and B height were set
      double xy_rad = Vec3d(verts[1][0], verts[1][1], 0).len();
      double root = hts[1] * hts[1] - verts[1][2] * verts[1][2];
      if (root < 0 || xy_rad == 0)
        return false;
      double scale = sqrt(root) / xy_rad;
      geom.transform(Trans3d::scale(scale, scale, 1));
    }
  }
  else { // B only
    geom.transform(Trans3d::scale(hts[1] / verts[B_idx].len()));
  }

  kite.resize(4);
  kite[0] = verts[A_idx];
  kite[1] = verts[C_idx];
  kite[2] = verts[B_idx];
  kite[3] = Trans3d::reflection(vcross(Vec3d::Z, verts[B_idx])) * verts[C_idx];

  hts_used->resize(3);
  (*hts_used)[0] = kite[0][2];
  (*hts_used)[1] = kite[2].len();
  (*hts_used)[2] = kite[1].len();

  return true;
}

bool triangle_to_kite(const vector<Vec3d> &tri_verts, vector<Vec3d> &kite,
                      unsigned int hts_set, double hts[],
                      vector<double> *hts_used)
{
  Vec3d A = tri_verts[0].unit();
  if (hts_set & 1)
    A *= hts[0];
  else
    A *= 0.9;
  Vec3d B = tri_verts[1].unit();
  if (hts_set & 2)
    B *= hts[1];
  else
    B *= 1.3;
  Vec3d C = tri_verts[2].unit();
  if (hts_set & 4)
    C *= hts[2];

  Vec3d norm_AB = vcross(tri_verts[0], tri_verts[1]); // normal to A and B axes

  int ret = false;
  if (!(hts_set & 4)) // Use A and B
    ret = get_C(A, B, C);
  else if (hts_set == 4 || hts_set == 5) // C default A || C and A
    ret = get_B(A, B, C, norm_AB);
  else if (hts_set == 6) {                // C and B
    ret = get_B(B, A, C, norm_AB);        // switch A and B
    C = Trans3d::reflection(norm_AB) * C; // correct the kite orientation
  }

  kite.resize(4);
  kite[0] = A;
  kite[1] = C;
  kite[2] = B;
  kite[3] = Trans3d::reflection(norm_AB) * C;

  hts_used->resize(3);
  (*hts_used)[0] = vdot(kite[0], tri_verts[0]);
  (*hts_used)[1] = vdot(kite[2], tri_verts[1]);
  (*hts_used)[2] = vdot((hts_set == 6) ? kite[3] : kite[1], tri_verts[2]);

  return ret;
}

void repeat_kite_with_color(Geometry &out_geom, Geometry kite_geom,
                            Symmetry sym, const vector<int> &fracs,
                            int color_type)
{
  int vert_no = color_type / 2;
  int v_map[] = {0, 2, 1}; // B is at 3 and C is at 2
  Vec3d axis = kite_geom.verts(v_map[vert_no]);
  Symmetry c_sym;
  c_sym.init(Symmetry::C, fracs[2 * vert_no], Trans3d::rotate(axis, Vec3d::Z));

  Geometry vert_kites_geom;
  sym_repeat(vert_kites_geom, kite_geom, c_sym);
  Coloring clrngs[3];
  if (is_even(color_type)) // value letter is followed by index letter
    clrngs[FACES].add_cmap(colormap_from_name("spread"));

  Transformations min_ts;
  min_ts.min_set(sym.get_trans(), c_sym.get_trans());

  Geometry comp_geom;
  sym_repeat(out_geom, vert_kites_geom, min_ts, ELEM_FACES, clrngs);
}

// Flip the triangle to completes the lune opposite vertex
void flip_triangle(vector<int> &fracs, int frac)
{
  fracs[((2 * frac) + 3) % 6] =
      fracs[((2 * frac) + 2) % 6] - fracs[((2 * frac) + 3) % 6];
  fracs[((2 * frac) + 5) % 6] =
      fracs[((2 * frac) + 4) % 6] - fracs[((2 * frac) + 5) % 6];
}

void get_list(vector<vector<string>> &list, string model, int N)
{
  vector<int> fracs(6);
  int num_fracs = args2fracs(model.c_str(), fracs);
  string trap2;         // trapezohedron pair
  string schwarz_equiv; // equivalend base Schwarz model
  if (num_fracs == 1) { // trapezohedron
    trap2 = msg_str("%d/%d", 2 * fracs[0], 2 * fracs[1]);
    vector<string> desc;
    desc.push_back(model);
    list.push_back(desc);
  }
  else if (num_fracs ==
           3) { // check if fracs correspond to a base Schwarz model
    for (int i = 0; strcmp(base_models[i], "end") != 0; i++) {
      vector<int> b_fracs;
      if (tri_from_str(base_models[i], b_fracs) == 0 &&
          fracs[4] == b_fracs[4]) { // C
        if (fracs[5] != b_fracs[5])
          flip_triangle(fracs, 0);    // covert C to suplement
        if (fracs[5] == b_fracs[5]) { // C compatible
          // check the four A, B combinations, one at a time for simplicity
          if (fracs == b_fracs)
            schwarz_equiv = base_models[i];
          else {
            flip_triangle(fracs, 2);
            if (fracs == b_fracs)
              schwarz_equiv = base_models[i];
            else {
              std::swap(fracs[0], fracs[2]);
              std::swap(fracs[1], fracs[3]);
              if (fracs == b_fracs)
                schwarz_equiv = base_models[i];
              else {
                flip_triangle(fracs, 2);
                if (fracs == b_fracs)
                  schwarz_equiv = base_models[i];
              }
            }
          }
        }
      }
      if (schwarz_equiv != "") { // found a compatible equivalent base model
        vector<string> desc;
        desc.push_back(model);
        list.push_back(desc);
        break;
      }
    }
  }
  else if (model == "") { // all base models match
    for (int i = 0; strcmp(base_models[i], "end") != 0; i++) {
      vector<string> desc;
      desc.push_back(base_models[i]);
      list.push_back(desc);
    }
  }

  for (int i = 0; strcmp(compounds[i][0], "end") != 0; i++) {
    if (model == "" || model == compounds[i][0] ||
        schwarz_equiv == compounds[i][0] || trap2 == compounds[i][0]) {
      vector<string> desc;
      if (schwarz_equiv == compounds[i][0])
        desc.push_back(model);
      else
        desc.push_back(compounds[i][0]);
      desc.push_back(compounds[i][1]);
      desc.push_back(compounds[i][2]);
      desc.push_back(compounds[i][3]);
      list.push_back(desc);
    }
    if (num_fracs == 1 && compounds[i][0][0] == 'n') {
      //   { "n/d_odd" , "DN_odd*nv" , "S2n", ""},    // uc22
      //   { "n/d_odd" , "DN_even*nh", "S2n", ""},    // uc22
      //   { "n/d_odd" , "DN_odd*nv" , "Dnv", ""},    // uc23
      //   { "n/d_odd" , "DN_even*nh", "Dnv", ""},    // uc23
      //   { "n/d_even", "DN*nh"     , "Cnh", ""},    // uc24
      //   { "n/d_even", "DN*nh"     , "Dnh", ""},    // uc25
      int denom_needed = 0;
      if (strcmp(compounds[i][0] + 4, "odd") == 0)
        denom_needed = 1;
      else if (strcmp(compounds[i][0] + 4, "even") == 0)
        denom_needed = 2;

      int N_needed = 0;
      if (strncmp(compounds[i][1] + 3, "odd", 3) == 0)
        N_needed = 1;
      else if (strncmp(compounds[i][1] + 3, "even", 4) == 0)
        N_needed = 2;

      bool denom_ok = is_even(fracs[1]) == is_even(denom_needed);
      bool N_ok = !N || is_even(fracs[1]) || is_even(N) == is_even(N_needed);
      if (denom_ok && N_ok) {
        vector<string> desc;
        desc.push_back(model);

        desc.push_back(compounds[i][1]);
        desc.push_back(compounds[i][2]);
        list.push_back(desc);
      }
    }
  }
}

void print_list(vector<vector<string>> &list)
{
  for (unsigned int i = 0; i < list.size(); i++) {
    printf("%2d: %12s", i + 1, list[i][0].c_str());
    for (unsigned int j = 1; j < list[i].size(); j++) {
      printf(" %12s", list[i][j].c_str());
    }
    printf("\n");
  }
}

void print_base_model_report(const string &model_name,
                             const string &orig_model_name,
                             const vector<double> &hts_used,
                             const string &model_sym, const string &sym_used)
{
  fprintf(stderr, "\n");
  fprintf(stderr, "Base Model\n");
  fprintf(stderr, "   model:               %s", model_name.c_str());
  if (model_name != orig_model_name)
    fprintf(stderr, " (from %s)", orig_model_name.c_str());
  fprintf(stderr, "\n");
  fprintf(stderr, "   A:                   %-16.14g\n", hts_used[0]);
  fprintf(stderr, "   B:                   %-16.14g\n", hts_used[1]);
  fprintf(stderr, "   C:                   %-16.14g\n", hts_used[2]);
  fprintf(stderr, "   symmetry:            %s", model_sym.c_str());
  if (model_sym != sym_used)
    fprintf(stderr, " (expected %s)", sym_used.c_str());
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
}

bool make_base_model(Geometry &geom, const vector<int> &fracs, int num_fracs,
                     unsigned int heights_set, double heights[],
                     char color_type, vector<double> *hts_used,
                     string *sym_used, bool kite_only, char *errmsg)
{
  *errmsg = '\0';
  Geometry kite_geom;
  vector<Vec3d> schwarz_verts;
  Symmetry sym;
  bool is_schwarz = (num_fracs > 1);
  if (is_schwarz) {
    if (!get_schwarz_tri_verts(fracs, schwarz_verts, &sym)) {
      strcpy_msg(errmsg, "fractions do not correspond to a Schwarz "
                         "triangle");
      return false;
    }
    if (!triangle_to_kite(schwarz_verts, kite_geom.raw_verts(), heights_set,
                          heights, hts_used)) {
      strcpy_msg(errmsg, "could not calculate Schwarz model with "
                         "these heights");
      return false;
    }
  }
  else {
    if (!fraction_to_trap_kite(fracs, sym, kite_geom.raw_verts(), heights_set,
                               heights, hts_used)) {
      strcpy_msg(errmsg, "could not calculate trapezohedron with these "
                         "heights (try increasing B height)");
      return false;
    }
  }

  kite_geom.add_face(0, 1, 2, 3, -1);

  *sym_used = sym.get_symbol();
  if (sym.get_symbol()[0] == 'D')
    *sym_used += is_even(fracs[1]) ? 'h' : 'v';
  else if (sym.get_symbol()[0] == 'T')
    *sym_used += 'd';
  else
    *sym_used += 'h';

  char errmsg2[MSG_SZ];
  if (is_degenerate_model(kite_geom.verts(), errmsg2))
    sprintf(errmsg, "degenerate model: %s", errmsg2);

  if (color_type == 0 || color_type == 1 ||             // All Aa
      (is_schwarz && color_type > 1 && color_type < 6)) // Schwarz BbCc
    repeat_kite_with_color(geom, kite_geom, sym, fracs, color_type);
  else
    sym_repeat(geom, kite_geom, sym);

  if (!kite_only) // avoid reordering faces if single kite is to be output
    merge_coincident_elements(geom, "v", epsilon);
  return true;
}

bool make_list_model(Geometry &geom, const string model_str,
                     vector<string> &model_desc, unsigned int heights_set,
                     double heights[], int num_parts, double angle,
                     char color_type, bool kite_only, char *errmsg,
                     vector<string> &warnings, int verb)
{
  warnings.clear();
  string base_model =
      (model_str.size() && model_desc.size() <= 1) ? model_str : model_desc[0];
  string orig_base_model = base_model;
  if (base_model == "n/d" || base_model == "n/d_odd")
    base_model = "7"; // example model
  else if (base_model == "n/d_even")
    base_model = "7/2"; // example model

  vector<int> fracs(6);
  int num_fracs = args2fracs(base_model.c_str(), fracs, errmsg);

  Geometry base_geom;
  vector<double> hts_used;
  string sym_used;
  if (!make_base_model(base_geom, fracs, num_fracs, heights_set, heights,
                       color_type, &hts_used, &sym_used, kite_only, errmsg))
    return false;
  if (*errmsg)
    warnings.push_back(errmsg);

  Coloring clrngs[3];
  if (is_even(color_type)) // value letter is followed by index letter
    clrngs[FACES].add_cmap(colormap_from_name("spread"));

  if (model_desc.size() == 1) {
    geom = base_geom;
    if (color_type == 6 || color_type == 7) { // Kk
      clrngs[FACES].set_geom(&geom);
      clrngs[FACES].f_one_col(clrngs[2].get_col(0));
    }
    Symmetry base_sym(base_geom);
    if (verb) {
      print_base_model_report(base_model, orig_base_model, hts_used,
                              base_sym.get_symbol(), sym_used);
    }
    return true;
  }

  string full_sym_str = model_desc[1];
  string part_sym_str = model_desc[2];

  int realign = 0;
  int sz = part_sym_str.size();
  if (sz >= 2 && part_sym_str.substr(sz - 2) == ",1") {
    realign = 1;
    part_sym_str.resize(sz - 2);
  }

  bool N_unneeded = full_sym_str.substr(0, 2) != "DN";
  if (num_parts && N_unneeded)
    warnings.push_back("option -N: given, but not used for this model");

  if (num_fracs == 1) { // trapezium
    if (num_parts == 0)
      num_parts = (full_sym_str == "DN_odd*nv") ? 3 : 2;
    if (full_sym_str == "DN_even*nh" || full_sym_str == "DN_odd*nv" ||
        full_sym_str == "DN*nh")
      full_sym_str =
          msg_str("D%d%c", num_parts * fracs[0], *full_sym_str.rbegin());

    if (part_sym_str == "S2n")
      part_sym_str = msg_str("S%d", 2 * fracs[0]);
    if (part_sym_str[1] == 'n')
      part_sym_str =
          msg_str("%c%d%c", part_sym_str[0], fracs[0], *part_sym_str.rbegin());
    if (heights_set & 4)
      warnings.push_back(
          "option -C: C height cannot be set for a trapezohedron\n");
  }

  Symmetry full_sym(full_sym_str);
  Symmetry part_sym(part_sym_str);

  Symmetry final_sym;
  full_sym.get_sub_sym(part_sym, &final_sym);
  final_sym.get_autos().set_fixed_type(realign);

  double ang = !std::isnan(angle) ? angle : 3; // example (small) angle
  bool angle_unneeded = !final_sym.get_autos().set_rot_principal(ang);
  if (angle_unneeded && !std::isnan(angle)) {
    warnings.push_back("option -a: given, but not used for this model");
    angle_unneeded = true;
  }

  Trans3d trans =
      (final_sym.get_autos().get_realignment() * final_sym.get_to_std())
          .inverse(); // to_std not needed, but leave for reference
  base_geom.transform(trans);

  Symmetry base_sym(base_geom);
  if (model_desc.size() > 3 && model_desc[3] != "") {
    int conj_type = atoi(model_desc[3].c_str());
    Symmetry sym;
    base_sym.get_sub_sym(part_sym, &sym, conj_type);
    Trans3d trans = sym.get_to_std().inverse() * base_sym.get_to_std();
    base_geom.transform(trans);
  }

  Transformations min_ts;
  min_ts.min_set(full_sym.get_trans(), base_sym.get_trans());

  if (color_type == 6 || color_type == 6) { // Kk
    sym_repeat(geom, base_geom, min_ts, ELEM_FACES, clrngs);
  }
  else
    sym_repeat(geom, base_geom, min_ts);

  if (verb) {
    print_base_model_report(base_model, orig_base_model, hts_used,
                            base_sym.get_symbol(), sym_used);
    fprintf(stderr, "\n");
    fprintf(stderr, "Compound\n");
    fprintf(stderr, "   full symmetry:       %s", full_sym_str.c_str());
    if (full_sym_str != model_desc[1])
      fprintf(stderr, " (from %s)", model_desc[1].c_str());
    fprintf(stderr, "\n");
    fprintf(stderr, "   alignment symmetry:  %s", part_sym_str.c_str());
    if (part_sym_str != model_desc[2])
      fprintf(stderr, " (from %s)", model_desc[2].c_str());
    fprintf(stderr, "\n");
    fprintf(stderr, "   N:                   %s\n",
            (N_unneeded) ? "n/a" : msg_str("%d", num_parts).c_str());
    fprintf(stderr, "   angle:               %s\n",
            (angle_unneeded) ? "n/a" : msg_str("%-16.14g", ang).c_str());
    fprintf(stderr, "   components:          %d", (int)min_ts.size());
    int expected_num_comps =
        full_sym.get_trans().size() / part_sym.get_trans().size();
    if ((int)min_ts.size() != expected_num_comps)
      fprintf(stderr, " (expected %d)", expected_num_comps);
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
  }

  return true;
}

int main(int argc, char *argv[])
{
  kt_opts opts;
  opts.process_command_line(argc, argv);

  char errmsg[MSG_SZ];
  Geometry out_geom;
  if (opts.list_idx < 0) {
    vector<double> hts_used;
    string sym_used;
    if (!make_base_model(out_geom, opts.triangle, opts.num_fracs,
                         opts.heights_set, opts.heights, opts.color_type,
                         &hts_used, &sym_used, opts.kite_only, errmsg))
      opts.error(errmsg);
    if (*errmsg)
      opts.warning(errmsg);
    if (opts.verb)
      print_base_model_report(opts.model_name, opts.model_name, hts_used,
                              Symmetry(out_geom).get_symbol(), sym_used);
  }
  else {
    vector<vector<string>> list;
    get_list(list, opts.model_name, opts.num_parts);
    if (!list.size())
      opts.error(
          "list cannot be used with degenerate Schwarz models \n"
          "(see supported models by running the program without options)\n",
          'l');
    if (opts.list_idx == 0) {
      print_list(list);
      return 0;
    }
    else if (opts.list_idx <= (int)list.size()) {
      vector<string> warnings;
      if (!make_list_model(out_geom, opts.model_name, list[opts.list_idx - 1],
                           opts.heights_set, opts.heights, opts.num_parts,
                           opts.angle, opts.color_type, opts.kite_only, errmsg,
                           warnings, opts.verb))
        opts.error(errmsg);
      for (auto &warning : warnings)
        opts.warning(warning.c_str());
    }
    else
      opts.error(msg_str("list number too large (maximum %d)", list.size()),
                 'l');
  }

  if (opts.kite_only) {
    Geometry kite;
    for (int i = 0; i < 4; i++)
      kite.add_vert(out_geom.face_v(0, i));
    kite.add_face(0, 1, 2, 3, -1);
    out_geom = kite;
  }

  opts.write_or_error(out_geom, opts.ofile);

  return 0;
}
