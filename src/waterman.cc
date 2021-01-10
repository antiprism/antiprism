/*
   Copyright (c) 2008-2021, Roger Kaufman, Adrian Rossiter

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
   Name: waterman.cc
   Description: sphere-ray intersection for producing waterman polyhedra
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"
#include "lattice_grid.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib> // avoid ambiguities with std::abs(long) on OSX
#include <string>
#include <vector>

using std::make_pair;
using std::pair;
using std::string;
using std::vector;

using namespace anti;

int get_num_decs(const char *str)
{
  const char *p = strchr(str, '.');
  if (!p)
    return 0;
  // Find the number of digits after the decimal point
  size_t num_decs = strspn(p + 1, "1234567890");
  // Disregard trailing 0's
  while (*(p + num_decs) == '0')
    num_decs--;
  return num_decs;
}

int get_max_num_decs(const char *str, const int num_parts)
{
  Split parts(str, ",");
  int max_num_decs = 0;
  for (int i = 0; i < (int)parts.size() && i < num_parts; i++) {
    int num_decs = get_num_decs(parts[i]);
    if (num_decs > max_num_decs)
      max_num_decs = num_decs;
  }
  return max_num_decs;
}

class waterman_opts : public ProgramOpts {
public:
  string ofile;

  int lattice_type;
  double radius;
  double R_squared;
  bool origin_based;
  Vec3d center;
  int method;
  bool fill;
  bool verbose;
  long scale;
  bool tester_defeat;
  bool convex_hull;
  bool add_hull;
  Color vert_col;
  Color edge_col;
  Color face_col;
  Color fill_col;
  char color_method;
  int face_opacity;
  double eps;

  waterman_opts()
      : ProgramOpts("waterman"), lattice_type(-1), radius(0), R_squared(0),
        origin_based(true), method(1), fill(false), verbose(false), scale(0),
        tester_defeat(false), convex_hull(true), add_hull(false),
        color_method('\0'), face_opacity(-1), eps(anti::epsilon)
  {
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

void waterman_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] lattice

Use sphere-ray intersection for producing Waterman Polyhedra. Lattice can be
SC, FCC, or BCC

Options
%s
  -r <r,n>  clip radius. r is radius taken to optional root n. n = 2 is sqrt
  -q <cent> center of lattice, in form "x_val,y_val,z_val" (default: origin)
  -m <mthd> 1 - sphere-ray intersection  2 - z guess (default: 1)
  -C <opt>  c - convex hull only, i - keep interior, s - suppress (default: c)
  -f        fill interior points (not for -C c)
  -t        defeat computational error testing for sphere-ray method
  -v        verbose output (on computational errors)
  -l <lim>  minimum distance for unique vertex locations as negative exponent
               (default: %d giving %.0e)
  -o <file> write output to file (default: write to standard output)

Coloring Options (run 'off_util -H color' for help on color formats)
  -V <col>  model vertex color (default: none)
  -E <col>  edge color (for convex hull, default: none)
  -F <col>  face color (for convex hull, default: none)
               lower case outputs map indexes. upper case outputs color values
               keyword: s,S color by symmetry using face normals
               keyword: c,C color by symmetry using face normals (chiral)
  -Z <col>  fill vertex color (default: model vertex color)
  -T <tran> face transparency. valid range from 0 (invisible) to 255 (opaque)

)",
          prog_name(), help_ver_text, int(-log(anti::epsilon) / log(10) + 0.5),
          anti::epsilon);
}

void waterman_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  vector<double> double_parms;
  double root = 1;
  int cent_num_decs = 0;
  int R_num_decs = 0;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hr:q:m:ftvC:V:E:F:Z:T:l:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'r': {
      R_num_decs = get_max_num_decs(optarg, 1);
      bool sqrt_found = (!strncmp(optarg, "sqrt", 4)) ? true : false;
      print_status_or_exit(read_double_list(optarg, double_parms, 2), c);
      if (double_parms.size() == 2 || sqrt_found) {
        root = (sqrt_found) ? 2 : double_parms[1];
        if (!root)
          error("root for radius must be non-zero", c);
        if (sqrt_found)
          radius = double_parms[0];
        else
          radius = pow(double_parms[0], 1 / root);
        R_squared = radius * radius;
        if (root == 2) { // more precise radius squared for sphere-ray method
          if (sqrt_found)
            R_squared = atof(optarg + 4);
          else
            R_squared = double_parms[0];
          R_num_decs = (R_num_decs + 1) / 2; // change in case of root 2
        }
        else if (root != 1)
          R_num_decs = std::numeric_limits<int>::max(); // not suitable for
                                                        // integer calculations
      }
      else { // radius given as explicit number
        radius = double_parms[0];
        R_squared = radius * radius;
      }
      if (radius <= 0)
        error("sphere radius cannot be negative", c);
      break;
    }

    case 'q':
      cent_num_decs = get_max_num_decs(optarg, 3);
      print_status_or_exit(center.read(optarg), c);
      if (compare(center, Vec3d(0, 0, 0), anti::epsilon))
        origin_based = false;
      break;

    case 'm':
      print_status_or_exit(read_int(optarg, &method), c);
      if (method < 1 || method > 2) {
        error("method must be 1 or 2", c);
      }
      break;

    case 'f':
      fill = true;
      break;

    case 't':
      tester_defeat = true;
      break;

    case 'v':
      verbose = true;
      break;

    case 'C':
      if (strlen(optarg) > 1 || !strchr("cis", *optarg))
        error("convex hull arg is '" + string(optarg) + "' must be c, i or s",
              c);
      if (strchr("s", *optarg))
        convex_hull = false;
      if (strchr("i", *optarg))
        add_hull = true;
      break;

    case 'V':
      print_status_or_exit(vert_col.read(optarg), c);
      break;

    case 'E':
      print_status_or_exit(edge_col.read(optarg), c);
      break;

    case 'F':
      if (strchr("SsCc", *optarg))
        color_method = *optarg;
      else
        print_status_or_exit(face_col.read(optarg), c);
      break;

    case 'Z':
      print_status_or_exit(fill_col.read(optarg), c);
      break;

    case 'T':
      print_status_or_exit(read_int(optarg, &face_opacity), c);
      if (face_opacity < 0 || face_opacity > 255) {
        error("face transparency must be between 0 and 255", c);
      }
      break;

    case 'l':
      int sig_compare;
      print_status_or_exit(read_int(optarg, &sig_compare), c);
      if (sig_compare > DEF_SIG_DGTS)
        warning("limit is very small, may not be attainable", c);
      eps = pow(10, -sig_compare);
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (argc == optind)
    error("no lattice specified", "lattice");

  if (argc - optind > 1)
    error("too many arguments");

  string lattice_type_str;
  if (optind < argc)
    lattice_type_str = argv[optind++];

  transform(lattice_type_str.begin(), lattice_type_str.end(),
            lattice_type_str.begin(), ::tolower);

  if (lattice_type_str == "sc")
    lattice_type = 0;
  else if (lattice_type_str == "fcc")
    lattice_type = 1;
  else if (lattice_type_str == "bcc")
    lattice_type = 2;

  if (lattice_type == -1)
    error("lattice must be sc, fcc, or bcc");

  if (!center.is_set())
    center = Vec3d(0, 0, 0);

  if (face_opacity > -1) {
    // if no face color is set but there is transparency set, use white
    if (!face_col.is_set())
      face_col = Color(255, 255, 255);

    // face color can only be made transparent if not index and not invisible
    if (!face_col.set_alpha(face_opacity))
      warning("transparency has no effect on map indexes or invisible", "T");

    if (color_method == 's' || color_method == 'c')
      warning("when writing indexes transparency setting ignored", "T");
  }

  if (fill && convex_hull && !add_hull) {
    warning("fill cannot be used with -C c", 'f');
    fill = false;
  }

  // fill color vertex default is that of outer vertex color
  if (vert_col.is_set() && !fill_col.is_set())
    fill_col = vert_col;

  if (tester_defeat) {
    if (method == 1)
      warning("computational error testing has been disabled!");
    else
      warning("for z-guess method -t has no effect");
  }

  // Choose scale to clear decimals
  int max_num_decs = (R_num_decs > cent_num_decs) ? R_num_decs : cent_num_decs;
  if (max_num_decs < 50) {
    double test_scale = pow(10, max_num_decs);
    // if((center.mag() + radius + 1)*test_scale > LONG_MAX/sqrt(3))
    if ((center.len() + radius + 1) * test_scale >
        sqrt(double(std::numeric_limits<long long>::max())) / sqrt(3))
      scale = 0; // don't use integer calculations
    else
      scale = (long)(test_scale + 0.5);
  }
  else
    scale = 0;

  if (!scale) {
    warning("scale is set to zero! computational errors may occur");
    if (method == 2)
      error("z-guess method cannot be used in this case");
  }
}

// returns:
// false ray misses sphere
// true  ray strikes sphere
//    ray tangent to sphere, z_near == z_far
//    ray passes through two points of sphere, z_near != z_far
// z_near and z_far are modified
bool sphere_ray_z_intersect_points(
    long &z_near, long &z_far, const double x0,
    const double y0,         // ray's xy (dz = 0)
    const bool origin_based, // if true, then cx, cy, cz all equal 0
    const double cx, const double cy, const double cz, // center
    const double R_squared,                            // radius squared
    const double eps)
{
  // original formulas taken from
  // http://www.ccs.neu.edu/home/fell/CSU540/programs/RayTracingFormulas.htm
  // double dx = x1 - x0;
  // double dy = y1 - y0;
  // double dz = z1 - z0;

  // double a = dx*dx + dy*dy + dz*dz;
  // double b = 2*dx*(x0-cx) + 2*dy*(y0-cy) + 2*dz*(z0-cz);
  // double c = cx*cx + cy*cy + cz*cz + x0*x0 + y0*y0 + z0*z0 + -2*(cx*x0 +
  // cy*y0 + cz*z0) - R*R;
  // end of paste

  // optimization: dx = dy = 0 dz = 1 z0 = 0 (reverses near and far point) by
  // Adrian Rossiter
  double b;
  double c;
  double discriminant;
  if (origin_based) {
    // optimization: if true, then cx, cy, cz all equal 0
    b = 0;
    c = x0 * x0 + y0 * y0 - R_squared;
    discriminant = -4 * c;
  }
  else {
    b = 2 * -cz;
    c = cx * cx + cy * cy + cz * cz + x0 * x0 + y0 * y0 +
        -2 * (cx * x0 + cy * y0) - R_squared;
    discriminant = b * b - 4 * c;
  }

  // explicitly set discriminant to 0 if in -eps to +eps
  // it is set to zero so that sqrt(disciminant) will not result in -nan
  if (double_eq(discriminant, 0.0, eps))
    discriminant = 0;
  if (discriminant < 0.0)
    return false;
  double disc_rt = sqrt(discriminant);

  vector<double> z_intersect(2);

  // near point
  z_intersect[0] = (-b + disc_rt) / 2.0;
  z_near = (long)floor(z_intersect[0] + eps);
  //   if (double_eq(discriminant, 0.0, eps)) {
  if (z_near == (long)floor(cz + eps)) {
    z_far = z_near; // prevent z_far from being stale
  }
  else {
    // far point
    z_intersect[1] = (-b - disc_rt) / 2.0;
    z_far = (long)ceil(z_intersect[1] - eps);
  }
  return true;
}

// combine lattice test for fcc(type = 1) and bcc(type = 2)
bool valid_point(const int lattice_type, const long x, const long y,
                 const long z)
{
  return lattice_type == 1
             ? (x + y + z) % 2 == 0
             : (x % 2 && y % 2 && z % 2) || !(x % 2 || y % 2 || z % 2);
}

bool inside_exact(const long z, const long z_cent,
                  const long long xy_contribution, const long long i_R2)
{
  // fprintf(stderr,"z = %ld, z_cent = %ld   xy_c = %I64d, (z-z_cent)^2 = %I64d,
  // i_R2 = %I64d\n",z,z_cent,xy_contribution,(long
  // long)(z-z_cent)*(z-z_cent),i_R2);
  return (xy_contribution + (long long)(z - z_cent) * (z - z_cent)) <= i_R2;
}

// tester code furnished by Adrian Rossiter
// z_near and z_far are modified
void refine_z_vals(long &z_near, long &z_far, const long x, const long y,
                   const int lattice_type, const long scale,
                   const vector<long> i_center, const long long i_R2)
{
  long long xy_contribution =
      ((long long)x * scale - i_center[0]) * (x * scale - i_center[0]) +
      ((long long)y * scale - i_center[1]) * (y * scale - i_center[1]);
  // fprintf(stderr, "x=%ld, y=%ld, z_near=%ld, z_far=%ld, i_center=(%ld, %ld,
  // %ld), i_R2=%lld\n", x, y, z_near, z_far, i_center[0], i_center[1],
  // i_center[2], i_R2);
  // fprintf(stderr, "xy_cont=%lld, scale=%ld\n", xy_contribution, scale);

  long z_nr;
  // Make sure z_start is a valid lattice point on or above z_centre
  long z_start = (z_near * scale > i_center[2])
                     ? z_near
                     : (long)ceil((double)i_center[2] / scale);
  while (lattice_type && !valid_point(lattice_type, x, y, z_start))
    z_start++;

  if (inside_exact(z_start * scale, i_center[2], xy_contribution, i_R2)) {
    z_nr = z_start;
    // Inside, so test for valid points above
    for (long z_val = z_start + 1;; ++z_val) {
      // fprintf(stderr, "zn in z_val=%ld\n", z_val);
      if (lattice_type && !valid_point(lattice_type, x, y, z_val))
        continue;
      if (inside_exact(z_val * scale, i_center[2], xy_contribution, i_R2))
        z_nr = z_val;
      else
        break;
    }
  }
  else {
    // Outside, so test for valid points below
    z_nr = std::numeric_limits<long>::max(); // an invalid value
    for (long z_val = z_start - 1;; --z_val) {
      // fprintf(stderr, "zn out z_val=%ld\n", z_val);
      if (z_val * scale < i_center[2]) {
        break;
      }
      if (lattice_type && !valid_point(lattice_type, x, y, z_val))
        continue;
      if (inside_exact(z_val * scale, i_center[2], xy_contribution, i_R2)) {
        z_nr = z_val;
        break;
      }
    }
  }
  z_near = z_nr;

  long z_fr;
  // Make sure z_start is a valid lattice point on or below z_centre
  z_start = (z_far * scale < i_center[2])
                ? z_far
                : (long)floor((double)i_center[2] / scale);
  while (lattice_type && !valid_point(lattice_type, x, y, z_start))
    z_start--;

  if (inside_exact(z_start * scale, i_center[2], xy_contribution, i_R2)) {
    z_fr = z_start;
    // Inside, so test for valid points below
    for (long z_val = z_start - 1;; --z_val) {
      // fprintf(stderr, "zf in z_val=%ld\n", z_val);
      if (lattice_type && !valid_point(lattice_type, x, y, z_val))
        continue;
      if (inside_exact(z_val * scale, i_center[2], xy_contribution, i_R2))
        z_fr = z_val;
      else
        break;
    }
  }
  else {
    // Outside, so test for valid points above
    z_fr = std::numeric_limits<long>::max(); // an invalid value
    for (long z_val = z_start + 1;; ++z_val) {
      // fprintf(stderr, "zf out z_val=%ld\n", z_val);
      if (z_val * scale > i_center[2]) {
        break;
      }
      if (lattice_type && !valid_point(lattice_type, x, y, z_val))
        continue;
      if (inside_exact(z_val * scale, i_center[2], xy_contribution, i_R2)) {
        z_fr = z_val;
        break;
      }
    }
  }
  z_far = z_fr;
}

// Separate function to contain protability problems with abs(long)
long long_abs(long val) { return std::abs((long)val); }

void sphere_ray_waterman(Geometry &geom, const int lattice_type,
                         const bool origin_based, const Vec3d &center,
                         const double radius, const double R_squared,
                         const long scale, const bool verbose,
                         const bool tester_defeat, const double eps)
{
  vector<Vec3d> &verts = geom.raw_verts();

  // check if z of center is on integer value
  bool cent_z_int = true;
  if (!origin_based) {
    double int_part;
    double fract_part = modf(center[2], &int_part);
    cent_z_int = double_eq(fract_part, 0.0, eps);
  }

  long rad_left_x = (long)ceil(center[0] - radius);
  long rad_right_x = (long)floor(center[0] + radius);
  long rad_bottom_y = (long)ceil(center[1] - radius);
  long rad_top_y = (long)floor(center[1] + radius);

  vector<long> i_center(3);
  for (int i = 0; i < 3; i++)
    i_center[i] = (long)floor(center[i] * scale + 0.5);

  long long i_R2 = (long long)floor(radius * radius * scale * scale + 0.5);

  long z_near = 0;
  long z_far = 0;

  long total_errors = 0;
  long total_misses = 0;

  for (long y = rad_bottom_y; y <= rad_top_y; y++) {
    for (long x = rad_left_x; x <= rad_right_x; x++) {
      // faster miss determination, but using for false miss detection
      bool miss = true;
      long long xy_contribution =
          ((long long)x * scale - i_center[0]) * (x * scale - i_center[0]) +
          ((long long)y * scale - i_center[1]) * (y * scale - i_center[1]);
      if (inside_exact(i_center[2], i_center[2], xy_contribution, i_R2))
        miss = false;
      // continue;

      if (!sphere_ray_z_intersect_points(z_near, z_far, x, y, origin_based,
                                         center[0], center[1], center[2],
                                         R_squared, eps)) {
        // fprintf(stderr,"Ray missed the Sphere\n");
        if (!miss) {
          // if (verbose)
          //   fprintf(stderr,"error: at x = %ld, y = %ld, a false miss
          //   happened\n",x,y);
          total_misses++;
        }
        continue;
      }

      // ray tangent points are never on integer when z of center is not on
      // integer value
      // NEEDS MORE TESTING
      if (!cent_z_int && z_near == z_far)
        continue;

      if (lattice_type != 0) { // lattice type is not equal to SC (type = 0)
        // if z_near is not on the lattice then find if a point 1 layer deeper
        // is on the lattice
        if (!valid_point(lattice_type, long_abs(x), long_abs(y),
                         long_abs(z_near))) {
          // if it is a tangent point, there is no valid deeper coordinate. It
          // was on "zero" already.
          // if bcc and z_near-1 is invalid then there is no valid z point
          // (z_far+1 will be invalid as well)
          if (z_near == z_far ||
              (lattice_type == 2 &&
               !valid_point(lattice_type, long_abs(x), long_abs(y),
                            long_abs(z_near - 1))))
            continue;
          else
            z_near--;
        }
        // if still in the loop, z_far is only advanced if on invalid point
        if (!valid_point(lattice_type, long_abs(x), long_abs(y),
                         long_abs(z_far)))
          z_far++;
      }

      // uncommenting next 2 lines forces errors
      // z_near += 5;
      // z_far += 5;
      if (!tester_defeat && scale) {
        long z_near2 = z_near;
        long z_far2 = z_far;
        refine_z_vals(z_near2, z_far2, x, y, lattice_type, scale, i_center,
                      i_R2);

        if (z_near2 != z_near) {
          total_errors++;
          // if (verbose)
          //   fprintf(stderr, "(%ld, %ld) z_near %ld -> %s\n", x, y, z_near,
          //          (z_near2!=LONG_MAX) ? itostr(z_near2).c_str() :
          //          "invalid");
          z_near = z_near2;
        }

        if (z_far2 != z_far) {
          total_errors++;
          // if (verbose)
          //   fprintf(stderr, "(%ld, %ld) z_far %ld -> %s\n", x, y, z_far,
          //          (z_far2!=LONG_MAX) ? itostr(z_far2).c_str() : "invalid");
          z_far = z_far2;
        }
      }

      // don't write invalid points
      if (z_near != std::numeric_limits<long>::max())
        verts.push_back(Vec3d(x, y, z_near));
      if (z_far != std::numeric_limits<long>::max() &&
          z_near != z_far) // don't rewrite tangent point
        verts.push_back(Vec3d(x, y, z_far));
    }
  }

  if (verbose && !tester_defeat)
    fprintf(stderr, "Total computational errors found and corrected: %ld\n",
            total_errors);
  if (verbose && total_misses)
    fprintf(stderr, "Total number of false misses: %ld\n", total_misses);
}

void z_guess_waterman(Geometry &geom, const int lattice_type,
                      const Vec3d &center, const double radius,
                      const long scale, const bool verbose)
{
  vector<Vec3d> &verts = geom.raw_verts();

  long rad_left_x = (long)ceil(center[0] - radius);
  long rad_right_x = (long)floor(center[0] + radius);
  long rad_bottom_y = (long)ceil(center[1] - radius);
  long rad_top_y = (long)floor(center[1] + radius);

  vector<long> i_center(3);
  for (int i = 0; i < 3; i++)
    i_center[i] = (long)floor(center[i] * scale + 0.5);

  long long i_R2 = (long long)floor(radius * radius * scale * scale + 0.5);

  long z_near = 0;
  long z_far = 0;

  long total_errors = 0;
  // long total_amount = 0;

  for (long y = rad_bottom_y; y <= rad_top_y; y++) {
    for (long x = rad_left_x; x <= rad_right_x; x++) {
      // see if some z point on this x,y is inside the radius
      long long xy_contribution =
          ((long long)x * scale - i_center[0]) * (x * scale - i_center[0]) +
          ((long long)y * scale - i_center[1]) * (y * scale - i_center[1]);
      if (!inside_exact(i_center[2], i_center[2], xy_contribution, i_R2)) {
        // reset z_near and z_far for next guess
        z_near = 0;
        z_far = 0;
        continue; // miss
      }
      else {
        // if we are on a bcc "tunnel" skip this x,y
        if (lattice_type == 2 &&
            !valid_point(lattice_type, long_abs(x), long_abs(y), 0) &&
            !valid_point(lattice_type, long_abs(x), long_abs(y), 1))
          continue;

        long z_near2 = z_near;
        long z_far2 = z_far;
        refine_z_vals(z_near2, z_far2, x, y, lattice_type, scale, i_center,
                      i_R2);

        if (z_near2 != z_near) {
          total_errors++;
          // total_amount+=long_abs(z_near2-z_near);
          // if (verbose)
          //   fprintf(stderr, "(%ld, %ld) z_near %ld -> %s\n", x, y, z_near,
          //          (z_near2!=LONG_MAX) ? itostr(z_near2).c_str() :
          //          "invalid");
          z_near = z_near2;
        }

        if (z_far2 != z_far) {
          total_errors++;
          // total_amount+=long_abs(z_far2-z_far);
          // if (verbose)
          //   fprintf(stderr, "(%ld, %ld) z_far %ld -> %s\n", x, y, z_far,
          //          (z_far2!=LONG_MAX) ? itostr(z_far2).c_str() : "invalid");
          z_far = z_far2;
        }

        // don't write invalid points
        if (z_near != std::numeric_limits<long>::max())
          verts.push_back(Vec3d(x, y, z_near));
        else
          z_near = 0; // when invalid, reset z_far for next guess

        if (z_far != std::numeric_limits<long>::max() &&
            z_near != z_far) // don't rewrite tangent point
          verts.push_back(Vec3d(x, y, z_far));
        else
          z_far = 0; // when invalid, reset z_far for next guess
      }
    }
  }

  if (verbose)
    fprintf(stderr, "Total computational errors found and corrected: %ld\n",
            total_errors);
  // fprintf(stderr,"Total errors amount: %ld\n",total_amount);
}

// fill interior points
Geometry fill_interior(const Geometry &geom, const int lattice_type)
{
  const vector<Vec3d> &verts = geom.verts();

  map<pair<long, long>, long> min_z_vert;
  map<pair<long, long>, long> max_z_vert;

  map<pair<long, long>, long>::iterator it;
  for (long i = 0; i < (long)verts.size(); i++) {
    long x = lround(verts[i][0]);
    long y = lround(verts[i][1]);
    long z = lround(verts[i][2]);

    pair<long, long> key = make_pair(x, y);

    // initialize values
    it = min_z_vert.find(key);
    if (it == min_z_vert.end())
      min_z_vert[key] = std::numeric_limits<long>::max();
    it = max_z_vert.find(key);
    if (it == max_z_vert.end())
      max_z_vert[key] = std::numeric_limits<long>::min();

    // find minimum z and maximum z for an x,y
    if (z < min_z_vert[key])
      min_z_vert[key] = z;
    if (z > max_z_vert[key])
      max_z_vert[key] = z;
  }

  Geometry fill;
  vector<Vec3d> &fill_verts = fill.raw_verts();

  for (auto const &key : min_z_vert) {
    long x = key.first.first;
    long y = key.first.second;
    long z_min = key.second;
    long z_max = max_z_vert[make_pair(x, y)];
    for (long i = z_min + 1; i < z_max; i++) {
      if (!lattice_type || valid_point(lattice_type, x, y, i))
        fill_verts.push_back(Vec3d(x, y, i));
    }
  }

  return fill;
}

int main(int argc, char *argv[])
{
  waterman_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;

  if (opts.verbose)
    fprintf(stderr, "calculating points\n");

  if (opts.method == 1)
    sphere_ray_waterman(geom, opts.lattice_type, opts.origin_based, opts.center,
                        opts.radius, opts.R_squared, opts.scale, opts.verbose,
                        opts.tester_defeat, opts.eps);
  else
    z_guess_waterman(geom, opts.lattice_type, opts.center, opts.radius,
                     opts.scale, opts.verbose);

  // interior filling
  Geometry fill_verts;
  if (opts.fill) {
    if (opts.verbose)
      fprintf(stderr, "filling interior\n");

    fill_verts = fill_interior(geom, opts.lattice_type);
    if (opts.fill_col.is_set())
      Coloring(&fill_verts).v_one_col(opts.fill_col);
  }

  // convex hull and coloring
  if (opts.convex_hull) {
    if (opts.verbose)
      fprintf(stderr, "performing convex hull\n");

    Status stat = (opts.add_hull) ? geom.add_hull() : geom.set_hull();
    if (stat.is_error()) {
      if (opts.verbose)
        fprintf(stderr, "%s\n", stat.c_msg());
    }
    else {
      geom.orient();
      if (opts.verbose)
        convex_hull_report(geom, opts.add_hull);

      if (opts.color_method)
        color_by_symmetry_normals(geom, opts.color_method, opts.face_opacity,
                                  opts.eps);
      else
        Coloring(&geom).f_one_col(opts.face_col);

      if (opts.edge_col.is_set()) {
        geom.add_missing_impl_edges(opts.edge_col);
      }
    }
  }

  // color vertices
  if (opts.vert_col.is_set())
    Coloring(&geom).v_one_col(opts.vert_col);

  // append fill points
  if (opts.fill)
    geom.append(fill_verts);

  if (opts.verbose)
    fprintf(stderr, "writing output\n");

  opts.write_or_error(geom, opts.ofile);

  if (opts.verbose)
    fprintf(stderr, "done!\n");

  return 0;
}
