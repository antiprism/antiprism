/*
   Copyright (c) 2007-2021, Roger Kaufman

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
   Name: off2txt.cc
   Description: Convert files in OFF format to Hedron format
   Project: Antiprism - http://www.antiprism.com
*/

#include <cstdio>
#include <cstdlib>

#include "../base/antiprism.h"

#include <cctype>
#include <string>
#include <vector>

using std::string;
using std::vector;

using namespace anti;

class o2t_opts : public ProgramOpts {
public:
  string ifile;
  string ofile;

  bool estimate_colors = false;      // determine estimated hedron color
  bool detect_rhombi = false;        // set detect rhombi color
  bool detect_star_polygons = false; // set detect star polygons
  bool exclude_coordinates = false;  // don't add coordinates to output
  bool force_transparent = false;    // set transparent symbol

  int sig_digits = DEF_SIG_DGTS; // significant digits output (system default)

  double eps = anti::epsilon;

  o2t_opts() : ProgramOpts("off2txt") {}

  void process_command_line(int argc, char **argv);
  void usage();
};

void o2t_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] [input_file]

Convert files in OFF format Hedron text file format. If
input_file is not given the program reads from standard input.

Options
%s
  -l <lim>  minimum distance for unique vertex locations as negative exponent
               (default: %d giving %.0e)
  -d <dgts> number of significant digits (default %d) or if negative
            then the number of digits after the decimal point
  -o <file> file name for output (otherwise prints to stdout)

Program Options
  -c        estimates colors from OFF file
  -r        detect rhombi. D parameter added if found
  -p        detect star polygons
  -x        exclude coordinates
  -t        force all faces transparent

)",
          prog_name(), help_ver_text, int(-log(anti::epsilon) / log(10) + 0.5),
          anti::epsilon, DEF_SIG_DGTS);
}

void o2t_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hcrpxtl:d:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'c':
      estimate_colors = true;
      break;

    case 'r':
      detect_rhombi = true;
      break;

    case 'p':
      detect_star_polygons = true;
      break;

    case 'x':
      exclude_coordinates = true;
      break;

    case 't':
      force_transparent = true;
      break;

    case 'l':
      int sig_compare;
      print_status_or_exit(read_int(optarg, &sig_compare), c);
      if (sig_compare > DEF_SIG_DGTS)
        warning("limit is very small, may not be attainable", c);
      eps = pow(10, -sig_compare);
      break;

    case 'd':
      print_status_or_exit(read_int(optarg, &sig_digits), c);
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (argc - optind > 1)
    error("too many arguments");

  if (argc - optind == 1)
    ifile = argv[optind];
}

string estimated_color(Color col)
{
  string color;
  int v[3];

  for (unsigned int i = 0; i < 3; i++)
    v[i] = (col[i] <= 128) ? 0 : 255;

  // Special case for orange
  if (v[0] == 255 && (col[1] >= 64 && col[1] <= 192) && v[2] == 0)
    color = "a";
  else if (v[0] == 0 && v[1] == 0 && v[2] == 255)
    color = "b";
  else
      // Special case for green. "0,128,0" is exact green would become 0 0 0
      // Make more dark greens be g instead of using 0 255 0.
      if (v[0] == 0 && col[1] >= 64 && v[2] == 0)
    color = "g";
  else
      // Hedron has no black so specify cyan
      if (v[0] == 0 && v[1] == 0 && v[2] == 0)
    color = "c";
  else if (v[0] == 0 && v[1] == 255 && v[2] == 255)
    color = "c";
  else if (v[0] == 255 && v[1] == 0 && v[2] == 0)
    color = "r";
  else if (v[0] == 255 && v[1] == 0 && v[2] == 255)
    color = "m";
  else if (v[0] == 255 && v[1] == 255 && v[2] == 0)
    color = "y";
  else if (v[0] == 255 && v[1] == 255 && v[2] == 255)
    color = "w";
  else
    color = "\0";

  return (color);
}

bool is_square(Geometry &geom, const int &face_idx, const double &eps)
{
  const vector<Vec3d> &verts = geom.verts();
  const vector<vector<int>> &faces = geom.faces();
  vector<int> face = faces[face_idx];

  Vec3d P1 = verts[face[0]];
  Vec3d P2 = verts[face[2]];
  double diag1 = (P2 - P1).len();

  P1 = verts[face[1]];
  P2 = verts[face[3]];
  double diag2 = (P2 - P1).len();

  return (double_eq(diag1, diag2, eps));
}

int detect_star_polygon(Geometry &geom, const int &face_idx)
{
  const vector<Vec3d> &verts = geom.verts();
  const vector<vector<int>> &faces = geom.faces();
  vector<int> face = faces[face_idx];

  Vec3d v1 = verts[face[1]] - verts[face[0]];
  Vec3d v2 = verts[face[1]] - verts[face[2]];
  double angle = acos(safe_for_trig(vdot(v1, v2) / (v1.len() * v2.len())));
  angle = rad2deg(angle);

  int star = 1;
  double m = face.size();
  for (double n = 2; n < m / 2; n++) {
    if ((int)(m / n) * n != m && gcd((int)m, (int)n) == 1) {
      double sp = 180 * (1 - 2 * n / m);
      // angle allowed plus/minus 2 degrees slop
      if (fabs(sp - angle) <= 2) {
        star = (int)n;
        break;
      }
    }
  }
  return (star);
}

string Vtxt(const Vec3d &v, const int &dgts)
{
  char buf[128];
  if (dgts > 0)
    snprintf(buf, 128, "%.*g, %.*g, %.*g,", dgts, v[0], dgts, v[1], dgts, v[2]);
  else
    snprintf(buf, 128, "%.*f, %.*f, %.*f,", -dgts, v[0], -dgts, v[1], -dgts,
             v[2]);
  return buf;
}

void print_hedron_txt(FILE *ofile, Geometry &geom, const int &sig_digits,
                      const bool &estimate_colors, const bool &detect_rhombi,
                      const bool &detect_star_polygons,
                      const bool &exclude_coordinates,
                      const bool &force_transparent, const double &eps)
{
  const vector<vector<int>> &faces = geom.faces();
  const vector<Vec3d> &verts = geom.verts();

  fprintf(ofile, "{\n");

  for (unsigned int i = 0; i < faces.size(); i++) {
    // fprintf(ofile, "");

    if (detect_rhombi && faces[i].size() == 4 && !is_square(geom, (int)i, eps))
      fprintf(ofile, "D");

    if (estimate_colors) {
      Color col = geom.colors(FACES).get((int)i);
      if (col.is_value())
        fprintf(ofile, "%s", estimated_color(col).c_str());
    }

    if (force_transparent) {
      fprintf(ofile, "t");
    }

    for (unsigned int j = 0; j < faces[i].size(); j++)
      fprintf(ofile, "%d,", faces[i][j]);

    int star = 1;
    if (detect_star_polygons && faces[i].size() > 4)
      star = detect_star_polygon(geom, (int)i);
    fprintf(ofile, "-%d,\n", star);
  }

  if (!exclude_coordinates) {
    fprintf(ofile, "(\n");

    for (unsigned int i = 0; i < verts.size(); i++)
      fprintf(ofile, "%s\n", Vtxt(verts[i], sig_digits).c_str());

    fprintf(ofile, ")\n");
  }

  fprintf(ofile, "}\n");
}

int main(int argc, char *argv[])
{
  o2t_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);

  FILE *ofile = stdout; // write to stdout by default
  if (opts.ofile != "") {
    ofile = fopen(opts.ofile.c_str(), "w");
    if (ofile == 0)
      opts.error("could not open output file \'" + opts.ofile + "\'");
  }

  print_hedron_txt(ofile, geom, opts.sig_digits, opts.estimate_colors,
                   opts.detect_rhombi, opts.detect_star_polygons,
                   opts.exclude_coordinates, opts.force_transparent, opts.eps);

  if (opts.ofile != "")
    fclose(ofile);

  return 0;
}
