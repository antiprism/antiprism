/*
   Copyright (c) 2003-2016, Adrian Rossiter

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
   Name: sph_saff.cc
   Description:
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

using std::string;

using namespace anti;

class sphc_opts : public ProgramOpts {
public:
  int num_rings;
  int num_divs;
  bool stagger;
  bool pts_at_poles;

  string ofile;

  sphc_opts()
      : ProgramOpts("sph_ring"), num_rings(10), num_divs(0), stagger(false),
        pts_at_poles(true)
  {
  }
  void process_command_line(int argc, char **argv);
  void usage();
};

void sphc_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] num_rings num_divs

Distribute points on rings a sphere. num_rings is the number of
rings between the poles. num_divs is the number of points on each
ring, or if not given then it is as many points as will fit in each
ring at least a band's width apart.

Options
%s
  -s        stagger placement of balls between cirles
  -x        don't place points at the two poles
  -o <file> write output to file (default: write to standard output)

)",
          prog_name(), help_ver_text);
}

void sphc_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hsxo:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 's':
      stagger = true;
      break;

    case 'x':
      pts_at_poles = false;
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (argc - optind > 0) {
    print_status_or_exit(read_int(argv[optind], &num_rings), "num_rings");
    if (num_rings < 1)
      error("number of rings must be a positive integer");
  }

  if (argc - optind > 1) {
    print_status_or_exit(read_int(argv[optind + 1], &num_divs), "num_pts");
    if (num_divs < 1)
      error("number of divisions must be a positive integer");
  }

  if (argc - optind > 2)
    error("too many arguments");
}

void make_sph_rings(Geometry &geom, int rings, int divs, bool stagger,
                    bool pts_at_poles)
{
  if (pts_at_poles) {
    geom.add_vert(Vec3d(0, 1, 0));
    geom.add_vert(Vec3d(0, -1, 0));
  }
  double horz_stagger = 0;
  double vert_ang_inc = M_PI / (rings + 1);
  double ball_dia = sqrt(2 - 2 * cos(vert_ang_inc));
  for (int i = 1; i < rings + 1; i++) {
    double num_pts = divs;
    double vert_ang = i * vert_ang_inc;
    double rad = sin(vert_ang);
    if (!divs) {
      double val = (2 * rad * rad - ball_dia * ball_dia) / (2 * rad * rad);
      if (val < -1 || val > 1)
        num_pts = 1;
      else
        num_pts = floor(2 * M_PI / acos(safe_for_trig(val)));
    }

    double horz_ang_inc = 2 * M_PI / num_pts;
    for (int j = 0; j < num_pts; j++) {
      double horz_ang = j * horz_ang_inc + horz_stagger;
      Vec3d P(rad * sin(horz_ang), cos(vert_ang), rad * cos(horz_ang));
      geom.add_vert(P);
    }

    if (stagger)
      horz_stagger = horz_ang_inc * (i % 2) / 2;
  }
}

int main(int argc, char **argv)
{
  sphc_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  make_sph_rings(geom, opts.num_rings, opts.num_divs, opts.stagger,
                 opts.pts_at_poles);

  opts.write_or_error(geom, opts.ofile);

  return 0;
}
