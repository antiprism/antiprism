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

/*
   Name: kcycle.cc
   Description: make a kaliedocycle with a polyhedron
   Project: Antiprism - http://www.antiprism.com
*/

#include <algorithm>
#include <cctype>
#include <map>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>

#include "../base/antiprism.h"

using std::string;
using std::vector;

using namespace anti;

class kc_opts : public ProgramOpts {
private:
public:
  vector<int> edge_idxs;
  double angle;
  int num_prs;

  string ifile;
  string ofile;

  kc_opts() : ProgramOpts("kcycle"), angle(0.0), num_prs(6) {}

  void process_command_line(int argc, char **argv);
  void usage();
};

// clang-format off
void kc_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] hinge_vertex_idxs\n"
"\n"
"Make a kaleidocyle from a polyhedron. hinge_vertex_idxs is four integers\n"
"separated by commas (e.g. 0,1,2,3). These are the index numbers of the\n"
"vertices for the two eges that will be hinges in the cycle.\n"
"\n"
"Options\n"
"%s"
"  -a <ang>  angle in degrees to rotate the first hinge from\n"
"            horizontal (default: 0.0)\n"
"  -n <num>  number of pairs of polyhedra in cycle (default: 6)\n"
"  -i <file> input file in OFF format. If '-' then read file from stdin.\n"
"            (default: a tetrahedron with hinge_vertex_idxs 0,1,2,3)\n" 
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}
// clang-format on

void kc_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":ha:n:i:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'a':
      print_status_or_exit(read_double(optarg, &angle), c);
      angle = deg2rad(angle);
      break;

    case 'n':
      print_status_or_exit(read_int(optarg, &num_prs), c);
      if (num_prs < 1)
        error("must be at least two pairs", c);
      break;

    case 'i':
      ifile = optarg;
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (argc - optind > 1) {
    error("too many arguments (leave no spaces between the edge indexes)");
    exit(1);
  }

  if (argc - optind < 1) {
    if (ifile == "") {
      for (int i = 0; i < 4; i++)
        edge_idxs.push_back(i);
    }
    else
      error("edge index numbers not given");
  }

  if (argc - optind == 1)
    print_status_or_exit(read_int_list(argv[optind], edge_idxs, true),
                         "edge index numbers");

  if (edge_idxs.size() != 4)
    error(msg_str("exacly four numbers must be specified (%lu given)",
                  (unsigned long)edge_idxs.size()),
          "edge index numbers");

  if (edge_idxs[0] == edge_idxs[1] || edge_idxs[2] == edge_idxs[3])
    error(msg_str("%s pair has repeated index number",
                  (edge_idxs[0] == edge_idxs[1]) ? "first" : "second"),
          "edge index numbers");
}

Vec3d lines_intersection(Vec3d P, Vec3d P_dir, Vec3d Q, Vec3d Q_dir)
{
  Vec3d N1, N2;
  if (lines_nearest_points(P, P + P_dir, Q, Q + Q_dir, N1, N2))
    return (N1 + N2) / 2.0;
  else
    return Vec3d();
}

void kcycle(Geometry geom, Geometry &cycle, int num_prs, vector<int> idxs,
            double angle)
{

  const vector<Vec3d> &verts = geom.verts();
  Vec3d P[] = {verts[idxs[0]], verts[idxs[1]], verts[idxs[2]], verts[idxs[3]]};
  Vec3d e1(cos(angle), sin(angle), 0); // direction of first edge

  Trans3d trans1 = Trans3d::rotate(P[1] - P[0], e1) * Trans3d::translate(-P[0]);
  for (auto &i : P)
    i = trans1 * i;

  Vec3d e2 = (P[3] - P[2]).unit(); // current direction of second edge

  double wedge_angle = -M_PI / num_prs;
  Vec3d n2(sin(wedge_angle), 0, cos(wedge_angle)); // normal to other slice

  // Use spherical trig to find the angle to rotate e2 about e1 until
  // it is normal to n2
  int sign1 = vtriple(e1, e2, n2) > 0;
  double ang1 =
      acos(safe_for_trig((vdot(e2, n2) - vdot(e1, n2) * vdot(e1, e2)) /
                         (vcross(e1, n2).len() * vcross(e1, e2).len())));
  int sign2 = vdot(e1, e2) > -epsilon;
  double ang2 =
      acos(safe_for_trig((-vdot(e1, n2) * vdot(e1, e2)) /
                         (vcross(e1, n2).len() * vcross(e1, e2).len())));
  double ang = (2 * sign1 - 1) * ang1 - (2 * sign2 - 1) * ang2;

  Trans3d trans2 = Trans3d::rotate(e1, ang);
  for (auto &i : P)
    i = trans2 * i;

  Vec3d N0, N1;
  lines_nearest_points(P[0], P[1], P[2], P[3], N0, N1);
  Vec3d perp = N1 - N0;
  double dist = -perp[0] - perp[2] / tan(wedge_angle);

  Trans3d trans3 = Trans3d::translate(Vec3d(dist, 0, 0) - N0) * trans2 * trans1;

  // align for symmetry exp
  trans3 = Trans3d::rotate(Vec3d(0, 1, 0), Vec3d(0, 0, 1)) *
           Trans3d::rotate(Vec3d(0, 1, 0), M_PI / num_prs) * trans3;
  for (unsigned int i = 0; i < verts.size(); i++)
    geom.verts(i) = trans3 * geom.verts(i);

  // align centroid with z=0
  Vec3d cent = centroid(verts);
  Trans3d trans_c = Trans3d::translate(Vec3d(0, 0, -cent[2]));
  for (unsigned int i = 0; i < verts.size(); i++)
    geom.verts(i) = trans_c * geom.verts(i);

  // make the cycle
  sym_repeat(cycle, geom, Symmetry(Symmetry::Cv, num_prs));
}

int main(int argc, char *argv[])
{
  kc_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  if (opts.ifile != "")
    opts.read_or_error(geom, opts.ifile);
  else
    geom.read_resource("std_tet");

  Geometry cycle;
  kcycle(geom, cycle, opts.num_prs, opts.edge_idxs, opts.angle);

  opts.write_or_error(cycle, opts.ofile);

  return 0;
}
