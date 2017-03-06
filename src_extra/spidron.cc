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
   Name: spidron.cc
   Description: convert polyhedron to one built from spidron units
   Project: Antiprism - http://www.antiprism.com
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <string>
#include <vector>

#include "../base/antiprism.h"

using std::string;
using std::vector;
using std::map;
using std::swap;

using namespace anti;

class spid_opts : public ProgramOpts {
public:
  int unit_len;
  int type;
  double ang;
  double ang2;
  string pattern;
  string ifile;
  string ofile;

  spid_opts()
      : ProgramOpts("spidron"), unit_len(5), type(1), ang(M_PI / 6), ang2(0),
        pattern("01")
  {
  }
  void process_command_line(int argc, char **argv);
  void usage();
};

// clang-format off
void spid_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Read a file in OFF format and convert to model built from spidron units.\n"
"If input_file is not given the program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"  -l <len>  unit length, half the number of triangles from edge to centre\n"
"  -t <type> folding type, 1 - spidronmyd (default)\n"
"  -a <ang>  angle of first triangle (default: 30)\n"
"  -b <ang>  angle of second triangle (default: two times angle of option -a)\n"
"  -p <pat>  folding pattern, a series of 0's and 1's, to old in or out\n"
"            (default: 01)\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}
// clang-format on

void spid_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":ht:a:b:l:p:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'o':
      ofile = optarg;
      break;

    case 't':
      print_status_or_exit(read_int(optarg, &type), c);
      if (type < 1 || type > 1)
        error("type must be 1", c);
      break;

    case 'a':
      print_status_or_exit(read_double(optarg, &ang), c);
      ang *= M_PI / 180;
      break;

    case 'b':
      print_status_or_exit(read_double(optarg, &ang2), c);
      if (ang2 == 0.0)
        error("angle cannot be 0.0", c);
      ang2 *= M_PI / 180;
      break;

    case 'l':
      print_status_or_exit(read_int(optarg, &unit_len), c);
      if (unit_len < 1)
        error("unit length must be 1 or greater", c);
      break;

    case 'p':
      if (strspn(optarg, "01") != strlen(optarg))
        error("pattern can contain only 0 and 1", c);
      pattern = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (ang2 == 0.0)
    ang2 = 2 * ang;

  if (argc - optind > 1)
    error("too many arguments");

  if (argc - optind == 1)
    ifile = argv[optind];
}

int add_spidron_unit(Geometry &spid, double ang, double ang2, Vec3d &cent,
                     bool fold, Color col = Color(0))
{
  int v0idx = spid.verts().size() - 2;
  Vec3d v0 = spid.verts(v0idx);
  Vec3d v1 = spid.verts(v0idx + 1);
  Vec3d mid = (v0 + v1) / 2.0;               // Midpoint of edge
  double R = (cent - v0).len();              // Radius of polygon
  double L = (mid - v0).len();               // Half edge length
  double H = (cent - mid).len();             // Height from edge to centre
  double theta = acos(safe_for_trig(H / R)); // Half angle of edge at centre
  double l = cos(ang2) * L / cos(ang);       // Half inner edge length
  // double l = L/(2*cos(ang));             // Half inner edge length
  double r = l / sin(theta);   // Radius of inner edge
  double h_tri = L * tan(ang); // Height of base triangle
  double delta_x = H - r;

  // Are the angles suitable for this polygon
  double test_val = h_tri * h_tri - delta_x * delta_x;
  if (test_val < -epsilon)
    return false;
  else if (test_val < 0)
    test_val = 0;

  double delta_y = (1.0 - 2 * fold) * sqrt(test_val);
  Vec3d norm = vcross(v1 - v0, cent - v0).unit();
  Vec3d perp = (cent - mid).unit();
  Vec3d P = mid + perp * delta_x + norm * delta_y;
  Trans3d trans = Trans3d::translate(cent) * Trans3d::rotate(norm, theta * 2) *
                  Trans3d::translate(-cent);
  Vec3d Q = trans * P;

  /*
  fprintf(stderr, "delta_x=%g, delta_y=%g, theta=%g, H=%g, l=%g\n", delta_x,
  delta_y, theta, H, L);
  cent.dump("cent");
  norm.dump("norm");
  perp.dump("perp");
  mid.dump("mid");
  v0.dump("v0");
  v1.dump("v1");
  P.dump("P");
  Q.dump("Q");
  */

  spid.add_vert(P);
  spid.add_vert(Q);
  int f = spid.add_face(vector<int>());
  spid.colors(FACES).set(f, col);
  spid.faces(f).push_back(v0idx);
  spid.faces(f).push_back(v0idx + 1);
  spid.faces(f).push_back(v0idx + 2);
  f = spid.add_face(vector<int>());
  spid.colors(FACES).set(f, col);
  spid.faces(f).push_back(v0idx + 1);
  spid.faces(f).push_back(v0idx + 3);
  spid.faces(f).push_back(v0idx + 2);
  cent += norm * delta_y;
  return true;
}

int add_spidron_unit2(Geometry &spid, double ang, double ang2, Vec3d &cent,
                      bool fold, Color col = Color(0))
{
  int v0idx = spid.verts().size() - 2;
  Vec3d v0 = spid.verts(v0idx);
  Vec3d v1 = spid.verts(v0idx + 1);
  Vec3d norm = vcross(v1 - v0, cent - v0).unit();
  Vec3d mid = (v0 + v1) / 2.0;               // Midpoint of edge
  double R = (cent - v0).len();              // Radius of polygon
  double L = (mid - v0).len();               // Half edge length
  double H = (cent - mid).len();             // Hieght from edge to centre
  double theta = acos(safe_for_trig(H / R)); // Half angle of edge at centre

  double h = sin(ang2) * L / cos(ang); // 2nd triangle height
  double mid2r = R - h;                // radius to mid inner edgd
  // double r = mid2r/cos(theta);           // radius to inner edge
  double l = cos(ang2) * L / cos(ang); // Half inner edge length
  double delta_x = mid2r * tan(theta); // length proj of half inner edge
  double delta_y = (1.0 - 2 * fold) * sqrt(l * l - delta_x * delta_x);
  // double rot_ang = acos(lp/l);           // angle to rot 2nd triangle

  Vec3d mid2 = v1 + h * (cent - v1).unit(); // mid-point inner edge
  Vec3d perp = vcross(cent - v1, norm).unit();
  Vec3d P = mid2 - perp * delta_x + norm * delta_y;
  Vec3d Q = mid2 + perp * delta_x - norm * delta_y;

  perp = (cent - mid).unit();
  delta_x = L * tan(ang); // Height of base triangle
  delta_y = 0;
  //   Vec3d P = mid - perp*delta_x + norm*delta_y;
  //   Vec3d Q = mid + perp*delta_x - norm*delta_y;

  fprintf(stderr, "\ndelta_x=%g, delta_y=%g, theta=%g, H=%g, l=%g\n", delta_x,
          delta_y, theta, H, l);
  cent.dump("cent");
  norm.dump("norm");
  perp.dump("perp");
  mid.dump("mid");
  mid2.dump("mid2");
  v0.dump("v0");
  v1.dump("v1");
  P.dump("P");
  Q.dump("Q");

  spid.add_vert(P);
  spid.add_vert(Q);
  int f = spid.add_face(vector<int>());
  spid.colors(FACES).set(f, col);
  spid.faces(f).push_back(v0idx);
  spid.faces(f).push_back(v0idx + 1);
  spid.faces(f).push_back(v0idx + 2);
  f = spid.add_face(vector<int>());
  spid.colors(FACES).set(f, col);
  spid.faces(f).push_back(v0idx + 1);
  spid.faces(f).push_back(v0idx + 3);
  spid.faces(f).push_back(v0idx + 2);

  Trans3d trans = Trans3d::translate(cent) * Trans3d::rotate(norm, theta * 4) *
                  Trans3d::translate(-cent);
  Vec3d P2 = trans * P;
  Vec3d V0 = (P - Q).unit();
  Vec3d V1 = (P2 - Q).unit();
  double beta = acos(safe_for_trig(vdot(V0, V1)));
  double slant = l / cos(beta / 2);
  Vec3d rad_dir = ((P + P2) / 2.0 - Q).unit();
  cent = Q + rad_dir * slant;
  fprintf(stderr, "\nbeta=%g, slant=%g\n", beta * 180 / M_PI, slant);
  cent.dump("cent");
  P2.dump("P2");
  // spid.add_vert(cent);
  return 1;
}

int make_spidron(Geometry &spid, Geometry &base, double ang, double ang2,
                 int len, vector<bool> folds, int type)
{
  for (unsigned int i = 0; i < base.faces().size(); i++) {
    int f_sz = base.faces(i).size();
    if (f_sz < 3)
      continue;
    Vec3d cent = base.face_cent(i);
    for (int j = 0; j < f_sz; j++) {
      int v0 = base.faces(i, j);
      int v1 = base.faces(i, (j + 1) % f_sz);
      spid.add_vert(base.verts(v0));
      spid.add_vert(base.verts(v1));
      Vec3d edge_cent = cent;
      Color col(0, 0, 0, 0);
      vector<int> edge(2);
      edge[0] = v0;
      edge[1] = v1;
      if (edge[0] > edge[1])
        swap(edge[0], edge[1]);
      vector<vector<int>>::const_iterator ei;
      const vector<vector<int>> &edges = base.edges();
      ei = find(edges.begin(), edges.end(), edge);
      if (ei != edges.end()) {
        unsigned int e_idx = ei - edges.begin();
        col = base.colors(EDGES).get(e_idx);
      }
      for (int l = 0; l < len; l++) {
        if (type == 1) {
          if (!add_spidron_unit(spid, ang, ang2, edge_cent,
                                folds[l % folds.size()], col))
            return false;
        }
        if (type == 2) {
          bool fold = folds[l % folds.size()];
          if (!is_even(j))
            fold = !fold;
          add_spidron_unit2(spid, ang, ang2, edge_cent, fold, col);
        }
      }
    }
  }
  return 1;
}

double v_ang_at_ax(const Vec3d &v0, const Vec3d &v1, const Vec3d &ax)
{
  Vec3d n0 = vcross(v0, ax).unit();
  Vec3d n1 = vcross(v1, ax).unit();
  double ang = acos(safe_for_trig(vdot(n0, n1)));
  if (vdot(ax, vcross(n0, n1)) < 0)
    ang = 2 * M_PI - ang;
  return ang;
}

int main(int argc, char *argv[])
{
  spid_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);
  geom.add_missing_impl_edges();

  vector<bool> fold(opts.pattern.size());
  for (unsigned int f = 0; f < opts.pattern.size(); f++)
    fold[f] = opts.pattern[f] == '0';

  Geometry spid;
  if (!make_spidron(spid, geom, opts.ang, opts.ang2, opts.unit_len, fold,
                    opts.type))
    opts.error("unable to make spidron units with these angles for at least "
               "one of the polygon faces");

  opts.write_or_error(spid, opts.ofile);

  return 0;
}
