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
   Name: canonical.cc
   Description: canonicalize a polyhedron
                Uses George Hart's two canonicalization algorithm
                http://library.wolfram.com/infocenter/Articles/2012/
                http://www.georgehart.com/virtual-polyhedra/conway_notation.html
   Project: Antiprism - http://www.antiprism.com
*/

#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>

#include "../base/antiprism.h"

using std::string;
using std::vector;

using namespace anti;

class cn_opts : public ProgramOpts {
public:
  string ifile;
  string ofile;

  int num_iters;
  double edge_factor;
  double plane_factor;
  char method;
  char cent_type;
  int num_iters_preplanar;
  int rep_count;
  int divergence_test;

  double epsilon;

  cn_opts()
      : ProgramOpts("canonical"), num_iters(-1), edge_factor(50),
        plane_factor(20), method('n'), cent_type('\0'),
        num_iters_preplanar(1000), rep_count(50), divergence_test(10),
        epsilon(0)
  {
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

void cn_opts::usage()
{
  fprintf(
      stdout,
      "\n"
      "Usage: %s [options] [input_file]\n"
      "\n"
      "Read a polyhedron from a file in OFF format. Canonicalize or planarize "
      "it.\n"
      "Uses algorithms by George W. Hart, http://www.georgehart.com/\n"
      "http://www.georgehart.com/virtual-polyhedra/conway_notation.html\n"
      "http://www.georgehart.com/virtual-polyhedra/canonical.html\n"
      "If input_file is not given the program reads from standard input.\n"
      "\n"
      "Options\n"
      "%s"
      "  -n <itrs> maximum number of iterations (default: no limit)\n"
      "  -l <lim>  minimum distance change to terminate, as negative exponent\n"
      "               (default: %d giving %.0e)\n"
      "  -d <int>  divergence test. 0 for no test. (default 10)\n"
      "  -M <mthd> canonicalizing method,\n"
      "            m - mathematica version of canonicalization (default)\n"
      "            n - conway notation version of canonicalization\n"
      "            l - mathematica planarize portion only\n"
      "            p - conway notation planarize (face centroids reciprocal)\n"
      "            q - conway notation planarize (face centroids magnitude "
      "reciprocal)\n"
      "            x - face centroids only (no reciprocal) planarize method\n"
      "  -C <cent> initial 'centering'\n"
      "            x - none, c - centroid (-M p and -M l default)\n"
      "            s - centroid and project vertices onto a sphere (-M m "
      "default)\n"
      "            p - centroid and pre-planarized (-M n default)\n"
      "            q - centroid and pre-planarized with magnitude reciprocal\n"
      "  -z <n>    status reporting every n lines. -1 for no status. (default "
      "50)\n"
      "  -o <file> write output to file (default: write to standard output)\n"
      "\n"
      "Mathematica Canonicalize Options (-M m and -M l)\n"
      "  -e <perc> percentage to scale the edge tangency error (default: 50)\n"
      "  -p <perc> percentage to scale the face planarity error (default: 20)\n"
      "\n"
      "Pre-planarization Options (-C p and -C q)\n"
      "  -i <itrs> maximum number of pre-planarize iterations (default: 1000)\n"
      "\n"
      "\n",
      prog_name(), help_ver_text, int(-log(::epsilon) / log(10) + 0.5),
      ::epsilon);
}

void cn_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  int sig_compare = INT_MAX;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hn:e:p:l:d:M:C:i:z:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'o':
      ofile = optarg;
      break;

    case 'n':
      print_status_or_exit(read_int(optarg, &num_iters), c);
      if (num_iters < 0)
        error("number of iterations must be 0 or greater", c);
      break;

    case 'd':
      print_status_or_exit(read_int(optarg, &divergence_test), c);
      if (divergence_test < 0)
        error("divergence test must be 0 or greater", c);
      break;

    case 'l':
      print_status_or_exit(read_int(optarg, &sig_compare), c);
      if (sig_compare < 0) {
        warning("limit is negative, and so ignored", c);
      }
      if (sig_compare > DEF_SIG_DGTS) {
        warning("limit is very small, may not be attainable", c);
      }
      break;

    case 'e':
      print_status_or_exit(read_double(optarg, &edge_factor), c);
      if (edge_factor <= 0 || edge_factor >= 100)
        warning("not inside range 0 to 100", c);
      break;

    case 'p':
      print_status_or_exit(read_double(optarg, &plane_factor), c);
      if (plane_factor <= 0 || plane_factor >= 100) {
        warning("not inside range 0 to 100", c);
      }
      break;

    case 'M':
      if (strlen(optarg) == 1 && strchr("mnlpqx", int(*optarg)))
        method = *optarg;
      else
        error("method type must be m, n, l, p, q or x", c);
      break;

    case 'C':
      if (strlen(optarg) != 1)
        error("initial centering must be exactly one character", c);
      if (!strchr("xcspq", *optarg))
        error("initial centering must be x, c, s, p or q", c);
      cent_type = *optarg;
      break;

    case 'i':
      print_status_or_exit(read_int(optarg, &num_iters_preplanar), c);
      if (num_iters_preplanar <= 0)
        error(
            "number of iterations for preplanarization must be greater than 0",
            c);
      break;

    case 'z':
      print_status_or_exit(read_int(optarg, &rep_count), c);
      if (rep_count < -1)
        error("number of iterations must be -1 or greater", c);
      break;

    default:
      error("unknown command line error");
    }
  }

  if (!cent_type) {
    if (method == 'c')
      cent_type = 's';
    else if (method == 'n')
      cent_type = 'p';
    else if (method == 'l' || method == 'p' || method == 'q')
      cent_type = 'c';
    else
      cent_type = 'x';
  }

  if (argc - optind > 1)
    error("too many arguments");

  if (argc - optind == 1)
    ifile = argv[optind];

  epsilon = (sig_compare != INT_MAX) ? pow(10, -sig_compare) : ::epsilon;
}

void centroid_to_origin(Geometry &geom)
{
  geom.transform(Trans3d::transl(-centroid(geom.verts())));
}

int main(int argc, char *argv[])
{
  cn_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);

  if (opts.cent_type == 'c' || opts.cent_type == 'p' || opts.cent_type == 'q' ||
      opts.cent_type == 's')
    centroid_to_origin(geom);

  if (opts.cent_type == 'p' || opts.cent_type == 'q')
    canonicalize_cn(&geom, opts.num_iters_preplanar, opts.cent_type,
                    opts.divergence_test, opts.rep_count, opts.epsilon);
  else if (opts.cent_type == 's')
    project_onto_sphere(&geom);

  if (opts.method == 'm' || opts.method == 'l') {
    bool planarize_only = (opts.method == 'l') ? true : false;
    canonicalize_mm(&geom, opts.edge_factor / 100, opts.plane_factor / 100,
                    opts.num_iters, opts.divergence_test, opts.rep_count,
                    planarize_only, opts.epsilon);
  }
  else {
    // conway notation canonicalize expects 'c' for canonicalize
    canonicalize_cn(&geom, opts.num_iters, opts.method, opts.divergence_test,
                    opts.rep_count, opts.epsilon);
  }

  opts.write_or_error(geom, opts.ofile);

  return 0;
}
