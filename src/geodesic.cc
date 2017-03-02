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
   Name: geodesic.cc
   Description: program to make geodesic spheres and polyhedra in OFF format
   Project: Antiprism - http://www.antiprism.com
*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>

#include "../base/antiprism.h"

using std::string;

using namespace anti;

class geo_opts : public ProgramOpts {
public:
  Vec3d centre;
  double radius;
  int m;
  int n;
  int pat_freq;
  bool use_step_freq;
  char method;
  bool keep_flat;
  bool equal_len_div;
  string ifile;
  string ofile;

  geo_opts()
      : ProgramOpts("geodesic"), centre(Vec3d(0, 0, 0)), m(1), n(0),
        pat_freq(1), use_step_freq(false), method('s')
  {
  }
  void process_command_line(int argc, char **argv);
  void usage();
};

// clang-format off
void geo_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Read a file in OFF format and make a higher frequency, plane-faced\n"
"polyhedron or geodesic sphere. If input_file is not given the program\n"
"reads from standard input\n"

"\n"
"Options\n"
"%s"
"  -f <freq> pattern frequency, a positive integer (default: 1) giving the\n"
"            number of repeats of the specified pattern along an edge\n"
"  -F <freq> final step frequency, minimum number of edges to move between\n"
"            base vertices in the geodesic model. For a pattern m,n the\n"
"            step frequency is pattern_frequency/(m+n)\n"
"  -c <clss> face division pattern,  1 (Class I, default), 2 (Class II), or\n"
"            two numbers separated by a comma to determine the pattern\n"
"            (Class III, but n,0 or 0,n is Class I, and n,n is Class II)\n"
"  -M <mthd> Method of applying the frequency pattern:\n"
"            s - geodesic sphere (default). The pattern grid is formed\n"
"                from divisions along each edge that make an equal angle\n"
"                at the centre. The geodesic vertices are centred at the\n"
"                origin and projected on to a unit sphere.\n"
"            p - planar. The pattern grid is formed from equal length\n"
"                divisions along each edge, the new vertices lie on the\n"
"                surface of the original polyhedron.\n"
"  -C <cent> centre of points, in form \"x_val,y_val,z_val\" (default: 0,0,0)\n"
"            used for geodesic spheres\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}
// clang-format on

void geo_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hf:F:c:M:C:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'f':
    case 'F':
      if (!read_int(optarg, &pat_freq) || pat_freq < 1)
        error(
            msg_str("frequency is '%s', should be a positive integer", optarg),
            c);
      use_step_freq = (c == 'F');
      break;

    case 'c':
      char *pos;
      if (!(pos = strchr(optarg, ','))) {
        if (!(strcmp(optarg, "1") == 0 || strcmp(optarg, "2") == 0))
          error("class type must be 1 or 2, or two integers "
                "separated by a comma\n",
                c);
        m = 1;
        n = *optarg == '2' ? 1 : 0;
      }
      else {
        *pos = 0;
        char *pat[] = {optarg, pos + 1};
        int nums[2];
        for (int i = 0; i < 2; i++) {
          if (!read_int(pat[i], &nums[i]) || nums[i] < 0)
            error(msg_str("%s value is '%s' should be a positive "
                          "integer or zero",
                          (i == 0) ? "first" : "second", pat[i]),
                  c);
        }
        if (nums[0] == 0 && nums[1] == 0)
          error("pattern must include at least one positive integer", c);

        m = nums[0];
        n = nums[1];
      }
      break;

    case 'M':
      if (strlen(optarg) == 1 && strchr("sp", int(*optarg)))
        method = *optarg;
      else
        error("method type must be s or p", c);
      break;

    case 'C':
      if (!centre.read(optarg))
        error(msg_str("centre is '%s', must be three numbers separated "
                      "by commas",
                      optarg),
              c);
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

  if (use_step_freq) {
    if (pat_freq % (m + n)) {
      if (m == 1 && n == 1)
        error("frequency must be even with Class II pattern", "F");
      else
        error("frequency must be divisible by (m+n) with general "
              "Class III pattern",
              "F");
    }
    pat_freq /= (m + n);
  }

  m *= pat_freq;
  n *= pat_freq;

  return;
}

int main(int argc, char **argv)
{
  geo_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);

  Geometry geo;
  if (opts.method == 's')
    make_geodesic_sphere(&geo, geom, opts.m, opts.n, opts.centre);
  else if (opts.method == 'p')
    make_geodesic_planar(&geo, geom, opts.m, opts.n);

  opts.write_or_error(geo, opts.ofile);

  return 0;
}
