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
   Name: conv_hull.cc
   Description: convex hulls (wrapper for qhull)
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <string>

using std::string;

using namespace anti;

class ch_opts : public ProgramOpts {
public:
  bool append_flg;
  string ifile;
  string ofile;
  string qh_args;

  ch_opts() : ProgramOpts("conv_hull"), append_flg(false) {}
  void process_command_line(int argc, char **argv);
  void usage();
};

void ch_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] [input_file]

Read a file in OFF format and make a convex hull (using Qhull). If
input_file is not given the program reads from standard input.

Options
%s
  -a        append the convex hull to the input file
  -Q <args> additional arguments to pass to qhull (unsupported, may not
            work, check output)
  -o <file> write output to file (default: write to standard output)

)",
          prog_name(), help_ver_text);
}

void ch_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":haQ:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'a':
      append_flg = true;
      break;

    case 'o':
      ofile = optarg;
      break;

    case 'Q':
      qh_args = optarg;
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

int main(int argc, char *argv[])
{
  ch_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);

  int dimension;
  Status stat = (opts.append_flg) ? geom.add_hull(opts.qh_args, &dimension)
                                  : geom.set_hull(opts.qh_args, &dimension);

  if (stat.is_error())
    opts.error(stat.msg());

  const char *dimension_desc[] = {"point", "line segment", "polygon"};
  if (dimension < 3)
    opts.warning(msg_str("result is a %s", dimension_desc[dimension]));

  opts.write_or_error(geom, opts.ofile);

  return dimension;
}
