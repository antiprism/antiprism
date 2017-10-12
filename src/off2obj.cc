/*
   Copyright (c) 2014-2016, Roger Kaufman

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
   Name: off2obj.cc
   Description: extract coordinates from an OFF file
   Project: Antiprism - http://www.antiprism.com
*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <vector>

#include "../../antiprism_rk/base/antiprism.h"

using std::string;
using std::vector;

using namespace anti;

class o2o_opts : public ProgramOpts {
public:
  string ifile;
  string ofile;

  string mtl_file;
  string sep;
  int sig_digits;

  o2o_opts() : ProgramOpts("off2obj"), sep(" "), sig_digits(DEF_SIG_DGTS) {}
  void process_command_line(int argc, char **argv);
  void usage();
};

// clang-format off
void o2o_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Convert an OFF file to Wavefront OBJ file format.\n"
"\n"
"Options\n"
"%s"
"  -m <file> generate mtl file. file name is hardcoded into obj file\n"
"  -d <dgts> number of significant digits (default %d) or if negative\n"
"            then the number of digits after the decimal point\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text, DEF_SIG_DGTS);
}
// clang-format on

void o2o_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hm:o:d:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'm':
      mtl_file = optarg;
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

int main(int argc, char *argv[])
{
  o2o_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);

  opts.print_status_or_exit(geom.write_obj(opts.ofile, opts.mtl_file,
                                           opts.sep.c_str(), opts.sig_digits));

  return 0;
}
