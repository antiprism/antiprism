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
   Name: off2vrml.cc
   Description: convert an OFF file to VRML format
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"
#include <math.h>
#include <string.h>
#include <string>
#include <vector>

using std::string;
using std::vector;

using namespace anti;

class o2v_opts : public ViewOpts {
public:
  int sig_dgts;
  string ofile;

  o2v_opts() : ViewOpts("off2vrml"), sig_dgts(DEF_SIG_DGTS) {}

  void process_command_line(int argc, char **argv);
  void usage();
};

// clang-format off
void o2v_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] input_files\n"
"\n"
"Convert files in OFF format to VRML format. If input_files are not\n"
"given the program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"%s"
"  -l        use lines for edges, points for vertices, in default colours\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"  Scene options\n"
"%s"
"\n"
"  Precision options\n"
"%s"
"\n"
"\n", prog_name(), help_ver_text, help_view_text,
      help_scene_text, help_prec_text);
}
// clang-format on

void o2v_opts::process_command_line(int argc, char **argv)
{
  Status stat;
  opterr = 0;
  int c;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv,
                     ":hv:e:V:E:F:m:x:n:s:lo:D:C:L:R:P:I:B:d:t:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'o':
      ofile = optarg;
      break;

    case 'l':
      get_geom_defs().set_use_lines(true);
      break;

    case 'd':
      print_status_or_exit(read_int(optarg, &sig_dgts), c);
      if (sig_dgts < 1)
        error("number of significant digits must be a positive integer", c);
      break;

    default:
      if (!(stat = read_disp_option(c, optarg))) {
        if (stat.is_warning())
          warning(stat.msg(), c);
        else
          error(stat.msg(), c);
      }
    }
  }

  if (argc - optind >= 1)
    while (argc - optind >= 1)
      ifiles.push_back(argv[optind++]);
  else
    ifiles.push_back("");
}

int main(int argc, char *argv[])
{
  o2v_opts opts;
  opts.process_command_line(argc, argv);
  Scene scen = opts.scen_defs;
  opts.set_view_vals(scen);

  FILE *ofile = stdout; // write to stdout by default
  if (opts.ofile != "") {
    ofile = fopen(opts.ofile.c_str(), "w");
    if (ofile == nullptr)
      opts.error("could not open output file \'" + opts.ofile + "\'");
  }

  VrmlWriter vrml;
  vrml.write(ofile, scen, opts.sig_dgts);

  if (opts.ofile != "")
    fclose(ofile);

  return 0;
}
