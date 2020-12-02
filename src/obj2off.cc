/*
   Copyright (c) 2014-2020, Roger Kaufman

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
   Name: obj2off.cc
   Description: Convert files in OBJ format to OFF file format
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"

#include <cctype>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using std::string;
using std::vector;

using namespace anti;

class obj2off_opts : public ProgramOpts {
public:
  int sig_digits;
  string ifile;
  string ofile;

  obj2off_opts() : ProgramOpts("obj2off"), sig_digits(DEF_SIG_DGTS) {}

  void process_command_line(int argc, char **argv);
  void usage();
};

void obj2off_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] [input_file]

Convert files in OBJ format to OFF format. Only v, e and f statements are
processed. If input_file is not given the program reads from standard input.

Options
%s
  -d <dgts> number of significant digits (default %d) or if negative
            then the number of digits after the decimal point
  -o <file> write output to file (default: write to standard output)

)",
          prog_name(), help_ver_text, DEF_SIG_DGTS);
}

void obj2off_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hd:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {

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

bool convert_obj_to_off(string &file_name, Geometry &geom,
                        string *error_msg = nullptr)
{
  FILE *ifile;
  if (file_name == "" || file_name == "-") {
    ifile = stdin;
    file_name = "stdin";
  }
  else {
    ifile = fopen(file_name.c_str(), "r");
    if (!ifile) {
      if (error_msg)
        *error_msg =
            msg_str("could not open input file '%s'", file_name.c_str());
      return false;
    }
  }

  int offset = 1; // obj files start indexes from 1

  char *line = nullptr;
  while (read_line(ifile, &line) == 0) {
    for (char *p = line; *p; p++) // convert whitespace to spaces
      if (isspace(*p))
        *p = ' ';

    Split parts(line, " ");
    unsigned int parts_sz = parts.size();

    char key = '\0';
    if (strlen(parts[0]))
      key = parts[0][0];

    if (key == 'v') {
      // only use x y z
      double coord[3];
      for (unsigned int i = 1; i < parts_sz; i++) {
        if (!read_double(parts[i], &coord[i - 1])) {
          if (error_msg)
            *error_msg = msg_str("invalid coordinate '%s'", parts[i]);
          return false;
        }
      }
      geom.add_vert(Vec3d(coord[0], coord[1], coord[2]));
    }
    else if (key == 'f' || key == 'l') {
      int idx;
      vector<int> indexes;
      for (unsigned int i = 1; i < parts_sz; i++) {
        if (!read_int(parts[i], &idx)) {
          if (error_msg)
            *error_msg = msg_str("invalid face or edge index '%s'", parts[i]);
          return false;
        }
        indexes.push_back(idx - offset);
      }
      if (key == 'f')
        geom.add_face(indexes);
      else if (key == 'l')
        geom.add_edge(indexes);
    }
    /* RK - The p key only outputs sequential vertex color map numbers
        else if (key == 'p') {
          int idx;
          if (!read_int(parts[1], &idx)) {
            if (error_msg)
              *error_msg = msg_str("invalid vertex color index '%s'", parts[1]);
            return false;
          }
          geom.colors(VERTS).set(idx - offset, idx);
        }
    */

    free(line);
  }

  if (file_name != "stdin")
    fclose(ifile);

  return true;
}

int main(int argc, char *argv[])
{
  obj2off_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;

  string error_msg;

  // obj is enough like OFF that it can be parsed and converted in line
  if (!convert_obj_to_off(opts.ifile, geom, &error_msg))
    if (!error_msg.empty())
      opts.error(error_msg);

  opts.write_or_error(geom, opts.ofile, opts.sig_digits);

  return 0;
}
