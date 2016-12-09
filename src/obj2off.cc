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
   Name: obj2off.cc
   Description: Convert files in OBJ format to OFF file format
   Project: Antiprism - http://www.antiprism.com
*/

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include <ctype.h>

#include <string>
#include <vector>

#include "../base/antiprism.h"

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
  fprintf(
      stdout,
      "\n"
      "Usage: %s [options] [input_file]\n"
      "\n"
      "Convert files in OBJ format to OFF format. Only v, e and f statements "
      "are\n"
      "processed. If input_file is not given the program reads from standard "
      "input.\n"
      "\n"
      "Options\n"
      "%s"
      "  -d <dgts> number of significant digits (default %d) or if negative\n"
      "            then the number of digits after the decimal point\n"
      "  -o <file> write output to file (default: write to standard output)\n"
      "\n"
      "\n",
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

void convert_obj_to_off(string &file_name, Geometry &geom, char *errmsg)
{
  if (errmsg)
    *errmsg = '\0';

  FILE *ifile;
  if (file_name == "" || file_name == "-") {
    ifile = stdin;
    file_name = "stdin";
  }
  else {
    ifile = fopen(file_name.c_str(), "r");
    if (!ifile) {
      if (errmsg)
        snprintf(errmsg, MSG_SZ, "could not open input file \'%s\'",
                 file_name.c_str());
    }
  }

  char parse_key[] = " ";

  int offset = 1; // obj files start indexes from 1

  char *line = nullptr;
  while (read_line(ifile, &line) == 0) {
    for (char *p = line; *p; p++) // convert whitespace to spaces
      if (isspace(*p))
        *p = ' ';

    char *ptok = strtok(line, parse_key);
    char key = ptok ? *ptok : '\0';
    int idx;

    // only use x y z
    if (key == 'v') {
      double coord[3];
      for (double &j : coord) {
        ptok = strtok(nullptr, parse_key);
        sscanf(ptok, "%lf", &j);
      }
      geom.add_vert(Vec3d(coord[0], coord[1], coord[2]));
    }
    else if (key == 'f' || key == 'l') {
      vector<int> indexes;
      ptok = strtok(nullptr, parse_key);
      while (ptok != nullptr) {
        sscanf(ptok, "%d", &idx);
        indexes.push_back(idx - offset);
        ptok = strtok(nullptr, parse_key);
      }
      if (key == 'f')
        geom.add_face(indexes);
      else if (key == 'l')
        geom.add_edge(indexes);
    }
    else if (key == 'p') {
      ptok = strtok(nullptr, parse_key);
      sscanf(ptok, "%d", &idx);
      geom.colors(VERTS).set(idx - offset, Color());
    }

    free(line);
  }

  if (file_name != "stdin")
    fclose(ifile);
}

int main(int argc, char *argv[])
{
  obj2off_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;

  char errmsg[MSG_SZ];
  // obj is enough like OFF that it can be parsed and converted in line
  convert_obj_to_off(opts.ifile, geom, errmsg);
  if (*errmsg)
    opts.error(errmsg);

  opts.write_or_error(geom, opts.ofile, opts.sig_digits);

  return 0;
}
