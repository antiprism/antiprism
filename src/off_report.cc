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
   Name: off_report.cc
   Description: analyse an off_file and print a report
   Project: Antiprism - http://www.antiprism.com
*/

#include <ctype.h>
#include <math.h>
#include <set>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>

#include "../base/antiprism.h"
#include "rep_print.h"

using std::string;
using std::vector;
using std::set;

using namespace anti;

class or_opts : public ProgramOpts {
public:
  Vec3d center;
  bool center_is_centroid;
  int sig_digits;
  string sections;
  string counts;
  bool orient;
  bool detect_symmetry;
  string sub_sym;
  char edge_type;
  string ifile;
  string ofile;

  or_opts()
      : ProgramOpts("off_report"), center(Vec3d(0, 0, 0)),
        center_is_centroid(false), sig_digits(17), orient(true),
        detect_symmetry(false), edge_type('a')
  {
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

// clang-format off
void or_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Read a file in OFF format and generate a report\n"
"\n"
"Options\n"
"%s"
"  -c <cent> centre of shape in form 'X,Y,Z', or C to use\n"
"            centroid (default, '0,0,0')\n"
"  -S <secs> Print values by sections, given as a list of letters\n"
"            A - all                     G - general\n"
"            F - faces                   E - edges\n"
"            S - solid Angles            a - plane angles\n"
"            D - distances (min/max)     s - symmetry\n"
"  -C <vals> Print counts of values, given as a list of letters\n"
"            A - All                     F - faces type by angles\n"
"            E - edge lengths            w - windings\n"
"            S - solid angles            D - dihedral angles\n"
"            s - face sides              o - vertex orders\n"
"            h - vertex heights (z-crds) O - symmetry orbits\n"
"  -k        keep orientation, don't try to orient the faces\n"
"  -y        subsymmetry for orbits: symmetry subgroup (Schoenflies notation)\n"
"            optionally followed by a comma and conjugation type (integer)\n"
"  -E <type> edges for report, e - explicit edges, i - implicit edges\n"
"            a - explicit and implicit (default)\n"
"  -o <file> write output to file (default: write to standard output)\n"
"  -d <dgts> number of significant digits (default 17) or if negative\n"
"            then the number of digits after the decimal point\n"
"\n"
"\n", prog_name(), help_ver_text);
}
// clang-format on

void or_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hc:S:C:kE:y:o:d:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'c':
      if (strcmp(optarg, "C") == 0)
        center_is_centroid = true;
      else
        print_status_or_exit(center.read(optarg), c);
      break;

    case 'S': {
      const char *all_section_letters = "AGFEaSsD";
      size_t len;
      if ((len = strspn(optarg, all_section_letters)) == strlen(optarg)) {
        if (strchr(optarg, 'A'))
          sections = all_section_letters; // keep the A
        else
          sections = optarg;
        // Record whether any sections require detecting symmetry
        if (strcspn(sections.c_str(), "s") != sections.size())
          detect_symmetry = true;
      }
      else
        error(msg_str("contains incorrect section type letter '%c'",
                      *(optarg + len)),
              c);
      break;
    }

    case 'C': {
      const char *all_count_letters = "AFwEDSsohO";
      size_t len;
      if ((len = strspn(optarg, all_count_letters)) == strlen(optarg)) {
        if (strchr(optarg, 'A'))
          counts = all_count_letters;
        else
          counts = optarg;
      }
      else
        error(msg_str("contains incorrect section type letter '%c'",
                      *(optarg + len)),
              c);
      // Record whether any counts require detecting symmetry
      if (strcspn(counts.c_str(), "O") != counts.size())
        detect_symmetry = true;
      break;
    }

    case 'k':
      orient = false;
      break;

    case 'y': {
      int n = 0;
      for (const char *p = optarg; *p; p++) {
        n += (*p == ',');
        if (n > 2)
          error("too many comma separated parts", c);
      }
      sub_sym = optarg;
      break;
    }

    case 'E':
      if (strlen(optarg) != 1 || !strchr("eia", *optarg))
        error("reporting edge type must be e, i or a");
      edge_type = *optarg;
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

  if (sections == "" && counts == "") {
    warning("no print options set, setting option -S G");
    sections = "G";
  }
}

void print_sections(rep_printer &rep, const char *sections)
{
  for (const char *c = sections; *c; c++) {
    switch (*c) {
    case 'A': // ignore
      break;
    case 'G':
      rep.general_sec();
      break;
    case 'F':
      rep.faces_sec();
      break;
    case 'E':
      rep.edges_sec();
      break;
    case 'a':
      rep.angles_sec();
      break;
    case 'S':
      rep.solid_angles_sec();
      break;
    case 'D':
      rep.distances_sec();
      break;
    case 's':
      rep.symmetry();
      break;
    }
  }
}

void print_counts(rep_printer &rep, const char *counts)
{
  for (const char *c = counts; *c; c++) {
    switch (*c) {
    case 'A': // ignore
      break;
    case 'F':
      rep.face_angles_cnts();
      break;
    case 'w':
      rep.windings();
      break;
    case 'E':
      rep.edge_lengths_cnts();
      break;
    case 'D':
      rep.dihedral_angles_cnts();
      break;
    case 'S':
      rep.solid_angles_cnts();
      break;
    case 's':
      rep.face_sides_cnts();
      break;
    case 'o':
      rep.vert_order_cnts();
      break;
    case 'h':
      rep.vert_heights_cnts();
      break;
    case 'O':
      rep.sym_orbit_cnts();
      break;
    }
  }
}

int main(int argc, char *argv[])
{
  or_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);

  if (opts.edge_type == 'a')
    geom.add_missing_impl_edges();
  else if (opts.edge_type == 'i') {
    geom.clear(EDGES);
    geom.add_missing_impl_edges();
  }

  if (opts.center_is_centroid)
    opts.center = geom.centroid();

  FILE *ofile = stdout; // write to stdout by default
  if (opts.ofile != "") {
    ofile = fopen(opts.ofile.c_str(), "w");
    if (ofile == nullptr)
      opts.error("could not open output file '" + opts.ofile + "'");
  }

  rep_printer rep(geom, ofile);
  rep.set_sig_dgts(opts.sig_digits);
  rep.set_center(opts.center);

  if (opts.detect_symmetry && !rep.set_sub_symmetry(opts.sub_sym))
    opts.error(("could not set subsymmetry: " + opts.sub_sym).c_str(), 'y');

  rep.is_oriented(); // set oriented value before orienting
  if (opts.orient) {
    geom.orient();
    if (GeometryInfo(geom).volume() < 0) // inefficient
      geom.orient_reverse();
  }

  print_sections(rep, opts.sections.c_str());
  print_counts(rep, opts.counts.c_str());

  if (opts.ofile == "")
    fclose(ofile);

  return 0;
}
