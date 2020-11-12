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

/* \file off2pov.cc
   \brief Convert an OFF file to POV format
*/

#include "../base/antiprism.h"

#include <cctype>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>

using std::string;
using std::vector;

using namespace anti;

class o2p_opts : public ViewOpts {
public:
  bool shadow;
  int stereo_type;
  int sig_dgts;
  char o_type;
  vector<string> scene_incs;
  vector<string> obj_incs;
  vector<string> geom_incs;

  string ofile;

  o2p_opts()
      : ViewOpts("off2pov"), shadow(false), stereo_type(-1),
        sig_dgts(DEF_SIG_DGTS), o_type('a')
  {
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

void o2p_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] input_files

Convert files in OFF format to POV format for display in POV-ray. If
input_files are not given the program reads from standard input.

Options
%s
%s
  -O <type> output type, can be: 'a' all in one POV file (default),
            's' separate files, 'o' objects only, 't' template only
  -i <fils> include files (separated by commas) for every POV geometry
  -j <fils> include files (separated by commas) for the POV scene file
  -J <fils> include files (separated by commas) containing additional POV
            objects for the POV scene file
  -o <file> write output to file (default: write to standard output)

  Scene options
%s
  -P <pers> narrow the angle of perspective (range 0-100,
            default: 2, recommend 4 for stereo option -S 1)
  -W        use lighting with shado
  -S <type> produce stereo output, type is 0 (default) mono, 1 stereo
            with one image file, 2 stereo with two image files (use
            the POV-Ray +KFF2 option for output)

  Precision options
%s

)",
          prog_name(), help_ver_text, help_view_text, help_scene_text,
          help_prec_text);
}

void o2p_opts::process_command_line(int argc, char **argv)
{
  Status stat;
  opterr = 0;
  int c;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv,
                     ":hv:e:V:E:F:m:x:s:n:o:D:C:L:R:P:W:S:B:d:t:I:j:J:i:O:")) !=
         -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'o':
      ofile = optarg;
      break;

    case 'W':
      shadow = true;
      break;

    case 'S':
      print_status_or_exit(read_int(optarg, &stereo_type), c);
      if (stereo_type < 0 || stereo_type > 3)
        error("stereo type must be 0, 1, 2 or 3", c);
      break;

    case 'd':
      print_status_or_exit(read_int(optarg, &sig_dgts), c);
      break;

    case 'j':
      scene_incs.push_back(optarg);
      break;

    case 'J':
      obj_incs.push_back(optarg);
      break;

    case 'i':
      geom_incs.push_back(optarg);
      break;

    case 'O':
      if (strlen(optarg) != 1)
        error("output type must be exactly one character", c);
      if (!strchr("asot", *optarg))
        error("output type must be 'a' all in one, 's' separate files, 'o' "
              "objects only 't' template only",
              c);
      o_type = *optarg;
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
  else if (o_type != 't')
    ifiles.push_back("");
}

void set_geom_includes(Scene &scen, const vector<string> &includes)
{
  vector<SceneGeometry> &sc_geoms = scen.get_geoms();
  vector<SceneGeometry>::iterator sc_i;
  for (sc_i = sc_geoms.begin(); sc_i != sc_geoms.end(); ++sc_i) {
    vector<GeometryDisplay *> &disps = sc_i->get_disps();
    vector<GeometryDisplay *>::iterator disp_i;
    for (disp_i = disps.begin(); disp_i != disps.end(); ++disp_i) {
      DisplayPoly *disp = dynamic_cast<DisplayPoly *>(*disp_i);
      if (disp)
        disp->set_includes(includes);
    }
  }
}

int main(int argc, char *argv[])
{
  o2p_opts opts;
  opts.process_command_line(argc, argv);
  Scene scen = opts.scen_defs;
  opts.set_view_vals(scen);
  set_geom_includes(scen, opts.geom_incs);

  PovWriter pov;
  pov.set_o_type(opts.o_type);

  if (opts.ofile != "")
    pov.set_file_name(basename2(opts.ofile.c_str()));
  else
    pov.set_file_name("stdout");

  pov.set_includes(opts.scene_incs);
  pov.set_obj_includes(opts.obj_incs);

  if (opts.stereo_type >= 0)
    pov.set_stereo_type(opts.stereo_type);

  if (opts.shadow)
    pov.set_shadow(opts.shadow);

  FILE *ofile = stdout; // write to stdout by default
  if (opts.ofile != "") {
    ofile = fopen(opts.ofile.c_str(), "w");
    if (ofile == nullptr)
      opts.error("could not open output file \'" + opts.ofile + "\'");
  }

  pov.write(ofile, scen, opts.sig_dgts);

  if (opts.ofile != "")
    fclose(ofile);

  return 0;
}
