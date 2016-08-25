/*
   Copyright (c) 2003-2009, Adrian Rossiter

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

#include <string.h>
#include <math.h>
#include <ctype.h>
#include <string>
#include <vector>
#include "../base/antiprism.h"

using std::string;
using std::vector;


class o2p_opts: public view_opts {
   public:
      bool shadow;
      int stereo_type;
      int sig_dgts;
      char o_type;
      vector<string> scene_incs;
      vector<string> obj_incs;
      vector<string> geom_incs;
      
      string ofile;

      o2p_opts(): view_opts("off2pov"),
                  shadow(false),
                  stereo_type(-1),
                  sig_dgts(DEF_SIG_DGTS),
                  o_type('a')
                  {}

      void process_command_line(int argc, char **argv);
      void usage();
};


void o2p_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] input_files\n"
"\n"
"Convert files in OFF format to POV format for display in POV-ray. If\n"
"input_files are not given the program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"%s"
"  -O <type> output type, can be: 'a' all in one POV file (default),\n"
"            's' separate files, 'o' objects only, 't' template only\n"
"  -i <fils> include files (separated by commas) for every POV geometry\n"
"  -j <fils> include files (separated by commas) for the POV scene file\n"
"  -J <fils> include files (separated by commas) containing additional POV\n"
"            objects for the POV scene file\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"  Scene options\n"
"%s"
"  -P <pers> narrow the angle of perspective (range 0-100,\n"
"            default: 2, recommend 4 for stereo option -S 1)\n"
"  -W        use lighting with shadows"
"  -S <type> produce stereo output, type is 0 (default) mono, 1 stereo\n"
"            with one image file, 2 stereo with two image files (use\n"
"            the POV-Ray +KFF2 option for output)\n"
"\n"
"  Precision options\n"
"%s"
"\n"
"\n", prog_name(), help_ver_text, help_view_text,
      help_scene_text, help_prec_text);
}


void o2p_opts::process_command_line(int argc, char **argv)
{
   char errmsg[MSG_SZ];
   vector<string> warnings;
   opterr = 0;
   int c;
  
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hv:e:V:E:F:m:x:s:n:o:D:C:L:R:P:W:S:B:d:t:I:j:J:i:O:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         case 'o':
            ofile = optarg;
            break;

         case 'W':
            shadow = true;
            break;

         case 'S':
            if(!read_int(optarg, &stereo_type, errmsg))
               error(errmsg, c);
            if(stereo_type<0 || stereo_type>3)
               error("stereo type must be 0, 1, 2 or 3", c);
            break;
         
         case 'd':
            if(!read_int(optarg, &sig_dgts, errmsg))
               error(errmsg, c);
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
            if(strlen(optarg)!=1)
               error("output type must be exactly one character", c);
            if(!strchr("asot", *optarg))
               error("output type must be 'a' all in one, 's' separate files, 'o' objects only 't' template only", c);
            o_type = *optarg;
            break;
         
         default:
            if(read_disp_option(c, optarg, errmsg, warnings)) {
               if(*errmsg)
                  error(errmsg, c);
               for(unsigned int i=0; i<warnings.size(); i++)
                  warning(warnings[i], c);
            }
            else 
               error("unknown command line error");
      }
   }

   if(argc-optind >= 1)
      while(argc-optind >= 1)
         ifiles.push_back(argv[optind++]);
   else if(o_type!='t')
      ifiles.push_back("");
}

void set_geom_includes(scene &scen, const vector<string> &includes)
{
   vector<scene_geom> &sc_geoms = scen.get_geoms();
   vector<scene_geom>::iterator sc_i;
   for(sc_i=sc_geoms.begin(); sc_i!=sc_geoms.end(); ++sc_i) {
      vector<geom_disp *> &disps = sc_i->get_disps();
      vector<geom_disp *>::iterator disp_i;
      for(disp_i=disps.begin(); disp_i!=disps.end(); ++disp_i) {
         disp_poly *disp = dynamic_cast<disp_poly *>(*disp_i);
         if(disp)
            disp->set_includes(includes);
      }
   }
}
 

int main(int argc, char *argv[])
{
   o2p_opts opts;
   opts.process_command_line(argc, argv);
   scene scen = opts.scen_defs;
   opts.set_view_vals(scen);
   set_geom_includes(scen, opts.geom_incs);
 
   pov_writer pov;
   pov.set_o_type(opts.o_type);
  
   if(opts.ofile!="")
      pov.set_file_name(basename2(opts.ofile.c_str()));
   else
      pov.set_file_name("stdout");
   
   pov.set_includes(opts.scene_incs);
   pov.set_obj_includes(opts.obj_incs);
  
   if(opts.stereo_type >= 0)
      pov.set_stereo_type(opts.stereo_type);
 
   if(opts.shadow)
      pov.set_shadow(opts.shadow);

   FILE *ofile = stdout;  // write to stdout by default
   if(opts.ofile != "") {
      ofile = fopen(opts.ofile.c_str(), "w");
      if(ofile == 0)
         opts.error("could not open output file \'"+opts.ofile+"\'");
   }

   pov.write(ofile, scen, opts.sig_dgts);
   
   if(opts.ofile!="")
      fclose(ofile);

   return 0;
}
   

