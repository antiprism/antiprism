/*
   Copyright (c) 2003-2009, Adrian Rossiter

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



#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>

#include "../base/antiprism.h"


class geo_opts: public prog_opts {
   public:
      vec3d centre;
      double radius;
      int m;
      int n;
      int pat_freq;
      bool trad_freq;
      char method;
      bool keep_flat;
      bool equal_len_div;
      string ifile;
      string ofile;

      geo_opts(): prog_opts("geodesic"),
                  centre(vec3d(0,0,0)),
                  m(1), n(0),
                  pat_freq(1), trad_freq(false),
                  method('s')
                  {}
      void process_command_line(int argc, char **argv);
      void usage();
};



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
"  -f <freq> A positive integer (default: 1), the number of repeats of\n"
"            the specified pattern along an edge\n"
"  -F <freq> A positive integer. As -f, but if pattern resolves to Class II\n"
"            then freq must be even and the repeat is taken as -f freq/2\n"
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


void geo_opts::process_command_line(int argc, char **argv)
{
   char errmsg[MSG_SZ];
   opterr = 0;
   char c;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hf:F:c:M:C:o:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         case 'f':
         case 'F':
            if(!read_int(optarg, &pat_freq, errmsg) || pat_freq<1) {
               snprintf(errmsg, MSG_SZ, "frequency is '%s', should be a positive integer", optarg);
               error(errmsg, c);
            }
            trad_freq = (c=='F');
            break;

         case 'c':
            char *pos;
            if(!(pos = strchr(optarg, ','))) {
               if( !(strcmp(optarg, "1")==0 || strcmp(optarg, "2")==0) )
                  error("class type must be 1 or 2, or two integers separated by a comma\n", c);
               m = 1;
               n = *optarg=='2' ? 1 : 0;
            }
            else {
               *pos = 0;
               char *pat[] = {optarg, pos+1};
               int nums[2];
               for(int i=0; i<2; i++) {
                  if(!read_int(pat[i], &nums[i], errmsg) || nums[i] < 0) {
                     snprintf(errmsg, MSG_SZ, "%s value is '%s' should be a positive integer or zero", (i==0) ? "first" : "second", pat[i]);
                     error(errmsg, c);
                  }
               }
               if(nums[0]==0 && nums[1]==0)
                  error("pattern must include at least one positive integer",c);

               m = nums[0];
               n = nums[1];
            }
            break;

         case 'M':
            if(strlen(optarg)==1 && strchr("sp", int(*optarg)))
               method = *optarg;
            else
               error("method type must be s or p\n", c);
            break;
         
         case 'C':
            if(!centre.read(optarg)) {
               snprintf(errmsg, MSG_SZ, "centre is '%s', must be three numbers separated by commas", optarg);
               error(errmsg, c);
            }
            break;

         case 'o':
            ofile = optarg;
            break;

         default:
            error("unknown command line error");
      }
   }

   if(argc-optind > 1)
      error("too many arguments");
   
   if(argc-optind == 1)
      ifile = argv[optind];
   

   if(m==n && trad_freq) {
      if(!is_even(pat_freq)) {
         snprintf(errmsg, MSG_SZ, "frequency must be even with basic Class II pattern\n");
         error(errmsg, "F");
      }
      pat_freq /= 2;
   }

   m *= pat_freq;
   n *= pat_freq;
   
   return;
}


int main(int argc, char **argv)
{
   geo_opts opts;
   opts.process_command_line(argc, argv);

   col_geom_v geom;
   char errmsg[MSG_SZ];
   if(!geom.read(opts.ifile, errmsg))
      opts.error(errmsg);
   if(*errmsg)
      opts.warning(errmsg);

   col_geom_v geo;
   if(opts.method=='s')
      geo.set_geodesic_sphere(geom, opts.m, opts.n, opts.centre);
   else if(opts.method=='p')
      geo.set_geodesic_planar(geom, opts.m, opts.n);
   
   if(!geo.write(opts.ofile, errmsg))
      opts.error(errmsg);

   return 0;
}


  

