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
   Name: zono.cc
   Description: make zonohedra 
   Project: Antiprism - http://www.antiprism.com
*/



#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include <string>
#include <vector>

#include "../base/antiprism.h"

using std::string;
using std::vector;


class zo_opts: public prog_opts{
   public:
      char method;
      vec3d centre;
      bool centroid;
      bool out_star;
      bool unit_len;
      int pol_num;
      bool non_polar_opt;
      string ifile;
      string ofile;

      zo_opts(): prog_opts("zono"), method('v'),
                 centre(vec3d(0, 0, 0)), centroid(false),
                 out_star(false), unit_len(false),
                 pol_num(0), non_polar_opt(false) {}
      void process_command_line(int argc, char **argv);
      void usage();
};


void zo_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Make a polar zonohedron, or a zonohedron based on the vectors from an\n"
"OFF file. If input_file is not given the program reads from standard input\n"
"\n"
"Options\n"
"%s"
"  -m <mthd> method to create star from input, can be\n"
"               v - centre to vertices are vectors (default)\n"
"               a - all vertex to vertex are vectors\n"
"               i - implicit edges (face sides) are vectors\n"
"               e - explicit edges are vectors\n"
"  -c <cent> centre of points for method v, C for centroid (default: 0,0,0)\n"
"  -s        output the star (instead of the zonohedron)\n"
"  -u        make vectors unit length\n"
"  -P <num>  polar zonohedron with given number of vectors\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}

void zo_opts::process_command_line(int argc, char **argv)
{
   char errmsg[MSG_SZ];
   opterr = 0;
   char c;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hm:c:suP:o:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         case 'm':
            if(!(strlen(optarg)==1 && strchr("vaei", *optarg)))
               error("unknown method '"+string(optarg)+"'", c);
            method = *optarg;
            non_polar_opt = true;
            break;

         case 'c':
            if((strlen(optarg)==1 && strchr("Cc", *optarg)))
               centroid = true;
            else if(!centre.read(optarg))
               error("invalid centre '"+string(optarg)+"'", c);
            non_polar_opt = true;
            break;
 
         case 's':
            out_star = true;
            break;
         
         case 'u':
            unit_len = true;
            break;
         
         case 'P':
            if(!read_int(optarg, &pol_num, errmsg))
               error(errmsg, c);
            if(pol_num < 3)
               error("must be 3 or greater", c);
            break;
         
         case 'o':
            ofile = optarg;
            break;

         default:
            error("unknown command line error");
      }
   }

   if(pol_num && (non_polar_opt||argc-optind))
      warning("using -P, ignoring options m, a, c or an input file");
   
   if(argc-optind > 1)
      error("too many arguments");
   
   if(argc-optind == 1)
      ifile=argv[optind];

}




bool make_star_and_zono(geom_if &zono, zo_opts &opts, char *errmsg)
{
   zono.clear_all();

   col_geom_v geom;
   if(!geom.read(opts.ifile, errmsg))
      opts.error(errmsg);
      
   if(opts.centroid)
      opts.centre = centroid(*geom.get_verts());
   
   vector<vec3d> star = geom.get_star(opts.method, opts.centre);

   if(star.size()<3) {
      snprintf(errmsg, MSG_SZ, "star contains %lu vectors, needs at least 3",
            (unsigned long)star.size());
      return false;
   }

   if(opts.unit_len) {
      for(unsigned int i=0; i<star.size(); i++)
         star[i].to_unit();
   }
         
   if(opts.out_star) {
      zono.add_verts(star);
      return true;
   }
   else
      return zono.set_zono(star, errmsg);
}


int main(int argc, char **argv)
{
   zo_opts opts;
   opts.process_command_line(argc, argv);
   geom_v zono;

   char errmsg[MSG_SZ] = "";
   if(opts.pol_num) 
      make_polar_zono(zono, opts.pol_num, opts.out_star);
   else
      if(!make_star_and_zono(zono, opts, errmsg))
         opts.error(errmsg);
   
   if(!zono.write(opts.ofile, errmsg))
      opts.error(errmsg);

   return 0;
   
}


