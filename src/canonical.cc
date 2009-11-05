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
   Name: canonical.cc
   Description: canonicalize a polyhedron
                Uses George Hart's two canonicalization algorithm
                http://library.wolfram.com/infocenter/Articles/2012/
                http://www.georgehart.com/virtual-polyhedra/conway_notation.html
   Project: Antiprism - http://www.antiprism.com
*/

#include <string.h>
#include <unistd.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <limits.h>
#include <string>
#include <vector>

#include "../base/antiprism.h"

using std::string;
using std::vector;


class cn_opts: public prog_opts
{
   public:
      int num_iters;
      double edge_factor;
      double plane_factor;
      int lim_exp_canonical;
      char method;
      char cent_type;
      int num_iters_preplanar;
      int lim_exp_preplanar;
      int rep_count;
      int divergence_test;
      
      string ifile;
      string ofile;
      
      int epsilon_num;

      cn_opts(): prog_opts("canonical"),
                 num_iters(-1),
                 edge_factor(50),
                 plane_factor(20),
                 lim_exp_canonical(INT_MAX),
                 method('m'),
                 cent_type('\0'),
                 num_iters_preplanar(-1),
                 lim_exp_preplanar(INT_MAX),
                 rep_count(50),
                 divergence_test(10)
                 {}

      void process_command_line(int argc, char **argv);
      void usage();
};

void cn_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Read a polyhedron from a file in OFF format. Canonicalize or planarize it.\n"
"Uses algorithms by George W. Hart, http://www.georgehart.com/\n"
"http://www.georgehart.com/virtual-polyhedra/conway_notation.html\n"
"http://www.georgehart.com/virtual-polyhedra/canonical.html\n"
"If input_file is not given the program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"  -n <itrs> maximum number of iterations (default: no limit)\n"
"  -l <lim>  minimum distance change to terminate, as negative exponent\n"
"               (default: %d giving %.0e)\n"
"  -d <int>  divergence test. 0 for no test. (default 10)\n"
"  -M <mthd> canonicalizing method,\n"
"            m - mathematica version of canonicalization (default)\n"
"            n - conway notation version of canonicalization\n"
"            l - mathematica planarize portion only\n"
"            p - conway notation planarize (face centroids reciprocal)\n"
"            q - conway notation planarize (face centroids magnitude reciprocal)\n"
"            x - face centroids only (no reciprocal) planarize method\n"
"  -C <cent> initial 'centering'\n"
"            x - none, c - centroid (-M p and -M l default)\n"
"            s - centroid and project vertices onto a sphere (-M m default)\n"
"            p - centroid and pre-planarized (-M n default)\n"
"            q - centroid and pre-planarized with magnitude reciprocal\n"
"  -z <n>    status reporting every n lines. -1 for no status. (default 50)\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"Mathematica Canonicalize Options (-M m and -M l)\n"
"  -e <perc> percentage to scale the edge tangency error (default: 50)\n" 
"  -p <perc> percentage to scale the face planarity error (default: 20)\n" 
"\n"
"Pre-planarization Options (-C p and -C q)\n"
"  -i <itrs> maximum number of pre-planarize iterations (default: no limit)\n"
"  -j <lim>  minimum distance change to terminate pre-planarize, as negative\n"
"               exponent (default: %d giving %.0e)\n"
"\n"
"\n",prog_name(), help_ver_text, epsilon_num,::epsilon,epsilon_num,::epsilon);
}

void cn_opts::process_command_line(int argc, char **argv)
{
   char errmsg[MSG_SZ];
   opterr = 0;
   char c;
   
   epsilon_num = int(-log(::epsilon)/log(10) + 0.5);
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hn:e:p:l:d:M:C:i:j:z:o:")) != -1 ) {
      if(common_opts(c))
         continue;

      switch(c) {
         case 'o':
            ofile = optarg;
            break;

         case 'n':
            if(!read_int(optarg, &num_iters, errmsg))
               error(errmsg, c);
            if(num_iters < 0)
               error("number of iterations must be 0 or greater", c);
            break;

         case 'd':
            if(!read_int(optarg, &divergence_test, errmsg))
               error(errmsg, c);
            if(divergence_test < 0)
               error("divergence test must be 0 or greater", c);
            break;

         case 'l':
            if(!read_int(optarg, &lim_exp_canonical, errmsg))
               error(errmsg, c);
            if(lim_exp_canonical < 0) {
               warning("limit is negative, and so ignored", c);
            }
            if(lim_exp_canonical > 16) {
               warning("limit is very small, may not be attainable", c);
            }
            break;

         case 'e':
            if(!read_double(optarg, &edge_factor, errmsg))
               error(errmsg, c);
            if(edge_factor <= 0 || edge_factor >=100)
               warning("not inside range 0 to 100", c);
            break;

         case 'p':
            if(!read_double(optarg, &plane_factor, errmsg))
               error(errmsg, c);
            if(plane_factor <= 0 || plane_factor >=100) {
               warning("not inside range 0 to 100", c);
            }
            break;

         case 'M':
            if(strlen(optarg)==1 && strchr("mnlpqx", int(*optarg)))
               method = *optarg;
            else
               error("method type must be m, n, l, p, q or x\n", c);
            break;

          case 'C':
            if(strlen(optarg)!=1)
               error("initial centering must be exactly one character", c);
            if(!strchr("xcspq", *optarg))
               error("initial centering must be x, c, s, p or q\n", c);
            cent_type = *optarg;
            break;

         case 'i':
            if(!read_int(optarg, &num_iters_preplanar, errmsg))
               error(errmsg, c);
            if(num_iters_preplanar <= 0)
               error("number of iterations for preplanarization must be greater than 0", c);
            break;

         case 'j':
            if(!read_int(optarg, &lim_exp_preplanar, errmsg))
               error(errmsg, c);
            if(lim_exp_preplanar < 0) {
               warning("limit is negative, and so ignored", c);
            }
            if(lim_exp_preplanar > 16) {
               warning("limit is very small, may not be attainable", c);
            }
            break;

         case 'z':
            if(!read_int(optarg, &rep_count, errmsg))
               error(errmsg, c);
            if(rep_count < -1)
               error("number of iterations must be -1 or greater", c);
            break;

         default:
            error("unknown command line error");
      }
   }

   if(!cent_type) {
      if(method=='c')
         cent_type = 's';
      else
      if(method=='n')
         cent_type = 'p';
      else   
      if(method=='l' || method=='p' || method=='q')
         cent_type = 'c';
      else
         cent_type = 'x';
   }
   
   if(argc-optind > 1)
      error("too many arguments");
   
   if(argc-optind == 1)
      ifile=argv[optind];

   lim_exp_preplanar = (lim_exp_preplanar != INT_MAX) ? lim_exp_preplanar : epsilon_num;
   lim_exp_canonical = (lim_exp_canonical != INT_MAX) ? lim_exp_canonical : epsilon_num;
}

void centroid_to_origin(geom_if &geom)
{
   geom.transform(mat3d::transl(-centroid(geom.verts())));
}

void project_onto_sphere(geom_if &geom)
{
   vector<vec3d> &verts = geom.raw_verts();
   for(unsigned int i=0; i<verts.size(); i++)
      verts[i].to_unit();
}

int main(int argc, char *argv[])
{
   cn_opts opts;
   opts.process_command_line(argc, argv);

   char errmsg[MSG_SZ] = "";
   col_geom_v geom;
   if(!geom.read(opts.ifile, errmsg))
      opts.error(errmsg);
   if(*errmsg)
      opts.warning(errmsg);

   if(opts.cent_type=='c' || opts.cent_type=='p' || opts.cent_type=='q' || opts.cent_type=='s')
      centroid_to_origin(geom);

   if(opts.cent_type=='p' || opts.cent_type=='q')
      canonicalize_cn(geom, opts.num_iters_preplanar, pow(10, -opts.lim_exp_preplanar),
                      opts.cent_type, opts.divergence_test, opts.rep_count);
   else
   if(opts.cent_type=='s')
      project_onto_sphere(geom);

   if(opts.method=='m' || opts.method=='l') {
      bool planarize_only = (opts.method=='l') ? true : false;
      canonicalize_mm(geom, opts.edge_factor/100, opts.plane_factor/100, opts.num_iters, pow(10, -opts.lim_exp_canonical),
                       opts.divergence_test, opts.rep_count, planarize_only);
   }
   else {
      // conway notation canonicalize expects 'c' for canonicalize
      canonicalize_cn(geom, opts.num_iters, pow(10, -opts.lim_exp_canonical), opts.method, opts.divergence_test, opts.rep_count);
   }

   if(!geom.write(opts.ofile, errmsg))
      opts.error(errmsg);

   return 0;
}
