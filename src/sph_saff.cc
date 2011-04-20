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
   Name: sph_saff.cc
   Description: 
   Project: Antiprism - http://www.antiprism.com
*/



#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "../base/antiprism.h"


class saff_opts: public prog_opts {
   public:
      int num_pts;
      double ang;
      bool pts_at_poles;
      
      string ofile;

      saff_opts(): prog_opts("sph_saff"),
                   ang(-1), pts_at_poles(false)
                   {}
      void process_command_line(int argc, char **argv);
      void usage();
};



void saff_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] num_points\n"
"\n"
"Distribute a number of points on a sphere using the algorithm\n"
"from \"Distributing many points on a sphere\" by E.B. Saff and\n"
"A.B.J. Kuijlaars, Mathematical Intelligencer 19.1 (1997) 5--11.\n"
"\n"
"Options\n"
"%s"
"  -a <ang>  increment each point placement by a fixed angle instead of\n"
"            using Saff and Kuiljaars' placement method\n"
"  -p        place points at the two poles\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}


void saff_opts::process_command_line(int argc, char **argv)
{
   char errmsg[MSG_SZ];
   opterr = 0;
   char c;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":ha:po:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         case 'a':
            if(!read_double(optarg, &ang, errmsg))
               error(errmsg, c);
            ang = fmod(ang, 360);
            if(ang<0)
               ang += 360;
            break;

         case 'p':
            pts_at_poles = true;
            break;

         case 'o':
            ofile = optarg;
            break;

         default:
            error("unknown command line error");
      }
   }

   if(argc-optind < 1)
      error("must specify the number of points");
   
   if(argc-optind > 1)
      error("too many arguments");
  
   if(!read_int(argv[optind], &num_pts, errmsg))
      error(errmsg, "num_pts");
   if(num_pts < 1)
      error("number of points must be a positive integer");
}

void make_saff(geom_if &geom, int num, double ang, bool pts_at_poles)
{
   if(!pts_at_poles)
      num += 2;
   double phi = 0;
   double theta = 0;
   for(int i=0; i<num; i++) {
      //double h = -1 + 2.0*(i-1)/(num-1);
      double h = -1 + (2.0*i)/(num-1);
      if(i==0 || i==num-1) {
         if(pts_at_poles) {
            phi = 0;
            theta = acos(-1+2*(h>0));
         }
         else
            continue;
      }
      else {
         theta = acos(h);
         if(ang<0)
            phi += 3.6/sqrt(num*(1-h*h));
         else
            phi += ang;
      }
      phi = fmod(phi, 2*M_PI);

      geom.get_verts()->push_back(
            vec3d(sin(phi)*sin(theta), cos(theta), cos(phi)*sin(theta) ) );
   }
}

   
int main(int argc, char **argv)
{
   saff_opts opts;
   opts.process_command_line(argc, argv);

   col_geom_v geom;
   make_saff(geom, opts.num_pts, deg2rad(opts.ang), opts.pts_at_poles);

   char errmsg[MSG_SZ];
   if(!geom.write(opts.ofile, errmsg))
      opts.error(errmsg);

   return 0;
}


  

