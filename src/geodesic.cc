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
      int clas;
      int m;
      int n;
      char method;
      int pat_freq;
      bool keep_flat;
      bool equal_len_div;
      char poly;
      string ifile;
      string ofile;

      geo_opts(): prog_opts("geodesic"),
                  centre(vec3d(0,0,0)),
                  clas(1), m(1), n(0),
                  method('s'),
                  poly('x') {}
      void process_command_line(int argc, char **argv);
      void usage();
};



void geo_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [freq]\n"
"\n"
"Make higher frequency, plane-faced polyhedra or geodesic spheres.\n"
"For Class I and II patterns freq (default: 1, or 2 if -d is 2) is the\n"
"number of sections along an edge, for Class III patterns (and those\n"
"specified by two numbers) freq is the number of times the pattern is\n"
"repeated along an edge. Default is a frequency 1 icosahedral geodesic\n"
"sphere\n"
"\n"
"Options\n"
"%s"
"  -p <poly> type of poly: i - icosahedron (default), o - octahedron,\n"
"            t - tetrahedron, T - triangle\n"
"  -M <mthd> Method of applying the frequency pattern:\n"
"            s - geodesic sphere (default). The pattern grid is formed\n"
"                from divisions along each edge that make an equal angle\n"
"                at the centre. The geodesic vertices are centred at the\n"
"                origin and projected on to a unit sphere.\n"
"            p - planar. The pattern grid is formed from equal length\n"
"                divisions along each edge, the new vertices lie on the\n"
"                surface of the original polyhedron.\n"
"  -c <clss> face division pattern,  1 (Class I, default) or 2 (Class II),\n"
"            or two numbers separated by a comma to determine the pattern\n"
"            (Class III generally, but 1,0 is Class I, 1,1 is Class II, etc)\n"
"  -C <cent> centre of points, in form \"x_val,y_val,z_val\" (default: 0,0,0)\n"
"            used for geodesic spheres\n"
"  -i <file> input file in OFF format containing the base triangle-faced\n"
"            polyhedron to be divided. If '-' then read file from stdin\n" 
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

   while( (c = getopt(argc, argv, ":hp:c:M:C:i:o:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         case 'p':
            if(strlen(optarg)==1 && strchr("toiT", int(*optarg)))
               poly = *optarg;
            else
               error("polyhedron type must be t, o, i or T\n", c);
            break;

         case 'M':
            if(strlen(optarg)==1 && strchr("sp", int(*optarg)))
               method = *optarg;
            else
               error("method type must be s or p\n", c);
            break;
         
         case 'c':
            char *pos;
            if(!(pos = strchr(optarg, ','))) {
               if( !(strcmp(optarg, "1")==0 || strcmp(optarg, "2")==0) )
                  error("class type must be 1 or 2, or two integers separated by a comma\n", c);
               clas = atoi(optarg);
               m = 1;
               n = *optarg=='2' ? 1 : 0;
            }
            else {
               clas = 3;
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

         case 'C':
            if(!centre.read(optarg)) {
               snprintf(errmsg, MSG_SZ, "centre is '%s', must be three numbers separated by commas", optarg);
               error(errmsg, c);
            }
            break;

         case 'i':
            ifile = optarg;
            break;
            
         case 'o':
            ofile = optarg;
            break;

         default:
            error("unknown command line error");
      }
   }

   if(ifile !="" && poly!='x') {
      snprintf(errmsg, MSG_SZ, "option ignored, reading poly from '%s'", (ifile!="-")?ifile.c_str():"stdin");
      warning(errmsg, 'p');
   }

   if(ifile =="" && poly=='x')
      poly='i';

   if(argc-optind > 1)
      error("too many arguments");
   
   pat_freq = 1;
   if(argc-optind == 1) {
      if(!read_int(argv[optind], &pat_freq, errmsg)) {
         snprintf(errmsg, MSG_SZ, "value is '%s', should be a positive integer", argv[optind]);
         error(errmsg, "frequency");
      }

      if(clas==2) {
         if(!is_even(pat_freq)) {
            snprintf(errmsg, MSG_SZ, "value is '%d', for Class II divisions it must be an even integer\n", pat_freq);
            error(errmsg, "frequency");
         }
         pat_freq /= 2;
      }
   }

   m *= pat_freq;
   n *= pat_freq;
   
   return;
}


col_geom_v get_poly(char p_type)
{
   col_geom_v geom;
   vector<vec3d> &verts = *geom.get_verts();
   if(p_type=='i')
      icosahedron(geom);
   else if(p_type=='o')
      octahedron(geom);
   else if(p_type=='t')
      tetrahedron(geom);
   else { //if(p_type=='T') {
      double X = 0.25;
      double Y = sqrt(3.0)/12;
      double Z = 0.8;
      verts.push_back(vec3d(-X, -Y,  Z)); // 0
      verts.push_back(vec3d( X, -Y,  Z)); // 1
      verts.push_back(vec3d( 0,2*Y,  Z)); // 3

      vector<vector<int> > &faces = *geom.get_faces();
      int f[] = { 0, 1, 2,   0, 2, 1};
      for(int i=0; i<1; i++) {
         faces.push_back(vector<int>());
         for(int j=0; j<3; j++)
            faces.back().push_back(f[i*3+j]);
      }
   }
   
   geom.transform(mat3d::scale(1/verts[0].mag()));
   return geom;
}


int main(int argc, char **argv)
{
   geo_opts opts;
   opts.process_command_line(argc, argv);

   col_geom_v geom;
   char errmsg[MSG_SZ];
   if(opts.ifile != "") {
      if(!geom.read(opts.ifile, errmsg))
         opts.error(errmsg);
      if(*errmsg)
         opts.warning(errmsg);
   }
   else {
      geom = get_poly(opts.poly);
   }

   col_geom_v geo;
   if(opts.method=='s')
      geo.set_geodesic_sphere(geom, opts.m, opts.n, opts.centre);
   else if(opts.method=='p')
      geo.set_geodesic_planar(geom, opts.m, opts.n);
   
   if(!geo.write(opts.ofile, errmsg))
      opts.error(errmsg);

   return 0;
}


  

