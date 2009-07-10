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
   Name: kcycle.cc
   Description: make a kaliedocycle with a polyhedron
   Project: Antiprism - http://www.antiprism.com
*/

#include <string.h>
#include <unistd.h>
#include <math.h>
#include <ctype.h>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include "../base/antiprism.h"

using std::string;
using std::vector;

bool read_ints(char *ints_str, vector<int> &ints, char *errmsg)
{
   int i_idx;
   char *i_str = strtok(ints_str, ",");
   int i=0;
   while(i_str) {
      i++;
      if(!read_int(i_str, &i_idx, errmsg) || i_idx < 0) {
         snprintf(errmsg, MSG_SZ, "index no. %d, \"%s\" is not a positive integer", i, i_str);
         return false;
      }
      ints.push_back(i_idx);
      i_str = strtok(NULL, ",");
   }
     
   return true;
}




class kc_opts: public prog_opts
{
   private:
   
   public:
      vector<int> edge_idxs;
      double angle;
      int num_prs;

      string ifile;
      string ofile;

      kc_opts(): prog_opts("kcycle"), angle(0.0), num_prs(6) {}

      void process_command_line(int argc, char **argv);
      void usage();
};

void kc_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] hinge_vertex_idxs\n"
"\n"
"Make a kaleidocyle from a polyhedron. hinge_vertex_idxs is four integers\n"
"separated by commas (e.g. 0,1,2,3). These are the index numbers of the\n"
"vertices for the two eges that will be hinges in the cycle.\n"
"\n"
"Options\n"
"%s"
"  -a <ang>  angle in degrees to rotate the first hinge from\n"
"            horizontal (default: 0.0)\n"
"  -n <num>  number of pairs of polyhedra in cycle (default: 6)\n"
"  -i <file> input file in OFF format. If '-' then read file from stdin.\n"
"            (default: a tetrahedron with hinge_vertex_idxs 0,1,2,3)\n" 
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}


     
void kc_opts::process_command_line(int argc, char **argv)
{
   char errmsg[MSG_SZ];
   extern char *optarg;
   extern int optind, opterr;
   opterr = 0;
   char c;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":ha:n:i:o:")) != -1 ) {
      if(common_opts(c))
         continue;

      switch(c) {
         case 'a':
            if(!read_double(optarg, &angle, errmsg))
               error(errmsg, c);
            angle = deg2rad(angle);
            break;

         case 'n':
            if(!read_int(optarg, &num_prs, errmsg))
               error(errmsg, c);
            if(num_prs < 1)
               error("must be at least two pairs", c);
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

   if(argc-optind > 1) {
      error("too many arguments (leave no spaces between the edge indexes)");
      exit(1);
   }
   
   if(argc-optind < 1) {
      if(ifile == "") {
         for(int i=0; i<4; i++)
            edge_idxs.push_back(i);
      }
      else
         error("edge index numbers not given");
   }
   
   if(argc-optind==1 && !read_ints(argv[optind], edge_idxs, errmsg))
         error(errmsg, "edge index numbers");


}


void kcycle(col_geom_v geom, col_geom_v &cycle, int num_prs, vector<int> idxs, double angle)
{
   
   vector<vec3d> &verts = *geom.get_verts();
   vec3d P[] = {verts[idxs[0]], verts[idxs[1]], verts[idxs[2]], verts[idxs[3]]};
   vec3d e1(cos(angle), sin(angle), 0);  // direction of first edge
   
   mat3d trans1 = mat3d::rot(P[1]-P[0], e1)
                * mat3d::transl(-P[0]);
   for(int i=0; i<4; i++)
      P[i] = trans1 * P[i];
   
   vec3d e2=(P[3]-P[2]).unit();          // current direction of second edge
   
   double wedge_angle = -M_PI/num_prs;
   vec3d n2(sin(wedge_angle), 0, cos(wedge_angle)); // normal to other slice

   // Use spherical trig to find the angle to rotate e2 about e1 until
   // it is normal to n2 
   int sign1 = vtriple(e1, e2, n2) > 0;
   double ang1 =  acos( (vdot(e2, n2) - vdot(e1, n2)*vdot(e1, e2) ) /
                        (vcross(e1, n2).mag()*vcross(e1, e2).mag())   );
   int sign2 = vdot(e1, e2) > 0;
   double ang2 =  acos( (- vdot(e1, n2)*vdot(e1, e2) ) /
                        (vcross(e1, n2).mag()*vcross(e1, e2).mag())   );
   double ang =  (2*sign1 - 1)*ang1 - (2*sign2 - 1)*ang2;
   
   mat3d trans2 = mat3d::rot(e1, ang);
   for(int i=0; i<4; i++)
      P[i] = trans2 * P[i];
   
   vec3d N0, N1;
   lines_nearest_points(P[0], P[1], P[2], P[3], N0, N1);
   vec3d perp = N1-N0;
   double dist = -perp[0] - perp[2]/tan(wedge_angle);
   
   mat3d trans3 = mat3d::transl(vec3d(dist,0,0) - N0) * trans2 * trans1;
   
   // align for symmetry exp
   trans3 = mat3d::rot(vec3d(0,1,0), vec3d(0,0,1)) *
            mat3d::rot(vec3d(0,1,0), M_PI/num_prs) *
            trans3;
   for(unsigned int i=0; i<verts.size(); i++)
      verts[i] = trans3 * verts[i];

   // align centroid with z=0
   vec3d cent = centroid(verts);
   mat3d trans_c = mat3d::transl(vec3d(0,0,-cent[2]));
   for(unsigned int i=0; i<verts.size(); i++)
      verts[i] = trans_c * verts[i];
   
   // make the cycle
   sym_repeat(cycle, geom, sch_sym(sch_sym::Cv, num_prs));
}


int main(int argc, char *argv[])
{
   kc_opts opts;
   opts.process_command_line(argc, argv);

   char errmsg[MSG_SZ];
   col_geom_v geom;
   if(opts.ifile != "") {
      if(!geom.read(opts.ifile, errmsg))
         opts.error(errmsg);
      if(*errmsg)
         opts.warning(errmsg);
   }
   else
      tetrahedron(geom);

   col_geom_v cycle;
   kcycle(geom, cycle, opts.num_prs, opts.edge_idxs, opts.angle);

   if(!cycle.write(opts.ofile, errmsg))
      opts.error(errmsg);

   return 0;
}



