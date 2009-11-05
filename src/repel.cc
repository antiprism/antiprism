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
   Name: repel.cc
   Description: points repelling (on a sphere)
   Project: Antiprism - http://www.antiprism.com
*/

#include <string.h>
#include <unistd.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <string>
#include <vector>

#include "../base/antiprism.h"

using std::string;
using std::vector;



class rep_opts: public prog_opts
{
   public:
      int num_iters;
      int num_pts;
      int rep_form;
      double shorten_by;
      int lim_exp;
      
      string ifile;
      string ofile;

      rep_opts(): prog_opts("repel"), num_iters(-1), num_pts(-1),
                  rep_form(2), shorten_by(-1), lim_exp(13)
                 {}

      void process_command_line(int argc, char **argv);
      void usage();
};


void rep_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"An equilibrium position is found for a set of points which repel each\n"
"other. The initial coordinates are read from input_file if given (or\n"
"from standard input), otherwise use -N to generate a random set.\n"
"\n"
"Options\n"
"%s"
"  -N <num>  initialise with a number of randomly placed points\n"
"  -n <itrs> maximum number of iterations (default: no limit)\n" 
"  -s <perc> percentage to shorten the travel distance (default: adaptive)\n" 
"  -l <lim>  minimum distance change to terminate, as negative\n"
"            exponent 1e-lim (default: 13 giving 1e-13)\n" 
"  -r <rep>  repelling formula\n"
"              1 - inverse of distance\n"
"              2 - inverse square of distance (default)\n"
"              3 - inverse cube of distance\n"
"              4 - inverse square root of distance\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}


     
void rep_opts::process_command_line(int argc, char **argv)
{
   char errmsg[MSG_SZ];
   opterr = 0;
   char c;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hn:N:s:l:r:o:")) != -1 ) {
      if(common_opts(c))
         continue;

      switch(c) {
         case 'n':
            if(!read_int(optarg, &num_iters, errmsg))
               error(errmsg, c);
            if(num_iters < 0)
               error("number of iterations must be greater than 0", c);
            break;

         case 'N':
            if(!read_int(optarg, &num_pts, errmsg))
               error(errmsg, c);
            if(num_pts < 2)
               error("number of points must be 2 or more", c);
            break;

         case 's':
            if(!read_double(optarg, &shorten_by, errmsg))
               error(errmsg, c);
            if(shorten_by <= 0 || shorten_by >=100)
               warning("not inside range 0 to 100", c);
            break;

         case 'l':
            if(!read_int(optarg, &lim_exp, errmsg))
               error(errmsg, c);
            if(lim_exp < 0) {
               warning("limit is negative, and so ignored", c);
            }
            if(lim_exp > 16) {
               warning("limit is very small, may not be attainable", c);
            }
            break;

         case 'r':
            if(!read_int(optarg, &rep_form, errmsg))
               error(errmsg, c);
            if(rep_form < 1 || rep_form > 4)
               error("formula is given by its number, 1 - 4", c);
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
   
   if(argc-optind == 1) {
      if(num_pts>0)
         error("cannot specify a number of points and also to read them from a file", 'N');
      else
         ifile=argv[optind];
   }

}

typedef vec3d (*REPEL_FN)(vec3d, vec3d);

vec3d rep_inv_dist1(vec3d v1, vec3d v2)
{ double len = pow((v2-v1).mag2(), -0.5); return (v2-v1).unit()*len; }

vec3d rep_inv_dist2(vec3d v1, vec3d v2)
{ double len = 1/(v2-v1).mag2(); return (v2-v1).unit()*len; }

vec3d rep_inv_dist3(vec3d v1, vec3d v2)
{ double len = pow((v2-v1).mag2(), -1.5); return (v2-v1).unit()*len; }

vec3d rep_inv_dist05(vec3d v1, vec3d v2)
{ double len = pow((v2-v1).mag2(), -0.25); return (v2-v1).unit()*len; }


void random_placement(geom_v &geom, int n)
{
   vector<vec3d> &verts = *geom.get_verts();
   verts.clear();
   time_t t;
   time(&t);
   srand(t);
   for(int i=0; i<n; i++)
      verts.push_back(vec3d::random().unit());
}

void repel(col_geom_v &geom, REPEL_FN rep_fn, double shorten_factor, double limit, int n)
{
   vector<vec3d> &verts = *geom.get_verts();
   vector<int> wts(verts.size());
   for(unsigned int i=0; i<verts.size(); i++) {
      col_val col = geom.get_v_col(i);
      wts[i] = col.is_idx() ? col.get_idx() : 1;
   }
   vector<vec3d> offsets(verts.size());
   double dist2, max_dist2=0;
   double last_av_max_dist2=0, max_dist2_sum=0;
   bool adaptive = false;
   int chng_cnt=0;
   int converge = 0;
   const int anum = 50;
   if(shorten_factor<0) {
      adaptive = true;
      shorten_factor = 0.001;
   }
   
   fprintf(stderr, "\n   ");
   
   unsigned int cnt;
   for(cnt=0; cnt< (unsigned int)n; cnt++) {
      std::fill(offsets.begin(), offsets.end(), vec3d(0,0,0));
      max_dist2 = 0;
      
      for(unsigned int i=0; i<verts.size()-1; i++) {
         for(unsigned int j=i+1; j<verts.size(); j++) {
            vec3d offset = rep_fn(verts[i], verts[j])*(wts[i]*wts[j]);
            offsets[i] -= offset/wts[i];
            offsets[j] += offset/wts[j];
         }
      }

      for(unsigned int i=0; i<verts.size(); i++) {
         //verts[i].dump("vert");
         //offsets[i].dump("   offset");
         vec3d new_pos = (verts[i]+offsets[i]*shorten_factor).unit();
         dist2 = (new_pos - verts[i]).mag2();
         if(dist2 > max_dist2)
            max_dist2 = dist2;
         if(dist2 > 0.001)
            new_pos = (new_pos+verts[i]).unit();
         verts[i] = new_pos;
      }
      
      if(sqrt(max_dist2) < limit)
         break;

      if(adaptive) {
         max_dist2_sum += max_dist2;
         if(max_dist2<last_av_max_dist2)
            converge +=1;
         if(!(cnt%anum)) {
            if(converge>anum/1.5) {
               chng_cnt = chng_cnt>0 ? chng_cnt+1 : 1;
               shorten_factor *= 1+0.005*chng_cnt;
               if(shorten_factor > 0.5)
                  shorten_factor = 0.5;
            }
            else if(converge<anum/2) {
               chng_cnt = chng_cnt<0 ? chng_cnt-1 : -2;
               shorten_factor *= 0.995+0.005*chng_cnt;
               if(shorten_factor < 0)
                  shorten_factor = 0.00001;
            }
            else
               chng_cnt = 0;

            last_av_max_dist2 = max_dist2_sum / anum;
            fprintf(stderr, "%2d ", converge/5);
            //fprintf(stderr, "\n%d %g %g (%g)\n", converge, shorten_factor, last_av_max_dist2, max_dist2);
            converge = 0;
            max_dist2_sum = 0;
         }
      }
      
      //if((cnt+1)%100 == 0)
      //   fprintf(stderr, ".");
      if((cnt+1)%1000 == 0) {
         double offset_sum=0;
         for(unsigned int i=0; i<offsets.size(); i++)
            offset_sum += offsets[i].mag();
         //fprintf(stderr, "\n%-15d %12.10g %g\n   ", cnt+1, sqrt(max_dist2), shorten_factor);
         fprintf(stderr, "\n%-13d  movement=%13.10g  s=%7.6g  F-sum=%.10g\n   ", cnt+1, sqrt(max_dist2), shorten_factor, offset_sum);
      }
   }

   if((cnt)%1000 != 0) {
      double offset_sum=0;
      for(unsigned int i=0; i<offsets.size(); i++)
         offset_sum += offsets[i].mag();
      fprintf(stderr, "\n%-13d  movement=%13.10g  s=%7.6g  F-sum=%.10g\n   ", cnt, sqrt(max_dist2), shorten_factor, offset_sum);
   }
   fprintf(stderr, "\n");
}
            
      

int main(int argc, char *argv[])
{
   rep_opts opts;
   opts.process_command_line(argc, argv);

   char errmsg[MSG_SZ] = "";
   col_geom_v geom;
   if(opts.num_pts>0)
      random_placement(geom, opts.num_pts);
   else {
      if(!geom.read(opts.ifile, errmsg))
         opts.error(errmsg);
      if(*errmsg)
         opts.warning(errmsg);
   }
   

   REPEL_FN fn[] = {rep_inv_dist1, rep_inv_dist2,
      rep_inv_dist3, rep_inv_dist05};
   repel(geom, fn[opts.rep_form-1], opts.shorten_by/100,
         pow(10, -opts.lim_exp), opts.num_iters);

   if(!geom.write(opts.ofile, errmsg))
      opts.error(errmsg);

   return 0;
}
   

