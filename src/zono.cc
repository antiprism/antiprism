/*
   Copyright (c) 2003-2014, Adrian Rossiter

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
#include <math.h>
#include <ctype.h>

#include <string>
#include <vector>
#include <algorithm>

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
      col_val zone_col;
      int pol_num;
      bool non_polar_opt;
      string ifile;
      col_geom_v seed_geom;
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
"Usage: %s [options] [star_file] \n"
"\n"
"Make a zonohedron or add zones to a convex seed polyhdron. The zones are\n"
"created from a star of vectors, which can be based on a polyhdron (input\n"
"model and option -m) or initialised (option -P) to make a polar zonohedron.\n"
"If input_file is not given the program reads from standard input\n"
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
"  -S        seed model to add zones to, must be convex\n"
"  -u        make vectors unit length\n"
"  -C <col>  colour for new zone faces\n"
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

   while( (c = getopt(argc, argv, ":hm:c:S:suC:P:o:")) != -1 ) {
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

         case 'S':
         {
            geom_read_or_error(seed_geom, optarg, *this);
            col_geom_v convex_chk = seed_geom;
            convex_chk.set_hull();
            if(!check_congruence(seed_geom, convex_chk))
               error("seed geometry is not convex", c);
            break;
         }

         case 'u':
            unit_len = true;
            break;

         case 'C':
            if(!zone_col.read(optarg, errmsg))
               error(errmsg, c);
            break;

         case 'P':
            if(!read_int(optarg, &pol_num, errmsg))
               error(errmsg, c);
            if(pol_num < 2)
               error("must be 2 or greater", c);
            break;

         case 'o':
            ofile = optarg;
            break;

         default:
            error("unknown command line error");
      }
   }

   if(pol_num && (non_polar_opt||argc-optind))
      error("option -P, is not compatible with options m, c or an input file");

   if(argc-optind > 1)
      error("too many arguments");

   if(argc-optind == 1)
      ifile=argv[optind];
}


struct vec_less {
   bool operator() (const vec3d &v1, const vec3d &v2) {
      return compare(v1, v2, epsilon)==-1;
   }
};


bool make_zonohedron_with_seed(geom_if &zono, const vector<vec3d> &star,
      const geom_if &seed, col_val col, char *errmsg)
{
   zono.clear_all();
   zono.append(seed);

   // Store original face colours by normal
   map<vec3d, col_val, vec_less> orig_cols;
   col_geom *cg = dynamic_cast<col_geom *>(&zono);
   if(cg) {
      for(unsigned int i=0; i<zono.faces().size(); i++) {
         orig_cols[zono.face_norm(i).to_unit()] = cg->get_f_col(i);
      }
   }

   for(unsigned int i=0; i<star.size(); i++) {
      int v_sz = zono.verts().size();
      zono.raw_verts().resize(v_sz*2);
      for(int j=0; j<v_sz; j++) {
         zono.raw_verts()[j+v_sz] = zono.verts(j) + star[i];
      }
      if(!zono.set_hull("", errmsg))
         return false;
   }

   // Restore original face colours by normal
   if(cg) {
      for(unsigned int i=0; i<zono.faces().size(); i++) {
         map<vec3d, col_val, vec_less>::iterator mi =
               orig_cols.find(zono.face_norm(i).to_unit());
         cg->set_f_col(i, (mi != orig_cols.end()) ? mi->second : col);
      }
   }

   return true;
}



bool make_zonohedron(geom_if &zono, const vector<vec3d> &star, char *errmsg)
{
   zono.clear_all();

   if(star.size()<2) {
      snprintf(errmsg, MSG_SZ, "star contains %lu vectors, needs at least 2",
            (unsigned long)star.size());
      return false;
   }

   else
      return zono.set_zono(star, errmsg);
}


int main(int argc, char **argv)
{
   zo_opts opts;
   opts.process_command_line(argc, argv);
   col_geom_v zono;

   vector<vec3d> star;
   if(opts.pol_num) {
      for(int n=0; n<opts.pol_num; n++) {
         double ang = 2*M_PI*n/opts.pol_num;
         star.push_back(vec3d(cos(ang), sin(ang), 1) * sqrt(0.5));
      }
   }
   else {
      col_geom_v geom;
      geom_read_or_error(geom, opts.ifile, opts);

      if(opts.centroid)
         opts.centre = centroid(*geom.get_verts());

      star = geom.get_star(opts.method, opts.centre);
   }

   // Set star vectors to unit, and remove any parallel vectors
   if(opts.unit_len) {
      for(unsigned int i=0; i<star.size(); i++)
         star[i].to_unit();
      int star_final_sz = star.size();
      for(int i=0; i<star_final_sz; i++)
         for(int j=i+1; j<star_final_sz; j++) {
            if(vcross(star[i], star[j]).mag()<epsilon)
               std::swap(star[j], star[--star_final_sz]);
         }
      star.resize(star_final_sz);
   }

   char errmsg[MSG_SZ] = "";
   if(opts.out_star)
      zono.add_verts(star);
   else if(opts.seed_geom) {
      if(!make_zonohedron_with_seed(zono, star, opts.seed_geom,
               opts.zone_col, errmsg))
         opts.error(errmsg);
   }
   else {
      if(!make_zonohedron(zono, star, errmsg))
         opts.error(errmsg);
      if(opts.zone_col.is_set())
         coloring(&zono).f_one_col(opts.zone_col);
   }

   geom_write_or_error(zono, opts.ofile, opts);

   return 0;
}


