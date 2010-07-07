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
   Name: polygon.cc
   Description: program to generate polyhedra based on polygons
   Project: Antiprism - http://www.antiprism.com
*/


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include <string>

#include "../base/antiprism.h"


using std::string;

class pg_opts: public prog_opts {
   public:
      enum {t_polygon, t_prism, t_antiprism, t_pyramid,
         t_dipyramid, t_cupola, t_orthobicupola, t_gyrobicupola,
         t_dihedron, t_snub_antiprism};
      int type;
      int subtype;
      double edge;
      double circumrad;
      double inrad;
      double height;
      double edge2;
      int num_sides;
      int fraction;
      double twist_ang;
      bool make_trapezo;

      string ofile;

      bool e_given;
      bool h_given;
      bool a_given;

      pg_opts(): prog_opts("polygon"), type(t_prism), subtype(-1),
                 edge(1.0), circumrad(0.0), inrad(0.0),
                 height(-1), edge2(-1),
                 num_sides(5), fraction(1), twist_ang(0.0),
                 make_trapezo(false),
                 e_given(false), h_given(false), a_given(false) {}
      void process_command_line(int argc, char **argv);
      void usage();
};



void pg_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] num_sides\n"
"\n"
"Make polyhedra based on polygons: polygons, prisms, antiprisms,\n"
"pyramids dipyramids, cupolas, orthobicupolas, gyrobicupolas and\n"
"snub-antiprisms. num_sides is a number (N) optionally followed by / and\n"
"a second number (N/M). N is the number of vertices spaced equally on a\n"
"circle, and each vertex is joined to the Mth (default: 1) vertex moving\n"
"around the circle.\n"
"\n"
"Options\n"
"%s"
"  -t <type> type can be polygon, prism (default), antiprism, pyramid\n"
"            dipyramid, cupola, orthobicupola, gyrobicupola, dihedron or\n"
"            snub-antiprism. Type can be abreviated to three letters.\n"
"  -s <subt> a number indicting the subtype of a polyhedron (snub-antiprisms\n"
"            only, 0 or 1)\n"
"  -e <len>  length of polygon edges (default: 1)\n"
"  -R <rad>  circumradius of polygon (default: use -e)\n"
"  -r <rad>  inradius of polygon (default: use -e)\n"
"  -l <hgt>  height of upper vertices above base polygon\n"
"            (default: circumradius of polygon\n"
"  -E <len>  length of non-polygon edges (default: use -l)\n"
"  -a <ang>  angle (degrees) to rotate antiprism polygons, in opposite\n"
"            directions\n"
"  -T        for a given pyramid or antiprism make a trapezohedron that\n"
"            includes it as a part\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}

void pg_opts::process_command_line(int argc, char **argv)
{
   char errmsg[MSG_SZ];
   opterr = 0;
   char c;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":ht:e:r:R:E:l:s:o:a:T")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         case 't':
            if(strlen(optarg)<3)
               error("solid type is '"+string(optarg)+"', must be at least 3 letters");
            if(strncmp("polygon", optarg, 3)==0)
               type = t_polygon;
            else if(strncmp("prism", optarg, 3)==0)
               type = t_prism;
            else if(strncmp("antiprism", optarg, 3)==0)
               type = t_antiprism;
            else if(strncmp("pyramid", optarg, 3)==0)
               type = t_pyramid;
            else if(strncmp("dipyramid", optarg, 3)==0)
               type = t_dipyramid;
            else if(strncmp("cupola", optarg, 3)==0)
               type = t_cupola;
            else if(strncmp("orthobicupola", optarg, 3)==0)
               type = t_orthobicupola;
            else if(strncmp("gyrobicupola", optarg, 3)==0)
               type = t_gyrobicupola;
            else if(strncmp("snub-antiprism", optarg, 3)==0)
               type = t_snub_antiprism;
            else if(strncmp("dihedron", optarg, 3)==0)
               type = t_dihedron;
            else
               error("unknown solid type '"+string(optarg)+"'");
            break;

         case 's':
            if(!read_int(optarg, &subtype, errmsg))
               error(errmsg, c);
            if(subtype<0)
               error("subtype cannot be negative", c);
            break;
            
         
         case 'l':
            if(h_given)
               error("can only use one of options -E, -l", "non-polygon edge");
            if(!read_double(optarg, &height, errmsg))
               error(errmsg, c);
            if(height<0)
               error("height cannot be negative", c);
            break;
            
         case 'E':
            if(h_given)
               error("can only use one of options -E, -l", "non-polygon edge");
            if(!read_double(optarg, &edge2, errmsg))
               error(errmsg, c);
            if(edge2<0)
               error("edge length cannot be negative", c);
            break;
            
          case 'e':
            if(e_given)
               error("can only use one of options -e, -R, -r", "polygon edge");
            if(!read_double(optarg, &edge, errmsg))
               error(errmsg, c);
            if(edge<0)
               error("edge length cannot be negative", c);
            break;
            
         case 'R':
            if(e_given)
               error("can only use one of options -e, -R, -r", "polygon edge");
            if(!read_double(optarg, &circumrad, errmsg))
               error(errmsg, c);
            if(circumrad<0)
               error("radius cannot be negative", c);
            break;
            
         case 'r':
            if(e_given)
               error("can only use one of options -e, -R, -r", "polygon edge");
            if(!read_double(optarg, &inrad, errmsg))
               error(errmsg, c);
            if(inrad<0)
               error("radius cannot be negative", c);
            break;

         case 'a':
            a_given = true;
            if(!read_double(optarg, &twist_ang, errmsg))
               error(errmsg, c);
            twist_ang = deg2rad(twist_ang);
            break;

         case 'T':
            make_trapezo = true;
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
  
   if(argc-optind > 0) {
      char *p = strchr(argv[optind], '/');
      if(p!=0) {
         *p++='\0';
         if(!read_int(p, &fraction, errmsg))
            error(errmsg, "number of sides (fractional part)");
      }
              
      if(!read_int(argv[optind], &num_sides, errmsg))
         error(errmsg, "number of sides");
      if(num_sides<2)
         error("must be an integer 2 or greater", "number of sides");
      if(fraction < 1)
         error("fractional part must be 1 or greater", "number of sides");
      if(fraction >= num_sides)
         error("fractional part must be less than number of sides", "number of sides");
   }

   if(a_given && type!=t_antiprism)
      error("angle option can only be used with antiprisms", "a");
   
   if(make_trapezo && !(type==t_antiprism || type==t_pyramid))
      error("only a pyramid or an antiprisms can be made into a trapezohedron"
            , "T");
}



int main(int argc, char *argv[])
{
   pg_opts opts;
   opts.process_command_line(argc, argv);

   polygon pgon(opts.num_sides, opts.fraction);

   // make n-polygon in xz plane with centre at origin
   if(opts.circumrad!=0)
      pgon.set_radius(opts.circumrad);
   else if(opts.inrad!=0)
      pgon.set_inradius(opts.inrad);
   else
      pgon.set_edge(opts.edge);

   polygon *poly=0;
   switch(opts.type) {
      case pg_opts::t_polygon:
         poly = new polygon(pgon);
         break;
      case pg_opts::t_prism:
         poly = new prism(pgon);
         break;
      case pg_opts::t_antiprism:
         {
            antiprism *ant = new antiprism(pgon);
            ant->set_twist_angle(opts.twist_ang);
            ant->set_output_trapezohedron(opts.make_trapezo);
            poly = ant;
            break;
         }
      case pg_opts::t_pyramid:
         {
            pyramid *pyr = new pyramid(pgon);
            pyr->set_output_trapezohedron(opts.make_trapezo);
            poly = pyr;
            break;
         }
      case pg_opts::t_dipyramid:
         poly = new dipyramid(pgon);
         break;
      case pg_opts::t_cupola:
         poly = new cupola(pgon);
         break;
      case pg_opts::t_orthobicupola:
         poly = new orthobicupola(pgon);
         break;
      case pg_opts::t_gyrobicupola:
         poly = new gyrobicupola(pgon);
         break;
      case pg_opts::t_dihedron:
         poly = new dihedron(pgon);
         break;
      case pg_opts::t_snub_antiprism:
         poly = new snub_antiprism(pgon);
         break;
   }

   char errmsg[MSG_SZ];
   if(opts.subtype>-1) {
      if(!poly->set_subtype(opts.subtype, errmsg))
         opts.warning(errmsg, 's');
   }
   if(opts.height>=0) {
      if(!poly->set_height(opts.height, errmsg))
         opts.warning(errmsg, 'l');
   }
   else if(opts.edge2>=0)
      if(!poly->set_edge2(opts.edge2, errmsg))
         opts.warning(errmsg, 'E');
         
   geom_v geom;
   poly->make_poly(geom);
   delete poly;
   if(!geom)
      opts.error("a polyhedron with these parameters could not be made");

   if(!geom.write(opts.ofile, errmsg))
      opts.error(errmsg);

   return 0;
}


