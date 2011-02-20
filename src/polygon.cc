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
#include <math.h>


#include <string>

#include "../base/antiprism.h"


using std::string;

class pg_opts: public prog_opts {
   public:
      enum {t_prism=1, t_antiprism, t_pyramid,
         t_dipyramid, t_cupola, t_orthobicupola, t_gyrobicupola,
         t_snub_antiprism, t_dihedron };
      int type;
      string subtype_str;
      int subtype;
      double edge;
      double radius;
      double radius2;
      double height;
      double height2;
      double edge2;
      int num_sides;
      int step;
      double twist_ang;

      string ofile;
      bool e_given;
      bool h_given;

      pg_opts(): prog_opts("polygon"), type(t_prism), subtype(-1),
                 edge(NAN), radius(NAN), radius2(NAN),
                 height(NAN), height2(NAN), edge2(NAN), step(1),
                 twist_ang(NAN),
                 e_given(false), h_given(false)
                 {}
      void process_command_line(int argc, char **argv);
      void usage();
};



void pg_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] type num_sides\n"
"\n"
"Make polyhedra based on polygons. The polyhedron type and subtype may be\n"
"given by (partial) name or number\n"
"   1. prism (angle: polygons twist, forming an antiprism)\n"
"        subtypes: 1. antiprism, from triangulating side faces\n"
"                  2. trapezohedron, with this antiprism-ised belt\n"
"   2. antiprism (angle: polygons twist)\n"
"        subtypes: 1. trapezohedron, with this antiprism belt\n"
"                  2. antihermaphrodite, with this antiprism base\n"
"                  3. scalenohedron (-L for apex height)\n"
"                  4. subdivided_scalenohedron (-L for apex height)\n"
"   3. pyramid (angle: base separates, polygons twist to antihermaphrodite)\n"
"        subtypes: 1. antihermaphrodite, with this antiprism base\n"
"                  2. elongated (-L for prism height)\n"
"                  3. gyroelongated (-L for antiprism height)\n"
"   4. dipyramid (angle: base separates, polygons twist to trapezohedron)\n"
"        subtypes: 1. trapezohedron, with this pyramid apex\n"
"                  2. elongated (-L for prism height)\n"
"                  3. gyroelongated (-L for antiprism height)\n"
"                  4. dipyramid_scalenohedron (-R for alternate vertex radius)\n"
"   5. cupola\n"
"        subtypes: 1. elongated (-L for prism height)\n"
"                  2. gyroelongated (-L for antiprism height)\n"
"                  3. cupoloid\n"
"   6. orthobicupola\n"
"        subtypes: 1. elongated (-L for prism height)\n"
"                  2. gyroelongated (-L for antiprism height)\n"
"   7. gyrobicupola\n"
"        subtypes: 1. elongated (-L for prism height)\n"
"                  2. gyroelongated (-L for antiprism height)\n"
"   8. snub-antiprism\n"
"        subtypes: 1. inverted, triangle band inverted\n"
"   9. dihedron\n"
"        subtypes: 1. polygon\n"
"\n"
"num_sides is a number (N) optionally followed by / and a second\n"
"number (N/M). N is the number of vertices spaced equally on a\n"
"circle, and each vertex is joined to the Mth (default: 1) vertex\n"
"moving around the circle.\n"
"\n"
"Options\n"
"%s"
"  -s <subt> a number or name (see type list above) indicting a subtype\n"
"            or modification of a polyhedron\n"
"  -e <len>  length of polygon edges (default: 1)\n"
"  -E <len>  length of non-polygon edges (default: calculate from -l)\n"
"  -r <rad>  circumradius of polygon (default: calculate from -e)\n"
"  -R <rad>  a second radius value\n"
"  -l <hgt>  height of upper vertices above base polygon\n"
"            (default: circumradius of polygon)\n"
"  -L <hgt>  a second height or length value\n"
"  -a <ang>  twist angle (degrees)\n"
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

   while( (c = getopt(argc, argv, ":he:E:r:R:l:L:s:o:a:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         case 's':
            subtype_str = optarg;
            break;
         
         case 'l':
            if(h_given)
               error("can only use one of options -E, -l", "non-polygon edge");
            h_given = true;
            if(!read_double(optarg, &height, errmsg))
               error(errmsg, c);
            break;

          case 'L':
            if(!read_double(optarg, &height2, errmsg))
               error(errmsg, c);
            break;
            
         case 'E':
            if(h_given)
               error("can only use one of options -E, -l", "non-polygon edge");
            h_given = true;
            if(!read_double(optarg, &edge2, errmsg))
               error(errmsg, c);
            if(edge2<0)
               error("edge length cannot be negative", c);
            break;
            
          case 'e':
            if(e_given)
               error("can only use one of options -e, -r", "polygon edge");
            e_given = true;
            if(!read_double(optarg, &edge, errmsg))
               error(errmsg, c);
            if(edge<0)
               error("edge length cannot be negative", c);
            break;
            
         case 'r':
            if(e_given)
               error("can only use one of options -e, -r", "polygon edge");
            e_given = true;
            if(!read_double(optarg, &radius, errmsg))
               error(errmsg, c);
            break;

         case 'R':
            if(!read_double(optarg, &radius2, errmsg))
               error(errmsg, c);
            break;
            
         case 'a':
            if(!read_double(optarg, &twist_ang, errmsg))
               error(errmsg, c);
            twist_ang = deg2rad(twist_ang);
            break;

         case 'o':
            ofile = optarg;
            break;

         default:
            error("unknown command line error");
      }
   }
   
   if(argc-optind != 2)
      error("must give two arguments");

   // Determine type
   string params = "|prism|antiprism|pyramid|dipyramid|cupola|orthobicupola"
                   "|gyrobicupola|snub-antiprism|dihedron";
             
   string arg_id = get_arg_id(argv[optind], params.c_str(),
         argmatch_default | argmatch_add_id_maps, errmsg);
   if(arg_id=="")
      error(errmsg, "polyhedron type");

   type = atoi(arg_id.c_str());
       
   // Determine subtype
   if(subtype_str != "") {
      params="0";
      if(type==t_prism)
         params += "|antiprism|trapezohedron";
      else if(type==t_antiprism)
         params += "|trapezohedron|antihermaphrodite|scalenohedron"
                   "|subdivided_scalenohedron";
      else if(type==t_pyramid)
         params += "|antihermaphrodite|elongated|gyroelongated";
      else if(type==t_dipyramid)
         params += "|trapezohedron|elongated|gyroelongated"
                   "|dipyramid_scalenohedron";
      else if(type==t_orthobicupola ||type==t_gyrobicupola)
         params += "|elongated|gyroelongated";
      else if(type==t_cupola)
         params += "|elongated|gyroelongated|cupoloid";
      else if(type==t_snub_antiprism)
         params += "|inverted";
      else if(type==t_dihedron)
         params += "|polygon";

      arg_id = get_arg_id(subtype_str.c_str(), params.c_str(),
            argmatch_default | argmatch_add_id_maps, errmsg);
      if(arg_id=="")
         error(errmsg, "polyhedron subtype");

         subtype = atoi(arg_id.c_str());
   }

 
   // read polygon   
   char *p = strchr(argv[optind+1], '/');
   if(p!=0) {
      *p++='\0';
      if(!read_int(p, &step, errmsg))
         error(errmsg, "number of sides (denominator of polygon fraction)");
   }

   if(!read_int(argv[optind+1], &num_sides, errmsg))
      error(errmsg, "number of sides");
   if(num_sides<2)
      error("must be an integer 2 or greater", "number of sides");
   if(step < 1)
      error("denominator of polygon fraction must be 1 or greater", "number of sides");
   if(step%num_sides==0)
      error("denominator of polygon fraction cannot be a multiple of the number of sides", "number of sides");

}



int main(int argc, char *argv[])
{
   pg_opts opts;
   opts.process_command_line(argc, argv);

   polygon pgon(opts.num_sides, opts.step);

   // make n-polygon in xz plane with centre at origin
   if(!isnan(opts.radius))
      pgon.set_radius(opts.radius);
   else if(!isnan(opts.edge))
      pgon.set_edge(opts.edge);
   else
      pgon.set_edge(1.0);
   
   polygon *poly=0;
   switch(opts.type) {
      case pg_opts::t_prism:
         poly = new prism(pgon);
         break;
      case pg_opts::t_antiprism:
         poly = new antiprism(pgon);
         break;
      case pg_opts::t_pyramid:
         poly = new pyramid(pgon);
         break;
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
         opts.error(errmsg, 's');
   }
   
   if(opts.type==pg_opts::t_cupola && opts.subtype==3 &&
         ((poly->get_step()*poly->get_parts())%2) )
      opts.error("cannot make a cupoloid unless the specified polygon fraction (in lowest form) has an even denominator", "polygon type");
         
   if(poly->get_num_sides()==2 && (opts.subtype==1||!isnan(opts.twist_ang)) ) {
      if(opts.type==pg_opts::t_pyramid)
         opts.error("cannot make an antihermaphrodite from a digonal pyramid",
               "subtype or twist angle");
      if(opts.type==pg_opts::t_dipyramid)
         opts.error("cannot make a trapezohedron from a digonal dipyramid",
               "subtype or twist angle");
   }
         
    
   if(!isnan(opts.twist_ang)) {
      if(!poly->set_twist_angle(opts.twist_ang, errmsg))
         opts.error(errmsg, 'a');
   }
   
   if(!isnan(opts.height)) {
      if(!poly->set_height(opts.height, errmsg))
         opts.error(errmsg, 'l');
   }
   else if(!isnan(opts.edge2))
      if(!poly->set_edge2(opts.edge2, errmsg))
         opts.error(errmsg, 'E');

   if(!poly->set_height2(opts.height2, errmsg))
         opts.error(errmsg, 'L');
   if(!poly->set_radius2(opts.radius2, errmsg))
         opts.error(errmsg, 'R');

   geom_v geom;
   poly->make_poly(geom);
   delete poly;
   if(!geom)
      opts.error("a polyhedron with these parameters could not be made");

   if(!geom.write(opts.ofile, errmsg))
      opts.error(errmsg);

   return 0;
}


