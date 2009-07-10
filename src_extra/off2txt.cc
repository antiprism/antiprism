/*
   Copyright (c) 2007-2009, Roger Kaufman

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
   Name: off2txt.cc
   Description: Convert files in OFF format to Hedron format
   Project: Antiprism - http://www.antiprism.com
*/



#include <stdio.h>
#include <stdlib.h>

#include <ctype.h>
#include <unistd.h>

#include <string>
#include <vector>

#include "../base/antiprism.h"


using std::string;
using std::vector;


class o2t_opts: public prog_opts {
   public:
      bool estimate_colors;
      bool detect_rhombi;
      bool detect_star_polygons;
      bool exclude_coordinates;
      bool force_transparent;
      int sig_compare;
      int sig_digits;
      string ifile;
      string ofile;

      o2t_opts(): prog_opts("off2txt"),
             estimate_colors(false),
             detect_rhombi(false),
             detect_star_polygons(false),
             exclude_coordinates(false),
             force_transparent(false),
             sig_compare(8),
             sig_digits(17)
             {}

      void process_command_line(int argc, char **argv);
      void usage();
};



void o2t_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Convert files in OFF format Hedron text file format. If\n"
"input_file is not given the program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"  -c        estimates colors from OFF file\n"
"  -r        detect rhombi. D parameter added if found\n"
"  -p        detect star polygons\n"
"  -x        exclude coordinates\n"
"  -t        force all faces transparent\n"
"  -l <lim>  minimum distance for rhombi detection\n"
"            exponent 1e-lim (default: 8 giving 1e-8)\n" 
"  -o <file> file name for output (otherwise prints to stdout)\n"
"\n"
"  Precision options\n"
"  -d <dgts> number of significant digits (default 17) or if negative\n"
"            then the number of digits after the decimal point\n"
"\n"
"\n", prog_name(), help_ver_text);
}


void o2t_opts::process_command_line(int argc, char **argv)
{
   extern char *optarg;
   extern int optind, opterr;
   opterr = 0;
   char c;
   char errmsg[MSG_SZ];
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hcrpxtl:d:o:")) != -1 ) {
      if(common_opts(c))
         continue;

      switch(c) {
         case 'c':
            estimate_colors = true;
            break;

         case 'r':
            detect_rhombi = true;
            break;

         case 'p':
            detect_star_polygons = true;
            break;

         case 'x':
            exclude_coordinates = true;
            break;

         case 't':
            force_transparent = true;
            break;

         case 'l':
            if(!read_int(optarg, &sig_compare, errmsg))
               error(errmsg, c);
            if(sig_compare < 0) {
               warning("limit is negative, and so ignored", c);
            }
            if(sig_compare > 16) {
               warning("limit is very small, may not be attainable", c);
            }
            break;

         case 'd':
            if(!read_int(optarg, &sig_digits, errmsg))
               error(errmsg, c);
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
      ifile=argv[optind];
   
}

string estimated_color(col_val col)
{
   string color;
   int v[3];

   for(unsigned int i=0; i<3; i++)
      v[i] = ( col[i] <= 128 ) ? 0 : 255;

   // Special case for orange
	if ( v[0] == 255 && ( col[1] >= 64 && col[1] <= 192 ) && v[2] == 0 )
		color = "a";
	else
	if ( v[0] == 0 && v[1] == 0 && v[2] == 255 )
		color = "b";
	else
   // Special case for green. "0,128,0" is exact green would become 0 0 0
   // Make more dark greens be g instead of using 0 255 0.
	if ( v[0] == 0 && col[1] >= 64 && v[2] == 0 )
		color = "g";
	else
   // Hedron has no black so specify cyan
	if ( v[0] == 0 && v[1] == 0 && v[2] == 0 )
		color = "c";
   else
	if ( v[0] == 0 && v[1] == 255 && v[2] == 255 )
		color = "c";
	else
	if ( v[0] == 255 && v[1] == 0 && v[2] == 0 )
		color = "r";
	else
	if ( v[0] == 255 && v[1] == 0 && v[2] == 255 )
		color = "m";
	else
	if ( v[0] == 255 && v[1] == 255 && v[2] == 0 )
		color = "y";
   else
	if ( v[0] == 255 && v[1] == 255 && v[2] == 255 )
		color = "w";
   else
      color = "\0";

   return(color);
}

bool doubleEquality(const double &d1, const double &d2, const double &epsilon)
{
   const double diff = d1 - d2;
   return diff < epsilon && diff > -epsilon;
}

bool is_square(geom_if *geom, int face_idx, double epsilon)
{
   vector<vec3d> &verts = *geom->get_verts();
   vector<vector<int> > &faces = *geom->get_faces();
   vector<int> face = faces[face_idx];

   vec3d P1 = verts[face[0]];
   vec3d P2 = verts[face[2]];
   double diag1 = (P2-P1).mag();

   P1 = verts[face[1]];
   P2 = verts[face[3]];
   double diag2 = (P2-P1).mag();

   return(doubleEquality(diag1, diag2, epsilon));
}

int detect_star_polygon(geom_if *geom, int face_idx)
{
   vector<vec3d> &verts = *geom->get_verts();
   vector<vector<int> > &faces = *geom->get_faces();
   vector<int> face = faces[face_idx];
   
   int star = 1;
   static const double radian = 57.295779513082320876798;

   vec3d v1 = verts[face[1]] - verts[face[0]];
   vec3d v2 = verts[face[1]] - verts[face[2]];
   double angle = acos(vdot(v1, v2)/(v1.mag()*v2.mag()));
   angle *= radian;

   double m = face.size();
   for(double n=2;n<m/2;n++) {
      if((int)(m/n)*n!=m && gcd((int)m,(int)n)==1){
         double sp=180*(1-2*n/m);
         // angle allowed plus/minus 2 degrees slop
         if (fabs(sp-angle)<=2) {
            star = (int)n;
            break;
         }
      }
   }
   return(star);
}

string Vtxt(vec3d v, int dgts)
{
   char buf[128];
   if(dgts>0)
      snprintf(buf, 128, "%.*g, %.*g, %.*g,", dgts, v[0], dgts, v[1], dgts,v[2]);
   else
      snprintf(buf, 128, "%.*f, %.*f, %.*f,", -dgts, v[0], -dgts, v[1], -dgts, v[2]);
   return buf;
}

void print_hedron_txt(FILE *ofile, geom_if *geom, int sig_digits, bool estimate_colors, bool detect_rhombi, 
                      bool detect_star_polygons, bool exclude_coordinates, bool force_transparent, double epsilon)
{
   vector<vector<int> > &faces = *geom->get_faces();
   vector<vec3d> &verts = *geom->get_verts();
   col_geom *cg = dynamic_cast<col_geom *>(geom);

	fprintf(ofile, "{\n");
            
   for(unsigned int i=0; i<faces.size(); i++)
   {
      //fprintf(ofile, "");

      if ( detect_rhombi && faces[i].size()==4 && !is_square(geom, (int)i, epsilon) )
         fprintf(ofile, "D");

      if ( estimate_colors ) {
         col_val col = cg->get_f_col((int)i);
         if ( col.is_val() )
		      fprintf(ofile, "%s", estimated_color(col).c_str());
      }

      if ( force_transparent ) {
		      fprintf(ofile, "t");
      }

      for(unsigned int j=0; j<faces[i].size(); j++)
		   fprintf(ofile, "%d,", faces[i][j]);

      int star = 1;
      if ( detect_star_polygons && faces[i].size()>4 )
         star = detect_star_polygon(geom, (int)i);
      fprintf(ofile, "-%d,\n",star );
    }

   if ( !exclude_coordinates ) {
      fprintf(ofile, "(\n");

      for(unsigned int i=0; i<verts.size(); i++)
         fprintf(ofile, "%s\n", Vtxt(verts[i], sig_digits).c_str());

      fprintf(ofile, ")\n");
   }

   fprintf(ofile, "}\n");
}

int main(int argc, char *argv[])
{
   o2t_opts opts;
   opts.process_command_line(argc, argv);

   char errmsg[MSG_SZ];
   col_geom_v geom;
   if(!geom.read(opts.ifile, errmsg))
      opts.error(errmsg);
   if(*errmsg)
      opts.warning(errmsg);
      
   FILE *ofile = stdout;  // write to stdout by default
   if(opts.ofile != "") {
      ofile = fopen(opts.ofile.c_str(), "w");
      if(ofile == 0)
         opts.error("could not open output file \'"+opts.ofile+"\'");
   }

   print_hedron_txt(ofile, &geom, opts.sig_digits, opts.estimate_colors, opts.detect_rhombi, 
      opts.detect_star_polygons, opts.exclude_coordinates, opts.force_transparent, pow(10, -opts.sig_compare));

   if(opts.ofile!="")
      fclose(ofile);

   return 0;
}
