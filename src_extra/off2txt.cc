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

#include <string>
#include <vector>

#include "../base/antiprism.h"


using std::string;
using std::vector;


class o2t_opts: public prog_opts {
   public:
      string ifile;
      string ofile;

      bool estimate_colors;
      bool detect_rhombi;
      bool detect_star_polygons;
      bool exclude_coordinates;
      bool force_transparent;
      double epsilon;
      int sig_digits;

      o2t_opts(): prog_opts("off2txt"),
             estimate_colors(false),
             detect_rhombi(false),
             detect_star_polygons(false),
             exclude_coordinates(false),
             force_transparent(false),
             epsilon(0),
             sig_digits(DEF_SIG_DGTS)
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
"  -l <lim>  minimum distance for unique vertex locations as negative exponent\n"
"               (default: %d giving %.0e)\n"
"  -d <dgts> number of significant digits (default %d) or if negative\n"
"            then the number of digits after the decimal point\n"
"  -o <file> file name for output (otherwise prints to stdout)\n"
"\n"
"\n", prog_name(), help_ver_text, int(-log(::epsilon)/log(10) + 0.5),::epsilon, DEF_SIG_DGTS);
}


void o2t_opts::process_command_line(int argc, char **argv)
{
   opterr = 0;
   char c;
   char errmsg[MSG_SZ];

   int sig_compare = INT_MAX;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hcrpxtl:d:o:")) != -1 ) {
      if(common_opts(c, optopt))
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

   epsilon = (sig_compare != INT_MAX) ? pow(10, -sig_compare) : ::epsilon;   
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

bool is_square(geom_if &geom, const int &face_idx, const double &eps)
{
   const vector<vec3d> &verts = geom.verts();
   const vector<vector<int> > &faces = geom.faces();
   vector<int> face = faces[face_idx];

   vec3d P1 = verts[face[0]];
   vec3d P2 = verts[face[2]];
   double diag1 = (P2-P1).mag();

   P1 = verts[face[1]];
   P2 = verts[face[3]];
   double diag2 = (P2-P1).mag();

   return(double_eq(diag1, diag2, eps));
}

int detect_star_polygon(geom_if &geom, const int &face_idx)
{
   const vector<vec3d> &verts = geom.verts();
   const vector<vector<int> > &faces = geom.faces();
   vector<int> face = faces[face_idx];

   vec3d v1 = verts[face[1]] - verts[face[0]];
   vec3d v2 = verts[face[1]] - verts[face[2]];
   double angle = acos(safe_for_trig(vdot(v1, v2)/(v1.mag()*v2.mag())));
   angle = rad2deg(angle);

   int star = 1;
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

string Vtxt(const vec3d &v, const int &dgts)
{
   char buf[128];
   if(dgts>0)
      snprintf(buf, 128, "%.*g, %.*g, %.*g,", dgts, v[0], dgts, v[1], dgts,v[2]);
   else
      snprintf(buf, 128, "%.*f, %.*f, %.*f,", -dgts, v[0], -dgts, v[1], -dgts, v[2]);
   return buf;
}

void print_hedron_txt(FILE *ofile, col_geom_v &geom, const int &sig_digits, const bool &estimate_colors, const bool &detect_rhombi, 
                      const bool &detect_star_polygons, const bool &exclude_coordinates, const bool &force_transparent, const double &eps)
{
   const vector<vector<int> > &faces = geom.faces();
   const vector<vec3d> &verts = geom.verts();

	fprintf(ofile, "{\n");
            
   for(unsigned int i=0; i<faces.size(); i++)
   {
      //fprintf(ofile, "");

      if ( detect_rhombi && faces[i].size()==4 && !is_square(geom, (int)i, eps) )
         fprintf(ofile, "D");

      if ( estimate_colors ) {
         col_val col = geom.get_f_col((int)i);
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

   print_hedron_txt(ofile, geom, opts.sig_digits, opts.estimate_colors, opts.detect_rhombi, 
      opts.detect_star_polygons, opts.exclude_coordinates, opts.force_transparent, opts.epsilon);

   if(opts.ofile!="")
      fclose(ofile);

   return 0;
}
