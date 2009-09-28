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
   Name: m2off.cc
   Description: Convert files in 'm' format for LiveGraphics3D to OFF file format
   Project: Antiprism - http://www.antiprism.com
*/


#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include <ctype.h>
#include <unistd.h>

#include <string>
#include <vector>

#include "../base/antiprism.h"

using std::string;
using std::vector;
 

class m2off_opts: public prog_opts {
   public:
      string ifile;
      string ofile;

      string hide_elems;
      int sig_compare;
      bool disallow_back_faces;
      bool live3D_do_viewpoint;
      string lights_geom_file;

      m2off_opts(): prog_opts("m2off"),
                    sig_compare(INT_MAX),
                    disallow_back_faces(true),
                    live3D_do_viewpoint(false)
             {}

      void process_command_line(int argc, char **argv);
      void usage();
};

void m2off_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Convert files in 'm' format, used by LiveGraphics3D, to OFF format. If\n"
"input_file is not given the program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"  -x <elms> hide elements. The element string can include v, e and f\n"
"               to hide, respectively, vertices, edges and faces\n"
"  -l <lim>  minimum distance for unique vertex locations as negative exponent\n"
"               (default: %d giving %.0e)\n"
"  -b        use back face colors instead of front ones, if available\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\nScene Options\n"
"  -v        turn model to LiveGraphics3D viewpoint\n"
"  -C <file> dump color lights into OFF file\n"
"\n"
"\n",prog_name(), help_ver_text, int(-log(::epsilon)/log(10) + 0.5),::epsilon);
}

void m2off_opts::process_command_line(int argc, char **argv)
{
   extern char *optarg;
   extern int optind, opterr;
   opterr = 0;
   char c;
   char errmsg[MSG_SZ];
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hx:l:bvC:o:")) != -1 ) {
      if(common_opts(c))
         continue;

      switch(c) {
         case 'x':
            if(strspn(optarg, "vef") != strlen(optarg)) {
               snprintf(errmsg, MSG_SZ, "elements to hide are %s must be v, e, or f\n", optarg);
               error(errmsg, c);
            }
            if(strlen(optarg)==3) {
               snprintf(errmsg, MSG_SZ, "cannot hide all v, e, and f. arg was: %s\n", optarg);
               error(errmsg, c);
            }
            hide_elems=optarg;
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

         case 'b':
            disallow_back_faces = false;
            break;

         case 'v':
            live3D_do_viewpoint = true;
            break;

         case 'C':
            lights_geom_file = optarg;
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
      
   sig_compare = (sig_compare != INT_MAX) ? sig_compare : int(-log(::epsilon)/log(10) + 0.5);
}

string replaceAllOccurrances(string s, string f, string r)
{
   size_t found = s.find(f);
   while(found != string::npos) {
      s.replace(found, f.length(), r);
      found = s.find(f);
   }
   return s;
}

string read_file_to_str(string file_name, char *errmsg)
{
   string str;

   if(errmsg)
      *errmsg='\0';

   FILE *ifile;
   if(file_name == "" || file_name == "-") {
      ifile = stdin;
      file_name = "stdin";
   }
   else {
      ifile = fopen(file_name.c_str(), "r");
      if(!ifile) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "could not open input file \'%s\'", file_name.c_str());
         return str;
      }
   }

   char *line=0;
   while(read_line(ifile, &line)==0) {
      for(char *p=line; *p; p++)   // convert whitespace to spaces
         if(isspace(*p))
            *p = ' ';

      str.append(line);
      str.append(" ");
      free(line);
   }
   
   if(file_name!="stdin")
      fclose(ifile);

   str = replaceAllOccurrances(str, "->", "  ");

   return str;
}

bool is_numeric(char *number)
{
   char tmp[1];
   float x;
   if (sscanf(number, "%f%c", &x, tmp) == 1)
      return true;
   return false;
}

void m_parse(string &m_txt, col_geom_v &geom, col_geom_v &geom_cv, string hide_elems, bool disallow_back_faces,
             col_geom_v &lights_geom, vec3d &view_point, vec3d &view_vertical, char *errmsg)
{
   const vector<vector<int> > &faces = geom.faces();
   const vector<vec3d> &verts = geom.verts();

   // seperate geom for colored verts
   const vector<vec3d> &verts_cv = geom_cv.verts();

   char parse_key[] = " ,{}[]";
   char parse_key_polygon[] = " ,{}[";
   bool delay_read = false;

   vector<char *> col_read;
   col_val current_col = 0;
   bool face_form = false;
   bool poly_backface_ignore = false;
   int color_case = 0;

   double coord[3];
   vector<int> face;

   size_t line_size = strlen((m_txt).c_str());

   char *buffer = (char*) malloc (line_size+1);
   if (!buffer) {
      snprintf(errmsg, MSG_SZ, "out of memory");
      return;
   }

   strncpy(buffer,(m_txt).c_str(),line_size);
   char *ptok = strtok(buffer,parse_key);

   // valid Live3D file?
   if ( strcmp(ptok,"Graphics3D") ) {
      snprintf(errmsg, MSG_SZ, "Input does not appear to be a Graphics3D file");
      return;
   }

   while( ptok != NULL ) {
//fprintf(stderr,"ptok = %s \n",ptok);
      if ( !strcmp(ptok,"FaceForm") ) {
         face_form = true;
      }
      else
      if ( !strcmp(ptok,"GrayLevel") || !strcmp(ptok,"RGBColor") || !strcmp(ptok,"Hue") || !strcmp(ptok,"CMYKColor") ) {
         if ( !strcmp(ptok,"GrayLevel") )
            color_case = 0;
         else
         if ( !strcmp(ptok,"RGBColor") )
            color_case = 1;
         else
         if ( !strcmp(ptok,"Hue") )
            color_case = 2;
         else
         if ( !strcmp(ptok,"CMYKColor") )
            color_case = 3;

         if ( color_case == 0 ) {
            ptok = strtok(NULL,parse_key);
            for(unsigned int j=0;j<3;j++) {
               col_read.push_back(ptok);
            }
         }
         else {
            for(unsigned int j=0;j<3;j++) {
               ptok = strtok(NULL,parse_key);
               col_read.push_back(ptok);
            }
         }

         if ( !poly_backface_ignore ) {
            current_col.read_decvals(col_read,errmsg);
            if(*errmsg)
               break;

            if ( color_case == 2 ) {
               vec3d hsv = current_col.get_vec3d();
               current_col.set_hsva(hsv[0], hsv[1], hsv[2]);
            }
            else
            if ( color_case == 3 ) {
               vec3d col = current_col.get_vec3d();
               col[0] = 1 - col[0];
               col[1] = 1 - col[1];
               col[2] = 1 - col[2];
               current_col = col;
            }
         }
         else {
            poly_backface_ignore = false;
            face_form = false;
         }

         if (face_form)
            poly_backface_ignore = true && disallow_back_faces;

         col_read.clear();
      }
      else
      if ( !strcmp(ptok,"Point") && (!strchr(hide_elems.c_str(), 'v')) ) {
         for(unsigned int j=0;j<3;j++) {
            ptok = strtok(NULL,parse_key);
            sscanf(ptok, "%lf", &coord[j]);
         }
         geom.add_vert(vec3d(coord[0], coord[1], coord[2]));
         if (current_col.is_set()) {
            // collect colored vertices in seperate geom here
            geom_cv.add_vert(vec3d(coord[0], coord[1], coord[2]));
            geom.set_v_col(verts_cv.size()-1, current_col);
            current_col = 0;
         }
      }
      else
      if ( (!strcmp(ptok,"Polygon") && !strchr(hide_elems.c_str(), 'f')) ||
           (!strcmp(ptok,"Line") && !strchr(hide_elems.c_str(), 'e')) ) {
         ptok = strtok(NULL,parse_key_polygon);
         int j = 0;
         while ( strcmp(ptok,"]") ) {
            // Sometimes find strange exponential format for very small numbers
//               if ( strstr( ptok, "*^-" ) )
            if ( strstr( ptok, "^" ) )
               coord[j++] = 0;
            else
//               sscanf(ptok, "%lf", &coord[j++]);
               coord[j++] = atof(ptok);
            if ( j > 2 ) {
               geom.add_vert(vec3d(coord[0], coord[1], coord[2]));
               face.push_back(verts.size()-1);
               j = 0;
            }
            ptok = strtok(NULL,parse_key_polygon);
         }
         geom.add_face(face);
         if (current_col.is_set()) {
            geom.set_f_col(faces.size()-1, current_col);
            current_col = 0;
         }

         face.clear();
      }
      else
      if ( !strcmp(ptok,"ViewPoint") ) {
         for(unsigned int j=0;j<3;j++) {
            ptok = strtok(NULL,parse_key);
            sscanf(ptok, "%lf", &coord[j]);
         }
         view_point = vec3d(coord[0], coord[1], coord[2]);
      }
      else
         if ( !strcmp(ptok,"ViewVertical") ) {
         for(unsigned int j=0;j<3;j++) {
            ptok = strtok(NULL,parse_key);
            sscanf(ptok, "%lf", &coord[j]);
         }
         view_vertical = vec3d(coord[0], coord[1], coord[2]);
      }
      else
      if ( !strcmp(ptok,"LightSources") ) {
         ptok = strtok(NULL,parse_key);
         while ( is_numeric(ptok) || !strcmp(ptok,"RGBColor") || !strcmp(ptok,"GrayLevel") ||
                  !strcmp(ptok,"Hue") || !strcmp(ptok,"CMYKColor") ) {
            if ( is_numeric(ptok) ) {
               sscanf(ptok, "%lf", &coord[0]);
               for(unsigned int j=1;j<3;j++) {
                  ptok = strtok(NULL,parse_key);
                  sscanf(ptok, "%lf", &coord[j]);
               }
               lights_geom.add_vert(vec3d(coord[0], coord[1], coord[2]));
            }
            else
            if ( !strcmp(ptok,"GrayLevel") || !strcmp(ptok,"RGBColor") || !strcmp(ptok,"Hue") || !strcmp(ptok,"CMYKColor") ) {
               if ( !strcmp(ptok,"GrayLevel") )
                  color_case = 0;
               else
               if ( !strcmp(ptok,"RGBColor") )
                  color_case = 1;
               else
               if ( !strcmp(ptok,"Hue") )
                  color_case = 2;
               else
               if ( !strcmp(ptok,"CMYKColor") )
                  color_case = 3;
               
               if ( color_case == 0 ) {
                  ptok = strtok(NULL,parse_key);
                  for(unsigned int j=0;j<3;j++) {
                     col_read.push_back(ptok);
                  }
               }
               else {
                  for(unsigned int j=0;j<3;j++) {
                     ptok = strtok(NULL,parse_key);
                     col_read.push_back(ptok);
                  }
               }

               current_col.read_decvals(col_read,errmsg);
               if(*errmsg)
                  break;

               if ( color_case == 2 ) {
                  vec3d hsv = current_col.get_vec3d();
                  current_col.set_hsva(hsv[0], hsv[1], hsv[2]);
               }
               else
               if ( color_case == 3 ) {
                  vec3d col = current_col.get_vec3d();
                  col[0] = 1 - col[0];
                  col[1] = 1 - col[1];
                  col[2] = 1 - col[2];
                  current_col = col;
               }

               lights_geom.set_v_col(lights_geom.verts().size()-1, current_col);
               col_read.clear();
            }
            ptok = strtok(NULL,parse_key_polygon);
         }
         delay_read = true;
      }

      if (!delay_read)
         ptok = strtok(NULL,parse_key);
      delay_read = false;
   }

   free(buffer);
}

void live3D_check_values(col_geom_v &lights_geom, vec3d &view_point, vec3d &view_vertical)
{
   if (!view_point.is_set())
      view_point = vec3d(1.3,-2.4,2);

   if (!view_vertical.is_set())
      view_vertical = vec3d(0,0,1);

   if (!lights_geom.verts().size()) {
      lights_geom.add_col_vert(vec3d(1,0,1), col_val(1.0,0.0,0.0));
      lights_geom.add_col_vert(vec3d(1,1,1), col_val(0.0,1.0,0.0));
      lights_geom.add_col_vert(vec3d(0,1,1), col_val(0.0,0.0,1.0));
      lights_geom.add_col_vert(vec3d(-1,-1,-1), col_val(1.0,1.0,1.0));
   }
}

void live3D_derotate(col_geom_v &geom, double &angle, vec3d &view_point)
{
   mat3d trans = mat3d::rot(0,0,-angle);
   trans = mat3d::rot(vec3d(0,0,1),view_point) * trans;
   geom.transform(trans);
}

void live3D_viewpoint(col_geom_v &geom, bool live3D_do_viewpoint, double &angle, vec3d &view_point, vec3d &view_vertical)
{
   mat3d trans = mat3d::rot(view_point,vec3d(0,0,1));
   vec3d rotated_view_vertical = trans * view_vertical;
   angle = atan2(rotated_view_vertical[0], rotated_view_vertical[1]);
   trans = mat3d::rot(0,0,angle) * trans;
   if (live3D_do_viewpoint)
      geom.transform(trans);
}

void live3D_dump_lights_geom(col_geom_v &lights_geom, string lights_geom_file, char *errmsg)
{
   if (lights_geom_file.find(".off",0) == string::npos)
      lights_geom_file += ".off";

   FILE *ofile = fopen(lights_geom_file.c_str(), "w");
   if(!lights_geom.write(lights_geom_file, errmsg)) {
      if(*errmsg)
         snprintf(errmsg, MSG_SZ, "could not open output file for color table \'%s\'", lights_geom_file.c_str());
      return;
   }
   fclose(ofile);
}

int main(int argc, char *argv[])
{
   m2off_opts opts;
   opts.process_command_line(argc, argv);
   
   vec3d view_point;
   vec3d view_vertical;
   double angle;
   col_geom_v lights_geom;

   char errmsg[MSG_SZ];
   // read mtxt into one long string
   string m_txt = read_file_to_str(opts.ifile, errmsg);
   if(*errmsg)
      opts.error(errmsg);

   // put colored vertices in seperate geom so the can be prepended to the list
   col_geom_v geom, geom_cv;
   m_parse(m_txt, geom, geom_cv, opts.hide_elems, opts.disallow_back_faces, 
           lights_geom, view_point, view_vertical, errmsg);
   if(*errmsg)
      opts.error(errmsg);
   m_txt.clear();
      
   geom_cv.append(geom);
   geom = geom_cv;
   geom_cv.clear_all();

   live3D_check_values(lights_geom, view_point, view_vertical);
   live3D_viewpoint(geom, opts.live3D_do_viewpoint, angle, view_point, view_vertical);
   // only if -v wasn't specified "de-rotate" lights
   if (!opts.live3D_do_viewpoint)
      live3D_derotate(lights_geom, angle, view_point);

   if (opts.lights_geom_file.size())
      live3D_dump_lights_geom(lights_geom, opts.lights_geom_file, errmsg);
   if(*errmsg)
      opts.error(errmsg);

   // sort/merge all and orient faces
   sort_merge_elems(geom, "vef", pow(10, -opts.sig_compare));
   geom.orient();

   if(!geom.write(opts.ofile, errmsg))
      opts.error(errmsg);

   return 0;
}
