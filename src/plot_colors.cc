/*
   Copyright (c) 2010, Roger Kaufman

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
   Name: plot_colors.cc
   Description: Creates models of colors. Color wheel, RGB, HSV, HSL
   Project: Antiprism - http://www.antiprism.com
*/

#include <stdio.h>
#include <stdlib.h>

#include <ctype.h>
#include <unistd.h>

#include <string>
#include <vector>
#include <algorithm>

#include "../base/antiprism.h"

using std::string;
using std::vector;
using std::min;
using std::max;


class plot_colors_opts: public prog_opts {
   public:
      string ifile;
      string ofile;
      string gfile;
      string mfile;
      
      int display_type;
      int color_system_mode;
      int show_container;
      bool unique_colors;
      int upright_view;
      int chroma_level;
      int hsl_height;
      bool plot_centroid;
      double sat_threshold;
      double value_power;
      double value_advance;
      bool cmy_mode;
      bool ryb_mode;
      bool seven_mode;
      bool antiprism_map;
      int map_maximum;
      bool verbose;

      vector<double> sat_powers;
      vector<double> value_powers;
      
      vector<col_val> cols;
      coloring clrngs[3];

      plot_colors_opts(): prog_opts("plot_colors"),
                     display_type(2),
                     color_system_mode(3),
                     show_container(2),
                     unique_colors(false),
                     upright_view(0),
                     chroma_level(0),
                     hsl_height(0),
                     plot_centroid(false),
                     sat_threshold(1.0),
                     value_advance(0.0),
                     cmy_mode(false),
                     ryb_mode(false),
                     seven_mode(false),
                     antiprism_map(false),
                     map_maximum(256),
                     verbose(false)
                     {}

      void process_command_line(int argc, char **argv);
      void usage();
};

void plot_colors_opts::usage()
{
   fprintf(stderr,
"\n"
"Usage: %s [options]\n"
"\n"
"various plots of colors and blendings\n"
"\n"
"Options\n"
"%s"
"  -f <file> also output the colors to a simple RGB map file\n"
"  -g        output map file in -f as antiprism map\n"
"  -v        verbose\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\nScene Options\n"
"  -d <int>  display type. 1 - color wheel  2 - plot (default)\n"
"  -m <mode> color system mode. 1 - HSV  2 - HSL  3 - RGB (default)\n"
"               if display is color wheel and mode is mode is RGB, then outer\n"
"               ring in RGB and inner rings are HSV\n"
"  -r <int>  HSV/HSL chroma  0 - none (cylinder)  1 - conic  2 - hexagonal\n"
"  -S        HSV/HSL distribute colors heptagonally, 7 ways, as in the rainbow\n"
"  -k <int>  display type 2: container control\n"
"               0 - suppress  1 - no facets on chroma-2  2 - full (default)\n"
"  -V <int>  display type 2: 1 - show upright  2 - don\'t rotate cube\n"
"  -q <int>  height of the HSL cones  1 - double  2 - make 90 degree dihedral\n"
"\nColor Blending Options\n"
"  -b        display type 2: plot color centroid(s)\n"
"  -s <sat>  HSV/HSL saturation curve. Greather than 0  (default: 1)\n"
"               1.0 - no curve. lower than 1.0 makes blends more pastel\n"
"               4 numbers can be entered seperated by commas\n"
"  -t <val>  HSV/HSL threshold to use average saturation. (default: 1)\n"
"               between 0.0 (all averaging) and 1.0 (no averaging)\n"
"  -l <val>  HSV/HSL value curve (default: 0)\n"
"               simulates subtractive coloring for blending 3 or more colors\n"
"               RGB: Red+Green+Blue = White   Cyan+Magenta+Yellow = Black\n"
"               RYB: Red+Yellow+Blue = Black  Green+Magenta+Orange = White\n"
"               1.0 - no curve. lower than 1.0 number makes blends lighter\n"
"               0.0 - use average value instead\n"
"               4 numbers can be entered seperated by commas\n"
"  -u <val>  HSV/HSL value advance. Rotates meaning of white and black\n"
"               valid values 0.0 to 120.0 degrees (default: 0)\n"
"  -y        RYB mode. Blend colors as in Red-Yellow-Blue color wheel\n"
"  -c        CMY mode. Complimentary colors.  RGB->(RYB/GMO)->CMY->blend\n"
"\nColor Source Options\n"
"\n"
"A colour value a value in form 'R,G,B,A' (3 or 4 values 0.0-1.0, or 0-255)\n"
"or hex 'xFFFFFF'\n"
"HSV color format may also be used in the form of \'(h or H),S,V,A\'\n"
"  note that the prefix of \'h\' or \'H\' must be used\n"
"  h is any angle. e.g h30.0  or  H is any real number. e.g H0.0833\n"
"  S (saturation), V (value) and A (alpha) are optional real\n"
"  numbers from 0.0 to 1.0 (default: S,V,A = 1.0,1.0,1.0)\n"
"\n"
"note: -M -O and -C may all be used together\n"
"\n"
"  -M <map>  get colors from a color map, or multiple maps seperated by commas\n"
"  -O <file> get colors from an OFF file\n"
"  -C <elms> one or more color values separated by \":\"\n"
"               note: 'invisible', 'none' or color indexes are not accepted\n"
"  -U        allow only unique colors (sorts by hue)\n"
"  -Z <int>  maximum entries to read from open ended maps (default: 256)\n"
"\n"
"\n",prog_name(), help_ver_text);
}

void plot_colors_opts::process_command_line(int argc, char **argv)
{
   extern char *optarg;
   extern int optind, opterr, optopt;
   opterr = 0;
   char c;
   char errmsg[MSG_SZ];
   
   vector<char *> color_strings;
   vector<double> double_parms;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hd:m:k:q:V:f:gvr:s:t:u:l:bcySl:M:O:C:UZ:o:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         case 'd':
            if(!read_int(optarg, &display_type, errmsg))
               error(errmsg, c);
            if(display_type < 1 || display_type > 2)
               error("display type must be between 1 and 2", c);
            break;
            
         case 'm':
            if(!read_int(optarg, &color_system_mode, errmsg))
               error(errmsg, c);
            if(color_system_mode < 1 || color_system_mode > 3)
               error("color system mode must be between 1 and 3", c);
            break;
            
         case 'k':
            if(!read_int(optarg, &show_container, errmsg))
               error(errmsg, c);
            if(show_container < 0 || show_container > 2)
               error("container mode must be between 0 and 2", c);
            break;
            
         case 'q':
            if(!read_int(optarg, &hsl_height, errmsg))
               error(errmsg, c);
            if(hsl_height < 1 || hsl_height > 2)
               error("hsl height must be between 1 and 2", c);
            break;
            
         case 'V':
            if(!read_int(optarg, &upright_view, errmsg))
               error(errmsg, c);
            if(upright_view < 1 || upright_view > 2)
               error("view type must be between 1 and 2", c);
            break;
            
         case 'f':
            mfile = optarg;
            break;
            
         case 'g':
            antiprism_map = true;
            break;
            
         case 'v':
            verbose = true;
            break;
            
         case 'r':
            if(!read_int(optarg, &chroma_level, errmsg))
               error(errmsg, c);
            if(chroma_level < 0 || chroma_level > 3)
               error("display type must be between 0 and 3", c);
            break;
            
         case 's':
            if(!read_double_list(optarg, double_parms, errmsg, 4))
               error(errmsg, c);
            for (unsigned int i=0;i<double_parms.size();i++) {
               if(double_parms[i] <= 0.0)
                  error("color centroid saturation curve must be greater than zero", c);
               sat_powers.push_back(double_parms[i]);
            }
            break;
            
         case 't':
            if(!read_double(optarg, &sat_threshold, errmsg))
               error(errmsg, c);
            if(sat_threshold < 0.0 || sat_threshold > 1.0)
               error("HSV/HSL threshold must be between 0 and 1", c);
            break;
            
         case 'u':
            if(!read_double(optarg, &value_advance, errmsg))
               error(errmsg, c);
            if(value_advance < 0.0 || value_advance > 120.0)
               error("HSV/HSL value advance must be between 0 and 120", c);
            break;
            
         case 'l':
            if(!read_double_list(optarg, double_parms, errmsg, 4))
               error(errmsg, c);
            for (unsigned int i=0;i<double_parms.size();i++) {
               if(double_parms[i] < 0.0)
                  error("value curve must be greater than or equal to zero", c);
               value_powers.push_back(double_parms[i]);
            }
            break;
            
         case 'b':
            plot_centroid = true;
            break;
            
         case 'c':
            cmy_mode = true;
            break;
            
         case 'y':
            ryb_mode = true;
            break;
            
         case 'S':
            seven_mode = true;
            break;
            
         case 'M':
            if(!read_colorings(clrngs, optarg, errmsg))
               error(errmsg, c);
            break;
            
         case 'O':
            gfile = optarg;
            break;
            
         case 'C':
            split_line(optarg, color_strings, ":");
            
            for (unsigned int i=0;i<color_strings.size();i++) {
//fprintf(stderr,"color_strings[%d] = %s\n",i,color_strings[i]);
               col_val col;
               if (!col.read(color_strings[i], errmsg))
                  error(errmsg, c);
               if(!col.is_set() || !col.is_val() || col.is_inv()) {
                  snprintf(errmsg, MSG_SZ, "color value not accepted here: %s\n",color_strings[i]);
                  error(errmsg, c);
               }
               cols.push_back(col);
            }
            break;
            
         case 'U':
            unique_colors = true;
            break;
            
         case 'Z':
            if(!read_int(optarg, &map_maximum, errmsg))
               error(errmsg, c);
            if(map_maximum < 0)
               error("maximum map elements to read in must be greater than 0", c);
            break;
            
         case 'o':
            ofile = optarg;
            break;

         case '?':
            error("unknown option", string("-")+(char)optopt);

         case ':':
            error("missing argument", string("-")+(char)optopt);

         default:
            error("unknown command line error");
      }
   }

   if(argc-optind > 0)
      error("too many arguments");

   // fill in missing sat_powers with 1.0, meaning use centroid saturation
   for (int i=sat_powers.size();i<4;i++)
      sat_powers.push_back(1.0);
      
   // fill in missing value_powers with -1.0, meaning use average values
   for (int i=value_powers.size();i<4;i++)
      value_powers.push_back(-1.0);
      
   if (antiprism_map && !mfile.length())
      error("not valid without -f","g");
      
   if (hsl_height && color_system_mode != 2)
      warning("HSL height adjustment has no effect in RGB or HSV mode","q");

   if (show_container == 1 && (color_system_mode == 3 || chroma_level < 2))
      warning("facets can only be shown in HSV or HSL chroma 2","k");
      
   if (upright_view == 2 && color_system_mode != 3)
      warning("view 2 has no effect when not in RGB mode","V");
      
   if (color_system_mode == 3 && chroma_level)
      warning("chroma has no effect in RGB mode","r");
      
   if (color_system_mode == 3 && seven_mode)
      warning("heptagonal view has no effect in RGB mode","S");
}

// http://en.wikipedia.org/wiki/HSL_and_HSV#Hue_and_chroma
// hue and chroma output are from 0 to 1
void get_chroma(const col_val &col, int chroma_level, double &hue, double &chroma)
{
   int R = col[0];
   int G = col[1];
   int B = col[2];
   
   int M = max(R,max(G,B));
   int m = min(R,min(G,B));
   int C = M - m;
   
   double H = 0.0;
   if (!C)
      H = 0.0;
   else
   if (M == R) // red
      H = fmod(((double)(G-B)/C),6);  // (g-b/c) mod 6;
   else
   if (M == G) // green
      H = ((double)(B-R)/C)+2;        // (b-r/c)+2;
   else
   if (M == B) // blue
      H = ((double)(R-G)/C)+4;        // (r-g/c)+4;
   H *= 60.0;
   if (H<0) H+=360.0;
   
   if (chroma_level == 1) {
      hue = H/360.0;
      chroma = C/255.0;
   }
   else
   if (chroma_level == 2) {
      double a = (2*R-G-B)*0.5;
      double b = sqrt(3)/2.0*(G-B);
   
      hue = atan2(b,a)/(2*M_PI);
      if (hue < 0.0)
         hue += 1.0;
      chroma = sqrt(a*a+b*b)/255.0;
   }
}

/*
RGB->RYB
0 >= 60 degrees, multiply by 2
60 >= 120 degrees, add 60
120 >= 180 maps to 180 to 210
180 >= 240 maps to 210 to 240

So, compression happens
at 120 it is adding 60
at 180 it is adding 30
at 240 it is adding 0

Then
120 >= 240 degrees, N+((240-N)/2)


RYB->RGB
0 >= 120 degrees, divide by 2
120 >= 180 degrees, subtract 60
180 maps to 120
210 maps to 180
240 maps to itself

So, decompression happens
at 180 it is subtracting 60
at 210 it is subtracting 30
at 240 it is subtracting 0

Then
180 >= 240 degrees, N-(240-N)
*/

// angle represented by 0 to 360 degrees
// input: HSV/HSL angle
// output: angle adjusted for RYB mode
double hsx_to_ryb(double angle)
{
   if (angle > 0.0 && angle <= 60.0)
      angle *= 2.0;
   else
   if (angle > 60.0 && angle <= 120.0)
      angle += 60.0;
   else
   if (angle > 120.0 && angle <= 240.0)
      angle += (240.0-angle)/2;
   return angle;
}

// angle represented by 0 to 360 degrees
// input: angle adjusted for RYB mode
// output: HSV/HSL angle
double ryb_to_hsx(double angle)
{
   if (angle > 0.0 && angle <= 120.0)
      angle /= 2.0;
   else
   if (angle > 120.0 && angle <= 180.0)
      angle -= 60.0;
   else
   if (angle > 180.0 && angle <= 240.0)
      angle -= (240.0-angle);
   return angle;
}

col_val rgb_simple_compliment(const col_val &col)
{
   int rc = 255 - col[0];
   int gc = 255 - col[1];
   int bc = 255 - col[2];
   return(col_val(rc,gc,bc,col[3]));
}

col_val rgb_compliment(const col_val &col, bool ryb_mode)
{
   // only need hue so algorithm doesn't matter
   vec4d hsxa = col.get_hsva();
   double angle = rad2deg(2*M_PI*hsxa[0]);
 
   if (ryb_mode)
      angle = hsx_to_ryb(angle);

   angle += 180.0;
   if (angle >= 360.0)
      angle -= 360.0;
      
   if (ryb_mode)
      angle = ryb_to_hsx(angle);
   angle /= 360.0;
   
   col_val rcol;
   rcol.set_hsva(angle,hsxa[1],hsxa[2],hsxa[3]);
   return rcol;
}

// wrapper for multiple HSV/HSL algorithms
vec4d get_hsxa(const col_val &col, int color_system_mode)
{
   return (color_system_mode == 1) ? col.get_hsva() : col.get_hsla();
}

col_val set_hsxa(double hue, double sat, double val, double alpha, int color_system_mode)
{
   col_val col;
   if (color_system_mode == 1)
      col.set_hsva(hue,sat,val,alpha);
   else
   if (color_system_mode == 2)
      col.set_hsla(hue,sat,val,alpha);
   return col;
}

// core code furnished by Adrian Rossiter
col_val blend_HSX_centroid(const vector<col_val> &cols, int color_system_mode, double sat_power, double sat_threshold, double value_power, double value_advance, bool ryb_mode, bool verbose)
{
   // safety conditions
   if (!cols.size())
      return col_val();
   else
   if (cols.size() == 1)
      return cols[0];
   
   // can't blend two or less colors to black
   if (cols.size() < 3)
      value_power = 0.0;
   
   double saturation_sum = 0.0;
 
   vec4d sum(0.0, 0.0, 0.0, 0.0);
   for(unsigned int i=0; i<cols.size(); i++) {
      vec4d hsxa = get_hsxa(cols[i], color_system_mode);
         
      double S = pow(hsxa[1], sat_power); // map onto distorted disc
      double angle = 2*M_PI*hsxa[0];
      
      // RYB mode
      if (ryb_mode)
         angle = deg2rad(hsx_to_ryb(rad2deg(angle)));

      // if value_power is set, simulate subtractive coloring for 3 or more colors
      double V = (value_power <= 0.0) ? hsxa[2] : pow(fabs(60.0-fmod(rad2deg(angle)+(ryb_mode ? 60.0 : 0.0)+value_advance,120.0))/60, value_power);

      sum += vec4d(S*cos(angle), S*sin(angle), V, hsxa[3]); // point in cylinder

      if (sat_threshold < 1.0)
         saturation_sum+=hsxa[1];
   }
   // average
   sum /= cols.size();

   double H = 0.0;
   double S = 0.0;
   
   // saturation
   S = pow(sum[0]*sum[0]+sum[1]*sum[1], 0.5/sat_power); // map back from distorted disc
   
   // saturations less than 1/255 will happen due to inaccuracy of HSx->RGB conversion
   // note that 255,255,254 has a saturation of 1/255 = 0.00392157..., which is the smallest valid saturation
   if (S<1/255.0)
      S = 0.0;

   if (verbose)
      fprintf(stderr,"blend_HSX_centroid: pre-threshold saturation = %g\n",S);

   // saturation of color centroid is higher than sat_threshold, use average saturation         
   if (S > sat_threshold)
      S = saturation_sum/cols.size();
   
   // hue
   // if saturation is 0, no need to calculate hue
   //if (!double_equality(sum[0],0.0,epsilon) && !double_equality(sum[1],0.0,epsilon)) { // old error, made nice pattern, but wrong
   if (S != 0.0) {
      H = atan2(sum[1], sum[0])/(2*M_PI);
      if(H<0)
         H += 1.0;

      // RYB mode
      if (ryb_mode)
         H = deg2rad(ryb_to_hsx(rad2deg(H*2*M_PI)))/(2*M_PI);
   }

   if (verbose)
      fprintf(stderr, "blend_HSX_centroid: angle = %g (hue = %g, saturation = %g)\n", rad2deg(H*2*M_PI), H, S);
      
   col_val col = set_hsxa(H, S, sum[2], sum[3], color_system_mode);
      
   if (verbose)
      fprintf(stderr,"blend_HSX_centroid: RGBA result = %d %d %d %d\n\n",col[0],col[1],col[2],col[3]);
   return col;
}

col_val blend_RGB_centroid(const vector<col_val> &cols, bool ryb_mode)
{
   vec4d col(0.0, 0.0, 0.0, 0.0);
   for(unsigned int i=0; i<cols.size(); i++) {
      col_val rcol = cols[i];
      if (ryb_mode) {
         // only need hue so algorithm doesn't matter
         vec4d hsxa = rcol.get_hsva();
         hsxa[0] = hsx_to_ryb(rad2deg(2*M_PI*hsxa[0]));
         rcol.set_hsva(hsxa[0]/360.0,hsxa[1],hsxa[2],hsxa[3]);
      }
      col += rcol.get_vec4d();
   }
   col /= cols.size();
   return col;
}

// furnished by Adrian Rossiter
void color_wheel(col_geom_v &geom, const vector<col_val> &cols, int color_system_mode,
                 vector<double> &sat_powers, double sat_threshold, vector<double> &value_powers, double value_advance, bool ryb_mode, bool verbose)
{
   // RK - dynamic polygon size
   unsigned int sz = cols.size();
   int polygon_sz = (sz < 120) ? 120 : sz;

   int num_points = polygon_sz/sz;
   int face_sz = 2*(num_points+1);
   double sect_ang = 2*M_PI/sz;
   double ang_inc = sect_ang/num_points;
   for(unsigned int i=0; i<sz; i++) {
      vector<int> face(face_sz);
      for(int j=0; j<num_points+1; j++) {
         double ang = i*sect_ang + j*ang_inc;
         geom.add_vert(vec3d(cos(ang), sin(ang), 0.0));
         face[j] = geom.verts().size()-1;
         geom.add_vert(vec3d(cos(ang)/2, sin(ang)/2, 0.0));
         face[face_sz-1-j] = geom.verts().size()-1;
      }
      geom.add_col_face(face, cols[i]);
   }
   
   // RK - user gets multiple level control
   double sat_power = sat_powers[0];
   double value_power = value_powers[0];
   
   for(int lvl=0; lvl<4; lvl++) {
      vector<int> face(polygon_sz);
      for(int i=0; i<polygon_sz; i++) {
         double ang = 2*M_PI*i/polygon_sz;
         double rad = 0.5*(4-lvl)/4;
         geom.add_vert(vec3d(rad*cos(ang), rad*sin(ang), 0.01*lvl));
         face[i] = geom.verts().size()-1;
      }
      if (color_system_mode == 3)
         // RK - if RGB mode, only show that blend in all four levels
         geom.add_col_face(face, blend_RGB_centroid(cols, ryb_mode));
      else {
         // RK - only change the powers if there is a new one not the default
         // keeps the center of the bullseye consistent with the last valid power
         if (sat_powers[lvl] != 1.0)
            sat_power = sat_powers[lvl];
         if (value_powers[lvl] > -1.0)
            value_power = value_powers[lvl];
               
         geom.add_col_face(face, blend_HSX_centroid(cols, color_system_mode, sat_power, sat_threshold, value_power, value_advance, ryb_mode, verbose));
      }
   }
   
   
   coloring clrng(&geom);
   geom.add_missing_impl_edges();
   clrng.e_one_col(col_val::invisible);
   clrng.v_one_col(col_val::invisible);
}

geom_v make_unit_circle(int polygon_size)
{
   geom_v geom;
   double arc = deg2rad(360.0/(double)polygon_size);

   double angle = 0.0;
   for (int i=0;i<polygon_size;i++) {
      geom.add_vert(vec3d(cos(angle), sin(angle), 0.0));
      geom.add_edge(make_edge(i, (i+1)%polygon_size));
      angle += arc;
   }
   
   return(geom);
}

void make_chroma2_container(col_geom_v &geom, int polygon_size, int color_system_mode, int show_container)
{
   geom.append(make_unit_circle(polygon_size));
   
   // if show_container is 1 do not show facets
   if (show_container < 2)
      return;
      
   int new_vert = 0;
   if (color_system_mode == 2) {
      geom.add_vert(vec3d(0.0,0.0,0.5));
      new_vert = geom.verts().size()-1;
      for (int i=0;i<polygon_size;i++)
         geom.add_edge(make_edge(i, new_vert));
   }
   geom.add_vert(vec3d(0.0,0.0,(color_system_mode == 2) ? -0.5 : -1.0));
   new_vert = geom.verts().size()-1;
   for (int i=0;i<polygon_size;i++)
      geom.add_edge(make_edge(i, new_vert));
}

col_geom_v make_hsx_container(int chroma_level, int color_system_mode, bool seven_mode, int show_container)
{
   col_geom_v geom;
   if (chroma_level == 2) {
      make_chroma2_container(geom, (seven_mode ? 7 : 6), color_system_mode, show_container);
      geom.transform(mat3d::transl(color_system_mode == 2 ? vec3d(0.0,0.0,0.5) : vec3d(0.0,0.0,1.0)));
   }
   else
   {
      geom = make_unit_circle(60);
      geom.transform(mat3d::transl((color_system_mode == 1 || !chroma_level) ? vec3d(0.0,0.0,1.0) : vec3d(0.0,0.0,0.5)));
      if (!chroma_level)
         geom.append(make_unit_circle(60));
   }
   
   coloring clrng(&geom);
   clrng.v_one_col(col_val::invisible);
   clrng.e_one_col(col_val(1.0,1.0,1.0,0.1));
   
   return(geom);
}

// angle represented by 0 to 360 degrees
// input: HSV/HSL angle
// output: angle adjusted for heptagon
double hsx_to_7gon(double angle)
{
//fprintf(stderr,"angle in %g\n",angle);
   double heptagon_angle = 360.0/7.0; // 51.42857143...
   
   // take 0 to 60 range to 0 to 102.8571429...
   if (angle > 0.0 && angle <= 60.0)
      angle *= heptagon_angle/30.0;
   else
   if (angle > 60.0)
      angle = ((angle-60.0) * heptagon_angle/60.0) + heptagon_angle*2;

//fprintf(stderr,"angle out %g\n",angle);
   return angle;
}

// code to draw cone into heptagon furnished by Adrian Rossiter
double saturation_correction_7gon(double hue, double sat)
{
   double angle = (2*M_PI/7)*(0.5-fmod(7.0*hue, 1.0));
   double scale = cos(M_PI/7)/cos(angle);
   return (sat * scale);
}

void plot_hsx_point(col_geom_v &geom, const col_val &col, int color_system_mode, int chroma_level, bool ryb_mode, bool seven_mode)
{
   if (!col.is_set())
      return;

   vec4d hsxa = get_hsxa(col, color_system_mode);

   // heptagonal display overrides RYB distribution
   bool seven_mode_chroma_level2 = false;
   if (seven_mode) {
      ryb_mode = false;
      if (chroma_level == 2) {
         chroma_level = 1;
         seven_mode_chroma_level2 = true;
      }
   }
      
   double H = (ryb_mode) ? hsx_to_ryb(rad2deg(2*M_PI*hsxa[0]))/360.0 : hsxa[0];
   double S = hsxa[1];

   if (chroma_level) {
      col_val ccol = col;
      if (ryb_mode)
         ccol = set_hsxa(H, S, hsxa[2], hsxa[3], color_system_mode);
      get_chroma(ccol,chroma_level,H,S);
   }
   
   if (seven_mode) {
      // distribute hex oriented hue into 7-way mode
      double angle = rad2deg(2*M_PI*H);
      H = hsx_to_7gon(angle)/360.0;
      // draw cone into heptagon
      if (seven_mode_chroma_level2)
         S = saturation_correction_7gon(H,S);
   }

   double angle = 2*M_PI*H;

   geom.add_col_vert(vec3d(S*cos(angle), S*sin(angle), hsxa[2]), col);
}

void plot_hsx(col_geom_v &geom, const vector<col_val> &cols, int color_system_mode, int chroma_level, bool ryb_mode, bool seven_mode)
{
   for(unsigned int i=0; i<cols.size(); i++)
      plot_hsx_point(geom, cols[i], color_system_mode, chroma_level, ryb_mode, seven_mode);
}

col_geom_v make_cube()
{
   col_geom_v geom;
      
   geom.add_vert(vec3d(1, 1, 1)); // 0
   geom.add_vert(vec3d(1, 1, 0)); // 1
   geom.add_vert(vec3d(1, 0, 1)); // 2
   geom.add_vert(vec3d(1, 0, 0)); // 3
   geom.add_vert(vec3d(0, 1, 1)); // 4
   geom.add_vert(vec3d(0, 1, 0)); // 5
   geom.add_vert(vec3d(0, 0, 1)); // 6
   geom.add_vert(vec3d(0, 0, 0)); // 7
   
   geom.add_edge(make_edge(0,1));
   geom.add_edge(make_edge(0,2));
   geom.add_edge(make_edge(0,4));
   geom.add_edge(make_edge(1,3));
   geom.add_edge(make_edge(2,3));
   geom.add_edge(make_edge(2,6));
   geom.add_edge(make_edge(4,6));
   geom.add_edge(make_edge(1,5));
   geom.add_edge(make_edge(4,5));
   geom.add_edge(make_edge(7,3));
   geom.add_edge(make_edge(7,5));
   geom.add_edge(make_edge(7,6));
   
   coloring clrng(&geom);
   clrng.v_one_col(col_val::invisible);
   clrng.e_one_col(col_val(1.0,1.0,1.0,0.1));
   
   return geom;
}

void plot_rgb_point(col_geom_v &geom, const col_val &col, bool ryb_mode)
{
   if (!col.is_set())
      return;
  
   col_val rcol = col;
   if (ryb_mode) {
      // only need hue so algorithm doesn't matter
      vec4d hsxa = rcol.get_hsva();
      hsxa[0] = hsx_to_ryb(rad2deg(2*M_PI*hsxa[0]));
      rcol.set_hsva(hsxa[0]/360.0,hsxa[1],hsxa[2],hsxa[3]);
   }

   geom.add_col_vert(vec3d(rcol[0]/255.0,rcol[1]/255.0,rcol[2]/255.0),col);
}

void plot_rgb_cube(col_geom_v &geom, const vector<col_val> &cols, bool ryb_mode)
{
   for(unsigned int i=0; i<cols.size(); i++)
      plot_rgb_point(geom, cols[i], ryb_mode);
}

bool cmp_hue(const col_val &a, const col_val &b, double eps)
{
   return (a.get_hsva()[0] < b.get_hsva()[0]);
}

class hue_cmp
{
public:
   double eps;
   hue_cmp(double ep): eps(ep) {}
   bool operator() (const col_val &a, const col_val &b) { return cmp_hue(a, b, eps); }
};

void collect_col(vector<col_val> &cols, const col_val &col)
{
   if (!col.is_set())
      return;
   cols.push_back(col);
}

void collect_cols_from_geom(const col_geom_v &geom, vector<col_val> &cols)
{
   for(unsigned int i=0; i<geom.verts().size(); i++)
      collect_col(cols,geom.get_v_col(i));
   for(unsigned int i=0; i<geom.edges().size(); i++)
      collect_col(cols,geom.get_e_col(i));
   for(unsigned int i=0; i<geom.faces().size(); i++)
      collect_col(cols,geom.get_f_col(i));
}

void collect_cols(vector<col_val> &cols, plot_colors_opts &opts)
{
  // maps
  if (opts.clrngs[0].get_cmaps().size()) {
      vector<color_map *> maps = opts.clrngs[0].get_cmaps();
      for(unsigned int i=0; i<maps.size(); i++) {
         unsigned int map_sz = maps[i]->effective_size();
//fprintf(stderr,"map size = %u\n",map_sz);
         bool open_ended_map = (map_sz >= INT_MAX);
         unsigned int max_map_sz = (open_ended_map) ? (unsigned int)opts.map_maximum : map_sz;
//fprintf(stderr,"max map size = %u\n",max_map_sz);
         if (open_ended_map) {
            char errmsg[MSG_SZ];
            snprintf(errmsg, MSG_SZ, "map entry %d: only %d out of %u map entries read in\n",i+1,opts.map_maximum,map_sz);
            opts.warning(errmsg, 'Z');
         }
         for(unsigned int j=0; j<max_map_sz; j++)
            collect_col(cols,maps[i]->get_col(j));
      }
   }
   
   // file
   if(opts.gfile.length()) {
      col_geom_v geom;
      char errmsg[MSG_SZ];
      if(!geom.read(opts.gfile, errmsg))
         opts.error(errmsg);
      if(*errmsg)
         opts.warning(errmsg);
      
      collect_cols_from_geom(geom, cols);
   }
   
   // command line
   for(unsigned int i=0; i<opts.cols.size(); i++)
      collect_col(cols,opts.cols[i]);

   if (opts.unique_colors) {
      sort( cols.begin(), cols.end(), hue_cmp(epsilon) );    
      vector<col_val>::iterator ci = unique(cols.begin(), cols.end());
      cols.resize( ci - cols.begin() );
   }

/* debug. output hsv values  
   for(unsigned int i=0; i<cols.size(); i++) {
      vec4d hcol;
      hcol = cols[i].get_hsva();
      fprintf(stderr,"%g %g %g %g\n",hcol[0],hcol[1],hcol[2],hcol[3]);
   }
*/
   
   if (opts.cmy_mode)
      for(unsigned int i=0; i<cols.size(); i++)
         cols[i] = rgb_compliment(cols[i], opts.ryb_mode);
         
   for(unsigned int i=0; i<cols.size(); i++) {
      if (cols[i][3] < 32) {
         opts.warning("colors with a very small alpha values may not be seen");
         break;
      }
   }
}

int main(int argc, char *argv[])
{
   plot_colors_opts opts;
   opts.process_command_line(argc, argv);
   
   col_geom_v geom;
   char errmsg[MSG_SZ];
      
   vector<col_val> cols;
   collect_cols(cols, opts);
   
   if (!cols.size())
      opts.error("found no colors to plot");
      
   if (cols.size() < 3) {
      for(unsigned int i=0; i<4; i++) {
         if (opts.value_powers[i] > -1.0) {
            opts.warning("subtractive colors only works on three colors or more","v");
            break;
         }
      }
   }
             
   if (opts.display_type == 1)
      color_wheel(geom, cols, opts.color_system_mode, opts.sat_powers, opts.sat_threshold, opts.value_powers, opts.value_advance, opts.ryb_mode, opts.verbose);
   else
   if (opts.display_type == 2) {
      if (opts.color_system_mode == 1 || opts.color_system_mode == 2) {
         plot_hsx(geom, cols, opts.color_system_mode, opts.chroma_level, opts.ryb_mode, opts.seven_mode);
         
         if (opts.plot_centroid) {
            for(unsigned int i=0; i<4; i++) {
               if (!i || opts.sat_powers[i] != 1.0 || opts.value_powers[i] > -1.0) {
                  col_val col = blend_HSX_centroid(cols, opts.color_system_mode, opts.sat_powers[i], opts.sat_threshold, opts.value_powers[i], opts.value_advance, opts.ryb_mode, opts.verbose);
                  plot_hsx_point(geom, col, opts.color_system_mode, opts.chroma_level, opts.ryb_mode, opts.seven_mode);
               }
            }
         }
      }
      else
      if (opts.color_system_mode == 3) {
         plot_rgb_cube(geom, cols, opts.ryb_mode);
         
         if (opts.plot_centroid) {
            col_val col = blend_RGB_centroid(cols, opts.ryb_mode);
            plot_rgb_point(geom, col, opts.ryb_mode);
         }
      }
   }
   
   if (opts.display_type == 2) {
      // container
      if (opts.show_container) {
         if (opts.color_system_mode == 1 || opts.color_system_mode == 2)
            geom.append(make_hsx_container(opts.chroma_level, opts.color_system_mode, opts.seven_mode, opts.show_container));
         else
         if (opts.color_system_mode == 3)
            geom.append(make_cube());
      }
   
      // view for cube
      if (opts.color_system_mode == 3) {
         if (opts.upright_view != 2) {
            geom.transform(mat3d::rot(deg2rad(45),0,0));
            geom.transform(mat3d::rot(0,-asin(1/sqrt(3)),0));
         }
      }
   
      // HSL height, before turning upright
      if (opts.hsl_height && opts.color_system_mode == 2) {
         geom.transform(mat3d::scale(1.0,1.0,2.0));
         // if option to make 90 degree dihedral permeters
         if (opts.hsl_height == 2 && opts.chroma_level) {
            double scale = 1/(2*tan(M_PI*(1.0 / (opts.seven_mode ? 7.0 : 6.0))));
            // scale of heptagon to 1 is .8677674782351162 = cos(Pi*5/14)+cos(Pi*5/14) or cos(Pi*3/14)/sin(Pi*5/14)
            if (opts.seven_mode)
               scale *= cos(M_PI*5.0/14.0)*2.0;
//fprintf(stderr,"scale = %g\n",scale);
            geom.transform(mat3d::scale(1.0,1.0,fabs(scale)));
         }
      }
      
      // view upright for all models if option 1    
      if (opts.upright_view == 1)
         geom.transform(mat3d::rot(deg2rad(-90),0,0));
   }

   if(!geom.write(opts.ofile, errmsg))
      opts.error(errmsg);
      
   // output map colors as an RGB file
   if(opts.mfile.length()) {
      FILE *ofile = fopen(opts.mfile.c_str(), "w");
      if(!ofile)
         opts.error("could not open output file \'"+opts.ofile+"\'");
      for(unsigned int i=0; i<cols.size(); i++) {
         //fprintf(stderr,"%3d %3d %3d %3d\n",cols[i][0],cols[i][1],cols[i][2],cols[i][3]);
         if (opts.antiprism_map)
            fprintf(ofile,"%-5d = ",i);
         fprintf(ofile,"%3d %3d %3d\n",cols[i][0],cols[i][1],cols[i][2]);
      }
   }
      
   return 0;
}
