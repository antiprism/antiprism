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
   Name: col_util.cc
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


class col_util_opts: public prog_opts {
   public:
      string ifile;
      string ofile;
      string gfile;
      
      int display_type;
      int color_system_mode;
      int map_type;
      int show_container;
      bool unique_colors;
      int upright_view;
      int chroma_level;
      int hsl_height;
      bool plot_centroid;
      double sat_threshold;
      double value_power;
      double value_advance;
      int alpha_mode;
      bool cmy_mode;
      bool ryb_mode;
      bool seven_mode;
      double brightness_adj;
      int map_maximum;

      vector<double> sat_powers;
      vector<double> value_powers;
      
      color_map_multi map;

      col_util_opts(): prog_opts("col_util"),
                     display_type(1),
                     color_system_mode(3),
                     map_type(1),
                     show_container(2),
                     unique_colors(false),
                     upright_view(0),
                     chroma_level(0),
                     hsl_height(0),
                     plot_centroid(false),
                     sat_threshold(1.0),
                     value_advance(0.0),
                     alpha_mode(3),
                     cmy_mode(false),
                     ryb_mode(false),
                     seven_mode(false),
                     brightness_adj(0),
                     map_maximum(0)
                     {}

      void process_command_line(int argc, char **argv);
      void usage();
};

void col_util_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options]\n"
"\n"
"various plots of maps, colors and blendings\n"
"\n"
"Options\n"
"%s"
"  -o <file> write output to file (default: write to standard output)\n"
"  -f <type> for output type 4, map is output instead of OFF file\n"
"               map type: rgb=1  antiprism=2  decimal=3 (default: 1)\n"
"\nScene Options\n"
"  -M <mode> color system mode. HSV=1  HSL=2  RGB=3 (default: 3)\n"
"               plot as cone, bicone or cube respectively\n"
"  -d <int>  output type. plot=1  wheel=2  grid=3  map=4 (default: 1)\n"
"               if color wheel and RGB, only one color blend is possible\n"
"\nOutput type 1 options\n" 
"  -r <int>  HSV/HSL chroma  0 - none (cylinder)  1 - conic  2 - hexagonal\n"
"  -S        HSV/HSL distribute colors heptagonally, 7 ways, as in the rainbow\n"
"  -k <int>  container control (default: 2)\n"
"               0 - suppress  1 - no facets on chroma-2  2 - full\n"
"  -z <int>  1 - show model upright  2 - upright but don\'t rotate RGB cube\n"
"  -q <int>  height of the HSL cones (default: no change)\n"
"               1 - double  2 - make 90 degree dihedral\n"
"  -p        plot color centroid(s)\n"
"\nColor Blending Options\n"
"  -s <sat>  HSV/HSL saturation curve. Greather than 0 (default: 1)\n"
"               1.0 - no curve. lower than 1.0 makes blends more pastel\n"
"               4 numbers can be entered seperated by commas\n"
"  -t <val>  HSV/HSL threshold to use average saturation (default: 1)\n"
"               between 0.0 (all averaging) and 1.0 (no averaging)\n"
"  -v <val>  HSV/HSL value curve (default: 0)\n"
"               simulates subtractive coloring for blending 3 or more colors\n"
"               RGB: Red+Green+Blue = White   Cyan+Magenta+Yellow = Black\n"
"               RYB: Red+Yellow+Blue = Black  Green+Magenta+Orange = White\n"
"               1.0 - no curve. lower than 1.0 number makes blends lighter\n"
"               0.0 - use average value instead\n"
"               4 numbers can be entered seperated by commas\n"
"  -u <val>  HSV/HSL value advance. Rotates meaning of blend to white and black\n"
"               valid values 0.0 to 120.0 degrees (default: 0)\n"
"  -a <int>  alpha to use for blend. average=1  minimum=2  maximum=3 (default: 3)\n"
"  -y        RYB mode. Blend colors as in Red-Yellow-Blue color wheel\n"
"  -c        CMY mode. Complementary colors.  RGB->(RYB/GMO)->CMY->blend\n"
"  -b <val>  All color sytem modes. brightness adjustment for final blended color\n"
"               valid values -1.0 to +1.0 (default: 0)\n"
"               negative for darker, positive for lighter, 0 for no change\n"
"\nColoring Options (run 'off_util -H color' for help on color formats)\n"
"  -m <map>  get colors from a color map, or multiple maps seperated by commas\n"
"  -O <file> get colors from an OFF file\n"
"               note: -m and -O may be used together\n"
"  -U        allow only unique colors (sorts by color)\n"
"  -Z <int>  maximum entries to read from open ended maps (default: 256)\n"
"\n"
"\n",prog_name(), help_ver_text);
}

void col_util_opts::process_command_line(int argc, char **argv)
{
   opterr = 0;
   char c;
   char errmsg[MSG_SZ];
   
   vector<double> double_parms;
   string id;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hd:m:k:q:z:f:r:b:s:t:u:v:pa:cySl:M:O:UZ:o:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         case 'd':
            id = get_arg_id(optarg, "plot=1|wheel=2|grid=3|map=4", argmatch_add_id_maps, errmsg);
            if(id=="")
               error(errmsg);
            display_type = atoi(id.c_str());
            break;
            
         case 'M':
            id = get_arg_id(optarg, "hsv=1|hsl=2|rgb=3", argmatch_add_id_maps, errmsg);
            if(id=="")
               error(errmsg);
            color_system_mode = atoi(id.c_str());
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
            
         case 'z':
            if(!read_int(optarg, &upright_view, errmsg))
               error(errmsg, c);
            if(upright_view < 1 || upright_view > 2)
               error("view type must be between 1 and 2", c);
            break;
            
         case 'f':
            id = get_arg_id(optarg, "rgb=1|antiprism=2|decimal=3", argmatch_add_id_maps, errmsg);
            if(id=="")
               error(errmsg);
            map_type = atoi(id.c_str());
            break;
            
         case 'r':
            if(!read_int(optarg, &chroma_level, errmsg))
               error(errmsg, c);
            if(chroma_level < 0 || chroma_level > 3)
               error("display type must be between 0 and 3", c);
            break;

        case 'b':
            if(!read_double(optarg, &brightness_adj, errmsg))
               error(errmsg, c);
            if(brightness_adj < -1.0 || brightness_adj > 1.0)
               error("brightness adjustment must be between -1.0 and 1.0", c);
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
            
         case 'v':
            if(!read_double_list(optarg, double_parms, errmsg, 4))
               error(errmsg, c);
            for (unsigned int i=0;i<double_parms.size();i++) {
               if(double_parms[i] < 0.0)
                  error("value curve must be greater than or equal to zero", c);
               value_powers.push_back(double_parms[i]);
            }
            break;
            
         case 'p':
            plot_centroid = true;
            break;

         case 'a':
            id = get_arg_id(optarg, "average=1|minimum=2|maximum=3", argmatch_add_id_maps, errmsg);
            if(id=="")
               error(errmsg);
            alpha_mode = atoi(id.c_str());
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
            
         case 'm':
            if(!map.init(optarg, errmsg))
               error(errmsg, c);
            break;
            
         case 'O':
            gfile = optarg;
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

   if (display_type == 3 || display_type == 4) {
      if (plot_centroid)
         warning("plotting centroid colors not valid in this output mode","b");
      if (sat_powers.size())
         warning("saturation entries are not valid in this output mode","s");
      if (value_powers.size())
         warning("value entries are not valid in this output mode","v");
   }

   // fill in missing sat_powers with -1.0, meaning use centroid saturation
   for (int i=sat_powers.size();i<4;i++)
      sat_powers.push_back(-1.0);
      
   // fill in missing value_powers with -1.0, meaning use average values
   for (int i=value_powers.size();i<4;i++)
      value_powers.push_back(-1.0);
      
   if (hsl_height && color_system_mode != 2)
      warning("HSL height adjustment has no effect in RGB or HSV mode","q");

   if (show_container == 1 && (color_system_mode == 3 || chroma_level < 2))
      warning("facets can only be shown in HSV or HSL chroma 2","k");
      
   if (upright_view == 2 && color_system_mode != 3)
      warning("view 2 has no effect when not in RGB mode","z");
      
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

// furnished by Adrian Rossiter
void color_wheel(col_geom_v &geom, const vector<col_val> &cols, int color_system_mode,
                 vector<double> &sat_powers, double sat_threshold, vector<double> &value_powers, double value_advance, int alpha_mode, bool ryb_mode, double brightness_adj)
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
      if (color_system_mode == 3) {
         // RK - if RGB mode, only show that blend in all four levels
         col_val col = blend_RGB_centroid(cols, alpha_mode, ryb_mode);
         if (brightness_adj)
            col.set_brightness(brightness_adj);
         geom.add_col_face(face, col);
      }
      else {
         // RK - only change the powers if there is a new one not the default
         // keeps the center of the bullseye consistent with the last valid power
         if (sat_powers[lvl] > -1.0)
            sat_power = sat_powers[lvl];
         if (value_powers[lvl] > -1.0)
            value_power = value_powers[lvl];

         col_val col = blend_HSX_centroid(cols, color_system_mode, sat_power, sat_threshold, value_power, value_advance, alpha_mode, ryb_mode);
         if (brightness_adj)
            col.set_brightness(brightness_adj);
         geom.add_col_face(face, col);
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

void plot_hsx_point(col_geom_v &geom, col_val &col, int color_system_mode, int chroma_level, bool ryb_mode, bool seven_mode, double brightness_adj)
{
   if (!col.is_val())
      return;

   if (brightness_adj)
      col.set_brightness(brightness_adj);

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
   for(unsigned int i=0; i<cols.size(); i++) {
      col_val col = cols[i]; // can't be const
      plot_hsx_point(geom, col, color_system_mode, chroma_level, ryb_mode, seven_mode, 0);
   }
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

void plot_rgb_point(col_geom_v &geom, col_val &col, bool ryb_mode, double brightness_adj)
{
   if (!col.is_val())
      return;

   if (brightness_adj)
      col.set_brightness(brightness_adj);

   col_val rcol = col;
   if (ryb_mode) {
      // only need hue so algorithm doesn't matter
      vec4d hsxa = rcol.get_hsva();
      hsxa[0] = hsx_to_ryb(rad2deg(2*M_PI*hsxa[0]))/360.0;
      rcol.set_hsva(hsxa);
   }

   geom.add_col_vert(vec3d(rcol[0]/255.0,rcol[1]/255.0,rcol[2]/255.0),col);
}

void plot_rgb_cube(col_geom_v &geom, const vector<col_val> &cols, bool ryb_mode)
{
   for(unsigned int i=0; i<cols.size(); i++) {
      col_val col = cols[i]; // can't be const
      plot_rgb_point(geom, col, ryb_mode, 0);
   }
}

col_geom_v make_square()
{
   col_geom_v geom;
      
   geom.add_vert(vec3d(0, 0, 0)); // 0
   geom.add_vert(vec3d(0, 1, 0)); // 1
   geom.add_vert(vec3d(1, 1, 0)); // 2
   geom.add_vert(vec3d(1, 0, 0)); // 3
   
   geom.add_edge(make_edge(0,1));
   geom.add_edge(make_edge(1,2));
   geom.add_edge(make_edge(2,3));
   geom.add_edge(make_edge(3,0));

   vector<int> face;
   face.push_back(0);
   face.push_back(1);
   face.push_back(2);
   face.push_back(3);
   geom.add_face(face);
   
   coloring clrng(&geom);
   clrng.v_one_col(col_val::invisible);
   clrng.e_one_col(col_val::invisible);
   
   return geom;
}

void color_grid(col_geom_v &geom, const vector<col_val> &cols)
{
   col_geom_v sgeom;
   sgeom.append(make_square());

   int cols_sz = cols.size();
   int dim = (int)ceil(sqrt(cols_sz));

   int k = 0;
   for(int i=0; i<dim; i++) {
      for(int j=0; j<dim; j++) {
         col_geom_v tgeom = sgeom;
         if (k<cols_sz && cols[k].is_idx()) {
            tgeom.add_col_vert(vec3d(0.5,0.45,0.0),col_val(0.0,0.0,0.0));
            tgeom.add_col_vert(vec3d(0.5,0.55,0.0),col_val(1.0,1.0,1.0));
            tgeom.add_col_edge(make_edge(4,5),col_val(0.5,0.5,0.5));
         }
         tgeom.transform(mat3d::transl(vec3d(i,j,0)));
         col_val c = (k>=cols_sz ? col_val(col_val::invisible) : (cols[k].is_idx() ? cols[k].get_idx() : cols[k]));
         k++;
         tgeom.set_f_col(0,c);
         geom.append(tgeom);
         tgeom.clear_all();
      }
   }
   sort_merge_elems(geom, "ve", epsilon);
   geom.transform(mat3d::rot(vec3d(0.0,0.0,deg2rad(-90.0))));
}

bool cmp_col(const col_val &a, const col_val &b)
{
   bool ret = false;

   vec4d hsva_a = a.get_hsva();
   vec4d hsva_b = b.get_hsva();
   for(unsigned int i=0; i<4; i++) {
      if (hsva_a[i] != hsva_b[i]) {
         ret = (hsva_a[i] < hsva_b[i]);
         break;
      }
   }

   return ret;
}

class col_cmp
{
public:
   col_cmp() {}
   bool operator() (const col_val &a, const col_val &b) { return cmp_col(a, b); }
};

void collect_col(vector<col_val> &cols, const col_val &col, bool no_indexes)
{
   if (!col.is_set())
      return;
   else
   if (no_indexes && col.is_idx())
      return;
   cols.push_back(col);
}

void collect_cols_from_geom(const col_geom_v &geom, vector<col_val> &cols, bool no_indexes)
{
   for(unsigned int i=0; i<geom.verts().size(); i++)
      collect_col(cols, geom.get_v_col(i), no_indexes);
   for(unsigned int i=0; i<geom.edges().size(); i++)
      collect_col(cols, geom.get_e_col(i), no_indexes);
   for(unsigned int i=0; i<geom.faces().size(); i++)
      collect_col(cols, geom.get_f_col(i), no_indexes);
}

void collect_cols(vector<col_val> &cols, col_util_opts &opts)
{
   int map_sz = opts.map.effective_size();
   bool open_ended_map = (map_sz >= INT_MAX);
   // map size priority:
   // -Z given
   // default for unlimited map (256)
   // use actual map size
   int max_map_sz = (opts.map_maximum) ? opts.map_maximum : ((open_ended_map) ? 256 : map_sz);

   if (open_ended_map)
      opts.warning(msg_str("map list: only %d out of %d map entries read in",max_map_sz,map_sz), 'Z');
   for(int j=0; j<max_map_sz; j++)
      collect_col(cols, opts.map.get_col(j), (opts.display_type==1 || opts.display_type==2));
   
   // file
   if(opts.gfile.length()) {
      col_geom_v geom;
      char errmsg[MSG_SZ];
      if(!geom.read(opts.gfile, errmsg))
         opts.error(errmsg);
      if(*errmsg)
         opts.warning(errmsg);
      
      collect_cols_from_geom(geom, cols, (opts.display_type==1 || opts.display_type==2));
   }

   if (opts.unique_colors) {
      sort( cols.begin(), cols.end(), col_cmp() );    
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
         cols[i] = rgb_complement(cols[i], opts.ryb_mode);
         
   for(unsigned int i=0; i<cols.size(); i++) {
      if (cols[i][3] < 32) {
         opts.warning("colors with a very small alpha values may not be seen");
         break;
      }
   }

   // grid may have indexes
   if (opts.display_type == 3) {
      for(unsigned int i=0; i<cols.size(); i++) {
         if (cols[i].is_idx()) {
            opts.warning("color indexes detected. unmapped cells will result");
            break;
         }
      }
   }
}

int main(int argc, char *argv[])
{
   col_util_opts opts;
   opts.process_command_line(argc, argv);
   
   col_geom_v geom;
   char errmsg[MSG_SZ];
      
   vector<col_val> cols;
   collect_cols(cols, opts);
   
   if (!cols.size())
      opts.error("found no color values to plot");
      
   if (cols.size() < 3) {
      for(unsigned int i=0; i<4; i++) {
         if (opts.value_powers[i] > -1.0) {
            opts.warning("subtractive colors only works on three colors or more","v");
            break;
         }
      }
   }

   // output map file
   if(opts.display_type == 4) {
      FILE *ofile = stdout;  // write to stdout by default
      if(opts.ofile != "") {
         ofile = fopen(opts.ofile.c_str(), "w");
         if(ofile == 0)
            opts.error("could not open output file \'"+opts.ofile+"\'");
      }
      for(unsigned int i=0; i<cols.size(); i++) {
         if (opts.map_type == 3) {
            if (cols[i].is_val()) {
               vec4d c = cols[i].get_vec4d();
               fprintf(ofile, "%g%s %g%s %g%s %g%s\n",
                  c[0], (c[0]==1.0 || c[0]==0.0) ? ".0" : "",
                  c[1], (c[1]==1.0 || c[1]==0.0) ? ".0" : "",
                  c[2], (c[2]==1.0 || c[2]==0.0) ? ".0" : "",
                  c[3], (c[3]==1.0 || c[3]==0.0) ? ".0" : "");
            }
         }
         else {
            if (opts.map_type == 2)
               fprintf(ofile,"%-5d = ",i);
            char buffer[MSG_SZ];
            buffer[0] = '\0';
            if (cols[i][3] != 255)
               sprintf(buffer,"%3d",cols[i][3]);
            if (cols[i].is_val())
               fprintf(ofile,"%3d %3d %3d %s\n",cols[i][0],cols[i][1],cols[i][2],buffer);
            else
               fprintf(ofile,"%3d\n",cols[i].get_idx());
         }
      }
      if(opts.ofile!="")
         fclose(ofile);
   }
   // else it is a plot/model
   else {
      // a plot
      if (opts.display_type == 1) {
         if (opts.color_system_mode == 1 || opts.color_system_mode == 2) {
            plot_hsx(geom, cols, opts.color_system_mode, opts.chroma_level, opts.ryb_mode, opts.seven_mode);
            
            if (opts.plot_centroid) {
               for(unsigned int i=0; i<4; i++) {
                  if (!i || opts.sat_powers[i] > -1.0 || opts.value_powers[i] > -1.0) {
                     col_val col = blend_HSX_centroid(cols, opts.color_system_mode, opts.sat_powers[i], opts.sat_threshold, opts.value_powers[i], opts.value_advance, opts.alpha_mode, opts.ryb_mode);
                     plot_hsx_point(geom, col, opts.color_system_mode, opts.chroma_level, opts.ryb_mode, opts.seven_mode, opts.brightness_adj);
                  }
               }
            }
         }
         else
         if (opts.color_system_mode == 3) {
            plot_rgb_cube(geom, cols, opts.ryb_mode);
            
            if (opts.plot_centroid) {
               col_val col = blend_RGB_centroid(cols, opts.alpha_mode, opts.ryb_mode);
               plot_rgb_point(geom, col, opts.ryb_mode, opts.brightness_adj);
            }
         }

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
         if (opts.color_system_mode == 2 && opts.hsl_height) {
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
      else
      if (opts.display_type == 2)
         color_wheel(geom, cols, opts.color_system_mode, opts.sat_powers, opts.sat_threshold, opts.value_powers, opts.value_advance, opts.alpha_mode, opts.ryb_mode, opts.brightness_adj);
      else
      if (opts.display_type == 3)
         color_grid(geom, cols);

      if(!geom.write(opts.ofile, errmsg))
         opts.error(errmsg);
   }
      
   return 0;
}

