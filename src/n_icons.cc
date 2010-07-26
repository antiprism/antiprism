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
   Name: n_icons.cc
   Description: Creates Sphericon like Polyhedra
   Project: Antiprism - http://www.antiprism.com
*/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include <ctype.h>
#include <unistd.h>

#include <string>
#include <vector>
#include <algorithm>

#include "../base/antiprism.h"
#include "n_icons.h"

using std::string;
using std::vector;
using std::swap;
using std::map;
using std::pair;

// default color needs to be a grey, but one not often chosen
#define DEFAULT_COLOR 193,192,191,255


// Need these here since used in opts
bool half_model(vector<int> longitudes)
{
   return(2*longitudes.back() == longitudes.front());
}

bool full_model(vector<int> longitudes)
{
   return(longitudes.front() == longitudes.back());
}

color_map_map * alloc_default_map()
{
   color_map_map *col_map = new color_map_map;

   col_map->set_col(0, col_val(1.0,0.0,0.0));      // red
   col_map->set_col(1, col_val(1.0,0.49804,0.0));  // darkoranage1
   col_map->set_col(2, col_val(1.0,1.0,0.0));      // yellow
   col_map->set_col(3, col_val(0.0,0.39216,0.0));  // darkgreen
   col_map->set_col(4, col_val(0.0,1.0,1.0));      // cyan
   col_map->set_col(5, col_val(0.0,0.0,1.0));      // blue
   col_map->set_col(6, col_val(1.0,0.0,1.0));      // magenta
   col_map->set_col(7, col_val(1.0,1.0,1.0));      // white
   col_map->set_col(8, col_val(0.5,0.5,0.5));      // grey
   col_map->set_col(9, col_val(0.0,0.0,0.0));      // black

   col_map->set_wrap();

   return col_map;
}

color_map_map * alloc_no_color_map()
{
   color_map_map *col_map = new color_map_map;

   col_map->set_col(0, col_val());  // unset edge color

   col_map->set_wrap();

   return col_map;
}

class ncon_opts: public prog_opts {
   public:
      string ofile;

      int ncon_order;
      int d;
      bool nm_shell;
      bool hide_indent;
      double inner_radius;
      double outer_radius;
      bool point_cut;
      bool hybrid;
      bool add_poles;
      int twist;
      bool info;
      char face_coloring_method;
      int face_opacity;
      string face_pattern;
      bool face_sequential_colors;
      char edge_coloring_method;
      int edge_opacity;
      string edge_pattern;
      bool edge_set_no_color;
      col_val unused_edge_color;
      bool edge_sequential_colors;
      bool symmetric_coloring;
      string closure;
      vector<int> longitudes;
      string hide_elems;
      string ncon_surf;
      vector<int> ncon_range;
      bool long_form;
      bool filter_case2;

      color_map_multi face_map;
      color_map_multi edge_map;
      
      // former global variable
      bool split;
      int half_model_marker;

      ncon_opts(): prog_opts("n_icons"),
                   ncon_order(4),
                   d(1),
                   nm_shell(true),
                   hide_indent(true),
                   inner_radius(FLT_MAX),
                   outer_radius(FLT_MAX),
                   point_cut(true),
                   hybrid(false),
                   add_poles(false),
                   twist(1),
                   info(false),
                   face_coloring_method('s'),
                   face_opacity(-1),
                   face_pattern("1"),
                   face_sequential_colors(true),
                   edge_coloring_method('\0'),
                   edge_opacity(-1),
                   edge_pattern("1"),
                   edge_set_no_color(false),
                   unused_edge_color(col_val::invisible),
                   edge_sequential_colors(false),
                   symmetric_coloring(false),
                   long_form(false),
                   filter_case2(false),
                   split(false),
                   half_model_marker(0)
             {}

      void process_command_line(int argc, char **argv);
      void usage();
};

void ncon_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options]\n"
"\n"
"Creates Sphericon like Polyhedra.\n"
"\n"
"Options\n"
"%s"
"  -n <n/d>  n-icon of order n. n must be 3 or greater (default: 4)\n"
"               use d to make star n-icon. d less than n\n"
"  -t <twst> number of twists. Can be negative, positive or 0 (default: 1)\n"
"  -s        side-cut of even order n-icon (default is point-cut)\n"
"  -H        hybrid of even order n-icon\n"
"               -a -c and m2 have no effect with hybrids\n"
"  -M <m,m2> longitudes of model of m sides with optional m2 of m sides showing\n"
"               m must be even and 3 or greater (default: 36,36)\n"
"  -x <elms> t and b to exclude top and/or bottom polygons if they exist\n"
"               v, e and f to remove OFF faces with one vertex (vertices),\n"
"               two-vertices (edges) and three or more vertices (faces)\n"
"  -a        place a north and south pole in top and bottom if they exist\n"
"                only valid if m2<m. Not valid with -c h\n"
"  -c <clse> close open model if m2<m. Valid values h or v\n"
"               h = horizontal closure  v = vertical closure\n"
"  -z        for n/d, make non-shell model. n/d must be co-prime\n"
"  -r        for n/d shell model, override inner radius. greater than 0\n"
"  -R        for n/d shell model, override outer radius. greater than 0\n"
"  -I        info on current n-icon\n"     
"  -o <file> write output to file (default: write to standard output)\n"
"\nColoring Options\n"
"  -f <mthd> mthd is face coloring method. The coloring is done before twist\n"
"               key word: none - sets no color (default: s)\n"
"               lower case outputs map indexes. upper case outputs color values\n"
"               s - color circuits with colors using list sequentially\n"
"               t - color circuits with colors using circuit numbers\n"
"               l - color latitudinally\n"
"               m - color longitudinally\n"
"               c - checkerboard with first two colors in face color list\n"
"               n - use each color in succession\n"
"               x - first two colors based on sign of x\n"
"               y - first two colors based on sign of y\n"
"               z - first two colors based on sign of z\n"
"                      note: z is also the twist plane\n"
"               o - use first eight colors per xyz octants\n"
"  -S        color circuits symmetrically when using coloring method s or t\n"
"  -T <tran> face transparency. valid range from 0 to 255\n"
"               0 - invisible  255 - opaque (default: 255)\n"
"  -O <strg> face transparency pattern string. valid values\n"
"               0 - T value suppressed  1 - T value applied  (default: '1')\n"
"  -e <mthd> mthd is edge coloring method. The coloring is done before twist\n"
"               key word: none - sets no color (default: Q)\n"
"               key word: Q - defer coloring all edges to option Q\n"
"                  or use the same letter options specified in -f\n"
"  -U <tran> edge transparency. valid range from 0 to 255\n"
"               0 - invisible  255 - opaque (default: 255)\n"
"  -P <strg> edge transparency pattern string. valid values\n"
"               0 - U value suppressed  1 - U value applied  (default: '1')\n"
"  -Q <col>  color given to uncolored edges and vertices of final model\n"
"               key word: none - sets no color (default: invisible)\n"
"  -Y        for n/d shells, when showing edges, show indented edges\n"
"  -m <maps> color maps to be tried in turn. (default: map_red:darkorange1:\n"
"               yellow:darkgreen:cyan:blue:magenta:white:grey:black%%) optionally\n"
"               followed by elements to map from v, e or f (default: vef)\n"
"\nSurface Count Reporting (options above igonored)\n"
"  -J <type> list n-icons with more than one surface. Valid values for type\n"
"               n = point cut even order n_icons\n"
"               s = side cut even order n-icons (surfaces > 2)\n"
"               o = odd order n_icons\n"
"               h = hybrids (all)\n"
"               i = hybrids (where N/2 is even)\n"
"               j = hybrids (where N/2 is odd)\n"
"               k = hybrids (where N/4 is even)\n"
"               l = hybrids (where N/4 is odd)\n"
"  -K <k,k2> range of n-icons to list for multiple surfaces\n"
"  -L        long form report\n"
"  -Z        filter out case 2 types\n"
"\n"
"\n", prog_name(), help_ver_text);
}

void ncon_opts::process_command_line(int argc, char **argv)
{
   opterr = 0;
   char c;
   char errmsg[MSG_SZ];

   coloring clrngs[3];
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hn:t:sHM:x:ac:zr:R:IJ:K:LZm:f:ST:O:e:U:P:Q:Yo:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         case 'n': {
            char *p;
            p = strchr(optarg, '/');
            if(p!=0) {
               *p++='\0';
               if(!read_int(p, &d, errmsg))
                  error(errmsg, "n/d (d part)");
            }

            if(!read_int(optarg, &ncon_order, errmsg))
               error(errmsg, "n/d (n part)");
            if(ncon_order<3)
               error("n must be an integer 3 or greater", "n/d (n part)");
            if(d < 1)
               error("d must be 1 or greater", "n/d (d part)");

            if(d >= ncon_order)
               error("d must be less than n", "n/d (d part)");
            break;
         }

         case 't':
            if(!read_int(optarg, &twist, errmsg))
               error(errmsg, c);
            break;
            
         case 's':
            point_cut = false;
            break;
            
         case 'H':
            hybrid = true;
            break;
            
         case 'M':
            if(!read_int_list(optarg, longitudes, errmsg, true, 2))
               error(errmsg, c);
            if(longitudes.front()<3)
               error("m must be 3 or greater", c);
            if(longitudes.size() == 1)
               longitudes.push_back(longitudes.front());
            else
            if(longitudes.size() == 2) {
               if(longitudes.front() < longitudes.back())
                  error("sides shown (m2) must be less than or equal to longitudes (m)", c);
               if(longitudes.back() <= 0)
                  error("sides shown (m2) must be at least 1", c);
            }
            break;

         case 'x':
            if(strspn(optarg, "tbvef") != strlen(optarg))
               error(msg_str("elements to hide are '%s' must be from "
                        "t, b, v, e and f", optarg), c);
            hide_elems=optarg;
            break;

         case 'a':
            add_poles = true;
            break;

         case 'c':
            if(strspn(optarg, "hv") != strlen(optarg) || strlen(optarg)==2)
               error(msg_str("closure is '%s', must be h or v (not both)",
                        optarg), c);
            closure=optarg;
            break;

         case 'z':
            nm_shell = false;
            break;

         case 'r':
            if(!read_double(optarg, &inner_radius, errmsg))
               error(errmsg, c);
           if (inner_radius <= 0.0)
               error("inner radius must be greater than 0", c);
            break;

         case 'R':
            if(!read_double(optarg, &outer_radius, errmsg))
               error(errmsg, c);
            if (outer_radius <= 0.0)
               error("outer radius must be greater than 0", c);
            break;

         case 'I':
            info = true;
            break;

         case 'J':
            if(strspn(optarg, "nsohijkl") != strlen(optarg) || strlen(optarg)>1)
               error(msg_str("n-icon type is '%s', must be only one of n, s, "
                        "o, h, i, j, k, or l\n", optarg), c);
            ncon_surf=optarg;
            break;

         case 'K':
            if(!read_int_list(optarg, ncon_range, errmsg, true, 2))
               error(errmsg, c);
            if(ncon_range.front()<3)
               error("k must be 3 or greater", c);
            if(ncon_range.size() == 1)
               ncon_range.push_back(ncon_range.front());
            else
            if((ncon_range.size() == 2) && (ncon_range.front() > ncon_range.back()))
               error("k2 shown must be greater than or equal to k", c);
            break;

         case 'L':
            long_form = true;
            break;

         case 'Z':
            filter_case2 = true;
            break;

         case 'm':
            if(!read_colorings(clrngs, optarg, errmsg))
               error(errmsg, c);
            break;
            
         case 'f':
            if(!strcasecmp(optarg,"none"))
               face_coloring_method = '\0';
            else
            if(strspn(optarg, "stlmcnxyzo") != strlen(optarg) || strlen(optarg)>1)
               error(msg_str("invalid face coloring method '%c'", *optarg), c);
            else {
               face_coloring_method = *optarg;
               if ((face_coloring_method == 's') || (face_coloring_method == 't')) {
                  face_sequential_colors = (face_coloring_method == 's') ? true : false;
                  face_coloring_method = 's';
               }
            }
            break;

         case 'S':
            symmetric_coloring = true;
            break;

         case 'T':
            if(!read_int(optarg, &face_opacity, errmsg))
               error(errmsg, c);
            if(face_opacity < 0 || face_opacity > 255) {
               error("face transparency must be between 0 and 255", c);
            }
            break;

         case 'O':
            if(strspn(optarg, "01") != strlen(optarg))
               error(msg_str("transparency string is '%s', must consist of "
                        "0 and 1's", optarg), c);
            face_pattern=optarg;
            break;

         case 'e':
            if(!strcasecmp(optarg,"none")) {
               edge_set_no_color = true;
               edge_coloring_method = 'n';
            }
            else
            if(!strcmp(optarg,"Q"))
              edge_coloring_method = '\0';
            else
            if(strspn(optarg, "stlmcnxyzo") != strlen(optarg) || strlen(optarg)>1)
               error(msg_str("invalid edge coloring method '%s'", optarg), c);
            else {
               edge_coloring_method = *optarg;
               if ((edge_coloring_method == 's') || (edge_coloring_method == 't')) {
                  edge_sequential_colors = (edge_coloring_method == 's') ? true : false;
                  edge_coloring_method = 's';
               }
            }
            break;

         case 'U':
            if(!read_int(optarg, &edge_opacity, errmsg))
               error(errmsg, c);
            if(edge_opacity < 0 || edge_opacity > 255) {
               error("edge transparency must be between 0 and 255", c);
            }
            break;
            
         case 'P':
            if(strspn(optarg, "01") != strlen(optarg))
               error(msg_str("transparency string %s must consist of 0 and 1's",
                        optarg), c);
            edge_pattern=optarg;
            break;
            
         case 'Q':
            if(!unused_edge_color.read(optarg, errmsg))
               error(errmsg, c);
            break;

        case 'Y':
            hide_indent = false;
            break;

         case 'o':
            ofile = optarg;
            break;

         default:
            error("unknown command line error");
      }
   }

   if(argc-optind > 0)
      error("too many arguments");

   // surfaces subsystem
   if (ncon_surf.length() > 0 ) {
      if (ncon_range.size() == 0)
         error("for surfaces reporting -K must be specified","J");
      else
      if (!is_even(ncon_range.front()) && ncon_range.front() == ncon_range.back()) {
         if (strchr(ncon_surf.c_str(), 'h') || strchr(ncon_surf.c_str(), 'i') || strchr(ncon_surf.c_str(), 'j'))
            error("for hybrid surfaces k must be even","K");
         else
         if (strchr(ncon_surf.c_str(), 'n'))
            error("for even order n-icon surfaces k must be even","K");
         else
         if (strchr(ncon_surf.c_str(), 's'))
            error("for side cut n-icon surfaces k must be even","K");
      }
      else {
      if (is_even(ncon_range.front()) && ncon_range.front() == ncon_range.back())
         if (strchr(ncon_surf.c_str(), 'o'))
            error("-K: for odd order n-icons surfaces k must be odd");
      }
   }
   // n_icons option processing
   else {
      // Default longitudes to use is 36
      if (longitudes.size() == 0 ) {
         longitudes.push_back(36);
         longitudes.push_back(36);
      }

     if (d == 1 || (ncon_order-d) == 1)
         nm_shell = false;

      if (!nm_shell)
         if (ncon_order>0 && d>0 && gcd(ncon_order,d) != 1)
            error("when not making shells, n and d must be co-prime","z");

      if (!hide_indent && !nm_shell)
         warning("show indented edges only valid in n/d shells","y");

      if (!hide_indent && !edge_coloring_method)
         warning("indented edges will not be shown unless an edge coloring is used (-e)","y");
      
      if (ncon_range.size() > 0)
         error("not valid without -J","K");

      if (long_form)
         error("not valid without -J","L");
 
      if (!point_cut) {
         if (!is_even(ncon_order))
            error("side cut is not valid with n-icons of odd n","s");
         if (hybrid)
            error("side cut is not valid with hybrid","s");
      }

      if (hybrid) {
         if (!is_even(ncon_order))
            error("hybrids must be even order","n");

         if (twist == 0)
            error("hybrids have no twist 0","t");

         if (longitudes.front() != longitudes.back())
            warning("for hybrids m2 has no effect. Full models only","M");
         longitudes.back() = longitudes.front()/2;

         if (closure.length() > 0) {
            warning("closure options not valid with -H","c");
            closure.clear();
         }

         if (add_poles) {
            warning("poles not valid with -H","a");
            add_poles = false;
         }
         
         if (hybrid && symmetric_coloring)
            warning("symmetric coloring is the same an non-symmetric coloring for hybrids","S");
      }

      // Let us allow globes (Twist = 0) to have uneven number of longitudes
      if((twist != 0) && (!is_even(longitudes.front())))
         error("if twist -t is not 0 then -M m must be even","M");

      // poles are only for odd order or even order side cuts
      if (is_even(ncon_order) && point_cut) {
         if ( add_poles )
            warning("poles not valid with even point cut n-icons","a");
         add_poles = false;
      }

      // Don't have poles on half or whole models - Cover longitudinal doesn't need poles
      if ((half_model(longitudes) || full_model(longitudes)) || strchr(closure.c_str(), 'h')) {
         if ( add_poles )
            warning("poles not valid with half or full model or when -c h","a");
         add_poles = false;
      }

      // Twist needs split. Can't have poles. (except for hybrids)
      if (twist != 0 && !hybrid) {
         split = true;
         if ( add_poles )
            warning("poles not valid when twist is not 0","a");
         add_poles = false;
      }

      // Full models have no need for closure
      if (full_model(longitudes)) {
         if ( closure.length() > 0 )
            warning("closure not needed for full model","c");
         closure.clear();
      }
   }
   
   if (((face_coloring_method != 's') && (face_coloring_method != 't') &&
        (edge_coloring_method != 's') && (edge_coloring_method != 't')) && symmetric_coloring)
      error("symmetric coloring is only for coloring methods s or t","S");
      
   if (!is_even(ncon_order) && symmetric_coloring)
      warning("symmetric coloring is the same an non-symmetric coloring for odd order n-icons","S");

   if((clrngs[2].get_cmaps()).size())
      face_map = clrngs[2];
   else
      face_map.add_cmap(alloc_default_map());

   // patch for setting edges with no color
   if (edge_set_no_color)
      edge_map.add_cmap(alloc_no_color_map());
   else
   // process as was done with face map
   if((clrngs[1].get_cmaps()).size())
      edge_map = clrngs[1];
   else
      edge_map.add_cmap(alloc_default_map());

   // vertex color map is the same as the edge color map
   //if((clrngs[0].get_cmaps()).size())
   //   warning("vertex map has no effect","m");
}

int longitudinal_faces(int ncon_order, bool point_cut)
{
   int lf = (int)floor((double)ncon_order/2);
   if ( is_even(ncon_order) && !point_cut )
      lf--;
   return(lf);
}

int num_lats(int ncon_order, bool point_cut)
{
   int lats = 0;
   if (is_even(ncon_order) && point_cut)
      lats = (ncon_order/2);
   else if (is_even(ncon_order) && !point_cut)
      lats = (ncon_order/2)+1;
   else if (!is_even(ncon_order))
      lats = (int)ceil((double)ncon_order/2);
   return lats;
}

void add_coord(col_geom_v &geom, vector<coordList *> &coordinates, vec3d vert)
{
   coordinates.push_back(new coordList(geom.add_vert(vert)));
}

void clear_coord(vector<coordList *> &coordinates)
{
   for (unsigned int i=0;i<coordinates.size();i++)
      delete coordinates[i];
   coordinates.clear();
}

void add_face(col_geom_v &geom, vector<faceList *> &face_list, vector<int> face, int lat, int lon)
{
   face_list.push_back(new faceList(geom.add_face(face),lat,lon));
}

void clear_faces(vector<faceList *> &face_list)
{
   for (unsigned int i=0;i<face_list.size();i++)
      delete face_list[i];
   face_list.clear();
}

void add_edge(col_geom_v &geom, vector<edgeList *> &edge_list, vector<int> edge, int lat, int lon)
{
   if(edge[0]>edge[1])
      swap(edge[0], edge[1]);
   edge_list.push_back(new edgeList(geom.add_edge_raw(edge),lat,lon));
}

void clear_edges(vector<edgeList *> &edge_list)
{
   for (unsigned int i=0;i<edge_list.size();i++)
      delete edge_list[i];
   edge_list.clear();
}

void build_prime_meridian(col_geom_v &geom, vector<int> &prime_meridian, vector<coordList *> &coordinates, int ncon_order, int d, bool nm_shell,
                          double &inner_radius, double &outer_radius, bool &point_cut)
{
   double arc = (nm_shell) ? 360.0/ncon_order : 360.0/ncon_order*d;
   double interior_angle = (180.0-arc)/2.0;
   double outer_radius_calc = sin(deg2rad(interior_angle))/sin(deg2rad(arc));
   if (outer_radius == FLT_MAX)
      outer_radius = outer_radius_calc;

   // formula furnished by Adrian Rossiter
   //r = R * cos(pi*m/n) / cos(pi*(m-1)/n)
   if (nm_shell) {
      if (inner_radius == FLT_MAX) {
         int n = ncon_order/2;
         if (2*d>n)
            d = n-d;
         inner_radius = outer_radius_calc * cos(M_PI*d/n) / cos(M_PI*(d-1)/n);
      }
   }

   bool radii_swapped = false;
   double angle = -90.0;
   if ( is_even(ncon_order) && !point_cut ) {
      if (nm_shell) {
         swap(outer_radius,inner_radius);
         radii_swapped = true;
         // now treat it like a point cut
         point_cut = true;
      }
      else
         angle += ( arc / 2.0 );
   }

   int num_vertices = longitudinal_faces(ncon_order, point_cut) + 1;
   for (int i=0;i<num_vertices;i++) {
      prime_meridian.push_back(i);
      double radius = (!nm_shell || is_even(i)) ? outer_radius : inner_radius;
      add_coord(geom, coordinates, vec3d(cos(deg2rad(angle))*radius, sin(deg2rad(angle))*radius, 0));
      angle += arc;
   }

   // swap these back for future reference
   if ( radii_swapped )
      swap(outer_radius,inner_radius);
}

void form_globe(col_geom_v &geom, vector<int> &prime_meridian, vector<coordList *> &coordinates, vector<faceList *> &face_list, vector<edgeList *> &edge_list,
                char edge_coloring_method, int ncon_order, bool point_cut, vector<int> longitudes, string closure, int &half_model_marker)
{
   const vector<vector<int> > &faces = geom.faces();
   const vector<vec3d> &verts = geom.verts();

   double arc = 360.0/longitudes.front();

   int lon_faces = longitudinal_faces(ncon_order, point_cut);
   // used to coordinate latitudes with bands
   int inc1 = (is_even(ncon_order) && point_cut) ? 0 : 1;
   
   // We close an open model longitudinally only if is desired, not a full model, or only 1 zone is shown
   int longitudinal_cover = 0;
   if (strchr(closure.c_str(), 'h') && !full_model(longitudes) && longitudes.back() != 1 )
      longitudinal_cover++;

   for (int i=1;i<=longitudes.back()+longitudinal_cover;i++) {
      for (int j=0;j<(int)prime_meridian.size();j++) {
         if (((is_even(ncon_order) && point_cut) && (j != 0) && (j < prime_meridian.back())) ||
              (is_even(ncon_order) && !point_cut) ||
              (!is_even(ncon_order) && (j > 0))) {
            if ((i != longitudes.front()) && (i != longitudes.back() + 1)) {
               // Rotate Point Counter-Clockwise about Y Axis (looking down through Y Axis)
               vec3d v = mat3d::rot(0,deg2rad(arc*i),0) *
                  verts[prime_meridian[j]];
               add_coord(geom, coordinates, v);
            } 
            if ((i == ((longitudes.front()/2) + 1)) && (half_model_marker == 0))
               half_model_marker = verts.size() - 2;
         }

         bool use_prime;
         int k = 0;
         int l = 0;

         // store in face list for coloring faces later
         int lat = lon_faces - j + inc1;
         int lon = i-1;

         if ((i == 1) || (i == longitudes.front()) || (i == longitudes.back() + 1)) {
            use_prime = true;
            l = 0;

            // Reset coordinate count for last row of faces
            if ((i == longitudes.front()) || (i == longitudes.back()+1)) {
               k = lon_faces - j;
               if (is_even(ncon_order) && point_cut)
                  k--;
               // no color for longitudes that span more than one meridian
               if (longitudes.front() - longitudes.back() > 1) {
                  lat = -1;
                  lon = -1;
               }
            }
         }
         else {
            use_prime = false;
            l = faces.size() - lon_faces;
         }

         // when j = 0 "prime the system"
         if ( j > 0 ) {
            if (((is_even(ncon_order) && point_cut) || !is_even(ncon_order)) && (j == 1)) {
               //fprintf(stderr,"South Triangle\n");
               vector<int> face;
               face.push_back(verts.size()-1-k);
               face.push_back(prime_meridian.front());

               if ( use_prime )
                  face.push_back(prime_meridian[1]);
               else
                  face.push_back(faces[l][0]);
               add_face(geom, face_list, face, lat, lon);

               if (edge_coloring_method) {
                  vector<int> edge;
                  edge.push_back(face[0]);
                  edge.push_back(face[2]);
                  add_edge(geom, edge_list, edge, lat, lon);
               }
            }
            else 
            if ((is_even(ncon_order) && point_cut) && (j == (int)prime_meridian.size()-1)) {
               //fprintf(stderr,"North Triangle\n");
               vector<int> face;
               face.push_back(verts.size()-1);
               face.push_back(prime_meridian.back());

               if ( use_prime )
                  face.push_back(prime_meridian.back()-1);
               else
                  face.push_back(faces[l][0]);
               add_face(geom, face_list, face, lat, lon);
            }
            else {
               //fprintf(stderr,"Square\n");
               vector<int> face;
               face.push_back(verts.size()-1-k);
               face.push_back(verts.size()-2-k);

               if ( use_prime ) {
                  face.push_back(prime_meridian[prime_meridian[j-1]]);
                  face.push_back(prime_meridian[prime_meridian[j]]);
               }
               else {
                  face.push_back(faces[l][1]);
                  face.push_back(faces[l][0]);
               }
               add_face(geom, face_list, face, lat, lon);

               if ( edge_coloring_method ) {
                  vector<int> edge;
                  edge.push_back(face[0]);
                  edge.push_back(face[3]);
                  add_edge(geom, edge_list, edge, lat, lon);
                  
                  // patch with the extra edge below, one rotate flag gets missed
                  if (half_model_marker != 0)
                     edge_list.back()->rotate = true;

                  // the last square has two edges if there is no bottom triangle
                  if (lat == lon_faces && (is_even(ncon_order) && !point_cut)) {
                     edge.clear();
                     edge.push_back(face[1]);
                     edge.push_back(face[2]);
                     add_edge(geom, edge_list, edge, lat+1, lon);
                  }
               }
            }

            if (half_model_marker != 0) {
               face_list.back()->rotate = true;
               if ( edge_coloring_method )
                  edge_list.back()->rotate = true;
            }
         }
      }
   }
}

void add_caps(col_geom_v &geom, vector<coordList *> &coordinates, vector<faceList *> &face_list, vector<poleList *> &pole, int ncon_order, bool point_cut, bool hybrid,
              vector<int> longitudes, bool split, bool add_poles, string hide_elems)
{
   const vector<vector<int> > &faces = geom.faces();
   const vector<vec3d> &verts = geom.verts();

   int lats = num_lats(ncon_order,point_cut);
         
   pole[0]->idx = -1;
   pole[0]->lat = 0;
   pole[1]->idx = -1;
   pole[1]->lat = lats; 

   // Even order and point cut always have poles (and south pole = 0)
   if (is_even(ncon_order) && point_cut) {
      pole[0]->idx = longitudinal_faces(ncon_order, point_cut);
      pole[1]->idx = 0;
   }
   else
   if (!is_even(ncon_order) || add_poles ) {
      pole[0]->idx = 0;
      if (!half_model(longitudes) || add_poles)
         pole[1]->idx = 0;
   }

   if (add_poles) {
      if (!strchr(hide_elems.c_str(), 't') && ((is_even(ncon_order) && !point_cut) || !is_even(ncon_order))) {
         vec3d v;
         v[0]=0;
         v[1]=verts.back()[1];
         v[2]=0;
         add_coord(geom, coordinates, v);
         pole[0]->idx = verts.size()-1;

         if (!strchr(hide_elems.c_str(), 'b') && ((is_even(ncon_order) && !point_cut))) { 
            vec3d v;
            v[0]=0;
            v[1]=verts.front()[1];
            v[2]=0;
            add_coord(geom, coordinates, v);
            pole[1]->idx = verts.size()-1;
         }
      }
   }

   // Top
   if (!strchr(hide_elems.c_str(), 't') && ((is_even(ncon_order) && !point_cut) || !is_even(ncon_order))) {
      int j = 1;
      bool split_done = false;
      vector<int> face;
      for (int i=longitudinal_faces(ncon_order, point_cut)-1; j<=longitudes.back(); i+=longitudinal_faces(ncon_order, point_cut)) {
         if (j == longitudes.back()) {
            if (!full_model(longitudes))
               face.push_back(faces[i].back());
            if (split && !split_done && add_poles && (longitudes.back() > longitudes.front()/2))
               face.push_back(pole[0]->idx);
            else
            if (!split || (split && split_done) || half_model(longitudes) || (!half_model(longitudes) && (longitudes.back() != (longitudes.front()/2) + 1)))
               face.push_back(faces[i].front());
         }
         else
            face.push_back(faces[i].back());

         if (split && !split_done && (j == (longitudes.front()/2) + 1)) {
            add_face(geom, face_list, face, 0, -1);
            face.clear();
            i -= longitudinal_faces(ncon_order, point_cut);
            j--;
            split_done = true;
         }
         j++;
      }

      if (split_done && !add_poles)
         face.push_back((faces.back()).front());

      if (add_poles)
         face.push_back(pole[0]->idx);

      add_face(geom, face_list, face, 0, -1);
      if (!hybrid)
         face_list.back()->rotate = true;
   }

   // Bottom
   if (!strchr(hide_elems.c_str(), 'b') && ((is_even(ncon_order) && !point_cut))) {
      int j = 1;
      bool split_done = false;
      vector<int> face;

      for (int i=0; j<=longitudes.back(); i+=longitudinal_faces(ncon_order, point_cut)) {
         if (j == longitudes.back()) {
            if (!full_model(longitudes))
               face.push_back(faces[i][2]);
            if (split && !split_done && add_poles && (longitudes.back() > longitudes.front()/2))
               face.push_back(pole[1]->idx);
            else
            if (!split || (split && split_done) || half_model(longitudes) || (!half_model(longitudes) && (longitudes.back() != (longitudes.front()/2) + 1)))
               face.push_back(faces[i][1]);
         }
         else
            face.push_back(faces[i][2]);

         if (split && !split_done && (j == (longitudes.front()/2) + 1)) { 
            add_face(geom, face_list, face, lats-1, -1);
            face.clear();
            i -= longitudinal_faces(ncon_order, point_cut);
            j--;
            split_done = true;
         }

         j++;
      }

      if (split_done && !add_poles)
         face.push_back((faces.back()).front());

      if ( add_poles )
         face.push_back(pole[1]->idx);

      add_face(geom, face_list, face, lats-1, -1);
      if (!hybrid)
         face_list.back()->rotate = true;
   }
}

// Reuse this code to find the twist plane
void close_latitudinal_or_find_twist_plane(col_geom_v &geom, vector<polarOrb *> &polar_orbit, vector<faceList *> &face_list, vector<poleList *> &pole,
                                           int ncon_order, bool point_cut, vector<int> longitudes, bool add_poles, string closure, int half_model_marker)
{
   //const vector<vector<int> > &faces = geom.faces();   
   vector<vector<int> > face(3);

   // Cover one side            
   if (add_poles && (is_even(ncon_order) && !point_cut))
      face[0].push_back(pole[1]->idx);

   int j = 0;
   for (int i=1;i<=longitudinal_faces(ncon_order, point_cut)+1; i++) {
      polar_orbit.push_back(new polarOrb(j));
      face[0].push_back(j);
      j++;
   }

   // If half model it is all one closure, continue face[0] as face[1]
   if (half_model(longitudes))
      face[1] = face[0];
   else {
      if (add_poles && ((is_even(ncon_order) && !point_cut) || !is_even(ncon_order)))
         face[0].push_back(pole[0]->idx);

      if ( strchr(closure.c_str(), 'v') && (face[0].size() > 2)) {
         add_face(geom, face_list, face[0], -2, -2);
         face_list.back()->rotate = true;
      }
   }

   // Cover other side
   j = geom.verts().size()-1;
   int m = half_model_marker;       

   if (add_poles) {
      if (!is_even(ncon_order))
         j--;
      else
      if (is_even(ncon_order) && !point_cut)
         j -= 2;
   }

   int k = longitudinal_faces(ncon_order, point_cut) + 1;
   if (add_poles || (is_even(ncon_order) && point_cut)) {
      if (!half_model(longitudes))
         face[1].push_back(pole[0]->idx);
      if (is_even(ncon_order) && point_cut)
         k--;
   }

   for (int i=1;i<k;i++) {
      polar_orbit.push_back(new polarOrb(m--));
      face[1].push_back(j--);
   }

   if ( is_even(ncon_order) && !point_cut ) {
      polar_orbit.push_back(new polarOrb(m));
      face[1].push_back(j);
   }

   if (add_poles || (((is_even(ncon_order) && point_cut) || !is_even(ncon_order)) && !half_model(longitudes)))
      face[1].push_back(pole[1]->idx);

   if ( strchr(closure.c_str(), 'v') && (face[1].size() > 2)) {  
      add_face(geom, face_list, face[1], -2, -2);
      face_list.back()->rotate = true;
   }  

   // In cases of odd longitudes sides or even longitudes side cuts which are not half models this generates a third side to cover    
   if (!add_poles && 
      ((is_even(ncon_order) && !point_cut) || !is_even(ncon_order)) &&
      (!half_model(longitudes)) &&
      (strchr(closure.c_str(), 'v'))) {
      face[2].push_back(0);
      if (is_even(ncon_order) && !point_cut)
         face[2].push_back(face[1].back());
      face[2].push_back(face[1].front());
      face[2].push_back(face[0].back());
      add_face(geom, face_list, face[2], -2, -2);
      face_list.back()->rotate = true;
   }
}

// untangle polar orbit. patch for n/m
void sort_polar_orbit(col_geom_v &geom, vector<polarOrb *> &polar_orbit)
{
   const vector<vec3d> &verts = geom.verts();

   vec3d v0 = verts[polar_orbit[0]->coord_no];
   int sz = polar_orbit.size();
   vector<pair<double, int> > angles(sz);
   for(int i=0; i<sz; i++) {
      int j = polar_orbit[i]->coord_no;
      angles[i].second = j;
      angles[i].first = rad2deg(angle_around_axis(v0,verts[j],vec3d(0,0,1)));
//fprintf(stderr,"%d = %d  %g\n",i,polar_orbit[i]->coord_no,angles[i].first);
   }
   
   // sort on angles
   sort( angles.begin(), angles.end() );

   for (int i=0; i<sz; i++)
      polar_orbit[i]->coord_no = angles[i].second;
}

void ncon_twist(col_geom_v &geom, vector<polarOrb *> &polar_orbit, vector<coordList *> &coordinates, vector<faceList *> &face_list, 
              vector<edgeList *> &edge_list, int twist, int ncon_order, vector<int> longitudes)
{
   // this function wasn't designed for twist 0
   if (twist == 0)
      return;

   // can't twist when half or less of model is showing
   if (2*longitudes.back() <= longitudes.front())
      return;

   vector<vector<int> > &faces = geom.raw_faces();
   vector<vector<int> > &edges = geom.raw_edges();
   vector<vec3d> &verts = geom.raw_verts();

   double arc = (360.0 / ncon_order) * twist * (+1);
   mat3d rot = mat3d::rot(0, 0, deg2rad(-arc));

   for (unsigned int i=0;i<face_list.size();i++) {
      if ( face_list[i]->rotate ) {
         for (unsigned int j=0;j<faces[i].size();j++) {
            int coord_no = faces[i][j];
            if ( !coordinates[coord_no]->rotated ) {
               verts[coord_no] = rot * verts[coord_no];
               coordinates[coord_no]->rotated = true;
            }
         }
      }
   }

   // Patch - If an open model, some of the polar circle coordinates are not yet rotated
   if ( longitudes.front() > longitudes.back() )
      for (int i=0;i<(int)polar_orbit.size();i++) {
         int j = polar_orbit[i]->coord_no;
         if ( !coordinates[j]->rotated ) {
               verts[j] = rot * verts[j];
               coordinates[j]->rotated = true;
         }
      }

   // patch for n/m models
   sort_polar_orbit(geom, polar_orbit);

   // Create Doubly Circularly Linked List
   for (int i=0;i<(int)polar_orbit.size();i++)
      polar_orbit[i]->forward = i + 1;
   polar_orbit.back()->forward = 0;
   for (int i=(int)polar_orbit.size()-1;i>0;i--)
      polar_orbit[i]->backward = i - 1;
   polar_orbit.front()->backward = polar_orbit.size() - 1;

   for (int k=0;k<(int)polar_orbit.size();k++)
   {
      int p,q;
      p = polar_orbit[k]->coord_no;
      if ( twist > 0 )
         q = polar_orbit[k]->forward;
      else
         q = polar_orbit[k]->backward;
      for (int n=1;n<abs(twist);n++) {
         if ( twist > 0 )
            q = polar_orbit[q]->forward;
         else
            q = polar_orbit[q]->backward;
      }
      q = polar_orbit[q]->coord_no;

      for (unsigned int i=0;i<faces.size();i++)
         if ( !face_list[i]->rotate )
            for (unsigned int j=0;j<faces[i].size();j++)
               if ( faces[i][j] == p )
                  faces[i][j] = q * (-1);

      for (unsigned int i=0;i<edges.size();i++)
         if ( !edge_list[i]->rotate )
            for (unsigned int j=0;j<edges[i].size();j++)
               if ( edges[i][j] == p )
                  edges[i][j] = q * (-1);
   }

   for (unsigned int i=0; i<faces.size(); i++)
      for (unsigned int j=0;j<faces[i].size();j++)
         if ( faces[i][j] < 0 )
            faces[i][j] = abs( faces[i][j] );

   for (unsigned int i=0; i<edges.size(); i++)
      for (unsigned int j=0;j<edges[i].size();j++)
         if ( edges[i][j] < 0 )
            edges[i][j] = abs( edges[i][j] );
}

void find_surface_count(vector<surfaceTable *> &surface_table, int twist, int &surfaces, bool &case2, int &case1_twist)
{
   surfaces = 0;
   case2 = false;
   case1_twist = 0;
   for (unsigned int i=0;i<surface_table.size();i++) {
      if (twist == surface_table[i]->twist) {
         surfaces = surface_table[i]->surfaces;
         case2 = surface_table[i]->case2;
         case1_twist = surface_table[i]->case1_twist;
         break;
      }
   }
}

void build_surface_table(vector<surfaceTable *> &surface_table, int max_twist, int ncon_order, bool point_cut, bool hybrid)
{
   // coding idea furnished by Adrian Rossiter
   int axis_edges = 0;
   int n = ncon_order;
   int t_mod = 0;

   if (!is_even(ncon_order) || hybrid) {
      axis_edges = 1;
      if (hybrid) {
         n*=2;
         t_mod = 1;
      }
   }
   else
   if (is_even(ncon_order) && !point_cut)
      axis_edges = 2;

   for (int twist=2;twist<=max_twist;twist++) {
      int total_surfaces = (int)((gcd(n, 2*twist-t_mod) + axis_edges) / 2);
      //fprintf(stderr,"twist = %d surfaces = %d\n",twist,total_surfaces);
      
      // subtract out discontinuous surfaces for table
      int continuous_surfaces = total_surfaces - axis_edges;
      bool factor = false;
      if (hybrid)
         factor = (n%(2*twist-1) == 0) ? true : false;
      else
         factor = (n%twist == 0) ? true : false;
      int case1_twist = total_surfaces;
      if (is_even(ncon_order) && !point_cut)
         case1_twist--; 

      // only store those with more than minimum surface counts
      if ((continuous_surfaces > 1) || (continuous_surfaces > 0 && (!(is_even(ncon_order) && point_cut) || hybrid))) {
         surface_table.push_back(new surfaceTable(twist,continuous_surfaces,(factor ? twist : case1_twist),(factor ? false : true)));
      }
   }

/*
   for (int twist=2;twist<=max_twist;twist++) {
      int surfaces = 0;
      if (hybrid) {
         if ((ncon_order*2)%(twist*2-1) == 0)
            surfaces = twist - 1;
      }
      else {
         if (ncon_order%twist == 0) {
            if (!is_even(ncon_order))
               surfaces = (int)floor((double)twist/2);
            else if (point_cut) {
               surfaces = twist;
               if (!is_even(ncon_order/twist))
                  surfaces /= 2;
            }
            else if (!point_cut) {
               surfaces = twist - 1;
               if (!is_even(ncon_order/twist))
                  surfaces = twist/2 - 1;
            }
         }
      }
      if (surfaces != 0)
         surface_table.push_back(new surfaceTable(twist,surfaces,twist,false));
   }

   int s = surface_table.size();
   for (int i=s-1;i>=0;i--) {
      int base_twist=surface_table[i]->twist;
      int base_surfaces=surface_table[i]->surfaces;
      int increment = 0;
      if (hybrid)
         increment = 2*base_twist - 1;
      else
         increment = base_twist;
      for (int case2_twist=base_twist+increment;case2_twist<=max_twist;case2_twist+=increment)
      {
         int s;
         bool dummy;
         int dummy2;
         find_surface_count(*surface_table,case2_twist,&s,&dummy,&dummy2);
         if (s==0)
            surface_table.push_back(new surfaceTable(case2_twist,base_surfaces,base_twist,true));
      }
   }
*/
}

void model_info(col_geom_v &geom, bool info)
{
   if (info) {
      unsigned long fsz = geom.faces().size();
      unsigned long vsz = geom.verts().size();
      fprintf(stderr,"The graphical model shown has %lu faces, %lu vertices, and %lu edges\n",
         fsz,vsz,(fsz + vsz - 2));
      fprintf(stderr,"\n");
   }
}

void ncon_info(int ncon_order, bool point_cut, int twist, bool hybrid, bool info, vector<surfaceTable *> &surface_table, surfaceData &sd)
{
   int first, last, forms, chiral, nonchiral, unique;

   if (info) { 
      fprintf(stderr,"\n");
      fprintf(stderr,"This is %s n-icon of order %d twisted %d time%s\n",
         (hybrid ? "a hybrid" : (point_cut ? "point cut" : "side cut")),
         ncon_order,twist,((twist == 1) ? "" : "s"));

      fprintf(stderr,"\n");
      if (hybrid) {
         fprintf(stderr,"This is an order %d, hybrid n-icon. By hybrid, it means it is one half a\n",ncon_order);
         fprintf(stderr,"point cut n-icon joined to half of a side cut n-icon. Hybrid n-icons are\n");
         fprintf(stderr,"of even order and have at least one discontinuous surface and one\n");
         fprintf(stderr,"discontinuous edge. As a non-faceted smooth model it would be self dual.\n");
         fprintf(stderr,"Circuit patterns are self dual.\n");
      }
      else
      if (is_even(ncon_order) && point_cut) {
         fprintf(stderr,"The order of this n-icon, %d, is even. It is, what is termed, a Point Cut\n",ncon_order);
         fprintf(stderr,"Even order, point cut n-icons have at least one continuous surface and\n");
         fprintf(stderr,"two discontinuous edges. As a non-faceted smooth model it would be dual\n");
         fprintf(stderr,"to the Side Cut n-icon. Circuit patterns are dual to Side Cut n-icons.\n");
      }
      else
      if (is_even(ncon_order) && !point_cut) {
         fprintf(stderr,"The order of this n-icon, %d, is even. It is, what is termed, a Side Cut\n",ncon_order);
         fprintf(stderr,"Even order, side cut n-icons have at least two discontinuous surfaces and\n");
         fprintf(stderr,"one continuous edge. As a non-faceted smooth model it would be dual to the\n");
         fprintf(stderr,"Point Cut n-icon. Circuit patterns are dual to Point Cut n-icons.\n");
      }
      else
      if (!is_even(ncon_order)) {
         fprintf(stderr,"The order of this n-icon, %d, is odd. It could be termed both a Point Cut\n",ncon_order);
         fprintf(stderr,"and a Side Cut depending on the poles. Odd order n-icons have at least one\n");
         fprintf(stderr,"discontinuous surface and one discontinuous edge. As a non-faceted smooth\n");
         fprintf(stderr,"model it would be self dual. Circuit patterns are self dual.\n");
      }
   }

   if (info && twist == 0) {
      fprintf(stderr,"\n");
      fprintf(stderr,"Since it is twisted 0 times, it is as a globe with %d latitudes\n", num_lats(ncon_order, point_cut));
   
      if (is_even(ncon_order) && point_cut) {
         fprintf(stderr,"It has a north and south pole\n");
         fprintf(stderr,"It is the polyhedral dual of the twist 0 side cut order %d n-icon\n",ncon_order);
      }
      else
      if (is_even(ncon_order) && !point_cut) {
         fprintf(stderr,"It has two polar caps\n");
         fprintf(stderr,"It is the polyhedral dual of the twist 0 point cut order %d n-icon\n",ncon_order);
      }
      else
      if (!is_even(ncon_order)) {
         fprintf(stderr,"It has a north polar cap and south pole\n");
         fprintf(stderr,"As a polyhedra it is self dual\n");
      }
   }
   
   if (info && hybrid) {
      fprintf(stderr,"\n");
      fprintf(stderr,"Note: Hybrids have no twist 0\n");
   }

   int base_form = 1;
   if (hybrid)
      base_form--;

   last = (int)floor((double)ncon_order/2);
   first = -last;
   if (is_even(ncon_order) && !hybrid)
      first++;

   if (info) {
      fprintf(stderr,"\n");
      fprintf(stderr,"Using coloring, there are %d distinct twists ranging from %d through +%d\n",
         ncon_order,first,last);
   }
   
   int mod_twist = twist%ncon_order;
   if (hybrid && mod_twist == 0)
      mod_twist = -1;
   int color_twist = mod_twist%last;
   if (is_even(ncon_order) && color_twist == -last)
      color_twist = abs(color_twist);
   if (mod_twist != color_twist) {
      if (mod_twist > last)
         color_twist = first + color_twist - 1;
      else
      if (mod_twist < first)
         color_twist = last + color_twist;
      else
         color_twist = mod_twist;
   }
   
   if (info) {
      if (color_twist != twist)
         fprintf(stderr,"This twist is the same as color twist %d\n", color_twist);
   }

   nonchiral = 0;
   if (is_even(ncon_order)) {
      if (hybrid)
         last = (int)floor((double)(ncon_order+2)/4);
      else
         last = (int)floor((double)ncon_order/4);
      forms = (last*2)+1;
      first = -last;
      if ((!hybrid && ncon_order%4 == 0) || (hybrid && (ncon_order+2)%4 == 0)) {
         forms--;
         first++;
      }
      // no twist 0
      if (hybrid)
         forms--;
   }
   else {
      forms = ncon_order;
      last = (int)floor((double)forms/2);
      first = -last;
   }

   chiral = last;
   if ((!hybrid && ncon_order%4 == 0) || (hybrid && (ncon_order+2)%4 == 0)) {
      nonchiral = 1;
      chiral--;
   }
   unique = forms-last-base_form+nonchiral;

   if (info) {
      fprintf(stderr,"\n");
      fprintf(stderr,"Without coloring, there are %d distinct twists ranging from %d through +%d\n",forms,first,last);
      fprintf(stderr,"   %d base model(twist 0)\n",base_form);
      fprintf(stderr,"   %d form%s chiral and %s mirror image%s\n",
         chiral,((chiral == 1) ? " is" : "s are"),((chiral == 1) ? "has a" : "have"),((chiral == 1) ? "" : "s"));
      fprintf(stderr,"   %d form%s non-chiral\n",
         nonchiral,((nonchiral == 1) ? " is" : "s are"));
      fprintf(stderr,"   %d unique form%s of order %d twisted n-icon%s\n\n",
         unique,((unique == 1) ? "" : "s"),ncon_order,((unique == 1) ? "" : "s"));
   }

   int tmp_forms = forms;
   if (is_even(ncon_order))
      tmp_forms = ncon_order/2;

   // also correct for no base twist 0 in hybrids
   int base_twist = twist;
   if (base_twist > 0)
      while(base_twist > last) {
         base_twist-=tmp_forms;
         if (hybrid && (base_twist <= 0 && base_twist+tmp_forms > 0))
            base_twist--;
      }
   else
   if (base_twist < 0)
      while(base_twist < first) {
         base_twist+=tmp_forms;
         if (hybrid && (base_twist >= 0 && base_twist-tmp_forms < 0))
            base_twist++;
      }

   // using abs(base_twist) so much it might as well be stored
   int posi_twist = abs(base_twist);

   sd.nonchiral = (((chiral == 0 || base_twist == 0 ||
      (!hybrid && is_even(ncon_order) && (ncon_order%4 == 0) && base_twist == last) ||
      ((hybrid && (ncon_order+2)%4 == 0) && base_twist == last))) ? true : false);

   if (info) {
      if (twist < first || twist > last)
         fprintf(stderr,"It is the same as an n-icon of order %d twisted %d %s\n",
            ncon_order,base_twist,((base_twist == 1) ? "time" : "times"));

      if (sd.nonchiral)
         fprintf(stderr,"It is not chiral\n");
      else
         fprintf(stderr,"It is the mirror image of an n-icon of order %d twisted %d time%s\n",
            ncon_order,-base_twist,((-base_twist == 1) ? "" : "s"));
   }

   if (info) {
      fprintf(stderr,"\n");
      double twist_angle = ((double)360/ncon_order)*twist;
      if (hybrid) {
         if ( twist_angle > 0 )
            twist_angle -= (double)360/(ncon_order*2);
         else
            twist_angle += (double)360/(ncon_order*2);
      }
      fprintf(stderr,"It has an absolute twist of %lf degrees, %lf radians\n",
         twist_angle, deg2rad(twist_angle));

      if ( twist != base_twist ) {
         twist_angle = ((double)360/ncon_order)*base_twist;
         if (hybrid) {
            if ( twist_angle > 0 )
               twist_angle -= (double)360/(ncon_order*2);
            else
               twist_angle += (double)360/(ncon_order*2);
         }
         fprintf(stderr,"which is the same as one at %lf degrees, %lf radians\n",
            twist_angle, deg2rad(twist_angle));
      }
   }

   // if ncon_info() is called multiple times, this keeps surface table from being called only once
   if (!surface_table.size())
      build_surface_table(surface_table, last, ncon_order, point_cut, hybrid);

   sd.c_surfaces = 0;
   sd.c_edges = 0;
   sd.d_surfaces = 0;
   sd.d_edges = 0;
   sd.total_surfaces = 0;
   int total_edges = 0;

   sd.ncon_case2 = false;
   int case1_twist = 0;

   if (posi_twist == 0) {
      sd.c_surfaces = (int)floor((double)ncon_order+1)/2;
      if (is_even(ncon_order) && !point_cut)
         sd.c_surfaces++;
      sd.c_edges = sd.c_surfaces-1;

      if (!is_even(ncon_order) || !point_cut) {
         sd.c_surfaces--;
         sd.d_surfaces++;
         if (is_even(ncon_order)) {
            sd.c_surfaces--;
            sd.d_surfaces++;
         }
      }
      sd.d_edges = 0;
   }
   else
   if (hybrid) {
      find_surface_count(surface_table, posi_twist, sd.c_surfaces, sd.ncon_case2, case1_twist);      
      sd.d_surfaces = 1;
      sd.c_edges = sd.c_surfaces + sd.d_surfaces - 1;
      sd.d_edges = 1;
   }
   else
   if (is_even(ncon_order) && point_cut) {
      find_surface_count(surface_table, posi_twist, sd.c_surfaces, sd.ncon_case2, case1_twist);
      if (sd.c_surfaces == 0)
         sd.c_surfaces = 1;
      sd.d_surfaces = 0;
      sd.c_edges = sd.c_surfaces + sd.d_surfaces - 1;
      sd.d_edges = 2;
   }
   else
   if (is_even(ncon_order) && !point_cut) {
      find_surface_count(surface_table, posi_twist, sd.c_surfaces, sd.ncon_case2, case1_twist);
      sd.d_surfaces = 2;
      sd.c_edges = sd.c_surfaces + sd.d_surfaces - 1;
      sd.d_edges = 0;
   }
   else
   if (!is_even(ncon_order)) {
      find_surface_count(surface_table, posi_twist, sd.c_surfaces, sd.ncon_case2, case1_twist);
      sd.c_edges = sd.c_surfaces;
      sd.d_surfaces = 1;
      sd.d_edges = 1;
   }
   sd.total_surfaces = sd.c_surfaces + sd.d_surfaces;
   total_edges = sd.c_edges + sd.d_edges;

   if (info) {
      fprintf(stderr,"\n");
      fprintf(stderr,"Treated as conic, it has a total of %d surface%s and a total of %d edge%s\n",
         sd.total_surfaces,((sd.total_surfaces == 1) ? "" : "s"),
         sd.c_edges+sd.d_edges,((sd.c_edges+sd.d_edges == 1) ? "" : "s"));
      fprintf(stderr,"   %d surface%s continuous\n",
         sd.c_surfaces,((sd.c_surfaces == 1) ? " is" : "s are"));
      fprintf(stderr,"   %d surface%s discontinuous\n",
         sd.d_surfaces,((sd.d_surfaces == 1) ? " is" : "s are"));
      fprintf(stderr,"   %d edge%s continuous\n",
         sd.c_edges,((sd.c_edges == 1) ? " is" : "s are"));
      fprintf(stderr,"   %d edge%s discontinuous\n",
         sd.d_edges,((sd.d_edges == 1) ? " is" : "s are"));
   }

   if (info) {
      if ((sd.c_surfaces > 1) || (sd.c_surfaces > 0 && (!(is_even(ncon_order) && point_cut) || hybrid))) {
         fprintf(stderr,"\n");
         if (!sd.ncon_case2) {
            fprintf(stderr,"This is a Case 1 n-icon. Surfaces cannot be colored based on an earlier twist\n");
         }
         else {
            char ntype = 'p';
            if (hybrid)
               ntype = 'h';
            else
            if (is_even(ncon_order) && !point_cut)
               ntype = 's';
            fprintf(stderr,"This is a Case 2 n-icon. Surfaces can be colored based on an earlier twist\n");
            char sign = (twist > 0) ? '+' : '-';
            fprintf(stderr,"It can be derived by twisting N%d%cT%d%c %c%d increments\n",
               ncon_order,sign,case1_twist,ntype,sign,posi_twist - case1_twist);
         }
      }
      fprintf(stderr,"\n");
   }
}

void color_uncolored_faces(col_geom_v &geom, col_val default_color)
{
   for (unsigned int i=0;i<geom.faces().size();i++) {
      if (!(geom.get_f_col(i)).is_set())
         geom.set_f_col(i,default_color);
   }
}

void color_uncolored_edges(col_geom_v &geom, vector<edgeList *> &edge_list, vector<poleList *> &pole, col_val default_color)
{
   const vector<vector<int> > &edges = geom.edges();
   
   for (unsigned int i=0;i<edge_list.size();i++) {
      int j = edge_list[i]->edge_no;
      if (!(geom.get_e_col(j)).is_set()) {
         geom.set_e_col(j,default_color);
         int v1 = edges[j][0];
         int v2 = edges[j][1];
         if (!(geom.get_v_col(v1)).is_set())
            geom.set_v_col(v1,default_color);
         if (!(geom.get_v_col(v2)).is_set())
            geom.set_v_col(v2,default_color);
      }
   }

   for (unsigned int i=0;i<2;i++) {
      if (pole[i]->idx > -1) {
         if (!(geom.get_v_col(pole[i]->idx)).is_set())
            geom.set_v_col(pole[i]->idx,default_color);
      }
   }
}

void color_unused_edges(col_geom_v &geom, col_val unused_edge_color)
{
   geom.add_missing_impl_edges();   
   for (unsigned int i=0;i<geom.edges().size();i++) {
      if (!(geom.get_e_col(i)).is_set())
         geom.set_e_col(i,unused_edge_color);
   }
   
   for (unsigned int i=0;i<geom.verts().size();i++) {
      if (!(geom.get_v_col(i)).is_set())
         geom.set_v_col(i,unused_edge_color);
   }
}

void reassert_colored_edges(col_geom_v &geom, col_val default_color)
{
   const vector<vector<int> > &edges = geom.edges();
   
   for (unsigned int i=0;i<geom.edges().size();i++) {
      col_val c = geom.get_e_col(i);
      // if its index or it is a value and is completely opaque
      if (c.is_set() && (c.is_idx() || !c.get_trans())) {
         int v1 = edges[i][0];
         int v2 = edges[i][1];
         col_val cv1 = geom.get_v_col(v1);
         col_val cv2 = geom.get_v_col(v2);
         // let default color be dominant
         if (cv1.is_set() && cv1 != default_color)
            geom.set_v_col(v1,c);
         if (cv2.is_set() && cv2 != default_color)
            geom.set_v_col(v2,c);
      }
   }
}

col_val set_alpha(col_val c, int a)
{  
   return col_val(c[0],c[1],c[2],a);
}

// safe transparency setting. don't allow alpha set if
// color is index
// color is invisible
// -T or -U have not been set
void set_vert_color(col_geom_v &geom, int i, col_val c, int opacity)
{
   if (c.is_val() && !c.is_inv() && opacity != -1)
      c = set_alpha(c,opacity);
   geom.set_v_col(i,c);
}

void set_edge_col(col_geom_v &geom, int i, col_val c, int opacity)
{
   if (c.is_val() && !c.is_inv() && opacity != -1)
      c = set_alpha(c,opacity);
   geom.set_e_col(i,c);
}

void set_face_color(col_geom_v &geom, int i, col_val c, int opacity)
{
   if (c.is_val() && !c.is_inv() && opacity != -1)
      c = set_alpha(c,opacity);
   geom.set_f_col(i,c);
}

void set_edge_and_verts_col(col_geom_v &geom, int i, col_val c, int opacity)
{
   set_edge_col(geom,i,c,opacity);

   set_vert_color(geom,geom.edges()[i][0],c,opacity);
   set_vert_color(geom,geom.edges()[i][1],c,opacity);
}

void set_edge_color(col_geom_v &geom, int i, col_val c, int opacity)
{
   set_edge_and_verts_col(geom,i,c,opacity);
}

void set_shell_indent_edges_invisible(col_geom_v &geom, vector<edgeList *> &edge_list, vector<poleList *> &pole, int ncon_order, bool point_cut, bool radius_reverse)
{
   int n = ncon_order/2;

   for (unsigned int i=0;i<edge_list.size();i++) {
      int j = edge_list[i]->edge_no;
      int lat = edge_list[i]->lat;
      bool set_invisible = (is_even(n) && point_cut && !is_even(lat)) || (is_even(n) && !point_cut && is_even(lat)) || (!is_even(n) && is_even(lat));
      if (radius_reverse) {
         set_invisible = (set_invisible) ? false : true;
fprintf(stderr,"doing it\n");
}
      if (set_invisible)
         set_edge_color(geom,j,col_val::invisible,255);
   }

   for (unsigned int i=0;i<2;i++) {
      if (pole[i]->idx > -1) {
         int lat = pole[i]->lat;
         bool set_invisible = (is_even(n) && point_cut && !is_even(lat)) || (is_even(n) && !point_cut && is_even(lat)) || (!is_even(n) && is_even(lat));
         if (radius_reverse)
            set_invisible = (set_invisible) ? false : true;
         if (set_invisible)
            set_vert_color(geom,pole[i]->idx,col_val::invisible,255);
      }
   }
}

void ncon_edge_coloring(col_geom_v &geom, vector<edgeList *> &edge_list, vector<poleList *> &pole, char edge_coloring_method,
                        color_map_multi edge_map, int edge_opacity, string edge_pattern,
                        map<int, pair<int, int> > &edge_color_table, bool point_cut, bool hybrid)
{
   const vector<vector<int> > &edges = geom.edges();
   const vector<vec3d> &verts = geom.verts();

   int opq = 255;

   if ( edge_coloring_method == 's' ) {
      for (unsigned int i=0;i<edge_list.size();i++) {
         int j = edge_list[i]->edge_no;
         int lat = edge_list[i]->lat;

         if ( edge_list[i]->lat < 0 ) {
            set_edge_color(geom,j,col_val(),edge_opacity);
         }
         else {
            int col_idx = 0;
            if ( edge_list[i]->rotate || (hybrid && point_cut)) // front side
               col_idx = edge_color_table[lat].second;
            else
               col_idx = edge_color_table[lat].first;

            opq = edge_pattern[col_idx%edge_pattern.size()] == '1' ? edge_opacity : 255;
            set_edge_color(geom,j,edge_map.get_col(col_idx),opq);
         }
      }
      
      for (unsigned int i=0;i<2;i++) {
         if (pole[i]->idx > -1) {
            int lat = pole[i]->lat;
            int col_idx = edge_color_table[lat].second;
            opq = edge_pattern[col_idx%edge_pattern.size()] == '1' ? edge_opacity : 255;;
            set_vert_color(geom,pole[i]->idx,edge_map.get_col(col_idx),opq);
         }
      }
   }
   else
   if ( edge_coloring_method == 'l' ) {
      for (unsigned int i=0;i<edge_list.size();i++) {
         int j = edge_list[i]->edge_no; 
         if ( edge_list[i]->lat < 0 ) {
            set_edge_color(geom,j,col_val(),edge_opacity);
         }
         else {
            int lat = edge_list[i]->lat;
            opq = edge_pattern[lat%edge_pattern.size()] == '1' ? edge_opacity : 255;
            set_edge_color(geom,j,edge_map.get_col(lat),opq);
         }
      }
           
      for (unsigned int i=0;i<2;i++) {
         if (pole[i]->idx > -1) {
            int lat = pole[i]->lat;
            opq = edge_pattern[lat%edge_pattern.size()] == '1' ? edge_opacity : 255;
            set_vert_color(geom,pole[i]->idx,edge_map.get_col(lat),opq);
         }
      }
   }
   else
   if ( edge_coloring_method == 'm' ) {
      for (unsigned int i=0;i<edge_list.size();i++) {
         int j = edge_list[i]->edge_no;
         if ( edge_list[i]->lon < 0 ) {
            set_edge_color(geom,j,col_val(),edge_opacity);
         }
         else {
            int lon = edge_list[i]->lon;
            opq = edge_pattern[lon%edge_pattern.size()] == '1' ? edge_opacity : 255;
            set_edge_color(geom,j,edge_map.get_col(lon),opq);
         }
      }
      
      // poles don't have any longitude
      for (unsigned int i=0;i<2;i++) {
         if (pole[i]->idx > -1) {
            set_vert_color(geom,pole[i]->idx,col_val(DEFAULT_COLOR),edge_opacity);
         }
      }
   }
   else
   if ( edge_coloring_method == 'c' ) {
      for (unsigned int i=0;i<edge_list.size();i++) {
         int j = edge_list[i]->edge_no;
         if ( edge_list[i]->lat < 0 )
            set_edge_color(geom,j,col_val(),edge_opacity);
         else {
            int n = -1;
            if ( (is_even(edge_list[i]->lat) && is_even(edge_list[i]->lon)) ||
                 (!is_even(edge_list[i]->lat) && !is_even(edge_list[i]->lon)) )
               n = 0;
            else
               n = 1;

            opq = edge_pattern[n%edge_pattern.size()] == '1' ? edge_opacity : 255;
            set_edge_color(geom,j,edge_map.get_col(n),opq);
         }
      }
      
      // poles will be colored based North/South
      for (unsigned int i=0;i<2;i++) {
         if (pole[i]->idx > -1) {
            opq = edge_pattern[i%edge_pattern.size()] == '1' ? edge_opacity : 255;
            set_vert_color(geom,pole[i]->idx,edge_map.get_col(i),opq);
         }
      }
   }
   else
   if ( edge_coloring_method == 'n' ) {
      // keep track of index beyond loop
      unsigned int k = 0;
      for (unsigned int i=0;i<edge_list.size();i++) {
         int j = edge_list[i]->edge_no;
         opq = edge_pattern[i%edge_pattern.size()] == '1' ? edge_opacity : 255;
         set_edge_color(geom,j,edge_map.get_col(i),opq);
         k++;
      }
      
      for (unsigned int i=0;i<2;i++) {
         int col_idx = k;
         if (pole[i]->idx > -1) {
            opq = edge_pattern[k%edge_pattern.size()] == '1' ? edge_opacity : 255;
            set_vert_color(geom,pole[i]->idx,edge_map.get_col(col_idx),opq);
         }
         k++;
      }
   }
   else
   if ( edge_coloring_method == 'x' ||  edge_coloring_method == 'y' || edge_coloring_method == 'z' ) {
      for (unsigned int i=0;i<edge_list.size();i++) {
         int j = edge_list[i]->edge_no;
         double d = 0.0;
         for (unsigned int k=0;k<edges[j].size();k++) {
            if ( edge_coloring_method == 'x' )
               d += verts[edges[j][k]][0];
            else
            if ( edge_coloring_method == 'y' )
               d += verts[edges[j][k]][1];
            else
            if ( edge_coloring_method == 'z' )
               d += verts[edges[j][k]][2];
            }

         // The hybrid base portion will end up in +Z
         if (hybrid && point_cut)
            d *= -1;

         col_val c;
         int n = -1;
         if ( d > 0.0 )
            n = 0;
         else
         if ( d < 0.0 )
            n = 1;
         else {
            n = -1;
            opq = edge_opacity;
            c = col_val(DEFAULT_COLOR);
         }

         if (n > -1) {
            opq = edge_pattern[n%edge_pattern.size()] == '1' ? edge_opacity : 255;
            c = edge_map.get_col(n);
         }
         set_edge_color(geom,j,c,opq);
      }
      
      // poles are on the axis
      for (unsigned int i=0;i<2;i++) {
         if (pole[i]->idx > -1) {
            set_vert_color(geom,pole[i]->idx,col_val(DEFAULT_COLOR),edge_opacity);
         }
      }
   }
   else
   if ( edge_coloring_method == 'o' ) {
      for (unsigned int i=0;i<edge_list.size();i++) {
         int j = edge_list[i]->edge_no;
         double dx = 0.0;
         double dy = 0.0;
         double dz = 0.0;
         for (unsigned int k=0;k<edges[j].size();k++) {
            dx += verts[edges[j][k]][0];
            dy += verts[edges[j][k]][1];
            dz += verts[edges[j][k]][2];
         }

         // The hybrid base portion will end up in +Z
         if (hybrid && point_cut)
            dz *= -1;

         col_val c;
         int n = -1;
         // by octant number 1 to 8
         if (( dx > 0.0 ) && ( dy > 0.0 ) && ( dz > 0.0 ))
            n = 0;
         else
         if (( dx < 0.0 ) && ( dy > 0.0 ) && ( dz > 0.0 ))
            n = 1;
         else
         if (( dx < 0.0 ) && ( dy < 0.0 ) && ( dz > 0.0 ))
            n = 2;
         else
         if (( dx > 0.0 ) && ( dy < 0.0 ) && ( dz > 0.0 ))
            n = 3;
         else
         if (( dx > 0.0 ) && ( dy > 0.0 ) && ( dz < 0.0 ))
            n = 4;
         else
         if (( dx < 0.0 ) && ( dy > 0.0 ) && ( dz < 0.0 ))
            n = 5;
         else
         if (( dx < 0.0 ) && ( dy < 0.0 ) && ( dz < 0.0 ))
            n = 6;
         else
         if (( dx > 0.0 ) && ( dy < 0.0 ) && ( dz < 0.0 ))
            n = 7;
         else {
            n = -1;
            opq = edge_opacity;
            c = col_val(DEFAULT_COLOR);
         }

         if (n > -1) {
            opq = edge_pattern[n%edge_pattern.size()] == '1' ? edge_opacity : 255;
            c = edge_map.get_col(n);
         }
         set_edge_color(geom,j,c,opq);
      }

      // poles are on the axis
      for (unsigned int i=0;i<2;i++) {
         if (pole[i]->idx > -1) {
            set_vert_color(geom,pole[i]->idx,col_val(DEFAULT_COLOR),edge_opacity);
         }
      }
   }
}

void ncon_face_coloring(col_geom_v &geom, vector<faceList *> &face_list, char face_coloring_method,
                        color_map_multi face_map, int face_opacity, string face_pattern,
                        map<int, pair<int, int> > &face_color_table, bool point_cut, bool hybrid)
{
   const vector<vector<int> > &faces = geom.faces();
   const vector<vec3d> &verts = geom.verts();

   int opq = 255;
   
   if ( face_coloring_method == 's' ) {
      for (unsigned int i=0;i<face_list.size();i++) {
         int j = face_list[i]->face_no;
         int lat = face_list[i]->lat;

         if ( face_list[i]->lat < 0 ) {
            set_face_color(geom,j,col_val(),face_opacity);
         }
         else {
            int col_idx = 0;
            if ( face_list[i]->rotate || (hybrid && point_cut)) // front side
               col_idx = face_color_table[lat].second;
            else
               col_idx = face_color_table[lat].first;

            opq = face_pattern[col_idx%face_pattern.size()] == '1' ? face_opacity : 255;
            set_face_color(geom,j,face_map.get_col(col_idx),opq);
         }
      }
   }
   else
   if ( face_coloring_method == 'l' ) {
      for (unsigned int i=0;i<face_list.size();i++) {
         int j = face_list[i]->face_no; 
         if ( face_list[i]->lat < 0 ) {
            set_face_color(geom,j,col_val(),face_opacity);
         }
         else {
            int lat = face_list[i]->lat;;
            opq = face_pattern[lat%face_pattern.size()] == '1' ? face_opacity : 255;
            set_face_color(geom,j,face_map.get_col(lat),opq);
         }
      }
   }
   else
   if ( face_coloring_method == 'm' ) {
      for (unsigned int i=0;i<face_list.size();i++) {
         int j = face_list[i]->face_no;
         if ( face_list[i]->lon < 0 ) {
            set_face_color(geom,j,col_val(),face_opacity);
         }
         else {
            int lon = face_list[i]->lon;
            opq = face_pattern[lon%face_pattern.size()] == '1' ? face_opacity : 255;
            set_face_color(geom,j,face_map.get_col(lon),opq);
         }
      }
   }
   else
   if ( face_coloring_method == 'c' ) {
      for (unsigned int i=0;i<face_list.size();i++) {
         int j = face_list[i]->face_no;
         if ( face_list[i]->lat < 0 )
            set_face_color(geom,j,col_val(),face_opacity);
         else {
            int n = -1;
            if ( (is_even(face_list[i]->lat) && is_even(face_list[i]->lon)) ||
                 (!is_even(face_list[i]->lat) && !is_even(face_list[i]->lon)) )
               n = 0;
            else
               n = 1;

            opq = face_pattern[n%face_pattern.size()] == '1' ? face_opacity : 255;
            set_face_color(geom,j,face_map.get_col(n),opq);
         }
      }
   }
   else
   if ( face_coloring_method == 'n' ) {
      for (unsigned int i=0;i<face_list.size();i++) {
         int j = face_list[i]->face_no;
         opq = face_pattern[i%face_pattern.size()] == '1' ? face_opacity : 255;
         set_face_color(geom,j,face_map.get_col(i),opq);
      }
   }
   else
   if ( face_coloring_method == 'x' ||  face_coloring_method == 'y' || face_coloring_method == 'z' ) {
      for (unsigned int i=0;i<face_list.size();i++) {
         int j = face_list[i]->face_no;
         double d = 0.0;
         for (unsigned int k=0;k<faces[j].size();k++) {
            if ( face_coloring_method == 'x' )
               d += verts[faces[j][k]][0];
            else
            if ( face_coloring_method == 'y' )
               d += verts[faces[j][k]][1];
            else
            if ( face_coloring_method == 'z' )
               d += verts[faces[j][k]][2];
            }

         // The hybrid base portion will end up in +Z
         if (hybrid && point_cut)
            d *= -1;

         col_val c;
         int n = -1;
         if ( d > 0.0 )
            n = 0;
         else
         if ( d < 0.0 )
            n = 1;
         else {
            n = -1;
            opq = face_opacity;
            c = col_val(DEFAULT_COLOR);
         }

         if (n > -1) {
            opq = face_pattern[n%face_pattern.size()] == '1' ? face_opacity : 255;
            c = face_map.get_col(n);
         }
         set_face_color(geom,j,c,opq);
      }
   }
   else
   if ( face_coloring_method == 'o' ) {
      for (unsigned int i=0;i<face_list.size();i++) {
         int j = face_list[i]->face_no;
         double dx = 0.0;
         double dy = 0.0;
         double dz = 0.0;
         for (unsigned int k=0;k<faces[j].size();k++) {
            dx += verts[faces[j][k]][0];
            dy += verts[faces[j][k]][1];
            dz += verts[faces[j][k]][2];
         }

         // The hybrid base portion will end up in +Z
         if (hybrid && point_cut)
            dz *= -1;

         // by octant number 1 to 8
         col_val c;
         int n = -1;
         if (( dx > 0.0 ) && ( dy > 0.0 ) && ( dz > 0.0 ))
            n = 0;
         else
         if (( dx < 0.0 ) && ( dy > 0.0 ) && ( dz > 0.0 ))
            n = 1;
         else
         if (( dx < 0.0 ) && ( dy < 0.0 ) && ( dz > 0.0 ))
            n = 2;
         else
         if (( dx > 0.0 ) && ( dy < 0.0 ) && ( dz > 0.0 ))
            n = 3;
         else
         if (( dx > 0.0 ) && ( dy > 0.0 ) && ( dz < 0.0 ))
            n = 4;
         else
         if (( dx < 0.0 ) && ( dy > 0.0 ) && ( dz < 0.0 ))
            n = 5;
         else
         if (( dx < 0.0 ) && ( dy < 0.0 ) && ( dz < 0.0 ))
            n = 6;
         else
         if (( dx > 0.0 ) && ( dy < 0.0 ) && ( dz < 0.0 ))
            n = 7;
         else {
            n = -1;
            opq = face_opacity;
            c = col_val(DEFAULT_COLOR);
         }

         if (n > -1) {
            opq = face_pattern[n%face_pattern.size()] == '1' ? face_opacity : 255;
            c = face_map.get_col(n);
         }
         set_face_color(geom,j,c,opq);
      }
   }
}

void make_sequential_map(map<int, pair<int, int> > &color_table)
{
   unsigned int sz = color_table.size();
   
   // add two to avoid zeroes and make all elements negative.
   for (unsigned int i=0;i<sz;i++) {
      color_table[i].first = -(color_table[i].first + 2);
      color_table[i].second = -(color_table[i].second + 2);
   }
   
   int k = 0;
   for (unsigned int i=0;i<sz;i++) {
      int element = color_table[i].first;
      if (element < 0) {
         for (unsigned int j=i;j<sz;j++) {
            if (color_table[j].first == element)
               color_table[j].first = k;
         }
         for (unsigned int j=0;j<sz;j++) {
            if (color_table[j].second == element)
               color_table[j].second = k;
         }
         k++;
      }
   }
}

// coding idea for circuit coloring furnished by Adrian Rossiter
void build_circuit_table(int n, int twist, bool hybrid, bool symmetric_coloring, map<int,int> &circuit_table)
{
   // use a double size polygon
   n *= 2;
   int t = 2*twist;
 
   if (hybrid)
      t--;
   
   int t_mult = symmetric_coloring ? 1 : 2; 

   int d = gcd(n,t_mult*t);

   for(int i=0; i<=d/2; i++) {
      for(int j=i; j<n; j+=d)
         circuit_table[j] = i;
      for(int j=d-i; j<n; j+=d)
         circuit_table[j] = i;
   }
}

void build_color_tables(map<int, pair<int, int> > &edge_color_table, map<int, pair<int, int> > &face_color_table,
                        int n, bool point_cut, bool hybrid, int twist,
                        bool symmetric_coloring, bool face_sequential_colors, bool edge_sequential_colors)
{
   // for the color tables twist is made positive. The positive twist colors work for negative twists.
   twist = abs(twist);
   
   map<int,int> circuit_table;
   build_circuit_table(n, twist, hybrid, symmetric_coloring, circuit_table);

   // first is back half.
   // second is front half.
   
   // build edge table
   int inc = ((is_even(n) && point_cut) & !hybrid) ? 0 : -1;
      
   for (int i=0;i<n;i++) {
      // when inc = -1, optional north pole will not be correct color if pos < 0. Make it 0 in that case
      int pos = 2*i+inc;
      pos = (pos < 0) ? 0 : pos;
      int cir_idx = pos%circuit_table.size();
      edge_color_table[i].first = circuit_table[cir_idx];
   }
   
   // build face table
   inc = ((is_even(n) && point_cut) & !hybrid) ? 1 : 0;

   for (int i=0;i<n;i++) {
      int cir_idx = (2*i+inc)%circuit_table.size();
      face_color_table[i].first = circuit_table[cir_idx];
   }

   for (int i=0;i<n;i++) {
      edge_color_table[i].second = edge_color_table[(twist+i)%n].first;
      face_color_table[i].second = face_color_table[(twist+i)%n].first;
   }
   
   if (edge_sequential_colors)
      make_sequential_map(edge_color_table);
   if (face_sequential_colors)
      make_sequential_map(face_color_table);
}

void ncon_coloring(col_geom_v &geom, vector<faceList *> &face_list, vector<edgeList *> &edge_list, vector<poleList *> &pole, ncon_opts &opts)
{
   int twist = (opts.nm_shell && opts.hybrid) ? opts.twist/2 : opts.twist;

   map<int, pair<int, int> > edge_color_table;
   map<int, pair<int, int> > face_color_table;
   if (opts.face_coloring_method == 's' || opts.edge_coloring_method == 's')
      build_color_tables(edge_color_table, face_color_table,
                         opts.ncon_order, opts.point_cut, opts.hybrid, twist,
                         opts.symmetric_coloring, opts.face_sequential_colors, opts.edge_sequential_colors);
   
   if (opts.face_coloring_method)
      ncon_face_coloring(geom, face_list, opts.face_coloring_method,
                         opts.face_map, opts.face_opacity, opts.face_pattern,
                         face_color_table, opts.point_cut, opts.hybrid);

   if (opts.edge_coloring_method)
      ncon_edge_coloring(geom, edge_list, pole, opts.edge_coloring_method,
                         opts.edge_map, opts.edge_opacity, opts.edge_pattern,
                         edge_color_table, opts.point_cut, opts.hybrid);

   if (opts.face_coloring_method)
      color_uncolored_faces(geom, col_val(DEFAULT_COLOR));
   if (opts.edge_coloring_method)
      color_uncolored_edges(geom, edge_list, pole, col_val(DEFAULT_COLOR));
}
      
double hybrid_twist_angle(int n, int t, bool nm_shell)
{
   if (nm_shell) {
      n /= 2;
      t /= 2;
   }

   // there is no twist 0 so positive twists have to be adjusted
   int twist = t;
   if ( twist > 0 )
      twist -= 1;

   // twist 1, 2, 3 is really twist 0.5, 1.5, 2.5 so add half twist
   double half_twist = (double)360/(n*2);
   double angle = half_twist + half_twist*twist*2;

   return angle;
}

void build_globe(col_geom_v &geom, vector<coordList *> &coordinates, vector<faceList *> &face_list, vector<edgeList *> &edge_list, vector<poleList *> &pole, ncon_opts &opts)
{
   bool point_cut_save = opts.point_cut;
   vector<int> prime_meridian;
   build_prime_meridian(geom, prime_meridian, coordinates, opts.ncon_order, opts.d, opts.nm_shell, opts.inner_radius, opts.outer_radius, opts.point_cut);
   form_globe(geom, prime_meridian, coordinates, face_list, edge_list, opts.edge_coloring_method,
              opts.ncon_order, opts.point_cut, opts.longitudes, opts.closure, opts.half_model_marker);
   add_caps(geom, coordinates, face_list, pole, opts.ncon_order, opts.point_cut, opts.hybrid, opts.longitudes,
      opts.split, opts.add_poles, opts.hide_elems);
   ncon_coloring(geom, face_list, edge_list, pole, opts);
   if (opts.nm_shell && opts.hide_indent)
      set_shell_indent_edges_invisible(geom, edge_list, pole, opts.ncon_order, point_cut_save, (opts.inner_radius > opts.outer_radius));
}

void process_hybrid(col_geom_v &geom, ncon_opts &opts)
{
   vector<coordList *> coordinates;
   vector<faceList *> face_list;
   vector<edgeList *> edge_list;

   // create memory for poles 0 - North Pole 1 - South Pole
   vector<poleList *> pole;
   pole.push_back(new poleList);
   pole.push_back(new poleList);

   col_geom_v geom_d;

   // build side cut half first
   bool point_cut_save = opts.point_cut;
   double inner_radius_save = opts.inner_radius;
   double outer_radius_save = opts.outer_radius;
   opts.point_cut = false;
   build_globe(geom_d, coordinates, face_list, edge_list, pole, opts);

   // start over. build base part second, then rotate it
   clear_coord(coordinates);
   clear_faces(face_list);
   clear_edges(edge_list);

   opts.inner_radius = inner_radius_save;
   opts.outer_radius = outer_radius_save;
   opts.point_cut = true;
   build_globe(geom, coordinates, face_list, edge_list, pole, opts);
   opts.point_cut = point_cut_save;

   // adjust longitude values for second half
   for (unsigned int i=0;i<face_list.size();i++)
      face_list[i]->lon += opts.longitudes.front()/2;

   // we can do the twist by transforming just one part.
   // negative angle because z reflection
   mat3d trans = mat3d::rot(0, 0, deg2rad(-hybrid_twist_angle(opts.ncon_order, opts.twist, opts.nm_shell)) ) * mat3d::refl(vec3d(0,0,1));
   geom.transform(trans);

   // finally merge the two halves and merge the vertices
   geom.append(geom_d);
   sort_merge_elems(geom, "v", pow(10, -abs(8)));

   // clean up
   clear_coord(coordinates);
   clear_faces(face_list);
   clear_edges(edge_list);
}

void process_normal(col_geom_v &geom, ncon_opts &opts)
{
   vector<coordList *> coordinates;
   vector<faceList *> face_list;
   vector<edgeList *> edge_list;

   // create memory for poles 0 - North Pole 1 - South Pole
   vector<poleList *> pole;
   pole.push_back(new poleList);
   pole.push_back(new poleList);

   build_globe(geom, coordinates, face_list, edge_list, pole, opts);

   // now we do the twisting
   vector<polarOrb *> polar_orbit;
   if (strchr(opts.closure.c_str(), 'v') || (opts.twist != 0))
      close_latitudinal_or_find_twist_plane(geom, polar_orbit, face_list, pole, opts.ncon_order, opts.point_cut,
                                            opts.longitudes, opts.add_poles, opts.closure, opts.half_model_marker);
 
   ncon_twist(geom, polar_orbit, coordinates, face_list, edge_list, opts.twist, opts.ncon_order, opts.longitudes);
   polar_orbit.clear();
   
   // clean up
   clear_coord(coordinates);
   clear_faces(face_list);
   clear_edges(edge_list);
}

void filter(col_geom_v &geom, const char *elems)
{
   if(strchr(elems, 'v'))
      geom.clear_v_cols();
   if(strchr(elems, 'e'))
      geom.clear_edges();
   if(strchr(elems, 'f'))
      geom.clear_faces();
}

void ncon_subsystem(ncon_opts opts)
{
   col_geom_v geom;
   char errmsg[MSG_SZ];

   bool point_cut_save = opts.point_cut;
   if (opts.nm_shell) {
      opts.ncon_order *= 2;
      opts.twist *= 2;
   }

   if (opts.hybrid)
      process_hybrid(geom, opts);
   else
      process_normal(geom, opts);

   if (opts.nm_shell) {
      opts.ncon_order /= 2;
      opts.twist /= 2;
      opts.point_cut = point_cut_save;
   }

   if (opts.info) {
      vector<surfaceTable *> surface_table;
      if (opts.d > 1) {
         fprintf(stderr,"surface info not available for star n_icons\n");
         if (opts.nm_shell)
            fprintf(stderr,"shell model: outer radius = %.17lf  inner radius = %.17lf\n",opts.outer_radius,opts.inner_radius);
      }
      else {
         surfaceData sd;
         ncon_info(opts.ncon_order, opts.point_cut, opts.twist, opts.hybrid, opts.info, surface_table, sd);
         surface_table.clear();
      }
      model_info(geom, opts.info);
   }
   
   // Color post-processing
   // process edges with no color
   if (opts.unused_edge_color.is_set())
      color_unused_edges(geom, opts.unused_edge_color);
   // vertices can get out of sync with their edges.
   if (!opts.edge_set_no_color)
      reassert_colored_edges(geom, col_val(DEFAULT_COLOR));
   // some faces can remain uncolored
   if ( opts.face_coloring_method )
      color_uncolored_faces(geom, col_val(DEFAULT_COLOR));

   geom.orient();
   
   filter(geom,opts.hide_elems.c_str());

   if(!geom.write(opts.ofile, errmsg))
      opts.error(errmsg);
}

void surface_subsystem(ncon_opts &opts)
{
   vector<surfaceTable *> surface_table;
   surfaceData sd;

   char form = opts.ncon_surf[0];

   fprintf(stderr,"\n");

   if (!opts.filter_case2)
      fprintf(stderr,"Note: case 2 n-icons are depicted with curly brackets\n\n");

   if ((form != 'o') && (form != 'i'))
      fprintf(stderr,"Note: non-chiral n-icons are depicted with square brackets\n\n");

   if (form == 's')
      fprintf(stderr,"Note: all even order side cut n-icons have at least two surfaces\n\n");
   else
   if (form == 'i')
      fprintf(stderr,"Note: only hybrids such that N/2 is even are shown\n\n");
   else
   if (form == 'j')
      fprintf(stderr,"Note: only hybrids such that N/2 is odd are shown\n\n");
   else
   if (form == 'k')
      fprintf(stderr,"Note: only hybrids such that N/4 is even are shown\n\n");
   else
   if (form == 'l')
      fprintf(stderr,"Note: only hybrids such that N/4 is odd are shown\n\n");

   int inc = 2;
   if (form == 'o') {
      if (is_even(opts.ncon_range.front()))
         opts.ncon_range.front()++;
   }
   else {
      if (!is_even(opts.ncon_range.front()))
         opts.ncon_range.front()++;
   }

   if ((form == 'i') || (form == 'j')) {
      if (((form == 'i') && !is_even(opts.ncon_range.front()/2)) ||
          ((form == 'j') &&  is_even(opts.ncon_range.front()/2)))
            opts.ncon_range.front()+=2;

      inc += 2;
      form = 'h';
   }
   else
   if ((form == 'k') || (form == 'l')) {
      if (((form == 'k') && !is_even(opts.ncon_range.front()/4)) ||
          ((form == 'l') &&  is_even(opts.ncon_range.front()/4)))
            opts.ncon_range.front()+=4;

      inc += 6;
      form = 'h';
   }

   if (opts.long_form) {
      fprintf(stderr,"                       Surfaces                       Edges\n");
      fprintf(stderr,"                       ------------------------------ ------------------------\n");
      fprintf(stderr,"%s  %s          %s %s %s %s %s\n\n",
         "Order","n-icon","Total","Continuous","Discontinuous","Continuous","Discontinuous");
   }

   int last = 0;
   for (int ncon_order=opts.ncon_range.front(); ncon_order<=opts.ncon_range.back(); ncon_order+=inc) {
      if (is_even(ncon_order)) {
         if (form == 'h')
            last = (int)floor((double)(ncon_order+2)/4);
         else
            last = (int)floor((double)ncon_order/4);
      }
      else
         last = (int)floor((double)ncon_order/2);

      bool point_cut = false;
      bool hybrid = false;
      bool info = false;

      if (form == 'n' || form == 'o' )
         point_cut = true;
      else
      if (form == 's')
         point_cut = false;
      else
      if (form == 'h')
         hybrid = true;

      bool none = true;

      if (opts.long_form)
         fprintf(stderr,"%-5d: ",ncon_order);
      else
         fprintf(stderr,"%d: ",ncon_order);

      for (int twist=2;twist<=last;twist++) {
         ncon_info(ncon_order, point_cut, twist, hybrid, info, surface_table, sd);

         if ((!is_even(ncon_order) && sd.total_surfaces > 1) ||
            (form != 's' && sd.total_surfaces > 1) || (form == 's' && sd.total_surfaces > 2)) {
            if (!sd.ncon_case2 || (sd.ncon_case2 && !opts.filter_case2)) {
               if (!none) {
                  if (opts.long_form)
                     fprintf(stderr,"%-5d: ",ncon_order);
                  else
                     fprintf(stderr,", ");
               }
               char buffer[80];
               if (sd.nonchiral)
                  sprintf(buffer,"[%d+%d]",ncon_order,twist);
               else
               if (sd.ncon_case2)
                  sprintf(buffer,"{%d+%d}",ncon_order,twist);
               else
                  sprintf(buffer,"(%d+%d)",ncon_order,twist);
               if (opts.long_form)
                     fprintf(stderr,"%-15s %5d %10d %13d %10d %13d\n",
                        buffer, sd.total_surfaces, sd.c_surfaces, sd.d_surfaces, sd.c_edges, sd.d_edges);
               else
                  fprintf(stderr,"%s",buffer);
               none = false;
            }
         }
      }
      surface_table.clear();

      if (none) {
         fprintf(stderr,"none");
         if (opts.long_form)
            fprintf(stderr,"\n");
      }
      fprintf(stderr,"\n");
   }
}

int main(int argc, char *argv[])
{
   ncon_opts opts;
   opts.process_command_line(argc, argv);

   if (!opts.ncon_surf.length())
      ncon_subsystem(opts);
   else
      surface_subsystem(opts);

   return 0;
}
