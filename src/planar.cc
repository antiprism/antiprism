/*
   Copyright (c) 2010-2011, Roger Kaufman

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
   Name: planar.cc
   Description: Various utilities for polyhedra with planar faces
   Project: Antiprism - http://www.antiprism.com
*/


#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include <ctype.h>

#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include "../base/antiprism.h"

using std::string;
using std::vector;
using std::pair;
using std::make_pair;
using std::map;
using std::min;
using std::max;


class planar_opts: public prog_opts {
   public:
      string ifile;
      string ofile;
      
      char face_color_method;
      int planar_merge_type;
      int polygon_fill_type;
      int winding_rule;
      int winding_rule_mode;
      char color_by_winding_number;
      bool winding_div2; // not implemented
      bool find_direction;
      bool verbose;
      int orient;
      bool hole_detection;
      vec3d center;
      bool stitch_faces;
      bool simplify_face_edges;
      bool delete_invisible_faces;
      char edge_blending;
      char special_edge_processing;
      col_val zero_density_color;
      bool zero_density_force_blend;
      double brightness_adj;
      int color_system_mode;
      bool cmy_mode;
      bool ryb_mode;
      double sat_power;
      double value_power;
      double sat_threshold;
      double value_advance;
      int alpha_mode;
      int face_opacity;

      color_map_multi map;
      color_map_multi map_negative;

      double epsilon;

      planar_opts(): prog_opts("planar"),
                        face_color_method('\0'),
                        planar_merge_type(0),
                        polygon_fill_type(0),
                        winding_rule(INT_MAX),
                        winding_rule_mode(INT_MAX),
                        color_by_winding_number('\0'),
                        winding_div2(false), // not implemented
                        find_direction(false),
                        verbose(false),
                        orient(0),
                        hole_detection(true),
                        stitch_faces(false),
                        simplify_face_edges(false),
                        delete_invisible_faces(false),
                        edge_blending('\0'),
                        special_edge_processing('\0'),
                        zero_density_color(col_val::invisible),
                        zero_density_force_blend(false),
                        brightness_adj(-2.0),
                        color_system_mode(3),
                        cmy_mode(false),
                        ryb_mode(false),
                        sat_power(0.0),
                        value_power(0.0),
                        sat_threshold(1.0),
                        value_advance(0.0),
                        alpha_mode(3),
                        face_opacity(255),
                        epsilon(0)
                        {}

      void process_command_line(int argc, char **argv);
      void usage();
};

void planar_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Read a file in OFF format and convert overlapping coplanar polygons into\n"
"non-overlapping polygon tiles. If input_file is not given the program \n"
"reads from standard input.\n"
"\n"
"Options\n"
"%s"
"  -d <opt>  blend overlapping (tile) or adjacent (merge) planar faces\n"  
"               tile=1, merge=2 (default: none)\n"
"  -p <opt>  polygon fill algorithm.  angular=1, modulo2=2 (pnpoly)\n"
"               triangulation=3, even_overlap=4, alt_modulo2=5 (default: 1)\n"
"  -w <opt>  winding rule, include face parts according to winding number\n"
"               odd, even, positive, negative, nonzero, zero (default: none)\n"
"               zodd, zeven, zpositive, znegative (includes zero)\n"
"               or symbol proceeding integer: eq, ne, gt, ge, lt, le  meaning\n"
"               equal, not equal, greater than, greater than or equal, less\n"
"               than, less than or equal. use 'a' for absolute value  e.g. gea2\n"
"  -z        use direction of normals for hemispherical winding numbers\n"
"  -V        verbose output (of minimum and maximum winding numbers)\n"
"  -H        turn off hole detection\n"
"  -S        stitch seams created by tiling or merging\n"
"  -I        rid faces of extra in-line vertices\n"
"  -e <opt>  blend existing explicit edges and/or vertices using blend options\n"
"               e - edges, v - vertices, b - both (-d forces b, otherwise none)\n"
"  -E <opt>  remove explicit edges, blend new ones using face colors (sets -S)\n"
"               e - edges only, v - also blend vertices (default: none)\n"
"               V - also blend invisible vertices, s - strip edges and vertices\n"
"  -D        delete invisible faces created by winding rule\n"
"  -O <opt>  orient the faces first (if possible) then for volume\n"
"               positive=1, negative=2, reverse=3, or use flip=4\n"
"               which reverses the orientation of the model as it was input\n"
"  -C <xyz>  center of model, in form 'X,Y,Z' (default: centroid)\n"
"  -l <lim>  minimum distance for unique vertex locations as negative exponent\n"
"               (default: %d giving %.0e)\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\nColor Blending Options (for option -d)\n"
"  -M <mode> color blending mode. HSV=1, HSL=2, RGB=3 (default: 3)\n"
"  -s <sat>  HSV/HSL saturation curve. Greather than 0 (default: 1)\n"
"               1.0 - no curve. lower than 1.0 makes blends more pastel\n"
"  -t <val>  HSV/HSL threshold to use average saturation (default: 1)\n"
"               between 0.0 (all averaging) and 1.0 (no averaging)\n"
"  -v <val>  HSV/HSL value curve (default: 0)\n"
"               simulates subtractive coloring for blending 3 or more colors\n"
"               RGB: Red+Green+Blue = White   Cyan+Magenta+Yellow = Black\n"
"               RYB: Red+Yellow+Blue = Black  Green+Magenta+Orange = White\n"
"               1.0 - no curve. lower than 1.0 number makes blends lighter\n"
"               0.0 - use average value instead\n"
"  -u <val>  HSV/HSL value advance. Rotates meaning of white and black\n"
"               valid values 0.0 to 120.0 degrees (default: 0)\n"
"  -a <int>  alpha for blend. average=1, minimum=2, maximum=3 (default: 3)\n"
"  -y        RYB mode. Blend colors as in Red-Yellow-Blue color wheel\n"
"  -c        CMY mode. Complementary colors.  RGB->(RYB/GMO)->CMY->blend\n"
"  -b <val>  brightness adjustment for areas of blended colors\n"
"               valid values -1.0 to +1.0 (default: no adjustment)\n"
"               negative for darker, positive for lighter\n"
"               at 0, an area with 2 blended colors will not change\n"
"               but areas with 3 or more will become slightly darker\n"
"\nColoring Options (run 'off_util -H color' for help on color formats)\n"
"  -f <opt>  take face colors from map (processed before -d)\n"
"               n - unique color for faces with same normals\n"
"               p - unique color for faces on same planes only\n"
"               o - unique color for faces on same and opposite normals\n"
"  -T <tran> face transparency. valid range from 0 (invisible) to 255 (opaque)\n"
"  -Z <col>  color for areas found colorless by winding (default: invisible)\n"
"               key word: b - force a color blend\n"
"  -W <opt>  color by winding number, using maps (overrides option -f)\n"
"               w - use actual winding number\n"
"               a - absolute value of winding number\n"
"               n - negative of absolute value of winding number\n"
"  -m <maps> color maps for faces to be tried in turn (default: compound)\n"
"  -n <maps> maps for negative winding numbers (default: rng17_S0V0.5:0)\n"
"               (map position zero not used. default is 16 gradients)\n"
"\n",prog_name(), help_ver_text, int(-log(::epsilon)/log(10) + 0.5), ::epsilon);
}

void planar_opts::process_command_line(int argc, char **argv)
{
   opterr = 0;
   char c;
   char errmsg[MSG_SZ];

   int sig_compare = INT_MAX;
   string id;
   string map_file;
   string map_file_negative;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hd:p:w:zVO:HC:SIe:E:Db:M:s:t:v:u:a:cyf:T:m:Z:W:n:l:o:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         case 'd':
            id = get_arg_id(optarg, "tile=1|merge=2", argmatch_add_id_maps, errmsg);
            if(id=="")
               error(errmsg);
            planar_merge_type = atoi(id.c_str());
            break;

         case 'p':
            id = get_arg_id(optarg, "angular=1|modulo2=2|triangulation=3|even_overlap=4|alt_modulo2=5", argmatch_add_id_maps, errmsg);
            if(id=="")
               error(errmsg);
            polygon_fill_type = atoi(id.c_str());
            break;

         case 'w': {
            id = optarg;
            // symbol proceeding integer. Symbols: eq, ne, gt, ge, lt, le"
            // correspond to 1,2,3,4,5,6
            if(id.substr(0,2)=="eq")
               winding_rule_mode = 1;
            else
            if((id.substr(0,2)=="ne") && (id.substr(0,3)!="neg"))
               winding_rule_mode = 2;
            else
            if(id.substr(0,2)=="gt")
               winding_rule_mode = 3;
            else
            if(id.substr(0,2)=="ge")
               winding_rule_mode = 4;
            else
            if(id.substr(0,2)=="lt")
               winding_rule_mode = 5;
            else
            if(id.substr(0,2)=="le")
               winding_rule_mode = 6;

            if (winding_rule_mode != INT_MAX) {
               int substr_start = 2;
               // negative winding_rule_mode will trigger absolute value evaluation later
               if((id.length() > 2) && (id.substr(2,1)=="a")) {
                  winding_rule_mode = -winding_rule_mode;
                  substr_start = 3;
               }

               int tmp;
               char buff;
               if(sscanf((id.substr(substr_start)).c_str(), " %d %c", &tmp, &buff) != 1)
                  error("expecting integer number after eq,ne,gt,ge,lt,le for winding number", c);
               winding_rule = tmp;
            }
            else {
               winding_rule_mode = 0;
               // built in winding rules become 1 through 10 (when winding_rule_mode is zero)
               id = get_arg_id(optarg, "odd|even|positive|negative|nonzero|zodd|zeven|zpositive|znegative|zero", argmatch_default, errmsg);
               if(id=="")
                  error(errmsg);
               winding_rule = (atoi(id.c_str()) + 1);
            }
            break;
         }

         case 'z':
            find_direction = true;
            break;

         case 'V':
            verbose = true;
            break;

         case 'O':
            id = get_arg_id(optarg, "positive=1|negative=2|reverse=3|flip=4", argmatch_add_id_maps, errmsg);
            if(id=="")
               error(errmsg);
            orient = atoi(id.c_str());
            break;

         case 'H':
            hole_detection = false;
            break;

         case 'C':
            if(!center.read(optarg, errmsg))
               error(errmsg, c);
            break;

         case 'S':
            stitch_faces = true;
            break;

         case 'I':
            simplify_face_edges = true;
            break;

         case 'e':
            if(strspn(optarg, "evb") != strlen(optarg) || strlen(optarg)>1)
               error(msg_str("pre edge blending is '%s', must be e, v, or b", optarg), c);
            edge_blending = *optarg;
            break;

         case 'E':
            if(strspn(optarg, "evVs") != strlen(optarg) || strlen(optarg)>1)
               error(msg_str("post edge blending is '%s', must be e, v, V, or s", optarg), c);
            special_edge_processing = *optarg;

            stitch_faces = true;
            break;

         case 'D':
            delete_invisible_faces = true;
            break;

         case 'b':
            if(!read_double(optarg, &brightness_adj, errmsg))
               error(errmsg, c);
            if(brightness_adj < -1.0 || brightness_adj > 1.0)
               error("brightness adjustment must be between -1.0 and 1.0", c);
            break;

         case 'M':
            id = get_arg_id(optarg, "hsv=1|hsl=2|rgb=3", argmatch_add_id_maps, errmsg);
            if(id=="")
               error(errmsg);
            color_system_mode = atoi(id.c_str());
            break;

         case 's':
            if(!read_double(optarg, &sat_power, errmsg))
               error(errmsg, c);
            if(sat_power <= 0.0)
               error("color centroid saturation curve must be greater than zero", c);
            break;
            
         case 't':
            if(!read_double(optarg, &sat_threshold, errmsg))
               error(errmsg, c);
            if(sat_threshold < 0.0 || sat_threshold > 1.0)
               error("HSV/HSL threshold must be between 0 and 1", c);
            break;
            
         case 'v':
            if(!read_double(optarg, &value_power, errmsg))
               error(errmsg, c);
            if(value_power < 0.0)
               error("value curve must be greater than or equal to zero", c);
            break;

        case 'u':
            if(!read_double(optarg, &value_advance, errmsg))
               error(errmsg, c);
            if(value_advance < 0.0 || value_advance > 120.0)
               error("HSV/HSL value advance must be between 0 and 120", c);
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

         case 'f':
            if(strlen(optarg)!=1 || !strchr("npo", *optarg))
               error("color method must be n, p or o");
            face_color_method = *optarg;
            break;
         
         case 'T':
            if(!read_int(optarg, &face_opacity, errmsg))
               error(errmsg, c);
            if(face_opacity < 0 || face_opacity > 255)
               error("face transparency must be between 0 and 255", c);
            break;

         case 'm':
            map_file = optarg;
            break;

         case 'Z':
            if(strlen(optarg)==1 && strchr("b", *optarg))
               zero_density_force_blend = true;
            else
            if(!zero_density_color.read(optarg, errmsg))
               error(errmsg, c);
            break;

         case 'W':
            if(strspn(optarg, "wan") != strlen(optarg) || strlen(optarg)>1)
               error(msg_str("color by winding number is '%s', must be w, a, or n", optarg), c);
            color_by_winding_number = *optarg;
            break;

         case 'n':
            map_file_negative = optarg;
            break;

         case 'l':
            if(!read_int(optarg, &sig_compare, errmsg))
               error(errmsg, c);
            if(sig_compare < 0) {
               warning("limit is negative, and so ignored", c);
            }
            if(sig_compare > DEF_SIG_DGTS) {
               warning("limit is very small, may not be attainable", c);
            }
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

   //if (edge_blending && special_edge_processing)
   //   error("edge blending and special edge processing cannot be used together","e");
      
   if (!planar_merge_type) {
      if (polygon_fill_type)
         warning("polygon fill has no effect if tile or merge is not selected","p");
      if (winding_rule != INT_MAX)
         warning("winding rule has no effect if tile or merge is not selected","w");
      // zero density area cannot be colored with not tile or merge
      if (color_by_winding_number && (!zero_density_color.is_inv() || zero_density_force_blend))
         warning("zero density areas cannot be colored or blended if tile or merge is not selected","Z");
   }
   else {
      // set default polygon fill type here
      if (!polygon_fill_type)
         polygon_fill_type = 1;
      // this is only true when tile/merge
      if (color_by_winding_number && polygon_fill_type != 1 && polygon_fill_type != 3)
         error("when tile or merge, color by winding number polygon fill type must be 1 or 3","p");
      // when tiling or merging, edges must be blended
      if (edge_blending && (edge_blending != 'b'))
         warning("when tile or merge, both edges and vertices are always blended","e");
      edge_blending = 'b';
   }

   if (!map_file.size())
      map_file = "compound";
   if(!map.init(map_file.c_str(), errmsg))
      error(errmsg, 'm');

   if (!map_file_negative.size())
      map_file_negative = "rng17_S0V0.5:0";
   if(!map_negative.init(map_file_negative.c_str(), errmsg))
      error(errmsg, 'n');

   epsilon = (sig_compare != INT_MAX) ? pow(10, -sig_compare) : ::epsilon;
}


bool cmp_verts(const pair<vec3d, int> &a, const pair<vec3d, int> &b, const double &eps)
{
   return (compare(a.first,b.first,eps)<0);
}

class vert_cmp
{
public:
   double eps;
   vert_cmp(double ep): eps(ep) {}
   bool operator() (const pair<vec3d, int> &a, const pair<vec3d, int> &b) { return cmp_verts(a, b, eps); }
};

// second check: if polygons are not actually on the same plane, split them up
// RK - this can only happen when the normals are not forced outward
vector<vector<int> > on_same_plane_filter(const geom_if &geom, const vec3d &normal, const vector<int> &coplanar_faces, const bool &filtered, const double &eps)
{
   const vector<vector<int> > &faces = geom.faces();
   const vector<vec3d> &verts = geom.verts();

   vector<vector<int> > coplanar_faces_filtered;

   if (!filtered) {
      coplanar_faces_filtered.push_back(coplanar_faces);
      return coplanar_faces_filtered;
   }
      
   vector<int> coplanar_faces_actual;
   vector<int> written;

   int sz = coplanar_faces.size();
   for(int i=0; i<sz; i++) {
      int face_idx1 = coplanar_faces[i];
      if(find(written.begin(), written.end(), face_idx1) != written.end())
         continue;
      coplanar_faces_actual.push_back(face_idx1);
      for(int j=i+1; j<sz; j++) {
         if (i==j)
            continue;
         int face_idx2 = coplanar_faces[j];
         vec3d v0 = verts[faces[face_idx1][0]];
         vec3d P = verts[faces[face_idx2][0]];
         if (double_eq(vdot(v0-P, normal), 0.0, eps)) {
            coplanar_faces_actual.push_back(face_idx2);
            written.push_back(face_idx2);
         }
      }

      coplanar_faces_filtered.push_back(coplanar_faces_actual);
      coplanar_faces_actual.clear();
   }

   return coplanar_faces_filtered;
}

void build_coplanar_faces_list(const geom_if &geom, vector<vector<int> > &coplanar_faces_list, vector<xnormal> &coplanar_normals, const fnormals &fnormals,
                               const bool &point_outward, const bool &fold_normals, const bool &fold_normals_hemispherical, const bool &filtered, const double &eps)
{
   // seperate out hemispherical because they often get mixed in
   vector<pair<vec3d, int> > hemispherical_table;
   vector<pair<vec3d, int> > face_normal_table;

   for(unsigned int i=0; i<geom.faces().size(); i++) {
      bool hemi = fnormals[i].is_hemispherical();

      vec3d normal = fnormals[i].raw();
      // point outward option here
      if (point_outward && !hemi)
         normal = fnormals[i].outward();
      normal = normal.unit();

      pair<vec3d, int> face_normal_pair;
      face_normal_pair.first = normal;
      face_normal_pair.second = i;

      if (hemi)
         hemispherical_table.push_back(face_normal_pair);
      else
         face_normal_table.push_back(face_normal_pair);
   }

   vector<int> coplanar_faces;

   // hemispherical normals are folded only with specific option
   // if folded, this is what associates them on the same plane
   int sz = hemispherical_table.size();
   if (fold_normals_hemispherical && (sz > 1)) {
      for(int i=0; i<sz-1; i++) {
         for(int j=i+1; j<sz; j++) {
            if (!compare(hemispherical_table[i].first,-hemispherical_table[j].first,eps)) {
               hemispherical_table[j].first = -hemispherical_table[j].first;
            }
         }
      }
   }

   // collect hemispherical which are coplanar
   if (sz) {
      stable_sort( hemispherical_table.begin(), hemispherical_table.end(), vert_cmp(eps) );

      // at least one face is on a plane
      coplanar_faces.push_back(hemispherical_table[0].second);
      // find faces on same plane and group them together
      for(int i=1; i<sz; i++) {
         if (compare(hemispherical_table[i].first,hemispherical_table[i-1].first,eps)) {
            vector<vector<int> > coplanar_faces_filtered = on_same_plane_filter(geom, hemispherical_table[i-1].first, coplanar_faces, filtered, eps);
            for(unsigned int j=0; j<coplanar_faces_filtered.size(); j++) {
               coplanar_faces_list.push_back(coplanar_faces_filtered[j]);
               coplanar_normals.push_back(fnormals[hemispherical_table[i-1].second]);
            }
            coplanar_faces.clear();
         }
         coplanar_faces.push_back(hemispherical_table[i].second);
      }
      vector<vector<int> > coplanar_faces_filtered = on_same_plane_filter(geom, hemispherical_table[sz-1].first, coplanar_faces, filtered, eps);
      for(unsigned int j=0; j<coplanar_faces_filtered.size(); j++) {
         coplanar_faces_list.push_back(coplanar_faces_filtered[j]);
         coplanar_normals.push_back(fnormals[hemispherical_table[sz-1].second]);
      }
      coplanar_faces.clear();
   }

   
   // non-hemispherical normals are folded only with specific option
   sz = face_normal_table.size();
   if (fold_normals && (sz > 1)) {
      for(int i=0; i<sz-1; i++) {
         for(int j=i+1; j<sz; j++) {
            if (!compare(face_normal_table[i].first,-face_normal_table[j].first,eps)) {
               face_normal_table[j].first = -face_normal_table[j].first;
            }
         }
      }
   }

   // collect non-hemispherical which are coplanar
   if (sz) { 
      stable_sort( face_normal_table.begin(), face_normal_table.end(), vert_cmp(eps) );

      // at least one face is on a plane
      coplanar_faces.push_back(face_normal_table[0].second);
      // find faces on same plane and group them together
      for(int i=1; i<sz; i++) {
         if (compare(face_normal_table[i].first,face_normal_table[i-1].first,eps)) {
            vector<vector<int> > coplanar_faces_filtered = on_same_plane_filter(geom, face_normal_table[i-1].first, coplanar_faces, filtered, eps);
            for(unsigned int j=0; j<coplanar_faces_filtered.size(); j++) {
               coplanar_faces_list.push_back(coplanar_faces_filtered[j]);
               coplanar_normals.push_back(fnormals[face_normal_table[i-1].second]);
            }
            coplanar_faces.clear();
         }
         coplanar_faces.push_back(face_normal_table[i].second);
      }
      vector<vector<int> > coplanar_faces_filtered = on_same_plane_filter(geom, face_normal_table[sz-1].first, coplanar_faces, filtered, eps);
      for(unsigned int j=0; j<coplanar_faces_filtered.size(); j++) {
         coplanar_faces_list.push_back(coplanar_faces_filtered[j]);
         coplanar_normals.push_back(fnormals[face_normal_table[sz-1].second]);
      }
      coplanar_faces.clear();
   }

   hemispherical_table.clear();
   face_normal_table.clear();
}

// idx will be the dimension not included so that 3D -> 2D projection occurs
void project_using_normal(const vec3d &normal, int &idx, int &sign)
{
   vector<pair<double, int> > normal_info(3);

   for(unsigned int i=0; i<3; i++) {
      normal_info[i].second = i;
      normal_info[i].first = fabs(normal[i]);   // absolute value for sorting
   }

   stable_sort( normal_info.begin(), normal_info.end() );
   
   idx = normal_info[normal_info.size()-1].second;
   sign = (normal[idx] >= 0.0) ? 1 : -1;
}

// if intersection point does not exist, insert a new one. Otherwise only return the index of the existing one
int vertex_into_geom(col_geom_v &geom, const vec3d &P, col_val vcol, const double &eps)
{
   int v_idx = find_vertex_by_coordinate(geom, P, eps);
   if (v_idx == -1) {
      geom.add_col_vert(P, vcol);
      v_idx = geom.verts().size()-1;
   }

   return v_idx;
}

// if edge already exists, do not create another one and return false. return true if new edge created
// check if edge1 and edge2 indexes are equal. If so do not allow an edge length of 0 to occur
bool edge_into_geom(col_geom_v &geom, const int &v_idx1, const int &v_idx2, col_val ecol)
{
   // no edge length 0 allowed
   if (v_idx1 == v_idx2)
      return false;

   vector<int> new_edge = make_edge(v_idx1, v_idx2);

   int answer = find_edge_in_edge_list(geom.edges(), new_edge);
   if (answer < 0)
      geom.add_col_edge(new_edge, ecol);

   return ((answer < 0) ? true : false);
}

// save free edges into a geom
col_geom_v free_edges_into_geom(const col_geom_v &geom)
{
   const vector<vector<int> > &faces = geom.faces();
   const vector<vector<int> > &edges = geom.edges();

   col_geom_v egeom;
   egeom.add_verts(geom.verts());
   for (unsigned int i=0;i<edges.size();i++) {
      vector<int> face_idx = find_faces_with_edge(faces, edges[i]);
      if (!face_idx.size())
         egeom.add_col_edge(edges[i],geom.get_e_col(i));
   }
   egeom.delete_verts(egeom.get_info().get_free_verts());

   return egeom;
}

bool is_point_on_polygon_edges(const geom_if &polygon, const vec3d &P, const double &eps)
{
   const vector<int> &face = polygon.faces()[0];
   const vector<vec3d> &verts = polygon.verts();

   bool answer = false;

   int fsz = face.size();
   for(int i=0; i<fsz; i++) {
      vec3d v1 = verts[face[i]];
      vec3d v2 = verts[face[(i+1)%fsz]];
      if ((point_in_segment(P, v1, v2, eps)).is_set()) {
         answer = true;
         break;
      }
   }

   return answer;
}

bool is_point_in_bbox_2D(const geom_if &geom, const vec3d &P, const int &idx, const double &eps)
{
   bound_box bb(geom.verts());
   vec3d min = bb.get_min();
   vec3d max = bb.get_max();

   int idx1 = (idx+1)%3;
   int idx2 = (idx+2)%3;

   //return !(P[idx1] < min[idx1] || P[idx1] > max[idx1] || P[idx2] < min[idx2] || P[idx2] > max[idx2]);
   return !(double_lt(P[idx1],min[idx1],eps) || double_ge(P[idx1],max[idx1],eps) || double_le(P[idx2],min[idx2],eps) || double_gt(P[idx2],max[idx2],eps));
}

// furshished by Adrian Rossiter
int side_of_line(const vec3d &P, const vec3d &A, const vec3d &B, const int &idx1, const int &idx2, const double &eps)
{
   double a = B[idx1] - A[idx1];
   double b = B[idx2] - A[idx2];
   double c = P[idx1] - A[idx1];
   double d = P[idx2] - A[idx2];
   double det = a*d - b*c;

/*
   if(det<-eps)
      return -1;
   else
   if(det>eps)
      return 1;
   else
      return 0;
*/
   return double_compare(det,0,eps);
}

// edge is between -eps and eps. '>-eps' will include edge. >-eps same as >=0
bool is_point_inside_triangle_or_edge(const vec3d &P, const vec3d &A, const vec3d &B, const vec3d &C,
                                      const int &idx1, const int &idx2, const int &sign, const double &eps)                          
{
   return (side_of_line(P,A,B,idx1,idx2,eps)*sign>-eps && side_of_line(P,B,C,idx1,idx2,eps)*sign>-eps && side_of_line(P,C,A,idx1,idx2,eps)*sign>-eps);
}

// polygon sent to this function MUST be triangulated
// includes edges and vertices as inside (bounds are considered inside points)
bool InsidePolygon_triangulated(const geom_if &polygon, const vec3d &P, const int &idx, const int &sign, const bool &include_edges, const double &eps)
{
   const vector<vector<int> > &faces = polygon.faces();
   const vector<vec3d> &verts = polygon.verts();

   int idx1 = (idx+1)%3;
   int idx2 = (idx+2)%3;

   bool answer = false;
   for(unsigned int i=0; i<faces.size(); i++) {
      // speed up. don't measure if point is clearly outside of triangles min max area
      vector<int> vec(1);
      vec[0] = i;
      col_geom_v triangle = faces_to_geom(polygon, vec);
      if (!is_point_in_bbox_2D(triangle, P, idx, eps))
         continue;

      bool found = false;
      if (is_point_inside_triangle_or_edge(P,verts[faces[i][0]],verts[faces[i][1]],verts[faces[i][2]],idx1,idx2,sign,eps)) {
         answer = true;
         found = true;
      }

      // if point IS inside polygon and bounds are not considered inside, reject if on polygon edges
      if (answer && !include_edges) {
         if (is_point_on_polygon_edges(triangle, P, eps))
            answer = false;
      }

      if (found)
         break;
      }

   return answer;
}

// solution 1 function ported from http://paulbourke.net/geometry/insidepoly/
// code is by Alexander Motrichuk, it returns true for interior points and false for exterior points
// RK - modification: bound is true if on bound (edge) or vertex. return true for inside polygon
// RK - does modulo-2 implicitly

//SOLUTION #1 (2D) - Redesigned
//horizintal left cross over direction algorithm
//-----------------------------------------------
//  bound   |   value that will be returned only if (p) lies on the bound or vertex
bool InsidePolygon_solution1(const geom_if &polygon, const vec3d &P, bool &bound, const int &idx, const double &eps)
{
   const vector<vec3d> &verts = polygon.verts();

   int idx1 = (idx+1)%3;
   int idx2 = (idx+2)%3;

   double p_x = P[idx1];
   double p_y = P[idx2];

   const vector<int> &face = polygon.faces()[0];
   int N = face.size();

   //left vertex
   double p1_x = verts[face[0]][idx1];
   double p1_y = verts[face[0]][idx2];

   //cross points count of x
   int __count = 0;

   bound = false;

   //check all rays
   for(int i = 1; i <= N; ++i) {

      //point is a vertex
      if (double_eq(p_x, p1_x, eps) && double_eq(p_y, p1_y, eps)) {
         bound = true;
         return true;
      }

      //right vertex
      double p2_x = verts[face[i%N]][idx1];
      double p2_y = verts[face[i%N]][idx2];

      //ray is outside of our interests
      //if (p_y < min(p1_y, p2_y) || p_y > max(p1_y, p2_y)) {
      if (double_lt(p_y,min(p1_y, p2_y),eps) || double_gt(p_y,max(p1_y, p2_y),eps)) {
         //next ray left point
         p1_x = p2_x;
         p1_y = p2_y;
         continue;
      }

      //ray is crossing over by the algorithm (common part of)
      //if (p_y > min(p1_y, p2_y) && p_y < max(p1_y, p2_y)) {
      if (double_gt(p_y,min(p1_y, p2_y),eps) && double_lt(p_y,max(p1_y, p2_y),eps)) {

         //x is before of ray
         //if (p_x <= max(p1_x, p2_x)) {
         if (double_le(p_x, max(p1_x, p2_x), eps)) {

            //overlies on a horizontal ray
            //if (double_eq(p1_y, p2_y, eps) && p_x >= min(p1_x, p2_x)) {
            if (double_eq(p1_y, p2_y, eps) && (double_ge(p_x, min(p1_x, p2_x), eps))) {
               bound = true;
               return true;
            }

            //ray is vertical
            if (double_eq(p1_x, p2_x, eps)) {

               //overlies on a ray
               if (double_eq(p1_x, p_x, eps)) {
                  bound = true;
                  return true;
               }
               //before ray
               else
                   ++__count;
            }

            //cross point on the left side
            else {

               //cross point of x
               double xinters = (p_y - p1_y) * (p2_x - p1_x) / (p2_y - p1_y) + p1_x;

               //overlies on a ray
               //if (fabs(p_x - xinters) < eps) {
               if (double_lt(fabs(p_x - xinters), 0, eps)) {
                 bound = true;
                 return true;
               }

               //before ray
               //if (p_x < xinters)
               if (double_le(p_x,xinters,eps))
                  ++__count;
            }
         }
      }

      //special case when ray is crossing through the vertex
      else
      {
         //p crossing over p2
         //if (p_y == p2_y && p_x <= p2_x) {
         if (double_eq(p_y, p2_y, eps) && double_le(p_x, p2_x, eps)) {

            //next vertex
            double p3_y = verts[face[(i+1)%N]][idx2];

            //p_y lies between p1_y & p3_y
            //if (p_y >= min(p1_y, p3_y) && p_y <= max(p1_y, p3_y)) {
            if (double_ge(p_y, min(p1_y, p3_y), eps) && double_le(p_y, max(p1_y, p3_y), eps)) {
               ++__count;
            }
            else {
               __count += 2;
            }
         }
      }

      //next ray left point
      p1_x = p2_x;
      p1_y = p2_y;
   }

   //EVEN
   if (__count % 2 == 0)
      return(false);
   //ODD
   else
      return(true);
}

// wrapper
bool InsidePolygon_solution1(const geom_if &polygon, const vec3d &P, const int &idx, const bool &include_edges, const double &eps)
{
   bool bound = false;
   bool answer = InsidePolygon_solution1(polygon, P, bound, idx, eps);

   // if point IS inside polygon and bounds are not considered inside, reject if on polygon edges
   if (answer && bound && !include_edges)
      answer = false;

   return answer;
}

// RK - a second solution 1 function found on http://paulbourke.net/geometry/insidepoly/
// RK - function ported from http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
// RK - does modulo-2 implicitly

// by Randolph Franklin, it returns 1 for interior points and 0 for exterior points

bool pnpoly(const geom_if &polygon, const vec3d &P, const int &idx, const double &eps)
{
   const vector<vec3d> &verts = polygon.verts();

   int idx1 = (idx+1)%3;
   int idx2 = (idx+2)%3;

   double testx = P[idx1];
   double testy = P[idx2];

   const vector<int> &face = polygon.faces()[0];
   int nvert = face.size();

   int i, j, c = 0;
   for (i = 0, j = nvert-1; i < nvert; j = i++) {
      double vertx_i = verts[face[i]][idx1];
      double verty_i = verts[face[i]][idx2];
      double vertx_j = verts[face[j]][idx1];
      double verty_j = verts[face[j]][idx2];

      //if ( ((verty_i>testy) != (verty_j>testy)) &&
	   //     (testx < (vertx_j-vertx_i) * (testy-verty_i) / (verty_j-verty_i) + vertx_i) )
      if ( (double_gt(verty_i,testy,eps) != double_gt(verty_j,testy,eps)) &&
	        double_lt(testx,(vertx_j-vertx_i) * (testy-verty_i) / (verty_j-verty_i) + vertx_i,eps) )
         c = !c;
   }

   return (c ? true : false);
}

// wrapper
bool pnpoly(const geom_if &polygon, const vec3d &P, const int &idx, const bool &include_edges, const double &eps)
{
   bool answer = pnpoly(polygon, P, idx, eps);

   // have to check bounds conditions if we might have to negate the answer
   if ((answer && !include_edges) || (!answer && include_edges)) {
      if (is_point_on_polygon_edges(polygon, P, eps))
         answer = !answer;
   }

   return answer;
}

//===================================================================

// ported from http://paulbourke.net/geometry/insidepoly/

// Return the angle between two vectors on a plane
// The angle is from vector 1 to vector 2, positive anticlockwise
// The result is between -pi -> pi

double Angle2D(const double &x1, const double &y1, const double &x2, const double &y2, const double &eps)
{
   double dtheta,theta1,theta2;

   theta1 = atan2(y1,x1);
   theta2 = atan2(y2,x2);
   dtheta = theta2 - theta1;
   //while (dtheta > M_PI)
   while (double_gt(dtheta,M_PI,eps))
      dtheta -= M_PI * 2.0;
   //while (dtheta < -M_PI)
   while (double_lt(dtheta,-M_PI,eps))
      dtheta += M_PI * 2.0;

   return(dtheta);
}

// method 2 function ported from http://paulbourke.net/geometry/insidepoly/
// RK - does not do modulo-2 implicitly

// Another solution forwarded by Philippe Reverdy is to compute the sum of the angles made between the test point and each pair of points making up the polygon.
// If this sum is 2pi then the point is an interior point, if 0 then the point is an exterior point.
// This also works for polygons with holes given the polygon is defined with a path made up of coincident edges into and out of the hole as is common practice in many CAD packages.

bool InsidePolygon_solution2(const geom_if &polygon, const vec3d &P, const int &idx, const double &eps)
{
   const vector<vec3d> &verts = polygon.verts();

   int idx1 = (idx+1)%3;
   int idx2 = (idx+2)%3;

   double p_x = P[idx1];
   double p_y = P[idx2];

   const vector<int> &face = polygon.faces()[0];
   int n = face.size();

   double angle=0;
   for (int i=0;i<n;i++) {
      double p1_x = verts[face[i]][idx1] - p_x;
      double p1_y = verts[face[i]][idx2] - p_y;
      double p2_x = verts[face[(i+1)%n]][idx1] - p_x;
      double p2_y = verts[face[(i+1)%n]][idx2] - p_y;
      angle += Angle2D(p1_x,p1_y,p2_x,p2_y,eps);
   }

   //if (fabs(angle) < M_PI)
   if (double_le(fabs(angle),M_PI,eps))
      return(false);
   else
      return(true);
}

// wrapper
bool InsidePolygon_solution2(const geom_if &polygon, const vec3d &P, const int &idx, const bool &include_edges, const double &eps)
{
   bool answer = InsidePolygon_solution2(polygon, P, idx, eps);

   // have to check bounds conditions if we have to negate the answer
   if ((answer && !include_edges) || (!answer && include_edges)) {
      if (is_point_on_polygon_edges(polygon, P, eps))
         answer = !answer;
   }

   return answer;
}

// find if a point is in a polygon
// the polygon is input as the only face in a geom
// set include_edges = true if point is allowed to land on outer edges or vertices
// the normal is expected to be a unit normal
bool is_point_inside_polygon(const col_geom_v &polygon, const vec3d &P, const vec3d &normal, const bool &include_edges, const bool &on_same_plane_check,
                             const int &polygon_fill_type, const double &eps)
{  
   // quick check to see if P is actually on the polygon
   // funished by Adrian Rossiter; is there any distance from P to the unit normal?
   if (on_same_plane_check)
      if (double_ne(vdot(polygon.verts()[0]-P, normal), 0.0, eps))
         return false;

   int idx = 0;
   int sign = 0;
   project_using_normal(normal, idx, sign);

   // this does not filter much. more overhead than it is worth
   //if (!is_point_in_bbox_2D(polygon, P, idx))
   //   return false;

   bool answer = false;

   if (polygon_fill_type == 1)
      answer = InsidePolygon_solution2(polygon, P, idx, include_edges, eps);
   else
   if (polygon_fill_type == 2 || polygon_fill_type == 4)
      answer = pnpoly(polygon, P, idx, include_edges, eps);
   else
   if (polygon_fill_type == 3)
      // include_edges set to true because point can fall on edge of triangulation
      answer = InsidePolygon_triangulated(polygon, P, idx, sign, true, eps);
   else
   if (polygon_fill_type == 5)
      answer = InsidePolygon_solution1(polygon, P, idx, include_edges, eps);

   return answer;
}

int intersection_is_end_point(const vec3d &intersection_point, const vec3d &P0, const vec3d &P1, const double &eps)
{
   int which_one = 0;
   if (!compare(intersection_point,P0,eps))
      which_one = 1;
   else
   if (!compare(intersection_point,P1,eps))
      which_one = 2;
   return which_one;
}

// anywhere a vertex is on an edge, split that edge
bool mesh_verts(col_geom_v &geom, const double &eps)
{
   const vector<vec3d> &verts = geom.verts();
   const vector<vector<int> > &edges = geom.edges();

   // remember original sizes as the geom will be changing size
   int vsz = verts.size();
   int esz = edges.size();
   
   vector<int> deleted_edges;
   
   // compare only existing edges and verts
   for(int i=0; i<esz; i++) {
      vector<pair<double, int> > line_intersections;
      vec3d Q0 = verts[edges[i][0]];
      vec3d Q1 = verts[edges[i][1]];
      for(int v_idx=0; v_idx<vsz; v_idx++) {
         // don't compare to self
         if (edges[i][0]==v_idx || edges[i][1]==v_idx)
            continue;

         // find if P is in Q0,Q1
         vec3d P = verts[v_idx];
         vec3d intersection_point = point_in_segment(P, Q0, Q1, eps);
         if ((intersection_point.is_set()) && (!intersection_is_end_point(intersection_point, Q0, Q1, eps)))
            // store distance from Q0 to intersection and also the intersection vertex index
            line_intersections.push_back(make_pair((Q0-P).mag(),v_idx));
      }

      if (line_intersections.size()) {
         // edge i will be replaced. mark it for deletion
         deleted_edges.push_back(i);
         // sort based on distance from Q0
         sort( line_intersections.begin(), line_intersections.end() );
         // create edgelets from Q0 through intersection points to Q1 (using indexes)
         edge_into_geom(geom, edges[i][0], line_intersections[0].second, col_val::invisible);
         for(unsigned int k=0; k<line_intersections.size()-1; k++)
            edge_into_geom(geom, line_intersections[k].second, line_intersections[k+1].second, col_val::invisible);
         edge_into_geom(geom, line_intersections[line_intersections.size()-1].second, edges[i][1], col_val::invisible);
      }
   }
   
   // delete the replaced edges
   geom.delete_edges(deleted_edges);

   return (deleted_edges.size() ? true : false);
}

// of geom, is face i inside face j
bool is_face_inside_face(const geom_if &geom, const int &i, const int &j, const double &eps)
{
   const vector<vector<int> > &faces = geom.faces();
   const vector<vec3d> &verts = geom.verts();

   vector<int> face_idxs(1);
   face_idxs[0] = i;
   col_geom_v polygon1 = faces_to_geom(geom, face_idxs);
   const vector<int> &polygon2 = faces[j];

   bool answer = true;

   for(unsigned int i=0; i<polygon2.size(); i++) {
      vec3d P = verts[polygon2[i]];
      vec3d normal = face_norm(polygon1.verts(), polygon1.faces()[0]).unit();
      answer = is_point_inside_polygon(polygon1, P, normal, false, false, 1, eps);
      if (!answer)
         break;
   }

   return answer;
}

void check_for_holes(geom_if &geom, vector<pair<int, int> > &polygon_hierarchy, const double &eps)
{
   const vector<vector<int> > &faces = geom.faces();

   // there has to be at least 2 faces for there to be a hole
   if (faces.size() < 2)
      return;

   for(unsigned int i=0; i<faces.size(); i++) {
      for(unsigned int j=0; j<faces.size(); j++) {
         if (i == j)
            continue;
         if (is_face_inside_face(geom, i, j, eps)) {
            // store as child/parent
            polygon_hierarchy.push_back(make_pair(j,i));
         }
      }
   }

   sort( polygon_hierarchy.begin(), polygon_hierarchy.end() );
}

void make_hole_connectors(geom_if &geom, vector<vector<int> > &connectors, vector<pair<vec3d, vec3d> > &connectors_verts, const vector<pair<int, int> > &polygon_hierarchy)
{
   const vector<vector<int> > &faces = geom.faces();
   const vector<vec3d> &verts = geom.verts();
 
   int last_idx = -1;
   for(unsigned int i=0; i<polygon_hierarchy.size(); i++) {
      pair<int, int> hierarchy = polygon_hierarchy[i];
      // only allow child connection once
      if (hierarchy.first == last_idx)
         continue;
      int v_idx1 = faces[hierarchy.first][0];  // child
      int v_idx2 = faces[hierarchy.second][0]; // parent
      connectors.push_back(make_edge(v_idx1, v_idx2));
      connectors_verts.push_back(make_pair(verts[v_idx1], verts[v_idx2]));
      last_idx = hierarchy.first;
   }
}

void add_hole_connectors(col_geom_v &geom, const vector<vector<int> > &connectors)
{
   coloring clrng(&geom);
   
   col_val ecol = col_val::invisible;
   // uncomment for testing
   //ecol = col_val(1.0,1.0,1.0);
   
   for(unsigned int i=0; i<connectors.size(); i++) {
      geom.add_col_edge(connectors[i], ecol);
   }
}

// mark hole connectors for later edge processing. use color index of INT_MAX
void mark_hole_connectors(col_geom_v &geom, const vector<pair<vec3d, vec3d> > &connectors_verts, const double &eps)
{
   const vector<vec3d> &verts = geom.verts();
   const vector<vector<int> > &edges = geom.edges();

   for(unsigned int i=0; i<connectors_verts.size(); i++) {
      vec3d v1 = connectors_verts[i].first;
      vec3d v2 = connectors_verts[i].second;
      for(unsigned int j=0; j<edges.size(); j++) {
         // if edge is within connector, mark it connector color
         vec3d P1 = verts[edges[j][0]];
         vec3d P2 = verts[edges[j][1]];
         if ((point_in_segment(P1, v1, v2, eps)).is_set() && (point_in_segment(P2, v1, v2, eps)).is_set()) {
            col_val ecol;
            ecol.set_idx(INT_MAX);
            geom.set_e_col(j, ecol);
         }
      }
   }
}

// input seperate networks of overlapping edges and merge them into one network
bool mesh_edges(col_geom_v &geom, const double &eps)
{
   const vector<vec3d> &verts = geom.verts();
   const vector<vector<int> > &edges = geom.edges();

   // remember original sizes as the geom will be changing size
   int vsz = verts.size();
   int esz = edges.size();
   
   vector<int> deleted_edges;
   vector<pair<pair<int, int>, int> > new_verts;
   
   // compare only existing edges
   for(int i=0; i<esz; i++) {
      vector<pair<double, int> > line_intersections;
      for(int j=0; j<esz; j++) {
         // don't compare to self
         if (i==j)
            continue;

         // see if the new vertex was already created
         int v_idx = -1;
         for(unsigned int k=0; k<new_verts.size(); k++) {
            if (new_verts[k].first.first == i && new_verts[k].first.second == j) {
               v_idx = new_verts[k].second;
               break;
            }
         }
         
         // if it doesn't already exist, see if it needs to be created
         if (v_idx == -1) {
            // compare segements P0,P1 with Q0,Q1
            vec3d intersection_point = lines_intersection_in_segments(verts[edges[i][0]], verts[edges[i][1]], verts[edges[j][0]], verts[edges[j][1]], eps);
            if (intersection_point.is_set()) {
               // find (or create) index of this vertex
               v_idx = vertex_into_geom(geom, intersection_point, col_val::invisible, eps);
               // don't include existing vertices
               if (v_idx < vsz)
                  v_idx = -1;
               else {
                  // store index of vert at i,j. Reverse index i,j so it will be found when encountering edges j,i
                  pair<pair<int, int>, int> new_vert;
                  new_vert.first.first = j;
                  new_vert.first.second = i;
                  new_vert.second = v_idx;
                  new_verts.push_back(new_vert);
               }
            }
         }
         
         // if ultimately an intersection was found
         if (v_idx != -1)
            // store distance from P0 to intersection and also the intersection vertex index
            line_intersections.push_back(make_pair((verts[edges[i][0]]-verts[v_idx]).mag(),v_idx));
      }

      if (line_intersections.size()) {
         // edge i will be replaced. mark it for deletion
         deleted_edges.push_back(i);
         // sort based on distance from P0
         sort( line_intersections.begin(), line_intersections.end() );
         // create edgelets from P0 through intersection points to P1 (using indexes)
         edge_into_geom(geom, edges[i][0], line_intersections[0].second, col_val::invisible);
         for(unsigned int k=0; k<line_intersections.size()-1; k++)
            edge_into_geom(geom, line_intersections[k].second, line_intersections[k+1].second, col_val::invisible);
         edge_into_geom(geom, line_intersections[line_intersections.size()-1].second, edges[i][1], col_val::invisible);
      }
   }
   
   // delete the replaced edges
   geom.delete_edges(deleted_edges);

   return (deleted_edges.size() ? true : false);
}

void delete_duplicate_index_edges(col_geom_v &geom)
{
   const vector<vector<int> > &edges = geom.edges();

   // delete edges in the form of 2 x x (where edge indexes are equal)
   vector<int> deleted_edges;
   for(unsigned int i=0; i<edges.size(); i++) {
      if (edges[i][0] == edges[i][1])
         deleted_edges.push_back(i);
   }
   geom.delete_edges(deleted_edges);
}

// can control color of skeleton. Here it is made invisible
void make_skeleton(col_geom_v &geom)
{
   geom.add_missing_impl_edges();
   geom.clear_faces();

   // make skeleton invisible
   coloring clrng(&geom);
   
   col_val ecol = col_val::invisible;
   // uncomment for testing
   //ecol = col_val(1.0,1.0,1.0);
   
   clrng.e_one_col(ecol);
   clrng.v_one_col(ecol);
}

// find connections from a vertex v_idx
void find_connections(const geom_if &geom, vector<int> &vcons, const int &v_idx)
{
   const vector<vector<int> > &edges = geom.edges();

   vcons.clear();
   for(unsigned int i=0; i<edges.size(); i++) {
      if (edges[i][0] == v_idx )
         vcons.push_back(edges[i][1]);
      else
      if (edges[i][1] == v_idx )
         vcons.push_back(edges[i][0]);
   }
}

void build_angle_map(const geom_if &geom, map<pair<int, int>, double> &angle_map, const vec3d &normal, const double &eps)
{
   const vector<vec3d> &verts = geom.verts();
   
   int idx;
   int sign;

   project_using_normal(normal, idx, sign);

   int idx0 = (idx+1)%3;
   int idx1 = (idx+2)%3;
   
   for(unsigned int i=0; i<verts.size(); i++) {
      vector<int> vcons;
      find_connections(geom,vcons,i);
      for(unsigned int j=0; j<vcons.size(); j++) {
         int k = vcons[j];
         double y = verts[k][idx1] - verts[i][idx1];
         double x = verts[k][idx0] - verts[i][idx0];
         double angle = rad2deg(atan2(y,x));
         // eliminate small error
         if (double_eq(angle, 0.0, eps))
            angle = 0.0;
         // make all angles positive
         if (angle < 0.0)
            angle += 360;
         angle_map[make_pair(i,k)] = angle;
      }
   }
}

void build_turn_map(const geom_if &geom, map<pair<int, int>, int> &turn_map, map<pair<int, int>, double> &angle_map)
{
   const vector<vector<int> > &edges = geom.edges();
 
   for(unsigned int i=0; i<edges.size(); i++) {
      for(unsigned int j=0; j<2; j++) { 
         vector<int> vcons;
         int a = edges[i][!j ? 0 : 1];
         int b = edges[i][!j ? 1 : 0];

         double base_angle = angle_map[make_pair(b,a)];
         find_connections(geom,vcons,b);
         vector<pair<double, int> > angles;
         for(unsigned int k=0; k<vcons.size(); k++) {
            int c = vcons[k];
            // can't go backwards unless it is the only possible path
            if (c == a && vcons.size() != 1)
               continue;
            double angle = angle_map[make_pair(b,c)];
            if (angle >= 0 && angle < base_angle)
               angle += 360;
            angles.push_back(make_pair(angle,c));
         }
         // the minimum angle is now the next higher angle than base angle. it will be at the top of the sort
         // map that turn for AB -> C
         sort( angles.begin(), angles.end() );
         turn_map[make_pair(a, b)] = angles[0].second;
      }
   }
}

void get_first_unvisited_triad(const geom_if &geom, int &first_unvisited, map<pair<int, int>, bool > &visited, map<pair<int, int>, int> &turn_map, vector<int> &face)
{
   const vector<vector<int> > &edges = geom.edges();
 
   face.clear();
   for(unsigned int i=first_unvisited; i<edges.size(); i++) {
      for(unsigned int j=0; j<2; j++) { 
         int a = edges[i][!j ? 0 : 1];
         int b = edges[i][!j ? 1 : 0];
         if (visited[make_pair(a,b)])
            continue;
         pair<int, int> edge_pair = make_pair(a,b);
         int c = turn_map[edge_pair];
         if (!visited[make_pair(b,c)]) {
            face.push_back(a);
            face.push_back(b);
            face.push_back(c);
            first_unvisited = i;
            // if j is 1 this edge has been traversed in both directions
            if (j)
               first_unvisited++;
            return;
         }
      }
   }
   return;
}

void construct_faces(geom_if &geom, map<pair<int, int>, int> &turn_map)
{
   map<pair<int, int>, bool > visited;
   
   int first_unvisited = 0;
   vector<int> face;
   get_first_unvisited_triad(geom,first_unvisited,visited,turn_map,face);
   while(face.size()) {
      int fsz = face.size();
      int a = face[fsz-2];
      int b = face[fsz-1];
      int c = turn_map[make_pair(a,b)];
      // if c == face[0], still consider face[0] is part of complex polygon if next turn does not lead to face[1]
      if (c != face[0] || turn_map[make_pair(b,c)] != face[1])
         face.push_back(c);
      else {
         // RK: patch. Do not allow faces to have sequential duplicate indexes
         vector<int>::iterator fi = unique(face.begin(), face.end());
         face.resize( fi - face.begin() );

         geom.add_face(face);
         // mark turns visited and start a new face
         fsz = face.size();
         for(int i=0; i<fsz; i++)
            visited[make_pair(face[i],face[(i+1)%fsz])] = true;
         get_first_unvisited_triad(geom,first_unvisited,visited,turn_map,face);
      }
   }
}

void analyze_faces(geom_if &geom, const int &planar_merge_type, vector<int> &nonconvex_faces, map<pair<int, int>, double> &angle_map)
{
   const vector<vector<int> > &faces = geom.faces();
   vector<int> deleted_faces;
   int deleted_faces_count = 0; // should be only one
   
   for(unsigned int i=0; i<faces.size(); i++) {
      double angle_sum = 0;
      bool all_negative_turns = true;
      int fsz = faces[i].size();
      for(int j=0; j<fsz; j++) {
         int a = faces[i][j];
         int v = faces[i][(j+1)%fsz];
         int b = faces[i][(j+2)%fsz];
         double angle = angle_map[make_pair(v,b)]-angle_map[make_pair(v,a)]-180.0;
         if (angle < -180.0)
            angle += 360.0;
         else
         if (angle > 180.0)
            angle -= 360.0;
         if (angle > 0.0)
            all_negative_turns = false;
         angle_sum += angle;
      }
      int sum = (int)floorf(angle_sum + 0.5);
      // old statement: if ((sum == 360 && planar_merge_type == 1) || (sum != 360 && planar_merge_type == 2)) {
      // ideally the interior angle sum of the tiles is -360 degrees. However it is known to be otherwise. Check less than 0
      // the angle sum of the perimeter face is ideally 360. However it is known to sometimes be otherwise. Check greater than or equal to 0
      if ((sum >= 0 && planar_merge_type == 1) || (sum < 0 && planar_merge_type == 2)) {
         deleted_faces.push_back(i);
         deleted_faces_count++;
      }
      else
      if (!all_negative_turns)
         nonconvex_faces.push_back(i-deleted_faces_count);
   }
   geom.delete_faces(deleted_faces);
}

void fill_in_faces(col_geom_v &geom, const int &planar_merge_type, vector<int> &nonconvex_faces, const vec3d &normal, const double &eps)
{
   map<pair<int, int>, double> angle_map;
   map<pair<int, int>, int> turn_map;
   build_angle_map(geom, angle_map, normal, eps);
   build_turn_map(geom, turn_map, angle_map);
   construct_faces(geom, turn_map);
   analyze_faces(geom, planar_merge_type, nonconvex_faces, angle_map);
}

col_val average_color(const vector<col_val> &cols, const int &color_system_mode,
                      const double &sat_power, const double &sat_threshold, const double &value_power, const double &value_advance, const int &alpha_mode, const bool &ryb_mode)
{
   col_val avg_col;

   if (color_system_mode < 3)
      avg_col = blend_HSX_centroid(cols, color_system_mode, sat_power, sat_threshold, value_power, value_advance, alpha_mode, ryb_mode);
   else
   if (color_system_mode == 3) {
      avg_col = blend_RGB_centroid(cols, alpha_mode, ryb_mode);
   }

   return avg_col;
}

// all values can be locally altered
bool winding_rule_filter(int winding_rule_mode, int winding_rule, int winding_number)
{
   bool answer = true;

   // default is accept all winding numbers
   if (winding_rule_mode == INT_MAX)
      return answer;

   // negative winding rule mode means do absolute value comparison
   if (winding_rule_mode < 0) {
      winding_rule_mode = abs(winding_rule_mode);
      winding_rule = abs(winding_rule);
      winding_number = abs(winding_number);
   }

   if (winding_rule_mode == 0) {
      // odd=1  even=2  positive=3  negative=4  nonzero=5
      // zodd=6  zeven=7  zpositive=8  znegative=9  zero=10
      // rule 5 implicitly done
      // when rule 5 or less zero is not included (probably couldn't get here)
      if ((winding_rule < 6) && !winding_number)
         answer = false;
      else
      if ((winding_rule == 1) && (is_even(winding_number)))
         answer = false;
      else
      if ((winding_rule == 2 || winding_rule == 7) && (!is_even(winding_number)))
         answer = false;
      else
      if ((winding_rule == 3 || winding_rule == 8) && (winding_number < 0))
         answer = false;
      else
      if ((winding_rule == 4 || winding_rule == 9) && (winding_number > 0))
         answer = false;
      else
      if ((winding_rule == 6) && (winding_number && is_even(winding_number)))
         answer = false;
      else
      if ((winding_rule == 10) && (winding_number))
         answer = false;
   }
   else
   // symbol proceeding integer. Symbols: eq, ne, gt, ge, lt, le"
   // correspond to 1,2,3,4,5,6 (negative if absolute value)
   // equal to
   if (winding_rule_mode == 1)
      answer = (winding_number == winding_rule) ? true : false;
   else
   // not equal to
   if (winding_rule_mode == 2)
      answer = (winding_number != winding_rule) ? true : false;
   else
   // greater than
   if (winding_rule_mode == 3)
      answer = (winding_number > winding_rule) ? true : false;
   else
   // greater than or equal to
   if (winding_rule_mode == 4)
       answer = (winding_number >= winding_rule) ? true : false;
   else
   // less than
   if (winding_rule_mode == 5)
       answer = (winding_number < winding_rule) ? true : false;
   else
   // less than or equal to
   if (winding_rule_mode == 6)
       answer = (winding_number <= winding_rule) ? true : false;

   return answer;
}

void sample_colors(col_geom_v &sgeom, const col_geom_v &cgeom, const int &planar_merge_type, const vector<xnormal> &original_normals, const vector<int> &nonconvex_faces, const vec3d &center,
   const int &polygon_fill_type, const int &winding_rule_mode, const int &winding_rule, const char &color_by_winding_number, const bool &find_direction, const bool &winding_div2,
   col_val &zero_density_color, const bool &zero_density_force_blend, const double &brightness_adj,
   const int &color_system_mode, const double &sat_power, const double &sat_threshold, const double &value_power, const double &value_advance, const int &alpha_mode, const bool &ryb_mode,
   int &winding_total_min, int &winding_total_max, const double &eps)
{
   const vector<vector<int> > &sfaces = sgeom.faces();
   const vector<vector<int> > &cfaces = cgeom.faces();

   vector<col_val> cols;
   col_val average_color_all_faces;
   for(unsigned int i=0; i<cfaces.size(); i++)
      cols.push_back(cgeom.get_f_col(i));
   average_color_all_faces = average_color(cols, color_system_mode, sat_power, sat_threshold, value_power, value_advance, alpha_mode, ryb_mode);
   cols.clear();

   if (zero_density_force_blend)
      zero_density_color = average_color_all_faces;

   for(unsigned int i=0; i<sfaces.size(); i++) {
      vector<vec3d> points;

      // if the sampling face is convex
      if(find(nonconvex_faces.begin(), nonconvex_faces.end(), (int)i) == nonconvex_faces.end())
         points.push_back(sgeom.face_cent(i));
      else {
      // else non-convex. need to triangulate the non-convex sample polygon to sample on centroid(s) of the triangles (until one is hit)
         vector<int> sface_idxs;
         sface_idxs.push_back(i);
         col_geom_v spolygon = faces_to_geom(sgeom, sface_idxs);
         spolygon.triangulate();
         if (!spolygon.faces().size()) {
            //trangulation of a polygon of zero density leaves no faces
            points.push_back(centroid(spolygon.verts()));
         }
         else {
            // when merging and it is a nonconvex face, then have to sample all the centroids. else only sample one of them
            int sz = (planar_merge_type==1) ? 1 : spolygon.faces().size();
            for(int k=0; k<sz; k++)
               points.push_back(spolygon.face_cent(k));
         }
      }

      // accumulate winding numbers
      int winding_total = 0;

      // put colored faces to sample (one at a time) into geom named polygon
      for(unsigned int j=0; j<cfaces.size(); j++) {
         vector<int> face_idxs;
         face_idxs.push_back(j);
         col_geom_v polygon = faces_to_geom(cgeom, face_idxs);

         col_geom_v tpolygon; // needed for triangulation method. the sample polygon needs to be triangulated
         if (polygon_fill_type == 3) {
            tpolygon = polygon;
            tpolygon.triangulate();
         }

         xnormal original_normal =  original_normals[j];
         vec3d normal = original_normal.unit();

         // if merging we have to sample all the centroids until there is a hit. otherwise k will begin and end at 0
         for(unsigned int k=0; k<points.size(); k++) {
            bool answer = is_point_inside_polygon((polygon_fill_type == 3 ? tpolygon : polygon), points[k], normal, true, false, polygon_fill_type, eps);
            if (answer) {
               vector<vec3d> one_point;
               one_point.push_back(points[k]);
               int winding_number = get_winding_number(polygon, one_point, original_normal, find_direction && original_normal.is_hemispherical(), eps);

               // if merging, find largest magnitude of winding number
               // if they are -W and +W, chose the positive
               if (planar_merge_type==2) {
                  if ((abs(winding_number) > abs(winding_total)) || ((winding_number > 0) && (winding_number == -winding_total)))
                     winding_total = winding_number;
               }
               // otherwise just total
               else
                  winding_total += winding_number;

               cols.push_back(cgeom.get_f_col(j));
               break;
            }
         }
      }

      // correct the winding number of the new face. (reverse face if necessary)
      vector<int> face_idxs;
      face_idxs.push_back(i);
      col_geom_v polygon = faces_to_geom(sgeom, face_idxs);
      vector<vec3d> one_point;
      one_point.push_back(sgeom.face_cent(i));
      xnormal sface_normal(sgeom, i, center, eps);
      int winding_number = get_winding_number(polygon, one_point, sface_normal, find_direction && sface_normal.is_hemispherical(), eps);
      if ((winding_number < 0 && winding_total > 0) || (winding_number > 0 && winding_total < 0))
         reverse(sgeom.raw_faces()[i].begin(), sgeom.raw_faces()[i].end());

      if (winding_div2)
         winding_total = (winding_total+1)/2;

      if (winding_total < winding_total_min)
         winding_total_min = winding_total;
      if (winding_total > winding_total_max)
         winding_total_max = winding_total;

      if ((winding_rule != INT_MAX)) {
         // if cols.size() is not zero then there were hits
         if (cols.size()) {
            if (!winding_rule_filter(winding_rule_mode, winding_rule, winding_total))
               cols.clear();
         }
         // the magic. if there are no hits
         // if winding rule is less than 5 then zero is included so color it average
         else
         if ((winding_rule_mode == 0 && winding_rule > 5) || (winding_rule_mode > 0 && winding_rule == 0))
            cols.push_back(average_color_all_faces);
      }

      int sz = cols.size();

      col_val col;

      if (color_by_winding_number && sz) {
         int wtotal = winding_total;
         // if absolute value of winding numbers (or their negative) are colored
         if (color_by_winding_number != 'w') {
            wtotal = abs(wtotal);
            if (color_by_winding_number == 'n')
               wtotal = -wtotal;
         }

         if (wtotal >= 0)
            col.set_idx(wtotal);
         else
            // negative winding numbers are set up to near INT_MAX to be subtracted out later
            col.set_idx(wtotal + INT_MAX);
      }
      else
      // if there is no hit, then that patch is of zero density color. if file_type 4 then all even numbered patches are also
      if (!sz || (polygon_fill_type == 4 && !(sz%2)))
         col = col_val(zero_density_color);
      else
         col = average_color(cols, color_system_mode, sat_power, sat_threshold, value_power, value_advance, alpha_mode, ryb_mode);

      if ((brightness_adj > -2.0) && col.is_set() && !col.is_inv() && (sz > 1)) {
         // add 0.5 to brightness so it is in the range of -0.5 to 1.5
         // if sz = 2, 1.0/2 = 0.5 so when -B is 0, a patch with 2 colors blended will not be changed 
         double brightness = (1.0/sz)+(brightness_adj + 0.5);
         col.set_brightness(-1.0+brightness);
      }

      sgeom.set_f_col(i,col);
      cols.clear();
   }
}

void collect_original_normals(vector<xnormal> &original_normals, const vector<int> &coplanar_face_list, const fnormals &fnormals)
{
   for(unsigned int i=0; i<coplanar_face_list.size(); i++)
      original_normals.push_back(fnormals[coplanar_face_list[i]]);
}

void blend_overlapping_faces(col_geom_v &geom, const vector<vector<int> > &coplanar_faces_list, const vector<xnormal> &coplanar_normals, const fnormals &fnormals, const vec3d &center,
   const int &planar_merge_type, const int &polygon_fill_type,
   const bool &hole_detection, const int &winding_rule_mode, const int &winding_rule, const char &color_by_winding_number, const bool &find_direction, const bool &winding_div2, const bool &verbose,
   col_val &zero_density_color, const bool &zero_density_force_blend, const double &brightness_adj,
   const int &color_system_mode, const double &sat_power, const double &sat_threshold, const double &value_power, const double &value_advance, const int &alpha_mode, const bool &ryb_mode, const double &eps)
{
   col_geom_v bgeom;
   vector<int> deleted_faces;

   int winding_number_min = INT_MAX;
   int winding_number_max = INT_MIN;
   
   for(unsigned int i=0; i<coplanar_faces_list.size(); i++) {
      // load a geom with color faces. keep it and copy it.
      col_geom_v cgeom = faces_to_geom(geom, coplanar_faces_list[i]);
      col_geom_v sgeom = cgeom;

      // check here for polygons within polygons
      vector<vector<int> > connectors;
      vector<pair<vec3d, vec3d> > connectors_verts;
      if (hole_detection) {
         vector<pair<int, int> > polygon_hierarchy;
         check_for_holes(sgeom, polygon_hierarchy, eps);
         if (polygon_hierarchy.size())
            make_hole_connectors(sgeom, connectors, connectors_verts, polygon_hierarchy);
      }

      make_skeleton(sgeom);

      // edges with duplicate indexes can happen if faces have duplicate sequential indexes
      delete_duplicate_index_edges(geom);

      if (connectors.size())
         add_hole_connectors(sgeom, connectors);
      connectors.clear();

      // duplicate vertices and edges can cause problems
      sort_merge_elems(sgeom, "ve", 0, eps);

      // sort merge can destill more duplicate indexes
      delete_duplicate_index_edges(sgeom);

      mesh_verts(sgeom, eps);
      mesh_edges(sgeom, eps);

      // have to use vertex location for marking because indexes have been scrambled
      if (connectors_verts.size())
         mark_hole_connectors(sgeom, connectors_verts, eps);
      connectors_verts.clear();

      vector<int> nonconvex_faces;
      fill_in_faces(sgeom, planar_merge_type, nonconvex_faces, coplanar_normals[i].outward().unit(), eps);
  
      // original normals are needed for sampling colors
      vector<xnormal> original_normals;
      collect_original_normals(original_normals, coplanar_faces_list[i], fnormals);

      int winding_total_min = INT_MAX;
      int winding_total_max = INT_MIN;
      sample_colors(sgeom, cgeom, planar_merge_type, original_normals, nonconvex_faces, center,
                    polygon_fill_type, winding_rule_mode, winding_rule, color_by_winding_number, find_direction, winding_div2,
                    zero_density_color, zero_density_force_blend, brightness_adj,
                    color_system_mode, sat_power, sat_threshold, value_power, value_advance, alpha_mode, ryb_mode,
                    winding_total_min, winding_total_max, eps);

      if (winding_total_min < winding_number_min)
         winding_number_min = winding_total_min;
      if (winding_total_max > winding_number_max)
         winding_number_max = winding_total_max;
 
      bgeom.append(sgeom);

      // mark the faces in the cluster for deletion at the end
      deleted_faces.insert(deleted_faces.end(), coplanar_faces_list[i].begin(), coplanar_faces_list[i].end());
   }

   geom.delete_faces(deleted_faces);
   geom.append(bgeom);

   if (verbose) {
      fprintf(stderr,"minimum winding number = %d\n",winding_number_min);
      fprintf(stderr,"maximum winding number = %d\n",winding_number_max);
   }
}

void planar_merge(col_geom_v &geom, const int &planar_merge_type, const int &polygon_fill_type, const bool &hole_detection, const vec3d &center,
   const int &winding_rule_mode, const int &winding_rule, const char &color_by_winding_number, const bool &find_direction, const bool &winding_div2, const bool &verbose,
   col_val &zero_density_color, const bool &zero_density_force_blend, const double &brightness_adj,
   const int &color_system_mode, const double &sat_power, const double &sat_threshold, const double &value_power, const double &value_advance, const int &alpha_mode, const bool &ryb_mode, const double &eps)
{
   fnormals face_normals(geom, center, eps);

   vector<vector<int> > coplanar_faces_list;
   vector<xnormal> coplanar_normals;
   
   bool point_outward = true;
   bool fold_normals = false;
   bool fold_normals_hemispherical = true;
   bool filtered = true;
   build_coplanar_faces_list(geom, coplanar_faces_list, coplanar_normals, face_normals, point_outward, fold_normals, fold_normals_hemispherical, filtered, eps);

   blend_overlapping_faces(geom, coplanar_faces_list, coplanar_normals, face_normals, center,
                           planar_merge_type, polygon_fill_type,
                           hole_detection, winding_rule_mode, winding_rule, color_by_winding_number, find_direction, winding_div2, verbose,
                           zero_density_color, zero_density_force_blend, brightness_adj,
                           color_system_mode, sat_power, sat_threshold, value_power, value_advance, alpha_mode, ryb_mode, eps);
}

void color_by_plane(col_geom_v &geom, const char &face_color_method, const color_map_multi &map, const vec3d &center, const double &eps)
{
   fnormals fnormals(geom, center, eps);

   vector<vector<int> > coplanar_faces_list;
   vector<xnormal> coplanar_normals;
  
   //bool point_outward = true; // always do this
   bool point_outward = (face_color_method == 'n') ? false : true;
   bool fold_normals_hemispherical = (face_color_method == 'n') ? false : true;
   bool filtered = (face_color_method == 'p') ? true : false;
   bool fold_normals = (face_color_method == 'o') ? true : false;
   build_coplanar_faces_list(geom, coplanar_faces_list, coplanar_normals, fnormals, point_outward, fold_normals, fold_normals_hemispherical, filtered, eps);
   
   for(unsigned int i=0; i<coplanar_faces_list.size(); i++) {
      for(unsigned int j=0; j<coplanar_faces_list[i].size(); j++) {
         int k = coplanar_faces_list[i][j];
         col_val col = map.get_col(i);
         geom.set_f_col(k,col);
      }
   }
}

void build_colinear_edge_list(const geom_if &geom, const vector<vector<int> > &edges, vector<vector<int> > &colinear_edge_list, const double &eps)
{
   rand_gen rnd;
   rnd.time_seed();
   vec3d C = vec3d(vec3d::random(rnd)).unit();

   vector<pair<vec3d, int> > edge_unit_nearpoints;
   for(unsigned int i=0;i<edges.size();i++) {
      vector<int> edge = edges[i];
      //vec3d P = nearest_point(C, verts[edge[0]], verts[edge[1]]).unit();
      vec3d P = (geom.edge_nearpt(edge, C)).unit();
      edge_unit_nearpoints.push_back(make_pair(P,i));
   }

   if (!edge_unit_nearpoints.size())
      return;

   stable_sort( edge_unit_nearpoints.begin(), edge_unit_nearpoints.end(), vert_cmp(eps) ); 

   vector<int> colinear_edges;
   // at least one edge is on a line
   colinear_edges.push_back(edge_unit_nearpoints[0].second);
   // find edges on same line and group them together
   for(unsigned int i=1; i<edge_unit_nearpoints.size(); i++) {
      if (compare(edge_unit_nearpoints[i].first,edge_unit_nearpoints[i-1].first,eps)) {
         colinear_edge_list.push_back(colinear_edges);
         colinear_edges.clear();
      }
      colinear_edges.push_back(edge_unit_nearpoints[i].second);
   }
   colinear_edge_list.push_back(colinear_edges);
   colinear_edges.clear();

   edge_unit_nearpoints.clear();
}

vector<int> find_end_points_vertex(const geom_if &geom, const vector<int> &vert_indexes, const double &eps)
{
   const vector<vec3d> &verts = geom.verts();

   col_geom_v vgeom;
   for(unsigned int i=0;i<vert_indexes.size();i++)
      vertex_into_geom(vgeom, verts[vert_indexes[i]], col_val::invisible, eps);
   vgeom.set_hull();

   const vector<vec3d> &gverts = vgeom.verts();
   vec3d end_point1 = gverts[0];
   vec3d end_point2 = gverts[1];
   // patch: unfortunately the verts() size can somtimes, due to math error, come out greater than 2
   // then have to take a vertex from longest edge. vert of longest edge will be at the bottom of the list
   int sz = gverts.size();
   if (sz > 2) {
      vector<pair<double, int> > distance_table;
      for(unsigned int i=0;i<gverts.size();i++) {
         double dist = (gverts[i]-gverts[(i+1)%sz]).mag();
         distance_table.push_back(make_pair(dist,i));
      }
      sort( distance_table.begin(), distance_table.end() );
      end_point1 = gverts[distance_table[sz-1].second];
      end_point2 = gverts[distance_table[sz-2].second];
   }
   vgeom.clear_all();

   vector<int> end_indexes;
   for(unsigned int i=0;i<vert_indexes.size();i++) {
      if (!compare(end_point1, verts[vert_indexes[i]], eps) || !compare(end_point2, verts[vert_indexes[i]], eps)) {
         end_indexes.push_back(vert_indexes[i]);
      }
   }

   return end_indexes;   
}

void collect_ordered_vert_indexes(const geom_if &geom, const vector<vector<int> > &edges, const vector<int> &colinear_edges, vector<int> &colinear_verts, const double &eps)
{
   const vector<vec3d> &verts = geom.verts();

   vector<int> vert_indexes;
   for(unsigned int i=0;i<colinear_edges.size();i++) {
      vector<int> edge = edges[colinear_edges[i]];
      for(unsigned int j=0;j<edge.size();j++)
         vert_indexes.push_back(edge[j]);
   }

   sort( vert_indexes.begin(), vert_indexes.end() );    
   vector<int>::iterator vi = unique(vert_indexes.begin(), vert_indexes.end());
   vert_indexes.resize( vi - vert_indexes.begin() );

   vector<int> end_indexes = find_end_points_vertex(geom, vert_indexes, eps);
   int end_idx = end_indexes[0];
   vec3d end_vertex = verts[end_idx];

   vector<pair<double, int> > distance_table;
   distance_table.push_back(make_pair(0.0, end_idx));
   for(unsigned int i=0;i<vert_indexes.size();i++) {
      if (vert_indexes[i] == end_idx)
         continue;
      double dist = (verts[vert_indexes[i]]-end_vertex).mag();
      distance_table.push_back(make_pair(dist,vert_indexes[i]));
   }
   vert_indexes.clear();

   sort( distance_table.begin(), distance_table.end() );

   for(unsigned int i=0;i<distance_table.size();i++)
      colinear_verts.push_back(distance_table[i].second);
}

void build_colinear_vertex_list(const geom_if &geom, vector<vector<int> > &colinear_vertex_list, const double &eps)
{
   vector<vector<int> > implicit_edges;
   geom.get_impl_edges(implicit_edges);

   vector<vector<int> > colinear_edge_list;
   build_colinear_edge_list(geom, implicit_edges, colinear_edge_list, eps);

   for(unsigned int i=0;i<colinear_edge_list.size();i++) {
      vector<int> colinear_verts;
      collect_ordered_vert_indexes(geom, implicit_edges, colinear_edge_list[i], colinear_verts, eps);
      colinear_vertex_list.push_back(colinear_verts);
      colinear_verts.clear();
   }
}

// inserts added_vertices into face[face_idx] between start_v_idx and end_v_idx.
// if end_v_idx comes before start_v_idx, reverse added_vertices
void insert_verts_into_face(geom_if &geom, vector<int> &added_vertices, const int &face_idx, const int &start_v_idx, const int &end_v_idx)
{
   vector<vector<int> > &faces = geom.raw_faces();
   vector<int>::iterator iter1 = find(faces[face_idx].begin(), faces[face_idx].end(), start_v_idx);
   vector<int>::iterator iter2 = find(faces[face_idx].begin(), faces[face_idx].end(), end_v_idx);

   // this is how the first and last element is designated
   vector<int>::iterator first = faces[face_idx].begin();
   vector<int>::iterator last = faces[face_idx].end();
   last--;

   if (!(*first == end_v_idx && *last == start_v_idx) && ((iter2 < iter1) || (*first == start_v_idx && *last == end_v_idx))) {
      iter1 = iter2;
      reverse(added_vertices.begin(), added_vertices.end());
   }

   faces[face_idx].insert((iter1+1), added_vertices.begin(), added_vertices.end());
}

// adds vertices to faces on seams where there is a mismatch in vertices from other face's edges
// may close the polyhedron
void stitch_faces_on_seams(col_geom_v &geom, const double &eps)
{
   const vector<vector<int> > &faces = geom.faces();
   const vector<vec3d> &verts = geom.verts();

   vector<vector<int> > colinear_vertex_list;
   build_colinear_vertex_list(geom, colinear_vertex_list, eps);

   // find what vertices are in faces
   vector<vector<int> > vert_has_faces;
   for(unsigned int i=0;i<verts.size();i++)
      vert_has_faces.push_back(find_faces_with_vertex(faces, i));

   for(unsigned int i=0;i<colinear_vertex_list.size();i++) {
      vector<int> colinear_verts = colinear_vertex_list[i];
      int sz = colinear_verts.size();
      // if there are only two verts there can be no face between them
      for(int j=0;j<sz-2;j++) {
         int start_v_idx = colinear_verts[j];
         vector<int> vert_has_faces_start = vert_has_faces[start_v_idx];
         for(unsigned int k=0;k<vert_has_faces_start.size();k++) {
            int face_idx = vert_has_faces_start[k];
            vector<int> added_vertices;
            for(int l=j+1;l<sz;l++) {
               int test_v_idx = colinear_verts[l];
               vector<int> vert_has_faces_test = vert_has_faces[test_v_idx];
               if(find(vert_has_faces_test.begin(), vert_has_faces_test.end(), face_idx) == vert_has_faces_test.end()) {
                  added_vertices.push_back(test_v_idx);
               }
               else {
                  if (added_vertices.size())
                     insert_verts_into_face(geom, added_vertices, face_idx, start_v_idx, test_v_idx);
                  break;
               }
            }
            added_vertices.clear();
         }
      }
   }
}

void remove_in_line_vertices(col_geom_v &geom, const double &eps)
{
  const vector<vec3d> &verts = geom.verts();

   vector<vector<int> > implicit_edges;
   geom.get_impl_edges(implicit_edges);

   vector<vector<int> > colinear_vertex_list;
   build_colinear_vertex_list(geom, colinear_vertex_list, eps);

   vector<int> end_points;
   for(unsigned int i=0;i<colinear_vertex_list.size();i++) {
      vector<int> end_indexes = find_end_points_vertex(geom, colinear_vertex_list[i], eps);
      end_points.push_back(end_indexes[0]);
      end_points.push_back(end_indexes[1]);
   }

   sort( end_points.begin(), end_points.end() );    
   vector<int>::iterator ei = unique(end_points.begin(), end_points.end());
   end_points.resize( ei - end_points.begin() );

   vector<int> deleted_verts;
   for(unsigned int i=0;i<verts.size();i++) {
      if(find(end_points.begin(), end_points.end(), i) == end_points.end())
         deleted_verts.push_back(i);
   }
   end_points.clear();

   geom.delete_verts(deleted_verts);
}

bool compare_edge_verts(col_geom_v &geom, const int &i, const int &j, const double &eps)
{
   const vector<vector<int> > &edges = geom.edges();
   const vector<vec3d> &verts = geom.verts();

   vec3d v_i0 = verts[edges[i][0]];
   vec3d v_i1 = verts[edges[i][1]];
   vec3d v_j0 = verts[edges[j][0]];
   vec3d v_j1 = verts[edges[j][1]];

   return ((!compare(v_i0,v_j0,eps) && !compare(v_i1,v_j1,eps)) || (!compare(v_i0,v_j1,eps) && !compare(v_i1,v_j0,eps)));
}

// set first color to the blend
// for sort_merge to reduce to first color for edges and/or verts
string pre_edge_blend(col_geom_v &geom, const char &edge_blending, const int &color_system_mode,
                      const double &sat_power, const double &sat_threshold, const double &value_power, const double &value_advance, const int &alpha_mode, const bool &ryb_mode, const double &eps)
{
   const vector<vector<int> > &edges = geom.edges();
   const vector<vec3d> &verts = geom.verts();

   string elems = "";

   if (edge_blending == 'e' || edge_blending == 'b') {
      int sz = edges.size();
      vector<bool> used(sz);
      for(int i=0; i<sz; i++) {
         if (used[i])
            continue;
         vector<col_val> cols;
         cols.push_back(geom.get_e_col(i));
         for(int j=i+1; j<sz; j++) {
            if (used[j])
               continue;
            if (compare_edge_verts(geom,i,j,eps)) {
               cols.push_back(geom.get_e_col(j));
               used[j] = true;
            }
         }
         col_val col = average_color(cols, color_system_mode, sat_power, sat_threshold, value_power, value_advance, alpha_mode, ryb_mode);
         geom.set_e_col(i,col);
      }

      elems += 'e';
   }

   if (edge_blending == 'v' || edge_blending == 'b') {
      int sz = verts.size();
      vector<bool> used(sz);
      for(int i=0; i<sz; i++) {
         if (used[i])
            continue;
         vector<col_val> cols;
         cols.push_back(geom.get_v_col(i));
         for(int j=i+1; j<sz; j++) {
            if (used[j])
               continue;
            if (!compare(verts[i],verts[j],eps)) {
               cols.push_back(geom.get_v_col(j));
               used[j] = true;
            }
         }
         col_val col = average_color(cols, color_system_mode, sat_power, sat_threshold, value_power, value_advance, alpha_mode, ryb_mode);
         geom.set_v_col(i,col);
      }

      elems += 'v';
   }

   return elems;
}

// similar to edge_blending, but wait till planar_merge to get new invisible edgelets to color blend
string post_edge_blend(col_geom_v &geom, const char &edge_blending, const int &original_edges_size, const int &color_system_mode,
                       const double &sat_power, const double &sat_threshold, const double &value_power, const double &value_advance, const int &alpha_mode, const bool &ryb_mode, const double &eps)
{
   const vector<vector<int> > &edges = geom.edges();
   const vector<vec3d> &verts = geom.verts();

   string elems = "";

   // Note: after tiling or merging faces, it makes no sense not to merge the edges and vertices
   // so the program forces 'b' for edge_blending

   // Patch: if edge_blending == 'v' then just use first encounted edge color but sort_merge anyway
   // this is because antiview displays duplicate edges with the last color
   // another solution would be to cut the original edges and paste them to the end and use e,v,b as normal
   if (true) {
   //if (edge_blending == 'e' || edge_blending == 'b') {
      vector<vector<int> > colinear_edge_list;
      build_colinear_edge_list(geom, edges, colinear_edge_list, eps);

      for(unsigned int i=0;i<colinear_edge_list.size();i++) {
         vector<int> original_edge_idx;
         vector<int> added_edge_idx;
         for(unsigned int j=0;j<colinear_edge_list[i].size();j++) {
            int edge_idx = colinear_edge_list[i][j];
            if (edge_idx<original_edges_size)
               original_edge_idx.push_back(edge_idx);
            else
               added_edge_idx.push_back(edge_idx);
         }

         for(unsigned int j=0;j<added_edge_idx.size();j++) {
            int added_edge = added_edge_idx[j];
            col_val col = geom.get_e_col(added_edge);
            if (col.is_idx() && col == INT_MAX)
               continue;

            vec3d P1 = verts[edges[added_edge][0]];
            vec3d P2 = verts[edges[added_edge][1]];

            vector<col_val> cols;
            for(unsigned int k=0;k<original_edge_idx.size();k++) {
               int original_edge = original_edge_idx[k];
               vec3d v1 = verts[edges[original_edge][0]];
               vec3d v2 = verts[edges[original_edge][1]];
               if ((point_in_segment(P1, v1, v2, eps)).is_set() && (point_in_segment(P2, v1, v2, eps)).is_set()) {
                  col = geom.get_e_col(original_edge);
                  cols.push_back(col);
               }
            }

            if (edge_blending == 'v')
               col = cols[0];
            else
               col = average_color(cols, color_system_mode, sat_power, sat_threshold, value_power, value_advance, alpha_mode, ryb_mode);
            geom.set_e_col(added_edge,col);
         }
      }

      // delete the original edges
      vector<int> deleted_edges;
      for(int i=0; i<original_edges_size; i++)
         deleted_edges.push_back(i);
      geom.delete_edges(deleted_edges);

      elems += 'e';
   }

   // for verts, the same processing as before merging can be done
   if (edge_blending == 'v' || edge_blending == 'b') {
      pre_edge_blend(geom, 'v',
                     color_system_mode, sat_power, sat_threshold, value_power, value_advance, alpha_mode, ryb_mode, eps);
      elems += 'v';
   }

   return elems;
}

void special_edge_process(col_geom_v &geom, const char &special_edge_processing, const int &color_system_mode,
                          const double &sat_power, const double &sat_threshold, const double &value_power, const double &value_advance, const int &alpha_mode, const bool &ryb_mode)
{
   const vector<vector<int> > &faces = geom.faces();
   const vector<vector<int> > &edges = geom.edges();
   const vector<vec3d> &verts = geom.verts();

   // except for hole connectors ...
   // strip visible edges and vertices
   vector<int> deleted_edges;
   for(unsigned int i=0; i<edges.size(); i++) {
      col_val col = geom.get_e_col(i);
      if (!(col.is_idx() && col == INT_MAX))
         deleted_edges.push_back(i);
   }
   geom.delete_edges(deleted_edges);

   for(unsigned int i=0; i<verts.size(); i++) {
      col_val col = geom.get_v_col(i);
      if (!col.is_inv())
         geom.set_v_col(i,col_val());
   }

   // and strip any standed vertices
   geom.delete_verts(geom.get_info().get_free_verts());

   // if stripping only
   if (special_edge_processing == 's')
      return;

   geom.add_missing_impl_edges();

   // special_edge_processing == 'e' will do at least this
   for(unsigned int i=0; i<edges.size(); i++) {
      col_val col = geom.get_e_col(i);
      if (col.is_idx() && col == INT_MAX)
         continue;
      vector<int> face_idx = find_faces_with_edge(faces, edges[i]);
      vector<col_val> cols;
      for(unsigned int j=0; j<face_idx.size(); j++)
         cols.push_back(geom.get_f_col(face_idx[j]));
      col = average_color(cols, color_system_mode, sat_power, sat_threshold, value_power, value_advance, alpha_mode, ryb_mode);
      geom.set_e_col(i,col);
      cols.clear();
   }

   if (special_edge_processing == 'v' || special_edge_processing == 'V') {
      for(unsigned int i=0; i<verts.size(); i++) {
         col_val col;
         // invisible vertices were the new ones for the tiles. skip them
         if (special_edge_processing == 'v') {
            col = geom.get_v_col(i);
            if (col.is_inv())
               continue;
         }
         vector<int> edge_idx = find_edges_with_vertex(edges, i);
         vector<col_val> cols;
         for(unsigned int j=0; j<edge_idx.size(); j++) {
            col_val col = geom.get_e_col(edge_idx[j]);
            if (col.is_idx() && col == INT_MAX)
               continue;
            cols.push_back(col);
         }
         col = average_color(cols, color_system_mode, sat_power, sat_threshold, value_power, value_advance, alpha_mode, ryb_mode);
         geom.set_v_col(i,col);
         cols.clear();
      }
   }
}

void delete_invisible_faces(col_geom_v &geom, const bool &hole_detection)
{
   const vector<vector<int> > &faces = geom.faces();
   const vector<vector<int> > &edges = geom.edges();

   // if passed a geom with only verts, keep
   if (!faces.size() && !edges.size())
      return;

   vector<int> deleted_elems;
   for(unsigned int i=0; i<faces.size(); i++) {
      col_val col = geom.get_f_col(i);
      if (col.is_inv())
         deleted_elems.push_back(i);
   }
   geom.delete_faces(deleted_elems);
   deleted_elems.clear();

   if (hole_detection) {
      for(unsigned int i=0; i<edges.size(); i++) {
         col_val col = geom.get_e_col(i);
         if (!(col.is_idx() && col == INT_MAX))
            continue;      
         if (!find_faces_with_edge(faces,edges[i]).size())
            deleted_elems.push_back(i);
      }
      geom.delete_edges(deleted_elems);
      deleted_elems.clear();
   }
}

void make_hole_connectors_invisible(col_geom_v &geom)
{
   const vector<vector<int> > &edges = geom.edges();

   for(unsigned int i=0; i<edges.size(); i++) {
      col_val col = geom.get_e_col(i);
      if (col.is_idx() && col == INT_MAX)
         geom.set_e_col(i,col_val::invisible);
   }
}

void resolve_winding_number_indexes(col_geom_v &geom, const color_map_multi &map, const color_map_multi &map_negative)
{
   const vector<vector<int> > &faces = geom.faces();

   for(unsigned int i=0; i<faces.size(); i++) {
      col_val col = geom.get_f_col(i);
      if (col.is_idx()) {
         int c_idx = col.get_idx();
         if (c_idx > INT_MAX/2)
            c_idx -= INT_MAX;

         if (c_idx >= 0)
            geom.set_f_col(i,map.get_col(c_idx));
         else
            geom.set_f_col(i,map_negative.get_col(abs(c_idx)));
      }
   }
}

void do_cmy_mode(col_geom_v &geom, const bool &ryb_mode, const char &edge_blending)
{
   const vector<vector<int> > &faces = geom.faces();
   const vector<vector<int> > &edges = geom.edges();
   const vector<vec3d> &verts = geom.verts();

   for(unsigned int i=0; i<faces.size(); i++)
      geom.set_f_col(i,rgb_complement(geom.get_f_col(i), ryb_mode));
   if (edge_blending == 'e' || edge_blending == 'b')
      for(unsigned int i=0; i<edges.size(); i++)
         geom.set_e_col(i,rgb_complement(geom.get_e_col(i), ryb_mode));
   if (edge_blending == 'v' || edge_blending == 'b')
      for(unsigned int i=0; i<verts.size(); i++)
         geom.set_v_col(i,rgb_complement(geom.get_v_col(i), ryb_mode));
}

void apply_transparency(col_geom_v &geom, int face_opacity)
{
   for (unsigned int i=0;i<geom.faces().size();i++) {
      col_val col = geom.get_f_col(i);
      if (col.is_val() && !col.is_inv())
         col = col_val(col[0],col[1],col[2],face_opacity);
      geom.set_f_col(i,col);
   }
}

void color_by_winding_number_raw(col_geom_v &geom, const char &color_by_winding_number, const double &eps)
{
   const vector<vector<int> > &faces = geom.faces();
   
   // get signed winding numbers
   for(unsigned int i=0; i<faces.size(); i++) {
      int wtotal = find_polygon_denominator_signed(geom, i, eps);
      
      int fsz = (int)faces[i].size();
      if (wtotal > fsz/2)
         wtotal -= fsz;
      
      // if absolute value of winding numbers (or their negative) are colored
      if (color_by_winding_number != 'w') {
         wtotal = abs(wtotal);
         if (color_by_winding_number == 'n')
            wtotal = -wtotal;
      }

      col_val col;
      if (wtotal >= 0)
         col.set_idx(wtotal);
      else
         // negative winding numbers are set up to near INT_MAX to be subtracted out later
         col.set_idx(wtotal + INT_MAX);
      geom.set_f_col(i,col);
   }
}

int main(int argc, char *argv[])
{
   planar_opts opts;
   opts.process_command_line(argc, argv);

   col_geom_v geom;
   geom_read_or_error(geom, opts.ifile, opts);

   char errmsg[MSG_SZ];
   if (opts.orient) {
      geom.orient(opts.orient, errmsg);
      if (*errmsg)
         opts.warning(errmsg, 'O');
   }

   // default center is centroid
   if (!opts.center.is_set())
      opts.center = centroid(geom.verts());

   // collect free edges get deleted later. save them if they exist
   // color merge
   // not much control over this feature
   col_geom_v egeom = free_edges_into_geom(geom);
   if (egeom.verts().size()) {
      string elems = pre_edge_blend(egeom, 'b',
                                    opts.color_system_mode, opts.sat_power, opts.sat_threshold, opts.value_power, opts.value_advance, opts.alpha_mode, opts.ryb_mode, opts.epsilon);
      sort_merge_elems(egeom, elems, 1, opts.epsilon);
   }

   if (opts.cmy_mode)
      do_cmy_mode(geom, opts.ryb_mode, opts.edge_blending);

   // blend edges using blending parameters without doing merging
   if (opts.edge_blending && !opts.planar_merge_type) {
      string elems = pre_edge_blend(geom, opts.edge_blending,
                                    opts.color_system_mode, opts.sat_power, opts.sat_threshold, opts.value_power, opts.value_advance, opts.alpha_mode, opts.ryb_mode, opts.epsilon);
      sort_merge_elems(geom, elems, 1, opts.epsilon);
   }

   if (opts.face_color_method)
      color_by_plane(geom, opts.face_color_method, opts.map, opts.center, opts.epsilon);

   int original_edges_size = 0;
   if (opts.planar_merge_type) {
      // missing edges are added so that they won't blend to invisible
      geom.add_missing_impl_edges();
      original_edges_size = geom.edges().size();

      planar_merge(geom, opts.planar_merge_type, opts.polygon_fill_type, opts.hole_detection, opts.center,
                   opts.winding_rule_mode, opts.winding_rule, opts.color_by_winding_number, opts.find_direction, opts.winding_div2, opts.verbose,
                   opts.zero_density_color, opts.zero_density_force_blend, opts.brightness_adj,
                   opts.color_system_mode, opts.sat_power, opts.sat_threshold, opts.value_power, opts.value_advance, opts.alpha_mode, opts.ryb_mode, opts.epsilon);
   }

   string elems;
   if (opts.edge_blending && opts.planar_merge_type) {
      elems = post_edge_blend(geom, opts.edge_blending, original_edges_size, opts.color_system_mode,
                              opts.sat_power, opts.sat_threshold, opts.value_power, opts.value_advance, opts.alpha_mode, opts.ryb_mode, opts.epsilon);
      // add free edges back
      geom.append(egeom);
   }

   // delay sort_merge so multiple edgelets still exist for post_edge_blend
   if (opts.planar_merge_type && elems!="")
      sort_merge_elems(geom, elems, 1, opts.epsilon);

   // resolve indexes of winding numbers to positive and negative color maps
   if (opts.color_by_winding_number) {
      if (!opts.planar_merge_type)
         color_by_winding_number_raw(geom, opts.color_by_winding_number, opts.epsilon);
      resolve_winding_number_indexes(geom, opts.map, opts.map_negative);
   }

   // add vertices to 'stitch' faces with dangling edges due to tiling or merging
   if (opts.stitch_faces)
      stitch_faces_on_seams(geom, opts.epsilon);

   if (opts.special_edge_processing) {
      special_edge_process(geom, opts.special_edge_processing, opts.color_system_mode,
                           opts.sat_power, opts.sat_threshold, opts.value_power, opts.value_advance, opts.alpha_mode, opts.ryb_mode);
      // add free edges back
      geom.append(egeom);
   }

   if (opts.delete_invisible_faces)
      delete_invisible_faces(geom, opts.hole_detection);

   // remove extra in line vertices from faces
   if (opts.simplify_face_edges)
      remove_in_line_vertices(geom, opts.epsilon);

   if (opts.hole_detection)
      make_hole_connectors_invisible(geom);

   // transparency
   if (opts.face_opacity != 255)
      apply_transparency(geom, opts.face_opacity);

   geom_write_or_error(geom, opts.ofile, opts);

   return 0;
}
