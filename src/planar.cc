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
#include <unistd.h>

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
      bool color_by_winding_number;
      bool orient;
      bool hole_detection;
      bool use_unit_vectors;
      bool stitch_faces;
      bool delete_invisible_faces;
      char blend_edges;
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

      int sig_compare;
      double epsilon;

      planar_opts(): prog_opts("planar"),
                        face_color_method('\0'),
                        planar_merge_type(0),
                        polygon_fill_type(0),
                        winding_rule(0),
                        color_by_winding_number(false),
                        orient(false),
                        hole_detection(false),
                        use_unit_vectors(false),
                        stitch_faces(false),
                        delete_invisible_faces(false),
                        blend_edges('\0'),
                        zero_density_color(col_val::invisible),
                        zero_density_force_blend(false),
                        brightness_adj(-2.0),
                        color_system_mode(2),
                        cmy_mode(false),
                        ryb_mode(false),
                        sat_power(0.0),
                        value_power(0.0),
                        sat_threshold(1.0),
                        value_advance(0.0),
                        alpha_mode(3),
                        face_opacity(255),
                        sig_compare(INT_MAX),
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
"               tile=1  merge=2 (default: none)\n"
"  -p <opt>  polygon fill algorithm.  angular=1  modulo2=2 (pnpoly)\n"
"               triangulation=3  even_overlap=4  alt_modulo2=5 (default: 1)\n"
"  -w <opt>  winding rule, include face parts according to winding number\n"
"               odd=1  even=2  positive=3  negative=4  nonzero=5 (default: none)\n"
"               zero=6  zodd=7  zpositive=8  znegative=9 (7,8,9 include zero)\n"
"  -O        orient the faces (if possible), flip orientation if oriented\n"
"  -H        hole detection. attempt to place connectors to holes\n"
"  -U <opt>  unit vectors. (default: on for hole detection or -f, off otherwise)\n"
"               off=0  on=1\n"
"  -S        stitch seams created by tiling or merging\n"
"  -E <opt>  remove explicit edges and blend new ones (sets -S)\n"
"               e - edges only  v - also blend vertices (default: none)\n"
"               V - also blend invisible vertices  s - strip edges and vertices\n"
"  -D        delete invisible faces created by modulo\n"
"  -l <lim>  minimum distance for unique vertex locations as negative exponent\n"
"               (default: %d giving %.0e)\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\nColor Blending Options (for option -d)\n"
"  -M <mode> color blending mode. HSV=1  HSL=2  RGB=3 (default: 2)\n"
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
"  -a <int>  alpha to use for blend. average=1  minimum=2  maximum=3 (default: 3)\n"
"  -y        RYB mode. Blend colors as in Red-Yellow-Blue color wheel\n"
"  -c        CMY mode. Complementary colors.  RGB->(RYB/GMO)->CMY->blend\n"
"  -b <val>  brightness adjustment for areas of blended colors\n"
"               valid values -1.0 to +1.0 (default: no adjustment)\n"
"               negative for darker, positive for lighter\n"
"               at 0, an area with 2 blended colors will not change\n"
"               but areas with 3 or more will become slightly darker\n"
"\nColoring Options (run 'off_util -H color' for help on color formats)\n"
"  -f <opt>  face color option (processed before -d)\n"
"               d - unique color for faces with opposite normals\n"
"               f - unique color for faces on same and opposite planes\n"
"               o - unique color for faces on same planes only\n"
"  -T <tran> face opacity for -f options. valid range from 0 to 255\n"
"  -m <maps> color maps for faces to be tried in turn (default: compound)\n"
"  -Z <col>  color for areas found colorless by modulo (default: invisible)\n"
"               key word: b - force a color blend\n"
"  -W        color by winding number. Overrides all other color options\n"
"  -n <maps> maps for negative winding numbers (default: rng6_S0V0.5:0)\n"
"\n",prog_name(), help_ver_text, int(-log(::epsilon)/log(10) + 0.5), ::epsilon);
}

void planar_opts::process_command_line(int argc, char **argv)
{
   opterr = 0;
   char c;
   char errmsg[MSG_SZ];

   string id;
   string map_file;
   string map_file_negative;
   int use_unit_vectors_int = -1;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hl:d:p:w:OHU:SE:Db:M:s:t:v:u:a:cyf:T:m:Z:Wn:o:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
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

         case 'w':
            id = get_arg_id(optarg, "odd=1|even=2|positive=3|negative=4|nonzero=5|zero=6|zodd=7|zpositive=8|znegative=9", argmatch_add_id_maps, errmsg);
            if(id=="")
               error(errmsg);
            winding_rule = atoi(id.c_str());
            break;

         case 'O':
            orient = true;
            break;

         case 'H':
            hole_detection = true;
            break;

         case 'U':
            id = get_arg_id(optarg, "off=0|on=1", argmatch_add_id_maps, errmsg);
            if(id=="")
               error(errmsg);
            use_unit_vectors_int = atoi(id.c_str());
            use_unit_vectors = (use_unit_vectors_int) ? true : false;
            break;

         case 'S':
            stitch_faces = true;
            break;

         case 'E':
            if(strspn(optarg, "evVs") != strlen(optarg) || strlen(optarg)>1)
               error(msg_str("edge blending is '%s', must be e, v, V or s", optarg), c);
            blend_edges = *optarg;

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
            if(!strlen(optarg)==1 || !strchr("dfo", *optarg))
               error("color method must be d, f or o");
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
            color_by_winding_number = true;
            break;

         case 'n':
            map_file_negative = optarg;
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

   if (!planar_merge_type) {
      if (polygon_fill_type)
         warning("polygon fill has no effect if tile or merge is not selected","p");
      if (winding_rule)
         warning("winding rule has no effect if tile of merge is not selected","w");
      if (hole_detection)
         warning("hole detection has no effect if tile or merge is not selected","H");
   }
   else {
      // set default polygon fill type here
      if (!polygon_fill_type)
         polygon_fill_type = 1;
   }

   if (use_unit_vectors_int == -1) {
      if (face_color_method || hole_detection)
         use_unit_vectors = true;
      else
         use_unit_vectors = false;
   }
//fprintf(stderr,"use_unit_vectors = %s\n",(use_unit_vectors ? "true" : "false"));

   if (!map_file.size())
      map_file = "compound";
   if(!map.init(map_file.c_str(), errmsg))
      error(errmsg, 'm');

   if (!map_file_negative.size())
      map_file_negative = "rng6_S0V0.5:0";
   if(!map_negative.init(map_file_negative.c_str(), errmsg))
      error(errmsg, 'n');

   epsilon = (sig_compare != INT_MAX) ? pow(10, -sig_compare) : ::epsilon;
}

void geom_dump(const col_geom_v &geom, string s)
{
   char errmsg[MSG_SZ];
   char filename[MSG_SZ];
   sprintf(filename,"geom_dump_%s.off",s.c_str());
   geom.write(filename, errmsg);
}

bool cmp_verts(const pair<vec3d, int> &a, const pair<vec3d, int> &b, double eps)
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

void build_coplanar_faces_list(const geom_if &geom, vector<vector<int> > &coplanar_faces_list, vector<xnormal> &coplanar_normals, const fnormals &fnormals,
                               bool point_outward, bool fold_normals, bool use_unit_vectors, double eps)
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
      if (use_unit_vectors)
         normal = normal.unit();

      pair<vec3d, int> face_normal_pair;
      face_normal_pair.first = normal;
      face_normal_pair.second = i;

      if (hemi)
         hemispherical_table.push_back(face_normal_pair);
      else
         face_normal_table.push_back(face_normal_pair);
   }

   vector<int> coplanar_face;

   // fold in hemispherical normals. this is what associates them on the same plane
   int sz = hemispherical_table.size();
   if (sz > 1) {
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
      coplanar_face.push_back(hemispherical_table[0].second);
      // find faces on same plane and group them together
      for(int i=1; i<sz; i++) {
         if (compare(hemispherical_table[i].first,hemispherical_table[i-1].first,eps)) {
            coplanar_faces_list.push_back(coplanar_face);
            coplanar_normals.push_back(fnormals[hemispherical_table[i-1].second]);
            coplanar_face.clear();
         }
         coplanar_face.push_back(hemispherical_table[i].second);
      }
      coplanar_faces_list.push_back(coplanar_face);
      coplanar_normals.push_back(fnormals[hemispherical_table[sz-1].second]);
      coplanar_face.clear();
   }

   
   // if non-hemispherical normals are folded only with specific option
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
      coplanar_face.push_back(face_normal_table[0].second);
      // find faces on same plane and group them together
      for(int i=1; i<sz; i++) {
         if (compare(face_normal_table[i].first,face_normal_table[i-1].first,eps)) {
            coplanar_faces_list.push_back(coplanar_face);
            coplanar_normals.push_back(fnormals[face_normal_table[i-1].second]);
            coplanar_face.clear();
         }
         coplanar_face.push_back(face_normal_table[i].second);
      }
      coplanar_faces_list.push_back(coplanar_face);
      coplanar_normals.push_back(fnormals[face_normal_table[sz-1].second]);
      coplanar_face.clear();
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
int vertex_into_geom(col_geom_v &geom, const vec3d P, double eps)
{
   col_val vcol = col_val::invisible;
   // uncomment for testing
   //vcol = col_val(0.0,0.0,1.0);
   
   int v_idx = find_vertex_by_coordinate(geom, P, eps);
   if (v_idx == -1) {
      geom.add_col_vert(P, vcol);
      v_idx = geom.verts().size()-1;
   }

   return v_idx;
}

// if edge already exists, do not create another one and return false. return true if new edge created
bool edge_into_geom(col_geom_v &geom, const int v_idx1, const int v_idx2)
{
   col_val ecol = col_val::invisible;
   // uncomment for testing
   //ecol = col_val(1.0,0.0,0.0);
   
   vector<int> new_edge = make_edge(v_idx1, v_idx2);

   int answer = find_edge_in_edge_list(geom.edges(), new_edge);
   if (answer < 0)
      geom.add_col_edge(new_edge, ecol);

   return ((answer < 0) ? true : false);
}

// put faces numbers in face_idxs into fgeom
col_geom_v faces_to_geom(const col_geom_v &geom, const vector<int> &face_idxs)
{
   col_geom_v fgeom;
   fgeom.add_verts(geom.verts());
   for(unsigned int i=0; i<face_idxs.size(); i++) {
      unsigned int j = face_idxs[i];
      fgeom.add_col_face(geom.faces()[j],geom.get_f_col(j));
   }
   fgeom.delete_verts(fgeom.get_info().get_free_verts());
   return fgeom;
}

bool is_point_on_polygon_edges(const geom_if &polygon, const vec3d &P, double eps)
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

bool is_point_in_bbox_2D(const geom_if &geom, const vec3d &P, const int idx)
{
   bound_box bb(geom.verts());
   vec3d min = bb.get_min();
   vec3d max = bb.get_max();

   int idx1 = (idx+1)%3;
   int idx2 = (idx+2)%3;

   return !(P[idx1] < min[idx1] || P[idx1] > max[idx1] || P[idx2] < min[idx2] || P[idx2] > max[idx2]);
}

// furshished by Adrian Rossiter
int side_of_line(const vec3d &P, const vec3d &A, const vec3d &B, const int idx1, const int idx2, double eps)
{
   double a = B[idx1] - A[idx1];
   double b = B[idx2] - A[idx2];
   double c = P[idx1] - A[idx1];
   double d = P[idx2] - A[idx2];
   double det = a*d - b*c;

   if(det<-eps)
      return -1;
   else
   if(det>eps)
      return 1;
   else
      return 0;
}

// edge is between -eps and eps. '>-eps' will include edge
bool is_point_inside_triangle_or_edge(const vec3d &P, const vec3d &A, const vec3d &B, const vec3d &C,
                                      const int idx1, const int idx2, const int sign, double eps)                          
{
   return (side_of_line(P,A,B,idx1,idx2,eps)*sign>-eps && side_of_line(P,B,C,idx1,idx2,eps)*sign>-eps && side_of_line(P,C,A,idx1,idx2,eps)*sign>-eps);
}

// polygon sent to this function MUST be triangulated
// includes edges and vertices as inside (bounds are considered inside points)
bool InsidePolygon_triangulated(const geom_if &polygon, const vec3d &P, const int idx, const int sign, const bool include_edges, double eps)
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
      if (!is_point_in_bbox_2D(triangle, P, idx))
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
bool InsidePolygon_solution1(const geom_if &polygon, const vec3d &P, bool &bound, const int idx, double eps)
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
      if (double_equality(p_x, p1_x, eps) && double_equality(p_y, p1_y, eps)) {
         bound = true;
         return true;
      }

      //right vertex
      double p2_x = verts[face[i%N]][idx1];
      double p2_y = verts[face[i%N]][idx2];

      //ray is outside of our interests
      if (p_y < min(p1_y, p2_y) || p_y > max(p1_y, p2_y)) {
         //next ray left point
         p1_x = p2_x;
         p1_y = p2_y;
         continue;
      }

      //ray is crossing over by the algorithm (common part of)
      if (p_y > min(p1_y, p2_y) && p_y < max(p1_y, p2_y)) {

         //x is before of ray
         if (p_x <= max(p1_x, p2_x)+eps) {

            //overlies on a horizontal ray
            if (double_equality(p1_y, p2_y, eps) && p_x+eps >= min(p1_x, p2_x)) {
               bound = true;
               return true;
            }

            //ray is vertical
            if (double_equality(p1_x, p2_x, eps)) {

               //overlies on a ray
               if (double_equality(p1_x, p_x, eps)) {
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
               if (fabs(p_x - xinters) < eps) {
                 bound = true;
                 return true;
               }

               //before ray
               if (p_x < xinters)
                  ++__count;
            }
         }
      }

      //special case when ray is crossing through the vertex
      else
      {
         //p crossing over p2
         if (double_equality(p_y, p2_y, eps) && p_x < p2_x+eps) {

            //next vertex
            double p3_y = verts[face[(i+1)%N]][idx2];

            //p_y lies between p1_y & p3_y
            if ((p_y+eps) >= min(p1_y, p3_y) && p_y <= max(p1_y, p3_y)+eps) {
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
bool InsidePolygon_solution1(const geom_if &polygon, const vec3d &P, const int idx, const bool include_edges, double eps)
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

bool pnpoly(const geom_if &polygon, const vec3d &P, const int idx)
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

      if ( ((verty_i>testy) != (verty_j>testy)) &&
	        (testx < (vertx_j-vertx_i) * (testy-verty_i) / (verty_j-verty_i) + vertx_i) )
         c = !c;
   }

   return (c ? true : false);
}

// wrapper
bool pnpoly(const geom_if &polygon, const vec3d &P, int idx, bool include_edges, double eps)
{
   bool answer = pnpoly(polygon, P, idx);

   // have to check bounds conditions if we might have to negate the answer
   if ((answer && !include_edges) || (!answer && include_edges)) {
      if (is_point_on_polygon_edges(polygon, P, eps))
         answer = !answer;
   }

   return answer;
}

// Copyright 2001, softSurfer (www.softsurfer.com)
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.

//    a Point is defined by its coordinates {int x, y;}
//===================================================================

// isLeft(): tests if a point is Left|On|Right of an infinite line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 and P1
//            =0 for P2 on the line
//            <0 for P2 right of the line
//    See: the January 2001 Algorithm "Area of 2D and 3D Triangles and Polygons"

//inline int
//isLeft( Point P0, Point P1, Point P2 )

// RK - this works with wn_PnPoly better if a double is passed back

double isLeft(const vec3d &P0, const vec3d &P1, const vec3d &P2, const int idx)
{
   int idx1 = (idx+1)%3;
   int idx2 = (idx+2)%3;

   double P0_x = P0[idx1];
   double P0_y = P0[idx2];

   double P1_x = P1[idx1];
   double P1_y = P1[idx2];

   double P2_x = P2[idx1];
   double P2_y = P2[idx2];

   return ( (P1_x - P0_x) * (P2_y - P0_y)
          - (P2_x - P0_x) * (P1_y - P0_y) );
}

//===================================================================

// cn_PnPoly(): crossing number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  0 = outside, 1 = inside
// This code is patterned after [Franklin, 2000]

//int
//cn_PnPoly( Point P, Point* V, int n )

// RK - try to use same variable names as original pnpoly

/* RK - Not currently used
bool cn_PnPoly(const geom_if &polygon, const vec3d &P, const int idx, int &crossing_number, double eps)
{
   const vector<vec3d> &verts = polygon.verts();

   int idx1 = (idx+1)%3;
   int idx2 = (idx+2)%3;

   double testx = P[idx1]; // P.x
   double testy = P[idx2]; // P.y

   const vector<int> &face = polygon.faces()[0];
   int n = face.size();

   int cn = 0;    // the crossing number counter

   // loop through all edges of the polygon
   for (int i=0; i<n; i++) {    // edge from face[i] to face[i+1]
      int j = (i+1)%n;

      double vertx_i = verts[face[i]][idx1]; // V[i].x
      double verty_i = verts[face[i]][idx2]; // V[i].y
      double vertx_j = verts[face[j]][idx1]; // V[i+1].x
      double verty_j = verts[face[j]][idx2]; // V[i+1].y

      if (((verty_i <= testy) && (verty_j > testy))          // an upward crossing
       || ((verty_i > testy) && (verty_j <= testy))) {       // a downward crossing
         // compute the actual edge-ray intersect x-coordinate
         double vt = (testy - verty_i) / (verty_j - verty_i);
         if (testx < vertx_i + vt * (vertx_j - vertx_i))     // testx < intersect
            ++cn;   // a valid crossing of y=testy right of testx
       }
   }

   crossing_number = cn;

   //return (cn&1);    // 0 if even (out), and 1 if odd (in)
   return (cn&1 ? true : false); // 0 if even (out), and 1 if odd (in)
}
*/
//===================================================================

// wn_PnPoly(): winding number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  wn = the winding number (=0 only if P is outside V[])

//int
//wn_PnPoly( Point P, Point* V, int n )

// RK - try to use same variable names as original pnpoly
// RK - wn is n time too large. winding_number = wn/n

int wn_PnPoly(const geom_if &polygon, const vec3d &P, const int idx, int &winding_number, double eps)
{
   const vector<vec3d> &verts = polygon.verts();

   int idx2 = (idx+2)%3;

   double testy = P[idx2]; // P.y

   const vector<int> &face = polygon.faces()[0];
   int n = face.size();

   int wn = 0;    // the winding number counter

   // loop through all edges of the polygon
   for (int i=0; i<n; i++) {   // edge from face[i] to face[i+1]
      int j = (i+1)%n;

      vec3d vert_i = verts[face[i]]; // V[i]
      vec3d vert_j = verts[face[j]]; // V[[i+1]

      double verty_i = vert_i[idx2]; // V[i].y
      double verty_j = vert_j[idx2]; // V[i+1].y

      for (int i=0; i<n; i++) {      // edge from face[i] to face[i+1]
         if (verty_i <= testy) {     // start y <= P.y
            if (verty_j > testy)     // an upward crossing
               if (isLeft( vert_i, vert_j, P, idx ) > eps)  // P left of edge
                  ++wn;              // have a valid up intersect
         }
         else {                      // start y > P.y (no test needed)
            if (verty_j <= testy)    // a downward crossing
               if (isLeft( vert_i, vert_j, P, idx ) < -eps)  // P right of edge
                  --wn;              // have a valid down intersect
         }
      }
   }

   winding_number = wn / n;

   return (wn ? true : false);
}

//===================================================================

// ported from http://paulbourke.net/geometry/insidepoly/

// Return the angle between two vectors on a plane
// The angle is from vector 1 to vector 2, positive anticlockwise
// The result is between -pi -> pi

double Angle2D(const double &x1, const double &y1, const double &x2, const double &y2)
{
   double dtheta,theta1,theta2;

   theta1 = atan2(y1,x1);
   theta2 = atan2(y2,x2);
   dtheta = theta2 - theta1;
   while (dtheta > M_PI)
      dtheta -= M_PI * 2.0;
   while (dtheta < -M_PI)
      dtheta += M_PI * 2.0;

   return(dtheta);
}

// method 2 function ported from http://paulbourke.net/geometry/insidepoly/
// RK - does not do modulo-2 implicitly

// Another solution forwarded by Philippe Reverdy is to compute the sum of the angles made between the test point and each pair of points making up the polygon.
// If this sum is 2pi then the point is an interior point, if 0 then the point is an exterior point.
// This also works for polygons with holes given the polygon is defined with a path made up of coincident edges into and out of the hole as is common practice in many CAD packages.

bool InsidePolygon_solution2(const geom_if &polygon, const vec3d &P, const int idx)
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
      angle += Angle2D(p1_x,p1_y,p2_x,p2_y);
   }

   if (fabs(angle) < M_PI)
      return(false);
   else
      return(true);
}

// wrapper
bool InsidePolygon_solution2(const geom_if &polygon, const vec3d &P, const int idx, const bool include_edges, double eps)
{
   bool answer = InsidePolygon_solution2(polygon, P, idx);

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
bool is_point_inside_polygon(const col_geom_v &polygon, const vec3d &P, const vec3d &normal, const bool include_edges, const bool on_same_plane_check,
                             const int polygon_fill_type, double eps)
{  
   // quick check to see if P is actually on the polygon
   // funished by Adrian Rossiter; is there any distance from P to the unit normal?
   if (on_same_plane_check)
      if (!double_equality(vdot(polygon.verts()[0]-P, normal), 0.0, eps))
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

int intersection_is_end_point(const vec3d &intersection_point, const vec3d &P0, const vec3d &P1, double eps)
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
bool mesh_verts(col_geom_v &geom, double eps)
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
         edge_into_geom(geom, edges[i][0], line_intersections[0].second);
         for(unsigned int k=0; k<line_intersections.size()-1; k++)
            edge_into_geom(geom, line_intersections[k].second, line_intersections[k+1].second);
         edge_into_geom(geom, line_intersections[line_intersections.size()-1].second, edges[i][1]);
      }
   }
   
   // delete the replaced edges
   geom.delete_edges(deleted_edges);

   return (deleted_edges.size() ? true : false);
}

// of geom, is face i inside face j
bool is_face_inside_face(geom_if &geom, const int i, const int j, double eps)
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

void check_for_holes(geom_if &geom, vector<pair<int, int> > &polygon_hierarchy, double eps)
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

void make_connectors(geom_if &geom, vector<vector<int> > &connectors, vector<pair<int, int> > &polygon_hierarchy)
{
   const vector<vector<int> > &faces = geom.faces();
 
   int last_idx = -1;
   for(unsigned int i=0; i<polygon_hierarchy.size(); i++) {
      pair<int, int> hierarchy = polygon_hierarchy[i];
      // only allow child connection once
      if (hierarchy.first == last_idx)
         continue;
      int v_idx1 = faces[hierarchy.first][0];  // child
      int v_idx2 = faces[hierarchy.second][0]; // parent
      connectors.push_back(make_edge(v_idx1, v_idx2));
      last_idx = hierarchy.first;
   }
}

void add_connectors(col_geom_v &geom, vector<vector<int> > &connectors)
{
   coloring clrng(&geom);
   
   col_val ecol = col_val::invisible;
   // uncomment for testing
   //ecol = col_val(1.0,1.0,1.0);
   
   for(unsigned int i=0; i<connectors.size(); i++) {
      geom.add_col_edge(connectors[i], ecol);
   }
}

// input seperate networks of overlapping edges and merge them into one network
bool mesh_edges(col_geom_v &geom, double eps)
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
               v_idx = vertex_into_geom(geom, intersection_point, eps);
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
         edge_into_geom(geom, edges[i][0], line_intersections[0].second);
         for(unsigned int k=0; k<line_intersections.size()-1; k++)
            edge_into_geom(geom, line_intersections[k].second, line_intersections[k+1].second);
         edge_into_geom(geom, line_intersections[line_intersections.size()-1].second, edges[i][1]);
      }
   }
   
   // delete the replaced edges
   geom.delete_edges(deleted_edges);

   return (deleted_edges.size() ? true : false);
}

void make_skeleton(col_geom_v &geom)
{
   geom.add_missing_impl_edges();
   geom.clear_faces();

   coloring clrng(&geom);
   
   col_val ecol = col_val::invisible;
   // uncomment for testing
   //ecol = col_val(1.0,1.0,1.0);
   
   clrng.e_one_col(ecol);
   clrng.v_one_col(ecol);
}

// use faces to find vertices instead of edge numbers
// count it as an intersection it is an edge end point
bool check_if_faces_intersect(const col_geom_v &geom, const int face_idx1, const int face_idx2, double eps)
{
   const vector<vector<int> > &faces = geom.faces();
   const vector<vec3d> &verts = geom.verts();
   
   vector<int> face1 = faces[face_idx1];
   vector<int> face2 = faces[face_idx2];
   
   int sz1 = face1.size();
   int sz2 = face2.size();
   for(int i=0; i<sz1; i++) {
      vec3d P0 = verts[face1[i]];
      vec3d P1 = verts[face1[(i+1)%sz1]];
      for(int j=0; j<sz2; j++) {
         vec3d Q0 = verts[face2[j]];
         vec3d Q1 = verts[face2[(j+1)%sz2]];
         
         vec3d intersection_point = lines_intersection_in_segments(P0, P1, Q0, Q1, eps);
         if (intersection_point.is_set())
            return true;
      }
   }
  
   return false;
}

// of a set of coplanar faces, find those that actually overlap and group them
void find_coplanar_clusters(const col_geom_v &geom, vector<vector<int> > &coplanar_clusters, const vector<int> &coplanar_face_list, bool force, double eps)
{
   vector<pair<int, int> > overlap_pairs;
  
   for(unsigned int i=0; i<coplanar_face_list.size()-1; i++) {
      for(unsigned int j=i+1; j<coplanar_face_list.size(); j++) {
         // if triangulating, there will always be overlaps 
         if (force || check_if_faces_intersect(geom, coplanar_face_list[i], coplanar_face_list[j], eps)) {
            overlap_pairs.push_back(make_pair(coplanar_face_list[i], coplanar_face_list[j]));
         }
      }
   }

   vector<int> already_listed;
   for(unsigned int i=0; i<coplanar_face_list.size(); i++) {
      if(find(already_listed.begin(), already_listed.end(), coplanar_face_list[i]) != already_listed.end())
         continue;
      vector<int> cluster;
      cluster.push_back(coplanar_face_list[i]);
      for(unsigned int j=0; j<cluster.size(); j++) {
         for(unsigned int k=0; k<overlap_pairs.size(); k++) {
            int l = -1;
            if (overlap_pairs[k].first == cluster[j])
               l = overlap_pairs[k].second;
            else
            if (overlap_pairs[k].second == cluster[j])
               l = overlap_pairs[k].first;
         
            if (l>-1) {
               if(find(cluster.begin(), cluster.end(), l) == cluster.end()) {
                  cluster.push_back(l);
                  already_listed.push_back(l);
               }
            }
         }
      }
      coplanar_clusters.push_back(cluster);
      cluster.clear();
   }
}

// find connections from a vertex v_idx
void find_connections(const geom_if &geom, vector<int> &vcons, const int v_idx)
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

void build_angle_map(const geom_if &geom, map<pair<int, int>, double> &angle_map, const vec3d normal, double eps)
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
         if (double_equality(angle, 0.0, eps))
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
            // can't go backwards
            if (c == a)
               continue;
            double angle = angle_map[make_pair(b,c)];;
            if (angle >= 0 && angle < base_angle)
               angle +=360;
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
         geom.add_face(face);
         // mark turns visited and start a new face
         fsz = face.size();
         for(int i=0; i<fsz; i++)
            visited[make_pair(face[i],face[(i+1)%fsz])] = true;
         get_first_unvisited_triad(geom,first_unvisited,visited,turn_map,face);
      }
   }
}

void analyze_faces(geom_if &geom, int planar_merge_type, vector<int> &nonconvex_faces, map<pair<int, int>, double> &angle_map)
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
      if ((sum == 360 && planar_merge_type == 1) || (sum != 360 && planar_merge_type == 2)) {
         deleted_faces.push_back(i);
         deleted_faces_count++;
      }
      else
      if (!all_negative_turns)
         nonconvex_faces.push_back(i-deleted_faces_count);
   }
   geom.delete_faces(deleted_faces);
}

void fill_in_faces(col_geom_v &geom, int planar_merge_type, vector<int> &nonconvex_faces, const vec3d normal, double eps)
{
   map<pair<int, int>, double> angle_map;
   map<pair<int, int>, int> turn_map;
   build_angle_map(geom, angle_map, normal, eps);
   build_turn_map(geom, turn_map, angle_map);
   construct_faces(geom, turn_map);
   analyze_faces(geom, planar_merge_type, nonconvex_faces, angle_map);
}

col_val average_color(const vector<col_val> &cols, int color_system_mode, double sat_power, double sat_threshold, double value_power, double value_advance, int alpha_mode, bool ryb_mode)
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

// odd=1  even=2  positive=3  negative=4  nonzero=5
// zero=6  zodd=7  zpositive=8  znegative=9 (7,8,9 include zero)
bool winding_rule_filter(int winding_rule, int winding_number)
{
   bool answer = true;

   if (winding_rule < 6 && !winding_number)
      answer = false;
   else
   if ((winding_rule == 1) && (is_even(winding_number)))
      answer = false;
   else
   if ((winding_rule == 2) && (!is_even(winding_number)))
      answer = false;
   else
   if ((winding_rule == 3 || winding_rule == 8) && (winding_number < 0))
      answer = false;
   else
   if ((winding_rule == 4 || winding_rule == 9) && (winding_number > 0))
      answer = false;
   else
   if ((winding_rule == 5) && (!winding_number))
      answer = false;
   else
   if ((winding_rule == 6) && (winding_number))
      answer = false;
   else
   if ((winding_rule == 7) && (winding_number && is_even(winding_number)))
      answer = false;

   return answer;
}

// the winding number is correct when the faces are placed on the xy-plane
// if 2D, this decision cannot be made
// assumes point is already a hit so answer is not checked
// geom and P are copies because they are going to be rotated
int get_winding_number(col_geom_v geom, vec3d P, xnormal &pnormal, double eps)
{
   if (!pnormal.is_hemispherical()) {
      mat3d trans = mat3d::rot(pnormal.outward().unit(), vec3d(0,0,1));
      geom.transform(trans);

      col_geom_v vgeom;
      int v_idx = vertex_into_geom(vgeom, P, eps);
      vgeom.transform(trans);
      P = vgeom.verts()[v_idx];
      vgeom.clear_all();
   }

   int winding_number = 0;
   // idx = 2;
   wn_PnPoly(geom, P, 2, winding_number, eps);
   return winding_number;
}

void sample_colors(col_geom_v &sgeom, const col_geom_v &cgeom, int planar_merge_type, vector<xnormal> &original_normals, const vector<int> &nonconvex_faces,
                   int polygon_fill_type, int winding_rule, bool color_by_winding_number, col_val zero_density_color, bool zero_density_force_blend, double brightness_adj,
                   int color_system_mode, double sat_power, double sat_threshold, double value_power, double value_advance, int alpha_mode, bool ryb_mode, double eps)
{
   const vector<vector<int> > &sfaces = sgeom.faces();
   const vector<vector<int> > &cfaces = cgeom.faces();

   vector<col_val> cols;
   if (zero_density_force_blend) {
      for(unsigned int i=0; i<cfaces.size(); i++)
         cols.push_back(cgeom.get_f_col(i));
         zero_density_color = average_color(cols, color_system_mode, sat_power, sat_threshold, value_power, value_advance, alpha_mode, ryb_mode);
   }
   cols.clear();

   for(unsigned int i=0; i<sfaces.size(); i++) {
      vector<vec3d> points;

      // if the sampling face is convex
      if(find(nonconvex_faces.begin(), nonconvex_faces.end(), (int)i) == nonconvex_faces.end())
         points.push_back(centroid(sgeom.verts(), sfaces[i]));
      else {
      // else non-convex. need to triangulate the non-convex sample polygon to sample on centroid(s) of the triangles (until one is hit)
         vector<int> sface_idxs;
         sface_idxs.push_back(i);
         col_geom_v spolygon = faces_to_geom(sgeom, sface_idxs);
         spolygon.triangulate();
         const vector<vector<int> > &spfaces = spolygon.faces();
         // when merging and it is a nonconvex face, then have to sample all the centroids. else only sample one of them
         int sz = (planar_merge_type==1) ? 1 : spfaces.size();
         for(int k=0; k<sz; k++)
            points.push_back(centroid(spolygon.verts(), spfaces[k]));
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
        vec3d normal = original_normal.raw().unit();

         // if merging we have to sample all the centroids until there is a hit. otherwise k will begin and end at 0
         // checking to see if points landing on edges is not necessary since centroids should always be away from edges
         for(unsigned int k=0; k<points.size(); k++) {
            bool answer = is_point_inside_polygon((polygon_fill_type == 3 ? tpolygon : polygon), points[k], normal, false, false, polygon_fill_type, eps);
            if (answer) {
               if (winding_rule || color_by_winding_number) {
                  int winding_number = get_winding_number(polygon, points[k], original_normal, eps);
                  winding_total += winding_number;
               }

               cols.push_back(cgeom.get_f_col(j));
               break;
            }
         }
      }

      if (!color_by_winding_number) { 
         if (cols.size() && winding_rule)
            if (!winding_rule_filter(winding_rule, winding_total))
               cols.clear();
      }

      int sz = cols.size();

      col_val c;

      if (color_by_winding_number) {
         if (winding_total >= 0)
            c.set_idx(winding_total);
         else
            // negative winding numbers are set up to near INT_MAX to be subtracted out later
            c.set_idx(winding_total + INT_MAX);
      }
      else
      // if there is no hit, then that patch is invisible. if file_type 4 then all even numbered patches are made invisible
      if (!sz || (polygon_fill_type == 4 && !(sz%2)))
         c = col_val(zero_density_color);
      else
         c = average_color(cols, color_system_mode, sat_power, sat_threshold, value_power, value_advance, alpha_mode, ryb_mode);

      if ((brightness_adj > -2.0) && c.is_set() && !c.is_inv() && (sz > 1)) {
         // add 0.5 to brightness so it is in the range of -0.5 to 1.5
         // if sz = 2, 1.0/2 = 0.5 so when -B is 0, a patch with 2 colors blended will not be changed 
         double brightness = (1.0/sz)+(brightness_adj + 0.5);
         c.set_brightness(-1.0+brightness);
      }

      sgeom.set_f_col(i,c);
      cols.clear();
   }
}

void collect_original_normals(vector<xnormal> &original_normals, const vector<int> &coplanar_cluster, const fnormals &fnormals)
{
   for(unsigned int i=0; i<coplanar_cluster.size(); i++)
      original_normals.push_back(fnormals[coplanar_cluster[i]]);
}

void blend_overlapping_colors(col_geom_v &geom, const vector<vector<int> > &coplanar_faces_list, vector<xnormal> &coplanar_normals, fnormals &fnormals,
                              int planar_merge_type, int polygon_fill_type,
                              bool hole_detection, int winding_rule, bool color_by_winding_number, col_val zero_density_color, bool zero_density_force_blend, double brightness_adj,
                              int color_system_mode, double sat_power, double sat_threshold, double value_power, double value_advance, int alpha_mode, bool ryb_mode, double eps)
{
   col_geom_v bgeom;
   vector<int> deleted_faces;
   
   for(unsigned int i=0; i<coplanar_faces_list.size(); i++) {
      vector<vector<int> > coplanar_clusters;
      find_coplanar_clusters(geom, coplanar_clusters, coplanar_faces_list[i], (hole_detection||winding_rule||color_by_winding_number), eps);
      for(unsigned int j=0; j<coplanar_clusters.size(); j++) {
         // load a geom with color faces. keep it and copy it.
         col_geom_v cgeom = faces_to_geom(geom, coplanar_clusters[j]);
         col_geom_v sgeom = cgeom;

         // check here for polygons within polygons
         vector<vector<int> > connectors;
         if (hole_detection) {
            vector<pair<int, int> > polygon_hierarchy;
            check_for_holes(sgeom, polygon_hierarchy, eps);
            if (polygon_hierarchy.size())
               make_connectors(sgeom, connectors, polygon_hierarchy);
         }

         make_skeleton(sgeom);

         if (connectors.size())
            add_connectors(sgeom, connectors);

         // duplicate vertices and edges can cause problems
         sort_merge_elems(sgeom, "ve", eps);
         bool meshed_verts = mesh_verts(sgeom, eps);
         bool meshed_edges = mesh_edges(sgeom, eps);

         // if nothing meshed then there is nothing to do
         if (!meshed_verts && !meshed_edges && !connectors.size() && !winding_rule && !color_by_winding_number)
            continue;

         vector<int> nonconvex_faces;
         fill_in_faces(sgeom, planar_merge_type, nonconvex_faces, coplanar_normals[i].outward().unit(), eps);
     
         // original normals are needed for sampling colors
         vector<xnormal> original_normals;
         collect_original_normals(original_normals, coplanar_clusters[j], fnormals);
         sample_colors(sgeom, cgeom, planar_merge_type, original_normals, nonconvex_faces,
                       polygon_fill_type, winding_rule, color_by_winding_number, zero_density_color, zero_density_force_blend, brightness_adj,
                       color_system_mode, sat_power, sat_threshold, value_power, value_advance, alpha_mode, ryb_mode, eps);
    
         bgeom.append(sgeom);

         // mark the faces in the cluster for deletion at the end
         deleted_faces.insert(deleted_faces.end(), coplanar_clusters[j].begin(), coplanar_clusters[j].end());
      }
   }

   geom.delete_faces(deleted_faces);
   geom.append(bgeom);

   // there are lots of extra vertices and edges
   sort_merge_elems(geom, "ve", eps);
}

void planar_merge(col_geom_v &geom, int planar_merge_type, int polygon_fill_type, bool hole_detection, bool use_unit_vectors,
                  int winding_rule, bool color_by_winding_number, col_val zero_density_color, bool zero_density_force_blend, double brightness_adj,
                  int color_system_mode, double sat_power, double sat_threshold, double value_power, double value_advance, int alpha_mode, bool ryb_mode, double eps)
{
   fnormals face_normals(geom, vec3d(), eps);

   vector<vector<int> > coplanar_faces_list;
   vector<xnormal> coplanar_normals;
   
   bool point_outward = true;
   bool fold_normals = false;
   build_coplanar_faces_list(geom, coplanar_faces_list, coplanar_normals, face_normals, point_outward, fold_normals, use_unit_vectors, eps);

   // missing edges are added so that they won't blend to invisible
   geom.add_missing_impl_edges();
   blend_overlapping_colors(geom, coplanar_faces_list, coplanar_normals, face_normals,
                            planar_merge_type, polygon_fill_type,
                            hole_detection, winding_rule, color_by_winding_number, zero_density_color, zero_density_force_blend, brightness_adj,
                            color_system_mode, sat_power, sat_threshold, value_power, value_advance, alpha_mode, ryb_mode, eps);
}

void color_by_plane(col_geom_v &geom, char face_color_method, int face_opacity, color_map_multi &map, int use_unit_vectors, double eps)
{
   fnormals fnormals(geom, vec3d(), eps);

   vector<vector<int> > coplanar_faces_list;
   vector<xnormal> coplanar_normals;
   
   // pointed outward unless the d,D options was chosen
   bool point_outward = (face_color_method == 'd' || face_color_method == 'D') ? false : true;
   bool fold_normals = (face_color_method == 'f' || face_color_method == 'F') ? true : false;
   build_coplanar_faces_list(geom, coplanar_faces_list, coplanar_normals, fnormals, point_outward, fold_normals, use_unit_vectors, eps);
   
   for(unsigned int i=0; i<coplanar_faces_list.size(); i++) {
      for(unsigned int j=0; j<coplanar_faces_list[i].size(); j++) {
         int k = coplanar_faces_list[i][j];
         col_val col = map.get_col(i);
         geom.set_f_col(k,col);
      }
   }

   // transparency
   if (face_opacity != 255) {
      for (unsigned int i=0;i<geom.faces().size();i++) {
         col_val col = geom.get_f_col(i);
         if (col.is_val())
            col = col_val(col[0],col[1],col[2],face_opacity);
         geom.set_f_col(i,col);
      }
   }
}

void build_colinear_edge_list(const geom_if &geom, const vector<vector<int> > &edges, vector<vector<int> > &colinear_edge_list, double eps)
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

int find_end_point_vertex(const geom_if &geom, vector<int> vert_indexes, double eps)
{
   const vector<vec3d> &verts = geom.verts();

   col_geom_v vgeom;
   for(unsigned int i=0;i<vert_indexes.size();i++)
      vertex_into_geom(vgeom, verts[vert_indexes[i]], eps);
   vgeom.set_hull();

   const vector<vec3d> &gverts = vgeom.verts();
   vec3d end_point = gverts[0];
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
      end_point = gverts[distance_table[sz-1].second];
   }
   vgeom.clear_all();

   int end_idx = -1;
   for(unsigned int i=0;i<vert_indexes.size();i++) {
      if (!compare(end_point, verts[vert_indexes[i]], eps)) {
         end_idx = vert_indexes[i];
         break;
      }
   }

   return end_idx;   
}

void collect_ordered_vert_indexes(const geom_if &geom, const vector<vector<int> > &edges, const vector<int> &colinear_edges, vector<int> &colinear_verts, double eps)
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

   int end_idx = find_end_point_vertex(geom, vert_indexes, eps);
   vec3d end_vertex = verts[end_idx];

   vector<pair<double, int> > distance_table;
   distance_table.push_back(make_pair(0.0, end_idx));
   for(unsigned int i=0;i<vert_indexes.size();i++) {
      if (vert_indexes[i] == end_idx)
         continue;;
      double dist = (verts[vert_indexes[i]]-end_vertex).mag();
      distance_table.push_back(make_pair(dist,vert_indexes[i]));
   }
   vert_indexes.clear();

   sort( distance_table.begin(), distance_table.end() );

   for(unsigned int i=0;i<distance_table.size();i++)
      colinear_verts.push_back(distance_table[i].second);
}

void build_colinear_vertex_list(const geom_if &geom, vector<vector<int> > &colinear_vertex_list, double eps)
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
void insert_verts_into_face(geom_if &geom, vector<int> &added_vertices, int face_idx, int start_v_idx, int end_v_idx)
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
void stitch_faces_on_seams(col_geom_v &geom, double eps)
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

void blend_edges(col_geom_v &geom, char blend_edges,
                 int color_system_mode, double sat_power, double sat_threshold, double value_power, double value_advance, int alpha_mode, bool ryb_mode)
{
   const vector<vector<int> > &faces = geom.faces();
   const vector<vector<int> > &edges = geom.edges();
   const vector<vec3d> &verts = geom.verts();

   // clear old edges
   geom.clear_edges();

   // if stripping only, also strip vertex colors that are not invisible
   if (blend_edges == 's') {
      for(unsigned int i=0; i<verts.size(); i++) {
         col_val c = geom.get_v_col(i);
         if (!c.is_inv())
            geom.set_v_col(i,col_val());
      }
      return;
   }

   geom.add_missing_impl_edges();

   for(unsigned int i=0; i<edges.size(); i++) {
      vector<int> face_idx = find_faces_with_edge(faces, edges[i]);
      vector<col_val> cols;
      for(unsigned int j=0; j<face_idx.size(); j++)
         cols.push_back(geom.get_f_col(face_idx[j]));
      col_val c = average_color(cols, color_system_mode, sat_power, sat_threshold, value_power, value_advance, alpha_mode, ryb_mode);
      geom.set_e_col(i,c);
      cols.clear();
   }

   if (blend_edges == 'v' || blend_edges == 'V') {
      for(unsigned int i=0; i<verts.size(); i++) {
         col_val c;
         // invisible vertices were the new ones for the tiles. skip them
         if (blend_edges == 'v') {
            c = geom.get_v_col(i);
            if (c.is_inv())
               continue;
         }
         vector<int> edge_idx = find_edges_with_vertex(edges, i);
         vector<col_val> cols;
         for(unsigned int j=0; j<edge_idx.size(); j++)
            cols.push_back(geom.get_e_col(edge_idx[j]));
         c = average_color(cols, color_system_mode, sat_power, sat_threshold, value_power, value_advance, alpha_mode, ryb_mode);
         geom.set_v_col(i,c);
         cols.clear();
      }
   }
}

void resolve_winding_number_indexes(col_geom_v &geom, color_map_multi &map, color_map_multi &map_negative)
{
   const vector<vector<int> > &faces = geom.faces();

   for(unsigned int i=0; i<faces.size(); i++) {
      col_val c = geom.get_f_col(i);
      if (c.is_idx()) {
         int c_idx = c.get_idx();
         if (c_idx > INT_MAX/2)
            c_idx -= INT_MAX;

         if (c_idx >= 0)
            geom.set_f_col(i,map.get_col(c_idx));
         else
            geom.set_f_col(i,map_negative.get_col(abs(c_idx)));
      }
   }
}

void delete_invisible_faces(col_geom_v &geom)
{
   const vector<vector<int> > &faces = geom.faces();
   vector<int> deleted_faces;

   for(unsigned int i=0; i<faces.size(); i++) {
      col_val c = geom.get_f_col(i);
      if (c.is_inv())
         deleted_faces.push_back(i);
   }

   geom.delete_faces(deleted_faces);
}

int main(int argc, char *argv[])
{
   planar_opts opts;
   opts.process_command_line(argc, argv);

   char errmsg[MSG_SZ];
   col_geom_v geom;
   if(!geom.read(opts.ifile, errmsg))
      opts.error(errmsg);
   if(*errmsg)
      opts.warning(errmsg);

   if (opts.orient) {
      geom.orient_reverse();
      geom.orient();
      geom_info rep(geom);
      if(!rep.is_oriented()) {
         snprintf(errmsg, MSG_SZ, "input file contains a non-orientable geometry");
         opts.warning(errmsg, 'O');
      }
   }

   if (opts.cmy_mode)
      for(unsigned int i=0; i<geom.faces().size(); i++)
         geom.set_f_col(i,rgb_complement(geom.get_f_col(i), opts.ryb_mode));

   if (opts.face_color_method)
      color_by_plane(geom, opts.face_color_method, opts.face_opacity, opts.map, opts.use_unit_vectors, opts.epsilon);

   if (opts.planar_merge_type)
      planar_merge(geom, opts.planar_merge_type, opts.polygon_fill_type, opts.hole_detection, opts.use_unit_vectors,
                   opts.winding_rule, opts.color_by_winding_number, opts.zero_density_color, opts.zero_density_force_blend, opts.brightness_adj,
                   opts.color_system_mode, opts.sat_power, opts.sat_threshold, opts.value_power, opts.value_advance, opts.alpha_mode, opts.ryb_mode, opts.epsilon);

   // resolve 'indexes' of winding numbers
   if (opts.color_by_winding_number)
      resolve_winding_number_indexes(geom, opts.map, opts.map_negative);

   // add vertices to 'stitch' faces with dangling edges due to tiling or merging
   if (opts.stitch_faces)
      stitch_faces_on_seams(geom, opts.epsilon);

   if (opts.blend_edges)
      blend_edges(geom, opts.blend_edges,
                  opts.color_system_mode, opts.sat_power, opts.sat_threshold, opts.value_power, opts.value_advance, opts.alpha_mode, opts.ryb_mode);

   if (opts.delete_invisible_faces)
      delete_invisible_faces(geom);

   if(!geom.write(opts.ofile, errmsg))
      opts.error(errmsg);

   return 0;
}
