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
   Name: planar.cc
   Description: Various utilities for polyhedra with planar faces
   Project: Antiprism - http://www.antiprism.com
*/


#include <stdio.h>
#include <stdlib.h>

#include <ctype.h>
#include <unistd.h>

#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include "../base/antiprism.h"
#include "../base/bbox.h"

using std::string;
using std::vector;
using std::pair;
using std::make_pair;
using std::map;
using std::swap;
using std::min;
using std::max;
 

class planar_opts: public prog_opts {
   public:
      string ifile;
      string ofile;
      
      bool how_planar;
      char face_color_method;
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

      int sig_compare;
      double epsilon;

      planar_opts(): prog_opts("planar"),
                        how_planar(false),
                        face_color_method('\0'),
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
"  -x        measure planar-ness (0 = planar)\n"
"  -l <lim>  minimum distance for unique vertex locations as negative exponent\n"
"               (default: %d giving %.0e)\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\nColor Blending Options (for option -F B)\n"
"  -M <mode> color blending mode. HSV=1  HSL=2 (default)  RGB=3\n"
"  -s <sat>  HSV/HSL saturation curve. Greather than 0  (default: 1)\n"
"               1.0 - no curve. lower than 1.0 makes blends more pastel\n"
"  -t <val>  HSV/HSL threshold to use average saturation. (default: 1)\n"
"               between 0.0 (all averaging) and 1.0 (no averaging)\n"
"  -v <val>  HSV/HSL value curve (default: 0)\n"
"               simulates subtractive coloring for blending 3 or more colors\n"
"               RGB: Red+Green+Blue = White   Cyan+Magenta+Yellow = Black\n"
"               RYB: Red+Yellow+Blue = Black  Green+Magenta+Orange = White\n"
"               1.0 - no curve. lower than 1.0 number makes blends lighter\n"
"               0.0 - use average value instead\n"
"  -u <val>  HSV/HSL value advance. Rotates meaning of white and black\n"
"               valid values 0.0 to 120.0 degrees (default: 0)\n"
"  -a <int>  alpha to use for blend. average=1  minimum=2  maximum=3 (default)\n"
"  -y        RYB mode. Blend colors as in Red-Yellow-Blue color wheel\n"
"  -c        CMY mode. Complementary colors.  RGB->(RYB/GMO)->CMY->blend\n"
"\nColoring Options\n"
"  -f <opt>  face color option\n"
"               d - unique color for faces on the same plane\n"
"               f - unique color for faces on opposite planes\n"
"               b - blend existing colors of overlaping planes\n"
"                 note: -T and -m do not apply to -f b\n"
"  -T <tran> face opacity for color by symmetry. valid range from 0 to 255\n"
"  -m <maps> color maps for faces to be tried in turn (default: compound)\n"
"\n",prog_name(), help_ver_text, int(-log(::epsilon)/log(10) + 0.5), ::epsilon);
}

void planar_opts::process_command_line(int argc, char **argv)
{
   opterr = 0;
   char c;
   char errmsg[MSG_SZ];

   string id;
   string map_file;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hl:xM:s:t:v:u:a:cyf:T:m:o:")) != -1 ) {
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
            
         case 'x':
            how_planar = true;
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
            if(!strlen(optarg)==1 || !strchr("dfb", *optarg))
               error("color method must be d, f or b");
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

   if (!map_file.size())
      map_file = "compound";
   if(!map.init(map_file.c_str(), errmsg))
      error(errmsg, c);

   epsilon = (sig_compare != INT_MAX) ? pow(10, -sig_compare) : ::epsilon;
}

void geom_dump(const col_geom_v &geom, string s)
{
   char errmsg[MSG_SZ];
   char filename[MSG_SZ];
   sprintf(filename,"geom_dump_%s.off",s.c_str());
   geom.write(filename, errmsg);
}

double get_iq(const geom_if &geom)
{
   double iq = 0;
   geom_info rep(geom);
   if (rep.num_verts() > 2) {
      // make sure area is non-zero
      if (rep.face_areas().sum)
         iq = rep.isoperimetric_quotient();
   }
   return iq;
}

double measure_planarity(const geom_if &geom)
{
   double planarity = 0;
   
   for(unsigned int i=0; i<geom.faces().size(); i++) {
      geom_v tgeom;
      vector<int> face = geom.faces(i);
      for(unsigned int j=0; j<face.size(); j++)
         tgeom.add_vert(geom.verts(face[j]));
      tgeom.add_face(face);
      int ret = tgeom.set_hull();
      double iq = (ret > 2) ? get_iq(tgeom) : 0;
      fprintf(stderr,"face(%d) iq = %.17g\n",i,iq);
      planarity += iq;
      tgeom.clear_all();
   }
   
   return planarity;
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

class face_normal
{
   private:
      vec3d normal;
      bool out;
      bool hemi;
      
   public:
      face_normal() { out = true; hemi = false; }
      face_normal(const geom_if &geom, int idx, vec3d C=vec3d(), double eps=epsilon);

      virtual ~face_normal() {}

      bool is_set() const { return normal.is_set(); }

      vec3d raw() const { return normal; }      
      vec3d outward() const { return (out ? normal : -normal); }
      vec3d inward() const { return (!out ? normal : -normal); }
      
      bool is_hemispherical() const { return hemi; }
      bool is_outward() const { return out; }
      bool is_inward() const { return !out; }
};

face_normal::face_normal(const geom_if &geom, int idx, vec3d C, double eps)
{
   const vector<vector<int> > &faces = geom.faces();
   const vector<vec3d> &verts = geom.verts();
   
   normal = face_norm(verts, faces[idx]);
//normal.dump("normal");
   
   if (!C.is_set())
      C = centroid(verts);
//(verts[faces[idx][0]]).dump("verts[faces[idx][0]]");
//C.dump("C");
//(verts[faces[idx][0]]-C).dump("verts[faces[idx][0]]-C)");


//vec3d fcent = centroid(verts, faces[idx]);
   //double D = vdot(fcent-C, normal);
   double D = vdot(verts[faces[idx][0]]-C, normal);

   out = true;
   hemi = false;
   //if (double_equality(D, 0.0, eps) || !compare(centroid(verts, faces[idx]), C, eps))
   //if (double_equality(D, 0.0, eps) || double_equality(vdot(verts[0]-C, normal), 0.0, eps))
   if (double_equality(D, 0.0, eps) || double_equality(vdot(centroid(verts, faces[idx])-C, normal), 0.0, eps))
      hemi = true;
   else
   if (D < -eps)
      out = false;
}

class face_normals
{
   private:
      vector<face_normal> normals;
      
   public:
      face_normals() {};
      face_normals(const geom_if &geom, double eps=epsilon) { refresh(geom,eps); }
      
      virtual ~face_normals() {}
      
      int size() const { return (normals.size()); }
      bool in_range(unsigned int idx) const { return idx < (unsigned int)size(); }
      bool is_set() const { return (size() > 0); }
      
      void refresh(const geom_if &, double eps=epsilon);
      face_normal operator [](unsigned int idx) const { return ((is_set() && in_range(idx)) ? normals[idx] : face_normal()); }
};

void face_normals::refresh(const geom_if &geom, double eps)
{
   vec3d C = centroid(geom.verts());
   
   normals.clear();
   for(unsigned int i=0;i<geom.faces().size();i++)
      normals.push_back(face_normal(geom,i,C,eps));
}

void build_coplanar_faces_list(const geom_if &geom, vector<vector<int> > &coplanar_faces_list, vector<face_normal> &coplanar_normals, const face_normals &fnormals,
                               bool point_outward, bool fold_normals, double eps)
{
   vector<pair<vec3d, int> > face_normal_table;
   
   vector<int> hemispherical_faces;
   for(unsigned int i=0; i<geom.faces().size(); i++) {
      pair<vec3d, int> face_normal_pair;
      vec3d normal = point_outward ? fnormals[i].outward() : fnormals[i].raw();
      if (fnormals[i].is_hemispherical())
         hemispherical_faces.push_back(i);
      face_normal_pair.first = normal.unit();
      face_normal_pair.second = i;
      face_normal_table.push_back(face_normal_pair);
   }

   // fold in hemispherical normals
   if (hemispherical_faces.size() > 1) {
      for(unsigned int i=0; i<hemispherical_faces.size()-1; i++) {
         for(unsigned int j=i+1; j<hemispherical_faces.size(); j++) {
            if (!compare(face_normal_table[hemispherical_faces[i]].first,-face_normal_table[hemispherical_faces[j]].first,eps)) {
               face_normal_table[hemispherical_faces[j]].first = -face_normal_table[hemispherical_faces[j]].first;
            }
         }
      }
   }
   
   // if rest of normals are folded with specific option. redundent for hemi's, they get passed over this time
   if (fold_normals) {
      for(unsigned int i=0; i<face_normal_table.size()-1; i++) {
         for(unsigned int j=i+1; j<face_normal_table.size(); j++) {
            if (!compare(face_normal_table[i].first,-face_normal_table[j].first,eps)) {
               face_normal_table[j].first = -face_normal_table[j].first;
            }
         }
      }
   }
   
   stable_sort( face_normal_table.begin(), face_normal_table.end(), vert_cmp(eps) );

// debug
//for(unsigned int i=0; i<face_normal_table.size(); i++) {
//   vec3d v = face_normal_table[i].first;
//   fprintf(stderr,"%d : %.17lf %.17lf %.17lf\n",i,v[0],v[1],v[2]);
//}

   vector<int> coplanar_face;
   // at least one face is on a plane
   coplanar_face.push_back(face_normal_table[0].second);
   // find faces on same plane and group them together
   for(unsigned int i=1; i<face_normal_table.size(); i++) {
      if (compare(face_normal_table[i].first,face_normal_table[i-1].first,eps)) {
         coplanar_faces_list.push_back(coplanar_face);
         //coplanar_normals.push_back(face_normal_table[i-1].first);
         coplanar_normals.push_back(fnormals[face_normal_table[i-1].second]);
         coplanar_face.clear();
      }
// debug
//else
//fprintf(stderr,"%d and %d are equal\n",i-1,i);
      coplanar_face.push_back(face_normal_table[i].second);
   }
   coplanar_faces_list.push_back(coplanar_face);
   //coplanar_normals.push_back(face_normal_table[face_normal_table.size()-1].first);
   coplanar_normals.push_back(fnormals[face_normal_table[face_normal_table.size()-1].second]);
   coplanar_face.clear();

   face_normal_table.clear();
}

// http://geometryalgorithms.com/Archive/algorithm_0106/algorithm_0106.htm
bool lines_nearest_points(const vec3d &P0, const vec3d &P1, const vec3d &Q0, const vec3d &Q1,
                          vec3d &P, vec3d &Q, double eps)
{
    vec3d u = (P1-P0);
    vec3d v = (Q1-Q0);
    vec3d w = P0-Q0;
    double a = vdot(u, u);
    double b = vdot(u, v);
    double c = vdot(v, v);
    double d = vdot(u, w);
    double e = vdot(v, w);
    double D = a*c-b*b;
    double sc = (b*e-c*d)/D;
    double tc = (a*e-b*d)/D;
    P = P0 + sc*u;
    Q = Q0 + tc*v;
    return (D > eps);
}

vec3d lines_intersection(const vec3d &P0, const vec3d &P1, const vec3d &Q0, const vec3d &Q1, double eps)
{
   vec3d N1, N2;
   // lines might not be parallel and still miss so check if nearest points is not zero
   if(lines_nearest_points(P0, P1, Q0, Q1, N1, N2, eps) && ((N1-N2).mag() < eps))
      return (N1+N2)/2.0;
   else
      return vec3d();
}

bool in_segment(const vec3d &P, const vec3d &Q0, const vec3d &Q1, double eps)
{
   bound_box bb;
   vector<vec3d> points(2);
   
   points[0] = Q0;
   points[1] = Q1;

   bb.add_points(points);
   vec3d min = bb.get_min();
   vec3d max = bb.get_max();

   return (compare(min,P,eps)<=0 && compare(max,P,eps)>=0);
}

vec3d lines_intersection_in_segments(const vec3d &P0, const vec3d &P1, const vec3d &Q0, const vec3d &Q1, double eps)
{
   vec3d intersection = lines_intersection(P0, P1, Q0, Q1, eps);
   
   if (intersection.is_set() && (!in_segment(intersection, P0, P1, eps) || !in_segment(intersection, Q0, Q1, eps)))
      intersection = vec3d();

   return intersection;
}

vec3d point_in_segment(const vec3d &P, const vec3d &Q0, const vec3d &Q1, double eps)
{
   vec3d intersection = nearest_point(P, Q0, Q1);

   if (intersection.is_set() && (compare(P,intersection,eps) || !in_segment(intersection, Q0, Q1, eps)))
      intersection = vec3d();

   return intersection;
   //return (intersection.is_set() && !compare(P,intersection,eps) && in_segment(intersection, Q0, Q1, eps));
}

// furshished by Adrian Rossiter
int side_of_line(const vec3d &P, const vec3d &A, const vec3d &B, const int idx0, const int idx1, double eps)
{
   double a = B[idx0] - A[idx0];
   double b = B[idx1] - A[idx1];
   double c = P[idx0] - A[idx0];
   double d = P[idx1] - A[idx1];
   double det = a*d - b*c;

   if(det<-eps)
      return -1;
   else if(det>eps)
      return 1;
   else
      return 0;
}
   
bool is_point_inside_triangle_or_edge(const vec3d &P, const vec3d &A, const vec3d &B, const vec3d &C,
                                  const int idx0, const int idx1, const int sign, double eps)                          
{
   return (side_of_line(P,A,B,idx0,idx1,eps)*sign>=-eps && side_of_line(P,B,C,idx0,idx1,eps)*sign>=-eps && side_of_line(P,C,A,idx0,idx1,eps)*sign>=-eps);
}

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

// find if a point is in a polygon
// the polygon is input as the only face in a geom
// the polygon is triangulated so it is changed from the original
// set include_edges = true if point is allowed to land on outer edges or vertices
// the normal is expected to be a unit normal
bool is_point_inside_polygon(col_geom_v &polygon, const vec3d &P, const vec3d &normal, bool include_edges, bool fast_on_plane_check, double eps)
{
   const vector<vector<int> > &faces = polygon.faces();
   const vector<vec3d> &verts = polygon.verts();
   
   // quick check to see if P is actually on the polygon
   // funished by Adrian Rossiter; is there any distance from P to the unit normal?
   if (fast_on_plane_check)
      if (!double_equality(vdot(verts[0]-P, normal), 0.0, eps))
         return false;
   
   int idx;
   int sign;

   project_using_normal(normal, idx, sign);
//fprintf(stderr,"is_point_inside_polygon: idx = %d  sign = %d\n",idx,sign);
   
   // before triangulation collect outer edge vertices for edge testing after face testing
   vector<pair<int,int> > outer_edge_vert_pairs;
   for(unsigned int i=0; i<faces[0].size(); i++)
      outer_edge_vert_pairs.push_back(make_pair(faces[0][i],faces[0][(i+1)%faces[0].size()]));
   
   polygon.triangulate(col_val::invisible);
   
   // test faces
   bool answer = false;
   for(unsigned int i=0; i<faces.size(); i++) {
      if (is_point_inside_triangle_or_edge(P,verts[faces[i][0]],verts[faces[i][1]],verts[faces[i][2]],(idx+1)%3,(idx+2)%3,sign,eps)) {
         answer = true;
         break;
      }
   }
   
   // test edges
   if (!include_edges && answer) {
      for(unsigned int i=0; i<outer_edge_vert_pairs.size(); i++) {
         vec3d v1 = verts[outer_edge_vert_pairs[i].first];
         vec3d v2 = verts[outer_edge_vert_pairs[i].second];
         if ((point_in_segment(P, v1, v2, eps)).is_set()) {
            answer = false;
            break;
         }
      }
   }
   
   return answer;
}

int intersection_is_P(const vec3d &intersection_point, const vec3d &P0, const vec3d &P1, double eps)
{
   int which_one = 0;
   if (!compare(intersection_point,P0,eps))
      which_one = 1;
   else
   if (!compare(intersection_point,P1,eps))
      which_one = 2;
   return which_one;
}

int intersection_is_Q(const vec3d &intersection_point, const vec3d &Q0, const vec3d &Q1, double eps)
{
   int which_one = 0;
   if (!compare(intersection_point,Q0,eps))
      which_one = 1;
   else
   if (!compare(intersection_point,Q1,eps))
      which_one = 2;
   return which_one;
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
 
   const vector<vector<int> > &edges = geom.edges();
   
   vector<int> new_edge = make_edge(v_idx1, v_idx2);
   
   bool edge_inserted = true;
   for(unsigned int i=0; i<edges.size(); i++) {
      if (new_edge[0] == edges[i][0] && new_edge[1] == edges[i][1]) {
         edge_inserted = false;
         break;
      }
   }
   
//fprintf(stderr,"inserting edge %d %d\n", v_idx1, v_idx2);
   if (edge_inserted)
      geom.add_col_edge(new_edge, ecol);

   return edge_inserted;
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
         if ((intersection_point.is_set()) && (!intersection_is_Q(intersection_point, Q0, Q1, eps)))
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
   return true;
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
   return true;
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
   
   for(unsigned int i=0; i<face1.size(); i++) {
      vec3d P0 = verts[face1[i]];
      vec3d P1 = verts[face1[(i+1)%face1.size()]];
      for(unsigned int j=0; j<face2.size(); j++) {
         vec3d Q0 = verts[face2[j]];
         vec3d Q1 = verts[face2[(j+1)%face2.size()]];
         
         vec3d intersection_point = lines_intersection_in_segments(P0, P1, Q0, Q1, eps);
         if (intersection_point.is_set())
            return true;
      }
   }
  
   return false;
}

// of a set of coplanar faces, find those that actually overlap and group them
void find_coplanar_clusters(const col_geom_v &geom, vector<vector<int> > &coplanar_clusters, const vector<int> &coplanar_face_list, double eps)
{
   vector<pair<int, int> > overlap_pairs;
   
   for(unsigned int i=0; i<coplanar_face_list.size()-1; i++) {
      for(unsigned int j=i+1; j<coplanar_face_list.size(); j++) {
         if (check_if_faces_intersect(geom, coplanar_face_list[i], coplanar_face_list[j], eps)) {
            overlap_pairs.push_back(make_pair(coplanar_face_list[i], coplanar_face_list[j]));
         }
      }
   }
   
//for(unsigned int i=0; i<overlap_pairs.size(); i++) {
//   fprintf(stderr,"overlap: %d %d\n",overlap_pairs[i].first,overlap_pairs[i].second);
//}
//fprintf(stderr,"\n");

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
//fprintf(stderr,"build_angle_map: idx = %d  sign = %d\n",idx,sign);

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
//fprintf(stderr,"angle at %d %d = %lf\n",i,vcons[j],angle);
      }
   }
}

void build_turn_map(const geom_if &geom, map<pair<int, int>, int> &turn_map, map<pair<int, int>, double> &angle_map)
{
   const vector<vector<int> > &edges = geom.edges();
 
   for(unsigned int i=0; i<edges.size(); i++) {
//fprintf(stderr,"edge = %d\n",i);
      for(unsigned int j=0; j<2; j++) { 
         vector<int> vcons;
         int a = edges[i][!j ? 0 : 1];
         int b = edges[i][!j ? 1 : 0];
//fprintf(stderr,"a b = %d %d\n",a,b);
         double base_angle = angle_map[make_pair(b,a)];
         find_connections(geom,vcons,b);
         vector<pair<double, int> > angles;
         for(unsigned int k=0; k<vcons.size(); k++) {
            int c = vcons[k];
            // can't go backwards
            if (c == a)
               continue;
            double angle = angle_map[make_pair(b,c)];
//fprintf(stderr,"angle at %d %d %d = %lf\n", a, b, c, angle);
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

/* debug   
   for(unsigned int i=0; i<edges.size(); i++) {
      for(unsigned int j=0; j<2; j++) { 
         int a = edges[i][!j ? 0 : 1];
         int b = edges[i][!j ? 1 : 0];
         fprintf(stderr,"turn at %d %d = %d\n", a, b, turn_map[make_pair(a,b)]);
      }
   }
*/
}

void get_first_unvisited_triad(const geom_if &geom, map<pair<pair<int, int>, int>, bool > &visited, map<pair<int, int>, int> &turn_map, vector<int> &face)
{
   const vector<vector<int> > &edges = geom.edges();
 
   face.clear();
   for(unsigned int i=0; i<edges.size(); i++) {
      for(unsigned int j=0; j<2; j++) { 
         int a = edges[i][!j ? 0 : 1];
         int b = edges[i][!j ? 1 : 0];
         int c = turn_map[make_pair(a,b)];
         if (!visited[make_pair(make_pair(a,b),c)]) {
            face.push_back(a);
            face.push_back(b);
            face.push_back(c);
            return;
         }
      }
   }
   return;
}

void construct_faces(geom_if &geom, map<pair<int, int>, int> &turn_map)
{
   map<pair<pair<int, int>, int>, bool > visited;
   
   vector<int> face;
   get_first_unvisited_triad(geom,visited,turn_map,face);
   while(face.size()) {
      int fsz = face.size();
      int a = face[fsz-2];
      int b = face[fsz-1];
      int c = turn_map[make_pair(a,b)];
      if (c != face[0])
         face.push_back(c);
      else {
         geom.add_face(face);
/*
fprintf(stderr,"COMPLETED FACE = ");
for(unsigned int i=0; i<face.size(); i++)
   fprintf(stderr,"%d ",face[i]);
fprintf(stderr,"\n");
*/
         // mark turns visited and start a new face
         for(unsigned int i=0; i<face.size(); i++)
            visited[make_pair(make_pair(face[i],face[(i+1)%face.size()]),face[(i+2)%face.size()])] = true;
         get_first_unvisited_triad(geom,visited,turn_map,face);
      }
/*
static int count = 0;
count++;
if (count == 30) {
   fprintf(stderr,"construct_faces: bailout reached\n");
   return;
}
*/
   }
}

void analyze_faces(geom_if &geom, vector<int> &nonconvex_faces, map<pair<int, int>, double> &angle_map)
{
   const vector<vector<int> > &faces = geom.faces();
   vector<int> deleted_faces;
   int deleted_faces_count = 0; // should be only one
   
   for(unsigned int i=0; i<faces.size(); i++) {
      double angle_sum = 0;
      bool all_negative_turns = true;
      for(unsigned int j=0; j<faces[i].size(); j++) {
         int a = faces[i][j];
         int v = faces[i][(j+1)%faces[i].size()];
         int b = faces[i][(j+2)%faces[i].size()];
         double angle = angle_map[make_pair(v,b)]-angle_map[make_pair(v,a)]-180.0;
         if (angle < -180.0)
            angle += 360.0;
         else
         if (angle > 180.0)
            angle -= 360.0;
         if (angle > 0.0)
            all_negative_turns = false;
         angle_sum += angle;
//fprintf(stderr,"angle = %lf\n",angle);
      }
//fprintf(stderr,"for %d: angle_sum = %.15lf\n",i,angle_sum);
//fprintf(stderr,"for %d: all_negative_turns = %s\n",i,(all_negative_turns ? "true" : "false"));
      if ((int)floorf(angle_sum + 0.5) == 360) {
         deleted_faces.push_back(i);
         deleted_faces_count++;
      }
      else
      if (!all_negative_turns)
         nonconvex_faces.push_back(i-deleted_faces_count);
   }
/*
fprintf(stderr,"nonconvex faces are\n");
for(unsigned int i=0; i<nonconvex_faces.size(); i++) {
   fprintf(stderr,"%d ",nonconvex_faces[i]);
}
fprintf(stderr,"\n");
fprintf(stderr,"deleted faces are\n");
for(unsigned int i=0; i<deleted_faces.size(); i++) {
   fprintf(stderr,"%d ",deleted_faces[i]);
}
fprintf(stderr,"\n\n");
*/
   geom.delete_faces(deleted_faces);
}

void fill_in_faces(col_geom_v &geom, vector<int> &nonconvex_faces, const vec3d normal, double eps)
{
   map<pair<int, int>, double> angle_map;
   map<pair<int, int>, int> turn_map;
   build_angle_map(geom, angle_map, normal, eps);
   build_turn_map(geom, turn_map, angle_map);
   construct_faces(geom, turn_map);
   analyze_faces(geom, nonconvex_faces, angle_map);
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
      
//fprintf(stderr,"avg_col = %d %d %d %d\n\n",avg_col[0],avg_col[1],avg_col[2],avg_col[3]);
   return avg_col;
}

void sample_colors(col_geom_v &sgeom, const col_geom_v &cgeom, vector<face_normal> &original_normals, const vector<int> &nonconvex_faces,
      int color_system_mode, double sat_power, double sat_threshold, double value_power, double value_advance, int alpha_mode, bool ryb_mode, double eps)
{
   const vector<vector<int> > &sfaces = sgeom.faces();
   const vector<vector<int> > &cfaces = cgeom.faces();

   vector<col_val> cols;
   for(unsigned int i=0; i<sfaces.size(); i++) {
//fprintf(stderr,"sample_colors: i = %d\n",i);
      vec3d P;
      if(find(nonconvex_faces.begin(), nonconvex_faces.end(), (int)i) == nonconvex_faces.end())
         P = centroid(sgeom.verts(), sfaces[i]);
      else {
         vector<int> sface_idxs;
         sface_idxs.push_back(i);
         col_geom_v spolygon = faces_to_geom(sgeom, sface_idxs);
         spolygon.triangulate(col_val::invisible);
         P = centroid(spolygon.verts(), spolygon.faces()[0]);
      }
      for(unsigned int j=0; j<cfaces.size(); j++) {
//fprintf(stderr,"sample_color: j = %d\n",j);
         vector<int> face_idxs;
         face_idxs.push_back(j);
         col_geom_v polygon = faces_to_geom(cgeom, face_idxs);

/* debug
vector<int> sface_idxs;
sface_idxs.push_back(j);
col_geom_v spolygon = faces_to_geom(cgeom, sface_idxs);
spolygon.add_vert(P);
char filename[MSG_SZ];
sprintf(filename,"spolygon_%d_%d",i,j);
geom_dump(spolygon,filename);
*/

         if (is_point_inside_polygon(polygon, P, original_normals[j].raw().unit(), false, false, eps))
            cols.push_back(cgeom.get_f_col(j));
/* debug            
char filename[MSG_SZ];
sprintf(filename,"polygon_%d",j);
polygon.add_col_vert(P,col_val(1.0,1.0,1.0));
geom_dump(polygon,filename);
*/
      }
      sgeom.set_f_col(i,average_color(cols, color_system_mode, sat_power, sat_threshold, value_power, value_advance, alpha_mode, ryb_mode));
      cols.clear();
   }
}

void collect_original_normals(vector<face_normal> &original_normals, const vector<int> &coplanar_cluster, const face_normals &fnormals)
{
   for(unsigned int i=0; i<coplanar_cluster.size(); i++) {
      original_normals.push_back(fnormals[coplanar_cluster[i]]);
//fprintf(stderr,"face %d has normal\n",coplanar_cluster[i]);
//original_normals[i].raw().dump("raw");
//original_normals[i].raw().unit().dump("raw unit");
//fprintf(stderr,"normal is %s\n",original_normals[i].is_outward() ? "outward" : "inward");
   }
}

void blend_overlapping_colors(col_geom_v &geom, const vector<vector<int> > &coplanar_faces_list, vector<face_normal> &coplanar_normals, face_normals &fnormals,
                              int color_system_mode, double sat_power, double sat_threshold, double value_power, double value_advance, int alpha_mode, bool ryb_mode, double eps)
{
   col_geom_v bgeom;
   vector<int> deleted_faces;
   
   for(unsigned int i=0; i<coplanar_faces_list.size(); i++) {
      vector<vector<int> > coplanar_clusters;
      find_coplanar_clusters(geom, coplanar_clusters, coplanar_faces_list[i], eps);
      for(unsigned int j=0; j<coplanar_clusters.size(); j++) {
         if (coplanar_clusters[j].size() == 1)
            continue;
/*
bool stop_next = false;
fprintf(stderr,"\n");
for(unsigned int k=0; k<coplanar_clusters[j].size(); k++) {
   fprintf(stderr,"cluster %d contains %d\n",j,coplanar_clusters[j][k]);
   //if (coplanar_clusters[j][k] == 21)
   //   stop_next = true;
}
*/
         // load a geom with color faces. keep it and copy it.
         col_geom_v cgeom = faces_to_geom(geom, coplanar_clusters[j]);
         col_geom_v sgeom;
         sgeom.append(cgeom);

         make_skeleton(sgeom);
         // duplicate vertices and edges can cause problems
         sort_merge_elems(sgeom, "ve", eps);
         mesh_verts(sgeom, eps);
         mesh_edges(sgeom, eps);

         vector<int> nonconvex_faces;
         fill_in_faces(sgeom, nonconvex_faces, coplanar_normals[i].outward().unit(), eps);
         
         // original normals are needed for sampling colors
         vector<face_normal> original_normals;
         collect_original_normals(original_normals, coplanar_clusters[j], fnormals);        
         sample_colors(sgeom, cgeom, original_normals, nonconvex_faces, color_system_mode, sat_power, sat_threshold, value_power, value_advance, alpha_mode, ryb_mode, eps);

//if (i==0) {
//geom.clear_all();
//geom.append(sgeom);
//return;
//}
/*
if (cgeom.get_f_col(0) == col_val(0.0,0.39216,0.0)) {
   char filename[MSG_SZ];
   sprintf(filename,"darkgreen_%d",i);
   geom_dump(sgeom,filename);
}
*/

/*
if (stop_next) {         
geom.clear_all();
geom.append(sgeom);
return;
}
*/       
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

void color_by_plane(col_geom_v &geom, char face_color_method, int face_opacity, color_map_multi &map,
                    int color_system_mode, double sat_power, double sat_threshold, double value_power, double value_advance, int alpha_mode, bool ryb_mode, double eps)
{
   face_normals fnormals(geom,eps);

   vector<vector<int> > coplanar_faces_list;
   vector<face_normal> coplanar_normals;
   
   // pointed outward unless the d,D options was chosen
   bool point_outward = (face_color_method == 'd' || face_color_method == 'D') ? false : true;
   bool fold_normals = (face_color_method == 'f' || face_color_method == 'F') ? true : false;
   build_coplanar_faces_list(geom, coplanar_faces_list, coplanar_normals, fnormals, point_outward, fold_normals, eps);
   
   if (face_color_method == 'b') {
      geom.add_missing_impl_edges();
      blend_overlapping_colors(geom, coplanar_faces_list, coplanar_normals, fnormals,
                               color_system_mode, sat_power, sat_threshold, value_power, value_advance, alpha_mode, ryb_mode, eps);
   }
   else
   if (face_color_method == 'd' || face_color_method == 'f') {
      for(unsigned int i=0; i<coplanar_faces_list.size(); i++) {
         for(unsigned int j=0; j<coplanar_faces_list[i].size(); j++) {
            //fprintf(stderr,"plane = %d  face = %d\n", i, coplanar_faces_list[i][j]);
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

   if(opts.how_planar)
      fprintf(stderr,"\nthe total planarity measures to: %.17g\n", measure_planarity(geom));

   if (opts.cmy_mode)
      for(unsigned int i=0; i<geom.faces().size(); i++)
         geom.set_f_col(i,rgb_complement(geom.get_f_col(i), opts.ryb_mode));

   // face color by symmetry normals
   if (opts.face_color_method)
      color_by_plane(geom, opts.face_color_method, opts.face_opacity, opts.map,
                     opts.color_system_mode, opts.sat_power, opts.sat_threshold, opts.value_power, opts.value_advance, opts.alpha_mode, opts.ryb_mode, opts.epsilon);

   if(!geom.write(opts.ofile, errmsg))
      opts.error(errmsg);

   return 0;
}
