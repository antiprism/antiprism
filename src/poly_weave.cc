/*
   Copyright (c) 2003-20011, Adrian Rossiter

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
   Name: weave.cc
   Description: make a weave based on edge mid-points
   Project: Antiprism - http://www.antiprism.com
*/

#include <string.h>
#include <math.h>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <functional>
#include "../base/antiprism.h"

using std::string;
using std::vector;
using std::swap;
using std::not_equal_to;

void make_meta(const geom_if &geom, col_geom_v &meta)
{
   meta.clear_all();
   meta.add_verts(geom.verts());
   int f_start = meta.verts().size();
   for(unsigned int f=0; f<geom.faces().size(); f++)
      meta.add_vert(geom.face_cent(f));

   col_val light(1.0,0.8,0.6,0.5);
   col_val dark(0.1,0.3,0.6,0.5);
   map<vector<int>, vector<int> > ef_pairs;
   geom.get_edge_face_pairs(ef_pairs);
   map<vector<int>, vector<int> >::iterator ef_i = ef_pairs.begin();
   for(ef_i = ef_pairs.begin(); ef_i != ef_pairs.end(); ++ef_i) {
      // The edge and face pair make a quadrilateral
      // Add the centre to this quadrilateral
      int e_idx = meta.add_vert(geom.edge_cent(ef_i->first));
      // Add four triangles
      int idx;
      idx = meta.add_face(ef_i->first[0], e_idx, ef_i->second[0]+f_start, -1);
      meta.set_f_col(idx, light);
      idx = meta.add_face(ef_i->first[0], e_idx, ef_i->second[1]+f_start, -1);
      meta.set_f_col(idx, dark);
      idx = meta.add_face(ef_i->first[1], e_idx, ef_i->second[1]+f_start, -1);
      meta.set_f_col(idx, light);
      idx = meta.add_face(ef_i->first[1], e_idx, ef_i->second[0]+f_start, -1);
      meta.set_f_col(idx, dark);
   }
   meta.add_missing_impl_edges();
   coloring (&meta).e_one_col(col_val::invisible);
   coloring (&meta).v_one_col(col_val::invisible);
}

class weave_path
{
   private:
      vec3d point;
      vector<vector<double> > path_points;
      int path_type;
   public:
      enum {LINE, CURVE};
      bool init(const string &path, char *errmsg=0);
      void set_path_type(int typ) {path_type = typ;}
      vec3d get_point() const { return point; }
      void add_points(vector<vec3d> *pts,
            const vector<vec3d> &start, const vector<vec3d> &end) const;
};


bool weave_path::init(const string &path, char *errmsg)
{
   char errmsg2[MSG_SZ];
   path_points.clear();
   vector<char> path_str(path.size()+1);
   strcpy(&path_str[0], path.c_str());
   vector<char *> parts;
   int parts_sz = split_line(&path_str[0], parts, ":", true);
   if(!point.read(parts[0], errmsg2)) {
      if(errmsg)
         sprintf(errmsg, "path '%s': first (join) point: '%s'",
               path.c_str(), errmsg2);
      return false;
   }
   point /= point[0] + point[1] + point[2];
   for(int i=1; i<parts_sz; i++) {
      vector<double> nums;
      if(!read_double_list(parts[i], nums, errmsg2)) {
         if(errmsg)
            sprintf(errmsg, "path '%s': path point %d: '%s'",
                  path.c_str(), i+1, errmsg2);
         return false;
      }
      if(nums.size()>3) {
         if(errmsg)
            sprintf(errmsg, "path '%s': path point %d: more than three numbers",
                  path.c_str(), i+1);
         return false;
      }
      path_points.push_back(nums);
      //fprintf(stderr, "parts[%d]='%s'\n", i, parts[i]);
   }

   return true;
}

static vec3d get_control_point(vec3d P, vec3d Q, vec3d P_up)
{
   vec3d along = Q-P;
   vec3d axis = vcross(along, P_up);
   double turn_ang = angle_around_axis(along, P_up, axis) - M_PI/2;
   mat3d trans = mat3d::transl(P) *
                 mat3d::rot(axis, turn_ang) *
                 mat3d::transl(-P);
   vec3d pt = (2*P + Q)/3;    // point 1/3 along, chosen experimentally
   return trans * pt;
}

   

void weave_path::add_points(vector<vec3d> *pts,
      const vector<vec3d> &start, const vector<vec3d> &end) const
{
   vec3d along = end[0] - start[0];
   const int path_sz = path_points.size();
   vec3d Q0 = get_control_point(start[0], end[0], start[1]); // for CURVE
   vec3d Q1 = get_control_point(end[0], start[0], end[1]);   // for CURVE
   pts->push_back(start[0]);
   for(int i=0; i<path_sz; i++) {
      double X = (path_points[i].size()>0) ? path_points[i][0] : 0.0;
      double Y = (path_points[i].size()>1) ? path_points[i][1] : 0.0;
      double Z = (path_points[i].size()>2) ? path_points[i][2]
                                           : (i+1.0)/(path_sz+1);
      vec3d up = (Z*end[1] + (1-Z)*start[1]).unit(); // rotation would be better
      vec3d side = vcross(along, up).unit();
      vec3d origin;
      if(path_type==LINE)
         origin = (Z*end[0] + (1-Z)*start[0]);
      else if(path_type==CURVE) {
         const double z = 1-Z;
         origin = z*z*z*start[0] + 3*z*z*Z*Q0 + 3*z*Z*Z*Q1 + Z*Z*Z*end[0];
      }

      pts->push_back(origin + up*X + side*Y);
   }
}

class weave_pattern
{
   private:
      vector<int> ops;
      vector<weave_path> paths;
      unsigned char start_faces;

      mutable vector<int>::const_iterator ops_i;
      mutable vector<weave_path>::const_iterator paths_i;

   public:
      enum {END=-1, V=0, E, F, P};
      bool set_pattern(const string &pat, char *errmsg=0);
      //string get_pattern() const;
      unsigned char get_start_faces() const { return start_faces; }

      void start_op() const { ops_i=ops.begin(); paths_i=paths.begin(); }
      void next_op() const
         { if(ops_i!=ops.end()) { ++ops_i; if(*ops_i==P) ++paths_i; } }
      int get_op() const { if(ops_i==ops.end()) return END; else return *ops_i;}
      const weave_path &get_path() const {return *paths_i; }
};


bool weave_pattern::set_pattern(const string &pat, char *errmsg)
{
   ops.clear();
   paths.clear();
   start_faces = 5;                    // 'left'/'dark' faces 
   int path_type = weave_path::CURVE;  // path points to follow a curve
   
   bool reverse = false;
   int pat_sz = pat.size();
   int pos=0;
   while(pos<pat_sz) {
      int len;
      // path points
      bool add_default_point = !paths.size() && strchr("vefVEF", pat[pos]);
      if((len=strspn(pat.substr(pos).c_str(), "0123456789.,-+:")) ||
            add_default_point) {
         ops.push_back(P);
         weave_path path;
         if(add_default_point)
            path.init("0,1,0", errmsg);
         else if(!path.init(pat.substr(pos, len).c_str(), errmsg))
            return false;
         path.set_path_type(path_type);
         paths.push_back(path);
         if(add_default_point)
            continue; // reprocess char that triggered adding of default point
      }
      else if(strchr("t", pat[pos])) {
         if((int)pat.size() == pos+1) {
            sprintf(errmsg, "'t' is last character, must be followed by l, r or b\n");
            return false;
         }
         char tris = pat[pos+1];
         if(!strchr("lrb", tris) ) {
            sprintf(errmsg, "'t' followed by '%c', must be followed by l, r or b\n", pat[pos]);
            return false;
         }
         start_faces = 5*(tris=='l' || tris=='b') + 10*(tris=='r' || tris=='b');
         len = 2;
      }

      // path shapes
      else if('L' == pat[pos]) {
         path_type = weave_path::LINE;
         if(paths.size())
            paths.back().set_path_type(path_type);
      }
      else if('C' == pat[pos]) {
         path_type = weave_path::CURVE;
         if(paths.size())
            paths.back().set_path_type(path_type);
      }
      
      // mirrors
      else if('v' == pat[pos])
         ops.push_back(V);
      else if('e' == pat[pos])
         ops.push_back(E);
      else if('f' == pat[pos])
         ops.push_back(F);
      
      // rotations
      else if('R' == pat[pos])
         reverse = !reverse;
      else if('V' == pat[pos]) {
         if(!reverse) {
            ops.push_back(E);
            ops.push_back(F);
         }
         else {
            ops.push_back(F);
            ops.push_back(E);
         }
      }
      else if('E' == pat[pos]) {
         ops.push_back(F);
         ops.push_back(V);
      }
      else if('F' == pat[pos]) {
         if(!reverse) {
            ops.push_back(V);
            ops.push_back(E);
         }
         else {
            ops.push_back(E);
            ops.push_back(V);
         }
      }
      
      else {
         if(errmsg)
            sprintf(errmsg, "invalid character '%c' in position %d",
                  pat[pos], pos+1);
         return false;
      }


      if(len)
         pos += len;
      else
         pos++;
   }
   return true;
}

/*
string weave_pattern::get_pattern() const
{
   string ret;
   for(unsigned int i=0; i<ops.size(); i++) {
      if(ops[i] == P)
         ret += "P";
         //ret += vtostr(paths[pt_idx++].get_path().get_point(), ",", 5);
      else if(ops[i] == V)
         ret += "v";
      else if(ops[i] == E)
         ret += "e";
      else if(ops[i] == F)
         ret += "f";
      else {
         ret = "X";
         break;
      }
   }
   ret += "t";
   if(start_faces==5)
      ret += 'l';
   if(start_faces==10)
      ret += 'r';
   if(start_faces==15)
      ret += 'b';

   return ret;
}
*/

class weave
{
   private:
      vector<weave_pattern> pats;

      col_geom_v meta;
      vector<vector<int> > nbrs;
      vector<vec3d> vert_norms;

      bool find_nbrs();
      vector<vec3d> point_on_face(int f_idx, const vec3d &crds) const;
      void add_circuit(col_geom_v &wv, int start_idx, const weave_pattern &pat,
            vector<bool> &seen) const;
      const vector<weave_pattern> &get_pats() const {return pats; }

   public:
      bool set_geom(const geom_if &geom);
      void add_pattern(const weave_pattern &pattern) { pats.push_back(pattern);}
      bool add_pattern(const string &pat, char *errmsg=0);
      void make_weave(col_geom_v &wv) const;
      
      const geom_if &get_meta() const {return meta; }
};

bool weave::find_nbrs()
{
   map<vector<int>, vector<int> > ef_pairs;
   meta.get_edge_face_pairs(ef_pairs, false);
   map<vector<int>, vector<int> >::iterator ef_i = ef_pairs.begin();
   
   // Find the neighbour face opposite each VEF vertex 
   nbrs.resize(meta.faces().size(), vector<int>(3));
   for(unsigned int f=0; f<meta.faces().size(); f++)
      for(int i=0; i<3; i++) {
         vector<int> e(2);
         e[0] = meta.faces_mod(f, i+1);
         e[1] = meta.faces_mod(f, i+2);
         if(e[0]>e[1])
            swap(e[0], e[1]);
         map<vector<int>, vector<int> >::iterator ef_i = ef_pairs.find(e);
         if(ef_i==ef_pairs.end())
            return false;
         nbrs[f][i] = (ef_i->second[0]!=(int)f)?ef_i->second[0]:ef_i->second[1];
      }
   return true;
}


inline vector<vec3d> weave::point_on_face(int f_idx, const vec3d &crds) const
{
   vector<vec3d> ret(2);
   // point coordinates
   ret[0]=crds[weave_pattern::V]*meta.face_v(f_idx, weave_pattern::V) +
          crds[weave_pattern::E]*meta.face_v(f_idx, weave_pattern::E) + 
          crds[weave_pattern::F]*meta.face_v(f_idx, weave_pattern::F);

   // point normal
   ret[1]=crds[weave_pattern::V]*vert_norms[meta.faces(f_idx,weave_pattern::V)]+
          crds[weave_pattern::E]*vert_norms[meta.faces(f_idx,weave_pattern::E)]+
          crds[weave_pattern::F]*vert_norms[meta.faces(f_idx,weave_pattern::F)];
   ret[1].to_unit();
   return ret;
}

void weave::add_circuit(col_geom_v &wv, int start_idx, const weave_pattern &pat,
      vector<bool> &seen) const
{
   // Apply pattern until circuit completes
   vector<vec3d> prev_pt;
   bool finish = false;
   int start_v_sz = wv.verts().size();
   int idx = start_idx;
   while(true) {        
      seen[idx] = true;
      pat.start_op();
      while(pat.get_op()!=weave_pattern::END) {
         //fprintf(stderr, "op=%d\n", pat.get_op());
         //fprintf(stderr, "nbrs[%d] = %d, %d, %d\n", idx, nbrs[idx][0],nbrs[idx][1],nbrs[idx][2]);
         //fprintf(stderr, "nbrs[%d][%d] = %d, idx=%d, start_idx=%d\n", idx, op, nbrs[idx][op], idx, start_idx);
         if(pat.get_op()==weave_pattern::P) {
            vector<vec3d> pt = point_on_face(idx, pat.get_path().get_point());
            if(prev_pt.size())
               pat.get_path().add_points(wv.get_verts(), prev_pt, pt);
            prev_pt = pt;
         }
         else {
            idx = nbrs[idx][pat.get_op()]; // move to next triangle
         }
         pat.next_op();
         //if(idx==start_idx && prev_pt.size())   // circuit complete
         //   finish = true;
      }
      if(finish)
         break;
      if(idx==start_idx && prev_pt.size())   // circuit complete
         finish=true;
   }

   const int f_sz = wv.verts().size()-start_v_sz;
   vector<int> face(f_sz);
   for(int i=0; i<f_sz; i++)
      face[i] = start_v_sz+i;
   wv.add_face(face);
}

static void reverse_odd_faces(geom_if &geom)
{
   const int f_sz = geom.faces().size();
   for(int i=0; i<f_sz; i++)
      if(i%2)
        reverse((*geom.get_faces())[i].begin(), (*geom.get_faces())[i].end());
}


bool weave::set_geom(const geom_if &geom)
{
   make_meta(geom, meta);
   find_nbrs();
   reverse_odd_faces(meta);
   vert_norms = meta.get_info().get_vert_norms();
   reverse_odd_faces(meta);
   return true;
}


bool weave::add_pattern(const string &pat, char *errmsg)
{
   weave_pattern pattern;
   if(pattern.set_pattern(pat, errmsg)) {
      add_pattern(pattern);
      return true;
   }
   else
      return false;
}

void weave::make_weave(col_geom_v &wv) const
{
   //for(unsigned int i=0; i<nbrs.size(); i++)
   //   fprintf(stderr, "nbrs[%d] = %d, %d, %d\n", i, nbrs[i][0],nbrs[i][1],nbrs[i][2]);

   wv.clear_all();
   int faces_sz = meta.faces().size();
   for(unsigned int p=0; p<pats.size(); p++) {
      //fprintf(stderr, "pattern is '%s'\n", pats[p].get_pattern().c_str());
      vector<bool> seen(faces_sz, false);
      unsigned char start_faces = pats[p].get_start_faces();
      for(int i=0; i<faces_sz; i++) {
         if(!seen[i] && (start_faces & (1<<(i%4))) )
            add_circuit(wv, i, pats[p], seen);
      }
   }
}









class wv_opts: public prog_opts
{
   private:
   
   public:
      bool add_meta;
      weave wv;
      bool use_default_pattern;
      string ifile;
      string ofile;

      wv_opts() : prog_opts("poly_weave"),
                  add_meta(false),
                  use_default_pattern(true)
                  {}

      void process_command_line(int argc, char **argv);
      void usage();
};

void wv_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Read a file in OFF format and make a weave following a specified pattern.\n"
"The polyhedron faces are divided by a 'meta' operation into triangles\n"
"each having vertices which are a vertex V, edge centre E and face centre F.\n"
"A start point is positioned on one of these triangles, the next point is\n"
"found by using the pattern to step between triangles, leading to a circuit.\n"
"Intermediate points may be added along the steps to help with the weave.\n"
"If input_file is not given the program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"  -p <pat>  weave pattern (default 'VEF'), a series of one or more paths.\n"
"            A path is:\n"
"               C,L - path of intermediate points is curve (default),line\n"
"                     (affects following paths, set as first char of pattern)\n"
"               Barycentric coords V,E,F of initial point (default 0,1,0)\n"
"               Optional intermediate points on the path - ':' followed\n"
"                  by 0 to 3 coordinates for 'up' (def: 0), 'side' (def: 0),\n"
"                  'along' (def: equal spacing)\n"
"               Series of letters to step between triangles\n"
"                  V,E,F - step two triangles rotating about V,E,F\n"
"                  R     - R reverse direction of following rotations\n"
"                  v,e,f - step over side opposite V,E,F\n"
"               t followed by l,r,b (def: tl), start circuits from 'left',\n"
"                  'right', both triangles\n"
"  -a        add the 'meta'-transformed base\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}

     
void wv_opts::process_command_line(int argc, char **argv)
{
   char errmsg[MSG_SZ];
   opterr = 0;
   char c;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":ho:p:a")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         case 'p':
            if(!wv.add_pattern(optarg, errmsg))
               error(errmsg, 'p');
            use_default_pattern = false;
            break;
         
         case 'a':
            add_meta = true;
            break;

         case 'o':
            ofile = optarg;
            break;

         default:
            error("unknown command line error");
      }
   }

   if(use_default_pattern)
      wv.add_pattern("FEV");

   if(argc-optind > 1) {
      error("too many arguments");
      exit(1);
   }
   
   if(argc-optind == 1)
      ifile=argv[optind];

}




int main(int argc, char *argv[])
{
   wv_opts opts;
   opts.process_command_line(argc, argv);

   char errmsg[MSG_SZ];
   col_geom_v geom;
   if(!geom.read(opts.ifile, errmsg))
      opts.error(errmsg);
   if(*errmsg)
      opts.warning(errmsg);

   geom_info info(geom);
   if(!info.is_polyhedron())
      opts.error("base geometry is not a polyhedron");
   if(!info.is_orientable())
      opts.error("base polyhedron is not orientable");
   if(!info.is_oriented()) {
      opts.warning("base polyhedron is not oriented; it will be oriented.");
      geom.orient();
   }

   weave &wv = opts.wv;

   wv.set_geom(geom);
   col_geom_v wv_geom;
   wv.make_weave(wv_geom);
   if(opts.add_meta)
      wv_geom.append(wv.get_meta());

   if(!wv_geom.write(opts.ofile, errmsg))
      opts.error(errmsg);

   return 0;
}
   

