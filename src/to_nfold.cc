/*
   Copyright (c) 2013, Adrian Rossiter

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
   Name: twist.cc
   Description: twist edges of a polyhedron
   Project: Antiprism - http://www.antiprism.com
*/

#include <string.h>
#include <math.h>
#include <ctype.h>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include "../base/antiprism.h"

using std::string;
using std::vector;
using std::swap;



class nfold_opts: public prog_opts
{
   private:

   public:
      int num;
      int denom;

      string ifile;
      string ofile;

      nfold_opts(): prog_opts("to_nfold")
                 {}

      void process_command_line(int argc, char **argv);
      void usage();
};

void nfold_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] fraction [input_file]\n"
"\n"
"Generalise an axial model by changing its rotational symmetry. Read a model,\n"
"in OFF format, with an m-fold rotational axis on the z-axis, and create a\n"
"new model, generally non-planar, with the same relative connections, but\n"
"with an n-fold axis instead. fraction is given as n, or n/d (n and d\n"
"integers). Vertices of a face originally separated by x/m of a turn around\n"
"the z-axis will be separated by xd/n of a turn in the final model. If\n"
"input_file is not given the program reads from standard input. %s is based \n"
"on an idea by Bruce R. Gilson.\n"
"\n"
"Options\n"
"%s"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), prog_name(), help_ver_text);
}


void nfold_opts::process_command_line(int argc, char **argv)
{
   char errmsg[MSG_SZ];
   opterr = 0;
   char c;

   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":ho:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         case 'o':
            ofile = optarg;
            break;

         default:
            error("unknown command line error");
      }
   }

   if(argc-optind < 1) {
      error("not enough arguments");
      exit(1);
   }
   if(argc-optind > 2) {
      error("too many arguments");
      exit(1);
   }

   // read fraction
   char *p = strchr(argv[optind], '/');
   denom = 1;
   if(p!=0) {
      *p++='\0';
      if(!read_int(p, &denom, errmsg))
         error(errmsg, "fraction n/d, denominator");
   }


   if(!read_int(argv[optind], &num, errmsg))
      error(errmsg, "fraction n/d, numerator");
   if(num<2)
      error("must be an integer 2 or greater", "fraction n/d, numerator");
   if(denom < 1)
      error("must be 1 or greater", "fraction n/d, denominator");
   if(denom%num==0)
      error("numerator cannot divide the denominator", "fraction n/d, numerator");


   if(argc-optind == 2)
      ifile=argv[optind+1];

}

class cyl_vertex
{
   public:
      double radius;
      double height;
      double angle;

      cyl_vertex(double radius, double height, double angle) :
         radius(radius), height(height), angle(angle) {}

      cyl_vertex(vec3d v)
      {
         radius = vec3d(v[0], v[1], 0).mag();
         height = v[2];
         angle = angle_around_axis(vec3d::X, v, vec3d::Z);
      }

      vec3d to_vec3d() const
      { return vec3d(radius*cos(angle), radius*sin(angle), height); }
};

class cyl_face
{
   public:
      double start_angle;
      vector<cyl_vertex> cyl_verts;

      cyl_face(const col_geom_v &geom, const vector<int> &face);
      void add_to_geom(col_geom_v &geom, int from_n, int to_n, int to_d);
      void add_chain_to_geom(col_geom_v &geom, int from_n, int to_n, int to_d);
};


cyl_face::cyl_face(const col_geom_v &geom, const vector<int> &face_orig)
{
   double eps=10e-8;
   vector<int> face = face_orig;
   if(face.size()>2) {
      vec3d N0, N1;
      lines_nearest_points(vec3d(0,0,0), vec3d(0,0,1),
            geom.verts(face[0]), geom.verts(face[1]), N0, N1);
      if((N0-N1).mag()<eps) {
         //fprintf(stderr, "cyclic face starts with axis intersecting edge\n");
         rotate(face.begin(), face.begin()+2, face.end());
      }
   }



   //cyl_vertex v_start(geom.face_cent(face[0]));
   cyl_vertex v_start(geom.verts(face[0]));
   start_angle = v_start.angle;
   //fprintf(stderr, "\nstart_angle=%.17g\n", rad2deg(start_angle));

   int sign = 0;
   for(unsigned int i=0; i<face.size(); i++) {
      cyl_vertex cyl_v(geom.verts(face[i]));

      // convert angle into offset from start_angle, in range -M_PI<=ang<=M_PI
      //fprintf(stderr, "initial cyl_v.angle=%g (%.17g)\n", rad2deg(cyl_v.angle), cyl_v.angle/M_PI);
      cyl_v.angle -= start_angle;
      // make sure that if an offset is at 180 degrees then all offsets have
      // the same sign
      //fprintf(stderr, "offset cyl_v.angle=%g (%.17g)\n", rad2deg(cyl_v.angle), cyl_v.angle/M_PI);
      if(cyl_v.angle>M_PI)
         cyl_v.angle -= 2*M_PI;
      //if(double_le(cyl_v.angle, -M_PI, eps))
      if(cyl_v.angle < -M_PI)
         cyl_v.angle += 2*M_PI;
      //fprintf(stderr, "\tfinal cyl_v.angle=%g (%.17g)\n", rad2deg(cyl_v.angle), cyl_v.angle/M_PI);


      if(!sign && !double_eq(cyl_v.angle, 0, eps) &&
                  !double_eq(fabs(cyl_v.angle), M_PI, eps) ) {
         sign = 1-2*(cyl_v.angle<0);
         //fprintf(stderr, "sign set to %d\n", sign);
      }

      cyl_verts.push_back(cyl_v);
   }

   if(sign) {
      for(unsigned int i=0; i<cyl_verts.size(); i++) {
         cyl_vertex &cyl_v = cyl_verts[i];
         if(double_ge(fabs(cyl_v.angle), M_PI, eps) && cyl_v.angle*sign<0) 
         cyl_v.angle = sign*M_PI;
      }
   }
}


void cyl_face::add_to_geom(col_geom_v &geom, int from_n, int to_n, int to_d)
{
   double new_start_angle = start_angle * from_n * to_d/to_n;
   vector<int> face;
   int orig_v_sz = geom.verts().size();
   for(unsigned int i=0; i<cyl_verts.size(); i++) {
      cyl_vertex cyl_v = cyl_verts[i];
      cyl_v.angle = new_start_angle + cyl_v.angle*from_n*to_d/to_n;
      geom.add_vert(cyl_v.to_vec3d());
      face.push_back(orig_v_sz+i);
   }
   geom.add_face(face);
}

void cyl_face::add_chain_to_geom(col_geom_v &geom,
      int from_n, int to_n, int to_d)
{
   //fprintf(stderr, "chain.verts.size()=%d\n", cyl_verts.size());
   double chain_ang = cyl_verts.back().angle;

   double new_start_angle = start_angle * from_n * to_d/to_n;
   vector<int> face;
   for(int n=0; n<to_n; n++) {
   //for(int n=0; n<1; n++) {
      double ang_inc = n*chain_ang*from_n*to_d/to_n;
      for(unsigned int i=0; i<cyl_verts.size()-1; i++) {
         cyl_vertex cyl_v = cyl_verts[i];
         //fprintf(stderr, "\thorz cyl_v.angle=%g (%.17g)\n", rad2deg(cyl_v.angle), cyl_v.angle/M_PI);
         cyl_v.angle = ang_inc + new_start_angle + cyl_v.angle*from_n*to_d/to_n;
         face.push_back(geom.verts().size());
         geom.add_vert(cyl_v.to_vec3d());
      }
   }
   geom.add_face(face);
}

void get_horz_face_chains(vector<vector<int> > &chains,
      const vector<vector<int> > &horz_faces, const col_geom_v &geom,
      const map<vector<int>, vector<int> > &e2f,
      const vector<set<int> > &f_equivs)
{
   chains.clear();
   vector<int> f2type(geom.faces().size(), -1);
   for(unsigned int type=0; type<f_equivs.size(); type++) {
      if(f_equivs[type].size()>1) { // only cycling faces
         for(set<int>::const_iterator si=f_equivs[type].begin();
               si!=f_equivs[type].end(); ++si) {
            f2type[*si] = type;
            //fprintf(stderr, "face %d -> type %d\n", f2type[*si], type);
         }
      }
   }

   for(unsigned int i=0; i<horz_faces.size(); i++) {
      const vector<int> &face = horz_faces[i];
      const int f_sz = face.size();
      //fprintf(stderr, "horz_face[%u].size()=%d\n", i, f_sz);
      vector<int> chain;
      chain.push_back(face[0]);
      set<int> type_seen;
      for(int j=1; j<f_sz; j++) {
         vector<int> edge = make_edge(face[j-1], face[j]);
         map<vector<int>, vector<int> >::const_iterator mi =
            e2f.find(edge);
         //fprintf(stderr, "(%d, %d)\n", edge[0], edge[1]);
         if(mi!=e2f.end()) {
            for(unsigned int f=0; f<mi->second.size(); f++) {
               int f_idx = mi->second[f];
               //fprintf(stderr, "face idx = %d\n", f_idx);
               int type = (f_idx>=0) ? f2type[f_idx] : -1;
               if(type>=0 && type_seen.find(type)==type_seen.end()) {
                  //fprintf(stderr, "not seen\n");
                  chain.push_back(face[j]);
                  type_seen.insert(type);
               }
               else {
                  //fprintf(stderr, "%s\n", type>=0 ? "already seen" : "not considered");
               }

            }
         }
      }
      if(chain.size()!=face.size()) {
         fprintf(stderr, "horizontal polygon, full chain not found\n");
      }

      //fprintf(stderr, "chain: ");
      //for(unsigned int c=0; c<chain.size(); c++)
      //   fprintf(stderr, " %d", chain[c]);
      //fprintf(stderr, "\n");

      chains.push_back(chain);
   }
}


bool to_nfold(const col_geom_v &geom, col_geom_v &out_geom,
      int from_n, int to_n, int to_d, vector<string> &msgs)
{
   msgs.clear();
   sch_sym sym(sch_sym::C, from_n);
   vector<vector<set<int> > > sym_equivs;
   get_equiv_elems(geom, sym.get_trans(), &sym_equivs);
   const vector<set<int> > &f_equivs = sym_equivs[2];
   geom_info info(geom);

   col_geom_v face_geom;
   vector<vector<int> > horz_faces;
   for(unsigned int i=0; i<f_equivs.size(); i++) {
      const int f_idx = *f_equivs[i].begin();
      //fprintf(stderr, "\n\nNEW FACE face_idx=%d\n", f_idx);
      const vector<int> &face = geom.faces(f_idx);
      if(f_equivs[i].size()>1) {   // cycle of faces exists
         //fprintf(stderr, "cyclic face\n");
         cyl_face cyl_f(geom, face);
         cyl_f.add_to_geom(face_geom, from_n, to_n, to_d);
      }
      else {    // single face of its class
         int f_sz = face.size();
         double dir_sign = geom.face_norm(face).unit()[2];
         if(!double_eq(fabs(dir_sign), 1.0)) {
            msgs.push_back(msg_str(
                  "no faces generated corresponding to face index %d,"
                  "single vertical face through axis", f_idx) );
         }
         else if(face.size()%from_n) {    // error from coincident vertices?
            msgs.push_back(msg_str(
                  "no faces generated corresponding to face index %d,"
                  "axis is %d-fold but horizontal face has %d vertices",
                  f_idx, from_n, f_sz) );
         }
         else {
            horz_faces.push_back(face);
            horz_faces.back().resize(1+face.size()/from_n); // min edge chain
         }
      }
   }

   // infer "horizontal" axial digons
   map<vector<int>, vector<int> > e2f;
   geom.get_edge_face_pairs(e2f, false);
   map<vector<int>, vector<int> >::iterator mi;
   for(mi=e2f.begin(); mi!=e2f.end(); ++mi) {
      //fprintf(stderr, "(%d, %d) = ", mi->first[0], mi->first[1]);
      //for(unsigned int z=0; z<mi->second.size(); z++)
      //   fprintf(stderr, "%d, ", mi->second[z]);
      //fprintf(stderr, "\n");

      // ADD TEST FOR FACE NORMALS SUMMING TO ZERO
      vec3d c = geom.edge_cent(mi->first);
      if(is_even(mi->second.size()) && cyl_vertex(c).radius < epsilon) {
         horz_faces.push_back(mi->first);
      }
   }


   vector<vector<int> > chains;
   get_horz_face_chains(chains, horz_faces, geom, e2f, f_equivs);
   for(unsigned int i=0; i<chains.size(); i++) {
      cyl_face cyl_ch(geom, chains[i]);
      cyl_ch.add_chain_to_geom(face_geom, from_n, to_n, to_d);
   }

   //sym_repeat(out_geom, face_geom, sch_sym(sch_sym::C, 1));
   sym_repeat(out_geom, face_geom, sch_sym(sch_sym::C, to_n));
   sort_merge_elems(out_geom, "vef", epsilon);

   return true;
}



int main(int argc, char *argv[])
{
   nfold_opts opts;
   opts.process_command_line(argc, argv);

   col_geom_v geom;
   geom_read_or_error(geom, opts.ifile, opts);

   sch_sym sym(geom);

   const set<sch_axis> &axes = sym.get_axes();
   int n_fold = 0;
   set<sch_axis>::const_iterator ax;
   for(ax=axes.begin(); ax!=axes.end(); ++ax) {
      if(double_eq(ax->get_axis()[2], 1.0)) {   // aligned with z-axis
         n_fold = ax->get_nfold();
         if(ax->get_sym_type()==sch_sym::S)
            n_fold /= 2;
      }
   }

   if(!n_fold)
      opts.error("model does not have cyclic symmetry round the z-axis",
            "input");

   vector<string> msgs;
   col_geom_v o_geom;
   to_nfold(geom, o_geom, n_fold, opts.num, opts.denom, msgs);
   geom.orient();

   for(unsigned int i=0; i<msgs.size(); i++)
      opts.warning(msgs[i]);

   geom_write_or_error(o_geom, opts.ofile, opts);

   return 0;
}


