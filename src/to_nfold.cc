/*
   Copyright (c) 2013-2016, Adrian Rossiter

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
   Name: to_nfold.cc
   Description: Generalise an axial model by changing its rotational symmetry.
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
   int c;

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


// ----------------------------------------------------------------------

double to_range_0_2pi(double ang)
{
   return ang - 2*M_PI*floor(ang/(2*M_PI));
}

double to_range_plusminus_pi(double ang)
{
   double ret = to_range_0_2pi(ang);
   return (ret>M_PI) ? ret-2*M_PI : ret;
}


class cyc_vert
{
   public:
      double radius;
      double height;
      double cyc_ang;

      cyc_vert(double radius=0.0, double height=0.0, double angle=0.0) :
         radius(radius), height(height), cyc_ang(angle) {}

      cyc_vert(vec3d v)
      {
         radius = vec3d(v[0], v[1], 0).mag();
         height = v[2];
         cyc_ang = angle_around_axis(vec3d::X, v, vec3d::Z);
      }

      void to_range_0_2pi() { cyc_ang = ::to_range_0_2pi(cyc_ang); }
      void to_range_plusminus_pi() {cyc_ang = ::to_range_plusminus_pi(cyc_ang);}
      vec3d to_vec3d() const
         { return vec3d(radius*cos(cyc_ang), radius*sin(cyc_ang), height); }

      // convert absolute angle to relative angle, range (-M_PI,M_PI]
      void angle_abs_to_rel(double from)
         { cyc_ang -= from; to_range_plusminus_pi(); }
      bool on_axis(double eps) { return double_eq(radius, 0.0, eps); }

      void dump() const
      {
         vec3d v = to_vec3d();
         fprintf(stderr, "r=%g, h=%g, ang=%g (%g) -> (%g, %g, %g)\n",
               radius, height, cyc_ang, rad2deg(cyc_ang), v[0], v[1], v[2]);
      }

};

// ----------------------------------------------------------------------

struct cyc_vert_cmp {
   int n_fold;
   double eps;
   cyc_vert_cmp(int n_fold=1, double eps=epsilon): n_fold(n_fold), eps(eps) {}
   bool operator()(const cyc_vert &from, const cyc_vert &to) const
   {
      int ret = double_compare(from.height, to.height, eps);
      if(!ret) {
         ret = double_compare(from.radius, to.radius, eps);
         if(!ret && !double_eq(from.radius, 0.0, eps))
            ret = double_compare(from.cyc_ang, to.cyc_ang, eps);
      }
      return ret < 0;
   }
};


class cyc_verts
{
   private:
      int n_fold;
      double eps;
      map<cyc_vert, int, cyc_vert_cmp> v2type;

      int get_sector(double ang) const;
      double first_sector(double ang) const;
      double sector_width() const { return 2*M_PI/n_fold; }

   public:
      cyc_verts(int n_fold, double eps): n_fold(n_fold), eps(eps)
      { v2type = map<cyc_vert, int, cyc_vert_cmp>(cyc_vert_cmp(n_fold, eps));}

      pair<cyc_vert, int> get_vert_type(cyc_vert cyc_v) const;
      int get_vert_idx(cyc_vert cyc_v) const;
      void add_vert(cyc_vert cyc_v);
      void assign_types_in_order();
      void to_vector(vector<cyc_vert> *v_types) const;
};

pair<cyc_vert, int> cyc_verts::get_vert_type(cyc_vert cyc_v) const
{
   cyc_v.to_range_plusminus_pi();
   double ang_orig = cyc_v.cyc_ang;
   cyc_v.cyc_ang = first_sector(ang_orig);
   const map<cyc_vert, int, cyc_vert_cmp>::const_iterator mi =
      v2type.find(cyc_v);
   if(mi == v2type.end())
      return std::make_pair(cyc_v, -1);
   else
      return *mi;
}

int cyc_verts::get_vert_idx(cyc_vert cyc_v) const
{
   pair<cyc_vert, int> v = get_vert_type(cyc_v);
   if(v.second < 0)
      return -1;      // vertex not found
   int idx_off = floor((cyc_v.cyc_ang-v.first.cyc_ang)/sector_width()+0.5);
   idx_off = (2*n_fold + idx_off)%n_fold;
   if(idx_off<0)
      idx_off += n_fold;
   return (v.second*n_fold) + idx_off;
}


void cyc_verts::add_vert(cyc_vert cyc_v)
{
   double ang_orig = cyc_v.cyc_ang;
   cyc_v.cyc_ang = first_sector(ang_orig);
   pair<cyc_vert, int> v = get_vert_type(cyc_v);
   if(v.second < 0) {
      int sz = v2type.size();
      v2type[cyc_v] = sz;
   }
}

int cyc_verts::get_sector(double ang) const
{
   ang = to_range_0_2pi(ang);
   if(first_sector(ang) == 0.0)
      ang += 0.5*sector_width();   // sector boundary! Move to middle of sector
   return floor(ang/sector_width());
}

double cyc_verts::first_sector(double ang) const
{
   ang = to_range_0_2pi(ang);
   double first_ang = fmod(ang, sector_width());
   if(ang<0 || double_eq(first_ang, sector_width(), eps))
      first_ang = 0.0;
   return first_ang;
}

void cyc_verts::assign_types_in_order()
{
   int i=0;
   map<cyc_vert, int, cyc_vert_cmp>::iterator mi;
   for(mi=v2type.begin(); mi!=v2type.end(); ++mi)
      mi->second = i++;
}

void cyc_verts::to_vector(vector<cyc_vert> *v_types) const
{
   v_types->resize(v2type.size());
   map<cyc_vert, int, cyc_vert_cmp>::const_iterator mi;
   for(mi=v2type.begin(); mi!=v2type.end(); ++mi)
      (*v_types)[mi->second] = mi->first;
}


// ----------------------------------------------------------------------

enum { WARN_NONE=0, WARN_EDGE180=1, WARN_EDGE180_CONSECUTIVE=2,
       WARN_EDGE180_AMBIGUOUS=4, WARN_NOCHAIN=8, WARN_N_OVER_N=16,
       WARN_N_OVER_N2=32, WARN_VERTEX_NOT_FOUND=64};


class cyc_chain
{
   public:
      double start_angle;
      vector<cyc_vert> cyc_vts;
      unsigned char warnings;
      cyc_chain(const col_geom_v &geom, const vector<int> &face, int len,
            double eps);
};


cyc_chain::cyc_chain(const col_geom_v &geom, const vector<int> &face_orig,
                     int len, double eps): warnings(WARN_NONE)
{
   // ensure that face starts with a non-axial vertex
   vector<int> face = face_orig;
   int fsz = face.size();
   for(int i=0; i<fsz; i++) {
      cyc_vert cyc_v = geom.verts(face[i]);
      if(!cyc_v.on_axis(eps)) {
         if(i)
            rotate(face.begin(), face.begin()+i, face.end());
         break;
      }
   }

   double prev_abs_ang = cyc_vert(geom.verts(face.back())).cyc_ang;
   cyc_vts.resize(fsz);
   vector<int> half_turns;
   for(int i=0; i<=fsz; i++) {
      int idx = i%fsz;
      //fprintf(stderr, "face[%d/%d] = %d", idx, fsz, face[idx]);
      //geom.verts(face[idx]).dump();
      cyc_vert cyc_v(geom.verts(face[idx]));
      bool on_axis = cyc_v.on_axis(eps);
      double cur_abs_ang = on_axis ? prev_abs_ang: cyc_v.cyc_ang;
      if(i==0) {
         start_angle = cyc_v.cyc_ang;
         //fprintf(stderr, "\nstart_angle=%.17g\n", rad2deg(start_angle));
      }
      else {
         if(on_axis)
            cyc_v.cyc_ang = 0;
         else {
            cyc_v.angle_abs_to_rel(prev_abs_ang);
            if(fsz>2 && double_eq(fabs(cyc_v.cyc_ang), M_PI, eps))
               half_turns.push_back(idx);
         }
         cyc_vts[idx] = cyc_v;
      }

      prev_abs_ang = cur_abs_ang;
   }

   int hts_sz = half_turns.size();
   if(fsz>2 && hts_sz)
      warnings |= WARN_EDGE180;

   for(int i=0; i<hts_sz; i++) {
      int idx = half_turns[i];
      int idx_before = (idx - 1 + fsz)%fsz;
      int idx_after = (idx + 1)%fsz;
      if(hts_sz>1 && abs(half_turns[i]-half_turns[(i+1)%fsz])==1)
         warnings |= WARN_EDGE180_CONSECUTIVE;
      // Find sum of adjacent angles
      double test_ang = cyc_vts[idx_before].cyc_ang+cyc_vts[idx_after].cyc_ang;
      if(double_eq(test_ang, 0.0))
         warnings |= WARN_EDGE180_AMBIGUOUS;

      int sign = 1 - 2*(test_ang > 0);
      cyc_vts[idx].cyc_ang = sign*M_PI;
   }
   cyc_vts.resize(len);
}


// ----------------------------------------------------------------------


class cyc_geom
{
   private:
      double eps;
      const col_geom_v *geom;
      col_geom_v chains;
      int from_n;
      unsigned char warnings;
      string error_msg;

      int set_from_n(const sch_sym &sym);
      void add_chain_to_geom(cyc_chain &chain, col_geom_v *o_geom, int to_n,
            int to_d, const cyc_verts &cyc_vs);

   public:
      cyc_geom(double eps) : eps(eps) {}
      int set_geom(const col_geom_v *geo);
      int make_geom(col_geom_v *o_geom, int to_n, int to_d=1);
      const col_geom_v *get_geom() const { return geom; }
      const string &get_error_msg() const { return error_msg; }
      vector<string> get_input_warnings();
      vector<string> get_output_warnings();
};


int cyc_geom::set_geom(const col_geom_v *geo)
{
   from_n = 0;

   sch_sym sym(*geo);
   const set<sch_axis> &axes = sym.get_axes();
   set<sch_axis>::const_iterator ax;
   for(ax=axes.begin(); ax!=axes.end(); ++ax) {
      if(double_eq(ax->get_axis()[2], 1.0)) {   // aligned with z-axis
         from_n = ax->get_nfold();
         if(ax->get_sym_type()==sch_sym::S)
            from_n /= 2;
      }
   }

   if(!from_n) {
      geom = 0;
      error_msg = "model does not have cyclic symmetry around the z-axis";
      return false;
   }

   geom = geo;
   return true;
}

// Is edge is horizontal with centre on the axis
bool is_axial_edge(const vector<int> &edge, const geom_if &geom, double eps)
{
   //  (x, y, z) == (-x, -y, z)
   const vec3d &v0 = geom.verts(edge[0]);
   const vec3d &v1 = geom.verts(edge[1]);
   return double_eq(v0[0], -v1[0], eps) &&
          double_eq(v0[1], -v1[1], eps) &&
          double_eq(v0[2],  v1[2], eps);
}

// Create a rotated face by stepping through vertex indices
vector<int> get_next_face(vector<int> face, int step, int to_n)
{
   vector<int> f(face.size());
   for(unsigned int i=0; i<face.size(); i++) {
      f[i] = (face[i]/to_n)*to_n + (face[i]+step)%to_n;
   }
   return f;
}

vector<string> cyc_geom::get_input_warnings()
{
   vector<string> msgs;
   if(warnings & WARN_EDGE180)
      msgs.push_back("includes 180 degree edges");
   if(warnings & WARN_EDGE180_CONSECUTIVE)
      msgs.push_back("includes consecutive 180 degree edges");
   if(warnings & WARN_EDGE180_AMBIGUOUS)
      msgs.push_back("includes a 180 degree edge and its sign cannot be "
            "unambiguously assigned");
   if(warnings & WARN_NOCHAIN)
      msgs.push_back("includes an edge chain that does not repeat "
            "symmetrically (invalid result)");
   return msgs;
}

vector<string> cyc_geom::get_output_warnings()
{
   vector<string> msgs;
   if(warnings & WARN_N_OVER_N)
      msgs.push_back("includes a face with a vertex adjacent to itself "
            "(invalid N/N result?)");
   if(warnings & WARN_N_OVER_N2)
      msgs.push_back("includes a non-axial two-vertex face"
            "(invalid N/N result?)");
   if(warnings & WARN_VERTEX_NOT_FOUND)
      msgs.push_back("could not determine all vertex positions"
            "(invalid result)");
   return msgs;
}


void cyc_geom::add_chain_to_geom(cyc_chain &chain, col_geom_v *o_geom,
      int to_n, int to_d, const cyc_verts &cyc_vs)
{
   //fprintf(stderr, "\nADDING CHAIN\n");

   double first_sector_start_angle = fmod(chain.start_angle, 2*M_PI/from_n);
   double new_start_angle = first_sector_start_angle * from_n * to_d/to_n;
   vector<int> face;
   //fprintf(stderr, "\nstart_angle = %g\n", chain.start_angle);
   //fprintf(stderr, "new_start_angle = %g\n", new_start_angle);
   double prev_angle = new_start_angle;
   bool cycle_completed = false;
   bool all_vertices_valid = true;
   int num_chains_in_face = 0;
   for(int n=0; n<to_n; n++) {
      //fprintf(stderr, "n = %d\n", n);
      for(unsigned int i=0; i<chain.cyc_vts.size(); i++) {
         cyc_vert cyc_v = chain.cyc_vts[i];
         if(n==0 && i==0)
            cyc_v.cyc_ang = 0.0; // set offset from start angle, not last vert
         //fprintf(stderr, "i=%d, a = %g * pi\n", i, cyc_v.cyc_ang/(M_PI));
         cyc_v.cyc_ang = prev_angle + cyc_v.cyc_ang*from_n*to_d/to_n;
         //cyc_v.dump();
         prev_angle = cyc_v.cyc_ang;
         int idx = cyc_vs.get_vert_idx(cyc_v);
         //fprintf(stderr, "idx = %d\n\n", idx);
         if(i==0 && face.size() && idx==face[0]) {
            cycle_completed = true;
            break;
         }
         if(face.size() && face.back()==idx)
            warnings |= WARN_N_OVER_N;
         if(idx == -1) {
            warnings |= WARN_VERTEX_NOT_FOUND;
            all_vertices_valid = false;
         }
         face.push_back(idx);
      }
      if(cycle_completed)
         break;
      else
         num_chains_in_face++;
   }

   if(face.size()==2 && !is_axial_edge(face, *o_geom, eps))
      warnings |= WARN_N_OVER_N2;

   if(all_vertices_valid) {
      for(int i=0; i<to_n/num_chains_in_face; i++)
         o_geom->add_face(get_next_face(face, i, to_n));
   }
}
int cyc_geom::make_geom(col_geom_v *o_geom, int to_n, int to_d)
{
   warnings = WARN_NONE;
   sch_sym sym_Cn(sch_sym::C, from_n);
   vector<vector<set<int> > > sym_equivs;
   get_equiv_elems(*geom, sym_Cn.get_trans(), &sym_equivs);
   const vector<set<int> > &f_equivs = sym_equivs[2];

   // Set up vertices, as types
   cyc_verts cyc_vs(to_n, sym_eps);         // vertex types in final model.
   cyc_verts cyc_vs_from(from_n, sym_eps);  // vertex types in original model
   for(unsigned int i=0; i<sym_equivs[0].size(); i++) {
      // get equiv vertex in first from_n sector
      cyc_vert cyc_v = cyc_vs_from.get_vert_type(
            geom->verts(*sym_equivs[0][i].begin())).first;
      // 'scale' to corresponding position in first to_n sector
      cyc_v.cyc_ang *= (double)from_n * to_d / to_n;
      cyc_vs.add_vert(cyc_v);
   }
   cyc_vs.assign_types_in_order(); // number the vertices nicely

   // Set up final vertices
   vector<cyc_vert> vs;
   cyc_vs.to_vector(&vs);
   for(unsigned int i=0; i<vs.size(); i++) {
      cyc_vert v = vs[i];
      for(int step=0; step<to_n; step++) {
         cyc_vert v_rot = v;
         v_rot.cyc_ang += step * 2*M_PI/to_n;
         o_geom->add_vert(v_rot.to_vec3d());
      }
   }

   vector<cyc_chain> chains;
   for(unsigned int i=0; i<f_equivs.size(); i++) {
      const int f_idx = *f_equivs[i].begin();
      //fprintf(stderr, "\n\nNEW CHAIN face_idx=%d\n", f_idx);
      //fprintf(stderr, "Num faces=%d, from_n=%d\n", (int)f_equivs[i].size(), from_n);
      vector<int> face = geom->faces(f_idx);
      int fsz = face.size();

      // Find the minimum length of chain of edges that make a repeat unit
      int chain_sz = fsz * f_equivs[i].size() / from_n;
      if(chain_sz<fsz) {
         double steps = angle_around_axis(geom->verts(face[0]),
               geom->verts(face[chain_sz%fsz]), vec3d::Z) *
               from_n / (2*M_PI);
         double diff = steps - round(steps);
         if(double_ne(diff, 0.0, eps))
            warnings |= WARN_NOCHAIN;
      }

      cyc_chain cyc_c(*geom, face, chain_sz, eps);
      chains.push_back(cyc_c);
      warnings |= cyc_c.warnings;
   }

   // infer "horizontal" axial digons as a horizontal edge shared by two
   // equivalent faces
   map<int, int> f2equiv;
   for(unsigned int i=0; i<f_equivs.size(); i++)
      for(set<int>::iterator si = f_equivs[i].begin(); si != f_equivs[i].end();
            ++si)
         f2equiv[*si] = i;

   map<vector<int>, vector<int> > e2f;
   geom->get_edge_face_pairs(e2f, false);
   map<vector<int>, vector<int> >::iterator mi;
   for(mi=e2f.begin(); mi!=e2f.end(); ++mi) {
      if(is_axial_edge(mi->first, *geom, eps) &&
            mi->second.size()==2 &&
            f2equiv[mi->second[0]] == f2equiv[mi->second[1]]) {
         cyc_chain cyc_c(*geom, mi->first, 1, eps);
         chains.push_back(cyc_c);
      }
   }

   for(vector<cyc_chain>::iterator vi=chains.begin(); vi!=chains.end(); ++vi)
      add_chain_to_geom(*vi, o_geom, to_n, to_d, cyc_vs);
   return true;
}

// ----------------------------------------------------------------------

int main(int argc, char *argv[])
{
   const double eps=10e-8; // hardcoded
   nfold_opts opts;
   opts.process_command_line(argc, argv);

   col_geom_v geom;
   geom_read_or_error(geom, opts.ifile, opts);
   sort_merge_elems(geom, "vef", 0, eps);

   cyc_geom cyc(eps);
   if(!cyc.set_geom(&geom))
      opts.error(cyc.get_error_msg(), "input geometry");

   col_geom_v o_geom;
   cyc.make_geom(&o_geom, opts.num, opts.denom);

   o_geom.orient();

   vector<string> msgs = cyc.get_input_warnings();
   for(unsigned int i=0; i<msgs.size(); i++)
      opts.warning(msgs[i], "input geometry");

   msgs = cyc.get_output_warnings();
   for(unsigned int i=0; i<msgs.size(); i++)
      opts.warning(msgs[i], "output geometry");

   sort_merge_elems(o_geom, "vef", 0, eps);
   geom_write_or_error(o_geom, opts.ofile, opts);

   return 0;
}



