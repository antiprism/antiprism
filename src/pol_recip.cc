/*
   Copyright (c) 2003-2009, Adrian Rossiter

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
   Name: polar_recip.cc
   Description: make a polar reciprocal
   Project: Antiprism - http://www.antiprism.com
*/

#include <string.h>
#include <unistd.h>
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


class pr_opts: public prog_opts
{
   private:
   
   public:
      vector<int> space_verts;
      double recip_rad;
      double init_rad;
      char recip_rad_type;
      vec3d centre;
      vec3d init_cent;
      char recip_cent_type;
      char init_cent_type;
      double inf;
      double extra_ideal_elems;
      int num_iters;
      int lim_exp;
      bool append;

      string ifile;
      string ofile;

      pr_opts(): prog_opts("pol_recip"),
                 recip_rad(0), init_rad(0), recip_rad_type('x'),
                 recip_cent_type('x'), init_cent_type('x'),
                 inf(1e15), extra_ideal_elems(true),
                 num_iters(100), lim_exp(13), append(false)
                 {}

      void process_command_line(int argc, char **argv);
      void usage();
};

void pr_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Read a file in OFF format and make a polar reciprocal from the face\n"
"planes. If input_file is not given the program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"  -c <cent> reciprocation centre (default: C)\n"
"              X,Y,Z - centre with these coordinates\n"
"              C - vertex centroid,    M - mid-sphere (approximate)\n"
"              e - edge balanced,      E - edge balanced with inversion\n"
"              v - vert/face balanced, V - vert/face balanced with inversion\n"
"              a - v/f/e balanced      A - v/f/e balanced with inversion\n"
"  -C <init> initial value for a centre calculation, in form 'X,Y,Z',\n"
"            or C to use centroid (default), M to calculate approx mid-sphere\n"
"  -r <rad>  reciprocation radius (default: calculated)\n"
"              radius value - use this value as the radius\n"
"              v - nearest vertex distance, V - furthest vertex distance\n"
"              e - nearest edge distance,   E - furthest edge distance\n"
"              f - nearest face distance,   F - furthest face distance\n"
"              X - topological dual using face centroids for dual vertices\n"  
"            or a comma separated list of vertex indices (starting from 0)\n"
"            and the distance is to the space containing those vertices\n"
"  -R <rad>  initial value for a radius calculation (default: calculated)\n"
"  -I <dist> maximum distance to project any normal or infinite dual vertex\n"
"            (default: 1e15), if 0 then use actual distances and delete\n"
"            infinite points\n"
"  -x        exclude extra elements added to duals with ideal vertices\n"
"  -n <itrs> maximum number of iterations (default: 10000)\n"
"  -l <lim>  minimum distance change to terminate, as negative\n"
"            exponent 1e-lim (default: 13 giving 1e-13)\n"
"  -a        append dual to original polyhedron\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}



     
void pr_opts::process_command_line(int argc, char **argv)
{
   char errmsg[MSG_SZ];
   extern char *optarg;
   extern int optind, opterr;
   opterr = 0;
   char c;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hc:C:r:R:xao:n:l:I:")) != -1 ) {
      if(common_opts(c))
         continue;

      switch(c) {
         case 'C':
            if(strlen(optarg)==1 && strchr("CM", *optarg))
               init_cent_type = *optarg;
            else if(init_cent.read(optarg, errmsg))
               init_cent_type = 'c';
            else
               error("initial centre must be three coordinates, C, or M", c);
            break;
            
         case 'c':
            if(strlen(optarg)==1 && strchr("CMeEvVaA", *optarg))
               recip_cent_type = *optarg;
            else if(centre.read(optarg, errmsg)) 
               recip_cent_type = 'c';
            else
               error("centre type must be three coordinates, C,M,e,E,v,V,a, or A", c);
            break;

         case 'R':
            if(!read_double(optarg, &init_rad, errmsg))
               error(errmsg, c);
            if(init_rad < epsilon)
               error("radius cannot be zero (or too close) or negative", c);
            break;

         case 'r':
            if(strlen(optarg)==1 && strchr("VvEeFfX", *optarg)) {
               recip_rad_type = *optarg;
            }
            else if(read_double(optarg, &recip_rad, errmsg)) {
               if(recip_rad < epsilon)
                  error("radius cannot be zero (or too close) or negative", c);
               recip_rad_type = 'r';
            }
            else if(read_int_list(optarg, space_verts, errmsg, true))
               recip_rad_type='S';
            else
               error("radius must be a radius value, or r,V,v,E,e,F,f,X or a list of index numbers", c);
            break;

         case 'x':
            extra_ideal_elems = false;
            break;
         
         case 'a':
            append = true;
            break;
         
         case 'o':
            ofile = optarg;
            break;
         
         case 'I':
            if(!read_double(optarg, &inf, errmsg))
               error(errmsg, c);
            if(inf < 0.0)
               error("distance must be positive, or 0 to disable", c);
            break;
         
         case 'n':
            if(!read_int(optarg, &num_iters, errmsg))
               error(errmsg, c);
            if(num_iters < 0)
               error("number of iterations must be greater than 0", c);
            break;

         case 'l':
            if(!read_int(optarg, &lim_exp, errmsg))
               error(errmsg, c);
            if(lim_exp < 0) {
               warning("limit is negative, the exponent used will be positive", c);
            }
            if(lim_exp > 16) {
               warning("limit is very small, may not be attainable", c);
            }
            break;

         default:
            error("unknown command line error");
      }
   }

   if(argc-optind > 1) {
      error("too many arguments");
      exit(1);
   }
   
   if(argc-optind == 1)
      ifile=argv[optind];

}

int find_mid_centre(geom_if &geom, double &rad, vec3d &cent, int n, double lim)
{
   vector<vector<int> > g_edges;
   geom.get_impl_edges(g_edges);
   int e_sz = g_edges.size();
   if(!cent.is_set())
      cent = centroid(*geom.get_verts());
   if(fabs(rad)<epsilon) {
      geom_info rep(geom);
      rep.set_center(cent);
      rad = rep.impl_edge_dists().sum/e_sz;
   }
   vec3d cur_cent = cent;
   double cur_rad = rad;
   double cent_test=-1, rad_test=-1;
   int cnt;
   for(cnt=0; cnt<n; cnt++) {
      vec3d c_diff(0,0,0);
      double rad_sum_g=0;
      for(int e=0; e<e_sz; ++e) {
         vec3d Pg = nearest_point(cur_cent, *geom.get_verts(), g_edges[e]);
         double r = (Pg-cur_cent).mag();
         rad_sum_g += r;
         c_diff += (1-r/cur_rad)*(Pg-cur_cent);
      }
      c_diff = c_diff/(double)e_sz;
      cur_cent = cent;
      cent = cur_cent - c_diff/(double)e_sz;
      cur_rad = rad;
      rad = rad_sum_g/e_sz;
      cent_test = (cur_cent - cent).mag()/rad;
      rad_test = (cur_rad - rad)/rad;
      //cent.dump("cent");
      //fprintf(stderr, "rad=%g, delt_rad=%g,delta_cent=%g\n", rad, cur_rad-rad,
      //      (cur_cent-cent).mag());
      if(fabs(cent_test)<lim && fabs(rad_test)<lim)
         break;
   }
   char str[MSG_SZ];
   fprintf(stderr, "[n=%d, limit=%sachieved, r_test=%g, c_test=%g]\n"
         "centre=(%s), radius=%.16g\n",
         cnt, cnt==n?"not ":"", cent_test, rad_test,
         vtostr(str, cent, " "), rad);
   return 1;
}
 
int find_can_centre(geom_if &geom, char type, double &rad, vec3d &cent, int n, double lim)
{
   bool use_e = false;
   bool use_vf = false;
   if(type=='e')
      use_e = true;
   else if(type=='v')
      use_vf = true;
   else {
      use_e = true;
      use_vf = true;
   }
      
   vector<vector<int> > g_edges, d_edges;
   geom.get_impl_edges(g_edges);
   int e_sz = g_edges.size();
   if(!cent.is_set())
      cent = centroid(*geom.get_verts());
   if(fabs(rad)<epsilon) {
      geom_info rep(geom);
      rep.set_center(cent);
      rad = (1-2*(rad<0))*rep.impl_edge_dists().sum/e_sz;
   }
   geom_v dual;
   get_dual(geom, dual, rad, cent);
   dual.get_impl_edges(d_edges);
   vec3d cur_cent = cent;
   double cur_rad = rad;
   double cent_test=-1, rad_test=-1;
   int cnt;
   for(cnt=0; cnt<n; cnt++) {
      vec3d cent_calc(0,0,0);
      double rad_g=0, rad_d=0;
      if(use_e) {
         vec3d e_sum_g(0,0,0);
         vec3d e_sum_d(0,0,0);
         double rad_sum_g=0, rad_sum_d=0;
         for(int e=0; e<e_sz; ++e) {
            vec3d Pg = nearest_point(cur_cent, *geom.get_verts(), g_edges[e]);
            rad_sum_g += (Pg-cur_cent).mag();
            e_sum_g += Pg;
            vec3d Pd = nearest_point(cur_cent, *dual.get_verts(), d_edges[e]);
            rad_sum_d += (Pd-cur_cent).mag();
            e_sum_d += Pd;
         }
         cent_calc += (e_sum_g+e_sum_d) / (2.0*e_sz);
         rad_g += rad_sum_g;
         rad_d += rad_sum_d;
      }
      if(use_vf) {
         vector<vec3d> &g_verts = *geom.get_verts();
         vector<vec3d> &d_verts = *dual.get_verts();
         int v_sz = g_verts.size();
         vector<vector<int> > &g_faces = *geom.get_faces();
         vector<vector<int> > &d_faces = *dual.get_faces();
         int f_sz = g_faces.size();
         vec3d f_sum_g(0,0,0);
         vec3d v_sum_d(0,0,0);
         vec3d v_sum_g(0,0,0);
         vec3d f_sum_d(0,0,0);
         double v_rad_sum_g=0, v_rad_sum_d=0;
         double f_rad_sum_g=0, f_rad_sum_d=0;
         for(int f=0; f<f_sz; ++f) {
            vec3d Pg = nearest_point(cur_cent, *geom.get_verts(), g_faces[f]);
            f_rad_sum_g += (Pg-cur_cent).mag();
            f_sum_g += Pg;
            vec3d Pd = d_verts[f];
            v_rad_sum_d += (Pd-cur_cent).mag();
            v_sum_d += Pd;
         }
      
         for(int v=0; v<v_sz; ++v) {
            vec3d Pg = g_verts[v];
            v_rad_sum_g += (Pg-cur_cent).mag();
            v_sum_g += Pg;
            vec3d Pd = nearest_point(cur_cent, *dual.get_verts(), d_faces[v]);
            f_rad_sum_d += (Pd-cur_cent).mag();
            f_sum_d += Pd;
         }
      
         vec3d vf_avg_g = (v_sum_g + f_sum_g)/double(v_sz+f_sz);
         vec3d vf_avg_d = (v_sum_d + f_sum_d)/double(v_sz+f_sz);
         vec3d vf_cent = (vf_avg_g + vf_avg_d)/2.0;
         cent_calc += vf_cent;
         rad_g += v_rad_sum_g * f_rad_sum_g;
         rad_d += v_rad_sum_d * f_rad_sum_d;
      }
      if(use_e && use_vf) {
         cent_calc /= 2.0;
      }
 
         
      if(is_even(cnt)) {
         cur_rad = rad;
         //rad = sqrt(rad_g/rad_d)*cur_rad;
         rad = 0.5*cur_rad + 0.5*sqrt(rad_g/rad_d)*cur_rad;
      }
      else {
         cur_cent = cent;
         cent = 0.5*cur_cent + 0.5*cent_calc;
      }
      
      cent_test = (cur_cent - cent).mag()/rad;
      rad_test = (cur_rad - rad)/rad;
      //cent.dump("cent");
      //fprintf(stderr, "rad=%g, delt_rad=%g,delta_cent=%g\n", rad, cur_rad-rad,
      //      (cur_cent-cent).mag());
      if(fabs(cent_test)<lim && fabs(rad_test)<lim)
         break;
      get_pol_recip_verts(geom, dual, rad, cur_cent);
   }
   char str[MSG_SZ];
   fprintf(stderr, "[n=%d, limit=%sachieved, r_test=%g, c_test=%g]\n"
         "centre=(%s), radius=%.16g\n",
         cnt, cnt==n?"not ":"", cent_test, rad_test,
         vtostr(str, cent, " "), rad);
   if(rad<0) {
      mat3d inv = mat3d::transl(cent) *
                  mat3d::inversion() *
                  mat3d::transl(-cent);
      geom.transform(inv);
   }

return 1;
}


void find_recip_centre(geom_if &geom, char type, double &rad, vec3d &cent, int n, double lim)
{
   if(fabs(rad)<epsilon)
      rad = epsilon/2.0;
   switch(type) {
      case 'x': // fefault is C
      case 'C':
         cent = centroid(*geom.get_verts());
         break;
      case 'M':
         find_mid_centre(geom, rad, cent, n, lim);
         break;
      case 'e':
         find_can_centre(geom, 'e', rad, cent, n, lim);
         break;
      case 'E':
         rad = -rad;
         find_can_centre(geom, 'e', rad, cent, n, lim);
         break;
      case 'v':
         find_can_centre(geom, 'v', rad, cent, n, lim);
         break;
      case 'V':
         rad = -rad;
         find_can_centre(geom, 'v', rad, cent, n, lim);
         break;
      case 'a':
         find_can_centre(geom, 'a', rad, cent, n, lim);
         break;
      case 'A':
         rad = -rad;
         find_can_centre(geom, 'a', rad, cent, n, lim);
         break;
   }   
}

double find_recip_rad(geom_if &geom, char type, vec3d cent,
      vector<int> vidxs=vector<int>())
{
   double rad = 1;
   geom_info rep(geom);
   rep.set_center(cent);
   switch(type) {
      case 'v':
         rad = rep.vert_dists().min;
         break;
      case 'V':
         rad = rep.vert_dists().max;
         break;
      case 'e':
         rad = rep.impl_edge_dists().min;
         break;
      case 'E':
         rad = rep.impl_edge_dists().max;
         break;
      case 'f':
         rad = rep.face_dists().min;
         break;
      case 'F':
         rad = rep.face_dists().max;
         break;
      case 'S':
         rad = (nearest_point(cent, *geom.get_verts(), vidxs) - cent).mag();
         break;
      case 'X':
         rad = 1;
         break;
   }
   return rad;
}


int is_polyhedron(geom_if &geom, char *errmsg=0)
{
   map<pair<int, int>, int> edge_cnts;
   map<pair<int, int>, int>::iterator ei;
   vector<vector<int> > &faces = *geom.get_faces();
   if(!faces.size()) {
      if(errmsg)
         strcpy(errmsg, "not a polyhedron, has no faces");
      return 0;
   }
   pair<int, int> edge;
   for(unsigned int i=0; i<faces.size(); ++i) {
      for(unsigned int j=0; j<faces[i].size(); ++j) {
         edge.first=faces[i][j];
         edge.second=faces[i][(j+1)%faces[i].size()];
         if(edge.first > edge.second)
            swap(edge.first, edge.second);
         ei = edge_cnts.find(edge);
         if(ei==edge_cnts.end())
            edge_cnts[edge] = 1;
         else
            ei->second++;
      }
   }

   for(ei=edge_cnts.begin(); ei!=edge_cnts.end(); ++ei) {
      if(ei->second!=2) {
         if(errmsg)
            sprintf(errmsg, "not a polyhedron, has an edge joined by"
                  " %d face(s)", ei->second);
         return ei->second;
      }
   }
   return 2;                        // each edge shared by exactly 2 faces
}


int main(int argc, char *argv[])
{
   pr_opts opts;
   opts.process_command_line(argc, argv);

   char errmsg[MSG_SZ];
   col_geom_v geom;
   if(!geom.read(opts.ifile, errmsg))
      opts.error(errmsg);
   if(*errmsg)
      opts.warning(errmsg);

   int e_cons = is_polyhedron(geom, errmsg);
   if(e_cons!=2)
      opts.error(errmsg, "input_file");

   double lim = pow(10, -opts.lim_exp);
   // Set up init_centre
   
   double tmp_rad = 0;
   if(opts.init_cent_type=='C')
      opts.init_cent = centroid(*geom.get_verts());
   else if(opts.init_cent_type=='M')
      find_mid_centre(geom, tmp_rad, opts.init_cent, opts.num_iters, lim);

   if(strchr("CMx", opts.recip_cent_type) && opts.recip_rad_type=='x') 
      opts.recip_rad_type='e';
   

   vec3d centre;
   if(opts.recip_cent_type=='c') {
      centre = opts.centre;
      opts.init_rad = 1.0;
   }
   else {
      centre = opts.init_cent;
      find_recip_centre(geom, opts.recip_cent_type, opts.init_rad, centre,
            opts.num_iters, lim);
   }
   
   double radius;
   if(opts.recip_rad_type=='r')
      radius = opts.recip_rad;
   else if(opts.recip_rad_type=='x')
      radius = opts.init_rad;
   else
      radius = find_recip_rad(geom, opts.recip_rad_type, centre,
            opts.space_verts);

   //fprintf(stderr, "radius=%g\n", radius);
   //centre.dump("cent");
   col_geom_v dual;
   get_dual(geom, dual, radius, centre, 1.01*opts.inf);

   
   vector<int> invalid_verts;
   for(unsigned int i=0; i<dual.get_verts()->size(); i++) {
      if(!(*dual.get_verts())[i].is_set()) {
         dual.delete_vert(i);
         int idx = invalid_verts.size()+i;
         invalid_verts.push_back(idx);
         i--;
      }
   }
   if(invalid_verts.size()) {
      string msg("removed invalid vertices (and associated faces) with indices - ");
      for(unsigned int i=0; i<invalid_verts.size()-1; i++) {
         snprintf(errmsg, MSG_SZ, "%d,", invalid_verts[i]);
         msg += string(errmsg);
      }
      snprintf(errmsg, MSG_SZ, "%d", invalid_verts.back());
      msg += string(errmsg);
      opts.warning(msg);
   }

   if(opts.recip_rad_type=='X')
      geom.face_cents(dual.raw_verts());

   if(opts.extra_ideal_elems)
      add_extra_ideal_elems(dual, centre, 1.005*opts.inf);

   dual.orient();
   if(opts.append)
      geom.append(dual);

   const col_geom_v &geom_out = (opts.append) ? geom : dual;
   if(!geom_out.write(opts.ofile, errmsg))
      opts.error(errmsg);

   return 0;
}
   

