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
   Name: spidron.cc
   Description: convert polyhedron to one built from spidron units
   Project: Antiprism - http://www.antiprism.com
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <string>
#include <vector>
#include <algorithm>

#include "../base/antiprism.h"

using std::string;
using std::vector;
using std::map;
using std::swap;

class spid_opts : public prog_opts {
   public:
      
      int unit_len;
      int type;
      double ang;
      double ang2;
      string pattern;
      string ifile;
      string ofile;

      spid_opts(): prog_opts("spidron"),
                   unit_len(5),
                   type(1),
                   ang(M_PI/6),
                   ang2(0),
                   pattern("01")
                   {}
      void process_command_line(int argc, char **argv);
      void usage();
};


void spid_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Read a file in OFF format and convert to model built from spidron units.\n"
"If input_file is not given the program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"  -l <len>  unit length, half the number of triangles from edge to centre\n"
"  -t <type> folding type, 1 - spidronmyd (default)\n"
"  -a <ang>  angle of first triangle (default: 30)\n"
"  -b <ang>  angle of second triangle (default: two times angle of option -a)\n"
"  -p <pat>  folding pattern, a series of 0's and 1's, to old in or out\n"
"            (default: 01)\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}


void spid_opts::process_command_line(int argc, char **argv)
{
   char errmsg[MSG_SZ];
   opterr = 0;
   char c;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":ht:a:b:l:p:o:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         case 'o':
            ofile = optarg;
            break;

         case 't':
            if(!read_int(optarg, &type, errmsg))
               error(errmsg, c);
            if(type<1 || type>1)
               error("type must be 1", c);
            break;

         case 'a':
            if(!read_double(optarg, &ang, errmsg))
               error(errmsg, c);
            ang *= M_PI/180;
            break;

         case 'b':
            if(!read_double(optarg, &ang2, errmsg))
               error(errmsg, c);
            if(ang2==0.0)
               error("angle cannot be 0.0", c);
            ang2 *= M_PI/180;
            break;

         case 'l':
            if(!read_int(optarg, &unit_len, errmsg))
               error(errmsg, c);
            if(unit_len<1)
               error("unit length must be 1 or greater", c);
            break;

         case 'p':
            if(strspn(optarg, "01")!=strlen(optarg))
               error("pattern can contain only 0 and 1", c);
            pattern = optarg;
            break;
            
         default:
            error("unknown command line error");
      }
   }

   if(ang2==0.0)
      ang2 = 2*ang;

   if(argc-optind > 1)
      error("too many arguments");
   
   if(argc-optind == 1)
      ifile=argv[optind];

}

int add_spidron_unit(geom_if &spid, double ang, double ang2, vec3d &cent, bool fold, col_val col=col_val(0))
{
   int v0idx = spid.get_verts()->size()-2;
   vec3d v0 = (*spid.get_verts())[v0idx];
   vec3d v1 = (*spid.get_verts())[v0idx+1];
   vec3d mid = (v0 + v1)/2.0;               // Midpoint of edge
   double R = (cent-v0).mag();              // Radius of polygon
   double L = (mid-v0).mag();               // Half edge length 
   double H = (cent-mid).mag();             // Height from edge to centre
   double theta = acos(safe_for_trig(H/R)); // Half angle of edge at centre
   double l = cos(ang2)*L/cos(ang);         // Half inner edge length
   //double l = L/(2*cos(ang));             // Half inner edge length
   double r = l/sin(theta);                 // Radius of inner edge
   double h_tri = L*tan(ang);               // Height of base triangle
   double delta_x = H-r;
   double delta_y = (1.0-2*fold)*sqrt(h_tri*h_tri - delta_x*delta_x);
   vec3d norm = vcross(v1-v0, cent-v0).unit();
   vec3d perp = (cent-mid).unit();
   vec3d P = mid + perp*delta_x + norm*delta_y;
   mat3d trans = mat3d::transl(cent) *
                 mat3d::rot(norm, theta*2) *
                 mat3d::transl(-cent);
   vec3d Q = trans*P;
   
   /*
   fprintf(stderr, "delta_x=%g, delta_y=%g, theta=%g, H=%g, l=%g\n", delta_x, delta_y, theta, H, L);
   cent.dump("cent");
   norm.dump("norm");
   perp.dump("perp");
   mid.dump("mid");
   v0.dump("v0");
   v1.dump("v1");
   P.dump("P");
   Q.dump("Q");
   */
   
   col_geom *cg_spid = dynamic_cast<col_geom *>(&spid);
   vector<vector<int> > &faces = *spid.get_faces();
   spid.add_vert(P);
   spid.add_vert(Q);
   int f = spid.add_face(vector<int>());
   if(cg_spid)
      cg_spid->set_f_col(f, col);
   faces[f].push_back(v0idx);
   faces[f].push_back(v0idx+1);
   faces[f].push_back(v0idx+2);
   f = spid.add_face(vector<int>());
   if(cg_spid)
      cg_spid->set_f_col(f, col);
   faces[f].push_back(v0idx+1);
   faces[f].push_back(v0idx+3);
   faces[f].push_back(v0idx+2);
   cent += norm*delta_y;
   return 1;
}
      

int add_spidron_unit2(geom_if &spid, double ang, double ang2, vec3d &cent, bool fold, col_val col=col_val(0))
{
   int v0idx = spid.get_verts()->size()-2;
   vec3d v0 = (*spid.get_verts())[v0idx];
   vec3d v1 = (*spid.get_verts())[v0idx+1];
   vec3d norm = vcross(v1-v0, cent-v0).unit();
   vec3d mid = (v0 + v1)/2.0;               // Midpoint of edge
   double R = (cent-v0).mag();              // Radius of polygon
   double L = (mid-v0).mag();               // Half edge length 
   double H = (cent-mid).mag();             // Hieght from edge to centre
   double theta = acos(safe_for_trig(H/R)); // Half angle of edge at centre
  
   double h = sin(ang2)*L/cos(ang);         // 2nd triangle height
   double mid2r = R - h;                    // radius to mid inner edgd
   //double r = mid2r/cos(theta);           // radius to inner edge
   double l = cos(ang2)*L/cos(ang);         // Half inner edge length
   double delta_x = mid2r*tan(theta);       // length proj of half inner edge
   double delta_y = (1.0-2*fold)*sqrt(l*l - delta_x*delta_x);
   //double rot_ang = acos(lp/l);           // angle to rot 2nd triangle
   
   vec3d mid2 = v1 + h*(cent-v1).unit();    // mid-point inner edge
   vec3d perp = vcross(cent-v1, norm).unit();
   vec3d P = mid2 - perp*delta_x + norm*delta_y;
   vec3d Q = mid2 + perp*delta_x - norm*delta_y;

   perp = (cent - mid).unit();
   delta_x = L*tan(ang);              // Height of base triangle
   delta_y=0;
   //   vec3d P = mid - perp*delta_x + norm*delta_y;
   //   vec3d Q = mid + perp*delta_x - norm*delta_y;
   
   fprintf(stderr, "\ndelta_x=%g, delta_y=%g, theta=%g, H=%g, l=%g\n", delta_x, delta_y, theta, H, l);
   cent.dump("cent");
   norm.dump("norm");
   perp.dump("perp");
   mid.dump("mid");
   mid2.dump("mid2");
   v0.dump("v0");
   v1.dump("v1");
   P.dump("P");
   Q.dump("Q");
   
   col_geom *cg_spid = dynamic_cast<col_geom *>(&spid);
   vector<vector<int> > &faces = *spid.get_faces();
   spid.add_vert(P);
   spid.add_vert(Q);
   int f = spid.add_face(vector<int>());
   if(cg_spid)
      cg_spid->set_f_col(f, col);
   faces[f].push_back(v0idx);
   faces[f].push_back(v0idx+1);
   faces[f].push_back(v0idx+2);
   f = spid.add_face(vector<int>());
   if(cg_spid)
      cg_spid->set_f_col(f, col);
   faces[f].push_back(v0idx+1);
   faces[f].push_back(v0idx+3);
   faces[f].push_back(v0idx+2);
   
   mat3d trans = mat3d::transl(cent) *
                 mat3d::rot(norm, theta*4) *
                 mat3d::transl(-cent);
   vec3d P2 = trans*P;
   vec3d V0 = (P-Q).unit();
   vec3d V1 = (P2-Q).unit();
   double beta = acos(safe_for_trig(vdot(V0, V1)));
   double slant = l/cos(beta/2);
   vec3d rad_dir = ((P+P2)/2.0 - Q).unit();
   cent = Q + rad_dir*slant;
   fprintf(stderr, "\nbeta=%g, slant=%g\n", beta*180/M_PI, slant);
   cent.dump("cent");
   P2.dump("P2");
   //spid.add_vert(cent);
   return 1;
}
      


int make_spidron(geom_if &spid, geom_if &base, double ang, double ang2, int len, vector<bool> folds, int type)
{
   vector<vec3d> &verts = *base.get_verts();
   vector<vector<int> > &edges = *base.get_edges();
   vector<vector<int> > &faces = *base.get_faces();
   col_geom *cg_spid = dynamic_cast<col_geom *>(&spid);
   col_geom *cg_base = dynamic_cast<col_geom *>(&base);
   for(unsigned int i=0; i<faces.size(); i++) {
      vec3d cent = centroid(verts, faces[i]);
      if(faces[i].size()<3)
         continue;
      for(unsigned int j=0; j<faces[i].size(); j++) {
         int v0 = faces[i][j];
         int v1 = faces[i][(j+1)%faces[i].size()];
         spid.add_vert(verts[v0]);
         spid.add_vert(verts[v1]);
         vec3d edge_cent = cent;
         col_val col(0,0,0,0);
         if(cg_spid && cg_base) {
            vector<int> edge(2);
            edge[0] = v0;
            edge[1] = v1;
            if(edge[0] > edge[1])
               swap(edge[0], edge[1]);
            vector<vector<int> >::iterator ei;
            ei = find(edges.begin(), edges.end(), edge);
            if(ei != edges.end()) {
               unsigned int e_idx = ei - edges.begin();
               col = cg_base->get_e_col(e_idx);
            }
         }
         for(int l=0; l<len; l++) {
            if(type==1)
               add_spidron_unit(spid, ang, ang2, edge_cent,
                     folds[l%folds.size()], col);
            if(type==2) {
               bool fold = folds[l%folds.size()];
               if(!is_even(j))
                  fold = !fold;
               add_spidron_unit2(spid, ang, ang2, edge_cent, fold, col);
            }
         }
      }
   }
   return 1;
}
   

double v_ang_at_ax(const vec3d &v0, const vec3d &v1, const vec3d &ax)
{
   vec3d n0 = vcross(v0, ax).unit();
   vec3d n1 = vcross(v1, ax).unit();
   double ang = acos(safe_for_trig(vdot(n0, n1)));
   if(vdot(ax, vcross(n0, n1))<0)
      ang = 2*M_PI - ang;
   return ang;
}



int main(int argc, char *argv[])
{
   spid_opts opts;
   opts.process_command_line(argc, argv);
   char errmsg[MSG_SZ];
   col_geom_v geom;
   if(!geom.read(opts.ifile, errmsg))
      opts.error(errmsg);
   if(*errmsg)
      opts.warning(errmsg);

   geom.add_missing_impl_edges();
   
   vector<bool> fold(opts.pattern.size());
   for(unsigned int f=0; f<opts.pattern.size(); f++)
      fold[f] = opts.pattern[f]=='0';
   
   col_geom_v spid;
   make_spidron(spid, geom, opts.ang, opts.ang2, opts.unit_len, fold, opts.type);

   if(!spid.write(opts.ofile, errmsg))
      opts.error(errmsg);
   
   return 0;
}


