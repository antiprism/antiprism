/*
   Copyright (c) 2006-2009, Adrian Rossiter

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
   Name: jitterbug.cc
   Description: jitterbug transformation
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
#include "../base/vec_utils.h"

using std::string;
using std::vector;


enum face_types{f_none = 0, f_equ = 1, f_oth = 2, f_quad = 4};

class jb_opts: public prog_opts
{
   private:
   
   public:
      double stage;
      unsigned int faces;
      char cycle_type;
      int lat_size;
      bool tri_fix;
      bool tri_top;
      bool fix_connect;

      string ofile;

      jb_opts(): prog_opts("jitterbug"), stage(0.0), faces(f_equ|f_oth),
                 cycle_type('f'), lat_size(1), tri_fix(false),
                 tri_top(false), fix_connect(false)
                 {}

      void process_command_line(int argc, char **argv);
      void usage();
};

void jb_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [cycle_stage]\n"
"\n"
"Create a stage of the jitterbug transformation in OFF format. The\n"
"cycles run from 0.0 to 1.0. The fractional part of cycle_stage\n"
"is tied to a stage in a full rotation jitterbug. This stage is\n"
"measured by the distance travelled by a point moving around a square\n"
"which halves a cube of edge 0.25, and which determines the jitterbug\n"
"shape.\n"
"\n"
"Options\n"
"%s"
"  -f <fces> faces to include:"
"              x - none\n"
"              a - all\n"
"              e - equilateral triangles (default)\n"
"              o - other triangles\n"
"  -c <cyc>  cycle type:\n"
"              f - full rotation\n"
"              o - alternating between left and_ang(cyc_val) right octahedron stages\n"
"              r - rolling between left and right octahedron stages\n"
"              i - alternating between left and right icosahedron stages\n"
"  -l <size> lattice size, nummber of jitterbugs along an edge (default 1)\n"
"  -r        allow lattice reversal by fixing connections between\n"
"            jitterbug neighbours\n"
"  -t        fix an opposing pair of triangles, stopping them from rotating\n"
"  -T        rotate so equilateral triangle is on top (suports -t)\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}


     
void jb_opts::process_command_line(int argc, char **argv)
{
   char errmsg[MSG_SZ];
   extern char *optarg;
   extern int optind, opterr;
   opterr = 0;
   char c;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hf:c:l:tTro:")) != -1 ) {
      if(common_opts(c))
         continue;

      switch(c) {
         case 'f':
            if(strspn(optarg, "xaAeo") != strlen(optarg)) {
               snprintf(errmsg, MSG_SZ, "faces to show are %s, must be x, a, e, or o\n", optarg);
               error(errmsg, c);
            }
            switch(*optarg) {
               case 'a':
                  faces = f_equ | f_oth;
                  break;
               case 'A':
                  faces = f_equ | f_quad;
                  break;
               case 'e':
                  faces = f_equ;
                  break;
               case 'o':
                  faces = f_oth;
                  break;
               default:
                  faces = f_none;
            }
            break;
            
         case 'c':
            if(strspn(optarg, "fori") != strlen(optarg)) {
               snprintf(errmsg, MSG_SZ, "cycle type is %s, must be f, o, r, or i\n", optarg);
               error(errmsg, c);
            }
            cycle_type = *optarg;
            break;
            
         case 'l':
            if(!read_int(optarg, &lat_size, errmsg))
               error(errmsg, c);
            if(lat_size <= 0)
               error("lattice size must be a positive integer", c);
            break;

         case 't':
            tri_fix = true;
            break;

         case 'T':
            tri_top = true;
            break;

         case 'r':
            fix_connect = true;
            break;

         case 'o':
            ofile = optarg;
            break;

         default:
            error("unknown command line error");
      }
   }

   if(argc-optind > 1) {
      error("too many arguments");
      exit(1);
   }
  
   stage = 0.0;
   if(argc-optind == 1) {
       if(!read_double(argv[optind], &stage, errmsg))
            error(errmsg, c);
   }
   stage = stage - floor(stage);
}


double get_tri_ang(double cyc_val)
{
   double tri_ang;
   double c = fmod(cyc_val, 0.25)*4;
   if((cyc_val>0.25&&cyc_val<0.5) || cyc_val>0.75)
      c=1-c;
   double ang = asin(sqrt(3*c*c/(4*c*c-2*c+1)));
   if(cyc_val<0.25)
      tri_ang = ang;
   else if(cyc_val<0.5)
      tri_ang = (M_PI - ang)*1;
   else if(cyc_val<0.75)
      tri_ang = (M_PI + ang)*1;
   else
      tri_ang = (2*M_PI - ang)*1;

   return tri_ang;
}


double jb_make_verts(col_geom_v &geom, double t)
{
   geom = col_geom_v();
   vector<vec3d> &verts = *geom.get_verts();

   t = t - floor(t);
   double x, y, z;
   x = 0;
   if(t>=0 && t<0.25) {
      y = 1;
      z = 1 - 8*t;
   }
   else if(t>=0.25 && t<0.5) {
      y = 3 - 8*t;
      z = -1;
   }
   else if(t>=0.5 && t<0.75) {
      y = -1;
      z = 8*t - 5;
   }
   else {
      y = 8*t - 7;
      z = 1;
   }

    verts.push_back(vec3d(x, y, z));
    verts.push_back(vec3d(x, -y, z));
    verts.push_back(vec3d(x, -y, -z));
    verts.push_back(vec3d(x, y, -z));
    verts.push_back(vec3d(z, x, y));
    verts.push_back(vec3d(z, x, -y));
    verts.push_back(vec3d(-z, x, -y));
    verts.push_back(vec3d(-z, x, y));
    verts.push_back(vec3d(y, z, x));
    verts.push_back(vec3d(-y, z, x));
    verts.push_back(vec3d(-y, -z, x));
    verts.push_back(vec3d(y, -z, x));

    double scale = 1.0/(verts[4] - verts[0]).mag();
    for(unsigned int i=0; i<12; i++)
        verts[i] *= scale;
    
    return scale;
}

void jb_add_equ_tris(col_geom_v &jb)
{
   int faces[] = {4,8,0, 5,3,8, 6,9,3, 7,0,9, 7,10,1, 6,2,10, 5,11,2, 4,1,11};
   vector<int> face(3);
   for(int i=0; i<8; i++) {
      for(int j=0; j<3; j++)
         face[j] = faces[3*i + j];
      jb.get_faces()->push_back(face);
      jb.set_f_col(jb.get_faces()->size()-1, 1);
   }
}

void jb_add_oth_quads(col_geom_v &jb)
{
   int quads[] = {7,1,4,0, 11,5,8,4, 6,3,5,2, 9,6,10,7, 3,9,0,8, 2,11,1,10};
   vector<int> face(4);
   for(int i=0; i<6; i++) {
      for(int j=0; j<4; j++)
         face[j] = (quads[4*i + j]);
      jb.get_faces()->push_back(face);
      jb.set_f_col(jb.get_faces()->size()-1, 2);
   }
}

void jb_add_oth_tris(col_geom_v &jb)
{
   int quads[] = {7,1,4,0, 11,5,8,4, 6,3,5,2, 9,6,10,7, 3,9,0,8, 2,11,1,10};
   vector<int> face(3);
   for(int i=0; i<6; i++) {
      for(int j=0; j<3; j++)
         face[j] = (quads[4*i + j]);
      jb.get_faces()->push_back(face);
      jb.set_f_col(jb.get_faces()->size()-1, 2);
      for(int j=0; j<3; j++)
         face[j] = quads[4*i + (j+2)%4];
      jb.get_faces()->push_back(face);
      jb.set_f_col(jb.get_faces()->size()-1, 2);
   }
}

double jb_make(col_geom_v &jb, double t, int f_types)
{
   double scale = jb_make_verts(jb, t);
   if(f_types & f_equ)
      jb_add_equ_tris(jb);
   if(f_types & f_oth)
      jb_add_oth_tris(jb);
   if(f_types & f_quad)
      jb_add_oth_quads(jb);
   return scale;
}


double cycle_value(double fract, char cycle_type)
{
   double val;
   double ico_fract = (3-sqrt(5))/16;
   switch(cycle_type) {
      case 'f':     /* full cycle */
         return fract;
      case 'o':     /* 1/8 cycle to octahedron */
         if(fract<0.25)
            val = fract/2;
         else if(fract<0.75)
            val = 0.125 - (fract - 0.25)/2;
         else
            val = -0.125 + (fract - 0.75)/2;
         if(val<0)
            val = val+1;
         return val;
      case 'r':     /* 1/8 cycle to octahedron with rotation*/
         if(fract<0.5)
            val = fract/4;
         else
            val = -0.125 + (fract - 0.5)/4;
         if(val<0)
            val = val+1;
         return val;
      case 'i':     /* 1/8 cycle to octahedron with rotation*/
         if(fract<0.25)
            val = 4*fract*ico_fract;
         else if(fract<0.75)
            val = ico_fract - 4*(fract - 0.25)*ico_fract;
         else
            val = -ico_fract + 4*(fract - 0.75)*ico_fract;
         if(val<0)
            val = val+1;
         return val;
      default:
         return 0;
   }
}

double get_r_scale(double c_val)
{
   double r_scale;
   if(c_val < 0.25)
      r_scale = 1;
   else if(c_val < 0.5)
      r_scale = 8*(0.375 - c_val);
   else if(c_val < 0.75)
      r_scale = -1;
   else
      r_scale = 8*(c_val - 0.875);

   return r_scale;
}

void lattice_make(col_geom_v &lat_jb, col_geom_v &jb, double cycle_val, double scale, int lat_sz, bool fix_connect)
{
   for(int i=1-lat_sz; i<lat_sz; i+=2)
      for(int j=1-lat_sz; j<lat_sz; j+=2)
         for(int k=1-lat_sz; k<lat_sz; k+=2) {
            col_geom_v jb_elem = jb;
            double scl = scale;
            if(fix_connect)
               scl *= get_r_scale(cycle_val);
            
            mat3d trans = mat3d::transl(vec3d(i*scl, j*scl, k*scl));
            for(unsigned int v=0; v<12; v++)
               (*jb_elem.get_verts())[v] = trans * (*jb_elem.get_verts())[v];
            lat_jb.append(jb_elem);
         }
}

void lattice_trans(col_geom_v &lat_jb, double cyc_val, bool tri_fix, bool tri_top)
{
   mat3d trans;
   if(tri_fix)
      trans = mat3d::rot(vec3d(1,1,1), get_tri_ang(cyc_val))
            * trans;
   if(tri_top)
      trans = mat3d::rot(vec3d(1,0,0), vec3d(1,0,1))
            * mat3d::rot(vec3d(1,1,1), vec3d(0,1,0))
            * trans;
   for(unsigned int v=0; v<lat_jb.get_verts()->size(); v++)
      (*lat_jb.get_verts())[v] = trans * (*lat_jb.get_verts())[v];
}


int main(int argc, char *argv[])
{
   jb_opts opts;
   opts.process_command_line(argc, argv);

   char errmsg[MSG_SZ];

   double cyc_val = cycle_value(opts.stage, opts.cycle_type);
   col_geom_v jb, lat_jb;
   double scale = jb_make(jb, cyc_val, opts.faces);
   lattice_make(lat_jb, jb, cyc_val, scale, opts.lat_size, opts.fix_connect);
   lattice_trans(lat_jb, cyc_val, opts.tri_fix, opts.tri_top); 
   
   if(!lat_jb.write(opts.ofile, errmsg))
      opts.error(errmsg);

   return 0;
}
   

