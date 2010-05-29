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
   Name: off2m.cc
   Description: Convert files in OFF format to 'm' format for LiveGraphics3D
   Project: Antiprism - http://www.antiprism.com
*/



#include <stdio.h>
#include <stdlib.h>

#include <ctype.h>
#include <unistd.h>

#include <string>
#include <vector>

#include "../base/antiprism.h"


using std::string;
using std::vector;


class o2m_opts: public prog_opts {
   public:
      double edge_size;
      double point_size;
      bool lighting;
      string exclude_elems;
      col_val face_col;
      col_val edge_col;
      col_val vert_col;
      int sig_digits;
      int f_dtype;
      col_val bg_col;
      vec3d view_point;
      string ifile;
      string ofile;

      o2m_opts(): prog_opts("off2m"),
             edge_size(0.0),
             point_size(0.0),
             lighting(false),
             exclude_elems(""),
             // face colors need to be explictly set since LG3D defaults to black
             face_col(vec3d(0.8,0.9,0.9)),
             // no longer set edge and vert colors, let them default to missing
             //edge_col(vec3d(0,0,0)),
             //vert_col(vec3d(0,0,0)),
             sig_digits(17),
             f_dtype(1),
             bg_col(vec3d(1,1,1)),
             view_point(vec3d(0,0,0))
             {}

      void process_command_line(int argc, char **argv);
      void usage();
};



void o2m_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Convert files in OFF format to 'm' format for display in LiveGraphics3D. If\n"
"input_file is not given the program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"  -v <size> vertex sphere size (default: 0.02 of bounding box diagonal)\n"
"  -e <size> frame model edge thickness size (default: 0.001 of bounding\n"
"            box diagonal)\n"
"  -x <elms> hide elements. The element string can include v, e and f\n"
"            to hide vertices, edges and faces\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"Coloring options\n"
"  Note: alpha values are ignored. color name \"invisible\" not allowed\n"
"  -V <col>  default vertex colour, in form 'R,G,B' (3 values\n"
"               0.0-1.0, or 0-255) or hex 'xFFFFFF'\n"
"  -E <col>  default edge colour, in form 'R,G,B' (3 values\n"
"               0.0-1.0, or 0-255) or hex 'xFFFFFF'\n"
"  -F <col>  default face colour, in form 'R,G,B' (3 values\n"
"               0.0-1.0, or 0-255) or hex 'xFFFFFF' (default: 0.8,0.8,0.9)\n"
"  -l        let LiveGraphics3D do the coloring itself\n"
"\n"
"Scene options\n"
"  -Y <view> specify the Live3D ViewPoint in form 'X,Y,Z'\n"
"  -B <col>  background colour, in form 'R,G,B' (3 values\n"
"               0.0-1.0, or 0-255) or hex 'xFFFFFF' (default: 1,1,1)\n"
"\n"
"Precision options\n"
"  -d <dgts> number of significant digits (default 17) or if negative\n"
"               then the number of digits after the decimal point\n"
"  -t <type> display type for faces 0 - polygons, 1 - triangulate\n"
"               polygons (default)\n"
"\n"
"\n", prog_name(), help_ver_text);
}


void o2m_opts::process_command_line(int argc, char **argv)
{
   opterr = 0;
   char c;
   char errmsg[MSG_SZ];
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":he:v:x:lF:E:V:d:t:o:Y:B:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         case 'e':
            if(!read_double(optarg, &edge_size, errmsg))
               error(errmsg, c);
            if(edge_size <= 0)
               error("frame width cannot be zero or negative", c);
            break;

         case 'v':
            if(!read_double(optarg, &point_size, errmsg))
               error(errmsg, c);
            if(point_size < 0)
               error("point size cannot be zero or negative", c);
            break;

         case 'l':
            lighting = true;
         break;

         case 'x':
            if(strspn(optarg, "vef") != strlen(optarg)) {
               snprintf(errmsg, MSG_SZ, "elements to exclude are %s must be v, e, or f\n", optarg);
               error(errmsg, c);
            }
            exclude_elems=optarg;
            break;

         case 'F':
            if(!face_col.read(optarg, errmsg))
               error(errmsg, c);
            if(face_col.is_inv())
               error("face color may not be invisible", c);
            break;

         case 'E':
            if(!edge_col.read(optarg, errmsg))
               error(errmsg, c);
            if(edge_col.is_inv())
               error("edge color may not be invisible", c);
            break;

         case 'V':
            if(!vert_col.read(optarg, errmsg))
               error(errmsg, c);
            if(vert_col.is_inv())
               error("vert color may not be invisible", c);
            break;

         case 'd':
            if(!read_int(optarg, &sig_digits, errmsg))
               error(errmsg, c);
            break;

         case 't':
            if(!read_int(optarg, &f_dtype, errmsg))
               error(errmsg, c);
            if(f_dtype !=0 && f_dtype !=1)
               error("display type for faces must be 0 or 1", c);
            break;

         case 'o':
            ofile = optarg;
            break;

         case 'Y':
            if(!view_point.read(optarg, errmsg))
               error(errmsg, c);
            break;

         case 'B':
            if(!bg_col.read(optarg, errmsg))
               error(errmsg, c);
            break;

         default:
            error("unknown command line error");
      }
   }

   if(argc-optind > 1)
      error("too many arguments");
   
   if(argc-optind == 1)
      ifile=argv[optind];
   
}

string RGBtxt(col_val col)
{
   char buf[128];
   vec3d cv = col.get_vec3d();
   snprintf(buf, 128, "RGBColor[%g, %g, %g]", cv[0], cv[1], cv[2]);
   return buf;
}

string Vtxt(vec3d v, int dgts)
{
   char buf[128];
   if(dgts>0)
      snprintf(buf, 128, "{%.*g, %.*g, %.*g}", dgts, v[0], dgts, v[1],
            dgts,v[2]);
   else
      snprintf(buf, 128, "{%.*f, %.*f, %.*f}", -dgts, v[0], -dgts, v[1],
            -dgts, v[2]);
   return buf;
}


int print_m_solid(FILE *ofile, col_geom_v &geom, int sig_digits,
      col_val face_col)
{
   const vector<vector<int> > &faces = geom.faces();
   const vector<vec3d> &verts = geom.verts();
   
   int last_f = -1;
   for(int i=faces.size()-1; i>=0; i--) {
      col_val fcol = geom.get_f_col((int)i);
      fcol = fcol.is_val() ? fcol : face_col;
      if(!fcol.is_inv()) {
         last_f = i;
         break;
      }
   }
   if (last_f == -1)
      return 0;

   for(int i=0; i<=last_f; i++) {
      fprintf(ofile,"{");

      col_val fcol = geom.get_f_col((int)i);
      fcol = fcol.is_val() ? fcol : face_col;
      if(fcol.is_inv())
         continue;
      string RGB = RGBtxt(fcol);
      fprintf(ofile, "FaceForm[ %s, %s ],\n", RGB.c_str(), RGB.c_str());

      fprintf(ofile,"EdgeForm[], Polygon[{");
      for(unsigned int j=0; j<faces[i].size(); j++) {
         fprintf(ofile, "%s", Vtxt(verts[faces[i][j]], sig_digits).c_str());
         if ((j+1)<faces[i].size())
            fprintf(ofile, ", ");
      }

      fprintf(ofile,"}]}");
      if(i==last_f)
         fprintf(ofile,"}");
      fprintf(ofile,",\n");
   }
   
   return 1;
} 


void print_m_frame_edge(FILE *ofile, vec3d v1, vec3d v2, int sig_digits,
      double edge_size, col_val edge_col)
{
   fprintf(ofile,"{");
   fprintf(ofile, "Thickness[%g], ", edge_size);
   if (edge_col.is_set())
      fprintf(ofile, "%s, ", RGBtxt(edge_col).c_str());
   fprintf(ofile, "Line[{%s, %s}", Vtxt(v1, sig_digits).c_str(),
         Vtxt(v2, sig_digits).c_str());
   fprintf(ofile, "]}");
} 

int print_m_frame(FILE *ofile, col_geom_v &geom, int sig_digits,
      double edge_size, col_val edge_col, bool more_to_print)
{
   const vector<vector<int> > &edges = geom.edges();
   const vector<vec3d> &verts = geom.verts();
   
   int last_e = -1;
   for(int i=edges.size()-1; i>=0; i--) {
      col_val ecol = geom.get_e_col((int)i);
      ecol = ecol.is_val() ? ecol : edge_col;
      if(!ecol.is_inv()) {
         last_e = i;
         break;
      }
   }
   if (last_e == -1)
      return 0;

   for(int i=0; i<=last_e; i++) {
      col_val ecol = geom.get_e_col((int)i);
      ecol = ecol.is_val() ? ecol : edge_col;
      if(ecol.is_inv())
         continue;
      print_m_frame_edge(ofile, verts[edges[i][0]], verts[edges[i][1]],
            sig_digits, edge_size, ecol);
      if (i==last_e && !more_to_print)
         fprintf(ofile, "}");
      fprintf(ofile,",\n");
   }
   
   return 1;
} 


int print_m_points(FILE *ofile, col_geom_v &geom, int sig_digits,
      double point_size, col_val vert_col, bool more_to_print)
{
   const vector<vec3d> &verts = geom.verts();
   vec3d vc;
   
   int last_v = -1;
   for(int i=verts.size()-1; i>=0; i--) {
      col_val vcol = geom.get_v_col((int)i);
      vcol = vcol.is_val() ? vcol : vert_col;
      if(!vcol.is_inv()) {
         last_v = i;
         break;
      }
   }
   if (last_v == -1)
      return 0;

   for(int i=0; i<=last_v; i++) {
      col_val vcol = geom.get_v_col(i);
      vcol = vcol.is_val() ? vcol : vert_col;
      if(vcol.is_inv())
         continue;
      fprintf(ofile,"{");
      if (vcol.is_set())
         fprintf(ofile, "%s, ", RGBtxt(vcol).c_str());
      fprintf(ofile, "PointSize[%g], ", point_size);
      fprintf(ofile, "Point[%s]}", Vtxt(verts[i], sig_digits).c_str());
      
      if (i==last_v && !more_to_print)
         fprintf(ofile, "}");
      fprintf(ofile,",\n");
   }
   
   return 1;
}

void print_m_head(FILE *ofile)
{
   fprintf(ofile,"Graphics3D[{\n");
}

void print_m_tail(FILE *ofile, col_geom_v &geom, bool lighting, col_val &bg,
      vec3d view_point)
{
   bound_sphere b_sph(geom.verts());
   vec3d center = b_sph.get_centre();
   double radius = b_sph.get_radius();

   // If ViewPoint isn't specified, Live3D sets an internal one of 1.3,-2.4,2 so the model would be tilted
   // ViewPoint(0,0,0) makes the model disappear in LG3D so don't allow it
   // The default is to let the Focal Length be normal such that ViewPoint(0,0,radius*100)
   // to set Z back far enough to avoid a "fish eye" view
   if (view_point[0] == 0 && view_point[1] == 0 && view_point[2] == 0)
      view_point[2] = radius * 100;

   fprintf(ofile,"ViewPoint -> {%g,%g,%g}, ",
         view_point[0], view_point[1], view_point[2]);

   // Setting ViewVertical with Y upright makes them model appear the same as in AntiView, off2pov, etc
   fprintf(ofile,"ViewVertical -> {0.0,1.0,0.0},\n");

   fprintf(ofile,"Background -> %s, ", RGBtxt(bg).c_str());
   if(lighting)
      fprintf(ofile,"Lighting->True, ");
   else
      fprintf(ofile,"Lighting->False, ");
   fprintf(ofile,"Boxed->False]");
}

int main(int argc, char *argv[])
{
   o2m_opts opts;
   opts.process_command_line(argc, argv);

   char errmsg[MSG_SZ];
   col_geom_v geom;
   if(!geom.read(opts.ifile, errmsg))
      opts.error(errmsg);
   if(*errmsg)
      opts.warning(errmsg);
   
   if(!opts.edge_col.is_inv())
      geom.add_missing_impl_edges();

   if(opts.f_dtype==1)
      geom.triangulate(col_val::invisible);
      
   bound_box bbox(geom.verts());
   double to_model_units = 2.0/(bbox.get_max()-bbox.get_min()).mag();
   if(opts.point_size==0.0)
      opts.point_size = 0.02/to_model_units;
   if(opts.edge_size==0.0)
      opts.edge_size = 0.001/to_model_units;
      

   FILE *ofile = stdout;  // write to stdout by default
   if(opts.ofile != "") {
      ofile = fopen(opts.ofile.c_str(), "w");
      if(ofile == 0)
         opts.error("could not open output file \'"+opts.ofile+"\'");
   }
   
   bool print_v = !strchr(opts.exclude_elems.c_str(), 'v') && geom.verts().size();
   bool print_e = !strchr(opts.exclude_elems.c_str(), 'e') && geom.edges().size();
   bool print_f = !strchr(opts.exclude_elems.c_str(), 'f') && geom.faces().size();
   
   if (!print_v && !print_e && !print_f)
      opts.error("there are no elements to output",'x');

   print_m_head(ofile);
   if(print_v) {
      if (!print_m_points(ofile, geom, opts.sig_digits, opts.point_size*to_model_units, opts.vert_col, (print_e || print_f)))
         opts.error("all vertices are invisible");
   }
   if(print_e) {
      if (!print_m_frame(ofile, geom, opts.sig_digits, opts.edge_size*to_model_units, opts.edge_col, print_f))
         opts.error("all edges are invisible");
   }
   if(print_f) {
      if (!print_m_solid(ofile, geom, opts.sig_digits, opts.face_col))
         opts.error("all faces are invisible");
   }
   print_m_tail(ofile, geom, opts.lighting, opts.bg_col, opts.view_point);

   if(opts.ofile!="")
      fclose(ofile);

   return 0;
}
