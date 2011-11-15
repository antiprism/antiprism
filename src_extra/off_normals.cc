/*
   Copyright (c) 2011, Roger Kaufman

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
   Name: off_normals.cc
   Description: Display normals of faces, implicit edges, or vertices
   Project: Antiprism - http://www.antiprism.com
*/


#include <stdio.h>
#include <stdlib.h>

#include <ctype.h>

#include <string>
#include <vector>
#include <algorithm>

#include "../base/antiprism.h"

using std::string;
using std::vector;
using std::pair;
using std::swap;
 

class off_normals_opts: public prog_opts {
   public:
      string ifile;
      string ofile;
      
      char plot_elem_normals_method;
      char normal_type;
      char exclude_normals_elems;
      char elem_normals_pointing;
      bool elem_normal_vecs;
      col_val outward_normal_col;
      col_val inward_normal_col;
      col_val edge_normal_col;
      col_val base_normal_col;
      char base_normal_method;
      string show_pointing;
      string show_elems;
      string average_pattern;
      vec3d center;

      int sig_compare;
      double epsilon;

      off_normals_opts(): prog_opts("off_normals"),
                          plot_elem_normals_method('n'),
                          normal_type('n'),
                          exclude_normals_elems('\0'),
                          elem_normals_pointing('\0'),
                          elem_normal_vecs(false),
                          base_normal_method('b'),
                          show_pointing("oih"),
                          show_elems("f"),
                          average_pattern("r"),
                          sig_compare(INT_MAX),
                          epsilon(0)
                          {}

      void process_command_line(int argc, char **argv);
      void usage();
};

void off_normals_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Display normals of faces, implicit edges, and vertices\n"
"input_file is not given the program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"  -n <opts> plot element normals\n"
"               n - raw normals (default),  u - unit normals\n"
"               o - set all outward,  i - set all inward\n"
"               r - reverse both inward and outward\n"
"               e - connect normal to element centroid\n"
"               d - delete faces of eliminated normals\n"
"               D - delete all of original model\n"
"  -t <opts> normal type\n"
"               n - positional, p - perpendicular to element (default: n)\n"
"  -i <elms> show pointing. The element string can include o, i and h\n"
"               to show, respectively, outward, inward and hemispherical\n"
"               note: exlusion occurs before o, i, r of -n (default: oih)\n"
"  -s <elms> show elements. The element string can include v, e and f\n"
"               to show, respectively, vertices, edges and faces (default: f)\n"
"  -p <opts> average pattern string for edge and vertex normals. Done before -n\n"
"               r - raw, o - outward, i - inward, u - unit  (default: r)\n"
"  -C <xyz>  center of model, in form 'X,Y,Z' (default: centroid)\n"
"  -l <lim>  minimum distance for unique vertex locations as negative exponent\n"
"               (default: %d giving %.0e)\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\nColoring Options (run 'off_util -H color' for help on color formats)\n"
"  -O <col>  outward normal vertex color\n"
"  -I <col>  inward normal vertex color\n"
"               default: vertex color is negative of outward col\n"
"  -E <col>  normal vector color. connected to element centroid\n"
"               default: color of normal vertex\n"
"  -B <col>  normal vector base color. color at element centroid\n"
"               key word: b take color of element (default)\n"
"               key word: n take color of normal vertex\n"
"\n"
"\n",prog_name(), help_ver_text, int(-log(::epsilon)/log(10) + 0.5),::epsilon);
}

void off_normals_opts::process_command_line(int argc, char **argv)
{
   opterr = 0;
   char c;
   char errmsg[MSG_SZ];
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hn:t:i:s:p:O:I:E:B:C:l:o:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) { 
         case 'n':
            if(strspn(optarg, "nuhHOIoiredD") != strlen(optarg)) {
               snprintf(errmsg, MSG_SZ, "plot elem normal options %s must be a combination of n, u, o, i, r, e, d, and D\n", optarg);
               error(errmsg, c);
            }
            plot_elem_normals_method = 'n'; //default            
            if(strchr(optarg, 'u'))
               plot_elem_normals_method = 'u';

            if(strchr(optarg, 'i'))
               elem_normals_pointing = 'i';
            else
            if(strchr(optarg, 'o'))
               elem_normals_pointing = 'o';
            else
            if(strchr(optarg, 'r'))
               elem_normals_pointing = 'r';
               
            if(strchr(optarg, 'e'))
               elem_normal_vecs = true;
               
            if(strchr(optarg, 'd'))
               exclude_normals_elems = 'd';
            else
            if(strchr(optarg, 'D'))
               exclude_normals_elems = 'D';
            break;

         case 't':
            if(strlen(optarg) > 1 || !strchr("np", *optarg))
               error("normal type is '"+string(optarg)+"' must be n or p", c);
            normal_type = *optarg;
            break;

         case 'i':
            if(strspn(optarg, "oih") != strlen(optarg))
               error(msg_str("pointing to include are '%s', must be from o, i, and h", optarg), c);
            show_pointing=optarg;
            break;

         case 's':
            if(strspn(optarg, "vef") != strlen(optarg))
               error(msg_str("elements to hide are '%s', must be from v, e, and f", optarg), c);
            show_elems=optarg;
            break;

         case 'p':
            if(strspn(optarg, "roiu") != strlen(optarg))
               error(msg_str("average string is '%s', must be from r, o, i, and u", optarg), c);
            average_pattern=optarg;
            break;

         case 'O':
            if(!outward_normal_col.read(optarg, errmsg))
               error(errmsg, c);
            break;
            
         case 'I':
            if(!inward_normal_col.read(optarg, errmsg))
               error(errmsg, c);
            break;
            
         case 'E':
            if(!edge_normal_col.read(optarg, errmsg))
               error(errmsg, c);
            break;
            
         case 'B':
            if(strchr("bn", *optarg))
               base_normal_method = *optarg;
            else
            if(!base_normal_col.read(optarg, errmsg))
               error(errmsg, c);
            break;

         case 'C':
            if(!center.read(optarg, errmsg))
               error(errmsg, c);
            break;

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

   epsilon = (sig_compare != INT_MAX) ? pow(10, -sig_compare) : ::epsilon;
}


void add_normals(col_geom_v &geom, const char &plot_elem_normals_method, const char &normal_type, const char &exclude_normals_elems, const char &elem_normals_pointing,
                 const bool &elem_normal_vecs,
                 const col_val &outward_normal_col, const col_val &inward_normal_col, const col_val &edge_normal_col, const col_val &base_normal_col, const char &base_normal_method,
                 const string &show_elems, const string &show_pointing, const string &average_pattern, const vec3d &center, const double &eps)
{
   col_val outward_col = outward_normal_col;
   col_val inward_col = inward_normal_col;
   if (!inward_col.is_set() && outward_normal_col.is_set()) {
      inward_col = outward_normal_col;
      inward_col.set_complement();
   }

   col_val col;
   
   col_geom_v ngeom;
   vector<int> deleted_faces;

   fnormals x_normals(geom, center, eps);

   if (strchr(show_elems.c_str(), 'f')) {
      for(unsigned int i=0;i<x_normals.size();i++) {
         xnormal x_normal = x_normals[i];
         
         if (x_normal.is_hemispherical()) {
            if (!strchr(show_pointing.c_str(), 'h')) {
               deleted_faces.push_back(i);
               continue;
            }
         }

         vec3d normal;
         if (x_normal.is_inward()) { // normal points inward
            if (!strchr(show_pointing.c_str(), 'i')) {
               deleted_faces.push_back(i);
               continue;
            }
            if ((elem_normals_pointing == 'o') || (elem_normals_pointing == 'r')) { // force it outwards
               normal = x_normal.outward();
               col = outward_col;
            }
            else {
               normal = x_normal.inward();
               col = inward_col;
            }
         }
         
         if (x_normal.is_outward()) { // normal points outward
            if (!strchr(show_pointing.c_str(), 'o')) {
               deleted_faces.push_back(i);
               continue; 
            }
            if ((elem_normals_pointing == 'i') || (elem_normals_pointing == 'r')) { // force it inwards
               normal = x_normal.inward();
               col = inward_col;
            }
            else {
               normal = x_normal.outward();
               col = outward_col;
            }
         }

         //if(find(deleted_faces.begin(), deleted_faces.end(), i) != deleted_faces.end())
         //   continue;
         
         if (!normal.is_set())
            normal = x_normal.raw();

         if (plot_elem_normals_method == 'u')
            normal = normal.unit();

         vec3d fc = geom.face_cent(i);
         if (normal_type == 'p')
            normal += fc;

         ngeom.add_col_vert(normal,col);

         if (elem_normal_vecs) {
            // get base color
            col_val bcol = (base_normal_col.is_set()) ? base_normal_col : ((base_normal_method == 'b') ? geom.get_f_col((int)i) : col);
            // add point at centroid
            ngeom.add_col_vert(fc, bcol);
            // get edge color
            col_val ecol = (edge_normal_col.is_set()) ? edge_normal_col : col;
            // edge from face centroid to normal
            ngeom.add_col_edge(make_edge(ngeom.verts().size()-1, ngeom.verts().size()-2), ecol);
         }
      }
   }

   if (strchr(show_elems.c_str(), 'e')) {
      vector<vector<int> > implicit_edges;
      geom.get_impl_edges(implicit_edges);

      for(unsigned int i=0;i<implicit_edges.size();i++) {
         vector<int> edge = implicit_edges[i];
         xnormal x_normal = x_normals.edge_normal(edge[0], edge[1], average_pattern);

         if (x_normal.is_hemispherical()) {
           if (!strchr(show_pointing.c_str(), 'h')) {
               continue;
            }
         }

         vec3d normal;
         if (x_normal.is_inward()) { // normal points inward
            if (!strchr(show_pointing.c_str(), 'i')) {
               continue; 
            }
            if ((elem_normals_pointing == 'o') || (elem_normals_pointing == 'r')) { // force it outwards
               normal = x_normal.outward();
               col = outward_col;
            }
            else {
               normal = x_normal.inward();
               col = inward_col;
            }
         }
         
         if (x_normal.is_outward()) { // normal points outward
            if (!strchr(show_pointing.c_str(), 'o')) {
               continue; 
            }
            if ((elem_normals_pointing == 'i') || (elem_normals_pointing == 'r')) { // force it inwards
               normal = x_normal.inward();
               col = inward_col;
            }
            else {
               normal = x_normal.outward();
               col = outward_col;
            }
         }
         
         if (!normal.is_set())
            normal = x_normal.raw();

         if (plot_elem_normals_method == 'u')
            normal = normal.unit();

         vec3d ec = centroid(geom.verts(), edge);
         if (normal_type == 'p')
            normal += ec;

         ngeom.add_col_vert(normal,col);

         if (elem_normal_vecs) {
            // get base color for edge. might not explicitly be colored
            int e_idx = find_edge_in_edge_list(geom.edges(), edge);
            col_val expl_col = (e_idx > -1) ? geom.get_e_col(e_idx) : col_val();
               
            col_val bcol = (base_normal_col.is_set()) ? base_normal_col : ((base_normal_method == 'b') ? expl_col : col);
            // add point at centroid
            ngeom.add_col_vert(ec, bcol);
            // get edge color
            col_val ecol = (edge_normal_col.is_set()) ? edge_normal_col : col;
            // edge from face centroid to normal
            ngeom.add_col_edge(make_edge(ngeom.verts().size()-1, ngeom.verts().size()-2), ecol);
         }
      }
   }

   if (strchr(show_elems.c_str(), 'v')) {
      const vector<vec3d> &verts = geom.verts();

      for(unsigned int i=0;i<verts.size();i++) {
         xnormal x_normal = x_normals.vertex_normal(i, average_pattern);

         if (x_normal.is_hemispherical()) {
           if (!strchr(show_pointing.c_str(), 'h')) {
               continue;
            }
         }

         vec3d normal;
         if (x_normal.is_inward()) { // normal points inward
            if (!strchr(show_pointing.c_str(), 'i')) {
               continue; 
            }
            if ((elem_normals_pointing == 'o') || (elem_normals_pointing == 'r')) { // force it outwards
               normal = x_normal.outward();
               col = outward_col;
            }
            else {
               normal = x_normal.inward();
               col = inward_col;
            }
         }
         
         if (x_normal.is_outward()) { // normal points outward
            if (!strchr(show_pointing.c_str(), 'o')) {
               continue; 
            }
            if ((elem_normals_pointing == 'i') || (elem_normals_pointing == 'r')) { // force it inwards
               normal = x_normal.inward();
               col = inward_col;
            }
            else {
               normal = x_normal.outward();
               col = outward_col;
            }
         }
         
         if (!normal.is_set())
            normal = x_normal.raw();

         if (plot_elem_normals_method == 'u')
            normal = normal.unit();

         if (normal_type == 'p')
            normal += verts[i];

         ngeom.add_col_vert(normal,col);

         if (elem_normal_vecs) {
            // get base color
            col_val bcol = (base_normal_col.is_set()) ? base_normal_col : ((base_normal_method == 'b') ? geom.get_v_col((int)i) : col);
            // add point at centroid
            ngeom.add_col_vert(verts[i], bcol);
            // get edge color
            col_val ecol = (edge_normal_col.is_set()) ? edge_normal_col : col;
            // edge from face centroid to normal
            ngeom.add_col_edge(make_edge(ngeom.verts().size()-1, ngeom.verts().size()-2), ecol);
         }
      }
   }
   
   if (exclude_normals_elems == 'd')
      geom.delete_faces(deleted_faces);
   else
   if (exclude_normals_elems == 'D')
      geom.clear_all();
   geom.append(ngeom);
}

/*
vec3d line_nearest_point(vec3d P, vec3d A, vec3d B)
{
   vec3d v1 = P-A;
   vec3d v2 = B-A;

   double v2m = v2.mag2();

   double D = vdot(v1,v2);

   double dist = D/v2m;

   return (A+(v2*dist));
}
*/

int main(int argc, char *argv[])
{
   off_normals_opts opts;
   opts.process_command_line(argc, argv);

   char errmsg[MSG_SZ];
   col_geom_v geom;
   if(!geom.read(opts.ifile, errmsg))
      opts.error(errmsg);
   if(*errmsg)
      opts.warning(errmsg);
      
   if(opts.plot_elem_normals_method)
      add_normals(geom, opts.plot_elem_normals_method, opts.normal_type, opts.exclude_normals_elems, opts.elem_normals_pointing,
                  opts.elem_normal_vecs, opts.outward_normal_col, opts.inward_normal_col, opts.edge_normal_col, opts.base_normal_col, opts.base_normal_method,
                  opts.show_elems, opts.show_pointing, opts.average_pattern, opts.center, opts.epsilon);

   if(!geom.write(opts.ofile, errmsg))
      opts.error(errmsg);

   return 0;
}
