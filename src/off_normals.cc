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
      
      bool unit_normals;
      char normal_type;
      char exclude_normals_elems;
      char force_normals_polarity;
      bool elem_normal_vecs;
      char base_normal_method;
      string show_pointing;
      string show_elems;
      string average_pattern;

      col_val outward_normal_col;
      col_val inward_normal_col;
      col_val hemispherical_normal_col;
      col_val edge_normal_col;
      col_val base_normal_col;

      vec3d center;

      int sig_compare;
      double epsilon;

      off_normals_opts(): prog_opts("off_normals"),
                          unit_normals(false),
                          normal_type('a'),
                          exclude_normals_elems('\0'),
                          force_normals_polarity('\0'),
                          elem_normal_vecs(false),
                          base_normal_method('b'),
                          show_pointing("oih"),
                          show_elems("f"),
                          average_pattern("r"),
                          hemispherical_normal_col(col_val(127,127,127)),
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
"If input_file is not given the program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"  -t <opt>  normal type.  a - added to element,  p - positional  (default: a)\n"
"  -u        unit normals  (raw normals otherwise)\n"
"  -e        connect normal to element centroid  (default for -t p)\n"
"  -p <opt>  force polarity. o - set all outward,  i - set all inward\n"
"               r - reverse both inward and outward\n"
"  -i <elms> include normals. The element string can include o, i and h\n"
"               to show, respectively, outward, inward and hemispherical\n"
"               note: exlusion occurs before -p  (default: oih)\n"
"  -s <elms> include elements. The element string can include v, e and f\n"
"               to show, respectively, vertices, edges and faces  (default: f)\n"
"  -d <opt>  delete elements.  f - delete faces of unincluded normals\n"
"               a - delete all of original model\n"
"  -c <opts> average pattern string for edge and vertex normals. Done before -p\n"
"               r - raw,  o - outward,  i - inward,  u - unit  (default: r)\n"
"  -C <xyz>  center of model, in form 'X,Y,Z'  (default: centroid)\n"
"  -l <lim>  minimum distance for unique vertex locations as negative exponent\n"
"               (default: %d giving %.0e)\n"
"  -o <file> write output to file  (default: write to standard output)\n"
"\nColoring Options (run 'off_util -H color' for help on color formats)\n"
"  -O <col>  outward normal vertex color\n"
"  -I <col>  inward normal vertex color\n"
"               default: vertex color is negative of outward col\n"
"  -H <col>  hemispherical normal vertex color  (default: gray50)\n"
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

   while( (c = getopt(argc, argv, ":ht:uep:i:s:d:c:O:I:H:E:B:C:l:o:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) { 
         case 't':
            if(strlen(optarg) > 1 || !strchr("ap", *optarg))
               error("normal type is '"+string(optarg)+"' must be a or p", c);
            normal_type = *optarg;
            break;

         case 'u':
            unit_normals = true;
            break;

         case 'e':
            elem_normal_vecs = true;
            break;

         case 'p':
             if(strlen(optarg) > 1 || !strchr("oir", *optarg))
               error(msg_str("force polarity is '%s', must be o, i, or r", optarg), c);
            force_normals_polarity = *optarg;
            break;

         case 'i':
            if(strspn(optarg, "oih") != strlen(optarg))
               error(msg_str("pointing to include are '%s', must be from o, i, and h", optarg), c);
            show_pointing = optarg;
            break;

         case 's':
            if(strspn(optarg, "vef") != strlen(optarg))
               error(msg_str("elements to hide are '%s', must be from v, e, and f", optarg), c);
            show_elems = optarg;
            break;

         case 'd':
            if(strlen(optarg) > 1 || !strchr("fa", *optarg))
               error(msg_str("delete elements is '%s', must be f or a", optarg), c);
            exclude_normals_elems = *optarg;
            break;

         case 'c':
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

        case 'H':
            if(!hemispherical_normal_col.read(optarg, errmsg))
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

   if (normal_type == 'a')
      elem_normal_vecs = true;

   epsilon = (sig_compare != INT_MAX) ? pow(10, -sig_compare) : ::epsilon;
}


void add_normals(col_geom_v &geom, const bool &unit_normals, const char &normal_type, const char &exclude_normals_elems, const char &force_normals_polarity,
                 const bool &elem_normal_vecs,
                 const col_val &outward_normal_col, const col_val &inward_normal_col, const col_val &hemispherical_normal_col, const col_val &edge_normal_col, const col_val &base_normal_col,
                 const char &base_normal_method, const string &show_elems, const string &show_pointing, const string &average_pattern, const vec3d &center, const double &eps)
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

         bool plotted = false;
         
         vec3d normal;
         if (x_normal.is_hemispherical()) {
            if (strchr(show_pointing.c_str(), 'h')) {
               plotted = true;

               normal = x_normal.raw();
               col = hemispherical_normal_col;
            }
         }
         else {
            if (x_normal.is_inward()) { // normal points inward
               if (strchr(show_pointing.c_str(), 'i')) {
                  plotted = true;

                  if ((force_normals_polarity == 'o') || (force_normals_polarity == 'r')) { // force it outwards
                     normal = x_normal.outward();
                     col = outward_col;
                  }
                  else {
                     normal = x_normal.inward();
                     col = inward_col;
                  }
               }
            }
            
            if (x_normal.is_outward()) { // normal points outward
               if (strchr(show_pointing.c_str(), 'o')) {
                  plotted = true;

                  if ((force_normals_polarity == 'i') || (force_normals_polarity == 'r')) { // force it inwards
                     normal = x_normal.inward();
                     col = inward_col;
                  }
                  else {
                     normal = x_normal.outward();
                     col = outward_col;
                  }
               }
            }
         }

         if (!plotted)
            deleted_faces.push_back(i); // if deleting faces of unplotted normals
         else {       
            if (!normal.is_set())
               normal = x_normal.raw();

            if (unit_normals)
               normal = normal.unit();

            vec3d fc = geom.face_cent(i);
            if (normal_type == 'a')
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
            if ((force_normals_polarity == 'o') || (force_normals_polarity == 'r')) { // force it outwards
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
            if ((force_normals_polarity == 'i') || (force_normals_polarity == 'r')) { // force it inwards
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

         if (unit_normals)
            normal = normal.unit();

         vec3d ec = centroid(geom.verts(), edge);
         if (normal_type == 'a')
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
            if ((force_normals_polarity == 'o') || (force_normals_polarity == 'r')) { // force it outwards
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
            if ((force_normals_polarity == 'i') || (force_normals_polarity == 'r')) { // force it inwards
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

         if (unit_normals)
            normal = normal.unit();

         if (normal_type == 'a')
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
   
   if (exclude_normals_elems == 'f')
      geom.delete_faces(deleted_faces);
   else
   if (exclude_normals_elems == 'a')
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
      
   add_normals(geom, opts.unit_normals, opts.normal_type, opts.exclude_normals_elems, opts.force_normals_polarity,
               opts.elem_normal_vecs, opts.outward_normal_col, opts.inward_normal_col, opts.hemispherical_normal_col, opts.edge_normal_col, opts.base_normal_col,
               opts.base_normal_method, opts.show_elems, opts.show_pointing, opts.average_pattern, opts.center, opts.epsilon);

   if(!geom.write(opts.ofile, errmsg))
      opts.error(errmsg);

   return 0;
}
