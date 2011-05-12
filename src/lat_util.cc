/*
   Copyright (c) 2008-2009, Roger Kaufman

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
   Name: lat_util.cc
   Description: Do various things to lattices
   Project: Antiprism - http://www.antiprism.com
*/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include <ctype.h>

#include <string>
#include <vector>

#include "../base/antiprism.h"
#include "lattice_grid.h"

using std::string;
using std::vector;


class lutil_opts: public prog_opts {
   public:
      vector<string> ifiles;
      string ofile;
      string cfile;
      string rfile;

      vector<double> strut_len;
      bool strip_faces;
      char color_edges_by_sqrt;
      bool is_vertex_color;
      char container;
      bool append_container;
      double radius;
      char radius_default;
      vec3d offset;
      bool voronoi_cells;
      bool voronoi_central_cell;
      bool convex_hull;
      bool add_hull;
      bool append_lattice;
      char color_method;
      int face_opacity;
      bool list_radii;
      bool list_struts;
      col_val cent_col;
      bool trans_to_origin;

      double epsilon;
      vector<col_val> remove_vertex_color_list;
      
      // 0 - lattice  1 - convex hull  2 - voronoi
      vector<col_val> vert_col;
      vector<col_val> edge_col;
      vector<col_val> face_col;

      lutil_opts(): prog_opts("lat_util"),
                    strip_faces(true),
                    color_edges_by_sqrt('\0'),
                    is_vertex_color(true),
                    container('c'),
                    append_container(false),
                    radius(0),
                    radius_default('s'),
                    voronoi_cells(false),
                    voronoi_central_cell(false),
                    convex_hull(false),
                    add_hull(false),
                    append_lattice(false),
                    color_method('\0'),
                    face_opacity(-1),
                    list_radii(false),
                    list_struts(false),
                    trans_to_origin(false),
                    epsilon(0) {}

      void process_command_line(int argc, char **argv);
      void usage();
};

void lutil_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_files]\n"
"\n"
"Read one or more files in OFF format, combine them into a single file and\n"
"process it. Operations take place in the order listed below. input_files is the\n"
"list of files to process. If they don't exist, explicit edges are created.\n"
"If the input possesses faces they are stripped by default.\n"
"\n"
"Options\n"
"%s"
"  -h        this help message\n"
"  -z        suppress stripping of faces\n"
"  -x <colr> remove every vertex of color colr\n"
"               use multiple -x parameters for multiple colors\n"
"  -X <colr> remove every vertex NOT of color colr\n"
"  -c <type> container, c - cube (default), s - sphere (uses radius)\n"
"  -k <file> container using convex polyhedron in off file (uses radius)\n"
"  -r <c,n>  radius. c is radius taken to optional root n. n = 2 is sqrt\n"
"               or  l - max insphere radius  s - min insphere radius (default: s)\n"
"               or  k - take radius from container specified by -k\n"
"  -q <vecs> center offset, in form \"a_val,b_val,c_val\" (default: none)\n"
"  -s <s,n>  create struts. s is strut length taken to optional root n\n"
"               use multiple -s parameters for multiple struts\n"
"  -D <opt>  Voronoi (a.k.a Dirichlet) cells (Brillouin zones for duals)\n"
"               c - cells only, i - cell(s) touching center only\n"
"  -C <opt>  c - convex hull only, i - keep interior\n"
"  -A        append the original lattice to the final product\n"
"  -O        translate center of final product to origin\n"
"  -Z <col>  add center vertex to final product in color col\n"
"  -R <file> repeat off file at every vertex in lattice\n"
"  -K        append cage of container of -k to final product\n"
"  -l <lim>  minimum distance for unique vertex locations as negative exponent\n"
"               (default: %d giving %.0e)\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\nListing Options\n"
"  -L        list unique radial distances of points from center (and offset)\n"
"  -S        list every possible strut value\n"
"\nColoring Options (run 'off_util -H color' for help on color formats)\n"
"  -V <col>  vertex color, (optional) elements, (optional) opacity\n"
"               elements to color are l - lattice  c - convex hull  v - voronoi\n"
"                  (default elements: lcv)\n"
"               opacity valid range from 0 to 255\n"
"                  0 - invisible  255 - opaque (default: 255)\n"
"  -E <col>  edge color (for struts, convex hulls, and voronoi)\n"
"               lower case outputs map indexes. upper case outputs color values\n"
"               key word: r,R for color edges by root value of final product\n"
"  -F <col>  face color (for convex hulls and voronoi)\n"
"               key word: s,S color by symmetry using face normals\n"
"               key word: c,C color by symmetry using face normals (chiral)\n"
"  -T <tran> face opacity for color by symmetry. valid range from 0 to 255\n"
"\n"
"\n",prog_name(), help_ver_text, int(-log(::epsilon)/log(10) + 0.5), ::epsilon);
}

void lutil_opts::process_command_line(int argc, char **argv)
{
   opterr = 0;
   char c;
   char errmsg[MSG_SZ];

   int sig_compare = INT_MAX;
   vector<double> double_parms;
   col_val col_tmp;
   int x_used = 0;
   
   vert_col.resize(3);
   edge_col.resize(3);
   face_col.resize(3);
   
   // set some default colors
   // 0 - lattice  1 - convex hull  2 - voronoi
   // voronoi cells  vcol = gold; ecol = lightgrey; fcol = transparent yellow
   vert_col[2] = col_val(255,215,0);
   edge_col[2] = col_val(211,211,211);
   face_col[2] = col_val(255,255,0,128);
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hzx:X:c:k:r:q:s:D:C:AV:E:F:T:Z:KOR:LSl:o:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         case 'z':
            strip_faces = false;
            break;
            
         case 'x':
            if (!x_used)
               x_used = 1;
            else
            if (x_used < 0)
               error("cannot mix -x with -X",c);
            is_vertex_color = true;
            if(!col_tmp.read(optarg, errmsg))
               error(errmsg, c);
            remove_vertex_color_list.push_back(col_tmp);
            break;

         // if X is called more than once, only the last one is honored
         case 'X':
            if (!x_used)
               x_used = -1;
            else
            if (x_used > 0)
               error("cannot mix -X with -x",c);
            is_vertex_color = false;
            if(!col_tmp.read(optarg, errmsg))
               error(errmsg, c);
            if (remove_vertex_color_list.size()) {
               warning("if more than one -X, only the last one is used");
               remove_vertex_color_list.clear();
            }
            remove_vertex_color_list.push_back(col_tmp);
            break;
            
         case 'c':
            if(strlen(optarg) > 1 || !strchr("cs", *optarg))
               error("method is '"+string(optarg)+"' must be c or s", c);
            container = *optarg;
            break;

         case 'k':
            cfile = optarg;
            break;
            
         case 'r':
            if(strlen(optarg) == 1 && strchr("lsk", *optarg))
               radius_default = *optarg;
            else {
               if(!read_double_list(optarg, double_parms, errmsg, 2))
                  error(errmsg, c);
               radius = double_parms[0];
               if(double_parms.size() == 2) {
                  if(double_parms[1] == 0)
                     error("root for radius must be non-zero", c);
                  radius = pow(radius, 1/double_parms[1]);
               }
               if(radius <= 0)
                  error("radius cannot be negative or zero", "s", c);
            }
            break;
            
         case 'q':
            if(!offset.read(optarg, errmsg))
               error(errmsg, c);
            break;

         case 's': {
            double strut_tmp;
            if(!read_double_list(optarg, double_parms, errmsg, 2))
               error(errmsg, c);
            strut_tmp = double_parms[0];
            if(double_parms.size() == 2)
               strut_tmp = pow(strut_tmp, 1/double_parms[1]);
            if(strut_tmp <= 0)
               error("strut lengths cannot be negative", "s");
            strut_len.push_back(strut_tmp);
            break;
         }

         case 'D':
            if(strlen(optarg) > 1 || !strchr("ci", *optarg))
               error("Voronoi cells arg is '"+string(optarg)+"' must be c, i", c);
            voronoi_cells = true;
            if(strchr("i", *optarg))
               voronoi_central_cell = true;
            break;

         case 'C':
            if(strlen(optarg) > 1 || !strchr("ci", *optarg))
               error("convex hull arg is '"+string(optarg)+"' must be c, i", c);
            convex_hull = true;
            if(strchr("i", *optarg))
               add_hull = true;
            break;
            
         case 'A':
            append_lattice = true;
            break;
            
         case 'V': {
            vector<char *> parts;
            int parts_sz = split_line(optarg, parts, ",");
            if(parts_sz>6)
               error("the argument has more than 6 parts", c);
            
            col_val col;
            bool valid_color = false;
            int next_parms_idx = 1;

            if (parts_sz >= 3) {
               char parts_test[80];
               parts_test[0] = '\0';
               strcat(parts_test,parts[0]);
               strcat(parts_test,",");
               strcat(parts_test,parts[1]);
               strcat(parts_test,",");
               strcat(parts_test,parts[2]);
               if(col.read(parts_test, errmsg)) {
                  valid_color = true;
                  next_parms_idx = 3;
               }
               col_val col_save = col;
               if (parts_sz >=4) {
                  strcat(parts_test,",");
                  strcat(parts_test,parts[3]);
                  if(col.read(parts_test, errmsg)) {
                     valid_color = true;
                     next_parms_idx = 4;
                  }
                  else
                     col = col_save;
               }
            }
            
            if (!valid_color) {
               if(!col.read(parts[0], errmsg))
                  error(errmsg, c);
            }

            unsigned int conv_elems = 15;
            if(parts_sz>next_parms_idx) {
               if(strspn(parts[next_parms_idx], "lcvh") !=
                                              strlen(parts[next_parms_idx]))
                  error(msg_str("elements to map are '%s' must be from "
                           "l, c, v or h", parts[next_parms_idx]), c);
               conv_elems = 8*(strchr(parts[next_parms_idx], 'h')!=0) +
                            4*(strchr(parts[next_parms_idx], 'v')!=0) +
                            2*(strchr(parts[next_parms_idx], 'c')!=0) +
                            1*(strchr(parts[next_parms_idx], 'l')!=0);
            }
            
            int opq = col[3];
            if(parts_sz>next_parms_idx+1)
               opq = atoi(parts[next_parms_idx+1]);
               
            if (opq < 0 || opq > 255)
               error("opacity value must be between 0 and 255", c);

            for(int i=0; i<3; i++) {
              if(conv_elems & (1<<i)) {
                  if (!col.is_set())
                     vert_col[i] = col_val();
                  else
                     vert_col[i] = col_val(col[0],col[1],col[2],opq);
               }
            }
            break;
         }
         
         case 'E': {
            vector<char *> parts;
            int parts_sz = split_line(optarg, parts, ",");
            if(parts_sz>6)
               error("the argument has more than 6 parts", c);
               
            if (!strcasecmp(parts[0],"r")) {
               color_edges_by_sqrt = parts[0][0];
               break;
            }
            
            col_val col;
            bool valid_color = false;
            int next_parms_idx = 1;

            if (parts_sz >= 3) {
               char parts_test[80];
               parts_test[0] = '\0';
               strcat(parts_test,parts[0]);
               strcat(parts_test,",");
               strcat(parts_test,parts[1]);
               strcat(parts_test,",");
               strcat(parts_test,parts[2]);
               if(col.read(parts_test, errmsg)) {
                  valid_color = true;
                  next_parms_idx = 3;
               }
               col_val col_save = col;
               if (parts_sz >=4) {
                  strcat(parts_test,",");
                  strcat(parts_test,parts[3]);
                  if(col.read(parts_test, errmsg)) {
                     valid_color = true;
                     next_parms_idx = 4;
                  }
                  else
                     col = col_save;
               }
            }
            
            if (!valid_color) {
               if(!col.read(parts[0], errmsg))
                  error(errmsg, c);
            }

            unsigned int conv_elems = 15;
            if(parts_sz>next_parms_idx) {
               if(strspn(parts[next_parms_idx], "lcvh") != 
                                            strlen(parts[next_parms_idx]))
                  error(msg_str("elements to map are '%s' must be from "
                           "l, c, v, h", parts[next_parms_idx]), c);
               conv_elems = 8*(strchr(parts[next_parms_idx], 'h')!=0) +
                            4*(strchr(parts[next_parms_idx], 'v')!=0) +
                            2*(strchr(parts[next_parms_idx], 'c')!=0) +
                            1*(strchr(parts[next_parms_idx], 'l')!=0);
            }
            
            int opq = col[3];
            if(parts_sz>next_parms_idx+1)
               opq = atoi(parts[next_parms_idx+1]);
               
            if (opq < 0 || opq > 255)
               error("opacity value must be between 0 and 255", c);

            for(int i=0; i<3; i++) {
              if(conv_elems & (1<<i)) {
                  if (!col.is_set())
                     edge_col[i] = col_val();
                  else
                     edge_col[i] = col_val(col[0],col[1],col[2],opq);
               }
            }
            break;
         }
         
         case 'F': {
            vector<char *> parts;
            int parts_sz = split_line(optarg, parts, ",");
            if(parts_sz>6)
               error("the argument has more than 6 parts", c);
               
            if ((!strcasecmp(parts[0],"s")) || (!strcasecmp(parts[0],"c"))) {
               color_method = parts[0][0];
               break;
            }
            
            col_val col;
            bool valid_color = false;
            int next_parms_idx = 1;

            if (parts_sz >= 3) {
               char parts_test[80];
               parts_test[0] = '\0';
               strcat(parts_test,parts[0]);
               strcat(parts_test,",");
               strcat(parts_test,parts[1]);
               strcat(parts_test,",");
               strcat(parts_test,parts[2]);
               if(col.read(parts_test, errmsg)) {
                  valid_color = true;
                  next_parms_idx = 3;
               }
               col_val col_save = col;
               if (parts_sz >=4) {
                  strcat(parts_test,",");
                  strcat(parts_test,parts[3]);
                  if(col.read(parts_test, errmsg)) {
                     valid_color = true;
                     next_parms_idx = 4;
                  }
                  else
                     col = col_save;
               }
            }
            
            if (!valid_color) {
               if(!col.read(parts[0], errmsg))
                  error(errmsg, c);
            }

            unsigned int conv_elems = 15;
            if(parts_sz>next_parms_idx) {
               if(strspn(parts[next_parms_idx], "lcvh") != 
                                           strlen(parts[next_parms_idx]))
                  error(msg_str("elements to map are '%s' must be from "
                           "l, c, v, h", parts[next_parms_idx]), c);
               conv_elems = 8*(strchr(parts[next_parms_idx], 'h')!=0) +
                            4*(strchr(parts[next_parms_idx], 'v')!=0) +
                            2*(strchr(parts[next_parms_idx], 'c')!=0) +
                            1*(strchr(parts[next_parms_idx], 'l')!=0);
            }
            
            int opq = col[3];
            if(parts_sz>next_parms_idx+1)
               opq = atoi(parts[next_parms_idx+1]);
               
            if (opq < 0 || opq > 255)
               error("opacity value must be between 0 and 255", c);

            for(int i=0; i<3; i++) {
              if(conv_elems & (1<<i)) {
                  if (!col.is_set())
                     face_col[i] = col_val();
                  else
                     face_col[i] = col_val(col[0],col[1],col[2],opq);
               }
            }
            break;
         }
         
         case 'T':
            if(!read_int(optarg, &face_opacity, errmsg))
               error(errmsg, c);
            if(face_opacity < 0 || face_opacity > 255) {
               error("face transparency must be between 0 and 255", c);
            }
            break;
         
         case 'Z':
            if(!cent_col.read(optarg, errmsg))
               error(errmsg, c);
            break;
            
         case 'K':
            append_container = true;
            break;

         case 'O':
            trans_to_origin = true;
            break;
            
         case 'R':
            rfile = optarg;
            break;

         case 'L':
            list_radii = true;
            break;

         case 'S':
            list_struts = true;
            break;
            
         case 'l':
            if(!read_int(optarg, &sig_compare, errmsg))
               error(errmsg, c);
            if(sig_compare < 0) {
               warning("limit is negative, and so ignored", c);
            }
            if(sig_compare > DEF_SIG_DGTS) {
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

   while(argc-optind)
      ifiles.push_back(string(argv[optind++]));

   if (radius_default == 'k' && !cfile.length())
      error("-r k can only be used if -k container is specified");
      
   if (container == 'c' && !cfile.length() && (radius != 0 || offset.is_set()))
      warning("cubic container in use. Radius and offset ignored");
      
   //if ((container == 's' || (cfile.length() && !radius_default)) && !radius)
   //   error("radius not set");
   
   epsilon = (sig_compare != INT_MAX) ? pow(10, -sig_compare) : ::epsilon;
}

void remove_vertex_by_color(col_geom_v &geom, const col_val &remove_vertex_color, const bool &is_vertex_color)
{
   const vector<vec3d> &verts = geom.verts();

   vector<int> del_verts;
   for(unsigned int i=0; i<verts.size(); i++) {
      if (( is_vertex_color && geom.get_v_col(i) == remove_vertex_color) ||
          (!is_vertex_color && geom.get_v_col(i) != remove_vertex_color))
         del_verts.push_back(i);
   }

   if (del_verts.size())
      geom.delete_verts(del_verts);

   if (!verts.size())
      fprintf(stderr,"remove_vertex_by_color: warning: all vertices were removed!\n");
}

void make_skeleton(col_geom_v &geom, const bool &strip_faces)
{
   geom.add_missing_impl_edges();
   if (strip_faces)
      geom.clear_faces();
}

void process_lattices(col_geom_v &geom, col_geom_v &container, const col_geom_v &repeater, lutil_opts &opts)
{
   // add explicit edges and remove faces if necessary
   make_skeleton(geom, opts.strip_faces);
   
   // save lattice in case if adding back in end
   col_geom_v tgeom;
   if (opts.append_lattice)
      tgeom = geom;
   
   geom.color_vef(opts.vert_col[0], opts.edge_col[0], opts.face_col[0]);

   for(unsigned int i=0; i<opts.remove_vertex_color_list.size(); i++)
      remove_vertex_by_color(geom, opts.remove_vertex_color_list[i], opts.is_vertex_color);

   for(unsigned int i=0; i<opts.strut_len.size(); i++)
      add_color_struts(geom, opts.strut_len[i]*opts.strut_len[i], opts.edge_col[0]);
      
   if (!opts.radius && opts.radius_default != 'k')
      opts.radius = lattice_radius(geom, opts.radius_default);
      
   if(opts.cfile.length())
      geom_container_clip(geom, container, (opts.radius_default == 'k') ? lattice_radius(container,opts.radius_default) : opts.radius, opts.offset, opts.epsilon);
   else
   if (opts.container == 's')
      geom_spherical_clip(geom, opts.radius, opts.offset, opts.epsilon);

   if(opts.voronoi_cells) {
      col_geom_v vgeom;
      if (get_voronoi_geom(geom, vgeom, opts.voronoi_central_cell, false, opts.epsilon)) {
         vgeom.color_vef(opts.vert_col[2], opts.edge_col[2], opts.face_col[2]);
         geom = vgeom;
      }
   }

   if (opts.convex_hull) {
      char errmsg[MSG_SZ]="";
      int ret = (opts.add_hull ? geom.add_hull("",errmsg) : geom.set_hull("",errmsg));
      if(ret < 0)
         fprintf(stderr,"%s\n",errmsg);
      else {
         geom.orient();
         if (true) // verbosity
            convex_hull_report(geom, opts.add_hull);
         geom.color_vef(opts.vert_col[1], opts.edge_col[1], opts.face_col[1]);
      }
   }
   
   if (opts.append_lattice) {
      geom.append(tgeom);
      tgeom.clear_all();
   }
   
   // face color by symmetry normals
   if (opts.color_method)
      color_by_symmetry_normals(geom, opts.color_method, opts.face_opacity, opts.epsilon);
   
   // if color by sqrt was used, override all edges of all structure
   if (opts.color_edges_by_sqrt)
      color_edges_by_sqrt(geom, opts.color_edges_by_sqrt);
   
   if (opts.trans_to_origin)
      geom.transform(mat3d::transl(-centroid(geom.verts())));

   if (opts.list_radii)
      list_grid_radii(geom, opts.offset, opts.epsilon);

   if (opts.list_struts)
      list_grid_struts(geom, opts.epsilon);

   // add central vertex last so not to alter listing outcomes
   if (opts.cent_col.is_set())
      color_centroid(geom, opts.cent_col, opts.epsilon);
      
   // place geom at every vertex in lattice
   if(opts.rfile.length()) {
      col_geom_v geom2;
      for(unsigned int i=0; i<geom.verts().size(); i++) {
         col_geom_v rep = repeater;
         rep.transform(mat3d::transl((geom.verts())[i]));
         geom2.append(rep);
      }
      geom = geom2;
      sort_merge_elems(geom, "vef", opts.epsilon);
   }
   
   if (opts.append_container) {
      container.add_missing_impl_edges();
      container.clear_faces();
      geom.append(container);
   }
}

int main(int argc, char *argv[])
{
   lutil_opts opts;
   opts.process_command_line(argc, argv);
   if(!opts.ifiles.size())
      opts.ifiles.push_back("");
      
   // read in container file if using. Check existance
   char errmsg[MSG_SZ]="";
   col_geom_v container;
   if(opts.cfile.length()) {
      if(!container.read(opts.cfile, errmsg))
         opts.error(errmsg);
      if(*errmsg)
         opts.warning(errmsg);
   }
   
   col_geom_v repeater;
   if(opts.rfile.length()) {
      if(!repeater.read(opts.rfile, errmsg))
         opts.error(errmsg);
      if(*errmsg)
         opts.warning(errmsg);
   }

   col_geom_v geoms;
   for(unsigned int i=0; i<opts.ifiles.size(); i++) {
      col_geom_v geom;
      if(!geom.read(opts.ifiles[i], errmsg))
         opts.error(errmsg);
      if(*errmsg)
         opts.warning(errmsg);
      geoms.append(geom);
   }
   
   process_lattices(geoms, container, repeater, opts);

   if(!geoms.write(opts.ofile, errmsg))
      opts.error(errmsg);
   
   return 0;
}
