/*
   Copyright (c) 2012, Adrian Rossiter

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
   Name: leonardo.cc
   Description: convert models to the form illustrated by Leonardo da Vinci,
                having thickened faces with a hole in the centre.
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"


class leo_opts: public prog_opts
{
   private:

   public:
      const double DEF_VAL;
      const double def_width;
      double width;
      bool width_is_perc;
      const double def_height;
      double height;
      bool height_is_perc;
      bool centre_height;
      char centre_type;       // v: vertex normal, c: set centre
      vec3d centre;
      bool hide_edges;
      bool col_from_edges;
      bool panels;

      string ifile;
      string ofile;

      leo_opts(): prog_opts("leonardo"),
                 DEF_VAL(1e100),
                 def_width(30), width(DEF_VAL), width_is_perc(true),
                 def_height(100), height(DEF_VAL), height_is_perc(true),
                 centre_height(false), centre_type('v'),
                 hide_edges(false),
                 col_from_edges(false),
                 panels(false)
                 {}

      void process_command_line(int argc, char **argv);
      void usage();
};

void leo_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Read a file in OFF format, and thicken the faces and cut a hole in their\n"
"centres to produce a model like those illustrated by Leonardo da Vinci. If\n"
"input_file is not given the program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"  -w <wdth> width of the perimeter border of the faces, follow by %% for\n"
"            percentage of maximum width without overlap (default: %g%%)\n"
"  -l <ht>   height to thicken faces, 0 for single polygon height, follow by\n"
"            %% for percentage of width value (default: %g%%)\n"
"  -c <cent> centre point to extrude faces towards, in form 'X,Y,Z', or C to\n"
"            use centroid (default: no centre, use calculated vertex normals)\n"
"  -p        faces converted to panels without holes (incompatible with\n"
"            -w, -x, -e)\n"
"  -m        distribute the height equally on both sides of the faces, so\n"
"            the original faces would lie in the middle of the new faces\n"
"            (use for non-orientable models)\n"
"  -x        hide the edges that join the outside of a face to the hole\n"
"  -e        take colours from the edge colours of the base polyhedron\n"
"            (default: use face colours)\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text, def_width, def_height);
}


void leo_opts::process_command_line(int argc, char **argv)
{
   opterr = 0;
   char c;
   char errmsg[MSG_SZ];
   int len;

   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hw:l:c:pmxieo:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         case 'w':
            len = strlen(optarg);
            if(len && optarg[len-1]=='%') {
               optarg[len-1] = '\0';
               width_is_perc = true;
            }
            else
               width_is_perc = false;
            if(!read_double(optarg, &width, errmsg))
               error(errmsg, c);
            break;

         case 'l':
            len = strlen(optarg);
            if(len && optarg[len-1]=='%') {
               optarg[len-1] = '\0';
               height_is_perc = true;
            }
            else
               height_is_perc = false;
            if(!read_double(optarg, &height, errmsg))
               error(errmsg, c);
            break;

         case 'c':
            centre_type = 'o';
            if(strlen(optarg)==1 && (*optarg=='C' || *optarg=='c'))
               break;   // leave unset centre to indicate centroid to be used
            if(!centre.read(optarg, errmsg))
               error(errmsg, c);
            break;

         case 'p':
            panels = true;
            break;

         case 'm':
            centre_height = true;
            break;

         case 'x':
            hide_edges = true;
            break;

         case 'e':
            col_from_edges = true;
            break;

         case 'o':
            ofile = optarg;
            break;

         default:
            error("unknown command line error");
      }
   }

   if(panels==true && (width!=DEF_VAL || col_from_edges || hide_edges))
      error("cannot set option -p and option -w, -x, -e", 'p');

   if(argc-optind > 1)
      error("too many arguments");

   if(argc-optind == 1)
      ifile=argv[optind];
}


void get_open_edges(const geom_if &geom, vector<vector<int> > &open_edges)
{
   open_edges.clear();
   map<vector<int>, vector<int> > efpairs;
   geom.get_edge_face_pairs(efpairs, false);

   map<vector<int>, vector<int> >::const_iterator ei;
   for(ei=efpairs.begin(); ei!=efpairs.end(); ++ei) {
      if(ei->second.size()==1)     // One faces at an edge
         open_edges.push_back(ei->first);
   }
}

vec3d force_line_plane_intersect(vec3d Q, vec3d n, vec3d P0, vec3d P1)
{
   vec3d pt = line_plane_intersect(Q, n, P0, P1);
   if(!pt.is_set())                     // probably coincident points
      pt = P0 - (P0 - Q).component(n);  // use nearest point to P0 on plane
   return pt;
}


bool leonardo_faces(geom_if &geom_out, const geom_if &geom, double width,
      double height, vec3d centre, double ht_offset, bool hide_edges,
      bool col_from_edges, bool panels, char *errmsg)
{
   if(!geom.faces().size()) {
      if(errmsg)
         strcpy(errmsg, "base polyhedron has no faces and cannot be converted");
      return false;
   }

   const col_geom *cg = dynamic_cast<const col_geom_v *>(&geom);
   col_geom_v *cg_out = dynamic_cast<col_geom_v *>(&geom_out);
   const bool take_colours = cg && cg_out;

   bool is_solid = (height != 0.0);
   const double height_below = height + ht_offset;    // height "below" surface
   const double height_above = ht_offset;             // height "above" surface

   geom_info info(geom);
   vector<vec3d> v_norms;
   if(!centre.is_set())
      v_norms = info.get_vert_norms();
   else {
      v_norms.resize(geom.verts().size());
      for(unsigned int i=0; i<geom.verts().size(); i++)
         v_norms[i] = (geom.verts(i) - centre).unit();
   }

   vector<vector<int> > open_edges;
   get_open_edges(geom, open_edges);

   for(unsigned int i=0; i<geom.faces().size(); i++) {
      const vector<int> &face = geom.faces(i);
      const int f_sz = face.size();
      const vec3d f_norm = geom.face_norm(i).unit();
      const vec3d F = geom.face_cent(i);
      col_val f_col;
      if(take_colours)
         f_col = cg->get_f_col(i);

      vector<vec3d> lines[2];
      for(int j=0; j<f_sz; j++) {
         int v_idx[2];
         vec3d v[8];
         for(int n=0; n<2; n++)
            v_idx[n] = face[(j+n)%f_sz];
         for(int n=0; n<2; n++)
            v[n] = geom.verts(v_idx[n]);
         vec3d perp = (F - nearest_point(F, v[0], v[1])).unit();
         for(int n=0; n<2; n++) {
            vec3d u = (F-v[n]).unit();
            u *= width / u.component(perp).mag();
            lines[n].push_back(v[n] + u);
         }
      }

      int start = geom_out.verts().size();
      if(panels) {
         vector<int> panels[2];
         for(int j=0; j<f_sz; j++)
            for(int p=0; p<1+is_solid; p++)
               panels[p].push_back(start+(1+is_solid)*j+p);
         for(int p=0; p<1+is_solid; p++) {
            int f_idx = geom_out.add_face(panels[p]);
            if(take_colours)
               cg_out->set_f_col(f_idx, f_col);
         }
      }

      for(int j=0; j<f_sz; j++) {
         vec3d v0 = geom.verts(face[j]);
         vec3d v1 = lines_intersection(
                  lines[0][(j-1+f_sz)%f_sz], lines[1][(j-1+f_sz)%f_sz],
                  lines[0][j],               lines[1][j], 0 );

         // vertex above "original" vertex
         geom_out.add_vert(force_line_plane_intersect(
                  F-height_above*f_norm, f_norm, v0, v0+v_norms[face[j]]));
         // vertex above vertex on face at 'width' from neighbouring edges
         if(!panels)
            geom_out.add_vert(force_line_plane_intersect(
                  F-height_above*f_norm, f_norm, v1, v1+f_norm));

         if(is_solid) {
            // vertex below "original" vertex
            geom_out.add_vert(force_line_plane_intersect(
                     F-height_below*f_norm, f_norm, v0, v0+v_norms[face[j]]));
            // vertex below vertex on face at 'width' from neighbouring edges
            if(!panels)
               geom_out.add_vert(force_line_plane_intersect(
                     F-height_below*f_norm, f_norm, v1, v1+f_norm));
         }

         int N = is_solid ? 4 : 2;     // number of vertices added per unit
         if(panels)
            N /= 2;
         bool is_open = false;         // polygon edge is open

         // top face
         if(!panels)
            geom_out.add_face(start+ N*j, start+ N*((j+1)%f_sz),
                              start+ N*((j+1)%f_sz)+1, start+ N*j+1, -1);

         if(is_solid) {
            // bottom face
            if(!panels) {
               geom_out.add_face(start+ N*j+2, start+ N*((j+1)%f_sz)+2,
                                 start+ N*((j+1)%f_sz)+3, start+ N*j+3, -1);
            // inner vertical face
               geom_out.add_face(start+ N*j+1, start+ N*((j+1)%f_sz)+1,
                                 start+ N*((j+1)%f_sz)+3, start+ N*j+3, -1);
            }

            // outer vertical face, if original edge is open
            is_open = find(open_edges.begin(), open_edges.end(),
                  make_edge(face[j], face[(j+1)%f_sz])) != open_edges.end();
            if(is_open) {
               int f_idx = geom_out.add_face(start+ N*j, start+ N*((j+1)%f_sz),
                      start+ N*((j+1)%f_sz)+2-panels, start+ N*j+2-panels, -1);
               if(panels && take_colours)  // easier to colour here
                  cg_out->set_f_col(f_idx, f_col);
            }
         }

         if(!panels) {
            col_val col;       // the colour to use for the faces
            if(take_colours) {
               if(col_from_edges) {
                  vector<vector<int> >::const_iterator ei =
                     find(geom.edges().begin(), geom.edges().end(),
                           make_edge(face[j], face[(j+1)%f_sz]));
                  if(ei != geom.edges().end()) {
                     unsigned int e_idx = ei - geom.edges().begin();
                     col = cg->get_e_col(e_idx);
                  }
               }
               else             // take colour from face
                  col = f_col;
            }
            if(cg_out && hide_edges) {
               cg_out->add_col_edge(start+N*j, start+N*j+1, col_val::invisible);
               if(is_solid)
                  cg_out->add_col_edge(start+N*j+2, start+N*j+3,
                        col_val::invisible);
            }
            if(take_colours && col.is_set()) {
               for(int i=0; i<1 + 2*is_solid + is_open; i++)
                  cg_out->set_f_col(geom_out.faces().size()-1-i, col);
            }
         }
      }
   }

   return true;
}

double min_dist_edge_to_face_cent(geom_if &geom)
{
   double min_dist2 = 1e100;  // large value
   for(unsigned int i=0; i<geom.faces().size(); i++) {
      const vector<int> &face = geom.faces(i);
      const int f_sz = face.size();
      const vec3d F = geom.face_cent(i);
      for(int j=0; j<f_sz; j++) {
         double dist2 = (F - nearest_point(F, geom.verts(face[j]),
                     geom.verts(face[(j+1)%f_sz]))).mag2();
         if(dist2 < min_dist2)
            min_dist2 = dist2;
      }
   }

   return sqrt(min_dist2);
}


int main(int argc, char *argv[])
{
   leo_opts opts;
   opts.process_command_line(argc, argv);

   char errmsg[MSG_SZ];
   col_geom_v geom;
   if(!geom.read(opts.ifile, errmsg))
      opts.error(errmsg);
   if(*errmsg)
      opts.warning(errmsg);

   double width = (opts.width==opts.DEF_VAL) ? opts.def_width : opts.width;
   if(opts.width_is_perc)
      width *= min_dist_edge_to_face_cent(geom) / 100;
   if(fabs(width)<epsilon)
      opts.warning("width is very small, if using default may need to set a value with -w");

   double height = (opts.height==opts.DEF_VAL) ? opts.def_height : opts.height;
   if(opts.height_is_perc)
      height *= width / 100;

   vec3d centre;
   if(opts.centre_type=='o') {
      if(opts.centre.is_set())  // coordinates specified
         centre = opts.centre;
      else             // centroid specified
         centre = geom.centroid();
   }

   col_geom_v geom_out;
   if(!leonardo_faces(geom_out, geom, width, height, centre,
            -height*0.5*opts.centre_height,
            opts.hide_edges, opts.col_from_edges, opts.panels, errmsg))
      opts.error(errmsg);

   if(!geom_out.write(opts.ofile, errmsg))
      opts.error(errmsg);

   return 0;
}




