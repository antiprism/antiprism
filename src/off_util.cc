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
   Name: off_util.cc
   Description: utility processing for OFF file, e.g. merge, orient
   Project: Antiprism - http://www.antiprism.com
*/

#include <string.h>
#include <unistd.h>
#include <math.h>
#include <ctype.h>
#include <string>
#include <vector>
#include <stack>
#include <algorithm>

#include "../base/antiprism.h"

#include "help.h"

using std::string;
using std::vector;
using std::stack;



class pr_opts: public prog_opts {
   public:
      vector<string> ifiles;
      bool orient;
      bool triangulate;
      bool skeleton;
      bool sph_proj;
      bool trunc;
      double trunc_ratio;
      int trunc_v_ord;
      bool edges_to_faces;
      string filt_elems;
      string merge_elems;
      int sig_compare;
      int sig_digits;
      double unzip_frac;
      int unzip_root;
      char unzip_centre;
      
      string ofile;

      pr_opts(): prog_opts("off_util"), orient(false),
                 triangulate(false), skeleton(false),
                 sph_proj(false), trunc(false), edges_to_faces(false),
                 sig_compare(-1), sig_digits(DEF_SIG_DGTS),
                 unzip_frac(100.0), unzip_root(0), unzip_centre('x')
                 {}
      void process_command_line(int argc, char **argv);
      void usage();
};


void pr_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] input_files\n"
"\n"
"Read one or more files in OFF format, combine them into a single file and\n"
"process it. Operations take place in the order listed below. input_files is\n"
"the list of files to process.\n"
"\n"
"Options\n"
"%s"
"  -M <elms> Sort and merge elements whose coordinates are the same to\n"
"            the number of decimal places given by option -l, elems can\n"
"            include v, e, or f to merge vertices, edges and faces, or\n"
"            s to sort without merging (default: vef if option -l set)\n"
"  -l <lim>  minimum distance for unique vertex locations as exponent 1e-lim\n"
"            (default: 8 giving 1e-8, if option -M set)\n"
"  -O        orient the faces (if possible), flip orientation if oriented\n"
"  -T <rat>  truncate vertices by cutting edges at a ratio from each vertex,\n"
"            can also be 'rat,num' to truncate only vertices of order num\n"
"  -E        turn edges into (non-planar) faces\n"
"  -s        skeleton, write the face edges and remove the faces\n"
"  -t        triangulate, divide faces triangles (may add new vertices)\n"
"  -x <elms> remove OFF face elements. The element string is processed in\n"
"            order and can include v, e and f to remove OFF faces with one\n"
"            vertex (vertices), two-vertices (edges) and three or more\n"
"            vertices (faces), and V to remove vertices that are not part\n"
"            of any face or edge.\n"
"  -S        project onto unit sphere centred at origin\n"
"  -u <args> unfold a polyhedron into a net, takes up to three comma separated\n"
"            values for base face index, dihedral fraction (normally 1.0 to\n"
"            -1.0, default: 0.0 flat), a final 'f' centres on centroid of\n"
"            face centres\n"
"  -d <dgts> number of significant digits (default %d) or if negative\n"
"            then the number of digits after the decimal point\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text, DEF_SIG_DGTS);
}

const char *get_help(const char* name)
{
   map<string, const char *> help;
   help["help"] = help_help;
   help["models"] = help_models;
   help["common_polys"] = help_common_polys;
   help["common_polyhedra"] = help_common_polys;
   help["uniform"] = help_uniform;
   help["archimedean"] = help_uniform;
   help["platonic"] = help_uniform;
   help["uniform_duals"] = help_uniform_duals;
   help["ud"] = help_uniform_duals;
   help["catalans"] = help_uniform_duals;
   help["std_polys"] = help_std_polys;
   help["std_polyhedra"] = help_std_polys;
   help["johnson"] = help_johnson;
   help["polygon"] = help_polygon;
   help["pyramid"] = help_polygon;
   help["prism"] = help_polygon;
   help["antiprism"] = help_polygon;
   help["uc"] = help_uniform_compounds;
   help["geo"] = help_geodesic;
   help["geodesic"] = help_geodesic;
   help["sym"] = help_sym;
   help["col_val"] = help_color_val;
   help["col_names"] = help_color_names;
   help["col_map"] = help_color_map;
   char hname[MSG_SZ];
   to_resource_name(hname, name);
   map<string, const char *>::iterator mi = help.find(hname);
   if(mi != help.end())
      return mi->second;
   else
      return 0;
}


   


void pr_opts::process_command_line(int argc, char **argv)
{
   char errmsg[MSG_SZ];
   opterr = 0;
   char c;
   vector<char *> parts;

   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hH:stOd:x:T:ESM:l:u:o:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         case 'H': {
            const char *help_text = get_help(optarg);
            if(help_text) {
               fprintf(stdout, "\n%s\n", help_text);
               exit(0);
            }
            else {
               snprintf(errmsg, MSG_SZ, "no help for '%s' "
                     "(try 'off_util -H help')", optarg);
               error(errmsg, c);
            }
            break;
         }

         case 'O':
            orient = true;
            break;

         case 'T': {
            vector<char *> parts;
            int parts_sz = split_line(optarg, parts, ",");
            if(parts_sz>2)
               error("truncation must be 'ratio' or 'ratio,vert_order'", c);
            if(!read_double(parts[0], &trunc_ratio, errmsg))
               error("truncation ratio is not a number", c);
            trunc_v_ord = 0; // truncate all vertices
            if(parts_sz>1 && !read_int(parts[1], &trunc_v_ord, errmsg)) {
               error("truncation vertex order is not an integer", c);
               if(trunc_v_ord<1)
                  error("truncation vertex order is not positive", c);
            }
            trunc = true;
            break;
         }

         case 'E':
            edges_to_faces = true;
            break;

         case 's':
            skeleton = true;
            break;

         case 't':
            triangulate = true;
            break;

         case 'd':
            if(!read_int(optarg, &sig_digits, errmsg))
               error(errmsg, c);
            break;

         case 'x':
            if(strspn(optarg, "Vvef") != strlen(optarg)) {
               snprintf(errmsg, MSG_SZ, "elements to hide are %s must be V, v, e, or f\n", optarg);
               error(errmsg, c);
            }
            filt_elems=optarg;
            break;

         case 'S':
            sph_proj = true;
            break;

         case 'M':
            if(strspn(optarg, "svef") != strlen(optarg)) {
               snprintf(errmsg, MSG_SZ, "elements to delete are %s must be v, e, f or s\n", optarg);
               error(errmsg, c);
            }
            if(strchr(optarg, 's') && strlen(optarg)>1) {
               error("s is for sorting only, cannot be used with v, e, or f", c);
            }
            if(strspn(optarg, "ef") && !strchr(optarg, 'v')) {
               warning("without v, some orphan vertices may result", c);
            }
            merge_elems=optarg;
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


         case 'u':
            split_line(optarg, parts, ",");
            if(!read_int(parts[0], &unzip_root, errmsg))
               error(errmsg, c);

            unzip_frac = 0.0;
            if(parts.size()>1 && !read_double(parts[1], &unzip_frac, errmsg))
               error(errmsg, c);
            
            unzip_centre = 'x';
            if(parts.size()>2) {
               if(strlen(parts[2])==1 && strchr("f", *parts[2]))
                   unzip_centre = *parts[2];
               else
                  error("invalid unzip centring type, can only be f", c);
            }
            break;

         case 'o':
            ofile = optarg;
            break;

         default:
            error("unknown command line error");
      }
   }

   if(merge_elems=="" && sig_compare>=0)
      merge_elems = "vef";
   if(sig_compare<0 and merge_elems!="")
      sig_compare = 8;

   while(argc-optind)
      ifiles.push_back(string(argv[optind++]));

}


struct tree_face
{
   int idx;
   int start;
   int cur;
   vector<int> orig_cons;
   vector<int> cons;
   tree_face(int index, vector<int> original_cons, int con_idx=-1);
   int get_next_idx();
};

tree_face::tree_face(int index, vector<int> original_cons, int con_idx):
   idx(index), orig_cons(original_cons)
{
   if(con_idx==-1) {
      cur = -1;
      start = 0;
   }
   else {
      vector<int>::iterator vi =
         std::find(orig_cons.begin(), orig_cons.end(), con_idx);
      start = cur = vi-orig_cons.begin();
   }
}

int tree_face::get_next_idx()
{
   if(cur-start+1 >= (int)orig_cons.size())       // all connections seen
      return -1;
   else
      return orig_cons[++cur % orig_cons.size()]; // wrap at end of cons list
}

class unzip_tree
{
   private:
      int root;
      map<int, vector<int> > tree;

   public:
      unzip_tree(): root(-1) {}
      int init_basic(col_geom_v &geom, int first_face);
      int flatten(const col_geom_v &geom, col_geom_v &net_geom, double fract);
};


int unzip_tree::init_basic(col_geom_v &geom, int first_face)
{
   tree.clear();
   root = first_face;
   
   geom_v dual;
   get_dual(geom, dual);
   geom_info info(dual);
   const vector<vector<int> > &f_cons = info.get_vert_cons();
   vector<tree_face> face_list;
   vector<bool> seen(geom.faces().size(), false);

   //start at first vertex of first face
   face_list.push_back(tree_face(root, f_cons[root]));

   int list_sz = -1;
   while(face_list.size()<geom.faces().size()) {  // quit when all faces seen
   //while((int)face_list.size()>list_sz) {  // quit when all faces seen
      list_sz = face_list.size();
      for(int i=0; i<list_sz; i++) {  // faces currently in tree
         int idx;
         while((idx=face_list[i].get_next_idx())>=0) { // face conections
            if(!seen[idx]) {  // unseen face
               seen[idx] = true;
               face_list.push_back(
                     tree_face(idx, f_cons[idx], face_list[i].idx) );
               face_list[i].cons.push_back(idx);
               //geom.set_f_col(idx, level);
               continue;
            }
         }
      }
   }

   for(unsigned int i=0; i<face_list.size(); i++)
      if(face_list[i].cons.size())
         tree[face_list[i].idx] = face_list[i].cons;

   return 1;
}



int unzip_tree::flatten(const col_geom_v &geom, col_geom_v &net_geom,
      double fract=0.0)
{
   net_geom = geom;
   int level=0;
   int cur_face = root;
   stack<int> face_stack;
   face_stack.push(cur_face);
   unsigned int cur_con_num = 0;
   stack<int> con_num_stack;
   con_num_stack.push(cur_con_num);
   map<int, int> cur_idx_map;
   stack<map<int, int> > idx_map_stack;
   idx_map_stack.push(cur_idx_map);
   stack<mat3d> trans_stack;
   trans_stack.push(mat3d());
   stack<vec3d> norm_stack;
   
   while(level>=0) {
      //fprintf(stderr, "new loop: L=%d, cur_con_num=%d, cur_face=%d, prev=%d\n", level, cur_con_num, cur_face, level==0 ? -1 : face_stack.top());
      
      if(cur_con_num==0) { // process new face
         //fprintf(stderr, "\tprocess face %d\n", cur_face);
         vector<int> &face = net_geom.raw_faces()[cur_face];
         map<int, int> new_idx_map;
         vector<int> join_edge;
         int first_mapped_i=-1;
         //mat3d inv_prev_trans = mat3d::inverse(trans_stack.top());
         for(unsigned int i=0; i<face.size(); i++) {
            int idx = face[i];
            int new_idx;
            if(cur_face == root) {  // new vertices to help deletion later
               net_geom.add_vert(net_geom.verts(idx));
               new_idx = net_geom.verts().size()-1;
            }
            else {
               // Find vertices common to current and previous faces, and use
               // duplicated vertices for the non-join vertices in the cur face.
               map<int, int>::const_iterator mi = 
                  idx_map_stack.top().find(face[i]);
               if(mi == idx_map_stack.top().end()) {
                  net_geom.add_vert(trans_stack.top()*net_geom.verts(idx));
                  new_idx = net_geom.verts().size()-1;
               }
               else {
                  new_idx = mi->second;
                  join_edge.push_back(new_idx);
                  // if the join edges don't follow each other they
                  // are in the first and last positions, and so are
                  // sequential in reverse order
                  if(join_edge.size()==1)
                     first_mapped_i = i;
                  else if(join_edge.size()==2) {
                     if(i-first_mapped_i>1)
                        std::swap(join_edge[0], join_edge[1]);
                  }
               }
            }

            face[i] = new_idx;
            new_idx_map[idx] = new_idx;
         }
         idx_map_stack.push(new_idx_map);
       
         mat3d trans;
         vec3d norm = net_geom.face_norm(cur_face);
         if(cur_face!=root) {
            // set join_edge so it is in order for previous face
            const vector<int> &prev_face = net_geom.faces(face_stack.top());
            unsigned int j;
            for(j=0; j<prev_face.size(); j++)
               if(prev_face[j] == join_edge[0])
                  break;
            if(prev_face[(j+1)%prev_face.size()] == join_edge[1]) {
               std::reverse(face.begin(), face.end());
               norm *= -1.0;
            }
            else
               std::swap(join_edge[0], join_edge[1]);

            
            vec3d axis = net_geom.edge_vec(join_edge).unit();
            double ang = angle_around_axis(norm, norm_stack.top(), axis);
            ang *= (1-fract);
            mat3d rot = mat3d::rot(axis, ang);
            vec3d offset = net_geom.verts(join_edge[0]);
            trans = mat3d::transl(offset) * rot * mat3d::transl(-offset);

            for(unsigned int i=0; i<face.size(); i++) {
               if(face[i]!=join_edge[0] && face[i]!=join_edge[1])
                  net_geom.raw_verts()[face[i]] = 
                     trans * net_geom.verts(face[i]);
            }
            norm = rot * norm;  // normal needs rotating because face rotated
         }
         norm_stack.push(norm);
         trans_stack.push(trans * trans_stack.top());
      }

      if(cur_con_num == tree[cur_face].size()) {  // finish at this level
         //fprintf(stderr, "\ttree.size()=%d, tree[%d].size()=%d\n", tree.size(), cur_face, tree[cur_face].size());
         //fprintf(stderr, "\tfinish at level\n");
         if(level) {
            cur_face = face_stack.top();
            face_stack.pop();
            cur_con_num = con_num_stack.top() + 1;  // reset and increment
            con_num_stack.pop();
            idx_map_stack.pop();
            norm_stack.pop();
            trans_stack.pop();
         }
         level--;
         continue;
      }

      unsigned int next_idx = tree[cur_face][cur_con_num];
      //fprintf(stderr, "\tgo out a level\n");
      level++;                         
      face_stack.push(cur_face);  
      cur_face = next_idx;
      con_num_stack.push(cur_con_num);
      cur_con_num = 0;
   }

   vector<int> del_verts(geom.verts().size());
   for(unsigned int i=0; i<geom.verts().size(); i++)
      del_verts[i] = i;
   net_geom.delete_verts(del_verts);
   return 1;
}


int unzip_poly(col_geom_v &geom, int root, double fract, char centring,
      char *errmsg)
{
   geom_info info(geom);
   if(!info.is_polyhedron()) {
      strcpy(errmsg, "input not a polyhedron (temporary restriction)");
      return 0;
   }
   if(info.num_parts()>1) {
      strcpy(errmsg, "input not connected (temporary restriction)");
      return 0;
   }
   if(root<0 || root>=(int)geom.faces().size()) {
      sprintf(errmsg, "root face '%d' is not a valid face index number", root);
      return 0;
   }
   if(!strchr("fx", centring)) {
      sprintf(errmsg, "invalid centring type '%c', must be f or x", centring);
      return 0;
   }

   unzip_tree tree;
   tree.init_basic(geom, root);
   col_geom_v net_geom;
   tree.flatten(geom, net_geom, fract);

   if(centring=='f') {
      vector<vec3d> f_cents;
      net_geom.face_cents(f_cents);
      net_geom.transform(mat3d::transl(-centroid(f_cents)));
   }

   geom = net_geom;
   return 1;
}


void make_edges_to_faces(geom_if &geom)
{
   col_geom_v egeom;
   edges_to_faces(geom, egeom, true);
   geom.clear_all();
   geom.append(egeom);
}

void triangulate_faces(geom_if &geom)
{
   geom.add_missing_impl_edges();
   geom.triangulate(col_val::invisible);
}

void make_skeleton(geom_if &geom)
{
   geom.add_missing_impl_edges();
   geom.clear_faces();
}


void filter(geom_if &geom, const char *elems)
{
   col_geom_v *cgv = dynamic_cast<col_geom_v *>(&geom);
   for(const char *p=elems; *p; p++) {
      switch(*p) {
         case 'V':
            geom.delete_verts(geom.get_info().get_free_verts());
            break;
         case 'v':
            cgv->clear_v_cols();
            break;
         case 'e':
            geom.clear_edges();
            break;
         case 'f':
            geom.clear_faces();
            break;
      }
   }
}

void proj_onto_sphere(geom_if &geom)
{
   vector<vec3d> &verts = *geom.get_verts();
   for(unsigned int i=0; i<verts.size(); i++)
      verts[i].to_unit();
}


void process_file(col_geom_v &geom, pr_opts opts)
{
   char errmsg[MSG_SZ]="";
   if(opts.merge_elems!="") {
      sort_merge_elems(geom, opts.merge_elems.c_str(),
            pow(10, -opts.sig_compare));
   }
      
   if(opts.orient) {
      geom.orient_reverse();
      //orient_faces(geom);
      geom.orient();
      geom_info rep(geom);
      if(!rep.is_oriented()) {
         snprintf(errmsg, MSG_SZ, "input file contains a non-orientable geometry");
         opts.warning(errmsg, 'O');
      }
   }
   if(opts.trunc)
      truncate_verts(geom, opts.trunc_ratio, opts.trunc_v_ord);
   if(opts.unzip_frac!=100.0) {
      if(!unzip_poly(geom, opts.unzip_root, opts.unzip_frac, 
            opts.unzip_centre, errmsg))
         opts.error(errmsg, 'u');
   }

   if(opts.edges_to_faces)
      make_edges_to_faces(geom);
   if(opts.triangulate)
      triangulate_faces(geom);
   if(opts.skeleton)
      make_skeleton(geom);
   if(opts.filt_elems!="")
      filter(geom, opts.filt_elems.c_str());
   if(opts.sph_proj)
      proj_onto_sphere(geom);
}



int main(int argc, char *argv[])
{
   pr_opts opts;
   opts.process_command_line(argc, argv);
   if(!opts.ifiles.size())
      opts.ifiles.push_back("");

   char errmsg[MSG_SZ]="";
   col_geom_v geoms;
   for(unsigned int i=0; i<opts.ifiles.size(); i++) {
      col_geom_v geom;
      if(!geom.read(opts.ifiles[i], errmsg))
         opts.error(errmsg);
      if(*errmsg)
         opts.warning(errmsg);
      geoms.append(geom);
   }
   
   process_file(geoms, opts);

   if(!geoms.write(opts.ofile, errmsg, opts.sig_digits))
      opts.error(errmsg);

   return 0;
}
   

