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
#include <math.h>
#include <ctype.h>
#include <string>
#include <vector>
#include <stack>
#include <algorithm>
#include <memory>

#include "../base/antiprism.h"

#include "help.h"

using std::string;
using std::vector;
using std::stack;
using std::auto_ptr;



class pr_opts: public prog_opts {
   public:
      vector<string> ifiles;
      int orient;
      unsigned int triangulate_rule;
      bool skeleton;
      bool sph_proj;
      bool trunc;
      double trunc_ratio;
      int trunc_v_ord;
      bool edges_to_faces;
      bool geometry_only;
      string filt_elems;
      vector<string> del_elems;
      vector<string> add_elems;
      bool close;
      col_val close_col;
      string merge_elems;
      int blend_type;
      double epsilon;
      int sig_digits;
      double unzip_frac;
      int unzip_root;
      char unzip_centre;
      char unzip_z_align;
      
      string ofile;

      pr_opts(): prog_opts("off_util"), orient(0),
                 triangulate_rule(0), skeleton(false),
                 sph_proj(false), trunc(false), edges_to_faces(false),
                 geometry_only(false), close(false), blend_type(3),
                 epsilon(0), sig_digits(DEF_SIG_DGTS),
                 unzip_frac(100.0), unzip_root(0), unzip_centre('x'),
                 unzip_z_align(false)
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
"            the number of decimal places given by option -l.  elems can\n"
"            include: v - vertices, e - edges, f - faces,  a - all (vef),\n"
"            b - bond (merge 've' and delete any face coincident with another),\n"
"            s - sort without merging\n"
"  -b <opt>  merge blend color. first=1, last=2, rgb=3, ryb=4 (default: 3)\n"
"  -l <lim>  minimum distance for unique vertex locations as negative exponent\n"
"               (default: %d giving %.0e)\n"
"  -O <opt>  orient the faces first (if possible) then for volume\n"
"            positive, negative, reverse, or flip which reverses the\n"
"            orientation of the model as it was input\n"
"  -T <rat>  truncate vertices by cutting edges at a ratio from each vertex,\n"
"            can also be 'rat,num' to truncate only vertices of order num\n"
"  -E        turn edges into (non-planar) faces\n"
"  -s        skeleton, write the face edges and remove the faces\n"
"  -t <disp> triangulate, include face parts according to winding number\n"
"            from: odd, nonzero, positive, negative, triangulate (synonym\n"
"            for nonzero)\n"
"  -g        geometry only, remove all colours, remove all two-vertex faces\n"
"            (edges) that are also a face edge\n"
"  -x <elms> remove OFF face elements. The element string is processed in\n"
"            order and can include v, e, f to remove OFF faces with one\n"
"            vertex (vertices), two-vertices (edges) and three or more\n"
"            vertices (faces), V to remove vertices that are not part\n"
"            of any face or edge, E to remove two-vertex faces (edges) that\n"
"            are also a face edge.\n"
"  -D <list> delete a list of elements, list starts with element letter\n"
"            (f,e, v, deleted in that order, only one list per element),\n"
"            followed by an index number list, given as index ranges\n"
"            separated by commas, range can be one number or two numbers\n"
"            separated by a hyphen (default range numbers: 0 and largest index)\n"
"  -A <elem> add element, elem is element letter (v, e, f), followed by\n"
"            element data, optionally followed by ':' and a colour. Data is\n"
"               v: three comma separated coordinates\n"
"               e: a comma separated list of index numbers, joined as a ring\n"
"               f: a comma separated list of index numbers\n"
"            negative index numbers are relative to the end of the vertex\n"
"            list, last vertex is -1 (useful to refer to added vertices.)\n"
"  -c <col>  close polyhedron, each hole converted to a face with colour col,\n"
"            holes having a vertex with more than two open edges are not filled\n"
"  -S        project onto unit sphere centred at origin\n"
"  -u <args> unfold a polyhedron into a net, takes up to three comma separated\n"
"            values for base face index, dihedral fraction (normally 1.0 to\n"
"            -1.0, default: 0.0 flat), and final option letters: 'f' centre\n"
"            on centroid of face centres, 'z' align base face normal to z_axis.\n"
"  -d <dgts> number of significant digits (default %d) or if negative\n"
"            then the number of digits after the decimal point\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text, int(-log(::epsilon)/log(10) + 0.5), ::epsilon, DEF_SIG_DGTS);
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
   help["sym_models"] = help_sym_models;
   help["color"] = help_color;
   help["colour"] = help_color;
   help["col_val"] = help_color_val;
   help["col_names"] = help_color_names;
   help["col_map"] = help_color_map;
   help["symmetry"] = help_symmetry;
   help["bowers"] = help_bowers;
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
   col_geom_v adding;
   string arg_id;

   int sig_compare = INT_MAX;

   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hH:st:O:d:x:D:A:c:gT:ESM:b:l:u:o:")) != -1 ) {
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
            arg_id = get_arg_id(optarg,"positive=1|negative=2|reverse=3|flip=4",
                  argmatch_add_id_maps, errmsg);
            if(arg_id=="")
               error(errmsg);
            orient = atoi(arg_id.c_str());
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
            arg_id = get_arg_id(optarg,
                  "odd|nonzero|positive|negative|triangulate=1",
                  argmatch_default, errmsg);
            if(arg_id=="")
               error(msg_str("invalid winding rule '%s'", optarg).c_str(), c);
            triangulate_rule = TESS_WINDING_ODD + atoi(arg_id.c_str());
            break;

         case 'd':
            if(!read_int(optarg, &sig_digits, errmsg))
               error(errmsg, c);
            break;

         case 'x':
            if(strspn(optarg, "vefVE") != strlen(optarg)) {
               snprintf(errmsg, MSG_SZ, "elements to hide are %s, can only include vefVE\n", optarg);
               error(errmsg, c);
            }
            filt_elems=optarg;
            break;

         case 'D':
            del_elems.push_back(optarg);
            break;

         case 'A':
            add_elems.push_back(optarg);
            break;

         case 'g':
            geometry_only = true;
            break;

         case 'c':
            if(!close_col.read(optarg, errmsg))
               error(errmsg, c);
            close = true;
            break;

         case 'S':
            sph_proj = true;
            break;

         case 'M':
            if(strspn(optarg, "svefab") != strlen(optarg)) {
               snprintf(errmsg, MSG_SZ, "elements to merge are %s must be v, e, f, a, b or s\n", optarg);
               error(errmsg, c);
            }
            if(strchr(optarg, 's') && strlen(optarg)>1) {
               error("s is for sorting only, cannot be used with v, e, or f", c);
            }
            if(strchr(optarg, 'a') && strlen(optarg)>1) {
               error("a includes vef, and must be used alone", c);
            }
            if(strchr(optarg, 'b') && strlen(optarg)>1) {
               error("b includes vef, and must be used alone", c);
            }
            if(strspn(optarg, "ef") && !strchr(optarg, 'v')) {
               warning("without v, some orphan vertices may result", c);
            }
            merge_elems=optarg;
            if(merge_elems == "a")
               merge_elems = "vef";
            break;

         case 'b':
         {
            string id = get_arg_id(optarg, "first=1|last=2|rgb=3|ryb=4", argmatch_add_id_maps, errmsg);
            if(id=="")
               error(errmsg);
            blend_type = atoi(id.c_str());
            break;
         }

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


         case 'u':
            split_line(optarg, parts, ",");
            if(!read_int(parts[0], &unzip_root, errmsg))
               error(errmsg, c);

            unzip_frac = 0.0;
            if(parts.size()>1 && !read_double(parts[1], &unzip_frac, errmsg))
               error(errmsg, c);
            
            unzip_centre = 'x';
            if(parts.size()>2) {
               if(strspn(parts[2], "zf") != strlen(parts[2])) {
                  snprintf(errmsg, MSG_SZ, "unzip options are '%s' must include"
                        "only f, u\n", parts[2]);
                  error(errmsg, c);
               }
               if(strchr(parts[2], 'f'))
                   unzip_centre = 'f';
               if(strchr(parts[2], 'z'))
                   unzip_z_align = true;
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

   epsilon = (sig_compare != INT_MAX) ? pow(10, -sig_compare) : ::epsilon;
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
      bool unzip_z_align, char *errmsg)
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

   if(unzip_z_align)
      geom.transform(mat3d::rot(geom.face_norm(root), vec3d::Z));
      
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

void triangulate_faces(geom_if &geom, unsigned int winding_rule)
{
   geom.triangulate(col_val::invisible, winding_rule);
}

void make_skeleton(geom_if &geom)
{
   geom.add_missing_impl_edges();
   geom.clear_faces();
}

// Roger Kaufman
void clear_unneeded_explicit_edges(geom_if &geom)
{
   vector<vector<int> > implicit_edges;
   geom.get_impl_edges(implicit_edges);

   const vector<vector<int> > &edges = geom.edges();

   vector<int> deleted_edges;

   for(unsigned int i=0; i<edges.size(); i++) {
      int answer = find_edge_in_edge_list(implicit_edges, edges[i]);
      if (answer > -1)
         deleted_edges.push_back(i);
   }

   geom.delete_edges(deleted_edges);
}

void geometry_only(col_geom_v &geom)
{
   geom.clear_v_cols();
   geom.clear_e_cols();
   geom.clear_f_cols();
   clear_unneeded_explicit_edges(geom);
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
         case 'E':
            clear_unneeded_explicit_edges(geom);
            break;
         case 'f':
            geom.clear_faces();
            break;
      }
   }
}


bool get_del_element_list(geom_if &geom, const string &elem, vector<vector<int> > &elem_lists, char *errmsg)
{
   *errmsg = '\0';
   if(!elem.size())
      return true;

   auto_ptr<char> elem_str(copy_str(elem.c_str()));
   if(!elem_str.get()) {
      strcpy(errmsg, "could not allocate memory");
      return false;
   }

   const char *elem_type_strs[] = {"vertex", "edge", "face"};
   char elem_type_char = *elem_str;
   int elem_type;
   int elems_sz;
   if(elem_type_char=='v') {
      elem_type = 0;
      elems_sz = geom.verts().size();
   }
   else if(elem_type_char=='e') {
      elem_type = 1;
      elems_sz = geom.edges().size();
   }
   else if(elem_type_char=='f') {
      elem_type = 2;
      elems_sz = geom.faces().size();
   }
   else {
      strcpy(errmsg, msg_str("invalid element type '%c', "
               "should be v, e, or f", elem_type_char).c_str());
      return false;
   }

   if(elem_lists[elem_type].size()) {
      strcpy(errmsg, msg_str("list for %s elements already given",
               elem_type_strs[elem_type]).c_str());
      return false;
   }

   char errmsg2[MSG_SZ];
   if(!read_idx_list(elem_str.get()+1, elem_lists[elem_type], elems_sz,
            false, errmsg2)) {
      strcpy(errmsg, msg_str("list for %s elements: %s",
               elem_type_strs[elem_type], errmsg2).c_str());
      return false;
   }

   return true;
}

bool delete_elements(col_geom_v &geom, vector<string> del_elems, char *errmsg)
{
   vector<vector <int> > elem_lists(3);
   vector<string>::const_iterator vi;
   for(vi=del_elems.begin(); vi!=del_elems.end(); ++vi) {
      if(!get_del_element_list(geom, *vi, elem_lists, errmsg))
         return false;
   }

   geom.delete_faces(elem_lists[2]);
   geom.delete_edges(elem_lists[1]);
   geom.delete_verts(elem_lists[0]);
   
   return true;
}


bool add_element(col_geom_v &geom, const string &elem, char *errmsg)
{
   *errmsg = '\0';
   char errmsg2[MSG_SZ];
   if(!elem.size())
      return true;

   auto_ptr<char> elem_str(copy_str(elem.c_str()));
   if(!elem_str.get()) {
      strcpy(errmsg, "could not allocate memory");
      return false;
   }

   char elem_type_char = *elem_str;
   const char *elem_type_str;
   if(elem_type_char=='v')
      elem_type_str = "vertex";
   else if(elem_type_char=='e')
      elem_type_str = "edge";
   else if(elem_type_char=='f')
      elem_type_str = "face";
   else {
      strcpy(errmsg, msg_str("invalid element type '%c', "
               "should be v, e, or f", elem_type_char).c_str());
      return false;
   }

   vector<char *> parts;
   int num_parts = split_line(elem_str.get()+1, parts, ":");
   if(num_parts>2) {
      strcpy(errmsg, "more than one ':' used");
      return false;
   }
   if(num_parts==0 || parts[0]=='\0') {
      strcpy(errmsg, msg_str("'%s': no element data", elem.c_str()).c_str());
      return false;
   }

   col_val col;
   if(num_parts>1 && parts[1]!='\0') {
      if(!col.read(parts[1], errmsg2)) {
        strcpy(errmsg, msg_str("%s colour '%s': %s",
                 elem_type_str, parts[1], errmsg2).c_str());
        return false;
      }
   }

   if(elem_type_char=='v') {
      vec3d vec;
      if(!vec.read(parts[0], errmsg2)) {
         strcpy(errmsg, msg_str("vertex coordinates '%s': %s", parts[0],
                  errmsg2).c_str());
         return false;
      }
      geom.add_col_vert(vec, col);
   }
   else if(elem_type_char=='e' || elem_type_char=='f') {
      vector<int> idx_list;
      if(!read_int_list(parts[0], idx_list, errmsg2)) {
         strcpy(errmsg, msg_str("%s index numbers: '%s': %s",
                              elem_type_str, parts[0], errmsg2).c_str());
         return false;
      }
      int list_sz = idx_list.size();
      for(int i=0; i<list_sz; i++) {
         int idx = idx_list[i];
         if(idx<0)
            idx += geom.verts().size();
         if(idx<0 || idx>=(int)geom.verts().size()) {
            strcpy(errmsg, msg_str("%s index numbers: '%s': %d out of range",
                              elem_type_str, parts[0], idx_list[i]).c_str());
            return false;
         }
         idx_list[i] = idx;
      }
      if(elem_type_char=='e') {
         if(list_sz<2) {
            strcpy(errmsg, msg_str("%s index numbers: '%s': "
                     "must give at least two numbers",
                              elem_type_str, parts[0]).c_str());
            return false;
         }
         if(list_sz==2)
            geom.add_col_face(idx_list, col);
         else {
            for(int i=0; i<list_sz; i++)
               geom.add_col_edge(
                     make_edge(idx_list[i], idx_list[(i+1)%list_sz]), col);
         }
      }
      if(elem_type_char=='f') {
         if(!list_sz) {
            strcpy(errmsg, msg_str("%s index numbers: '%s': "
                     "no index numbers given",
                              elem_type_str, parts[0]).c_str());
            return false;
         }
         geom.add_col_face(idx_list, col);
      }
   }
   return true;
}

bool add_elements(col_geom_v &geom, vector<string> add_elems, char *errmsg)
{
   vector<string>::const_iterator vi;
   for(vi=add_elems.begin(); vi!=add_elems.end(); ++vi) {
      if(!add_element(geom, *vi, errmsg))
         return false;
   }
   return true;
}

void close_poly(col_geom_v &geom, col_val col)
{
   int orig_faces_sz = geom.faces().size();
   close_poly_basic(geom);
   if(col.is_set()) {
      for(unsigned int i=orig_faces_sz; i<geom.faces().size(); i++)
         geom.set_f_col(i, col);
   }
}


void process_file(col_geom_v &geom, pr_opts opts)
{
   char errmsg[MSG_SZ]="";
   if(opts.merge_elems!="") {
      if(opts.merge_elems=="b") {
         sort_merge_elems(geom, "ve", opts.blend_type, opts.epsilon);
         col_geom_v tmp = geom;
         vector<map<int, set<int> > > equiv_elems;
         check_congruence(geom, tmp, &equiv_elems, opts.epsilon);
         vector<int> del_faces;
         map<int, set<int> >::iterator mi;
         for(mi=equiv_elems[2].begin(); mi!=equiv_elems[2].end(); ++mi) {
            set<int>::iterator si;
            si = mi->second.begin();
            ++si;
            if(si!=mi->second.end() && *si<(int)(geom.faces().size())) {
               for(si=mi->second.begin(); si!=mi->second.end(); ++si) {
                  del_faces.push_back(*si);
               }
            }
         }
         geom.delete_faces(del_faces);
      }
      else
         sort_merge_elems(geom, opts.merge_elems, opts.blend_type, opts.epsilon);
   }

   if(opts.orient) {
      if(!geom.orient(opts.orient, errmsg))
         opts.error(errmsg, 'O');
      if(*errmsg)
         opts.warning(errmsg, 'O');
   }
   if(opts.trunc)
      truncate_verts(geom, opts.trunc_ratio, opts.trunc_v_ord);
   if(opts.unzip_frac!=100.0) {
      if(!unzip_poly(geom, opts.unzip_root, opts.unzip_frac, 
            opts.unzip_centre, opts.unzip_z_align, errmsg))
         opts.error(errmsg, 'u');
   }

   if(opts.edges_to_faces)
      make_edges_to_faces(geom);
   if(opts.triangulate_rule)
      triangulate_faces(geom, opts.triangulate_rule);
   if(opts.skeleton)
      make_skeleton(geom);
   if(opts.geometry_only)
      geometry_only(geom);
   if(opts.filt_elems!="")
      filter(geom, opts.filt_elems.c_str());
   if(opts.add_elems.size())
      if(!add_elements(geom, opts.add_elems, errmsg))
         opts.error(errmsg, 'A');
   if(opts.del_elems.size())
      if(!delete_elements(geom, opts.del_elems, errmsg))
         opts.error(errmsg, 'D');
   if(opts.close)
      close_poly(geom, opts.close_col);
   if(opts.sph_proj)
      project_onto_sphere(geom);
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
   

