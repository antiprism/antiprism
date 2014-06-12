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
   Name: off_align.cc
   Description: position one OFF file relative to another
   Project: Antiprism - http://www.antiprism.com
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <string>
#include <vector>

#include "../base/antiprism.h"

using std::string;
using std::vector;
using std::map;

struct bond_brick
{
   col_geom_v geom;
   int align_type;
   vector<int> bond;
};


class bond_base
{
   private:
      col_geom_v base;
      sch_sym sym;
      vector<bond_brick> bricks;

   public:
      enum { align_verts, align_faces, align_faces_merge };
      enum { out_default=0, out_brick, out_base_brick, out_brick_base};
      bond_base(geom_if &bas): base(bas) {}
      bool set_sym(const char *sym_str, char *errmsg=0);
      bool add_brick(char type, const string &brick_str, char *errmsg=0);
      bool bond_all(geom_if &geom_out, int out_type, char *errmsg=0);
};

bool bond_base::set_sym(const char *sym_str, char *errmsg)
{
   char errmsg2[MSG_SZ];
   char sym_cpy[MSG_SZ]; // big enough for normal use
   strncpy(sym_cpy, sym_str, MSG_SZ);

   sch_sym full_sym(base);
   vector<char *> parts;
   split_line(sym_cpy, parts, ",");
   if(parts.size()==0 || parts.size()>2) {
      if(errmsg)
         sprintf(errmsg, "argument should have 1 or 2 comma separated parts");
      return false;
   }
           
   sch_sym sub_sym;
   if(strncmp(parts[0], "full", strlen(parts[0]))==0) {
      sub_sym = full_sym;
   }
   else if(!sub_sym.init(parts[0], mat3d(), errmsg2)) {
      if(errmsg)
         snprintf(errmsg, MSG_SZ, "sub-symmetry type: %s", errmsg2);
      return false;
   }

   int sub_sym_conj = 0;
   if(parts.size()>1) {
      if(!read_int(parts[1], &sub_sym_conj, errmsg2)) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "sub-symmetry conjugation number: %s",
                  errmsg2);
         return false;
      }
   }

   sym = full_sym.get_sub_sym(sub_sym, sub_sym_conj, errmsg2);
   if(sym.get_sym_type() == sch_sym::unknown) {
      if(errmsg)
         snprintf(errmsg, MSG_SZ, "sub-symmetry: %s", errmsg2);
      return false;
   }

   return true;
}


int get_aligns(const vector<double> &angs0, const vector<double> &angs1,
      vector<int> &aligns)
{
   vector<double> test = angs1;
   reverse(test.begin(), test.end());
   for(unsigned int i=0; i<angs0.size(); i++) {
      if(cmp_face_angles(angs0, test)==0)
         aligns.push_back(-1-i);
      rotate(test.begin(), test.begin()+1, test.end());
   }

   reverse(test.begin(), test.end());
   for(unsigned int i=0; i<angs0.size(); i++) {
      if(cmp_face_angles(angs0, test)==0)
         aligns.push_back(i);
      rotate(test.begin(), test.begin()+1, test.end());
   }
   return aligns.size();
}

 
bool bond_base::add_brick(char type, const string &brick_str, char *errmsg)
{
   if(errmsg)
      *errmsg = '\0';
   char errmsg2[MSG_SZ];
   char *str = copy_str(brick_str.c_str());
   char *first_comma = strchr(str, ',');
   if(first_comma)
      *first_comma = '\0';
   
   bricks.push_back(bond_brick());
   bond_brick &brick = bricks.back();
   
   if(!*str)
      brick.geom = base;
   else {
      if(!brick.geom.read(str, errmsg2)) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "brick geometry: %s", errmsg2);
         bricks.pop_back();
         free(str);
         return false;
      }
   }

   if(first_comma) {
      if(!read_int_list(first_comma+1, brick.bond, errmsg2, true)) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "bond values: %s", errmsg2);
         bricks.pop_back();
         free(str);
         return false;
      }
   }

   if(type=='v') {
      brick.align_type = align_verts;
      int n = brick.bond.size();
      if(n!=2 && n!=4 && n!=6) {
         if(errmsg)
            strcpy(errmsg, msg_str(
                     "must give 2, 4 or 6 vertices (%d were given)",n).c_str());
         bricks.pop_back();
         free(str);
         return false;
      }
      if(n>2) {
         for(int i=0; i<n/2; i++) {
            if(brick.bond[i]==brick.bond[(i+1)%(n/2)]) {
               if(errmsg)
                  strcpy(errmsg, "repeated vertex index in base");
               bricks.pop_back();
               free(str);
               return false;
            }
            if(brick.bond[n/2 + i]==brick.bond[n/2 + (i+1)%(n/2)]) {
               if(errmsg)
                  strcpy(errmsg, "repeated vertex index in brick");
               bricks.pop_back();
               free(str);
               return false;
            }
         }
      }
      for(unsigned int i=0; i<brick.bond.size(); i++) {
         const vector<vec3d> &vs = (i<brick.bond.size()/2) ?
            base.verts() : brick.geom.verts();
         if(brick.bond[i]<0 || brick.bond[i]>=(int)vs.size()) {
            if(errmsg)
               snprintf(errmsg, MSG_SZ,
                        "bond values: vertex %d (position %d) is out of bounds",
                        brick.bond[i], i+1);
            bricks.pop_back();
            free(str);
            return false;
         }
      }
   }
   else if(type=='f' || type=='F') {
      brick.align_type = (type=='f') ? align_faces : align_faces_merge;
      if(brick.bond.size()>3) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ,
                  "up to three arguments can be given (%lu were given)",
                  (unsigned long)brick.bond.size());
         bricks.pop_back();
         free(str);
         return false;
      }
      brick.bond.resize(3, 0);
      int f0 = brick.bond[0];
      int f1 = brick.bond[1];
      if(f0<0 || f0>=(int)base.faces().size()) {
         if(errmsg) {
            if(base.faces().size())
               snprintf(errmsg, MSG_SZ, "invalid base face '%d', "
                     "last face is %d", f0, (int)base.faces().size()-1);
            else
               snprintf(errmsg, MSG_SZ, "base has no faces");
         }
         bricks.pop_back();
         free(str);
         return false;
      }
      if(f1<0 || f1>=(int)brick.geom.faces().size()) {
         if(errmsg) {
            if(brick.geom.faces().size())
               strcpy(errmsg, msg_str(
                     "invalid brick face '%d', last face is %d",
                     f1, (int)brick.geom.faces().size()-1).c_str());
            else
               snprintf(errmsg, MSG_SZ, "brick has no faces");
         }
         bricks.pop_back();
         free(str);
         return false;
      }
      if(brick.align_type==align_faces_merge &&
            base.faces(f0).size()!= brick.geom.faces(f1).size() ) {
         if(errmsg)
            strcpy(errmsg, "faces to bond have different number of sides");
      }
      vector<int> aligns;
      vector<double> angs0, angs1;
      base.get_info().face_angles_lengths(f0, angs0);
      brick.geom.get_info().face_angles_lengths(f1, angs1);
      if(!get_aligns(angs0, angs1, aligns)) {
         if(errmsg)
            sprintf(errmsg,
                  "base and brick bonding faces are not the same shape\n");
            return false;
      }
      if(brick.bond[2]>=(int)aligns.size()) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ,
                  "bond selection number is %d, last selection is %d\n",
                  brick.bond[2], (int)aligns.size()-1);
            return false;
      }
      vector<int> &face = brick.geom.raw_faces()[f1];
      int sel = aligns[brick.bond[2]];
      if(sel>=0)
         rotate(face.begin(), face.begin()+sel, face.end());
      else {
         brick.geom.orient_reverse();
         rotate(face.begin(), face.begin()+(-sel-1), face.end());
      }
   }
   return true;
}


bool bond_base::bond_all(geom_if &geom_out, int out_type, char *errmsg)
{
   geom_out.clear_all();
   double base_rad = bound_sphere(base.verts()).get_radius();
   // check for errors before making any transfromations
   vector<bond_brick>::iterator bi;
   for(bi=bricks.begin(); bi!=bricks.end(); ++bi) {
      if(bi->align_type==align_faces_merge && 
            (out_type==out_brick || out_type==out_brick_base) ) {
         if(errmsg)
            sprintf(errmsg,
                  "output type is %s and not compatible with a merged brick",
                  (out_type==out_brick) ? "brick only" : "brick then base");
         return false;
      }
   }

   bool has_merged_brick = false;
   col_geom_v base_out;
   col_geom_v brick_out;
   vector<col_geom_v> merge_bricks;
   map<int, vector<int> > face_params;
   for(bi=bricks.begin(); bi!=bricks.end(); ++bi) {
      bond_brick &brick = *bi;
      if(brick.align_type==align_verts) {
         vector<vector<vec3d> > pts(2);
         for(unsigned int i=0; i<brick.bond.size(); i++) {
            const vector<vec3d> &vs = (i<brick.bond.size()/2) ?
               base.verts() : brick.geom.verts();
            pts[!(i<brick.bond.size()/2)].push_back(vs[brick.bond[i]]);
         }
         if(pts[0].size()>1) {
            mat3d scl = mat3d::scale(
                  (pts[0][0]-pts[0][1]).mag()/(pts[1][0]-pts[1][1]).mag()); 
            brick.geom.transform(scl);
            for(unsigned int i=0; i<pts[1].size(); i++)
               pts[1][i] = scl * pts[1][i];
         }
         brick.geom.transform(mat3d::alignment(pts[1],pts[0]));
         brick_out.append(brick.geom);
      }
      else if(brick.align_type==align_faces ||
              brick.align_type==align_faces_merge) {
         int f0 = brick.bond[0];
         int f1 = brick.bond[1];
           
         double e0 = base.edge_len(base.faces(f0)); // first edge
         double e1 = brick.geom.edge_len(
               make_edge(brick.geom.faces(f1, 0), brick.geom.faces(f1, 1)));
         brick.geom.transform(mat3d::scale(e0/e1) );
        
         if(bi->align_type==align_faces) {
            face_bond_direct(base, brick.geom, f0, f1, false);
            brick_out.append(brick.geom);
         }
         else {
            if(!has_merged_brick) {
               has_merged_brick = true;
               base_out = base;
            }
            double model_rad = base_rad +
               bound_sphere(bi->geom.verts()).get_radius();
            merge_bricks.push_back(brick.geom);
            t_set ts = sym.get_trans();
            if(!ts.size())        // set to unit if not set
               ts.add(mat3d());
            for(set<mat3d>::const_iterator si=ts.begin(); si!=ts.end(); ++si) {
               vector<vector<int> > elem_maps;
               get_congruence_maps(base, *si, elem_maps, model_rad*sym_eps);
               const int f0_map = elem_maps[2][f0];
               const vector<int> &f0_mapface = base.faces(f0_map);
               int f0_sz = base.faces(f0).size();
               map<int, vector<int> >::iterator vi = face_params.find(f0_map);
               bool direct = iso_type(*si).is_direct();
               if(vi==face_params.end() || (!(vi->second)[2] && direct) ) {
                  int map_offset;
                  const int idx0 = base.faces(f0)[0];
                  for(map_offset=0; map_offset<f0_sz; map_offset++) {
                     if(elem_maps[0][idx0]==f0_mapface[map_offset])
                        break;
                  }
                  const int idx1 = base.faces(f0)[1];
                  bool rev = (elem_maps[0][idx1]!=
                        f0_mapface[(map_offset+1)%f0_sz]); 
                  face_params[f0_map] = vector<int>(4);
                  face_params[f0_map][0] = f1;
                  face_params[f0_map][1] = map_offset;
                  face_params[f0_map][2] = merge_bricks.size()-1;
                  face_params[f0_map][3] = direct;
                  face_params[f0_map][4] = rev;
               }
            }
         }
      }
   }

   if(face_params.size()) {
      map<int, vector<int> >::reverse_iterator mi;  
      for(mi = face_params.rbegin(); mi!=face_params.rend(); ++mi) {
         col_geom_v bgeom = merge_bricks[mi->second[2]];
         vector<int> &brick_f = bgeom.raw_faces()[mi->second[0]];
         rotate(brick_f.begin(), brick_f.begin()+mi->second[1], brick_f.end());
         if(!mi->second[3])  // symmetry was indirect
            bgeom.transform(mat3d::inversion());
         if(mi->second[4])   // reverse was set
            bgeom.orient_reverse();
         face_bond_direct(base_out, bgeom, mi->first, mi->second[0], true);
      }
   }

   if(sym.is_set())
      sym_repeat(brick_out, brick_out, sym);

   if(out_type==out_default) {
      if(has_merged_brick)
         geom_out.append(base_out);
      geom_out.append(brick_out);
   }
   else if(out_type==out_brick)
      geom_out.append(brick_out);
   else if(out_type==out_base_brick) {
      if(has_merged_brick)
         geom_out.append(base_out);
      else
         geom_out.append(base);
      geom_out.append(brick_out);
   }
   else if(out_type==out_brick_base) {
      geom_out.append(brick_out);
      geom_out.append(base);
   }
   return true;
}


class align_opts: public prog_opts {
   public:
      vector<pair<char, string> > brick_args;
      string sym_str;
      int out_type;
      string ifile;
      string ofile;

      align_opts(): prog_opts("off_align"), out_type(0) {}

      void process_command_line(int argc, char **argv);
      void usage();
};


void align_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Read a base file and brick file in OFF format, and use vertices or faces\n"
"to position the brick with respect to the base. The brick may be repeated\n"
"symmetrically and/or merged with the base. If input_file is not given the\n"
"program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"  -v <arg>  align by vertices, arg is a comma separated list of a brick\n"
"            geometry (if empty use base) optionally followed by 'r' (reverse\n"
"            brick orientation) followed by two, four or six points given as\n"
"            the vertex index number in the OFF files (starting at 0). The base\n"
"            points are all given first and then the brick points.\n"
"            Formats:\n"
"               u1,v1 - point alignment\n"
"                  translation u1-v1 (so v1 of brick moves to u1 of base)\n"
"               u1,u2,v1,v2 - line alignment\n"
"                  point alignment as above followed by rotation through\n"
"                  u1 perpendicular to u1u2 and v1v2 to align u1u2 and v1v2.\n"
"               u1,u2,u3,v1,v2,v3 - face alignment\n"
"                  line alignment as above followed by a rotation\n"
"                  around u1,u2 so v3 lies in plane of u1u2u3.\n"
"  -f <arg>  align by face index, arg is a comma separated list of a brick\n"
"            geometry (if empty use base) followed by up to three numbers\n"
"            separated by commas: base face index, brick face index\n"
"            (default: 0), polygon alignment selection number (default: 0)\n"
"  -F <arg>  align and combine polyhedra by face index, arg is a comma\n"
"            separated list of a brick geometry (if empty use base) followed\n"
"            by up to three numbers: base face index, brick face index\n"
"            (default: 0), polygon alignment selection number (default: 0).\n"
"            Brick is after base, bond vertices are merged, bond faces are\n"
"            removed\n"
"  -M <val>  merge parts, select and order parts in the output, val may be:\n"
"               default (0):    any combined parts (-F) followed by any brick\n"
"                               parts (-v, -f)\n"
"               brick (1):      brick parts only\n"
"               base_brick (2): base part, possibly combined (-F). followed\n"
"                               any brick parts (-v, -f)\n"
"               brick_base (3): brick parts (-v, -f), followed by base part\n"
"  -y <sub>  repeat bricks according to symmetry of base. sub is symmetry\n"
"            subgroup (Schoenflies notation) or 'full' optionally followed\n"
"            by a ',' and conjugation type (integer)\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}


void align_opts::process_command_line(int argc, char **argv)
{
   char errmsg[MSG_SZ];
   opterr = 0;
   char c;
   string arg_id;
   mat3d trans_m2;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hv:f:F:M:y:o:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         case 'v':
         case 'F':
         case 'f':
            brick_args.push_back(pair<char, string>(c, optarg));
            break;

         case 'M':
            arg_id = get_arg_id(optarg,"default|brick|all_base_brick|brick_base",
                  argmatch_add_id_maps, errmsg);
            if(arg_id=="")
               error(errmsg);
            out_type = atoi(arg_id.c_str());
            break;
            
         case 'y':
            sym_str = optarg;
            break;
         
         case 'o':
            ofile = optarg;
            break;

         default:
            error("unknown command line error");
      }
   }

   if(brick_args.size()==0)
      error("no brick alignments given, must use option -v, -f or -F");
      
   if(argc-optind > 1)
      error("too many arguments");
   
   if(argc-optind == 1)
      ifile=argv[optind];

   // Cannot read more than one geometry from stdin
   bool base_from_stdin = (ifile=="" || ifile=="-");
   bool stdin_in_use = base_from_stdin;

   for(unsigned int i=0; i<brick_args.size(); i++) {
      const string &brick = brick_args[i].second;
      if(brick[0]=='-' && (brick.size()==1 || brick[1]==',')) {
         if(stdin_in_use) {
            if(base_from_stdin)
               error("cannot read both brick and base from standard input",
                     brick_args[i].first);
            else
               error("cannot read more than one brick from standard input",
                     brick_args[i].first);
         }
         else
            stdin_in_use = true;
      }
   }

}

int main(int argc, char *argv[])
{
   align_opts opts;
   opts.process_command_line(argc, argv);
   
   char errmsg[MSG_SZ];
   col_geom_v geom;
   geom_read_or_error(geom, opts.ifile, opts);

   bond_base base(geom);
   if(opts.sym_str!="" && !base.set_sym(opts.sym_str.c_str(), errmsg))
      opts.error(errmsg, 'y');
   
   vector<pair<char, string> >::iterator argi;
   for(argi=opts.brick_args.begin(); argi!=opts.brick_args.end(); ++argi) {
      if(!base.add_brick(argi->first, argi->second, errmsg))
         opts.error(errmsg, argi->first);
      else if(*errmsg)
         opts.warning(errmsg, argi->first);
   }

   col_geom_v geom_out;
   if(!base.bond_all(geom_out, opts.out_type, errmsg))
      opts.error(errmsg, 'M');
   geom_write_or_error(geom_out, opts.ofile, opts);

   return 0;
}



