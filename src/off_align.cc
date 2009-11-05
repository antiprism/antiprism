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
#include <unistd.h>
#include <math.h>

#include <string>
#include <vector>

#include "../base/antiprism.h"

using std::string;
using std::vector;
using std::map;


class align_opts: public prog_opts {
   public:
      vector<int> verts;
      vector<int> f_bond;
      bool rev_orient;
      int merge;
      bool merge_polys;
      bool scale;
      int align_cnt;
      string ifile;
      string brick_file;
      string ofile;

      align_opts(): prog_opts("off_align"), rev_orient(false), merge(-1),
                    merge_polys(false), scale(true), align_cnt(0) {}

      void process_command_line(int argc, char **argv);
      void usage();
};


void align_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Read a file in OFF format (the base) and use its vertices to position\n"
"copies of itself or another file in OFF format (the brick). If input_file\n"
"is not given the program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"  -p <pts>  two, four or six points given as the vertex index number\n"
"            in the OFF files (starting at 0), and separated by commas.\n"
"            The base points are all given first and then the brick points.\n"
"            Formats:\n"
"               u1,v1 - point alignment\n"
"                  translation u1-v1 (so v1 of brick moves to u1 of base)\n"
"               u1,u2,v1,v2 - line alignment\n"
"                  point alignment as above followed by rotation through\n"
"                  u1 perpendicular to u1u2 and v1v2 to align u1u2 and v1v2.\n"
"               u1,u2,u3,v1,v2,v3 - face alignment\n"
"                  line alignment as above followed by a rotation\n"
"                  around u1,u2 so v3 lies in plane of u1u2u3.\n"
"  -f <idxs> align by face index, up to three numbers separated by commas:\n"
"            base face index, brick face index (default: 0), polygon\n"
"            vertex offset (default: 0) to rotate the face\n"
"  -F <idxs> align and combine polyhedra by face index, up to three numbers\n"
"            separated by commas: base face index, brick face index\n"
"            (default: 0), polygon vertex offset (default: 0) to rotate the\n"
"            face. Brick after base, merge bond vertices, remove bond faces\n"
"  -b <file> brick file (default: base file)\n"
"  -M <val>  merge files. 0 (default) don't merge, 1 brick after\n"
"            base, 2 brick before base\n"
"  -R        reverse orientation of brick, make it bond on other side of face\n"
"  -x        don't scale the brick if the first base edge and brick edge are\n"
"            different lengths\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}


void align_opts::process_command_line(int argc, char **argv)
{
   char errmsg[MSG_SZ];
   opterr = 0;
   char c;
   int n;
   mat3d trans_m2;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hp:f:F:M:Rxb:o:")) != -1 ) {
      if(common_opts(c))
         continue;

      switch(c) {
         case 'p':
            align_cnt++;
            if(!read_int_list(optarg, verts, errmsg, true))
               error(errmsg, c);
            n = verts.size();
            if( n!=2 && n!=4 && n!=6) {
               snprintf(errmsg, MSG_SZ, "must give 2, 4 or 6 vertices (%d were given)",n);
               error(errmsg, c);
            }
            if(n>2) {
               for(int i=0; i<n/2; i++) {
                  if(verts[i]==verts[(i+1)%(n/2)])
                     error("repeated vertex index in base", c);
                  if(verts[n/2 + i]==verts[n/2 + (i+1)%(n/2)])
                     error("repeated vertex index in brick", c);
               }
            }
            break;

         case 'F':
            merge_polys = true;
            // fall through
         case 'f':
            align_cnt++;
            if(!read_int_list(optarg, f_bond, errmsg, true))
               error(errmsg, c);
            if(f_bond.size()>3) {
               snprintf(errmsg, MSG_SZ, "up to three arguments can be given (%lu were given)", (unsigned long)f_bond.size());
               error(errmsg, c);
            }
            f_bond.resize(3, 0);
            break;

         case 'M':
            if(!read_int(optarg, &merge, errmsg))
               error(errmsg, c);
            if(merge <0 || merge>3) {
               snprintf(errmsg, MSG_SZ, "merge value is %d, must be 0, 1, 2 or 3", merge);
               error(errmsg, c);
            }
            break;
            
         case 'R':
            rev_orient = true;
            break;
         
         case 'x':
            scale = false;
            break;
         
         case 'b':
            brick_file = optarg;
            break;
         
         case 'o':
            ofile = optarg;
            break;

         default:
            error("unknown command line error");
      }
   }

   if(align_cnt==0)
      error("no alignment given, must use option -p, -f or -F");
   if(align_cnt>1)
      error("more than one alignment given, must use only one -p, -f or -F");

   if(merge_polys && (merge==0 || merge==2))
      error("merge type 0 and 2 are not compatible with option -F", 'M');
   if(merge<0) // set default merge
      merge = 0;
      
   if(argc-optind > 1)
      error("too many arguments");
   
   if(argc-optind == 1)
      ifile=argv[optind];

}


int main(int argc, char *argv[])
{
   align_opts opts;
   opts.process_command_line(argc, argv);
   
   char errmsg[MSG_SZ];
   col_geom_v geom;
   if(!geom.read(opts.ifile, errmsg))
      opts.error(errmsg);

   col_geom_v brick_geom;
   if(opts.brick_file == "")
      brick_geom = geom;
   else {
      if(!brick_geom.read(opts.brick_file, errmsg))
         opts.error(errmsg);
   }
  
   if(opts.rev_orient)
      brick_geom.orient_reverse();

   vector<vec3d> &verts = *geom.get_verts();
   vector<vec3d> &brick_verts = *brick_geom.get_verts();

   vector<vector<vec3d> > pts(2);
   if(opts.verts.size()) {
      for(unsigned int i=0; i<opts.verts.size(); i++) {
         vector<vec3d> &vs = (i<opts.verts.size()/2) ? verts : brick_verts;
         if(opts.verts[i]<0 || opts.verts[i]>=(int)vs.size()) {
            snprintf(errmsg, MSG_SZ, "vertex %d (position %d) is out of bounds", opts.verts[i], i+1);
            opts.error(errmsg, 'p');
         }
         pts[!(i<opts.verts.size()/2)].push_back(vs[opts.verts[i]]);
         //fprintf(stderr, "pts[%d][%d] = ", (i<opts.verts.size()/2), pts[(i<opts.verts.size()/2)].size()-1);
         //pts[(i<opts.verts.size()/2)].back().dump();
      }
      if(opts.scale && pts[0].size()>1) {
         mat3d scl = mat3d::scale(
               (pts[0][0]-pts[0][1]).mag()/(pts[1][0]-pts[1][1]).mag()); 
         brick_geom.transform(scl);
         for(unsigned int i=0; i<pts[i].size(); i++)
            pts[1][i] = scl * pts[1][i];
      }
      brick_geom.transform(mat3d::alignment(pts[1],pts[0]));
   }
   if(opts.f_bond.size()) {
      vector<vector<int> > &fs = *geom.get_faces();
      vector<vector<int> > &bfs = *brick_geom.get_faces();
      char opt_c = opts.merge_polys ? 'F' : 'f';
      int f0 = opts.f_bond[0];
      int f1 = opts.f_bond[1];
      int f1_sz = bfs[f1].size(); 
      int offset = ((opts.f_bond[2]%f1_sz)+f1_sz)%f1_sz;
      if(f0<0 || f0>=(int)geom.get_faces()->size())
         opts.error("base face is out of bounds", opt_c);
      if(f1<0 || f1>=(int)brick_geom.get_faces()->size())
         opts.error("brick face is out of bounds", opt_c);
      if(opts.merge_polys &&
            (*geom.get_faces())[f0].size() != (*brick_geom.get_faces())[f1].size() )
         opts.warning("faces to bond are different sorts of polygon", 'F');
      if(opts.scale) {
         vec3d e0 = verts[fs[f0][1]] - verts[fs[f0][0]];
         vec3d e1 = brick_verts[bfs[f1][(offset+1)%f1_sz]] -
                                           brick_verts[bfs[f1][offset]];
         brick_geom.transform(mat3d::scale(e0.mag()/e1.mag()) );
      }
      face_bond(geom, brick_geom, f0, f1, offset, opts.merge_polys);
   }
      
   col_geom_v *out_geom = &geom;
   if(!opts.merge_polys) {
      if(opts.merge==0)            // Brick only
         out_geom = &brick_geom;
      else if(opts.merge==1) {     // Brick after base
         geom.append(brick_geom);
         out_geom = &geom;
      }
      else if(opts.merge==2) {     // Brick before base
         brick_geom.append(geom);
         out_geom = &brick_geom;
      }
   }

   /*col_geom_v *ogeom = &brick_geom;
   if(opts.merge==1) {         // Brick after base
      if(!opts.merge_polys)    // if -F then previously merged
         geom.append(brick_geom);
      ogeom = &geom;
   }
   else if(opts.merge==2)      // Brick before base
      brick_geom.append(geom);
   */

   if(!out_geom->write(opts.ofile, errmsg))
      opts.error(errmsg);

   return 0;
}


