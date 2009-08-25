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
   Name: poly_kscope.cc
   Description: linear transformations for OFF files
   Project: Antiprism - http://www.antiprism.com
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <unistd.h>
#include <math.h>
#include <limits.h>

#include <string>
#include <vector>
#include <algorithm>

#include "../base/antiprism.h"

using std::string;
using std::vector;
using std::map;
using std::swap;


class ksc_opts: public prog_opts
{
   public:
      mat3d trans_m;
      sch_sym sym;
      char col_elems;
      bool test;
      
      string sfile;
      string ifile;
      string ofile;

      ksc_opts(): prog_opts("poly_kscope"), col_elems('\0'), test(false) { }

      void process_command_line(int argc, char **argv);
      void usage();
};


void ksc_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"A polyhedral kaleidoscope. Read a file in OFF format and repeat it\n"
"in a symmetric arrangement like a kaleidoscope. If input_file is\n"
"not given the program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"  -s <symm> symmetry type. When a type contains 'n' this must be\n"
"            replaced by an integer (for S this integer must be even).\n"
"            Principal rotational axes are vertical. symm can be:\n"
"               Cs  - mirror\n"
"               Ci  - inversion\n"
"               Cn  - cyclic rotational\n"
"               Cnv - cyclic rotational with vertical mirror\n"
"               Cnh - cyclic rotational with horizontal mirror\n"
"               Dn  - dihedral rotational\n"
"               Dnv - dihedral rotational with vertical mirror\n"
"               Dnh - dihedral rotational with horizontal mirror\n"
"               Sn  - cyclic rotational (n/2-fold) with inversion\n"
"               T   - tetrahedral rotational\n"
"               Td  - tetrahedral rotational with mirror\n"
"               Th  - tetrahedral rotational with inversion\n"
"               O   - octahedral rotational\n"
"               Oh  - octahedral rotational with mirror\n"
"               I   - icosahedral rotational\n"
"               Ih  - icosahedral rotational with mirror\n"
"\n"
"  -S <file> use the symmetries of a polyhedron read from file\n"
"  -c <elms> color elements with a different index number for each part. The\n"
"            element string can include v, e and f to color, respectively,\n"
"            vertices, edges and faces\n"
"  -t        test the symmetry type by making a pattern of arrows\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}

void ksc_opts::process_command_line(int argc, char **argv)
{
   char errmsg[MSG_SZ];
   extern char *optarg;
   extern int optind, opterr;
   opterr = 0;
   char c;
   vector<double> nums;
   mat3d trans_m2;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hS:s:c:o:t")) != -1 ) {
      if(common_opts(c))
         continue;

      switch(c) {
         case 's':
            if(!sym.init(optarg, mat3d(), errmsg))
               error(errmsg, c);
            break;

         case 'S':
            sfile = optarg;
            break;

         case 'c':
            if(strspn(optarg, "vef") != strlen(optarg)) {
               snprintf(errmsg, MSG_SZ, "elements to color are %s must be v, e, or f\n", optarg);
               error(errmsg, c);
            }
            col_elems = (strchr(optarg, 'v')!=0)*ELEM_VERTS +
                        (strchr(optarg, 'e')!=0)*ELEM_EDGES +
                        (strchr(optarg, 'f')!=0)*ELEM_FACES;
            break;

         case 't':
            test = true;
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
   
   if(argc-optind == 1 && test)
      warning("input file ignored with -t");
   
   if(argc-optind == 1)
      ifile=argv[optind];

}

vector<int> args2face(int v1, ...)
{
   vector<int> face;
   va_list ap;
   va_start(ap, v1); 
   for(int i=v1; i!=-1; i=va_arg(ap, int))
      face.push_back(i);
   va_end(ap);
   return face;
}


void get_arrow(col_geom_v &geom, const sch_sym &sym)
{
   geom.clear_all();
   geom.add_vert(vec3d(0.0, 0.0, 0.0));
   geom.add_vert(vec3d(3.0, 0.0 ,0.0));
   geom.add_vert(vec3d(2.0, 1.0, 0.0));
   geom.add_vert(vec3d(2.0, 0.5, 0.0));
   geom.add_vert(vec3d(0.0, 0.5, 0.0));
   geom.add_vert(vec3d(0.0, 0.0, 0.5));
   geom.add_vert(vec3d(3.0, 0.0, 0.5));
   geom.add_vert(vec3d(2.0, 0.5, 0.25));
   geom.add_vert(vec3d(0.0, 0.5, 0.25));
   geom.add_col_face(args2face(1, 2, 3, 4, 0, -1), col_val(0.9, 0.9, 0.7));
   geom.add_col_face(args2face(6, 2, 7, 8, 5, -1), col_val(0.0, 0.2, 0.6));
   geom.add_col_face(args2face(0, 1, 6, 5, -1), col_val(0.9, 0.9, 0.7));
   geom.add_col_face(args2face(1, 2, 6, -1), col_val(0.9, 0.9, 0.7));
   geom.add_col_face(args2face(2, 3, 7, -1), col_val(0.9, 0.9, 0.7));
   geom.add_col_face(args2face(3, 4, 8, 7, -1), col_val(0.9, 0.9, 0.7));
   geom.add_col_face(args2face(4, 0, 5, 8, -1), col_val(0.9, 0.9, 0.7));

   mat3d trans = mat3d::transl(vec3d(0.1, 0.2, 0.4));
   int fold = sym.get_nfold();
   if(fold==1)
      fold = INT_MAX;

   switch(sym.get_sym_type()) {
      case sch_sym::C1:
      case sch_sym::Cs:
      case sch_sym::Ci:
         break;

      case sch_sym::C:
      case sch_sym::Ch:
      case sch_sym::S:
         trans = mat3d::transl(vec3d(1.6/sin(M_PI/fold), 0, 0)) *
                 mat3d::rot(vec3d::z, -M_PI/2) *
                 mat3d::transl(vec3d(-1.6,-0.2,0.0)) * trans;
         break;

      case sch_sym::Cv:
      case sch_sym::D:
      case sch_sym::Dh:
         trans = mat3d::transl(vec3d(3.2/tan(M_PI/fold), 0, 0)) *
                 mat3d::rot(vec3d::z, -M_PI/2) * trans;
         break;

      case sch_sym::Dv:
         trans = mat3d::rot(vec3d::z, 0.5*M_PI/fold) *
                 mat3d::transl(vec3d(3.2/tan(M_PI/fold), 0, 0)) *
                 mat3d::rot(vec3d::z, -M_PI/2) * trans;
         break;

      case sch_sym::T:
      case sch_sym::Td:
      case sch_sym::Th:
         trans = mat3d::alignment(vec3d::z, vec3d::x,
                                  vec3d(1,1,1), vec3d(0,-1,1)) *
                 mat3d::transl(vec3d(0, 0.5, 3)) *
                 mat3d::rot(vec3d::z, M_PI/2) * trans;
         break;

      case sch_sym::O:
      case sch_sym::Oh:
         //trans = mat3d::transl(vec3d(0, 3, 4)) * trans;
         trans = mat3d::transl(vec3d(1, 0, 4)) * trans;
         break;

      case sch_sym::I:
      case sch_sym::Ih:
         trans = mat3d::rot(vec3d::z, vec3d(0,1,1.618)) *
                 mat3d::transl(vec3d(0, 1, 5.5)) *
                 mat3d::rot(vec3d::z, M_PI/2) * trans;
         break;

   }

   geom.transform(trans);

}

int main(int argc, char *argv[])
{
   ksc_opts opts;
   opts.process_command_line(argc, argv);

   char errmsg[MSG_SZ];
   sch_sym final_sym;
   if(opts.sfile!="") {
      col_geom_v sgeom;
      if(!sgeom.read(opts.sfile, errmsg))
         opts.error(errmsg, 'S');
      if(*errmsg)
         opts.warning(errmsg, 'S');
      final_sym.init(sgeom);
   }
   else
      final_sym = opts.sym;

   col_geom_v geom;
   if(!opts.test) {
      if(!geom.read(opts.ifile, errmsg))
         opts.error(errmsg);
      if(*errmsg)
         opts.warning(errmsg);
   }
   else
      get_arrow(geom, final_sym);

   vector<vector<set<int> > > equivs;
   sch_sym part_sym(geom, &equivs);

   /*
      for(int i=0; i<3; i++) {
      const char *elems[] = { "verts", "edges", "faces" };
      fprintf(stderr, "\n%s\n", elems[i]);
      for(unsigned int j=0; j<equivs[i].size(); j++) {
      set<int>::iterator si;
      for(si=equivs[i][j].begin(); si!=equivs[i][j].end(); ++si)
      fprintf(stderr, "%d  ", *si);
      fprintf(stderr, "\n");
      }
      }
      */

   t_set min_ts;
   min_ts.min_set(final_sym.get_trans(), part_sym.get_trans());

   col_geom_v comp_geom;
   sym_repeat(comp_geom, geom, min_ts, opts.col_elems);

   if(!comp_geom.write(opts.ofile, errmsg))
      opts.error(errmsg);

   return 0;
}


