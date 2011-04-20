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
      string sfile;
      string ifile;
      string ofile;

      ksc_opts(): prog_opts("poly_kscope"), col_elems('\0') { }

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
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}

void ksc_opts::process_command_line(int argc, char **argv)
{
   char errmsg[MSG_SZ];
   opterr = 0;
   char c;
   vector<double> nums;
   mat3d trans_m2;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hS:s:c:o:")) != -1 ) {
      if(common_opts(c, optopt))
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
            if(strspn(optarg, "vef") != strlen(optarg))
               error(msg_str("elements to color are '%s' must be from v, e, f",
                        optarg), c);
            col_elems = (strchr(optarg, 'v')!=0)*ELEM_VERTS +
                        (strchr(optarg, 'e')!=0)*ELEM_EDGES +
                        (strchr(optarg, 'f')!=0)*ELEM_FACES;
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
   if(!geom.read(opts.ifile, errmsg))
      opts.error(errmsg);
   if(*errmsg)
      opts.warning(errmsg);

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


