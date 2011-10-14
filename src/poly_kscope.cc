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
#include <set>
#include <algorithm>

#include "../base/antiprism.h"

using std::string;
using std::vector;
using std::map;
using std::set;
using std::make_pair;
using std::swap;


class ksc_opts: public prog_opts
{
   public:
      mat3d trans_m;
      sch_sym sym;
      sch_sym sub_sym;
      int sub_sym_conj;
      char col_elems;
      bool consider_part_sym;

      bool compound_print_list;
      int compound_number;
      string compound_realignment;

      string sfile;
      string ifile;
      string ofile;

      ksc_opts(): prog_opts("poly_kscope"),
                  sub_sym_conj(0), col_elems('\0'),
                  consider_part_sym(true),
                  compound_print_list(false), compound_number(-1)
                  { }

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
"  -s <sym>  symmetry type for kaleidoscope, up to three comma separated\n"
"            parts: main symmetry (Schoenflies notation) or file name,\n"
"            subgroup (Schoenflies notation), and conjugation type (integer)\n"
"  -y <arg>  make a compound by aligning the component to match a subsymmetry\n"
"            with the kaleidoscope. Argument is either 'list' (to print the\n"
"            list of the compounds to standard output) or a number from the\n"
"            list optionally followed by a comma and realignment (colon\n"
"            separated list of an integer then decimal numbers)\n"
"  -c <elms> color elements with a different index number for each part. The\n"
"            element string can include v, e and f to color, respectively,\n"
"            vertices, edges and faces\n"
"  -I        ignore shared symmetries, full kaleidoscopic repetition of\n"
"            component\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}

void ksc_opts::process_command_line(int argc, char **argv)
{
   char errmsg[MSG_SZ];
   opterr = 0;
   char c;
   col_geom_v sgeom;
   vector<char *> parts;
   vector<double> nums;
   mat3d trans_m2;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hs:c:Iy:o:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         case 's':
            split_line(optarg, parts, ",");
            if(parts.size()==0 || parts.size()>3)
               error("argument should have 1 to 3 comma separated parts", c);
            
            if(sgeom.read(parts[0], errmsg))
               sym.init(sgeom);
            else if(!sym.init(parts[0], mat3d(), errmsg))
               error(msg_str("invalid filename or symmetry type name: %s",
                        errmsg), c);

            if(parts.size()>1) {
               if(!sub_sym.init(parts[1], mat3d(), errmsg))
                  error(msg_str("sub-symmetry type: %s", errmsg), c);

               if(parts.size()>2) {
                  if(!read_int(parts[2], &sub_sym_conj, errmsg))
                     error(msg_str("sub-symmetry conjugation number: %s",
                              errmsg), c);
               }

               sym = sym.get_sub_sym(sub_sym, sub_sym_conj, errmsg);
               if(sym.get_sym_type() == sch_sym::unknown)
                     error(msg_str("sub-symmetry: %s", errmsg), c);
            }

            break;

         case 'c':
            if(strspn(optarg, "vef") != strlen(optarg))
               error(msg_str("elements to color are '%s' must be from v, e, f",
                        optarg), c);
            col_elems = (strchr(optarg, 'v')!=0)*ELEM_VERTS +
                        (strchr(optarg, 'e')!=0)*ELEM_EDGES +
                        (strchr(optarg, 'f')!=0)*ELEM_FACES;
            break;

         case 'y':
            if(strncmp(optarg, "list", strlen(optarg))==0)
               compound_print_list = true;
            else {
               split_line(optarg, parts, ",");
               if(parts.size()==0 || parts.size()>2)
                  error("argument should have 1 or 2 comma separated parts", c);
               if(!read_int(parts[0], &compound_number, errmsg))
                  error(errmsg, c);
               if(parts.size()>1)
                  compound_realignment = parts[1];
            }
            break;

         case 'I':
            consider_part_sym = false;
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

struct compound_list_item
{
   string sub;
   int type_sub;
   int type_comp;
   compound_list_item(const string &su, int typ_sub=0, int typ_comp=0):
      sub(su), type_sub(typ_sub), type_comp(typ_comp) {}
   //bool operator < (const compound_list_item &item2)
   //   { return sub.get_symbol() < item2.sub.get_symbol(); }
};

void compound_get_list(vector<compound_list_item> &compound_list,
      const sch_sym &part_sym, const sch_sym &comp_sym)
{
   compound_list.clear();
   const set<sch_sym> &subs = part_sym.get_sub_syms();
   set<sch_sym>::iterator si;
   map<string, pair<int, int> > sub_cnts;
   for(si=subs.begin(); si!=subs.end(); ++si)
      sub_cnts[si->get_symbol()].first++;
   const set<sch_sym> &subs2 = comp_sym.get_sub_syms();
   for(si=subs2.begin(); si!=subs2.end(); ++si)
      sub_cnts[si->get_symbol()].second++;

   map<string, pair<int, int> >::iterator item_i;
   for(item_i=sub_cnts.begin(); item_i!=sub_cnts.end(); ++item_i)
      for(int i=0; i<item_i->second.first; i++)
         for(int j=0; j<item_i->second.second; j++)
            compound_list.push_back(compound_list_item(item_i->first, i, j));
 
}

void compound_print_list(const sch_sym &part_sym, const sch_sym &comp_sym)
{
   FILE *ofile = stdout;
   vector<compound_list_item> compound_list;
   compound_get_list(compound_list, part_sym, comp_sym);
   for(int i=0; i<(int)compound_list.size(); ++i) {
      sch_sym sub(compound_list[i].sub);
      fprintf(ofile, "%3d: %5s", i, compound_list[i].sub.c_str());
      fprintf(ofile, " (%2d, %2d) ", 
            compound_list[i].type_sub, compound_list[i].type_comp);
      fprintf(ofile, " + %2d fixed",
            sub.get_autos().get_fixed().size());
      int free_rots = sub.get_autos().num_free_rots();
      if(free_rots==1)
         fprintf(ofile, " x axial rotation   ");
      else if(free_rots==3)
         fprintf(ofile, " x full  rotation   ");
      int free_transls = sub.get_autos().num_free_transls();
      if(free_transls==1)
         fprintf(ofile, " x axial translation");
      else if(free_transls==2)
         fprintf(ofile, " x plane translation");
      else if(free_transls==3)
         fprintf(ofile, " x space translation");
      fprintf(ofile, "\n");
   }
}

bool compound_get_component_trans(mat3d &trans,
      sch_sym part_sym, const sch_sym &comp_sym,
      int compound_number, const string &realignment, char *errmsg)
{
   vector<compound_list_item> compound_list;
   compound_get_list(compound_list, part_sym, comp_sym);
   if(compound_number<0 || compound_number>=(int)compound_list.size()) {
      if(errmsg)
         sprintf(errmsg, "compound number '%d' is not in range 0 to %d\n",
               compound_number, compound_list.size()-1);
      return false;
   }

   compound_list_item compound = compound_list[compound_number];
   char errmsg2[MSG_SZ];
   sch_sym part_sub_sym = part_sym.get_sub_sym(
         sch_sym(compound.sub), compound.type_sub);
   sch_sym comp_sub_sym = comp_sym.get_sub_sym(
         sch_sym(compound.sub), compound.type_comp);

   if(!part_sub_sym.get_autos().set_realignment(realignment.c_str(), errmsg2)) {
      if(errmsg)
         snprintf(errmsg, MSG_SZ, "sub-symmetry realignment: %s", errmsg2);
      return false;
   }

   trans = mat3d::inverse(comp_sub_sym.get_to_std()) *
           part_sub_sym.get_autos().get_realignment() *
           part_sub_sym.get_to_std(); 

   return true;
}




int main(int argc, char *argv[])
{
   ksc_opts opts;
   opts.process_command_line(argc, argv);

   char errmsg[MSG_SZ];

   col_geom_v geom;
   if(!geom.read(opts.ifile, errmsg))
      opts.error(errmsg);
   if(*errmsg)
      opts.warning(errmsg);



   vector<vector<set<int> > > equivs;
   sch_sym part_sym;
   if(opts.consider_part_sym)
      part_sym.init(geom, &equivs);

   if(opts.compound_print_list) {
      compound_print_list(part_sym, opts.sym);
      return 0;
   }
   else if(opts.compound_number>=0) {
      mat3d trans;
      if(!compound_get_component_trans(trans, part_sym, opts.sym,
            opts.compound_number, opts.compound_realignment, errmsg))
         opts.error(errmsg, 'y');
      geom.transform(trans);
      part_sym.init(geom);
   }


   t_set min_ts;
   min_ts.min_set(opts.sym.get_trans(), part_sym.get_trans());

   col_geom_v comp_geom;
   sym_repeat(comp_geom, geom, min_ts, opts.col_elems);

   if(!comp_geom.write(opts.ofile, errmsg))
      opts.error(errmsg);

   return 0;
}


