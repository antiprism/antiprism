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
   Name: off_col.cc
   Description: program to colour OFF files
   Project: Antiprism - http://www.antiprism.com
*/

#include <stdio.h>
#include <stdlib.h>

#include <string.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include "../base/antiprism.h"

using std::string;
using std::vector;
using std::map;


///Colour processing using an HSVA range 
class color_proc_torange_hsv: public color_map_range
{
   protected:
      bool check_and_add_range(int idx, const char *rngs, char *errmsg);

   public:
      ///Initialise from a string
      bool init(const char *name, char *errmsg=0);

      ///Apply processing to the colour values of a geometry
      /** \param elem_cols element colours to process */
      void proc_map(map<int, col_val> &elem_cols);
      
      ///Get the colour value for an index number.
      /**\param col the color.
       * \return The processed colour. */
      col_val get_col(col_val col);
};

bool chk_range(vector<double> &v, char *errmsg)
{
   if(v.size()==1) {
      if(v[0]<-epsilon || v[0]>1+epsilon) {
         sprintf(errmsg, "value, %g, is not in rage 0.0 to 1.0", v[0]);
         return false;
      }
      v.push_back(v[0]);
   }
   else if(v.size()==2) {
      if(v[0]<-epsilon || v[0]>1+epsilon) {
         sprintf(errmsg, "first value, %g, is not in rage 0.0 to 1.0", v[0]);
         return false;
      }
      if(v[1]<-epsilon || v[1]>1+epsilon) {
         sprintf(errmsg, "second value, %g, is not in rage 0.0 to 1.0", v[1]);
         return false;
      }
   }
   else {
      sprintf(errmsg, "range has %lu values, must have 1 or 2",
            (unsigned long)v.size());
      return false;
   }
   if(v[0]>v[1]+epsilon)
      v[1]++;
   return true;
}   

bool color_proc_torange_hsv::check_and_add_range(int idx, const char *range,
      char *errmsg)
{
   char str[MSG_SZ];
   strncpy(str, range, MSG_SZ);
   str[MSG_SZ-1] = '\0';
   if(read_double_list(str, ranges[idx], errmsg, 2, ":") &&
         chk_range(ranges[idx], errmsg))
      return true;
   else
      return false;
}


bool color_proc_torange_hsv::init(const char *range_name, char *errmsg)
{
   set_func = &col_val::set_hsva;
   for(int i=0; i<4; i++)
      ranges[i].clear();
   ranges[0].push_back(0);
   ranges[0].push_back(1);
   ranges[1].push_back(0);
   ranges[1].push_back(1);
   ranges[2].push_back(0);
   ranges[2].push_back(1);
   ranges[3].push_back(0);
   ranges[3].push_back(1);

   char errmsg2[MSG_SZ];
   int name_len = strlen(range_name);
   char rngs[MSG_SZ];
   char *q = rngs;
   int cur_idx = -1;
   char cur_comp = 'X';
   for(const char *p=range_name; p-range_name<name_len+1; p++) {
      //fprintf(stderr, "*p = %c\n", *p);
      if(strchr("HhSsVvAa", *p) || cur_idx<0 || *p=='\0') {
         *q = '\0';
         if(cur_idx>=0 && !check_and_add_range(cur_idx, rngs, errmsg2)) {
            if(errmsg)
               sprintf(errmsg, "component '%c': %s", cur_comp, errmsg2);
            return false;
         }
         if(strchr("Hh", *p))
            cur_idx = 0;
         else if(strchr("Ss", *p))
            cur_idx = 1;
         else if(strchr("Vv", *p))
            cur_idx = 2;
         else if(strchr("Aa", *p))
            cur_idx = 3;
         else {
            if(errmsg)
               sprintf(errmsg, "invalid component letter '%c'", *p);
            return false;
         }
         cur_comp = *p;
         q = rngs;
      }
      else if(!(isdigit(*p) || *p == '.' || *p == ':')) {
         if(errmsg)
            sprintf(errmsg, "invalid component letter '%c'",
                  (cur_idx<0) ? *rngs : *p);
         return false;
      }
      else if(!isspace(*p)) {
         *q++ = *p;
      }
   }

   //for(int i=0; i<4; i++)
   //   for(unsigned int j=0; j<ranges[i].size(); j++)
   //      fprintf(stderr, "ranges[%d][%u] = %g\n", i, j, ranges[i][j]);

   return true;
}

inline double fract(vector<double> &rng, double frac)
{
   return fmod(rng[0] + (rng[1]-rng[0])*frac, 1+epsilon);
}


col_val color_proc_torange_hsv::get_col(col_val col)
{
   col_val c = col;
   if(col.is_val()) {
      vec4d hsva_orig = col.get_hsva();
      vec4d hsva( fract(ranges[0], hsva_orig[0]),
                  fract(ranges[1], hsva_orig[1]),
                  fract(ranges[2], hsva_orig[2]),
                  fract(ranges[3], hsva_orig[3]) );
      c.set_hsva(hsva);
   }
   return c;
}

void color_proc_torange_hsv::proc_map(map<int, col_val> &elem_cols)
{
   map<int, col_val>::iterator mi;
   for(mi=elem_cols.begin(); mi!=elem_cols.end(); ++mi)
      mi->second = get_col(mi->second);
}


void color_vals_to_idxs(col_geom_v &geom, char elems=ELEM_ALL,
      color_map_map *cmap=0)
{
   if(cmap)
      cmap->clear();
   
   map<col_val, vector<vector<int> > > val2idxs;
   int first_idx = 0;
   map<int, col_val> *elem_cols[3] = {
      (elems & ELEM_VERTS) ? &geom.raw_vert_cols() : 0, 
      (elems & ELEM_EDGES) ? &geom.raw_edge_cols() : 0, 
      (elems & ELEM_FACES) ? &geom.raw_face_cols() : 0
   };
   for(int i=0; i<3; i++) {
      if(elem_cols[i]) {
         map<int, col_val>::const_iterator mi;
         for(mi=elem_cols[i]->begin(); mi!=elem_cols[i]->end(); ++mi) {
            const col_val &col = mi->second;
            if(col.is_idx()) {
               if(col.get_idx()>first_idx)
                  first_idx = col.get_idx()+1;
            }
            else if(col.is_val()) {
               map<col_val, vector<vector<int> > >::iterator v2i_it;
               v2i_it = val2idxs.find(col);
               if(v2i_it==val2idxs.end()) {
                  pair<map<col_val,
                     vector<vector<int> > >::iterator, bool> ins =
                        val2idxs.insert(make_pair(col,vector<vector<int> >(3)));
                  v2i_it = ins.first;
                  v2i_it->second.resize(3);
               }
               v2i_it->second[i].push_back(mi->first);
            }
         }
      }
   }
   
   int idx_inc = 0;
   map<col_val, vector<vector<int> > >::const_iterator vmi;
   for(vmi=val2idxs.begin(); vmi!=val2idxs.end(); ++vmi) {
      int idx_no = first_idx + idx_inc++;
      for(int i=0; i<3; i++)
         if(elem_cols[i])
            for(unsigned int j=0; j<vmi->second[i].size(); j++)
               (*elem_cols[i])[vmi->second[i][j]] = col_val(idx_no);
      if(cmap)
         cmap->set_col(idx_no, vmi->first);
   }
}

enum {CV_UNSET=1, CV_INDEX=2, CV_VALUE=4, CV_INVISIBLE=8};

class o_col_opts: public prog_opts {
   public:
      char v_col_op;
      sch_sym v_sub_sym;
      vector<set<int> > v_equivs;
      col_val v_col;

      char e_col_op;
      sch_sym e_sub_sym;
      vector<set<int> > e_equivs;
      col_val e_col;
       
      char f_col_op;
      sch_sym f_sub_sym;
      vector<set<int> > f_equivs;
      col_val f_col;
     
      char edge_type;
      unsigned int selection;

      coloring clrngs[3];

      char range_elems;
      color_proc_torange_hsv col_procs[3];
      char v2i_elems;
      
      string lfile;
      string ifile;
      string ofile;

      o_col_opts(): prog_opts("off_color"),
                    v_col_op(0), e_col_op(0), f_col_op(0),
                    edge_type('x'), selection(0),
                    range_elems(ELEM_NONE), v2i_elems(ELEM_NONE)
         {}

      void process_command_line(int argc, char **argv);
      void usage();
};


void o_col_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Read a file in OFF format and add colours to it. If input_file is\n"
"not given the program reads from standard input.\n"
"\n"
"A colour value can be a single integer (a colour index), a value in\n"
"form 'R,G,B,A' (3 or 4 values 0.0-1.0, or 0-255) or hex 'xFFFFFF', a\n"
"colour name from the X11 colour map, 'invisible' or 'none' (which sets\n"
"without any colour information)\n"
"\n"
"Lowercase letters (except l) colour using index numbers and uppercase\n"
"letters colour using colour values\n"
"\n"
"Options\n"
"%s"
"  -f <col>  colour the faces according to:\n"
"               a colour value - apply to all faces\n"
"               u,U - unique colour\n"
"               p,P - minimal proper colouring\n"
"               s,S - symmetric colouring\n"
"               n,N - colour by number of sides\n"
"               a,A - colour by average internal angle (to nearest degree)\n"
"               k,K - sets of faces connected by face edges\n"
"               g,G - gradient on Y coordinate of normal\n"
"               c,C - gradient on Y coordinate of centroid\n"
"               L   - lighting effect by normal (see option -l)\n"
"               l   - lighting effect by centroid (see option -l)\n"
"               M   - use colour map to convert existing colour index numbers\n"
"                     into to values\n"
"  -E <type> colour by edge type, e - explicit edges, i - implicit edges\n"
"            I - implicit edges which are not explicit edges (default: if -e\n"
"            explicit and implicit, else explicit for mapping only)\n"
"  -e <col>  colour the edges according to:\n"
"               a colour value - apply to all edges\n"
"               u,U - unique colour\n"
"               p,P - minimal proper colouring\n"
"               s,S - symmetric colouring\n"
"               k,K - sets of edges connected by edges\n"
"               F   - colour with average adjoining face colour\n"
"               g,G - gradient on Y coordinate of edge direction\n"
"               c,C - gradient on Y coordinate of centroid\n"
"               L   - lighting effect (see option -l)\n"
"               M   - use colour map to convert existing colour index numbers\n"
"                     into to values\n"
"  -v <col>  colour the vertices according to:\n"
"               a colour value - apply to all vertices\n"
"               u,U - unique colour\n"
"               p,P - minimal proper colouring\n"
"               s,S - symmetric colouring\n"
"               n,N - colour by order of vertex\n"
"               F   - colour with average adjoining face colour\n"
"               E   - colour with average adjoining edge colour\n"
"               c,C - gradient on Y coordinate\n"
"               L   - lighting effect (see option -l)\n"
"               M   - use colour map to convert existing colour index numbers\n"
"                     into to values\n"
"  -l <file> lights for colouring type L in a file in OFF format, each\n"
"            vertex and its colour gives a light direction and colour\n"
"            (default: a set of six lights giving a rainbow colouring\n"
"  -m <maps> a comma separated list of colour maps used to transform colour\n"
"            indexes (default: rand), a part consisting of letters from\n"
"            v, e, f, selects the element types to apply the map list to\n"
"            (default 'vef').\n"
"  -r <rnge> Map HSVA values onto the specified HSVA ranges after other\n"
"            processing (but before -I), component letters are followed by\n"
"            one or two values separated by a colon e.g H0.5:0.8 followed by\n"
"            elements to map from v, e and f (H0:1,S0:1,V0:1,A0:1,vef)\n"
"  -I <elms> map color values to index numbers (after other procesing)\n"
"            elements to map are from v, e and f (default none)\n"
"  -w <wdth> width of sphere containing points (default: calculated)\n"
"  -U <typs> colour only elements with particular current colour types:\n"
"            u - unset, i - indexed, v - visible colour value, x - invisible.\n"
"            ~ before the letter will select the opposite.\n" 
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}

void o_col_opts::process_command_line(int argc, char **argv)
{
   char errmsg[MSG_SZ];
   char errmsg2[MSG_SZ];
   opterr = 0;
   vector<char *> parts;
   bool prev_char_was_not;
   char c;
   
   handle_long_opts(argc, argv);

   while( (c = getopt(argc, argv, ":hv:f:e:E:s:m:M:c:l:U:o:r:I:")) != -1 ) {
      if(common_opts(c, optopt))
         continue;

      switch(c) {
         case 'v':
            if(v_col.read(optarg, errmsg)) {
               v_col_op = 'o';
               break;
            }
            split_line(optarg, parts, ",");
            if(strlen(parts[0])==1 && strchr("uUpPsSnNFEcCLM", *parts[0]))
               v_col_op = *parts[0];
            else
               error("invalid colouring", c);

            if(parts.size()>2 || (strchr("sS", (char)v_col_op)&&parts.size()>3))
               error("too many comma separated parts", c);
            
            if(strchr("sS", v_col_op)) {
               v_sub_sym = sch_sym();
               if(parts.size()==2 && !v_sub_sym.init(parts[1], mat3d(),errmsg2))
                  error(msg_str("invalid subsymmetry: %s", errmsg2), c);
            }
            break;

         case 'f':
            if(f_col.read(optarg, errmsg)) {
               f_col_op = 'o';
               break;
            }
            split_line(optarg, parts, ",");
            if(strlen(parts[0])==1 &&strchr("uUpPsSnNaAkKgGcCLlM",*parts[0]))
               f_col_op = *parts[0];
            else
               error("invalid colouring", c);

            if(parts.size()>2 || (strchr("sS", (char)f_col_op)&&parts.size()>3))
               error("too many comma separated parts", c);
            
            if(strchr("sS", f_col_op)) {
               f_sub_sym = sch_sym();
               if(parts.size()==2 && !f_sub_sym.init(parts[1], mat3d(),errmsg2))
                  error(msg_str("invalid subsymmetry: %s", errmsg2), c);
            }
            break;

         case 'e':
            if(e_col.read(optarg, errmsg)) {
               e_col_op = 'o';
               break;
            }
            split_line(optarg, parts, ",");
            if(strlen(parts[0])==1 && strchr("uUpPsSkKFgGcCLM", *parts[0]))
               e_col_op = *parts[0];
            else
               error("invalid colouring", c);

            if(parts.size()>2 || (strchr("sS", (char)e_col_op)&&parts.size()>3))
               error("too many comma separated parts", c);
            
            if(strchr("sS", e_col_op)) {
               e_sub_sym = sch_sym();
               if(parts.size()==2 && !e_sub_sym.init(parts[1], mat3d(),errmsg2))
                  error(msg_str("invalid subsymmetry: %s", errmsg2), c);
            }
            break;

         case 'E':
            if(!strlen(optarg)==1 || !strchr("eiI", *optarg))
               error("edge type to color must be e, i or I");
            edge_type = *optarg;
            break;

         case 'l':
            lfile = optarg;
            break;

         case 'm':
            *errmsg = '\0';
            if(!read_colorings(clrngs, optarg, errmsg))
               error(errmsg, c);
            if(*errmsg)
               warning(errmsg, c);
            break;

         case 'r': {
            char r_elems;
            color_proc_torange_hsv col_proc;
            if(split_line(optarg, parts, ",")>2)
               error("too many comma separated parts", c);
            if(!col_proc.init(parts[0], errmsg))
               error(errmsg, c);
            if(parts.size()>1) {
               if(strspn(parts[1], "vef") != strlen(parts[1]))
                  error(msg_str("elements for colour ranges are '%s' must be "
                     "from v, e, and f", optarg), c);
               r_elems = (strchr(parts[1], 'v')!=0)*ELEM_VERTS +
                             (strchr(parts[1], 'e')!=0)*ELEM_EDGES +
                             (strchr(parts[1], 'f')!=0)*ELEM_FACES;
            }
            else
               r_elems = ELEM_VERTS + ELEM_EDGES + ELEM_FACES;
            
            for(int i=0; i<3; i++)
               if((r_elems & (1<<i)))
                  col_procs[i] = col_proc;

            range_elems |= r_elems;
            break;
         }

         case 'I':
            if(strspn(optarg, "vef") != strlen(optarg))
               error(msg_str("elements to map are '%s' must be "
                        "from v, e, and f", optarg), c);
            v2i_elems = (strchr(optarg, 'v')!=0)*ELEM_VERTS +
                        (strchr(optarg, 'e')!=0)*ELEM_EDGES +
                        (strchr(optarg, 'f')!=0)*ELEM_FACES;
            break;

          case 'U':
            selection = 0;
            prev_char_was_not = false;
            for(const char *p = optarg; *p; p++) {
               unsigned int new_select = 0;
               if(*p == '~') {
                  if(prev_char_was_not)
                     error("types cannot include repeated '~'", c);
                  prev_char_was_not = true;
               }
               else if(*p == 'u')
                  new_select = CV_UNSET;
               else if(*p == 'i')
                  new_select = CV_INDEX;
               else if(*p == 'v')
                  new_select = CV_VALUE;
               else if(*p == 'x')
                  new_select = CV_INVISIBLE;
               else
                  error(msg_str("invalid type character '%c'", *p), c);

               if(new_select)
                  selection |= (prev_char_was_not) ? ~new_select : new_select;
            }
            break;

         case 'o':
            ofile = optarg;
            break;

         default:
            error("unknown command line error");
      }
   }


   // default: explicit for mapping only othrwise explicit and implicit
   if(edge_type=='x')
      edge_type = (!e_col_op || e_col_op=='M') ? 'e' : 'a';
   
   if(argc-optind > 1)
      error("too many arguments");
   
   if(argc-optind == 1)
      ifile=argv[optind];

}


bool lights_read(const string &fname, col_geom_v *lights, char *errmsg)
{
   *errmsg = '\0';
   FILE *lfile = open_sup_file(fname.c_str(), "/col_lights/");
   if(lfile==0) {
      if(errmsg)
         snprintf(errmsg, MSG_SZ, "could not open colour lights file \'%s\'",
               fname.c_str());
      return 0;
   }
  
   if(!lights->read(lfile, errmsg))
      return false;
   bool no_color = false;
   for(unsigned int i=0; i<lights->get_verts()->size(); i++) {
      if(!lights->get_v_col(i).is_val()) {
         lights->set_v_col(i, col_val(0.5,0.5,0.5));
         no_color = true;
      }
   }
   if(no_color)
      strcpy(errmsg, "one or more missing colours, set to grey");
   
   return true;
}

inline unsigned int col_type(const col_val &col)
{ return CV_UNSET*col.is_def() + CV_INDEX*col.is_idx() + 
         CV_VALUE*(col.is_val()&&!col.is_inv()) + CV_INVISIBLE*col.is_inv(); }

void restore_orig_cols(col_geom_v &geom, col_geom_v &restore_geom,
      unsigned int selection, unsigned int orig_edges_sz)
{
   for(unsigned int i=0; i<restore_geom.verts().size(); i++) {
      col_val col = restore_geom.get_v_col(i);
      if(!(col_type(col)&selection))  // restore colours of unselected elements
         geom.set_v_col(i, col);
   }

   vector<int> del_edges;
   for(unsigned int i=0; i<geom.edges().size(); i++) {
      col_val col;
      if(i<restore_geom.edges().size())
         col = restore_geom.get_e_col(i);
      if(!(col_type(col)&selection))  // restore cols of unselected elements
         geom.set_e_col(i, col);
      // implicit edges with unset colour were not selected to be coloured
      if(i>=orig_edges_sz && !geom.get_e_col(i).is_set())
         del_edges.push_back(i);
   }
   geom.delete_edges(del_edges);

   for(unsigned int i=0; i<restore_geom.faces().size(); i++) {
      col_val col = restore_geom.get_f_col(i);
      if(!(col_type(col)&selection))  // restore colours of unselected elements
         geom.set_f_col(i, col);
   }
}
   

int main(int argc, char *argv[])
{
   o_col_opts opts;
   opts.process_command_line(argc, argv);

   char errmsg[MSG_SZ];
   col_geom_v geom;
   if(!geom.read(opts.ifile, errmsg))
   if(!geom)
      opts.error(errmsg);
   if(*errmsg)
      opts.warning(errmsg);

   // read lights
   col_geom_v lights;
   if(opts.lfile!="") {
      if(!lights_read(opts.lfile, &lights, errmsg))
         opts.error(errmsg, 'l');
      if(*errmsg)
         opts.warning(errmsg, 'l');
   }
   

   // store original edges
   vector<vector<int> > &edges = *geom.get_edges();
   unsigned int orig_edges_sz = edges.size();
   map<vector<int>, col_val> expl_edges;
   if(opts.e_col_op && strchr("pP", opts.e_col_op) &&
         !strchr("iI", opts.edge_type))
      opts.edge_type = 'i';
   if(opts.edge_type=='a')
      geom.add_missing_impl_edges();
   else if(strchr("iI", opts.edge_type)) {
      for(unsigned int i=0; i<edges.size(); i++)
         expl_edges[edges[i]] = geom.get_e_col(i);
      geom.clear_edges();
      geom.add_missing_impl_edges();
   }

   col_geom_v store_geom;
   if(opts.selection)
      store_geom = geom;


   // Get symmetry if necessary
   sch_sym sym;
   vector<vector<set<int> > > sym_equivs;
   if( (opts.f_col_op && strchr("sS", opts.f_col_op)) ||
       (opts.e_col_op && strchr("sS", opts.e_col_op)) ||
       (opts.v_col_op && strchr("sS", opts.v_col_op))   ) {
      sym.init(geom, &sym_equivs);
      opts.v_equivs = sym_equivs[0];
      opts.e_equivs = sym_equivs[1];
      opts.f_equivs = sym_equivs[2];

      if(opts.v_sub_sym.get_sym_type() != sch_sym::unknown) {
         sch_sym sub = sym.get_sub_sym(opts.v_sub_sym.get_sym_type(),
               opts.v_sub_sym.get_nfold());
         get_equiv_elems(geom, sub.get_trans(), &sym_equivs);
         opts.v_equivs = sym_equivs[0];
      }
      if(opts.e_sub_sym.get_sym_type() != sch_sym::unknown) {
         sch_sym sub = sym.get_sub_sym(opts.e_sub_sym.get_sym_type(),
               opts.e_sub_sym.get_nfold());
         get_equiv_elems(geom, sub.get_trans(), &sym_equivs);
         opts.e_equivs = sym_equivs[1];
      }  
      if(opts.f_sub_sym.get_sym_type() != sch_sym::unknown) {
         sch_sym sub = sym.get_sub_sym(opts.f_sub_sym.get_sym_type(),
               opts.f_sub_sym.get_nfold());
         get_equiv_elems(geom, sub.get_trans(), &sym_equivs);
         opts.f_equivs = sym_equivs[2];
      }
   }

   //sch_sym sub = sym.get_sub_sym(sch_sym::C, 5);
   //const t_set &ts = sub.get_trans();
   //t_set::const_iterator si;
   //for(si=ts.begin(); si!=ts.end(); si++)
   //   si->dump();

   //get_equiv_elems(geom, sub.get_trans(), &sym_equivs);


   //sym_repeat(geom, geom, ts);

  
   coloring &fc = opts.clrngs[2];
   fc.set_geom(&geom);
   if(opts.f_col_op) {
      char op = opts.f_col_op;
      color_map *cmap = 0;
      if(fc.get_cmaps().size()==0) {
         if(strchr("GgCc", op))
            cmap = init_color_map("range");
         else
            cmap = init_color_map("spread");
      }
      if(cmap)
         fc.add_cmap(cmap);
      if(op=='o')
         fc.f_one_col(opts.f_col);
      else if(strchr("uU", op))
         fc.f_unique(op=='U');
      else if(strchr("pP", op))
         fc.f_proper(op=='P');
      else if(strchr("sS", op))
         fc.f_sets(sym_equivs[2], op=='S');
      else if(strchr("nN", op))
         fc.f_sides(op=='N');
      else if(strchr("aA", op))
         fc.f_avg_angle(op=='A');
      else if(strchr("kK", op))
         fc.f_parts(op=='K');
      else if(strchr("Gg", op))
         fc.f_normal(op=='G');
      else if(strchr("Cc", op))
         fc.f_centroid(op=='C');
      else if(strchr("L", op))
         fc.f_lights(lights);
      else if(strchr("l", op))
         fc.f_lights2(lights);
      else if(strchr("M", op))
         fc.f_apply_cmap();
   }

  
   coloring &ec = opts.clrngs[1];
   ec.set_geom(&geom);
   if(opts.e_col_op) {
      char op = opts.e_col_op;
      color_map *cmap = 0;
      if(ec.get_cmaps().size()==0) {
         if(strchr("GgCc", op))
            cmap = init_color_map("range");
         else
            cmap = init_color_map("spread");
      }
      if(cmap)
         ec.add_cmap(cmap);
      if(op=='o')
         ec.e_one_col(opts.e_col);
      else if(strchr("uU", op))
         ec.e_unique(op=='U');
      else if(strchr("pP", op))
         ec.e_proper(op=='P');
      else if(strchr("sS", op))
         ec.e_sets(sym_equivs[1], op=='S');
      else if(strchr("kK", op))
         ec.e_parts(op=='K');
      else if(strchr("Gg", op))
         ec.e_direction(op=='G');
      else if(strchr("Cc", op))
         ec.e_mid_point(op=='C');
      else if(strchr("L", op))
         ec.e_lights(lights);
      else if(strchr("M", op))
         ec.e_apply_cmap();

   }
   
   coloring &vc = opts.clrngs[0];
   vc.set_geom(&geom);
   if(opts.v_col_op) {
      char op = opts.v_col_op;
      color_map *cmap = 0;
      if(vc.get_cmaps().size()==0) {
         if(strchr("Cc", op))
            cmap = init_color_map("range");
         else
            cmap = init_color_map("spread");
      }
      if(cmap)
         vc.add_cmap(cmap);
      if(op=='o')
         vc.v_one_col(opts.v_col);
      else if(strchr("uU", op))
         vc.v_unique(op=='U');
      else if(strchr("pP", op))
         vc.v_proper(op=='P');
      else if(strchr("sS", op))
         vc.v_sets(sym_equivs[0], op=='S');
      else if(strchr("nN", op))
         vc.v_order(op=='N');
      else if(strchr("cC", op))
         vc.v_position(op=='C');
      else if(strchr("L", op))
         vc.v_lights(lights);
      else if(strchr("M", op))
         vc.v_apply_cmap();
   }

   /*
   // convert index numbers to values after other processing
   if(col_map.size()) {
      if(strchr(opts.cmap_elems.c_str(), 'f'))
         fc.apply_cmap();
      if(strchr(opts.cmap_elems.c_str(), 'e'))
         ec.apply_cmap();
      if(strchr(opts.cmap_elems.c_str(), 'v'))
         vc.apply_cmap();
   }
   */

   // value to value mappings
   if(opts.range_elems & (ELEM_VERTS))
      opts.col_procs[0].proc_map(geom.raw_vert_cols());
   if(opts.range_elems & (ELEM_EDGES))
      opts.col_procs[1].proc_map(geom.raw_edge_cols());
   if(opts.range_elems & (ELEM_FACES))
      opts.col_procs[2].proc_map(geom.raw_face_cols());

   
   // Average colour values from adjoining elements after converting
   // index numbers
   if(opts.e_col_op=='F')
      ec.e_face_color();
   if(opts.v_col_op=='F')
      vc.v_face_color();
   else if(opts.v_col_op=='E')
      vc.v_edge_color();

   // Finally convert to index numbers
   color_vals_to_idxs(geom, opts.v2i_elems);


   if(opts.selection)
      restore_orig_cols(geom, store_geom, opts.selection, orig_edges_sz);
   

   // restore original edges
   map<vector<int>, col_val>::iterator ei;
   if(opts.edge_type=='I') {
      for(ei=expl_edges.begin(); ei!=expl_edges.end(); ++ei)
         geom.add_col_edge(ei->first, ei->second);
   }
   else if(opts.edge_type=='i') {
      for(ei=expl_edges.begin(); ei!=expl_edges.end(); ++ei)
         if(find(edges.begin(), edges.end(), ei->first)==edges.end())
            geom.add_col_edge(ei->first, ei->second);
   }

   if(!geom.write(opts.ofile, errmsg))
      opts.error(errmsg);

   return 0;
}


