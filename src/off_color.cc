/*
   Copyright (c) 2003-2016, Adrian Rossiter

   Antiprism - http://www.antiprism.com

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

#include <algorithm>
#include <ctype.h>
#include <map>
#include <math.h>
#include <set>
#include <string.h>
#include <string>
#include <time.h>
#include <vector>

#include "../base/antiprism.h"

using std::map;
using std::pair;
using std::set;
using std::string;
using std::vector;

using namespace anti;

void color_vals_to_idxs(Geometry &geom, char elems = ELEM_ALL,
                        ColorMapMap *cmap = nullptr)
{
  if (cmap)
    cmap->clear();

  map<Color, vector<vector<int>>> val2idxs;
  int first_idx = 0;
  map<int, Color> *elem_cols[3] = {
      (elems & ELEM_VERTS) ? &geom.colors(VERTS).get_properties() : nullptr,
      (elems & ELEM_EDGES) ? &geom.colors(EDGES).get_properties() : nullptr,
      (elems & ELEM_FACES) ? &geom.colors(FACES).get_properties() : nullptr};
  for (int i = 0; i < 3; i++) {
    if (elem_cols[i]) {
      map<int, Color>::const_iterator mi;
      for (mi = elem_cols[i]->begin(); mi != elem_cols[i]->end(); ++mi) {
        const Color &col = mi->second;
        if (col.is_index()) {
          if (col.get_index() > first_idx)
            first_idx = col.get_index() + 1;
        }
        else if (col.is_value()) {
          map<Color, vector<vector<int>>>::iterator v2i_it;
          v2i_it = val2idxs.find(col);
          if (v2i_it == val2idxs.end()) {
            pair<map<Color, vector<vector<int>>>::iterator, bool> ins =
                val2idxs.insert(make_pair(col, vector<vector<int>>(3)));
            v2i_it = ins.first;
            v2i_it->second.resize(3);
          }
          v2i_it->second[i].push_back(mi->first);
        }
      }
    }
  }

  int idx_inc = 0;
  map<Color, vector<vector<int>>>::const_iterator vmi;
  for (vmi = val2idxs.begin(); vmi != val2idxs.end(); ++vmi) {
    int idx_no = first_idx + idx_inc++;
    for (int i = 0; i < 3; i++)
      if (elem_cols[i])
        for (unsigned int j = 0; j < vmi->second[i].size(); j++)
          (*elem_cols[i])[vmi->second[i][j]] = Color(idx_no);
    if (cmap)
      cmap->set_col(idx_no, vmi->first);
  }
}

enum { CV_UNSET = 1, CV_INDEX = 2, CV_VALUE = 4, CV_INVISIBLE = 8 };

class o_col_opts : public ProgramOpts {
public:
  char v_col_op;
  string v_sub_sym;
  vector<set<int>> v_equivs;
  Color v_col;

  char e_col_op;
  string e_sub_sym;
  vector<set<int>> e_equivs;
  Color e_col;

  char f_col_op;
  string f_sub_sym;
  vector<set<int>> f_equivs;
  Color f_col;

  char edge_type;
  unsigned int selection;

  Coloring clrngs[3];

  char range_elems;
  ColorValuesToRangeHsva col_procs[3];
  char v2i_elems;

  string lfile;
  string ifile;
  string ofile;

  o_col_opts()
      : ProgramOpts("off_color"), v_col_op(0), e_col_op(0), f_col_op(0),
        edge_type('x'), selection(0), range_elems(ELEM_NONE),
        v2i_elems(ELEM_NONE)
  {
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

// clang-format off
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
"letters colour using colour values. Symmetric colourings are optionally\n"
"followed, separated by commas, by a subsymmetry (Schoenflies notation)\n"
"and a conjugation type (integer)\n"
"\n"
"Options\n"
"%s"
"  -f <col>  colour the faces according to:\n"
"               a colour value - apply to all faces\n"
"               u,U - unique colour\n"
"               p,P - minimal proper colouring\n"
"               s,S - symmetric colouring [,sub_group,conj_type] (see above)\n"
"               n,N - colour by number of sides\n"
"               a,A - colour by average internal angle (to nearest degree)\n"
"               k,K - sets of faces connected by face edges\n"
"               g,G - gradient on z-coordinate of normal\n"
"               c,C - gradient on z-coordinate of centroid\n"
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
"               s,S - symmetric colouring [,sub_group,conj_type] (see above)\n"
"               k,K - sets of edges connected by edges\n"
"               F   - colour with average adjoining face colour\n"
"               d,D - colour by edge direction\n"
"               g,G - gradient on z-coordinate of edge direction\n"
"               c,C - gradient on z-coordinate of centroid\n"
"               L   - lighting effect (see option -l)\n"
"               l   - lighting effect on edge directions (see option -l)\n"
"               M   - use colour map to convert existing colour index numbers\n"
"                     into to values\n"
"  -v <col>  colour the vertices according to:\n"
"               a colour value - apply to all vertices\n"
"               u,U - unique colour\n"
"               p,P - minimal proper colouring\n"
"               s,S - symmetric colouring [,sub_group,conj_type] (see above)\n"
"               n,N - colour by order of vertex\n"
"               a,A - colour by avg internal ang of vert-fig (to nearest deg)\n"
"               F   - colour with average adjoining face colour\n"
"               E   - colour with average adjoining edge colour\n"
"               c,C - gradient on z-coordinate\n"
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
"            a comma and elements to map from v, e, f (H0:1S0:1V0:1A0:1,vef)\n"
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
// clang-format on

void o_col_opts::process_command_line(int argc, char **argv)
{
  Status stat;
  opterr = 0;
  vector<char *> parts;
  string optarg_orig;
  bool prev_char_was_not;
  int c;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hv:f:e:E:s:m:c:l:U:o:r:I:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'v':
      if (v_col.read(optarg)) {
        v_col_op = 'o';
        break;
      }
      optarg_orig = optarg;
      split_line(optarg, parts, ",");
      if (strlen(parts[0]) == 1 && strchr("uUpPsSnNaAFEcCLM", *parts[0]))
        v_col_op = *parts[0];
      else
        error("invalid colouring", c);

      if (!((strchr("sS", (char)v_col_op) && parts.size() < 4) ||
            parts.size() < 2))
        error("too many comma separated parts", c);

      if (strchr("sS", v_col_op))
        v_sub_sym = optarg_orig.size() > 2 ? optarg_orig.substr(2) : "";
      break;

    case 'f':
      if (f_col.read(optarg)) {
        f_col_op = 'o';
        break;
      }
      optarg_orig = optarg;
      split_line(optarg, parts, ",");
      if (strlen(parts[0]) == 1 && strchr("uUpPsSnNaAkKgGcCLlM", *parts[0]))
        f_col_op = *parts[0];
      else
        error("invalid colouring", c);

      if (!((strchr("sS", (char)f_col_op) && parts.size() < 4) ||
            parts.size() < 2))
        error("too many comma separated parts", c);

      if (strchr("sS", f_col_op))
        f_sub_sym = optarg_orig.size() > 2 ? optarg_orig.substr(2) : "";
      break;

    case 'e':
      if (e_col.read(optarg)) {
        e_col_op = 'o';
        break;
      }
      optarg_orig = optarg;
      split_line(optarg, parts, ",");
      if (strlen(parts[0]) == 1 && strchr("uUpPsSkKFgGcCLldDM", *parts[0]))
        e_col_op = *parts[0];
      else
        error("invalid colouring", c);

      if (!((strchr("sS", (char)e_col_op) && parts.size() < 4) ||
            parts.size() < 2))
        error("too many comma separated parts", c);

      if (strchr("sS", e_col_op))
        e_sub_sym = optarg_orig.size() > 2 ? optarg_orig.substr(2) : "";
      break;

    case 'E':
      if (strlen(optarg) != 1 || !strchr("eiI", *optarg))
        error("edge type to color must be e, i or I");
      edge_type = *optarg;
      break;

    case 'l':
      lfile = optarg;
      break;

    case 'm':
      print_status_or_exit(read_colorings(clrngs, optarg), c);
      break;

    case 'r': {
      char r_elems;
      ColorValuesToRangeHsva col_proc;
      if (split_line(optarg, parts, ",") > 2)
        error("too many comma separated parts", c);
      print_status_or_exit(col_proc.init(parts[0]), c);
      if (parts.size() > 1) {
        if (strspn(parts[1], "vef") != strlen(parts[1]))
          error(msg_str("elements for colour ranges are '%s' must be "
                        "from v, e, and f",
                        optarg),
                c);
        r_elems = (strchr(parts[1], 'v') != nullptr) * ELEM_VERTS +
                  (strchr(parts[1], 'e') != nullptr) * ELEM_EDGES +
                  (strchr(parts[1], 'f') != nullptr) * ELEM_FACES;
      }
      else
        r_elems = ELEM_VERTS + ELEM_EDGES + ELEM_FACES;

      for (int i = 0; i < 3; i++)
        if ((r_elems & (1 << i)))
          col_procs[i] = col_proc;

      range_elems |= r_elems;
      break;
    }

    case 'I':
      if (strspn(optarg, "vef") != strlen(optarg))
        error(msg_str("elements to map are '%s' must be "
                      "from v, e, and f",
                      optarg),
              c);
      v2i_elems = (strchr(optarg, 'v') != nullptr) * ELEM_VERTS +
                  (strchr(optarg, 'e') != nullptr) * ELEM_EDGES +
                  (strchr(optarg, 'f') != nullptr) * ELEM_FACES;
      break;

    case 'U':
      selection = 0;
      prev_char_was_not = false;
      for (const char *p = optarg; *p; p++) {
        unsigned int new_select = 0;
        if (*p == '~') {
          if (prev_char_was_not)
            error("types cannot include repeated '~'", c);
          prev_char_was_not = true;
        }
        else if (*p == 'u')
          new_select = CV_UNSET;
        else if (*p == 'i')
          new_select = CV_INDEX;
        else if (*p == 'v')
          new_select = CV_VALUE;
        else if (*p == 'x')
          new_select = CV_INVISIBLE;
        else
          error(msg_str("invalid type character '%c'", *p), c);

        if (new_select)
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
  if (edge_type == 'x')
    edge_type = (!e_col_op || e_col_op == 'M') ? 'e' : 'a';

  if (argc - optind > 1)
    error("too many arguments");

  if (argc - optind == 1)
    ifile = argv[optind];
}

bool lights_read(const string &fname, Geometry *lights, char *errmsg)
{
  *errmsg = '\0';
  FILE *lfile = open_sup_file(fname.c_str(), "/col_lights/");
  if (lfile == nullptr) {
    if (errmsg)
      snprintf(errmsg, MSG_SZ, "could not open colour lights file \'%s\'",
               fname.c_str());
    return 0;
  }

  if (!lights->read(lfile))
    return false;
  bool no_color = false;
  for (unsigned int i = 0; i < lights->verts().size(); i++) {
    if (!lights->colors(VERTS).get(i).is_value()) {
      lights->colors(VERTS).set(i, Color(0.5, 0.5, 0.5));
      no_color = true;
    }
  }
  if (no_color)
    strcpy_msg(errmsg, "one or more missing colours, set to grey");

  return true;
}

inline unsigned int col_type(const Color &col)
{
  return CV_UNSET * !col.is_set() + CV_INDEX * col.is_index() +
         CV_VALUE * (col.is_value() && !col.is_invisible()) +
         CV_INVISIBLE * col.is_invisible();
}

void restore_orig_cols(Geometry &geom, Geometry &restore_geom,
                       unsigned int selection, unsigned int orig_edges_sz)
{
  for (unsigned int i = 0; i < restore_geom.verts().size(); i++) {
    Color col = restore_geom.colors(VERTS).get(i);
    if (!(col_type(col) & selection)) // restore colours of unselected elements
      geom.colors(VERTS).set(i, col);
  }

  vector<int> del_edges;
  for (unsigned int i = 0; i < geom.edges().size(); i++) {
    Color col;
    if (i < restore_geom.edges().size())
      col = restore_geom.colors(EDGES).get(i);
    if (!(col_type(col) & selection)) // restore cols of unselected elements
      geom.colors(EDGES).set(i, col);
    // implicit edges with unset colour were not selected to be coloured
    if (i >= orig_edges_sz && !geom.colors(EDGES).get(i).is_set())
      del_edges.push_back(i);
  }
  geom.del(EDGES, del_edges);

  for (unsigned int i = 0; i < restore_geom.faces().size(); i++) {
    Color col = restore_geom.colors(FACES).get(i);
    if (!(col_type(col) & selection)) // restore colours of unselected elements
      geom.colors(FACES).set(i, col);
  }
}

bool warn_if_not_subgroup(const ProgramOpts &opts, char opt,
                          const Symmetry &whole, const Symmetry &part)
{
  Transformations min;
  if (whole.get_trans().size() % part.get_trans().size() || // Lagrange
      min.intersection(whole.get_trans(), part.get_trans()).size() !=
          part.get_trans().size()) {
    opts.warning("group is not subgroup of model", opt);
    return false; // is not subgroup
  }
  return true; // is subgroup
}

int main(int argc, char *argv[])
{
  o_col_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);

  // read lights
  char errmsg[MSG_SZ];
  Geometry lights;
  if (opts.lfile != "") {
    if (!lights_read(opts.lfile, &lights, errmsg))
      opts.error(errmsg, 'l');
    if (*errmsg)
      opts.warning(errmsg, 'l');
  }

  // store original edges
  unsigned int orig_edges_sz = geom.edges().size();
  map<vector<int>, Color> expl_edges;
  if (opts.e_col_op && strchr("pP", opts.e_col_op) &&
      !strchr("iI", opts.edge_type))
    opts.edge_type = 'i';
  if (opts.edge_type == 'a')
    geom.add_missing_impl_edges();
  else if (strchr("iI", opts.edge_type)) {
    for (unsigned int i = 0; i < geom.edges().size(); i++)
      expl_edges[geom.edges(i)] = geom.colors(EDGES).get(i);
    geom.clear(EDGES);
    geom.add_missing_impl_edges();
  }

  Geometry store_geom;
  if (opts.selection)
    store_geom = geom;

  // Get symmetry if necessary
  Symmetry sym;
  vector<vector<set<int>>> sym_equivs;
  if ((opts.f_col_op && strchr("sS", opts.f_col_op)) ||
      (opts.e_col_op && strchr("sS", opts.e_col_op)) ||
      (opts.v_col_op && strchr("sS", opts.v_col_op))) {
    sym.init(geom, &sym_equivs);
    opts.v_equivs = sym_equivs[0];
    opts.e_equivs = sym_equivs[1];
    opts.f_equivs = sym_equivs[2];

    Symmetry sub;
    if (opts.v_col_op && strchr("sS", opts.v_col_op)) {
      opts.print_status_or_exit(sym.get_sub_sym(opts.v_sub_sym, &sub), 'v');
      get_equiv_elems(geom, sub.get_trans(), &sym_equivs);
      opts.v_equivs = sym_equivs[0];
    }
    if (opts.e_col_op && strchr("sS", opts.e_col_op)) {
      opts.print_status_or_exit(sym.get_sub_sym(opts.e_sub_sym, &sub), 'e');
      get_equiv_elems(geom, sub.get_trans(), &sym_equivs);
      opts.e_equivs = sym_equivs[1];
    }
    if (opts.f_col_op && strchr("sS", opts.f_col_op)) {
      opts.print_status_or_exit(sym.get_sub_sym(opts.f_sub_sym, &sub), 'f');
      get_equiv_elems(geom, sub.get_trans(), &sym_equivs);
      opts.f_equivs = sym_equivs[2];
    }
  }

  Coloring &fc = opts.clrngs[FACES];
  fc.set_geom(&geom);
  if (opts.f_col_op) {
    char op = opts.f_col_op;
    ColorMap *cmap = nullptr;
    if (fc.get_cmaps().size() == 0) {
      if (strchr("GgCc", op))
        cmap = colormap_from_name("range");
      else
        cmap = colormap_from_name("spread");
    }
    if (cmap)
      fc.add_cmap(cmap);
    if (op == 'o')
      fc.f_one_col(opts.f_col);
    else if (strchr("uU", op))
      fc.f_unique(op == 'U');
    else if (strchr("pP", op))
      fc.f_proper(op == 'P');
    else if (strchr("sS", op))
      fc.f_sets(sym_equivs[2], op == 'S');
    else if (strchr("nN", op))
      fc.f_sides(op == 'N');
    else if (strchr("aA", op))
      fc.f_avg_angle(op == 'A');
    else if (strchr("kK", op))
      fc.f_parts(op == 'K');
    else if (strchr("Gg", op))
      fc.f_normal(op == 'G');
    else if (strchr("Cc", op))
      fc.f_centroid(op == 'C');
    else if (strchr("L", op))
      fc.f_lights(lights);
    else if (strchr("l", op))
      fc.f_lights2(lights);
    else if (strchr("M", op))
      fc.f_apply_cmap();
  }

  Coloring &ec = opts.clrngs[EDGES];
  ec.set_geom(&geom);
  if (opts.e_col_op) {
    char op = opts.e_col_op;
    ColorMap *cmap = nullptr;
    if (ec.get_cmaps().size() == 0) {
      if (strchr("GgCc", op))
        cmap = colormap_from_name("range");
      else
        cmap = colormap_from_name("spread");
    }
    if (cmap)
      ec.add_cmap(cmap);
    if (op == 'o')
      ec.e_one_col(opts.e_col);
    else if (strchr("uU", op))
      ec.e_unique(op == 'U');
    else if (strchr("pP", op))
      ec.e_proper(op == 'P');
    else if (strchr("sS", op))
      ec.e_sets(sym_equivs[1], op == 'S');
    else if (strchr("kK", op))
      ec.e_parts(op == 'K');
    else if (strchr("Gg", op))
      ec.e_direction(op == 'G');
    else if (strchr("Cc", op))
      ec.e_mid_point(op == 'C');
    else if (strchr("Dd", op))
      ec.e_vector(op == 'D');
    else if (strchr("L", op))
      ec.e_lights(lights);
    else if (strchr("l", op))
      ec.e_dir_lights(lights);
    else if (strchr("M", op))
      ec.e_apply_cmap();
  }

  Coloring &vc = opts.clrngs[VERTS];
  vc.set_geom(&geom);
  if (opts.v_col_op) {
    char op = opts.v_col_op;
    ColorMap *cmap = nullptr;
    if (vc.get_cmaps().size() == 0) {
      if (strchr("Cc", op))
        cmap = colormap_from_name("range");
      else
        cmap = colormap_from_name("spread");
    }
    if (cmap)
      vc.add_cmap(cmap);
    if (op == 'o')
      vc.v_one_col(opts.v_col);
    else if (strchr("uU", op))
      vc.v_unique(op == 'U');
    else if (strchr("pP", op))
      vc.v_proper(op == 'P');
    else if (strchr("sS", op))
      vc.v_sets(sym_equivs[0], op == 'S');
    else if (strchr("nN", op))
      vc.v_order(op == 'N');
    else if (strchr("aA", op))
      vc.v_avg_angle(op == 'A');
    else if (strchr("cC", op))
      vc.v_position(op == 'C');
    else if (strchr("L", op))
      vc.v_lights(lights);
    else if (strchr("M", op))
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
  if (opts.range_elems & (ELEM_VERTS))
    opts.col_procs[0].apply(geom.colors(VERTS).get_properties());
  if (opts.range_elems & (ELEM_EDGES))
    opts.col_procs[1].apply(geom.colors(EDGES).get_properties());
  if (opts.range_elems & (ELEM_FACES))
    opts.col_procs[2].apply(geom.colors(FACES).get_properties());

  // Average colour values from adjoining elements after converting
  // index numbers
  if (opts.e_col_op == 'F')
    ec.e_face_color();
  if (opts.v_col_op == 'F')
    vc.v_face_color();
  else if (opts.v_col_op == 'E')
    vc.v_edge_color();

  // Finally convert to index numbers
  color_vals_to_idxs(geom, opts.v2i_elems);

  if (opts.selection)
    restore_orig_cols(geom, store_geom, opts.selection, orig_edges_sz);

  // restore original edges
  map<vector<int>, Color>::iterator ei;
  if (opts.edge_type == 'I') {
    for (ei = expl_edges.begin(); ei != expl_edges.end(); ++ei)
      geom.add_edge(ei->first, ei->second);
  }
  else if (opts.edge_type == 'i') {
    for (ei = expl_edges.begin(); ei != expl_edges.end(); ++ei)
      if (find(geom.edges().begin(), geom.edges().end(), ei->first) ==
          geom.edges().end())
        geom.add_edge(ei->first, ei->second);
  }

  opts.write_or_error(geom, opts.ofile);

  return 0;
}
