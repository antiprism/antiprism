/*
   Copyright (c) 2017, Roger Kaufman

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
   Name: miller.cc
   Description: stellations test program
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"

#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>

using std::map;
using std::set;
using std::string;
using std::vector;

using namespace anti;

class miller_opts : public ProgramOpts {
public:
  string ifile;
  string ofile;

  string output_parts;
  bool merge_faces;
  bool rebuild_compound_model;
  bool list_polys;

  char vertex_coloring_method;
  char edge_coloring_method;
  char face_coloring_method;

  Color vertex_color;
  Color edge_color;
  Color face_color;

  string map_string;
  int face_opacity;

  double epsilon;

  miller_opts()
      : ProgramOpts("miller"), output_parts("s"), merge_faces(false),
        rebuild_compound_model(false), list_polys(false),
        vertex_coloring_method('\0'), edge_coloring_method('\0'),
        face_coloring_method('\0'), vertex_color(Color::invisible),
        edge_color(Color::invisible), face_color(Color()),
        map_string("compound"), face_opacity(-1), epsilon(0)
  {
  }
  void process_command_line(int argc, char **argv);
  void usage();
};

void miller_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] input

Millers 59 Icosahedra Stellations. Plus additional stellations since discovered
input may be Miller list number from 1 to 75. Or m_string where string consists
of one or more cell names: A,B,C,D,E,F,G,H,e1,f1,f1',g1,e2,f2,g2  e.g m_De1f1g1
model string can be followed by I or Ih symmetry. e.g. m_e1f1',I
std_ may precede input string to output a raw model 

Options
%s
  -L        list models only
  -M        merge stellation facelets
  -r        rebuild compound model to separate vertices
  -O <args> output s - stellation, d - diagram (default: s)
  -l <lim>  minimum distance for unique vertex locations as negative exponent
               (default: %d giving %.0e)
  -o <file> write output to file (default: write to standard output)

Coloring Options (run 'off_util -H color' for help on color formats)
  -F <opt>  face coloring method. d - from diagram, s - symmetry
               c - color by compound (default: 255,193,37, if compound then c)
               C - face/face connection count using map n colors
               keyword: none - sets no color
  -E <col>  edge color. f - from faces (default: invisible)
               C - edge/face connection count using map n colors
               keyword: none - sets no color
  -V <col>  vertex color.  e - from edges (default: invisible)
               keyword: none - sets no color
  -T <tran> face transparency. valid range from 0 (invisible) to 255 (opaque)
  -m <maps> color maps. stellation diagram or face symmetry (default: compound)

)",
          prog_name(), help_ver_text, int(-log(::epsilon) / log(10) + 0.5),
          ::epsilon);
}

void miller_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  int sig_compare = INT_MAX;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hLMrO:V:E:F:T:m:l:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'L':
      list_polys = true;
      break;

    case 'M':
      merge_faces = true;
      break;

    case 'r':
      rebuild_compound_model = true;
      break;

    case 'O':
      if (strspn(optarg, "sd") != strlen(optarg))
        error(msg_str("output parts are '%s' must be any or all from "
                      "s, d",
                      optarg),
              c);
      output_parts = optarg;
      break;

    case 'V':
      if (strchr("e", *optarg))
        vertex_coloring_method = *optarg;
      else
        print_status_or_exit(vertex_color.read(optarg), c);
      break;

    case 'E':
      if (strlen(optarg) == 1 && strchr("fC", int(*optarg)))
        edge_coloring_method = *optarg;
      else
        print_status_or_exit(edge_color.read(optarg), c);
      break;

    case 'F':
      if (strlen(optarg) == 1 && strchr("dscC", int(*optarg)))
        face_coloring_method = *optarg;
      else {
        print_status_or_exit(face_color.read(optarg), c);
        face_coloring_method = 'x';
      }
      break;

    case 'T':
      print_status_or_exit(read_int(optarg, &face_opacity), c);
      if (face_opacity < 0 || face_opacity > 255) {
        error("face transparency must be between 0 and 255", c);
      }
      break;

    case 'm':
      map_string = optarg;
      break;

    case 'l':
      print_status_or_exit(read_int(optarg, &sig_compare), c);
      if (sig_compare < 0) {
        warning("limit is negative, and so ignored", c);
      }
      if (sig_compare > DEF_SIG_DGTS) {
        warning("limit is very small, may not be attainable", c);
      }
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (argc - optind > 1)
    error("too many arguments");

  if (argc - optind == 1)
    ifile = argv[optind];

  epsilon = (sig_compare != INT_MAX) ? pow(10, -sig_compare) : ::epsilon;
}

// copy of code from stellate.cc
void apply_transparency(Geometry &geom, const miller_opts &opts,
                        int alternate_opacity = -1)
{
  int face_opacity =
      (alternate_opacity != -1) ? alternate_opacity : opts.face_opacity;
  if (face_opacity > -1) {
    ColorValuesToRangeHsva valmap(msg_str("A%g", (double)face_opacity / 255));
    valmap.apply(geom, FACES);

    for (const auto &kp : geom.colors(FACES).get_properties()) {
      if (kp.second.is_index()) {
        opts.warning("map indexes cannot be made transparent", 'T');
        break;
      }
    }

    // check if some faces are not set
    if (geom.colors(FACES).get_properties().size() < geom.faces().size())
      opts.warning("unset faces cannot be made transparent", 'T');
  }
}

void color_stellation(Geometry &stellation, const miller_opts &opts)
{
  // set color map
  Coloring clrng(&stellation);
  ColorMap *cmap = colormap_from_name(opts.map_string.c_str());
  clrng.add_cmap(cmap);

  // in case of face coloring C
  Geometry kis;

  // stellation is built with color from the diagram
  if (opts.face_coloring_method != 'd') {
    if (!opts.face_coloring_method)
      // if no color specified, clear faces
      stellation.colors(FACES).clear();
    else if (opts.face_coloring_method == 's') {
      // color faces by symmetry
      Symmetry sym;
      vector<vector<set<int>>> sym_equivs;
      sym.init(stellation, &sym_equivs);
      clrng.f_sets(sym_equivs[2], true);
    }
    else if (opts.face_coloring_method == 'c')
      clrng.f_parts(true);
    else if (opts.face_coloring_method == 'C') {
      wythoff_make_tiling(kis, stellation, "k", true, false);
      // remove digons
      vector<int> dels;
      for (unsigned int i = 0; i < kis.faces().size(); i++) {
        if (kis.faces(i).size() < 3)
          dels.push_back((int)i);
      }
      kis.del(FACES, dels);
      kis.orient(1);

      // make new verts and edges invisible
      kis.add_missing_impl_edges();
      for (unsigned int i = 0; i < kis.verts().size(); i++) {
        int v_idx = find_vert_by_coords(stellation, kis.verts()[i], epsilon);
        if (v_idx == -1) {
          kis.colors(VERTS).set(i, Color::invisible);
          vector<int> edge_idx = find_edges_with_vertex(kis.edges(), i);
          for (unsigned int j = 0; j < edge_idx.size(); j++)
            kis.colors(EDGES).set(edge_idx[j], Color::invisible);
        }
      }
      // the old faces are cleared and kis faces added
      stellation.clear(FACES);
      stellation.append(kis);
      int blend_type = 1; // first color, invisible edges stay
      merge_coincident_elements(stellation, "vef", blend_type, epsilon);

      for (unsigned int i = 0; i < stellation.faces().size(); i++) {
        vector<int> face = stellation.faces()[i];
        unsigned int fsz = face.size();
        // face to face
        // connections with invisible faces are ignored
        int connections = 0;
        for (unsigned int j = 0; j < fsz; j++) {
          int v1 = face[j];
          int v2 = face[(j + 1) % fsz];
          vector<int> edge = make_edge(v1, v2);
          vector<int> face_idx = find_faces_with_edge(stellation.faces(), edge);
          int edge_no = find_edge_in_edge_list(stellation.edges(), edge);
          if (!(stellation.colors(EDGES).get(edge_no)).is_invisible())
            connections += face_idx.size();
        }
        stellation.colors(FACES).set(i, cmap->get_col(connections));
      }
    }
    else
      // use color selected
      clrng.f_one_col(opts.face_color);
  }

  // if edges take color from faces (if faces none, clear edges)
  if (opts.edge_coloring_method == 'f') {
    // if face colors is none
    if (!opts.face_coloring_method)
      stellation.colors(EDGES).clear();
    else
      clrng.e_face_color();
  }
  else if (opts.edge_coloring_method == 'C') {
    auto efpairs = stellation.get_edge_face_pairs(false);
    for (const auto &edge : stellation.edges()) {
      vector<int> faces = efpairs[edge];
      int i = find_edge_in_edge_list(stellation.edges(), edge);
      if (i > -1) {
        if (!(stellation.colors(EDGES).get(i)).is_invisible()) {
          unsigned int connections = faces.size();
          stellation.colors(EDGES).set(i, cmap->get_col(connections));
        }
      }
    }
  }
  else
    // use color selected
    clrng.e_one_col(opts.edge_color);

  // vertices taking edge coloring from stellation is the default
  // if not using default, sample colors again
  if (opts.vertex_coloring_method == 'e') {
    // vertices take color from edges
    clrng.v_edge_color();
  }
  else
    // use color selected
    clrng.v_one_col(opts.vertex_color);

  // if using face connection coloring
  // vertices from kis must be made invisible in stellation
  if (opts.face_coloring_method == 'C') {
    for (unsigned int i = 0; i < kis.verts().size(); i++) {
      int v_idx = find_vert_by_coords(stellation, kis.verts()[i], epsilon);
      if (v_idx != -1) {
        if ((kis.colors(VERTS).get(i)).is_invisible())
          stellation.colors(VERTS).set(v_idx, Color::invisible);
      }
    }
  }

  // set transparency
  if (opts.face_opacity > -1)
    apply_transparency(stellation, opts);
}

struct MillerItem {
  const char *cell_string;
  const char *sub_sym;
  bool remove_inline_verts;
};

class Miller {
private:
  MillerItem *Miller_items;
  int last_M;

public:
  Miller();
  int get_poly(anti::Geometry &geom, int sym, string cell_str, string sym_str,
               map<int, Geometry> &diagrams, miller_opts &opts);
  int get_last_M() { return last_M; }
  void list_polys();
};

// Millers 59 Stellations of the Icosahedron

/* notes:
32 Full Icosahedral Symmetry
27 Enantiomeric (have mirror images, i.e. I symmetry)
59 Total

1   Icosahedron
2   Triakis Icosahedron
3   Compound of 5 Octahedra
7   Great Icosahedron
8   Echidnahedron
16  Only model which is disconnected
22  Compound of 10 Tetrahedra
26  Excavated Dodecahedron
30  Great Triambic Icosahedron
??  Deltahedron of 60 faces
47  Compound of 5 Tetrahedra
*/

// clang-format off
// Table of miller Models which are not Uniforms, 75 items
// Number, Cells, symmetry, remove inline vertices
// f1 left handed form represented by f3
MillerItem miller_item_list[] = {
   // Full symmetry
   { "A",          "Ih", true  }, // w4
   { "B",          "Ih", true  }, // w26
   { "C",          "Ih", true  }, // w23
   { "D",          "Ih", true  },
   { "E",          "Ih", true  },
   { "F",          "Ih", true  }, // w27
   { "G",          "Ih", true  }, // w41
   { "H",          "Ih", true  }, // w42
   { "e1",         "Ih", true  }, // w37
   { "f1",         "Ih", false },
   { "g1",         "Ih", true  }, // w29
   { "e1f1",       "Ih", false },
   { "e1f1g1",     "Ih", true  },
   { "f1g1",       "Ih", false },
   { "e2",         "Ih", true  }, // cells
   { "f2",         "Ih", true  },
   { "g2",         "Ih", false }, // cells, faces not merged
   { "e2f2",       "Ih", true  },
   { "e2f2g2",     "Ih", true  },
   { "f2g2",       "Ih", true  }, // w30
   { "De1",        "Ih", true  }, // w32
   { "Ef1",        "Ih", true  }, // w25
   { "Fg1",        "Ih", true  }, // w31
   { "De1f1",      "Ih", true  }, // cells, faces not merged
   { "De1f1g1",    "Ih", true  },
   { "Ef1g1",      "Ih", true  }, // w28
   { "De2",        "Ih", true  },
   { "Ef2",        "Ih", true  },
   { "Fg2",        "Ih", true  }, // w33
   { "De2f2",      "Ih", true  }, // w34
   { "De2f2g2",    "Ih", true  },
   { "Ef2g2",      "Ih", true  }, // cells
   // Enantiomeric
   { "f1'",        "I",  true  }, // w35
   { "e1f1'",      "I",  true  }, // w36
   { "De1f1'",     "I",  false },
   { "f1'g1",      "I",  false },
   { "e1f1'g1",    "I",  false }, // w39
   { "De1f1'g1",   "I",  false },
   { "f1'g2",      "I",  false },
   { "e1f1'g2",    "I",  true  }, // cells, faces not merged
   { "De1f1'g2",   "I",  true  }, // cells, faces not merged
   { "f1'f2g2",    "I",  true  },
   { "e1f1'f2g2",  "I",  true  },
   { "De1f1'f2g2", "I",  true  }, // cells
   { "e2f1'",      "I",  false }, // w40
   { "De2f1'",     "I",  true  },
   { "Ef1'",       "I",  true  }, // w24
   { "e2f1'g1",    "I",  false },
   { "De2f1'g1",   "I",  true  },
   { "Ef1'g1",     "I",  true  },
   { "e2f1'f2",    "I",  true  }, // w38
   { "De2f1'f2",   "I",  true  },
   { "Ef1'f2",     "I",  true  },
   { "e2f1'f2g1",  "I",  true  },
   { "De2f1'f2g1", "I",  true  },
   { "Ef1'f2g1",   "I",  true  },
   { "e2f1'f2g2",  "I",  true  },
   { "De2f1'f2g2", "I",  true  },
   { "Ef1'f2g2",   "I",  true  },
   // begin the true lost stellations
   { "Be1",        "Ih", true  },
   { "Ce2",        "Ih", true  },
   { "Df1",        "I",  true  },
   { "Af2",        "Ih", true  },
   { "Df2",        "Ih", true  },
   { "De1f1'f2",   "I",  true  },
   { "De1g1",      "Ih", true  },
   { "Af2g1",      "Ih", true  },
   { "De1f1'f2g1", "I",  true  },
   { "Af2g2",      "Ih", true  },
   // begin the candidate lost stellations
   { "Be2",        "Ih", true  },
   { "De2f2",      "Ih", true  },
   { "e1g1",       "Ih", true  },
   { "Ef1g1",      "Ih", true  },
   { "Cf2g1",      "Ih", true  },
   { "ACDf2g1",    "Ih", true  },
};
// clang-format on

void Miller::list_polys()
{
  int l = get_last_M();
  for (int i = 0; i < l; i++) {
    int j = i + 1;
    if (j == 1)
      fprintf(stderr, "Full Symmetry\n");
    else if (j == 33)
      fprintf(stderr, "\nEnantiomeric\n");
    else if (j == 60)
      fprintf(stderr, "\nLost Stellations\n");
    else if (j == 70)
      fprintf(stderr, "\nCandidate Lost Stellations\n");

    fprintf(stderr, "%2d) %-10s %-2s\n", j, Miller_items[i].cell_string,
            Miller_items[i].sub_sym);
  }
}

Miller::Miller()
{
  Miller_items = miller_item_list;
  last_M = sizeof(miller_item_list) / sizeof(miller_item_list[0]);
}

// if cell string will not decode, return empty result string
// f1 left handed form represented by f3
vector<string> decode_cell_string(string cell_str)
{
  // clang-format off
  map<string, string> cell;
  cell["A"]  = "0,18";             // Ih
  cell["B"]  = "0,46";
  cell["C"]  = "0,30";
  cell["D"]  = "0,17,29,44";       // must be Ih symmetric
  cell["E"]  = "0,14,16,31,43,62"; // must be Ih symmetric
  cell["F"]  = "0,32,34,57";
  cell["G"]  = "0,11,33"; 
  cell["H"]  = "0,8,12";
  cell["e1"] = "0,14,43,44";       // must be Ih symmetric
  cell["f1"] = "0,15,31,43,61";
  cell["g1"] = "0,11,32,61";
  cell["e2"] = "0,16,29,31,40,58"; // must be Ih symmetric
  cell["f2"] = "0,16,34";          // Ih
  cell["g2"] = "0,10,15,33,34,57"; // must be Ih symmetric

  // opposite hand. use prime since no italics is possible
  cell["f1'"] = "0,14,32,57,62";

  // right (normal) hand
  cell["A+"]  = "0,18";
  cell["B+"]  = "0,46";
  cell["C+"]  = "0,30";
  cell["D+"]  = "0,29,44";
  cell["E+"]  = "0,16,31,43";
  cell["F+"]  = "0,32,34,57";
  cell["G+"]  = "0,11,33";
  cell["H+"]  = "0,8,12";
  cell["e1+"] = "0,43,44";
  cell["f1+"] = "0,15,31,43,61";
  cell["g1+"] = "0,11,32";
  cell["e2+"] = "0,16,29,31";
  cell["f2+"] = "0,16,34";
  cell["g2+"] = "0,10,33,34,57";

  // left hand
  cell["A+"]  = "0,18"; // same
  cell["B-"]  = "0,28";
  cell["C-"]  = "0,45";
  cell["D-"]  = "0,17,44";
  cell["E-"]  = "0,14,16,62";
  cell["F-"]  = "0,15,34,61";
  cell["G-"]  = "0,10,11";
  cell["H-"]  = "0,9,12";
  cell["e1-"] = "0,14,44";
  cell["f1-"] = "0,14,32,57,62";
  cell["g1-"] = "0,11,61";
  cell["e2-"] = "0,20,51,52";
  cell["f2-"] = "0,16,34"; // same
  cell["g2-"] = "0,5,49,55,64";
// clang-format off

  // code
  vector<string> diagram_list_strings;

  int len = cell_str.length();
  int len_cnt = len;
  while (len_cnt > 0) {
    string cell_substr = cell_str.substr(len-len_cnt);
    string token;
    int tok_len = 0;
    size_t pos = cell_substr.find_first_of("ABCDEFGH");
    if (pos != string::npos) {
      if (pos+1 <= cell_substr.size())
        tok_len = (cell_substr.substr(pos+1,1).find_first_of("'+-") != string::npos) ? 2 : 1;
    }
    else {
      pos = cell_substr.find_first_of("efg");
      if (pos != string::npos) {
        if (pos+2 <= cell_substr.size())
          tok_len = (cell_substr.substr(pos+2,1).find_first_of("'+-") != string::npos) ? 3 : 2;
      }
    }
    if (tok_len) {
      if (pos <= cell_substr.size())
        token = cell_substr.substr(pos,tok_len);
      len_cnt-=tok_len;
    }

//fprintf(stderr,"token = '%s'\n", token.c_str());
//fprintf(stderr,"cell[token] = %s\n", cell[token].c_str());

    // if not found, return empty string list
    if (!cell[token].length()) {
      diagram_list_strings.clear();
      break;
    }
    else {
      // append face number
      diagram_list_strings.push_back(cell[token]);
    }
  }

  return (diagram_list_strings);
}

int Miller::get_poly(Geometry &geom, int sym, string cell_str, string sym_str, map<int, Geometry> &diagrams, miller_opts &opts)
{
  geom.read_resource("ico");

  // decode the cell string
  if (!cell_str.length()) {
    cell_str = Miller_items[sym].cell_string;
    sym_str = Miller_items[sym].sub_sym;
  }
  else {
    if (!sym_str.length())
      sym_str = "Ih";
  }
  fprintf(stderr,"%s,%s\n", cell_str.c_str(), sym_str.c_str()); 

  vector<string> diagram_list_strings = decode_cell_string(cell_str);
  // if list is empty, code string was invalid
  unsigned int sz = diagram_list_strings.size();
  if (!sz)
    return -1;

  vector<vector<int> > idx_lists(sz);

  // data sizes are verified
  for (unsigned int i = 0; i < diagram_list_strings.size(); i++) {
    if (!diagram_list_strings[i].length())
      continue;

    read_idx_list((char *)diagram_list_strings[i].c_str(), idx_lists[i], INT_MAX, false);

    // stellation face index is in the first position
    int stellation_face_idx = idx_lists[i][0];

    // construct the diagrams
    if (!diagrams[stellation_face_idx].verts().size())
      diagrams[stellation_face_idx] = make_stellation_diagram(geom, stellation_face_idx, sym_str);
  }

  bool merge_faces = opts.merge_faces;
  bool remove_inline_verts = Miller_items[sym].remove_inline_verts;
  bool split_pinched = true;
  bool resolve_faces = true;
  bool remove_multiples = true;

  geom = make_stellation(geom, diagrams, idx_lists, sym_str, merge_faces, remove_inline_verts,
                         split_pinched, resolve_faces, remove_multiples, opts.map_string, opts.epsilon);

  vector<vector<int>> idx_lists_full =
      lists_full(diagrams, idx_lists, remove_multiples);

  if (opts.output_parts.find_first_of("d") != string::npos) {
    // color diagrams
    for (auto const &key1 : diagrams) {
      Coloring clrng(&diagrams[key1.first]);
      ColorMap *cmap = colormap_from_name(opts.map_string.c_str());
      clrng.add_cmap(cmap);
      clrng.f_apply_cmap();
    }

    // set large diagram transparency to 1/4th
    for (auto const &key1 : diagrams)
      apply_transparency(diagrams[key1.first], opts, 64);

    for (unsigned int i = 0; i < idx_lists_full.size(); i++) {
      // stellation face index is in the first position
      int stellation_face_idx = idx_lists_full[i][0];

      // set selected faces back to full visibility
      for (unsigned int j = 1; j < idx_lists_full[i].size(); j++) {
        int face_no = idx_lists_full[i][j];
        Color c = diagrams[stellation_face_idx].colors(FACES).get(face_no);
        c.set_alpha(255);
        diagrams[stellation_face_idx].colors(FACES).set(face_no, c);
      }
    }
  }

  return 1;
}

int make_resource_miller(Geometry &geom, string name, bool is_std, miller_opts &opts,
                         char *errmsg)
{
  *errmsg = '\0';

  int sym_no = 0;
  // check if it is just the index number, if so format as m%d
  if (read_int(name.c_str(), &sym_no))
    name = "m" + std::to_string(sym_no);
  else
  if (name.size() < 2 || !strchr("mM", name[0]) ||
      name.find('.') != string::npos)
    return -1; // not miller name (the "." indicates a likely local file)
               // so the name is not handled

  string cell_str;
  string sym_str;

  Miller mill;
  if (read_int(name.c_str() + 1, &sym_no)) {
    sym_no--;
    if (sym_no < 0 || sym_no >= mill.get_last_M()) {
      strcpy_msg(errmsg, "miller stellation number out of range");
      return 1; // fail
    }
  }
  else
  // if the next character is an underscore, try it as cell names
  if (name[1] == '_') {
    cell_str = name.substr(2);
    size_t pos = cell_str.find_first_of(",");
    if (pos != string::npos) {
      sym_str = cell_str.substr(pos+1);
      if (sym_str.length()) {
        if (sym_str != "I" && sym_str != "Ih") {
          strcpy_msg(errmsg, "miller cell name symmetry must be I or Ih");
          return 1; // fail
        }
      }
      cell_str = cell_str.substr(0,pos);
    }
  }
  else
    return -1; // not a string of cell names

  map<int, Geometry> diagrams;
  int ret = mill.get_poly(geom, sym_no, cell_str, sym_str, diagrams, opts);
  if (ret < 1)
    return 1; // fail

  if (opts.rebuild_compound_model) {
    rebuild_compound(geom);
    geom.add_missing_impl_edges();
  }

  // if is_std, have to strip built in color
  if (is_std) {
    geom.colors(VERTS).clear();
    geom.colors(EDGES).clear();
    geom.colors(FACES).clear();
  }
  else if (!opts.face_coloring_method) {
    Coloring clrng(&geom);
    ColorMap *cmap = colormap_from_name("compound");
    clrng.add_cmap(cmap);

    GeometryInfo info(geom);
    if (info.num_parts() > 1)
      // color by compound
      clrng.f_parts(true);
    else
      // faces to one color because model is not completely merged
      // goldenrod1: Color(255,193,37) ... from George Hart, Color(1.0,0.7,0.1)
      clrng.f_one_col(Color(255,193,37));

    // edges can be confusing. set them to invisible
    Coloring(&geom).vef_one_col(Color::invisible, Color::invisible, Color());

    // patch to keep color_stellation from firing
    opts.face_coloring_method = 'd';
  }

  color_stellation(geom, opts);

  // if only diagram is being output, clear geom
  if (opts.output_parts.find_first_of("s") == string::npos)
    geom.clear_all();

  // append diagrams if being output
  if (opts.output_parts.find_first_of("d") != string::npos) {
    for (auto const &key1 : diagrams)
      geom.append(diagrams[key1.first]);

    // if diagram only, move the diagram to the front
    if (opts.output_parts == "d") {
      Vec3d face_normal = face_norm(geom.verts(), geom.faces(0));
      if (vdot(face_normal, geom.face_cent(0)) < 0)
        face_normal *= -1.0;
      Trans3d trans = Trans3d::rotate(face_normal, Vec3d(0, 0, 1));
      geom.transform(trans);
    }
  }

  return 0; // name found
}

int try_miller(Geometry &geom, miller_opts &opts, char *errmsg)
{
  string name = opts.ifile;
  bool is_std = (name.size() > 3 && name.substr(0, 4) == "std_");
  if (is_std)
    name = name.substr(4);
  int idx = make_resource_miller(geom, name, is_std, opts, errmsg);
  return(idx);
}

int main(int argc, char *argv[])
{
  miller_opts opts;
  opts.process_command_line(argc, argv);

  if (opts.list_polys) {
    Miller mill;
    mill.list_polys();
    exit(0);
  }

  char errmsg[MSG_SZ] = {0};

  Geometry geom;
  if (try_miller(geom, opts, errmsg))
    if (*errmsg)
      opts.error(errmsg);

  opts.write_or_error(geom, opts.ofile);

  return 0;
}
