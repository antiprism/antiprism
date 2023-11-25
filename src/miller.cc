/*
   Copyright (c) 2017-2023, Roger Kaufman

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
#include "color_common.h"

#include <cstdio>
#include <map>
#include <string>
#include <vector>

using std::map;
using std::string;
using std::vector;

using namespace anti;

class miller_opts : public ProgramOpts {
public:
  string ifile;
  string ofile;

  string output_parts = "s";           // s - stellation d - diagram
  bool merge_faces = false;            // do a planar merge on facelet
  bool rebuild_compound_model = false; // rebuild with seperate constituents

  bool list_polys = false; // output the list of models to the screen

  double eps = anti::epsilon;

  // map name is used to pass to stellation function
  string map_string = "compound";

  // internal colorings are needed for local coloring
  OffColor off_color = OffColor(map_string);

  int opacity[3] = {-1, -1, -1}; // transparency from 0 to 255, for v,e,f

  miller_opts() : ProgramOpts("miller") {}

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

Options
%s
  -L        list models only
  -l <lim>  minimum distance for unique vertex locations as negative exponent
              (default: %d giving %.0e)
  -o <file> write output to file (default: write to standard output)

Program Options
  -M        merge stellation facelets (for cell name strings only)
  -r        rebuild compound model to separate vertices

Scene Options
  -O <args> output s - stellation, d - diagram (default: s)

Coloring Options (run 'off_util -H color' for help on color formats)
keyword: none - sets no color
  -F <col>  color the faces according to: (default: 255,193,37, if compound, k)
              a color value - apply to all faces
              k - sets of faces connected by face edges (compounds)
              s - symmetric colouring [,sub_group,conj_type]
              q - from stellation diagram
              h - face/face connection count (model is altered with kis)
  -E <col>  color the edges according to: (default: invisible)
              a color value - apply to all edges
              f - color with average adjacent face color
              n - colour by number of faces connected to each edge
  -V <col>  color the vertices according to: (default: invisible)
              a color value - apply to all vertices
              e - color with average adjacent edge color
              f - color with average adjacent face color
              n - color by order of vertex
  -T <t,e>  transparency. from 0 (invisible) to 255 (opaque). element is any
            or all of, v - vertices, e - edges, f - faces, a - all (default: f)
  -m <maps> a comma separated list of color maps used to transform color
            indexes (default: compound), a part consisting of letters from
            v, e, f, selects the element types to apply the map list to
            (default 'vef'). use map name of 'index' to output index numbers
              compound:   yellow,red,darkgreen,blue,magenta,cyan,darkorange1
            
)",
          prog_name(), help_ver_text, int(-log(anti::epsilon) / log(10) + 0.5),
          anti::epsilon);
}

void miller_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;
  int num;

  Split parts;
  Color col;

  off_color.set_e_col(Color::invisible);
  off_color.set_v_col(Color::invisible);

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
      if (col.read(optarg)) {
        off_color.set_v_col(col);
        break;
      }
      parts.init(optarg, ",");
      if (off_color.v_op_check((char *)parts[0], "efn"))
        off_color.set_v_col_op(*parts[0]);
      else
        error("invalid coloring", c);

      if (!((strchr("sS", off_color.get_v_col_op()) && parts.size() < 4) ||
            parts.size() < 2))
        error("too many comma separated parts", c);

      if (strchr("sS", off_color.get_v_col_op()))
        off_color.set_v_sub_sym(strlen(optarg) > 2 ? optarg + 2 : "");
      break;

    case 'E':
      if (col.read(optarg)) {
        off_color.set_e_col(col);
        break;
      }
      parts.init(optarg, ",");
      if (off_color.e_op_check((char *)parts[0], "fn"))
        off_color.set_e_col_op(*parts[0]);
      else
        error("invalid coloring", c);

      if (!((strchr("sS", off_color.get_e_col_op()) && parts.size() < 4) ||
            parts.size() < 2))
        error("too many comma separated parts", c);

      if (strchr("sS", off_color.get_e_col_op()))
        off_color.set_e_sub_sym(strlen(optarg) > 2 ? optarg + 2 : "");
      break;

    case 'F':
      if (col.read(optarg)) {
        off_color.set_f_col(col);
        break;
      }
      parts.init(optarg, ",");
      if (off_color.f_op_check((char *)parts[0], "ksqh"))
        off_color.set_f_col_op(*parts[0]);
      else
        error("invalid coloring", c);

      if (!((strchr("sS", off_color.get_f_col_op()) && parts.size() < 4) ||
            parts.size() < 2))
        error("too many comma separated parts", c);

      if (strchr("sS", off_color.get_f_col_op()))
        off_color.set_f_sub_sym(strlen(optarg) > 2 ? optarg + 2 : "");
      break;

    case 'T': {
      int parts_sz = parts.init(optarg, ",");
      if (parts_sz > 2)
        error("the argument has more than 2 parts", c);

      print_status_or_exit(read_int(parts[0], &num), c);
      if (num < 0 || num > 255)
        error("face transparency must be between 0 and 255", c);

      // if only one part, apply to faces as default
      if (parts_sz == 1) {
        opacity[FACES] = num;
      }
      else if (parts_sz > 1) {
        if (strspn(parts[1], "vefa") != strlen(parts[1]))
          error(msg_str("transparency elements are '%s' must be any or all "
                        "from  v, e, f, a",
                        parts[1]),
                c);

        string str = parts[1];
        if (str.find_first_of("va") != string::npos)
          opacity[VERTS] = num;
        if (str.find_first_of("ea") != string::npos)
          opacity[EDGES] = num;
        if (str.find_first_of("fa") != string::npos)
          opacity[FACES] = num;
      }
      break;
    }

    case 'm':
      // keep map name to pass to stellation function
      map_string = optarg;
      print_status_or_exit(read_colorings(off_color.clrngs, optarg), c);
      break;

    case 'l':
      int sig_compare;
      print_status_or_exit(read_int(optarg, &sig_compare), c);
      if (sig_compare > DEF_SIG_DGTS)
        warning("limit is very small, may not be attainable", c);
      eps = pow(10, -sig_compare);
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
}

struct MillerItem {
  const char *cell_string;
  const char *sub_sym;
  bool remove_inline_verts;
  bool merge_faces;
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
*/

// clang-format off
// Table of miller Models which are not Uniforms, 75 items
// Number, Cells, symmetry, remove inline vertices
// f1 left handed form represented by f3
MillerItem miller_item_list[] = {
   // Full symmetry
   { "A",          "Ih", true,  false }, // 01 w4  Icosahedron
   { "B",          "Ih", true,  false }, // 02 w26 Triakis Icosahedron (60-deltahedron)
   { "C",          "Ih", true,  false }, // 03 w23 Compound of 5 Octahedra
   { "D",          "Ih", true,  false }, // 04
   { "E",          "Ih", true,  false }, // 05
   { "F",          "Ih", true,  false }, // 06 w27
   { "G",          "Ih", true,  false }, // 07 w41 Great Icosahedron
   { "H",          "Ih", true,  false }, // 08 w42 Echidnahedron
   { "e1",         "Ih", true,  false }, // 09 w37
   { "f1",         "Ih", false, false }, // 10 *
   { "g1",         "Ih", true,  false }, // 11 w29
   { "e1f1",       "Ih", false, false }, // 12 *
   { "e1f1g1",     "Ih", true,  false }, // 13
   { "f1g1",       "Ih", false, false }, // 14 *
   { "e2",         "Ih", true,  false }, // 15 *
   { "f2",         "Ih", true,  false }, // 16 Only model which is disconnected
   { "g2",         "Ih", false, true  }, // 17
   { "e2f2",       "Ih", true,  false }, // 18
   { "e2f2g2",     "Ih", true,  false }, // 19
   { "f2g2",       "Ih", true,  false }, // 20 w30
   { "De1",        "Ih", true,  false }, // 21 w32
   { "Ef1",        "Ih", true,  false }, // 22 w25 Compound of 10 Tetrahedra
   { "Fg1",        "Ih", true,  false }, // 23 w31
   { "De1f1",      "Ih", true,  false }, // 24 *
   { "De1f1g1",    "Ih", true,  false }, // 25
   { "Ef1g1",      "Ih", true,  false }, // 26 w28 Excavated Dodecahedron
   { "De2",        "Ih", true,  false }, // 27
   { "Ef2",        "Ih", true,  false }, // 28
   { "Fg2",        "Ih", true,  false }, // 29 w33
   { "De2f2",      "Ih", true,  false }, // 30 w34 Great Triambic Icosahedron
   { "De2f2g2",    "Ih", true,  false }, // 31
   { "Ef2g2",      "Ih", true,  false }, // 32 *
   // Enantiomeric
   { "f1'",        "I",  true,  false }, // 33 w35
   { "e1f1'",      "I",  true,  false }, // 34 w36
   { "De1f1'",     "I",  false, false }, // 35
   { "f1'g1",      "I",  false, false }, // 36
   { "e1f1'g1",    "I",  false, false }, // 37 w39
   { "De1f1'g1",   "I",  false, false }, // 38
   { "f1'g2",      "I",  false, true  }, // 39
   { "e1f1'g2",    "I",  true,  false }, // 40 *
   { "De1f1'g2",   "I",  true,  false }, // 41 *
   { "f1'f2g2",    "I",  true,  false }, // 42
   { "e1f1'f2g2",  "I",  true,  false }, // 43
   { "De1f1'f2g2", "I",  true,  false }, // 44 *
   { "e2f1'",      "I",  false, true  }, // 45 w40
   { "De2f1'",     "I",  true,  false }, // 46
   { "Ef1'",       "I",  true,  false }, // 47 w24 Compound of 5 Tetrahedra
   { "e2f1'g1",    "I",  false, true  }, // 48
   { "De2f1'g1",   "I",  true,  false }, // 49
   { "Ef1'g1",     "I",  true,  false }, // 50
   { "e2f1'f2",    "I",  true,  false }, // 51 w38
   { "De2f1'f2",   "I",  true,  false }, // 52
   { "Ef1'f2",     "I",  true,  false }, // 53
   { "e2f1'f2g1",  "I",  true,  false }, // 54
   { "De2f1'f2g1", "I",  true,  false }, // 55
   { "Ef1'f2g1",   "I",  true,  false }, // 56
   { "e2f1'f2g2",  "I",  true,  false }, // 57
   { "De2f1'f2g2", "I",  true,  false }, // 58
   { "Ef1'f2g2",   "I",  true,  false }, // 59
   // begin the true lost stellations
   { "Be1",        "Ih", true,  false }, // 60
   { "Ce2",        "Ih", true,  false }, // 61 *
   { "Df1",        "I",  true,  false }, // 62 *
   { "Af2",        "Ih", true,  false }, // 63
   { "Df2",        "Ih", true,  false }, // 64 *
   { "De1f1'f2",   "I",  true,  false }, // 65 *
   { "De1g1",      "Ih", true,  false }, // 66 *
   { "Af2g1",      "Ih", true,  false }, // 67
   { "De1f1'f2g1", "I",  true,  false }, // 68 *
   { "Af2g2",      "Ih", true,  false }, // 69
   // begin the candidate lost stellations
   { "Be2",        "Ih", true,  false }, // 70 *
   { "De2f2",      "Ih", true,  false }, // 71
   { "e1g1",       "Ih", true,  false }, // 72 *
   { "Ef1g1",      "Ih", true,  false }, // 73
   { "Cf2g1",      "Ih", true,  false }, // 74
   { "ACDf2g1",    "Ih", true,  false }, // 75 *
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
  cell["A-"]  = "0,18"; // same
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
  // clang-format on

  // code
  vector<string> diagram_list_strings;

  int len = cell_str.length();
  int len_cnt = len;

  // quick check that only valid characters are in the string
  size_t pos = cell_str.find_first_not_of("ABCDEFGHefg12'+-");
  if (pos != string::npos) {
    // found an illegal character is in the string, return empty list
    diagram_list_strings.clear();
    len_cnt = 0; // prevent while loop
  }

  while (len_cnt > 0) {
    string cell_substr = cell_str.substr(len - len_cnt);
    string part;
    int part_len = 0;
    size_t pos = cell_substr.find_first_of("ABCDEFGH");
    if (pos != string::npos) {
      if (pos + 1 <= cell_substr.size())
        part_len = (cell_substr.substr(pos + 1, 1).find_first_of("'+-") !=
                    string::npos)
                       ? 2
                       : 1;
    }
    else {
      pos = cell_substr.find_first_of("efg");
      if (pos != string::npos) {
        if (pos + 2 <= cell_substr.size())
          part_len = (cell_substr.substr(pos + 2, 1).find_first_of("'+-") !=
                      string::npos)
                         ? 3
                         : 2;
      }
    }
    if (part_len) {
      if (pos <= cell_substr.size())
        part = cell_substr.substr(pos, part_len);
      len_cnt -= part_len;
    }

    // fprintf(stderr,"part = '%s'\n", part.c_str());
    // fprintf(stderr,"cell[part] = %s\n", cell[part].c_str());

    // if not found, return empty string list
    if (!cell[part].length()) {
      diagram_list_strings.clear();
      break;
    }
    else {
      // append face number
      diagram_list_strings.push_back(cell[part]);
    }
  }

  return (diagram_list_strings);
}

int Miller::get_poly(Geometry &geom, int sym, string cell_str, string sym_str,
                     map<int, Geometry> &diagrams, miller_opts &opts)
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
  fprintf(stderr, "cell string = %s symmetry = %s\n", cell_str.c_str(),
          sym_str.c_str());

  vector<string> diagram_list_strings = decode_cell_string(cell_str);
  // if list is empty, code string was invalid
  unsigned int sz = diagram_list_strings.size();
  if (!sz)
    return -1;

  vector<vector<int>> idx_lists(sz);

  // data sizes are verified
  for (unsigned int i = 0; i < diagram_list_strings.size(); i++) {
    if (!diagram_list_strings[i].length())
      continue;

    read_idx_list((char *)diagram_list_strings[i].c_str(), idx_lists[i],
                  std::numeric_limits<int>::max(), false);

    // stellation face index is in the first position
    int stellation_face_idx = idx_lists[i][0];

    // construct the diagrams
    if (!diagrams[stellation_face_idx].verts().size())
      diagrams[stellation_face_idx] =
          make_stellation_diagram(geom, stellation_face_idx, sym_str);
  }

  bool merge_faces =
      (sym > -1) ? Miller_items[sym].merge_faces : opts.merge_faces;
  bool remove_inline_verts = Miller_items[sym].remove_inline_verts;
  bool split_pinched = true;
  bool resolve_faces = true;
  bool remove_multiples = true;

  geom = make_stellation(geom, diagrams, idx_lists, sym_str, merge_faces,
                         remove_inline_verts, split_pinched, resolve_faces,
                         remove_multiples, opts.map_string, opts.eps);

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
      Status stat = Coloring(&diagrams[key1.first]).apply_transparency(64);

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

int make_resource_miller(Geometry &geom, string name, miller_opts &opts,
                         string *error_msg = nullptr)
{
  int sym_no = -1;
  // check if it is just the index number, if so format as m%d
  if (read_int(name.c_str(), &sym_no))
    name = "m" + std::to_string(sym_no);
  else if (name.size() < 2 || !strchr("mM", name[0])) {
    if (error_msg)
      *error_msg = "miller number and name designations start with m";
    return -1; // not miller name
  }

  string cell_str;
  string sym_str;

  Miller mill;
  if (read_int(name.c_str() + 1, &sym_no)) {
    sym_no--;
    if (sym_no < 0 || sym_no >= mill.get_last_M()) {
      if (error_msg)
        *error_msg = "miller stellation number out of range";
      return -1; // fail
    }
  }
  // if the next character is an underscore, try it as cell names
  else if (name[1] == '_') {
    cell_str = name.substr(2);
    size_t pos = cell_str.find_first_of(",");
    if (pos != string::npos) {
      sym_str = cell_str.substr(pos + 1);
      if (sym_str.length()) {
        if (sym_str != "I" && sym_str != "Ih") {
          if (error_msg)
            *error_msg = "miller cell name symmetry must be I or Ih";
          return -1; // fail
        }
      }
      cell_str = cell_str.substr(0, pos);
    }
  }
  else {
    if (error_msg)
      *error_msg = "not a string of cell names";
    return -1;
  }

  map<int, Geometry> diagrams;
  int ret = mill.get_poly(geom, sym_no, cell_str, sym_str, diagrams, opts);
  if (ret < 1) {
    if (error_msg)
      *error_msg = "code string was invalid";
    return -1;
  }

  // for items 1 to 75 for seperate compound items
  if (opts.rebuild_compound_model) {
    rebuild_compound(geom);
    geom.add_missing_impl_edges();
  }

  // if coloring method not set
  if (!opts.off_color.get_f_col_op()) {
    // check if it is a compound
    GeometryInfo info(geom);
    if (info.num_parts() > 1)
      opts.off_color.set_f_col_op('k');
    else
      // goldenrod1: Color(255,193,37) ... from George Hart, Color(1.0,0.7,0.1)
      opts.off_color.set_f_col(Color(255, 193, 37));
  }

  color_stellation(geom, opts.off_color);

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

  // apply all element transparencies
  apply_transparencies(geom, opts.opacity);

  return 0; // name found
}

int try_miller(Geometry &geom, miller_opts &opts, string *error_msg = nullptr)
{
  string name = opts.ifile;
  int idx = make_resource_miller(geom, name, opts, error_msg);
  return (idx);
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

  string error_msg;

  Geometry geom;
  if (try_miller(geom, opts, &error_msg))
    if (!error_msg.empty())
      opts.error(error_msg);

  opts.write_or_error(geom, opts.ofile);

  return 0;
}
