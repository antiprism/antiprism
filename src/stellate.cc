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
   Name: stellate.cc
   Description: stellate a polyhedron
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"
#include "color_common.h"

#include <cstdio>
#include <string>
#include <vector>

using std::map;
using std::string;
using std::vector;

using namespace anti;

class stellate_opts : public ProgramOpts {
public:
  string ifile;
  string ofile;

  string output_parts = "s";           // s, d, i, D, S, F, R
  bool merge_faces = false;            // do a planar merge on facelet
  bool remove_inline_vertices = true;  // if inline vertices exist, remove them
  bool split_pinched = true;           // split pinched bowties
  bool resolve_faces = false;          // resolves faces of same color
  bool remove_multiples = false;       // remove multiple resolved faces
  bool rebuild_compound_model = false; // rebuild with seperate constituents

  bool move_to_front = false; // move side with stellation to front
  int projection_width = 500; // magnification of diagram

  vector<string> diagram_list_strings; // face numbers of diagram to use
  string sym_str;                      // for sub-symmetry

  double eps = anti::epsilon;

  // map name is used to pass to stellation function
  string map_string = "compound";

  // internal colorings are needed for local coloring
  OffColor off_color = OffColor(map_string);

  int opacity[3] = {-1, -1, -1}; // transparency from 0 to 255, for v,e,f

  stellate_opts() : ProgramOpts("stellate") {}

  void process_command_line(int argc, char **argv);
  void usage();
};

void stellate_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] [input_file]

Stellate a polyhedron.

Options
%s
  -l <lim>  minimum distance for unique vertex locations as negative exponent
              (default: %d giving %.0e)
  -o <file> write output to file (default: write to standard output)

Program Options
  -f <fnos> face number of input file for stellation diagram (default: 0)
            followed by face numbers of stellation diagram for stellation
            model separated by commas. multiple -n parameters as needed
  -s <sym>  symmetry subgroup (Schoenflies notation)
  -M        do not merge stellation facelets
  -I        do not remove inline vertices (if not -M)
  -S        do not split pinched faces (if not -M)
  -R        resolve stellation facelets
  -D        remove multiples occurrences (sets -R)
  -r        rebuild compound model to separate vertices

Scene Options
  -O <args> output s - stellation, d - diagram, i - input model (default: s)
              (one of: D - diagram faces used highlighted, R - resolved faces
               used for stellation, S - expanded for symmetry)
              F - highlight used faces only (when using D, S or R)
  -z        move first diagram to face front (out of symmetry alignment)
  -w <int>  width to project stellation diagram (default: 500)
  
Coloring Options (run 'off_util -H color' for help on color formats)
keyword: none - sets no color
  -F <col>  color the faces according to: (default: q)
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

void stellate_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;
  int num;

  Split parts;
  Color col;

  off_color.set_e_col(Color::invisible);
  off_color.set_v_col(Color::invisible);

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hf:s:MISRDrzw:O:V:E:F:T:m:l:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'f':
      // have to check face list range later
      diagram_list_strings.push_back(optarg);
      break;

    case 's':
      sym_str = optarg;
      break;

    case 'M':
      merge_faces = false;
      break;

    case 'I':
      remove_inline_vertices = false;
      break;

    case 'S':
      split_pinched = false;
      break;

    case 'R':
      resolve_faces = true;
      break;

    case 'D':
      remove_multiples = true;
      resolve_faces = true;
      break;

    case 'r':
      rebuild_compound_model = true;
      break;

    case 'z':
      move_to_front = true;
      break;

    case 'w':
      print_status_or_exit(read_int(optarg, &projection_width), c);
      if (projection_width <= 0)
        error("projection width must be greater than zero", c);
      break;

    case 'O':
      if (strspn(optarg, "sdiDSFR") != strlen(optarg))
        error(msg_str("output parts are '%s' must be any or all from "
                      "s, d, i, D, S, F, R",
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

// idx_lists still contains stellation face number in position 0
Geometry construct_model(Geometry &geom, map<int, Geometry> &diagrams,
                         vector<vector<int>> &idx_lists, stellate_opts &opts)
{
  Geometry model;
  Geometry stellation;

  // construct the stellation, before coloring diagrams...
  // make_stellation needs diagrams to have color indexes
  string geom_sym_symbol = opts.sym_str;
  if (opts.output_parts.find("s") != string::npos) {
    // set the symmetry
    if (!geom_sym_symbol.length()) {
      Symmetry geom_full_sym(geom);
      geom_sym_symbol = geom_full_sym.get_symbol().c_str();
    }

    stellation = make_stellation(
        geom, diagrams, idx_lists, geom_sym_symbol, opts.merge_faces,
        opts.remove_inline_vertices, opts.split_pinched, opts.resolve_faces,
        opts.remove_multiples, opts.map_string, opts.eps);

    if (opts.rebuild_compound_model) {
      rebuild_compound(stellation);
      stellation.add_missing_impl_edges();
    }
  }

  // need to execute before coloring diagrams
  // these functions need diagrams to have color indexes
  vector<vector<int>> idx_lists_full =
      lists_full(diagrams, idx_lists, opts.remove_multiples);
  vector<vector<int>> idx_lists_resolved =
      lists_resolved(geom, geom_sym_symbol, diagrams, idx_lists, idx_lists_full,
                     opts.remove_multiples);

  // color diagrams
  for (auto const &key1 : diagrams) {
    Coloring clrng(&diagrams[key1.first]);
    ColorMap *cmap = colormap_from_name(opts.map_string.c_str());
    clrng.add_cmap(cmap);
    clrng.f_apply_cmap();
  }

  // diagrams
  if (opts.output_parts.find_first_of("dDSRF") != string::npos) {
    if (opts.output_parts.find_first_of("DSR") != string::npos) {
      // set large diagram transparency to 1/4th
      for (auto const &key1 : diagrams)
        Status stat = Coloring(&diagrams[key1.first]).apply_transparency(64);

      vector<vector<int>> *lists = &idx_lists;
      if (opts.output_parts.find_first_of("D") != string::npos)
        lists = &idx_lists;
      else if (opts.output_parts.find_first_of("S") != string::npos)
        lists = &idx_lists_full;
      else if (opts.output_parts.find_first_of("R") != string::npos)
        lists = &idx_lists_resolved;

      for (unsigned int i = 0; i < lists->size(); i++) {
        // stellation face index is in the first position
        int stellation_face_idx = (*lists)[i][0];

        // set selected faces back to full visibility
        for (unsigned int j = 1; j < (*lists)[i].size(); j++) {
          int face_no = (*lists)[i][j];
          Color c = diagrams[stellation_face_idx].colors(FACES).get(face_no);
          c.set_alpha(255);
          diagrams[stellation_face_idx].colors(FACES).set(face_no, c);
        }
      }

      if (opts.output_parts.find_first_of("F") != string::npos) {
        for (auto const &key1 : diagrams) {
          vector<int> del_faces;
          for (unsigned int i = 0; i < diagrams[key1.first].faces().size();
               i++) {
            Color c = diagrams[key1.first].colors(FACES).get(i);
            // if not fully opaque, delete faces
            if (c.get_transparency())
              del_faces.push_back(i);
          }

          if (del_faces.size()) {
            diagrams[key1.first].del(FACES, del_faces);
            // delete free edges
            vector<vector<int>> implicit_edges;
            diagrams[key1.first].get_impl_edges(implicit_edges);
            const vector<vector<int>> &edges = diagrams[key1.first].edges();
            vector<int> del_edges;
            for (unsigned int i = 0; i < edges.size(); i++) {
              if (find_edge_in_edge_list(implicit_edges, edges[i]) == -1)
                del_edges.push_back(i);
            }
            diagrams[key1.first].del(EDGES, del_edges);
            diagrams[key1.first].del(
                VERTS, diagrams[key1.first].get_info().get_free_verts());
          }
        }
      }
    }

    for (auto const &key1 : diagrams)
      model.append(diagrams[key1.first]);
  }

  // stellation model
  if (opts.output_parts.find("s") != string::npos) {
    color_stellation(stellation, opts.off_color);
    model.append(stellation);
  }

  // input model
  if (opts.output_parts.find("i") != string::npos) {
    // color faces by symmetry
    Symmetry geom_full_sym(geom);
    Symmetry geom_sub_sym = geom_full_sym;

    vector<vector<set<int>>> sym_equivs;
    geom_sub_sym.init(geom, &sym_equivs);

    // color by sub symmetry
    geom_full_sym.get_sub_sym(opts.sym_str, &geom_sub_sym);
    get_equiv_elems(geom, geom_sub_sym.get_trans(), &sym_equivs);

    Coloring clrng(&geom);
    ColorMap *cmap = colormap_from_name(opts.map_string.c_str());
    clrng.add_cmap(cmap);
    clrng.f_sets(sym_equivs[2], true);

    // apply transparency
    Status stat = Coloring(&geom).apply_transparency(opts.opacity[FACES]);
    if (stat.is_warning())
      opts.warning(stat.msg(), 'T');

    model.append(geom);
  }

  // move the first diagram to the front
  if (opts.move_to_front) {
    Geometry diagram = diagrams[idx_lists[0][0]];
    Vec3d face_normal = face_norm(diagram.verts(), diagram.faces(0));
    if (vdot(face_normal, diagram.face_cent(0)) < 0)
      face_normal *= -1.0;
    Trans3d trans = Trans3d::rotate(face_normal, Vec3d(0, 0, 1));
    model.transform(trans);
  }

  // apply all element transparencies
  apply_transparencies(model, opts.opacity);

  return model;
}

int main(int argc, char *argv[])
{
  stellate_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);

  // process symmetry symbol
  Symmetry sym(geom);
  if (!opts.sym_str.length())
    opts.sym_str = sym.get_symbol().c_str();
  else {
    // if specified, check validity
    sym.init(geom);
    Symmetry full_sym = sym;
    Status stat = full_sym.get_sub_sym(opts.sym_str, &sym);
    if (stat.is_error())
      opts.error(msg_str("invalid subsymmetry '%s': %s", opts.sym_str.c_str(),
                         stat.c_msg()),
                 's');
  }

  int fsz = geom.faces().size() - 1;

  map<int, Geometry> diagrams;

  int sz = opts.diagram_list_strings.size();
  vector<vector<int>> idx_lists(sz);

  bool display_diagrams_only = true;
  for (int i = 0; i < sz; i++) {
    opts.print_status_or_exit(
        read_idx_list((char *)opts.diagram_list_strings[i].c_str(),
                      idx_lists[i], std::numeric_limits<int>::max(), false),
        'f');
    if (!idx_lists[i].size())
      opts.error("no face number are input", 'f');

    // stellation face index is in the first position
    int stellation_face_idx = idx_lists[i][0];
    if (stellation_face_idx > fsz)
      opts.error(msg_str("stellation face(%d) number given: %d is larger than "
                         "maximum model face number: %d",
                         i + 1, stellation_face_idx, fsz),
                 'f');

    // will hold just 1 if it is just the stellation face
    if (idx_lists[i].size() > 1)
      display_diagrams_only = false;

    // construct the diagrams
    if (!diagrams[stellation_face_idx].verts().size())
      diagrams[stellation_face_idx] =
          make_stellation_diagram(geom, stellation_face_idx, opts.sym_str,
                                  opts.projection_width, opts.eps);

    // check face index range. start from 1 since 0 is a placeholder for
    // stellation face
    int max_idx = -1;
    for (unsigned int j = 1; j < idx_lists[i].size(); j++) {
      if (idx_lists[i][j] > max_idx)
        max_idx = idx_lists[i][j];
    }

    int dsz = (int)diagrams[stellation_face_idx].faces().size() - 1;
    if (max_idx > dsz)
      opts.error(msg_str("diagram(%d) number given: %d is larger than maximum "
                         "face number: %d",
                         i + 1, max_idx, dsz),
                 'f');
  }

  // if not enough information for the stellation, just show the diagrams
  if (display_diagrams_only && opts.output_parts != "i") {
    size_t found = opts.output_parts.find_first_of("s");
    if (found != string::npos)
      opts.output_parts.erase(found);
    opts.output_parts = "d" + opts.output_parts;
  }

  geom = construct_model(geom, diagrams, idx_lists, opts);
  // color_by_edge_usage(geom);

  opts.write_or_error(geom, opts.ofile);

  return 0;
}
