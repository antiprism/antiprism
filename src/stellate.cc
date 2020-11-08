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
   Name: stellate.cc
   Description: stellate a polyhedron
   Project: Antiprism - http://www.antiprism.com
*/

#include <cstdio>
#include <cstdlib>

#include <set>
#include <string>
#include <vector>

#include "../base/antiprism.h"

using std::map;
using std::set;
using std::string;
using std::vector;

using namespace anti;

class stellate_opts : public ProgramOpts {
public:
  string ifile;
  string ofile;

  vector<string> diagram_list_strings;
  string sym_str;
  bool merge_faces;
  bool remove_inline_vertices;
  bool split_pinched;
  bool resolve_faces;
  bool remove_multiples;
  bool rebuild_compound_model;

  string output_parts;
  bool move_to_front;
  int projection_width;

  char vertex_coloring_method;
  char edge_coloring_method;
  char face_coloring_method;

  Color vertex_color;
  Color edge_color;
  Color face_color;

  string map_string;
  int face_opacity;

  double epsilon;

  stellate_opts()
      : ProgramOpts("stellate"), merge_faces(true),
        remove_inline_vertices(true), split_pinched(true), resolve_faces(false),
        remove_multiples(false), rebuild_compound_model(false),
        output_parts("s"), move_to_front(false), projection_width(500),
        vertex_coloring_method('\0'), edge_coloring_method('\0'),
        face_coloring_method('d'), vertex_color(Color::invisible),
        edge_color(Color::invisible), face_color(Color()),
        map_string("compound"), face_opacity(-1), epsilon(0)
  {
  }
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
  -O <args> output s - stellation, d - diagram, i - input model (default: s)
               D - diagram faces used highlighted, S - with symmetry
               R - resolved faces used for stellation (when using D or S)
               F - highlighted faces only (when using D, S or R)
  -z        move first diagram to face front (out of symmetry alignment)
  -w <int>  width to project stellation diagram (default: 500)
  -l <lim>  minimum distance for unique vertex locations as negative exponent
               (default: %d giving %.0e)
  -o <file> write output to file (default: write to standard output)

Coloring Options (run 'off_util -H color' for help on color formats)
  -F <opt>  face coloring method. d - from diagram, s - symmetry (default: d)
               c - color by compound
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

void stellate_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  int sig_compare = INT_MAX;

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

void apply_transparency(Geometry &geom, const stellate_opts &opts,
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

void color_stellation(Geometry &stellation, const stellate_opts &opts)
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
      for (int i = 0; i < (int)kis.faces().size(); i++) {
        if (kis.faces(i).size() < 3)
          dels.push_back(i);
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
        int fsz = face.size();
        // face to face
        // connections with invisible faces are ignored
        int connections = 0;
        for (int j = 0; j < fsz; j++) {
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
          int connections = faces.size();
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

// idx_lists still contains stellation face number in position 0
Geometry construct_model(Geometry &geom, map<int, Geometry> &diagrams,
                         vector<vector<int>> &idx_lists,
                         const stellate_opts &opts)
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
        opts.remove_multiples, opts.map_string, opts.epsilon);

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
        apply_transparency(diagrams[key1.first], opts, 64);

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
    color_stellation(stellation, opts);
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

    if (opts.face_opacity > -1)
      apply_transparency(geom, opts);
    model.append(geom);
  }

  // move the first diagram to the front
  if (opts.move_to_front) {
    Geometry diagram = diagrams.begin()->second;
    Vec3d face_normal = face_norm(diagram.verts(), diagram.faces(0));
    if (vdot(face_normal, diagram.face_cent(0)) < 0)
      face_normal *= -1.0;
    Trans3d trans = Trans3d::rotate(face_normal, Vec3d(0, 0, 1));
    model.transform(trans);
  }

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
                      idx_lists[i], INT_MAX, false),
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
                                  opts.projection_width, opts.epsilon);

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
