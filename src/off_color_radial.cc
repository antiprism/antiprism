/*
   Copyright (c) 2017-2021, Roger Kaufman

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
   Name: off_color_radial.cc
   Description: Color in radial pattern based on symmetry
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <set>
#include <string>
#include <vector>

using std::map;
using std::pair;
using std::set;
using std::string;
using std::vector;

using namespace anti;

class radial_opts : public ProgramOpts {
public:
  string ifile;
  string ofile;

  vector<int> axis_orders;      // axis order priority
  vector<int> axis_percent;     // percent length of axes, 100 percent default
  bool axis_orders_set = false; // false if using element lists
  string sym_str;               // for sub-symmetry
  int show_axes = 0;            // if coloring axes

  // element index lists
  vector<pair<char, vector<int>>> elem_lists; // elements instead of axes

  double eps = anti::epsilon;

  char vertex_coloring_method = '\0'; // e
  char edge_coloring_method = '\0';   // f

  Color vertex_color = Color(Color::maximum_index); // unset
  Color edge_color = Color(Color::maximum_index);   // unset

  int coloring_method = 0;   // color radial or axes or both
  int axes_coloring = 1;     // color axes by nfold or order
  string map_string = "rng"; // default map name
  int face_opacity = -1;     // transparency from 0 to 255

  ColorMapMulti map;
  ColorMapMulti map_axes;

  radial_opts() : ProgramOpts("off_color_radial") {}

  void process_command_line(int argc, char **argv);
  void usage();
};

void radial_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] [input_file]

Color in radial pattern based on symmetry.

Options
%s
  -l <lim>  minimum distance change to terminate planarization, as negative
               exponent (default: %d giving %.0e)
  -o <file> write output to file (default: write to standard output)

Scene Options
  -d <opt>  coloring. radial=1, axes=2 (default: 1)
               (multiple -d as needed)
  -a <ax,p> axis order. primary=1, secondary=2, tertiary=3, all=4 (default: 1)
               p is percent length of the axis. (default: 100 percent)
               (multiple -a as needed)
  -f <list> specify a list of elements, list starts with element letter,
            followed by an index number list, given as index ranges separated
            by commas. range can be one number or two numbers separated by a
            hyphen (default range numbers: 0 and largest index).
            Index number list will be preceded by f, e, v for faces, edges and
            vertices. Elements resolve to connected faces to be staring point
            for radial coloring. special selector: s for number of face sides
               (multiple -f as needed, -f overrides -a)
  -s <sym>  symmetry subgroup (Schoenflies notation)

Coloring Options (run 'off_util -H color' for help on color formats)
  -E <col>  edge color (default: unchanged)
               keyword: none - sets no color
               f - color with average adjacent face color
  -V <col>  vertex color (default: unchanged)
               keyword: none - sets no color
               e - color with average adjacent edge color
               f - color with average adjacent face color
  -T <tran> face transparency. valid range from 0 (invisible) to 255 (opaque)
  -m <maps> color maps for all elements to be tried in turn (default: rng)
              rng and rainbow maps without entries specified are calculated 
  -n <maps> color maps for axes (A=1: calculated, A=2: map_red:blue:yellow)
  -A <opt>  color axes by nfold=1, order=2 (default: 1)

)",
          prog_name(), help_ver_text, int(-log(anti::epsilon) / log(10) + 0.5),
          anti::epsilon);
}

void radial_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  string arg_id;

  string map_string_axes;

  handle_long_opts(argc, argv);

  // load axis order 0, 1, and 2 places
  axis_orders.push_back(-1);
  axis_orders.push_back(-1);
  axis_orders.push_back(-1);

  // load percent for axis 0, 1, and 2
  axis_percent.push_back(100);
  axis_percent.push_back(100);
  axis_percent.push_back(100);

  while ((c = getopt(argc, argv, ":hd:f:a:s:E:V:T:m:n:A:l:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'd': {
      print_status_or_exit(
          get_arg_id(optarg, &arg_id, "radial=1|axes=2", argmatch_add_id_maps),
          c);
      int option = atoi(arg_id.c_str());
      if (option == 1)
        coloring_method = 1;
      else if (option == 2)
        show_axes = 2;
      break;
    }

    case 'f': {
      pair<char, vector<int>> elem_lst;
      string elements = "vefs";
      if (elements.find(optarg[0]) == string::npos)
        error("element type must be v, e, f, or s");
      else {
        elem_lst.first = optarg[0];
        optarg++;
      }
      vector<int> idx_lst;
      print_status_or_exit(read_idx_list(optarg, idx_lst,
                                         std::numeric_limits<int>::max(),
                                         false),
                           c);
      if (!idx_lst.size())
        error("no element numbers are in the list", c);
      elem_lst.second = idx_lst;
      elem_lists.push_back(elem_lst);
      break;
    }

    case 'a': {
      Split parts(optarg, ",");
      unsigned int parts_sz = parts.size();

      if (parts_sz > 2)
        error("excess entries for axis orders", c);

      int axis_order;
      for (unsigned int i = 0; i < parts_sz; i++) {
        if (i == 0) {
          print_status_or_exit(
              get_arg_id(parts[i], &arg_id,
                         "primary=1|secondary=2|tertiary=3|all=4",
                         argmatch_add_id_maps),
              c);
          axis_order = atoi(arg_id.c_str());
          if (axis_order == 4) {
            for (int j = 0; j < 3; j++)
              axis_orders[j] = j;
          }
          else {
            // axes start with 0
            axis_order--;
            axis_orders[axis_order] = axis_order;
          }
        }
        else if (i == 1) {
          int ax_pct;
          read_int(parts[i], &ax_pct);
          if (ax_pct < 1)
            error("axis percent must be 1 or greater", c);
          if (axis_order == 4) {
            for (int j = 0; j < 3; j++)
              axis_percent[j] = ax_pct;
          }
          else {
            axis_percent[axis_order] = ax_pct;
          }
        }
      }

      axis_orders_set = true;
      break;
    }

    case 's':
      sym_str = optarg;
      break;

    case 'V':
      if (strlen(optarg) == 1 && strchr("ef", int(*optarg)))
        vertex_coloring_method = *optarg;
      else
        print_status_or_exit(vertex_color.read(optarg), c);
      break;

    case 'E':
      if (strlen(optarg) == 1 && strchr("f", int(*optarg)))
        edge_coloring_method = *optarg;
      else
        print_status_or_exit(edge_color.read(optarg), c);
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

    case 'n':
      map_string_axes = optarg;
      break;

    case 'A':
      print_status_or_exit(
          get_arg_id(optarg, &arg_id, "nfold=1|order=2", argmatch_add_id_maps),
          c);
      axes_coloring = atoi(arg_id.c_str());
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

  // default coloring method
  if (show_axes == 0)
    coloring_method = 1;

  // if elements specified, don't use axes
  if (elem_lists.size()) {
    for (unsigned int i = 0; i < axis_orders.size(); i++)
      axis_orders[i] = -1;
    axis_orders_set = false;

    if (sym_str.length()) {
      sym_str = "";
      warning("symmetry specification has no effect when choosing faces", 's');
    }
  }
  // default axis order is 0
  else if (!axis_orders_set) {
    axis_orders[0] = 0;
    axis_orders_set = true;
  }

  print_status_or_exit(map.init(map_string.c_str()), 'm');

  // default axes map
  if (map_string_axes.size())
    print_status_or_exit(map_axes.init(map_string_axes.c_str()), 'n');
  else {
    if (axes_coloring == 1) {
      // nfold map is same as antiview
      auto *col_map0 = new ColorMapMap;
      col_map0->set_col(0, Color(0.6, 0.3, 0.0));
      col_map0->set_col(1, Color());
      col_map0->set_col(2, Color(0.8, 0.8, 0.2));
      col_map0->set_col(3, Color(0.3, 0.8, 0.3));
      col_map0->set_col(4, Color(0.6, 0.0, 0.0));
      col_map0->set_col(5, Color(0.0, 0.0, 0.6));
      map_axes.add_cmap(col_map0);

      // nfold 6 and higher same color
      auto *col_map = new ColorMapMap;
      col_map->set_col(0, Color(0.6, 0.3, 0.0));
      col_map->set_wrap();
      map_axes.add_cmap(col_map);
    }
    else if (axes_coloring == 2) {
      // match symmetro face colors
      map_string_axes = "map_red:blue:yellow";
      print_status_or_exit(map_axes.init(map_string_axes.c_str()), 'n');
    }
  }

  if (argc - optind > 1)
    error("too many arguments");

  if (argc - optind == 1)
    ifile = argv[optind];
}

// code to get axes by Adrian Rossiter
// sym.init() done before call
// color with map indexes by axis number
Geometry get_axes(const Geometry &geom, Symmetry &sym, int ax_idx,
                  const radial_opts &opts)
{
  if (opts.sym_str.length()) {
    Symmetry full_sym = sym;
    full_sym.get_sub_sym(opts.sym_str, &sym);
  }

  vector<int> n_folds;
  switch (sym.get_sym_type()) {
  case Symmetry::I:
  case Symmetry::Ih:
    n_folds = {5, 1, 3, 1, 2, 1};
    break;
  case Symmetry::O:
  case Symmetry::Oh:
    n_folds = {4, 1, 3, 1, 2, 1};
    break;
  case Symmetry::T:
  case Symmetry::Th:
  case Symmetry::Td:
    n_folds = {3, 1, 3, 1, 2, 1};
    break;
  case Symmetry::D:
  case Symmetry::Dv:
  case Symmetry::Dh:
    n_folds = {sym.get_nfold(), 1, 2, 1, 2, 1};
    break;
  case Symmetry::S:
  case Symmetry::C:
  case Symmetry::Cv:
  case Symmetry::Ch:
    n_folds = {sym.get_nfold()}; // Use principal axis vertex
    break;
  }

  Geometry axes;

  if (n_folds.size() && ((n_folds.size() == 6) || (ax_idx == 0))) {
    double max_dist = 0.0;
    for (const auto &vert : geom.verts())
      // Use sym.get_to_std() to centre on orig
      max_dist = std::max(max_dist, (sym.get_to_std() * vert).len());

    Vec3d axis_pt;             // A point on an axis
    if (n_folds.size() == 6) { // Three distinct axes
      vector<Vec3d> tri_verts;
      get_schwarz_tri_verts(n_folds, tri_verts);
      axis_pt = tri_verts[ax_idx];
    }
    else // Single axis
      axis_pt = Vec3d::Z;

    Geometry axis;

    axis.add_vert(
        sym.get_to_std().inverse() * // Carry from std pos'n to geom pos'n
        (max_dist * 1.1 * axis_pt)); // Scale the axis point first
    sym_repeat(axes, axis, sym);     // Repeat positioned axis point
    merge_coincident_elements(axes, "v", 1e-8);

    // set color map index to axis vertices
    int map_idx = 0;
    if (opts.axes_coloring == 1) {
      // reverse secondary and tertiary (why?)
      int ax_tweaked = (ax_idx == 1) ? 2 : ((ax_idx == 2) ? 1 : 0);
      map_idx = n_folds[ax_tweaked % n_folds.size()];
      if ((ax_tweaked == 0) && (sym.get_sym_type() == Symmetry::S))
        map_idx /= 2;
      // sometimes nfold is 1, needs to be 2
      if (map_idx == 1)
        map_idx++;
      // fprintf(stderr,"ax_tweaked = %d map_idx = %d\n", ax_tweaked, map_idx);
    }
    else if (opts.axes_coloring == 2) {
      map_idx = ax_idx;
    }
    for (unsigned int i = 0; i < axes.verts().size(); i++)
      axes.colors(VERTS).set(i, map_idx);
  }

  return axes;
}

// find starting faces for fronts
void find_fronts_from_axes(map<int, vector<int>> &fronts, const Geometry &axes,
                           const Geometry &geom, const Vec3d &cent,
                           const radial_opts &opts)
{
  const vector<Vec3d> &axes_verts = axes.verts();

  const vector<vector<int>> &faces = geom.faces();
  const vector<Vec3d> &verts = geom.verts();

  // fronts radiate from axes
  bool found = false;
  for (unsigned int i = 0; i < axes_verts.size(); i++) {
    // search vertices first
    for (unsigned int j = 0; j < verts.size(); j++) {
      // check if vertex is on axis because face might not be planar
      if (point_in_segment(verts[j], axes_verts[i], cent, opts.eps).is_set()) {
        // intersection is a vertex. find all faces at vertex
        fronts[i] = find_faces_with_vertex(faces, j);
        break;
      }
    }
    // if not found, check edges
    if (!fronts[i].size()) {
      // edges need to exist
      vector<vector<int>> edges;
      geom.get_impl_edges(edges);
      for (unsigned int j = 0; j < edges.size(); j++) {
        int v1 = edges[j][0];
        int v2 = edges[j][1];
        // if (lines_intersection(axes_verts[i], cent, verts[v1], verts[v2],
        // opts.eps).is_set()) {
        if (segments_intersection(axes_verts[i], cent, verts[v1], verts[v2],
                                  opts.eps)
                .is_set()) {
          // intersection is an edge
          fronts[i] = find_faces_with_edge(faces, make_edge(v1, v2));
          break;
        }
      }
    }
    // if not found, find if axis passes through a face centroid
    if (!fronts[i].size()) {
      for (unsigned int j = 0; j < faces.size(); j++) {
        // for (unsigned int j = 0; j < 1; j++) {
        Vec3d face_normal = geom.face_norm(j).unit();
        Vec3d face_centroid = geom.face_cent(j);
        // make sure face_normal points outward
        if (vdot(face_normal, face_centroid) < 0)
          face_normal *= -1.0;

        int where = 0;
        int *w = &where;
        Vec3d P = line_plane_intersect(face_centroid, face_normal,
                                       axes_verts[i], cent, w, opts.eps);

        // if point is on axis
        if (point_in_segment(P, axes_verts[i], cent, opts.eps).is_set()) {
          // get winding number, if not zero, point is on a polygon
          vector<int> face_idxs;
          face_idxs.push_back(j);
          Geometry polygon = faces_to_geom(geom, face_idxs);
          Normal normal(polygon, face_normal, 0, cent, opts.eps);
          vector<Vec3d> one_point;
          one_point.push_back(P);
          int winding_number = get_winding_number_polygon(
              polygon, one_point, normal, true, opts.eps);
          if (winding_number) {
            fronts[i].push_back(j);
            break;
          }
        }
      }
    }

    if (fronts[i].size())
      found = true;
    else
      opts.warning(msg_str("no element found at axis %d", i));
  }

  if (!found) {
    fronts[0].push_back(0);
    opts.warning("no axis intersections found. using face 0");
  }
}

// if index list is supplied, use it for fronts
void find_fronts_from_list(map<int, vector<int>> &fronts, const Geometry &geom,
                           const radial_opts &opts)
{
  // if index list is supplied, use it for fronts
  // edges need to exist
  vector<vector<int>> edges;
  geom.get_impl_edges(edges);
  int front_no = 0;
  for (unsigned int i = 0; i < opts.elem_lists.size(); i++) {
    pair<char, vector<int>> elem_lst = opts.elem_lists[i];
    char elem_type = elem_lst.first;
    vector<int> idx_list = elem_lst.second;
    for (unsigned int j = 0; j < idx_list.size(); j++) {
      vector<int> face_idx;
      if (elem_type == 'v') {
        face_idx = find_faces_with_vertex(geom.faces(), idx_list[j]);
        fronts[front_no].insert(fronts[front_no].end(), face_idx.begin(),
                                face_idx.end());
      }
      else if (elem_type == 'e') {
        int v1 = edges[idx_list[j]][0];
        int v2 = edges[idx_list[j]][1];
        face_idx = find_faces_with_vertex(geom.faces(), v1);
        fronts[front_no].insert(fronts[front_no].end(), face_idx.begin(),
                                face_idx.end());
        face_idx = find_faces_with_vertex(geom.faces(), v2);
        fronts[front_no].insert(fronts[front_no].end(), face_idx.begin(),
                                face_idx.end());
      }
      else if (elem_type == 'f') {
        fronts[front_no].push_back(idx_list[j]);
      }
      else if (elem_type == 's') {
        bool found = false;
        for (unsigned int k = 0; k < geom.faces().size(); k++) {
          if (idx_list[j] == (int)geom.faces(k).size()) {
            fronts[front_no].push_back(k);
            front_no++;
            found = true;
          }
        }
        // front number will be advance once too far
        if (found)
          front_no--;
      }
      front_no++;
    }
  }
}

// color with map indexes
void radial_coloring(Geometry &geom, map<int, vector<int>> &fronts)
{
  // clear all colors
  geom.colors(FACES).clear();

  int ridge = 0;
  bool found = true;

  while (found) {
    found = false;
    // starting from each radial point
    for (auto const &key1 : fronts) {
      vector<int> front = fronts[key1.first];
      vector<int> next_ridge;
      vector<int> adjacent_faces;

      for (unsigned int i = 0; i < front.size(); i++) {
        // color faces in radial ridge with indexes
        // it can be set in one final ridge but still listed in another so check
        if (!(geom.colors(FACES).get(front[i])).is_set())
          geom.colors(FACES).set(front[i], ridge);

        // get faces connected to this face via edge or vertex
        vector<int> face_idx;
        vector<int> face = geom.faces()[front[i]];
        unsigned int fsz = face.size();
        for (unsigned int j = 0; j < fsz; j++) {
          int v1 = face[j];
          int v2 = face[(j + 1) % fsz];
          // face_idx = find_faces_with_edge(geom.faces(), make_edge(v1, v2));
          // adjacent_faces.insert(adjacent_faces.end(), face_idx.begin(),
          //                      face_idx.end());
          face_idx = find_faces_with_vertex(geom.faces(), v1);
          adjacent_faces.insert(adjacent_faces.end(), face_idx.begin(),
                                face_idx.end());
          face_idx = find_faces_with_vertex(geom.faces(), v2);
          adjacent_faces.insert(adjacent_faces.end(), face_idx.begin(),
                                face_idx.end());
        }
      }

      // make list unique
      sort(adjacent_faces.begin(), adjacent_faces.end());
      auto vi = unique(adjacent_faces.begin(), adjacent_faces.end());
      adjacent_faces.resize(vi - adjacent_faces.begin());

      // find next ridge of faces not yet colored
      for (unsigned int i = 0; i < adjacent_faces.size(); i++) {
        if (!(geom.colors(FACES).get(adjacent_faces[i])).is_set()) {
          next_ridge.push_back(adjacent_faces[i]);
        }
      }
      fronts[key1.first] = next_ridge;

      if (next_ridge.size())
        found = true;
    }
    ridge++;
  }
}

// break down model into parts
void compound_parts(Geometry &geom, vector<Geometry> &parts)
{
  GeometryInfo info(geom);
  if (info.num_parts() == 1)
    parts.push_back(geom);
  else {
    vector<vector<int>> face_parts;
    geom.orient(&face_parts);

    for (unsigned int i = 0; i < face_parts.size(); i++) {
      Geometry part = faces_to_geom(geom, face_parts[i]);
      part.del(VERTS, part.get_info().get_free_verts());
      part.orient(1);
      parts.push_back(part);
    }
  }
}

void set_indexes_to_color(Geometry &geom, const radial_opts &opts)
{
  // find maximum index
  int max_ridge = 0;
  for (unsigned int i = 0; i < geom.faces().size(); i++) {
    int idx = geom.colors(FACES).get(i).get_index();
    if (idx > max_ridge)
      max_ridge = idx;
  }
  max_ridge++;

  // default is to make a rainbow map of number of ridges
  opts.message(msg_str("maximum ridges formed is %d", max_ridge));
  string map_name;
  if (opts.map_string == "rng" || opts.map_string == "rainbow") {
    map_name = opts.map_string + std::to_string(max_ridge);
    opts.message(msg_str("default map used is %s", map_name.c_str()),"option -m");
  }

  Coloring clrng;
  clrng.set_geom(&geom);
  if (map_name.length()) {
    ColorMap *cmap = colormap_from_name(map_name.c_str());
    clrng.add_cmap(cmap);
  }
  else
    clrng.add_cmap(opts.map.clone());
  clrng.f_apply_cmap();

  // apply transparency
  Status stat = Coloring(&geom).apply_transparency(opts.face_opacity);
  if (stat.is_warning())
    opts.warning(stat.msg(), 'T');
}

// for now show_axes is always 2 to display whole axis
void append_axes(Geometry &geom, Geometry &axes, const radial_opts &opts)
{
  if (opts.show_axes == 2) {
    // add axes
    const vector<Vec3d> &verts = axes.verts();
    unsigned int vsz = verts.size();

    // add a vertex to model at centroid
    Vec3d cent = centroid(geom.verts());
    axes.add_vert(cent, Color(255, 255, 255));
    int cent_vert = vsz;

    // connect centroid to axes points
    Color c;
    for (unsigned int i = 0; i < vsz; i++) {
      c = axes.colors(VERTS).get(i);
      // int v_mirror = find_vert_by_coords(axes, -verts[i], opts.eps);
      axes.add_edge(make_edge(i, cent_vert), c);
    }

    // for when there is only the base axis
    if (vsz == 1) {
      axes.transform(Trans3d::translate(-cent));
      axes.add_vert((-verts[0]), c);
      axes.transform(Trans3d::translate(cent));
      axes.add_edge(make_edge((int)verts.size() - 1, cent_vert), c);
    }
  }

  Coloring clrng;
  clrng.set_geom(&axes);
  clrng.add_cmap(opts.map_axes.clone());
  clrng.v_apply_cmap();
  clrng.e_apply_cmap();

  geom.append(axes);
}

void radial_coloring(Geometry &geom, const Geometry &axes,
                     const radial_opts &opts)
{
  // centroid needs to be preserved before possible compound break down
  Vec3d cent = centroid(geom.verts());

  // find fronts from axes or list
  map<int, vector<int>> fronts;
  // if fronts from elements, do this once
  if (opts.elem_lists.size())
    find_fronts_from_list(fronts, geom, opts);

  // process possible compound
  vector<Geometry> parts;
  if (opts.axis_orders_set)
    compound_parts(geom, parts);
  else
    // when front from elements, compound cannot be split
    parts.push_back(geom);

  // process each part, then rebuild geom
  geom.clear_all();
  for (unsigned int i = 0; i < parts.size(); i++) {
    if (opts.axis_orders_set) {
      fronts.clear();
      // fronts from axes needs to be done for each part
      find_fronts_from_axes(fronts, axes, parts[i], cent, opts);
    }
    radial_coloring(parts[i], fronts);
    geom.append(parts[i]);
  }

  // now turn all indexes to color
  set_indexes_to_color(geom, opts);
}

void ev_coloring(Geometry &geom, const radial_opts &opts)
{
  geom.add_missing_impl_edges();
  if (opts.edge_coloring_method == 'f') {
    // edges take colors from faces
    Coloring clrng(&geom);
    clrng.e_from_adjacent(FACES);
  }
  // if color hasn't been input let current color unchanged
  else if (!opts.edge_color.is_maximum_index())
    // use color selected
    Coloring(&geom).e_one_col(opts.edge_color);

  if (opts.vertex_coloring_method == 'f') {
    // edges take colors from edges
    Coloring clrng(&geom);
    clrng.v_from_adjacent(FACES);
  }
  else if (opts.vertex_coloring_method == 'e') {
    // edges take colors from edges
    Coloring clrng(&geom);
    clrng.v_from_adjacent(EDGES);
  }
  // if color hasn't been input let current color unchanged
  else if (!opts.vertex_color.is_maximum_index())
    // use color selected
    Coloring(&geom).v_one_col(opts.vertex_color);
}

int main(int argc, char *argv[])
{
  radial_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);

  if (opts.elem_lists.size()) {
    unsigned int vsz = geom.verts().size();
    vector<vector<int>> edges;
    geom.get_impl_edges(edges);
    unsigned int esz = edges.size();
    unsigned int fsz = geom.faces().size();
    for (unsigned int i = 0; i < opts.elem_lists.size(); i++) {
      pair<char, vector<int>> elem_lst = opts.elem_lists[i];
      char elem_type = elem_lst.first;
      vector<int> idx_list = elem_lst.second;
      if (elem_type == 's') {
        bool found = false;
        for (unsigned int j = 0; j < idx_list.size(); j++) {
          for (unsigned int k = 0; k < fsz; k++) {
            if (idx_list[j] == (int)geom.faces(k).size()) {
              found = true;
              break;
            }
          }
        }
        if (!found)
          opts.error("no faces found of size requested");
      }
      else {
        for (unsigned int j = 0; j < idx_list.size(); j++) {
          unsigned int sz = 0;
          string elem_name;
          if (elem_type == 'v') {
            sz = vsz;
            elem_name = "vertex";
          }
          else if (elem_type == 'e') {
            sz = esz;
            elem_name = "edge";
          }
          else if (elem_type == 'f') {
            sz = fsz;
            elem_name = "face";
          }
          if (idx_list[j] > (int)sz - 1) {
            opts.error(
                msg_str("term %d, element %d, %s number %d is larger than %s "
                        "numbers of input model: %d",
                        i + 1, j + 1, elem_name.c_str(), idx_list[i],
                        elem_name.c_str(), sz - 1),
                'f');
          }
        }
      }
    }
  }

  Symmetry sym;
  // only init if needed
  if (opts.axis_orders_set)
    sym.init(geom);

  // process symmetry symbol
  if (opts.sym_str.length()) {
    // if specified, check validity
    Symmetry full_sym = sym;
    Status stat = full_sym.get_sub_sym(opts.sym_str, &sym);
    if (stat.is_error())
      opts.error(msg_str("invalid subsymmetry '%s': %s", opts.sym_str.c_str(),
                         stat.c_msg()),
                 's');
  }

  // collect axes
  Geometry axes;
  if (opts.axis_orders_set) {
    for (unsigned int i = 0; i < opts.axis_orders.size(); i++) {
      if (opts.axis_orders[i] != -1) {
        Geometry ax = get_axes(geom, sym, opts.axis_orders[i], opts);
        if (!ax.verts().size()) {
          opts.warning(msg_str("no axes of order %d", i));
          continue;
        }
        if (opts.axis_percent[i] != 100) {
          double pct = (double)opts.axis_percent[i] / 100;
          ax.transform(Trans3d::scale(pct));
        }
        axes.append(ax);
      }
    }
    // fprintf(stderr, "axes.size = %d\n", (int)axes.verts().size());
  }

  if (opts.coloring_method == 1) {
    radial_coloring(geom, axes, opts);
    ev_coloring(geom, opts);
  }

  // append axes if set
  if (opts.show_axes)
    append_axes(geom, axes, opts);

  opts.write_or_error(geom, opts.ofile);

  return 0;
}
