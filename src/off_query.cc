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
   Name: off_query.cc
   Description: analyse an off_file and print information about elements
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"
#include "rep_print.h"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <set>
#include <string>
#include <utility>
#include <vector>

using std::pair;
using std::set;
using std::string;
using std::vector;

using namespace anti;

class oq_opts : public ProgramOpts {
public:
  Geometry geom;
  Vec3d center;
  bool center_is_centroid;
  vector<Vec3d> extra_verts;
  vector<vector<int>> extra_faces;
  vector<vector<int>> extra_edges;
  vector<int> idxs;
  int sig_digits;
  string query;
  bool orient;
  char edge_type;
  string ifile;
  string ofile;

  oq_opts()
      : ProgramOpts("off_query"), center(Vec3d(0, 0, 0)),
        center_is_centroid(false), sig_digits(17), orient(true), edge_type('a')
  {
  }

  void process_command_line(int argc, char **argv);
  void post_process_command_line(int argc, char **argv, Geometry &geom);
  void usage();
};

void oq_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] query [input_file]

Read a file in OFF format and list element data for specified elements.
Added elements of the query type are also added to the query list.
Added elements are referred to by xN, where N is the order the element
was added (starting with 0). Query is a list of letters, the first is the
query element type V, F, or E followed by the values to print
   V - Vertex:
       c - coordinates             d - distance from centre
       f - face index numbs        F - face angles
       a - solid angle             n - neighbours
       o - order                   g - vertex figure
       K - colour
   E - Edge:
       v - vertex index nums       d - distance from centre
       f - face index nums         D - direction
       a - diheral angle           l - length
       c - central angle           C - centroid
       K - colour
   F - Face:
       v - vertex index nums       d - distance from centre
       n - neighbours              A - area
       N - normal                  C - centroid
       a - angles                  l - lengths
       s - number of sides         p - max nonplanar distance
       P - perimeter               K - colour

Options
%s
  -c <cent> centre of shape in form 'X,Y,Z', or C to use
            centroid (default, '0,0,0')
  -I <idxs> list of elements by index number given as index ranges
            separated by commas, range can be one number or two
            numbers separated by a hyphen (default range numbers: 0 and
            largest index, default argument: any added elements else '-')
  -v <crds> add vertex, coordinates in form 'X,Y,Z'
  -f <vnos> add face, a list of vertex index numbers separated by commas
  -e <vnos> add edges, a list of vertex index numbers separated by commas
            and taken in pairs
  -k        keep orientation, don't try to orient the faces
  -E <type> edges for query, e - explicit edges, i - implicit edges
            a - explicit and implicit (default)
  -o <file> write output to file (default: write to standard output)
  -d <dgts> number of significant digits (default 17) or if negative
            then the number of digits after the decimal point

)",
          prog_name(), help_ver_text);
}

void add_to_list(vector<int> &list, vector<int> &to_add)
{
  list.insert(list.end(), to_add.begin(), to_add.end());
  sort(list.begin(), list.end());
  auto li = unique(list.begin(), list.end());
  list.erase(li, list.end());
}

void oq_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;
  vector<pair<char, char *>> args;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hc:I:v:f:e:kE:o:d:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'E':
      if (strlen(optarg) != 1 || !strchr("eia", *optarg))
        error("reporting edge type must be e, i or a");
      edge_type = *optarg;
      break;

    default:
      args.push_back(pair<char, char *>(c, optarg));
    }
  }

  if (argc - optind < 1)
    error("must give a query");

  if (!strchr("VEF", *argv[optind]))
    error("the first letter of the query must be V, E or F");
  query = argv[optind++];

  if (argc - optind == 1)
    ifile = argv[optind];

  read_or_error(geom, ifile);

  if (edge_type == 'a')
    geom.add_missing_impl_edges();
  else if (edge_type == 'i') {
    geom.clear(EDGES);
    geom.add_missing_impl_edges();
  }

  vector<int> idx_list;
  vector<int> edge(2);
  Vec3d vec;
  int max_elem_sz;
  if (query[0] == 'V')
    max_elem_sz = geom.verts().size();
  else if (query[0] == 'E')
    max_elem_sz = geom.edges().size();
  else
    max_elem_sz = geom.faces().size();

  for (auto &arg : args) {
    c = arg.first;
    char *optarg = arg.second;
    switch (c) {
    case 'c':
      if (strcmp(optarg, "C") == 0)
        center_is_centroid = true;
      else
        print_status_or_exit(center.read_maths(optarg), c);
      break;

    case 'I':
      print_status_or_exit(read_idx_list(optarg, idx_list, max_elem_sz, false),
                           c);
      if (!idx_list.size())
        warning("file contains no elements of this query type", c);
      add_to_list(idxs, idx_list);
      break;

    case 'v':
      print_status_or_exit(vec.read(optarg), c);
      extra_verts.push_back(vec);
      break;

    case 'f':
      print_status_or_exit(
          read_idx_list(optarg, idx_list, geom.verts().size(), true), c);
      extra_faces.push_back(idx_list);
      break;

    case 'e':
      print_status_or_exit(
          read_idx_list(optarg, idx_list, geom.verts().size(), true), c);
      if (!is_even(idx_list.size()))
        error("odd number of vertex index numbers", c);
      for (unsigned int i = 0; i < idx_list.size() / 2; i++) {
        edge[0] = idx_list[2 * i];
        edge[1] = idx_list[2 * i + 1];
        extra_edges.push_back(edge);
      }
      break;

    case 'k':
      orient = false;
      break;

    case 'd':
      print_status_or_exit(read_int(optarg, &sig_digits), c);
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  int query_elems_added;
  if (query[0] == 'V')
    query_elems_added = extra_verts.size();
  else if (query[0] == 'E')
    query_elems_added = extra_edges.size();
  else
    query_elems_added = extra_faces.size();

  idx_list.resize(query_elems_added);
  for (int i = 0; i < query_elems_added; i++)
    idx_list[i] = max_elem_sz + i;
  add_to_list(idxs, idx_list);

  int total_verts = geom.verts().size() + extra_verts.size();
  for (auto &extra_face : extra_faces)
    for (int j : extra_face)
      if (j > total_verts - 1)
        error(msg_str("extra face vertex x%d out of range",
                      j - (int)geom.verts().size()),
              'f');

  for (auto &extra_edge : extra_edges)
    for (int j : extra_edge)
      if (j > total_verts - 1)
        error(msg_str("extra edge vertex x%d out of range",
                      j - (int)geom.verts().size()),
              'e');

  int total_elems = max_elem_sz + query_elems_added;
  for (int idx : idxs)
    if (idx > total_elems - 1)
      error(msg_str("index x%d out of range", idx - max_elem_sz), 'I');

  if (!idxs.size()) { // default, process all if none specified
    read_idx_list("-", idx_list, max_elem_sz, false);
    add_to_list(idxs, idx_list);
  }
  if (!idxs.size()) // default, process all if none specified
    error("no elements of the query type are in the geometry or were added");
}

void vertex_query(FILE *ofile, rep_printer &rep, oq_opts &opts)
{
  vector<void (rep_printer::*)(int)> query_items;
  query_items.push_back(&rep_printer::v_index);
  for (unsigned int i = 1; i < opts.query.size(); i++) {
    switch (opts.query[i]) {
    case 'c':
      query_items.push_back(&rep_printer::v_coords);
      break;
    case 'd':
      query_items.push_back(&rep_printer::v_distance);
      break;
    case 'f':
      query_items.push_back(&rep_printer::v_face_idxs);
      break;
    case 'a':
      query_items.push_back(&rep_printer::v_solid_angle);
      break;
    case 'o':
      query_items.push_back(&rep_printer::v_order);
      break;
    case 'F':
      query_items.push_back(&rep_printer::v_angles);
      break;
    case 'n':
      query_items.push_back(&rep_printer::v_neighbours);
      break;
    case 'g':
      query_items.push_back(&rep_printer::v_figure);
      break;
    case 'K':
      query_items.push_back(&rep_printer::v_color);
      break;
    default:
      opts.error(msg_str("unknown query letter '%c'", opts.query[i]),
                 "V query");
    }
  }

  for (int idx : opts.idxs) {
    for (unsigned int j = 0; j < query_items.size(); j++) {
      if (j)
        fprintf(ofile, ",");
      (rep.*query_items[j])(idx);
    }
    fprintf(ofile, "\n");
  }
}

void edge_query(FILE *ofile, rep_printer &rep, oq_opts &opts)
{
  vector<void (rep_printer::*)(int)> query_items;
  query_items.push_back(&rep_printer::e_index);
  for (unsigned int i = 1; i < opts.query.size(); i++) {
    switch (opts.query[i]) {
    case 'v':
      query_items.push_back(&rep_printer::e_vert_idxs);
      break;
    case 'f':
      query_items.push_back(&rep_printer::e_face_idxs);
      break;
    case 'a':
      query_items.push_back(&rep_printer::e_dihedral_angle);
      break;
    case 'c':
      query_items.push_back(&rep_printer::e_central_angle);
      break;
    case 'd':
      query_items.push_back(&rep_printer::e_distance);
      break;
    case 'C':
      query_items.push_back(&rep_printer::e_centroid);
      break;
    case 'D':
      query_items.push_back(&rep_printer::e_direction);
      break;
    case 'l':
      query_items.push_back(&rep_printer::e_length);
      break;
    case 'K':
      query_items.push_back(&rep_printer::e_color);
      break;
    default:
      opts.error(msg_str("unknown query letter '%c'", opts.query[i]),
                 "E query");
    }
  }

  for (int idx : opts.idxs) {
    for (unsigned int j = 0; j < query_items.size(); j++) {
      if (j)
        fprintf(ofile, ",");
      (rep.*query_items[j])(idx);
    }
    fprintf(ofile, "\n");
  }
}

void face_query(FILE *ofile, rep_printer &rep, oq_opts &opts)
{
  vector<void (rep_printer::*)(int)> query_items;
  query_items.push_back(&rep_printer::f_index);
  for (unsigned int i = 1; i < opts.query.size(); i++) {
    switch (opts.query[i]) {
    case 'v':
      query_items.push_back(&rep_printer::f_vert_idxs);
      break;
    case 'n':
      query_items.push_back(&rep_printer::f_neighbours);
      break;
    case 'N':
      query_items.push_back(&rep_printer::f_normal);
      break;
    case 'a':
      query_items.push_back(&rep_printer::f_angles);
      break;
    case 's':
      query_items.push_back(&rep_printer::f_sides);
      break;
    case 'd':
      query_items.push_back(&rep_printer::f_distance);
      break;
    case 'A':
      query_items.push_back(&rep_printer::f_area);
      break;
    case 'C':
      query_items.push_back(&rep_printer::f_centroid);
      break;
    case 'l':
      query_items.push_back(&rep_printer::f_lengths);
      break;
    case 'p':
      query_items.push_back(&rep_printer::f_max_nonplanar);
      break;
    case 'P':
      query_items.push_back(&rep_printer::f_perimeter);
      break;
    case 'K':
      query_items.push_back(&rep_printer::f_color);
      break;
    default:
      opts.error(msg_str("unknown query letter '%c'", opts.query[i]),
                 "F query");
    }
  }

  for (int idx : opts.idxs) {
    for (unsigned int j = 0; j < query_items.size(); j++) {
      if (j)
        fprintf(ofile, ",");
      (rep.*query_items[j])(idx);
    }
    fprintf(ofile, "\n");
  }
}

int main(int argc, char *argv[])
{
  oq_opts opts;
  opts.process_command_line(argc, argv);
  Geometry &geom = opts.geom;

  if (opts.center_is_centroid)
    opts.center = geom.centroid();

  FILE *ofile = stdout; // write to stdout by default
  if (opts.ofile != "") {
    ofile = fopen(opts.ofile.c_str(), "w");
    if (ofile == nullptr)
      opts.error("could not open output file '%s'", opts.ofile.c_str());
  }

  rep_printer rep(geom, ofile);
  rep.set_sig_dgts(opts.sig_digits);
  rep.set_center(opts.center);
  rep.is_oriented(); // set oriented value before orienting
  if (opts.orient) {
    geom.orient();
    if (GeometryInfo(geom).volume() < 0) // inneficient
      geom.orient_reverse();
  }

  // print_ranges(rep, opts);

  for (auto &extra_vert : opts.extra_verts)
    geom.add_vert(extra_vert);
  vector<int> edge(2);
  for (auto &extra_edge : opts.extra_edges)
    geom.add_edge_raw(extra_edge);
  for (auto &extra_face : opts.extra_faces)
    geom.add_face(extra_face);
  rep.extra_elems_added(opts.extra_verts.size(), opts.extra_edges.size(),
                        opts.extra_faces.size());

  switch (opts.query[0]) {
  case 'V':
    vertex_query(ofile, rep, opts);
    break;
  case 'E':
    edge_query(ofile, rep, opts);
    break;
  case 'F':
    face_query(ofile, rep, opts);
    break;
  }

  if (opts.ofile == "")
    fclose(ofile);

  return 0;
}
