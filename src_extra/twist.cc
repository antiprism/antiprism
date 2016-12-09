/*
   Copyright (c) 2003-2009, Adrian Rossiter

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
   Name: twist.cc
   Description: twist edges of a polyhedron
   Project: Antiprism - http://www.antiprism.com
*/

#include <algorithm>
#include <ctype.h>
#include <map>
#include <math.h>
#include <string.h>
#include <string>
#include <utility>
#include <vector>

#include "../base/antiprism.h"
#include "../base/vec_utils.h"

using std::string;
using std::vector;
using std::map;
using std::pair;
using std::swap;

using namespace anti;

class tw_opts : public ProgramOpts {
private:
public:
  bool read_centroid_verts(char *optarg, char *errmsg);
  Vec3d centre;
  vector<int> centroid_verts;
  double twist_val;
  bool struts_only;

  string ifile;
  string ofile;

  tw_opts()
      : ProgramOpts("twist"), centre(Vec3d(0, 0, 0)), twist_val(-1),
        struts_only(false)
  {
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

void tw_opts::usage()
{
  fprintf(
      stdout,
      "\n"
      "Usage: %s [options] [input_file]\n"
      "\n"
      "Read a file in OFF format and twist the edges.\n"
      "If input_file is not given the program reads from standard input.\n"
      "\n"
      "Options\n"
      "%s"
      "  -c <cent> centre to twist about, in form \"x_val,y_val,z_val\"\n"
      "            (default, \"0,0,0\")\n"
      "  -t val    twist edges an amount between polyhedron (0.0) and dual "
      "(1.0)\n"
      "  -s        output struts only\n"
      "  -o <file> write output to file (default: write to standard output)\n"
      "\n"
      "\n",
      prog_name(), help_ver_text);
}

bool tw_opts::read_centroid_verts(char *optarg, char *errmsg)
{
  // vector<int> vecs;
  int vec_idx;
  char *v_str = strtok(optarg, ",");
  int i = 0;
  while (v_str) {
    i++;
    if (!read_int(v_str, &vec_idx) || vec_idx < 0) {
      snprintf(errmsg, MSG_SZ,
               "vertex index %d, '%s' is not a positive integer", i, v_str);
      return false;
    }
    centroid_verts.push_back(vec_idx);
    v_str = strtok(NULL, ",");
  }

  return true;
}

void tw_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hc:o:t:s")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'c':
      print_status_or_exit(centre.read(optarg), c);
      break;

    case 'o':
      ofile = optarg;
      break;

    case 't':
      print_status_or_exit(read_double(optarg, &twist_val), c);
      break;

    case 's':
      struts_only = true;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (argc - optind > 1) {
    error("too many arguments");
    exit(1);
  }

  if (argc - optind == 1)
    ifile = argv[optind];
}

Geometry twist(Geometry &poly, Geometry &dual, double twist_val, Vec3d centre,
               bool struts_only)
{

  map<vector<int>, vector<int>> edges;
  poly.get_edge_face_pairs(edges, true);
  map<vector<int>, vector<int>>::iterator mi, mi_next;
  Geometry twist;

  Vec3d v0, v1, vm; // strut ends and middle
  double ratio;
  Vec3d v0p, v1p, v0d, v1d; // vertices at ends 1 and 2 in poly and dual
  map<vector<int>, vector<int>> attach_edges;
  vector<int> tv_map;
  map<vector<int>, pair<Vec3d, Vec3d>> twist_edges;
  map<vector<int>, pair<Vec3d, Vec3d>>::iterator mi_tw, mi_tw2;
  for (mi = edges.begin(); mi != edges.end(); mi++) {
    v0p = poly.verts(mi->first[0]);
    v1p = poly.verts(mi->first[1]);
    v0d = dual.verts(mi->second[1]);
    v1d = dual.verts(mi->second[0]);
    // don't need this
    // ratio = (v1p - v0p).len() / (v1d - v0d).len();
    // v0d = centre + (v0d - centre)*ratio;
    // v1d = centre + (v1d - centre)*ratio;
    Vec3d norm = vcross(v1p - v0p, v1d - v0d);
    Trans3d trans = Trans3d::transl(-0.5 * (v0p + v1p));
    trans = Trans3d::rot(norm, twist_val * M_PI / 2) * trans;
    trans = Trans3d::transl(0.5 * (v0p + v1p)) * trans;
    trans = Trans3d::transl((v0d + v1d - v0p - v1p) * 0.5 * twist_val) * trans;

    v0 = trans * v0p;
    v1 = trans * v1p;

    // v0 = v0p + (v0d - v0p)*twist_val;
    // v1 = v1p + (v1d - v1p)*twist_val;

    twist_edges[mi->first] = pair<Vec3d, Vec3d>(v0, v1);
    vector<int> ve = mi->first;
    if (ve[1] < ve[0])
      swap(ve[1], ve[0]);
    attach_edges[ve] = vector<int>();
  }

  double edge_len = (v1p - v0p).len();
  vector<int>::const_iterator vi;
  for (mi_tw = twist_edges.begin(); mi_tw != twist_edges.end(); mi_tw++) {
    mi = edges.find(mi_tw->first);
    twist.add_face(vector<int>());
    for (int i = 0; i < 2; i++) {
      vector<int> e_next(2);
      e_next[0] = mi_tw->first[(i + 1) % 2];
      const vector<int> &face = poly.faces(mi->second[i]);
      vi = find(face.begin(), face.end(), e_next[0]);
      e_next[1] = ++vi != face.end() ? *vi : face.front();
      // fprintf(stderr, "twist_edge/next (%d, %d)/(%d, %d)\n", mi_tw->first[0],
      // mi_tw->first[1], e_next[0], e_next[1]);
      if (e_next[1] < e_next[0])
        swap(e_next[1], e_next[0]);
      mi_tw2 = twist_edges.find(e_next);
      attach_edges[e_next].push_back(twist.verts().size());
      int where;
      Vec3d v = line_plane_intersect(centre, mi_tw2->second.first,
                                     mi_tw2->second.second, mi_tw->second.first,
                                     mi_tw->second.second, &where);
      tv_map.push_back(mi_tw->first[i]);
      twist.add_vert(v);
      twist.faces(twist.faces().size() - 1).push_back(twist.verts().size() - 1);
    }

    Vec3d &v0 = twist.verts(twist.verts().size() - 2);
    Vec3d &v1 = twist.verts(twist.verts().size() - 1);
    ratio = edge_len / (v1 - v0).len();
    v0 = centre + (v0 - centre) * ratio;
    v1 = centre + (v1 - centre) * ratio;
  }

  if (struts_only) {
    for (unsigned int i = 0; i < twist.faces().size(); i++)
      twist.add_edge(make_edge(twist.faces(i, 0), twist.faces(i, 1)), 0);
    twist.clear(FACES);
  }
  else {
    for (unsigned int i = 0; i < twist.faces().size(); i++) {
      vector<int> edge(2);
      edge[0] = tv_map[twist.faces(i, 0)];
      edge[1] = tv_map[twist.faces(i, 1)];
      if (edge[1] < edge[0])
        swap(edge[1], edge[0]);
      // fprintf(stderr, "edge = (%d, %d) edge.size = %d\n", edge[0], edge[1],
      // (int)edge.size());
      vector<int> v_idxs = attach_edges.find(edge)->second;
      Vec3d P0 = twist.face_v(i, 0);
      Vec3d edge_vec = twist.face_v(i, 1) - P0;
      if (vdot(edge_vec, twist.verts(v_idxs[0]) - P0) <
          vdot(edge_vec, twist.verts(v_idxs[1]) - P0))
        swap(v_idxs[0], v_idxs[1]);
      twist.faces(i).push_back(v_idxs[0]);
      twist.faces(i).push_back(v_idxs[1]);

      for (int j = 0; j < 4; j++)
        twist.add_edge(
            make_edge(twist.faces(i, j), twist.faces(i, (j + 1) % 4)),
            Color(bool(j)));
    }
  }

  return twist;
}

int main(int argc, char *argv[])
{
  tw_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);
  geom.orient();

  Geometry dual;
  get_dual(&dual, geom, 1, opts.centre);

  vector<int> invalid_verts;
  for (unsigned int i = 0; i < dual.verts().size(); i++) {
    if (!dual.verts(i).is_set()) {
      dual.del(VERTS, i);
      int idx = invalid_verts.size() + i;
      invalid_verts.push_back(idx);
      i--;
    }
  }
  if (invalid_verts.size()) {
    string msg(
        "removed invalid vertices (and associated faces) with indices - ");
    char errmsg[MSG_SZ];
    for (unsigned int i = 0; i < invalid_verts.size() - 1; i++) {
      snprintf(errmsg, MSG_SZ, "%d,", invalid_verts[i]);
      msg += string(errmsg);
    }
    snprintf(errmsg, MSG_SZ, "%d", invalid_verts.back());
    msg += string(errmsg);
    opts.warning(msg);
  }

  Geometry twisted =
      twist(geom, dual, opts.twist_val, opts.centre, opts.struts_only);

  opts.write_or_error(twisted, opts.ofile);

  return 0;
}
