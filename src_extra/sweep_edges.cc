/*
   Copyright (c) 2021, Adrian Rossiter

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
   Name: sweep
   Description: sweep a polygon according to a transformation
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"

#include <cmath>
#include <ctype.h>
#include <map>
#include <set>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>

using std::map;
using std::string;
using std::vector;

using namespace anti;

class Sweep {
private:
  const string valid_types = "RSTMIAaXi"; // supported transformation types
  // transformations, as type character and specification string
  vector<std::pair<char, string>> transformations;
  int steps = 12;                   // total steps in sweep
  bool triangulate = false;         // convert quads to triangles
  bool join_coincident_ends = true; // join ends when coincident
  double local_epsilon = 1e-8;      // epsilon for checking coincidence

public:
  Status set_steps(int stps)
  {
    if (stps < 2)
      return Status::error("number of steps must be 2 or morei");
    steps = stps;
    return Status::ok();
  }
  int get_steps() const { return steps; }
  void set_triangulate(bool trianglate) { triangulate = trianglate; }
  bool get_triangulate() const { return triangulate; }
  void set_join_coincident_ends(bool join) { join_coincident_ends = join; }
  bool get_join_coincident_ends() const { return join_coincident_ends; }
  void set_local_epsilon(double eps) { local_epsilon = eps; }
  double get_local_epsilon() const { return local_epsilon; }
  Status add_transformation(char type, string trans_str);
  Status set_transformation(Trans3d &trans, int step);
  Status make_sweep(Geometry &geom, Geometry &polygon);
};

string replace_tag_all(const string &str, const string &tag,
                       const string &new_text)
{
  string new_str; // string with tags replaced
  size_t pos = 0; // current position in original string

  while (true) {
    auto tag_pos = str.find(tag, pos);

    // add any text before tag or string end (find_pos can be npos)
    new_str += str.substr(pos, tag_pos - pos);

    if (tag_pos == string::npos) // tag not found
      break;

    new_str += new_text;        // replace tag in new string
    pos = tag_pos + tag.size(); // skip past tag in original string
  }

  return new_str;
}

Status Sweep::add_transformation(char type, string trans_str)
{
  if (valid_types.find(type) == string::npos)
    return Status::error(msg_str("unknown type '%c'", type));

  transformations.push_back({type, trans_str});
  return Status::ok();
}

Status Sweep::set_transformation(Trans3d &trans, int step)
{
  Status stat;
  trans = Trans3d::unit();
  for (auto kp : transformations) {
    auto type = kp.first;
    const auto &trans_str_orig = kp.second;

    double frac = double(step) / steps;
    auto trans_str =
        replace_tag_all(trans_str_orig, "FRAC", std::to_string(frac));
    trans_str = replace_tag_all(trans_str, "ANG", std::to_string(frac * 360.0));
    auto msg_pre =
        msg_str("transformation type %c: %s: ", type, trans_str_orig.c_str());

    Trans3d trans_m;
    Trans3d trans_m2;
    vector<double> nums;
    if (type == 'R') {
      if (!(stat = read_double_list(trans_str.c_str(), nums)))
        return Status::error(msg_pre + stat.msg());
      if (nums.size() == 3)
        trans_m2 = Trans3d::rotate(deg2rad(nums[0]), deg2rad(nums[1]),
                                   deg2rad(nums[2]));
      else if (nums.size() == 4)
        trans_m2 =
            Trans3d::rotate(Vec3d(nums[0], nums[1], nums[2]), deg2rad(nums[3]));
      else if (nums.size() == 6)
        trans_m2 = Trans3d::rotate(Vec3d(nums[0], nums[1], nums[2]),
                                   Vec3d(nums[3], nums[4], nums[5]));
      else if (nums.size() == 12)
        trans_m2 = Trans3d::align(Vec3d(nums[0], nums[1], nums[2]),
                                  Vec3d(nums[3], nums[4], nums[5]),
                                  Vec3d(nums[6], nums[7], nums[8]),
                                  Vec3d(nums[9], nums[10], nums[11]));
      else
        return Status::error(
            msg_pre +
            msg_str("must give 3, 4, 6 of 12 numbers (%lu were given)",
                    (unsigned long)nums.size()));
      trans_m = trans_m2 * trans_m;
    }
    else if (type == 'S') {
      if (!(stat = read_double_list(trans_str.c_str(), nums)))
        return Status::error(msg_pre + stat.msg());
      if (nums.size() == 1)
        trans_m2 = Trans3d::scale(nums[0]);
      else if (nums.size() == 3)
        trans_m2 = Trans3d::scale(nums[0], nums[1], nums[2]);
      else if (nums.size() == 4)
        trans_m2 = Trans3d::scale(Vec3d(nums[0], nums[1], nums[2]), nums[3]);
      else
        return Status::error(
            msg_pre + msg_str("must give 1, 3 or 4 numbers (%lu were given)",
                              (unsigned long)nums.size()));

      trans_m = trans_m2 * trans_m;
    }
    else if (type == 'T') {
      if (!(stat = read_double_list(trans_str.c_str(), nums)))
        return Status::error(msg_pre + stat.msg());
      if (nums.size() != 3)
        return Status::error(
            msg_pre +
            msg_str("must give exactly three numbers (%lu were given)",
                    (unsigned long)nums.size()));

      trans_m2 = Trans3d::translate(Vec3d(nums[0], nums[1], nums[2]));
      trans_m = trans_m2 * trans_m;
    }
    else if (type == 'M') {
      if (!(stat = read_double_list(trans_str.c_str(), nums)))
        return Status::error(msg_pre + stat.msg());
      if (nums.size() != 3)
        return Status::error(
            msg_pre +
            msg_str("must give exactly three numbers (%lu were given)",
                    (unsigned long)nums.size()));

      trans_m2 = Trans3d::reflection(Vec3d(nums[0], nums[1], nums[2]));
      trans_m = trans_m2 * trans_m;
    }
    else if (type == 'I') {
      trans_m = Trans3d::inversion() * trans_m;
    }
    else if (type == 'A') {
      if (!(stat = read_double_list(trans_str.c_str(), nums)))
        return Status::error(msg_pre + stat.msg());
      if (nums.size() != 3)
        return Status::error(
            msg_pre + msg_str("must give exactly 18 numbers (%lu were given)",
                              (unsigned long)nums.size()));

      trans_m2 = Trans3d::align(
          Vec3d(nums[0], nums[1], nums[2]), Vec3d(nums[3], nums[4], nums[5]),
          Vec3d(nums[6], nums[7], nums[8]), Vec3d(nums[9], nums[10], nums[11]),
          Vec3d(nums[12], nums[13], nums[14]),
          Vec3d(nums[15], nums[16], nums[17]));
      trans_m = trans_m2 * trans_m;
    }
    else if (type == 'a') {
      if (!(stat = read_double_list(trans_str.c_str(), nums)))
        return Status::error(msg_pre + stat.msg());
      if (nums.size() != 3)
        return Status::error(
            msg_pre +
            msg_str("must give exactly three numbers (%lu were given)",
                    (unsigned long)nums.size()));
      bool valid;
      trans_m2 = Trans3d::angles_between_axes(
          deg2rad(nums[0]), deg2rad(nums[1]), deg2rad(nums[2]), &valid);
      if (!valid)
        return Status::error(
            msg_pre +
            "the sum of any two angles must be greater than the third");

      trans_m = trans_m2 * trans_m;
    }
    else if (type == 'X') {
      if (!(stat = read_double_list(trans_str.c_str(), nums)))
        return Status::error(msg_pre + stat.msg());
      trans_m2 = Trans3d::unit();
      if (nums.size() == 9) {
        trans_m2[0] = nums[0];
        trans_m2[1] = nums[1];
        trans_m2[2] = nums[2];
        trans_m2[4] = nums[3];
        trans_m2[5] = nums[4];
        trans_m2[6] = nums[5];
        trans_m2[8] = nums[6];
        trans_m2[9] = nums[7];
        trans_m2[10] = nums[8];
      }
      else if (nums.size() == 12) {
        for (int i = 0; i < 12; i++)
          trans_m2[i] = nums[i];
      }
      else
        return Status::error(
            msg_pre + msg_str("must give 9 or 12 numbers (%lu were given)",
                              (unsigned long)nums.size()));

      trans_m = trans_m2 * trans_m;
    }
    else if (type == 'i') {
      trans_m = trans_m.inverse();
    }
    else
      return Status::error(msg_str("option -%c: unhandled transformation"),
                           type);
    trans = trans_m * trans;
  }

  return stat;
}

Geometry color_by_join(const Geometry &poly, const map<int, int> &vmap,
                       const vector<vector<int>> &edges)
{
  Geometry colored_polygon = poly;

  // Each vertex is may be mapped to another, and by following the mappings
  // should form a circuit, which can be given a single colour.
  vector<bool> v_seen(poly.verts().size(), false);
  int col_idx = -1;
  for (int i = 0; i < (int)poly.verts().size(); i++) {
    if (v_seen[i])
      continue; // index was already included in a circuit, try next index
    col_idx++;
    int v_idx = i;
    do {
      v_seen[v_idx] = true;
      colored_polygon.colors(VERTS).set(v_idx, Color(col_idx));
      auto v_idx_it = vmap.find(v_idx);
      if (v_idx_it == vmap.end())
        break; // shouldn't happen
      v_idx = v_idx_it->second;
    } while (!v_seen[(v_idx)]); // check if circuit is complete
  }

  // reverse lookup for vector of edges
  map<vector<int>, int> e2idx;
  int idx = 0;
  for (auto e : edges)
    e2idx[e] = idx++;

  // Each edge may be mapped to another by the mappings of its vertex index
  // numbers, and by following the mappings should form a circuit, which can
  // be given a single colour.
  vector<bool> e_seen(poly.verts().size(), false);
  col_idx = -1;
  for (int i = 0; i < (int)poly.edges().size(); i++) {
    if (e_seen[i])
      continue; // index was already included in a circuit, try next index
    col_idx++;
    int e_idx = i;
    do {
      e_seen[e_idx] = true;
      colored_polygon.colors(EDGES).set(e_idx, Color(col_idx));
      auto v_to0_it = vmap.find(edges[e_idx][0]);
      auto v_to1_it = vmap.find(edges[e_idx][1]);
      if (v_to0_it == vmap.end() || v_to1_it == vmap.end())
        break; // shouldn't happen
      vector<int> edge_to = make_edge(v_to0_it->second, v_to1_it->second);
      auto e_to_it = e2idx.find(edge_to);
      if (e_to_it == e2idx.end())
        break; // shouldn't happen
      e_idx = e_to_it->second;
    } while (!e_seen[(e_idx)]); // check if circuit is complete
  }

  return colored_polygon;
}

Status Sweep::make_sweep(Geometry &geom, Geometry &polygon)
{
  Status stat;
  geom.clear_all();
  polygon.add_missing_impl_edges();
  auto color_polygon = polygon;
  auto edges = polygon.edges();

  Trans3d trans_first;
  if (!(stat = set_transformation(trans_first, 0)))
    return stat;
  Trans3d trans_last;
  if (!(stat = set_transformation(trans_last, steps)))
    return stat;

  map<int, int> vmap_geom; // vert to vert merge correspondence in geom
  if (join_coincident_ends) {
    // Check if first and last set of polygon vertices are coincident
    Geometry poly_first = polygon;
    poly_first.transform(trans_first);
    Geometry poly_last = polygon;
    poly_last.transform(trans_last);
    map<int, int> vmap_poly; // vert to vert correspondence in poly
    vector<map<int, std::set<int>>> equiv_elems;
    if (check_congruence(poly_first, poly_last, &equiv_elems, local_epsilon)) {
      // Set up the vertex correspondence
      const int v_from_offset = (steps - 1) * polygon.verts().size();
      for (auto &v_orbit : equiv_elems[VERTS]) {
        int v_to = *v_orbit.second.begin();
        for (auto v_from : v_orbit.second) {
          if (v_to != v_from) {
            vmap_poly[v_from - polygon.verts().size()] = v_to;
            vmap_geom[v_from + v_from_offset] = v_to;
          }
        }
      }

      color_polygon = color_by_join(polygon, vmap_poly, edges);
    }
  }

  Geometry pverts = color_polygon;
  pverts.clear(EDGES);
  pverts.clear(FACES);

  for (int i = 0; i < steps + 1; i++) {
    Trans3d trans;
    if (!(stat = set_transformation(trans, i)))
      return stat;

    if (i == 0)
      trans_first = trans;
    if (i == steps)
      trans_last = trans;

    auto next_verts = pverts;
    next_verts.transform(trans);
    geom.append(next_verts);
    if (i < steps) {
      int off1 = pverts.verts().size() * i;
      int off2 = pverts.verts().size() * (i + 1);
      for (int e_idx = 0; e_idx < (int)edges.size(); e_idx++) {
        const auto &edge = edges[e_idx];
        Color e_col = color_polygon.colors(EDGES).get(e_idx);
        if (triangulate) {
          geom.add_face({edge[0] + off1, edge[1] + off1, edge[1] + off2},
                        e_col);
          geom.add_face({edge[0] + off1, edge[1] + off2, edge[0] + off2},
                        e_col);
        }
        else
          geom.add_face(
              {edge[0] + off1, edge[1] + off1, edge[1] + off2, edge[0] + off2},
              e_col);
      }
      for (int v_idx = 0; v_idx < (int)pverts.verts().size(); v_idx++)
        geom.add_edge(v_idx + off1, v_idx + off2,
                      pverts.colors(VERTS).get(v_idx));
    }
  }
  if (vmap_geom.size())
    geom.verts_merge(vmap_geom); // merge first and last set of poly verts

  return stat;
}

class sw_opts : public ProgramOpts {
public:
  Sweep sweep;
  Coloring clrngs[3];
  string ifile;
  string ofile;

  sw_opts() : ProgramOpts("sweep")
  {
    for (int i = 0; i < 3; i++)
      clrngs[i].add_cmap(colormap_from_name("spread"));
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

void sw_opts::usage()
{
  const Sweep def_sweep;
  fprintf(stdout, R"(
Usage: %s [options] input

Read input model in OFF format and transform it multiple times to produce
a polygonal surface made of the quadrilaterals swept by the edges. At
each step the input model is transformed according to the specified
transformations, which can be defined using sweep progress variables FRAC
and ANG. FRAC runs from 0.0 at at step 0 to 1.0 on the final step.
  FRAC = current_step / total_steps
  ANG = 360.0 * FRAC
Colours are taken from the input model if the begining of the sweep does
not coincide with the end, or if option -O is set, otherwise sweep circuits
are found and the colours are set from the colour map (option -m). Vertex
colours are used for final vertices and swept vertex edges, edge colours
are used for swept edge faces.
If input is not given the program reads from standard input.

Options
%s
  -n <int>  number of steps (default: %d)
  -T <tran> translate, three numbers separated by commas which are
            used as the x, y and z displacements
  -R <rot>  rotate about an axis, three, four or six numbers separated by
            commas. If three numbers these are angles (degrees) to rotate
            about the x, y and z axes. If four numbers, the first three
            are a direction vector for the axis, the last number is the
            angle (degrees) to rotate. If six numbers, these are two
            vectors (from,to) and rotate to carry the first to the second.
            If twelve numbers these are four vectors (from1,from2,to1,to2)
            and rotate to carry the first onto the third then rotate around
            the third to carry the second onto the fourth
  -M <norm> reflect in a plane, three numbers separated by commas which
            give a vector normal to the plane of reflection.
  -S <scal> scale, one, three or four numbers separated by commas. If one
            number then scale by this factor in all directions. If three
            numbers these are the factors to scale along the x, y and
            z axes. If four numbers, the first three are a direction
            vector for the scaling, the last number is the factor to scale
  -I        inversion
  -A <crds> transformation that will align two sets of three points
            (18 numbers coordinates of from1,from2,from3,to1,to2,to3)
  -a <angs> transformation that makes particular angles between the
            mapped axes, angles in degrees in form yz_ang,zx_ang,xy_ang
            (corresponding to the angles opposite the x-, y- and z-axis)
  -X <mtrx> transformation matrix of 9 or 12 values, given left-to-right
            top-to-bottom, used to premultipy each coordinate
  -i        replace the current combined transformation by its inverse
  -t        triangulate
  -O        do not join coincident sweep ends
  -l <lim>  maximum distance between vertices that should be coincident, lim
            is an integer represeting the negative exponent of the distance
            (default: %d giving %.0e)
  -m <maps> a comma separated list of colour maps used to transform colour
            indexes (default: rand), a part consisting of letters from
            v, e, f, selects the element types to apply the map list to
            (default 'vef'). The 'compound' map should give useful results.
  -o <file> write output to file (default: write to standard output)
)",
          prog_name(), help_ver_text, def_sweep.get_steps(),
          int(-log(def_sweep.get_local_epsilon()) / log(10) + 0.5),
          def_sweep.get_local_epsilon());
}

void sw_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;
  int num;
  vector<double> nums;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hR:S:T:M:IA:a:X:in:tOl:m:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'R':
    case 'S':
    case 'T':
    case 'M':
    case 'I':
    case 'A':
    case 'a':
    case 'X':
    case 'i':
      print_status_or_exit(sweep.add_transformation(c, optarg));
      break;

    case 'n':
      print_status_or_exit(read_int(optarg, &num), c);
      print_status_or_exit(sweep.set_steps(num), c);
      break;

    case 't':
      sweep.set_triangulate(true);
      break;

    case 'O':
      sweep.set_join_coincident_ends(false);
      break;

    case 'l':
      int sig_compare;
      print_status_or_exit(read_int(optarg, &sig_compare), c);
      if (sig_compare > DEF_SIG_DGTS)
        warning("limit is very small, may not be attainable", c);
      sweep.set_local_epsilon(pow(10, -sig_compare));
      break;

    case 'o':
      ofile = optarg;
      break;

    case 'm':
      print_status_or_exit(read_colorings(clrngs, optarg), c);
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

int main(int argc, char *argv[])
{
  sw_opts opts;
  opts.process_command_line(argc, argv);

  Geometry polygon;
  opts.read_or_error(polygon, opts.ifile);

  Geometry geom;
  opts.print_status_or_exit(opts.sweep.make_sweep(geom, polygon));

  for (int i = 0; i < 3; i++)
    opts.clrngs[i].set_geom(&geom);

  opts.clrngs[0].v_apply_cmap();
  opts.clrngs[1].e_apply_cmap();
  opts.clrngs[2].f_apply_cmap();

  opts.write_or_error(geom, opts.ofile);

  return 0;
}
