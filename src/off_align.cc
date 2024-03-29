/*
   Copyright (c) 2003-2020, Adrian Rossiter

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
   Name: off_align.cc
   Description: position one OFF file relative to another
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

using std::map;
using std::pair;
using std::set;
using std::string;
using std::vector;

using namespace anti;

struct bond_brick {
  Geometry geom;
  int align_type;
  vector<int> bond;
};

class bond_base {
private:
  Geometry base;
  Symmetry sym;
  vector<bond_brick> bricks;

public:
  enum { align_verts, align_faces, align_faces_merge };
  enum { out_default = 0, out_brick, out_base_brick, out_brick_base };
  bond_base(Geometry &bas) : base(bas) {}
  Status set_sym(const char *sym_str);
  Status add_brick(char type, const string &brick_str);
  Status bond_all(Geometry &geom_out, int out_type);
};

Status bond_base::set_sym(const char *sym_str)
{
  Status stat;

  Symmetry full_sym(base);
  Split parts(sym_str, ",");
  if (parts.size() == 0 || parts.size() > 2)
    return Status::error("argument should have 1 or 2 comma separated parts");

  Symmetry sub_sym;
  if (strncmp(parts[0], "full", strlen(parts[0])) == 0)
    sub_sym = full_sym;
  else if (!(stat = sub_sym.init(parts[0], Trans3d())))
    return Status::error(msg_str("sub-symmetry type: %s", stat.c_msg()));

  int sub_sym_conj = 0;
  if (parts.size() > 1 && !(stat = read_int(parts[1], &sub_sym_conj)))
    return Status::error(
        msg_str("sub-symmetry conjugation number: %s", stat.c_msg()));

  if (!(stat = full_sym.get_sub_sym(sub_sym, &sym, sub_sym_conj)))
    return Status::error(msg_str("sub-symmetry: %s", stat.c_msg()));

  return Status::ok();
}

int get_aligns(const vector<double> &angs0, const vector<double> &angs1,
               vector<int> &aligns)
{
  vector<double> test = angs1;
  AngleVectLess cmp_face_angles;

  reverse(test.begin(), test.end());
  for (unsigned int i = 0; i < angs0.size(); i++) {
    if (cmp_face_angles(angs0, test) == 0)
      aligns.push_back(-1 - i);
    rotate(test.begin(), test.begin() + 1, test.end());
  }

  reverse(test.begin(), test.end());
  for (unsigned int i = 0; i < angs0.size(); i++) {
    if (cmp_face_angles(angs0, test) == 0)
      aligns.push_back(i);
    rotate(test.begin(), test.begin() + 1, test.end());
  }
  return aligns.size();
}

Status bond_base::add_brick(char type, const string &brick_str)
{
  auto first_comma_pos = brick_str.find(','); // split on first comma
  string geom_str = brick_str.substr(0, first_comma_pos); // first part
  string bond_idxs_str;
  if (first_comma_pos != string::npos)
    bond_idxs_str = brick_str.substr(first_comma_pos + 1); // all after comma

  bricks.push_back(bond_brick());
  bond_brick &brick = bricks.back();

  Status stat;
  if (geom_str.empty()) // brick geometry specifier is empty, use base geometry
    brick.geom = base;
  else {
    if (!(stat = brick.geom.read(geom_str)))
      return Status::error("brick geometry: " + stat.msg());
  }

  if (!(stat = read_int_list(bond_idxs_str.c_str(), brick.bond, true)))
    return Status::error("bond values: " + stat.msg());

  if (type == 'v') {
    brick.align_type = align_verts;
    int n = brick.bond.size();
    if (n != 2 && n != 4 && n != 6)
      return Status::error(
          msg_str("must give 2, 4 or 6 vertices (%d were given)", n));

    if (n > 2) {
      for (int i = 0; i < n / 2; i++) {
        if (brick.bond[i] == brick.bond[(i + 1) % (n / 2)])
          return Status::error("repeated vertex index in base");
        if (brick.bond[n / 2 + i] == brick.bond[n / 2 + (i + 1) % (n / 2)])
          return Status::error("repeated vertex index in base");
      }
    }
    for (unsigned int i = 0; i < brick.bond.size(); i++) {
      const vector<Vec3d> &vs =
          (i < brick.bond.size() / 2) ? base.verts() : brick.geom.verts();
      if (brick.bond[i] < 0 || brick.bond[i] >= (int)vs.size())
        return Status::error(
            msg_str("bond values: vertex %d (position %d) is out of bounds",
                    brick.bond[i], i + 1));
    }
  }
  else if (type == 'f' || type == 'F') {
    brick.align_type = (type == 'f') ? align_faces : align_faces_merge;
    if (brick.bond.size() > 3)
      return Status::error(
          msg_str("up to three arguments can be given (%d were given)",
                  (int)brick.bond.size()));
    brick.bond.resize(3, 0);
    int f0 = brick.bond[0];
    int f1 = brick.bond[1];
    if (f0 < 0 || f0 >= (int)base.faces().size()) {
      if (base.faces().size())
        return Status::error(msg_str("invalid base face '%d', last face is %d",
                                     f0, (int)base.faces().size() - 1));
      else
        return Status::error("base has no faces");
    }
    if (f1 < 0 || f1 >= (int)brick.geom.faces().size()) {
      if (brick.geom.faces().size())
        return Status::error(msg_str("invalid brick face '%d', last face is %d",
                                     f1, (int)brick.geom.faces().size() - 1));
      else
        return Status::error("brick has no faces");
    }
    if (brick.align_type == align_faces_merge &&
        base.faces(f0).size() != brick.geom.faces(f1).size())
      return Status::error("faces to bond have different number of sides");

    vector<int> aligns;
    vector<double> angs0, angs1;
    base.face_angles_lengths(f0, &angs0);
    brick.geom.face_angles_lengths(f1, &angs1);
    if (!get_aligns(angs0, angs1, aligns)) {
      return Status::error(
          "base and brick bonding faces are not the same shape\n");
    }
    if (brick.bond[2] >= (int)aligns.size()) {
      return Status::error(
          msg_str("bond selection number is %d, last selection is %d\n",
                  brick.bond[2], (int)aligns.size() - 1));
    }
    vector<int> &face = brick.geom.raw_faces()[f1];
    int sel = aligns[brick.bond[2]];
    if (sel >= 0)
      rotate(face.begin(), face.begin() + sel, face.end());
    else {
      brick.geom.orient_reverse();
      rotate(face.begin(), face.begin() + (-sel - 1), face.end());
    }
  }

  return Status::ok();
}

Status bond_base::bond_all(Geometry &geom_out, int out_type)
{
  geom_out.clear_all();
  double base_rad = BoundSphere(base.verts()).get_radius();
  // check for errors before making any transformations
  vector<bond_brick>::iterator bi;
  for (bi = bricks.begin(); bi != bricks.end(); ++bi) {
    if (bi->align_type == align_faces_merge &&
        (out_type == out_brick || out_type == out_brick_base))
      return Status::error(
          msg_str("output type is %s and not compatible with a merged brick",
                  (out_type == out_brick) ? "brick only" : "brick then base"));
  }

  bool has_merged_brick = false;
  bool face_used_more_than_once = false;
  Geometry base_out;
  Geometry brick_out;
  vector<Geometry> merge_bricks;
  map<int, vector<int>> face_params;
  for (bi = bricks.begin(); bi != bricks.end(); ++bi) {
    bond_brick &brick = *bi;
    if (brick.align_type == align_verts) {
      vector<vector<Vec3d>> pts(2);
      for (unsigned int i = 0; i < brick.bond.size(); i++) {
        const vector<Vec3d> &vs =
            (i < brick.bond.size() / 2) ? base.verts() : brick.geom.verts();
        pts[!(i < brick.bond.size() / 2)].push_back(vs[brick.bond[i]]);
      }
      if (pts[0].size() > 1) {
        Trans3d scl = Trans3d::scale((pts[0][0] - pts[0][1]).len() /
                                     (pts[1][0] - pts[1][1]).len());
        brick.geom.transform(scl);
        for (auto &i : pts[1])
          i = scl * i;
      }
      brick.geom.transform(Trans3d::align(pts[1], pts[0]));
      brick_out.append(brick.geom);
    }
    else if (brick.align_type == align_faces ||
             brick.align_type == align_faces_merge) {
      int f0 = brick.bond[0];
      int f1 = brick.bond[1];

      double e0 = base.edge_len(base.faces(f0)); // first edge
      double e1 = brick.geom.edge_len(
          make_edge(brick.geom.faces(f1, 0), brick.geom.faces(f1, 1)));
      brick.geom.transform(Trans3d::scale(e0 / e1));

      if (bi->align_type == align_faces) {
        face_bond(base, brick.geom, f0, f1, 0, false, true);
        brick_out.append(brick.geom);
      }
      else {
        if (!has_merged_brick) {
          has_merged_brick = true;
          base_out = base;
        }
        double model_rad =
            base_rad + BoundSphere(bi->geom.verts()).get_radius();
        merge_bricks.push_back(brick.geom);
        Transformations ts = sym.get_trans();
        if (!ts.size()) // set to unit if not set
          ts.add(Trans3d());
        for (const auto &t : ts) {
          vector<vector<int>> elem_maps;
          get_coincidence_maps(base, t, elem_maps, model_rad * sym_eps);
          const int f0_map = elem_maps[2][f0];
          const vector<int> &f0_mapface = base.faces(f0_map);
          int f0_sz = base.faces(f0).size();
          if (f0_sz == 0)
            return Status::error("a base bond face has no vertices");
          auto vi = face_params.find(f0_map);
          bool direct = Isometry(t).is_direct();
          if (vi == face_params.end() || (!(vi->second)[2] && direct)) {
            if (vi != face_params.end())
              face_used_more_than_once = true;
            int map_offset;
            const int idx0 = base.faces(f0)[0];
            for (map_offset = 0; map_offset < f0_sz; map_offset++) {
              if (elem_maps[0][idx0] == f0_mapface[map_offset])
                break;
            }
            const int idx1 = base.faces(f0)[1];
            bool rev =
                (elem_maps[0][idx1] != f0_mapface[(map_offset + 1) % f0_sz]);
            face_params[f0_map] = vector<int>(5);
            face_params[f0_map][0] = f1;
            face_params[f0_map][1] = map_offset;
            face_params[f0_map][2] = merge_bricks.size() - 1;
            face_params[f0_map][3] = direct;
            face_params[f0_map][4] = rev;
          }
        }
      }
    }
  }

  if (face_params.size()) {
    map<int, vector<int>>::reverse_iterator mi;
    for (mi = face_params.rbegin(); mi != face_params.rend(); ++mi) {
      Geometry bgeom = merge_bricks[mi->second[2]];
      vector<int> &brick_f = bgeom.raw_faces()[mi->second[0]];
      rotate(brick_f.begin(), brick_f.begin() + mi->second[1], brick_f.end());
      if (!mi->second[3]) // symmetry was indirect
        bgeom.transform(Trans3d::inversion());
      if (mi->second[4]) // reverse was set
        bgeom.orient_reverse();
      face_bond(base_out, bgeom, mi->first, mi->second[0], 0, true, true);
    }
  }

  if (sym.is_set())
    sym_repeat(brick_out, brick_out, sym);

  if (out_type == out_default) {
    if (has_merged_brick)
      geom_out.append(base_out);
    geom_out.append(brick_out);
  }
  else if (out_type == out_brick)
    geom_out.append(brick_out);
  else if (out_type == out_base_brick) {
    if (has_merged_brick)
      geom_out.append(base_out);
    else
      geom_out.append(base);
    geom_out.append(brick_out);
  }
  else if (out_type == out_brick_base) {
    geom_out.append(brick_out);
    geom_out.append(base);
  }

  if (face_used_more_than_once)
    return Status::warning(
        "face merge: the same bond face was specified for merging in more than "
        "one call, only the last specified bond will be used");
  else
    return Status::ok();
}

class align_opts : public ProgramOpts {
public:
  bool has_merge;
  vector<pair<char, string>> brick_args;
  string sym_str;
  int out_type;
  string ifile;
  string ofile;

  align_opts() : ProgramOpts("off_align"), has_merge(false), out_type(0) {}

  void process_command_line(int argc, char **argv);
  void usage();
};

void align_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] [input_file]

Read a base file and brick file in OFF format, and use vertices or faces
to position the brick with respect to the base. The brick may be repeated
symmetrically and/or merged with the base. If input_file is not given the
program reads from standard input.

Options
%s
  -v <arg>  align by vertices, arg is a comma separated list of a brick
            geometry (if empty use base) optionally followed by 'r' (reverse
            brick orientation) followed by two, four or six points given as
            the vertex index number in the OFF files (starting at 0). The base
            points are all given first and then the brick points.
            Formats:
               u1,v1 - point alignment
                  translation u1-v1 (so v1 of brick moves to u1 of base)
               u1,u2,v1,v2 - line alignment
                  point alignment as above followed by rotation through
                  u1 perpendicular to u1u2 and v1v2 to align u1u2 and v1v2.
               u1,u2,u3,v1,v2,v3 - face alignment
                  line alignment as above followed by a rotation
                  around u1,u2 so v3 lies in plane of u1u2u3.
  -f <arg>  align by face index, arg is a comma separated list of a brick
            geometry (if empty use base) followed by up to three numbers
            separated by commas: base face index, brick face index
            (default: 0), polygon alignment selection number (default: 0)
  -F <arg>  align and combine polyhedra by face index, arg is a comma
            separated list of a brick geometry (if empty use base) followed
            by up to three numbers: base face index, brick face index
            (default: 0), polygon alignment selection number (default: 0).
            Brick is after base, bond vertices are merged, bond faces are
            removed
  -M <val>  merge parts, select and order parts in the output, val may be:
               default (0):    any combined part (-F) followed by any brick
                               parts (-v, -f)
               brick (1):      brick parts only
               base_brick (2): base part, possibly combined (-F). followed
                               any brick parts (-v, -f)
               brick_base (3): brick parts (-v, -f), followed by base part
  -y <sub>  repeat bricks according to symmetry of base. sub is symmetry
            subgroup (Schoenflies notation) or 'full' optionally followed
            by a ',' and conjugation type (integer)
  -o <file> write output to file (default: write to standard output)

)",
          prog_name(), help_ver_text);
}

void align_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;
  string arg_id;
  Trans3d trans_m2;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hv:f:F:M:y:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'F':
      has_merge = true;
      // fall through
    case 'v':
    case 'f':
      brick_args.push_back(pair<char, string>(c, optarg));
      break;

    case 'M':
      print_status_or_exit(get_arg_id(optarg, &arg_id,
                                      "default|brick|base_brick|brick_base",
                                      argmatch_add_id_maps),
                           c);
      out_type = atoi(arg_id.c_str());
      break;

    case 'y':
      sym_str = optarg;
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (brick_args.size() == 0)
    error("no brick alignments given, must use option -v, -f or -F");

  if (argc - optind > 1)
    error("too many arguments");

  if (argc - optind == 1)
    ifile = argv[optind];

  // Cannot read more than one geometry from stdin
  bool base_from_stdin = (ifile == "" || ifile == "-");
  bool stdin_in_use = base_from_stdin;

  for (auto &brick_arg : brick_args) {
    const string &brick = brick_arg.second;
    if (brick[0] == '-' && (brick.size() == 1 || brick[1] == ',')) {
      if (stdin_in_use) {
        if (base_from_stdin)
          error("cannot read both brick and base from standard input",
                brick_arg.first);
        else
          error("cannot read more than one brick from standard input",
                brick_arg.first);
      }
      else
        stdin_in_use = true;
    }
  }
}

int main(int argc, char *argv[])
{
  align_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);

  if (opts.has_merge) {
    unsigned int num_verts = geom.verts().size();
    merge_coincident_elements(geom, "v");
    if (geom.verts().size() != num_verts)
      opts.warning("coincident vertices in base geometry have been merged",
                   'F');
  }

  bond_base base(geom);
  if (opts.sym_str != "")
    opts.print_status_or_exit(base.set_sym(opts.sym_str.c_str()), 'y');

  vector<pair<char, string>>::iterator argi;
  for (argi = opts.brick_args.begin(); argi != opts.brick_args.end(); ++argi)
    opts.print_status_or_exit(base.add_brick(argi->first, argi->second),
                              argi->first);

  Geometry geom_out;
  opts.print_status_or_exit(base.bond_all(geom_out, opts.out_type), 'M');

  opts.write_or_error(geom_out, opts.ofile);

  return 0;
}
