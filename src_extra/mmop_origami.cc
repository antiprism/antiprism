/*
   Copyright (c) 2017, Adrian Rossiter

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
   Name: mmop_origami.cc
   Description: solve multimodular origami polyhedra, described in
     Multimodular Origami Polyhedra: Archimedeans, Buckyballs and Duality
     by Rona Gurkewitz, Bennett Arnstein
     http://store.doverpublications.com/0486423174.html
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"
#include <cmath>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <vector>

using std::string;
using std::vector;

using namespace anti;

class iter_params {
public:
  int num_iters;
  int num_iters_status;
  int sig_digits;
  FILE *rep_file;
  iter_params()
      : num_iters(10000), num_iters_status(1000), sig_digits(16),
        rep_file(stderr)
  {
  }

  double get_test_val() { return pow(10, -sig_digits); }
  bool quiet() { return rep_file == nullptr; }
  bool check_status(int n)
  {
    return n == num_iters || n == 1 ||
           (num_iters_status > 0 && n % num_iters_status == 0);
  }
  bool checking_status() { return num_iters_status > 0; };
  bool print_progress_dot(int n)
  {
    if (!quiet() && num_iters_status >= 10)
      return (n % num_iters) % (num_iters_status / 10) == 0;
    else
      return false;
  }
};

class mmop_opts : public ProgramOpts {
public:
  iter_params it_params;
  double adjust_fact;
  double trunc_len;
  double keep_orient;
  double init_ht;
  char color_method;
  Coloring clrngs[3];

  string ifile;
  string ofile;

  mmop_opts()
      : ProgramOpts("mmop_origami"), adjust_fact(100.0), trunc_len(1.0),
        keep_orient(false), init_ht(-0.5), color_method('n')
  {
    clrngs[2].add_cmap(colormap_from_name("spread"));
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

// clang-format off
void mmop_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Read a file in OFF format containing an oriented polyhedron, and try to\n"
"convert to a multimodular origami polyhedron form with unit-edged polygons.\n"
"The form is described in Multimodular Origami Polyhedra: Archimedeans,\n"
"Buckyballs and Duality by Rona Gurkewitz, Bennett Arnstein\n"
"http://store.doverpublications.com/0486423174.html\n"
"If input_file is not given the program reads from standard input.\n"
"\n"
"Options\n"
"%s"
"  -n <itrs> number of iterations (default 10000)\n"
"  -s <perc> percentage to adjust corrections on iteration (default: 100)\n"
"  -t <val>  truncate polygon edge to this length (default: no truncation\n"
"  -k        keep orientation, affects face centre offset direction (default:\n"
"            set positive orientation)\n"
"  -p <ht>   set initial height of face centres (default: -0.5)\n"
"  -V        colour units with base model vertex colour\n"
"  -m <maps> a comma separated list of colour maps used to transform colour\n"
"            indexes (default: rand), a part consisting of letters from\n"
"            v, e, f, selects the element types to apply the map list to\n"
"            (default 'vef'). The 'compound' map should give useful results.\n"
"  -L <lim>  NOT USED: minimum change of distance/width_of_model to\n"
"               terminate, as negative exponent (default: %d giving %.0e)\n"
"  -z <n>    status checking and reporting every n iterations, -1 for no\n"
"            status (default: 1000)\n"
"  -q        quiet, do not print status messages\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text,
   it_params.sig_digits, it_params.get_test_val() );
}
// clang-format on

void mmop_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;
  vector<double> nums;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hn:s:t:kp:Vm:L:z:qo:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'o':
      ofile = optarg;
      break;

    case 'n':
      print_status_or_exit(read_int(optarg, &it_params.num_iters), c);
      if (it_params.num_iters < 0)
        error("number of iterations must be greater than 0", c);
      break;

    case 'z':
      print_status_or_exit(read_int(optarg, &it_params.num_iters_status), c);
      if (it_params.num_iters_status < -1)
        error("number of iterations must be -1 or greater", c);
      break;

    case 's':
      print_status_or_exit(read_double(optarg, &adjust_fact), c);
      if (adjust_fact < 0 || adjust_fact > 100)
        warning("not in range 0 to 100", c);
      break;

    case 't':
      print_status_or_exit(read_double(optarg, &trunc_len), c);
      break;

    case 'k':
      keep_orient = true;
      break;

    case 'p':
      print_status_or_exit(read_double(optarg, &init_ht), c);
      break;

    case 'V':
      color_method = 'v';
      break;

    case 'm':
      print_status_or_exit(read_colorings(clrngs, optarg), c);
      break;

    case 'L':
      print_status_or_exit(read_int(optarg, &it_params.sig_digits), c);
      if (it_params.sig_digits < 0) {
        warning("termination limit is negative, and so ignored", c);
      }
      if (it_params.sig_digits > DEF_SIG_DGTS) {
        warning("termination limit is very small, may not be attainable", c);
      }
      break;

    case 'q':
      it_params.rep_file = nullptr;
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

static void make_origami_faces(const Geometry &geom, Geometry &orig,
                               map<vector<int>, vector<double>> &lens,
                               double slant, double init_ht)
{
  orig.clear_all();
  for (const auto &vert : geom.verts())
    orig.add_vert(vert, Color(0));
  int f_start = orig.verts().size();
  for (unsigned int f = 0; f < geom.faces().size(); f++) {
    // Starting face centre is centroid dropped by the slant height
    Vec3d face_pt = geom.face_cent(f) + geom.face_norm(f).with_len(init_ht);
    orig.add_vert(face_pt, Color(2));
  }

  GeometryInfo info(geom);
  double total_slant = 0.0;
  map<vector<int>, int> e2v;
  for (const auto &e : info.get_impl_edges())
    e2v[e] = orig.add_vert(geom.edge_cent(e), Color(1));
  for (int f_idx = 0; f_idx < (int)geom.faces().size(); f_idx++) {
    int f_cent_idx = f_start + f_idx;
    for (int v = 0; v < (int)geom.faces(f_idx).size(); v++) {
      int v0 = geom.faces(f_idx, v);
      int v1 = geom.faces_mod(f_idx, v + 1);
      int v0_n = info.get_vert_cons()[v0].size();
      int v1_n = info.get_vert_cons()[v1].size();
      int e_cent_idx = e2v[make_edge(v0, v1)];
      orig.add_face({v0, e_cent_idx, f_cent_idx}, v0_n);
      orig.add_face({v1, e_cent_idx, f_cent_idx}, v1_n);
      lens[{v0, e_cent_idx, v1, f_cent_idx}] = {
          slant / tan(M_PI / v0_n),  // v0 to e_cent
          slant / tan(M_PI / v1_n),  // e_cent to v1
          slant / sin(M_PI / v1_n),  // v1 to f_cent
          slant / sin(M_PI / v0_n)}; // f_cent to v0
      total_slant += orig.edge_vec(e_cent_idx, f_cent_idx).len();
    }
  }
  orig.transform(Trans3d::scale(0.5 / (total_slant / (2 * e2v.size()))));
}

inline double adjust_edge_mid(Geometry &geom, int v0, int v1, int v2,
                              double len, double factor)
{
  // moving v1 on long edge v0 to v2
  const Vec3d edge_vec = geom.edge_vec(v2, v1);
  const double edge_len = edge_vec.len();
  double diff = edge_len - len;
  Vec3d u = geom.edge_vec(v2, v0).unit(); // long edge
  geom.verts(v1) = geom.verts(v2) + u * (len + diff * factor);
  return diff;
}

inline double adjust_edge(Geometry &geom, int v0, int v1, double len,
                          double factor, bool v1_fixed = false)
{
  const Vec3d edge_vec = geom.edge_vec(v0, v1);
  const double edge_len = edge_vec.len();
  double diff = edge_len - len;
  Vec3d u = edge_vec / edge_len;
  geom.verts(v0) += u * diff * factor;
  if (!v1_fixed)
    geom.verts(v1) -= u * diff * factor;
  return diff;
}

Status make_origami(const Geometry &geom, Geometry &orig, iter_params it_params,
                    double factor, double init_ht)
{
  double max_diff = 0;
  double slant = 0.5; // hardcoded for a model of reasonable and known size
  map<vector<int>, vector<double>> lens;
  make_origami_faces(geom, orig, lens, slant, init_ht);
  for (int cnt = 1; cnt <= it_params.num_iters; cnt++) {
    double fact = factor;
    // Start gently seems like good idea, but is commented out as it affects
    // final symmetry
    // if(cnt < 100)
    //  fact /= 100;
    // else if(cnt < 200)
    //  fact /= 10;

    max_diff = 0.0;
    for (const auto &elem : lens) {
      const auto &elem_vs = elem.first;
      const auto &elem_lens = elem.second;
      // adjust slant height, adjust face centre vertex only
      adjust_edge(orig, elem_vs[3], elem_vs[1], slant, fact, true);
      for (int i = 0; i < 4; i++) {
        double diff = 0.0;
        if (i != 1)
          diff = adjust_edge(orig, elem_vs[i], elem_vs[(i + 1) % 4],
                             elem_lens[i], fact);
        else
          diff = adjust_edge_mid(orig, elem_vs[0], elem_vs[1], elem_vs[2],
                                 elem_lens[1], fact);
        if (diff > max_diff)
          max_diff = diff;
      }
    }

    if (!it_params.quiet() && it_params.check_status(cnt))
      fprintf(it_params.rep_file, "\niter:%-15d max_diff:%17.15f ", cnt,
              max_diff);
    else if (it_params.print_progress_dot(cnt))
      fprintf(it_params.rep_file, ".");
  }
  if (!it_params.quiet() && it_params.checking_status())
    fprintf(it_params.rep_file, "\nFinal:\niter:%-15d max:%17.15f\n",
            it_params.num_iters, max_diff);
  return Status::ok();
}

void truncate_faces(Geometry &orig, double trunc_len)
{
  map< vector<int>, int > e2v; // edge to truncation vertex
  for (int i = 0; i < (int)orig.faces().size(); i++) {
    auto face = orig.faces(i);
    // Truncate the vertex at F.
    int N = orig.colors(FACES).get(i).get_index();
    double ang = M_PI / N;  // Original vertex angle
    double ang_b = ang - atan(2 * trunc_len * tan(ang));
    double len2 = sqrt(pow(trunc_len, 2) + pow(0.5 / tan(ang), 2)) * cos(ang_b);

    // Avoid duplicates. Vertices may be coincident so cannot merge
    int idx12;
    auto mi = e2v.find(make_edge(face[1], face[2]));
    if (mi == e2v.end()) {
      idx12 =
          orig.add_vert(orig.verts(face[1]) +
                        orig.edge_vec(face[1], face[2]).with_len(trunc_len));
      e2v[make_edge(face[1], face[2])] = idx12;
    }
    else
      idx12 = mi->second;

    int idx02;
    mi = e2v.find(make_edge(face[0], face[2]));
    if (mi == e2v.end()) {
      idx02 = orig.add_vert(orig.verts(face[0]) +
                            orig.edge_vec(face[0], face[2]).with_len(len2));
      e2v[make_edge(face[0], face[2])] = idx02;
    }
    else
      idx02 = mi->second;

    orig.faces(i) = {face[0], face[1], idx12, idx02};
  }
  orig.del(VERTS, orig.get_info().get_free_verts()); // delete F vertices
}


int main(int argc, char *argv[])
{
  mmop_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);

  if (!geom.faces().size())
    opts.error("input file contains no faces");

  if(!opts.keep_orient && !geom.orient(1)) // positive orientation
      opts.error("base polyhedron cannot be oriented: override with option -k");

  Geometry origami;
  make_origami(geom, origami, opts.it_params, opts.adjust_fact / 200,
               opts.init_ht);
  if (opts.trunc_len != 1)
    truncate_faces(origami, opts.trunc_len / 2);

  GeometryInfo info(origami);
  for(const auto &face : origami.faces()) {
    int hub_idx = face[0];
    origami.colors(VERTS).set(hub_idx, Color::invisible);
    for(int con : info.get_vert_cons()[hub_idx])
      origami.add_edge(hub_idx, con, Color::invisible);
  }

  if(opts.color_method == 'v')
    for (unsigned int i = 0; i < origami.faces().size(); i++)
      origami.colors(FACES).set(i, geom.colors(VERTS).get(origami.faces(i, 0)));

  for (int i = 0; i < 3; i++)
      opts.clrngs[i].set_geom(&origami);

  opts.clrngs[0].v_apply_cmap();
  opts.clrngs[1].e_apply_cmap();
  opts.clrngs[2].f_apply_cmap();

  opts.write_or_error(origami, opts.ofile);

  return 0;
}
