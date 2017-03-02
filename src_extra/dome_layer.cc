/*
   Copyright (c) 2009-2016, Adrian Rossiter

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
   Name: dome_layer.cc
   Description: add layers to a single layer dome
   Project: Antiprism - http://www.antiprism.com
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <algorithm>
#include <string>
#include <vector>

#include "../base/antiprism.h"

using std::string;
using std::vector;
using std::map;
using std::swap;

using namespace anti;

class dome_opts : public ProgramOpts {
public:
  string type;
  double radius;
  bool use_Colorues;
  string ifile;
  string ofile;

  dome_opts()
      : ProgramOpts("dome_layer"), type("eden"), radius(0.85),
        use_Colorues(true)
  {
  }
  void process_command_line(int argc, char **argv);
  void usage();
};

// clang-format off
void dome_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Read a dome polyhedron in OFF format, project onto unit sphere and add\n"
"additional layering. If input_file is not given the program reads from\n"
"standard input.\n"
"\n"
"Options\n"
"%s"
"  -r <rad>  radius of sphere for second layer vertices (default: 0.85)\n"
"  -t <type> layer type (default: dual)\n"
"               dual:      base vertices to dual vertices\n"
"               eden:      base vertices to edge centre vertices\n"
"               asm:       base vertices to inner base, vertices to cell centres\n"
"               honeycomb: base vertices to dual of 'kis' form\n"
"               prism:     repeat base and connect faces as prism tensegrities\n"
"  -i        write colours as index numbers (struts 0-4, faces 10)\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text);
}
// clang-format on

void dome_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;
  char name[MSG_SZ];

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hr:t:io:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'o':
      ofile = optarg;
      break;

    case 't':
      type = to_resource_name(name, optarg);
      if (type != "eden" && type != "dual" && type != "asm" &&
          type != "honeycomb" && type != "prism")
        error("unknown layering type", c);
      break;

    case 'r':
      print_status_or_exit(read_double(optarg, &radius), c);
      break;

    case 'i':
      use_Colorues = false;
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

void proj_onto_sphere(Geometry &geom, double radius = 1.0)
{
  for (unsigned int i = 0; i < geom.verts().size(); i++) {
    geom.raw_verts()[i].to_unit();
    geom.raw_verts()[i] *= radius;
  }
}

int get_edge_idx(const Geometry &geom, int v0, int v1)
{
  vector<int> edge(2);
  edge[0] = v0;
  edge[1] = v1;
  if (edge[0] > edge[1])
    std::swap(edge[0], edge[1]);
  vector<vector<int>>::const_iterator ei =
      find(geom.edges().begin(), geom.edges().end(), edge);
  return (ei != geom.edges().end()) ? ei - geom.edges().begin() : -1;
}

void make_dome_eden(const Geometry &geom, Geometry &dome, double radius)
{
  dome = geom;
  GeometryInfo info(dome);

  /*
   vector<vector<int> > v_cons(info.num_verts());
   for(int i=0; i<info.num_verts(); i++)
      if(info.get_vert_figs()[i].size())
         v_cons[i] = info.get_vert_figs()[i][0];
   */

  proj_onto_sphere(dome);

  dome.add_missing_impl_edges();
  Coloring d_clrng(&dome);
  d_clrng.f_one_col(Color(10));
  d_clrng.e_one_col(Color(0));

  int orig_num_faces = dome.faces().size();
  int orig_num_verts = dome.verts().size();
  for (unsigned int i = 0; i < dome.edges().size(); i++)
    dome.add_vert(dome.edge_cent(i).unit() * radius);

  for (int i = 0; i < orig_num_faces; i++) {
    int f_sz = dome.faces(i).size();
    vector<int> face(f_sz);
    for (int j = 0; j < f_sz; j++)
      face[j] = orig_num_verts + get_edge_idx(dome, dome.faces(i, j),
                                              dome.faces(i, (j + 1) % f_sz));
    dome.add_face(face, Color(10));
    for (int j = 0; j < f_sz; j++)
      dome.add_edge(face[j], face[(j + 1) % f_sz], Color(1));
  }

  for (int i = 0; i < orig_num_verts; i++) {
    const vector<vector<int>> &v_fig = info.get_vert_figs()[i];
    for (unsigned int j = 0; j < v_fig.size(); j++) {
      const vector<int> &cons = v_fig[j];
      int f_sz = cons.size();
      vector<int> face(f_sz);
      for (int k = 0; k < f_sz; k++)
        face[k] = orig_num_verts + get_edge_idx(dome, i, cons[k]);
      dome.add_face(face, Color(10));
      for (int k = 0; k < f_sz; k++)
        dome.add_edge(i, face[k], Color(2));
    }
  }
}

void make_dome_dual(const Geometry &geom, Geometry &dome, double radius)
{
  dome = geom;
  proj_onto_sphere(dome);

  Geometry dual;
  get_dual(&dual, dome);
  dome.face_cents(dual.raw_verts());
  proj_onto_sphere(dual, radius);

  dome.add_missing_impl_edges();
  Coloring d_clrng(&dome);
  d_clrng.f_one_col(Color(10));
  d_clrng.e_one_col(Color(0));

  dual.add_missing_impl_edges();
  Coloring dl_clrng(&dual);
  dl_clrng.f_one_col(Color(10));
  dl_clrng.e_one_col(Color(1));

  int orig_num_faces = dome.faces().size();
  int orig_num_verts = dome.verts().size();
  dome.append(dual);
  for (int i = 0; i < orig_num_faces; i++)
    for (unsigned int j = 0; j < dome.faces(i).size(); j++)
      dome.add_edge(dome.faces(i, j), orig_num_verts + i, Color(2));
}

void make_dome_asm(const Geometry &geom, Geometry &dome, double radius)
{
  dome = geom;
  proj_onto_sphere(dome);

  Geometry dual;
  get_dual(&dual, dome);
  dome.face_cents(dual.raw_verts());
  proj_onto_sphere(dual, sqrt(radius));
  dual.clear(FACES);

  Geometry inner_dome = geom;
  proj_onto_sphere(inner_dome, radius);

  dome.add_missing_impl_edges();
  Coloring d_clrng(&dome);
  d_clrng.f_one_col(Color(10));
  d_clrng.e_one_col(Color(0));

  inner_dome.add_missing_impl_edges();
  Coloring in_clrng(&inner_dome);
  in_clrng.f_one_col(Color(10));
  in_clrng.e_one_col(Color(1));

  int orig_num_faces = dome.faces().size();
  int orig_num_verts = dome.verts().size();
  dome.append(dual);
  dome.append(inner_dome);
  for (int i = 0; i < orig_num_faces; i++)
    for (unsigned int j = 0; j < dome.faces(i).size(); j++) {
      dome.add_edge(dome.faces(i, j), orig_num_verts + i, Color(2));
      dome.add_edge(dome.faces(i, j) + orig_num_verts + orig_num_faces,
                    orig_num_verts + i, Color(2));
    }

  for (int i = 0; i < orig_num_verts; i++)
    dome.add_edge(i, orig_num_verts + orig_num_faces + i, Color(3));
}

void make_dome_honeycomb(const Geometry &geom, Geometry &dome, double radius)
{
  dome = geom;
  proj_onto_sphere(dome);

  Geometry dual;
  get_dual(&dual, dome);
  dome.face_cents(dual.raw_verts());
  proj_onto_sphere(dual);
  truncate_verts(&dual, 1 / 3.0);
  proj_onto_sphere(dual, radius);

  dome.add_missing_impl_edges();
  Coloring d_clrng(&dome);
  d_clrng.f_one_col(Color(10));
  d_clrng.e_one_col(Color(0));

  dual.add_missing_impl_edges();
  Coloring dl_clrng(&dual);
  dl_clrng.f_one_col(Color(10));
  dl_clrng.e_one_col(Color(1));

  int orig_num_faces = dome.faces().size();
  int orig_num_verts = dome.verts().size();
  dome.append(dual);
  for (int i = 0; i < orig_num_verts; i++) {
    int f_idx = i + orig_num_faces;
    for (unsigned int j = 0; j < dome.faces(f_idx).size(); j++)
      dome.add_edge(dome.faces(f_idx, j), i, Color(2));
  }
}

void make_dome_prism(const Geometry &geom, Geometry &dome, double radius)
{
  for (unsigned int f_idx = 0; f_idx < geom.faces().size(); f_idx++) {
    int orig_vsz = dome.verts().size();
    const vector<int> &face = geom.faces(f_idx);
    const int fsz = face.size();
    for (unsigned int i = 0; i < face.size(); i++) {
      dome.add_vert(
          geom.edge_cent(make_edge(face[i], face[(i + 1) % fsz])).unit());
      dome.add_vert(geom.edge_cent(make_edge(face[(i + 2) % fsz],
                                             face[(i + 2 + 1) % fsz]))
                        .unit() *
                    radius);
      int vsz = dome.verts().size();
      dome.add_edge(vsz - 2, vsz - 1, 0);

      dome.add_edge(orig_vsz + (2 * i), orig_vsz + (2 * i + 2) % (2 * fsz), 1);
      dome.add_edge(orig_vsz + (2 * i + 1), orig_vsz + (2 * i + 3) % (2 * fsz),
                    2);
      dome.add_edge(orig_vsz + (2 * i),
                    orig_vsz + (2 * i + 2 * fsz - 1) % (2 * fsz), 3);
    }
  }
}

void set_color_values(Geometry &geom)
{
  ColorMapMap *cmap = new ColorMapMap;
  cmap->set_col(0, Color(1.0, 0.0, 0.0));
  cmap->set_col(1, Color(0.0, 0.0, 1.0));
  cmap->set_col(2, Color(0.0, 1.0, 0.0));
  cmap->set_col(3, Color(0.9, 0.9, 0.0));
  cmap->set_col(10, Color(1.0, 1.0, 1.0, 0.3));
  Coloring clrng(&geom);
  clrng.add_cmap(cmap);
  clrng.e_apply_cmap();
  clrng.f_apply_cmap();
}

int main(int argc, char *argv[])
{
  dome_opts opts;
  opts.process_command_line(argc, argv);
  Geometry geom;
  opts.read_or_error(geom, opts.ifile);

  Geometry dome;
  if (opts.type == "eden")
    make_dome_eden(geom, dome, opts.radius);
  else if (opts.type == "dual")
    make_dome_dual(geom, dome, opts.radius);
  else if (opts.type == "asm")
    make_dome_asm(geom, dome, opts.radius);
  else if (opts.type == "honeycomb")
    make_dome_honeycomb(geom, dome, opts.radius);
  else if (opts.type == "prism")
    make_dome_prism(geom, dome, opts.radius);

  if (opts.use_Colorues)
    set_color_values(dome);

  opts.write_or_error(dome, opts.ofile);

  return 0;
}
