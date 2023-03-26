/*
   Copyright (c) 2014-2023, Roger Kaufman

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
   Name: obj2off.cc
   Description: Convert files in OBJ format to OFF file format
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"

#define TINYOBJLOADER_IMPLEMENTATION
#include "./tiny_obj_loader.h"

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using std::string;
using std::vector;

using namespace anti;

class obj2off_opts : public ProgramOpts {
public:
  string ifile;
  string ofile;

  int sig_digits = DEF_SIG_DGTS; // significant digits output (system default)

  obj2off_opts() : ProgramOpts("obj2off") {}

  void process_command_line(int argc, char **argv);
  void usage();
};

void obj2off_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] [input_file]

Convert files in OBJ format to OFF format. Only v, e and l statements are
used. Face colors are take from diffuse coloring (Kd) in the mtl file. OBJ does
not support edge colors. If input_file is not given the program reads from
standard input.

Options
%s
  -d <dgts> number of significant digits (default %d) or if negative
            then the number of digits after the decimal point
  -o <file> write output to file (default: write to standard output)

)",
          prog_name(), help_ver_text, DEF_SIG_DGTS);
}

void obj2off_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hd:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {

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

  if (argc - optind > 1)
    error("too many arguments");

  if (argc - optind == 1)
    ifile = argv[optind];
}

// convert obj file to off  https://github.com/tinyobjloader/tinyobjloader
bool convert_obj_to_off(string &file_name, Geometry &geom,
                        const obj2off_opts &opts)
{
  // required parameters
  const char *filename = file_name.c_str();
  const char *basepath = NULL;
  bool triangulate = false;

  tinyobj::attrib_t attrib;
  std::vector<tinyobj::shape_t> shapes;
  std::vector<tinyobj::material_t> materials;

  std::string warn;
  std::string err;
  bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err,
                              filename, basepath, triangulate);

  if (!warn.empty())
    opts.warning(warn);

  if (!err.empty())
    opts.error(err);

  if (!ret) {
    opts.warning(msg_str("Failed to load/parse %s", filename));
    return false;
  }

  // vertices
  for (size_t v = 0; v < attrib.vertices.size() / 3; v++) {
    Vec3d vert = Vec3d(attrib.vertices[3 * v + 0], attrib.vertices[3 * v + 1],
                       attrib.vertices[3 * v + 2]);
    Color c = Color(attrib.colors[3 * v + 0], attrib.colors[3 * v + 1],
                    attrib.colors[3 * v + 2]);
    geom.add_vert(vert, c);
  }

  // for each shape
  for (size_t i = 0; i < shapes.size(); i++) {
    // for each face
    size_t index_offset = 0;
    for (size_t f = 0; f < shapes[i].mesh.num_face_vertices.size(); f++) {
      size_t fnum = shapes[i].mesh.num_face_vertices[f];
      // For each vertex in the face
      vector<int> face;
      for (size_t v = 0; v < fnum; v++) {
        tinyobj::index_t idx = shapes[i].mesh.indices[index_offset + v];
        face.push_back(idx.vertex_index);
      }
      // face color
      int material_id = shapes[i].mesh.material_ids[f];
      Color c = Color();
      if (material_id > -1)
        c = Color(materials[material_id].diffuse[0],
                  materials[material_id].diffuse[1],
                  materials[material_id].diffuse[2]);

      geom.add_face(face, c);

      index_offset += fnum;
    }

    // for each edge
    index_offset = 0;
    for (size_t e = 0; e < shapes[i].lines.num_line_vertices.size(); e++) {
      size_t eno = shapes[i].lines.num_line_vertices[e];
      //  For each vertex in the edge
      vector<int> edge;
      for (size_t v = 0; v < eno; v++) {
        tinyobj::index_t idx = shapes[i].lines.indices[index_offset + v];
        edge.push_back(idx.vertex_index);
      }

      geom.add_edge(make_edge(edge[0], edge[1]));

      index_offset += eno;
    }
  }

  return true;
}

int main(int argc, char *argv[])
{
  obj2off_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  convert_obj_to_off(opts.ifile, geom, opts);

  opts.write_or_error(geom, opts.ofile, opts.sig_digits);

  return 0;
}
