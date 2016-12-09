/*
   Copyright (c) 2007-2016, Roger Kaufman

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
   Name: off2m.cc
   Description: Convert files in OFF format to 'm' format for LiveGraphics3D
   Project: Antiprism - http://www.antiprism.com
*/

#include <stdio.h>
#include <stdlib.h>

#include <ctype.h>

#include <string>
#include <vector>

#include "../base/antiprism.h"

using std::string;
using std::vector;

using namespace anti;

class o2m_opts : public ProgramOpts {
public:
  string ifile;
  string ofile;

  double edge_size;
  double vert_size;
  bool lighting;
  string exclude_elems;
  Color face_col;
  Color edge_col;
  Color vert_col;
  Color bg_col;
  int sig_digits;
  int f_dtype;
  Vec3d view_point;

  o2m_opts()
      : ProgramOpts("off2m"), edge_size(0.0), vert_size(0.0), lighting(false),
        exclude_elems(""),
        // face colors need to be explictly set since LG3D defaults to black
        face_col(Color(0.8, 0.9, 0.9)), // the system default colors
        edge_col(Color(0.8, 0.6, 0.8)), vert_col(Color(1.0, 0.5, 0.0)),
        bg_col(Color(0.9, 0.9, 0.9)), sig_digits(DEF_SIG_DGTS), f_dtype(1),
        view_point(Vec3d(0, 0, 0))
  {
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

void o2m_opts::usage()
{
  fprintf(
      stdout,
      "\n"
      "Usage: %s [options] [input_file]\n"
      "\n"
      "Convert files in OFF format to 'm' format for display in "
      "LiveGraphics3D. If\n"
      "input_file is not given the program reads from standard input.\n"
      "\n"
      "Options\n"
      "%s"
      "  -v <size> vertex sphere size (default: 0.02 of bounding box "
      "diagonal)\n"
      "  -e <size> frame model edge thickness size (default: 0.01 of bounding\n"
      "            box diagonal if vertex sphere is 0, else vertex_rad/1.5)\n"
      "  -x <elms> hide elements. The element string can include v, e and f\n"
      "            to hide vertices, edges and faces\n"
      "  -o <file> write output to file (default: write to standard output)\n"
      "\n"
      "\nColoring Options (run 'off_util -H color' for help on color formats)\n"
      "Note: transparency (alpha) is ignored. color name \"invisible\" not "
      "allowed\n"
      "\n"
      "  -V <col>  vertex color (default: 1.0,0.5,0.0)\n"
      "  -E <col>  edge color (default: 0.8,0.6,0.8)\n"
      "  -F <col>  face color (default: 0.8,0.8,0.9)\n"
      "  -l        let LiveGraphics3D do the Coloring itself\n"
      "\n"
      "Scene options\n"
      "  -Y <view> specify the Live3D ViewPoint in form 'X,Y,Z'\n"
      "  -B <col>  background color (default: 0.9,0.9,0.9)\n"
      "\n"
      "Precision options\n"
      "  -d <dgts> number of significant digits (default %d) or if negative\n"
      "            then the number of digits after the decimal point\n"
      "  -t <type> display type for faces 0 - polygons, 1 - triangulate\n"
      "               polygons (default)\n"
      "\n"
      "\n",
      prog_name(), help_ver_text, DEF_SIG_DGTS);
}

void o2m_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":he:v:x:lF:E:V:d:t:o:Y:B:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'e':
      print_status_or_exit(read_double(optarg, &edge_size), c);
      if (edge_size <= 0)
        error("frame width cannot be zero or negative", c);
      break;

    case 'v':
      print_status_or_exit(read_double(optarg, &vert_size), c);
      if (vert_size < 0)
        error("point size cannot be zero or negative", c);
      break;

    case 'l':
      lighting = true;
      break;

    case 'x':
      if (strspn(optarg, "vef") != strlen(optarg))
        error(msg_str("elements to exclude are '%s' must be from "
                      "v, e, and f",
                      optarg),
              c);
      exclude_elems = optarg;
      break;

    case 'F':
      print_status_or_exit(face_col.read(optarg), c);
      if (face_col.is_inv())
        error("face color may not be invisible", c);
      break;

    case 'E':
      print_status_or_exit(edge_col.read(optarg), c);
      if (edge_col.is_inv())
        error("edge color may not be invisible", c);
      break;

    case 'V':
      print_status_or_exit(vert_col.read(optarg), c);
      if (vert_col.is_inv())
        error("vert color may not be invisible", c);
      break;

    case 'd':
      print_status_or_exit(read_int(optarg, &sig_digits), c);
      break;

    case 't':
      print_status_or_exit(read_int(optarg, &f_dtype), c);
      if (f_dtype != 0 && f_dtype != 1)
        error("display type for faces must be 0 or 1", c);
      break;

    case 'o':
      ofile = optarg;
      break;

    case 'Y':
      print_status_or_exit(view_point.read(optarg), c);
      break;

    case 'B':
      print_status_or_exit(bg_col.read(optarg), c);
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

string RGBtxt(const Color &col)
{
  char buf[128];
  Vec3d cv = col.get_Vec3d();
  snprintf(buf, 128, "RGBColor[%g, %g, %g]", cv[0], cv[1], cv[2]);
  return buf;
}

string Vtxt(const Vec3d &v, const int &dgts)
{
  char buf[128];
  if (dgts > 0)
    snprintf(buf, 128, "{%.*g, %.*g, %.*g}", dgts, v[0], dgts, v[1], dgts,
             v[2]);
  else
    snprintf(buf, 128, "{%.*f, %.*f, %.*f}", -dgts, v[0], -dgts, v[1], -dgts,
             v[2]);
  return buf;
}

// elem_type 1=vertex, 2=edge, 3=face
int get_last_visible(const Geometry &geom, const int &elem_type,
                     const Color &def_col)
{
  int end = (elem_type == 1 ? geom.verts().size()
                            : (elem_type == 2 ? geom.edges().size()
                                              : geom.faces().size()));

  int last = -1;
  for (int i = end - 1; i >= 0; i--) {
    Color col = (elem_type == 1 ? geom.colors(VERTS).get(i)
                                : (elem_type == 2 ? geom.colors(EDGES).get(i)
                                                  : geom.colors(FACES).get(i)));

    col = col.is_val() ? col : def_col;
    if (!col.is_inv()) {
      last = i;
      break;
    }
  }

  return last;
}

int print_m_solid(FILE *ofile, const Geometry &geom, const int &sig_digits,
                  const Color &face_col)
{
  const vector<vector<int>> &faces = geom.faces();
  const vector<Vec3d> &verts = geom.verts();

  int last = get_last_visible(geom, 3, face_col);
  if (last == -1)
    return 0;

  for (int i = 0; i <= last; i++) {
    fprintf(ofile, "{");

    Color fcol = geom.colors(FACES).get(i);
    fcol = fcol.is_val() ? fcol : face_col;
    if (fcol.is_inv())
      continue;
    string RGB = RGBtxt(fcol);
    fprintf(ofile, "FaceForm[ %s, %s ],\n", RGB.c_str(), RGB.c_str());

    fprintf(ofile, "EdgeForm[], Polygon[{");
    for (unsigned int j = 0; j < faces[i].size(); j++) {
      fprintf(ofile, "%s", Vtxt(verts[faces[i][j]], sig_digits).c_str());
      if ((j + 1) < faces[i].size())
        fprintf(ofile, ", ");
    }

    fprintf(ofile, "}]}");
    if (i == last)
      fprintf(ofile, "}");
    fprintf(ofile, ",\n");
  }

  return 1;
}

void print_m_frame_edge(FILE *ofile, const Vec3d &v1, const Vec3d &v2,
                        const int &sig_digits, const double &edge_size,
                        const Color &edge_col)
{
  fprintf(ofile, "{");
  fprintf(ofile, "Thickness[%g], ", edge_size);
  if (edge_col.is_set())
    fprintf(ofile, "%s, ", RGBtxt(edge_col).c_str());
  fprintf(ofile, "Line[{%s, %s}", Vtxt(v1, sig_digits).c_str(),
          Vtxt(v2, sig_digits).c_str());
  fprintf(ofile, "]}");
}

int print_m_frame(FILE *ofile, const Geometry &geom, const int &sig_digits,
                  const double &edge_size, const Color &edge_col,
                  const bool &more_to_print)
{
  const vector<vector<int>> &edges = geom.edges();
  const vector<Vec3d> &verts = geom.verts();

  int last = get_last_visible(geom, 2, edge_col);
  if (last == -1)
    return 0;

  for (int i = 0; i <= last; i++) {
    Color ecol = geom.colors(EDGES).get(i);
    ecol = ecol.is_val() ? ecol : edge_col;
    if (ecol.is_inv())
      continue;
    print_m_frame_edge(ofile, verts[edges[i][0]], verts[edges[i][1]],
                       sig_digits, edge_size, ecol);
    if (i == last && !more_to_print)
      fprintf(ofile, "}");
    fprintf(ofile, ",\n");
  }

  return 1;
}

int print_m_points(FILE *ofile, const Geometry &geom, const int &sig_digits,
                   const double &vert_size, const Color &vert_col,
                   const bool &more_to_print)
{
  const vector<Vec3d> &verts = geom.verts();

  int last = get_last_visible(geom, 1, vert_col);
  if (last == -1)
    return 0;

  for (int i = 0; i <= last; i++) {
    Color vcol = geom.colors(VERTS).get(i);
    vcol = vcol.is_val() ? vcol : vert_col;
    if (vcol.is_inv())
      continue;
    fprintf(ofile, "{");
    if (vcol.is_set())
      fprintf(ofile, "%s, ", RGBtxt(vcol).c_str());
    fprintf(ofile, "PointSize[%g], ", vert_size);
    fprintf(ofile, "Point[%s]}", Vtxt(verts[i], sig_digits).c_str());

    if (i == last && !more_to_print)
      fprintf(ofile, "}");
    fprintf(ofile, ",\n");
  }

  return 1;
}

void print_m_head(FILE *ofile) { fprintf(ofile, "Graphics3D[{\n"); }

// view_point is not changed
void print_m_tail(FILE *ofile, const Geometry &geom, const bool &lighting,
                  const Color &bg, Vec3d view_point)
{
  BoundSphere b_sph(geom.verts());
  double radius = b_sph.get_radius();

  // If ViewPoint isn't specified, Live3D sets an internal one of 1.3,-2.4,2 so
  // the model would be tilted
  // ViewPoint(0,0,0) makes the model disappear in LG3D so don't allow it
  // The default is to let the Focal Length be normal such that
  // ViewPoint(0,0,radius*100)
  // to set Z back far enough to avoid a "fish eye" view
  if (view_point[0] == 0 && view_point[1] == 0 && view_point[2] == 0)
    view_point[2] = radius * 100;

  fprintf(ofile, "ViewPoint -> {%g,%g,%g}, ", view_point[0], view_point[1],
          view_point[2]);

  // Setting ViewVertical with Y upright makes them model appear the same as in
  // AntiView, off2pov, etc
  fprintf(ofile, "ViewVertical -> {0.0,1.0,0.0},\n");

  fprintf(ofile, "Background -> %s, ", RGBtxt(bg).c_str());
  if (lighting)
    fprintf(ofile, "Lighting->True, ");
  else
    fprintf(ofile, "Lighting->False, ");
  fprintf(ofile, "Boxed->False]");
}

int main(int argc, char *argv[])
{
  o2m_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  opts.read_or_error(geom, opts.ifile);

  if (!opts.edge_col.is_inv())
    geom.add_missing_impl_edges();

  if (opts.f_dtype == 1)
    geom.triangulate(Color::invisible);

  BoundBox bbox(geom.verts());
  double to_model_units = 2.0 / (bbox.get_max() - bbox.get_min()).len();
  if (opts.vert_size == 0.0)
    opts.vert_size = 0.02 / to_model_units;
  if (opts.edge_size == 0.0) {
    if (opts.vert_size == 0.0)
      opts.edge_size = 0.01 / to_model_units;
    else
      opts.edge_size = opts.vert_size / 1.5;
  }

  FILE *ofile = stdout; // write to stdout by default
  if (opts.ofile != "") {
    ofile = fopen(opts.ofile.c_str(), "w");
    if (ofile == 0)
      opts.error("could not open output file \'" + opts.ofile + "\'");
  }

  bool print_v =
      !strchr(opts.exclude_elems.c_str(), 'v') && geom.verts().size();
  bool print_e =
      !strchr(opts.exclude_elems.c_str(), 'e') && geom.edges().size();
  bool print_f =
      !strchr(opts.exclude_elems.c_str(), 'f') && geom.faces().size();

  if (!print_v && !print_e && !print_f)
    opts.error("there are no elements to output", 'x');

  print_m_head(ofile);
  if (print_v) {
    if (!print_m_points(ofile, geom, opts.sig_digits,
                        opts.vert_size * to_model_units, opts.vert_col,
                        (print_e || print_f)))
      opts.error("all vertices are invisible. try excluding by -x v");
  }
  if (print_e) {
    if (!print_m_frame(ofile, geom, opts.sig_digits,
                       opts.edge_size * to_model_units, opts.edge_col, print_f))
      opts.error("all edges are invisible. try excluding by -x e");
  }
  if (print_f) {
    if (!print_m_solid(ofile, geom, opts.sig_digits, opts.face_col))
      opts.error("all faces are invisible. try excluding by -x f");
  }
  print_m_tail(ofile, geom, opts.lighting, opts.bg_col, opts.view_point);

  if (opts.ofile != "")
    fclose(ofile);

  return 0;
}
