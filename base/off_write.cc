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

/* \file off_file.cc
   \brief Write OFF files
*/

#include "private_off_file.h"
#include "utils.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

using std::map;
using std::string;
using std::vector;

FILE *file_open_w(string file_name, string &error_msg)
{
  error_msg.clear();
  FILE *ofile = stdout; // write to stdout by default
  if (file_name != "") {
    ofile = fopen(file_name.c_str(), "w");
    if (!ofile)
      error_msg = "could not open file for writing '" + file_name +
                  "': " + strerror(errno);
  }
  return ofile;
}

void file_close_w(FILE *ofile)
{
  if (ofile != stdout)
    fclose(ofile);
}

void crds_write(FILE *ofile, const Geometry &geom, const char *sep,
                int sig_dgts)
{
  for (unsigned int i = 0; i < geom.verts().size(); i++)
    fprintf(ofile, "%s\n", geom.verts(i).to_str(sep, sig_dgts).c_str());
}

Status crds_write(string file_name, const Geometry &geom, const char *sep,
                  int sig_dgts)
{
  string error_msg;
  FILE *ofile = file_open_w(file_name, error_msg);
  if (!ofile)
    return Status::error(error_msg);

  crds_write(ofile, geom, sep, sig_dgts);
  file_close_w(ofile);
  return Status::ok();
}

// RK - color sorting functions
bool cmp_col(const Color &a, const Color &b)
{
  bool ret = false;

  Vec4d hsva_a = a.get_hsva();
  Vec4d hsva_b = b.get_hsva();
  for (unsigned int i = 0; i < 4; i++) {
    if (hsva_a[i] != hsva_b[i]) {
      ret = (hsva_a[i] < hsva_b[i]);
      break;
    }
  }

  return ret;
}

class col_cmp {
public:
  col_cmp() = default;
  bool operator()(const Color &a, const Color &b) const
  {
    return cmp_col(a, b);
  }
};

string off_col(Color col)
{
  if (col.is_index())
    return msg_str(" %d", col.get_index());
  else if (col.is_value()) {
    if (col.get_transparency())
      return col.get_vec4d().to_str(" ", -5);
    else
      return col.get_vec3d().to_str(" ", -5);
  }
  else
    return string();
}

void write_mtl_color(FILE *mfile, Color c)
{
  fprintf(mfile, "newmtl color_%02x%02x%02x%02x\n", c[0], c[1], c[2], c[3]);
  fprintf(mfile, "Kd %g %g %g\n", c[0] / 255.0, c[1] / 255.0, c[2] / 255.0);
  fprintf(mfile, "illum 2\n");
  fprintf(mfile, "\n");
}

void write_mtl_file(FILE *mfile, vector<Color> &cols)
{
  // add default colors to mtl file
  fprintf(mfile, "newmtl color_vert_default\n");
  fprintf(mfile, "Kd 1.0 0.5 0.0\n");
  fprintf(mfile, "illum 2\n");
  fprintf(mfile, "\n");
  fprintf(mfile, "newmtl color_edge_default\n");
  fprintf(mfile, "Kd 0.8 0.6 0.8\n");
  fprintf(mfile, "illum 2\n");
  fprintf(mfile, "\n");
  fprintf(mfile, "newmtl color_face_default\n");
  fprintf(mfile, "Kd 0.8 0.9 0.9\n");
  fprintf(mfile, "illum 2\n");
  fprintf(mfile, "\n");

  // sort and find unique colors
  sort(cols.begin(), cols.end(), col_cmp());
  auto ci = unique(cols.begin(), cols.end());
  cols.resize(ci - cols.begin());

  for (auto &col : cols)
    write_mtl_color(mfile, col);
}

// RK - write OBJ file type for Meshlab and other
void obj_write(FILE *ofile, FILE *mfile, string mtl_file, const Geometry &geom,
               const char *sep, int sig_dgts)
{
  int offset = 1; // obj files start indexes from 1

  vector<Color> cols;

  fprintf(ofile, "# File type: ASCII OBJ\n");

  // materials file reference as string
  if (mfile)
    fprintf(ofile, "mtllib %s\n", mtl_file.c_str());

  // v entries
  for (unsigned int i = 0; i < geom.verts().size(); i++) {
    Color c = geom.colors(VERTS).get(i);
    // only color values can be used
    if (!c.is_value())
      c = Color();
    else if (!c.is_invisible()) {
      c.set_alpha(255); // future transparency possible?
    }
    fprintf(ofile, "v %s %s\n", geom.verts(i).to_str(sep, sig_dgts).c_str(),
            off_col(c).c_str());
  }

  Color last_color = Color();

  // f entries
  for (unsigned int i = 0; i < geom.faces().size(); i++) {
    // if materials, color logic
    if (mfile) {
      Color c = geom.colors(FACES).get(i);
      if (c.is_value() && !c.is_invisible()) {
        c.set_alpha(255); // future transparency possible?
        cols.push_back(c);
      }
      // first color might be unset
      if (c != last_color || i == 0) {
        if (c.is_value() && !c.is_invisible())
          fprintf(ofile, "usemtl color_%02x%02x%02x%02x\n", c[0], c[1], c[2],
                  c[3]);
        else
          fprintf(ofile, "usemtl color_face_default\n");
      }
      last_color = c;
    }
    fprintf(ofile, "f");
    for (unsigned int j = 0; j < geom.faces(i).size(); j++)
      fprintf(ofile, " %d", geom.faces(i, j) + offset);
    fprintf(ofile, "\n");
  }

  last_color = Color();

  // l entries
  for (unsigned int i = 0; i < geom.edges().size(); i++) {
    /* edges cannot currently have colors
        // if materials, color logic
        if (mfile) {
          Color c = geom.colors(EDGES).get(i);
          if (c.is_value() && !c.is_invisible()) {
            c.set_alpha(255); // future transparency possible?
            cols.push_back(c);
          }
          // first color might be unset
          if (c != last_color || i == 0) {
            if (c.is_value() && !c.is_invisible())
              fprintf(ofile, "usemtl color_%02x%02x%02x%02x\n", c[0], c[1],
       c[2], c[3]); else fprintf(ofile, "usemtl color_edge_default\n");
          }
          last_color = c;
        }
    */
    fprintf(ofile, "l %d %d\n", geom.edges(i, 0) + offset,
            geom.edges(i, 1) + offset);
  }

  if (mfile)
    write_mtl_file(mfile, cols);
}

Status obj_write(string file_name, string mtl_file, const Geometry &geom,
                 const char *sep, int sig_dgts)
{
  string error_msg;
  FILE *ofile = file_open_w(file_name, error_msg);
  if (!ofile)
    return Status::error(error_msg);

  FILE *mfile = nullptr;
  if (mtl_file.length()) {
    mfile = file_open_w(mtl_file, error_msg);
    if (!mfile)
      return Status::error(error_msg);
  }

  obj_write(ofile, mfile, mtl_file, geom, sep, sig_dgts);
  file_close_w(ofile);

  if (mfile)
    file_close_w(mfile);
  return Status::ok();
}

Status off_file_write(string file_name, const Geometry &geom, int sig_dgts)
{
  vector<const Geometry *> vg;
  vg.push_back(&geom);
  return off_file_write(file_name, vg, sig_dgts);
}

Status off_file_write(string file_name, const vector<const Geometry *> &geoms,
                      int sig_dgts)
{
  string error_msg;
  FILE *ofile = file_open_w(file_name, error_msg);
  if (!ofile)
    return Status::error(error_msg);

  off_file_write(ofile, geoms, sig_dgts);
  file_close_w(ofile);
  return Status::ok();
}

void off_polys_write(FILE *ofile, const Geometry &geom, int offset)
{
  for (unsigned int i = 0; i < geom.faces().size(); i++) {
    fprintf(ofile, "%lu", (unsigned long)geom.faces(i).size());
    for (unsigned int j = 0; j < geom.faces(i).size(); j++) {
      fprintf(ofile, " %d", geom.faces(i, j) + offset);
    }
    fprintf(ofile, " %s", off_col(geom.colors(FACES).get(i)).c_str());
    fprintf(ofile, "\n");
  }

  for (unsigned int i = 0; i < geom.edges().size(); i++) {
    fprintf(ofile, "2 %d %d", geom.edges(i, 0) + offset,
            geom.edges(i, 1) + offset);
    fprintf(ofile, " %s", off_col(geom.colors(EDGES).get(i)).c_str());
    fprintf(ofile, "\n");
  }
  // print coloured vertex elements
  map<int, Color>::const_iterator mi;
  for (mi = geom.colors(VERTS).get_properties().begin();
       mi != geom.colors(VERTS).get_properties().end(); mi++) {
    fprintf(ofile, "1 %d %s\n", (mi->first) + offset,
            off_col(mi->second).c_str());
  }
}

void off_file_write(FILE *ofile, const vector<const Geometry *> &geoms,
                    int sig_dgts)
{
  int vert_cnt = 0, face_cnt = 0, edge_cnt = 0;
  for (auto geom : geoms) {
    int num_v_col_elems = geom->colors(VERTS).get_properties().size();
    vert_cnt += geom->verts().size();
    edge_cnt += geom->edges().size();
    face_cnt += geom->faces().size() + num_v_col_elems + edge_cnt;
  }

  fprintf(ofile, "OFF\n%d %d 0\n", vert_cnt, face_cnt);

  for (auto geom : geoms)
    crds_write(ofile, *geom, " ", sig_dgts);

  int last_offset = 0;
  vert_cnt = 0;
  for (auto geom : geoms) {
    off_polys_write(ofile, *geom,
                    geom->verts().size() ? vert_cnt : last_offset);
    last_offset = vert_cnt;
    vert_cnt += geom->verts().size();
  }
}

void off_file_write(FILE *ofile, const Geometry &geom, int sig_dgts)
{
  vector<const Geometry *> vg;
  vg.push_back(&geom);
  off_file_write(ofile, vg, sig_dgts);
}
