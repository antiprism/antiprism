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

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <vector>

#include "private_off_file.h"
#include "utils.h"

using std::map;
using std::string;
using std::vector;

FILE *file_open_w(string file_name, char *errmsg)
{
  if (errmsg)
    *errmsg = '\0';
  FILE *ofile = stdout; // write to stdout by default
  if (file_name != "") {
    ofile = fopen(file_name.c_str(), "w");
    if (!ofile && errmsg)
      snprintf(errmsg, MSG_SZ, "could not output file \'%s\'",
               file_name.c_str());
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
  char line[MSG_SZ];
  for (unsigned int i = 0; i < geom.verts().size(); i++)
    fprintf(ofile, "%s\n", vtostr(line, geom.verts(i), sep, sig_dgts));
}

bool crds_write(string file_name, const Geometry &geom, char *errmsg,
                const char *sep, int sig_dgts)
{
  FILE *ofile = file_open_w(file_name, errmsg);
  if (!ofile)
    return false;

  crds_write(ofile, geom, sep, sig_dgts);
  file_close_w(ofile);
  return true;
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
  char line[MSG_SZ];
  for (unsigned int i = 0; i < geom.verts().size(); i++)
    fprintf(ofile, "v %s\n", vtostr(line, geom.verts(i), sep, sig_dgts));

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
          fprintf(ofile, "usemtl color_%02x%02x%02x%02x\n", c[0], c[1], c[2],
                  c[3]);
        else
          fprintf(ofile, "usemtl color_edge_default\n");
      }
      last_color = c;
    }
    fprintf(ofile, "l %d %d\n", geom.edges(i, 0) + offset,
            geom.edges(i, 1) + offset);
  }

  last_color = Color();

  // p entries
  for (unsigned int i = 0; i < geom.verts().size(); i++) {
    // if materials, color logic
    if (mfile) {
      Color c = geom.colors(VERTS).get(i);
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
          fprintf(ofile, "usemtl color_vert_default\n");
      }
      last_color = c;
    }
    fprintf(ofile, "p %d\n", i + offset);
  }

  if (mfile)
    write_mtl_file(mfile, cols);
}

bool obj_write(string file_name, string mtl_file, const Geometry &geom,
               char *errmsg, const char *sep, int sig_dgts)
{
  FILE *ofile = file_open_w(file_name, errmsg);
  if (!ofile)
    return false;

  FILE *mfile = nullptr;
  if (mtl_file.length()) {
    mfile = file_open_w(mtl_file, errmsg);
    if (!mfile)
      return false;
  }

  obj_write(ofile, mfile, mtl_file, geom, sep, sig_dgts);
  file_close_w(ofile);

  if (mfile)
    file_close_w(mfile);
  return true;
}

bool off_file_write(string file_name, const Geometry &geom, char *errmsg,
                    int sig_dgts)
{
  vector<const Geometry *> vg;
  vg.push_back(&geom);
  return off_file_write(file_name, vg, errmsg, sig_dgts);
}

bool off_file_write(string file_name, const vector<const Geometry *> &geoms,
                    char *errmsg, int sig_dgts)
{
  if (errmsg)
    *errmsg = '\0';
  FILE *ofile = file_open_w(file_name, errmsg);
  if (!ofile)
    return false;

  off_file_write(ofile, geoms, sig_dgts);
  file_close_w(ofile);
  return true;
}

char *off_col(char *str, Color col)
{
  if (col.is_index())
    snprintf(str, MSG_SZ - 1, " %d", col.get_index());
  else if (col.is_value()) {
    if (col.get_transparency())
      vtostr(str, col.get_vec4d(), " ", -5);
    else
      vtostr(str, col.get_vec3d(), " ", -5);
  }
  else
    *str = '\0';

  return str;
}

void off_polys_write(FILE *ofile, const Geometry &geom, int offset)
{
  char col_str[MSG_SZ];
  for (unsigned int i = 0; i < geom.faces().size(); i++) {
    fprintf(ofile, "%lu", (unsigned long)geom.faces(i).size());
    for (unsigned int j = 0; j < geom.faces(i).size(); j++) {
      fprintf(ofile, " %d", geom.faces(i, j) + offset);
    }
    fprintf(ofile, " %s", off_col(col_str, geom.colors(FACES).get(i)));
    fprintf(ofile, "\n");
  }

  for (unsigned int i = 0; i < geom.edges().size(); i++) {
    fprintf(ofile, "2 %d %d", geom.edges(i, 0) + offset,
            geom.edges(i, 1) + offset);
    fprintf(ofile, " %s", off_col(col_str, geom.colors(EDGES).get(i)));
    fprintf(ofile, "\n");
  }
  // print coloured vertex elements
  map<int, Color>::const_iterator mi;
  for (mi = geom.colors(VERTS).get_properties().begin();
       mi != geom.colors(VERTS).get_properties().end(); mi++) {
    fprintf(ofile, "1 %d %s\n", (mi->first) + offset,
            off_col(col_str, mi->second));
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
