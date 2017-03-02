/*
   Copyright (c) 2016, Adrian Rossiter

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
   Name: off2dae.cc
   Description: convert from OFF to Collada DAE
   Project: Antiprism - http://www.antiprism.com
*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <stack>
#include <string>
#include <utility>
#include <vector>

#include "../base/antiprism.h"

using std::string;
using std::vector;
using std::stack;
using std::map;
using std::pair;

using namespace anti;

// --------------------------------------------------------------------------
// XML classes

class XML_elem {
private:
  string name;
  vector<pair<string, string>> attrs;

public:
  XML_elem(string nam = "") { name = nam; }
  XML_elem &add_attr(string att, string val)
  {
    attrs.push_back(pair<string, string>(att, val));
    return *this;
  }
  XML_elem &add_id(string val) { return add_attr("id", val); }
  const string &get_name() const { return name; }
  const vector<pair<string, string>> &get_attrs() const { return attrs; }
};

class XML_writer {
protected:
  FILE *ofile;
  stack<string> tag_names;
  int in_line_level;
  void open_tag(const XML_elem &tag, bool empty = false);
  void in_line() { in_line_level = tag_names.size(); }
  bool is_in_line(int off = 0) const
  {
    return in_line_level >= 0 && in_line_level < (int)tag_names.size() + off;
  }
  void in_line_update()
  {
    if (!is_in_line())
      in_line_level = -1;
  }

public:
  XML_writer(FILE *file) : in_line_level(-1) { ofile = file; };
  FILE *get_ofile() const { return ofile; }
  void print_open(const XML_elem &tag) { open_tag(tag, false); }
  void print_open(string name) { open_tag(XML_elem(name), false); }
  void print_open_empty(const XML_elem &tag) { open_tag(tag, true); }
  void print_open_empty(string name) { open_tag(XML_elem(name), true); }
  void print_open_close(const XML_elem &tag, string data);
  void print_open_close(string name, string data)
  {
    print_open_close(XML_elem(name), data);
  }

  void print_close();
  string indent(int off = 0) const
  {
    return string((off + (int)tag_names.size()) * 2, ' ');
  }
  void print_indent() const { fprintf(ofile, "%s", indent().c_str()); }
  void print_newline() const { fprintf(ofile, "\n"); }
};

// --------------------------------------------------------------------------
// XML implementation

void XML_writer::open_tag(const XML_elem &tag, bool empty)
{
  fprintf(ofile, "%s<%s",
          is_in_line() ? "" : indent().c_str(), // if in_line, suppress indent
          tag.get_name().c_str());
  vector<pair<string, string>>::const_iterator vi;
  for (vi = tag.get_attrs().begin(); vi != tag.get_attrs().end(); ++vi)
    fprintf(ofile, " %s=\"%s\"", vi->first.c_str(), vi->second.c_str());
  if (!empty)
    tag_names.push(tag.get_name());
  fprintf(ofile, "%s>%s", empty ? " /" : "",
          is_in_line() ? "" : "\n"); // if in_line, suppress newline
  in_line_update();
}

void XML_writer::print_close()
{
  fprintf(ofile, "%s</%s>%s",
          is_in_line() ? "" : indent(-1).c_str(), // if in_line, suppress indent
          tag_names.top().c_str(),
          is_in_line(-1) ? "" : "\n"); // if final close of in_line, add newline
  tag_names.pop();
  in_line_update();
}

void XML_writer::print_open_close(const XML_elem &tag, string data)
{
  if (!is_in_line())
    in_line();
  print_open(tag);
  fprintf(ofile, "%s", data.c_str());
  print_close();
}

// --------------------------------------------------------------------------
// Collada class

class Collada_writer : public XML_writer {
private:
  int smooth_level; // controls number faces in sphere and cylinder elems

  static void get_unique_cols(const map<Color, vector<int>> *col2elems,
                              map<Color, string> *cols);
  static void add_col2elems(DisplayPoly &disp,
                            map<Color, vector<int>> *col2elems);
  static void get_col2elems(const Scene &scen,
                            map<Color, vector<int>> *col2elems);

  static string col2hex(const Color &col);
  static string format_vec(const Vec3d &v, int sig_digits);
  static string id_tag(const string &name, int g, int d);

  void print_asset();

  void print_effect(const Color &col, const string &col_str);
  void print_material(const string &col_str);
  void print_library_effects_and_materials(const Scene &scen);

  void print_vertex_coords(const vector<Vec3d> &vecs, int sig_digits);
  void print_geometry_vecs(const string &id, const vector<Vec3d> &vecs,
                           int sig_digits);
  void print_face_sizes(const Geometry &geom, const vector<int> &face_idxs);
  void print_face_indices(const Geometry &geom, const vector<int> &face_idxs,
                          bool with_normals);
  void print_geometry(const string &id, const Geometry &geom,
                      const map<Color, vector<int>> &g_col2elems,
                      bool triangulate, int sig_digits,
                      const vector<Vec3d> *normals = nullptr,
                      const string &mat = "");
  void print_geometry_sphere();
  void print_geometry_cylinder();
  void print_library_geometries(const Scene &scen, int sig_digits);

  void print_library_visual_scene_faces(DisplayPoly &disp, const string &id);
  void print_library_visual_scene_verts(DisplayPoly &disp, const string &id);
  void print_library_visual_scene_edges(DisplayPoly &disp, const string &id);
  void print_library_visual_scenes(const Scene &scen);

  void print_scene();

public:
  Collada_writer(FILE *file, int smooth) : XML_writer(file)
  {
    set_smooth_level(smooth);
  }
  void set_smooth_level(int smooth) { smooth_level = smooth < 1 ? 1 : smooth; }
  bool write(const Scene &scen, int sig_digits = 16);
};

// --------------------------------------------------------------------------
// Collada class implementation

void Collada_writer::get_unique_cols(const map<Color, vector<int>> *col2elems,
                                     map<Color, string> *cols)
{
  // Copy the colors from all elements into a single set
  map<Color, vector<int>>::const_iterator mi;
  for (int elem = 0; elem < 3; elem++)
    for (mi = col2elems[elem].begin(); mi != col2elems[elem].end(); ++mi)
      (*cols)[mi->first] = col2hex(mi->first);
}

void Collada_writer::add_col2elems(DisplayPoly &disp,
                                   map<Color, vector<int>> *col2elems)
{
  if (disp.elem(VERTS).get_show()) {
    const vector<Vec3d> &vs = disp.get_disp_geom().verts();
    for (unsigned int i = 0; i < vs.size(); i++) {
      Color col = disp.get_disp_geom().colors(VERTS).get((int)i);
      if (col.is_idx())
        col = disp.clrng(VERTS).get_col(col.get_idx());
      if (!col.is_val())
        col = disp.def_col(VERTS); // use default
      if (col.is_val() && !col.is_inv())
        col2elems[0][col].push_back(i);
    }
  }
  if (disp.elem(EDGES).get_show()) {
    const vector<vector<int>> &es = disp.get_disp_geom().edges();
    for (unsigned int i = 0; i < es.size(); i++) {
      Color col = disp.get_disp_geom().colors(EDGES).get((int)i);
      if (col.is_idx())
        col = disp.clrng(EDGES).get_col(col.get_idx());
      if (!col.is_val())
        col = disp.def_col(EDGES); // use default
      if (col.is_val() && !col.is_inv())
        col2elems[1][col].push_back(i);
    }
  }
  if (disp.elem(FACES).get_show()) {
    const vector<vector<int>> &fs = disp.get_disp_geom().faces();
    for (unsigned int i = 0; i < fs.size(); i++) {
      if (fs[i].size() < 3) // skip degenerate polygons
        continue;
      Color col = disp.get_disp_geom().colors(FACES).get((int)i);
      if (col.is_idx())
        col = disp.clrng(FACES).get_col(col.get_idx());
      if (!col.is_val())
        col = disp.def_col(FACES); // use default
      if (col.is_val() && !col.is_inv())
        col2elems[2][col].push_back(i);
    }
  }
}

void Collada_writer::get_col2elems(const Scene &scen,
                                   map<Color, vector<int>> col2elems[3])
{
  for (int elem = 0; elem < 3; elem++)
    col2elems[elem].clear();

  vector<SceneGeometry>::const_iterator geo;
  for (geo = scen.get_geoms().begin(); geo != scen.get_geoms().end(); ++geo) {
    vector<GeometryDisplay *>::const_iterator dsp;
    for (dsp = geo->get_disps().begin(); dsp != geo->get_disps().end(); ++dsp) {
      if (DisplayPoly *disp_p = dynamic_cast<DisplayPoly *>(*dsp))
        add_col2elems(*disp_p, col2elems);
    }
  }
}

string Collada_writer::col2hex(const Color &col)
{
  char val[9];
  sprintf(val, "%02x%02x%02x%02x", col[0], col[1], col[2], col[3]);
  return val;
}

string Collada_writer::format_vec(const Vec3d &v, int sig_digits)
{
  char buf[MSG_SZ];
  if (sig_digits > 0)
    snprintf(buf, MSG_SZ, "%.*g %.*g %.*g", sig_digits, v[0], sig_digits, v[1],
             sig_digits, v[2]);
  else
    snprintf(buf, MSG_SZ, "%.*f %.*f %.*f", -sig_digits, v[0], -sig_digits,
             v[1], -sig_digits, v[2]);
  return buf;
}

string Collada_writer::id_tag(const string &name, int g, int d)
{
  return itostr(g) + "_" + itostr(d) + "_" + name;
}

void Collada_writer::print_asset()
{
  print_open("asset");

  print_open("contributor");
  print_open_close("authoring_tool",
                   "off2dae (Antiprism - http://www.antiprism.com)");
  print_close(); // contributor

  time_t now;
  time(&now);
  char datetime_iso8601[MSG_SZ];
  strftime(datetime_iso8601, MSG_SZ, "%FT%TZ", gmtime(&now));
  print_open_close("created", datetime_iso8601);
  print_open_close("modified", datetime_iso8601);

  print_open_empty(
      XML_elem("unit").add_attr("meter", "1").add_attr("name", "meter"));
  print_open_close("up_axis", "Z_UP");

  print_close(); // asset
}

void Collada_writer::print_effect(const Color &col, const string &col_str)
{
  print_open(XML_elem("effect").add_id("c_" + col_str));
  print_open("profile_COMMON");
  print_open(XML_elem("technique").add_attr("sid", "common"));
  print_open("phong");

  Vec4d cv = col.get_Vec4d();
  in_line();
  print_open("diffuse");
  print_open("color");
  fprintf(ofile, "%.4f %.4f %.4f %.4f", cv[0], cv[1], cv[2], cv[3]);
  print_close(); // color
  print_close(); // diffuse

  in_line();
  print_open("transparency");
  print_open("float");
  fprintf(ofile, "%4f", cv[3]); // alpha from color
  print_close();                // float
  print_close();                // transparency

  print_close(); // phong
  print_close(); // technique
  print_close(); // profile_COMMON
  print_close(); // effect
}

void Collada_writer::print_material(const string &col_str)
{
  in_line();
  print_open(XML_elem("material").add_id("m_" + col_str));
  print_open_empty(
      XML_elem("instance_effect").add_attr("url", "#c_" + col_str));
  print_close(); // material
}

void Collada_writer::print_library_effects_and_materials(const Scene &scen)
{
  print_open("library_effects");

  map<Color, vector<int>> col2elems[3];
  get_col2elems(scen, col2elems);
  map<Color, string> cols;
  get_unique_cols(col2elems, &cols);

  map<Color, string>::const_iterator mi;
  for (mi = cols.begin(); mi != cols.end(); ++mi)
    print_effect(mi->first, mi->second);

  print_close(); // library_effects

  print_open("library_materials");
  for (mi = cols.begin(); mi != cols.end(); ++mi)
    print_material(mi->second);

  print_close(); // library_materials
}

void Collada_writer::print_vertex_coords(const vector<Vec3d> &vecs,
                                         int sig_digits)
{
  print_indent();
  for (const auto &vec : vecs)
    fprintf(ofile, "%s ", format_vec(vec, sig_digits).c_str());
  print_newline();
}

void Collada_writer::print_geometry_vecs(const string &id,
                                         const vector<Vec3d> &vecs,
                                         int sig_digits)
{
  print_open(XML_elem("source").add_id("source_" + id));

  print_open(XML_elem("float_array")
                 .add_id("float_array_" + id)
                 .add_attr("count", itostr(3 * vecs.size())));
  print_vertex_coords(vecs, sig_digits);
  print_close(); // float_array

  print_open("technique_common");
  print_open(XML_elem("accessor")
                 .add_attr("count", itostr(vecs.size()))
                 .add_attr("source", "#float_array_" + id)
                 .add_attr("offset", "0")
                 .add_attr("stride", "3"));

  print_open_empty(
      XML_elem("param").add_attr("name", "X").add_attr("type", "float"));
  print_open_empty(
      XML_elem("param").add_attr("name", "Y").add_attr("type", "float"));
  print_open_empty(
      XML_elem("param").add_attr("name", "Z").add_attr("type", "float"));

  print_close(); // accessor
  print_close(); // technique_common
  print_close(); // source
}

void Collada_writer::print_face_sizes(const Geometry &geom,
                                      const vector<int> &face_idxs)
{
  for (int f_idx : face_idxs) {
    fprintf(ofile, "%lu ", (unsigned long)geom.faces(f_idx).size());
  }
}

void Collada_writer::print_face_indices(const Geometry &geom,
                                        const vector<int> &face_idxs,
                                        bool with_normals)
{
  for (int f_idx : face_idxs) {
    for (unsigned int j = 0; j < geom.faces(f_idx).size(); j++)
      for (int cnt = 0; cnt < 1 + with_normals; cnt++)
        fprintf(ofile, "%d ", geom.faces(f_idx, j));
  }
}

void Collada_writer::print_geometry(const string &id, const Geometry &geom,
                                    const map<Color, vector<int>> &g_col2elems,
                                    bool triangulate, int sig_digits,
                                    const vector<Vec3d> *normals,
                                    const string &mat)
{
  print_open(XML_elem("geometry").add_id("geometry_" + id));
  print_open("mesh");

  if (normals)
    print_geometry_vecs("normals_" + id, *normals, 6); // low precisionn

  print_geometry_vecs(id, geom.verts(), sig_digits);

  print_open(XML_elem("vertices").add_id("vertices_" + id));
  print_open_empty(XML_elem("input")
                       .add_attr("semantic", "POSITION")
                       .add_attr("source", "#source_" + id));
  print_close(); // vertices

  map<Color, vector<int>>::const_iterator mi;
  for (mi = g_col2elems.begin(); mi != g_col2elems.end(); ++mi) {
    const Color &face_col = mi->first;
    const vector<int> face_idxs = mi->second;
    string imat = "im_" + ((mat != "") ? mat : col2hex(face_col));
    print_open(XML_elem(triangulate ? "triangles" : "polylist")
                   .add_attr("material", imat)
                   .add_attr("count", itostr(face_idxs.size())));

    print_open_empty(XML_elem("input")
                         .add_attr("semantic", "VERTEX")
                         .add_attr("offset", "0")
                         .add_attr("source", "#vertices_" + id));

    if (normals)
      print_open_empty(XML_elem("input")
                           .add_attr("semantic", "NORMAL")
                           .add_attr("offset", "1")
                           .add_attr("source", "#source_normals_" + id));

    if (!triangulate) {
      in_line();
      print_open("vcount");
      print_face_sizes(geom, face_idxs);
      print_close(); // vcount
    }

    in_line();
    print_open("p");
    print_face_indices(geom, face_idxs, normals);
    print_close(); // p

    print_close(); // polylist
  }
  print_close(); // mesh
  print_close(); // geometry
}

void Collada_writer::print_geometry_sphere()
{
  Geometry geom;
  geom.read_resource("std_geo_" + itostr(smooth_level));
  map<Color, vector<int>> g_col2fs;
  for (int i = 0; i < (int)geom.faces().size(); i++)
    g_col2fs[Color(0)].push_back(i);

  print_geometry("vertex_sphere", geom, g_col2fs, true, 6, &geom.verts(),
                 "vertex_sphere");
}

void Collada_writer::print_geometry_cylinder()
{
  Polygon pri(6 * smooth_level, 1, Polygon::prism);
  pri.set_radius(1, 1.0);
  pri.set_edge(2, 1.0);
  Geometry geom;
  pri.make_poly(geom);
  vector<Vec3d> norms = geom.verts();
  for (auto &norm : norms)
    norm[2] = 0; // normals are perpendicular to the axis

  map<Color, vector<int>> g_col2fs;
  for (int i = 0; i < (int)geom.faces().size(); i++)
    g_col2fs[Color(0)].push_back(i);

  print_geometry("edge_cylinder", geom, g_col2fs, false, 6, &norms,
                 "edge_cylinder");
}

void Collada_writer::print_library_geometries(const Scene &scen, int sig_digits)
{
  print_open("library_geometries");
  bool include_vert_geometry = false;
  bool include_edge_geometry = false;
  int g = 0;
  vector<SceneGeometry>::const_iterator geo;
  for (geo = scen.get_geoms().begin(); geo != scen.get_geoms().end(); ++geo) {
    vector<GeometryDisplay *>::const_iterator dsp;
    int d = 0;
    for (dsp = geo->get_disps().begin(); dsp != geo->get_disps().end(); ++dsp) {
      if (DisplayPoly *disp_p = dynamic_cast<DisplayPoly *>(*dsp)) {
        if (disp_p->elem(FACES).get_show()) {
          map<Color, vector<int>> g_col2elems[3];
          add_col2elems(*disp_p, g_col2elems);
          print_geometry(id_tag(geo->get_name(), g, d), disp_p->get_disp_geom(),
                         g_col2elems[2], disp_p->get_triangulate(), sig_digits);
        }
        if (disp_p->elem(VERTS).get_show())
          include_vert_geometry = true;
        if (disp_p->elem(EDGES).get_show())
          include_edge_geometry = true;
      }
    }
    g++;
  }

  if (include_vert_geometry)
    print_geometry_sphere();
  if (include_edge_geometry)
    print_geometry_cylinder();

  print_close(); // library_geometries
}

void Collada_writer::print_library_visual_scene_faces(DisplayPoly &disp,
                                                      const string &id)
{
  print_open(XML_elem("instance_geometry").add_attr("url", "#geometry_" + id));
  print_open("bind_material");
  print_open("technique_common");

  map<Color, vector<int>> g_col2elems[3];
  add_col2elems(disp, g_col2elems);
  map<Color, string> cols;
  get_unique_cols(g_col2elems, &cols);
  map<Color, string>::const_iterator mi;
  for (mi = cols.begin(); mi != cols.end(); ++mi)
    print_open_empty(XML_elem("instance_material")
                         .add_attr("symbol", "im_" + mi->second)
                         .add_attr("target", "#m_" + mi->second));

  print_close(); // technique_common
  print_close(); // bind_materials
  print_close(); // instance_geometry
}

void Collada_writer::print_library_visual_scene_verts(DisplayPoly &disp,
                                                      const string &id)
{
  print_open(XML_elem("node").add_id("geometry_verts" + id));
  const Geometry &dgeom = disp.get_disp_geom();
  for (unsigned int i = 0; i < dgeom.verts().size(); i++) {
    Color col = disp.get_disp_geom().colors(VERTS).get((int)i);
    if (col.is_idx())
      col = disp.clrng(VERTS).get_col(col.get_idx());
    if (!col.is_val())
      col = disp.def_col(VERTS); // use default
    if (col.is_val() && !col.is_inv()) {
      print_open("node");

      print_open_close("translate", format_vec(dgeom.verts(i), 8));
      print_open_close("scale",
                       format_vec(Vec3d(1, 1, 1) * disp.get_vert_rad(), 8));
      print_open(XML_elem("instance_geometry")
                     .add_attr("url", "#geometry_vertex_sphere"));
      print_open("bind_material");
      print_open("technique_common");
      print_open_empty(XML_elem("instance_material")
                           .add_attr("symbol", "im_vertex_sphere")
                           .add_attr("target", "#m_" + col2hex(col)));
      print_close(); // technique_common
      print_close(); // bind_materials
      print_close(); // instance_geometry
      print_close(); // node
    }
  }
  print_close(); // node
}

void Collada_writer::print_library_visual_scene_edges(DisplayPoly &disp,
                                                      const string &id)
{
  print_open(XML_elem("node").add_id("geometry_edges" + id));
  const Geometry &dgeom = disp.get_disp_geom();
  for (unsigned int i = 0; i < dgeom.edges().size(); i++) {
    Color col = disp.get_disp_geom().colors(EDGES).get((int)i);
    if (col.is_idx())
      col = disp.clrng(EDGES).get_col(col.get_idx());
    if (!col.is_val())
      col = disp.def_col(EDGES); // use default
    if (col.is_val() && !col.is_inv()) {
      print_open("node");

      Vec3d mid = dgeom.edge_cent(i);
      Vec3d dir = dgeom.edge_vec(i);
      double ht = dir.len();
      dir /= ht;                                 // to unit
      double ang = -acos(safe_for_trig(dir[2])); // ang betwn dir and y-axis
      Vec3d axis = vcross(dir, Vec3d(0, 0, 2)).unit(); // axis
      print_open_close("translate", format_vec(mid, 8));
      print_open_close("rotate",
                       format_vec(axis, 8) + " " + dtostr(rad2deg(ang), 8));
      print_open_close(
          "scale",
          format_vec(Vec3d(disp.get_edge_rad(), disp.get_edge_rad(), ht), 8));

      print_open(XML_elem("instance_geometry")
                     .add_attr("url", "#geometry_edge_cylinder"));
      print_open("bind_material");
      print_open("technique_common");
      print_open_empty(XML_elem("instance_material")
                           .add_attr("symbol", "im_edge_cylinder")
                           .add_attr("target", "#m_" + col2hex(col)));
      print_close(); // technique_common
      print_close(); // bind_materials
      print_close(); // instance_geometry
      print_close(); // node
    }
  }
  print_close(); // node
}

void Collada_writer::print_library_visual_scenes(const Scene &scen)
{
  print_open("library_visual_scenes");
  print_open(
      XML_elem("visual_scene").add_id("scene_0").add_attr("name", "Scene_0"));
  print_open(XML_elem("node").add_id("visgeom_0").add_attr("type", "NODE"));

  int g = 0; // geometry order index
  vector<SceneGeometry>::const_iterator geo;
  for (geo = scen.get_geoms().begin(); geo != scen.get_geoms().end(); ++geo) {
    int d = 0; // geometry display type order index
    vector<GeometryDisplay *>::const_iterator dsp;
    for (dsp = geo->get_disps().begin(); dsp != geo->get_disps().end(); ++dsp) {
      if (DisplayPoly *disp_p = dynamic_cast<DisplayPoly *>(*dsp)) {
        string id = id_tag(geo->get_name(), g, d);
        if (disp_p->elem(FACES).get_show())
          print_library_visual_scene_faces(*disp_p, id);
        if (disp_p->elem(VERTS).get_show())
          print_library_visual_scene_verts(*disp_p, id);
        if (disp_p->elem(EDGES).get_show())
          print_library_visual_scene_edges(*disp_p, id);
      }
      d++;
    }
    g++;
  }
  print_close(); // node

  print_close(); // visual_scenes
  print_close(); // library_visual_scenes
}

void Collada_writer::print_scene()
{
  print_open("scene");
  print_open_empty(
      XML_elem("instance_visual_scene").add_attr("url", "#scene_0"));
  print_close(); // scene
}

bool Collada_writer::write(const Scene &scen, int sig_digits)
{
  fprintf(ofile, "<?xml version=\"1.0\"?>\n");
  print_open(
      XML_elem("COLLADA")
          .add_attr("xmlns", "http://www.collada.org/2005/11/COLLADASchema")
          .add_attr("version", "1.4.1"));

  print_asset();
  print_library_effects_and_materials(scen);
  print_library_geometries(scen, sig_digits);
  print_library_visual_scenes(scen);
  print_scene();

  print_close(); // collada
  return true;
}

// --------------------------------------------------------------------------
// Options class

class o2d_opts : public ViewOpts {
public:
  int smooth_lvl;
  int sig_digits;
  string ofile;

  o2d_opts() : ViewOpts("off2dae"), smooth_lvl(3), sig_digits(DEF_SIG_DGTS) {}
  void usage();
  void process_command_line(int argc, char **argv);
};

// --------------------------------------------------------------------------
// Options implementation

// clang-format off
void o2d_opts::usage()
{
   fprintf(stdout,
"\n"
"Usage: %s [options] [input_file]\n"
"\n"
"Convert an OFF file to Collada DAE file format\n"
"\n"
"Options\n"
"%s"
"  -v <rad>  radius of vertex spheres, or 'b' to have radius of balls\n"
"            of the maximum size without overlap (default: ball_rad/15)\n"
"  -e <rad>  radius of edge cylinders (default: vertex_rad/1.5)\n"
"  -V <col>  default vertex colour, in form 'R,G,B,A' (3 or 4 values\n"
"            0.0-1.0, or 0-255) or hex 'xFFFFFF' (default: 1.0,0.5,0.0)\n"
"  -E <col>  default edge colour, in form 'R,G,B,A' (3 or 4 values\n"
"            0.0-1.0, or 0-255) or hex 'xFFFFFF', 'x' to hide implicit edges\n"
"            (default: 0.8,0.6,0.8)\n"
"  -F <col>  default face colour, in form 'R,G,B,A' (3 or 4 values\n"
"            0.0-1.0, or 0-255) or hex 'xFFFFFF' (default: 0.8,0.9,0.9)\n"
"  -x <elms> hide elements. The element string can include v, e and f\n"
"            to hide vertices, edges and faces\n"
"  -t <disp> select face parts to display according to winding number from:\n"
"            odd, nonzero (default), positive, negative, no_triangulation\n"
"            (use native polygon display)\n"
"  -m <maps> a comma separated list of colour maps used to transform colour\n"
"            indexes, a part consisting of letters from v, e, f, selects \n"
"            the element types to apply the map list to (default 'vef').\n"
"  -K <lvl>  a positive integer to specify the level of smoothing (using more\n"
"            polygons) of the vertex spheres and edge cylinders (default: %d)\n"
"  -d <dgts> number of significant digits (default %d) or if negative\n"
"            then the number of digits after the decimal point\n"
"  -o <file> write output to file (default: write to standard output)\n"
"\n"
"\n", prog_name(), help_ver_text, smooth_lvl, DEF_SIG_DGTS);
}
// clang-format on

void o2d_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;
  Status stat;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hv:e:V:E:F:x:t:m:K:o:d:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {

    case 'K':
      print_status_or_exit(read_int(optarg, &smooth_lvl), c);
      if (smooth_lvl < 1)
        error("smooth level cannot be less than 1", c);
      if (smooth_lvl > 7)
        warning("higher values increase the model size and may not be "
                "visually distinguishable from a lower value",
                c);
      break;

    case 'd':
      print_status_or_exit(read_int(optarg, &sig_digits), c);
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      if (!(stat = read_disp_option(c, optarg))) {
        if (stat.is_warning())
          warning(stat.msg(), c);
        else
          error(stat.msg(), c);
      }
    }
  }

  if (!get_geom_defs().elem(VERTS).get_show() &&
      !get_geom_defs().elem(EDGES).get_show() &&
      !get_geom_defs().elem(FACES).get_show())
    error("cannot hide all elements", 'x');

  if (argc - optind >= 1)
    while (argc - optind >= 1)
      ifiles.push_back(argv[optind++]);
  else
    ifiles.push_back("");
}

// --------------------------------------------------------------------------
// Main

int main(int argc, char *argv[])
{
  o2d_opts opts;
  opts.process_command_line(argc, argv);
  Scene scen = opts.scen_defs;
  opts.set_view_vals(scen);

  FILE *ofile = stdout; // write to stdout by default
  if (opts.ofile != "") {
    ofile = fopen(opts.ofile.c_str(), "w");
    if (ofile == nullptr)
      opts.error("could not open output file \'" + opts.ofile + "\'");
  }

  Collada_writer collada(ofile, opts.smooth_lvl);
  collada.write(scen, opts.sig_digits);

  if (opts.ofile != "")
    fclose(ofile);

  return 0;
}
