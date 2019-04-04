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

/* \file DisplayPoly.cc
   \brief display a polyhedron as plane faces, edge rods and vertex balls and
   with element number labels
*/

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#ifdef HAVE_GL_GL_H
#include <GL/gl.h>
#elif defined HAVE_OPENGL_GL_H
#include <OpenGL/gl.h>
#endif

#ifdef HAVE_GL_GLU_H
#include <GL/glu.h>
#elif defined HAVE_OPENGL_GLU_H
#include <OpenGL/glu.h>
#endif

#include <string.h>
#include <vector>

#include "../base/antiprism.h"
#include "displaypoly_gl.h"
#include "gl_writer.h"

using std::vector;

using namespace anti;

extern unsigned char stippleMask[17][128];

void gl_set_material(Color col = Vec3d(.5, .5, .5), bool trans = true,
                     int sides = GL_FRONT_AND_BACK);

void gl_set_material(Color col, bool trans, int sides)
{
  Vec4d cv = col.get_vec4d();
  GLfloat f_specular[] = {(GLfloat)cv[0] / 2, (GLfloat)cv[1] / 2,
                          (GLfloat)cv[2] / 2, 1.0};
  GLfloat f_diffuse[] = {(GLfloat)cv[0], (GLfloat)cv[1], (GLfloat)cv[2], 1.0};
  GLfloat f_shininess[] = {100.0};
  glMaterialfv(sides, GL_SPECULAR, f_specular);
  glMaterialfv(sides, GL_DIFFUSE, f_diffuse);
  glMaterialfv(sides, GL_SHININESS, f_shininess);
  if (trans)
    glPolygonStipple(stippleMask[int(cv[3] * 16 + 0.5)]);
}

void DisplayPoly_gl::gl_verts(const Scene &scen)
{
  double v_rad = get_vert_rad();
  double extra = 10 * v_rad / scen.get_width();
  int long_div = int(11 * (1 + extra));
  int lat_div = int(7 * (1 + extra));
  GLUquadric *quad = gluNewQuadric();
  GLuint sph = glGenLists(1);
  glNewList(sph, GL_COMPILE);
  gluSphere(quad, v_rad, long_div, lat_div);
  glEndList();

  const vector<Vec3d> &verts = disp_geom.verts();
  for (unsigned int i = 0; i < verts.size(); i++) {
    Color col = disp_geom.colors(VERTS).get((int)i);
    if (col.is_index())
      col = clrng(VERTS).get_col(col.get_index());
    if (!col.is_value())
      col = def_col(VERTS); // use default
    if (col.is_invisible())
      continue;

    gl_set_material(col, get_elem_trans(), GL_FRONT);
    glPushMatrix();
    glTranslated(verts[i][0], verts[i][1], verts[i][2]);
    glCallList(sph);
    glPopMatrix();
  }
  glDeleteLists(sph, 1);
  gluDeleteQuadric(quad);
}

static void edge_cyl_trans(Vec3d p1, Vec3d p2)
{
  Vec3d p1to2 = p2 - p1;
  glTranslated(p1[0], p1[1], p1[2]);
  Trans3d rot = Trans3d::rotate(Vec3d(0, 0, 1), p1to2);
  glMultMatrixd(rot.transpose().get_m());
  glScaled(1, 1, p1to2.len());
}

void DisplayPoly_gl::gl_edges(const Scene &scen)
{
  double e_rad = get_edge_rad();
  double extra = 10 * e_rad / scen.get_width();
  int long_div = int(11 * (1 + extra));
  int lat_div = 1;

  GLUquadric *quad = gluNewQuadric();
  GLuint cyl = glGenLists(1);
  glNewList(cyl, GL_COMPILE);
  gluCylinder(quad, e_rad, e_rad, 1, long_div, lat_div);
  glEndList();

  const vector<Vec3d> &verts = disp_geom.verts();
  const vector<vector<int>> &edges = disp_geom.edges();
  for (unsigned int i = 0; i < edges.size(); i++) {
    Color col = disp_geom.colors(EDGES).get((int)i);
    if (col.is_index())
      col = clrng(EDGES).get_col(col.get_index());
    if (!col.is_value())
      col = def_col(EDGES); // use default
    if (col.is_invisible())
      continue;

    gl_set_material(col, get_elem_trans(), GL_FRONT);
    glPushMatrix();
    edge_cyl_trans(verts[edges[i][0]], verts[edges[i][1]]);
    glCallList(cyl);
    glPopMatrix();
  }
  glDeleteLists(cyl, 1);
  gluDeleteQuadric(quad);
}

void DisplayPoly_gl::gl_faces(const Scene &)
{
  const vector<Vec3d> &verts = disp_geom.verts();
  const vector<vector<int>> &faces = disp_geom.faces();
  for (unsigned int i = 0; i < faces.size(); i++) {
    if (faces[i].size() < 3)
      continue;
    Color col = disp_geom.colors(FACES).get((int)i);
    if (col.is_index())
      col = clrng(FACES).get_col(col.get_index());
    if (!col.is_value())
      col = def_col(FACES); // use default
    if (col.is_invisible())
      continue;
    if (show_orientation) {
      gl_set_material(Color(1.0, 1.0, 1.0), get_elem_trans(), GL_FRONT);
      gl_set_material(Color(0.0, 0.0, 0.0), get_elem_trans(), GL_BACK);
    }
    else {
      if (get_transparency_type() == trans_50pc)
        col.set_rgba(col[0], col[1], col[2], 128);
      else if (get_transparency_type() == trans_0pc)
        col.set_rgba(col[0], col[1], col[2], 255);

      gl_set_material(col, get_elem_trans());
    }

    glBegin(GL_POLYGON);
    Vec3d norm = face_norm(verts, faces[i]);
    glNormal3dv(norm.get_v());
    for (unsigned int j = 0; j < faces[i].size(); j++)
      glVertex3dv(verts[faces[i][j]].get_v());
    glEnd();
  }
}

DisplayPoly_gl::DisplayPoly_gl()
    : DisplayPoly(), show_orientation(false), transparency_type(0)
{
}

void DisplayPoly_gl::gl_geom(const Scene &scen)
{
  if (elem(VERTS).get_show())
    gl_verts(scen);
  if (elem(FACES).get_show())
    gl_faces(scen);
  if (elem(EDGES).get_show())
    gl_edges(scen);
}

static void draw_text(char *str, double font_sz, Vec3d pos,
                      Vec3d norm = Vec3d(), int halign = 1, int valign = 1,
                      const anti_StrokeFont *font = ANTI_STROKE_ROMAN)
{
  const float drop = 33.33;
  const float height = 152.38;
  float scale = font_sz / height;
  float width = antiStrokeLength(font, (unsigned char *)str);

  float off_x = 0, off_y = 0;
  switch (halign) {
  case 0: // left
    off_x = 0;
    break;
  case 1: // centre
    off_x -= 0.5 * width;
    break;
  case 2: // right
    off_x -= width;
    break;
  }

  switch (valign) {
  case 0: // bottom
    off_y -= -drop;
    break;
  case 1: // centre
    off_y -= height / 2 - drop;
    break;
  case 2: // right
    off_y -= height - drop;
    break;
  case 3: // right
    off_y = 0;
    break;
  }

  glLineWidth(1.5);
  glPushMatrix();
  glTranslated(pos[0], pos[1], pos[2]);

  float m[16];
  glGetFloatv(GL_MODELVIEW_MATRIX, m);
  if (!norm.is_set()) {
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
        if (i == j)
          m[i * 4 + j] = 1;
        else
          m[i * 4 + j] = 0;
      }
    glLoadMatrixf(m);
  }
  else {
    Trans3d trans_cur; // Model coords to View coords
    // from_gl_matrixf(m, trans_cur);
    for (int i = 0; i < 16; i++)
      trans_cur[i] = m[i];
    trans_cur = trans_cur.transpose();
    trans_cur[3] = trans_cur[7] = trans_cur[11] = 0;

    // Work in View coords. Rotate label into position by rotating
    // about Y, then towards Y, to use View up direction for labels
    Vec3d n = trans_cur * norm.to_unit(); // norm in View coords
    Vec3d n_xz(n[0], 0.0, n[2]);
    auto ang = angle_around_axis(Vec3d::Z, n_xz, Vec3d::Y);
    const auto trans_tr = (trans_cur.transpose() * Trans3d::rotate(n_xz, n) *
                           Trans3d::rotate(Vec3d::Y, ang))
                              .transpose();
    for (int i = 0; i < 16; i++)
      m[i] = (float)trans_tr[i];
    glMultMatrixf(m);

    // Don't view at number labels from behind
    if (n[2] < 0)
      glRotatef(180, 0, 1, 0);
  }

  glScalef(scale, scale, scale);
  glTranslatef(off_x, off_y, 10);

  glNormal3f(0, 0, 1);
  for (char *c = str; *c; c++)
    antiStrokeCharacter(font, *c);
  glPopMatrix();
}

void DisplayNumLabels_gl::write_label(const Scene &scen, char *label, Vec3d pos,
                                      Vec3d norm)
{
  const float font_size = scen.cur_camera().get_text_sz(pos) * text_scale;
  draw_text(label, font_size, pos, norm);
}

void DisplayNumLabels_gl::gl_verts(const Scene &scen)
{
  Geometry &geom = sc_geom->get_geom();
  gl_set_material(get_label_col(elem(VERTS).get_col()), get_elem_trans());
  char label[64];
  const vector<Vec3d> &verts = geom.verts();
  for (unsigned int i = 0; i < verts.size(); i++) {
    if (geom.colors(VERTS).get((int)i).is_invisible())
      continue;
    sprintf(label, "%u", i);
    write_label(scen, label, sc_geom->get_v_label_pos(i));
  }
}

void DisplayNumLabels_gl::gl_edges(const Scene &scen)
{
  Geometry &geom = sc_geom->get_geom();
  gl_set_material(get_label_col(elem(EDGES).get_col()), get_elem_trans());
  char label[64];
  const vector<vector<int>> &edges = geom.edges();
  for (unsigned int i = 0; i < edges.size(); i++) {
    if (geom.colors(EDGES).get((int)i).is_invisible())
      continue;
    sprintf(label, "%u", i);
    write_label(scen, label, sc_geom->get_e_label_pos(i));
  }
}

void DisplayNumLabels_gl::gl_faces(const Scene &scen)
{
  Geometry &geom = sc_geom->get_geom();
  gl_set_material(get_label_col(elem(FACES).get_col()), get_elem_trans());
  char label[64];
  const vector<vector<int>> &faces = geom.faces();
  for (unsigned int i = 0; i < faces.size(); i++) {
    if (geom.colors(FACES).get((int)i).is_invisible())
      continue;
    sprintf(label, "%u", i);
    if (use_alt_labels)
      write_label(scen, label, sc_geom->get_geom().face_cent(i),
                  sc_geom->get_geom().face_norm(i));
    else
      write_label(scen, label, sc_geom->get_f_label_pos(i));
  }
}

void DisplayNumLabels_gl::gl_geom(const Scene &scen)
{
  if (elem(VERTS).get_show())
    gl_verts(scen);
  if (elem(FACES).get_show())
    gl_faces(scen);
  if (elem(EDGES).get_show())
    gl_edges(scen);
}

void DisplaySymmetry_gl::gl_geom(const Scene &scen)
{
  DisplayPoly_gl::gl_geom(scen);
}
