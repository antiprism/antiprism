/*
   Copyright (c) 2006-2016, Adrian Rossiter

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

/* \file vw_glut.cc
   \brief antiview - interface with GLUT
*/

#include <vector>

#include "displaypoly_gl.h"
#include "vw_glut.h"

using std::vector;

using namespace anti;

extern glut_state glut_s;

void glut_state::make_menu()
{
  int show_menu = glutCreateMenu(menu_cb);
  glutAddMenuEntry("Vertices 'V'", 'V');
  glutAddMenuEntry("Edges 'E'", 'E');
  glutAddMenuEntry("Faces 'F'", 'F');
  glutAddMenuEntry("Vertex Nums 'n'", 'n');
  glutAddMenuEntry("Face Nums 'N'", 'N');
  glutAddMenuEntry("Edge Nums 'M'", 'M');
  glutAddMenuEntry("Transparency 'T'", 'T');
  glutAddMenuEntry("Sym Elems 'Y'", 'Y');
  glutAddMenuEntry("Orientation 'O'", 'O');

  glutCreateMenu(menu_cb);
  glutAddMenuEntry("Rotate  'r'", 'r');
  glutAddMenuEntry("Drag  'd'", 'd');
  glutAddMenuEntry("Zoom  'z'", 'z');
  glutAddMenuEntry("Spin Control 's'", 's');
  glutAddMenuEntry("Slice  'S'", 'S');
  glutAddMenuEntry("Proj Type  'P'", 'P');
  glutAddSubMenu("Show/Hide", show_menu);
  glutAddMenuEntry("Reset  'R'", 'R');
  glutAddMenuEntry("Quit  'Q'", 'Q');
  glutAttachMenu(GLUT_RIGHT_BUTTON);
}

void glut_state::set_op(int op)
{
  cur_op = op;
  if (cur_op == tr_spin) {
    scen.cur_camera().stop_spinning();
    // glutIdleFunc(glut_idle_cb);
  }
  // else {
  //   if(!gl_cam.is_spinning())
  //      glutIdleFunc(0);
  //}
}

void glut_state::change_sym_disp()
{
  sym_disp_type = (sym_disp_type + 1) % 4;
  vector<SceneGeometry>::iterator geo;
  for (geo = scen.get_geoms().begin(); geo != scen.get_geoms().end(); ++geo) {
    if (DisplaySymmetry *sym =
            dynamic_cast<DisplaySymmetry *>(geo->get_sym())) {
      sym->set_show_axes(sym_disp_type);
      sym->set_show_mirrors(sym_disp_type > 1);
      sym->set_show_rotrefls(sym_disp_type > 2);
    }
  }
}

void glut_state::change_transparency_disp()
{
  transparency_disp_type = (transparency_disp_type + 1) % 3;
  vector<SceneGeometry>::const_iterator geo;
  for (geo = glut_s.scen.get_geoms().begin();
       geo != glut_s.scen.get_geoms().end(); ++geo) {
    vector<GeometryDisplay *>::const_iterator disp;
    for (disp = geo->get_disps().begin(); disp != geo->get_disps().end();
         ++disp) {
      if (DisplayPoly_gl *disp_p = dynamic_cast<DisplayPoly_gl *>(*disp)) {
        disp_p->set_transparency_type(transparency_disp_type);
      }
    }
  }
}

static void toggle(char elem)
{
  vector<SceneGeometry>::const_iterator geo;
  for (geo = glut_s.scen.get_geoms().begin();
       geo != glut_s.scen.get_geoms().end(); ++geo) {
    vector<GeometryDisplay *>::const_iterator disp;
    for (disp = geo->get_disps().begin(); disp != geo->get_disps().end();
         ++disp) {
      if (DisplayPoly *disp_p = dynamic_cast<DisplayPoly *>(*disp)) {
        if (elem == 'v')
          disp_p->elem(VERTS).set_show(!disp_p->elem(VERTS).get_show());
        else if (elem == 'e')
          disp_p->elem(EDGES).set_show(!disp_p->elem(EDGES).get_show());
        else if (elem == 'f')
          disp_p->elem(FACES).set_show(!disp_p->elem(FACES).get_show());
      }
      if (DisplayPoly_gl *disp_p = dynamic_cast<DisplayPoly_gl *>(*disp)) {
        if (elem == 'O')
          disp_p->set_show_orientation(!disp_p->get_show_orientation());
      }
    }
  }
}

static void toggle_lab(char elem)
{
  vector<SceneGeometry>::const_iterator geo;
  for (geo = glut_s.scen.get_geoms().begin();
       geo != glut_s.scen.get_geoms().end(); ++geo) {
    if (GeometryDisplayLabel *lab = geo->get_label()) {
      if (elem == 'n')
        lab->elem(VERTS).set_show(!lab->elem(VERTS).get_show());
      else if (elem == 'm')
        lab->elem(EDGES).set_show(!lab->elem(EDGES).get_show());
      else if (elem == 'N')
        lab->elem(FACES).set_show(!lab->elem(FACES).get_show());
      else if (elem == 'c')
        lab->set_label_light(!lab->get_label_light());
      else if (elem == 'C')
        lab->set_label_invert(!lab->get_label_invert());
    }
  }
}

static void scale(char elem, double scal)
{
  vector<SceneGeometry>::iterator geo;
  for (geo = glut_s.scen.get_geoms().begin();
       geo != glut_s.scen.get_geoms().end(); ++geo) {
    vector<GeometryDisplay *>::iterator disp;
    for (disp = geo->get_disps().begin(); disp != geo->get_disps().end();
         ++disp) {
      if (DisplayPoly *disp_p = dynamic_cast<DisplayPoly *>(*disp)) {
        if (elem == 'v')
          disp_p->elem(VERTS).set_size(disp_p->get_vert_rad() * scal);
        else if (elem == 'e')
          disp_p->elem(EDGES).set_size(disp_p->get_edge_rad() * scal);
      }
    }
  }
}

void keyboard_cb(unsigned char key, int /*x*/, int /*y*/)
{
  double elem_scale = 1.07;
  switch (key) {
  case 'Q':
  case 'q':
  case 27:
    exit(0);
    break;
  case 'P':
  case 'p':
    glut_s.toggle_projection();
    glut_s.reshape();
    break;
  case 'R':
    glut_s.restore_camera();
    glut_s.set_projection();
    glut_s.reshape();
    break;
  case 'r':
    glut_s.set_op(glut_state::tr_rot);
    break;
  case 'Z':
  case 'z':
    glut_s.set_op(glut_state::tr_z);
    break;
  case 'D':
  case 'd':
    glut_s.set_op(glut_state::tr_xy);
    break;
  case 'S':
    glut_s.set_op(glut_state::tr_slice);
    break;
  case 's':
    glut_s.set_op(glut_state::tr_spin);
    break;
  case 'Y':
  case 'y':
    glut_s.change_sym_disp();
    break;
  case 'O':
    toggle('O');
    break;
  case 'V':
  case 'v':
    toggle('v');
    break;
  case 'T':
  case 't':
    glut_s.change_transparency_disp();
    break;
  case 'E':
  case 'e':
    toggle('e');
    break;
  case 'F':
  case 'f':
    toggle('f');
    break;
  case 'N':
    toggle_lab('N');
    break;
  case 'n':
    toggle_lab('n');
    break;
  case 'M':
  case 'm':
    toggle_lab('m');
    break;
  case 'c':
    toggle_lab('c');
    break;
  case 'C':
    toggle_lab('C');
    break;
  case 'J':
    scale('v', elem_scale);
    scale('e', elem_scale);
    break;
  case 'j':
    scale('v', 1 / elem_scale);
    scale('e', 1 / elem_scale);
    break;
  case 'K':
    scale('v', elem_scale);
    break;
  case 'k':
    scale('v', 1 / elem_scale);
    break;
  case 'L':
    scale('e', elem_scale);
    break;
  case 'l':
    scale('e', 1 / elem_scale);
    break;
  default:
    break;
  }

  glutPostRedisplay();
}

void special_cb(int key, int /*x*/, int /*y*/)
{
  int x_inc = 0;
  int y_inc = 0;
  int inc = 3;
  switch (key) {
  case GLUT_KEY_UP:
    y_inc -= inc;
    break;
  case GLUT_KEY_DOWN:
    y_inc += inc;
    break;
  case GLUT_KEY_LEFT:
    x_inc -= inc;
    break;
  case GLUT_KEY_RIGHT:
    x_inc += inc;
    break;
  }

  if (x_inc || y_inc) {
    mouse_cb(GLUT_LEFT_BUTTON, GLUT_DOWN, 0, 0);
    mouse_cb(GLUT_LEFT_BUTTON, GLUT_UP, x_inc, y_inc);
  }
}

void motion_cb(int x, int y)
{
  if (glut_s.l_button == GLUT_DOWN) {
    double x_inc = x - glut_s.start_x;
    double y_inc = -y + glut_s.start_y;
    double scale;
    Camera &cam = glut_s.scen.cur_camera();
    switch (glut_s.get_op()) {
    case glut_state::tr_rot: { // scale to window size
      x_inc = fmod(x_inc, 360) * M_PI / 180;
      y_inc = fmod(y_inc, 360) * M_PI / 180;
      cam.inc_rotation(Trans3d::rot(Vec3d(-y_inc, x_inc, 0.0),
                                    sqrt(x_inc * x_inc + y_inc * y_inc)));
      break;
    }
    case glut_state::tr_xy:
      scale = glut_s.trans_scale * glut_s.scen.get_width();
      cam.inc_lookat(Vec3d(-x_inc, -y_inc, 0) * scale);
      break;
    case glut_state::tr_z:
      scale = glut_s.trans_scale * glut_s.scen.get_width();
      cam.inc_lookat(Vec3d(0, 0, y_inc) * scale);
      break;
    case glut_state::tr_slice:
      scale = glut_s.slice_scale * glut_s.scen.get_width();
      cam.inc_cut_dist(y_inc * scale);
      glut_s.reshape();
      break;
    case glut_state::tr_spin: { // scale to window size
      double x_off = 2.0 * x / glut_s.width - 1;
      double y_off = -2.0 * y / glut_s.height + 1;
      x_off = fmod(x_off, 100.0) * 2.0 / glut_s.get_fps();
      y_off = fmod(y_off, 100.0) * 2.0 / glut_s.get_fps();
      cam.set_spin_rot(Vec3d(-y_off, x_off, 0.0) * cam.get_rotation() *
                           cam.get_spin_rot(),
                       sqrt(x_off * x_off + y_off * y_off));
      break;
    }
    }
  }
  glut_s.start_x = x;
  glut_s.start_y = y;
  glutPostRedisplay();
}

void menu_cb(int item) { keyboard_cb(item, 0, 0); }

void glut_idle_cb()
{
  static Timer tmr;
  tmr.sleep_until_finished();
  tmr.set_timer(1 / glut_s.get_fps());
  if (glut_s.scen.animate())
    glutPostRedisplay();
}

void mouse_cb(int button, int state, int x, int y)
{
  switch (button) {
  case GLUT_LEFT_BUTTON:
    if (state == GLUT_DOWN) {
      glut_s.l_button = GLUT_DOWN;
      glut_s.start_x = x;
      glut_s.start_y = y;
    }
    // glutIdleFunc(spinDisplay);
    if (state == GLUT_UP) {
      motion_cb(x, y);
      glut_s.l_button = GLUT_UP;
    }
    break;
  default:
    break;
  }
}

void display_cb()
{
  glut_s.write();
  glutSwapBuffers();
}

void reshape_cb(int w, int h)
{
  glut_s.width = w;
  glut_s.height = h;
  glViewport(0, 0, (GLsizei)w, (GLsizei)h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  double wth = glut_s.scen.get_width();
  Camera &cam = glut_s.scen.cur_camera();
  if (glut_s.persp_proj)
    gluPerspective(30, (float)w / h, cam.get_cut_dist(), 100 * wth);
  else
    glOrtho(-wth / 2 * w / h, wth / 2 * w / h, -wth / 2, wth / 2,
            cam.get_cut_dist(), 100 * wth);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glutPostRedisplay();
}
