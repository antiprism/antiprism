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

/* \file vw_glut.h
   \brief antiview - interface with GLUT
*/

#ifndef VW_GLUT_H
#define VW_GLUT_H

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

#include "../base/antiprism.h"
#include "gl_writer.h"

#define FOUND_NO_GLUT 0
#define FOUND_GLUT 1
#define FOUND_OPENGLUT 2
#define FOUND_FREEGLUT 3
#define FOUND_FLTKGLUT 4

#if GLUT_TYPE == FOUND_GLUT
#ifdef HAVE_GL_GLUT_H
#include <GL/glut.h>
#elif defined HAVE_GLUT_GLUT_H
#include <GLUT/glut.h>
#endif
#elif GLUT_TYPE == FOUND_OPENGLUT
#include <GL/openglut.h>
#elif GLUT_TYPE == FOUND_FREEGLUT
#include <GL/freeglut.h>
#elif GLUT_TYPE == FOUND_FLTKGLUT
#include <FL/glut.H>
#undef Status
#else
#error "No GLUT library has been specified"
#endif

using namespace anti;

void keyboard_cb(unsigned char key, int x, int y);
void special_cb(int key, int x, int y);
void menu_cb(int value);
void mouse_cb(int button, int state, int x, int y);
void motion_cb(int x, int y);
void display_cb();
void reshape_cb(int w, int h);
void glut_idle_cb();

void draw_text(char *str, double font_sz, Vec3d pos, int halign = 1,
               int valign = 1, bool bill = true,
               void *font = GLUT_STROKE_ROMAN);

class glut_state {
private:
  int cur_op;
  double fps; // frames per second
  Camera saved_camera;
  gl_writer gl_wrtr;

public:
  enum { tr_none = 0, tr_rot, tr_xy, tr_z, tr_slice, tr_spin };
  int width;
  int height;
  int start_x;
  int start_y;
  int l_button;
  double trans_scale;
  double slice_scale;
  bool persp_proj;
  int sym_disp_type;
  int transparency_disp_type;
  Scene scen;

  glut_state()
      : cur_op(tr_rot), fps(50), start_x(0), start_y(0), l_button(GLUT_UP),
        trans_scale(0.005), slice_scale(0.003), persp_proj(true),
        sym_disp_type(0), transparency_disp_type(0)
  {
  }
  void make_menu();
  void toggle_projection() { persp_proj = !persp_proj; }
  void set_projection(bool per = true) { persp_proj = per; }
  void reshape() { reshape_cb(width, height); }
  int get_op() { return cur_op; }
  void set_op(int op);
  double get_fps() { return fps; }
  void init() { gl_wrtr.init(scen); }
  void write() { gl_wrtr.write(scen); }
  void save_camera() { saved_camera = scen.cur_camera(); }
  void restore_camera() { scen.cur_camera() = saved_camera; }
  void change_sym_disp();
  void change_transparency_disp();
};

#endif // VW_GLUT_H
