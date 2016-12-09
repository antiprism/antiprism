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

/*
   Name: gl_writer.h
   Description: write a OpenGL output
   Project: Antiprism - http://www.antiprism.com
*/

#include <vector>

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

#include "gl_writer.h"

using std::vector;

using namespace anti;

void gl_writer::init(const Scene &)
{
  GLfloat light_position[] = {-100.0, 200.0, 1000.0, 0.0};
  GLfloat light_scene_ambient[] = {0.3, 0.3, 0.3, 1.0};
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glClearDepth(1);
  glShadeModel(GL_SMOOTH);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, light_scene_ambient);

  glLightfv(GL_LIGHT0, GL_POSITION, light_position);
  // glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  // glEnable(GL_COLOR_MATERIAL);
  // glEnable(GL_POLYGON_OFFSET_FILL);   // remove overlap
  glEnable(GL_POLYGON_STIPPLE); // for transparency
  glEnable(GL_NORMALIZE);
  glEnable(GL_DEPTH_TEST);
}

void gl_writer::geometry_objects(const Scene &scen)
{
  vector<SceneGeometry>::const_iterator geo;
  for (geo = scen.get_geoms().begin(); geo != scen.get_geoms().end(); ++geo) {
    vector<GeometryDisplay *>::const_iterator disp;
    for (disp = geo->get_disps().begin(); disp != geo->get_disps().end();
         ++disp)
      (*disp)->gl_geom(scen);
    if (geo->get_label())
      geo->get_label()->gl_geom(scen);
    if (geo->get_sym())
      geo->get_sym()->gl_geom(scen);
  }
}

void gl_writer::cur_camera(const Camera &cam)
{
  glLoadIdentity();
  // adjust distance to give the same view as in other export formats
  glTranslatef(0.0, 0.0, -1.57 * cam.get_distance());
  Vec3d lookat = cam.get_lookat();
  glTranslated(-lookat[0], -lookat[1], -lookat[2]);
  Vec3d centre = cam.get_centre();
  glTranslated(centre[0], centre[1], centre[2]);
  double rot[16];
  to_gl_matrix(rot, cam.get_rotation() * cam.get_spin_rot());

  glMultMatrixd(rot);
  glTranslated(-centre[0], -centre[1], -centre[2]);
}

void gl_writer::write(const Scene &scen)
{
  Vec3d bg = scen.get_bg_col().get_Vec3d();
  glClearColor(bg[0], bg[1], bg[2], 1);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  cur_camera(scen.cur_camera());

  glPushMatrix();
  geometry_objects(scen);
  glPopMatrix();
}
