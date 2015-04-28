/*
   Copyright (c) 2003-2008, Adrian Rossiter

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

/* \file disp_poly.cc
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

#include "../base/antiprism.h"
#include "disp_poly_gl.h"
#include "gl_writer.h"



extern unsigned char stippleMask[17][128];

void gl_set_material(col_val col = vec3d(.5,.5,.5), bool trans=true,
      int sides=GL_FRONT_AND_BACK);
 
void gl_set_material(col_val col, bool trans, int sides)
{
   vec4d cv = col.get_vec4d();
   GLfloat f_specular[] = { (GLfloat)cv[0]/2, (GLfloat)cv[1]/2,
                            (GLfloat)cv[2]/2, 1.0 };
   GLfloat f_diffuse[] =  { (GLfloat)cv[0], (GLfloat)cv[1],
                            (GLfloat)cv[2], 1.0 };
   GLfloat f_shininess[] = { 100.0 };
   glMaterialfv(sides, GL_SPECULAR, f_specular);
   glMaterialfv(sides, GL_DIFFUSE, f_diffuse);
   glMaterialfv(sides, GL_SHININESS, f_shininess);
   if(trans)
      glPolygonStipple(stippleMask[int(cv[3]*16+0.5)]);
}

 

void disp_poly_gl::gl_verts(const scene &scen)
{
   double v_rad = get_vert_rad();
   double extra = 10*v_rad/scen.get_width();
   int long_div = int(11*(1+extra));
   int lat_div = int(7*(1+extra));
   GLUquadric *quad = gluNewQuadric();
   GLuint sph = glGenLists(1);
   glNewList(sph, GL_COMPILE);
   gluSphere(quad, v_rad, long_div, lat_div);
   glEndList();
   
   const vector<vec3d> &verts = disp_geom.verts();
   for(unsigned int i=0; i<verts.size(); i++) { 
      col_val col = disp_geom.get_v_col((int)i);
      if(col.is_idx())
         col = get_v_clrng().get_col(col.get_idx());
      if(!col.is_val())
         col = get_def_v_col(); // use default
      if(col.is_inv())
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

static void edge_cyl_trans(vec3d p1, vec3d p2)
{
   vec3d p1to2 = p2 - p1;
   glTranslated(p1[0], p1[1], p1[2]);
   mat3d rot = mat3d::rot(vec3d(0,0,1), p1to2);
   glMultMatrixd(rot.transpose().get_m());
   glScaled(1, 1, p1to2.mag());
} 


void disp_poly_gl::gl_edges(const scene &scen)
{
   double e_rad = get_edge_rad();
   double extra = 10*e_rad/scen.get_width();
   int long_div = int(11*(1+extra));
   int lat_div = 1;
   
   GLUquadric *quad = gluNewQuadric();
   GLuint cyl = glGenLists(1);
   glNewList(cyl, GL_COMPILE);
   gluCylinder(quad, e_rad, e_rad, 1, long_div, lat_div);
   glEndList();
   
   const vector<vec3d> &verts = disp_geom.verts();
   const vector<vector<int> > &edges = disp_geom.edges();
   for(unsigned int i=0; i<edges.size(); i++) { 
      col_val col = disp_geom.get_e_col((int)i);
      if(col.is_idx())
         col = get_e_clrng().get_col(col.get_idx());
      if(!col.is_val())
         col = get_def_e_col(); // use default
      if(col.is_inv())
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


void disp_poly_gl::gl_faces(const scene &)
{
   const vector<vec3d> &verts = disp_geom.verts();
   const vector<vector<int> > &faces = disp_geom.faces();
   for(unsigned int i=0; i<faces.size(); i++) { 
      if(faces[i].size()<3)
         continue;
      col_val col = disp_geom.get_f_col((int)i);
      if(col.is_idx())
         col = get_f_clrng().get_col(col.get_idx());
      if(!col.is_val())
         col = get_def_f_col(); // use default
      if(col.is_inv())
         continue;
      if(show_orientation) {
         gl_set_material(col_val(1.0,1.0,1.0), get_elem_trans(), GL_FRONT);
         gl_set_material(col_val(0.0,0.0,0.0), get_elem_trans(), GL_BACK);
      }
      else {
         gl_set_material(col, get_elem_trans());
      }

      glBegin(GL_POLYGON);
      vec3d norm = face_norm(verts,faces[i]);
      glNormal3dv(norm.get_v());
      for(unsigned int j=0; j<faces[i].size(); j++)
         glVertex3dv(verts[faces[i][j]].get_v());
      glEnd();
   }
}

disp_poly_gl::disp_poly_gl() : disp_poly(), show_orientation(false)
{
}

void disp_poly_gl::gl_geom(const scene &scen)
{
   if(v().get_show())
      gl_verts(scen);
   if(f().get_show())
      gl_faces(scen);
   if(e().get_show())
      gl_edges(scen);
}



static void draw_text(char *str, double font_sz, vec3d pos, int halign=1,
      int valign=1, bool bill=true,
      const anti_StrokeFont*font=ANTI_STROKE_ROMAN)
{  
   const float drop = 33.33;
   const float height = 152.38;
   float scale = font_sz/height;
   float width = antiStrokeLength(font, (unsigned char *)str);
  
   float off_x=0, off_y=0;
   switch(halign) {
      case 0: // left
         off_x=0;
         break;
      case 1: // centre
         off_x -= 0.5 * width;
         break;
      case 2: // right
         off_x -= width;
         break;
   }
   
   switch(valign) {
      case 0: // bottom
         off_y -= -drop;
         break;
      case 1: // centre
         off_y -= height/2 - drop;
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
   
   if(bill) {
      float m[16];
      glGetFloatv(GL_MODELVIEW_MATRIX , m);
      for(int i=0; i<3; i++ ) 
         for(int j=0; j<3; j++ ) {
            if ( i==j )
               m[i*4+j] = 1;
            else
               m[i*4+j] = 0;
         }
      glLoadMatrixf(m);
   }
  
   glScalef(scale, scale, scale);
   glTranslatef(off_x, off_y, 0);
   glNormal3f(0,0,1);
   for (char *c=str; *c ; c++)
      antiStrokeCharacter(font, *c);
   glPopMatrix();
}



static void gl_write_label(char *label, vec3d pos, const camera &cam)
{
   draw_text(label, cam.get_text_sz(pos), pos);
}


void disp_num_labels_gl::gl_verts(const scene &scen)
{
   col_geom_v &geom = sc_geom->get_geom();
   gl_set_material(get_label_col(v().get_col()), get_elem_trans());
   char label[64];
   const vector<vec3d> &verts = geom.verts();
   for(unsigned int i=0; i<verts.size(); i++) {
      if(geom.get_v_col((int)i).is_inv())
         continue;
      sprintf(label, "%u", i);
      gl_write_label(label, sc_geom->get_v_label_pos(i), scen.cur_camera());
   }
}

void disp_num_labels_gl::gl_edges(const scene &scen)
{
   col_geom_v &geom = sc_geom->get_geom();
   gl_set_material(get_label_col(e().get_col()), get_elem_trans());
   char label[64];
   const vector<vector<int> > &edges = geom.edges();
   for(unsigned int i=0; i<edges.size(); i++) { 
      if(geom.get_e_col((int)i).is_inv())
         continue;
      sprintf(label, "%u", i);
      gl_write_label(label, sc_geom->get_e_label_pos(i), scen.cur_camera());
   }
}


void disp_num_labels_gl::gl_faces(const scene &scen)
{
   col_geom_v &geom = sc_geom->get_geom();
   gl_set_material(get_label_col(f().get_col()), get_elem_trans());
   char label[64];
   const vector<vector<int> > &faces = geom.faces();
   for(unsigned int i=0; i<faces.size(); i++) { 
      if(geom.get_f_col((int)i).is_inv())
         continue;
      sprintf(label, "%u", i);
      gl_write_label(label, sc_geom->get_f_label_pos(i), scen.cur_camera());
   }
}

void disp_num_labels_gl::gl_geom(const scene &scen)
{
   if(v().get_show())
      gl_verts(scen);
   if(f().get_show())
      gl_faces(scen);
   if(e().get_show())
      gl_edges(scen);
}


void disp_sym_gl::gl_geom(const scene &scen)
{
   disp_poly_gl::gl_geom(scen);
}

 
