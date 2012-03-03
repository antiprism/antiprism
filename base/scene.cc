/*
   Copyright (c) 2008, Adrian Rossiter

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

/* \file scene.cc
   Scene classes for a camera and geometry.
*/


#include "scene.h"
#include "utils.h"
#include "vec_utils.h"
#include "geom_utils.h"


// -------------------------------------------------------------- 
// geom_disp

double geom_disp::get_vert_rad() const
{
   if(sc_geom && v().get_size()==rad_default)
      return sc_geom->get_v_ball_rad()/20;
   else if(sc_geom && v().get_size()==rad_ball)
      return sc_geom->get_v_ball_rad();
   else
      return v().get_size();
}

double geom_disp::get_edge_rad() const
{
   if(e().get_size()==rad_default)
      return get_vert_rad()/1.5;
   else
      return e().get_size();
}


vec3d geom_disp::get_label_pos(const vec3d &point, double elem_sz)
{
   vec3d off_set(point - sc_geom->get_centre());
   if(off_set.mag2()<1e-20)
      off_set = vec3d(0,0,1e-10);
   off_set *= 1.07 + elem_sz/off_set.mag();
   return sc_geom->get_centre() + off_set;
}

vec3d geom_disp::get_v_label_pos(int idx)
{
   return get_label_pos(sc_geom->get_geom().verts(idx), v().get_size());
}

vec3d geom_disp::get_e_label_pos(int idx)
{
   return get_label_pos(sc_geom->get_geom().edge_cent(idx), e().get_size());
}


vec3d geom_disp::get_f_label_pos(int idx)
{
   // use edge radius for face
   return get_label_pos(sc_geom->get_geom().face_cent(idx), e().get_size());
}


geom_disp_label::geom_disp_label() :
   label_light(false), label_invert(false)
{ 
   v().set_show(0);
   e().set_show(0);
   f().set_show(0);
}


col_val geom_disp_label::get_label_col(col_val col) const
{
   col_val lab_col = col;
   if(label_invert)
      lab_col.set_complement();
   if(label_light)
      lab_col.set_brightness(0.5-label_invert);
   return lab_col;
}
 

// -------------------------------------------------------------- 
// scene_geom

scene_geom::scene_geom(const geom_if &geo): label(0), sym(0)
{
   set_geom(geo);
}

scene_geom::~scene_geom()
{
   for(unsigned int i=0; i<disps.size(); i++)
      delete disps[i];
   delete label;
   delete sym;
}


scene_geom::scene_geom(const scene_geom &sc_geo): scene_item(sc_geo),
   bound_sph(sc_geo.bound_sph), width(sc_geo.width), centre(sc_geo.centre),
   v_ball_rad(sc_geo.v_ball_rad), geom(sc_geo.geom), label(0), sym(0)
{
   for(unsigned int i=0; i<sc_geo.disps.size(); i++)
      add_disp(*sc_geo.disps[i]);
   if(sc_geo.label)
      set_label(*sc_geo.label);
   if(sc_geo.sym)
      set_sym(*sc_geo.sym);
}


scene_geom &scene_geom::operator=(const scene_geom &sc_geo)
{
   if(this!=&sc_geo) {
      this->scene_item::operator =(sc_geo);

      bound_sph = sc_geo.bound_sph;
      width = sc_geo.width;
      centre = sc_geo.centre;
      v_ball_rad = sc_geo.v_ball_rad;
      geom = sc_geo.geom;
      for(unsigned int i=0; i<sc_geo.disps.size(); i++)
         add_disp(*sc_geo.disps[i]);
      if(sc_geo.label)
         set_label(*sc_geo.label);
      else {
         delete label;
         label=0;
      }
      if(sc_geo.sym)
         set_sym(*sc_geo.sym);
      else {
         delete sym;
         sym=0;
      }

   }
   return *this;
}


void scene_geom::geom_changed()
{
   bound_sph = bound_sphere(geom.verts(),
            (scen) ? scen->get_inf_dist() : -1);
   width = 2*bound_sph.get_radius();
   centre = bound_sph.get_centre();
   v_ball_rad = minimum_distance(geom, width/1000)/2.0;

   for(unsigned int i=0; i<disps.size(); i++)
      disps[i]->geom_changed();

   if(label)
      label->geom_changed();

   if(sym)
      sym->geom_changed();
}


void scene_geom::set_geom(const geom_if &geo)
{ 
   geom=geo;
   geom_changed();
}


void scene_geom::add_disp(geom_disp &disp)
{ 
   disps.push_back(disp.clone());
   disps.back()->set_scene_geom(this);
   disps.back()->geom_changed();
}


bool scene_geom::delete_disp(int idx)
{
   if(idx<0 || idx>=(int)disps.size())
      return false;
   delete disps[idx];
   disps.erase(disps.begin()+idx);
   return true;
}

geom_disp_label *scene_geom::get_label() const
{
   return label;
}


void scene_geom::set_label(const geom_disp_label &lab)
{ 
   delete label;
   label = dynamic_cast<geom_disp_label *>(lab.clone());
   label->set_scene_geom(this);
}

geom_disp *scene_geom::get_sym() const
{
   return sym;
}


void scene_geom::set_sym(const geom_disp &sy)
{ 
   delete sym;
   sym = dynamic_cast<geom_disp *>(sy.clone());
   sym->set_scene_geom(this);
}


vec3d scene_geom::get_v_label_pos(int idx) const
{
   if(disps.size())
      return(disps[0]->get_v_label_pos(idx));
   else
      return geom.verts(idx);
}


vec3d scene_geom::get_e_label_pos(int idx) const
{
   if(disps.size())
      return(disps[0]->get_e_label_pos(idx));
   else
      return geom.edge_cent(idx);
}

vec3d scene_geom::get_f_label_pos(int idx) const
{
   if(disps.size())
      return(disps[0]->get_f_label_pos(idx));
   else
      return geom.face_cent(idx);
} 

int scene_geom::animate()
{
   int num_changes = 0;
   for(unsigned int i=0; i<disps.size(); i++)
      num_changes += disps[i]->animate();
   return num_changes;
}



// ------------------------------------------------------------------- 
// camera

camera::camera(): distance(0), width(0), spin_inc(0)
{
   set_rotation();
   set_persp();
   set_cut_dist();
   stop_spinning();
}

double camera::get_width() const
{ 
   double w = width;
   if(width==0 && scen) {
      w = scen->get_width();
      if(scen->get_bound_sph().get_cut_off_cnt())  // add width when inf verts
         w *= 2;
   }

   return w>epsilon? w : epsilon;
}


vec3d camera::get_centre() const
{ 
   if(centre.is_set())
      return centre;
   else if (scen)
      return scen->get_bound_sph().get_centre();
   else
      return vec3d(0,0,0);
}

void camera::set_cut_dist(double dist)
{
   cut_dist = dist;
}

double camera::get_cut_dist() const
{
   double width = scen ? scen->get_bound_sph().get_width() : 1.0;
   return (cut_dist>=0.1*width) ? cut_dist : 0.1*width;
}

double camera::get_text_sz(vec3d pos) const
{
   double dist;
   if(persp)
      dist = 2*get_distance() - (label_rot*pos)[2] + get_lookat()[2];
   else
      dist = 2*scen->get_width();
   return dist/60;
}

void camera::set_label_rot()
{
   if(scen)
      label_rot = mat3d::transl(scen->get_centre()) *
                  get_rotation() *
                  mat3d::transl(-scen->get_centre());
}

int camera::animate()
{
   int num_changes = 0;
   if(is_spinning()) {
      inc_spin_rot();
      num_changes += 1;
   }

   return num_changes;
}

// ------------------------------------------------------------------- 
// scene

scene::scene()
{
   set_bg_col();
   set_inf_dist(DEF_CAMERA_INF_DIST);
   add_camera();
   cur_camera().set_name("default");
}

scene::scene(const scene &scen): scene_item(scen)
{
   anim_timer = scen.anim_timer;
   cycle_rate = scen.cycle_rate;
   bg_col = scen.bg_col;
   inf_dist = scen.inf_dist;
   bound_sph = scen.bound_sph;

   geoms = scen.geoms;
   for(unsigned int i=0; i<geoms.size(); i++)
      geoms[i].set_scene(this);
      
   cur_cam_num = scen.cur_cam_num;
   cams = scen.cams;
   for(unsigned int i=0; i<cams.size(); i++)
      cams[i].set_scene(this);
}


scene &scene::operator=(const scene &scen)
{
   if(this != &scen) {
      anim_timer = scen.anim_timer;
      cycle_rate = scen.cycle_rate;
      bg_col = scen.bg_col;
      inf_dist = scen.inf_dist;
      bound_sph = scen.bound_sph;
      
      geoms = scen.geoms;
      for(unsigned int i=0; i<geoms.size(); i++)
         geoms[i].set_scene(this);

      cur_cam_num = scen.cur_cam_num;
      cams = scen.cams;
      for(unsigned int i=0; i<cams.size(); i++)
         cams[i].set_scene(this);
   }
   return *this;
}


void scene::add_geom(const scene_geom sc_geom)
{
   geoms.push_back(sc_geom);
   geoms.back().set_scene(this);
   bound_sph.add_b_sphere(geoms.back().get_bound_sph());
}


void scene::add_geom(const geom_if &geom)
{ 
   geoms.push_back(scene_geom(geom));
   geoms.back().set_scene(this);
   bound_sph.add_b_sphere(geoms.back().get_bound_sph());
}


bool scene::delete_geom(int idx)
{
   if(idx<0 || idx>=(int)geoms.size())
      return false;
   geoms.erase(geoms.begin()+idx);

   bound_sph = bound_sphere();
   bound_sph.set_cut_off(inf_dist);
   vector<scene_geom>::iterator vi;
   for(vi=geoms.begin(); vi!=geoms.end(); ++vi)
      bound_sph.add_b_sphere(vi->get_bound_sph());
   return true;
}

const camera &scene::cur_camera() const
{
   return cams[cur_cam_num];
}

camera &scene::cur_camera()
{
   return cams[cur_cam_num];
}

string scene::get_camera_name(int idx) const
{
   char str[MSG_SZ];
   snprintf(str, MSG_SZ, "%s_camera_%d", cams[idx].get_name().c_str(), idx);
   return str;
}

void scene::add_camera(const camera &cam)
{
   cams.push_back(cam);
   cams.back().set_scene(this);
   cur_cam_num = cams.size()-1;
}

bool scene::delete_camera(int idx)
{
   if(idx<0 || idx>=(int)cams.size())
      return false; // camera out of range
   if(cams.size()==1)
      return false; // can't remove last camera
   // Adjust current camera number back if beyond idx
   if(cur_cam_num>0 && cur_cam_num > idx)
      cur_cam_num--;
   cams.erase(cams.begin()+idx);
   return true;
}
  
        
int scene::animate()
{
   int num_changes = 0;
   if(anim_timer.finished()) {
      vector<camera>::iterator cam;
      for(cam=cams.begin(); cam!=cams.end(); ++cam)
         num_changes += cam->animate();
      
      vector<scene_geom>::iterator sc_geo;
      for(sc_geo=geoms.begin(); sc_geo!=geoms.end(); ++sc_geo)
         num_changes += sc_geo->animate();
      
   }
   anim_timer.inc_timer(1.0/cycle_rate);

   return num_changes;
} 

 



