/*
   Copyright (c) 2008-2016, Adrian Rossiter

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

#include <string>
#include <vector>

#include "private_misc.h"
#include "scene.h"
#include "utils.h"
#include "vec_utils.h"

using std::string;
using std::vector;

namespace anti {

// --------------------------------------------------------------
// GeometryDisplay

double GeometryDisplay::get_vert_rad() const
{
  if (sc_geom && elem(VERTS).get_size() == rad_default)
    return sc_geom->get_v_ball_rad() / 20;
  else if (sc_geom && elem(VERTS).get_size() == rad_ball)
    return sc_geom->get_v_ball_rad();
  else
    return elem(VERTS).get_size();
}

double GeometryDisplay::get_edge_rad() const
{
  if (elem(EDGES).get_size() == rad_default)
    return get_vert_rad() / 1.5;
  else
    return elem(EDGES).get_size();
}

Vec3d GeometryDisplay::get_label_pos(const Vec3d &point, double elem_sz)
{
  Vec3d off_set(point - sc_geom->get_centre());
  if (off_set.len2() < 1e-20)
    off_set = Vec3d(0, 0, 1e-10);
  off_set *= 1.07 + elem_sz / off_set.len();
  return sc_geom->get_centre() + off_set;
}

Vec3d GeometryDisplay::get_v_label_pos(int idx)
{
  return get_label_pos(sc_geom->get_geom().verts(idx), elem(VERTS).get_size());
}

Vec3d GeometryDisplay::get_e_label_pos(int idx)
{
  return get_label_pos(sc_geom->get_geom().edge_cent(idx),
                       elem(EDGES).get_size());
}

Vec3d GeometryDisplay::get_f_label_pos(int idx)
{
  // use edge radius for face
  return get_label_pos(sc_geom->get_geom().face_cent(idx),
                       elem(EDGES).get_size());
}

GeometryDisplayLabel::GeometryDisplayLabel()
    : label_light(false), label_invert(false)
{
  elem(VERTS).set_show(0);
  elem(EDGES).set_show(0);
  elem(FACES).set_show(0);
}

Color GeometryDisplayLabel::get_label_col(Color col) const
{
  Color lab_col = col;
  if (label_invert)
    lab_col.set_complement();
  if (label_light)
    lab_col.set_brightness(0.5 - label_invert);
  return lab_col;
}

// --------------------------------------------------------------
// SceneGeometry

SceneGeometry::SceneGeometry(const Geometry &geo) : label(nullptr), sym(nullptr)
{
  set_geom(geo);
}

SceneGeometry::~SceneGeometry()
{
  for (auto &disp : disps)
    delete disp;
  delete label;
  delete sym;
}

SceneGeometry::SceneGeometry(const SceneGeometry &sc_geo)
    : SceneItem(sc_geo), bound_sph(sc_geo.bound_sph), width(sc_geo.width),
      centre(sc_geo.centre), v_ball_rad(sc_geo.v_ball_rad), geom(sc_geo.geom),
      label(nullptr), sym(nullptr)
{
  for (auto disp : sc_geo.disps)
    add_disp(*disp);
  if (sc_geo.label)
    set_label(*sc_geo.label);
  if (sc_geo.sym)
    set_sym(*sc_geo.sym);
}

SceneGeometry &SceneGeometry::operator=(const SceneGeometry &sc_geo)
{
  if (this != &sc_geo) {
    this->SceneItem::operator=(sc_geo);

    bound_sph = sc_geo.bound_sph;
    width = sc_geo.width;
    centre = sc_geo.centre;
    v_ball_rad = sc_geo.v_ball_rad;
    geom = sc_geo.geom;
    for (auto disp : sc_geo.disps)
      add_disp(*disp);
    if (sc_geo.label)
      set_label(*sc_geo.label);
    else {
      delete label;
      label = nullptr;
    }
    if (sc_geo.sym)
      set_sym(*sc_geo.sym);
    else {
      delete sym;
      sym = nullptr;
    }
  }
  return *this;
}

void SceneGeometry::geom_changed()
{
  bound_sph = BoundSphere(geom.verts(), (scen) ? scen->get_inf_dist() : -1);
  width = 2 * bound_sph.get_radius();
  centre = bound_sph.get_centre();
  v_ball_rad = get_min_vert_to_vert_dist(geom.verts(), width / 1000) / 2.0;

  for (auto &disp : disps)
    disp->geom_changed();

  if (label)
    label->geom_changed();

  if (sym)
    sym->geom_changed();
}

void SceneGeometry::set_geom(const Geometry &geo)
{
  geom = geo;
  geom_changed();
}

void SceneGeometry::add_disp(GeometryDisplay &disp)
{
  disps.push_back(disp.clone());
  disps.back()->set_scene_geom(this);
  disps.back()->geom_changed();
}

bool SceneGeometry::delete_disp(int idx)
{
  if (idx < 0 || idx >= (int)disps.size())
    return false;
  delete disps[idx];
  disps.erase(disps.begin() + idx);
  return true;
}

GeometryDisplayLabel *SceneGeometry::get_label() const { return label; }

void SceneGeometry::set_label(const GeometryDisplayLabel &lab)
{
  delete label;
  label = dynamic_cast<GeometryDisplayLabel *>(lab.clone());
  label->set_scene_geom(this);
}

GeometryDisplay *SceneGeometry::get_sym() const { return sym; }

void SceneGeometry::set_sym(const GeometryDisplay &sy)
{
  delete sym;
  sym = dynamic_cast<GeometryDisplay *>(sy.clone());
  sym->set_scene_geom(this);
}

Vec3d SceneGeometry::get_v_label_pos(int idx) const
{
  if (disps.size())
    return (disps[0]->get_v_label_pos(idx));
  else
    return geom.verts(idx);
}

Vec3d SceneGeometry::get_e_label_pos(int idx) const
{
  if (disps.size())
    return (disps[0]->get_e_label_pos(idx));
  else
    return geom.edge_cent(idx);
}

Vec3d SceneGeometry::get_f_label_pos(int idx) const
{
  if (disps.size())
    return (disps[0]->get_f_label_pos(idx));
  else
    return geom.face_cent(idx);
}

int SceneGeometry::animate()
{
  int num_changes = 0;
  for (auto &disp : disps)
    num_changes += disp->animate();
  return num_changes;
}

// -------------------------------------------------------------------
// Camera

Camera::Camera() : distance(0), width(0), spin_inc(0)
{
  set_rotation();
  set_persp();
  set_cut_dist();
  stop_spinning();
}

double Camera::get_width() const
{
  double w = width;
  if (width == 0 && scen) {
    w = scen->get_width();
    if (scen->get_bound_sph().get_cut_off_cnt()) // add width when inf verts
      w *= 2;
  }

  return w > epsilon ? w : epsilon;
}

Vec3d Camera::get_centre() const
{
  if (centre.is_set())
    return centre;
  else if (scen)
    return scen->get_bound_sph().get_centre();
  else
    return Vec3d(0, 0, 0);
}

void Camera::set_cut_dist(double dist) { cut_dist = dist; }

double Camera::get_cut_dist() const
{
  double width = scen ? scen->get_bound_sph().get_width() : 1.0;
  return (cut_dist >= 0.1 * width) ? cut_dist : 0.1 * width;
}

double Camera::get_text_sz(Vec3d pos) const
{
  double dist;
  if (persp)
    dist = 2 * get_distance() - (label_rot * pos)[2] + get_lookat()[2];
  else
    dist = 2 * scen->get_width();
  return dist / 60;
}

void Camera::set_label_rot()
{
  if (scen)
    label_rot = Trans3d::translate(scen->get_centre()) * get_rotation() *
                Trans3d::translate(-scen->get_centre());
}

int Camera::animate()
{
  int num_changes = 0;
  if (is_spinning()) {
    inc_spin_rot();
    num_changes += 1;
  }

  return num_changes;
}

// -------------------------------------------------------------------
// Scene

Scene::Scene(): cycle_rate(-1 /* instant updates */)
{
  set_bg_col();
  set_inf_dist(DEF_CAMERA_INF_DIST);
  add_camera();
  cur_camera().set_name("default");
}

Scene::Scene(const Scene &scen) : SceneItem(scen)
{
  anim_timer = scen.anim_timer;
  cycle_rate = scen.cycle_rate;
  bg_col = scen.bg_col;
  inf_dist = scen.inf_dist;
  bound_sph = scen.bound_sph;

  geoms = scen.geoms;
  for (auto &geom : geoms)
    geom.set_scene(this);

  cur_cam_num = scen.cur_cam_num;
  cams = scen.cams;
  for (auto &cam : cams)
    cam.set_scene(this);
}

Scene &Scene::operator=(const Scene &scen)
{
  if (this != &scen) {
    anim_timer = scen.anim_timer;
    cycle_rate = scen.cycle_rate;
    bg_col = scen.bg_col;
    inf_dist = scen.inf_dist;
    bound_sph = scen.bound_sph;

    geoms = scen.geoms;
    for (auto &geom : geoms)
      geom.set_scene(this);

    cur_cam_num = scen.cur_cam_num;
    cams = scen.cams;
    for (auto &cam : cams)
      cam.set_scene(this);
  }
  return *this;
}

void Scene::add_geom(const SceneGeometry sc_geom)
{
  geoms.push_back(sc_geom);
  geoms.back().set_scene(this);
  bound_sph.add_b_sphere(geoms.back().get_bound_sph());
}

void Scene::add_geom(const Geometry &geom)
{
  geoms.push_back(SceneGeometry(geom));
  geoms.back().set_scene(this);
  bound_sph.add_b_sphere(geoms.back().get_bound_sph());
}

bool Scene::delete_geom(int idx)
{
  if (idx < 0 || idx >= (int)geoms.size())
    return false;
  geoms.erase(geoms.begin() + idx);

  bound_sph = BoundSphere();
  bound_sph.set_cut_off(inf_dist);
  vector<SceneGeometry>::iterator vi;
  for (vi = geoms.begin(); vi != geoms.end(); ++vi)
    bound_sph.add_b_sphere(vi->get_bound_sph());
  return true;
}

const Camera &Scene::cur_camera() const { return cams[cur_cam_num]; }

Camera &Scene::cur_camera() { return cams[cur_cam_num]; }

string Scene::get_camera_name(int idx) const
{
  char str[MSG_SZ];
  snprintf(str, MSG_SZ, "%s_camera_%d", cams[idx].get_name().c_str(), idx);
  return str;
}

void Scene::add_camera(const Camera &cam)
{
  cams.push_back(cam);
  cams.back().set_scene(this);
  cur_cam_num = cams.size() - 1;
}

bool Scene::delete_camera(int idx)
{
  if (idx < 0 || idx >= (int)cams.size())
    return false; // camera out of range
  if (cams.size() == 1)
    return false; // can't remove last camera
  // Adjust current camera number back if beyond idx
  if (cur_cam_num > 0 && cur_cam_num > idx)
    cur_cam_num--;
  cams.erase(cams.begin() + idx);
  return true;
}

int Scene::animate()
{
  int num_changes = 0;
  if (anim_timer.finished()) {
    vector<Camera>::iterator cam;
    for (cam = cams.begin(); cam != cams.end(); ++cam)
      num_changes += cam->animate();

    vector<SceneGeometry>::iterator sc_geo;
    for (sc_geo = geoms.begin(); sc_geo != geoms.end(); ++sc_geo)
      num_changes += sc_geo->animate();
  }
  anim_timer.inc_timer(1.0 / cycle_rate);

  return num_changes;
}

} // namespace anti
