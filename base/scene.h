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

/**\file scene.h
   Classes for a scene, camera, geometry, labelling, etc.
*/


#ifndef SCENE_H
#define SCENE_H

#include <string>
#include "geom.h"
#include "timing.h"
#include "coloring.h"
#include "symmetry.h"
#include "bbox.h"

using std::string;

class geom_disp;

class scene;
class camera;
class scene_geom;

///Class for basic display properties of an element type.
class elem_disp
{
   private:
      bool show;
      double size;
      col_val col;

   public:
      ///Constructor
      elem_disp(): show(true), size(0) {}

      ///Set visibility for the element type.
      /**\param sho \c true for visible, \c false for hidden. */
      void set_show(bool sho=true) { show = sho; }

      ///Get visibility of the element type.
      /**\return \c true if visible, \c false if hidden. */
      bool get_show() const { return show;}

      ///Set the size for the element type.
      /**\param sz The size for the element type. */
      void set_size(double sz=0) { size = sz; }

      ///Get the size of the element type.
      /**\return The size for the element type. */
      double get_size() const { return size; }

      ///Set the colour of the element type.
      /**\param clr the color for the element type.*/
      void set_col(col_val clr) { col = clr; }

      ///Get the colour of the element type.
      /**\return The color of the element type.*/
      col_val get_col() const { return col; }
};

//class scene_geom;

///Base class for displaying a geometry.
class geom_disp
{
   private:
      elem_disp elems[3];
      bool elem_trans;

   protected:
      
      scene_geom *sc_geom; ///The scene geometry to be displaed.
   
   public:
      enum {
         rad_default=0,  ///Default radius
         rad_ball=-1     ///Greatest radius for non-overlapping balls
      };
      
      ///Constructor
      geom_disp():  elem_trans(true), sc_geom(0) {}

      //Destructor
      virtual ~geom_disp() {}
  
      ///Get a reference to the basic vertex properties.
      /**\return A reference to the proprties. */
      elem_disp &v() { return elems[0]; }
  
      ///Get a reference to the basic vertex properties.
      /**\return A reference to the proprties. */
      const elem_disp &v() const { return elems[0]; }
  
      ///Get a reference to the basic edge properties.
      /**\return A reference to the proprties. */
      elem_disp &e() { return elems[1]; }
  
      ///Get a reference to the basic edge properties.
      /**\return A reference to the proprties. */
      const elem_disp &e() const { return elems[1]; }
  
      ///Get a reference to the basic face properties.
      /**\return A reference to the proprties. */
      elem_disp &f() { return elems[2]; }
  
      ///Get a reference to the basic face properties.
      /**\return A reference to the proprties. */
      const elem_disp &f() const { return elems[2]; }
    
      ///Get the edge rod radius.
      /**\return The edge rod radius. */
      double get_edge_rad() const;
    
      ///Get the vertex ball radius.
      /**\return The vertex ball radius. */
      double get_vert_rad() const;

      ///Set display of transparency.
      /**\param trans \c true to display any transparency,
       * \c false to display as solid. */
      void set_elem_trans(bool trans) { elem_trans = trans; }

      ///Get display of transparency.
      /**\return \c true if transparency will be displayed,
       * \c false if transparency is displayed as solid. */
      bool get_elem_trans() { return elem_trans; }

      ///Position of a label for a point.
      /**\param point the point to be labelled.
       * \param elem_sz the size of the element.
       * \return The position of the label. */
      virtual vec3d get_label_pos(const vec3d &point, double elem_sz);
      
      ///Position of a vertex label.
      /**\param idx vertex index number.
       * \return position of the label */
      virtual vec3d get_v_label_pos(int idx);

      ///Position of an edge label.
      /**\param idx edge index number.
       * \return position of the label */
      virtual vec3d get_e_label_pos(int idx);

      ///Position of a face label.
      /**\param idx face index number.
       * \return position of the label */
      virtual vec3d get_f_label_pos(int idx);
     
      ///Set the scene geometry to display.
      /**\param sc_geo the scene geometry to display. */
      void set_scene_geom(scene_geom *sc_geo) { sc_geom = sc_geo; }

      ///Clone
      /**\return A dynamically allocated copy of the object. */
      virtual geom_disp *clone() const =0;

      ///Action performed when the geometry has changed.
      virtual void geom_changed()=0;

      ///Write geometry for inclusion in VRML.
      /**\param ofile file pointer to write to.
       * \param scene the scene that the display is part of.
       * \param sig_dgts the number of significant digits for the output. */
      virtual void vrml_geom(FILE *ofile, const scene &scene, int sig_dgts)=0;

      ///Write geometry for inclusion in a POV-Ray script.
      /**\param ofile file pointer to write to.
       * \param scene the scene that the display is part of.
       * \param sig_dgts the number of significant digits for the output. */
      virtual void pov_geom(FILE *ofile, const scene &scene, int sig_dgts)=0;

      ///Display geometry as OpenGL.
      /**\param scen the scene that the display is part of.*/
      virtual void gl_geom(const scene &scen)=0;

      ///Update animated properties.
      virtual void animate() {}
 
};

class geom_disp_label: public geom_disp
{
   private:
      bool label_light;
      bool label_invert;

   public:
      ///Constructor
      geom_disp_label();

      ///Set labels to be a light colour.
      /**\param light \c true for light labels, \c false for normal labels. */
      void set_label_light(bool light) {label_light=light;}

      ///Get setting for light colour labels.
      /**\return \c true for light labels, \c false for normal labels. */
      int get_label_light() const {return label_light;}

      ///Set labels to an inverted colour.
      /**\param inv \c true for inverted, \c false for normal colour. */
      void set_label_invert(bool inv) {label_invert=inv;}

      ///Get setting for inverted colour labels.
      /**\return \c true for inverted, \c false for normal colour. */
      int get_label_invert() const {return label_invert;}

      ///Convert a base label colour into a display label colour.
      /**\param col the base colour of a label.
       * \return The display colour for the label. */
      col_val get_label_col(col_val col) const ;
};


///Base class for items in a scene.
class scene_item
{
   protected:
      scene *scen;  ///The scene the item is part of
      string name;  ///The name of the item

   public:
      ///Constructor
      scene_item(): scen(0), name("unnamed") {}

      ///Set the name
      /**\param nam the name of the item. */
      void set_name(string nam) { name=nam; }

      ///Get the name
      /**\return The name of the item. */
      string get_name() const { return name; }

      ///Set the scene that the item is part of.
      /**\param sc the scene the item is part of. */
      void set_scene(scene *sc) { scen = sc; }

      ///Get the scene that the item is part of.
      /**\return The scene the item is part of. */
      scene *get_scene() { return scen; }
};


///Class for a geometry in a scene
class scene_geom: public scene_item
{
   private:
      bound_sphere bound_sph;
      double width;
      vec3d centre;
      double v_ball_rad;
      
      col_geom_v geom;

      vector<geom_disp *> disps;
      geom_disp_label *label;
      geom_disp *sym;

   public:
      ///Constructor
      /**\param geo the geometry to display in a scene. */
      scene_geom(const geom_if &geo=col_geom_v());
      
      ///Copy constructor
      /**\param sc_geo the scene geometry to copy from. */
      scene_geom(const scene_geom &sc_geo);
      
      ///Assignment operator
      /**\param sc_geo the scene geometry to assign from. */
      scene_geom &operator=(const scene_geom &sc_geo);

      ///Destructor
      ~scene_geom();

      ///Set the geometry.
      /**\param geo the geometry to display. */
      void set_geom(const geom_if &geo=col_geom_v());

      ///Get the geometry.
      /**\return the geometry to display. */
      const geom_if &get_geom() const { return geom; }

      ///Get the geometry.
      /**\return the geometry to display. */
      col_geom_v &get_geom() { return geom; }

      ///Indicate that the geometry has changed
      /**Call this if the disply geometry has been changed through a
       * reference to the data member. */
      void geom_changed();

      ///Get the vertex ball radius (the maximum without overlap.)
      /**\return The vertex ball radius. */
      double get_v_ball_rad() const { return v_ball_rad; }

      ///Get the minimum display radius
      /**\return The minimum radius. */
      double min_rad() const { return width/1000; }

      ///Get a sphere that bounds the geometry
      /**\return A sphere that bonds the geometry. */
      const bound_sphere &get_bound_sph() const { return bound_sph; }
     
      ///Set the width of the geometry
      /**\param wdth the width to set */
      void set_width(double wdth) { width=wdth; }
     
      ///Get the width of the geometry
      /**\return The width of the geometry */
      double get_width() const { return width; }
     
      ///Set the centre of the geometry
      /**\param cent the centre to set */
      void set_centre(vec3d cent) { centre=cent; }
     
      ///Get the centre of the geometry
      /**\return The centre of the geometry. */
      vec3d get_centre() const { return centre; }
 
      ///Add a new display
      /**\param disp display to add. */
      void add_disp(geom_disp &disp);
 
      ///Delete a display
      /**\param idx index number of the display to delete. */
      bool delete_disp(int idx);

      ///Get the displays
      /**\return The vector of displays. */
      const vector<geom_disp *> &get_disps() const { return disps; }

      ///Get the displays
      /**\return The vector of displays. */
      vector<geom_disp *> &get_disps() { return disps; }
      
      ///Get the geometry label display
      /**\return The geometry label display. */
      geom_disp_label *get_label() const;
      
      ///Set the geometry label display
      /**\param lab the geometry label display to set. */
      void set_label(const geom_disp_label &lab);
      
      ///Position of a vertex label.
      /**\param idx vertex index number.
       * \return position of the label */
      vec3d get_v_label_pos(int idx) const;

      ///Position of an edge label.
      /**\param idx edge index number.
       * \return position of the label */
      vec3d get_e_label_pos(int idx) const;

      ///Position of a face label.
      /**\param idx face index number.
       * \return position of the label */
      vec3d get_f_label_pos(int idx) const;
 
      ///Get the symmetry element display
      /**\return The symmetry element display. */
      geom_disp *get_sym() const;

      ///Set the symmetry element display
      /**\param sym the symmetry element display to set. */
      void set_sym(const geom_disp &sym);

      ///Update animated displays.
      void animate();
 
};


class camera : public scene_item
{
   private:
      vec3d centre;
      vec3d lookat;
      
      mat3d rotation;
      mat3d label_rot;
      double distance;
      double width;
      double persp;
      double cut_dist;

      vec3d spin_axis;
      double spin_inc;
      double spin_rot;

      //Set up the label rotation matrix
      void set_label_rot();

   public:

      ///Constructor
      camera();

      ///Set the rotation centre
      /**\param cent the rotation centre. */
      void set_centre(vec3d cent = vec3d()) { centre = cent; }

      ///Get the rotation centre
      /**\return The rotation centre. */
      vec3d get_centre() const;
      
      ///Set the look-at point
      /**\param look the point the camera looks at. */
      void set_lookat(vec3d look=vec3d()) { lookat = look; }
      
      ///Get the look-at point
      /**\return The point the camera looks at. */
      vec3d get_lookat() const
         { return (lookat.is_set()) ? lookat : get_centre(); }
      
      ///Increment the look-at point
      /**\param inc a vector to displace the look at point. */
      void inc_lookat(vec3d inc) { set_lookat(get_lookat() + inc); }
      
      ///Set the distance from the look-at point
      /**\param dist the distance of the camera from the look-at point. */
      void set_distance(double dist) { distance = dist; }
      
      ///Get the distance from the look-at point
      /**\return The distance of the camera from the look-at point. */
      double get_distance() const
         { return (distance>0)? distance: 1.4*get_width(); }
      
      ///Set the width to view
      /**\param wdth width to view. */
      void set_width(double wdth) { width=wdth; }
      
      ///Get the width to view
      /**\return Width to view. */
      double get_width() const;

      ///Set the perspective
      /**Default value is 2, lower values widen the perspective
       * while larger values narrow it.
       * \param per value for the perspective. */
      void set_persp(double per=2) { persp=per; }

      ///Get the perspective
      /**\return The value of the perspective. */
      double get_persp() const { return persp; }
     
      ///Set the rotation
      /**\param rot the rotation. */
      void set_rotation(mat3d rot = mat3d()) { rotation = rot; set_label_rot();}
     
      ///Get the rotation
      /**\return The rotation. */
      const mat3d &get_rotation() const { return rotation; }
     
      ///Increment the rotation
      /**\param rot_inc the rotation increment. */
      void inc_rotation(mat3d rot_inc) {set_rotation(rot_inc * get_rotation());}
      
      ///Set the distance of the display cut plane
      /**This plane, perpendicular to the viewing direction, may cut
       * through objects and expose a section. Some displays may
       * have a minimum distance, which will be used if the
       * value passed is too low.
       * \param dist the cut plane distance. */
      void set_cut_dist(double dist=-1);
      
      ///Get the distance of the display cut plane
      /**This plane, perpendicular to the viewing direction, may cut
       * through objects and expose a section.
       * \return The cut plane distance. */
      double get_cut_dist() const;

      ///Increment the distance of the display cut plane
      /**This plane, perpendicular to the viewing direction, may cut
       * through objects and expose a section. Some displays may
       * have a minimum distance, which will be used if the incremented
       * value is too low.
       * \param inc the distance to increment the cut plane. */
      void inc_cut_dist(double inc) { set_cut_dist(get_cut_dist() + inc); }
      
      ///Set the spin rotation
      /**\param axis the spin axis.
       * \param rot the angle to rotate on each update. */
      void set_spin_rot(vec3d axis=vec3d(0,1,0), double rot=0)
         { stop_spinning(); spin_axis=axis; spin_inc=rot; }
      
      ///Get the spin rotation matrix
      /**\return The rotation on each update. */
      mat3d get_spin_rot() const
         {  return is_spinning()?mat3d::rot(spin_axis, spin_rot):mat3d(); }

      ///Increment the spin rotation angle
      void inc_spin_rot()
         { spin_rot += spin_inc; }

      ///Stop the spinning
      void stop_spinning()
         { if(is_spinning()) set_rotation(get_rotation()*get_spin_rot());
           spin_rot=0, spin_inc=0; }

      ///Check if spinning
      bool is_spinning() const { return fabs(spin_inc)>epsilon; }
              
      ///Get the size of text for labels
      double get_text_sz(vec3d pos) const;

      ///Update animated properties
      void animate();
};


class scene : public scene_item {
   private:
      timer anim_timer;
      double cycle_rate;

      col_val bg_col;
      double inf_dist;
      bound_sphere bound_sph;
      vector<scene_geom> geoms;
 
      int cur_cam_num;
      vector<camera> cams;

   public:
      ///Constructor
      scene();

      ///Copy constructor
      /**\param scen the scene to copy from. */
      scene(const scene &scen);

      ///Assignment operator
      /**\param scen the scene to assign from. */
      scene &operator=(const scene &scen);

      ///Set the background colour
      /**\param bg the background colour. */
      void set_bg_col(col_val bg=col_val(vec3d(0.9, 0.9, 0.9))) { bg_col=bg; }

      ///Get the background colour
      /**\return The background colour. */
      col_val get_bg_col() const { return bg_col; }
      
      ///Set the infinity distance
      /**Points further away than this are considered to be at infinity.
       * A non-positive distance indicates that all points are finite.
       * \param dist the infinity distance. */
      void set_inf_dist(double dist=-1) { inf_dist = dist; }
      
      ///Get the infinity distance
      /**Points further away than this are considered to be at infinity.
       * A non-positive distance indicates that all points are finite.
       * \return The infinity distance. */
      double get_inf_dist() const { return inf_dist; }

      ///Get a bounding sphere for the scene
      /**\return A bounding sphere. */
      const bound_sphere &get_bound_sph() const { return bound_sph; }
      
      ///Add a geometry to the scene
      /**\param sc_geom the geometry to add. */
      void add_geom(const scene_geom sc_geom);
      
      ///Add a geometry to the scene
      /**\param geom the geometry to add. */
      void add_geom(const geom_if &geom);

      ///Delete a geometry from the scene
      /**\param idx the index of the geometry to delete. */
      bool delete_geom(int idx);

      ///Get the scene geometries
      /**\return The scene geometries. */
      const vector<scene_geom> &get_geoms() const { return geoms; }

      ///Get the scene geometries
      /**\return The scene geometries. */
      vector<scene_geom> &get_geoms() { return geoms; }
   
      ///Add a camera to the scene
      /**\param cam the camera to add. */
      void add_camera(const camera &cam=camera());
   
      ///Delete a camera from the scene
      /**\param idx the index of the camera to delete. */
      bool delete_camera(int idx);

      ///Get the cameras
      /**\return The cameras. */
      const vector<camera> &get_cameras() const { return cams; }

      ///Get the cameras
      /**\return The cameras. */
      vector<camera> &get_cameras() { return cams; }

      ///Get the current camera
      /**\return The current camera. */
      const camera &cur_camera() const;

      ///Get the current camera
      /**\return The current camera. */
      camera &cur_camera();

      ///Get the camera name
      /**\param idx the index of the camera.
       * \return The camera name. */
      string get_camera_name(int idx) const;

      ///Get the width of the scene
      /**The width is the diameter of a bounding sphere.
       * \return The width of the scene. */
      double get_width() const { return bound_sph.get_width(); }

      ///Get the centre of the scene
      /**The centre is the centre of a bounding sphere.
       * \return The centre of the scene. */
      vec3d get_centre() const { return bound_sph.get_centre(); }
      
      ///Set cycle rate
      /**The rate at which animation updates occur.
       * \param cps The rate, in cycles per second. */
      void set_cycle_rate(double cps) {cycle_rate=(cps<0)?0:cps;}
      
      ///Update animated items.
      void animate();
};

#endif //SCENE



