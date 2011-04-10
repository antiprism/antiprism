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

/*!\file bbox.h
   \brief Bounding box and sphere
*/

#ifndef BBOX_H
#define BBOX_H

#include "geom.h"

/// Bounding box
/** A bounding box aligned with the coordinate planes. */ 
class bound_box
{
   private:
      vec3d min_coords; /**< minimum coordinates, at left bottom back corner*/
      vec3d max_coords; /**< maximum coordinates, at right top front corner*/

   public:
      /// Constructer
      /**\param points points to find the bounding box for
       * \param cutoff ignore points beyond this distance from the
       * origin. A negative value indicates there is no cut off
       * distance.  */
      bound_box(const vector<vec3d> &points=vector<vec3d>(), double cutoff=-1);
      
      /// Destructor
      virtual ~bound_box() {}

      /// Add points and calculate the new bounding box
      /**\param points points to add
       * \param cutoff ignore points beyond this distance from the
       * origin. A negative value indicates there is no cut off distance.*/
      void add_points(const vector<vec3d> &points, double cutoff=-1);
      
      /// Add a bounding box and calculate the new bounding box
      /**\param b_box bounding box to add */
      void add_b_box(const bound_box &b_box);
      
      /// Maximum width of the box
      /**\return The length of a diagonal */
      double max_width() const { return (max_coords-min_coords).mag(); }
      
      /// Centre of the box
      /**\return The centre coordinates */
      virtual vec3d get_centre() const { return (max_coords+min_coords)*0.5; }
      
      /// Minimum coordinates
      /**\return The minimum coordinates, at left bottom back corner */
      const vec3d &get_min() const { return min_coords; }
      
      /// Maximum coordinates
      /**\return The maximum coordinates, at right top front corner */
      const vec3d &get_max() const { return max_coords; }
};


/// Bounding Sphere
/** An approximate bounding sphere.
 */
class bound_sphere
{
   private:
      double radius;
      vec3d centre;
      double cut_off;   // ignore points beyond this distance, if positive
      int cut_off_cnt;  // number of points that were cut off

      void find_radius_centre(const vector<vec3d> &pts);
      void add_b_sphere(vec3d cent, double rad);

   public:
      /// Constructer
      /**\param points points to find the bounding sphere for
       * \param cutoff ignore points beyond this distance from the
       * origin. A negative value indicates there is no cut off
       * distance.  */
      bound_sphere(const vector<vec3d> &points=vector<vec3d>(),
            double cutoff=-1);
      
      ///Destructor
      virtual ~bound_sphere() {}
      
      /// Add points and calculate the new bounding sphere
      /**\param points points to add
       * \param cutoff ignore points beyond this distance from the
       * origin. A negative value indicates there is no cut off distance.*/
      void add_points(const vector<vec3d> &points, double cutoff=-1);
      
      
      /// Add a bounding sphere and calculate the new bounding sphere
      /**\param b_sphere bounding sphere to add */
      void add_b_sphere(const bound_sphere &b_sphere);
      
      /// Set the cut off distance
      /**\param cutoff ignore points beyond this distance from the
       * origin. A negative value indicates there is no cut off distance.*/
      void set_cut_off(double cutoff);

      /// Centre of the sphere
      /**\return The centre coordinates */
      virtual vec3d get_centre() const { return centre; }

      /// Radius of the sphere
      /**\return The radius */
      double get_radius() const { return radius; }
      
      /// Width of the points (diameter of the sphere)
      /**\return The width */
      double get_width() const { return 2*radius; }
      
      /// The number of points excluded by the cut-off distance
      /**\return The number of points cut off */
      double get_cut_off_cnt() const { return cut_off_cnt; }
};


#endif // BBOX_H

