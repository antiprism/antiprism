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

/*!\file coloring.h
 * \brief Classes to color all element of a type.
 */

#ifndef COLORING_H
#define COLORING_H

#include <set>

#include "geom.h"
#include "col_map.h"

using std::set;

///Class for colouring the elements of a geometry 
class coloring {
   private:
      //The geometry to colour.
      col_geom_v *geom;
      
      //The colour maps to be tried sequentially.
      vector<color_map *> cmaps;

      unsigned int cycle_msecs;

      void face_edge_color(const vector<vector<int> > &elems,
            const map<int,col_val> &cmap);
      
      void edge_color_and_branch(int idx, int part, bool as_values, 
            vector<vector<int> > &vcons, vector<bool> &seen);

   protected:

      ///Colour index for a unit vector by its y value
      /**Index numbers are assigned in the range of 0 up to the last
       * entry in the first colour map.
       * \param vec the unit vector
       * \param cent the centre of the geometry.
       * \param height the height of the geometry.
       * \param def_sz allocate index numbers upto def_sz if no map
       * has been set.
       * size of the first map if one is set.
       * \return The index. */
      int y_gradient(vec3d vec, vec3d cent=vec3d(0,0,0), double height=2,
            int def_sz=256);
      

      ///Set up the lights
      /**\param lts the lights to prepare for use, or if there are no
       * lights then a default set is provided. */
      void setup_lights(col_geom_v &lts);

      ///Colour for a position vector from a set of lights.
      /**\param vec the position vector.
       * \param lts the lights (prepared for use by setup_lights().)
       * \return The final RGBA colour. */
      col_val light(vec3d vec, col_geom_v &lts);
      
      ///Convert all colour index numbers into colour values.
      /**\param cols the colours of the elements, by element index. */
      void set_all_idx_to_val(map<int, col_val> &cols);

      ///Get the geometry that is being coloured.
      /**\return A pointer to the geometry. */
      col_geom_v *get_geom() { return geom; }
      
   public:
      ///Constructor
      /**\param geo the geometry to be coloured. */
      coloring(col_geom_v *geo=0);

      ///Copy Constructor
      /**\param clrng the coloring to copy from. */
      coloring(const coloring &clrng);
      
      ///Copy Assignment
      /**\param clrng the coloring to copy from. */
      coloring &operator =(const coloring &clrng);
      
      ///Destructor
      ~coloring();

      ///Add a colour map.
      /**\param col_map the colour map.
       * \param pos the position to add it, or at the end if pos is
       * greater then the current size */
      void add_cmap(color_map *col_map, unsigned int pos=UINT_MAX);

      ///Delete a colour map.
      /**\param pos the position of the colour map to delete, or delete
       * the last colour map if \c pos is greater than or equal to the
       * current size */
      void del_cmap(unsigned int pos=UINT_MAX);

      ///Get a the colour maps.
      /**\return The colour maps. */
      const vector<color_map *> &get_cmaps() {return cmaps;}

      ///Set the geometry to colour.
      /**\param geo the geometry to colour. */
      void set_geom(col_geom_v *geo) { geom = geo; }

      ///Convert a colour index into a colour value.
      /** \param idx the index to convert.
       * \return The final RGBA colour. */
      col_val idx_to_val(int idx) const;

      ///Cycle the map colours
      /**Each colour index number is mapped to the previous colour
       * value in the colour map. */
      void cycle_map_cols();

      ///Set the time between colour map cycles
      /**An interval of 0 indicates no cycling.
       * \param interval the number of milliseconds between colour map cycles.*/
      void set_cycle_msecs(unsigned int interval) { cycle_msecs = interval; }

      ///Get the time between colour map cycles
      /**An interval of 0 indicates no cycling.
       * \return The number of milliseconds between colour map cycles. */
      unsigned int get_cycle_msecs() const { return cycle_msecs; }

      
      ///Colour vertices with a single colour
      /**\param col colour for the vertices. */
      void v_one_col(col_val col);

      ///Colour each vertex set with a single colour
      /**\param equivs the indexes in each set are given the same colour.
       * \param as_values if true then colour with values, if false then
       * colour with index numbers. */
      void v_sets(const vector<set<int> > &equivs, bool as_values=true);

      ///Colour each vertex with a different colour.
      /**\param as_values if true then colour with values, if false then
       * colour with index numbers. */
      void v_unique(bool as_values=true);

      ///Colour vertices with a proper colouring.
      /**Colour so that no two adjacent vertices on a face have the same colour.
       * \param as_values if true then colour with values, if false then
       * colour with index numbers. */
      void v_proper(bool as_values=true);

      ///Colour vertices by face colour.
      /**Colour the vertices with the colour of a face they are on. */
      void v_face_color();

      ///Colour vertices by edge colour.
      /**Colour the vertices with the colour of an edge they are on. */
      void v_edge_color();
      
      ///Colour each vertex by the number of faces it lies on.
      /**\param as_values if true then colour with values, if false then
       * colour with index numbers. */
      void v_order(bool as_values=true);

      ///Colour each vertex by its y-coordinate.
      /**\param as_values if true then colour with values, if false then
       * colour with index numbers. */
      void v_position(bool as_values=true);
      
      ///Colour each vertex using a set of lights.
      /**\param lts a geometry holding coloured vertices to use as lights.*/
      void v_lights(col_geom_v lts);
      
      ///Apply the colour map to turn vertex index numbers into colour values.
      /**Uses the value set by \c handle_no_map() for converting index
       * numbers not in the map . The HSVA ranges are respected
       * if a random colour is allocated. */
      void v_apply_cmap();
      
      
      ///Colour faces with a single colour.
      /**\param col colour for the faces. */
      void f_one_col(col_val col);

      ///Colour each face set with a single colour
      /**\param equivs the indexes in each set are given the same colour.
       * \param as_values if true then colour with values, if false then
       * colour with index numbers. */
      void f_sets(const vector<set<int> > &equivs, bool as_values=true);

      ///Colour each face with a different colour.
      /**\param as_values if true then colour with values, if false then
       * colour with index numbers. */
      void f_unique(bool as_values=true);

      ///Colour faces with a proper colouring.
      /**Colour so that no two adjoining faces have the same colour.
       * \param as_values if true then colour with values, if false then
       * colour with index numbers. */
      void f_proper(bool as_values=true);
      
      ///Colour each face by the number of sides it has.
      /**\param as_values if true then colour with values, if false then
       * colour with index numbers. */
      void f_sides(bool as_values=true);

      ///Colour each face by the average internal angle (to nearest degree.)
      /**\param as_values if true then colour with values, if false then
       * colour with index numbers. */
      void f_avg_angle(bool as_values=true);

      ///Colour each face by the set of connected faces it is part of.
      /**To be connected two faces must share an edge.
       * \param as_values if true then colour with values, if false then
       * colour with index numbers. */
      void f_parts(bool as_values=true);
      
      ///Colour each face by the y-component of the normal.
      /**\param as_values if true then colour with values, if false then
       * colour with index numbers. */
      void f_normal(bool as_values=true);
      
      ///Colour each face by the y-coordinate of the centroid.
      /**\param as_values if true then colour with values, if false then
       * colour with index numbers. */
      void f_centroid(bool as_values=true);
      
      ///Colour each face by normal using a set of lights 
      /**\param lts a geometry holding coloured vertices to use as lights. */
      void f_lights(col_geom_v lts);

      ///Colour each face by centroid using a set of lights 
      /**\param lts a geometry holding coloured vertices to use as lights.*/
      void f_lights2(col_geom_v lts);

      ///Apply the colour map to turn face index numbers into colour values.
      /**Uses the value set by \c handle_no_map() for converting index
       * numbers not in the map . The HSVA ranges are respected
       * if a random colour is allocated. */
      void f_apply_cmap();
      
      
      ///Colour edges with a single colour
      /**\param col colour for the edges. */
      void e_one_col(col_val col);

      ///Colour each edge set with a single colour
      /**\param equivs the indexes in each set are given the same colour.
       * \param as_values if true then colour with values, if false then
       * colour with index numbers. */
      void e_sets(const vector<set<int> > &equivs, bool as_values=true);

      ///Colour each edge with a different colour.
      /**\param as_values if true then colour with values, if false then
       * colour with index numbers. */
      void e_unique(bool as_values=true);

      ///Proper edge colouring
      /**Colour so that no two adjacent edges on a face have the same colour.
       * \param as_values if true then colour with values, if false then
       * colour with index numbers. */
      void e_proper(bool as_values=true);

      ///Colour by face colour
      /**Colour the edges with the colour of a face they are part of. */
      void e_face_color();

      ///Colour each edge by the set of connected edges it is part of.
      /**To be connected two edges must share a vertex.
       * \param as_values if true then colour with values, if false then
       * colour with index numbers. */
      void e_parts(bool as_values=true);

      ///Colour each edge by the y-component of its direction.*/
      /**\param as_values if true then colour with values, if false then
       * colour with index numbers. */
      void e_direction(bool as_values=true);
      
      ///Colour each edge by the y-coordinate of its mid-point.
      /**\param as_values if true then colour with values, if false then
       * colour with index numbers. */
      void e_mid_point(bool as_values=true);
      
      ///Colour each edge using a set of lights.
      /**\param lts a geometry holding coloured vertices to use as lights.*/
      void e_lights(col_geom_v lts);
      
      ///Apply the colour map to turn edge index numbers into colour values.
      /**Uses the value set by \c handle_no_map() for converting index
       * numbers not in the map . The HSVA ranges are respected
       * if a random colour is allocated. */
      void e_apply_cmap();
};


bool read_colorings(coloring clrng[], const char *line, char *errmsg=0,
      int max_parts=4);

#endif // COLORING_H


