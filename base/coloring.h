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

/*!\file coloring.h
 * \brief Classes to color all element of a type.
 */

#ifndef COLORING_H
#define COLORING_H

#include "colormap.h"
#include "geometry.h"

#include <set>

namespace anti {

/// Class for colouring the elements of a geometry
class Coloring : public ColorMapMulti {
private:
  // The geometry to colour.
  Geometry *geom;

  unsigned int cycle_msecs;

  void face_edge_color(const std::vector<std::vector<int>> &elems,
                       const ElemProps<Color> &cols);

  void edge_color_and_branch(int idx, int part, bool apply_map,
                             std::vector<std::vector<int>> &vcons,
                             std::vector<bool> &seen);

protected:
  /// Colour index for a unit vector by its z-value
  /** Index numbers are assigned in the range of 0 up to the last
   *  entry in the first colour map.
   * \param vec the unit vector
   * \param cent the centre of the geometry.
   * \param height the height of the geometry along the z-axis.
   * \param def_sz allocate index numbers up to def_sz if no map
   *  has been set, size of the first map if one is set.
   * \return The index. */
  int z_gradient(Vec3d vec, Vec3d cent = Vec3d(0, 0, 0), double height = 2,
                 int def_sz = 256);

  /// Set up the lights
  /**\param lts the lights to prepare for use, or if there are no
   *  lights then a default set is provided. */
  void setup_lights(Geometry &lts);

  /// Colour for a position vector from a set of lights.
  /**\param vec the position vector.
   * \param lts the lights (prepared for use by setup_lights().)
   * \return The final RGBA colour. */
  Color light(Vec3d vec, Geometry &lts);

  /// Convert all colour index numbers into colour values.
  /**\param cols the colours of the elements, by element index. */
  void set_all_idx_to_val(std::map<int, Color> &cols);

  /// Get the geometry that is being coloured.
  /**\return A pointer to the geometry. */
  Geometry *get_geom() { return geom; }

public:
  /// Constructor
  /**\param geo the geometry to be coloured. */
  Coloring(Geometry *geo = nullptr);

  /// Copy Constructor
  /**\param clrng the Coloring to copy from. */
  Coloring(const Coloring &clrng);

  /// Copy Assignment
  /**\param clrng the Coloring to copy from.
   * \return A reference to this object.*/
  Coloring &operator=(const Coloring &clrng);

  /// Destructor
  ~Coloring();

  /// Set the geometry to colour.
  /**\param geo the geometry to colour. */
  void set_geom(Geometry *geo) { geom = geo; }

  /// Cycle the map colours
  /** Each colour index number is mapped to the previous colour
   *  value in the colour map. */
  void cycle_map_cols();

  /// Set the time between colour map cycles
  /** An interval of 0 indicates no cycling.
   * \param interval the number of milliseconds between colour map cycles.*/
  void set_cycle_msecs(unsigned int interval) { cycle_msecs = interval; }

  /// Get the time between colour map cycles
  /** An interval of 0 indicates no cycling.
   * \return The number of milliseconds between colour map cycles. */
  unsigned int get_cycle_msecs() const { return cycle_msecs; }

  /// Colour vertices with a single colour
  /**\param col colour for the vertices. */
  void v_one_col(Color col);

  /// Colour each vertex set with a single colour
  /**\param equivs the indexes in each set are given the same colour.
   * \param apply_map if \c false, colour with index numbers,
   *  if \c true, convert
   *  these index numbers using the colour maps */
  void v_sets(const std::vector<std::set<int>> &equivs, bool apply_map = true);

  /// Colour each vertex with a different colour.
  /**\param apply_map if \c false, colour with index numbers,
   *  if \c true, convert these index numbers using the colour maps */
  void v_unique(bool apply_map = true);

  /// Colour vertices with a proper colouring.
  /** Colour so that no two adjacent vertices on a face have the
   *  same colour.
   * \param apply_map if \c false, colour with index numbers,
   *  if \c true, convert these index numbers using the colour maps */
  void v_proper(bool apply_map = true);

  /// Colour vertices by face colour.
  /** Colour the vertices with the colour of a face they are on. */
  void v_face_color();

  /// Colour vertices by edge colour.
  /** Colour the vertices with the colour of an edge they are on. */
  void v_edge_color();

  /// Colour each vertex by the number of faces it lies on.
  /**\param apply_map if \c false, colour with index numbers,
   *  if \c true, convert these index numbers using the colour maps */
  void v_order(bool apply_map = true);

  /// Colour each vertex by its y-coordinate.
  /**\param apply_map if \c false, colour with index numbers,
   *  if \c true, convert these index numbers using the colour maps */
  void v_position(bool apply_map = true);

  /// Colour each vertex by the average internal angle of a "vertex figure" (to
  /// nearest degree.)
  /**\param apply_map if \c false, colour with index numbers,
   *  if \c true, convert these index numbers using the colour maps */
  void v_avg_angle(bool apply_map = true);

  /// Colour each vertex using a set of lights.
  /**\param lts a geometry holding coloured vertices to use as lights.*/
  void v_lights(Geometry lts);

  /// Apply the colour map to turn vertex index numbers into colour values.
  /** Uses the value set by \c handle_no_map() for converting index
   *  numbers not in the map . The HSVA ranges are respected
   *  if a random colour is allocated. */
  void v_apply_cmap();

  /// Colour faces with a single colour.
  /**\param col colour for the faces. */
  void f_one_col(Color col);

  /// Colour each face set with a single colour
  /**\param equivs the indexes in each set are given the same colour.
   * \param apply_map if \c false, colour with index numbers,
   *  if \c true, convert these index numbers using the colour maps */
  void f_sets(const std::vector<std::set<int>> &equivs, bool apply_map = true);

  /// Colour each face with a different colour.
  /**\param apply_map if \c false, colour with index numbers,
   *  if \c true, convert these index numbers using the colour maps */
  void f_unique(bool apply_map = true);

  /// Colour faces with a proper colouring.
  /** Colour so that no two adjoining faces have the same colour.
   * \param apply_map if \c false, colour with index numbers,
   *  if \c true, convert these index numbers using the colour maps */
  void f_proper(bool apply_map = true);

  /// Colour each face by the number of sides it has.
  /**\param apply_map if \c false, colour with index numbers,
   * if \c true, convert these index numbers using the colour maps */
  void f_sides(bool apply_map = true);

  /// Colour each face by the average internal angle (to nearest degree.)
  /**\param apply_map if \c false, colour with index numbers,
   *  if \c true, convert these index numbers using the colour maps */
  void f_avg_angle(bool apply_map = true);

  /// Colour each face by the set of connected faces it is part of.
  /** To be connected two faces must share an edge.
   * \param apply_map if \c false, colour with index numbers,
   *  if \c true, convert these index numbers using the colour maps */
  void f_parts(bool apply_map = true);

  /// Colour each face by the y-component of the normal.
  /**\param apply_map if \c false, colour with index numbers,
   * if \c true, convert these index numbers using the colour maps */
  void f_normal(bool apply_map = true);

  /// Colour each face by the y-coordinate of the centroid.
  /**\param apply_map if \c false, colour with index numbers,
   *  if \c true, convert these index numbers using the colour maps */
  void f_centroid(bool apply_map = true);

  /// Colour each face by normal using a set of lights
  /**\param lts a geometry holding coloured vertices to use as lights. */
  void f_lights(Geometry lts);

  /// Colour each face by centroid using a set of lights
  /**\param lts a geometry holding coloured vertices to use as lights.*/
  void f_lights2(Geometry lts);

  /// Apply the colour map to turn face index numbers into colour values.
  /** Uses the value set by \c handle_no_map() for converting index
   *  numbers not in the map . The HSVA ranges are respected
   *  if a random colour is allocated. */
  void f_apply_cmap();

  /// Colour edges with a single colour
  /**\param col colour for the edges. */
  void e_one_col(Color col);

  /// Colour each edge by the number of edges connected to it.
  /**\param apply_map if \c false, colour with index numbers,
   *  if \c true, convert these index numbers using the colour maps */
  void e_order(bool apply_map = true);

  /// Colour each edge set with a single colour
  /**\param equivs the indexes in each set are given the same colour.
   * \param apply_map if \c false, colour with index numbers,
   *  if \c true, convert these index numbers using the colour maps */
  void e_sets(const std::vector<std::set<int>> &equivs, bool apply_map = true);

  /// Colour each edge with a different colour.
  /**\param apply_map if \c false, colour with index numbers,
   *  if \c true, convert these index numbers using the colour maps */
  void e_unique(bool apply_map = true);

  /// Proper edge colouring
  /** Colour so that no two adjacent edges on a face have the same colour.
   * \param apply_map if \c false, colour with index numbers,
   *  if \c true, convert these index numbers using the colour maps */
  void e_proper(bool apply_map = true);

  /// Colour by face colour
  /** Colour the edges with the colour of a face they are part of. */
  void e_face_color();

  /// Colour each edge by the set of connected edges it is part of.
  /** To be connected two edges must share a vertex.
   * \param apply_map if \c false, colour with index numbers,
   *  if \c true, convert these index numbers using the colour maps */
  void e_parts(bool apply_map = true);

  /// Colour each edge by the y-component of its direction.*/
  /**\param apply_map if \c false, colour with index numbers,
   *  if \c true, convert these index numbers using the colour maps */
  void e_direction(bool apply_map = true);

  /// Colour each edge by the y-coordinate of its mid-point.
  /**\param apply_map if \c false, colour with index numbers,
   *  if \c true, convert these index numbers using the colour maps */
  void e_mid_point(bool apply_map = true);

  /// Colour each edge by its unsigned direction.*/
  /**\param apply_map if false, colour with index numbers, if true, convert
   * these index numbers using the colour maps */
  void e_vector(bool apply_map = true);

  /// Colour each edge according to unique lengths
  /** Edges are considered in length order until an edge is found
   *  that is is greater than the minimum difference from the first
   *  edge in the sequence
   * \param min_diff the minimum difference
   * \param apply_map if \c false, colour with index numbers,
   *  if \c true, convert these index numbers using the colour maps */
  void e_lengths(double min_diff, bool apply_map = true);

  /// Colour each edge using a set of lights.
  /**\param lts a geometry holding coloured vertices to use as lights.*/
  void e_lights(Geometry lts);

  /// Colour each edge by edge direction using a set of lights.
  /**\param lts a geometry holding coloured vertices to use as lights.*/
  void e_dir_lights(Geometry lts);

  /// Apply the colour map to turn edge index numbers into colour values.
  /**Uses the value set by \c handle_no_map() for converting index
   * numbers not in the map . The HSVA ranges are respected
   * if a random colour is allocated. */
  void e_apply_cmap();

  /// Colour all elements with a colour for each type
  /**Original colours for a type are left if the %Color value is not set
   * \param vert_col vertex colour.
   * \param edge_col edge colour.
   * \param face_col face colour. */
  void vef_one_col(Color vert_col, Color edge_col, Color face_col);

  /// Apply transparancy to faces of the geom
  /** unset color and map indexes cannot be made transparent
   * \param face_opacity value from 1 to 255 */
  Status apply_transparency(const int);
};

/// Read colourings for vertex, edge and face elements from a string.
/**\param clrngs an array of three colorings for the vertex, edge and
 *  face elements.
 * \param line the string containing the map specifier, applying to all
 *  element types by default, but letters from vef may by specified after
 *  a comma to indicate the element types that the coloring will be set for.
 * \return status, which evaluates to \c true if a valid colourings line was
 *  read, otherwise \c false */
Status read_colorings(Coloring clrngs[], const char *line);

} // namespace anti

#endif // COLORING_H
