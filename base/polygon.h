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

/*!\file polygon.h
   \brief Generate polyhedra based on polygons.
*/

#ifndef POLYGON_H
#define POLYGON_H

#include <cmath>
#include <string.h>

#include "geometry.h"
#include "geometryutils.h"
#include "mathutils.h"

namespace anti {

/// Make a uniform polygon based model
/**Model types, subtypes and parameters.
 * 1.  prism (A0: polygons twist, forming an antiprism)
 *        subtypes: 1. antiprism, from triangulating side faces
 *                  2. trapezohedron, with this antiprism-ised belt
 *                  3. crown (A1 integer to specify a second polygon step)
 * 2.  antiprism (A0: polygons twist)
 *        subtypes: 1. trapezohedron, with this antiprism belt
 *                  2. antihermaphrodite, with this antiprism base
 *                  3. scalenohedron (L1 for apex height)
 *                  4. subdivided_scalenohedron (L1 for apex height)
 *                  5. crown (A1 integer to specify a second polygon step)
 * 3.  pyramid (A0: base separates, polygons twist to antihermaphrodite)
 *        subtypes: 1. antihermaphrodite, with this antiprism base
 *                  2. elongated (L1 for prism height)
 *                  3. gyroelongated (L1 for antiprism height)
 * 4.  dipyramid (A0: base separates, polygons twist to trapezohedron)
 *        subtypes: 1. trapezohedron, with this pyramid apex
 *                  2. elongated (L1 for prism height)
 *                  3. gyroelongated (L1 for antiprism height)
 *                  4. dipyramid_scalenohedron (R1 for alternate vertex radius)
 * 5.  cupola
 *        subtypes: 1. elongated (L1 for prism height)
 *                  2. gyroelongated (L1 for antiprism height)
 *                  3. cupoloid
 * 6.  orthobicupola
 *        subtypes: 1. elongated (L1 for prism height)
 *                  2. gyroelongated (L1 for antiprism height)
 * 7.  gyrobicupola
 *        subtypes: 1. elongated (L1 for prism height)
 *                  2. gyroelongated (L1 for antiprism height)
 * 8.  snub-antiprism
 *        subtypes: 1. inverted, triangle band inverted
 * 9.  dihedron
 *        subtypes: 1. polygon
 * 10. crown polyhedron (A1 integer to specify a second polygon step
 *                                    or A0 angle to specify a twist)
 */

class Polygon {
private:
  Status make_dihedron_part(Geometry &geom);
  Status make_prism_part(Geometry &geom);
  Status make_antiprism_part(Geometry &geom);
  Status make_pyramid_part(Geometry &geom);
  Status make_dipyramid_part(Geometry &geom);
  Status make_cupola_part(Geometry &geom);
  Status make_orthobicupola_part(Geometry &geom);
  Status make_gyrobicupola_part(Geometry &geom);
  Status make_snub_antiprism_part(Geometry &geom);
  Status make_dipyramid_scal_part(Geometry &geom);
  Status make_antiprism_scal_part(Geometry &geom);
  Status make_antiprism_subscal_part(Geometry &geom);
  Status make_prism_crown_part(Geometry &geom);
  Status make_antiprism_crown_part(Geometry &geom);
  Status make_crown_full(Geometry &geom);
  Status make_crown_part(Geometry &geom);
  double get_antiprism_height();

protected:
  int num_sides; ///< The number of sides of the polygon (n of {n/d})
  int step;      ///< The number of verts stepped by a side (d of {n/d})
  int parts;     ///< The number of parts (polygon may be compound).
  std::vector<double> radius;      ///< radius values
  std::vector<double> height;      ///< height values
  std::vector<double> edge;        ///< edge values
  std::vector<double> twist_angle; ///< twist angle values
  int type;                        ///< The type of a %Polygon based polyhedron.
  int subtype; ///< The subtype of a %Polygon based polyhedron.

  /// Set a parameter
  /**\param param the parameter type.
   * \param idx the parameter index.
   * \param val the value to set the parameter
   * \return true if parameter could be set, otherwise false. */
  bool set_param(std::vector<double> &param, int idx, double val);

  /// Get a parameter
  /**\param param the parameter type.
   * \param idx the parameter index.
   * \return The circumradius (or NAN if radius not set or index is out
   *  of range). */
  double get_param(std::vector<double> &param, int idx);

public:
  enum {
    unknown = 0,
    prism,
    antiprism,
    pyramid,
    dipyramid,
    cupola,
    orthobicupola,
    gyrobicupola,
    snub_antiprism,
    dihedron,
    crown,
    types_end
  };
  enum { sub_default = 0 };
  enum { sub_dihedron_polygon = 1 };
  enum { sub_prism_antiprism = 1, sub_prism_trapezohedron, sub_prism_crown };
  enum {
    sub_antiprism_trapezohedron = 1,
    sub_antiprism_antihermaphrodite,
    sub_antiprism_scalenohedron,
    sub_antiprism_subdivided_scalenohedron,
    sub_antiprism_crown
  };
  enum {
    sub_pyramid_antihermaphrodite = 1,
    sub_pyramid_elongated,
    sub_pyramid_gyroelongated
  };
  enum {
    sub_dipyramid_trapezohedron = 1,
    sub_dipyramid_elongated,
    sub_dipyramid_gyroelongated,
    sub_dipyramid_scalenohedron
  };
  enum {
    sub_cupola_elongated = 1,
    sub_cupola_gyroelongated,
    sub_cupola_cuploid
  };

  enum { R0 = 1, R1 = 2, H0 = 4, H1 = 8, E0 = 16, E1 = 32, A0 = 64, A1 = 128 };

  /// Constructor
  /** Polygon in form {N/D} (with N/D not necessarily in lowest form.)
   * \param N number of sides to the (compound) polygon.
   * \param D the number of vertices stepped by an edge (default 1)
   * \param type polyhedron type
   * \param subtype polyhedron subtype */
  Polygon(int N = 3, int D = 1, int type = 0, int subtype = 0);

  /// Set polygon fraction
  /** Polygon in form {N/D} (with N/D not necessarily in lowest form.)
   * \param N number of sides to the (compound) polygon.
   * \param D the number of vertices stepped by an edge (default 1)
   * \return status, evaluates to \c true if the fraction was valid,
   *  otherwise \c false.*/
  Status set_fraction(int N, int D = 1);

  /// Return a polygon intialised as a particular type
  /** Polygon in form {N/D} (with N/D not necessarily in lowest form.)
   * \param type the type of polygon-based model.
   * \param N number of sides to the (compound) polygon.
   * \param D the number of vertices stepped by an edge (default 1)
   * \return The polygon. */
  static Polygon as_type(int type, int N, int D = 1);

  /// Destructor
  virtual ~Polygon() = default;

  /// Set a radius
  /**\param idx the radius index number
   * \param val value for the radius with this index number
   * \return true if radius could be set, otherwise false. */
  bool set_radius(int idx, double val);

  /// Get a radius
  /**\param idx the radius index number
   * \return The radius (or NAN if radius not set or index is out
   *  of range). */
  double get_radius(int idx);

  /// Set a height
  /**\param idx the height index number
   * \param val value for the height with this index number
   * \return true if height could be set, otherwise false. */
  bool set_height(int idx, double val);

  /// Get a height
  /**\param idx the height index number
   * \return The height (or NAN if height not set or index is out
   *  of range). */
  double get_height(int idx);

  /// Set an edge length
  /**\param idx the edge length index number
   * \param val value for the edge length with this index number
   * \return true if edge length could be set, otherwise false. */
  bool set_edge(int idx, double val);

  /// Get an edge length
  /**\param idx the edge length index number
   * \return The edge length (or NAN if edge length not set or index is out
   *  of range). */
  double get_edge(int idx);

  /// Set a twist angle
  /**\param idx the twist angle index number
   * \param val value for the twist angle (degrees) with this index number
   * \return true if the angle could be set, otherwise false. */
  bool set_twist_angle(int idx, double val);

  /// Get a twist angle
  /**\param idx the twist angle index number
   * \return The twist angle (or NAN if angle not set or index is out
   *  of range). */
  double get_twist_angle(int idx);

  /// Get the angle that an edge makes at the centre
  /**\return the angle, in radians. */
  double angle() { return 2 * M_PI * step / num_sides; }

  /// Get the number of sides of the (component) polygon.
  /**\return the number of sides. */
  int get_num_sides() { return num_sides; }

  /// Get the number of vertices stepped by a side of the (component) polygon
  /**\return the number of vertices stepped. */
  int get_step() { return step; }

  /// Get the number of component parts of a polygon
  /** If the polygon step is not in lowest form the polygon will
   *  be compound and have more than one part.
   * \return the number of sides. */
  int get_parts() { return parts; }

  /// Get the polygon radius
  /**\return The radius */
  double get_polygon_radius();

  /// Get the polygon radius
  /**\return The radius */
  double get_polygon_edge();

  /// Add a polygon aligned with the xy plane and at a given z-height
  /**\param geom the geometry to add the polygon to.
   * \param ht the height on the z-axis to place the polygon. */
  void add_polygon(Geometry &geom, double ht = 0);

  /// Check if parameter value was set
  /**\param val the parameter value to check
   * \return \c true if set, otherwise \c false. */
  bool value_is_set(double val);

  /// Set the type of the %Polygon based polyhedron.
  /**\param typ the type number, may be \c Polygon::prism,
   * \c Polygon::antiprism, \c Polygon::pyramid, \c Polygon::dipyramid,
   * \c Polygon::cupola \c Polygon::orthobicupola, \c Polygon::gyrobicupola,
   * \c Polygon::snub-antiprism, \c Polygon::dihedron.
   * \return status, evaluates to \c true if the type was valid,
   *  otherwise \c false. */
  virtual Status set_type(int typ);

  /// Get the type of the %Polygon based polyhedron.
  /**\return the type, an integer that may be tested against
   *  \c Polygon::antiprism, \c Polygon::pyramid, \c Polygon::dipyramid,
   *  \c Polygon::cupola \c Polygon::orthobicupola,
   *  \c Polygon::gyrobicupola, \c Polygon::snub-antiprism,
   *  \c Polygon::dihedron.*/
  virtual int get_type();

  /// Set the subtype of the %Polygon based polyhedron.
  /** Some %Polygon based polyhedra come in several forms, and setting
   *  the sub-type can select a particular form.
   * \param subtyp the sub-type number.
   * \return status, evaluates to \c true if the sub-type was valid,
   *  otherwise \c false. */
  virtual Status set_subtype(int subtyp);

  /// Get the subtype of the %Polygon based polyhedron.
  /**\return the subtype, an integer that may be tested against
   *  the valid subtypes for type.*/
  virtual int get_subtype();

  /// Get the parameters compatible with the model type and subtype
  /**\param params flags for the accepatable parameters. Test against
   *  R0, R1, H0, etc.
   * \param conflicts flags indicting conflicting parameters. Test against
   *  R0, R1, etc.
   * \return the maximum subtype for the type (or -1 if type is invalid).*/
  int get_acceptable_params(unsigned int &params,
                            std::vector<unsigned int> &conflicts);

  /// Check whether a model can use particular parameters
  /**\param params flags for R0, R1, etc to check against acceptable
   *  parameters for the model type and subtype.
   * \return \c true if the parameters can be set, otherwise \c false */
  bool has_params(unsigned int params);

  /// Get flags for the parameters which have been set
  /** This can be used for checking against conflicts flags
   * \return The flags for set parameters */
  unsigned int get_params_set();

  /// Repeat a component of a Polygon-based polyhedron to make the compound
  /** Polygon::make_poly will repeat the part \c parts times.
   * \param geom a geometry to return the polyhedron.
   * \param part the component to repeat */
  void repeat_part(Geometry &geom, const Geometry &part);

  /// Make a polygon based polyhedron.
  /**\param geom a geometry to return the polyhedron.
   * \return status, evaluates to \c true if the construction succeeded,
   *  otherwise \c false.*/
  virtual Status make_poly(Geometry &geom);

  /// Dump polygon data to stderr
  void dump();
};

/// Make a uniform model of a polygon-based polyhedron.
/** The model has all its edges set to one and its faces are regular. For
 *  prism's, antiprisms and pyramids the resulting model wil be uniform.
 * \param geom a geometry to return the model.
 * \param pgon a Polygon-derived object of the polyhedron type required;
 * \return true if the edge types could be set to length \c 1.0,
 *  otherwise \c false. */
template <class T> bool uni_pgon(Geometry &geom, T pgon)
{
  pgon.set_edge(0, 1.0);
  bool ret = pgon.set_edge(1, 1.0);
  pgon.make_poly(geom);
  return ret;
}

// Inline functions
//
inline int Polygon::get_type() { return type; }

inline int Polygon::get_subtype() { return type; }

inline Polygon Polygon::as_type(int type, int N, int D)
{
  return Polygon(N, D, type);
}

inline bool Polygon::set_param(std::vector<double> &param, int idx, double val)
{
  if (idx >= 0 && idx < (int)param.size()) {
    param[idx] = val;
    return true;
  }
  else
    return false;
}

inline double Polygon::get_param(std::vector<double> &param, int idx)
{
  if (idx >= 0 && idx < (int)param.size())
    return param[idx];
  else
    return NAN;
}

inline bool Polygon::set_radius(int idx, double val)
{
  return set_param(radius, idx, val);
}

inline double Polygon::get_radius(int idx) { return get_param(radius, idx); }

inline bool Polygon::set_height(int idx, double val)
{
  return set_param(height, idx, val);
}

inline double Polygon::get_height(int idx) { return get_param(height, idx); }

inline bool Polygon::set_edge(int idx, double val)
{
  return set_param(edge, idx, val);
}

inline double Polygon::get_edge(int idx) { return get_param(edge, idx); }

inline bool Polygon::set_twist_angle(int idx, double val)
{
  return set_param(twist_angle, idx, val);
}

inline double Polygon::get_twist_angle(int idx)
{
  return get_param(twist_angle, idx);
}

inline bool Polygon::value_is_set(double val) { return !std::isnan(val); }

} // namespace anti

#endif // POLYGON_H
