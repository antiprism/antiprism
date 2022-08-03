/*
   Copyright (c) 2003-2022, Adrian Rossiter, Roger Kaufman

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
   Name: canonical_common.h
   Description: canonical code shared in /src
   Project: Antiprism - http://www.antiprism.com
*/

#ifndef CANONICAL_COMMON_H
#define CANONICAL_COMMON_H

#include <cstdio>
#include <string>
#include <vector>

using std::string;
using std::vector;

#include "../base/antiprism.h"

// for canonical.cc and conway.cc

/// find nearpoints radius, sets range minimum and maximum
/**\param geom geometry.
 * \param min returns the minimum nearpoints radius.
 * \param max returns the maximum nearpoints radius.
 * \param center returns the cente rof the nearpoints.
 * \returns the average radius of the nearpoints. */
double edge_nearpoints_radius(const anti::Geometry &geom, double &min,
                              double &max, anti::Vec3d &center);

/// wrapper for above.
/**\param geom geometry. */
double edge_nearpoints_radius(const anti::Geometry &geom);

/// sets radius of geom to average of edge near points radius
/**\param geom geometry. */
void unitize_nearpoints_radius(anti::Geometry &geom);

/// return true if maximum vertex radius is radius_range_percent (0.0 to ...)
/**greater than minimum vertex radius (visible for canonical.cc)
 * \param geom geometry to measure.
 * \param radius_range_percent limit to maximum radius over minimum radius */
bool canonical_radius_range_test(const anti::Geometry &geom,
                                 const double radius_range_percent);

/// returns the edge near points centroid
/**\param geom geometry to measure
 * \param cent centre from which to calculate nearpoints on edges
 * \return the centroid of the nearpoints. */
anti::Vec3d edge_nearpoints_centroid(anti::Geometry &geom,
                                     const anti::Vec3d cent = anti::Vec3d(0, 0,
                                                                          0));

/// Canonicalize (George Hart "Conway Notation" algorithm)
/**See http://www.georgehart.com/virtual-polyhedra/conway_notation.html
 * \param base geometry to canonicalise.
 * \param it_ctrl interation control.
 * \param radius_range_percent if the model outer radius increases this
 *  much over the inner radius then it is growing too much, terminate.
 * \param planarize_only planarise only.
 * \return \c true if success, otherwise \c false */
bool canonicalize_bd(anti::Geometry &base, anti::IterationControl it_ctrl,
                     double radius_range_percent, const bool planarize_only);

/// an abbreviated wrapper for planarization with the base/dual method
/**\param geom geometry to planarize.
 * \param it_ctrl interation control.
 * \return \c true if success, otherwise \c false */
bool planarize_bd(anti::Geometry &geom, anti::IterationControl it_ctrl);

#endif // CANONICAL_COMMON_H
