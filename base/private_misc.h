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
   Name: private_misc.h
   Description: miscellaneous functions
   Project: Antiprism - http://www.antiprism.com
*/

#ifndef PRIVATE_MISC_H
#define PRIVATE_MISC_H

#include <map>
#include <string>
#include <vector>

#include "geometry.h"

namespace anti {

// triangulate.cc
int triangulate(Geometry &geom, Color inv = Color(),
                unsigned int winding_rule = TESS_WINDING_NONZERO,
                std::vector<int> *fmap = nullptr);
void triangulate_basic(Geometry &geom, bool sq_diag = true, Color inv = Color(),
                       std::vector<int> *fmap = nullptr);

Status add_hull(Geometry &geom, std::string qh_args = "", int *dim = nullptr);
Status set_hull(Geometry &geom, std::string qh_args = "", int *dim = nullptr);

/// Get Voronoi cells.
/**Get all Voronoi cells of the vertex points which are finite polyhedra.
 * \param geom geometry to find the Voronoi cells for
 * \param cells to return the Voronoi cells
 * \param qh_args additional arguments to pass to qhull (unsupported,
 * may not work, check output.)
 * \return status, evaluates to \c true if the cells were calculated,
 * otherwise false.*/
Status get_voronoi_cells(Geometry &geom, std::vector<Geometry> &cells,
                         std::string qh_args = "");

/// Get Delaunay edges.
/**Get the edges of the Delaunay tesselation.
 * \param geom geometry to find the Dealunay edges for
 * \param edges the Delaunay edges, as a pairs of vertex index numbers,
 * mapped to the number of Delaunay tetrahedra the edge was part of
 * \param qh_args additional arguments to pass to qhull (unsupported,
 * may not work, check output.)
 * \return status, evaluates to \c true if the edges were calculated,
 * otherwise false.*/
Status get_delaunay_edges(const Geometry &geom,
                          std::map<std::pair<int, int>, int> &edges,
                          std::string qh_args);

double get_min_vert_to_vert_dist(const std::vector<Vec3d> &verts,
                                 double sig_dist);

int orient_geom(Geometry &geom, std::vector<std::vector<int>> *parts = nullptr);
bool orient_geom(Geometry &geom, int type, char *errmsg = nullptr);
void orient_reverse(Geometry &geom);
void get_star(const Geometry &geom, std::vector<Vec3d> &star, char type = 'v',
              Vec3d centre = Vec3d(0, 0, 0));
Status make_zono(Geometry &zono, const std::vector<Vec3d> &star);

void orient_face(std::vector<int> &face, int v0, int v1);

} // namespace anti

#endif // PRIVATE_MISC_H
