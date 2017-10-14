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

/*!\file liling.h
   \brief Generate polyhedra based on polygons.
*/

#ifndef TILING_H
#define TILING_H

#include "status.h"

#include <string>
#include <utility>
#include <vector>

namespace anti {

class Tile {
private:
  std::vector<int> idxs; // index numbers of points visited
  std::vector<int> ops;  // operations between visiting points
  char start_faces;

  mutable std::vector<int>::const_iterator ops_i;
  mutable std::vector<int>::const_iterator idxs_i;

public:
  enum { END = -1, V = 0, E, F, VE, EF, FV, VEF, P };

  Status read(const std::string &pat);
  unsigned char get_start_faces() const { return start_faces; }

  void start_op() const;
  void next_op() const;
  int get_op() const;
  int get_idx() const;

  void relabel(std::vector<int> relab);
  std::vector<int> check_index_range(int num_points) const;
  std::string tile_string();
};

class Tiling {
private:
  std::vector<std::pair<Vec3d, Color>> points;
  std::vector<Tile> pat_paths;

  Geometry meta;
  std::vector<std::vector<int>> nbrs;
  std::vector<Vec3d> vert_norms;

  bool one_of_each_tile;

  bool find_nbrs();
  int get_vert_idx(int tri, const std::pair<Vec3d, Color> &point) const;
  Vec3d point_on_face(int f_idx, const Vec3d &crds) const;
  void
  add_circuit(Geometry &geom, int start_idx, const Tile &pat,
              std::vector<bool> &seen, Color col,
              const std::vector<std::map<std::vector<int>, std::pair<int, int>>>
                  &index_order,
              const std::vector<int> &point_vertex_offsets) const;
  const std::vector<Tile> &get_pat_paths() const { return pat_paths; }

public:
  Tiling() : one_of_each_tile(false) {}
  Status set_geom(const Geometry &geom, bool is_meta = false,
                  double face_ht = 0.0);
  Status add_tile(const std::string &pat);
  Status make_tiling(Geometry &geom,
                     std::vector<int> *tile_counts = nullptr) const;
  Status read_pattern(const std::string &pat);
  Status relabel_pattern(const std::string &relabel);
  Status read_conway(const std::string &op);
  const Geometry &get_meta() const { return meta; }
  void set_one_of_each_tile(bool val = true) { one_of_each_tile = val; };

  std::string pattern_string();
  void print_conway_list(FILE *ofile = stdout);
};

/// Make a tiling of a base model using Wythoff constructive notation
/**The base model must be orientated. A constructive notation pattern
 * is applied to the base model to make a new model. This model is coloured
 * by index. The faces are coloured by tile number in pattern, the vertices
 * are coloured by elements used in weighting V=0, E, F, VE, EF, FV, VEF
 * \param tiled_geom, to return the final model
 * \param base geom, the oriented polyhedron (or tiling) to be processed
 * \param pat the pattern to apply, in Wythoff constructive notation (or
 *  specifying a Conway Notation operator if \c pat_is_conway_op is \c true)
 * \param pat_is_conway_op if \c true the pattern is a Conway Notation
 *  operator, typically a single letter possibly followed by an integer.
 * \return status, which evaluates to \c true if the geometry and
 *  pattern were valid, otherwise \c false to indicate an error. */
Status wythoff_make_tiling(Geometry &tiled_geom, const Geometry &base_geom,
                           const std::string &pat,
                           bool pat_is_conway_op = false);

/// Get vertex points of a Schwarz triangle, and its symmetry group
/**\param fracs six integers, taken in pairs as the angle fractions.
 * \param verts the three vertices of the triangle.
 * \param sym to return the symmetry group.
 * \return \c true if the triangle is valid, therwise \c false. */
bool get_schwarz_tri_verts(const std::vector<int> &fracs,
                           std::vector<Vec3d> &verts, Symmetry *sym = nullptr);

/// Get angle fractions for a Schwarz triangle by index number
/**\param tri_idx ranges from 0 - 43 (see base/wythoff.cc for the list).
 * \param fracs six integers, taken in pairs as the angle fractions.
 * \return \c true if the triangle index is valid, therwise \c false. */
bool get_schwarz_tri_fracs(int tri_idx, std::vector<int> &fracs);

}; // namespace anti

#endif // TILING_H
