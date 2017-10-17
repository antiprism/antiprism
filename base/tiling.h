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

/*!\file tiling.h
   \brief Generate polyhedra based on polygons.
*/

#ifndef TILING_H
#define TILING_H

#include "status.h"

#include <string>
#include <utility>
#include <vector>

namespace anti {

// Tile used in Wythoff constructive notation
class Tile {
private:
  std::vector<int> idxs; // index numbers of points visited
  std::vector<int> ops;  // operations between visiting points
  char start_faces;      // start tile on: '+', '-' or '*' (both) face types

  mutable std::vector<int>::const_iterator ops_i;
  mutable std::vector<int>::const_iterator idxs_i;

public:
  enum {
    END = -1, ///< Inidicates tile has been read
    V = 0,    ///< Vertex mirror operation, or point lying on vertex V
    E,        ///< Edge mirror operation, or point lying on vertex E
    F,        ///< Face mirror operation, or point lying on vertex F
    VE,       ///< Point lying in edge VE
    EF,       ///< Point lying in edge EF
    FV,       ///< Point lying in edge FE
    VEF,      ///< Point lying inside trinagle VEF
    P         ///< Plot point operation
  };

  /// Read a tile pattern string
  /**\param pat tile pattern
   * \return Status, which evaluates to \c true if the pattern was
   *  successfully read, otherwise \c false to indicate an error. */
  Status read(const std::string &pat);

  /// Get the start face types where the for the tile
  /**\return Start face types '+' positive, '-' negavie or '*' both. */
  unsigned char get_start_faces() const { return start_faces; }

  /// Start reading operations
  void start_op() const;

  /// Move to next operation
  void next_op() const;

  /// Get the operation
  int get_op() const;

  /// If current operation is P, get the point index number
  int get_idx() const;

  /// Relabel the tile pattern by a permutation of VEF
  /**\param relab relabelling permutation VEF -> . */
  void relabel(std::vector<int> relab);

  /// Check point index numbers are within rage
  /** \param num_points the number of points
   *  \return The index numbers which were out of range. */
  std::vector<int> check_index_range(int num_points) const;

  /// Convert the tile to a string representation
  std::string tile_string();
};

// Tiling using Wythoff constructive notation
class Tiling {
private:
  std::vector<std::pair<Vec3d, Color>> points; ///< Points
  std::vector<Tile> pat_paths;                 ///< Tile patterns

  Geometry meta;                      ///< Base triangle tiling
  std::vector<std::vector<int>> nbrs; ///< Base tiling face neighbours
  // std::vector<Vec3d> vert_norms;            ///< Base tiling vertex normals

  bool one_of_each_tile; ///< Only plot one tile per kind

  /// Find base tiling face neighbours
  /** \return Status, which evaluates to \c true if the nieghbours were
   *  found successfully, otherwise \c false to indicate an error. */
  bool find_nbrs();

  /// Add a circuit (face) for an individual tile pattern
  /**\param geom the geometry to add the circuit face to
   * \param start_index the base triangle to start the circuit
   * \param seen base triangles that have already been used
   * \param index order used to calculate final vertex index numbers
   * \param point_vertex_offsets used to calculate final vertex index numbers */
  void
  add_circuit(Geometry &geom, int start_idx, const Tile &pat,
              std::vector<bool> &seen, Color col,
              const std::vector<std::map<std::vector<int>, std::pair<int, int>>>
                  &index_order,
              const std::vector<int> &point_vertex_offsets) const;
  /// Get the tile patterns
  /** \return The tile patterns. */
  const std::vector<Tile> &get_pat_paths() const { return pat_paths; }

public:
  /// Constructor
  Tiling() : one_of_each_tile(false) {}

  /// Set the base geometry
  /**\param geom the base geometry
   * \param is_meta the base geometry is already a meta-like tiling
   *  and should be used as-is
   * \param face_ht rais the F vertex by this amount above the face
   * \return Status, which evaluates to \c true if the base geometry was
   *  successfully set, otherwise \c false to indicate an error. */
  Status set_geom(const Geometry &geom, bool is_meta = false,
                  double face_ht = 0.0);

  /// Add a tile
  /**\param pat the tile pattern string
   * \return Status, which evaluates to \c true if the pattern string was
   *  successfully added, otherwise \c false to indicate an error. */
  Status add_tile(const std::string &pat);

  /// Make the tiling
  /**\param geom the geometry to return the tiling
   * \param tile_counts to return how many tiles of each kind were created
   * \return Status, which evaluates to \c true if the tiling was
   *  successfully created, otherwise \c false to indicate an error. */
  Status make_tiling(Geometry &geom,
                     std::vector<int> *tile_counts = nullptr) const;
  /// Read tiling pattern
  /**\param pat the tiling pattern string
   * \return Status, which evaluates to \c true if the pattern string was
   *  successfully read, otherwise \c false to indicate an error. */
  Status read_pattern(const std::string &pat);

  /// Relabel pattern (permute VEF in pattern string)
  /**\param relabel a three character string permution of the letters VEF
   * \return Status, which evaluates to \c true if the permutation was
   *  valid, otherwise \c false to indicate an error. */
  Status relabel_pattern(const std::string &relabel);

  /// Read Conway operation string
  /**\param op a single Conway operation
   * \return Status, which evaluates to \c true if the operation was
   *  valid, otherwise \c false to indicate an error. */
  Status read_conway(const std::string &op);

  /// Get the base meta tiling
  /**\return The base meta tiling */
  const Geometry &get_meta() const { return meta; }

  /// Set that only tile of each kind should be plotted
  /**\param val \c true, only one of each kind, otherwise \c false, all */
  void set_one_of_each_tile(bool val = true) { one_of_each_tile = val; };

  /// Get a Wythoff constructive notation string representing the pattern
  /**\return The pattern string. */
  std::string pattern_string();

  /// Print a list of the Conway notation operators supported
  /**\param ofile print to this stream */
  void print_conway_list(FILE *ofile = stdout);
};

/// Make a tiling of a base model using Wythoff constructive notation
/**The base model must be orientated. A constructive notation pattern
 * is applied to the base model to make a new model. This model is coloured
 * by index. The faces are coloured by tile number in pattern, the vertices
 * are coloured by elements used in weighting V=0, E, F, VE, EF, FV, VEF
 * \param tiled_geom to return the final model
 * \param base_geom the oriented polyhedron (or tiling) to be processed
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
