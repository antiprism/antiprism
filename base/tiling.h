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

#include "color.h"
#include "colormap.h"
#include "geometry.h"
#include "symmetry.h"

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
    VEF,      ///< Point lying inside triangle VEF
    P         ///< Plot point operation
  };

  struct TileReport {
    int count;             // count of tiles generated for a path
    std::string step;      // steps to reach associated element (mirrors vef)
    std::string assoc;     // associated element type (mirrors vef)
    std::string step_back; // steps to return from assoc'd element (mirrors vef)
    int assoc_type;        // associated element type
  };

  /// Read a tile pattern string
  /**\param pat tile pattern
   * \return Status, which evaluates to \c true if the pattern was
   *  successfully read, otherwise \c false to indicate an error. */
  Status read(const std::string &pat);

  /// Get the start face types for the tile
  /**\return Start face types '+' positive, '-' negative or '*' both. */
  unsigned char get_start_faces() const { return start_faces; }

  /// Set the start face types for the tile
  /**\param start start face types: '+' positive, '-' negative or '*' both. */
  void set_start_faces(unsigned char start = '*') { start_faces = start; }

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

  // Flip start faces
  /**Flip +/-, * is left unchanged. */
  void flip_start_faces();

  /// Check point index numbers are within rage
  /** \param num_points the number of points
   *  \return The index numbers which were out of range. */
  std::vector<int> check_index_range(int num_points) const;

  /// Convert the tile to a string representation
  std::string tile_string() const;

  /// Get element association
  /**\return vector with four strings:
   *         1 - steps to reach associated element (mirrors vef)
   *         2 - associated element type (mirrors vef)
   *         3 - steps to return from associated element (mirrors vef)
   *         4 - associated element type letter (VEF or X for none) */
  TileReport get_element_association() const;

  /// Normalize pattern (deduplicate operators, remove trailing _)
  /**\param pattern pattern string to mormalize
   * \return normalized pattern string. */
  static std::string normalize_pattern(const std::string &pattern);
};

class TilingPoint {
public:
  TilingPoint(Vec3d crds = Vec3d(), int idx = -1) : coords(crds), index(idx)
  {
    // val {1, 2, 4, 3, 6, 5, 7, 0} = v, e, f, ve, ef, fv, vef, 0 }
    const int incl2idx[] = {7, 0, 1, 3, 2, 5, 4, 6}; // v,e,f,ve,ef,fv,vef,0
    inclusion = incl2idx[(coords[0] != 0) + (coords[1] != 0) * 2 +
                         (coords[2] != 0) * 4];
  }
  const Vec3d &get_coords() const { return coords; }
  int get_index() const { return index; }
  int get_inclusion() const { return inclusion; }

private:
  Vec3d coords;
  int inclusion;
  int index;
};

class TilingColoring {
public:
  enum class ColoringType {
    none,               ///< do not colour tile/point)
    index,              ///< colour tile/point by index number
    associated_element, ///< colour tile same as corresponding base element
    component,          ///< colour point by presence of VEF components
    weight,             ///< colour point by barycentric RGB colour
  };
  enum { TILES = 0, POINTS = 1, ELEMS_SZ = 2 };

  /// Constructor
  /**\param str colouring method for tiles, can be given as initial substring
   * none: do not colour tiles
   * index, value (default): use the path index
   * association: colour tiles using corresponding base element
   *   colour, optionally follow by comma and
   *      element type (from VEF): to colour local tiles by that element
   *                               type (default: F)
   *     colour: to colour all local tiles by that colour */
  TilingColoring(const std::string &clrng_str = std::string())
  {
    if (clrng_str.empty() || !read_coloring(clrng_str)) {
      for (int elem : {TILES, POINTS})
        set_coloring(elem, ColoringType::index);
    }
  }

  /// Read tile colouring
  /**\param clrng_str colouring method for tiles, can be given as
   * initial substring
   *  none: do not colour tiles
   *  index, value (default): use the path index
   *  association: colour tiles using corresponding base element
   *    colour, optionally follow by comma and
   *      element type (from VEF): to colour local tiles by that element
   *                               type (default: F)
   *      colour: to colour all local tiles by that colour
   * \return Status, which evaluates to \c true if the coloring was
   *  valid, otherwise \c false to indicate an error. */
  // Status read_tile_coloring(const std::string &clrng_str);
  Status read_coloring(const std::string &clrng_str);

  /// Read point colouring
  /**\param clrng_str colouring method for points, can be given as
   * initial substring
   *  none: do not colour points
   *  index, value: use the point index
   *  component: colour points using the presence of VEF components (default)
   *  weight: colour points using barycentric VEF components as RGBA
   * \return Status, which evaluates to \c true if the coloring was
   *  valid, otherwise \c false to indicate an error. */
  // Status read_point_coloring(const std::string &clrng_str);

  /// Get option help
  /**\param op_char the option letter
   * \return formatted help for option. */
  static std::string get_option_help(char op_char);

  /// Check if tiling colouring is by path index
  /**\param elem 0 - tile, 1 - point
   * \return \c true if tiling colouring is path index, otherwise \c false. */
  bool is_index(int elem) const
  {
    return coloring_types[elem] == ColoringType::index;
  }

  /// Check if tiling colouring is by associated element
  /**\param elem 0 - tile, 1 - point
   * \return \c true if tiling colouring is associated element,
   *  otherwise \c false. */
  bool is_associated_element(int elem) const
  {
    return coloring_types[elem] == ColoringType::associated_element;
  }

  /// Get associations to use for local elements
  /**\param elem 0 - tile, 1 - point
   * \return \c Tile::V - vertices, \c Tile::E - edges, \c Tile::F - faces,
   *  \c Tile::P - the tiling colouring \c color.*/
  int get_local_association(int elem) const { return local_associations[elem]; }

  /// Get the colour to use for local elements
  /**\param elem 0 - tile, 1 - point
   * \return the colour to use for local tiles with associated colouring,
   *  or test with \c is_index() or \c is_value() for tile index colouring. */
  const Color &get_local_color(int elem) const { return local_colors[elem]; }

  // ColoringType get_point_coloring_type() const { return point_clrng_type;}
  Color get_point_color(const TilingPoint &point) const;

  /// Get default element colormap
  /**\return the default colormap to use for tiles (created with new, nust
   * be deleted after use, or added to a \c Coloring which will manage the
   * deletion) */
  ColorMap *get_default_colormap(int elem) const;

private:
  // The set functions are private, as no external code uses them
  Status set_coloring(int elem, ColoringType clrng_type, int assoc = -1,
                      Color local = Color());
  /*  Status set_tile_coloring(ColoringType clrng_type, int elem = -1,
                           Color local = Color());

    Status set_point_coloring(ColoringType clrng_type, int elem = -1,
                              Color local = Color());
  */
  std::vector<ColoringType> coloring_types =
      std::vector<ColoringType>(ELEMS_SZ, ColoringType::index);
  std::vector<int> local_associations = std::vector<int>(ELEMS_SZ, Tile::F);
  std::vector<Color> local_colors = std::vector<Color>(ELEMS_SZ);
};

// Tiling using Wythoff constructive notation
class Tiling {
private:
  std::vector<TilingPoint> points; ///< Points
  ElemProps<Color> orig_colors;    ///< Original colours
  std::vector<Tile> pat_paths;     ///< Tile patterns
  TilingColoring coloring;         ///< Coloring for points, tiles

  Geometry meta;                      ///< Base triangle tiling
  std::vector<std::vector<int>> nbrs; ///< Base tiling face neighbours
  // std::vector<Vec3d> vert_norms;            ///< Base tiling vertex normals

  bool one_of_each_tile; ///< Only plot one tile per kind

  /// Find base tiling face neighbours
  /** \return Status, which evaluates to \c true if the nieghbours were
   *  found successfully, otherwise \c false to indicate an error. */
  bool find_nbrs();

  /// Find the colour of the element associated with a tiling vertex
  /**\param f_idx the index of the meta tiling face
   * \param incl the element types included in the pattern point specifier
   * \param coloring the tiling coloring
   * \return colour of the associated element. */
  Color get_associated_element_point_color(int f_idx, int incl) const;

  /// Find the element associated with a circuit
  /**\param start_idx the index of the start face
   * \param step the sequence of mirrors vef to step to the associated meta
   *  triangle
   * \param type the element type, from V, E, F or VEF for no element
   * \param local_color colour for interior points
   * \return index of the associated element, or -1 if no element. */
  int get_associated_element(int start_idx, const std::string &step,
                             int assoc_type) const;

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

  /// Set coloring
  /**\param clrng the coloring to use for points and tiles */
  void set_coloring(const TilingColoring &clrng) { coloring = clrng; }

  /// Get coloring
  /**\return the coloring to use for points and tiles */
  const TilingColoring &get_coloring() const { return coloring; }

  /// Get default tile colormap
  /**\return the default colormap to use for tiles (created with new, nust
   * be deleted after use, or added to a \c Coloring which will manage the
   * deletion) */
  ColorMap *get_default_tile_colormap() const
  {
    return coloring.get_default_colormap(TilingColoring::TILES);
  }

  /// Get default point colormap
  /**\return the default colormap to use for tiles (created with new, nust
   * be deleted after use, or added to a \c Coloring which will manage the
   * deletion) */
  ColorMap *get_default_point_colormap() const
  {
    return coloring.get_default_colormap(TilingColoring::POINTS);
  }

  /// Add a tile
  /**\param pat the tile pattern string
   * \return Status, which evaluates to \c true if the pattern string was
   *  successfully added, otherwise \c false to indicate an error. */
  Status add_tile(const std::string &pat);

  /// Make the tiling
  /**\param geom the geometry to return the tiling
   * \param tile_reports reports about the paths and the tiles produced
   * \return Status, which evaluates to \c true if the tiling was
   *  successfully created, otherwise \c false to indicate an error. */
  Status
  make_tiling(Geometry &geom,
              std::vector<Tile::TileReport> *tile_reports = nullptr) const;

  /// Read tiling pattern
  /**\param pat the tiling pattern string
   * \return Status, which evaluates to \c true if the pattern string was
   *  successfully read, otherwise \c false to indicate an error. */
  Status read_pattern(const std::string &pat);

  /// Relabel pattern (permute VEF in pattern string)
  /**\param relabel a three character string permution of the letters VEF
   * \return Status, which evaluates to \c true if the permutation was
   *  valid, otherwise \c false to indicate an error. */
  Status relabel_pattern(const std::string relabel);

  /// Read Conway operation string
  /**\param op a single Conway operation
   * \return Status, which evaluates to \c true if the operation was
   *  valid, otherwise \c false to indicate an error. */
  Status read_conway(const std::string &op);

  /// Get the base meta tiling
  /**\return The base meta tiling */
  const Geometry &get_meta() const { return meta; }

  /// Reverse pattern, switch current tile start faces
  void reverse_pattern();

  /// Set current tile start faces to both (*)
  void start_everywhere();

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
 * \param pat the pattern to apply, in Wythoff constructive notation or
 *  by specifying a single Conway Notation operator (typically a single
 *  letter possibly followed by an integer.)
 * \param oriented geom is oriented, otherwise tiles use all start triangles
 *  and result may not be a polyhedron.
 * \param reverse reverse the pattern, flip start triangles
 * \param col_type type of colouring
 * \return status, which evaluates to \c true if the geometry and
 *  pattern were valid, otherwise \c false to indicate an error. */
Status wythoff_make_tiling(Geometry &tiled_geom, const Geometry &base_geom,
                           const std::string &pat, bool oriented = true,
                           bool reverse = false,
                           TilingColoring col_type = TilingColoring());

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
