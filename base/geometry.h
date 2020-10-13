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

/*!\file geometry.h
 * \brief Classes to represent a geometry
 */

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <map>
#include <string>
#include <vector>

#include "elemprops.h"
#include "status.h"
#include "trans3d.h"
#include "vec_utils.h"

namespace anti {

class GeometryInfo;

/// Geometry Interface
class Geometry {
private:
  std::vector<Vec3d> vert_elems;
  std::vector<std::vector<int>> face_elems;
  std::vector<std::vector<int>> edge_elems;

  GeomElemProps<Color> cols;

public:
  /// Constructor
  Geometry() = default;

  /// Copy Constructor
  /** Initialise from another geometry that implements \c Geometry
   * \param geom geometry to copy from.*/
  Geometry(const Geometry &geom) = default;

  /// Copy Assignment
  /** Copy from another geometry that implements \c Geometry
   * \param geom geometry to copy from.
   * \return A reference to this object.*/
  Geometry &operator=(const Geometry &geom) = default;

  /// Destructor
  virtual ~Geometry() = default;

  /// Check whether geometry is set.
  /**\return \c true if the geometry is set, otherwise \c false. */
  bool is_set() const { return verts().size() > 0; }

  //-------------------------------------------
  // Element Access
  //-------------------------------------------

  /// Read access to the vertices.
  /**\return A reference to vertex coordinates. */
  virtual const std::vector<Vec3d> &verts() const;

  /// Read/Write access to the vertices.
  /**\return A reference to vertex coordinates. */
  virtual std::vector<Vec3d> &raw_verts();

  /// Read access to a vertex.
  /**\param v_idx index number of the vertex.
   * \return A reference to the vertex coordinates. */
  virtual const Vec3d &verts(int v_idx) const;

  /// Read/Write access to a vertex.
  /**\param v_idx index number of the vertex.
   * \return A reference to the vertex coordinates. */
  virtual Vec3d &verts(int v_idx);

  /// Read access to the edges.
  /**\return A reference to the edge data. */
  virtual const std::vector<std::vector<int>> &edges() const;

  /// Read/Write access to the edges.
  /**\return A reference to the edge data. */
  virtual std::vector<std::vector<int>> &raw_edges();

  /// Read access to an edge.
  /**\param e_idx index number of the edge.
   * \return A reference to the edge data. */
  virtual const std::vector<int> &edges(int e_idx) const;

  /// Get the index number of a vertex of an edge.
  /**\param e_idx edge index number.
   * \param v_no the position the vertex appears in the edge, \c 0 or \c 1
   * \return The vertex index number.
   */
  virtual int edges(int e_idx, int v_no) const;

  /// Get the coordinates of a vertex of an edge.
  /**\param e_idx edge index number.
   * \param v_no the position the vertex appears in the edge, \c 0 or \c 1
   * \return The vertex coordinates.
   */
  Vec3d edge_v(int e_idx, int v_no) const;

  /// Read access to the faces.
  /**\return A reference to the face data. */
  virtual const std::vector<std::vector<int>> &faces() const;

  /// Read/Write access to the faces.
  /**\return A reference to the face data. */
  virtual std::vector<std::vector<int>> &raw_faces();

  /// Read access to a face.
  /**\param f_idx index number of the face.
   * \return A reference to the face data. */
  virtual const std::vector<int> &faces(int f_idx) const;

  /// Read/Write access to a face.
  /**\param f_idx index number of the face.
   * \return A reference to the face data. */
  virtual std::vector<int> &faces(int f_idx);

  /// Get the vertex index number of a face vertex *face vertex in range).
  /**\param f_idx face index number.
   * \param v_no the position the vertex appears in the face,
   * \c 0, \c 1, \c 2, ...
   * \return The vertex index number. */
  virtual int faces(int f_idx, int v_no) const;

  /// Get the vertex index number of a face vertex (face vertex in or out of
  /// range).
  /**\param f_idx face index number.
   * \param v_no the position the vertex appears in the face,
   * \c 0, \c 1, \c 2, ...
   * \return The vertex index number. */
  virtual int faces_mod(int f_idx, int v_no) const;

  /// Get the coordinates of a face vertex (face vertex in range)
  /**\param f_idx face index number.
   * \param v_no the position the vertex appears in the face,
   * \c 0, \c 1, \c 2, ...
   * \return The vertex coordinates.
   */
  Vec3d face_v(int f_idx, int v_no) const;

  /// Get the coordinates of a face vertex (face vertex in or out of range).
  /**\param f_idx face index number.
   * \param v_no the position the vertex appears in the face,
   * \c 0, \c 1, \c 2, ...
   * \return The vertex coordinates.
   */
  Vec3d face_v_mod(int f_idx, int v_no) const;

  //-------------------------------------------
  // Add and Delete Elements
  //-------------------------------------------
  /// Add a vertex with colour
  /**\param vert vertex to add.
   * \param col colour of the vertex.
   * \return index number of newly added vertex. */
  virtual int add_vert(Vec3d vert, Color col = Color());

  /// Add several vertices without colour.
  /**\param vrts vertices to add.
   * \return index number of last added vertex. */
  virtual int add_verts(const std::vector<Vec3d> &vrts);

  /// Add an edge with colour
  /**\param edge edge to add.
   * \param col colour of the edge.
   * \return index number of newly added edge. */
  virtual int add_edge(std::vector<int> edge, Color col = Color());

  /// Add an edge from vertex index numbers, with colour
  /**\param v_idx1 index number of first vertex.
   * \param v_idx2 index number of second vertex.
   * \param col colour of the edge.
   * \return index number of newly added edge.*/
  virtual int add_edge(int v_idx1, int v_idx2, Color col = Color());

  /// Add several edges, without colour
  /** Each edge is added only if it is not already present in the edge list.
   * \param edgs edges to add.
   * \return index number of last added edge in edge list. */
  virtual int add_edges(const std::vector<std::vector<int>> &edgs);

  /// Add an edge without checking it is valid, with colour
  /** This function is much faster than \c add_edge().
   *  The edge should have the lowest vertex index first. The edge
   *  should not already be included in the edge list.
   * \param edge edge to add.
   * \param col colour of the edge.
   * \return index number of newly added edge. */
  virtual int add_edge_raw(const std::vector<int> &edge, Color col = Color());

  /// Add several edges without checking they are valid, without colour
  /** This function is much faster than \c add_edges().
   *  The edges should have the lowest vertex index first. The edges
   *  should not already be included in the edge list.
   * \param edgs edges to add.
   * \return index number of last added edge. */
  virtual int add_edges_raw(const std::vector<std::vector<int>> &edgs);

  /// Add a face, with colour
  /**\param face face to add.
   * \param col colour of the face.
   * \return index number of newly added face. */
  virtual int add_face(const std::vector<int> &face, Color col = Color());

  /// Add a face made from integer arguments, without colour
  /**The final argument must be a dummy value of \c -1.
   * \param v1 the first index
   * \param ... further indexes in the face, and a final -1 to terminate.
   * \return index number of newly added face. */
  virtual int add_face(int v1, ...);

  /// Add several faces, without colour
  /**\param fces faces to add.
   * \return index number of last added face. */
  virtual int add_faces(const std::vector<std::vector<int>> &fces);

  /// Delete an element
  /** When deleting several vertices is is more efficient to use
   *  \c delete_verts().
   * \param type from VERTS, EDGES, FACES.
   * \param idx vertex index number to delete
   * \param elem_map a map of old index numbers to new index numbers,
   *  deleted elements map to index \c -1.  */
  virtual void del(int type, int idx, std::map<int, int> *elem_map = nullptr);

  /// Delete several elements
  /**\param type from VERTS, EDGES, FACES.
   * \param idxs vertex index numbers to delete
   * \param elem_map a map of old index numbers to new index numbers,
   *  deleted elements map to index \c -1.  */
  virtual void del(int type, const std::vector<int> &idxs,
                   std::map<int, int> *elem_map = nullptr);

  /// Delete all elements of a type
  virtual void clear(int type);

  /// Delete all elements of all types
  virtual void clear_all();

  /// Get colours for an element type
  /**\param type element type from VERTS, EDGES, FACES.
   * \return the colour. */
  const ElemProps<Color> &colors(int type) const;
  ElemProps<Color> &colors(int type);

  /// Get colours
  /**\return the colour. */
  const GeomElemProps<Color> &get_cols() const;

  /// Get colours
  /**\return the colour. */
  GeomElemProps<Color> &get_cols();

  /// Clear all colours for all elements
  void clear_cols();

  /// Append a geometry
  /** Include the elements of a geometry after the current set of elements
   * \param geom geometry to append.*/
  virtual void append(const Geometry &geom);

  /// Merge some of the vertices
  /**\param vmap a map from each vertex index number to the index
   *  number of the vertex it will be replaced with */
  virtual void verts_merge(std::map<int, int> &vmap);

  /// Add missing implicit edges
  /** Add implicit edges (edges of faces) to the edge list if they are
   *  not already included.
   * \param col colour for the edges (default: no colour set). */
  virtual void add_missing_impl_edges(const Color &col = Color());

  //-------------------------------------------
  // Transformations
  //-------------------------------------------

  /// Add the faces of the convex hull to the geometry.
  /**\param qh_args additional arguments to pass to qhull (unsupported,
   *  may not work, check output.)
   * \param dim dimension of the hull 3, 2, 1 or 0.
   * \return status, which evaluates to \c true if qhull could
   *  calculate the hull(possibly with warnings), otherwise \c false
   *  to indicate an error. */
  Status add_hull(std::string qh_args = "", int *dim = nullptr);

  /// Set the geometry to its convex hull.
  /** If the convex hull could not be calculated the the geometry
   *  will be left empty.
   * \param qh_args additional arguments to pass to qhull (unsupported,
   *  may not work, check output.)
   * \param dim dimension of the hull 3, 2, 1 or 0.
   * \return status, which evaluates to \c true if qhull could
   *  calculate the hull(possibly with warnings), otherwise \c false
   *  to indicate an error. */
  Status set_hull(std::string qh_args = "", int *dim = nullptr);

  /// Orient the geometry (if possible.)
  /**\param parts used to return the index numbers of the faces
   *  in each of the disconnected parts (if \a parts is not \c 0 .)
   * \return the number of disconnected parts. */
  int orient(std::vector<std::vector<int>> *parts = nullptr);

  /// Orient the geometry (if possible.)
  /**\param type type of orientation, 1 - positive volume, 2 negative
   *  volume, 3 - reverse first face and orient, 4 - reverse faces without
   *  orienting.
   * \return status, which evaluates to \c true if the orientation was
   *  successful (possibly with warnings), otherwise \c false to indicate
   *  an error. */
  Status orient(int type);

  /// Reverse the orientation of all the faces.
  void orient_reverse();

  /// Apply a transformation matrix.
  /**\param trans the transformation matrix. */
  void transform(const Trans3d &trans);

  /// Align a polyhedron with the standard alignment for its symmetry type.
  void sym_align();

  /// Triangulate (tesselate) the faces
  /** Divide the faces into triangles, adding new vertices as necessary,
   *  and coloured edges if a colour is set for new elements.
   * \param col the colour for any new edges or vertices that were added.
   *  The default leaves them with default colours. If it is set to
   *  Color::invisible then the new elements are given a colour that
   *  indicates that they should not be displayed.
   * \param winding selects resulting faces acording to winding number
   *    TESS_WINDING_ODD
   *    TESS_WINDING_NONZERO (default)
   *    TESS_WINDING_POSITIVE
   *    TESS_WINDING_NEGATIVE
   *    TESS_WINDING_ABS_GEQ_TWO
   * \param fmap a vector to return the face mapping. Each old face
   *  index maps to the first index of faces it was converted to. A
   *  final index holds the total number of new faces */
  void triangulate(Color col = Color(),
                   unsigned int winding = TESS_WINDING_NONZERO,
                   std::vector<int> *fmap = nullptr);

  //-------------------------------------------
  // Geometric Utilities
  //-------------------------------------------

  /// Get the centroid of all the vertices.
  /**\return The coordinates of the centroid. */
  Vec3d centroid() const;

  /// Get all the face centroids.
  /**\param ctds vector used to return the coordinates of the centroids. */
  void face_cents(std::vector<Vec3d> &ctds) const;

  /// Get all the face normals.
  /**\param norms vector used to return the normals.
   * \param allow_zero if \c true then the magnitude of the normal
   *  is the area of the face, if \c false then for faces with near-zero
   *  area a normal with a more accurate direction is calculated. */
  void face_norms(std::vector<Vec3d> &norms, bool allow_zero = false) const;

  /// Get a face centroid, for a face index.
  /**\param f_idx face index number.
   * \return The coordinates of the centroid. */
  Vec3d face_cent(int f_idx) const;

  /// Get a face centroid, for a face.
  /**\param face contains the vertex index numbers.
   * \return The coordinates of the centroid. */
  Vec3d face_cent(const std::vector<int> &face) const;

  /// Get a face normal, for a face index.
  /**\param f_idx the face index number.
   * \param allow_zero if \c true then the magnitude of the normal
   *  is the area of the face, if \c false then for faces with near-zero
   *  area a normal with a more accurate direction is calculated.
   * \return The normal */
  Vec3d face_norm(int f_idx, bool allow_zero = false) const;

  /// Get a face normal, for a face.
  /**\param face contains the vertex index numbers in the face.
   * \param allow_zero if \c true then the magnitude of the normal
   *  is the area of the face, if \c false then for faces with near-zero
   *  area a normal with a more accurate direction is calculated.
   * \return The normal */
  Vec3d face_norm(const std::vector<int> &face, bool allow_zero = false) const;

  /// Get a face area, for a face index.
  /**\param f_idx the face index number.
   * \return The signed area */
  double face_area(int f_idx) const;

  /// Get a face area, for a face.
  /**\param face contains the vertex index numbers in the face.
   * \return The signed area */
  double face_area(const std::vector<int> &face) const;

  /// Get the nearest point on a face to another point, for a face index.
  /**\param f_idx face index number.
   * \param P the point
   * \return The coordinates of the nearest point to \c P on the
   *  face plane. */
  Vec3d face_nearpt(int f_idx, Vec3d P) const;

  /// Get the nearest point on a face to another point, for a face.
  /**\param face contains the vertex index numbers.
   * \param P the point
   * \return The coordinates of the nearest point to \c P on the
   *  face plane. */
  Vec3d face_nearpt(const std::vector<int> &face, Vec3d P) const;

  /// Get the angles and edge lengths around a face
  /**\param f_idx face index number.
   * \param angles to return the plane angles, starting with the angle
   *  at face vertex 0.
   * \param lengths to return the edge lengths, starting with the edge
   *  leading from face vertex 0.*/
  void face_angles_lengths(int f_idx, std::vector<double> *angles,
                           std::vector<double> *lengths = nullptr) const;

  /// Get an edge centre (the centroid), for an edge index.
  /**\param e_idx edge index number.
   * \return The coordinates of the centre. */
  Vec3d edge_cent(int e_idx) const;

  /// Get an edge centre (the centroid), for an edge.
  /**\param edge contains the two vertex index numbers.
   * \return The coordinates of the centre. */
  Vec3d edge_cent(const std::vector<int> &edge) const;

  /// Get the vector between the ends of an edge, for an edge index.
  /**\param e_idx edge index number.
   * \return The vector. */
  Vec3d edge_vec(int e_idx) const;

  /// Get the vector between the ends of an edge, for an edge.
  /**\param edge contains the two vertex index numbers.
   * \return The vector. */
  Vec3d edge_vec(const std::vector<int> &edge) const;

  /// Get the vector between two vertices by vertex index.
  /**\param v_idx0 start vertex index number.
   * \param v_idx1 end vertex index number.
   * \return The vector. */
  Vec3d edge_vec(int v_idx0, int v_idx1) const;

  /// Get the nearest point on an edge to another point, for a edge index.
  /**\param e_idx edge index number.
   * \param P the point
   * \return The coordinates of the nearest point to \c P on the
   * line of the edge. */
  Vec3d edge_nearpt(int e_idx, Vec3d P) const;

  /// Get the nearest point on an edge to another point, for an edge.
  /**\param edge contains the two vertex index numbers.
   * \param P the point
   * \return The coordinates of the nearest point to \c P on the
   *  line of the edge. */
  Vec3d edge_nearpt(const std::vector<int> &edge, Vec3d P) const;

  /// Get the length of an edge, for an edge index.
  /**\param e_idx edge index number.
   * \return The length of the edge. */
  double edge_len(int e_idx) const;

  /// Get the length of an edge, for an edge.
  /**\param edge contains the two vertex index numbers.
   * \return The length of the edge. */
  double edge_len(const std::vector<int> &edge) const;

  //-------------------------------------------
  //  Other Utilities
  //-------------------------------------------

  /// Read geometry from a file
  /** The file is first read as a normal OFF file, if that fails it will be
   *  read as a Qhull formatted OFF file, and if that fails the file will be
   *  read for any coordinates (lines that contains three numbers separated
   *  by commas and/or spaces will be taken as a set of coordinates.)
   * \param file_name the file name ("" or "-" for standard input).
   * \return status, which evaluates to \c true if the file could be read
   *  (possibly with warnings), otherwise \c false to indicate an error. */
  virtual Status read(std::string file_name = "");

  /// Read geometry from a file stream
  /** The file is first read as a normal OFF file, if that fails it will be
   *  read as a Qhull formatted OFF file, and if that fails the file will be
   *  read for any coordinates (lines that contains three numbers separated
   *  by commas and/or spaces will be taken as a set of coordinates.)
   * \param file the file stream.
   * \return status, which evaluates to \c true if the file could be read
   *  (possibly with warnings), otherwise \c false to indicate an error. */
  virtual Status read(FILE *file);

  /// Read resource geometry from resource model name
  /**\param res_name the resource model name.
   * \return status, which evaluates to \c true if the resource could be
   *  read (possibly with warnings), otherwise \c false to indicate an
   *  error. */
  virtual Status read_resource(std::string res_name = "");

  /// Write geometry to a file
  /**\param file_name the file name ("" for standard output.)
   * \param sig_dgts the number of significant digits to write,
   *  or if negative then the number of digits after the decimal point.
   * \return status, which evaluates to \c true if the file could be written
   *  (possibly with warnings), otherwise \c false to indicate an error. */
  virtual Status write(std::string file_name = "",
                       int sig_dgts = DEF_SIG_DGTS) const;

  /// Write geometry to a file stream
  /**\param file the file stream.
   * \param sig_dgts the number of significant digits to write,
   *  or if negative then the number of digits after the decimal point. */
  virtual void write(FILE *file, int sig_dgts = DEF_SIG_DGTS) const;

  /// Write coordinates to a file
  /**\param file_name the file name ("" for standard output.)
   * \param sep a string to use as the seperator between coordinates.
   * \param sig_dgts the number of significant digits to write,
   *  or if negative then the number of digits after the decimal point.
   * \return status, which evaluates to \c true if the file could be written
   *  (possibly with warnings), otherwise \c false to indicate an error. */
  virtual Status write_crds(std::string file_name = "", const char *sep = " ",
                            int sig_dgts = DEF_SIG_DGTS) const;

  /// Write coordinates to a file stream
  /**\param file the file stream.
   * \param sep a string to use as the seperator between coordinates.
   * \param sig_dgts the number of significant digits to write,
   *  or if negative then the number of digits after the decimal point. */
  virtual void write_crds(FILE *file, const char *sep = " ",
                          int sig_dgts = DEF_SIG_DGTS) const;

  /// Write coordinates to a file
  /**\param file_name the file name ("" for standard output.)
   * \param mtl_file is the materials file if specified.
   * \param sep a string to use as the seperator between coordinates.
   * \param sig_dgts the number of significant digits to write,
   * or if negative then the number of digits after the decimal point.
   * \return status, which evaluates to \c true if the file could be written
   *  (possibly with warnings), otherwise \c false to indicate an error. */
  virtual Status write_obj(std::string file_name = "",
                           std::string mtl_file = "", const char *sep = " ",
                           int sig_dgts = DEF_SIG_DGTS) const;

  /// Write coordinates to a file stream
  /**\param file the file stream.
   * \param mfile stream for materials file.
   * \param mtl_file is the materials file name.
   * \param sep a string to use as the seperator between coordinates.
   * \param sig_dgts the number of significant digits to write,
   * or if negative then the number of digits after the decimal point. */
  virtual void write_obj(FILE *file, FILE *mfile, std::string mtl_file = "",
                         const char *sep = " ",
                         int sig_dgts = DEF_SIG_DGTS) const;

  /// Check if geomtery is consistently oriented
  /**\return \c true if consistently oriented, otherwise \c false. */
  bool is_oriented() const;

  /// Get GeometryInfo object
  /**\return GeometryInfo object associated with this geometry. */
  GeometryInfo get_info() const;

  /// Get implicit edges
  /** Returns the edges of the polygon faces
   * \param edgs the edges are returned here */
  void get_impl_edges(std::vector<std::vector<int>> &edgs) const;

  /// Get faces lying on each side of an edge for all edges
  /**\param oriented \c true the geometry is oriented and the first face in
   *  the face pair has the vertices of the edge pair in order, -1 for
   *  a face index indicates the edge is open to that sid. \c false the
   * value is a list of all face sharing the edge,in any order.
   * \return edge2facepr return map of edges to face pairs.*/
  std::map<std::vector<int>, std::vector<int>>
  get_edge_face_pairs(bool oriented = true) const;
};

/// Make a face from vertex index numbers
/** The final argument must be a dummy value of \c -1.
 * \param v1 the first index
 * \param ... further indexes in the face, and a final -1 to terminate.
 * \return index number of newly added face. */
std::vector<int> make_face(int v1, ...);

/// Make an edge from vertex index numbers
/**\param v_idx1 index number of first vertex.
 * \param v_idx2 index number of second vertex.
 * \return edge containg index numbers in numerical order. */
std::vector<int> make_edge(int v_idx1, int v_idx2);

// Implementation of inline functions

// -------------------------------------------------------------------
// Geometry::

inline Vec3d Geometry::edge_v(int e_idx, int v_no) const
{
  return verts(edges(e_idx, v_no));
}

inline int Geometry::faces_mod(int f_idx, int v_no) const
{
  unsigned int f_sz = faces(f_idx).size();
  if (v_no < 0)
    v_no = f_sz - ((-v_no) % f_sz);
  else
    v_no = v_no % f_sz;
  return faces(f_idx, v_no);
}

inline Vec3d Geometry::face_v(int f_idx, int v_no) const
{
  return verts(faces(f_idx, v_no));
}

inline Vec3d Geometry::face_v_mod(int f_idx, int v_no) const
{
  return verts(faces_mod(f_idx, v_no));
}

inline Vec3d Geometry::centroid() const { return anti::centroid(verts()); }

inline void Geometry::face_cents(std::vector<Vec3d> &ctds) const
{
  ctds.resize(faces().size());
  for (unsigned int i = 0; i < faces().size(); i++)
    ctds[i] = face_cent(i);
}

inline void Geometry::face_norms(std::vector<Vec3d> &norms,
                                 bool allow_zero) const
{
  norms.resize(faces().size());
  for (unsigned int i = 0; i < faces().size(); i++)
    norms[i] = face_norm(i, allow_zero).unit();
}

inline Vec3d Geometry::face_cent(int f_idx) const
{
  return face_cent(faces(f_idx));
}

inline Vec3d Geometry::face_cent(const std::vector<int> &face) const
{
  return anti::centroid(verts(), face);
}

inline Vec3d Geometry::face_norm(int f_idx, bool allow_zero) const
{
  return face_norm(faces(f_idx), allow_zero);
}

inline Vec3d Geometry::face_norm(const std::vector<int> &face,
                                 bool allow_zero) const
{
  return anti::face_norm(verts(), face, allow_zero);
}

inline double Geometry::face_area(int f_idx) const
{
  return face_norm(f_idx, true).len();
}

inline double Geometry::face_area(const std::vector<int> &face) const
{
  return face_norm(face, true).len();
}

inline Vec3d Geometry::face_nearpt(int f_idx, Vec3d P) const
{
  return face_nearpt(faces(f_idx), P);
}

inline Vec3d Geometry::face_nearpt(const std::vector<int> &face, Vec3d P) const
{
  return nearest_point(P, verts(), face);
}

inline Vec3d Geometry::edge_cent(int e_idx) const
{
  return edge_cent(edges(e_idx));
}

inline Vec3d Geometry::edge_cent(const std::vector<int> &edge) const
{
  return 0.5 * (verts(edge[1]) + verts(edge[0]));
}

inline Vec3d Geometry::edge_vec(int e_idx) const
{
  return edge_vec(edges(e_idx));
}

inline Vec3d Geometry::edge_vec(const std::vector<int> &edge) const
{
  return verts(edge[1]) - verts(edge[0]);
}

inline Vec3d Geometry::edge_vec(int v_idx0, int v_idx1) const
{
  return verts(v_idx1) - verts(v_idx0);
}

inline Vec3d Geometry::edge_nearpt(int e_idx, Vec3d P) const
{
  return edge_nearpt(edges(e_idx), P);
}

inline Vec3d Geometry::edge_nearpt(const std::vector<int> &edge, Vec3d P) const
{
  return nearest_point(P, verts(), edge);
}

inline double Geometry::edge_len(int e_idx) const
{
  return edge_vec(e_idx).len();
}

inline double Geometry::edge_len(const std::vector<int> &edge) const
{
  return edge_vec(edge).len();
}

inline void Geometry::transform(const Trans3d &trans)
{
  anti::transform(raw_verts(), trans);
}

// -------------------------------------------------------------------
// Geometry::

inline std::vector<Vec3d> &Geometry::raw_verts() { return vert_elems; }

inline const std::vector<Vec3d> &Geometry::verts() const { return vert_elems; }

inline const Vec3d &Geometry::verts(int v_idx) const
{
  return vert_elems[v_idx];
}

inline Vec3d &Geometry::verts(int v_idx) { return vert_elems[v_idx]; }

inline std::vector<std::vector<int>> &Geometry::raw_edges()
{
  return edge_elems;
}

inline const std::vector<std::vector<int>> &Geometry::edges() const
{
  return edge_elems;
}

inline const std::vector<int> &Geometry::edges(int e_idx) const
{
  return edge_elems[e_idx];
}

inline int Geometry::edges(int e_idx, int v_no) const
{
  return edge_elems[e_idx][v_no];
}

inline std::vector<std::vector<int>> &Geometry::raw_faces()
{
  return face_elems;
}

inline const std::vector<std::vector<int>> &Geometry::faces() const
{
  return face_elems;
}

inline const std::vector<int> &Geometry::faces(int f_idx) const
{
  return face_elems[f_idx];
}

inline std::vector<int> &Geometry::faces(int f_idx)
{
  return face_elems[f_idx];
}

inline int Geometry::faces(int f_idx, int v_no) const
{
  return face_elems[f_idx][v_no];
}

inline const ElemProps<Color> &Geometry::colors(int type) const
{
  return cols[type];
}

inline ElemProps<Color> &Geometry::colors(int type) { return cols[type]; }

inline const GeomElemProps<Color> &Geometry::get_cols() const { return cols; }

inline GeomElemProps<Color> &Geometry::get_cols() { return cols; }

inline void Geometry::clear_cols() { cols.clear(); }

} // namespace anti

#endif // GEOMETRY_H
