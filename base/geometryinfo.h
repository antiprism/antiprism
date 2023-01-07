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

/*!\file geometryinfo.h
 * \brief Anayse and get information about a geometry
 */

#ifndef GEOMETRYINFO_H
#define GEOMETRYINFO_H

#include "geometry.h"
#include "geometryutils.h"

namespace anti {

/// Less, for angles
struct AngleLess {
  /// Less, for angles
  /**\param ang1 first angle to compare
   * \param ang2 second angle to compare
   * \return \c true if first angle is less than second angle,
   * otherwise \c false.*/
  bool operator()(const double &ang1, const double &ang2) const;
};

/// Less, for vectors of angles
struct AngleVectLess {
  /// Less, for vectors of angles
  /**\param angs1 first vector of angles to compare
   * \param angs2 second angle to compare
   * \return \c true if first angle is less than second angle,
   * otherwise \c false.*/
  bool operator()(const std::vector<double> &angs1,
                  const std::vector<double> &angs2) const;
};

/// Collection of limit values for some measure for elements of some type
/** This is simple store for the information. Processing is done externally.
 *  Not all limit types will be available for all measures. Edges are held by
 *  their two vertex index numbers, e.g. \c idx[IDX_MIN] and \c idx[IDX_MIN2].*/
class ElementLimits {
public:
  /// Index numbers to select limit elements
  enum {
    IDX_MIN = 0, ///> Select element with minimum value
    IDX_MIN2,    ///> Select element with minimum value (for edges)
    IDX_MAX,     ///> Select element with maximum value
    IDX_MAX2,    ///> Select element with maximum value (for edges)
    IDX_ZERO,    ///> Select element closest to zero value
    IDX_ZERO2    ///> Select element closest to zero value (for edges)
  };

  int idx[6];  ///> Limit element indexes, select with \c IDX_MIN, etc
  double max;  ///> Maximum value for the elements
  double min;  ///> Minimum value for the elements
  double zero; ///> Value nearest to zero for the elements
  double sum;  ///> Sum of the values from all the elements

  /// Constructor
  ElementLimits();

  /// Initialiase, if reusing the object
  void init();

  /// Check whether values have been set
  /**\return \c true if set, otherwise \c false.*/
  bool is_set() const;
};

struct double_range_cnt {
  int cnt;
  double min;
  double max;

  double_range_cnt() : cnt(0), min(1e100), max(-1e100) {}
  double_range_cnt update(double val)
  {
    cnt++;
    if (val < min)
      min = val;
    if (val > max)
      max = val;
    return *this;
  }
  double mid() const { return (min + max) / 2; } // middle of range
  double rad() const { return (max - min) / 2; } // radius of range
};

/// Find values and properties of a geometry
/** The properties and values are cached, including any intermediate values,
 *  (some calculations find several assoociated properties).
 *  Properties are returned by reference, and may be accessed directly
 *  using the calling function without them being recalculated each
 *  time they are accessed.*/
class GeometryInfo {
private:
  Vec3d cent;
  int oriented;
  int orientable;
  bool found_connectivity;
  bool closed;
  bool polyhedron;
  bool known_connectivity;
  bool even_connectivity;
  int number_parts;
  int genus_val;
  ElementLimits iedge_len;
  ElementLimits edge_len;
  ElementLimits so_angles;
  ElementLimits dih_angles;
  ElementLimits ang;
  int num_angs;
  ElementLimits area;
  double vol;
  Vec3d vol_cent;
  ElementLimits v_dists;
  ElementLimits e_dists;
  ElementLimits ie_dists;
  ElementLimits f_dists;

  std::vector<std::vector<int>> impl_edges;
  std::map<std::vector<int>, std::vector<int>> efpairs;
  std::vector<std::vector<int>> edge_parts;
  std::map<std::vector<double>, int, AngleVectLess> face_angles;
  std::map<std::vector<double>, int, AngleVectLess> vert_dihed;
  std::map<std::vector<int>, int> edge_index_numbers;
  std::map<double, double_range_cnt, AngleLess> dihedral_angles;
  std::map<double, double_range_cnt, AngleLess> e_lengths;
  std::map<double, double_range_cnt, AngleLess> ie_lengths;
  std::map<double, double_range_cnt, AngleLess> plane_angles;
  std::map<double, double_range_cnt, AngleLess> sol_angles;
  std::vector<double> vertex_angles;
  std::map<std::pair<int, int>, double> vf_plane_angles;
  std::vector<double> edge_dihedrals;
  std::vector<double> f_areas;
  std::vector<double> f_perimeters;
  std::vector<double> f_max_nonplanars;
  std::vector<std::vector<int>> vert_cons;
  std::vector<std::vector<int>> vert_cons_orig;
  std::vector<std::vector<int>> vert_faces;
  std::vector<std::vector<int>> vert_impl_edges;
  std::vector<std::vector<std::vector<int>>> face_cons;
  std::vector<std::vector<std::vector<int>>> vert_figs;
  std::vector<Vec3d> vert_norms;
  bool vert_norms_local_orient;
  std::vector<int> free_verts;
  bool found_free_verts;
  Geometry dual;
  Symmetry sym;

  void find_impl_edges();
  void find_edge_index_numbers();
  void find_edge_face_pairs();
  void find_edge_parts();
  void find_connectivity();
  void find_face_angles();
  void find_dihedral_angles();
  void find_vert_cons();
  void find_vert_cons_orig();
  void find_face_cons();
  void find_vert_figs();
  void find_vert_elems(const std::vector<std::vector<int>> &elems,
                       std::vector<std::vector<int>> &vert_elems);
  void find_vert_norms(bool local_orient = false);
  void find_free_verts();
  void find_solid_angles();
  void find_e_lengths(std::map<double, double_range_cnt, AngleLess> &e_lens,
                      const std::vector<std::vector<int>> &edges,
                      ElementLimits &lens);
  void find_f_areas();
  void find_f_perimeters();
  void find_f_max_nonplanars();
  void find_oriented();
  void find_v_dist_lims();
  void find_e_dist_lims();
  void find_ie_dist_lims();
  void find_f_dist_lims();
  void find_symmetry();

protected:
  const Geometry &geom;

public:
  /// Constructor
  /**\param geo geometry to get information about
   * \param center used for any properties that are relative to a centre */
  GeometryInfo(const Geometry &geo, Vec3d center = Vec3d(0, 0, 0));

  /// Reset, clear all setting
  void reset();

  /// Get the geometry being analysed
  /**\return The geometry.*/
  const Geometry &get_geom() const;

  /// Get the centre used for calculations
  /** This is the centre used for calculations that are relative to
   *  a centre, and is not a property of the geometry.
   * \return The centre.*/
  Vec3d get_center() const;

  /// Set the centre used for calculations
  /** This centre is used for calculations that are relative to
   *  a centre.
   * \param center the centre*/
  void set_center(Vec3d center);

  // elements
  // ----------------------------------------------------------------

  /// Get number of vertices
  /**\return The number of vertices.*/
  int num_verts() const;

  /// Get number of edges
  /**\return The number of edges.*/
  int num_edges() const;

  /// Get number of implicit edges
  /** Edges lying on some face
   * \return The number of implicit edges.*/
  int num_iedges();

  /// Get number of faces
  /**\return The number of faces.*/
  int num_faces() const;

  /// Get number of parts
  /** Parts consist of sets of faces joined by edges, where any two
   *  faces are connected by a chain of edge connected faces.
   * \return The number of face parts.*/
  int num_parts();

  // ----------------------------------------------------------------

  /// Check if oriented
  /**\return \c true if oriented, otherwise \c false.*/
  bool is_oriented();

  /// Check if orientable
  /**\return \c true if orientable, otherwise \c false.*/
  bool is_orientable();

  /// Check if closed
  /** Each edge is met by at least two faces
   * \return \c true if closed, otherwise \c false.*/
  bool is_closed();

  /// Check if polyhedron
  /** Each edge is met by exactly two faces.
   * \return \c true if polyhedron, otherwise \c false.*/
  bool is_polyhedron();

  /// Check if even connectivity
  /** Each edge is met by an even number of faces.
   * \return \c true if even connectivity, otherwise \c false.*/
  bool is_even_connectivity();

  /// Check if known connectivity
  /** Each edge is met by a maximum of two faces.
   * \return \c true if know connectivity, otherwise \c false.*/
  bool is_known_connectivity();

  /// Check if genus is known
  /**\return \c true if genus could be calculated, otherwise \c false.*/
  bool is_known_genus();

  /// Get genus
  /**\return The genus, if the geometry is orientable, otherwise
   *  the unorientable genus (demigenus).*/
  int genus();

  /// Get the volume centroid
  /** Note: the result will only be valid for models that enclose volume
   * \return The coordinates of the volume centroid. */
  Vec3d volume_centroid();

  // ----------------------------------------------------------------
  // Limits: largest and smallest by some measure

  /// Get face area limits
  /** Elements are faces.
   * \return The face area limits.*/
  ElementLimits face_areas();

  /// Get edge length limits
  /** Elements are edges, identified by two vertex index numbers.
   * \return The edge length limits.*/
  ElementLimits edge_length_lims();

  /// Get implicit edge length limits
  /** Elements are edges, identified by two vertex index numbers.
   * \return The implicit edge length limits.*/
  ElementLimits iedge_length_lims();

  /// Get dihedral angle limits
  /** Elements are edges, identified by two vertex index numbers.
   * \return The dihedral angle limits.*/
  ElementLimits dihed_angle_lims();

  /// Get solid angle limits
  /** Elements are vertices.
   * \return The solid angle limits.*/
  ElementLimits solid_angle_lims();

  /// Get vertex distance limits
  /** Distance from centre to vertex.
   *  Elements are vertices.
   * \return The vertex distance limits.*/
  ElementLimits &vert_dist_lims();

  /// Get edge distance limits
  /** Distance from centre to line through edge
   *  Elements are edges, identified by two vertex index numbers.
   * \return The edge distance limits.*/
  ElementLimits &edge_dist_lims();

  /// Get implicit edge distance limits
  /** Distance from centre to line through edge.
   *  Elements are edges, identified by two vertex index numbers.
   * \return The implicit edge distance limits.*/
  ElementLimits &iedge_dist_lims();

  /// Get facee distance limits
  /** Distance from centre to plane through face.
   *  Elements are faces.
   * \return The face distance limits.*/
  ElementLimits &face_dist_lims();

  /// Get face angle limits
  /** Angles between a vertex and its two neighbours on a face.
   *  No elements are returned, only values
   * \return The face angle limits.*/
  ElementLimits angle_lims();

  // ----------------------------------------------------------------

  /// Get signed volume
  /** The geometry should consist of oriented polyhedra for the result
   *  to be thought of as a volume.
   * \return The signed volume.*/
  double volume();

  /// Get isoperimetric quotient
  /** The isoperimetric quotient is a measure of sphericity, and is
   *  independant of scale. It is defined as 36*PI*Vol^2/Area^3, and
   *  has a maximum value of 1, correspomding to a perfect sphere.
   * \return The isoperimetric quotient.*/
  double isoperimetric_quotient();

  /// Get the number of face angles
  /** This is also half the number of implicit edges
   * \return The number of face angles.*/
  int num_angles();

  /// Get the angle defect
  /** Defined as  <tt>2PI * number_vertices - sum_of_face_angles</tt> \n
   *  This will be 4PI for a convex polyhedron (or one connected
   *  like a sphere).
   * \return The angle defect.*/
  double angle_defect();

  // ----------------------------------------------------------
  // Vertices

  /// Get vertex normals
  /** These are calculated from the faces that surround the vertex,
   *  as the average of the unit normals of these faces
   * \param local_orient \c true indicates that the faces that
   *  surround the vertex may not be consistently oriented, and that
   *  the calculation should correct for this. \c false indicates
   *  that the given face normals should be used, which is an
   *  optimisation if the geometry is known to be oriented.
   * \return Vertx normals.*/
  const std::vector<Vec3d> &get_vert_norms(bool local_orient = true);

  /// Get vertex connections
  /**Get, for each vertex, the vertices connected to it by
   * an edge, either implicit (face edhe) or explicit.
   * \return The vertices connected to each vertex.*/
  const std::vector<std::vector<int>> &get_vert_cons();

  /// Get vertex figures
  /**A vertex may lie on several faces, and some of these faces may be
   * joined to each other and surround the vertex, implying a circuit of
   * neigbouring vertices that surround the vertex. These circuits may
   * be informally thought of as vertex figures. There may be several
   * such circuits around a vertex. Note: the vertex figure can, in
   * general, be any sort of graph, and won't always break down into
   * disjoint circuits.
   * \return The vertex figure circuits for each vertex.*/
  const std::vector<std::vector<std::vector<int>>> &get_vert_figs();

  /// Get, vertex implicit edges
  /**Get, for each vertex, the implicit edges it is a part of
   * (in numeric order)
   * \return The implicit edges connected to each vertex.*/
  const std::vector<std::vector<int>> &get_vert_impl_edges();

  /// Get, vertex faces
  /**Get, for each vertex, the faces it is a part of (in numeric order)
   * \return The faces connected to each vertex.*/
  const std::vector<std::vector<int>> &get_vert_faces();

  /// Get free verts
  /** Free vertices are vertices that are not part of any face
   *  or explicit edge.
   * \return The free vertices.*/
  const std::vector<int> &get_free_verts();

  /// Get solid vertex angles
  /** \return The solid angle ant each vertex, in steradians.*/
  const std::vector<double> &get_vert_solid_angles();

  /// Get solid angles by size
  /** \return The sizes of solid angles (steradians) and the number
   * of vertices having an angle of that size.*/
  const std::map<double, double_range_cnt, AngleLess> &
  get_solid_angles_by_size();

  /// Get plane angles
  /** The plane angles are the angles between neighbouring edges on a
   * face.
   * \return The plane angles, indexed by a vertex face pair (in that
   * order).*/
  const std::map<std::pair<int, int>, double> &get_plane_angles();

  /// Get plane angles by size
  /** \return The sizes of plane angles and the number of angles of
   * each size.*/
  const std::map<std::vector<double>, int, AngleVectLess> &
  get_plane_angles_by_size();

  // ----------------------------------------------------------
  // Edges

  /// Get index number of edge between two vertices
  /** Look up the edge index in the edge list
   * \param v_idx0 first vertex index number
   * \param v_idx1 second vertex index number
   * \return The edge index number, or -1 if not found. */
  int get_edge_index(int v_idx0, int v_idx1);

  /// Get edge face pairs
  /** An edge has two vertices, but may be part of any number of faces,
   * \return A map of the vertex pair of an edge to the faces it lies on.*/
  const std::map<std::vector<int>, std::vector<int>> &get_edge_face_pairs();

  /// Get the dihedral angle at each edge
  /**\return The dihedral angles.*/
  const std::vector<double> &get_edge_dihedrals();

  /// Get edge parts
  /** Edge parts consist of sets of edges joined by vertices, where any two
   *  edges are connected by a chain of vertex connected edges.
   * \return The edge parts.*/
  const std::vector<std::vector<int>> &get_edge_parts();

  /// Get dihedral angles by size
  /**\return The sizes of dihedral angles, and the number of edges having
   * each size.*/
  const std::map<double, double_range_cnt, AngleLess> &
  get_dihedral_angles_by_size();

  /// Get explicit edge lengths by size
  /**\return The lengths of edges, and the number of edges having
   * each length.*/
  const std::map<double, double_range_cnt, AngleLess> &
  get_edge_lengths_by_size();

  /// Get implicit edges
  /** Implicit edges are edges that lie on faces.
   * \return The implicit edges.*/
  const std::vector<std::vector<int>> &get_impl_edges();

  /// Get explicit edge lengths by size
  /**\return The lengths of edges, and the number of edges having
   * each length.*/
  const std::map<double, double_range_cnt, AngleLess> &
  get_iedge_lengths_by_size();

  // ----------------------------------------------------------
  // Faces

  /// Get face areas
  /**\return The face areas.*/
  const std::vector<double> &get_f_areas();

  /// Get face perimeters
  /**\return The face perimeters.*/
  const std::vector<double> &get_f_perimeters();

  /// Get face non-planarity value
  /** The value is the maximum height of a vertex above or below
   *  the plane of the face.
   * \return The face non-planarity values.*/
  const std::vector<double> &get_f_max_nonplanars();

  /// Get face connections
  /**Get, for each face, the faces connected to it across an edge.
   * \return The faces connected to each face.*/
  const std::vector<std::vector<std::vector<int>>> &get_face_cons();

  /// Get a dual
  /**\return A dual.*/
  const Geometry &get_dual();

  // ----------------------------------------------------------
  // Symmetry

  /// Get the symmetry group
  /**\return The symmetry group.*/
  const Symmetry &get_symmetry();

  /// Get the name of the symmetry type
  /** This will be the Schoenflies name
   * \return The name of the symmetry type.*/
  std::string get_symmetry_type_name();

  /// Get the symmetry axes
  /** These are the all axes for a rotational symmetry subgroup
   * \return The symmetry axes.*/
  const std::set<SymmetryAxis> &get_symmetry_axes();

  /// Get the subsymmetries
  /**\return The subsymmetries.*/
  const std::set<Symmetry> &get_symmetry_subgroups();

  /// Get the symmetry automorphisms
  /**These leave the symmetry group unchanged, by conjugation.
   * \return The symmetry automorphism.*/
  const SymmetryAutos &get_symmetry_autos();

  /// Get the transformation to align with the standard symmetry group
  /**\return The transformation.*/
  Trans3d get_symmetry_alignment_to_std();
};

} // namespace anti

#endif // GEOMETRYINFO_H
