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

/*!\file geometryutils.h
 * \brief Utilities associated with geometries
*/

#ifndef GEOMETRYUTILS_H
#define GEOMETRYUTILS_H

#include "coloring.h"
#include "normal.h"
#include "symmetry.h"

namespace anti {
class GeometryInfo;

/// Get Voronoi cells.
/**Get all finite Voronoi cells of the vertex points of a geometry
 * \param verts the vertices to find the Voronoi cells for
 * \param cells to return the Voronoi cells
 * \param qh_args additional arguments to pass to qhull (unsupported,
 * may not work, check output.)
 * \return status, evaluates to \c true if the cells were calculated,
 * otherwise false.*/
Status get_voronoi_cells(const std::vector<Vec3d> &verts,
                         std::vector<Geometry> *cells,
                         std::string qh_args = "");

/// Get a star of vectors to use for making a zonohedron.
/**\param geom geometry to get the star from
 * \param type the type of star to make can be
 * <ul>
 * <li>v - centre to vertices are vectors (default)
 * <li>a - all vertex to vertex are vectors
 * <li>i - implicit edges (face sides) are vectors
 * <li>e - explicit edges are vectors
 * </ul>
 * \param centre the centre to use when the type is v.
 * \return The star of vectors. */
std::vector<Vec3d> get_star(const Geometry &geom, char type = 'v',
                            Vec3d centre = Vec3d(0, 0, 0));

/// Make a zonohedron from a star
/**\param geom to return the zonohedron
 * \param star the star of vectors.
 * \return status, which evaluates to true if the zonohedron could be
 *  calculated (possibly with warnings), otherwise \c false to indicate
 *  an error. */
Status make_zonohedron(Geometry &geom, const std::vector<Vec3d> &star);

/// Make a zonohedrified polyhedron from a seed polyhedron and a star
/**\param geom to return the zonohedrified polyhedron
 * \param seed the seed polyhedron
 * \param star the star of vectors.
 * \param col the colour for the faces in the generated zones
 * \return status, which evaluates to true if the zonohedron could be
 *  calculated (possibly with warnings), otherwise \c false to indicate
 *  an error. */
Status make_zonohedrified_polyhedron(Geometry &geom, const Geometry &seed,
                                     const std::vector<Vec3d> &star, Color col);

/// Make a polar zonohedron from an ordered star
/**\param geom to return the polar zonohedron
 * \param star the ordered star of vectors.
 * \param step step this many places to get to the next vector in the star
 * \return status, which evaluates to true if the zonohedron could be
 *  calculated (possibly with warnings), otherwise \c false to indicate
 *  an error. */
Status make_polar_zonohedron(Geometry &geom, const std::vector<Vec3d> &star,
                             int step = 1);

/// Set planar geodesic division.
/** A Class I pattern is made with m=0,n=1. A Class II pattern
 *  is made with m=1,n=1. A Class III pattern is made with n>=1,m>1.
 * \param geom to return the planar geodesic polyhedron
 * \param base the base polyhedron
 * \param m the first pattern specifier.
 * \param n the second pattern specifier.
 * \return \c true if the pattern was valid, otherwise \c false. */
bool make_geodesic_planar(Geometry &geom, const Geometry &base, int m,
                          int n = 0);

/// Set spherical geodesic division.
/** A Class I pattern is made with m=0,n=1. A Class II pattern
 *  is made with m=1,n=1. A Class III pattern is made with n>=1,m>1.
 * \param geom to return the geodesic sphere
 * \param base the base polyhedron
 * \param m the first pattern specifier.
 * \param n the second pattern specifier.
 * \param cent the centre of projection.
 * \return \c true if the pattern was valid, otherwise \c false. */
bool make_geodesic_sphere(Geometry &geom, const Geometry &base, int m,
                          int n = 0, Vec3d cent = Vec3d(0, 0, 0));

/// Project the vertices onto a sphere
/**\param geom whose vertices will be projected
 * \param centre the centre of the sphere.
 * \param radius the radius of the sphere. */
void project_onto_sphere(Geometry &geom, Vec3d centre = Vec3d(0, 0, 0),
                         double radius = 1.0);

/// Combine face circuits into a single face by bridging between the circuits
/// with double edges
/**\param contours a set of circuits that comprise a single composite face
 * \param new_edges to return the bridges between the circuits
 * \return the face, as a single circuit. */
std::vector<int>
make_face_from_contours(const std::vector<std::vector<int>> &contours,
                        std::vector<std::vector<int>> *new_edges = nullptr);

/// Repeat a part by a set of symmetry transformations
/**\param geom geometry to return the final model.
 * \param part geometry to be repeated.
 * \param ts transformations to be used for the repeats.
 * \param col_part_elems element types, combining flags ELEM_VERTS,
 *   ELEM_EDGES an ELEM_FACES, to be coloured, based on the
 *   order position of the transformation that produced them.
 * \param clrngs an array of three Colorings applied, correspondingly, to the
 *  index coloured vertices, edges and faces.*/
void sym_repeat(Geometry &geom, const Geometry &part, const Transformations &ts,
                char col_part_elems = ELEM_NONE, Coloring *clrngs = nullptr);

/// Repeat a part by a set of symmetry transformations
/**\param geom geometry to return the final model.
 * \param part geometry to be repeated.
 * \param sym symmetry group transformations to be used for the repeats.
 * \param col_part_elems element types, combining flags ELEM_VERTS,
 *   ELEM_EDGES an ELEM_FACES, to be coloured, based on the
 *   order position of the transformation that produced them.
 * \param clrngs an array of three Colorings applied, correspondingly, to the
 *  index coloured vertices, edges and faces.
 * \return \c true if the symmetry group was valid, otherwise \c false. */
bool sym_repeat(Geometry &geom, const Geometry &part, const Symmetry &sym,
                char col_part_elems = ELEM_NONE, Coloring *clrngs = nullptr);

/// Repeat a part by a set of symmetry transformations
/**\param geom geometry to return the final model.
 * \param sym_to target symmetry.
 * \param sym_from initial symmetry.
 * \param pos to realign sym_from. */
void transform_and_repeat(Geometry &geom, std::string sym_to,
                          std::string sym_from, Trans3d pos = Trans3d());

/// Convert the face planes to vertices by polar reciprocation
/** \param dual output geometry containing polar vertices.
 * \param geom input geometry containing faces.
 * \param recip_rad radius of reciprocation sphere.
 * \param centre centre of reciprocation sphere.
 * \param inf maximum distance a vertex will be placed. */
void get_pol_recip_verts(Geometry &dual, const Geometry &geom, double recip_rad,
                         Vec3d centre, double inf = 1e15);

/// Convert the face planes to vertices by polar reciprocation in sphere
/**\param dual output geometry containing polar reciprocal.
 * \param geom input geometry containing faces.
 * \param recip_rad radius of reciprocation sphere.
 * \param centre centre of reciprocation sphere.
 * \param inf maximum distance a vertex will be placed. */
void get_dual(Geometry &dual, const Geometry &geom, double recip_rad = 0,
              Vec3d centre = Vec3d(0, 0, 0), double inf = 1e20);

/// Add extra elements when an element joined to an ideal point.
/**An ideal point is a which is at "infinity" in, arbitrarily, either
 * direction. To respect this symmetry, the points can be repeated
 * and also the elements they are joined to. Should be considered a
 * visual aid only.
 * \param geom input geometry, generally containing a dual.
 * \param centre centre of original reciprocation sphere.
 * \param inf distance at which a point is considered to be ideal. */
void add_extra_ideal_elems(Geometry &geom, Vec3d centre, double inf);

/// Make a geometry with a face for each edge in the original
/**Like the Conway 'join' operation.
 * \param geom geometry with edges to convert into faces.*/
void make_edges_to_faces(Geometry &geom);

/// Truncate specified vertices
/**\param geom geometry to truncate
 * \param v_idxs index numbers of vertices to truncate.
 * \param ratio truncation points divide edges in this ratio.
 * \param info use if passed, otherwise a local one will be used. */
void truncate_verts(Geometry &geom, std::vector<int> &v_idxs, double ratio,
                    GeometryInfo *info = nullptr);

/// Truncate vertices
/**\param geom geometry to truncate
 * \param ratio truncation points divide edges in this ratio.
 * \param order if positive, truncate only vertices with this vertex order.
 * \param info use if passed, otherwise a local one will be used. */
void truncate_verts(Geometry &geom, double ratio, int order = 0,
                    GeometryInfo *info = nullptr);

/// Merge coincident elements
/**\param geom geometry with elements to be merged
 * \param merge_elems a string containg letters v, e, and f to indicate
 *  the element types to be merged. An empty string will sort the elements.
 * \param equiv_elems vector (0:vertices, 1:edges, 2:faces) of maps from an
 *  original element index to a set of the original element indexes that
 *  were considered to be equivalent.
 * \param eps a small number, coordinates differing by less than eps are
 *  the same. */
void merge_coincident_elements(
    Geometry &geom, const std::string &merge_elems,
    std::vector<std::map<int, std::set<int>>> *equiv_elems,
    double eps = epsilon);

/// Merge coincident elements
/**\param geom geometry with elements to be merged
 * \param merge_elems a string containg letters v, e, and f to indicate
 *  the element types to be merged. An empty string will sort the elements.
 * \param blend_type how to blend the colours of merged elements - 1:first color
    2:last color, 3: RGB average, 4:RYB mode.
 * \param eps a small number, coordinates differing by less than eps are
 *  the same. */
void merge_coincident_elements(Geometry &geom, const std::string &merge_elems,
                               const int blend_type, double eps = epsilon);

/// Merge coincident elements
/**\param geom geometry with elements to be merged
 * \param merge_elems a string containg letters v, e, and f to indicate
 *  the element types to be merged. An empty string will sort the elements.
 * \param eps a small number, coordinates differing by less than eps are
 *  the same. */
void merge_coincident_elements(Geometry &geom, const std::string &merge_elems,
                               double eps = epsilon);

/// Check congruence of coincident polyhedra
/**\param geom1 first geometry
 * \param geom2 second geometry
 * \param equiv_elems vector (0:vertices, 1:edges, 2:faces) of maps from an
 *  original element index to a set of the original element indexes that
 *  were considered to be equivalent (after geom2 is appended to geom1).
 * \param eps a small number, coordinates differing by less than eps are
 *  the same.
 * \return \c true if the polyhedra are congruent, otherwise \c false. */
bool check_congruence(
    const Geometry &geom1, const Geometry &geom2,
    std::vector<std::map<int, std::set<int>>> *equiv_elems = nullptr,
    double eps = epsilon);

/// Get elements which are equivalent by a symmetry transformation
/**The input geometry cannot have coincident elements
 * \param geom geometry.
 * \param trans transformation that carries the geometry onto itself.
 * \param elem_maps vector (0:vertices, 1:edges, 2:faces) of vectors
 *  mapping an original element index to the index it is coincident with
 *  after the transformation.
 * \param eps a small number, coordinates differing by less than eps are
 *  the same. */
void get_congruence_maps(const Geometry &geom, Trans3d trans,
                         std::vector<std::vector<int>> &elem_maps,
                         double eps = epsilon);

/// return the unit normal of all perimeter triangles
/**\param geom geometry.
 * \param face contains the vertex index numbers in the face. */
Vec3d face_norm_nonplanar_triangles(const Geometry &geom,
                                    const std::vector<int> &face);

/// return the unit normal of all perimeter triangles
/**\param geom geometry.
 * \param f_idx the face index number. */
Vec3d face_norm_nonplanar_triangles(const Geometry &geom, const int f_idx);

/// return the unit normal of all quads in polygon
/**\param geom geometry.
 * \param face contains the vertex index numbers in the face. */
Vec3d face_norm_nonplanar_quads(const Geometry &geom,
                                const std::vector<int> &face);

/// return the unit normal of all quads in polygon
/**\param geom geometry.
 * \param f_idx the face index number.
 * \return unit normal */
Vec3d face_norm_nonplanar_quads(const Geometry &geom, const int f_idx);

/// select normal by type. Newell, triangles, or quads
/**\param geom geometry.
 * \param face contains the vertex index numbers in the face.
 * \param normal_type: n - Newell, t - triangular, q - quads */
Vec3d face_normal_by_type(const Geometry &geom, const std::vector<int> &face,
                          const char normal_type);

/// select normal by type. Newell, triangles, or quads
/**\param geom geometry.
 * \param f_idx the face index number.
 * \param normal_type: n - Newell, t - triangular, q - quads */
Vec3d face_normal_by_type(const Geometry &geom, const int f_idx,
                          const char normal_type);

/// return true if maximum vertex radius is radius_range_percent (0.0 to ...)
/**greater than minimum vertex radius (visible for canonical.cc)
 * \param geom geometry to measure.
 * \param radius_range_percent limit to maximum radius over minimum radius */
bool canonical_radius_range_test(const Geometry &geom,
                                 const double radius_range_percent);

/// Canonicalize (George Hart "Mathematica" algorithm)
/**See http://library.wolfram.com/infocenter/Articles/2012/
 * \param geom geometry to canonicalise.
 * \param edge_factor small number to scale edge adjustments.
 * \param plane_factor small number to scale plane adjustments.
 * \param num_iters maximumn number of iterations.
 * \param radius_range_percent if the model outer radius increases this
 *  much over the inner radius then it is growing too much, terminate.
 * \param rep_count report on propgress after this many iterations.
 * \param alternate_loop use alternate loop.
 * \param planar_only planarise only.
 * \param normal_type: n - Newell, t -triangles, q - quads (default n)
 * \param eps a small number, coordinates differing by less than eps are
 *  the same. */
bool canonicalize_mm(Geometry &geom, const double edge_factor,
                     const double plane_factor, const int num_iters,
                     const double radius_range_percent, const int rep_count,
                     const bool alternate_loop, const bool planar_only,
                     const char normal_type = 'n', const double eps = epsilon);

/// an abbreviated wrapper for canonicalization with mathematica
/**\param geom geometry to planarize.
 * \param num_iters maximumn number of iterations.
 * \param rep_count report on propgress after this many iterations.
 * \param eps a small number, coordinates differing by less than eps are
 *  the same. */
bool canonicalize_mm(Geometry &geom, const int num_iters,
                     const int rep_count = -1, const double eps = epsilon);

/// an abbreviated wrapper for planarize with mathematica
/**\param geom geometry to planarize.
 * \param num_iters maximumn number of iterations.
 * \param rep_count report on propgress after this many iterations.
 * \param eps a small number, coordinates differing by less than eps are
 *  the same.
 * \return \c true if success, otherwise \c false */
bool planarize_mm(Geometry &geom, const int num_iters, const int rep_count = -1,
                  const double eps = epsilon);

/// returns the edge near points centroid
/**\param geom geometry to measure
 * \param cent centre from which to calculate nearpoints on edges
 * \return the centroid of the nearpoints. */
Vec3d edge_nearpoints_centroid(Geometry &geom,
                               const Vec3d cent = Vec3d(0, 0, 0));

/// Canonicalize (George Hart "Conway Notation" algorithm)
/**See http://www.georgehart.com/virtual-polyhedra/conway_notation.html
 * \param base geometry to canonicalise.
 * \param canonical_method - 'b': base/dual, 'p': adjust vertices with
 * side effect of planarization (len2() version), 'q': adjust vertices with
 * side effect of planarization (len() version), case 'f': use face centres.
 * \param num_iters maximumn number of iterations.
 * \param radius_range_percent if the model outer radius increases this
 *  much over the inner radius then it is growing too much, terminate.
 * \param rep_count report on propgress after this many iterations.
 * \param centering passed from canonical program, when centering is not used
 * \param normal_type: n - Newell, t -triangles, q - quads (default n)
 * \param eps a small number, coordinates differing by less than eps are
 *  the same.
 * \return \c true if success, otherwise \c false */
bool canonicalize_bd(Geometry &base, const int num_iters,
                     const char canonical_method,
                     const double radius_range_percent, const int rep_count,
                     const char centering, const char normal_type = 'n',
                     const double eps = epsilon);

/// an abbreviated wrapper for canonicalization with the base/dual method
/**\param geom geometry to planarize.
 * \param num_iters maximumn number of iterations.
 * \param rep_count report on propgress after this many iterations.
 * \param eps a small number, coordinates differing by less than eps are
 *  the same. */
bool canonicalize_bd(Geometry &geom, const int num_iters,
                     const int rep_count = -1, const double eps = epsilon);

/// an abbreviated wrapper for planarize with the base/dual method
/**\param geom geometry to planarize.
 * \param num_iters maximumn number of iterations.
 * \param rep_count report on propgress after this many iterations.
 * \param eps a small number, coordinates differing by less than eps are
 *  the same. */
bool planarize_bd(Geometry &geom, const int num_iters, const int rep_count = -1,
                  const double eps = epsilon);

/// Close polyhedron (basic)
/**Each hole (open circuit of edges) is converted to a face with colour col
 * holes having a vertex with more than two open edges are not filled
 * \param geom geometry to close.
 * \param col colour for the new faces (default: no colour set).
 * \return \c true if all holes could be closed, otherwise \c false. */
bool close_poly_basic(Geometry &geom, const Color &col = Color());

/// Bond polyhedra at a face
/**\param base base geometry
 * \param brick brick geometry, an alignment will be applied to it
 * \param base_f_idx base face number
 * \param brick_f_idx brick face number
 * \param off how many vertex places to turn the bond.
 * \param merge if \c true, make a base a combined geometry with the common
 * face removed, otherwise \c false, only align brick.
 * \param flip if \c true bond brick face with reverse orientation, otherwise
 *  \c false original orientation.
 * \return \c true if operation succeeded, otherwise \c false. */
bool face_bond(Geometry &base, Geometry &brick, int base_f_idx = 0,
               int brick_f_idx = 0, int off = 0, bool merge = true,
               bool flip = false);

/// Are points in convex hull
/**\param points the points to test
 * \param hull geometry containing the convex hull
 * \param inclusion_test from ORing flags INCLUSION_IN, INCLUSION_ON
 *  and INCLUSION_OUT
 * \param eps a small number, coordinates differing by less than eps are
 *  the same.
 * \return \c true or \c false */
bool are_points_in_hull(const std::vector<Vec3d> &points, const Geometry &hull,
                        unsigned int inclusion_test, const double &eps);

/// Find the index number of a vertex with a set of coordinates
/**\param geom the geometry
 * \param coords the coordinates
 * \param eps a small number, coordinates differing by less than eps are
 *  the same.
 * \return The coincident vertex with lowest index number, otherwise -1 */
int find_vert_by_coords(const Geometry &geom, const Vec3d &coords,
                        double eps = epsilon);

/// Is an edge part of a face
/**\param face the face.
 * \param edge the edge to find.
 * \return \c true if the edge is part of the face, otherwise \c false. */
bool edge_exists_in_face(const std::vector<int> &face,
                         const std::vector<int> &edge);

/// Find faces that include a particular edge
/**\param faces the faces to searc.
 * \param edge the edge to find.
 * \return The index numbers of the faces that include the edge. */
std::vector<int>
find_faces_with_edge(const std::vector<std::vector<int>> &faces,
                     const std::vector<int> &edge);

/// Does vertex index number appear in a face
/**\param face the face.
 * \param v_idx the vertex index number to find.
 * \return \c true if the vertex is part of the face, otherwise \c false. */
bool vertex_exists_in_face(const std::vector<int> &face, int v_idx);

/// Does vertex index number appear in an edge
/**\param edge the edge.
 * \param v_idx the vertex index number to find.
 * \return \c true if the vertex is part of the edge, otherwise \c false. */
bool vertex_exists_in_edge(const std::vector<int> &edge, int v_idx);

/// Find faces that include a particular vertex
/**\param faces the faces to search.
 * \param v_idx the vertex index number to find.
 * \return The index numbers of the faces that include the vertex. */
std::vector<int>
find_faces_with_vertex(const std::vector<std::vector<int>> &faces, int v_idx);

/// Find edges that include a particular vertex
/**\param edges the edges to search.
 * \param v_idx the vertex index number to find.
 * \return The index numbers of the edges that include the vertex. */
std::vector<int>
find_edges_with_vertex(const std::vector<std::vector<int>> &edges, int v_idx);

/// Find edge in edge list
/**\param edges the edges to search.
 * \param edge the edge to find.
 * \return The corresponding edge with lowest index number, otherwise -1 */
int find_edge_in_edge_list(const std::vector<std::vector<int>> &edges,
                           const std::vector<int> &edge);

/// Find edges which do not correspond to the edge of a face
/**\param geom the geometry
 * \return the unmatched edges. */
std::vector<std::vector<int>> find_unmatched_edges(const Geometry &geom);

} // namespace anti

#endif // GEOMETRYUTILS_H
