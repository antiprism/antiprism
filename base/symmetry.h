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

/*!\file symmetry.h
 *\brief symmetry and transformation handling
 */

#ifndef SYMMETRY_H
#define SYMMETRY_H

#include "geometry.h"
#include "trans3d.h"

#include <cmath>
#include <set>
#include <string>
#include <vector>

namespace anti {

/// A set of transformations
class Transformations {
private:
  struct SymTransLess {
    bool operator()(const Trans3d &t1, const Trans3d &t2) const
    {
      return compare(t1, t2, sym_eps) < 0;
    }
  };
  std::set<Trans3d, SymTransLess> trans;

public:
  /// A constant iterator over the set
  typedef std::set<Trans3d, SymTransLess> SymTransSet;
  typedef SymTransSet::const_iterator const_iterator;

  /// Add a transformation
  /**\param tr transformation to add
   * \return A reference to this object with the transformation added.*/
  Transformations &add(const Trans3d &tr)
  {
    trans.insert(tr);
    return *this;
  }

  /// Clear all transformations
  /** Remove all the transformations
   * \return A reference to this object with the transformationd cleared.*/
  Transformations &clear()
  {
    trans.clear();
    return *this;
  }

  /// Set to direct product of two sets of transformations
  /** Each transformation of one set is combined with every
   *  transformation from the other set.
   * \param s1 one of the sets.
   * \param s2 the other set.
   * \return A reference to this object with the direct product
   *  of transformations set. */
  Transformations &product(const Transformations &s1,
                           const Transformations &s2);

  /// Set to direct product of current transformations with another set of
  /// transformations
  /** Each transformation of current set is combined with every
   *  transformation from the other set.
   * \param s the other set.
   * \return A reference to this object with the direct product
   *  of transformations set. */
  Transformations &product_with(const Transformations &s);

  /// Set transformations to their conjugates by a given transformation.
  /** If the given transformation is M then each transformation
   *  T is set to MTM^-1.
   * \param t the transformation to form the conjugates with.
   * \return A reference to this object with the conjugates set. */
  Transformations &conjugate(const Trans3d &t);

  /// Intersection of two sets of transformations.
  /** Find the transformations that are common to both sets of
   *  transformations.
   * \param s1 one of the sets.
   * \param s2 the other set.
   * \return A reference to this object with common transformations set. */
  Transformations &intersection(const Transformations &s1,
                                const Transformations &s2);

  /// Subtract from current transformations any common to a second set of
  /// transformations.
  /** After this operation the two sets will have no transformations
   *  in common.
   * \param s the second set.
   * \return A reference to this object with the common transformations
   *  removed. */
  Transformations &subtract(const Transformations &s);

  /// Form the left coset with a given transformation.
  /** If the given transformation is M then each transformation
   *  T is set to MT in the coset.
   * \param t the transformation to form the left coset with.
   * \return The coset with \arg m. */
  Transformations lcoset(const Trans3d &t) const;

  /// Form the left cosets for a given subgroup.
  /**\param sub a subgroup of this transformation group to form the
   * left cosets with.
   * \param lcosets the left cosets.
   * \return The number of cosets.*/
  int lcosets(const Transformations sub,
              std::vector<Transformations> &lcosets) const;

  /// The minimum set of transformations that will generate one group of
  /// transformations from another
  /**\param tr_whole the target group of transformations.
   * \param tr_part the group of transformations that will be used to
   *  generate the target group.
   * \param pos A transformation that will align \arg tr_part with
   *  \arg tr_whole.
   * \return A reference to this object containing a minimum set of
   *  transformations that when applied to a geometry with symmetry
   *  \arg tr_part will create a compound with symmetry \arg tr_whole. */
  Transformations &min_set(const Transformations &tr_whole,
                           const Transformations &tr_part,
                           const Trans3d &pos = Trans3d());

  /// Begin, for iterating through the set.
  /**\return A constant iterator pointing to the first transformation. */
  const_iterator begin() const { return trans.begin(); }

  /// End, for iterating through the set.
  /**\return A constant iterator pointing to one past the last
   *  transformation. */
  const_iterator end() const { return trans.end(); }

  /// Size, the number of transformations.
  /**\return The number of transformations. */
  size_t size() const { return trans.size(); }

  /// Get the underlying transformation set
  /**\return The underlying transformation set */
  const SymTransSet &get_trans() const { return trans; }

  /// Get the underlying transformation set
  /**\return The underlying transformation set */
  SymTransSet &get_trans() { return trans; }

  /// Check if there are transformations in the set.
  /**\return \c true if set, otherwise \c false */
  bool is_set() const { return trans.size() > 0; }
};

/// Form direct product of two sets of transformations
/** Each transformation from one set is combined with every transformation
 *  from the other set.
 * \param s1 one of the sets.
 * \param s2 the other set.
 * \return The direct productof the two sets. */
Transformations operator*(const Transformations &s1, const Transformations &s2);

/// Form direct product of two sets of transformations
/** Each transformation of the first set is combined with every
 *  transformation from the second set and the result added to the
 *  first set.
 * \param s1 the first set.
 * \param s2 the second set.
 * \return A reference to s1 with the direct product of transformations set. */
Transformations &operator*=(Transformations &s1, const Transformations &s2);

/// Add a transformation
/**\param s a set of transformations.
 * \param m the transformation to add.
 * \return A copy of \arg s with \arg m added. */
Transformations operator+(const Transformations &s, const Trans3d &m);

/// Add a transformation
/**\param s a set of transformations.
 * \param m the transformation to add.
 * \return A refernce to \arg s with \arg m added. */
Transformations &operator+=(Transformations &s, const Trans3d &m);

/// Compare two Transformationss for order
/** Returns -1, 0, or 1 to indicate less, equal or greater. Order by
 *  size, and if equal compare matrices sequentially until not equal,
 *  and if all equal return 0.
 * \param t0 first Transformations
 * \param t1 second Transformations
 * \return  -1, 0, or 1 to indicate t0 is less, equal or greater than t1 */

int compare(const Transformations &t0, const Transformations &t1);

/// Isometry type
class Isometry {
public:
  /// Rotation type
  enum { rt_none = 0, rt_unit, rt_rot, rt_inv, rt_refl, rt_rot_refl };

private:
  int rot_type;
  Vec3d axis;
  double ang;
  Vec3d transl;

  bool is_isometry(Trans3d m) const;
  void normalise();

public:
  /// Constructor
  Isometry() : rot_type(rt_none) {}

  /// Constructor
  /** Set up with the details for a particular transformation matrix.
   * \param m the transformation matrix. */
  Isometry(const Trans3d &m) { init(m); }

  /// Initialise with a transformation matrix.
  /** Set up with the details for a particular transformation matrix.
   * \param m the transformation matrix.
   * \return A reference to this symmetry axis. */
  Isometry &init(Trans3d m);

  /// Get the rotation type
  /** A return value of \c rt_none indicates that the object has not
   *  been initialised with an isometry.
   * \return The rotation type. */
  int get_rot_type() const { return rot_type; }

  /// Get principal axis.
  /** The axis will be unset if the rotation type doesn't have an axis.
   * \return The axis. */
  Vec3d get_axis() const { return axis; }

  /// Get rotation angle.
  /** The angle will be 0 if the rotation type doesn't have an angle.
   * \return The angle. */
  double get_ang() const { return ang; }

  /// Get translation.
  /** The translation will be unset if the object has not been
   *  initialised with an isometry.
   * \return The translation. */
  Vec3d get_transl() const { return transl; }

  /// Check if isometry is direct.
  /**\return \c true if the isometry is direct, otherwise \c false */
  bool is_direct() const { return rot_type > rt_none && rot_type < rt_inv; }

  /// Dump
  /** Print the object data to \c stdout for debugging. */
  void dump() const;
};

/// Compare two Transformationss for equality
/**\param t0 first Transformations
 * \param t1 second Transformations
 * \return  \c true id \c t0 and \c t1 contain the same set of transformations,
 * otherwise \c false.*/
inline bool operator==(const Transformations &t0, const Transformations &t1)
{
  return compare(t0, t1) == 0;
}

/// Get sets of elements that are equivalent under a set of transformations
/**\param geom the geometry.
 * \param ts the set of transfromations to apply.
 * \param equiv_sets vectors of sets of equivalent elements for
 *  vertices (0), edges (1) and faces (2). */
void get_equiv_elems(const Geometry &geom, const Transformations &ts,
                     std::vector<std::vector<std::set<int>>> *equiv_sets);
/// Subspace
class Subspace {
public:
  /// Type of subspace, corresponding to dimension
  enum class SubspaceType {
    none,  ///< Uninitialized
    point, ///< A point
    line,  ///< A line
    plane, ///< A plane
    space  ///< The whole space
  };
  /// Constructor
  /**\param type type of subspace
   * \param point a point fixed by the subspace (or uninitialized)
   * \param direction direction of axis or normal to plane (or uninitialized)*/
  Subspace(SubspaceType type = SubspaceType::none, Vec3d point = Vec3d(),
           Vec3d direction = Vec3d());

  /// Nearest point in the subspace
  /**\param P the point
   * \return nearest point to P in the subspace. */
  Vec3d nearest_point(Vec3d P) const;

  /// A string representation of the subspace
  std::string str() const;

private:
  SubspaceType type;
  Vec3d point;
  Vec3d direction;
};

class SymmetryAxis;
class Symmetry;

/// Automorphisms: transformations maintaining symmetry alignment
class SymmetryAutos {
private:
  std::vector<Trans3d> fixed_trans;
  unsigned int free_vars;
  int fixed_type;
  double rot[3];
  double transl[3];

  void init();

public:
  enum {
    FREE_NONE = 0,
    FREE_ROT_PRINCIPAL = 1,
    FREE_ROT_FULL = 2,
    FREE_TRANSL_PRINCIPAL = 4,
    FREE_TRANSL_PLANE = 8,
    FREE_TRANSL_SPACE = 16
  };

  /// Constructer
  SymmetryAutos() { init(); }

  /// Constructer
  /**\param sym symmetry group to transform */
  SymmetryAutos(const Symmetry &sym);

  /// Check if has been set
  bool is_set() const { return fixed_trans.size(); }

  /// Set fixed transformation type
  /**\param type the index number of the type (0 always the identity)
   * \return status, evaluates to true if the fixed transformation
   *  index is valid, otherwise false.*/
  Status set_fixed_type(int type);

  /// Set rotation about principal axis
  /**\param rot_ang rotation angle about prinicipal axis, in degrees
   * \return status, evaluates to \c true if the free variable could
   *  be set, otherwise \c false.*/
  Status set_rot_principal(double rot_ang);

  /// Set rotation about origin
  /**\param rot_x rotation about x-axis, in degrees
   * \param rot_y rotation about y-axis, in degrees
   * \param rot_z rotation about z-axis, in degrees
   * \return status, evaluates to \c true if the free variable could
   *  be set, otherwise \c false.*/
  Status set_rot_full(double rot_x, double rot_y, double rot_z);

  /// Set translation distance along principal direction
  /**\param transl0 translation distance
   * \return status, evaluates to \c true if the free variable could
   *  be set, otherwise \c false.*/
  Status set_transl_principal(double transl0);

  /// Set translation distances along directions which span a plane
  /**\param transl0 first translation distance
   * \param transl1 second translation distance
   * \return status, evaluates to \c true if the free variable could
   *  be set, otherwise \c false.*/
  Status set_transl_plane(double transl0, double transl1);

  /// Set translation distances along directions which span all space
  /**\param transl0 first translation distance
   * \param transl1 second translation distance
   * \param transl2 third translation distance
   * \return status, evaluates to \c true if the free variable could
   *  be set, otherwise \c false.*/
  Status set_transl_space(double transl0, double transl1, double transl2);

  /// Get  fixed realignment transformations
  /**\return The fixed realigment transformations.*/
  const std::vector<Trans3d> &get_fixed() const { return fixed_trans; }

  /// Set  fixed realignment transformations
  /**\param fixed the fixed realigment transformations.*/
  void set_fixed(const Transformations &fixed);

  /// Number of free rotation variables
  /**\return Number of free rotation variables (0, 1 or 3).*/
  int num_free_rots() const;

  /// Number of free translation variables
  /**\return Number of free translation variables (0 to 3).*/
  int num_free_transls() const;

  /// Get realignment transformtion from current settings
  /**\return Realignment transformation */
  Trans3d get_realignment() const;

  /// Set and get realignment transformation from a string of numbers
  /**\param realign fixed type number and free variable list as
   *  colon separated string
   * \return status, evaluates to \c true if the realignment was valid,
   *  otherwise false.*/
  Status set_realignment(const std::string &realign);
};

/// Symmetry group, using Schoenflies notation.
class Symmetry {
public:
  /// Schoenflies symmetry type identifiers
  enum {
    unknown,
    C1,
    Ci,
    Cs,
    C,
    Cv,
    Ch,
    D,
    Dv,
    Dh,
    S,
    T,
    Td,
    Th,
    O,
    Oh,
    I,
    Ih
  };

private:
  int sym_type;
  int nfold;
  mutable std::set<Symmetry> sub_syms;
  mutable std::set<SymmetryAxis> axes;
  mutable std::set<Vec3d> mirrors;
  SymmetryAutos autos;
  Trans3d to_std;

  void add_sub_axes(const Symmetry &sub) const;
  void find_full_sym_type(const std::set<SymmetryAxis> &full_sym);

  /// Set the symmetry type for the axis
  /**\param type the symmetry type of the axis as the Schoenflies
   *  identifier. */
  void set_sym_type(int type) { sym_type = type; }

  /// Set the n-fold order of the axis.
  /** If the symmetry type is \c sym_S then the axis has
   *  rotational n/2-fold symmetry.
   * \param n the n-fold order of the axis. */
  void set_nfold(int n) { nfold = n; }

public:
  /// Constructor
  /** Find the symmetry in Schoenflies notation.
   * \param geom geometry to find the symmetry type for.
   * \param equiv_sets for the vertices, edges and faces set up a
   *  vector of sets of equivalent elements. */
  Symmetry(const Geometry &geom,
           std::vector<std::vector<std::set<int>>> *equiv_sets = nullptr);

  /// Constructor
  /** Find the symmetry in Schoenflies notation.
   * \param ts group of transformations to find the symmetry type for.*/
  Symmetry(const Transformations &ts);

  /// Constructor
  /**\param type the symmetry type.
   * \param n n-fold order, for a principal axis, otherwise 0.
   * \param pos transformation that carries an object of this
   *  symmetry type onto the standard symetry type.
   * \param stat to return any error message. */
  Symmetry(int type = 0, int n = 0, const Trans3d &pos = Trans3d(),
           Status *stat = nullptr);

  /// Constructor
  /**\param name name of symmetry type to set up.
   * \param pos transformation that carries an object of this
   *  symmetry type onto the standard symetry type.
   * \param stat to return any error message. */
  Symmetry(const std::string &name, const Trans3d &pos = Trans3d(),
           Status *stat = nullptr);

  /// Constructor
  /**\param sym_axis symmetry type in symmetry axis form.
   * \param cent coordinates of symmetry centre */
  Symmetry(const SymmetryAxis &sym_axis, const Vec3d &cent = Vec3d::zero);

  /// Initialiser
  /** Find the symmetry in Schoenflies notation.
   * \param geom geometry to find the symmetry type for.
   * \param equiv_sets for the vertices, edges and faces set up a
   *  vector of sets of equivalent elements
   * \return status, evaluates to \c true if the symmetry type
   *  could be determined, otherwise \c false.*/
  Status init(const Geometry &geom,
              std::vector<std::vector<std::set<int>>> *equiv_sets = nullptr);

  /// Initialiser
  /**\param type the symmetry type.
   * \param n n-fold order, for a principal axis, otherwise 0.
   * \param pos transformation that carries an object of this
   *  symmetry type onto the standard symetry type.
   * \return status, evaluates to \c true if the symmetry specification
   *  was valid, otherwise false.*/
  Status init(int type, int n = 0, const Trans3d &pos = Trans3d());

  /// Initialiser
  /**\param name name of symmetry type to set up.
   * \param pos transformation that carries an object of this
   *  symmetry type onto the standard symetry type.
   * \return true if the symmetry type name was valid, otherwise false.*/
  Status init(const std::string &name, const Trans3d &pos = Trans3d());

  /// Initialiser
  /**\param sym_axis symmetry type in symmetry axis form.
   * \param cent coordinates of symmetry centre */
  void init(const SymmetryAxis &sym_axis, const Vec3d &cent = Vec3d::zero);

  /// Get the symmetry type
  /**\return the symmetry type of the axis as the Schoenflies
   *  identifier. */
  int get_sym_type() const { return sym_type; }

  /// Get the n-fold order of the axis with highest n.
  /** This value is used in conjunction with \c get_sym_type()
   *  to describe principal axis symmetry types. If the symmetry
   *  type is \c sym_S and the axis has rotational n-fold symmetry
   *  then 2n is returned.
   * \return the n-fold order of the axis. */
  int get_nfold() const { return nfold; }

  /// Get the Schoenflies symbol for the symmetry type
  /**\return The symbol. */
  std::string get_symbol() const;

  /// Get the tranformation to standard symmetry type.
  /**\return The transformation that carries an object with the symmetry
   *  type onto the standard set of symmetries for that type. */
  Trans3d get_to_std() const { return to_std; }

  /// Set the tranformation to standard symmetry type.
  /**\param trans the transformation that carries an object with the
   *  symmetry type onto the standard set of symmetries for that type. */
  void set_to_std(const Trans3d &trans) { to_std = trans; }

  /// Get the symmetry transformations
  /**\param ts to return the set of symmetry transformations for this
   *  symmetry type.
   * \return A reference to the symmetry transformations.*/
  Transformations &get_trans(Transformations &ts) const;

  /// Get the symmetry transformations
  /**\return The set of symmetry transformations for this symmetry type. */
  Transformations get_trans() const;

  /// Get the symmetry subgroups
  /** Only one example is included from each conjugacy class
   * \return The symmetry subgroups. */
  const std::set<Symmetry> &get_sub_syms() const;

  /// Get the Euclidean outer automorphisms
  /** These map the symmetry transformations onto themselves, with one
   *  example included from each conjugacy class, the first value
   *  in the returned transformations is always the identity.
   * \return The symmetry automorphisms. */
  SymmetryAutos &get_autos();

  /// Get a symmetry subgroup
  /**\param sub_sym the symmetry subgroup
   * \param sub to return the symmetry subgroup
   * \param conj_type use to select from inequivalent (non-conjugate)
   *  subgroups.
   * \return status, evaluates to \c true if the symmetry subgroup
   *  could be found, otherwise \c false.*/
  Status get_sub_sym(const Symmetry &sub_sym, Symmetry *sub,
                     int conj_type = 0) const;

  /// Get a symmetry subgroup
  /**\param sub_name the symmetry subgroup name and conjugation number,
   *  separated by a comma
   * \param sub to return the symmetry subgroup
   * \return status, evaluates to \c true if the symmetry subgroup
   *  could be found, otherwise \c false.*/
  Status get_sub_sym(const std::string &sub_name, Symmetry *sub) const;

  /// Get a maximal direct symmetry subgroup
  /**\return the symmetry group containing all direct transformations.**/
  Symmetry get_max_direct_sub_sym() const;

  /// Get the subspace fixed by the symmetry group
  /**\return The subspace fixed by the symmetry group. */
  Subspace get_fixed_subspace() const;

  /// Get the axes or mirrors.
  /**\return The symmetry axes or mirrors. */
  const std::set<SymmetryAxis> &get_axes() const;

  /// Get the mirrors.
  /**\return The mirror normals. */
  const std::set<Vec3d> &get_mirrors() const;

  /// Is an inversion a symmetry.
  /**\return \c true if an inversion is a symmetry, otherwise \c false */
  bool has_inversion_symmetry() const;

  /// Check if a valid symmetry type is set.
  /**\return \c true if set to valid type, otherwise \c false */
  bool is_set() const { return sym_type != unknown; }

  /// Less than.
  /** Only for containers that need it
   * \param s the symmetry axis to compare for less than.
   * \return \c true if this axis is less than \arg s, otherwise \c false.*/
  bool operator<(const Symmetry &s) const;
};

/// Principal axis or mirror symmetry, using Schoenflies notation.
class SymmetryAxis {
private:
  int sym_type;
  Vec3d axis;
  Vec3d perp;
  long nfold;

public:
  /// Constructor
  /** Set up with the smallest symmetry type that includes a
   *  particular transfomation.
   * \param m the transformation. */
  SymmetryAxis(const Trans3d &m = Trans3d());

  /// Less than.
  /** Only for containers that need it
   * \param s the symmetry axis to compare for less than.
   * \return \c true if this axis is less than \arg s, otherwise \c false.*/
  bool operator<(const SymmetryAxis &s) const;

  /// Set the principal axis
  /**\param ax the axis. */
  void set_axis(Vec3d ax) { axis = ax; }

  /// Get the principal axis
  /**\return the axis. */
  const Vec3d &get_axis() const { return axis; }

  /// Set the perpendicular axis
  /**\param perp_ax the axis. */
  void set_perp(Vec3d perp_ax) { perp = perp_ax; }

  /// Get the perpendicular axis
  /** Not all symmetry types have a perpendicular axis.
   * \return The perpendicular axis. */
  const Vec3d &get_perp() const { return perp; }

  /// Set the symmetry type for the axis
  /**\param type the symmetry type of the axis as the Schoenflies
   *  identifier. */
  void set_sym_type(int type) { sym_type = type; }

  /// Get the symmetry type for the axis
  /**\return the symmetry type of the axis as the Schoenflies
   *  identifier. */
  int get_sym_type() const { return sym_type; }

  /// Check for horizontal mirror
  /**\return \c true if there is a horizontal mirror,
   *  otherwise \c false. */
  bool has_horz_mirror() const
  {
    return (sym_type == Symmetry::Cs || sym_type == Symmetry::Ch ||
            sym_type == Symmetry::Dh);
  }

  /// Set the n-fold order of the axis.
  /** If the symmetry type is \c sym_S then the axis has
   *  rotational n/2-fold symmetry.
   * \param n the n-fold order of the axis. */
  void set_nfold(int n) { nfold = n; }

  /// Get the n-fold order of the axis.
  /** If the symmetry type is \c sym_S and the axis has
   *  rotational n-fold symmetry then 2n is returned.
   * \return the n-fold order of the axis. */
  int get_nfold() const { return nfold; }

  /// Dump
  /** Print the object data to \c stdout for debugging. */
  void dump() const;
};

/// Orbit relationship mapping
class ElemOrbitMapping {
public:
  /// Constructor
  /**\param orbit_no the orbit number that the element lies on
   * \param trans_to transformation that carries a particular element onto the
   *                 principal orbit element
   * \param trans_from transformation that carries the principal orbit element
   *                   onto a particular element*/
  ElemOrbitMapping(int orbit_no, std::set<Trans3d>::const_iterator trans_to,
                   std::set<Trans3d>::const_iterator trans_from)
      : orbit_no(orbit_no), trans_to(trans_to), trans_from(trans_from)
  {
  }
  /// Get the orbit number
  /**\return the orbit number*/
  int get_orbit_no() const { return orbit_no; }

  /// Get the transformation onto the principal orbit element
  /**\return the transformation that carries a particular element onto the
   *         principal orbit element*/
  const Trans3d &get_trans_to() const { return *trans_to; }

  /// Get the transformation from the principal orbit element
  /**\return the transformation that carries the principal orbit element
   *         onto a particular element*/
  const Trans3d &get_trans_from() const { return *trans_from; }

private:
  int orbit_no;
  std::set<Trans3d>::const_iterator trans_to;
  std::set<Trans3d>::const_iterator trans_from;
};

/// Symmetric propogation of movements of individual vertices
class SymmetricUpdater {
public:
  /// Constructor
  /**\param geom the geometry
   * \param sym symmetry group, or subgroup, of the geometry */
  SymmetricUpdater(const Geometry &geom, Symmetry sym);

  /// Get equivalent sets for an element type
  /**\param elem_type VERTS, EDGES or FACES
   * \return vector of sets of equivalent elements */
  const std::vector<std::set<int>> get_equiv_sets(int elem_type) const
  {
    return equiv_sets[elem_type];
  }

  /// Get principal (first) elements of an element type
  /**\param elem_type VERTS, EDGES or FACES
   * \return index numbers of principal elements of a type. */
  std::vector<int> get_principal(int elem_type);

  /// Remove duplicates from a list of index numbers
  /**\param idxs the list of index numbers to process */
  static void to_unique_index_list(std::vector<int> &idxs);

  /// Make a sequential list of all index numbers
  /**\param size the number of index numbers
   * \return the index numbers. */
  static std::vector<int> sequential_index_list(int size);

  /// Get a list of all included vertices
  /**\param elem_idxs the key elements that a vertex is part of
   * \param elem_verts the vertices that are part of each element
   * \return all associated vertices. */
  static std::vector<int>
  get_included_verts(const std::vector<int> &elem_idxs,
                     const std::vector<std::vector<int>> &elem_verts);

  /// Make a list of associated element index numbers for principal vertices
  /**\param elems element lists in principal vertex order
   * \return all included element index numbers. */
  std::vector<int>
  get_associated_elems(const std::vector<std::vector<int>> &elems);

  /// Update the principal orbit vertex location using any vertex on the
  /// orbit
  /**\param v_idx the index of the vertex whose location will be updated
   * \param point the new location of the vertex */
  void update_principal_vertex(int v_idx, Vec3d point);

  /// Update a vertex location using the principal orbit vertex
  /**\param v_idx the index of the vertex which will be updated
   * \return the updated vertex coordinates */
  Vec3d update_from_principal_vertex(int v_idx);

  /// Update all vertex locations using principal orbit vertices
  void update_all();

  /// Prepare for the next iteration
  /** If deferred, switch current geometry and updated geometry. */
  void prepare_for_next_iteration();

  /// Get the working geometry (for updating vertices of interest)
  /** Use \c get_geom_final() to get a fully updated geometry
   * \return the working geometry */
  const Geometry &get_geom_working() const { return geom; }

  /// Get the final geometry after all iteration has finished
  /** \return the final geometry */
  const Geometry &get_geom_final();

private:
  Geometry geom;
  Symmetry symmetry;
  Transformations transformations;
  std::vector<std::vector<std::set<int>>> equiv_sets;
  std::vector<int> orbit_vertex_idx;
  std::vector<Subspace> orbit_invariant_subspaces;
  std::vector<ElemOrbitMapping> orbit_mapping;

  void init_vert_orbit(int orbit_idx, const std::set<int> &orbit);
};

/// Get element equivalence transformations
/** Partition the elements into sets of equivalents. The stabilizer of
 * an element is a subgroup of a symmetry group consisting of all the
 * transformations that map the element onto itself.
 * \param geom the geometry
 * \param idx the element index number
 * \param elem_type type of element: ELEM_VERTS, ELEM_EDGES, ELEM_FACES
 * \param sym the symmetry group
 * \return The stabilizer group.*/
Symmetry get_elem_stabilizer(const Geometry &geom, int idx, int elem_type,
                             const Symmetry &sym);

/// Check if a set of transformations fixes an element
/**An element is fixed if all the transformations map the element onto itself.
 * \param geom the geometry
 * \param idx the element index number
 * \param elem_type type of element: ELEM_VERTS, ELEM_EDGES, ELEM_FACES
 * \param sym the symmetry group
 * \return \c true if the element was fixed, otherwise \c false.*/
bool elem_fixed(const Geometry &geom, int idx, int elem_type,
                const Symmetry &sym);

} // namespace anti

#endif // SYMMETRY_H
