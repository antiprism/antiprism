/*
   Copyright (c) 2010-2016, Roger Kaufman

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

/*!\file base/normal.h
 * \brief classes for working with normals
*/

#ifndef NORMAL_H
#define NORMAL_H

#include "geometry.h"

namespace anti {

// Single normal
class Normal {
private:
  Vec3d normal;
  int direction; // 1 = outward  0 = hemispherical  -1 = inward

public:
  /// Constructor
  /**Unset normal. */
  Normal();

  /// Constructor
  /**Face normal. */
  Normal(const Geometry &geom, const int face_idx, Vec3d C = Vec3d(),
         double eps = epsilon);

  /// Constructor
  /**For precalculated normals. */
  Normal(const Geometry &geom, const Vec3d &norm, const int v_idx,
         Vec3d C = Vec3d(), double eps = epsilon);

  /// Check if set
  /**\return \c true if set, othewise \c false. */
  bool is_set() const;

  /// Get the raw normal
  /**\return The raw normal. */
  Vec3d raw() const;

  /// Get the normal, converted to unit length
  /**\return The unit normal. */
  Vec3d unit() const;

  /// Get the normal, directed to point outwards
  /**\return The normal directed outwards. */
  Vec3d outward() const;

  /// Get the normal, directed to point outwards
  /**\return The normal directed inwards. */
  Vec3d inward() const;

  /// Check if normal is outward
  /**\return \c true if outward, othewise \c false. */
  bool is_outward() const;

  /// Check if normal is outward
  /**\return \c true if outward, othewise \c false. */
  bool is_inward() const;

  /// Check if normal is outward
  /**\return \c true if outward, othewise \c false. */
  bool is_hemispherical() const;
};

// Geometry face normals
class FaceNormals {
private:
  const Geometry *ngeom;
  std::vector<Normal> normals;

  Vec3d average_normals(std::vector<int> &face_idx,
                        const std::string &average_pattern);

public:
  /// Constructor
  FaceNormals();

  /// Constructor
  /**\param geom the geometry.
   * \param cent the centre of the model, if set, otherwise vertex centroid.
   * \param eps a small number, coordinates differing by less than eps are
   *  the same. */
  FaceNormals(const Geometry &geom, Vec3d cent = Vec3d(), double eps = epsilon);

  /// Refresh
  /**\param geom the geometry.
   * \param cent the centre of the model, if set, otherwise vertex centroid.
   * \param eps a small number, coordinates differing by less than eps are
   *  the same. */
  void refresh(const Geometry &geom, Vec3d cent = Vec3d(),
               double eps = epsilon);

  /// Get edge normal
  /**\param v_idx1 vertex index 1
   * \param v_idx2 vertex index 2
   * \param average_pattern type of normals to average - r:raw,
   *  o:outward, i:inward, u: unit.
   * \param cent the centre of the model, if set, otherwise vertex centroid.
   * \param eps a small number, coordinates differing by less than eps are
   *  the same. */
  Normal edge_normal(int v_idx1, int v_idx2, const std::string &average_pattern,
                     Vec3d cent = Vec3d(), double eps = epsilon);

  /// Get vertex normal
  /**\param v_idx vertex index
   * \param average_pattern type of normals to average - r:raw,
   *  o:outward, i:inward, u: unit.
   * \param cent the centre of the model, if set, otherwise vertex centroid.
   * \param eps a small number, coordinates differing by less than eps are
   *  the same. */
  Normal vert_normal(int v_idx, const std::string &average_pattern,
                     Vec3d cent = Vec3d(), double eps = epsilon);

  /// Size of normal vector
  /**\return Size of normal vector. */
  unsigned int size() const;

  /// Is normal index valid
  /**\param idx normal index.
   * \return \c true if valid, otherwise \c false invalid. */
  bool in_range(unsigned int idx) const;

  /// Check if set
  /**\return \c true if set, othewise \c false. */
  bool is_set() const;

  /// Access to normal vector
  /**\param idx index.
   * \return The normal at that index. */
  Normal operator[](unsigned int idx) const;
};

// inline functions

// -------------------------------------------------------------------
// Normal::

inline Normal::Normal() : direction(0) {}

inline bool Normal::is_set() const { return normal.is_set(); }

inline Vec3d Normal::raw() const { return normal; }

inline Vec3d Normal::unit() const
{
  return normal.unit();
} // unit() of raw normal

inline Vec3d Normal::outward() const
{
  return (is_inward() ? -normal : normal);
}

inline Vec3d Normal::inward() const
{
  return (is_outward() ? -normal : normal);
}

inline bool Normal::is_outward() const { return (direction == 1); }

inline bool Normal::is_inward() const { return (direction == -1); }

inline bool Normal::is_hemispherical() const { return (direction == 0); }

// -------------------------------------------------------------------
// FaceNormals::

inline FaceNormals::FaceNormals() {}

inline FaceNormals::FaceNormals(const Geometry &geom, Vec3d cent, double eps)
{
  refresh(geom, cent, eps);
}

inline unsigned int FaceNormals::size() const { return (normals.size()); }

inline bool FaceNormals::in_range(unsigned int idx) const
{
  return (idx < (unsigned int)size());
}

inline bool FaceNormals::is_set() const { return (size() > 0); }

inline Normal FaceNormals::operator[](unsigned int idx) const
{
  return ((is_set() && in_range(idx)) ? normals[idx] : Normal());
}

} // namespace anti

#endif // NORMAL_H
