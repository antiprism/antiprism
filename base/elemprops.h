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

/*!\file elemprops.h
   \brief A class to represent colours in geometries
*/

#ifndef ELEMPROPS_H
#define ELEMPROPS_H

#include "color.h"

#include <map>

namespace anti {

template <class T> class ElemProps {
private:
  // Element index to element propert mapping
  std::map<int, T> ElemProps;

public:
  /// Set an element property.
  /**\param idx the element index number.
   * \param prop the property to set. */
  void set(int idx, const T &prop);

  /// Delete property.
  /**\param idx the element index number. */
  void del(int idx);

  /// Get an element property.
  /**\param idx the element index number.
   * \return The property. */
  T get(int idx) const;

  /// Clear all element properties.
  void clear();

  /// Get the properties map
  /**\return The properties map. */
  const std::map<int, T> &get_properties() const;

  /// Get the properties map
  /**\return The properties map. */
  std::map<int, T> &get_properties();

  /// Map properties to different index numbers.
  /**Used to maintain properties when index numbers are changed. This
   * can happen after deletions.
   * \param chg_map a map of old index numbers to new index numbers.
   *           if the new index number is \c -1 then the element index
   *           has been deleted so the property is deleted. */
  void remap(const std::map<int, int> &chg_map);
};

/// Geometry property container
template <class T> class GeomElemProps {
private:
  /// Colours for vertices, edges and faces
  ElemProps<T> elem_props[3];

public:
  /// Get the element properties for an element type
  /**\param type from VERTS, EDGES, FACES.
   * \return The element properties. */
  const ElemProps<T> &operator[](int type) const;

  /// Get the element properties for an element type
  /**\param type from VERTS, EDGES, FACES.
   * \return The element properties. */
  ElemProps<T> &operator[](int type);

  /// Clear all properties for all elements
  void clear();

  /// Append a geometry property_container.
  /**\param geom_props geometry property holder to append.
   * \param v_size number of vertices in geometry associated with geom.
   * \param e_size number of edges in geometry associated with geom.
   * \param f_size number of faces in geometry associated with geomi. */
  void append(const GeomElemProps &geom_props, int v_size, int e_size,
              int f_size);
};

// Implementation

template <class T> void ElemProps<T>::set(int idx, const T &prop)
{
  if (prop.is_set())
    ElemProps[idx] = prop;
  else
    del(idx);
}

template <class T> void ElemProps<T>::del(int idx) { ElemProps.erase(idx); }

template <class T> T ElemProps<T>::get(int idx) const
{
  auto mi = ElemProps.find(idx);
  if (mi != ElemProps.end())
    return mi->second;
  else
    return T();
}

template <class T> void ElemProps<T>::clear() { ElemProps.clear(); }

template <class T> const std::map<int, T> &ElemProps<T>::get_properties() const
{
  return ElemProps;
}

template <class T> std::map<int, T> &ElemProps<T>::get_properties()
{
  return ElemProps;
}

template <class T> void ElemProps<T>::remap(const std::map<int, int> &chg_map)
{
  if (!chg_map.size())
    return;
  std::map<int, T> new_props;
  for (const auto &kp : chg_map) {
    if (kp.second != -1) {
      auto cmi = ElemProps.find(kp.first);
      if (cmi != ElemProps.end())
        new_props[kp.second] = cmi->second;
    }
  }

  ElemProps = new_props;
}

template <class T>
const ElemProps<T> &GeomElemProps<T>::operator[](int type) const
{
  return elem_props[type];
}

template <class T> ElemProps<T> &GeomElemProps<T>::operator[](int type)
{
  return elem_props[type];
}

template <class T> void GeomElemProps<T>::clear()
{
  for (int typ = 0; typ < 3; ++typ)
    elem_props[typ].clear();
}

template <class T>
void GeomElemProps<T>::append(const GeomElemProps &geom_props, int v_size,
                              int e_size, int f_size)
{
  int offs[] = {v_size, e_size, f_size};
  for (int i = 0; i < 3; i++) {
    std::map<int, Color>::const_iterator mi;
    for (const auto &pk : geom_props[i].get_properties())
      elem_props[i].set(pk.first + offs[i], pk.second);
  }
}

} // namespace anti

#endif // ELEMPROPS_H
