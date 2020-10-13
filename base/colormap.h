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

/*!\file colormap.h
   \brief A colour map class
*/

#ifndef COLORMAP_H
#define COLORMAP_H

#include "color.h"
#include "geometry.h"
#include "random.h"
#include "status.h"
#include <limits.h>
#include <map>
#include <vector>

namespace anti {

/// A colour map class
class ColorMap {
private:
  int shift; // the first entry starts this many places into the map.
  int step;  // step this many places with each index increment.
  int wrap;  // wrap index numbers back to zero at this index size.

protected:
  /// Copy parameters
  /**\param cmap the colour map to copy parameters from. */
  void copy_params(const ColorMap &cmap);

public:
  enum { max_map_sz = INT_MAX }; // The largest size of a map

  /// Constructor
  ColorMap() : shift(0), step(1), wrap(0) {}

  /// Destructor
  virtual ~ColorMap() = default;

  /// Initialise from a string
  /**\param params the map parameters.
   * \return status, evaluates to \c true if the map could be initialised
   *  (possibly with warnings) otherwise \c false. */
  virtual Status init(const char *params);

  /// Initialise from a string, stripping parameters from the end
  /**\param map_name the map name.
   * \return status, evaluates to \c true if the map could be initialised
   *  (possibly with warnings) otherwise \c false. */
  Status init_strip(char *map_name);

  /// Read parameters
  /**\param params the map name to read parameters from.
   * \return status, evaluates to \c true if the parameters could be read
   *  (possibly with warnings) otherwise \c false. */
  Status read_params(const char *params);

  /// Get a copy of the map with any with the values indexed by order
  /**\return a pointer to the dynamically allocated condensed copy,
   *  which must be freed by the caller with \c delete, 0 indicates
   *  that the copy failed. */
  virtual ColorMap *get_condensed() const { return clone(); }

  /// Get a copy of the map
  /**\return a pointer to the dynamically allocated copy,
   *  which must be freed by the caller with \c delete, 0 indicates
   *  that the clone failed. */
  virtual ColorMap *clone() const { return new ColorMap(*this); }

  /// The effective size of the map
  /** The effective size of a map is one greater than the highest
   *  index number in the map. It is the size of the smallest
   *  map (sequential, starting at 0) that will include all
   *  the entries of the map.
   * \return The effective size */
  virtual int effective_size() const { return max_map_sz; };

  /// Get the map shift
  /** Lookup of colour values is determined by (shift + index*step)%wrap
   * \return the shift */
  int get_shift() const { return shift; }

  /// Set the map shift
  /** Lookup of colour values is determined by (shift + index*step)%wrap
   * \param shft the shift */
  void set_shift(int shft) { shift = shft; }

  /// Get the map step
  /** Look up of colour values is determined by (shift + index*step)%wrap
   * \return the step */
  int get_step() const { return step; }

  /// Set the map shift
  /** Look up of colour values is determined by (shift + index*step)%wrap
   * \param stp the step */
  void set_step(int stp) { step = stp; }

  /// Get the map wrap size
  /** Look up of colour values is determined by (shift + index*step)%wrap
   *  a value of 0 indicates no wrapping.
   * \return the wrap */
  int get_wrap() const { return wrap; }

  /// Set the map wrap size
  /** Look up of colour values is determined by (shift + index*step)%wrap
   *  a value of 0 indicates no wrapping, a negative value indicates
   *  wrapping after the highest index number.
   * \param wrp the wrap */
  void set_wrap(int wrp = -1);

  /// Get effective index
  /** The effective index is (shift + index*step)%wrap
   * \param idx the index
   * \return The effective index. */
  int get_effective_index(int idx) const;

  /// Get the colour value for an index number.
  /**
   * \return The colour. */
  virtual Color get_col(int /*idx*/) const { return Color(); }

  /// Cycle the map colours
  /** Each colour index number is mapped to the previous colour
   *  value in the colour map. */
  void cycle_colors();
};

/// A colour map with the mappings held in a map
class ColorMapMap : public ColorMap {
private:
  std::map<int, Color> cmap; // color map

public:
  /// Constructor
  ColorMapMap() : ColorMap() {}

  /// Initialise from a file
  /** The colour map can be in the Antiprism, GIMP or Fractint format. If
   *  the filename isn't found then the name will be looked for in the
   *  Antiprism data directory colour map resources.
   * \param map_name the map file name.
   * \return status, evaluates to \c true if the map could be initialised
   *  (possibly with warnings) otherwise \c false. */
  virtual Status init(const char *map_name);

  /// Initialise from a formatted line
  /** The colour map line is converted to Antiprism format in
   *  the following way: ':' is converted to a newline, and '_'
   *  to a space. If the lines are bare colours then "idx=" is
   *  added.
   * \param map_line the map line text.
   * \return status, evaluates to \c true if the map line was valid
   *  (possibly with warnings) otherwise \c false. */
  Status init_from_line(const char *map_line);

  /// Get a copy of the map
  /**\return a pointer to the dynamically allocated copy,
   *  which must be freed by the caller with \c delete, 0 indicates
   *  that the clone failed. */
  ColorMap *clone() const { return new ColorMapMap(*this); }

  /// Get a copy of the map with any with the values indexed by order
  /**\return a pointer to the dynamically allocated condensed copy,
   *  which must be freed by the caller with \c delete, 0 indicates
   *  that the copy failed. */
  virtual ColorMap *get_condensed() const;

  /// Read the named colours into a colour map
  /** This is equivalent to reading in the x11 colour map from
   *  the resources. */
  void read_named_colors();

  /// Set the colour map
  /**\param col_map the colour map to use. */
  void set_map(const std::map<int, Color> col_map) { cmap = col_map; }

  /// Get the colour map
  /**\return the colour map */
  const std::map<int, Color> &get_map() const { return cmap; }

  /// Clear the colour map
  void clear() { cmap.clear(); }

  /// The size of the map
  /**\return The size */
  unsigned int size() const { return cmap.size(); }

  /// The effective size of the map
  /** The effective size of a map is one greater than the highest
   *  index number in the map. It is the size of the smallest
   *  map (sequential, starting at 0) that will include all
   *  the entries of the map.
   * \return The effective size */
  virtual int effective_size() const;

  /// Get the colour value for an index number.
  /**\param idx the index.
   * \return The colour. */
  virtual Color get_col(int idx) const;

  /// Set a colour value for an index number.
  /**\param idx the index.
   * \param col the colour to set. */
  void set_col(int idx, Color col);
};

/// A colour map that looks up in other colour maps in order
class ColorMapMulti : public ColorMap {
private:
  // The colour maps to be tried sequentially.
  std::vector<ColorMap *> cmaps;

  // Largest effective size
  int max_eff_map_sz;

  // Specified map size
  int map_sz;

  // Set map_sz to the largest effective size
  void set_max_eff_map_sz();

public:
  /// Constructor
  ColorMapMulti() : ColorMap(), max_eff_map_sz(0), map_sz(-1) {}

  /// Copy Constructor
  /**\param cmap the multiple colour map to copy from. */
  ColorMapMulti(const ColorMapMulti &cmap);

  /// Copy Assignment
  /**\param cmap the multiple colour map to copy from.
   * \return A reference to this object.*/
  ColorMapMulti &operator=(const ColorMapMulti &cmap);

  /// Destructor
  ~ColorMapMulti();

  /// Initialise from a string
  /** The colour map can be in the Antiprism, GIMP or Fractint format. If
   *  the filename isn't found then the name will be looked for in the
   *  Antiprism data directory colour map resources.
   * \param map_name the map name.
   * \return status, evaluates to \c true if the map could be initialised
   *  (possibly with warnings) otherwise \c false. */
  virtual Status init(const char *map_name);

  /// Get a copy of the map
  /**\return a pointer to the dynamically allocated copy,
   *  which must be freed by the caller with \c delete, 0 indicates
   *  that the clone failed. */
  virtual ColorMap *clone() const { return new ColorMapMulti(*this); }

  /// Set the map size
  /**\param sz the number of entries in the map. */
  void set_map_sz(int sz) { map_sz = sz; }

  /// The effective size of the map
  /** The effective size of a map is one greater than the highest
   *  index number in the map. It is the size of the smallest
   *  map (sequential, starting at 0) that will include all
   *  the entries of the map.
   * \return The effective size */
  virtual int effective_size() const
  {
    return (map_sz >= 0) ? map_sz : max_eff_map_sz;
  }

  /// Add a colour map.
  /** The ColorMap must be dynamically allocated, using \c new directly
   *  or through, for example, \c colormap_from_name(name). The calling
   *  program must not delete it. The \c ColorMapMulti object will
   *  delete the added map when it is no longer needed
   * \param col_map the colour map.
   * \param pos the position to add it, or at the end if pos is
   *  greater then the current size */
  void add_cmap(ColorMap *col_map, unsigned int pos = INT_MAX);

  /// Delete a colour map.
  /**\param pos the position of the colour map to delete, or delete
   *  the last colour map if \c pos is greater than or equal to the
   *  current size */
  void del_cmap(unsigned int pos = INT_MAX);

  /// Get a the colour maps.
  /**\return The colour maps. */
  const std::vector<ColorMap *> &get_cmaps() const { return cmaps; }

  /// Get the colour value for an index number.
  /**\param idx the index.
   * \return The colour. */
  virtual Color get_col(int idx) const;
};

/// Create a colour map from its name
/** The map may be read from a file or generated.
 * \param map_name the map name
 * \param stat status, evaluates to \c true if the map could be initialised
 *  (possibly with warnings) otherwise \c false
 * \return a pointer to the dynamically allocated map which must
 *  be freed by the caller with \c delete, \c nullptr is returned for an
 *  invalid map name and the error is detailed in \a stat. */
ColorMap *colormap_from_name(const char *map_name, Status *stat = nullptr);

/// Colour value map to an HSVA range
class ColorValuesToRangeHsva {
protected:
  Color default_color;
  std::vector<double> ranges[4];
  Status add_range(int idx, const char *rngs);

public:
  /// Constructor
  /**\param range_name string with ranges.
   * \param def_col default colour for unset colours (processed like others). */
  ColorValuesToRangeHsva(const std::string &range_name = "",
                         const Color &def_col = Color());

  /// Initialise
  /**\param range_name string with ranges
   * \param def_col default colour for unset colours (processed like others)
   * \return status, evaluates to \c true if the value map could be initialised
   *  (possibly with warnings) otherwise \c false. */
  Status init(const std::string &range_name, const Color &def_col = Color());

  /// Set the default colour
  /**\param def_col default colour for unset colours (processed like others).*/
  void set_default_color(const Color &def_col) { default_color = def_col; }

  /// Apply processing to the colour values of a geometry
  /**\param geom geometry to map colours for.
   * \param elem_type element type to map colours for. */
  void apply(anti::Geometry &geom, int elem_type);

  /// Apply processing to colour values in a map
  /**\param elem_cols element type to map colours for. */
  void apply(std::map<int, Color> &elem_cols);

  /// Get the processed colour
  /**\param col the color.
   * \return The processed colour. */
  virtual Color get_col(Color col);
};

// -------------------------------------------------------------------
// inline functions

inline int ColorMap::get_effective_index(int idx) const
{
  int eff_idx = shift + idx * step;
  if (wrap > 0)
    eff_idx %= wrap;
  return eff_idx;
}

inline void ColorMap::cycle_colors()
{
  step++;
  if (wrap)
    step %= wrap;
}

} // namespace anti

#endif // COLORMAP_H
