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

/*!\file color.h
 * \brief A colour value class
 */

#ifndef COLOR_H
#define COLOR_H

#include "vec3d.h"
#include "vec4d.h"

#include <vector>

namespace anti {

/// Colour holding RGBA or an index number
class Color {
private:
  int index;
  unsigned char rgba[4];

public:
  /// Constructor
  /** Initialise with an index number
   * \param idx a positive index number (default value \c -2
   *  indicates that no colour is set) */
  Color(int idx = -2);

  /// Constructor
  /** Initialise with integer RGBA values. Using out of range values will
   *  leave the colour unset.
   * \param red red value in range \c 0 - \c 255
   * \param green green value in range \c 0 - \c 255
   * \param blue blue value in range \c 0 - \c 255
   * \param alpha alpha value in range \c 0 (clear) - \c 255 (opaque) */
  Color(int red, int green, int blue, int alpha = 255);

  /// Constructor
  /** Initialise with RGB values in a Vec3d. Using out of range
   *  values will leave the colour unset.
   * \param col components hold in turn the red, green and blue values,
   *  in the range \c 0.0 to \c 1.0 */
  Color(Vec3d col);

  /// Constructor
  /** Initialise with RGBA values in a Vec4d. Using out of range
   *  values will leave the colour unset.
   * \param col components hold in turn the red, green, blue and
   *  alpha values, in the range \c 0.0 to \c 1.0 */
  Color(Vec4d col);

  /// Constructor
  /** Initialise with floating point RGBA values. Using out of range
   *  values will leave the colour unset.
   * \param red red value in range \c 0.0 - \c 1.0
   * \param green green value in range \c 0.0 - \c 1.0
   * \param blue blue value in range \c 0.0 - \c 1.0
   * \param alpha alpha value in range \c 0.0 (clear) - \c 1.0 (opaque) */
  Color(double red, double green, double blue, double alpha = 1.0);

  /// Unset the colour, indicating no colour is set.
  void unset();

  /// Set the colour with an index number
  /**\param idx a positive index number (default value \c -2
   *  indicates that no colour is set)
   * \return \c true if valid index, else \c false */
  bool set_index(int idx = -2);

  /// Set the colour with integer RGBA values.
  /** Using out of range values will leave the colour unset.
   * \param red red value in range \c 0 - \c 255
   * \param green green value in range \c 0 - \c 255
   * \param blue blue value in range \c 0 - \c 255
   * \param alpha alpha value in range \c 0 (clear) - \c 255 (opaque)
   * \return \c true if valid RGBA values, else \c false */
  bool set_rgba(int red, int green, int blue, int alpha = 255);

  /// Set the colour with floating point RGBA values.
  /** Using out of range values will leave the colour unset.
   * \param red red value in range \c 0.0 - \c 1.0
   * \param green green value in range \c 0.0 - \c 1.0
   * \param blue blue value in range \c 0.0 - \c 1.0
   * \param alpha alpha value in range \c 0.0 (clear) - \c 1.0 (opaque)
   * \return \c true if valid RGBA values, else \c false */
  bool set_rgba(double red, double green, double blue, double alpha = 1.0);

  /// Set the colour with floating point HSVA values.
  /** Using out of range values will leave the colour unset.
   * \param hue hue in range \c 0.0 - \c 1.0
   * \param sat saturation in range \c 0.0 - \c 1.0
   * \param val value in range \c 0.0 - \c 1.0
   * \param alpha alpha value in range \c 0.0 (clear) - \c 1.0 (opaque)
   * \return \c true if valid HSVA values, else \c false */
  bool set_hsva(double hue, double sat, double val, double alpha = 1.0);

  /// Set the colour with floating point HSVA values.
  /** The HSVA values are in the range 0.0-1.0. Using out of range
   *  values will leave the colour unset.
   * \param hsva the HSVA values.
   * \return \c true if valid HSVA values, else \c false */
  bool set_hsva(const Vec4d &hsva);

  /// Set the colour with floating point hsla values.
  /** Using out of range values will leave the colour unset.
   * \param hue hue in range \c 0.0 - \c 1.0
   * \param sat saturation in range \c 0.0 - \c 1.0
   * \param light lightness in range \c 0.0 - \c 1.0
   * \param alpha alpha value in range \c 0.0 (clear) - \c 1.0 (opaque)
   * \return \c true if valid HSLA values, else \c false */
  bool set_hsla(double hue, double sat, double light, double alpha = 1.0);

  /// Set the colour with floating point hsla values.
  /** The hsla values are in the range 0.0-1.0. Using out of range
   *  values will leave the colour unset.
   * \param hsla the HSLA values.
   * \return \c true if valid HSLA values, else \c false */
  bool set_hsla(const Vec4d &hsla);

  /// Set the alpha component of a color value
  /**\param alpha the alpha value in range 0 - 255.
   * \return \c true if the color was a valid color value and alpha
   * was valid, else \c false .*/
  bool set_alpha(int alpha);

  /// Set the alpha component of a color value
  /**\param alpha the alpha value in range 0.0 - 1.0.
   * \return \c true if the color was a valid color value and alpha
   * was valid, else \c false .*/
  bool set_alpha(double alpha);

  /// Set the colour to its complement, or the complement of another colour
  /**\param col the colour to take the complement from, or use current
   *  colour if \a col is unset (default). If the base colour is not a
   *  colour value the final colour will be unset.
   * \return \c true if the color was a valid color value, else \c false */
  bool set_complement(Color col = Color());

  /// Set the colour brightness, or from brightening another colour
  /**\param brt_val the value to brighten in the range -1.0 to 1.0.
   *  Positive values brighten towards white and negative values
   *  darken towards black
   * \param col the colour as the base to brighten, or use current
   *  colour if \a col is unset (default). If the base colour is not a
   *  colour value the final colour will be unset.
   * \return \c true if the color was a valid color value, else \c false */
  bool set_brightness(double brt_val, Color col = Color());

  /// Get the index number
  /**\return The index number. */
  int get_index() const;

  /// Get the transparency in integer format
  /**\return The transparency in the range \c 0 (opaque) - \c 255 (clear).*/
  int get_transparency() const;

  /// Get the transparency in floating point format
  /**\return The transparency in the range
   *  \c 0.0 (opaque) - \c 1.0 (clear).*/
  double get_transparency_d() const;

  /// Get the RGB values
  /**\return The RGB values as components in range \c 0.0 to \c 1.0. */
  Vec3d get_vec3d() const;

  /// Get the RGBA values
  /**\return The RGBA values as components in range \c 0.0 to \c 1.0. */
  Vec4d get_vec4d() const;

  /// Get the RGBA values as a long integer
  /**\return The RGBA values converted to a long using 2^12R+2^8G+2^4G+A */
  unsigned long get_long() const;

  /// Get the HSVA values
  /**\return The RGBA values as components in range \c 0.0 to \c 1.0. */
  Vec4d get_hsva() const;

  /// Get the HSLA values
  /**\return The RGBA values as components in range \c 0.0 to \c 1.0. */
  Vec4d get_hsla() const;

  /// Read access to the integer RGBA values
  /**\param idx index number of component
   * \return The component */
  int operator[](int idx) const;

  /// Check whether a colour is set
  /**\return \c true if the colour has been set to a valid value,
   *  otherwise \c false. */
  bool is_set() const;

  /// Check whether a colour is held as RGBA
  /**\return \c true if the colour holds an RGB value,
   *  otherwise \c false. */
  bool is_value() const;

  /// Check whether a colour is held as RGBA and is not the invisible value
  /**\return \c true if the colour holds an RGB value an is not the invisible
   *  value, otherwise \c false. */
  bool is_visible_value() const;

  /// Check whether a colour is held as an index number
  /**\return \c true if the colour holds an index number,
   *  otherwise \c false. */
  bool is_index() const;

  /// Check whether a colour is invisible
  /**\return \c true if the colour holds the invisible colour value,
   *  otherwise \c false. */
  bool is_invisible() const;

  /// Check for equality
  /**\param c the colour to compare with for equality
   * \return \c true if both have equal index numbers, or both have
   *  equal RGBA values, or both are not yet set, otherwise \c false. */
  bool operator==(Color c) const;

  /// Check for inequality
  /**\param c the colour to compare with for equality
   * \return \c true if colours have unequal index numbers, or they have
   *  unequal RGBA values, or only one is set, otherwise return \c false. */
  bool operator!=(Color c) const;

  /// Read a colour given as text, detecting the format.
  /** Will read any of the formats below that work with a string
   * \param col_str the color string.
   * \return status, which evaluates to \c true if a valid colour was read,
   *  otherwise \c false */
  Status read(const char *col_str);

  /// Read an OFF colour given as number strings.
  /** Read decimal values, integer values or an index number
   * \param vals numbers held as strings.
   * \param col_type used to return the colour type that was read
   *    - \c 0 - not valid
   *    - \c 1 - one integer (index number)
   *    - \c 2 - three integer values (RGB)
   *    - \c 3 - four integer values (RGBA)
   *    - \c 4 - three decimal values (RGB)
   *    - \c 5 - four decimal values (RGBA)
   * \return status, which evaluates to \c true if a valid colour was read,
   *  or there was no color to read, otherwise \c false */
  Status from_offvals(const std::vector<char *> &vals, int *col_type = nullptr);

  /// Read a colour given as floating point values.
  /**\param vals 3 (RGB) or 4 (RGBA) values in range \c 0.0 - \c 1.0.
   * \return status, which evaluates to \c true if a valid colour was read,
   *  otherwise \c false */
  Status from_decvals(const std::vector<double> &vals);

  /// Read a colour given as integers.
  /**\param vals 3 (RGB) or 4 (RGBA) values in range \c 0 - \c 255.
   * \return status, which evaluates to \c true if a valid colour was read,
   *  otherwise \c false */
  Status from_intvals(const std::vector<int> &vals);

  /// Read a colour given as decimal strings.
  /**\param vals 3 (RGB) or 4 (RGBA) numbers in range
   *  \c 0.0 - \c 1.0 held as strings.
   * \return status, which evaluates to \c true if a valid colour was read,
   *  otherwise \c false */
  Status read_decvals(const std::vector<char *> &vals);

  /// Read a colour given as integer strings.
  /**\param vals 3 (RGB) or 4 (RGBA) numbers in range
   *  \c 0 - \c 255 held as strings.
   * \return status, which evaluates to \c true if a valid colour was read,
   *  otherwise \c false */
  Status read_intvals(const std::vector<char *> &vals);

  /// Read a colour given as decimals in a string.
  /**\param str a string of 3 (RGB) or 4 (RGBA) comma separated
   *  decimals in range \c 0.0 - \c 1.0.
   * \return status, which evaluates to \c true if a valid colour was read,
   *  otherwise \c false */
  Status read_decvals(const char *str);

  /// Read a colour given as integers in a string.
  /**\param str a string of 3 (RGB) or 4 (RGBA) comma separated
   *  integers in range \c 0 - \c 255.
   * \return status, which evaluates to \c true if a valid colour was read,
   *  otherwise \c false */
  Status read_intvals(const char *str);

  /// Read a colour given as hex format in a string.
  /**\param str a string starting with X, x or # followed by
   *  3 (RGB) or 4 (RGBA) pairs of hexadecimal digits representing
   *  integers in range \c 0 - \c 255.
   * \return status, which evaluates to \c true if a valid colour was read,
   *  otherwise \c false */
  Status read_hexvals(const char *str);

  /// Read a colour given as HSVA format in a string.
  /**\param str a string starting with H or h followed by
   *  param str a string of 1 (h) value of plus/minus N degrees or
   *  1 (H) decimal value in range \c 0.0 - \c 1.0.
   *  and up to 3 (SVA) optional comma separated decimals in
   *  range \c 0.0 - \c 1.0.
   * \return status, which evaluates to \c true if a valid colour was read,
   *  otherwise \c false */
  Status read_hsva_vals(const char *str);

  /// Read a colour given by name.
  /**\param str the colour name to be looked up in the internal list
   *  of colours.
   * \param as_index if \c true then set the colour as the colour
   *  index number in the X11 color map.
   * \return status, which evaluates to \c true if a valid colour was read,
   *  otherwise \c false */
  Status read_colorname(const char *str, bool as_index = false);

  /// The invisible colour for non-display elements.
  static Color invisible;

  /// Convert from floating point to integer format.
  /**\param f a floating point number in range \c 0.0 - \c 1.0.
   * \return an integer in range \c 0 - \c 255 */
  static int f2i(double f);

  /// Convert from integer to floating point format.
  /**\param i an integer in range \c 0 - \c 255.
   * \return a floating point number in range \c 0.0 - \c 1.0 */
  static double i2f(int i);

  /// Check that an integer is in the valid colour range.
  /**\param i number to check.
   * \return \c true if \a i is in range \c 0 - \c 255,
   *  otherwise \c false. */
  static bool check(int i);

  /// Check that integer components are in the valid colour range.
  /**\param r the red component.
   * \param g the green component.
   * \param b the blue component.
   * \param a the alpha component.
   * \return \c true if all components are in the range \c 0 - \c 255,
   *  otherwise \c false. */
  static bool in_range(int r, int g, int b, int a);

  /// Debugging print of a vector variable
  /**\param var a string to identify the vector variable.
   * \param file file stream to print the variable. */
  void dump(const char *var = "", FILE *file = stderr) const;
};

/// Less than
/** Order by invalid first, then index number then RGBA component values.
 * \param c1 first colour to compare.
 * \param c2 second colour to compare.
 * \return \c true if <tt>c1 \c < \c c2</tt>, otherwise \c false. */
bool operator<(const Color &c1, const Color &c2);

namespace col_blend {

/// Convert HSX angle to RYB angle
/**\param angle HSX angle.
 * \return RYB angle. */
double hsx_to_ryb(double angle);

/// Convert RYB angle to HSV/HSL angle
/**\param angle HSX angle.
 * \return RYB angle. */
double ryb_to_hsx(double angle);

/// Get RGB complement
/**\param col base colour.
 * \param ryb_mode \c true treat as RYB color, otherwise \c false treat as RGB.
 * \return complement of colour. */
Color rgb_complement(const Color &col, bool ryb_mode);

/// Convert RGBA colour to HSVA/HSLA
/**\param col colour to convert
 * \param color_system_mode - 1:HSVA, 2:HSLA
 * \return Color component values. */
Vec4d get_hsxa(const Color &col, int color_system_mode);

/// Set the colour with floating point HSVA values.
/**Using out of range values will leave the colour unset.
 * \param hue hue in range \c 0.0 - \c 1.0
 * \param sat saturation in range \c 0.0 - \c 1.0
 * \param val value in range \c 0.0 - \c 1.0
 * \param alpha alpha value in range \c 0.0 (clear) - \c 1.0 (opaque)
 * \param color_system_mode - 1:HSVA, 2:HSLA.
 * \return The color. */
Color set_hsxa(double hue, double sat, double val, double alpha,
               int color_system_mode);

/// Blend colors as HSVA/HSLA centroid
/**\param cols colours to blend.
 * \param color_system_mode - 1:HSVA, 2:HSLA.
 * \param sat_power saturation power
 * \param sat_threshold saturation threshold.
 * \param value_power value power
 * \param value_advance value avance
 * \param alpha_mode alpha mode.
 * \param ryb_mode \c true treat as RYB color, otherwise \c false treat as RGB.
 * \return the color resulting from the blend. */
Color blend_HSX_centroid(const std::vector<Color> &cols,
                         int color_system_mode = 2, double sat_power = 0,
                         double sat_threshold = 1.0, double value_power = 0,
                         double value_advance = 0, int alpha_mode = 3,
                         bool ryb_mode = false);

/// Blend colors as RGB centroid
/**\param cols colours to blend.
 * \param alpha_mode alpha mode.
 * \param ryb_mode \c true treat as RYB color, otherwise \c false treat as RGB.
 * \return the color resulting from the blend. */
Color blend_RGB_centroid(const std::vector<Color> &cols, int alpha_mode = 3,
                         bool ryb_mode = false);
} // namespace col_blend

// ----------------------------------------------------------------------
// inline functions

inline Color::Color(int idx) { set_index(idx); }

inline Color::Color(int red, int green, int blue, int alpha)
{
  set_rgba(red, green, blue, alpha);
}

inline Color::Color(double red, double green, double blue, double alpha)
{
  set_rgba(red, green, blue, alpha);
}

inline Color::Color(Vec3d col) { set_rgba(col[0], col[1], col[2]); }

inline Color::Color(Vec4d col) { set_rgba(col[0], col[1], col[2], col[3]); }

inline void Color::unset() { set_index(); }

inline bool Color::set_index(int idx)
{
  for (auto &comp : rgba)
    comp = 0;
  index = (idx < 0) ? -2 : idx;
  return is_set();
}

inline bool Color::set_rgba(int red, int green, int blue, int alpha)
{
  rgba[0] = red;
  rgba[1] = green;
  rgba[2] = blue;
  rgba[3] = alpha;
  index = -2 + in_range(red, green, blue, alpha);
  return is_set();
}

inline bool Color::set_rgba(double red, double green, double blue, double alpha)
{
  return set_rgba(f2i(red), f2i(green), f2i(blue), f2i(alpha));
}

inline bool Color::is_value() const { return index == -1; }

inline bool Color::is_visible_value() const
{
  return is_value() && !is_invisible();
}

inline bool Color::is_index() const { return index >= 0; }

inline bool Color::is_invisible() const { return *this == invisible; }

inline bool Color::is_set() const { return index > -2; } // is_val()||is_idx()

inline bool Color::operator!=(Color c) const { return !(*this == c); }

inline bool Color::set_alpha(int alpha)
{
  bool valid = false;
  if (is_visible_value() && check(alpha)) {
    rgba[3] = alpha;
    valid = true;
  }
  return valid;
}

inline bool Color::set_alpha(double alpha) { return set_alpha(f2i(alpha)); }

inline int Color::get_index() const { return is_index() ? index : -2; }

inline int Color::get_transparency() const { return 255 - rgba[3]; }

inline double Color::get_transparency_d() const { return 1 - i2f(rgba[3]); }

inline Vec3d Color::get_vec3d() const
{
  return Vec3d(i2f(rgba[0]), i2f(rgba[1]), i2f(rgba[2]));
}

inline Vec4d Color::get_vec4d() const
{
  return Vec4d(i2f(rgba[0]), i2f(rgba[1]), i2f(rgba[2]), i2f(rgba[3]));
}

inline unsigned long Color::get_long() const
{
  return (1UL << 12) * rgba[0] + (1UL << 8) * rgba[1] + (1UL << 4) * rgba[2] +
         rgba[3];
}

inline int Color::operator[](int idx) const { return rgba[idx]; }

inline int Color::f2i(double f) { return (int)(f * 255 + 0.5); }

inline double Color::i2f(int i) { return i / 255.0; }

inline bool Color::check(int i) { return i >= 0 && i <= 255; }

inline bool Color::in_range(int r, int g, int b, int a)
{
  return check(r) && check(g) && check(b) && check(a);
}

} // namespace anti

#endif // COLOR_H
