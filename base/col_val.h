/*
   Copyright (c) 2003-2008, Adrian Rossiter
   
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

/*!\file col_val.h
 * \brief A colour value class
*/


#ifndef COL_VAL_H
#define COL_VAL_H

#include <vector>
#include <map>
#include <string>

#include "vec3d.h"
#include "vec4d.h"

using std::vector;
using std::map;
using std::string;

/// Colour holding RGBA or an index number
class col_val
{
   private:
      int index;
      unsigned char rgba[4];

   public:
      /// Constructor
      /**Initialise with an index number
       * \param idx a positive index number (default value \c -2
       * indicates that no colour is set) */
      col_val(int idx=-2);

      /// Constructor
      /**Initialise with integer RGBA values. Using out of range values will
       * leave the colour unset.
       *\param red red value in range \c 0 - \c 255
       *\param green green value in range \c 0 - \c 255
       *\param blue blue value in range \c 0 - \c 255
       *\param alpha alpha value in range \c 0 (clear) - \c 255 (opaque) */
      col_val(int red, int green, int blue, int alpha=255);

      /// Constructor
      /**Initialise with RGB values in a vec3d. Using out of range
       * values will leave the colour unset.
       * \param col components hold in turn the red, green and blue values,
       * in the range \c 0.0 to \c 1.0 */
      col_val(vec3d col);

      /// Constructor
      /**Initialise with RGBA values in a vec4d. Using out of range
       * values will leave the colour unset.
       * \param col components hold in turn the red, green, blue and
       * alpha values, in the range \c 0.0 to \c 1.0 */
      col_val(vec4d col);

      /// Constructor
      /**Initialise with floating point RGBA values. Using out of range
       * values will leave the colour unset.
       *\param red red value in range \c 0.0 - \c 1.0
       *\param green green value in range \c 0.0 - \c 1.0
       *\param blue blue value in range \c 0.0 - \c 1.0
       *\param alpha alpha value in range \c 0.0 (clear) - \c 1.0 (opaque) */
      col_val(double red, double green, double blue, double alpha=1.0);
     
      ///Unset the colour, indicating no colour is set.
      void unset();

      /// Set the colour with an index number
      /** \param idx a positive index number (default value \c -2
       * indicates that no colour is set) */
      void set_idx(int idx=-2);

      ///Set the colour with integer RGBA values.
      /**Using out of range values will leave the colour unset.
       *\param red red value in range \c 0 - \c 255
       *\param green green value in range \c 0 - \c 255
       *\param blue blue value in range \c 0 - \c 255
       *\param alpha alpha value in range \c 0 (clear) - \c 255 (opaque) */
      void set_rgba(int red, int green, int blue, int alpha=255);

      ///Set the colour with floating point RGBA values.
      /**Using out of range values will leave the colour unset.
       *\param red red value in range \c 0.0 - \c 1.0
       *\param green green value in range \c 0.0 - \c 1.0
       *\param blue blue value in range \c 0.0 - \c 1.0
       *\param alpha alpha value in range \c 0.0 (clear) - \c 1.0 (opaque) */
      void set_rgba(double red, double green, double blue, double alpha=1.0);
     
      ///Set the colour with floating point HSVA values.
      /**Using out of range values will leave the colour unset.
       *\param hue hue in range \c 0.0 - \c 1.0
       *\param sat saturation in range \c 0.0 - \c 1.0
       *\param val value in range \c 0.0 - \c 1.0
       *\param alpha alpha value in range \c 0.0 (clear) - \c 1.0 (opaque) */
      void set_hsva(double hue, double sat, double val, double alpha=1.0);
      
      ///Set the colour with floating point HSVA values.
      /**The HSVA values are in the range 0.0-1.0. Using out of range
       * values will leave the colour unset. */
      void set_hsva(const vec4d &hsva);
      
      ///Set the colour with floating point hsla values.
      /**Using out of range values will leave the colour unset.
       *\param hue hue in range \c 0.0 - \c 1.0
       *\param sat saturation in range \c 0.0 - \c 1.0
       *\param lum value in range \c 0.0 - \c 1.0
       *\param alpha alpha value in range \c 0.0 (clear) - \c 1.0 (opaque) */
      void set_hsla(double hue, double sat, double lum, double alpha=1.0);
      
      ///Set the colour with floating point hsla values.
      /**The hsla values are in the range 0.0-1.0. Using out of range
       * values will leave the colour unset. */
      void set_hsla(const vec4d &hsla);

      ///Get the index number
      /**\return The index number. */
      int get_idx() const;

      ///Get the transparency in integer format
      /**\return The transparency in the range \c 0 (opaque) - \c 255 (clear).*/
      int get_trans() const;

      ///Get the transparency in floating point format
      /**\return The transparency in the range
       * \c 0.0 (opaque) - \c 1.0 (clear).*/
      double get_transd() const;

      ///Get the RGB values
      /**\return The RGB values as components in range \c 0.0 to \c 1.0. */
      vec3d get_vec3d() const;

      ///Get the RGBA values
      /**\return The RGBA values as components in range \c 0.0 to \c 1.0. */
      vec4d get_vec4d() const;

      ///Get the RGBA values as a long integer
      /**\return The RGBA values converted to a long using 2^12R+2^8G+2^4G+A */
      unsigned long get_long() const;

      ///Get the HSVA values
      /**\return The RGBA values as components in range \c 0.0 to \c 1.0. */
      vec4d get_hsva() const;

     ///Get the HSLA values
      /**\return The RGBA values as components in range \c 0.0 to \c 1.0. */
      vec4d get_hsla() const;
      
      ///Read access to the integer RGBA values
      int operator [](int i) const;

      ///Check whether a colour is set
      /**\return \c true if the colour has been set to a valid value,
       * otherwise \c false. */
      bool is_set() const;
      
      ///Check whether a colour is held as RGBA
      /**\return \c true if the colour holds an RGB value,
       * otherwise \c false. */
      bool is_val() const;

      ///Check whether a colour is held as an index number
      /**\return \c true if the colour holds an index number,
       * otherwise \c false. */
      bool is_idx() const;
      
      ///Check whether a colour has the default unset colour state
      /**The opposite of \c is_set()
       * \return \c true if the colour has not been set to a valid value,
       * otherwise \c false. */
      bool is_def() const;

      ///Check whether a colour is invisible
      /**\return \c true if the colour holds the invisible colour value,
       * otherwise \c false. */
      bool is_inv() const;
      
      ///Check for equality
      /**\return true if both have equal index numbers, or both have
       * equal RGBA values, or both are not yet set, otherwise return false. */
      bool operator == (col_val c) const;
      
      ///Check for inequality
      /**\return true if colours have unequal index numbers, or they have
       * unequal RGBA values, or only one is set, * otherwise return false. */
      bool operator !=(col_val c) const;
     
    
      ///Read a colour given as text, detecting the format.
      /**Will read any of the formats below that work with a string
       * \param col_str the color string.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if a valid colour was read, otherwise false
       * and the error is detailed in \a errmsg. */
      bool read(char *col_str, char *errmsg=0);

      ///Read an OFF colour given as number strings.
      /**Read decimal values, integer values or an index number
       * \param vals numbers held as strings.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return an error message if \a col_type is \c 0.
       * \param col_type used to return the colour type that was read
       *    - \c 0 - not valid
       *    - \c 1 - one integer (index number)
       *    - \c 2 - three integer values (RGB)
       *    - \c 3 - four integer values (RGBA)
       *    - \c 4 - three decimal values (RGB)
       *    - \c 5 - four decimal values (RGBA)
       * \return true if a valid colour was read or there was no color to read,
          * otherwise false */
      bool from_offvals(vector<char *> &vals, char *errmsg=0, int *col_type=0);
      
      ///Read a colour given as floating point values.
      /**\param vals 3 (RGB) or 4 (RGBA) values in range \c 0.0 - \c 1.0.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if a valid colour was read, otherwise false
       * and the error is detailed in \a errmsg. */
      bool from_decvals(vector<double> &vals, char *errmsg=0);
      
      ///Read a colour given as integers.
      /**\param vals 3 (RGB) or 4 (RGBA) values in range \c 0 - \c 255.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if a valid colour was read, otherwise false
       * and the error is detailed in \a errmsg. */
      bool from_intvals(vector<int> &vals, char *errmsg=0);
      
      ///Read a colour given as decimal strings.
      /**\param vals 3 (RGB) or 4 (RGBA) numbers in range
       * \c 0.0 - \c 1.0 held as strings.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if a valid colour was read, otherwise false
       * and the error is detailed in \a errmsg. */
      bool read_decvals(vector<char *> &vals, char *errmsg=0);
      
      ///Read a colour given as integer strings.
      /**\param vals 3 (RGB) or 4 (RGBA) numbers in range
       * \c 0 - \c 255 held as strings.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if a valid colour was read, otherwise false
       * and the error is detailed in \a errmsg. */
      bool read_intvals(vector<char *> &vals, char *errmsg=0);

      ///Read a colour given as decimals in a string.
      /**\param str a string of 3 (RGB) or 4 (RGBA) comma separated
       * decimals in range \c 0.0 - \c 1.0.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if a valid colour was read, otherwise false
       * and the error is detailed in \a errmsg. */
      bool read_decvals(char *str, char *errmsg=0);

      ///Read a colour given as integers in a string.
      /**\param str a string of 3 (RGB) or 4 (RGBA) comma separated
       * integers in range \c 0 - \c 255.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if a valid colour was read, otherwise false
       * and the error is detailed in \a errmsg. */
      bool read_intvals(char *str, char *errmsg=0);

      ///Read a colour given as hex format in a string.
      /**\param str a string starting with X, x or # followed by
       * 3 (RGB) or 4 (RGBA) pairs of hexadecimal digits representing
       * integers in range \c 0 - \c 255.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if a valid colour was read, otherwise false
       * and the error is detailed in \a errmsg. */
      bool read_hexvals(char *str, char *errmsg=0);
      
      ///Read a colour given as HSVA format in a string.
      /**\param str a string starting with H or h followed by
       * param str a string of 1 (h) value of plus/minus N degrees or
       * 1 (H) decimal value in range \c 0.0 - \c 1.0.
       * and up to 3 (SVA) optional comma separated decimals in range \c 0.0 - \c 1.0.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \return true if a valid colour was read, otherwise false
       * and the error is detailed in \a errmsg. */
      bool read_hsva_vals(char *str, char *errmsg=0);

      ///Read a colour given by name.
      /**\param str the colour name to be looked up in the internal list
       * of colours.
       * \param errmsg an array at least \c MSG_SZ chars long to
       * return any error message.
       * \param as_index if true then set the colour as the colour
       * index number in the X11 color map.
       * \return true if a valid colour was read, otherwise false
       * and the error is detailed in \a errmsg. */
      bool read_colorname(char *str, char *errmsg=0, bool as_index=false);
     
      ///The invisible colour for non-display elements.
      static col_val invisible;

      ///Convert from floating point to integer format.
      /**\param f a floating point number in range \c 0.0 - \c 1.0.
       * \return an integer in range \c 0 - \c 255 */
      static int f2i(double f);

      ///Convert from integer to floating point format.
      /**\param i an integer in range \c 0 - \c 255.
       * \return a floating point number in range \c 0.0 - \c 1.0 */
      static double i2f(int i);
     
      ///Check that an integer is in the valid colour range.
      /**\param i number to check.
       * \return \c true if \a i is in range \c 0 - \c 255,
       * otherwise \c false. */
      static bool check(int i);

      ///Check that integer components are in the valid colour range.
      /**\param r the red component.
       * \param g the green component.
       * \param b the blue component.
       * \param a the alpha component.
       * \return \c true if all components are in the range \c 0 - \c 255,
       * otherwise \c false. */
      static bool in_range(int r, int g, int b, int a);

};

///Less than
/**Order by invalid first, then index number then RGBA component values.
 * \param c1 first colour to compare.
 * \param c2 second colour to compare.
 * \return true if \arg c1 \c < \arg c2, otherwise false. */
bool operator <(const col_val &c1, const col_val &c2);
 


// Implementation of inline functions

inline col_val::col_val(int idx)
{ 
   set_idx(idx);
}

inline col_val::col_val(int r, int g, int b, int a)
{ 
   set_rgba(r, g, b, a);
}

inline col_val::col_val(vec3d c)
{ 
   set_rgba(c[0], c[1], c[2]);
}

inline col_val::col_val(vec4d c)
{ 
   set_rgba(c[0], c[1], c[2], c[3]);
}

inline col_val::col_val(double r, double g, double b, double a)
{ 
   set_rgba(r, g, b, a); 
}
      
inline void col_val::unset()
{ 
   set_idx(); 
}

inline void col_val::set_idx(int idx)
{
   for(int i=0; i<4; i++)
      rgba[i] = 0;
   index = (idx<0) ? -2 : idx;
}

inline void col_val::set_rgba(int r, int g, int b, int a)
{ 
   rgba[0] = r;
   rgba[1] = g;
   rgba[2] = b;
   rgba[3] = a;
   index = -2 + in_range(r,g,b,a);
}

inline void col_val::set_rgba(double r, double g, double b, double a)
{ 
   set_rgba(f2i(r), f2i(g), f2i(b), f2i(a));
}
      
inline bool col_val::is_val() const
{
   return index==-1;
}

inline bool col_val::is_idx() const
{ 
   return index>=0;
}

inline bool col_val::is_def() const
{
   return !is_set();
}

inline bool col_val::is_inv() const
{ 
   return *this == invisible;
}

inline bool col_val::is_set() const
{ 
   return (is_val() || is_idx());
}

inline bool col_val::operator !=(col_val c) const
{
   return !(*this==c);
}

inline int col_val::get_idx() const
{
   return is_idx() ? index : -2;
}

inline int col_val::get_trans() const
{
   return 255 - rgba[3];
}

inline double col_val::get_transd() const
{
   return 1- i2f(rgba[3]);
}

inline vec3d col_val::get_vec3d() const
{
   return vec3d(i2f(rgba[0]), i2f(rgba[1]), i2f(rgba[2]));
}

inline vec4d col_val::get_vec4d() const
{
   return vec4d(i2f(rgba[0]), i2f(rgba[1]), i2f(rgba[2]), i2f(rgba[3]));
}

inline unsigned long col_val::get_long() const
{
   return (1UL<<12)*rgba[0] + (1UL<<8)*rgba[1] + (1UL<<4)*rgba[2] + rgba[3];
}
      
inline int col_val::operator [](int i) const
{
   return rgba[i];
}

inline int col_val::f2i(double f)
{
   return (int)(f*255+0.5);
}

inline double col_val::i2f(int i)
{
   return i/255.0;
}

inline bool col_val::check(int i)
{
   return i>=0 && i<=255;
}

inline bool col_val::in_range(int r, int g, int b, int a) 
{ 
   return check(r) && check(g) && check(b) && check(a);
}



#endif // COL_VAL_H

