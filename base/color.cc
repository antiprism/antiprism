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

/*
   Name: color.cc
   Description: representation of color values
   Project: Antiprism - http://www.antiprism.com
*/

#include "color.h"
#include "private_named_cols.h"
#include "utils.h"
#include <cstring>
#include <string>
#include <vector>

using std::max;
using std::min;
using std::string;
using std::vector;

namespace anti {

Color Color::invisible(0, 0, 0, 0);

bool Color::operator==(Color c) const
{
  if (!is_set() && !c.is_set())
    return true;
  if (is_index() && c.is_index() && get_index() == c.get_index())
    return true;
  if (is_value() && c.is_value() && memcmp(rgba, c.rgba, 4) == 0)
    return true;
  return false;
}

bool operator<(const Color &c1, const Color &c2)
{
  if (!c1.is_set())
    return c2.is_set();

  if (c1.is_index())
    return c2.is_value() || c1.get_index() < c2.get_index();

  return c1.get_long() < c2.get_long();
}

void Color::dump(const char *var, FILE *file) const
{
  if (var)
    fprintf(file, "%s=", var);
  if (!is_set())
    fprintf(file, "(not set)\n");
  else if (is_invisible())
    fprintf(file, "invisible\n");
  else if (is_index())
    fprintf(file, "%d (index)\n", get_index());
  else
    fprintf(file, "(%d,%d,%d,%d)\n", rgba[0], rgba[1], rgba[2], rgba[3]);
}

bool Color::set_complement(Color col)
{
  const Color &base_col = col.is_set() ? col : *this;
  if (base_col.is_value()) {
    for (int i = 0; i < 3; i++)
      rgba[i] = ~base_col.rgba[i];
  }
  else
    unset();

  return is_set();
}

bool Color::set_brightness(double brt_val, Color col)
{
  const Color &base_col = col.is_set() ? col : *this;
  if (base_col.is_value()) {
    if (brt_val > 1.0)
      brt_val = 1.0;
    else if (brt_val < -1.0)
      brt_val = -1.0;

    for (int i = 0; i < 3; i++) {
      if (brt_val > 0)
        rgba[i] = (unsigned char)(255 * brt_val +
                                  (1 - brt_val) * base_col.rgba[i]); // to white
      else
        rgba[i] = (unsigned char)((1 + brt_val) * base_col.rgba[i]); // to black
    }
  }
  else
    unset();

  return is_set();
}

// The following RGB / HSV functions are taken from
// http://www.cs.rit.edu/~ncs/color/t_convert.html

// r,g,b values are from 0 to 1
// h = [0,1], s = [0,1], v = [0,1]
//    if s == 0, then h = -1 (undefined)

static void RGBtoHSV(double r, double g, double b, double *h, double *s,
                     double *v)
{
  double *rgb[] = {&r, &g, &b};
  for (auto &i : rgb) {
    if (*i < 0)
      *i = 0;
    else if (*i > 1)
      *i = 1;
  }

  double min, max, delta;

  if (r <= g && r <= b)
    min = r;
  else if (g <= r && g <= b)
    min = g;
  else
    min = b;

  if (r >= g && r >= b)
    max = r;
  else if (g >= r && g >= b)
    max = g;
  else
    max = b;

  *v = max; // v

  delta = max - min;

  if (max > epsilon)  // avoid division by zero
    *s = delta / max; // s
  else {
    // r = g = b = 0     // s = 0, v and h are not important
    *s = 0;
    *h = 0;
    return;
  }

  if (delta < epsilon) // grey range, avoid division by zero
    *h = 0;
  else if (r == max)
    *h = (g - b) / delta; // between yellow & magenta
  else if (g == max)
    *h = 2 + (b - r) / delta; // between cyan & yellow
  else
    *h = 4 + (r - g) / delta; // between magenta & cyan

  *h *= 60; // degrees
  if (*h < 0)
    *h += 360;
  *h /= 360;
}

static void HSVtoRGB(double *r, double *g, double *b, double h, double s,
                     double v)
{
  // bring h into range 0-360
  h = fmod(h, 1.0) * 360;
  if (h < 0)
    h += 360;

  if (s == 0) {
    // achromatic (grey)
    *r = *g = *b = v;
    return;
  }

  if (h > 360 - 0.0001 || h < 0.0001)
    h = 0.01;
  h /= 60; // sector 0 to 5
  int i = (int)floor(h);
  double f = h - i; // factorial part of h
  double p = v * (1 - s);
  double q = v * (1 - s * f);
  double t = v * (1 - s * (1 - f));

  switch (i) {
  case 0:
    *r = v;
    *g = t;
    *b = p;
    break;
  case 1:
    *r = q;
    *g = v;
    *b = p;
    break;
  case 2:
    *r = p;
    *g = v;
    *b = t;
    break;
  case 3:
    *r = p;
    *g = q;
    *b = v;
    break;
  case 4:
    *r = t;
    *g = p;
    *b = v;
    break;
  default: // case 5:
    *r = v;
    *g = p;
    *b = q;
    break;
  }
}

bool Color::set_hsva(double hue, double sat, double val, double alpha)
{
  double *hsva[] = {&hue, &sat, &val, &alpha};
  for (int i = 1; i < 4; i++) { // skip i=0 as hue can wrap
    if (*hsva[i] < 0)
      *hsva[i] = 0;
    else if (*hsva[i] > 1)
      *hsva[i] = 1;
  }

  double r, g, b;
  HSVtoRGB(&r, &g, &b, hue, sat, val);

  return set_rgba(r, g, b, alpha);
}

bool Color::set_hsva(const Vec4d &hsva)
{
  return set_hsva(hsva[0], hsva[1], hsva[2], hsva[3]);
}

Vec4d Color::get_hsva() const
{
  Vec4d hsva;
  Vec4d rgba_d = get_vec4d();
  RGBtoHSV(rgba_d[0], rgba_d[1], rgba_d[2], &hsva[0], &hsva[1], &hsva[2]);
  hsva[3] = rgba_d[3];

  return hsva;
}

// HSL algorithms by Paul Bourke
// http://local.wasp.uwa.edu.au/~pbourke/texture_colour/convert/
/*
   Calculate HSL from RGB
   Hue is in degrees
   Lightness is between 0 and 1
   Saturation is between 0 and 1
*/
static Vec4d RGB2HSL(const Color &c)
{
  Vec4d c1(c[0], c[1], c[2], c[3]);

  double themin = 0.0;
  double themax = 0.0;
  double delta = 0.0;

  themin = min(c1[0], min(c1[1], c1[2]));
  themax = max(c1[0], max(c1[1], c1[2]));
  delta = themax - themin;
  double l = (themin + themax) / 2.0;
  l /= 255.0; // Antiprism
  double s = 0.0;
  if (l > 0.0 && l < 1.0)
    s = delta / (l < 0.5 ? (2.0 * l) : (2.0 - 2.0 * l));
  s /= 255.0; // Antiprism
  double h = 0.0;
  if (delta > 0.0) {
    if (themax == c1[0] && themax != c1[1])
      h += (c1[1] - c1[2]) / delta;
    if (themax == c1[1] && themax != c1[2])
      h += (2.0 + (c1[2] - c1[0]) / delta);
    if (themax == c1[2] && themax != c1[0])
      h += (4.0 + (c1[0] - c1[1]) / delta);
    h *= 60.0;
  }
  if (h < 0.0)
    h += 360.0; // Antiprism

  return (Vec4d(h / 360.0, s, l, c[3] / 255.0));
}

/*
   Calculate RGB from HSL, reverse of RGB2HSL()
   Hue is in degrees
   Lightness is between 0 and 1
   Saturation is between 0 and 1
*/
static Color HSL2RGB(Vec4d c1)
{
  c1[0] = fmod(c1[0], 1.0) * 360; // Antiprism
  if (c1[0] < 0)
    c1[0] += 360;

  Vec4d c2, sat, ctmp;

  while (c1[0] < 0.0)
    c1[0] += 360.0;
  while (c1[0] > 360.0)
    c1[0] -= 360.0;

  if (c1[0] < 120.0) {
    sat[0] = (120.0 - c1[0]) / 60.0;
    sat[1] = c1[0] / 60.0;
    sat[2] = 0.0;
  }
  else if (c1[0] < 240.0) {
    sat[0] = 0.0;
    sat[1] = (240.0 - c1[0]) / 60.0;
    sat[2] = (c1[0] - 120.0) / 60.0;
  }
  else {
    sat[0] = (c1[0] - 240.0) / 60.0;
    sat[1] = 0.0;
    sat[2] = (360.0 - c1[0]) / 60.0;
  }
  sat[0] = min(sat[0], 1.0);
  sat[1] = min(sat[1], 1.0);
  sat[2] = min(sat[2], 1.0);

  ctmp[0] = 2.0 * c1[1] * sat[0] + (1 - c1[1]);
  ctmp[1] = 2.0 * c1[1] * sat[1] + (1 - c1[1]);
  ctmp[2] = 2.0 * c1[1] * sat[2] + (1 - c1[1]);

  if (c1[2] < 0.5) {
    c2[0] = c1[2] * ctmp[0];
    c2[1] = c1[2] * ctmp[1];
    c2[2] = c1[2] * ctmp[2];
  }
  else {
    c2[0] = (1.0 - c1[2]) * ctmp[0] + 2.0 * c1[2] - 1.0;
    c2[1] = (1.0 - c1[2]) * ctmp[1] + 2.0 * c1[2] - 1.0;
    c2[2] = (1.0 - c1[2]) * ctmp[2] + 2.0 * c1[2] - 1.0;
  }

  return (Color(c2[0], c2[1], c2[2], c1[3]));
}

bool Color::set_hsla(double hue, double sat, double light, double alpha)
{
  double *hsla[] = {&hue, &sat, &light, &alpha};
  for (int i = 1; i < 4; i++) {
    if (*hsla[i] < 0)
      *hsla[i] = 0;
    else if (*hsla[i] > 1)
      *hsla[i] = 1;
  }

  *this = HSL2RGB(Vec4d(*hsla[0], *hsla[1], *hsla[2], *hsla[3]));

  return is_set();
}

bool Color::set_hsla(const Vec4d &hsla)
{
  return set_hsla(hsla[0], hsla[1], hsla[2], hsla[3]);
}

Vec4d Color::get_hsla() const { return RGB2HSL(*this); }

Status Color::from_offvals(const vector<char *> &vals, int *col_type)
{
  unset();

  int dummy_type;
  if (col_type == nullptr)
    col_type = &dummy_type;

  Status stat;
  *col_type = -1;
  int i;
  switch (vals.size()) {
  case 0:
    *col_type = 0;
    break;
  case 1:
    if (read_int(vals[0], &i) && i > -1) {
      index = i;
      *col_type = 1;
    }
    else
      stat.set_error(msg_str("colour index '%s' is not a valid "
                             "colour index",
                             vals[0]));
    break;
  case 3:
  case 4: // includes alpha value
    if ((stat = read_intvals(vals)))
      *col_type = 3 + (vals.size() == 4);
    else if ((stat = read_decvals(vals)))
      *col_type = 5 + (vals.size() == 4);
    break;
  default:
    stat.set_error("incorrect number of colour values");
  }

  return (*col_type == 0 || is_set()) ? Status::ok() : stat;
}

Status Color::read_intvals(const vector<char *> &vals)
{
  unset();

  vector<int> ivals;
  Status stat = read_int_list(vals, ivals, false);
  if (stat.is_error())
    return stat;

  return from_intvals(ivals);
}

Status Color::read_intvals(const char *str)
{
  unset();

  vector<int> vals;
  Status stat = read_int_list(str, vals, false);
  if (stat.is_error())
    return stat;

  return from_intvals(vals);
}

Status Color::from_intvals(const vector<int> &vals)
{
  unset();

  const auto vals_sz = vals.size();
  if (vals_sz < 3 || vals_sz > 4)
    return Status::error(msg_str("integer format, %lu numbers "
                                 "given, must give 3 or 4",
                                 (unsigned long)vals_sz));

  for (size_t i = 0; i < vals_sz; i++) {
    if (vals[i] < 0 || vals[i] > 255)
      return Status::error(msg_str("integer format, \"%d\" is "
                                   "not in range 0-255",
                                   vals[i]));
  }

  // default of 255 (opaque) for alpha component
  set_rgba(vals[0], vals[1], vals[2], (vals_sz > 3) ? vals[3] : 255);

  return Status::ok();
}

Status Color::read_decvals(const char *str)
{
  unset();

  vector<double> vals;
  Status stat = read_double_list_noparse(str, vals);
  if (stat.is_error())
    return stat;

  return from_decvals(vals);
}

Status Color::read_decvals(const vector<char *> &vals)
{
  unset();

  vector<double> dvals;
  Status stat = read_double_list_noparse(vals, dvals);
  if (stat.is_error())
    return stat;

  return from_decvals(dvals);
}

Status Color::from_decvals(const vector<double> &vals)
{
  unset();

  const auto vals_sz = vals.size();
  if (vals_sz < 3 || vals_sz > 4)
    return Status::error(msg_str("decimal format, %lu numbers given, "
                                 "must give 3 or 4",
                                 (unsigned long)vals_sz));

  for (size_t i = 0; i < vals_sz; i++) {
    if (vals[i] < 0 || vals[i] > 1)
      return Status::error(msg_str("decimal format, \"%g\" is not "
                                   "in range 0.0-1.0",
                                   vals[i]));
  }

  // default of 1.0 (opaque) for alpha component
  set_rgba(vals[0], vals[1], vals[2], (vals_sz > 3) ? vals[3] : 1.0);

  return Status::ok();
}

Status Color::read_hexvals(const char *str)
{
  unset();
  string hexstr(str); // work with a copy as it will be modified
  if (hexstr.empty() || !strchr("Xx#", hexstr[0]))
    return Status::error(
        msg_str("hex format, first character is not X, x, or #"));

  if (hexstr.size() == 1)
    hexstr += "00000000";
  else if (hexstr.size() == 7)
    hexstr += "FF";

  if (hexstr.size() != 9 ||
      strspn(hexstr.c_str() + 1, "0123456789aAbBcCdDeEfF") != 8)
    return Status::error(msg_str("hex format, %c is not followed by "
                                 "6 or 8 hexadecimal digits",
                                 hexstr[0]));

  unsigned int val;
  sscanf(hexstr.c_str(), "%*c%8x", &val);
  set_rgba(int(val / (256 * 256 * 256)), int((val / (256 * 256)) % 256),
           int((val / 256) % 256), int(val % 256));
  return Status::ok();
}

Status Color::read_hsva_vals(const char *str)
{
  if (strlen(str) == 0)
    return Status::error("hsva format: no format given");
  if (str[0] != 'h' && str[0] != 'H')
    return Status::error("hsva format: does not start with H or h");

  string hsva_str(str); // work with a copy
  vector<double> vals;
  bool hue_degrees = (str[0] == 'h');
  Status stat = read_double_list(&hsva_str[1], vals, 4);
  if (stat.is_error())
    return Status::error(msg_str("hsva format: components: %s", stat.c_msg()));

  int sz = vals.size();
  double hue = vals[0];
  if (hue_degrees)
    hue /= 360.0;
  for (int i = 1; i < sz; i++) {
    if (vals[i] < 0 || vals[i] > 1)
      return Status::error(msg_str("hsva format: components: \"%g\" is "
                                   "not in range 0.0-1.0",
                                   vals[i]));
  }

  set_hsva(hue, (sz < 2 ? 1.0 : vals[1]), (sz < 3 ? 1.0 : vals[2]),
           (sz < 4 ? 1.0 : vals[3]));

  return Status::ok();
}

Status Color::read_colorname(const char *str, bool as_index)
{
  unset();

  string str_cpy(str);         // copy, do not access as C++ string
  char *str_ptr = &str_cpy[0]; // may be used to modify character

  // Colour name 'none' leaves the colour unset
  if (strcmp(str_ptr, "none") == 0 && !as_index)
    return Status::ok();

  // strip whitespace, convert to lowercase
  char *p = str_ptr;
  char *p_to = p;
  while (*p) {
    if (!isspace(*p)) {
      if (isupper(*p))
        *p = tolower(*p);
      *p_to++ = *p;
    }
    p++;
  }
  *p_to = '\0';
  // change grey to gray
  char *pgrey = strstr(str_ptr, "grey");
  if (pgrey)
    pgrey[2] = 'a';

  if (strcmp(str_ptr, "invisible") == 0) {
    if (!as_index)
      *this = invisible;
  }
  else {
    // inefficient, not likely to be called much
    for (unsigned int i = 0; *named_colors[i].name; i++) {
      if (strcmp(str_ptr, named_colors[i].name) == 0) {
        if (as_index)
          set_index(i);
        else {
          set_rgba(named_colors[i].r, named_colors[i].g, named_colors[i].b);
          break;
        }
      }
    }
  }

  return (is_set())
             ? Status::ok()
             : Status::error(msg_str("unknown colour name '%s'", str_ptr));
}

Status Color::read(const char *col_str)
{
  unset();

  string col_str_trimmed(col_str);
  clear_extra_whitespace(col_str_trimmed);
  if (col_str_trimmed.empty()) // don't interpret whitespace only as index 0
    return Status::error("colour name is all whitespace characters");

  Split split;
  if (strchr(col_str_trimmed.c_str(), ',') != nullptr)
    split.init(col_str_trimmed, ",");
  else
    split.init(col_str_trimmed);

  from_offvals(split.get_parts());
  if (is_set())
    return Status::ok();

  if (read_hexvals(col_str_trimmed.c_str()))
    return Status::ok();

  if (read_hsva_vals(col_str_trimmed.c_str()))
    return Status::ok();

  if (read_colorname(col_str_trimmed.c_str()))
    return Status::ok();

  return Status::error("unknown colour name or format '" + col_str_trimmed +
                       "'");
}

namespace col_blend {

/*
RGB->RYB
0 >= 60 degrees, multiply by 2
60 >= 120 degrees, add 60
120 >= 180 maps to 180 to 210
180 >= 240 maps to 210 to 240

So, compression happens
at 120 it is adding 60
at 180 it is adding 30
at 240 it is adding 0

Then
120 >= 240 degrees, N+((240-N)/2)


RYB->RGB
0 >= 120 degrees, divide by 2
120 >= 180 degrees, subtract 60
180 maps to 120
210 maps to 180
240 maps to itself

So, decompression happens
at 180 it is subtracting 60
at 210 it is subtracting 30
at 240 it is subtracting 0

Then
180 >= 240 degrees, N-(240-N)
*/

// angle represented by 0 to 360 degrees
// input: HSV/HSL angle
// output: angle adjusted for RYB mode

double hsx_to_ryb(double angle)
{
  if (angle > 0.0 && angle <= 60.0)
    angle *= 2.0;
  else if (angle > 60.0 && angle <= 120.0)
    angle += 60.0;
  else if (angle > 120.0 && angle <= 240.0)
    angle += (240.0 - angle) / 2;

  return angle;
}

// angle represented by 0 to 360 degrees
// input: angle adjusted for RYB mode
// output: HSV/HSL angle
double ryb_to_hsx(double angle)
{
  if (angle > 0.0 && angle <= 120.0)
    angle /= 2.0;
  else if (angle > 120.0 && angle <= 180.0)
    angle -= 60.0;
  else if (angle > 180.0 && angle <= 240.0)
    angle -= (240.0 - angle);

  return angle;
}

Color rgb_complement(const Color &col, bool ryb_mode)
{
  if (!col.is_value() || col.is_invisible())
    return col;

  // only need hue so algorithm doesn't matter
  Vec4d hsxa = col.get_hsva();
  double angle = rad2deg(2 * M_PI * hsxa[0]);

  if (ryb_mode)
    angle = hsx_to_ryb(angle);

  angle += 180.0;
  if (angle >= 360.0)
    angle -= 360.0;

  if (ryb_mode)
    angle = ryb_to_hsx(angle);
  angle /= 360.0;

  Color rcol;
  rcol.set_hsva(angle, hsxa[1], hsxa[2], hsxa[3]);

  return rcol;
}

// wrapper for multiple HSV/HSL algorithms
Vec4d get_hsxa(const Color &col, int color_system_mode)
{
  return (color_system_mode == 1) ? col.get_hsva() : col.get_hsla();
}

Color set_hsxa(double hue, double sat, double val, double alpha,
               int color_system_mode)
{
  Color col;
  if (color_system_mode <= 1)
    col.set_hsva(hue, sat, val, alpha);
  else if (color_system_mode == 2)
    col.set_hsla(hue, sat, val, alpha);

  return col;
}

// core code furnished by Adrian Rossiter
Color blend_HSX_centroid(const vector<Color> &cols, int color_system_mode,
                         double sat_power, double sat_threshold,
                         double value_power, double value_advance,
                         int alpha_mode, bool ryb_mode)
{
  // no colors, return unset color. one color, return that color
  if (!cols.size())
    return Color();
  else if (cols.size() == 1)
    return cols[0];

  // saturation power can't be 0 or less
  if (sat_power <= 0.0)
    sat_power = 1.0;

  // can't blend two or less colors to black
  if (cols.size() < 3)
    value_power = 0.0;

  double saturation_sum = 0.0;

  double alpha_min = 1.1;
  double alpha_max = -0.1;

  // check for unset, map index, or invisible
  Color map_found;
  bool invisible_found = false;
  bool unset_found = false;

  int cols_sz = cols.size();

  Vec4d sum(0.0, 0.0, 0.0, 0.0);
  for (auto col : cols) {
    // indexes, invisible or unset are not averaged in
    if (col.is_index()) {
      if (!map_found.is_set())
        map_found = col;
      cols_sz--;
      continue;
    }
    else if (col.is_invisible()) {
      invisible_found = true;
      cols_sz--;
      continue;
    }
    else if (!col.is_set()) {
      unset_found = true;
      cols_sz--;
      continue;
    }

    Vec4d hsxa = get_hsxa(col, color_system_mode);

    double S = pow(hsxa[1], sat_power); // map onto distorted disc
    double angle = 2 * M_PI * hsxa[0];

    // RYB mode
    if (ryb_mode)
      angle = deg2rad(hsx_to_ryb(rad2deg(angle)));

    // if value_power is set, simulate subtractive Coloring for 3 or more colors
    double V =
        (value_power <= 0.0)
            ? hsxa[2]
            : pow(fabs(60.0 - fmod(rad2deg(angle) + (ryb_mode ? 60.0 : 0.0) +
                                       value_advance,
                                   120.0)) /
                      60,
                  value_power);

    alpha_min = (hsxa[3] < alpha_min) ? hsxa[3] : alpha_min;
    alpha_max = (hsxa[3] > alpha_max) ? hsxa[3] : alpha_max;

    sum +=
        Vec4d(S * cos(angle), S * sin(angle), V, hsxa[3]); // point in cylinder

    if (sat_threshold < 1.0)
      saturation_sum += hsxa[1];
  }

  // if no colors are being averaged, all are a mix of unset, indexes, and/or
  // invisible
  // heirarchy: first map index, unset, invisible
  // note: invisible should be after unset, else an invisible element placed
  // where an
  //       unset one is will cause it to disappear
  if (!cols_sz) {
    if (map_found.is_set())
      return (map_found);
    else if (unset_found)
      return (Color());
    else if (invisible_found)
      return (Color(Color::invisible));
  }

  // average
  sum /= cols_sz;

  double H = 0.0;
  double S = 0.0;

  // saturation
  S = pow(sum[0] * sum[0] + sum[1] * sum[1],
          0.5 / sat_power); // map back from distorted disc

  // saturations less than 1/255 will happen due to inaccuracy of HSx->RGB
  // conversion
  // note that 255,255,254 has a saturation of 1/255 = 0.00392157..., which is
  // the smallest valid saturation
  if (S < 1 / 255.0)
    S = 0.0;

  // saturation of color centroid is higher than sat_threshold, use average
  // saturation
  if (S > sat_threshold)
    S = saturation_sum / cols_sz;

  // hue
  // if saturation is 0, no need to calculate hue
  // if (double_ne(sum[0],0.0,epsilon) && double_ne(sum[1],0.0,epsilon)) { //
  // old error, made nice pattern, but wrong
  if (S != 0.0) {
    H = atan2(sum[1], sum[0]) / (2 * M_PI);
    if (H < 0)
      H += 1.0;

    // RYB mode
    if (ryb_mode)
      H = deg2rad(ryb_to_hsx(rad2deg(H * 2 * M_PI))) / (2 * M_PI);
  }

  double A =
      (alpha_mode <= 1) ? sum[3] : ((alpha_mode == 2) ? alpha_min : alpha_max);

  Color col = set_hsxa(H, S, sum[2], A, color_system_mode);

  return col;
}

Color blend_RGB_centroid(const vector<Color> &cols, int alpha_mode,
                         bool ryb_mode)
{
  // no colors, return unset color. one color, return that color
  if (!cols.size())
    return Color();
  else if (cols.size() == 1)
    return cols[0];

  double alpha_min = 1.1;
  double alpha_max = -0.1;

  // check for unset, map index, or invisible
  Color map_found;
  bool invisible_found = false;
  bool unset_found = false;

  int cols_sz = cols.size();

  Vec4d col(0.0, 0.0, 0.0, 0.0);
  for (auto i : cols) {
    // indexes, invisible or unset are not averaged in
    if (i.is_index()) {
      if (!map_found.is_set())
        map_found = i;
      cols_sz--;
      continue;
    }
    else if (i.is_invisible()) {
      invisible_found = true;
      cols_sz--;
      continue;
    }
    else if (!i.is_set()) {
      unset_found = true;
      cols_sz--;
      continue;
    }

    Color rcol = i;
    if (ryb_mode) {
      // only need hue so algorithm doesn't matter
      Vec4d hsxa = rcol.get_hsva();
      hsxa[0] = hsx_to_ryb(rad2deg(2 * M_PI * hsxa[0])) / 360.0;
      rcol.set_hsva(hsxa);
    }

    double A = 1.0 - rcol.get_transparency_d();
    alpha_min = (A < alpha_min) ? A : alpha_min;
    alpha_max = (A > alpha_max) ? A : alpha_max;

    col += rcol.get_vec4d();
  }

  // if no colors are being averaged, all are a mix of unset, indexes, and/or
  // invisible
  // heirarchy: first map index, unset, invisible
  // note: invisible should be after unset, else an invisible element placed
  // where an
  //       unset one is will cause it to disappear
  if (!cols_sz) {
    if (map_found.is_set())
      return (map_found);
    else if (unset_found)
      return (Color());
    else if (invisible_found)
      return (Color(Color::invisible));
  }

  col /= cols_sz;

  if (ryb_mode) {
    Color rcol(col);
    // only need hue so algorithm doesn't matter
    Vec4d hsxa = rcol.get_hsva();
    hsxa[0] = ryb_to_hsx(rad2deg(2 * M_PI * hsxa[0])) / 360.0;
    rcol.set_hsva(hsxa);
    col = rcol.get_vec4d();
  }

  if (alpha_mode > 1)
    col[3] = (alpha_mode == 2) ? alpha_min : alpha_max;

  return col;
}
} // namespace col_blend

} // namespace anti
