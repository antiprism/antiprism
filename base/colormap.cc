/*
   Copyright (c) 2008-2023, Adrian Rossiter

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

/*!\file col_map.h
   \brief A colour map class
*/

#include "coloring.h"
#include "geometry.h"
#include "private_named_cols.h"
#include "random.h"
#include "utils.h"

#include <algorithm>
#include <cctype>
#include <cstring>
#include <map>

using std::map;
using std::string;
using std::vector;

namespace anti {

/// Whitespace characters
const char WHITESPACE[] = " \t\r\n\f\v";
/*
ColorMap &ColorMap::operator=(const ColorMap &cmap)
{
   if(this!=&cmap) {
      copy_params(cmap);
      map_sz = cmap.map_sz;
      max_eff_map_sz = cmap.max_eff_map_sz;
   }
   return *this;
}
*/

//-----------------------------------------------------------------------
Status ColorMap::init(const char *params)
{
  wrap = 0;
  step = 1;
  shift = 0;

  if (!params) // treat null pointer like the empty string
    return Status::ok();

  // Make a copy for parsing
  string params_str = params;
  char *prms = &params_str[0];

  Status stat;

  char *p = strchr(prms, '%');
  if (p) {
    if (*(p + 1) == '\0')
      wrap = -1;
    else if (!(stat = read_int(p + 1, &wrap)))
      return Status::error(msg_str("wrap: %s", stat.c_msg()));
    *p = '\0';
  }

  p = strchr(prms, '*');
  if (p) {
    if (!(stat = read_int(p + 1, &step)))
      return Status::error(msg_str("step size: %s", stat.c_msg()));
    *p = '\0';
  }

  p = strchr(prms, '+');
  if (p) {
    if (!(stat = read_int(p + 1, &shift)))
      return Status::error(msg_str("shift size: %s", stat.c_msg()));
    *p = '\0';
  }

  return Status::ok();
}

void ColorMap::copy_params(const ColorMap &cmap)
{
  shift = cmap.shift;
  step = cmap.step;
  wrap = cmap.wrap;
}

Status ColorMap::read_params(const char *params)
{
  ColorMap cm;
  Status stat;
  if (!(stat = cm.init(params)))
    return stat;

  copy_params(cm);
  return Status::ok();
}

Status ColorMap::init_strip(char *map_name)
{
  Status stat;
  if (!(stat = ColorMap::init(map_name)))
    return stat;
  size_t name_len = strcspn(map_name, "+*%");
  map_name[name_len] = '\0';
  return Status::ok();
}

void ColorMap::set_wrap(int wrp)
{
  if (wrp < 0)
    wrp = effective_size();
  wrap = (wrp < 0) ? 0 : wrp;
}

ColorMap *ColorMapMap::get_condensed() const
{
  auto *colmap = new ColorMapMap();
  if (!colmap)
    return nullptr;

  int idx = 0;
  map<int, Color>::const_iterator mi;
  for (mi = cmap.begin(); mi != cmap.end(); mi++)
    colmap->cmap[idx++] = (mi->second);

  return colmap;
}

//-----------------------------------------------------------------------

/// A colour map that remaps index numbers
class ColorMapRemap : public ColorMap {
public:
  /// Get a copy of the map
  /**\return a pointer to the dynamically allocated copy,
   *  which must be freed by the caller with \c delete, 0 indicates
   *  that the clone failed. */
  ColorMap *clone() const { return new ColorMapRemap(*this); }

  /// Get the colour value for an index number.
  /**\param idx the index.
   * \return The colour. */
  virtual Color get_col(int idx) const { return get_effective_index(idx); }
};

class UniformRandom : public Random {
public:
  UniformRandom() { time_seed(); }
  using result_type = long unsigned int;
  static constexpr result_type min() { return 0; }
  static constexpr result_type max() { return 2147483647; } // 2*31 - 1
  result_type operator()() { return ranlui(); }
};

/// A colour map that maps index numbers to shuffled packs of numbers
class ColorMapDeal : public ColorMap {
private:
  int map_sz;
  int pack_sz;
  std::vector<int> map_vals;
  UniformRandom urnd;

public:
  /// Initialise from a string
  /**\param map_name the map name.
   * \return status, evaluates to \c true if the map could be initialised
   *  (possibly with warnings) otherwise \c false. */
  virtual Status init(const char *map_name);

  /// Get a copy of the map
  /**\return a pointer to the dynamically allocated copy,
   *  which must be freed by the caller with \c delete, 0 indicates
   *  that the clone failed. */
  ColorMap *clone() const { return new ColorMapDeal(*this); }

  /// The effective size of the map
  /** The effective size of a map is one greater than the highest
   *  index number in the map. It is the size of the smallest
   *  map (sequential, starting at 0) that will include all
   *  the entries of the map.
   * \return The effective size */
  virtual int effective_size() const { return map_vals.size(); };

  /// Get the colour value for an index number.
  /**\param idx the index.
   * \return The colour. */
  virtual Color get_col(int idx) const;

  /// Shuffle the mapping
  void shuffle();
};

Status ColorMapDeal::init(const char *map_name)
{
  // Make a copy for parsing
  string name_str(map_name);
  char *name = &name_str[0];

  Status stat = init_strip(name);
  if (stat.is_error())
    return stat;

  Split vals(name, "_", true);

  if (vals.size() > 2)
    return Status::error("map_name contains more than one '_'");

  // Get the map size
  map_sz = 256;
  if (*vals[0]) {
    if (!(stat = read_int(vals[0], &map_sz)))
      return Status::error(msg_str("map size: '%s' %s", vals[0], stat.c_msg()));
    if (map_sz < 0)
      return Status::error("map size: cannot be negative");
  }

  pack_sz = map_sz;
  if (vals.size() > 1 && *vals[1]) { // pack size given
    if (!(stat = read_int(vals[1], &pack_sz)))
      return Status::error(
          msg_str("range size: '%s' %s", vals[1], stat.c_msg()));
    if (pack_sz < 1)
      return Status::error("range size: cannot be less than 1");
  }

  shuffle();
  return Status::ok();
}

void ColorMapDeal::shuffle()
{
  if (!map_sz) {
    map_vals.clear();
    return;
  }
  int num_packs = map_sz / pack_sz;
  if (num_packs * pack_sz < map_sz)
    num_packs += 1;
  map_vals.resize(num_packs * pack_sz);
  for (int i = 0; i < num_packs; i++) {
    int off = i * pack_sz;
    for (int j = 0; j < pack_sz; j++)
      map_vals[off + j] = j;
    std::shuffle(map_vals.begin() + off, map_vals.begin() + off + pack_sz,
                 urnd);
  }
  map_vals.resize(map_sz);
}

Color ColorMapDeal::get_col(int idx) const
{
  int eff_idx = get_effective_index(idx);
  if (eff_idx < (int)map_vals.size())
    return Color(map_vals[eff_idx]);
  else
    return Color();
}

//-----------------------------------------------------------------------

/// A colour map that reverses index numbers
class ColorMapReverse : public ColorMap {
private:
  vector<int> switch_idxs = {0, 256};

public:
  /// Initialise from a string
  /**\param map_name the map name.
   * \return status, evaluates to \c true if the map could be initialised
   *  (possibly with warnings) otherwise \c false. */
  virtual Status init(const char *map_name);

  /// Get a copy of the map
  /**\return a pointer to the dynamically allocated copy,
   *  which must be freed by the caller with \c delete, 0 indicates
   *  that the clone failed. */
  ColorMap *clone() const { return new ColorMapReverse(*this); }

  /// Get the colour value for an index number.
  /**\param idx the index.
   * \return The colour. */
  virtual Color get_col(int idx) const;
};

Status ColorMapReverse::init(const char *map_name)
{
  // Make a copy for parsing
  string name_str(map_name);
  char *name = &name_str[0];

  Status stat = init_strip(name);
  if (stat.is_error())
    return stat;

  Split vals(name, ":", true);

  if (vals.size() > 2)
    return Status::error("map_name contains more than one ':'");

  // if no N
  // default
  switch_idxs = {0, 256};

  if (*vals[0]) {
    // Formats:
    // N - Reverse to switch N-1 and 0
    // I:J - Reverse to switch I and J
    vector<int> params;
    int param;
    for (size_t i = 0; i < vals.size(); i++) {
      if (vals.size() > i && *vals[i]) {
        if (!(stat = read_int(vals[i], &param))) {
          const vector<string> param_names = {"first", "second"};
          return Status::error(msg_str("%s parameter: '%s' %s",
                                       param_names[i].c_str(), vals[i],
                                       stat.c_msg()));
        }
        else
          params.push_back(param);
      }
    }

    if (vals.size() == 2) // format I:J
      switch_idxs = params;
    else if (vals.size() == 1) // format N
      switch_idxs = {0, params[0] - 1};
  }

  return Status::ok();
}

Color ColorMapReverse::get_col(int idx) const
{
  int eff_idx = get_effective_index(idx);
  int new_idx = switch_idxs[0] + switch_idxs[1] - eff_idx;
  return Color(new_idx);
}

//-----------------------------------------------------------------------

/// A colour map using a range
class ColorMapRange : public ColorMap {
private:
  int map_sz;

protected:
  std::vector<double> ranges[4];
  bool (Color::*set_func)(double, double, double, double);

public:
  /// Initialise from a string
  /**\param map_name the map name.
   * \return status, evaluates to \c true if the map could be initialised
   *  (possibly with warnings) otherwise \c false. */
  virtual Status init(const char *map_name);

  /// Get a copy of the map
  /**\return a pointer to the dynamically allocated copy,
   *  which must be freed by the caller with \c delete, 0 indicates
   *  that the clone failed. */
  ColorMap *clone() const { return new ColorMapRange(*this); }

  /// Set a range
  /**\param idx the index number of the component (0-3 for RGBA or HSVA)
   * \param range the range to set.
   * \return \c true if the range was valid, else \c false and the
   *  range was not changed. */
  virtual bool set_range(int idx, std::vector<double> range);

  /// Get the colour value for an index number.
  /**\param idx the index.
   * \return The colour. */
  virtual Color get_col(int idx) const;

  /// Set the map size
  /**\param sz the number of entries in the map. */
  void set_map_sz(int sz) { map_sz = sz; }

  /// Get the map size
  /**\return the number of entries in the map. */
  int get_map_sz() const { return map_sz; }

  /// The effective size of the map
  /** The effective size of a map is one greater than the highest
   *  index number in the map. It is the size of the smallest
   *  map (sequential, starting at 0) that will include all
   *  the entries of the map.
   * \return The effective size */
  virtual int effective_size() const { return map_sz; };
};

static double interpolate(int num, int map_sz, vector<double> vals)
{
  int interval_sz = (int)vals.size() - 1;
  if (interval_sz < 1)
    return vals[0];

  double pos = (double)num / map_sz;
  int low_idx = (int)floor(interval_sz * pos);
  int high_idx = (int)ceil(interval_sz * pos);
  double low_frac = (double)low_idx / interval_sz;
  double high_frac = (double)high_idx / interval_sz;
  double frac =
      (high_idx == low_idx) ? 0 : (pos - low_frac) / (high_frac - low_frac);
  double val = vals[low_idx] + frac * (vals[high_idx] - vals[low_idx]);
  return val;
}

Color ColorMapRange::get_col(int idx) const
{
  Color col;
  int eff_idx = get_effective_index(idx);
  if (eff_idx < map_sz) {
    double comp[4];
    for (int i = 0; i < 4; i++)
      comp[i] = interpolate(get_effective_index(idx), map_sz, ranges[i]);
    (col.*set_func)(comp[0], comp[1], comp[2], comp[3]);
  }
  return col;
}

Status ColorMapRange::init(const char *map_name)
{
  // Make a copy for parsing
  string name_str(map_name);
  char *name = &name_str[0];

  const char *p = strchr(map_name, '_');
  if (p && *(p + 1) == '\0')
    return Status::error("map_name contains trailing '_'");

  Status stat = init_strip(name);
  if (stat.is_error())
    return stat;

  Split vals(name, "_");
  // for(unsigned int i=0; i<vals.size(); i++)
  //   fprintf(stderr, "vals[%d] = '%s'\n", i, vals[i]);

  if (vals.size() > 2)
    return Status::error("map_name contains more than one '_'");
  // Get the map size
  if (*map_name != '_') {
    if (vals.size()) {
      if (!(stat = read_int(vals[0], &map_sz)))
        return Status::error(msg_str("map size: %s", stat.c_msg()));
      if (map_sz < 0)
        return Status::error("map size: cannot be negative");
    }
  }

  if (get_wrap() == -1)
    set_wrap(effective_size());

  if (*map_name != '_' &&
      vals.size() < 2) // A size was given but no comp ranges
    return Status::ok();

  const char *vals_back = (vals.size()) ? vals[vals.size() - 1] : "";
  if (strpbrk(vals_back, "HhSsVv") && strpbrk(vals_back, "RrGgBb"))
    return Status::error(
        "component letters include both 'HhSsVv' "
        "and 'RrGgBb', can only include letters from one set\n");

  int rng_len = strlen(vals_back);
  vector<char> rngs(rng_len + 1, 0);
  char *q = rngs.data();
  int cur_idx = -1;
  char cur_comp = 'X';
  for (const char *p = vals_back; p - vals_back < rng_len + 1; p++) {
    // fprintf(stderr, "*p = %c\n", *p);
    if (strchr("HhSsVvAaRrGgBb", *p) || cur_idx < 0 || *p == '\0') {
      *q = '\0';
      if (cur_idx >= 0) {
        if ((stat = read_double_list(rngs.data(), ranges[cur_idx], 0, ":"))) {
          if (ranges[cur_idx].size() == 0)
            return Status::error(msg_str("component letter '%c' "
                                         "isn't followed by any values",
                                         cur_comp));

          for (double &j : ranges[cur_idx]) {
            if (j < 0)
              return Status::error(msg_str("component letter '%c' "
                                           "contains a negative value",
                                           cur_comp));
            if (cur_comp == 'h')
              j /= 360; // h is hue in range 0-360
          }
        }
        else {
          return Status::error(
              msg_str("component letter '%c': %s", cur_comp, stat.c_msg()));
        }
      }
      if (strchr("HhRr", *p))
        cur_idx = 0;
      else if (strchr("SsGg", *p))
        cur_idx = 1;
      else if (strchr("VvBb", *p))
        cur_idx = 2;
      else if (strchr("Aa", *p))
        cur_idx = 3;
      else
        return Status::error(msg_str("invalid component letter '%c'", *p));
      cur_comp = *p;
      q = rngs.data();
    }
    else if (!(isdigit(*p) || *p == '.' || *p == ':')) {
      return Status::error(msg_str("invalid component letter '%c'", *p));
    }
    else if (!isspace(*p)) {
      *q++ = *p;
    }
  }

  // for(int i=0; i<4; i++)
  //   for(unsigned int j=0; j<ranges[i].size(); j++)
  //      fprintf(stderr, "ranges[%d][%u] = %g\n", i, j, ranges[i][j]);

  return Status::ok();
}

bool ColorMapRange::set_range(int idx, vector<double> range)
{
  if (idx < 0 || idx > 3 || !range.size())
    return false;
  ranges[idx] = range;
  return true;
}

//-----------------------------------------------------------------------

/// A colour map using values in an HSVA range
class ColorMapRangeHsv : public ColorMapRange {
public:
  /// Initialise from a string
  /**\param map_name the map name.
   * \return status, evaluates to \c true if the map could be initialised
   *  (possibly with warnings) otherwise \c false. */
  virtual Status init(const char *map_name);

  /// Get a copy of the map
  /**\return a pointer to the dynamically allocated copy,
   *  which must be freed by the caller with \c delete, 0 indicates
   *  that the clone failed. */
  ColorMap *clone() const { return new ColorMapRangeHsv(*this); }
};

Status ColorMapRangeHsv::init(const char *map_name)
{
  set_func = &Color::set_hsva;
  ranges[0].push_back(0);
  ranges[0].push_back(1);
  ranges[1].push_back(0.9);
  ranges[2].push_back(0.9);
  ranges[3].push_back(1);
  set_map_sz(256);

  return ColorMapRange::init(map_name);
}

static double rand_in_range(vector<double> rng, int seed)
{
  if (rng.size() > 1) {
    Random rnd((seed + 130));
    rnd.seedi((rnd.ranlui() * rnd.ranlui() * rnd.ranlui()) & 0xFFFFFFFF);
    const double val =
        fmod(rng[0] + (rng[1] - rng[0]) * rnd.ranf(), 1 + anti::epsilon);
    return val;
  }
  else
    return rng[0];
}

//-----------------------------------------------------------------------

/// A colour map using values in an RGBA range
class ColorMapRangeRgb : public ColorMapRange {
public:
  /// Initialise from a string
  /**\param map_name the map name.
   * \return status, evaluates to \c true if the map could be initialised
   *  (possibly with warnings) otherwise \c false. */
  virtual Status init(const char *map_name);

  /// Get a copy of the map
  /**\return a pointer to the dynamically allocated copy,
   *  which must be freed by the caller with \c delete, 0 indicates
   *  that the clone failed. */
  ColorMap *clone() const { return new ColorMapRangeRgb(*this); }
};

Status ColorMapRangeRgb::init(const char *map_name)
{
  set_func = &Color::set_rgba;
  ranges[0].push_back(0.3);
  ranges[0].push_back(1);
  ranges[1].push_back(0.3);
  ranges[1].push_back(1);
  ranges[2].push_back(0.3);
  ranges[2].push_back(1);
  ranges[3].push_back(1);
  set_map_sz(256);

  return ColorMapRange::init(map_name);
}

//-----------------------------------------------------------------------

/// A colour map using random values in a range
class ColorMapRangeRand : public ColorMapRange {
public:
  // Check the component ranges have 1 or 2 values
  Status check_range_sizes(const char *components);

  /// Get the colour value for an index number.
  /**\param idx the index.
   * \return The colour. */
  virtual Color get_col(int idx) const;

  /// Get a copy of the map
  /**\return a pointer to the dynamically allocated copy,
   *  which must be freed by the caller with \c delete, 0 indicates
   *  that the clone failed. */
  ColorMap *clone() const { return new ColorMapRangeRand(*this); }
};

Status ColorMapRangeRand::check_range_sizes(const char *components)
{
  for (int i = 0; i < 4; i++)
    if (ranges[i].size() == 0 || ranges[i].size() > 2)
      return Status::error(
          msg_str("component %c or %c: %s values given", components[i],
                  tolower(components[i]),
                  (ranges[i].size() == 0) ? "no" : "more than two"));
  return Status::ok();
}

Color ColorMapRangeRand::get_col(int idx) const
{
  Color col;
  idx = get_effective_index(idx);
  if (get_wrap() || idx < get_map_sz())
    (col.*set_func)(
        rand_in_range(ranges[0], idx * 1), rand_in_range(ranges[1], idx * 2),
        rand_in_range(ranges[2], idx * 3), rand_in_range(ranges[3], idx * 4));
  return col;
}

//-----------------------------------------------------------------------

/// A colour map using random values in an HSVA range
class ColorMapRangeRandHsv : public ColorMapRangeRand {
public:
  /// Initialise from a string
  /**\param map_name the map name.
   * \return status, evaluates to \c true if the map could be initialised
   *  (possibly with warnings) otherwise \c false. */
  virtual Status init(const char *map_name);

  /// Get a copy of the map
  /**\return a pointer to the dynamically allocated copy,
   *  which must be freed by the caller with \c delete, 0 indicates
   *  that the clone failed. */
  ColorMap *clone() const { return new ColorMapRangeRandHsv(*this); }
};

Status ColorMapRangeRandHsv::init(const char *map_name)
{
  set_func = &Color::set_hsva;

  ranges[0] = {0, 1};   // H
  ranges[1] = {0.7, 1}; // S
  ranges[2] = {0.7, 1}; // V
  ranges[3] = {1};      // A

  set_map_sz(max_map_sz);

  Status stat = ColorMapRange::init(map_name);
  if (stat)
    stat = check_range_sizes("HSVA");

  return stat;
}

//-----------------------------------------------------------------------

/// A colour map using random values in an RGBA range
class ColorMapRangeRandRgb : public ColorMapRangeRand {
public:
  /// Initialise from a string
  /**\param map_name the map name.
   * \return status, evaluates to \c true if the map could be initialised
   *  (possibly with warnings) otherwise \c false. */
  virtual Status init(const char *map_name);

  /// Get a copy of the map
  /**\return a pointer to the dynamically allocated copy,
   *  which must be freed by the caller with \c delete, 0 indicates
   *  that the clone failed. */
  ColorMap *clone() const { return new ColorMapRangeRandRgb(*this); }
};

Status ColorMapRangeRandRgb::init(const char *map_name)
{
  set_func = &Color::set_rgba;

  ranges[0] = {0.3, 1}; // R
  ranges[1] = {0.3, 1}; // G
  ranges[2] = {0.4, 1}; // B
  ranges[3] = {1};      // A

  set_map_sz(max_map_sz);

  Status stat = ColorMapRange::init(map_name);
  if (stat)
    stat = check_range_sizes("RGBA");

  return stat;
}

//-----------------------------------------------------------------------

/// A colour map with a good spread of colours
class ColorMapSpread : public ColorMapRange {
public:
  /// Initialise from a string
  /**\param map_name the map name.
   * \return status, evaluates to \c true if the map could be initialised
   *  (possibly with warnings) otherwise \c false. */
  virtual Status init(const char *map_name);

  /// Get a copy of the map
  /**\return a pointer to the dynamically allocated copy,
   *  which must be freed by the caller with \c delete, 0 indicates
   *  that the clone failed. */
  ColorMap *clone() const { return new ColorMapSpread(*this); }

  /// Get the colour value for an index number.
  /**\param idx the index.
   * \return The colour. */
  virtual Color get_col(int idx) const;
};

Status ColorMapSpread::init(const char *map_name)
{
  if (strchr(map_name, '_'))
    return Status::error(
        "spread map cannot contain '_' (does not take range specifiers)");

  set_map_sz(max_map_sz);

  return ColorMapRange::init(map_name);
}

Color ColorMapSpread::get_col(int idx) const
{
  int eff_idx = get_effective_index(idx);
  if (eff_idx >= get_map_sz())
    return Color();

  int num_entries = 1024;
  int num_intervals = 4;
  // int step_by = 29+196;
  // int step_by = 37+128;
  int step_by = 53;
  int entries_per_interval = num_entries / num_intervals;

  eff_idx = ((long)eff_idx * step_by) % num_entries;
  int interval = eff_idx / entries_per_interval;
  float H = (float)eff_idx / entries_per_interval - interval;
  float S = 0.0;
  float V = 0.0;
  switch (interval) {
  case 0:
    S = 0.9;
    V = 1.0;
    break;
  case 1:
    S = 0.5;
    V = 1.0;
    break;
  case 2:
    S = 0.9;
    V = 0.5;
    break;
  case 3:
    S = 0.5;
    V = 0.6;
    break;
  }

  Color col;
  col.set_hsva(H, S, V);

  return col;
}

//-----------------------------------------------------------------------

// ColorMapMap

static Status parse_gimp_file(FILE *cfile, map<int, Color> *cmap)
{
  const int line_size = 1024;
  char line[line_size];
  char buf[line_size];

  int stage = 0;

  int idx_no = 0;
  int line_no = 0;
  while (fgets(line, line_size, cfile)) {
    line_no++;

    // ignore comments
    char *first_hash = strchr(line, '#');
    if (first_hash)
      *first_hash = '\0';

    // skip blank lines
    if (sscanf(line, " %s", buf) == EOF)
      continue;

    if (stage == 0) {
      stage++;
      // ignore header line
      if (strncasecmp(line, "GIMP Palette", 12) == 0)
        continue;
    }
    if (stage == 1) {
      stage++;
      // ignore header line
      if (strncasecmp(line, "Name", 4) == 0)
        continue;
    }

    if (stage == 2) {
      stage++;
      // ignore header line
      if (strncasecmp(line, "Columns", 7) == 0)
        continue;
    }

    // stage == 3
    char *r = strtok(line, WHITESPACE);
    char *g = (r) ? strtok(nullptr, WHITESPACE) : nullptr;
    char *b = (g) ? strtok(nullptr, WHITESPACE) : nullptr;
    // char *name = (b) ? strtok(NULL, WHITESPACE) : 0;

    if (!b)
      return Status::error(msg_str(
          "gimp colour map: line %d: not enough colour values", line_no));

    sprintf(buf, "%s %s %s", r, g, b);
    Color col;
    col.read(buf);
    if (!col.is_set())
      return Status::error(
          msg_str("gimp colour map: line %d: invalid colour '%.*s'", line_no,
                  150, buf));

    // col.get_val().dump();
    (*cmap)[idx_no++] = col;
  }

  return Status::ok();
}

static Status parse_file(FILE *cfile, map<int, Color> *cmap)
{
  string messages;
  const int line_size = 1024;
  char line[line_size];

  int line_no = 0;
  int next_idx = 0; // index to use if no map index is given with a colour
  while (fgets(line, line_size, cfile)) {
    line_no++;

    // copy the map entry string
    string entry_cpy(line);
    char *entry = &entry_cpy[0];

    // ignore comments
    char *first_hash = strchr(entry, '#');
    if (first_hash)
      *first_hash = '\0';

    // skip blank lines
    char c;
    if (sscanf(entry, " %c", &c) == EOF)
      continue;

    if (!cmap->size() && strncasecmp(entry, "GIMP Palette", 12) == 0) {
      rewind(cfile);
      return parse_gimp_file(cfile, cmap);
    }

    char *col_pos = entry;
    char *eq_pos = strchr(entry, '=');
    if (eq_pos) {
      col_pos = eq_pos + 1;
      if (strchr(col_pos, '='))
        return Status::error(
            msg_str("colour map: line %d: more than one =, '%.*s'", line_no,
                    150, line));

      *eq_pos = '\0';
      if (!read_int(entry, &next_idx) || next_idx < 0)
        return Status::error(
            msg_str("colour map: line %d: invalid index number, '%.*s'",
                    line_no, 150, entry));
    }

    Color col;
    col.read(col_pos);
    if (!col.is_set())
      return Status::error(
          msg_str("colour map: line %d: invalid colour, '%.*s'", line_no, 150,
                  col_pos));

    if (cmap->find(next_idx) != cmap->end())
      messages += msg_str(
          "colour map: line %d: mapping for index %d is being overwritten\n",
          line_no, next_idx);

    (*cmap)[next_idx++] = col;
  }

  return (messages.empty()) ? Status::ok() : Status::warning(messages);
}

static Status parse_map_from_line(const char *line, map<int, Color> *cmap)
{
  string messages;
  Split entries(line, ":");
  int next_idx = 0; // index to use if no map index is given with a colour
  for (unsigned int i = 0; i < entries.size(); i++) {
    // copy the map entry string
    string entry_cpy(entries[i]);
    char *entry = &entry_cpy[0];

    // ignore comments
    char *first_hash = strchr(entry, '#');
    if (first_hash)
      *first_hash = '\0';

    // skip blank lines
    char c;
    if (sscanf(entry, " %c", &c) == EOF)
      continue;

    char *col_pos = entry;
    char *eq_pos = strchr(entry, '=');
    if (eq_pos) {
      col_pos = eq_pos + 1;
      if (strchr(col_pos, '='))
        return Status::error(
            msg_str("entry %d: more than one =, '%s'", i + 1, entries[i]));

      *eq_pos = '\0';
      if (!read_int(entry, &next_idx) || next_idx < 0)
        return Status::error(msg_str("entry %d: invalid index number, '%.*s'",
                                     i + 1, 150, entry));
    }

    // Allow '' as a number separator
    for (char *p = col_pos; *p; p++)
      if (*p == '/')
        *p = ' ';

    Color col;
    col.read(col_pos);
    if (!col.is_set())
      return Status::error(
          msg_str("entry %d: invalid colour, '%.*s'", i + 1, 150, col_pos));

    if (cmap->find(next_idx) != cmap->end())
      messages += msg_str("entry %d: mapping for index %d is being overwritten",
                          i, next_idx);

    (*cmap)[next_idx++] = col;
  }

  return (messages.empty()) ? Status::ok() : Status::warning(messages);
}

Color ColorMapMap::get_col(int idx) const
{
  auto mi_idx = cmap.find(get_effective_index(idx));
  if (mi_idx != cmap.end()) // index is in colour map
    return mi_idx->second;
  else
    return Color();
}

void ColorMapMap::set_col(int idx, Color col)
{
  if (col.is_set())
    cmap[idx] = col;
  else
    cmap.erase(idx);
}

int ColorMapMap::effective_size() const
{
  return size() ? cmap.rbegin()->first + 1 : 0;
}

Status ColorMapMap::init(const char *map_name)
{
  // Make a copy for parsing
  string name_str(map_name);
  char *name = &name_str[0];

  Status stat = init_strip(name);
  if (stat.is_error())
    return stat;

  string alt_name;
  FILE *cfile = open_sup_file(name, "/col_maps/", &alt_name);
  if (cfile)
    stat = parse_file(cfile, &cmap);
  else
    stat.set_error("map file not found");

  if (get_wrap() == -1)
    set_wrap(effective_size());

  if (cfile)
    fclose(cfile);

  return stat;
}

Status ColorMapMap::init_from_line(const char *map_line)
{
  string map_line_cpy(map_line);
  char *name = &map_line_cpy[0];

  Status stat = init_strip(name);
  if (stat)
    stat = parse_map_from_line(name, &cmap);

  if (stat && get_wrap() == -1)
    set_wrap(effective_size());

  return stat;
}

void ColorMapMap::read_named_colors()
{
  Color col;
  for (int i = 0; *named_colors[i].name; i++) {
    col.set_rgba(named_colors[i].r, named_colors[i].g, named_colors[i].b);
    cmap[i] = col;
  }
}

//-----------------------------------------------------------------------

// ColorMapMulti

Status ColorMapMulti::init(const char *map_name)
{
  ColorMap::init(nullptr);
  map_sz = -1;
  while (cmaps.size())
    del_cmap();

  // Make a copy for parsing
  string names_str(map_name);
  char *names = &names_str[0];

  Status stat;
  Split parts(names, ",");
  int parts_sz = parts.size();
  for (int i = 0; i < parts_sz; i++) {
    ColorMap *col_map = colormap_from_name(parts[i], &stat);
    string msg = stat.msg();
    if (col_map) {
      add_cmap(col_map);
      if (msg != "")
        stat.set_warning(msg_str("map '%s': %s", parts[i], msg.c_str()));
    }
    else {
      return stat.set_error(msg_str("map '%s': %s", parts[i], msg.c_str()));
    }
  }

  return stat;
}

ColorMapMulti::~ColorMapMulti()
{
  for (auto &cmap : cmaps)
    delete cmap;
}

ColorMapMulti::ColorMapMulti(const ColorMapMulti &cmap) : ColorMap(cmap)
{
  copy_params(cmap);
  map_sz = cmap.map_sz;
  max_eff_map_sz = cmap.max_eff_map_sz;
  for (auto i : cmap.cmaps)
    add_cmap(i->clone());
}

ColorMapMulti &ColorMapMulti::operator=(const ColorMapMulti &cmap)
{
  if (this != &cmap) {
    copy_params(cmap);
    map_sz = cmap.map_sz;
    max_eff_map_sz = cmap.max_eff_map_sz;
    while (cmaps.size())
      del_cmap();
    for (auto i : cmap.cmaps)
      add_cmap(i->clone());
  }
  return *this;
}

void ColorMapMulti::set_max_eff_map_sz()
{
  max_eff_map_sz = 0;
  for (auto &cmap : cmaps) {
    if (cmap->effective_size() > max_eff_map_sz)
      max_eff_map_sz = cmap->effective_size();
    if (cmap->get_wrap() != 0)
      max_eff_map_sz = max_map_sz;
  }
}

void ColorMapMulti::add_cmap(ColorMap *col_map, unsigned int pos)
{
  vector<ColorMap *>::iterator mi;
  if (pos >= cmaps.size())
    mi = cmaps.end();
  else
    mi = cmaps.begin() + pos;
  cmaps.insert(mi, col_map);
  if (col_map->effective_size() > max_eff_map_sz)
    max_eff_map_sz = col_map->effective_size();
  if (col_map->get_wrap() != 0)
    max_eff_map_sz = max_map_sz;
}

void ColorMapMulti::del_cmap(unsigned int pos)
{
  if (cmaps.size()) {
    vector<ColorMap *>::iterator mi;
    if (pos >= cmaps.size())
      mi = cmaps.end() - 1;
    else
      mi = cmaps.begin() + pos;
    delete *mi;
    cmaps.erase(mi);
    set_max_eff_map_sz();
  }
}

Color ColorMapMulti::get_col(int idx) const
{
  Color col;
  int cur_idx = get_effective_index(idx);
  if (cur_idx < effective_size()) {
    for (auto cmap : cmaps) {
      col = cmap->get_col(cur_idx);
      if (col.is_value())
        break;
      if (col.is_index())
        cur_idx = col.get_index();
    }
  }

  return col.is_set() ? col : Color(idx);
}

//-----------------------------------------------------------------------

static ColorMap *colormap_from_name_generated(const char *map_name,
                                              Status *stat)
{
  stat->set_ok();

  // Make a copy for parsing
  string name_str(map_name);
  char *name = &name_str[0];
  size_t name_len = strspn(name, "abcdefghijklmnopqrstuvwxyz");
  name[name_len] = '\0';

  size_t num_len = strspn(map_name + name_len, "0123456789");
  bool extra_chars = *(map_name + name_len + num_len) != '\0';
  int map_size = -1;
  if (num_len)
    map_size = atoi(string(map_name + name_len, num_len).c_str());

  // fprintf(stderr, "map_size=%d, extra_chars=%d (%c)\n", map_size,
  // extra_chars, *(map_name+name_len+num_len));

  ColorMap *cmap = nullptr;

  if (*map_name == '\0' || strcmp(name, "null") == 0) { // A null map
    cmap = new ColorMapMap();
  }

  else if (strcmp(name, "rnd") == 0 || strcmp(name, "rand") == 0 ||
           strcmp(name, "random") == 0) {
    if (!strpbrk(map_name + name_len + num_len, "RrGgBb"))
      cmap = new ColorMapRangeRandHsv();
    else
      cmap = new ColorMapRangeRandRgb();
    if (cmap && !(*stat = cmap->init(map_name + name_len))) {
      delete cmap;
      cmap = nullptr;
    }
  }

  else if (strcmp(name, "rng") == 0 || strcmp(name, "range") == 0) {
    if (!strpbrk(map_name + name_len + num_len, "RrGgBb"))
      cmap = new ColorMapRangeHsv();
    else
      cmap = new ColorMapRangeRgb();
    if (cmap && !(*stat = cmap->init(map_name + name_len))) {
      delete cmap;
      cmap = nullptr;
    }
  }

  else if (strcmp(name, "spread") == 0) {
    cmap = new ColorMapSpread();
    if (cmap && !(*stat = cmap->init(map_name + name_len))) {
      delete cmap;
      cmap = nullptr;
    }
  }

  else if (strcmp(name, "grey") == 0 || strcmp(name, "greyw") == 0 ||
           strcmp(name, "gray") == 0 || strcmp(name, "grayw") == 0) {
    if (strspn(map_name + name_len, "0123456789-+*%") ==
        strlen(map_name + name_len)) {
      bool wrp = (name[strlen(name) - 1] == 'w');
      size_t num_dgts = strspn(map_name + name_len, "0123456789");

      string grey_spec_str(map_name + name_len, num_dgts);
      grey_spec_str += msg_str("_H0S0V0:1%s%s", wrp ? ":0" : "",
                               map_name + name_len + num_dgts);
      cmap = new ColorMapRangeHsv();
      if (cmap && !(*stat = cmap->init(grey_spec_str.c_str()))) {
        delete cmap;
        cmap = nullptr;
      }
    }
  }

  else if (strcmp(name, "rainbow") == 0) {
    double fract = 0; // light/dark value

    size_t to_params_len = strcspn(map_name, "+*%"); // length to +*%
    // value string follows base map name and map size, and runs to +*%
    string value_str(map_name + name_len + num_len,
                     to_params_len - (name_len + num_len));
    if (!value_str.empty()) { // there are _value characters to process
      if (value_str[0] != '_')
        stat->set_error(msg_str("rainbow map has '%c' instead of underscore",
                                value_str[0]));
      else if (value_str.size() == 1)
        stat->set_error("rainbow map needs value after underscore");
      else { // there is a value after '_'
        if ((*stat = read_double(&value_str[1], &fract))) { // true if number
          if (fabs(fract) > 1) { // number was read, check value
            stat->set_error(msg_str(
                "rainbow map value '%g', must be between -1 and 1", fract));
          }
        }
      }

      if (stat->is_error())
        return cmap;
    }

    double x = 0;
    double y = 1;
    if (fract > 0)
      x = fract;
    else
      y = fract + 1;

    double o = 0.5; // orange
    if (x > 0)
      o = (x + 1) / 2;
    else if (y < 1)
      o = y / 2;

    string rainbow_str;
    rainbow_str = msg_str("_R%g:%g:%g:%g:%g:%g:%g:%gG%g:%g:%g:%g:%g:%g:%g:%gB%"
                          "g:%g:%g:%g:%g:%g:%g:%g",
                          y, y, y, x, x, x, y, y, x, o, y, y, y, x, x, x, x, x,
                          x, x, y, y, y, x);

    cmap = new ColorMapRangeRgb();
    string map_size(map_name + name_len, num_len);
    const char *final_params = map_name + to_params_len;
    if (cmap && !(*stat = cmap->init(
                      (map_size + rainbow_str + final_params).c_str()))) {
      delete cmap;
      cmap = nullptr;
    }
  }

  else if (strcmp(name, "uniform") == 0) {

    auto *multi = new ColorMapMulti;
    auto *overrides = new ColorMapMap;
    ColorMap *spread_map = colormap_from_name("spread+53*12");

    if (multi && overrides && spread_map &&
        (*stat = multi->read_params(map_name))) {

      overrides->set_col(60, Color(0.9, 0.45, 0.0)); // triangle
      overrides->set_col(36, Color(0.7, 0.1, 0.2));  // pentagram
      multi->add_cmap(overrides);
      multi->add_cmap(spread_map);
      if (map_size >= 0)
        multi->set_map_sz(map_size);
      cmap = multi;
    }
    else {
      if (extra_chars)
        stat->set_error(msg_str("uniform map: trailing characters '%s'",
                                map_name + name_len + num_len));
      delete multi;
      delete overrides;
      delete spread_map;
      cmap = nullptr;
    }
  }

  // RK - for blending colors, make each two colors blend to different color
  // than adjacent ones. uses nicer color values for green and orange
  else if (strcmp(name, "compound") == 0) {
    auto *multi = new ColorMapMulti();
    auto *overrides = new ColorMapMap;
    ColorMap *spread_map = colormap_from_name("spread+2");
    if (multi && overrides && spread_map &&
        (*stat = multi->read_params(map_name))) {
      overrides->set_col(0, Color(255, 255, 0)); // yellow
      overrides->set_col(1, Color(255, 0, 0));   // red
      overrides->set_col(2, Color(0, 100, 0));   // darkgreen (X11)
      overrides->set_col(3, Color(0, 0, 255));   // blue
      overrides->set_col(4, Color(255, 0, 255)); // magnenta
      overrides->set_col(5, Color(0, 255, 255)); // cyan
      overrides->set_col(6, Color(255, 127, 0)); // darkorange1 (X11)
      multi->add_cmap(overrides);
      multi->add_cmap(spread_map);
      if (map_size >= 0)
        multi->set_map_sz(map_size);
      cmap = multi;
    }
    else {
      delete multi;
      delete overrides;
      delete spread_map;
      cmap = nullptr;
    }
  }

  // RK - a map used in RK's programs. Here for external use.
  // for using with coloring n faces shift as 'colorful+-3'
  // use nicer color values for green and orange
  else if (strcmp(name, "colorful") == 0) {
    auto *multi = new ColorMapMulti();
    auto *overrides = new ColorMapMap;
    ColorMap *spread_map = colormap_from_name("spread");
    if (multi && overrides && spread_map &&
        (*stat = multi->read_params(map_name))) {
      overrides->set_col(0, Color(255, 0, 0));     // 3-sided red
      overrides->set_col(1, Color(255, 127, 0));   // 4-sided darkoranage1 (X11)
      overrides->set_col(2, Color(255, 255, 0));   // 5-sided yellow
      overrides->set_col(3, Color(0, 100, 0));     // 6-sided darkgreen (X11)
      overrides->set_col(4, Color(0, 255, 255));   // 7-sided cyan
      overrides->set_col(5, Color(0, 0, 255));     // 8-sided blue
      overrides->set_col(6, Color(255, 0, 255));   // 9-sided magenta
      overrides->set_col(7, Color(255, 255, 255)); // 10-sided white
      overrides->set_col(8, Color(127, 127, 127)); // 11-sided gray50
      overrides->set_col(9, Color(0, 0, 0));       // 12-sided black
      multi->add_cmap(overrides);
      multi->add_cmap(spread_map);
      if (map_size >= 0)
        multi->set_map_sz(map_size);
      cmap = multi;
    }
    else {
      delete multi;
      delete overrides;
      delete spread_map;
      cmap = nullptr;
    }
  }

  // RK - a map used in RK's programs. Here for external use.
  // source was from from George Hart's original conway applet
  // for using with coloring n faces shift as 'ghart+-3'
  else if (strcmp(name, "ghart") == 0) {
    auto *multi = new ColorMapMulti();
    auto *overrides = new ColorMapMap;
    ColorMap *spread_map = colormap_from_name("spread");
    if (multi && overrides && spread_map &&
        (*stat = multi->read_params(map_name))) {
      overrides->set_col(0, Color(0.9, 0.3, 0.3));   // 3-sided red
      overrides->set_col(1, Color(0.4, 0.4, 1.0));   // 4-sided blue
      overrides->set_col(2, Color(0.2, 0.9, 0.3));   // 5-sided green
      overrides->set_col(3, Color(0.9, 0.9, 0.2));   // 6-sided yellow
      overrides->set_col(4, Color(0.5, 0.25, 0.25)); // 7-sided brown
      overrides->set_col(5, Color(0.8, 0.2, 0.8));   // 8-sided magenta
      overrides->set_col(6, Color(0.5, 0.2, 0.8));   // 9-sided purple
      overrides->set_col(7, Color(0.1, 0.9, 0.9));   // 10-sided grue
      overrides->set_col(8, Color(0.5, 0.5, 0.5));   // 11-sided gray
      overrides->set_col(9, Color(1.0, 0.6, 0.1));   // 12-sided orange
      multi->add_cmap(overrides);
      multi->add_cmap(spread_map);
      if (map_size >= 0)
        multi->set_map_sz(map_size);
      cmap = multi;
    }
    else {
      delete multi;
      delete overrides;
      delete spread_map;
      cmap = nullptr;
    }
  }

  // RK - a map used in RK's programs. Here for external use.
  // colors matching the models in Kaplan's symmetro paper
  // colors from PDF document measured from screen
  else if (strcmp(name, "kaplan") == 0) {
    auto *multi = new ColorMapMulti();
    auto *overrides = new ColorMapMap;
    // spread map replaced with a repeating color
    // Color(145, 160, 100)); // 9-sided faces and higher
    ColorMap *repeat_map = colormap_from_name("map_x91a064");
    if (multi && overrides && repeat_map &&
        (*stat = multi->read_params(map_name))) {
      overrides->set_col(3, Color(141, 108, 47));  // 3-sided faces
      overrides->set_col(4, Color(105, 123, 96));  // 4-sided faces
      overrides->set_col(5, Color(119, 176, 77));  // 5-sided faces
      overrides->set_col(6, Color(115, 128, 44));  // 6-sided faces
      overrides->set_col(7, Color(171, 164, 115)); // 7-sided faces (blend 6,8)
      overrides->set_col(8, Color(227, 200, 186)); // 8-sided faces
      multi->add_cmap(overrides);
      repeat_map->set_wrap();
      multi->add_cmap(repeat_map);
      if (map_size >= 0)
        multi->set_map_sz(map_size);
      cmap = multi;
    }
    else {
      delete multi;
      delete overrides;
      delete repeat_map;
      cmap = nullptr;
    }
  }

  // RK - a map used in RK's programs. Here for external use.
  // used to color by faces by axis number in the symmetro program
  // or axes in off_color_radial program
  else if (strcmp(name, "axes") == 0) {
    auto *multi = new ColorMapMulti();
    auto *overrides = new ColorMapMap;
    if (multi && overrides && (*stat = multi->read_params(map_name))) {
      overrides->set_col(0, Color(255, 0, 0));   // axis1 red
      overrides->set_col(1, Color(0, 0, 255));   // axis2 blue
      overrides->set_col(2, Color(255, 255, 0)); // axis3 yellow
      overrides->set_col(3, Color(0, 100, 0));   // convex hull darkgreen
      multi->add_cmap(overrides);
      if (map_size >= 0)
        multi->set_map_sz(map_size);
      multi->set_wrap();
      cmap = multi;
    }
    else {
      delete multi;
      delete overrides;
      cmap = nullptr;
    }
  }

  // RK - a map used in RK's programs. Here for external use.
  // for circuits in n_icons
  else if (strcmp(name, "circuits") == 0) {
    auto *multi = new ColorMapMulti();
    auto *overrides = new ColorMapMap;
    // required for circuits so that the original coloring has enough colors
    ColorMap *spread_map = colormap_from_name("spread");
    if (multi && overrides && (*stat = multi->read_params(map_name))) {
      overrides->set_col(0, Color(255, 255, 255)); // continuous
      overrides->set_col(1, Color(64, 64, 64));    // discontinuous
      multi->add_cmap(overrides);
      multi->add_cmap(spread_map);
      if (map_size >= 0)
        multi->set_map_sz(map_size);
      cmap = multi;
    }
    else {
      delete multi;
      delete overrides;
      delete spread_map;
      cmap = nullptr;
    }
  }

  // RK - a map used in RK's programs. Here for external use.
  // used to color by convexity of faces and edges
  else if (strcmp(name, "convexity") == 0) {
    auto *multi = new ColorMapMulti();
    auto *overrides = new ColorMapMap;
    if (multi && overrides && (*stat = multi->read_params(map_name))) {
      overrides->set_col(0, Color(255, 255, 255)); // convex
      overrides->set_col(1, Color(127, 127, 127)); // coplanar
      overrides->set_col(2, Color(64, 64, 64));    // nonconvex
      multi->add_cmap(overrides);
      if (map_size >= 0)
        multi->set_map_sz(map_size);
      multi->set_wrap();
      cmap = multi;
    }
    else {
      delete multi;
      delete overrides;
      cmap = nullptr;
    }
  }

  else if (strcmp(name, "index") == 0 || strcmp(name, "remap") == 0) {
    auto *cmr = new ColorMapRemap;
    if (cmr) {
      if ((*stat = cmr->init(map_name)))
        cmap = cmr;
      else
        delete cmr;
    }
  }

  else if (strcmp(name, "deal") == 0) {
    auto *cmr = new ColorMapDeal;
    if (cmr) {
      if ((*stat = cmr->init(map_name + name_len)))
        cmap = cmr;
      else
        delete cmr;
    }
  }

  else if (strcmp(name, "reverse") == 0) {
    auto *cmr = new ColorMapReverse;
    if (cmr) {
      if ((*stat = cmr->init(map_name + name_len)))
        cmap = cmr;
      else
        delete cmr;
    }
  }

  else if (strcmp(name, "map") == 0 && map_name[name_len] == '_') {
    auto *cmm = new ColorMapMap;
    if (cmm) {
      if ((*stat = cmm->init_from_line(map_name + name_len + 1)))
        cmap = cmm;
      else
        delete cmm;
    }
  }
  else if (!stat->is_error()) {
    stat->set_error("name not found");
  }

  return cmap;
}

ColorMap *colormap_from_name(const char *map_name, Status *stat)
{
  ColorMap *cmap = nullptr;

  // Make a copy for parsing
  string name_str(map_name);
  char *name = &name_str[0];
  size_t name_len = strcspn(name, "+*%");
  name[name_len] = '\0';

  Status stat2;
  string alt_name;
  FILE *cfile = open_sup_file(name, "/col_maps/", &alt_name);
  if (alt_name != "") { // an alt name found before a file with the name
    cmap = new ColorMapMulti;
    if (cmap) {
      if (!(stat2 = cmap->init(alt_name.c_str())) || // init first
          !(stat2 = cmap->read_params(map_name))) {  // then read shft-stp-wrp
        if (stat)
          stat->set_error(msg_str("could not open colour map file '%s=%s': %s",
                                  map_name, alt_name.c_str(), stat2.c_msg()));
        delete cmap;
        cmap = nullptr;
      }
    }
  }
  else if (cfile) {
    cmap = new ColorMapMap();
    if (cmap && !(stat2 = cmap->init(map_name))) {
      delete cmap;
      cmap = nullptr;
    }
    if (stat)
      *stat = stat2;
  }
  else {
    cmap = colormap_from_name_generated(map_name, &stat2);
    if (stat) {
      if (cmap)
        *stat = stat2;
      else
        stat->set_error(msg_str("could not open colour map file \'%s\': %s",
                                map_name, stat2.c_msg()));
    }
  }

  return cmap;
}

//-----------------------------------------------------------------------
// ColorValuesToRangeHsva

static Status chk_range(vector<double> &v)
{
  if (v.size() == 1) {
    if (v[0] < -anti::epsilon || v[0] > 1 + anti::epsilon)
      return Status::error(
          msg_str("value, %g, is not in rage 0.0 to 1.0", v[0]));

    v.push_back(v[0]);
  }
  else if (v.size() == 2) {
    if (v[0] < -anti::epsilon || v[0] > 1 + anti::epsilon)
      return Status::error(
          msg_str("first value, %g, is not in rage 0.0 to 1.0", v[0]));

    if (v[1] < -anti::epsilon || v[1] > 1 + anti::epsilon)
      return Status::error(
          msg_str("second value, %g, is not in rage 0.0 to 1.0", v[1]));
  }
  else
    return Status::error(msg_str("range has %lu values, must have 1 or 2",
                                 (unsigned long)v.size()));

  if (v[0] > v[1] + anti::epsilon)
    v[1]++;

  return Status::ok();
}

Status ColorValuesToRangeHsva::add_range(int idx, const char *rngs)
{
  Status stat = read_double_list(rngs, ranges[idx], 2, ":");
  if (stat.is_error())
    return stat;
  if (!(stat = chk_range(ranges[idx])))
    return stat;

  return Status::ok();
}

ColorValuesToRangeHsva::ColorValuesToRangeHsva(const string &range_name,
                                               const Color &def_col)
{
  init(range_name, def_col);
}

Status ColorValuesToRangeHsva::init(const string &range_name,
                                    const Color &def_col)
{
  for (auto &range : ranges)
    range.clear();
  ranges[0].push_back(0);
  ranges[0].push_back(1);
  ranges[1].push_back(0);
  ranges[1].push_back(1);
  ranges[2].push_back(0);
  ranges[2].push_back(1);
  ranges[3].push_back(0);
  ranges[3].push_back(1);

  Status stat;
  int name_len = range_name.size();
  // Make a copy for parsing
  string name_str(range_name);
  char *name = &name_str[0];
  vector<char> rngs(name_len + 1, 0);
  char *q = rngs.data();
  int cur_idx = -1;
  char cur_comp = 'X';
  for (const char *p = name; p - name < name_len + 1; p++) {
    if (strchr("HhSsVvAa", *p) || cur_idx < 0 || *p == '\0') {
      *q = '\0';
      if (cur_idx >= 0 && !(stat = add_range(cur_idx, rngs.data())))
        return Status::error(
            msg_str("component '%c': %s", cur_comp, stat.c_msg()));
      if (strchr("Hh", *p))
        cur_idx = 0;
      else if (strchr("Ss", *p))
        cur_idx = 1;
      else if (strchr("Vv", *p))
        cur_idx = 2;
      else if (strchr("Aa", *p))
        cur_idx = 3;
      else
        return Status::error(msg_str("invalid component letter '%c'", *p));

      cur_comp = *p;
      q = rngs.data();
    }
    else if (!(isdigit(*p) || *p == '.' || *p == ':'))
      return Status::error(msg_str("invalid component letter '%c'",
                                   (cur_idx < 0) ? rngs[0] : *p));
    else if (!isspace(*p)) {
      *q++ = *p;
    }
  }

  set_default_color(def_col);

  return Status::ok();
}

static inline double fract(const vector<double> &rng, double frac)
{
  return fmod(rng[0] + (rng[1] - rng[0]) * frac, 1 + anti::epsilon);
}

Color ColorValuesToRangeHsva::get_col(Color col)
{
  Color c = col;
  if (col.is_visible_value()) {
    Vec4d hsva_orig = col.get_hsva();
    Vec4d hsva(fract(ranges[0], hsva_orig[0]), fract(ranges[1], hsva_orig[1]),
               fract(ranges[2], hsva_orig[2]), fract(ranges[3], hsva_orig[3]));
    c.set_hsva(hsva);
  }
  return c;
}

void ColorValuesToRangeHsva::apply(map<int, Color> &elem_cols)
{
  map<int, Color>::iterator mi;
  for (mi = elem_cols.begin(); mi != elem_cols.end(); ++mi)
    mi->second = get_col(mi->second);
}

void ColorValuesToRangeHsva::apply(Geometry &geom, int elem_type)
{
  if (default_color.is_set()) {
    int num_elems;
    if (elem_type == VERTS)
      num_elems = geom.verts().size();
    else if (elem_type == EDGES)
      num_elems = geom.edges().size();
    else // (elem_type == FACES)
      num_elems = geom.faces().size();

    for (int i = 0; i < num_elems; i++) {
      auto col = geom.colors(elem_type).get(i);
      if (!col.is_set())
        col = default_color;
      geom.colors(elem_type).set(i, get_col(col));
    }
  }
  else {
    apply(geom.colors(elem_type).get_properties());
  }
}

} // namespace anti
