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

/*!\file utils.h
   \brief utility routines for text operations, I/O conversions, etc
*/

#ifndef UTILS_H
#define UTILS_H

#include "geometry.h"

#include "trans3d.h"
#include "vec3d.h"
#include "vec4d.h"
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace anti {

/// Get the Basename
/** Get the name of the file from a path. This is taken as the last part
 *  of the path after any '%\' or '/', otherwise the path is the file name.
 * \param path the path to the file.
 * \return A pointer to the start of the file name. */
const char *basename2(const char *path);

/// Convert dots to underscores.
/**\param str the string to convert.
 * \return The converted string. */
std::string dots2underscores(std::string str);

/// Read a floating point number, or mathematical expression, from a string.
/** The string should only hold the floating point number or expression, but may
 *  have leading and trailing whitespace.
 * \param str the string holding the floating point number or expressions.
 * \param f used to return the floating point number.
 * \return status, evaluates to \c true if a valid floating point number
 *  was read, otherwise \c false.*/
Status read_double(const char *str, double *f);

/// Read a floating point number from a string.
/** The string should only hold the floating point number, but may
 *  have leading and trailing whitespace.
 * \param str the string holding the floating point number.
 * \param f used to return the floating point number.
 * \return status, evaluates to \c true if a valid floating point number
 *  was read, otherwise \c false.*/
Status read_double_noparse(const char *str, double *f);

/// Read an integer from a string.
/** The string should only hold the integer, but may
 *  have leading and trailing whitespace.
 * \param str the string holding the integer.
 * \param i used to return the integer.
 * \return status, evaluates to \c true if a valid integer
 *  was read, otherwise \c false.*/
Status read_int(const char *str, int *i);

/// Read a fraction from a string.
/**The string should only hold an integer, or two integers seperated by '/',
 * but may have leading and trailing whitespace.
 * \param frac_str the string holding the fraction.
 * \param num used to return the numerator.
 * \param denom used to return the denominator.
 * \return status, evaluates to \c true if a valid fraction
 * was read, otherwise \c false.*/
Status read_fraction(const char *str, int *num, int *denom);

/// Read floating point numbers, or mathematical expressions, from a list of
/// strings.
/** The strings should only hold the floating point number, but may
 *  have leading and trailing whitespace.
 * \param vals the strings holding the floating point numbers.
 * \param nums used to return the floating point numbers.
 * \return status, evaluates to \c true if only valid floating point numbers
 *  were read, otherwise \c false.*/
Status read_double_list(const std::vector<char *> &vals,
                        std::vector<double> &nums);

/// Read floating point numbers from a list of strings.
/** The strings should only hold the floating point number, but may
 *  have leading and trailing whitespace.
 * \param vals the strings holding the floating point numbers.
 * \param nums used to return the floating point numbers.
 * \return status, evaluates to \c true if only valid floating point numbers
 *  were read, otherwise \c false.*/
Status read_double_list_noparse(const std::vector<char *> &vals,
                                std::vector<double> &nums);

/// Read floating point numbers, or mathematical expressions, listed in a single
/// string.
/** The numbers in the string should be comma separated, and may
 *  have leading and trailing whitespace.
 * \param str the string holding the comma-separated floating point numbers.
 * \param nums used to return the floating point numbers.
 * \param len the maximum total of numbers that should be in \a str, or if
 *  it is \c 0 then there is no maximum.
 * \param sep the characters that can separate the numbers.
 * \return status, evaluates to \c true if only valid floating point numbers
 *  were read, , and no more than \a len (if \c len>0),otherwise \c false.*/
Status read_double_list(const char *str, std::vector<double> &nums, int len = 0,
                        const char *sep = ",");

/// Read floating point numbers listed in a single string.
/** The numbers in the string should be comma separated, and may
 *  have leading and trailing whitespace.
 * \param str the string holding the comma-separated floating point numbers.
 * \param nums used to return the floating point numbers.
 * \param len the maximum total of numbers that should be in \a str, or if
 *  it is \c 0 then there is no maximum.
 * \param sep the characters that can separate the numbers.
 * \return status, evaluates to \c true if only valid floating point numbers
 *  were read, , and no more than \a len (if \c len>0),otherwise \c false.*/
Status read_double_list_noparse(const char *str, std::vector<double> &nums,
                                int len = 0, const char *sep = ",");

/// Read integers from a list of strings.
/** The strings should only hold the integers, but may
 *  have leading and trailing whitespace.
 * \param vals the strings holding the integers.
 * \param nums used to return the integers.
 * \param is_index if true then the integers cannot be negative.
 * \return status, evaluates to \c true if only valid integers
 *  were read, otherwise \c false.*/
Status read_int_list(const std::vector<char *> &vals, std::vector<int> &nums,
                     bool is_index = false);

/// Read integers listed in a single string.
/** The numbers in the string should be comma separated, and may
 *  have leading and trailing whitespace.
 * \param str the string holding the comma-separated integers.
 * \param nums used to return the integers.
 * \param is_index if true then the integers cannot be negative.
 * \param len the maximum total of numbers that should be in \a str, or if
 *  it is \c 0 then there is no maximum.
 * \param sep the characters that can separate the numbers.
 * \return status, evaluates to \c true if only valid integers
 *  were read, , and no more than \a len (if \c len>0),otherwise \c false.*/
Status read_int_list(const char *str, std::vector<int> &nums,
                     bool is_index = false, int len = 0, const char *sep = ",");

/// Read index numbers listed in a single string.
/** The string consists of comma separated index number ranges, and may
 *  have leading and trailing whitespace. A number range is either a single
 *  number, or a sequential list indicated by two numbers seperated by '-'.
 *  If the numbers are not given they default to the first and last index
 *  number respectively.
 * \param str the string holding the comma-separated integers.
 * \param nums used to return the integers.
 * \param num_idxs index numbers with this value or higher are out of range
 * \param allow_extra allow out-of-range index numbers, prefixed by an x or X,
 *  which are indexed relative to num_idxs, i.e. X0 = num_idxs.
 * \return status, evaluates to \c true if only valid index numbers
 *  were read, otherwise \c false.*/
Status read_idx_list(const char *str, std::vector<int> &nums, int num_idxs,
                     bool allow_extra = false);

/// Read a line of arbitrary length
/** The caller is responsible for freeing the memory allocated to line
 *  after each read.
 * \param file the file stream to read from.
 * \param line where the line is returned
 * \return <ul>
 *    <li>\c 0 if the line was read correctly.
 *    <li>\c -1 if memory for \a line could not be allocated.
 *    <li>\c 1 if an final unterminated empty line was read.
 * </ul> */
int read_line(FILE *file, char **line);

/// Remove leading and trailing space, convert any whitespace to a single space
/**\param str the string to convert. */
void clear_extra_whitespace(std::string &str);

/// Convert to a normalised resource name
/** Remove leading and trailing space, convert any whitespace to a
 *  single space, make lowercase
 * \param to the string to convert (up to \c MSG_SZ-1 characters used.)
 * \param from the string to convert (length \c MSG_SZ.)
 * \return A pointer to the \a to string. */
char *to_resource_name(char *to, const char *from);

/// Open a support file
/** Tries to open a file by its name, then tries to open it in
 *  \c $ANTIPRISM_DATA/sub_dir, finally tries to open it in
 *  \c sub_dir in the installation data directory.
 * \param fname the name of the file to open.
 * \param subdir the data directory subdirectory to search in.
 * \param alt_name a name that is found in an alt_names.txt file before
 *  a file with the name is found.
 * \param where used to return the place that the file was found <ul>
 *  <li>\c 0 - locally
 *  <li>\c 1 - in the ANTIPRISM_DATA directory
 *  <li>\c 2 - in the installation data directory
 *  </ul>
 * \param fpath used to return the full path to the file that was found.
 * \return A pointer to the opened file stream. */
FILE *open_sup_file(const char *fname, const char *subdir,
                    std::string *alt_name = nullptr, int *where = nullptr,
                    std::string *fpath = nullptr);

/// Convert a C formated message string to a C++ string
/** Converts the first MSG_SZ-1 characters of the C format string
 * \param fmt the formatted string
 * \param ... the values for the format
 * \return The converted string. */
std::string msg_str(const char *fmt, ...);

/// Convert an integer to a string
/**\param i the integer.
 * \return The string. */
std::string itostr(int i);

/// Convert a floating point number to a string
/**\param f the floating point number.
 * \param sig_dgts the number of significant digits in the conversion,
 *  or if negative then the number of digits after the decimal point.
 * \return The string. */
std::string dtostr(double f, int sig_dgts = 17);

/// Convert a coordinate vector to a string
/**\param v the vector.
 * \param sep the separator between the numbers.
 * \param sig_dgts the number of significant digits in the conversion,
 *  or if negative then the number of digits after the decimal point.
 * \return The string. */
std::string vtostr(Vec3d v, const char *sep = ", ", int sig_dgts = 17);

/// Convert a coordinate vector to a string
/**\param v the vector.
 * \param sep the separator between the numbers.
 * \param sig_dgts the number of significant digits in the conversion,
 *  or if negative then the number of digits after the decimal point.
 * \return The string. */
std::string vtostr(Vec4d v, const char *sep = ", ", int sig_dgts = 17);

/// Class to split a line into delimited parts
class Split {
private:
  std::string data;
  std::vector<char *> parts;

public:
  /// Constructor
  Split() = default;

  /// Constructor
  /**\param line the line to split
   * \param delims the characters to use as delimiters, if \c nullptr then use
   *  whitespace characters.
   * \param strict if true then treat every delimiter as a separator, returning
   *  null strings between adjacent delimiters, always returning at least
   *  one part.*/
  Split(const char *line, const char *delims = nullptr, bool strict = false)
  {
    init(line, delims, strict);
  }

  /// Constructor
  /**\param line_str the line to split
   * \param delims the characters to use as delimiters, if \c nullptr then use
   *  whitespace characters.
   * \param strict if true then treat every delimiter as a separator, returning
   *  null strings between adjacent delimiters, always returning at least
   *  one part.*/
  Split(const std::string &line_str, const char *delims = nullptr,
        bool strict = false)
  {
    init(line_str, delims, strict);
  }

  /// Init
  /**\param line the line to split
   * \param delims the characters to use as delimiters, if \c nullptr then use
   *  whitespace characters.
   * \param strict if true then treat every delimiter as a separator, returning
   *  null strings between adjacent delimiters, always returning at least
   *  one part.
   * \return The number of parts. */
  int init(const char *line, const char *delims = nullptr, bool strict = false);

  /// Init
  /**\param line_str the line to split
   * \param delims the characters to use as delimiters, if \c nullptr then use
   *  whitespace characters.
   * \param strict if true then treat every delimiter as a separator, returning
   *  null strings between adjacent delimiters, always returning at least
   *  one part.
   * \return The number of parts. */
  int init(const std::string &line_str, const char *delims = nullptr,
           bool strict = false)
  {
    return init(line_str.c_str(), delims, strict);
  }

  /// Size, the number of parts
  /**\return The number of parts. */
  size_t size() const { return parts.size(); }

  /// Get a part
  /**\param idx the part index.
   * \return The part. */
  const char *operator[](int idx) const { return parts[idx]; }

  /// Get the parts
  /**\return The parts. */
  const std::vector<char *> &get_parts() const { return parts; }

  // Copy disallowed
  Split(const Split &) = delete;
  void operator=(const Split &) = delete;
};

char *strcpy_msg(char *dest, const char *src);

// inline function definitions

inline std::string itostr(int i)
{
  char buf[MSG_SZ];
  buf[MSG_SZ - 1] = 0;
  sprintf(buf, "%d", i);
  return buf;
}

inline std::string dtostr(double f, int sig_dgts)
{
  char buf[MSG_SZ];
  buf[MSG_SZ - 1] = '\0';
  if (sig_dgts > 0)
    snprintf(buf, MSG_SZ - 1, "%.*g", sig_dgts, f);
  else
    snprintf(buf, MSG_SZ - 1, "%.*f", -sig_dgts, f);
  return buf;
}

inline std::string vtostr(Vec3d v, const char *sep, int sig_dgts)
{
  char buf[MSG_SZ];
  buf[MSG_SZ - 1] = 0;
  if (sig_dgts > 0)
    snprintf(buf, MSG_SZ - 1, "%.*g%s%.*g%s%.*g", sig_dgts, v[0], sep, sig_dgts,
             v[1], sep, sig_dgts, v[2]);
  else
    snprintf(buf, MSG_SZ - 1, "%.*f%s%.*f%s%.*f", -sig_dgts, v[0], sep,
             -sig_dgts, v[1], sep, -sig_dgts, v[2]);
  return buf;
}

inline std::string vtostr(Vec4d v, const char *sep, int sig_dgts)
{
  char buf[MSG_SZ];
  buf[MSG_SZ - 1] = 0;
  if (sig_dgts > 0)
    snprintf(buf, MSG_SZ - 1, "%.*g%s%.*g%s%.*g%s%.*g", sig_dgts, v[0], sep,
             sig_dgts, v[1], sep, sig_dgts, v[2], sep, sig_dgts, v[3]);
  else
    snprintf(buf, MSG_SZ - 1, "%.*f%s%.*f%s%.*f%s%.*f", -sig_dgts, v[0], sep,
             -sig_dgts, v[1], sep, -sig_dgts, v[2], sep, -sig_dgts, v[3]);
  return buf;
}

// for alternate name look up
std::string find_alt_name(FILE *, const char *);
std::string find_alt_name(const char *fname, const char *subdir);

} // namespace anti

#endif // UTILS_H
