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

/*!\file math_utils.h
   \brief utility routines for text operations, I/O conversions, etc
*/


#ifndef UTILS_H
#define UTILS_H


#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include "vec3d.h"
#include "vec4d.h"
#include "mat3d.h"
#include "geom.h"
#include "utils_ultragetopt.h"

using std::string;
using std::vector;

///Command line processing
class prog_opts: public ultra_getopt
{
   private:
      string program_name;

   public:
      enum {   argmatch_default=0,
               argmatch_case_sensitive=1,
               argmatch_no_partial=2,
               argmatch_add_id_maps=4
      };

      static const char *help_ver_text;

      ///Constructor
      /**\param prog_name the name of the program. */
      prog_opts(string prog_name): program_name(prog_name) {}

      ///Destructor
      virtual ~prog_opts() {}

      ///Process the command line
      /**In the derived class this will process the program options
       * and arguments, probably using \c getopt. */
      virtual void process_command_line(int /*argc*/, char ** /*argv*/) {};

      ///Usage message
      /**In the derived class this will print a program usage help message*/
      virtual void usage() {};

      ///Usage message
      /**print a version message*/
      void version();

      ///Get the program name
      /**\return a pointer to the program name. */
      const char *prog_name() const { return program_name.c_str(); }
      
      ///Print a message (to standard error).
      /**The message will be preceded by the program name, the
       * message type (if given), and the option letter or
       * argument name (if given).
       * \param msg the message to print.
       * \param msg_type the message type (e.g. 'warning').
       * \param opt the option letter or argument name. */
      void message(string msg, const char *msg_type=0, string opt="") const;


      ///Print an error message (to standard error) and exit.
      /**The message will be preceded by the program name, and the
       * option letter or argument name (if given).
       * \param msg the message to print.
       * \param opt the option letter or argument name.
       * \param exit_num The value to return when the program exits. */
      void error(string msg, string opt="", int exit_num=1) const
         { message(msg, "error", opt); exit(exit_num); }
      
      ///Print an error message (to standard error) and exit.
      /**The message will be preceded by the program name, and the
       * option letter (if given).
       * \param msg the message to print.
       * \param opt the option letter.
       * \param exit_num The value to return when the program exits. */
      void error(string msg, char opt, int exit_num=1) const
         { message(msg, "error", string()+opt); exit(exit_num); }

      ///Print a warning message (to standard error).
      /**The message will be preceded by the program name, and the
       * option letter or argument name (if given).
       * \param msg the message to print.
       * \param opt the option letter or argument name. */
      void warning(string msg, string opt="") const
         { message(msg, "warning", opt); }

      ///Print a warning message (to standard error).
      /**The message will be preceded by the program name, and the
       * option letter (if given).
       * \param msg the message to print.
       * \param opt the option letter. */
      void warning(string msg, char opt) const
         { message(msg, "warning", string()+opt); }

      ///Process long options
      /**\param argc the number of arguments.
       * \param argv pointers to the argument strings. */
      void handle_long_opts(int argc, char *argv[]);

      ///Process common options
      /**\param c the character returned by getopt.
       * \param opt the option character being considered by getopt.
       * \return whether the option was handled. */
      bool common_opts(char c, char opt);

      ///Map option arguments to identifiers using matching
      /**\param arg the option argument
       * \param maps a set of maps from argument strings to identifiers
       * separated by '|', e.g. 'string1=id1|sting2=id2|string3=id3'
       * \param match_flags the default is a icase insesetive match of
       * \a arg to a string or failing that to the start of exactly one string.
       * \c argmatch_case_sensitive distinguishes case, \c argmatch_no_partial
       * disallows partial matches, \c argmatch_add_id_maps add extra maps
       * so each identifiers maps to itself.
       * \param errmsg error message
       * \return The identifier corresponding to the matched string, or
       * if the argument is not matched "" is returned and the error
       * message is copied to \a errmsg.*/
      string get_arg_id(const char *arg, const char *maps,
            unsigned int match_flags=argmatch_default, char *errmsg=0);

};

///Read a geometry from a name passed as a program argument
/**Read geometry from a name, print any messages, and error out if necessary
 * \param geom to hold the model geometry read
 * \param name file name or resource name of the model
 * \param opts program options */
void geom_read_or_error(geom_if &geom, const string &name,
      const prog_opts &opts);

///Write a geometry to a file name passed as a program argument
/**Write geometry to a file name, print any messages, and error out if necessary
 * \param geom the model geometry
 * \param name file name or resource name of the model
 * \param opts program options
 * \param sig_dgts the number of significant digits to write,
 * or if negative then the number of digits after the decimal point. */
void geom_write_or_error(const geom_if &geom, const string &name,
      const prog_opts &opts, int sig_dgts=DEF_SIG_DGTS);


//Text utilities

///Get the Basename
/**Get the name of the file from a path. This is taken as the last part
 * of the path after any '%\' or '/', otherwise the path is the file name.
 * \param path the path to the file.
 * \return A pointer to the start of the file name. */
const char *basename2(const char *path);

///Convert dots to underscores.
/**\param str the string to convert.
 * \return The converted string. */
string dots2underscores(string str);

///Read a floating point number, which may be a mathematical expression, from a string.
/**The string should only hold the floating point number or expression, but may
 * have leading and trailing whitespace.
 * \param str the string holding the floating point number or expressions.
 * \param f used to return the floating point number.
 * \param errmsg an array at least \c MSG_SZ chars long to
 * return any error message.
 * \return true if a valid floating point number was read, otherwise false
 * and the error is detailed in \a errmsg. */
bool read_double(const char *str, double *f, char *errmsg=0);

///Read a floating point number from a string.
/**The string should only hold the floating point number, but may
 * have leading and trailing whitespace.
 * \param str the string holding the floating point number.
 * \param f used to return the floating point number.
 * \param errmsg an array at least \c MSG_SZ chars long to
 * return any error message.
 * \return true if a valid floating point number was read, otherwise false
 * and the error is detailed in \a errmsg. */
bool read_double_noparse(const char *str, double *f, char *errmsg=0);

///Read an integer from a string.
/**The string should only hold the integer, but may
 * have leading and trailing whitespace.
 * \param str the string holding the integer.
 * \param i used to return the integer.
 * \param errmsg an array at least \c MSG_SZ chars long to
 * return any error message.
 * \return true if a valid integer was read, otherwise false
 * and the error is detailed in \a errmsg. */
bool read_int(const char *str, int *i, char *errmsg=0);

///Read a fraction from a string.
/**The string should only hold an integer, or two integers seperated by '/',
 * but may have leading and trailing whitespace.
 * \param str the string holding the fraction.
 * \param num used to return the numerator.
 * \param denom used to return the denominator.
 * \param errmsg an array at least \c MSG_SZ chars long to
 * return any error message.
 * \return true if a valid fraction was read, otherwise false
 * and the error is detailed in \a errmsg. */
bool read_fraction(const char *str, int *num, int *denom, char *errmsg=0);

///Read floating point numbers, which may be maththematical expressions, from a list of strings.
/**The strings should only hold the floating point number, but may
 * have leading and trailing whitespace.
 * \param vals the strings holding the floating point numbers.
 * \param nums used to return the floating point numbers.
 * \param errmsg an array at least \c MSG_SZ chars long to
 * return any error message.
 * \return true if only valid floating point numbers were read, otherwise false
 * and the error is detailed in \a errmsg. */
bool read_double_list(vector<char *> &vals, vector<double> &nums,
      char *errmsg=0);

///Read floating point numbers from a list of strings.
/**The strings should only hold the floating point number, but may
 * have leading and trailing whitespace.
 * \param vals the strings holding the floating point numbers.
 * \param nums used to return the floating point numbers.
 * \param errmsg an array at least \c MSG_SZ chars long to
 * return any error message.
 * \return true if only valid floating point numbers were read, otherwise false
 * and the error is detailed in \a errmsg. */
bool read_double_list_noparse(vector<char *> &vals, vector<double> &nums,
      char *errmsg=0);

///Read floating point numbers, which may be mathematical expressions, listed in a single string.
/**The numbers in the string should be comma separated, and may
 * have leading and trailing whitespace.
 * \param str the string holding the comma-separated floating point numbers.
 * \param nums used to return the floating point numbers.
 * \param errmsg an array at least \c MSG_SZ chars long to
 * return any error message.
 * \param len the maximum total of numbers that should be in \a str, or if
 * it is \c 0 then there is no maximum.
 * \param sep the characters that can separate the numbers.
 * \return true if only valid floating point numbers were read, and no more
 * than \a len (if \c len>0), otherwise false.
 * and the error is detailed in \a errmsg. */
bool read_double_list(char *str, vector<double> &nums, char *errmsg=0,
      int len=0, const char *sep=",");

///Read floating point numbers listed in a single string.
/**The numbers in the string should be comma separated, and may
 * have leading and trailing whitespace.
 * \param str the string holding the comma-separated floating point numbers.
 * \param nums used to return the floating point numbers.
 * \param errmsg an array at least \c MSG_SZ chars long to
 * return any error message.
 * \param len the maximum total of numbers that should be in \a str, or if
 * it is \c 0 then there is no maximum.
 * \param sep the characters that can separate the numbers.
 * \return true if only valid floating point numbers were read, and no more
 * than \a len (if \c len>0), otherwise false.
 * and the error is detailed in \a errmsg. */
bool read_double_list_noparse(char *str, vector<double> &nums, char *errmsg=0,
      int len=0, const char *sep=",");

///Read integers from a list of strings.
/**The strings should only hold the integers, but may
 * have leading and trailing whitespace.
 * \param vals the strings holding the integers.
 * \param nums used to return the integers.
 * \param errmsg an array at least \c MSG_SZ chars long to
 * return any error message.
 * \param is_index if true then the integers cannot be negative.
 * \return true if only valid integers were read, otherwise false
 * and the error is detailed in \a errmsg. */
bool read_int_list(vector<char *> &vals, vector<int> &nums, char *errmsg=0,
      bool is_index=false);

///Read integers listed in a single string.
/**The numbers in the string should be comma separated, and may
 * have leading and trailing whitespace.
 * \param str the string holding the comma-separated integers.
 * \param nums used to return the integers.
 * \param errmsg an array at least \c MSG_SZ chars long to
 * return any error message.
 * \param is_index if true then the integers cannot be negative.
 * \param len the maximum total of numbers that should be in \a str, or if
 * it is \c 0 then there is no maximum.
 * \param sep the characters that can separate the numbers.
 * \return true if only valid integers were read, and no more
 * than \a len (if \c len>0), otherwise false.
 * and the error is detailed in \a errmsg. */
bool read_int_list(char *str, vector<int> &nums, char *errmsg=0,
      bool is_index=false, int len=0, const char *sep=",");

///Read index numbers listed in a single string.
/**The string consists of comma separated index number ranges, and may
 * have leading and trailing whitespace. A number range is either a single
 * number, or a sequential list indicated by two numbers seperated by '-'.
 * If the numbers are not given they default to the first and last index
 * number respectively.
 * \param str the string holding the comma-separated integers.
 * \param nums used to return the integers.
 * \param num_idxs index numbers with this value or higher are out of range
 * \param allow_extra allow out-of-range index numbers, prefixed by an x or X,
 * which are indexed relative to num_idxs, i.e. X0 = num_idxs.
 * \param errmsg an array at least \c MSG_SZ chars long to return any
 * error message.
 * \return true if only valid index numbers were read, otherwise false.
 * and the error is detailed in \a errmsg. */
bool read_idx_list(char *str, vector<int> &nums, int num_idxs,
      bool allow_extra=false, char *errmsg=0);

/// Read a line of arbitrary length
/**The caller is responsible for freeing the memory allocated to line
 * after each read.
 * \param file the file stream to read from.
 * \param line where the line is returned
 * \return <ul>
 *    <li>\c 0 if the line was read correctly.
 *    <li>\c -1 if memory for \a line could not be allocated.
 *    <li>\c 1 if an final unterminated empty line was read.
 * </ul> */
int read_line(FILE *file, char **line);

///Split a line into delimited parts
/**\param line the line to split (this will be modified).
 * \param parts the parts of the split line.
 * \param delims the characters to use as delimiters, if \c 0 then use
 * whitespace characters.
 * \param strict if true then treat every delimiter as a separator, returning
 * null strings between adjacent delimiters, always returning at least one part.
 * \return The number of parts. */
int split_line(char *line, vector<char *> &parts,const char *delims=0,
      bool strict=false);


///Remove leading and trailing space, convert any whitespace to a single space
/**\param str the string to convert.
 * \return A pointer to the string. */
char *clear_extra_whitespace(char *str);

///Convert to a normalised resource name
/**Remove leading and trailing space, convert any whitespace to a single space, make lowercase
 * \param to the string to convert (up to \c MSG_SZ-1 characters used.)
 * \param from the string to convert (length \c MSG_SZ.)
 * \return A pointer to the \a to string. */
char *to_resource_name(char *to, const char *from);


///Open a support file
/**Tries to open a file by its name, then tries to open it in
 * \c $ANTIPRISM_DATA/sub_dir, finally tries to open it in
 * \c sub_dir in the installation data directory.
 * \param fname the name of the file to open.
 * \param subdir the data directory subdirectory to search in.
 * \param alt_name a name that is found in an alt_names.txt file before
 * a file with the name is found.
 * \param where used to return the place that the file was found <ul>
 * <li>\c 0 - locally
 * <li>\c 1 - in the ANTIPRISM_DATA directory
 * <li>\c 2 - in the installation data directory
 * </ul>
 * \param fpath used to return the full path to the file that was found.
 * \return A pointer to the opened file stream. */
FILE *open_sup_file(const char *fname, const char *subdir,
      string *alt_name=0, int *where=0, string *fpath=0);

//Convert a C formated message string to a C++ string
/** Converts the first MSG_SZ-1 characters of the C format string
 * \param fmt the formatted string
 * \param ... the values for the format
 * \return The converted string. */
string msg_str(const char *fmt, ...);

//Copy a C string
/** The copy is dynamically allocated and must be freed with free()
 * \param str the string to copy
 * \return A pointer to the newly allocated string, or 0*/
char *copy_str(const char *str);

///Convert an integer to a string
/**\param i the integer.
 * \return A pointer to \a buf, which holds the string. */
string itostr(int i);

///Convert an integer to a string
/**\param buf a buffer to return the string.
 * \param i the integer.
 * \return The string. */
char *itostr(char *buf, int i);

///Convert a floating point number to a string
/**\param f the floating point number.
 * \param sig_dgts the number of significant digits in the conversion,
 * or if negative then the number of digits after the decimal point.
 * \return The string. */
string dtostr(double f, int sig_dgts=17);

///Convert a floating point number to a string
/**\param buf a buffer to return the string.
 * \param f the floating point number.
 * \param sig_dgts the number of significant digits in the conversion,
 * or if negative then the number of digits after the decimal point.
 * \return A pointer to \a buf, which holds the string. */
char *dtostr(char *buf, double f, int sig_dgts=17);


///Convert a coordinate vector to a string
/**\param v the vector.
 * \param sep the separator between the numbers.
 * \param sig_dgts the number of significant digits in the conversion,
 * or if negative then the number of digits after the decimal point.
 * \return The string. */
string vtostr(vec3d v, const char *sep=", ", int sig_dgts=17);

///Convert a coordinate vector to a string
/**\param buf a buffer to return the string.
 *\param v the vector.
 * \param sep the separator between the numbers.
 * \param sig_dgts the number of significant digits in the conversion,
 * or if negative then the number of digits after the decimal point.
 * \return A pointer to \a buf, which holds the string. */
char* vtostr(char *buf, vec3d v, const char *sep=", ", int sig_dgts=17);


///Convert a coordinate vector to a string
/**\param buf a buffer to return the string.
 *\param v the vector.
 * \param sep the separator between the numbers.
 * \param sig_dgts the number of significant digits in the conversion,
 * or if negative then the number of digits after the decimal point.
 * \return A pointer to \a buf, which holds the string. */
inline char* vtostr(char *buf, vec4d v, const char *sep=", ", int sig_dgts=17);



// inline function definitions

inline string itostr(int i)
{ 
   char buf[MSG_SZ];
   return string(itostr(buf, i));
}

inline char *itostr(char *buf, int i)
{ 
   buf[MSG_SZ-1] = 0;
   sprintf(buf, "%d", i);
   return buf;
}


inline string dtostr(double f, int sig_dgts)
{ 
   char buf[MSG_SZ];
   return string(dtostr(buf, f, sig_dgts));
}

inline char *dtostr(char *buf, double f, int sig_dgts)
{ 
   buf[MSG_SZ-1] = 0;
   if(sig_dgts>0)
       snprintf(buf, MSG_SZ-1, "%.*g", sig_dgts, f);
   else 
       snprintf(buf, MSG_SZ-1, "%.*f", -sig_dgts, f);
   return buf;
}

inline string vtostr(vec3d v, const char *sep, int sig_dgts)
{ 
   char buf[MSG_SZ];
   return string(vtostr(buf, v, sep, sig_dgts));
}

inline char* vtostr(char *buf, vec3d v, const char *sep, int sig_dgts)
{ 
   buf[MSG_SZ-1] = 0;
   if(sig_dgts>0)
      snprintf(buf, MSG_SZ-1, "%.*g%s%.*g%s%.*g",
            sig_dgts, v[0], sep, sig_dgts, v[1], sep, sig_dgts, v[2]);
   else
      snprintf(buf, MSG_SZ-1, "%.*f%s%.*f%s%.*f",
            -sig_dgts, v[0], sep, -sig_dgts, v[1], sep, -sig_dgts, v[2]);
   return buf;
}

inline char* vtostr(char *buf, vec4d v, const char *sep, int sig_dgts)
{ 
   buf[MSG_SZ-1] = 0;
   if(sig_dgts>0)
      snprintf(buf, MSG_SZ-1, "%.*g%s%.*g%s%.*g%s%.*g",
            sig_dgts, v[0], sep, sig_dgts, v[1], sep, sig_dgts, v[2],
            sep, sig_dgts, v[3]);
   else
      snprintf(buf, MSG_SZ-1, "%.*f%s%.*f%s%.*f%s%.*f",
            -sig_dgts, v[0], sep, -sig_dgts, v[1], sep, -sig_dgts, v[2],
            sep, -sig_dgts, v[3]);
   return buf;
}

// for alternate name look up
string find_alt_name(FILE *, const char *);
string find_alt_name(const char *fname, const char *subdir);

#endif // UTILS_H

