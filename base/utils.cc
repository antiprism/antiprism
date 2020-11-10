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

/* \file utils.cc
   \brief utility routines for maths operations, text operations,
   I/O conversions, etc
*/

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <cctype>
#include <climits>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

#include "utils.h"

using std::map;
using std::string;
using std::vector;

namespace anti {

/// Whitespace characters
const char WHITESPACE[] = " \t\r\n\f\v";

const char *basename2(const char *path) // basename - forward and back slashes
{
  const char *fpart = path;
  while (*path) {
    if (*path == '/' || *path == '\\')
      fpart = path + 1;
    path++;
  }
  return fpart;
}

string dots2underscores(string str)
{
  for (char &i : str)
    if (i == '.')
      i = '_';
  return str;
}

Status read_double_noparse(const char *str, double *f)
{
  bool to_sqrt;
  char buff;
  if (sscanf(str, " sqrt%lf %c", f, &buff) == 1)
    to_sqrt = true;
  else if (sscanf(str, " %lf %c", f, &buff) == 1)
    to_sqrt = false;
  else
    return Status::error("not a number");

  if (std::isinf(*f))
    return Status::error("number too large\n");

  if (std::isnan(*f))
    return Status::error("not a number\n");

  if (to_sqrt)
    *f = sqrt(*f);

  return Status::ok();
}

Status read_int(const char *str, int *i)
{
  char buff;
  if (sscanf(str, " %d %c", i, &buff) != 1)
    return Status::error("not an integer");

  if (*i == INT_MAX)
    return Status::error("integer too large\n");

  return Status::ok();
}

Status read_int_list(const vector<char *> &vals, vector<int> &nums,
                     bool is_index)
{
  nums.clear();
  int num;
  for (auto v_str : vals) {
    if (!read_int(v_str, &num))
      return Status::error(msg_str("'%s' is not an integer", v_str));

    if (is_index && num < 0)
      return Status::error(msg_str("'%s' is not a positive integer", v_str));

    nums.push_back(num);
  }

  return Status::ok();
}

Status read_int_list(const char *str, vector<int> &nums, bool is_index, int len,
                     const char *sep)
{
  nums.clear();
  int vec_idx;

  string str_cpy(str);         // copy, do not access as C++ string
  char *str_ptr = &str_cpy[0]; // may be used to modify characters

  char *v_str = strtok(str_ptr, sep);
  int i = 0;
  while (v_str) {
    i++;
    if (!read_int(v_str, &vec_idx))
      return Status::error(msg_str("'%s' is not an integer", v_str));

    if (is_index && vec_idx < 0)
      return Status::error(msg_str("'%s' is not a positive integer", v_str));

    if (len && i > len)
      return Status::error(msg_str("more than %d integers given", len));

    nums.push_back(vec_idx);
    v_str = strtok(nullptr, sep);
  }

  return Status::ok();
}

static Status read_idx(const char *str, int *idx, int num_idxs)
{
  if (!read_int(str, idx))
    return Status::error(msg_str("'%s' is not an integer", str));

  if (*idx < 0)
    return Status::error(msg_str("'%s' is not a positive integer", str));

  if (*idx >= num_idxs) {
    if (num_idxs == 0)
      return Status::error(msg_str("'%s' is out of range (no elements of "
                                   "this type)",
                                   str));
    else
      return Status::error(msg_str("'%s' is out of range (last index is %d)",
                                   str, num_idxs - 1));
  }

  return Status::ok();
}

Status read_idx_list(const char *str, vector<int> &nums, int num_idxs,
                     bool allow_extra)
{
  Status stat;
  nums.clear();
  int idx, idx2;
  char *p;

  string str_cpy(str);         // copy, do not access as C++ string
  char *str_ptr = &str_cpy[0]; // may be used to modify characters

  char *v_str = strtok(str_ptr, ",");
  while (v_str) {
    if ((p = strchr(v_str, '-'))) { // process a range
      *p = '\0';                    // terminate first index
      if (*v_str) {
        if (!(stat = read_idx(v_str, &idx, num_idxs)))
          return stat;
      }
      else
        idx = 0;
      if (*(p + 1)) {
        if (!(stat = read_idx(p + 1, &idx2, num_idxs)))
          return stat;
      }
      else
        idx2 = num_idxs - 1;
      if (*v_str && *(p + 1) && idx > idx2)
        return Status::error(
            msg_str("index range, %s is greater than %s", v_str, p + 1));

      if ((*v_str || *(p + 1)) && !num_idxs)
        return Status::error(
            "invalid range, there are no elements of the query type");

      for (int i = idx; i <= idx2; i++)
        nums.push_back(i);
    }
    else {
      int extra = false;
      if (allow_extra && (*v_str == 'x' || *v_str == 'X')) {
        extra = true;
        v_str++;
      }
      if (v_str) {
        if (!(stat = read_idx(v_str, &idx, num_idxs)))
          return stat;
      }
      nums.push_back(idx + extra * num_idxs);
    }
    v_str = strtok(nullptr, ",");
  }

  return Status::ok();
}

static Status read_double_list(const vector<char *> &vals, vector<double> &nums,
                               bool parse)
{
  nums.clear();
  double num;
  for (auto &val : vals) {
    Status stat =
        (parse) ? read_double(val, &num) : read_double_noparse(val, &num);
    if (stat.is_error())
      return Status::error(msg_str("%s: '%s'", stat.c_msg(), val));

    nums.push_back(num);
  }
  return Status::ok();
}

Status read_double_list(const vector<char *> &vals, vector<double> &nums)
{
  return read_double_list(vals, nums, true);
}

Status read_double_list_noparse(const vector<char *> &vals,
                                vector<double> &nums)
{
  return read_double_list(vals, nums, false);
}

static Status read_double_list(const char *str, vector<double> &nums, int len,
                               const char *sep, bool parse)
{
  nums.clear();

  string str_cpy(str);         // copy, do not access as C++ string
  char *str_ptr = &str_cpy[0]; // may be used to modify characters

  double num;
  char *num_str = strtok(str_ptr, sep);
  int i = 0;
  while (num_str) {
    i++;
    Status stat = (parse) ? read_double(num_str, &num)
                          : read_double_noparse(num_str, &num);
    if (stat.is_error())
      return Status::error(msg_str("%s: '%s'", stat.c_msg(), num_str));

    if (len && i > len)
      return Status::error(msg_str("more than %d numbers given", len));

    nums.push_back(num);
    num_str = strtok(nullptr, sep);
  }

  return Status::ok();
}

Status read_double_list(const char *str, vector<double> &nums, int len,
                        const char *sep)
{
  return read_double_list(str, nums, len, sep, true);
}

Status read_double_list_noparse(const char *str, vector<double> &nums, int len,
                                const char *sep)
{
  return read_double_list(str, nums, len, sep, false);
}

Status read_fraction(const char *str, int *num, int *denom)
{
  Status stat;
  Status stat2;

  string str_cpy(str);         // copy, do not access as C++ string
  char *str_ptr = &str_cpy[0]; // may be used to modify characters

  *denom = 1;
  char *p = strchr(str_ptr, '/');
  if (p != nullptr) {
    *p++ = '\0';
    if (!(stat2 = read_int(p, denom)))
      return stat.set_error("denominator, " + stat2.msg());
  }
  if (!(stat2 = read_int(str_ptr, num)))
    return stat.set_error("numerator, " + stat2.msg());

  return stat;
}

int read_line(FILE *file, char **line)
{

  int linesize = 128;
  *line = (char *)malloc(linesize);
  if (!*line)
    return -1;

  int offset = 0;
  while (true) {
    if (!fgets(*line + offset, linesize - offset, file)) {
      if (offset != 0)
        return 0;
      else {
        *(*line + offset) = '\0'; // terminate the line
        return (ferror(file)) ? -1 : 1;
      }
    }
    int len = offset + strlen(*line + offset);
    if ((*line)[len - 1] == '\n') {
      (*line)[len - 1] = 0;
      return 0;
    }
    offset = len;

    char *newline = (char *)realloc(*line, linesize * 2);
    if (!newline)
      return -1;
    *line = newline;
    linesize *= 2;
  }
}

// remove leading and trailing space, covert whitepspace to single space
char *clear_extra_whitespace(char *str)
{
  char *p = str;
  int cnt = 0;
  bool prev_is_space = true;
  while (*p) {
    if (strchr(WHITESPACE, *p)) {
      if (prev_is_space) {
        p++;
        continue;
      }
      prev_is_space = true;
      *p = ' ';
    }
    else
      prev_is_space = false;
    str[cnt++] = *p++;
  }
  str[cnt] = '\0';
  if (cnt && prev_is_space)
    str[cnt - 1] = '\0';

  return str;
}

void clear_extra_whitespace(string &str)
{
  clear_extra_whitespace(&str[0]); // clear whitespace inplace
  str.resize(strlen(&str[0]));     // resize to the length of the C string
}

char *to_resource_name(char *to, const char *from)
{
  strcpy_msg(to, from);
  clear_extra_whitespace(to);
  for (char *p = to; *p; p++)
    *p = tolower(*p);
  return to;
}

void backslash_to_forward(string &path)
{
  for (char &si : path)
    if (si == '\\')
      si = '/';
}

FILE *fopen_file(string &fpath)
{
  backslash_to_forward(fpath);
  FILE *file = fopen(fpath.c_str(), "r");
  if (file) {
    struct stat st;
    fstat(fileno(file), &st);
    if (S_ISDIR(st.st_mode)) {
      fclose(file);
      file = nullptr;
    }
  }
  return file;
}

string find_alt_name(FILE *afile, const char *a_name)
{
  const int line_size = 1024;
  char line[line_size];
  char aname[line_size];
  strncpy(aname, a_name, line_size - 1);
  aname[line_size - 1] = '\0';
  clear_extra_whitespace(aname);

  int line_no = 0;
  while (fgets(line, line_size, afile)) {
    line_no++;

    // ignore comments
    char *first_hash = strchr(line, '#');
    if (first_hash)
      *first_hash = '\0';

    char *altname, *name;
    // skip lines without =
    if (!(altname = strtok(line, "=")))
      continue;

    if ((name = strtok(nullptr, "\n"))) {
      if (strcasecmp(clear_extra_whitespace(altname), aname) == 0) {
        clear_extra_whitespace(name);
        for (char *p = name; *p; p++)
          *p = tolower(*p);
        return string((name));
      }
    }
  }

  return string("");
}

string find_alt_name(const char *fname, const char *subdir)
{
  string aname;
  FILE *alt = open_sup_file("alt_names.txt", subdir);
  if (alt) {
    char f_name[MSG_SZ];
    strcpy_msg(f_name, fname);
    aname = find_alt_name(alt, fname);
    fclose(alt);
  }
  return aname;
}

FILE *open_file_data(const string &dir, const string &fname,
                     string *aname = nullptr)
{
  if (aname)
    *aname = "";

  FILE *file = nullptr;
  // convert fname to lowercase
  char f_name[MSG_SZ];
  unsigned int i;
  for (i = 0; i < fname.size() && i < MSG_SZ - 1; i++)
    f_name[i] = tolower(fname[i]);
  f_name[i] = '\0';

  string fpath;
  // don't allow escape from data directory
  if (!strchr(f_name, '\\') && !strchr(f_name, '/')) {
    fpath = dir + f_name;
    if ((file = fopen_file(fpath)))
      return file;
  }
  string alt_names_file = "alt_names.txt";
  fpath = dir + alt_names_file;
  FILE *alt = fopen_file(fpath);
  if (alt) {
    string name = find_alt_name(alt, f_name);
    fclose(alt);
    if (name != "" && !strchr(name.c_str(), '\\') &&
        !strchr(name.c_str(), '/')) {
      fpath = dir + name;
      if ((file = fopen_file(fpath)))
        return file;
    }
    *aname = name; // may be used to open an internal resource
  }

  return nullptr;
}

FILE *open_sup_file(const char *fname, const char *subdir, string *alt_name,
                    int *where, string *fpath)
{
  string alt_nam;
  if (!alt_name)
    alt_name = &alt_nam;

  int whr;
  if (!where)
    where = &whr;

  string fpth;
  if (!fpath)
    fpath = &fpth;

  string aname;
  FILE *file;
  // try to open fname alone
  *fpath = fname;
  *where = 0; // local
  if ((file = fopen_file(*fpath)))
    return file;

  // The file open will fail, but will read alt_name.txt in local directory
  if ((file = open_file_data("", *fpath, alt_name)))
    return file;
  if (*alt_name != "")
    return nullptr;

  // try environment variable for data directory
  char *sup_dir = getenv("ANTIPRISM_DATA");
  if (sup_dir) {
    *where = 1; // environment
    string fdir = string(sup_dir) + subdir;
    if ((file = open_file_data(fdir, *fpath, alt_name)))
      return file;
    if (*alt_name != "")
      return nullptr;
  }

  // try hardcoded install path for data directory
  *where = 2; // installed
  string fdir = string(SUPDIR) + subdir;
  if ((file = open_file_data(fdir, *fpath, alt_name)))
    return file;
  if (*alt_name != "")
    return nullptr;

  return nullptr;
}

char *strcpy_msg(char *dest, const char *src)
{
  strncpy(dest, src, MSG_SZ - 1);
  dest[MSG_SZ - 1] = '\0';
  return dest;
}

string msg_str(const char *fmt, ...)
{
  char message[MSG_SZ];
  va_list args;
  va_start(args, fmt);
  vsnprintf(message, MSG_SZ - 1, fmt, args);
  return message;
}

namespace {
int split_line(char *line, vector<char *> &parts, const char *delims,
               bool strict)
{
  parts.clear();
  if (!delims)
    delims = WHITESPACE;

  if (strict) {
    char *cur = line;
    parts.push_back(cur);                    // always an entry, even if null
    while (*(cur += strcspn(cur, delims))) { // quit at end of string
      *cur = '\0';                           // terminate part
      cur++;                                 // start of next part
      parts.push_back(cur);                  // add even if final null
    }
  }
  else {
    char *val;
    if (!(val = strtok(line, delims)))
      return 0;

    parts.push_back(val);
    while ((val = strtok(nullptr, delims)))
      parts.push_back(val);
  }

  return parts.size();
}
}; // namespace

int Split::init(const char *line, const char *delims, bool strict)
{
  data = line;
  return split_line(&data[0], parts, delims, strict);
}

} // namespace anti
