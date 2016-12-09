/*
   Copyright (c) 2012-2016, Adrian Rossiter

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

/* \file utils_parser.cc
   \brief parse mathematical expressions (using muParser)
*/

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "mathutils.h"
#include "muparser/muParser.h"
#include "utils.h"

// The muParser header includes <cmath> which conflicts with <math.h>
// on, at least, Cygwin, regarding isnan() and isinf(). Include all
// of namespace std to make sure functions are found.
using namespace std;

using namespace anti;
using namespace mu;

const int num_vars = 10; // number of variables of type var0, var1, ...

class ExpParser : public Parser {
private:
  double vars[num_vars];

  static double deg_sin(double a) { return sin(deg2rad(a)); }
  static double deg_cos(double a) { return cos(deg2rad(a)); }
  static double deg_tan(double a) { return tan(deg2rad(a)); }
  static double deg_asin(double x) { return rad2deg(asin(x)); }
  static double deg_acos(double x) { return rad2deg(acos(x)); }
  static double deg_atan(double x) { return rad2deg(atan(x)); }
  static double deg_atan2(double x, double y) { return rad2deg(atan2(x, y)); }
  static double deg(double a) { return rad2deg(a); }
  static double rad(double a) { return deg2rad(a); }

  static int rt_tok(const char *tok, int *pos, double *val);

public:
  ExpParser();
};

// Determine if a token is of form rt2.2, and set value to sqrt(2.2)
int ExpParser::rt_tok(const char *tok, int *pos, double *val)
{
  int len = strlen(tok);
  int prefix_len = 2; // "rt"

  if (len <= prefix_len)
    return false;
  if (strncmp(tok, "rt", prefix_len))
    return false;

  // Read a decimal number string made of digits and up to one decimal point
  string num_str;
  int dec_point_cnt = 0;
  int p;
  for (p = prefix_len; p < len; p++) {
    if (isdigit(tok[p]) || (tok[p] == '.' && dec_point_cnt < 1)) {
      num_str += tok[p];
      if (tok[p] == '.')
        dec_point_cnt++;
    }
    else
      break;
  }

  double num;
  Status stat = read_double_noparse(num_str.c_str(), &num);
  if (stat.is_error())
    throw exception_type(string("rt: invalid number: ") + stat.msg());

  *val = sqrt(num);
  *pos += p;
  return true;
}

ExpParser::ExpParser() : Parser()
{
  // Seperator for expressions and function arguments
  SetArgSep(';');

  for (int i = 0; i < num_vars; i++) {
    vars[i] = 0.0;
    DefineVar(msg_str("var%d", i), &vars[i]);
  }

  // Replacement functions
  DefineFun(_T("sin"), deg_sin);
  DefineFun(_T("cos"), deg_cos);
  DefineFun(_T("tan"), deg_tan);
  DefineFun(_T("asin"), deg_asin);
  DefineFun(_T("acos"), deg_acos);
  DefineFun(_T("atan"), deg_atan);
  DefineFun(_T("atan2"), deg_atan2);

  // New functions
  DefineFun(_T("deg"), deg2rad);
  DefineFun(_T("rad"), rad2deg);

  // Clear existing constants
  ClearConst();

  // New constants
  DefineConst("phi", phi);
  DefineConst("pi", M_PI);

  // Tokens of form rt2.2 will return sqrt(2.2)
  AddValIdent(rt_tok);
}

namespace anti {

Status read_double(const char *str, double *f)
{
  char msg_type[] = "maths expression";
  Status stat;
  try {
    ExpParser p;

    p.SetExpr(str);

    *f = p.Eval();

    if (isnan(*f))
      stat.set_error(
          msg_str("%s: result is not a number (domain error, etc)", msg_type));
    else if (isinf(*f))
      stat.set_error(
          msg_str("%s: result is not a finite number (division by zero, etc)",
                  msg_type));
  }
  catch (Parser::exception_type &e) {
    stat.set_error(msg_str("%s: %s", msg_type, e.GetMsg().c_str()));
  }

  return stat;
}

} // namespace anti
