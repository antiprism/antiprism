/*
   Copyright (c) 2003-2012, Adrian Rossiter

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

/* \file antiparser.cc
   \brief parse mathematical expressions (using muParser)
*/

#ifdef HAVE_CONFIG_H
   #include "../config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <cmath>
#include "muparser/muParser.h"

const int MSG_SZ = 256;

using namespace mu;
using namespace std;

class AntiParser : public Parser
{
   private:
      double vars[10];
      static double d2r(double ang) { return ang * M_PI/180; }
      static double r2d(double ang) { return ang * 180/M_PI; }
      static double deg_sin(double a)  { return sin(d2r(a)); }
      static double deg_cos(double a)  { return cos(d2r(a)); }
      static double deg_tan(double a)  { return tan(d2r(a)); }
      static double deg_asin(double x) { return r2d(asin(x)); }
      static double deg_acos(double x) { return r2d(acos(x)); }
      static double deg_atan(double x) { return r2d(atan(x)); }
      static double deg_atan2(double x, double y) { return r2d(atan2(x, y)); }
      static double deg(double a)      { return r2d(a); }
      static double rad(double a)      { return d2r(a); }
   public:
      AntiParser();
};

AntiParser::AntiParser() : Parser()
{
   char var_name[16];
   for(int i=0; i<10; i++) {
      sprintf(var_name, "var%d", i);
      DefineVar(var_name, &vars[i]);
   }
   SetArgSep(';');

   // Replacement functions
   DefineFun(_T("sin"), deg_sin);
   DefineFun(_T("cos"), deg_cos);
   DefineFun(_T("tan"), deg_tan);
   DefineFun(_T("asin"), deg_asin);
   DefineFun(_T("acos"), deg_acos);
   DefineFun(_T("atan"), deg_atan);
   DefineFun(_T("atan2"), deg_atan2);
   //New functions
   DefineFun(_T("deg"), deg);
   DefineFun(_T("rad"), rad);

   // Clear existing constants
   ClearConst();

   //New constants
   DefineConst("phi", (sqrt(5)+1)/2);
   DefineConst("pi", M_PI);
}




bool read_double(const char *str, double *f, char *errmsg)
{
   char msg_type[] = "maths expression";
   bool exp_good = true;
   try {
      AntiParser p;

      p.SetExpr(str);

      *f = p.Eval();

      if(isnan(*f)) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ,
                  "%s: result is not a number (domain error, etc)", msg_type);
         exp_good = false;
      }
      else if(isinf(*f)) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ,
                  "%s: result is not a finite number (division by zero, etc)",
                  msg_type);
         exp_good = false;
      }
   }
   catch (Parser::exception_type &e) {
      if(errmsg)
         snprintf(errmsg, MSG_SZ, "%s: %s", msg_type, e.GetMsg().c_str());
      exp_good = false;
   }

   return exp_good;
}


