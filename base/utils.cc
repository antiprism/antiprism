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

/* \file utils.cc
   \brief utility routines for maths operations, text operations,
   I/O conversions, etc
*/

#ifdef HAVE_CONFIG_H
   #include "../config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>

#include "muparser/muParser.h"
#include "utils.h"

///Whitespace characters
const char WHITESPACE[] = " \t\r\n\f\v";


const char *basename2(const char *path)  // basename - forward and back slashes
{
  const char *fpart = path;
  while(*path) {
     if (*path == '/' || *path == '\\')
        fpart = path+1;
     path++;
  }
  return fpart;
}

string dots2underscores(string str)
{ 
   for(unsigned int i=0; i<str.size(); i++)
      if(str[i]=='.')
         str[i] = '_';
   return str;
}


const char *prog_opts::help_ver_text =
"  -h,--help this help message (run 'off_util -H help' for general help)\n"
"  --version version information\n";


void prog_opts::message(string msg, const char *msg_type, string opt) const
{
   fprintf(stderr, "%s: ", prog_name());
   if(msg_type)
      fprintf(stderr, "%s: ", msg_type);
   if(opt != "") {
     if(opt.size()==1 || opt[0]=='\0') 
        fprintf(stderr, "option -%s: ", opt.c_str());
     else
        fprintf(stderr, "%s: ", opt.c_str());
   }
   
   fprintf(stderr, "%s\n", msg.c_str());
}


bool prog_opts::common_opts(char c, char opt)
{
   switch(c) {
      case 'h':
         usage();
         exit(0);

      case '?':
         error("unknown option", string("-")+opt);

      case ':':
         error("missing argument", string("-")+opt);

      default:
         return false;
   }
   return true;
}

void prog_opts::handle_long_opts(int argc, char *argv[])
{
   for(int i=1; i<argc; i++) {
      if(strcmp(argv[i], "--help")==0) {
         usage();
         exit(0);
      }
      else if(strcmp(argv[i], "--version")==0) {
         version();
         exit(0);
      }
      else if(strncmp(argv[i], "--", 2)==0 && strlen(argv[i])>2)
         error("unknown option", argv[i]);
   }
}


string prog_opts::get_arg_id(const char *arg, const char *maps,
      unsigned int match_flags, char *errmsg)
{
   if(errmsg)
      *errmsg = '\0';

   char arg_lower[MSG_SZ];
   bool ignore_case = !(argmatch_case_sensitive & match_flags);
   if(ignore_case) {              // make lowercase copy of arg
      char *q=arg_lower;
      for(const char *p=arg; *p; ++p)
         *q++ = tolower(*p);
      *q = '\0';
   }
   const char *argu = (ignore_case) ? arg_lower : arg;

   // set up arg -> id map
   map<string, string> mps;
   string argument, id;
   int pos_no = 0;
   int in_arg = true;
   const char *p = maps;
   while(true) {
      if(*p=='|' || *p=='\0') {     // end of argument, id pair
         if(id == "")               // no '=' in pair so ID is position number
            id = itostr(pos_no);
         if(argument != "")         // don't stor empty pairs
            mps[argument] = id;
         if(!*p)                    // end of string, finish processing
            break;
         argument.clear();          // start of new argument, id pair
         id.clear();
         in_arg=true;
         pos_no++;
      }
      else if(*p=='=' && in_arg)    // end of arg, consider id next
         in_arg = false;
      else {
         if(in_arg) {
            // store argument as lowercase if 'no_case' compare
            char q = (argmatch_case_sensitive & match_flags) ? *p : tolower(*p);
            argument += q;
         }
         else
            id += *p;
      }

      p++;
   }

   map<string, string>::iterator mi;
   if(argmatch_add_id_maps & match_flags) {
      for(mi=mps.begin(); mi!=mps.end(); ++mi)
         if(mps.find(mi->second)==mps.end())
            mps[mi->second] = mi->second;
   }

   // look for match
   mi = mps.lower_bound(argu);

   if(mi!=mps.end() &&  mi->first==argu)   // exact match
      return mi->second;

   // look for valid partial match
   if( !(argmatch_no_partial & match_flags) &&  // partial matches allowed
         mi!=mps.end() && mi->first.substr(0, strlen(argu))==argu) {
      // is partial match, check for ambiguity
      map<string, string>::iterator mi2 = mi; 
      bool ambiguous = false;
      // run through subsequent partial matches
      while(++mi2!=mps.end() && mi2->first.substr(0, strlen(argu))==argu)
         if(mi2->second != mi->second) { // different values, ambiguous
            ambiguous = true;
            break;
         }

      if(!ambiguous)
         return mi->second;
      else {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "ambiguous argument '%s'", arg);
         return "";
      }
   }

   // no valid match found
   if(errmsg)
      snprintf(errmsg, MSG_SZ, "invalid argument '%s'", arg);
   return "";
}


void prog_opts::version()
{ 
   fprintf(stdout, "%s: Antiprism %s - http://www.antiprism.com\n",
         prog_name(), VERSION);
}

using namespace mu;

class AntiParser : public Parser
{
   private:
      double vars[10];
   protected:
      static double deg_sin(double a)  { return sin(deg2rad(a)); }
      static double deg_cos(double a)  { return cos(deg2rad(a)); }
      static double deg_tan(double a)  { return tan(deg2rad(a)); }
      static double deg_asin(double x) { return rad2deg(asin(x)); }
      static double deg_acos(double x) { return rad2deg(acos(x)); }
      static double deg_atan(double x) { return rad2deg(atan(x)); }
      static double deg_atan2(double x, double y)
         { return rad2deg(atan2(x, y)); }
      static double deg(double a)      { return rad2deg(a); }
      static double rad(double a)      { return deg2rad(a); }
   public:
      AntiParser();
};

AntiParser::AntiParser() : Parser()
{

   for(int i=0; i<10; i++)
      DefineVar(msg_str("var%d", i), &vars[i]);
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
   DefineConst("rt2", sqrt(2));
   DefineConst("rt3", sqrt(3));
   DefineConst("rt5", sqrt(5));
   DefineConst("phi", phi);
   DefineConst("pi", M_PI);
}




bool read_double(const char *str, double *f, char *errmsg)
{
   char msg_type[] = "maths expression:";
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


bool read_double_noparse(const char *str, double *f, char *errmsg)
{
   bool to_sqrt;
   char buff;
   if(sscanf(str, " sqrt%lf %c", f,  &buff) == 1)
      to_sqrt = true;
   else if(sscanf(str, " %lf %c", f,  &buff) == 1)
      to_sqrt = false;
   else {
      if(errmsg)
         strcpy(errmsg, "not a number");
      return false;
   }

   if(isinf(*f)) {
      if(errmsg)
         sprintf(errmsg, "number too large\n");
      return false;
   }

   if(isnan(*f)) {
      if(errmsg)
         sprintf(errmsg, "not a number\n");
      return false;
   }

   if(to_sqrt)
      *f = sqrt(*f);

   return true;
}


bool read_int(const char *str, int *i, char *errmsg)
{
   char buff;
   if( sscanf(str, " %d %c", i,  &buff) != 1) {
      if(errmsg)
         strcpy(errmsg, "not an integer");
      return false;
   }

   if(*i==INT_MAX) {
      if(errmsg)
         sprintf(errmsg, "integer too large\n");
      return false;
   }

   return true;
}



bool read_int_list(vector<char *> &vals, vector<int> &nums, char *errmsg, bool is_index)
{
   nums.clear();
   int num;
   for(unsigned int i=0; i<vals.size(); i++) {
      char *v_str = vals[i];
      if(!read_int(v_str, &num)) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "'%s' is not an integer", v_str);
         return false;
      }
      if(is_index && num<0) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "'%s' is not a positive integer", v_str);
         return false;
      }
      nums.push_back(num);
   }     
   return true;
}

bool read_int_list(char *str, vector<int> &nums, char *errmsg, bool is_index,
      int len, const char *sep)
{
   nums.clear();
   int vec_idx;
   char *v_str = strtok(str, sep);
   int i=0;
   while(v_str) {
      i++;
      if(!read_int(v_str, &vec_idx)) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "'%s' is not an integer", v_str);
         return false;
      }
      if(is_index && vec_idx<0) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "'%s' is not a positive integer", v_str);
         return false;
      }
      if(len && i>len) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "more than %d integers given", len);
         return false;
      }
      nums.push_back(vec_idx);
      v_str = strtok(NULL, sep);
   }
     
   return true;
}


static bool read_idx(char *str, int *idx, int num_idxs, char *errmsg)
{
   if(!read_int(str, idx, errmsg)) {
      snprintf(errmsg, MSG_SZ, "'%s' is not an integer", str);
      return false;
   }
   if(*idx<0) {
      snprintf(errmsg, MSG_SZ, "'%s' is not a positive integer", str);
      return false;
   }
   if(*idx>=num_idxs) {
      if(num_idxs==0)
         snprintf(errmsg, MSG_SZ, "'%s' is out of range (no elements of "
               "this type)", str);
      else
         snprintf(errmsg, MSG_SZ, "'%s' is out of range (too large)", str);
      return false;
   }
   return true;
}

bool read_idx_list(char *str, vector<int> &nums, int num_idxs,
      bool allow_extra, char *errmsg)
{
   nums.clear();
   int idx, idx2;
   char *p;
   char *v_str = strtok(str, ",");
   while(v_str) {
      if((p=strchr(v_str, '-'))) {  // process a range
         *p = '\0';  // terminate first index
         if(*v_str) {
            if(!read_idx(v_str, &idx, num_idxs, errmsg))
               return false;
         }
         else
            idx=0;
         if(*(p+1)) {
            if(!read_idx(p+1, &idx2, num_idxs, errmsg))
               return false;
         }
         else
            idx2=num_idxs-1;
         if(*v_str && *(p+1) && idx>idx2) {
            snprintf(errmsg, MSG_SZ, "index range, %s is greater than %s",
                  v_str, p+1);
            return false;
         }
         if((*v_str || *(p+1)) && !num_idxs) {
            snprintf(errmsg, MSG_SZ, "invalid range, there are no elements of the query type");
            return false;
         }
         for(int i=idx; i<=idx2; i++)
            nums.push_back(i);
      }
      else {
         int extra = false;
         if(allow_extra && (*v_str=='x'||*v_str=='X')) {
            extra = true;
            v_str++;
         }
         if(v_str) {
            if(!read_idx(v_str, &idx, num_idxs, errmsg))
               return false;
         }
         nums.push_back(idx + extra*num_idxs);
      }
      v_str = strtok(NULL, ",");
   }

   return true;
}


static bool read_double_list(vector<char *> &vals, vector<double> &nums,
      char *errmsg, bool parse)
{
   char errmsg2[MSG_SZ];
   nums.clear();
   double num;
   for(unsigned int i=0; i<vals.size(); i++) {
      int ret = (parse) ? read_double(vals[i], &num, errmsg2)
                        : read_double_noparse(vals[i], &num, errmsg2);
      if(!ret ) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "%s: '%s'", errmsg2, vals[i]);
         return false;
      }
      nums.push_back(num);
   }
   return true;
}

bool read_double_list(vector<char *> &vals, vector<double> &nums, char *errmsg)
{
   return read_double_list(vals, nums, errmsg, true);
}

bool read_double_list_noparse(vector<char *> &vals, vector<double> &nums,
      char *errmsg)
{
   return read_double_list(vals, nums, errmsg, false);
}


static bool read_double_list(char *str, vector<double> &nums, char *errmsg,
      int len, const char *sep, bool parse)
{
   char errmsg2[MSG_SZ];
   nums.clear();
   double num;
   char *num_str = strtok(str, sep);
   int i=0;
   while(num_str) {
      i++;
      int ret = (parse) ? read_double(num_str, &num, errmsg2)
                        : read_double_noparse(num_str, &num, errmsg2);
      if(!ret ) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "%s: '%s'", errmsg2, num_str);
         return false;
      }
      if(len && i>len) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "more than %d numbers given", len);
         return false;
      }
      nums.push_back(num);
      num_str = strtok(NULL, sep);
   }

   return true;
}

bool read_double_list(char *str, vector<double> &nums, char *errmsg,
      int len, const char *sep)
{
   return read_double_list(str, nums, errmsg, len, sep, true);
}

bool read_double_list_noparse(char *str, vector<double> &nums, char *errmsg,
      int len, const char *sep)
{
   return read_double_list(str, nums, errmsg, len, sep, false);
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
         if(offset != 0)
            return 0;
         else {
            *(*line+offset) = '\0'; // terminate the line
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

int split_line(char *line, vector<char *> &vals, const char *delims,bool strict)
{
   vals.clear();
   if(!delims)
      delims = WHITESPACE;

   if(strict) {
      char *cur = line;
      vals.push_back(cur);                     // always an entry, even if null
      while(*(cur += strcspn(cur, delims))) {  // quit at end of string
         *cur = '\0';                          // terminate part
         cur++;                                // start of next part
         vals.push_back(cur);                  // add even if final null
      }
      /*while(*cur) {                             // quit at end of string
         cur += strcspn(cur, delims);           // next delimiter
         if(*cur) {                             // part ended with delimiter
            *cur = '\0';                        // terminate part
            cur++;                              // start of next part
            vals.push_back(cur);                // add even if final null
         }
      }*/
   }
   else {
      char *val;
      if(!(val = strtok(line, delims)))
         return 0;

      vals.push_back(val);
      while((val=strtok(NULL, delims)))
         vals.push_back(val);
   }

   return vals.size();
}


void backslash_to_forward(string &path)
{
   for(string::iterator si=path.begin(); si!=path.end(); ++si)
      if(*si=='\\')
         *si = '/';
}

// remove leading and trailing space, covert whitepspace to single space
char *clear_extra_whitespace(char *str)
{
   char *p = str;
   int cnt=0;
   bool prev_is_space=true;
   while(*p) {
      if(strchr(WHITESPACE, *p)) {
         if(prev_is_space) {
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
   if(cnt && prev_is_space)
      str[cnt-1] = '\0';
   
   return str;
}

char *to_resource_name(char *to, const char *from)
{
   strncpy(to, from, MSG_SZ);
   to[MSG_SZ-1] = '\0';
   clear_extra_whitespace(to);
   for(char *p=to; *p; p++)
      *p = tolower(*p);
   return to;
}

char *copy_str(const char *str)
{
   char *p = (char *)malloc(strlen(str) + 1);
   if(p) {
      strcpy(p, str);
   }
   return p;
}



FILE *fopen_file(string &fpath)
{
   backslash_to_forward(fpath);
   FILE *file=fopen(fpath.c_str(), "r");
   if(file) {
      struct stat st;
      fstat(fileno(file), &st);
      if(S_ISDIR(st.st_mode)) {
         fclose(file);
         file = 0;
      }
   }
   return file;
}
         
string find_alt_name(FILE *afile, const char *a_name)
{
   const int line_size=1024;
   char line[line_size];
   char aname[line_size];
   strncpy(aname, a_name, line_size);
   aname[line_size-1] = '\0';
   clear_extra_whitespace(aname);
   
   int line_no = 0;
   while(fgets(line, line_size, afile)) {
      line_no++;
     
      // ignore comments
      char *first_hash = strchr(line, '#');
      if(first_hash)
         *first_hash = '\0';

      char *altname, *name;
      // skip lines without =
      if(!(altname = strtok(line, "=")))
         continue;

      if((name = strtok(NULL, "\n"))) {
         if(strcasecmp(clear_extra_whitespace(altname), aname)==0) {
            clear_extra_whitespace(name);
            for(char *p=name; *p; p++)
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
   if(alt) {
      char f_name[MSG_SZ];
      strcpy(f_name,fname);
      aname = find_alt_name(alt, fname);
      fclose(alt);
   }
   return aname;
}


FILE *open_file_data(const string &dir, const string &fname, string *aname=0)
{
   if(aname)
      *aname = "";

   FILE *file = 0;
   // convert fname to lowercase
   char f_name[MSG_SZ];
   unsigned int i;
   for(i=0; i<fname.size() && i<MSG_SZ-1; i++)
      f_name[i] = tolower(fname[i]);
   f_name[i] = '\0';
   
   string fpath;
   // don't allow escape from data directory
   if(!strchr(f_name, '\\') && !strchr(f_name, '/')) {
      fpath = dir + f_name;
      if((file=fopen_file(fpath)))
         return file;
   }
   string alt_names_file = "alt_names.txt";
   fpath = dir + alt_names_file;
   FILE *alt = fopen_file(fpath);
   if(alt) {
      string name = find_alt_name(alt, f_name);
      fclose(alt);
      if(name != "" && !strchr(name.c_str(), '\\') && !strchr(name.c_str(), '/')) {
         fpath = dir + name;
         if((file=fopen_file(fpath)))
            return file;
      }
      *aname = name;    // may be used to open an internal resource
   }
  
   return 0;
}


FILE *open_sup_file(const char *fname, const char *subdir, string *alt_name,
      int *where, string *fpath)
{
   string alt_nam;
   if(!alt_name)
      alt_name = &alt_nam;

   int whr;
   if(!where)
      where = &whr;
   
   string fpth;
   if(!fpath)
      fpath = &fpth;
   
   string aname;
   FILE *file;
   // try to open fname alone
   *fpath = fname;
   *where = 0; // local
   if((file=fopen_file(*fpath)))
      return file;
   
   // The file open will fail, but will read alt_name.txt in local directory
   if((file = open_file_data("", *fpath, alt_name)))
      return file;
   if(*alt_name != "")
      return 0;

   
   // try environment variable for data directory
   char *sup_dir = getenv("ANTIPRISM_DATA");
   if(sup_dir) {
      *where = 1; // environment
      string fdir = string(sup_dir) + subdir;
      if((file=open_file_data(fdir, *fpath, alt_name)))
         return file;
      if(*alt_name != "")
         return 0;
   }
   
   // try hardcoded install path for data directory
   *where = 2; // installed
   string fdir = string(SUPDIR) + subdir;
   if((file=open_file_data(fdir, *fpath, alt_name)))
      return file;
   if(*alt_name != "")
      return 0;
   
   return 0;
}


string msg_str(const char *fmt, ...)
{
   char message[MSG_SZ];
   va_list args;
   va_start(args, fmt);
   vsnprintf(message, MSG_SZ-1, fmt, args);
   return message;
}
