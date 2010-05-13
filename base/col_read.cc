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

/*\file col_val_read.cc
 * \brief Input processing for col_val
*/

#include <ctype.h>
#include <string.h>

#include "utils.h"
#include "col_val.h"
#include "named_cols.h"

bool col_val::from_offvals(vector<char *> &vals, char *errmsg, int *type)
{
   unset();
   int dummy_type;
   if(!type)
      type = &dummy_type;
   
   *type = -1;
   int i;
   switch(vals.size()) {
      case 0:
         *type = 0;
         break;
      case 1:
         if(read_int(vals[0], &i) && i>-1) {
            index = i;
            *type = 1;
         }
         else
            if(errmsg)
               snprintf(errmsg, MSG_SZ,
                     "colour index '%s' is not a valid colour index", vals[0]);
         break;
      case 3:
      case 4:  // includes alpha value
         if(read_intvals(vals, errmsg))
            *type=3+(vals.size()==4);
         else if(read_decvals(vals, errmsg))
            *type=5+(vals.size()==4);
         break;
      default:
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "incorrect number of colour values");
         
   }
   return (*type==0 || is_set());
}


bool col_val::read_intvals(vector<char *> &vals, char *errmsg)
{
   unset();
   vector<int> ivals;
   if(!read_int_list(vals, ivals, errmsg, false))
      return false;
   return from_intvals(ivals, errmsg);
}


bool col_val::read_intvals(char *str, char *errmsg)
{
   unset();
   vector<int> vals;
   if(!read_int_list(str, vals, errmsg, false))
      return false;
   return from_intvals(vals, errmsg);
}


bool col_val::from_intvals(vector<int> &vals, char *errmsg)
{
   unset();
   if(vals.size()<3 || vals.size()>4) {
      if(errmsg)
         snprintf(errmsg, MSG_SZ,
             "integer format, %lu numbers given, must give 3 or 4",
             (unsigned long)vals.size());
      return false;
      }
   
   if(vals.size()==3)
      vals.push_back(255);
   
   for(unsigned int i=0; i<4; i++) {
      if(vals[i]<0 || vals[i]>255) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ,
                  "integer format, \"%d\" is not in range 0-255", vals[i]);
         return false;
      }
   }
   
   set_rgba(vals[0], vals[1], vals[2], vals[3]);
   return true;
}
      
bool col_val::read_decvals(char *str, char *errmsg)
{
   unset();
   vector<double> vals;
   if(!read_double_list(str, vals, errmsg))
      return false;
   return from_decvals(vals, errmsg);
}


bool col_val::read_decvals(vector<char *> &vals, char *errmsg)
{
   unset();
   vector<double> dvals;
   if(!read_double_list(vals, dvals, errmsg))
      return false;
   return from_decvals(dvals, errmsg);
}



bool col_val::from_decvals(vector<double> &vals, char *errmsg)
{
   unset();
   if(vals.size()<3 || vals.size()>4) {
      if(errmsg)
         snprintf(errmsg, MSG_SZ,
               "decimal format, %lu numbers given, must give 3",
               (unsigned long) vals.size());
      return false;
   }
   
   if(vals.size()==3)
      vals.push_back(1.0);
   
   for(unsigned int i=0; i<4; i++) {
      if(vals[i]<0 || vals[i]>1) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ,
                  "decimal format, \"%g\" is not in range 0.0-1.0", vals[i]);
         return false;
      }
   }
   set_rgba(vals[0], vals[1], vals[2], vals[3]);
   return true;   
}


bool col_val::read_hexvals(char *str, char *errmsg)
{
   unset();
   char hexstr[MSG_SZ];
   strncpy(hexstr, str, MSG_SZ-1);
   hexstr[MSG_SZ-1] = '\0';
   if(!strchr("Xx#", *hexstr)) {
      if(errmsg)
         snprintf(errmsg, MSG_SZ,
               "hex format, first character is not X, x, or #");
      return false;
   }
   if(strlen(hexstr)==1)
      strcat(hexstr, "00000000");
   else if(strlen(hexstr)==7)
      strcat(hexstr, "FF");
   if(strlen(hexstr)!=9 || strspn(hexstr+1, "0123456789aAbBcCdDeEfF")!=8) {
      if(errmsg)
         snprintf(errmsg, MSG_SZ, "hex format, %c is not followed by "
               "6 or 8 hexadecimal digits", *hexstr);
      return false;
   }
   unsigned int val;
   sscanf(hexstr, "%*c%8x", &val);
   set_rgba(int(val/(256*256*256)), int((val/(256*256))%256),
         int((val/256)%256), int(val%256) );
   return true;
}



bool col_val::read_hsva_vals(char *str, char *errmsg)
{
   bool ret = false;
   if((str[0] == 'h' || str[0] == 'H') && strlen(str) > 1) {
      char hsva_str[MSG_SZ];
      strncpy(hsva_str, str+1, MSG_SZ-1);
      hsva_str[MSG_SZ-1] = '\0';
      vector<double> vals;
      bool hue_degrees = (str[0] == 'h');
      if(read_double_list(hsva_str, vals, errmsg, 4)) {
         int sz = vals.size();
         double hue = vals[0];
         if (hue_degrees)
            hue /= 360.0;
         for(int i=1;i<sz; i++) {
            if(vals[i]<0 || vals[i]>1) {
               if(errmsg)
                  snprintf(errmsg, MSG_SZ,
                        "hsva format, \"%g\" is not in range 0.0-1.0", vals[i]);
               return ret;
            }
         }
         set_hsva(hue,(sz < 2 ? 1.0 : vals[1]),(sz < 3 ? 1.0 : vals[2]),(sz < 4 ? 1.0 : vals[3]));
         ret = true;
      }
   }
   return ret;
}
   

bool col_val::read_colorname(char *str, char *errmsg, bool as_index)
{
   unset();
   // Colour name 'none' leaves the colour unset
   if(strcmp(str, "none")==0 && !as_index)
      return true;
   
   // strip whitespace, convert to lowercase
   char *p = str;
   char *p_to = str;
   while(*p) {
      if(!isspace(*p)) {
        if(isupper(*p))
         *p = tolower(*p);
        *p_to++ = *p;
      }
      *p++;
   }
   *p_to = '\0';
   // change grey to gray
   char *pgrey = strstr(str, "grey");
   if(pgrey)
      pgrey[2] = 'a';
  
   if(strcmp(str, "invisible")==0) {
      if(!as_index)
         *this = invisible;
   }
   else {
      // inefficient, not likely to be called much
      for(unsigned int i=0; *named_colors[i].name; i++) {
         if(strcmp(str, named_colors[i].name)==0) {
            if(as_index)
               set_idx(i);
            else {
               set_rgba(named_colors[i].r,named_colors[i].g,named_colors[i].b);
               break;
            }
         }
      }
   }
   
   if(errmsg && !is_set())
      snprintf(errmsg, MSG_SZ, "unknown colour name '%s'", str);

   return is_set();     
}
      

bool col_val::read(char *str, char *errmsg)
{
   unset();
   char str2[MSG_SZ];
   strncpy(str2, str, MSG_SZ);
   str2[MSG_SZ-1] = '\0';

   clear_extra_whitespace(str2);
   if(!*str2) // don't interpret whitespace only as index 0
      return false;
   
   string orig = str2;
   
   vector<char *> vals;
   int split_ret, typ;
   if(strchr(str2, ',')!=0)
      split_ret = split_line(str2, vals, ",");
   else
      split_ret = split_line(str2, vals);
   
   from_offvals(vals, errmsg, &typ);
   if(is_set())
      return true;
   
            
   strcpy(str2, orig.c_str());
   if(read_hexvals(str2, errmsg))
      return true;
      
   if(read_hsva_vals(str2, errmsg))
      return true;

   strcpy(str2, orig.c_str());
   if(read_colorname(str2, errmsg))
      return true;
 
   if(errmsg)
      snprintf(errmsg, MSG_SZ, "unknown colour name or format '%s'", str);
   return false;
}


