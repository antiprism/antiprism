/*
   Copyright (c) 2008, Adrian Rossiter

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

#include <ctype.h>
#include <string.h>
#include <map>

#include "geom.h"
#include "coloring.h"
#include "utils.h"
#include "rand_gen.h"

#include "named_cols.h"

using std::map;

///Whitespace characters
const char WHITESPACE[] = " \t\r\n\f\v";


bool color_map::init(const char *params, char *errmsg)
{
   wrap = 0;
   step = 1;
   shift = 0;

   if(!params)            //treat null pointer like the empty string
      return true;

   char prms[MSG_SZ];
   strncpy(prms, params, MSG_SZ);
   prms[MSG_SZ-1] = '\0';

   char errmsg2[MSG_SZ];
   
   char *p = strchr(prms, '%');
   if(p) {
      if(*(p+1)=='\0')
         wrap=-1;
      else if(!read_int(p+1, &wrap, errmsg2)) {
         snprintf(errmsg, MSG_SZ, "wrap: %s", errmsg2);
         return false;
      }
      *p = '\0';
   }

   p = strchr(prms, '*');
   if(p) {
      if(!read_int(p+1, &step, errmsg2)) {
         snprintf(errmsg, MSG_SZ, "step size: %s", errmsg2);
         return false;
      }
      *p = '\0';
   }

   p = strchr(prms, '+');
   if(p) {
      if(!read_int(p+1, &shift, errmsg2)) {
         snprintf(errmsg, MSG_SZ, "shift size: %s", errmsg2);
         return false;
      }
      *p = '\0';
   }

   return true;
}

void color_map::copy_params(const color_map &cmap)
{
   shift = cmap.shift;
   step = cmap.step;
   wrap = cmap.wrap;
}

bool color_map::init_strip(char *map_name, char *errmsg)
{
   if(!color_map::init(map_name, errmsg))
      return false;
   size_t name_len = strcspn(map_name, "+*%");
   map_name[name_len] = '\0';
   return true;
}

static color_map* init_color_map_generated(const char *map_name, char *errmsg=0)
{
   char errmsg_tmp[MSG_SZ];
   if(!errmsg)
      errmsg = errmsg_tmp;
   strcpy(errmsg, "name not found");  // default error message

   char name[MSG_SZ];
   strncpy(name, map_name, MSG_SZ);
   name[MSG_SZ-1] = '\0';
   size_t name_len = strspn(name, "abcdefghijklmnopqrstuvwxyz");
   name[name_len] = '\0';

   color_map *cmap = 0;
     
   if(*map_name=='\0' || strcmp(name, "null")==0) {  // A null map
      cmap = new color_map_map();
   }
   
   else if(strcmp(name, "rnd")==0 ||
           strcmp(name, "rand")==0 ||
           strcmp(name, "random")==0  ) {
      cmap = new color_map_range_rand_hsv();
      if(cmap && !cmap->init(map_name+name_len, errmsg)) {
         delete cmap;
         cmap = 0;
      }
      if(!cmap) {
         cmap = new color_map_range_rand_rgb();
         if(cmap && !cmap->init(map_name+name_len, errmsg)) {
            delete cmap;
            cmap = 0;
         }
      }
   }

   else if(strcmp(name, "rng")==0 ||
           strcmp(name, "range")==0 ) {
      cmap = new color_map_range_hsv();
      if(cmap && !cmap->init(map_name+name_len, errmsg)) {
         delete cmap;
         cmap = 0;
      }
      if(!cmap) {
         cmap = new color_map_range_rgb();
         if(cmap && !cmap->init(map_name+name_len, errmsg)) {
            delete cmap;
            cmap = 0;
         }
      }
   }

   else if(strcmp(name, "spread")==0) {
      cmap = new color_map_spread();
      if(cmap && !cmap->init(map_name+name_len, errmsg)) {
         delete cmap;
         cmap = 0;
      }
   }

   else if(strcmp(name, "grey")==0 ||
           strcmp(name, "greyw")==0 ||
           strcmp(name, "gray")==0  ||
           strcmp(name, "grayw")==0 ) {
      if(strspn(map_name+name_len,"0123456789-+*%")==strlen(map_name+name_len)){
         bool wrp = (name[strlen(name)-1]=='w');
         char grey_name[MSG_SZ];
         size_t num_dgts = strspn(map_name+name_len,"0123456789");
         strncpy(grey_name, map_name+name_len, num_dgts);
         snprintf(grey_name+num_dgts, MSG_SZ-num_dgts-1,
               "_H0S0V0:1%s%s", wrp?":0":"", map_name+name_len+num_dgts);
         cmap = new color_map_range_hsv();
         if(cmap && !cmap->init(grey_name, errmsg)) {
            delete cmap;
            cmap = 0;
         }
      }
   }

   else if(strcmp(name, "uniform")==0) {
      color_map_multi *multi = new color_map_multi;
      color_map_map *overrides = new color_map_map;
      color_map *spread_map = init_color_map("spread+53*12");
      if(multi && overrides && spread_map && multi->init(map_name, errmsg)) {
         overrides->set_col(60, col_val(0.9,0.45,0.0)); // triangle
         overrides->set_col(36, col_val(0.7,0.1,0.2));  // pentagram
         multi->add_cmap(overrides);
         multi->add_cmap(spread_map);
         cmap = multi;
      }
      else {
         delete cmap;
         delete overrides;
         delete spread_map;
         cmap = 0;
      }
   }

   else if(strcmp(name, "compound")==0) {
      color_map_multi *multi = new color_map_multi();
      color_map_map *overrides = new color_map_map;
      color_map *spread_map = init_color_map("spread+2");
      if(multi && overrides && spread_map && multi->init(map_name, errmsg)) {
         overrides->set_col(0, col_val(1.0,0.5,0.0)); // orange
         overrides->set_col(1, col_val(1.0,1.0,0.0)); // yellow
         overrides->set_col(2, col_val(1.0,0.0,0.0)); // red
         overrides->set_col(3, col_val(0.0,0.5,0.0)); // green
         overrides->set_col(4, col_val(0.0,0.0,1.0)); // blue
         overrides->set_col(5, col_val(1.0,0.0,1.0)); // magnenta
         overrides->set_col(6, col_val(0.0,1.0,1.0)); // cyan
         multi->add_cmap(overrides);
         multi->add_cmap(spread_map);
         cmap = multi;
      }
      else {
         delete cmap;
         delete overrides;
         delete spread_map;
         cmap = 0;
      }
   }

   else if(strcmp(name, "remap")==0) {
      cmap = new color_map_remap;
   }
   
   else if(strcmp(name, "mkmap")==0) {
      if(map_name[name_len]=='_') {
         color_map_map *cmm = new color_map_map;
         if(cmm && cmm->init_from_line(map_name+name_len+1, errmsg))
            cmap = cmm;
      }
   }

   return cmap;
}


color_map* init_color_map(const char *map_name, char *errmsg)
{
   color_map *cmap = 0;
   
   char name[MSG_SZ];
   strncpy(name, map_name, MSG_SZ);
   name[MSG_SZ-1] = '\0';
   size_t name_len = strcspn(name, "+*%");
   name[name_len] = '\0';

   char errmsg2[MSG_SZ];
   string alt_name;
   FILE *cfile = open_sup_file(name, "/col_maps/", &alt_name);
   if(alt_name!="") {  // an alt name found before a file with the name
      cmap = init_color_map_generated(alt_name.c_str(), errmsg2);
      if(errmsg) {
         if(!cmap)
            snprintf(errmsg, MSG_SZ, "could not open colour map file"
               " \'%s=%s\'", map_name, alt_name.c_str());
         else
            strcpy(errmsg, errmsg2);
      }
   }
   else if(cfile) {
      cmap = new color_map_map();
      if(cmap && !cmap->init(map_name, errmsg2)) {
         delete cmap;
         cmap = 0;
      }
      if(errmsg)
         strcpy(errmsg, errmsg2);
   }
   else {
      cmap = init_color_map_generated(map_name, errmsg2);
      if(errmsg) {
         if(cmap)
            strcpy(errmsg, errmsg2);
         else
            snprintf(errmsg, MSG_SZ,
                  "could not open colour map file \'%s\': %s",
                  map_name, errmsg2);
      }
   }

   return cmap;
}


void color_map::set_wrap(int wrp)
{
   if(wrp<0)
      wrp=effective_size();
   wrap = (wrp<0) ? 0 : wrp;
}


color_map *color_map_map::get_condensed()
{
   color_map_map *colmap = new color_map_map();
   if(!colmap)
      return 0;

   int idx=0;
   map<int, col_val>::iterator mi;
   for(mi=cmap.begin(); mi!=cmap.end(); mi++)
      colmap->cmap[idx++] = (mi->second);

   return colmap;
}

//-------------------------------------------------------------------
//color_map_range

static double interpolate(int num, int map_sz, vector<double> vals)
{
   int interval_sz = (int)vals.size()-1;
   if(interval_sz<1)
      return vals[0];

   double pos = (double)num/map_sz;
   int low_idx = (int)floor(interval_sz*pos);
   int high_idx = (int)ceil(interval_sz*pos);
   double low_frac = (double)low_idx/interval_sz;
   double high_frac = (double)high_idx/interval_sz;
   double frac = (high_idx==low_idx) ? 0 : (pos-low_frac)/(high_frac-low_frac);
   double val = vals[low_idx] + frac*(vals[high_idx]-vals[low_idx]);
   return val;
}


col_val color_map_range::get_col(int idx)
{
   col_val col;
   int eff_idx = get_effective_index(idx);
   if(eff_idx<map_sz) {
      double comp[4];
      for(int i=0; i<4; i++)
         comp[i] = interpolate(get_effective_index(idx), map_sz, ranges[i]);
      (col.*set_func)(comp[0], comp[1], comp[2], comp[3]);
   }
   return col;
}


bool color_map_range::init(const char *map_name, char *errmsg)
{
   if(errmsg)
      *errmsg = '\0';
   char name[MSG_SZ];
   strncpy(name, map_name, MSG_SZ-1);
   name[MSG_SZ-1] = '\0';

   if(!init_strip(name, errmsg))
      return false;

   vector<char *> vals;
   split_line(name, vals, "_");
   //for(unsigned int i=0; i<vals.size(); i++)
   //   fprintf(stderr, "vals[%d] = '%s'\n", i, vals[i]);

   if(vals.size() > 2) {
      if(errmsg)
         sprintf(errmsg, "map_name contains more than one '_'");
      return false;
   }

   // Get the map size
   char errmsg2[MSG_SZ];
   if(*map_name != '_') {
      if(vals.size() && !read_int(vals[0], &map_sz, errmsg2)) {
         if(errmsg)
            sprintf(errmsg, "map size: %s", errmsg2);
         return false;
      }
   }
   
   if(get_wrap()==-1)
      set_wrap(effective_size());
      
   if(*map_name != '_' && vals.size()<2) // A size was given but no comp ranges
         return true;

   if(strspn(vals.back(), "HhSsVv") && strspn(vals.back(), "RrGgBb")) {
      if(errmsg)
         sprintf(errmsg, "cannot include both RGB and HSV components");
      return false;
   }
   size_t char_cnt;
   bool is_hsv = (set_func==((void (col_val::*)(double, double, double, double))&col_val::set_hsva));
   if( ( is_hsv && (char_cnt=strspn(vals.back(), "RrGgBb"))) ||
       (!is_hsv && (char_cnt=strspn(vals.back(), "HhSsVv"))) ) {
      if(errmsg)
         sprintf(errmsg, "invalid component letter '%c'",
               *(vals.back()+char_cnt-1));
      return false;
   }
   
   int rng_len = strlen(vals.back());
   char rngs[MSG_SZ];
   char *q = rngs;
   int cur_idx = -1;
   char cur_comp = 'X';
   for(const char *p=vals.back(); p-vals.back()<rng_len+1; p++) {
      //fprintf(stderr, "*p = %c\n", *p);
      if(strchr("HhSsVvAaRrGgBb", *p) || cur_idx<0 || *p=='\0') {
         *q = '\0';
         if(cur_idx>=0) {
            if(read_double_list(rngs, ranges[cur_idx], errmsg2, 0, ":")) {
               if(ranges[cur_idx].size()==0) {
                  if(errmsg)
                     sprintf(errmsg, "component letter '%c' isn't followed "
                           "by any values", cur_comp);
                  return false;
               }

               for(unsigned int j=0; j<ranges[cur_idx].size(); j++) {
                  if(ranges[cur_idx][j]<0) {
                     if(errmsg)
                        sprintf(errmsg, "component letter '%c' contains a "
                              "negative value", cur_comp);
                     return false;
                  }
                  if(cur_comp=='h')
                     ranges[cur_idx][j] /= 360;   // h is hue in range 0-360 
               }
            }
            else {
               if(errmsg)
                  sprintf(errmsg, "component letter '%c': %s", cur_comp, errmsg2);
               return false;
            }
         }
         if(strchr("HhRr", *p))
            cur_idx = 0;
         else if(strchr("SsGg", *p))
            cur_idx = 1;
         else if(strchr("VvBb", *p))
            cur_idx = 2;
         else if(strchr("Aa", *p))
            cur_idx = 3;
         else {
            if(errmsg)
               sprintf(errmsg, "invalid component letter '%c'", *p);
            return false;
         }
         cur_comp = *p;
         q = rngs;
      }
      else if(!(isdigit(*p) || *p == '.' || *p == ':')) {
         if(errmsg)
            sprintf(errmsg, "invalid component letter '%c'", *p);
         return false;
      }
      else if(!isspace(*p)) {
         *q++ = *p;
      }
   }

   //for(int i=0; i<4; i++)
   //   for(unsigned int j=0; j<ranges[i].size(); j++)
   //      fprintf(stderr, "ranges[%d][%u] = %g\n", i, j, ranges[i][j]);
   
   return true;
}


bool color_map_range::set_range(int idx, vector<double> range)
{
   if(idx<0 || idx>3 || !range.size())
      return false;
   ranges[idx] = range;
   return true;
}



bool color_map_range_hsv::init(const char *map_name, char *errmsg)
{
   set_func = &col_val::set_hsva;
   ranges[0].push_back(0);
   ranges[0].push_back(1);
   ranges[1].push_back(0.9);
   ranges[2].push_back(0.9);
   ranges[3].push_back(1);
   
   set_map_sz(256);
   int ret = color_map_range::init(map_name, errmsg);
   
   return ret;
}


static double rand_in_range(vector<double> rng, int seed)
{
   if(rng.size()>1) {
      rand_gen rnd((seed+130));
      rnd.seedi((rnd.ranlui()*rnd.ranlui()*rnd.ranlui())&0xFFFFFFFF);
      const double val = fmod(rng[0] + (rng[1]-rng[0])*rnd.ranf(), 1+epsilon);
      return val;
   }
   else
      return rng[0];
}


bool color_map_range_rgb::init(const char *map_name, char *errmsg)
{
   set_func = &col_val::set_rgba;
   ranges[0].push_back(0.3);
   ranges[0].push_back(1);
   ranges[1].push_back(0.3);
   ranges[1].push_back(1);
   ranges[2].push_back(0.3);
   ranges[2].push_back(1);
   ranges[3].push_back(1);
   
   set_map_sz(256);
   int ret = color_map_range::init(map_name, errmsg);
 
   return ret;
}


bool color_map_range_rand_hsv::init(const char *map_name, char *errmsg)
{
   set_func = &col_val::set_hsva;
   ranges[0].push_back(0);
   ranges[0].push_back(1);
   ranges[1].push_back(0.7);
   ranges[1].push_back(1);
   ranges[2].push_back(0.7);
   ranges[2].push_back(1);
   ranges[3].push_back(1);
   
   set_map_sz(-1);
   int ret = color_map_range::init(map_name, errmsg);
   
   return ret;
}


bool color_map_range_rand_rgb::init(const char *map_name, char *errmsg)
{
   set_func = &col_val::set_rgba;
   ranges[0].push_back(0.3);
   ranges[0].push_back(1);
   ranges[1].push_back(0.3);
   ranges[1].push_back(1);
   ranges[2].push_back(0.3);
   ranges[2].push_back(1);
   ranges[3].push_back(1);
   
   set_map_sz(-1);
   int ret = color_map_range::init(map_name, errmsg);
   return ret;
}

col_val color_map_range_rand::get_col(int idx)
{
   col_val col;
   idx = get_effective_index(idx);
   if(get_wrap() || get_map_sz()==-1 || idx<get_map_sz())
      (col.*set_func)( rand_in_range(ranges[0], idx*1),
                       rand_in_range(ranges[1], idx*2),
                       rand_in_range(ranges[2], idx*3),
                       rand_in_range(ranges[3], idx*4) );
   return col;
}


col_val color_map_spread::get_col(int idx)
{
   int eff_idx = get_effective_index(idx);
   int num_entries = 1024;
   int num_intervals = 4;
   //int step_by = 29+196;
   //int step_by = 37+128;
   int step_by = 53;
   int entries_per_interval = num_entries/num_intervals;

   eff_idx = ((long)eff_idx*step_by) % num_entries;
   int interval = eff_idx/entries_per_interval;
   float H = (float)eff_idx/entries_per_interval - interval;
   float S = 0.0;
   float V = 0.0;
   switch(interval) {
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

   col_val col;
   col.set_hsva(H, S, V);

   return col;
}


//-------------------------------------------------------------------
//color_map_map



static bool parse_gimp_file(FILE *cfile, map<int, col_val> *cmap,char *errmsg=0)
{
   const int line_size=1024;
   char line[line_size];
   char buf[line_size];

   int stage = 0;

   if(errmsg)
      *errmsg='\0';

   int idx_no = 0;
   int line_no = 0;
   while(fgets(line, line_size, cfile)) {
      line_no++;
     
      // ignore comments
      char *first_hash = strchr(line, '#');
      if(first_hash)
         *first_hash = '\0';

      // skip blank lines
      if(sscanf(line, " %s", buf)==EOF)
         continue;

      if(stage == 0) {
         stage++;
         // ignore header line
         if(strncasecmp(line, "GIMP Palette", 12)==0)
            continue;
      }
      if(stage == 1) {
         stage++;
         // ignore header line
         if(strncasecmp(line, "Name", 4)==0)
            continue;
      }

      if(stage == 2) {
         stage++;
         // ignore header line
         if(strncasecmp(line, "Columns", 7)==0)
            continue;
      }

      // stage == 3
      char *r = strtok(line, WHITESPACE);
      char *g = (r) ? strtok(NULL, WHITESPACE) : 0;
      char *b = (g) ? strtok(NULL, WHITESPACE) : 0;
      //char *name = (b) ? strtok(NULL, WHITESPACE) : 0;

      if(!b) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "gimp colour map: line %d: not enough colour values", line_no);
         return false;
      }

      sprintf(buf, "%s %s %s", r, g, b);
      col_val col;
      col.read(buf, errmsg);
      if(!col.is_set()) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "gimp colour map: line %d: invalid colour '%s'", line_no, buf);
         return false;
      }

      //col.get_val().dump();
      (*cmap)[idx_no++] = col;
   }
   
   return true;
}


static bool parse_file(FILE *cfile, map<int, col_val> *cmap, char *errmsg=0)
{
   const int line_size=1024;
   char line[line_size];
   
   if(errmsg)
      *errmsg='\0';

   int line_no = 0;
   int next_idx = 0;   // index to use if no map index is given with a colour
   while(fgets(line, line_size, cfile)) {
      line_no++;
     
      // copy the map entry string
      char entry[MSG_SZ]; 
      strncpy(entry, line, MSG_SZ);
      entry[MSG_SZ-1] = '\0';
      
      // ignore comments
      char *first_hash = strchr(entry, '#');
      if(first_hash)
         *first_hash = '\0';

      // skip blank lines
      char c;
      if(sscanf(entry, " %c", &c)==EOF)
         continue;
 
      if(!cmap->size() && strncasecmp(entry, "GIMP Palette", 12)==0) {
         rewind(cfile);
         return parse_gimp_file(cfile, cmap, errmsg);
      }
            

      char *col_pos = entry;
      char *eq_pos = strchr(entry, '=');
      if(eq_pos) {
         col_pos = eq_pos + 1;
         if(strchr(col_pos, '=')) {
            if(errmsg)
               snprintf(errmsg, MSG_SZ, "colour map: line %d: more than one =, '%s'", line_no, line);
            return false;
         }
         *eq_pos = '\0';
         if(!read_int(entry, &next_idx) || next_idx<0) {
            if(errmsg)
               snprintf(errmsg, MSG_SZ, "colour map: line %d: invalid index number, '%s'", line_no, entry);
            return false;
         }
      }

      col_val col;
      col.read(col_pos);
      if(!col.is_set()) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "colour map: line %d: invalid colour, '%s'", line_no, col_pos);
         return false;
      }

      if(errmsg && !*errmsg && cmap->find(next_idx)!=cmap->end())
         snprintf(errmsg, MSG_SZ, "colour map: line %d: mapping for index %d is being overwritten", line_no, next_idx);

      (*cmap)[next_idx++] = col;
   }
   
   return true;

} 
      
static bool parse_map_from_line(const char *line, map<int, col_val> *cmap,
      char *errmsg=0)
{
   if(errmsg)
      *errmsg='\0';

   // copy the map string so the original will not be modified
   char *str = new char[strlen(line)+1];
   strcpy(str, line);
   
   vector<char *> entries;
   split_line(str, entries, ":");
   int next_idx = 0;   // index to use if no map index is given with a colour
   bool cmap_ok = true;
   for(unsigned int i=0; i<entries.size(); i++) {
      // copy the map entry string
      char entry[MSG_SZ]; 
      strncpy(entry, entries[i], MSG_SZ);
      entry[MSG_SZ-1] = '\0';

      // ignore comments
      char *first_hash = strchr(entry, '#');
      if(first_hash)
         *first_hash = '\0';

      // skip blank lines
      char c;
      if(sscanf(entry, " %c", &c)==EOF)
         continue;
      
      char *col_pos = entry;
      char *eq_pos = strchr(entry, '=');
      if(eq_pos) {
         col_pos = eq_pos + 1;
         if(strchr(col_pos, '=')) {
            if(errmsg)
               snprintf(errmsg, MSG_SZ, "entry %d: more than one =, '%s'", i+1, entries[i]);
            cmap_ok = false;
            break;
         }
         *eq_pos = '\0';
         if(!read_int(entry, &next_idx) || next_idx<0) {
            if(errmsg)
               snprintf(errmsg, MSG_SZ, "entry %d: invalid index number, '%s'", i+1, entry);
            cmap_ok = false;
            break;
         }
      }

      // Allow '' as a number separator
      for(char *p=col_pos; *p; p++)
         if(*p == '/')
            *p = ' ';

      col_val col;
      col.read(col_pos);
      if(!col.is_set()) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "entry %d: invalid colour, '%s'", i+1, col_pos);
         cmap_ok = false;
      }

      if(errmsg && !*errmsg && cmap->find(next_idx)!=cmap->end())
         snprintf(errmsg, MSG_SZ, "entry %d: mapping for index %d is being overwritten", i, next_idx);

      (*cmap)[next_idx++] = col;
   }

   delete[] str;
   return cmap_ok;
} 


col_val color_map_map::get_col(int idx)
{
   map<int, col_val>::iterator mi_idx = cmap.find(get_effective_index(idx));
   if(mi_idx!=cmap.end())  // index is in colour map
      return mi_idx->second;
   else
      return col_val();
}


void color_map_map::set_col(int idx, col_val col)
{
   if(col.is_set())
      cmap[idx] = col;
   else
      cmap.erase(idx);
}  



int color_map_map::effective_size() const
{
   return size() ? cmap.rbegin()->first+1 : 0;
}


bool color_map_map::init(const char *map_name, char *errmsg)
{
   char name[MSG_SZ];
   strncpy(name, map_name, MSG_SZ);
   name[MSG_SZ-1] = '\0';
   
   if(!init_strip(name, errmsg))
      return false;

   bool cmap_ok = false;
   string alt_name;
   FILE *cfile = open_sup_file(name, "/col_maps/", &alt_name);
   if(cfile) {
      char errmsg2[MSG_SZ];
      cmap_ok = parse_file(cfile, &cmap, errmsg2);
      if(/*!cmap_ok &&*/ errmsg)
         strcpy(errmsg, errmsg2);
   }
   if(get_wrap()==-1)
      set_wrap(effective_size());
   
   if(cfile)
      fclose(cfile);

   return cmap_ok;
}

bool color_map_map::init_from_line(const char *map_name, char *errmsg)
{
   char name[MSG_SZ];
   strncpy(name, map_name, MSG_SZ);
   name[MSG_SZ-1] = '\0';
   
   if(!init_strip(name, errmsg))
      return false;

   bool cmap_ok = parse_map_from_line(name, &cmap, errmsg);
   if(get_wrap()==-1)
      set_wrap(effective_size());

   return cmap_ok;
}

void color_map_map::read_named_colors()
{
   col_val col;
   for(int i=0; *named_colors[i].name; i++) {
      col.set_rgba(named_colors[i].r, named_colors[i].g, named_colors[i].b);
      cmap[i] = col;
   }
}



//-------------------------------------------------------------------
//color_map_multi

bool color_map_multi::init(const char *map_name, char *errmsg)
{
   if(!color_map::init(map_name, errmsg))
      return false;
   
   if(get_wrap()==-1)
      set_wrap(0);
   
   return true;
}


color_map_multi::~color_map_multi()
{
   for(unsigned int i=0; i<cmaps.size(); i++)
      delete cmaps[i];
}

color_map_multi::color_map_multi(const color_map_multi &cmap) : color_map(cmap)
{
   copy_params(cmap);
   map_sz = cmap.map_sz;
   for(unsigned int i=0; i<cmap.cmaps.size(); i++)
      add_cmap(cmap.cmaps[i]->clone());
}

color_map_multi &color_map_multi::operator=(const color_map_multi &cmap)
{
   if(this!=&cmap) {
      copy_params(cmap);
      map_sz = cmap.map_sz;
      while(cmaps.size())
         del_cmap();
      for(unsigned int i=0; i<cmap.cmaps.size(); i++)
         add_cmap(cmap.cmaps[i]->clone());
   }
   return *this;
}


void color_map_multi::set_map_sz()
{
   map_sz = 0;
   for(unsigned int i=0; i<cmaps.size(); i++)
      if(cmaps[i]->effective_size()>map_sz)
         map_sz = cmaps[i]->effective_size();
}


void color_map_multi::add_cmap(color_map *col_map, unsigned int pos)
{
   vector<color_map *>::iterator mi;
   if(pos>=cmaps.size())
      mi = cmaps.end();
   else
      mi = cmaps.begin()+pos;
   cmaps.insert(mi, col_map);
   if(col_map->effective_size()>map_sz)
      map_sz = col_map->effective_size();
}


void color_map_multi::del_cmap(unsigned int pos)
{
   if(cmaps.size()) {
      vector<color_map *>::iterator mi;
      if(pos>=cmaps.size())
         mi = cmaps.end()-1;
      else
         mi = cmaps.begin()+pos;
      delete *mi;
      cmaps.erase(mi);
      set_map_sz();
   }
}


col_val color_map_multi::get_col(int idx)
{
   col_val col;
   int eff_idx = get_effective_index(idx);
   for(unsigned int i=0; i<cmaps.size(); i++) {
      col = cmaps[i]->get_col(eff_idx);
      if(col.is_set())
         break;
   }

   return col.is_set() ? col : col_val(idx);
}


/*
bool color_map_range::init(const char *map_name, char *errmsg)
{
   char name[MSG_SZ];
   strncpy(name, map_name, MSG_SZ-1);
   name[MSG_SZ-1] = '\0';

   if(!init_strip(name, errmsg))
      return false;
   
   vector<char *> vals;
   split_line(name, vals, "_");
   if(vals.size() > 5)
      return false;

   // Followed by the map size
   int map_sz=-1;
   int first_comp_idx = 0;
   if(*map_name != '_') {
      if(vals.size() && !read_int(vals[0], &map_sz))
         return false;
      first_comp_idx = 1;
   }


   // Followed by lists of values for each of the four components
   for(unsigned int i=first_comp_idx; i<vals.size(); i++) {
      int comp_idx = i-first_comp_idx;
      if(!read_double_list(vals[i], ranges[comp_idx], errmsg, 0, ":"))
         return false;
      for(unsigned int j=0; j<ranges[comp_idx].size(); j++)
         if(ranges[comp_idx][j]<0) {
            if(errmsg)
               sprintf(errmsg, "component %d contains a negative value",
                     comp_idx+1);
            return false;
         }
   }
   set_map_sz((map_sz>0) ? map_sz : 256);
   
   return true;
}
*/

 
/*
bool color_proc_torange::init(const char *map_name, char *errmsg)
{
   char name[MSG_SZ];
   strncpy(name, map_name, MSG_SZ-1);
   name[MSG_SZ-1] = '\0';

   vector<char *> vals;
   split_line(name, vals, "_");
   if(vals.size() > 4)
      return false;

   // lists of values for each of the four components
   for(unsigned int i=0; i<vals.size(); i++) {
      if(!read_double_list(vals[i], ranges[i], errmsg, 0, ":"))
         return false;
      for(unsigned int j=0; j<ranges[i].size(); j++)
         if(ranges[i][j]<0) {
            if(errmsg)
               sprintf(errmsg, "component %d contains a negative value", i+1);
            return false;
         }
      if(ranges[i].size()>2) {
         if(errmsg)
            sprintf(errmsg, "component %d contains a negative value", i+1);
         return false;
      }
      if(ranges[i].size()==1)
         ranges[i].push_back(ranges[i][0]);
   }
   
   return true;
}

*/
