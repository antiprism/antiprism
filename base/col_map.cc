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

bool color_map::init_strip(char *map_name, char *errmsg)
{
   if(!color_map::init(map_name, errmsg))
      return false;
   size_t name_len = strcspn(map_name, "+*%");
   map_name[name_len] = '\0';
   return true;
}

static color_map* init_color_map_generated(const char *map_name)
{
   char name[MSG_SZ];
   strncpy(name, map_name, MSG_SZ);
   name[MSG_SZ-1] = '\0';
   size_t name_len = strspn(name, "abcdefghijklmnopqrstuvwxyz");
   name[name_len] = '\0';

   color_map *cmap = 0;
     
   if(*map_name=='\0' || strcmp(name, "null")==0)  // A null map
      cmap = new color_map_map();
   
   else if(strcmp(name, "rand")==0) {
      cmap = new color_map_range_rand_hsv();
      if(cmap && !cmap->init(map_name+4)) {
         delete cmap;
         cmap = 0;
      }
   }

   else if(strcmp(name, "randhsv")==0) {
      cmap = new color_map_range_rand_hsv();
      if(cmap && !cmap->init(map_name+7)) {
         delete cmap;
         cmap = 0;
      }
   }

   else if(strcmp(name, "randrgb")==0) {
      cmap = new color_map_range_rand_rgb();
      if(cmap && !cmap->init(map_name+7)) {
         delete cmap;
         cmap = 0;
      }
   }

   else if(strcmp(name, "hsv")==0) {
      cmap = new color_map_range_hsv();
      if(cmap && !cmap->init(map_name+3)) {
         delete cmap;
         cmap = 0;
      }
   }

   else if(strcmp(name, "rgb")==0) {
      cmap = new color_map_range_rgb();
      if(cmap && !cmap->init(map_name+3)) {
         delete cmap;
         cmap = 0;
      }
   }

   else if(strcmp(name, "spread")==0) {
      cmap = new color_map_spread();
      if(cmap && !cmap->init(map_name+6)) {
         delete cmap;
         cmap = 0;
      }
   }

   else if(strcmp(name, "remap")==0)
      cmap = new color_map_remap();

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

   string alt_name;
   FILE *cfile = open_sup_file(name, "/col_maps/", &alt_name);
   if(alt_name!="") {  // an alt name found before a file with the name
      cmap = init_color_map_generated(alt_name.c_str());
      if(!cmap && errmsg)
         snprintf(errmsg, MSG_SZ, "could not open colour map file"
               " \'%s=%s\'", map_name, alt_name.c_str());
   }
   else if(cfile) {
      char errmsg2[MSG_SZ];
      cmap = new color_map_map();
      if(cmap && !cmap->init(map_name, errmsg2)) {
         delete cmap;
         cmap = 0;
         if(errmsg)
            strcpy(errmsg, errmsg2);
      }
   }
   else {
      cmap = init_color_map_generated(map_name);
      if(errmsg)
         snprintf(errmsg, MSG_SZ, "could not open colour map file \'%s\'",
                  map_name);
   }

   return cmap;
}


void color_map::set_wrap(int wrp)
{
   if(wrp<0) wrp=max_index();
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
   
   int ret = color_map_range::init(map_name, errmsg);
   
   return ret;
}

col_val color_map_range_rand::get_col(int idx)
{
   col_val col;
   idx = get_effective_index(idx);
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
   float S, V;
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

      // stage == 2
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
            snprintf(errmsg, MSG_SZ, "colour map: line %d: invalid colour '%s'", line_no, buf);
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
   char buf[line_size];

   
   if(errmsg)
      *errmsg='\0';

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

      if(!cmap->size() && !strchr(line, '=')) {
         rewind(cfile);
         return parse_gimp_file(cfile, cmap, errmsg);
      }
            
      char *val;
      if(!(val = strtok(line, "="))) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "colour map: line %d: missing '='", line_no);
         return false;
      }

      int idx_no;
      if(!read_int(val, &idx_no, buf)) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "colour map: line %d: invalid index number '%s'", line_no, val);
         return false;
      }

      val=strtok(NULL, "\n");

      col_val col;
      col.read(val, errmsg);
      if(!col.is_set()) {
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "colour map: line %d: invalid colour '%s'", line_no, val);
         return false;
      }

      //col.get_val().dump();
      (*cmap)[idx_no] = col;
   }
   
   return true;

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



unsigned int color_map_map::max_index() const
{
   return size() ? cmap.rbegin()->first : 0;
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
      if(!cmap_ok && errmsg)
         strcpy(errmsg, errmsg2);
   }
   if(get_wrap()==-1)
      set_wrap(max_index());
   
   if(cfile)
      fclose(cfile);

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
