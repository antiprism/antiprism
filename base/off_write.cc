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

/* \file off_file.cc
   \brief Write OFF files
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <algorithm>
#include <vector>

#include "utils.h"
#include "off_file.h"

using std::vector;

FILE *file_open_w(string file_name, char *errmsg)
{
   FILE *ofile = stdout;  // write to stdout by default
   if(file_name != "") {
      ofile = fopen(file_name.c_str(), "w");
      if(!ofile && errmsg)
         snprintf(errmsg, MSG_SZ, "could not output file \'%s\'",
               file_name.c_str());
   }
   return ofile;
}

void file_close_w(FILE *ofile)
{
   if(ofile!=stdout)
      fclose(ofile);
}


void crds_write(FILE *ofile, const geom_if &geom, const char *sep, int sig_dgts)
{
   char line[MSG_SZ];
   for(unsigned int i=0; i<geom.verts().size(); i++)
      fprintf(ofile, "%s\n", vtostr(line, geom.verts(i), sep, sig_dgts));
}   


bool crds_write(string file_name, const geom_if &geom, char *errmsg,
      const char *sep, int sig_dgts)
{
   FILE *ofile = file_open_w(file_name, errmsg);
   if(!ofile)
      return false;

   crds_write(ofile, geom, sep, sig_dgts);
   file_close_w(ofile);
   return true;
}


bool off_file_write(string file_name, const geom_if &geom, char *errmsg,
      int sig_dgts)
{
   vector<const geom_if *> vg;
   vg.push_back(&geom);
   return off_file_write(file_name, vg, errmsg, sig_dgts);
}


bool off_file_write(string file_name, const vector<const geom_if *> &geoms,
      char *errmsg, int sig_dgts)
{
   FILE *ofile = file_open_w(file_name, errmsg);
   if(!ofile)
      return false;

   off_file_write(ofile, geoms, sig_dgts);
   file_close_w(ofile);
   return true;
}


char *off_col(char *str, col_val col)
{
    if(col.is_idx())
       snprintf(str, MSG_SZ-1, " %d", col.get_idx());
    else if(col.is_val()) {
       if(col.get_trans())
          vtostr(str, col.get_vec4d(), " ", -5);
       else 
          vtostr(str, col.get_vec3d(), " ", -5);
    }
    else
       *str = '\0';
    
    return str;
}




void off_polys_write(FILE *ofile, const geom_if &geom, int offset)
{
   char col_str[MSG_SZ];
   const col_geom *cg = dynamic_cast<const col_geom *>(&geom);
   for(unsigned int i=0; i<geom.faces().size(); i++) {
      fprintf(ofile, "%lu", (unsigned long)geom.faces(i).size());
      for(unsigned int j=0; j<geom.faces(i).size(); j++) {
         fprintf(ofile, " %d", geom.faces(i, j)+offset);
      }
      if(cg)
         fprintf(ofile, " %s", off_col(col_str, cg->get_f_col(i)));
      fprintf(ofile, "\n");
   }
   
   for(unsigned int i=0; i<geom.edges().size(); i++) {
      fprintf(ofile, "2 %d %d", geom.edges(i, 0)+offset,
            geom.edges(i, 1)+offset);
      if(cg)
         fprintf(ofile, " %s", off_col(col_str, cg->get_e_col(i)));
      fprintf(ofile, "\n");
   }
   // print coloured vertex elements
   if(cg) {
      map<int, col_val>::const_iterator mi;
      for(mi=cg->vert_cols().begin(); mi!=cg->vert_cols().end(); mi++){
         col_val col = mi->second;
         fprintf(ofile, "1 %d %s\n", mi->first, off_col(col_str, mi->second));
      }
   }
}   

void off_file_write(FILE *ofile, const vector<const geom_if *> &geoms, int sig_dgts)
{
   int vert_cnt=0, face_cnt=0, edge_cnt=0;
   for(unsigned int i=0; i<geoms.size(); i++) {
      const col_geom *cg = dynamic_cast<const col_geom *>(geoms[i]);
      int num_v_col_elems = cg ? cg->vert_cols().size(): 0;
      vert_cnt += geoms[i]->verts().size();
      edge_cnt += geoms[i]->edges().size();
      face_cnt += geoms[i]->faces().size() + num_v_col_elems + edge_cnt;
   }
   
   fprintf(ofile, "OFF\n%d %d 0\n", vert_cnt, face_cnt);
  
   for(unsigned int i=0; i<geoms.size(); i++)
      crds_write(ofile, *geoms[i], " ", sig_dgts);

   int last_offset=0;
   vert_cnt=0;
   for(unsigned int i=0; i<geoms.size(); i++) {
      off_polys_write(ofile, *geoms[i], geoms[i]->verts().size() ? vert_cnt : last_offset);
      last_offset = vert_cnt;
      vert_cnt += geoms[i]->verts().size();
   }
}

void off_file_write(FILE *ofile, const geom_if &geom, int sig_dgts)
{
   vector<const geom_if *> vg;
   vg.push_back(&geom);
   off_file_write(ofile, vg, sig_dgts);
}

