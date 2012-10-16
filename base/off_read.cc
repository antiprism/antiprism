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

/* \file off_read.cc
   \brief Read OFF files
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <algorithm>
#include <vector>

#include "utils.h"
#include "off_file.h"
#include "std_polys.h"
#include "polygons.h"

using std::vector;


int read_off_line(FILE *fp, char **line)
{
   int ret = read_line(fp, line);
   if(ret==0) {
      char *first_hash = strchr(*line, '#');
      if(first_hash)
         *first_hash = '\0';
   }
   return ret;
}


bool off_file_read(string file_name, geom_if &geom, char *errmsg)
{  
   if(errmsg)
      *errmsg = '\0';
   
   char errmsg2[MSG_SZ];
   bool geom_ok = false;
   string alt_name;
   FILE *ifile;
   if(file_name == "" || file_name == "-") {
      ifile = stdin;
      file_name = "stdin";
   }
   else
      ifile = open_sup_file(file_name.c_str(), "/models/", &alt_name);
      
   if(alt_name!="") {  // an alt name found before a file with the name
      if(make_resource_geom(geom, alt_name, errmsg2))
         geom_ok = true;
      else {
         if(errmsg)
            snprintf(errmsg, MSG_SZ-50, "could not open input"
                  "file \'%s=%s\': %s", file_name.c_str(), alt_name.c_str(),
                  errmsg2);
      }
   }
   else if(ifile) {       // the file name was found
      geom_ok = off_file_read(ifile, geom, errmsg);
      if(errmsg && *errmsg) {
         string msg("reading \'"+file_name+"\': "+ errmsg);
         snprintf(errmsg, MSG_SZ, "%s", msg.c_str());
      }
   }
   else {                 // try the name as an internal identifier
      if(make_resource_geom(geom, file_name, errmsg2))
         geom_ok = true;
      else {
         if(errmsg)
            snprintf(errmsg, MSG_SZ-50, "could not open input"
                  "file \'%s\': %s", file_name.c_str(), errmsg2);
      }
   }

   if(ifile && ifile != stdin)
      fclose(ifile);

   return geom_ok;
}


bool add_vert(geom_if &geom, vector<char *> vals, char *errmsg)
{
   char errmsg2[MSG_SZ];
   vec3d v;
   for(unsigned int i=0; (i<vals.size()&&i<3); i++) {
      if(!(read_double_noparse(vals[i], &v[i], errmsg2))) {
         sprintf(errmsg, "vertex coords: '%s' %s", vals[i], errmsg2);
         return false;
      }
   }
  
   if(vals.size() < 3) {
      snprintf(errmsg, MSG_SZ, "vertex coords: less than three coordinates");
      return false;
   }
   geom.add_vert(v);

   return true;
} 


bool add_face(geom_if &geom, col_geom *cg, vector<char *> vals,
      char *errmsg, col_geom_v &alt_cols, bool *contains_int_gt_1,
      bool *contains_adj_equal_idx)
{
   char errmsg2[MSG_SZ];
   int face_sz;
   if(!vals.size()) {
      sprintf(errmsg, "face: no face data");
      return false;
   }
   if(!read_int(vals[0], &face_sz, errmsg2)) {
      sprintf(errmsg, "face size: '%s' %s", vals[0], errmsg2);
      return 0;
   }
   if(face_sz<1) {
      sprintf(errmsg, "face size: '%d', must be 1 or more", face_sz);
      return 0;
   }
   *contains_adj_equal_idx = false;
   vector<int> face(face_sz);
   for(unsigned int i=1; (i<vals.size()&&(int)i<=face_sz); i++) {
      if(!(read_int(vals[i], &face[i-1], errmsg2))) {
         sprintf(errmsg, "face index: '%s' %s", vals[i], errmsg2);
         return false;
      }
      int last_vert = geom.verts().size()-1;
      if(face[i-1]<0 || face[i-1]>last_vert) {
         sprintf(errmsg, "face index: '%s' is not in range 0 to %d", vals[i],
               last_vert);
         return false;
      }
      if(i>1 && face[i-1]==face[i-2])
         *contains_adj_equal_idx = true;
   }
   if(face_sz>1 && face[0]==face[face_sz-1])
      *contains_adj_equal_idx = true;

   if((int)vals.size()-1 < face_sz) {
      snprintf(errmsg, MSG_SZ, "face: less than %d values", face_sz);
      return false;
   }
  
   vals.erase(vals.begin(), vals.begin()+face_sz+1);
   int col_type;
   col_val col, alt_col;
   if(!col.from_offvals(vals, errmsg2, &col_type)) {
      snprintf(errmsg, MSG_SZ,"face colour: invalid colour: %s", errmsg2);
      return false;
   }
   
   alt_col = col;
   if(col_type==3||col_type==4) { // read as integers
      if(col[0]>1||col[1]>1||col[2]>1||(col_type==4&&col[3]>1))
         *contains_int_gt_1 = true;
      else   // store alternative colour with integres taken as floats
         alt_col = col_val(col[0]*255, col[1]*255, col[2]*255,
                                         col_type==4?col[3]*255:255);
   }

   int idx;
   if(face_sz==1) {      // vertex element, only need to set colour
      if(cg) {
         cg->set_v_col(face[0], col);
         alt_cols.set_v_col(face[0], alt_col);
      }
   }
   else if(face_sz==2) { // digon edge element
      idx = geom.add_edge(face);
      if(cg) {
         cg->set_e_col(idx, col);
         alt_cols.set_e_col(idx, alt_col);
      }
   }
   else {               // face element
      idx = geom.add_face(face);
      if(cg) {
         cg->set_f_col(idx, col);
         alt_cols.set_f_col(idx, alt_col);
      }
   }
   
   return true;
} 

bool off_file_read(FILE *ifile, geom_if &geom, char *errmsg)
{   
   char errmsg2[MSG_SZ];
   col_geom *cg = dynamic_cast<col_geom *>(&geom);

   int file_line_no = 0; // line number in the file
   
   if(errmsg)
      *errmsg='\0';

   // read OFF type
   int read_ret;
   char *line=0;
   while((read_ret = read_off_line(ifile, &line))==0) {
      file_line_no++;
      if(sscanf(line, " %*s")!=EOF)
         break;
      else
         free(line);
   }

   if(!strstr(line, "OFF")) {
      if(*line == '3') {
         if(errmsg)
            strncpy(errmsg, "assuming file has Qhull OFF output format",MSG_SZ);
      }
      else {
         if(errmsg)
            strncpy(errmsg, "assuming file is list of coordinates", MSG_SZ);
         crds_file_read(ifile, geom, line);
         if(errmsg && !geom)
            strncat(errmsg, ": no coordinates found", MSG_SZ);
         return geom;
      }
   }

   free(line); // finished with file format line
      
   // read counts of coords, polys (and edges)
   while((read_ret = read_off_line(ifile, &line))==0) {
      file_line_no++;
      if(sscanf(line, " %*s")!=EOF)
         break;
      else
         free(line);
   }
 
   int num_pts, num_faces;
   int scan_ret = sscanf(line, " %d %d", &num_pts, &num_faces);
   free(line); // finished with vert and face count line
   
   if(scan_ret < 2) {
      if(errmsg)
         snprintf(errmsg, MSG_SZ, "line %d: didn't find face and vertex counts", file_line_no);
      return false;
   }

   int data_line_no = 2; // non blank lines
  
   // Variables so that if all integer color values
   // are 0 or 1, then they are all converted to decimals
   bool contains_int_gt_1 = false;
   col_geom_v alt_cols;

   // First few line numbers for faces with adjacent verts with equal indexs
   const unsigned int max_adj_equal_idx_lines = 6;
   vector<int> adj_equal_idx_lines;

   // read coords
   while((read_ret = read_off_line(ifile, &line))==0) {
      file_line_no++;
      
      vector<char *> vals;
      int split_ret = split_line(line, vals); 
      if(!split_ret)       // line was blank
            continue;      // skip the line

      data_line_no++;

      if(data_line_no <= 2+num_pts) {   // vertex line
         if(!add_vert(geom, vals, errmsg2)) {
            if(errmsg)
               snprintf(errmsg, MSG_SZ, "line %d: %s", file_line_no, errmsg2);
            geom.clear_all();
            break;
         }
      }
      else if(data_line_no <= 2+num_pts+num_faces) {   // face line
         bool contains_adj_equal_idx;
         if(!add_face(geom, cg, vals, errmsg2, alt_cols,
                  &contains_int_gt_1, &contains_adj_equal_idx)) {
            if(errmsg)
               snprintf(errmsg, MSG_SZ, "line %d: %s", file_line_no, errmsg2);
            geom.clear_all();
            break;
         }
         // only record the first few lines of elements with seq equal indexes
         if(contains_adj_equal_idx &&
               adj_equal_idx_lines.size() < max_adj_equal_idx_lines)
            adj_equal_idx_lines.push_back(file_line_no);

      }
      else {  // extra data at end
         if(errmsg)
            snprintf(errmsg, MSG_SZ, "line %d: data at end of file", file_line_no);
         geom.clear_all();
         break;
      }
      free(line);
   }
   
   free(line);

   if(cg && !contains_int_gt_1) {
      cg->raw_vert_cols() = alt_cols.vert_cols();
      cg->raw_edge_cols() = alt_cols.edge_cols();
      cg->raw_face_cols() = alt_cols.face_cols();
   }

   // create warning message for adjacent equal vertex numbers on faces
   if(errmsg && adj_equal_idx_lines.size()) {
      string msg("line");
      msg += ((adj_equal_idx_lines.size()>1)?"s ":" ");
      for(unsigned int i=0;
            i<adj_equal_idx_lines.size() && i<max_adj_equal_idx_lines-1; i++)
         msg += itostr(adj_equal_idx_lines[i])+", ";

      if(adj_equal_idx_lines.size() == max_adj_equal_idx_lines)
         msg += "...";              // the unmentioned last line and any others
      else
         msg.resize(msg.size()-2);  // the list was complete

      msg += ": face element has adjacent vertices with the same index number";
      if(*errmsg)  // already a message
         strncat(errmsg, ", and, ", MSG_SZ);
      strncat(errmsg, msg.c_str(), MSG_SZ);
   }

   if(errmsg && !bool(geom) && !*errmsg)  // no previous error message
      strncpy(errmsg, "no coordinate data", MSG_SZ);

   return bool(geom);
} 

