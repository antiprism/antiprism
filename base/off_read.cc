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

/* \file off_read.cc
   \brief Read OFF files
*/

#include "polygon.h"
#include "private_off_file.h"
#include "private_std_polys.h"
#include "utils.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using std::string;
using std::to_string;
using std::vector;

int read_off_line(FILE *fp, char **line)
{
  int ret = read_line(fp, line);
  if (ret == 0) {
    char *first_hash = strchr(*line, '#');
    if (first_hash)
      *first_hash = '\0';
  }
  return ret;
}

Status off_file_read(string file_name, Geometry &geom)
{
  Status stat;
  string alt_name;
  FILE *ifile;
  if (file_name == "" || file_name == "-") {
    ifile = stdin;
    file_name = "stdin";
  }
  else
    ifile = open_sup_file(file_name.c_str(), "/models/", &alt_name);

  if (alt_name != "") { // an alt name found before a file with the name
    if (!(stat = make_resource_geom(geom, alt_name)))
      stat.set_error("could not open input file from alternative name '" +
                     file_name + "' = '" + alt_name + "': " + stat.msg());
  }
  else if (ifile) { // the file name was found
    stat = off_file_read(ifile, geom);
    if (!stat.msg().empty()) {
      string msg("reading '" + file_name + "': " + stat.msg());
    }
  }
  else { // try the name as an internal identifier
    if (!(stat = make_resource_geom(geom, file_name)))
      stat.set_error("could not open input file '" + file_name +
                     "': " + stat.msg());
  }

  if (ifile && ifile != stdin)
    fclose(ifile);

  return stat;
}

Status add_vert(Geometry &geom, const vector<char *> &vals)
{
  Status stat;
  Vec3d v;
  for (unsigned int i = 0; (i < vals.size() && i < 3); i++) {
    if (!(stat = read_double_noparse(vals[i], &v[i])))
      return Status::error(
          msg_str("vertex coords: '%s' %s", vals[i], stat.c_msg()));
  }

  if (vals.size() < 3)
    return Status::error("vertex coords: less than three coordinates");

  geom.add_vert(v);

  return Status::ok();
}

Status add_face(Geometry &geom, vector<char *> vals, Geometry &alt_cols,
                bool *contains_int_gt_1, bool *contains_adj_equal_idx)
{
  Status stat;
  int face_sz;
  if (!vals.size())
    return Status::error("face: no face data");

  if (!(stat = read_int(vals[0], &face_sz)))
    return Status::error(msg_str("face size: '%s' %s", vals[0], stat.c_msg()));

  if (face_sz < 1)
    return Status::error(
        msg_str("face size: '%d', must be 1 or more", face_sz));

  *contains_adj_equal_idx = false;
  vector<int> face(face_sz);
  for (unsigned int i = 1; (i < vals.size() && (int)i <= face_sz); i++) {
    if (!(stat = read_int(vals[i], &face[i - 1])))
      return Status::error(
          msg_str("face index: '%s' %s", vals[i], stat.c_msg()));

    int last_vert = geom.verts().size() - 1;
    if (face[i - 1] < 0 || face[i - 1] > last_vert)
      return Status::error(msg_str("face index: '%s' is not in range 0 to %d",
                                   vals[i], last_vert));

    if (i > 1 && face[i - 1] == face[i - 2])
      *contains_adj_equal_idx = true;
  }
  if (face_sz > 1 && face[0] == face[face_sz - 1])
    *contains_adj_equal_idx = true;

  if ((int)vals.size() - 1 < face_sz)
    return Status::error(msg_str("face: less than %d values", face_sz));

  vals.erase(vals.begin(), vals.begin() + face_sz + 1);
  int col_type;
  Color col, alt_col;
  if (!(stat = col.from_offvals(vals, &col_type)))
    return Status::error(
        msg_str("face colour: invalid colour: %s", stat.c_msg()));

  alt_col = col;
  if (col_type == 3 || col_type == 4) { // read as integers
    if (col[0] > 1 || col[1] > 1 || col[2] > 1 || (col_type == 4 && col[3] > 1))
      *contains_int_gt_1 = true;
    else // store alternative colour with integres taken as floats
      alt_col = Color(col[0] * 255, col[1] * 255, col[2] * 255,
                      col_type == 4 ? col[3] * 255 : 255);
  }

  int idx;
  if (face_sz == 1) { // vertex element, only need to set colour
    geom.colors(VERTS).set(face[0], col);
    alt_cols.colors(VERTS).set(face[0], alt_col);
  }
  else if (face_sz == 2) { // digon edge element
    idx = geom.add_edge(face);
    geom.colors(EDGES).set(idx, col);
    alt_cols.colors(EDGES).set(idx, alt_col);
  }
  else { // face element
    idx = geom.add_face(face);
    geom.colors(FACES).set(idx, col);
    alt_cols.colors(FACES).set(idx, alt_col);
  }

  return Status::ok();
}

Status off_file_read(FILE *ifile, Geometry &geom)
{
  int file_line_no = 0; // line number in the file

  // read OFF type
  char *line = nullptr;
  while (read_off_line(ifile, &line) == 0) {
    file_line_no++;
    if (sscanf(line, " %*s") != EOF)
      break;
    else
      free(line);
  }

  string message;
  if (!strstr(line, "OFF")) {
    if (*line == '3')
      message = "assuming file has Qhull OFF output format";
    else {
      message = "assuming file is list of coordinates";
      crds_file_read(ifile, geom, line);
      return geom.is_set() ? Status::warning(message)
                           : Status::error(message + ": no coordinates found");
    }
  }

  free(line); // finished with file format line

  // read counts of coords, polys (and edges)
  while (read_off_line(ifile, &line) == 0) {
    file_line_no++;
    if (sscanf(line, " %*s") != EOF)
      break;
    else
      free(line);
  }

  int num_pts, num_faces;
  int scan_ret = sscanf(line, " %d %d", &num_pts, &num_faces);
  free(line); // finished with vert and face count line

  if (scan_ret < 2)
    return Status::error(
        msg_str("line %d: didn't find face and vertex counts", file_line_no));

  if (num_pts < 0 || num_faces < 0)
    return Status::error(
        msg_str("line %d: element counts: %s count is negative", file_line_no,
                (num_pts < 0) ? "vertex" : "face"));

  if (num_pts == 0 && num_faces != 0)
    return Status::error(msg_str("line %d: element counts: cannot have a "
                                 "positive face count if vertex count is zero ",
                                 file_line_no));

  int data_line_no = 2; // non blank lines

  // Variables so that if all integer color values
  // are 0 or 1, then they are all converted to decimals
  bool contains_int_gt_1 = false;
  Geometry alt_cols;

  // First few line numbers for faces with adjacent verts with equal indexs
  const unsigned int max_adj_equal_idx_lines = 6;
  vector<int> adj_equal_idx_lines;

  // read coords
  while (read_off_line(ifile, &line) == 0) {
    file_line_no++;

    Split vals(line);
    if (!vals.size()) // line was blank
      continue;       // skip the line

    data_line_no++;

    Status stat;
    if (data_line_no <= 2 + num_pts) { // vertex line
      if (!(stat = add_vert(geom, vals.get_parts()))) {
        message = msg_str("line %d: ", file_line_no) + stat.msg();
        geom.clear_all();
        break;
      }
    }
    else if (data_line_no <= 2 + num_pts + num_faces) { // face line
      bool contains_adj_equal_idx = false;
      if (!(stat = add_face(geom, vals.get_parts(), alt_cols,
                            &contains_int_gt_1, &contains_adj_equal_idx))) {
        message = msg_str("line %d: ", file_line_no) + stat.msg();
        geom.clear_all();
        break;
      }
      // only record the first few lines of elements with seq equal indexes
      if (contains_adj_equal_idx &&
          adj_equal_idx_lines.size() < max_adj_equal_idx_lines)
        adj_equal_idx_lines.push_back(file_line_no);
    }
    else { // extra data at end
      message = msg_str("line %d: data at end of file", file_line_no);
      geom.clear_all();
      break;
    }
    free(line);
  }

  free(line);

  if (!contains_int_gt_1)
    geom.get_cols() = alt_cols.get_cols();

  // create warning message for adjacent equal vertex numbers on faces
  if (adj_equal_idx_lines.size()) {
    string msg("line");
    msg += ((adj_equal_idx_lines.size() > 1) ? "s " : " ");
    for (unsigned int i = 0;
         i < adj_equal_idx_lines.size() && i < max_adj_equal_idx_lines - 1; i++)
      msg += std::to_string(adj_equal_idx_lines[i]) + ", ";

    if (adj_equal_idx_lines.size() == max_adj_equal_idx_lines)
      msg += "..."; // the unmentioned last line and any others
    else
      msg.resize(msg.size() - 2); // the list was complete

    msg += ": face element has adjacent vertices with the same index number";
    if (!message.empty()) // already a message
      message += ", and, ";
    message += msg;
  }

  if (!geom.is_set() && message.empty()) // no previous error message
    message = "no vertices (empty geometry)";

  if (geom.is_set())
    return message.empty() ? Status::ok() : Status::warning(message);
  else
    return Status::error(message);
}
