/*
   Copyright (c) 2006-2016, Adrian Rossiter

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

/*
   Name: string_art.cc
   Description: string art
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using std::map;
using std::string;
using std::vector;

using namespace anti;

class art_elem {
protected:
  int num_pins;
  Trans3d trans_m;

public:
  art_elem(int num = 1, Trans3d trans = Trans3d::unit())
      : num_pins(num), trans_m(trans)
  {
  }
  virtual ~art_elem() = default;
  void set_trans(Trans3d trans) { trans_m = trans; }
  int get_num_pins() { return num_pins; }
  Vec3d get_pin(int pin) { return trans_m * get_elem_pin(pin); }
  virtual Vec3d get_elem_pin(int pin) = 0;
};

class circle_elem : public art_elem {
public:
  circle_elem(int num = 1, Trans3d trans = Trans3d::unit())
      : art_elem(num, trans)
  {
  }

  Vec3d get_elem_pin(int pin)
  {
    double a = 2 * M_PI * (pin % num_pins) / num_pins;
    return Vec3d(cos(a), sin(a), 0);
  }
};

class line_elem : public art_elem {
private:
  Vec3d p1;
  Vec3d p2;

public:
  line_elem(int num = 1, Trans3d trans = Trans3d::unit()) : art_elem(num, trans)
  {
    set_coords(Vec3d(-1, 0, 0), Vec3d(1, 0, 0));
  }

  void set_coords(Vec3d v1, Vec3d v2)
  {
    p1 = v1;
    p2 = v2;
  }
  Vec3d get_elem_pin(int pin)
  {
    double d = (pin % num_pins) / (num_pins - 1.0);
    return p1 + (p2 - p1) * d;
  }
};

class string_art : public Geometry {
private:
  vector<art_elem *> elems;

public:
  string_art() {}
  ~string_art()
  {
    for (auto &elem : elems)
      delete elem;
  }
  void add_elem(art_elem *elem) { elems.push_back(elem); }
  void make_figure();
};

void string_art::make_figure()
{
  vector<Vec3d> &verts = raw_verts();
  vector<vector<int>> &faces = raw_faces();
  verts.clear();
  faces.clear();
  clear(EDGES);
  if (!elems.size())
    return;

  for (int i = 0; i < elems[0]->get_num_pins(); i++) {
    for (unsigned int e = 0; e < elems.size(); e++) {
      verts.push_back(elems[e]->get_pin(i));
      if (e) {
        faces.push_back(vector<int>());
        faces.back().push_back(verts.size() - 2);
        faces.back().push_back(verts.size() - 1);
      }
    }
    /*
    Vec3d &v0 = verts[verts.size()-2];
    Vec3d &v1 = verts[verts.size()-1];
    Vec3d v = v1 - v0;
    Vec3d v_mid = v0 + v/2.0;
    double len = (verts[1]-verts[0]).mag();
    double new_len = len/2 * (elems[0]->get_num_pins() - i) /
    elems[0]->get_num_pins();
    v = v.with_len(new_len);
    v0 = v_mid - v;
    v1 = v_mid + v;
    */
  }
}

class string_opts : public ProgramOpts {
public:
  string_art figure;
  string ifile;
  string ofile;

  string_opts() : ProgramOpts("string_art") {}
  void process_command_line(int argc, char **argv);
  void usage();
};

void string_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options]

Create a string art in OFF fomat. Each pin arrangement shape is
modified by any immediately following transformation options.

Options
%s
  -c <num>  a unit circle arrangement of num pins
  -l <num>  a 2 unit line arrangement of num pins
  -p <crds> line end points coordinates, six numbers separated by commas
  -T <tran> translate, three numbers separated by commas which are
            used as the x, y and z displacements
  -R <rot>  rotate about an axis, three, four or six numbers separated by
            commas. If three numbers these are angles (degrees) to rotate
            about the x, y and z axes. If four numbers, the first three
            are a direction vector for the axis, the last number is the
            angle (degrees) to rotate. If six numbers, these are two
            vectors and rotate to carry the first to the second
  -M <norm> reflect in a plane, three numbers separated by commas which
            give a vector normal to the plane of reflection.
  -S <scal> scale, one, three or four numbers separated by commas. If one
            number then scale by this factor in all directions. If three
            numbers these are the factors to scale along the x, y and
            z axes. If four numbers, the first three are a direction
            vector for the scaling, the last number is the factor to scale
  -I        invert
  -o <file> write output to file (default: write to standard output)

)",
          prog_name(), help_ver_text);
}

void string_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;
  vector<double> nums;
  vector<Vec3d> coords;
  art_elem *elem = nullptr;
  int num;
  Trans3d trans_m;
  Trans3d trans_m2;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hc:l:p:T:R:M:S:Io:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'c':
    case 'l':
      if (elem) {
        line_elem *le;
        if ((le = dynamic_cast<line_elem *>(elem)) && coords.size() == 2)
          le->set_coords(coords[0], coords[1]);
        elem->set_trans(trans_m);
        figure.add_elem(elem);
      }
      trans_m = Trans3d::unit();
      coords.clear();
      print_status_or_exit(read_int(optarg, &num), c);
      if (c == 'c')
        elem = new circle_elem(num);
      else
        elem = new line_elem(num);
      break;

    case 'p':
      print_status_or_exit(read_double_list(optarg, nums), c);
      if (nums.size() == 6) {
        coords.push_back(Vec3d(nums[0], nums[1], nums[2]));
        coords.push_back(Vec3d(nums[3], nums[4], nums[5]));
      }
      else
        error(msg_str("must give six numbers (%lu were given)",
                      (unsigned long)nums.size()),
              c);
      break;

    case 'R':
      print_status_or_exit(read_double_list(optarg, nums), c);
      if (nums.size() == 3)
        trans_m2 = Trans3d::rotate(deg2rad(nums[0]), deg2rad(nums[1]),
                                   deg2rad(nums[2]));
      else if (nums.size() == 4)
        trans_m2 =
            Trans3d::rotate(Vec3d(nums[0], nums[1], nums[2]), deg2rad(nums[3]));
      else if (nums.size() == 6)
        trans_m2 = Trans3d::rotate(Vec3d(nums[0], nums[1], nums[2]),
                                   Vec3d(nums[3], nums[4], nums[5]));
      else
        error(msg_str("must give three, four or six numbers "
                      "(%lu were given)",
                      (unsigned long)nums.size()),
              c);
      trans_m = trans_m2 * trans_m;
      break;

    case 'S':
      print_status_or_exit(read_double_list(optarg, nums), c);
      if (nums.size() == 1)
        trans_m2 = Trans3d::scale(nums[0]);
      else if (nums.size() == 3)
        trans_m2 = Trans3d::scale(nums[0], nums[1], nums[2]);
      else if (nums.size() == 4)
        trans_m2 = Trans3d::scale(Vec3d(nums[0], nums[1], nums[2]), nums[3]);
      else
        error(msg_str("must give 1, 3 or 4 numbers (%lu were given)",
                      (unsigned long)nums.size()),
              c);
      trans_m = trans_m2 * trans_m;
      break;

    case 'T':
      print_status_or_exit(read_double_list(optarg, nums), c);
      if (nums.size() != 3)
        error(msg_str("must give exactly three numbers (%lu were given)",
                      (unsigned long)nums.size()),
              c);
      trans_m2 = Trans3d::translate(Vec3d(nums[0], nums[1], nums[2]));
      trans_m = trans_m2 * trans_m;
      break;

    case 'M':
      print_status_or_exit(read_double_list(optarg, nums), c);
      if (nums.size() != 3)
        error(msg_str("must give exactly three numbers (%lu were given)",
                      (unsigned long)nums.size()),
              c);

      trans_m2 = Trans3d::reflection(Vec3d(nums[0], nums[1], nums[2]));
      trans_m = trans_m2 * trans_m;
      break;

    case 'I':
      trans_m = Trans3d::inversion() * trans_m;
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (elem) {
    line_elem *le;
    if ((le = dynamic_cast<line_elem *>(elem)) && coords.size() == 2)
      le->set_coords(coords[0], coords[1]);
    elem->set_trans(trans_m);
    figure.add_elem(elem);
  }
  else
    error("no art elements given with -c, -l");

  if (argc - optind > 0)
    error("too many arguments");
}

int main(int argc, char *argv[])
{
  string_opts opts;
  opts.process_command_line(argc, argv);

  opts.figure.make_figure();

  opts.write_or_error(opts.figure, opts.ofile);

  return 0;
}
