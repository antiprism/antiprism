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

/*
   Name: repel.cc
   Description: points repelling (on a sphere)
   Project: Antiprism - http://www.antiprism.com
*/

#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <time.h>
#include <vector>

#include "../base/antiprism.h"

using std::string;
using std::vector;

using namespace anti;

class rep_opts : public ProgramOpts {
public:
  IterationControl it_ctrl;
  int num_pts;
  int rep_form;
  double shorten_by;

  string ifile;
  string ofile;

  rep_opts() : ProgramOpts("repel"), num_pts(-1), rep_form(2), shorten_by(-1)
  {
    it_ctrl.set_max_iters_unlimited();
    it_ctrl.set_status_checks("1000,1");
    it_ctrl.set_sig_digits(12);
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

void rep_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] [input_file]

An equilibrium position is found for a set of points which repel each
other. The initial coordinates are read from input_file if given (or
from standard input), otherwise use -N to generate a random set.

Options
%s
  -N <num>  initialise with a number of randomly placed points
  -n <itrs> maximum number of iterations, -1 for unlimited (default: %d)
  -s <perc> percentage to shorten the travel distance (default: adaptive)
  -r <rep>  repelling formula
              1 - inverse of distance
              2 - inverse square of distance (default)
              3 - inverse cube of distance
              4 - inverse square root of distance
  -l <lim>  minimum change of distance/width_of_model to terminate, as 
               negative exponent (default: %d giving %.0e)
  -z <nums> number of iterations between status reports (implies termination
            check) (0 for final report only, -1 for no report), optionally
            followed by a comma and the number of iterations between
            termination checks (0 for report checks only) (default: %d,%d)
  -o <file> write output to file (default: write to standard output)

)",
          prog_name(), help_ver_text, it_ctrl.get_max_iters(),
          it_ctrl.get_sig_digits(), it_ctrl.get_test_val(),
          it_ctrl.get_status_check_and_report_iters(),
          it_ctrl.get_status_check_only_iters());
}
// clang-format on

void rep_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  int num;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hn:z:N:s:l:r:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'n':
      print_status_or_exit(read_int(optarg, &num), c);
      print_status_or_exit(it_ctrl.set_max_iters(num), c);
      break;

    case 'z':
      print_status_or_exit(it_ctrl.set_status_checks(optarg), c);
      break;

    case 'N':
      print_status_or_exit(read_int(optarg, &num_pts), c);
      if (num_pts < 2)
        error("number of points must be 2 or more", c);
      break;

    case 's':
      print_status_or_exit(read_double(optarg, &shorten_by), c);
      if (shorten_by <= 0 || shorten_by >= 100)
        warning("not inside range 0 to 100", c);
      break;

    case 'l':
      print_status_or_exit(read_int(optarg, &num), c);
      print_status_or_exit(it_ctrl.set_sig_digits(num), c);
      break;

    case 'r':
      print_status_or_exit(read_int(optarg, &rep_form), c);
      if (rep_form < 1 || rep_form > 4)
        error("formula is given by its number, 1 - 4", c);
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (argc - optind > 1)
    error("too many arguments");

  if (argc - optind == 1) {
    if (num_pts > 0)
      error(
          "cannot specify a number of points and also to read them from a file",
          'N');
    else
      ifile = argv[optind];
  }
}

typedef Vec3d (*REPEL_FN)(Vec3d, Vec3d);

Vec3d rep_inv_dist1(Vec3d v1, Vec3d v2)
{
  double len = pow((v2 - v1).len2(), -0.5);
  return (v2 - v1).with_len(len);
}

Vec3d rep_inv_dist2(Vec3d v1, Vec3d v2)
{
  double len = 1 / (v2 - v1).len2();
  return (v2 - v1).with_len(len);
}

Vec3d rep_inv_dist3(Vec3d v1, Vec3d v2)
{
  double len = pow((v2 - v1).len2(), -1.5);
  return (v2 - v1).with_len(len);
}

Vec3d rep_inv_dist05(Vec3d v1, Vec3d v2)
{
  double len = pow((v2 - v1).len2(), -0.25);
  return (v2 - v1).with_len(len);
}

void random_placement(Geometry &geom, int n)
{
  geom.clear_all();
  Random rnd;
  rnd.time_seed();
  for (int i = 0; i < n; i++)
    geom.add_vert(Vec3d::random(rnd).unit());
}

void repel(Geometry &geom, IterationControl it_ctrl, REPEL_FN rep_fn,
           double shorten_factor)
{
  const int v_sz = geom.verts().size();
  vector<int> wts(v_sz);
  for (int i = 0; i < v_sz; i++) {
    Color col = geom.colors(VERTS).get(i);
    wts[i] = col.is_index() ? col.get_index() : 1;
  }
  vector<Vec3d> offsets(v_sz);
  double dist2, max_dist2 = 0;
  double last_av_max_dist2 = 0, max_dist2_sum = 0;
  bool adaptive = false;
  int chng_cnt = 0;
  int converge = 0;
  const int anum = 50;
  if (shorten_factor < 0) {
    adaptive = true;
    shorten_factor = 0.001;
  }

  double test_val = it_ctrl.get_test_val();

  for (it_ctrl.start_iter(); !it_ctrl.is_done(); it_ctrl.next_iter()) {
    std::fill(offsets.begin(), offsets.end(), Vec3d::zero);
    max_dist2 = 0;

    for (int i = 0; i < v_sz - 1; i++) {
      for (int j = i + 1; j < v_sz; j++) {
        Vec3d offset = rep_fn(geom.verts(i), geom.verts(j)) * (wts[i] * wts[j]);
        offsets[i] -= offset / wts[i];
        offsets[j] += offset / wts[j];
      }
    }

    for (int i = 0; i < v_sz; i++) {
      Vec3d new_pos = (geom.verts(i) + offsets[i] * shorten_factor).unit();
      dist2 = (new_pos - geom.verts(i)).len2();
      if (dist2 > max_dist2)
        max_dist2 = dist2;
      if (dist2 > 0.001)
        new_pos = (new_pos + geom.verts(i)).unit();
      geom.verts(i) = new_pos;
    }

    if (adaptive) {
      max_dist2_sum += max_dist2;
      if (max_dist2 < last_av_max_dist2)
        converge += 1;
      if (!(it_ctrl.get_current_iter() % anum)) {
        if (converge > anum / 1.5) {
          chng_cnt = chng_cnt > 0 ? chng_cnt + 1 : 1;
          shorten_factor *= 1 + 0.005 * chng_cnt;
          if (shorten_factor > 0.5)
            shorten_factor = 0.5;
        }
        else if (converge < anum / 2) {
          chng_cnt = chng_cnt < 0 ? chng_cnt - 1 : -2;
          shorten_factor *= 0.995 + 0.005 * chng_cnt;
          if (shorten_factor < 0)
            shorten_factor = 0.00001;
        }
        else
          chng_cnt = 0;

        last_av_max_dist2 = max_dist2_sum / anum;
        //             fprintf(stderr, "%2d ", converge / 5);
        // fprintf(stderr, "\n%d %g %g (%g)\n", converge, shorten_factor,
        // last_av_max_dist2, max_dist2);
        converge = 0;
        max_dist2_sum = 0;
      }
    }

    if (it_ctrl.is_status_check_iter()) {
      string finish_reason;

      if (sqrt(max_dist2) < test_val) {
        it_ctrl.set_finished();
        finish_reason = "solved, test value achieved";
      }
      else if (it_ctrl.is_last_iter()) {
        // reached last iteration without solving
        it_ctrl.set_finished();
        finish_reason = "not solved, test value not achieved";
      }

      if (it_ctrl.is_status_report_iter()) {
        double offset_sum = 0.0;
        for (auto &offset : offsets)
          offset_sum += offset.len();

        if (it_ctrl.is_finished())
          it_ctrl.print("Final iteration (%s):\n", finish_reason.c_str());
        it_ctrl.print("%-12u  max_diff:%-16.10g  s:%-11.6g  F-sum:%-.16g\n",
                      it_ctrl.get_current_iter(), sqrt(max_dist2),
                      shorten_factor, offset_sum);
      }
    }
  }
}

int main(int argc, char *argv[])
{
  rep_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;
  if (opts.num_pts > 0)
    random_placement(geom, opts.num_pts);
  else
    opts.read_or_error(geom, opts.ifile);

  REPEL_FN fn[] = {rep_inv_dist1, rep_inv_dist2, rep_inv_dist3, rep_inv_dist05};
  repel(geom, opts.it_ctrl, fn[opts.rep_form - 1], opts.shorten_by / 100);

  opts.write_or_error(geom, opts.ofile);

  return 0;
}
