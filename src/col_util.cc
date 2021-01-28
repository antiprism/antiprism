/*
   Copyright (c) 2010-2021, Roger Kaufman

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
   Name: col_util.cc
   Description: Creates models of colors. Color wheel, RGB, HSV, HSL
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

using std::max;
using std::min;
using std::string;
using std::vector;

using namespace anti;

class col_util_opts : public ProgramOpts {
public:
  string ifile;
  string ofile;
  string gfile; // optional input off file to collect colors for graphing

  int display_type = 1;       // plot=1 wheel=2 grid=3 map=4 (default: 1)
  int grid_width = 0;         // width of grid (default : automatic)
  int collect_indexes = true; // collect index number from input
  int color_system_mode = 2;  // hsv=1 hsl=2 rgb=3 (default: 2)
  int map_type = 0;           // rgb=1 antiprism=2 decimal=3 (default: 1)
  int show_container = 2;     // 0 - suppress 1 - no facets  2 - full
  int show_reference = false; // show reference colors
  char sort_colors = '\0';    // s - simple,  u - unique
  int upright_view = 0;       // 1 - upright, 2 - upright cube (default : none)
  int container = 0;          // cylinder=1 conic=2 hexagonal=3 cube=4
  int hsl_height = 0;         // 1 - double 2 - angle 90 degrees dihedral
  bool plot_centroid = false; // plot centroid color
  double sat_threshold = 1.0; // saturation threshold to average saturation
  double value_advance = 0;   // value advance in degrees
  int alpha_mode = 3;         // average=1 minimum=2 maximum=3 (default: 3)
  bool cmy_mode = false;      // complementary colors
  bool ryb_mode = false;      // Red-Yellow-Blue color wheel
  bool seven_mode = false;    // use 7 divisions for color models
  double brightness_adj = 0;  // brightness adjustment. -1 - black, 1 - white
  int map_maximum = 0;        // maximum map entries to read for infinite maps
  bool verbose = false;       // output information to screen

  vector<double> sat_powers;   // saturation curve (4 values)
  vector<double> value_powers; // value curve (4 value)

  ColorMapMulti map;

  col_util_opts() : ProgramOpts("col_util") {}

  void process_command_line(int argc, char **argv);
  void usage();
};

void col_util_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options]

various plots of maps, colors and blendings

Options
%s
  -o <file> write output to file (default: write to standard output)
  -f <type> for output type 4, map is output instead of OFF file
               map type: rgb=1, antiprism=2, decimal=3 (default: 1)

Scene Options
  -M <mode> color system mode. HSV=1, HSL=2, RGB=3 (default: 2)
  -d <int>  output type. plot=1, wheel=2, grid=3, map=4 (default: 1)
  -y        RYB mode. Blend colors as in Red-Yellow-Blue color wheel
  -c        CMY mode. Complementary colors.  RGB->(RYB/GMO)->CMY->blend
  -b <val>  brightness adjustment. valid values -1.0 to +1.0 (default: 0)
               negative for darker, positive for lighter, 0 for no change
            (-y also affects location of color on plot. -c and -b affect
             all colors in all output types)

Output for grid (-d 3) Options
  -w <int>  width of grid. positive integer. (default: automatic)
              keyword: b - display as color bars
  -I        exclude map indexes on grid

Output for plot (-d 1) Options 
  -r <int>  HSV/HSL container (default: 4)
              cylinder=1, conic=2, hexagonal=3, cube=4
  -S        HSV/HSL distribute colors heptagonally, 7 ways, as in the rainbow
  -k <int>  container visibility (default: 2)
               0 - suppress, 1 - no facets (-r 3),  2 - full
  -z <int>  1 - show model upright, 2 - upright, don't rotate cube (RGB -M 3)
  -q <int>  height of cone or pyramid (-M 2, -r 2,3) (default: no change)
               1 - double, 2 - make perimeter angle 90 degrees dihedral
  -R        show reference primary and secondary colors on containers

Color Blending Options (-d 1, -M 1,2, -p)
  -p        plot color centroid (or multiple centroids -s, -v)
  -V        verbose output. show color blend RGBA components 
  -s <sat>  HSV/HSL saturation curve. Greater than 0 (default: 1)
               1.0 - no curve. lower than 1.0 makes blends more pastel
               4 numbers can be entered separated by commas
  -t <val>  HSV/HSL threshold to use average saturation (default: 1)
               between 0.0 (all averaging) and 1.0 (no averaging)
  -v <val>  HSV/HSL value curve. Positive integer (default: 0)
               simulates subtractive coloring for blending 3 or more colors
               RGB: Red+Green+Blue = White   Cyan+Magenta+Yellow = Black
               RYB: Red+Yellow+Blue = Black  Green+Magenta+Orange = White
               1.0 - no curve. lower than 1.0 number makes blends lighter
               0.0 - use average value instead
               4 numbers can be entered separated by commas
  -u <val>  HSV/HSL value advance. Rotates meaning of blend to white and black
               valid values 0.0 to 120.0 degrees (default: 0)
  -a <int>  alpha for blend. average=1, minimum=2, maximum=3 (default: 3)

Coloring Options (run 'off_util -H color' for help on color formats)
  -m <map>  get colors from a color map, or multiple maps separated by commas
  -O <file> get colors from an OFF file
               note: -m and -O may be used together
  -U <opt>  sort colors. s - simple,  u - unique  (by hue, saturation, value)
  -Z <int>  maximum entries to read from open ended maps (default: 256)

)",
          prog_name(), help_ver_text);
}

void col_util_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;

  vector<double> double_parms;
  string arg_id;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv,
                     ":hd:w:Im:k:Rq:z:f:r:b:s:t:u:v:pa:cySl:M:O:U:Z:Vo:")) !=
         -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'd':
      print_status_or_exit(get_arg_id(optarg, &arg_id,
                                      "plot=1|wheel=2|grid=3|map=4",
                                      argmatch_add_id_maps),
                           c);
      display_type = atoi(arg_id.c_str());
      break;

    case 'w':
      if (strlen(optarg) == 1 && strchr("b", *optarg))
        grid_width = -1; // color bars
      else {
        print_status_or_exit(read_int(optarg, &grid_width), c);
        if (grid_width < 1)
          error("grid width must be positive", c);
      }
      break;

    case 'I':
      collect_indexes = false;
      break;

    case 'M':
      print_status_or_exit(get_arg_id(optarg, &arg_id, "hsv=1|hsl=2|rgb=3",
                                      argmatch_add_id_maps),
                           c);
      color_system_mode = atoi(arg_id.c_str());
      break;

    case 'k':
      print_status_or_exit(read_int(optarg, &show_container), c);
      if (show_container < 0 || show_container > 2)
        error("container mode must be between 0 and 2", c);
      break;

    case 'R':
      show_reference = true;
      break;

    case 'q':
      print_status_or_exit(read_int(optarg, &hsl_height), c);
      if (hsl_height < 1 || hsl_height > 2)
        error("hsl height must be between 1 and 2", c);
      break;

    case 'z':
      print_status_or_exit(read_int(optarg, &upright_view), c);
      if (upright_view < 1 || upright_view > 2)
        error("view type must be between 1 and 2", c);
      break;

    case 'f':
      print_status_or_exit(get_arg_id(optarg, &arg_id,
                                      "rgb=1|antiprism=2|decimal=3",
                                      argmatch_add_id_maps),
                           c);
      map_type = atoi(arg_id.c_str());
      break;

    case 'r':
      print_status_or_exit(get_arg_id(optarg, &arg_id,
                                      "cylinder=1|conic=2|hexagonal=3|cube=4",
                                      argmatch_add_id_maps),
                           c);
      container = atoi(arg_id.c_str());
      break;

    case 'b':
      print_status_or_exit(read_double(optarg, &brightness_adj), c);
      if (brightness_adj < -1.0 || brightness_adj > 1.0)
        error("brightness adjustment must be between -1.0 and 1.0", c);
      break;

    case 's':
      print_status_or_exit(read_double_list(optarg, double_parms, 4), c);
      for (double &double_parm : double_parms) {
        if (double_parm <= 0.0)
          error("color centroid saturation curve must be greater than zero", c);
        sat_powers.push_back(double_parm);
      }
      break;

    case 't':
      print_status_or_exit(read_double(optarg, &sat_threshold), c);
      if (sat_threshold < 0.0 || sat_threshold > 1.0)
        error("HSV/HSL threshold must be between 0 and 1", c);
      break;

    case 'u':
      print_status_or_exit(read_double(optarg, &value_advance), c);
      if (value_advance < 0.0 || value_advance > 120.0)
        error("HSV/HSL value advance must be between 0 and 120", c);
      break;

    case 'v':
      print_status_or_exit(read_double_list(optarg, double_parms, 4), c);
      for (double &double_parm : double_parms) {
        if (double_parm < 0.0)
          error("value curve must be greater than or equal to zero", c);
        value_powers.push_back(double_parm);
      }
      break;

    case 'p':
      plot_centroid = true;
      break;

    case 'a':
      print_status_or_exit(get_arg_id(optarg, &arg_id,
                                      "average=1|minimum=2|maximum=3",
                                      argmatch_add_id_maps),
                           c);
      alpha_mode = atoi(arg_id.c_str());
      break;

    case 'c':
      cmy_mode = true;
      break;

    case 'y':
      ryb_mode = true;
      break;

    case 'S':
      seven_mode = true;
      break;

    case 'm':
      print_status_or_exit(map.init(optarg), c);
      break;

    case 'O':
      gfile = optarg;
      break;

    case 'U':
      if (strlen(optarg) > 1 || !strchr("su", *optarg))
        error(msg_str("sort colors is '%s', must be s or u", optarg), c);
      sort_colors = *optarg;
      break;

    case 'Z':
      print_status_or_exit(read_int(optarg, &map_maximum), c);
      if (map_maximum < 0)
        error("maximum map elements to read in must be greater than 0", c);
      break;

    case 'V':
      verbose = true;
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (argc - optind > 0)
    error("too many arguments");

  if (color_system_mode == 3) {
    if (display_type == 1 && plot_centroid)
      warning("only averge blend is displayed on plots in RGB mode", "p");
    else if (display_type == 2)
      warning("only averge blend is displayed on wheel in RGB mode", "p");

    if (sat_powers.size())
      warning("saturation entries are not used in RGB mode", "s");
    if (sat_threshold != 1.0)
      warning("saturation threshold is not used in RGB mode", "t");

    if (value_powers.size())
      warning("value entries are are not used in RGB mode", "v");
    if (value_advance > 0)
      warning("value advance is not used in RGB mode", "u");
  }

  // check these values before filling in defaults
  if (display_type == 3 || display_type == 4) {
    if (plot_centroid)
      warning("plotting centroid colors is only valid in plot or wheel mode",
              "p");
    if (sat_powers.size())
      warning("saturation entries are only valid in plot or wheel mode", "s");
    if (sat_threshold != 1.0)
      warning("saturation threshold is only valid in plot or wheel mode", "t");

    if (value_powers.size())
      warning("value entries are only valid in plot or wheel mode", "v");
    if (value_advance > 0)
      warning("value advance is only valid in plot or wheel mode", "u");
  }

  if (color_system_mode != 3) {
    if (sat_powers.size() == 0 && sat_threshold != 1.0)
      warning("saturation threshold has no effect if -s is not set", "t");
    if (value_powers.size() == 0 && value_advance > 0)
      warning("value advance has no effect if -v is not set", "u");
  }

  // fill in missing sat_powers with -1.0, meaning use centroid saturation
  for (unsigned int i = sat_powers.size(); i < 4; i++)
    sat_powers.push_back(-1.0);

  // fill in missing value_powers with -1.0, meaning use average values
  for (unsigned int i = value_powers.size(); i < 4; i++)
    value_powers.push_back(-1.0);

  if (map_type && display_type != 4) {
    warning("map type is only valid in map mode", "f");
    map_type = 0;
  }
  else {
    // set default map type here
    if (!map_type)
      map_type = 1;
  }

  if (!collect_indexes && display_type != 3) {
    warning("excluding map indexes only valid in grid mode", "I");
    collect_indexes = true;
  }

  if (grid_width && display_type != 3) {
    warning("grid width only valid in grid mode", "w");
    grid_width = 0;
  }

  if (container) {
    if (display_type > 1) {
      warning("container type is only valid in plot mode", "r");
      container = 0;
    }
    if (color_system_mode == 3) {
      warning("container type has no effect in RGB mode", "r");
      container = 0;
    }
  }
  else {
    // set default HSV/HSL container here (HSV/HSL cube)
    if (!container && color_system_mode != 3)
      container = 4;
  }

  if (hsl_height && color_system_mode != 2)
    warning("HSL height adjustment has no effect in RGB or HSV mode", "q");

  if (show_container == 1 &&
      (color_system_mode == 3 || (container == 1 || container == 4)))
    warning("facets are only in HSV or HSL conic or hexagonal container", "k");

  if (show_reference) {
    if (!show_container)
      warning("show reference has no effect when container is not shown", "R");
    if (display_type != 1)
      warning("show reference has no effect when not displaying plot", "R");
    show_reference = false;
  }

  if (upright_view && container == 4) {
    warning("upright view has no effect HSV/HSL cube container", "z");
    upright_view = 0;
  }

  if (upright_view == 2 && color_system_mode != 3)
    warning("view 2 has no effect when not in RGB mode", "z");

  if (seven_mode && (color_system_mode == 3 || container == 4))
    warning("heptagonal view has no effect in RGB mode or cube container", "S");
}

// furnished by Adrian Rossiter
void color_wheel(Geometry &geom, const vector<Color> &cols,
                 const col_util_opts &opts)
{
  // RK - dynamic polygon size
  unsigned int sz = cols.size();
  int polygon_sz = (sz < 120) ? 120 : sz;

  int num_points = polygon_sz / sz;
  int face_sz = 2 * (num_points + 1);
  double sect_ang = 2 * M_PI / sz;
  double ang_inc = sect_ang / num_points;
  for (unsigned int i = 0; i < sz; i++) {
    vector<int> face(face_sz);
    for (int j = 0; j < num_points + 1; j++) {
      double ang = i * sect_ang + j * ang_inc;
      geom.add_vert(Vec3d(cos(ang), sin(ang), 0.0));
      face[j] = geom.verts().size() - 1;
      geom.add_vert(Vec3d(cos(ang) / 2, sin(ang) / 2, 0.0));
      face[face_sz - 1 - j] = geom.verts().size() - 1;
    }
    geom.add_face(face, cols[i]);
  }

  // RK - user gets multiple level control
  double sat_power = opts.sat_powers[0];
  double value_power = opts.value_powers[0];

  for (int lvl = 0; lvl < 4; lvl++) {
    vector<int> face(polygon_sz);
    for (int i = 0; i < polygon_sz; i++) {
      double ang = 2 * M_PI * i / polygon_sz;
      double rad = 0.5 * (4 - lvl) / 4;
      geom.add_vert(Vec3d(rad * cos(ang), rad * sin(ang), 0.01 * lvl));
      face[i] = geom.verts().size() - 1;
    }
    if (opts.color_system_mode == 3) {
      // RK - if RGB mode, only show that blend in all four levels
      Color col =
          col_blend::blend_RGB_centroid(cols, opts.alpha_mode, opts.ryb_mode);
      if (opts.brightness_adj)
        col.set_brightness(opts.brightness_adj);

      if (opts.verbose && !lvl)
        col.dump("[blend]");

      geom.add_face(face, col);
    }
    else {
      // RK - only change the powers if there is a new one not the default
      // keeps the center of the bullseye consistent with the last valid power
      if (opts.sat_powers[lvl] > -1.0)
        sat_power = opts.sat_powers[lvl];
      if (opts.value_powers[lvl] > -1.0)
        value_power = opts.value_powers[lvl];

      Color col = col_blend::blend_HSX_centroid(
          cols, opts.color_system_mode, sat_power, opts.sat_threshold,
          value_power, opts.value_advance, opts.alpha_mode, opts.ryb_mode);
      if (opts.brightness_adj)
        col.set_brightness(opts.brightness_adj);

      if (opts.verbose) {
        string name;
        name = "[blend " + std::to_string(lvl + 1) + "]";
        col.dump(name.c_str());
      }

      geom.add_face(face, col);
    }
  }

  Coloring clrng(&geom);
  geom.add_missing_impl_edges();
  clrng.e_one_col(Color::invisible);
  clrng.v_one_col(Color::invisible);
}

// http://en.wikipedia.org/wiki/HSL_and_HSV#Hue_and_chroma
// hue and chroma output are from 0 to 1
void get_chroma(const Color &col, int container, double &hue, double &chroma)
{
  int R = col[0];
  int G = col[1];
  int B = col[2];

  int M = max(R, max(G, B));
  int m = min(R, min(G, B));
  int C = M - m;

  double H = 0.0;
  if (!C)
    H = 0.0;
  else if (M == R)                      // red
    H = fmod(((double)(G - B) / C), 6); // (g-b/c) mod 6;
  else if (M == G)                      // green
    H = ((double)(B - R) / C) + 2;      // (b-r/c)+2;
  else if (M == B)                      // blue
    H = ((double)(R - G) / C) + 4;      // (r-g/c)+4;
  H *= 60.0;
  if (H < 0)
    H += 360.0;

  if (container == 2) {
    hue = H / 360.0;
    chroma = C / 255.0;
  }
  else if (container == 3) {
    double a = (2 * R - G - B) * 0.5;
    double b = sqrt(3) / 2.0 * (G - B);

    hue = atan2(b, a) / (2 * M_PI);
    if (hue < 0.0)
      hue += 1.0;
    chroma = sqrt(a * a + b * b) / 255.0;
  }
}

// angle represented by 0 to 360 degrees
// input: HSV/HSL angle
// output: angle adjusted for heptagon
double hsx_to_7gon(double angle)
{
  double heptagon_angle = 360.0 / 7.0; // 51.42857143...

  // take 0 to 60 range to 0 to 102.8571429...
  if (angle > 0.0 && angle <= 60.0)
    angle *= heptagon_angle / 30.0;
  else if (angle > 60.0)
    angle = ((angle - 60.0) * heptagon_angle / 60.0) + heptagon_angle * 2;

  return angle;
}

// code to draw cone into heptagon furnished by Adrian Rossiter
double saturation_correction_7gon(double hue, double sat)
{
  double angle = (2 * M_PI / 7) * (0.5 - fmod(7.0 * hue, 1.0));
  double scale = cos(M_PI / 7) / cos(angle);
  return (sat * scale);
}

// return the color that was finally plotted if needed
Color plot_hsx_point(Geometry &geom, Color &col, bool bypass_brightness,
                     const col_util_opts &opts)
{
  if (!col.is_value())
    return col;

  if (opts.brightness_adj && !bypass_brightness)
    col.set_brightness(opts.brightness_adj);

  Vec4d hsxa = col_blend::get_hsxa(col, opts.color_system_mode);

  // heptagonal display overrides RYB distribution
  bool ryb = opts.ryb_mode;
  bool seven_mode_hsl = false;
  int chroma_calc = opts.container;
  if (opts.seven_mode) {
    ryb = false;
    if (chroma_calc == 3) {
      chroma_calc = 2;
      seven_mode_hsl = true;
    }
  }

  double H = (ryb) ? col_blend::hsx_to_ryb(rad2deg(2 * M_PI * hsxa[0])) / 360.0
                   : hsxa[0];
  double S = hsxa[1];

  if (chroma_calc > 1) {
    Color ccol = col;
    if (ryb)
      ccol =
          col_blend::set_hsxa(H, S, hsxa[2], hsxa[3], opts.color_system_mode);
    get_chroma(ccol, chroma_calc, H, S);
  }

  if (opts.seven_mode) {
    // distribute hex oriented hue into 7-way mode
    double angle = rad2deg(2 * M_PI * H);
    H = hsx_to_7gon(angle) / 360.0;
    // draw cone into heptagon
    if (seven_mode_hsl)
      S = saturation_correction_7gon(H, S);
  }

  double angle = 2 * M_PI * H;

  geom.add_vert(Vec3d(S * cos(angle), S * sin(angle), hsxa[2]), col);

  return col;
}

void plot_hsx(Geometry &geom, const vector<Color> &cols,
              const col_util_opts &opts)
{
  for (auto col : cols) {
    // can't be const
    plot_hsx_point(geom, col, false, opts);
  }
}

Geometry make_unit_circle(int polygon_size)
{
  Geometry geom;
  double arc = deg2rad(360.0 / (double)polygon_size);

  double angle = 0.0;
  for (int i = 0; i < polygon_size; i++) {
    geom.add_vert(Vec3d(cos(angle), sin(angle), 0.0));
    geom.add_edge(make_edge(i, (i + 1) % polygon_size));
    angle += arc;
  }

  return (geom);
}

void make_hsv_container(Geometry &geom, const col_util_opts &opts)
{
  int polygon_size = (opts.seven_mode ? 7 : 6);
  geom.append(make_unit_circle(polygon_size));

  // if show_container is 1 do not show facets
  if (opts.show_container < 2)
    return;

  int new_vert = 0;
  if (opts.color_system_mode == 2) {
    geom.add_vert(Vec3d(0.0, 0.0, 0.5));
    new_vert = geom.verts().size() - 1;
    for (int i = 0; i < polygon_size; i++)
      geom.add_edge(make_edge(i, new_vert));
  }
  geom.add_vert(Vec3d(0.0, 0.0, (opts.color_system_mode == 2) ? -0.5 : -1.0));
  new_vert = geom.verts().size() - 1;
  for (int i = 0; i < polygon_size; i++)
    geom.add_edge(make_edge(i, new_vert));
}

Geometry make_hsx_container(const col_util_opts &opts)
{
  Geometry geom;
  if (opts.container == 3) {
    make_hsv_container(geom, opts);
    geom.transform(Trans3d::translate(opts.color_system_mode == 2
                                          ? Vec3d(0.0, 0.0, 0.5)
                                          : Vec3d(0.0, 0.0, 1.0)));
  }
  else {
    geom = make_unit_circle(60);
    geom.transform(
        Trans3d::translate((opts.color_system_mode == 1 || opts.container == 1)
                               ? Vec3d(0.0, 0.0, 1.0)
                               : Vec3d(0.0, 0.0, 0.5)));
    if (opts.container == 1)
      geom.append(make_unit_circle(60));
  }

  Coloring clrng(&geom);
  clrng.e_one_col(Color(1.0, 1.0, 1.0, 0.1));

  clrng.v_one_col(Color::invisible);
  if (opts.show_reference) {
    vector<Color> col(9);
    col[0] = Color(255, 0, 0, 224);
    col[1] = Color(255, 127, 0, 224); // darkorange1
    col[2] = Color(255, 255, 0, 224);
    col[3] = Color(0, 255, 0, 224);
    col[4] = Color(0, 255, 255, 224);
    col[5] = Color(0, 0, 255, 224);
    col[6] = Color(255, 0, 255, 224);
    col[7] = Color(255, 255, 255, 224);
    col[8] = Color(0, 0, 0, 224);
    for (unsigned int i = 0; i < col.size(); i++) {
      if (i == 1 && !(opts.ryb_mode || opts.seven_mode))
        continue;
      // reference points can't be altered by brightness adjust
      plot_hsx_point(geom, col[i], true, opts);
    }
  }

  return (geom);
}

Geometry make_cube(const col_util_opts &opts)
{
  Geometry geom;

  // vertex reference colors
  Color col[8] = {Color()};
  if (opts.color_system_mode == 3) {
    // color the corner vertices as rgb values
    col[0] = Color(255, 255, 255, 224);
    // color in ryb_mode is darkorange1
    col[1] =
        (opts.ryb_mode ? Color(255, 127, 0, 224) : Color(255, 255, 0, 224));
    col[2] = Color(255, 0, 255, 224);
    col[3] = Color(255, 0, 0, 224);
    col[4] = (opts.ryb_mode ? Color(0, 255, 0, 244) : Color(0, 255, 255, 224));
    col[5] = (opts.ryb_mode ? Color(255, 255, 0, 244) : Color(0, 255, 0, 224));
    col[6] = Color(0, 0, 255, 224);
    col[7] = Color(0, 0, 0, 224);
  }
  else {
    // color the corner vertices as hsva values
    col[0].set_hsva(1.0, 1.0, 1.0, 0.875);
    col[1].set_hsva(1.0, 1.0, 0.0, 0.875);
    col[2].set_hsva(1.0, 0.0, 1.0, 0.875);
    col[3].set_hsva(1.0, 0.0, 0.0, 0.875);
    col[4].set_hsva(0.0, 1.0, 1.0, 0.875);
    col[5].set_hsva(0.0, 1.0, 0.0, 0.875);
    col[6].set_hsva(0.0, 0.0, 1.0, 0.875);
    col[7].set_hsva(0.0, 0.0, 0.0, 0.875);
  }

  geom.add_vert(Vec3d(1, 1, 1), col[0]); // 0
  geom.add_vert(Vec3d(1, 1, 0), col[1]); // 1
  geom.add_vert(Vec3d(1, 0, 1), col[2]); // 2
  geom.add_vert(Vec3d(1, 0, 0), col[3]); // 3
  geom.add_vert(Vec3d(0, 1, 1), col[4]); // 4
  geom.add_vert(Vec3d(0, 1, 0), col[5]); // 5
  geom.add_vert(Vec3d(0, 0, 1), col[6]); // 6
  geom.add_vert(Vec3d(0, 0, 0), col[7]); // 7

  geom.add_edge(make_edge(0, 1));
  geom.add_edge(make_edge(0, 2));
  geom.add_edge(make_edge(0, 4));
  geom.add_edge(make_edge(1, 3));
  geom.add_edge(make_edge(2, 3));
  geom.add_edge(make_edge(2, 6));
  geom.add_edge(make_edge(4, 6));
  geom.add_edge(make_edge(1, 5));
  geom.add_edge(make_edge(4, 5));
  geom.add_edge(make_edge(7, 3));
  geom.add_edge(make_edge(7, 5));
  geom.add_edge(make_edge(7, 6));

  Coloring clrng(&geom);
  if (!opts.show_reference)
    clrng.v_one_col(Color::invisible);
  else {
    // if ryb mode, have to add cyan in on mid-edge
    if (opts.ryb_mode && opts.color_system_mode == 3)
      geom.add_vert(Vec3d(0, 0.5, 1), Color(0, 255, 255, 224));
  }
  clrng.e_one_col(Color(1.0, 1.0, 1.0, 0.1));

  return geom;
}

// return the color that was finally plotted if needed
Color plot_rgb_point(Geometry &geom, Color &col, const col_util_opts &opts)
{
  if (!col.is_value())
    return col;

  if (opts.brightness_adj)
    col.set_brightness(opts.brightness_adj);

  // rcol is relocated color, not to be displayed
  Color rcol = col;
  if (opts.ryb_mode) {
    // only need hue so algorithm doesn't matter
    Vec4d hsxa = rcol.get_hsva();
    hsxa[0] = col_blend::hsx_to_ryb(rad2deg(2 * M_PI * hsxa[0])) / 360.0;
    rcol.set_hsva(hsxa);
  }

  if (opts.color_system_mode == 3)
    geom.add_vert(Vec3d(rcol[0] / 255.0, rcol[1] / 255.0, rcol[2] / 255.0),
                  col);
  else {
    // borrow function for plotting hsv cube
    // hsxa is relocated color, not to be displayed
    Vec4d hsxa = col_blend::get_hsxa(col, opts.color_system_mode);
    geom.add_vert(Vec3d(hsxa[0], hsxa[1], hsxa[2]), col);
  }

  return col;
}

void plot_rgb_cube(Geometry &geom, const vector<Color> &cols,
                   const col_util_opts &opts)
{
  for (auto col : cols) {
    // can't be const
    plot_rgb_point(geom, col, opts);
  }
}

Geometry make_square(const int &width)
{
  Geometry geom;

  geom.add_vert(Vec3d(0, 0, 0));     // 0
  geom.add_vert(Vec3d(0, 1, 0));     // 1
  geom.add_vert(Vec3d(width, 1, 0)); // 2
  geom.add_vert(Vec3d(width, 0, 0)); // 3

  geom.add_edge(make_edge(0, 1));
  geom.add_edge(make_edge(1, 2));
  geom.add_edge(make_edge(2, 3));
  geom.add_edge(make_edge(3, 0));

  vector<int> face;
  face.push_back(0);
  face.push_back(1);
  face.push_back(2);
  face.push_back(3);
  geom.add_face(face);

  Coloring clrng(&geom);
  clrng.v_one_col(Color::invisible);
  clrng.e_one_col(Color::invisible);

  return geom;
}

void color_grid(Geometry &geom, const vector<Color> &cols,
                const col_util_opts &opts)
{
  bool color_bars = (opts.grid_width == -1 ? true : false);
  int width = opts.grid_width;
  if (color_bars)
    width = cols.size();

  Geometry sgeom;
  sgeom.append(make_square((color_bars ? width : 1)));

  unsigned int cols_sz = cols.size();
  unsigned int dim1;
  unsigned int dim2;

  if (!width) {
    dim2 = (int)ceil(sqrt(cols_sz));
    dim1 = dim2;
  }
  else {
    dim2 = width;
    dim1 = cols_sz / width;
    if (cols_sz - (dim1 * width) > 0)
      dim1++;
  }

  unsigned int k = 0;
  for (unsigned int i = 0; i < dim1; i++) {
    for (unsigned int j = 0; j < dim2; j++) {
      Geometry tgeom = sgeom;
      if (k < cols_sz && cols[k].is_index()) {
        tgeom.add_vert(Vec3d(0.5, 0.45, 0.0), Color(0.0, 0.0, 0.0));
        tgeom.add_vert(Vec3d(0.5, 0.55, 0.0), Color(1.0, 1.0, 1.0));
        tgeom.add_edge(make_edge(4, 5), Color(0.5, 0.5, 0.5));
      }
      tgeom.transform(Trans3d::translate(Vec3d(i, j, 0)));
      Color c =
          (k >= cols_sz ? Color(Color::invisible)
                        : (cols[k].is_index() ? cols[k].get_index() : cols[k]));
      k++;
      tgeom.colors(FACES).set(0, c);
      geom.append(tgeom);
      tgeom.clear_all();
    }
  }
  // use anti::epsilon since there is no opts epsilon
  merge_coincident_elements(geom, "ve", anti::epsilon);
  geom.transform(Trans3d::rotate(Vec3d(0.0, 0.0, deg2rad(-90.0))));
}

bool cmp_col(const Color &a, const Color &b)
{
  bool ret = false;

  Vec4d hsva_a = a.get_hsva();
  Vec4d hsva_b = b.get_hsva();
  for (unsigned int i = 0; i < 4; i++) {
    if (hsva_a[i] != hsva_b[i]) {
      ret = (hsva_a[i] < hsva_b[i]);
      break;
    }
  }

  return ret;
}

class col_cmp {
public:
  col_cmp() = default;
  bool operator()(const Color &a, const Color &b) const
  {
    return cmp_col(a, b);
  }
};

void collect_col(vector<Color> &cols, const Color &col, bool collect_indexes)
{
  if (!col.is_set())
    return;
  else if (!collect_indexes && col.is_index())
    return;
  cols.push_back(col);
}

void collect_cols_from_geom(const Geometry &geom, vector<Color> &cols,
                            bool collect_indexes)
{
  for (unsigned int i = 0; i < geom.verts().size(); i++)
    collect_col(cols, geom.colors(VERTS).get(i), collect_indexes);
  for (unsigned int i = 0; i < geom.edges().size(); i++)
    collect_col(cols, geom.colors(EDGES).get(i), collect_indexes);
  for (unsigned int i = 0; i < geom.faces().size(); i++)
    collect_col(cols, geom.colors(FACES).get(i), collect_indexes);
}

void collect_cols(vector<Color> &cols, col_util_opts &opts)
{
  int map_sz = opts.map.effective_size();
  bool open_ended_map = (map_sz >= std::numeric_limits<int>::max());
  // map size priority:
  // -Z given
  // default for unlimited map (256)
  // use actual map size
  int max_map_sz =
      (opts.map_maximum) ? opts.map_maximum : ((open_ended_map) ? 256 : map_sz);

  if (open_ended_map)
    opts.warning(msg_str("map list: only %d out of %d map entries read in",
                         max_map_sz, map_sz),
                 'Z');
  for (int j = 0; j < max_map_sz; j++)
    collect_col(cols, opts.map.get_col(j), opts.collect_indexes);

  // file
  if (opts.gfile.length()) {
    Geometry geom;
    opts.read_or_error(geom, opts.gfile);
    collect_cols_from_geom(geom, cols, opts.collect_indexes);
  }

  if (opts.sort_colors) {
    sort(cols.begin(), cols.end(), col_cmp());
    if (opts.sort_colors == 'u') {
      auto ci = unique(cols.begin(), cols.end());
      cols.resize(ci - cols.begin());
    }
  }

  if (opts.cmy_mode)
    for (auto &col : cols)
      col = col_blend::rgb_complement(col, opts.ryb_mode);

  for (auto &col : cols) {
    if (!col.is_index() && col[3] < 32) {
      opts.warning("colors with a very small alpha values may not be seen");
      break;
    }
  }

  // grid may have indexes
  if (opts.display_type == 3) {
    for (auto &col : cols) {
      if (col.is_index()) {
        opts.warning("color indexes detected. unmapped cells will result");
        break;
      }
    }
  }
}

int main(int argc, char *argv[])
{
  col_util_opts opts;
  opts.process_command_line(argc, argv);

  Geometry geom;

  vector<Color> cols;
  collect_cols(cols, opts);

  if (!cols.size())
    opts.error("found no color values to plot");

  if (cols.size() < 3) {
    for (unsigned int i = 0; i < 4; i++) {
      if (opts.value_powers[i] > -1.0) {
        opts.warning("subtractive colors only works on three colors or more",
                     "v");
        break;
      }
    }
  }

  // output map file
  if (opts.display_type == 4) {
    FILE *ofile = stdout; // write to stdout by default
    if (opts.ofile != "") {
      ofile = fopen(opts.ofile.c_str(), "w");
      if (ofile == nullptr)
        opts.error("could not open output file \'" + opts.ofile + "\'");
    }
    for (unsigned int i = 0; i < cols.size(); i++) {
      if (opts.map_type == 3) {
        if (cols[i].is_value()) {
          Vec4d c = cols[i].get_vec4d();
          fprintf(ofile, "%g%s %g%s %g%s %g%s\n", c[0],
                  (c[0] == 1.0 || c[0] == 0.0) ? ".0" : "", c[1],
                  (c[1] == 1.0 || c[1] == 0.0) ? ".0" : "", c[2],
                  (c[2] == 1.0 || c[2] == 0.0) ? ".0" : "", c[3],
                  (c[3] == 1.0 || c[3] == 0.0) ? ".0" : "");
        }
        else
          fprintf(ofile, "%3d\n", cols[i].get_index());
      }
      else {
        if (opts.map_type == 2)
          fprintf(ofile, "%-5d = ", i);
        string buffer = "";
        if (cols[i][3] != 255)
          buffer = msg_str("%3d", cols[i][3]);
        if (cols[i].is_value())
          fprintf(ofile, "%3d %3d %3d %s\n", cols[i][0], cols[i][1], cols[i][2],
                  buffer.c_str());
        else
          fprintf(ofile, "%3d\n", cols[i].get_index());
      }
    }
    if (opts.ofile != "")
      fclose(ofile);
  }
  // else it is a plot/model
  else {
    // a plot
    if (opts.display_type == 1) {
      if (opts.color_system_mode == 1 || opts.color_system_mode == 2) {
        if (opts.container < 4)
          plot_hsx(geom, cols, opts);
        else
          // borrow rgb cube for plotting hsv cube container == 4
          plot_rgb_cube(geom, cols, opts);

        if (opts.plot_centroid) {
          for (unsigned int i = 0; i < 4; i++) {
            if (!i || opts.sat_powers[i] > -1.0 ||
                opts.value_powers[i] > -1.0) {
              Color col = col_blend::blend_HSX_centroid(
                  cols, opts.color_system_mode, opts.sat_powers[i],
                  opts.sat_threshold, opts.value_powers[i], opts.value_advance,
                  opts.alpha_mode, opts.ryb_mode);
              Color blend;
              if (opts.container < 4)
                blend = plot_hsx_point(geom, col, false, opts);
              else
                blend = plot_rgb_point(geom, col, opts);
              if (opts.verbose) {
                string name;
                name = "[blend " + std::to_string(i + 1) + "]";
                blend.dump(name.c_str());
              }
            }
          }
        }
      }
      else if (opts.color_system_mode == 3) {
        plot_rgb_cube(geom, cols, opts);

        if (opts.plot_centroid) {
          Color col = col_blend::blend_RGB_centroid(cols, opts.alpha_mode,
                                                    opts.ryb_mode);
          Color blend = plot_rgb_point(geom, col, opts);
          if (opts.verbose)
            blend.dump("[blend]");
        }
      }

      // container
      if (opts.show_container) {
        if (opts.color_system_mode == 3 || opts.container == 4)
          geom.append(make_cube(opts));
        else if (opts.color_system_mode == 1 || opts.color_system_mode == 2)
          geom.append(make_hsx_container(opts));
      }

      // view for cube
      if (opts.color_system_mode == 3) {
        if (opts.upright_view != 2) {
          geom.transform(Trans3d::rotate(deg2rad(45), 0, 0));
          geom.transform(Trans3d::rotate(0, -asin(1 / sqrt(3)), 0));
        }
      }

      // HSL height, before turning upright
      if (opts.color_system_mode == 2 && opts.hsl_height) {
        geom.transform(Trans3d::scale(1.0, 1.0, 2.0));
        // if option to make 90 degree dihedral permeters
        if (opts.hsl_height == 2 && opts.container) {
          double scale =
              1 / (2 * tan(M_PI * (1.0 / (opts.seven_mode ? 7.0 : 6.0))));
          // scale of heptagon to 1 is .8677674782351162 =
          // cos(Pi*5/14)+cos(Pi*5/14) or cos(Pi*3/14)/sin(Pi*5/14)
          if (opts.seven_mode)
            scale *= cos(M_PI * 5.0 / 14.0) * 2.0;
          geom.transform(Trans3d::scale(1.0, 1.0, fabs(scale)));
        }
      }

      // view upright for all models if option 1
      if (opts.upright_view == 1)
        geom.transform(Trans3d::rotate(deg2rad(-90), 0, 0));
    }
    else if (opts.display_type == 2)
      color_wheel(geom, cols, opts);
    else if (opts.display_type == 3)
      color_grid(geom, cols, opts);

    opts.write_or_error(geom, opts.ofile);
  }

  return 0;
}
