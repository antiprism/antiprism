/*
   Copyright (c) 2017-2020, Roger Kaufman

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

/* $Id: kaleido.h,v 3.27 2002-01-06 16:21:30+02 rl Exp $ */
/*
 *****************************************************************************
 * kaleido
 *
 * Kaleidoscopic construction of uniform polyhedra
 * Copyright © 1991-2002 Dr. Zvi Har'El <rl@math.technion.ac.il>
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * 3. The end-user documentation included with the redistribution,
 *    if any, must include the following acknowledgment:
 *  "This product includes software developed by
 *   Dr. Zvi Har'El (http://www.math.technion.ac.il/~rl/)."
 *    Alternately, this acknowledgment may appear in the software itself,
 *    if and wherever such third-party acknowledgments normally appear.
 *
 * This software is provided 'as-is', without any express or implied
 * warranty.  In no event will the author be held liable for any
 * damages arising from the use of this software.
 *
 * Author:
 *  Dr. Zvi Har'El,
 *  Technion, Israel Institue of Technology,
 *  Haifa 32000, Israel.
 *  E-Mail: rl@math.technion.ac.il
 *****************************************************************************
 */

#include "../base/antiprism.h"

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>

using std::set;
using std::swap;
using std::vector;

using namespace anti;

// definitions and structures used to be in kaleido.h
#ifndef MAXLONG
#define MAXLONG 0x7FFFFFFF
#endif

#ifndef MAXDIGITS
#define MAXDIGITS 10 /* (int)log10((double)MAXLONG) + 1 */
#endif

#ifndef DBL_EPSILON
#define DBL_EPSILON 2.2204460492503131e-16
#endif

#define BIG_EPSILON 3e-2
#define DEG (180 / M_PI)
#define AZ 0 // was M_PI / 7  /* axis azimuth */
#define EL 0 // was M_PI / 17 /* axis elevation */

#define Free(lvalue)                                                           \
  {                                                                            \
    if (lvalue) {                                                              \
      free((char *)lvalue);                                                    \
      lvalue = 0;                                                              \
    }                                                                          \
  }

#define Malloc(lvalue, n, type)                                                \
  {                                                                            \
    if (!(lvalue = (type *)malloc((n) * sizeof(type))))                        \
      Err(0);                                                                  \
  }

#define Realloc(lvalue, n, type)                                               \
  {                                                                            \
    if (!(lvalue = (type *)realloc(lvalue, (n) * sizeof(type))))               \
      Err(0);                                                                  \
  }

#define Calloc(lvalue, n, type)                                                \
  {                                                                            \
    if (!(lvalue = (type *)calloc(n, sizeof(type))))                           \
      Err(0);                                                                  \
  }

#define Matalloc(lvalue, n, m, type)                                           \
  {                                                                            \
    if (!(lvalue = (type **)matalloc(n, (m) * sizeof(type))))                  \
      Err(0);                                                                  \
  }

#define Sprintfrac(lvalue, x)                                                  \
  {                                                                            \
    if (!(lvalue = sprintfrac(x)))                                             \
      return 0;                                                                \
  }

#define numerator(x) (frac(x), frax.n)
#define denominator(x) (frac(x), frax.d)
#define compli(x) (frac(x), (double)frax.n / (frax.n - frax.d))

#define GR_DIRH_IDX 74

typedef struct {
  double x, y, z;
} Vector;

Vector x, y, z;

typedef struct {
  long n, d;
} Fraction;

Fraction frax;

typedef struct {
  /* NOTE: some of the int's can be replaced by short's, char's,
   or even bit fields, at the expense of readability!!!*/
  int index;       /* index to the standard list, the array uniform[] */
  int N;           /* number of faces types (atmost 5)*/
  int M;           /* vertex valency  (may be big for dihedral polyhedra) */
  int V;           /* vertex count */
  int E;           /* edge count */
  int F;           /* face count */
  int D;           /* density */
  int chi;         /* Euler characteristic */
  int g;           /* order of symmetry group */
  int K;           /* symmetry type: D=2, T=3, O=4, I=5 */
  int hemi;        /* flag hemi polyhedron */
  int onesided;    /* flag onesided polyhedron */
  int even;        /* removed face in pqr| */
  int *Fi;         /* face counts by type (array N)*/
  int *rot;        /* vertex configuration (array M of 0..N-1) */
  int *snub;       /* snub triangle configuration (array M of 0..1) */
  int *firstrot;   /* temporary for vertex generation (array V) */
  int *anti;       /* temporary for direction of ideal vertices (array E) */
  int *ftype;      /* face types (array F) */
  int **e;         /* edges (matrix 2 x E of 0..V-1)*/
  int **dual_e;    /* dual edges (matrix 2 x E of 0..F-1)*/
  int **incid;     /* vertex-face incidence (matrix M x V of 0..F-1)*/
  int **adj;       /* vertex-vertex adjacency (matrix M x V of 0..V-1)*/
  double p[4];     /* p, q and r; |=0 */
  double minr;     /* smallest nonzero inradius */
  double gon;      /* basis type for dihedral polyhedra */
  double *n;       /* number of side of a face of each type (array N) */
  double *m;       /* number of faces at a vertex of each type (array N) */
  double *gamma;   /* fundamental angles in radians (array N) */
  char *polyform;  /* printable Wythoff symbol */
  char *config;    /* printable vertex configuration */
  char *name;      /* name, standard or manifuctured */
  char *dual_name; /* dual name, standard or manifuctured */
  Vector *v;       /* vertex coordinates (array V) */
  Vector *f;       /* face coordinates (array F)*/
} Polyhedron;

typedef struct { /* See uniform.h for explanation of the fields */
  const char *Wythoff, *name, *dual;
  short Coxeter, Wenninger;
} Uniform;

/*
 *****************************************************************************
 * List of Uniform Polyhedra and Their Kaleidoscopic Formulae
 * ==========================================================
 *
 * Each entry contains the following items:
 *
 * 1) Wythoff symbol.
 * 2) Polyhedron name.
 * 3) Dual name.
 * 4) Coxeter &al. reference figure.
 * 5) Wenninger reference figure.
 *
 * Notes:
 *
 * (1) Cundy&Roulette's trapezohedron has been renamed to
 *  deltohedron, as its faces are deltoids, not trapezoids.
 * (2) The names of the non-dihedral polyhedra are those
 *  which appear in Wenninger (1984). Some of them are
 *  slightly modified versions of those in Wenninger (1971).
 *
 * References:
 *
 * Coxeter, H.S.M., Longuet-Higgins, M.S. & Miller, J.C.P.,
 *  Uniform polyhedra, Phil. Trans. Royal Soc. London, Ser. A,
 *  246 (1953), 401-409.
 * Cundy, H.M. & Rollett, A.P.,
 *  "Mathematical Models", 3rd Ed., Tarquin, 1981.
 * Har'El, Z.
 *  Unifom solution for uniform polyhedra, Geometriae Dedicata,
 *  47 (1993), 57-110.
 * Wenninger, M.J.,
 *  "Polyhedron Models", Cambridge University Press, 1971.
 *  "Dual Models", Cambridge University Press, 1984.
 *
 *****************************************************************************
 */

// clang-format off
Uniform uniform_list[] = {

/*******************************************************************************
*  Tetrahedral Schwarz Triangles
*******************************************************************************/

/* (2 3 3) (T1) */

/*1*/ {"3|2 3","tetrahedron",
   "tetrahedron",15,1},
/*2*/ {"2 3|3","truncated tetrahedron",
   "triakistetrahedron",16,6},

/* (3/2 3 3) (T2) */

/*3*/ {"3/2 3|3","octahemioctahedron",
   "octahemioctacron",37,68},

/* (3/2 2 3) (T3) */

/*4*/ {"3/2 3|2","tetrahemihexahedron",
   "tetrahemihexacron",36,67},

/******************************************************************************
*  Octahedral Schwarz Triangles
******************************************************************************/

/* (2 3 4) (O1) */

/*5*/ {"4|2 3","octahedron",
   "cube",17,2},
/*6*/ {"3|2 4","cube",
   "octahedron",18,3},
/*7*/ {"2|3 4","cuboctahedron",
   "rhombic dodecahedron",19,11},
/*8*/ {"2 4|3","truncated octahedron",
   "tetrakishexahedron",20,7},
/*9*/ {"2 3|4","truncated cube",
   "triakisoctahedron",21,8},
/*10*/ {"3 4|2","rhombicuboctahedron",
   "deltoidal icositetrahedron",22,13},
/*11*/ {"2 3 4|","truncated cuboctahedron",
   "disdyakisdodecahedron",23,15},
/*12*/ {"|2 3 4","snub cube",
   "pentagonal icositetrahedron",24,17},

/* (3/2 4 4) (O2b) */

/*13*/ {"3/2 4|4","small cubicuboctahedron",
   "small hexacronic icositetrahedron",38,69},

/* (4/3 3 4) (O4) */

/*14*/ {"3 4|4/3","great cubicuboctahedron",
   "great hexacronic icositetrahedron",50,77},
/*15*/ {"4/3 4|3","cubohemioctahedron",
   "hexahemioctacron",51,78},
/*16*/ {"4/3 3 4|","cubitruncated cuboctahedron",
   "tetradyakishexahedron",52,79},

/* (3/2 2 4) (O5) */

/*17*/ {"3/2 4|2","great rhombicuboctahedron",
   "great deltoidal icositetrahedron",59,85},
/*18*/ {"3/2 2 4|","small rhombihexahedron",
   "small rhombihexacron",60,86},

/* (4/3 2 3) (O7) */

/*19*/ {"2 3|4/3","stellated truncated hexahedron",
   "great triakisoctahedron",66,92},
/*20*/ {"4/3 2 3|","great truncated cuboctahedron",
   "great disdyakisdodecahedron",67,93},

/* (4/3 3/2 2) (O11) */

/*21*/ {"4/3 3/2 2|","great rhombihexahedron",
   "great rhombihexacron",82,103},

/******************************************************************************
*  Icosahedral Schwarz Triangles
******************************************************************************/

/* (2 3 5) (I1) */

/*22*/ {"5|2 3","icosahedron",
   "dodecahedron",25,4},
/*23*/ {"3|2 5","dodecahedron",
   "icosahedron",26,5},
/*24*/ {"2|3 5","icosidodecahedron",
   "rhombic triacontahedron",28,12},
/*25*/ {"2 5|3","truncated icosahedron",
   "pentakisdodecahedron",27,9},
/*26*/ {"2 3|5","truncated dodecahedron",
   "triakisicosahedron",29,10},
/*27*/ {"3 5|2","rhombicosidodecahedron",
   "deltoidal hexecontahedron",30,14},
/*28*/ {"2 3 5|","truncated icosidodecahedron",
   "disdyakistriacontahedron",31,16},
/*29*/ {"|2 3 5","snub dodecahedron",
   "pentagonal hexecontahedron",32,18},

/* (5/2 3 3) (I2a) */

/*30*/ {"3|5/2 3","small ditrigonal icosidodecahedron",
   "small triambic icosahedron",39,70},
/*31*/ {"5/2 3|3","small icosicosidodecahedron",
   "small icosacronic hexecontahedron",40,71},
/*32*/ {"|5/2 3 3","small snub icosicosidodecahedron",
   "small hexagonal hexecontahedron",41,110},

/* (3/2 5 5) (I2b) */

/*33*/ {"3/2 5|5","small dodecicosidodecahedron",
   "small dodecacronic hexecontahedron",42,72},

/* (2 5/2 5) (I3) */

/*34*/ {"5|2 5/2","small stellated dodecahedron",
   "great dodecahedron",43,20},
/*35*/ {"5/2|2 5","great dodecahedron",
   "small stellated dodecahedron",44,21},
/*36*/ {"2|5/2 5","great dodecadodecahedron",
   "medial rhombic triacontahedron",45,73},
/*37*/ {"2 5/2|5","truncated great dodecahedron",
   "small stellapentakisdodecahedron",47,75},
/*38*/ {"5/2 5|2","rhombidodecadodecahedron",
   "medial deltoidal hexecontahedron",48,76},
/*39*/ {"2 5/2 5|","small rhombidodecahedron",
   "small rhombidodecacron",46,74},
/*40*/ {"|2 5/2 5","snub dodecadodecahedron",
   "medial pentagonal hexecontahedron",49,111},

/* (5/3 3 5) (I4) */

/*41*/ {"3|5/3 5","ditrigonal dodecadodecahedron",
   "medial triambic icosahedron",53,80},
/*42*/ {"3 5|5/3","great ditrigonal dodecicosidodecahedron",
   "great ditrigonal dodecacronic hexecontahedron",54,81},
/*43*/ {"5/3 3|5","small ditrigonal dodecicosidodecahedron",
   "small ditrigonal dodecacronic hexecontahedron",55,82},
/*44*/ {"5/3 5|3","icosidodecadodecahedron",
   "medial icosacronic hexecontahedron",56,83},
/*45*/ {"5/3 3 5|","icositruncated dodecadodecahedron",
   "tridyakisicosahedron",57,84},
/*46*/ {"|5/3 3 5","snub icosidodecadodecahedron",
   "medial hexagonal hexecontahedron",58,112},

/* (3/2 3 5) (I6b) */

/*47*/ {"3/2|3 5","great ditrigonal icosidodecahedron",
   "great triambic icosahedron",61,87},
/*48*/ {"3/2 5|3","great icosicosidodecahedron",
   "great icosacronic hexecontahedron",62,88},
/*49*/ {"3/2 3|5","small icosihemidodecahedron",
   "small icosihemidodecacron",63,89},
/*50*/ {"3/2 3 5|","small dodecicosahedron",
   "small dodecicosacron",64,90},

/* (5/4 5 5) (I6c) */

/*51*/ {"5/4 5|5","small dodecahemidodecahedron",
   "small dodecahemidodecacron",65,91},

/* (2 5/2 3) (I7) */

/*52*/ {"3|2 5/2","great stellated dodecahedron",
   "great icosahedron",68,22},
/*53*/ {"5/2|2 3","great icosahedron",
   "great stellated dodecahedron",69,41},
/*54*/ {"2|5/2 3","great icosidodecahedron",
   "great rhombic triacontahedron",70,94},
/*55*/ {"2 5/2|3","great truncated icosahedron",
   "great stellapentakisdodecahedron",71,95},
/*56*/ {"2 5/2 3|","rhombicosahedron",
   "rhombicosacron",72,96},
/*57*/ {"|2 5/2 3","great snub icosidodecahedron",
   "great pentagonal hexecontahedron",73,113},

/* (5/3 2 5) (I9) */

/*58*/ {"2 5|5/3","small stellated truncated dodecahedron",
   "great pentakisdodekahedron",74,97},
/*59*/ {"5/3 2 5|","truncated dodecadodecahedron",
   "medial disdyakistriacontahedron",75,98},
/*60*/ {"|5/3 2 5","inverted snub dodecadodecahedron",
   "medial inverted pentagonal hexecontahedron",76,114},

/* (5/3 5/2 3) (I10a) */

/*61*/ {"5/2 3|5/3","great dodecicosidodecahedron",
   "great dodecacronic hexecontahedron",77,99},
/*62*/ {"5/3 5/2|3","small dodecahemicosahedron",
   "small dodecahemicosacron",78,100},
/*63*/ {"5/3 5/2 3|","great dodecicosahedron",
   "great dodecicosacron",79,101},
/*64*/ {"|5/3 5/2 3","great snub dodecicosidodecahedron",
   "great hexagonal hexecontahedron",80,115},

/* (5/4 3 5) (I10b) */

/*65*/ {"5/4 5|3","great dodecahemicosahedron",
   "great dodecahemicosacron",81,102},

/* (5/3 2 3) (I13) */

/*66*/ {"2 3|5/3","great stellated truncated dodecahedron",
   "great triakisicosahedron",83,104},
/*67*/ {"5/3 3|2","great rhombicosidodecahedron",
   "great deltoidal hexecontahedron",84,105},
/*68*/ {"5/3 2 3|","great truncated icosidodecahedron",
   "great disdyakistriacontahedron",87,108},
/*69*/ {"|5/3 2 3","great inverted snub icosidodecahedron",
   "great inverted pentagonal hexecontahedron",88,116},

/* (5/3 5/3 5/2) (I18a) */

/*70*/ {"5/3 5/2|5/3","great dodecahemidodecahedron",
   "great dodecahemidodecacron",86,107},

/* (3/2 5/3 3) (I18b) */

/*71*/ {"3/2 3|5/3","great icosihemidodecahedron",
   "great icosihemidodecacron",85,106},

/* (3/2 3/2 5/3) (I22) */

/*72*/ {"|3/2 3/2 5/2","small retrosnub icosicosidodecahedron",
   "small hexagrammic hexecontahedron",91,118},

/* (3/2 5/3 2) (I23) */

/*73*/ {"3/2 5/3 2|","great rhombidodecahedron",
   "great rhombidodecacron",89,109},
/*74*/ {"|3/2 5/3 2","great retrosnub icosidodecahedron",
   "great pentagrammic hexecontahedron",90,117},

/******************************************************************************
*  Last But Not Least
******************************************************************************/

/*75*/ {"3/2 5/3 3 5/2","great dirhombicosidodecahedron",
   "great dirhombicosidodecacron",92,119},

/******************************************************************************
*  Dihedral Schwarz Triangles (D5 only)
******************************************************************************/

/* (2 2 5) (D1/5) */

/*76*/ {"2 5|2","pentagonal prism",
   "pentagonal dipyramid",0,0},
/*77*/ {"|2 2 5","pentagonal antiprism",
   "pentagonal deltohedron",0,0},

/* (2 2 5/2) (D2/5) */

/*78*/ {"2 5/2|2","pentagrammic prism",
   "pentagrammic dipyramid",0,0},
/*79*/ {"|2 2 5/2","pentagrammic antiprism",
   "pentagrammic deltohedron",0,0},

/* (5/3 2 2) (D3/5) */

/*80*/ {"|2 2 5/3","pentagrammic crossed antiprism",
   "pentagrammic concave deltohedron",0,0}

};
// clang-format on

int get_uniform_list(Uniform **uniform)
{
  *uniform = uniform_list;
  return sizeof(uniform_list) / sizeof(uniform_list[0]);
}

class kaleido_opts : public ProgramOpts {
public:
  string ofile;

  string symbol;
  int just_list;
  int need_coordinates;
  int need_approx;
  int model;
  int base;
  double azimuth;
  double elevation;
  double freeze;
  int sig_digits;

  kaleido_opts()
      : ProgramOpts("kaleido"), symbol(""), just_list(0), need_coordinates(0),
        need_approx(0), model(1), base(1), azimuth(AZ), elevation(EL),
        freeze(0), sig_digits(DEF_SIG_DGTS)
  {
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

void kaleido_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] [Symbol]

Antiprism port of Kaleido that lists Uniform information, or outputs models
Kaleidoscopic Construction of Uniform Polyhedra, $Revision: 3.27 $
Copyright © 1991-2002 Dr. Zvi Har'El <rl@math.technion.ac.il>

Symbol may be a number Meader index (e.g u75), Kaleido index (e.g k75)
Coxeter index (e.g c75), Wenninger index (e.g w75). There are 80 models
which can be generated by indexes. Or the symbol may be a Wythoff symbol
in quotes (e.g "3|2 3" , "5/3 5/2|3", "|2 2 5/3", "|2 2 9/4")
Entering a listing option with symbol 0 will list all

Options
%s
  -w <opt>  off=1, vrml=2  (default: off)
  -b <opt>  base=1, dual=2 (default: base)
  -d <dgts> number of significant digits (default %d) or if negative
            then the number of digits after the decimal point
  -o <file> write output to file (default: write to standard output)

Listings (use 0 for Symbol to list all)
  -l        list polyhedron names, symbols and reference figures only
  -v        print vertex and face coordinates
  -x        print successive approximations

)",
          prog_name(), help_ver_text, DEF_SIG_DGTS);
}
// RK - rotation doesn't work well
//  -a <deg>  rotation axis azimouth  (default: 0)
//  -e <deg>  rotation axis elevation (default: 0)
//  -f <deg>  angle of rotation until freezing (default: 0)

void kaleido_opts::process_command_line(int argc, char **argv)
{
  char c;
  char errmsg[MSG_SZ] = {0};
  bool sig_digits_set = false;

  string arg_id;

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv, ":hw:b:a:e:f:lvxd:o:")) != -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'w':
      print_status_or_exit(
          get_arg_id(optarg, &arg_id, "off=1|vrml=2", argmatch_add_id_maps), c);
      if (arg_id == "")
        error(errmsg);
      model = atoi(arg_id.c_str());
      break;

    case 'b':
      print_status_or_exit(
          get_arg_id(optarg, &arg_id, "base=1|dual=2", argmatch_add_id_maps),
          c);
      if (arg_id == "")
        error(errmsg);
      base = atoi(arg_id.c_str());
      break;

    case 'a':
      azimuth = atof(optarg) / DEG;
      break;

    case 'e':
      elevation = atof(optarg) / DEG;
      break;

    case 'f':
      freeze = atof(optarg) / DEG;
      break;

    case 'l':
      just_list = 1;
      break;

    case 'v':
      need_coordinates = 1;
      break;

    case 'x':
      need_approx = 1;
      break;

    case 'd':
      sig_digits_set = true;
      print_status_or_exit(read_int(optarg, &sig_digits), c);
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  if (argc - optind == 0)
    error("no argments");
  else if (argc - optind > 1)
    error("too many arguments");
  else if (argc - optind == 1)
    symbol = argv[optind];

  // just_list suppresses need_coordinates. override.
  if (just_list && need_coordinates)
    just_list = false;

  // original code had significant digits at 6 for listings
  if (just_list || need_coordinates || need_approx) {
    if (!sig_digits_set) {
      sig_digits = 6;
      warning("for listings, default significant digits is 6");
    }
  }
  else {
    if (symbol == "0")
      error("can only output one model at at time\n", 'w');
  }
}

// from various libraries

// RK - putchar() replace with fprintf(fp);

// RK - more() becomes a noop
void more() {}

// RK - rewrote Err()
void Err(const char *errmsg)
{
  fprintf(stderr, "kaleido: error: %s\n", errmsg);
  exit(0);
}

void printErrorMessage(int error_number, string func_name)
{
  char buffer[MSG_SZ] = {0};
  snprintf(buffer, sizeof(buffer), "(%s) errno: %d %s", func_name.c_str(),
           error_number, strerror(error_number));
  Err(buffer);
}

int mod(int i, int j) { return (i %= j) >= 0 ? i : j < 0 ? i - j : i + j; }

double dot(Vector a, Vector b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

Vector scale(double k, Vector a)
{
  a.x *= k;
  a.y *= k;
  a.z *= k;
  return a;
}

Vector diff(Vector a, Vector b)
{
  a.x -= b.x;
  a.y -= b.y;
  a.z -= b.z;
  return a;
}

Vector cross(Vector a, Vector b)
{
  Vector p;
  p.x = a.y * b.z - a.z * b.y;
  p.y = a.z * b.x - a.x * b.z;
  p.z = a.x * b.y - a.y * b.x;
  return p;
}

Vector sum(Vector a, Vector b)
{
  a.x += b.x;
  a.y += b.y;
  a.z += b.z;
  return a;
}

Vector sum3(Vector a, Vector b, Vector c)
{
  a.x += b.x + c.x;
  a.y += b.y + c.y;
  a.z += b.z + c.z;
  return a;
}

Vector rotate(Vector vertex, Vector axis, double angle)
{
  Vector p;
  p = scale(dot(axis, vertex), axis);
  return sum3(p, scale(cos(angle), diff(vertex, p)),
              scale(sin(angle), cross(axis, vertex)));
}

int same(Vector a, Vector b, double epsilon)
{
  return fabs(a.x - b.x) < epsilon && fabs(a.y - b.y) < epsilon &&
         fabs(a.z - b.z) < epsilon;
}

/*
 * Compute the polar reciprocal of the plane containing a, b and c:
 *
 * If this plane does not contain the origin, return p such that
 * dot(p,a) = dot(p,b) = dot(p,b) = r.
 *
 * Otherwise, return p such that
 * dot(p,a) = dot(p,b) = dot(p,c) = 0
 * and
 * dot(p,p) = 1.
 */
Vector pole(double r, Vector a, Vector b, Vector c)
{
  Vector p;
  double k;
  p = cross(diff(b, a), diff(c, a));
  k = dot(p, a);
  if (fabs(k) < 1e-6)
    return scale(1 / sqrt(dot(p, p)), p);
  else
    return scale(r / k, p);
}

/*
 * Find the numerator and the denominator using the Euclidean algorithm.
 */
void frac(double x)
{
  static Fraction zero = {0, 1}, inf = {1, 0};
  Fraction r;
  double s = x;
  r = zero;
  frax = inf;
  for (;;) {
    Fraction q;
    long f;
    // double floor();
    if (fabs(s) > (double)MAXLONG)
      return;
    f = (long)floor(s);
    q = r;
    r = frax;
    frax.n = frax.n * f + q.n;
    frax.d = frax.d * f + q.d;
    if (x == (double)frax.n / (double)frax.d)
      return;
    s = 1 / (s - f);
  }
}

// RK - add fp for file support
void printfrac(double x, FILE *fp)
{
  frac(x);
  fprintf(fp, "%ld", frax.n);
  if (frax.d != 1)
    fprintf(fp, "/%ld", frax.d);
}

char *sprintfrac(double x)
{
  char *s;
  frac(x);
  if (!frax.d) {
    Malloc(s, sizeof("infinity"), char);
    strcpy(s, "infinity");
  }
  else if (frax.d == 1) {
    char n[MAXDIGITS + 1];
    sprintf(n, "%ld", frax.n);
    Malloc(s, strlen(n) + 1, char);
    strcpy(s, n);
  }
  else {
    char n[MAXDIGITS + 1], d[MAXDIGITS + 1];
    sprintf(n, "%ld", frax.n);
    sprintf(d, "%ld", frax.d);
    Malloc(s, strlen(n) + strlen(d) + 2, char);
    sprintf(s, "%s/%s", n, d);
  }
  return s;
}

Polyhedron *polyalloc()
{
  Polyhedron *P;
  Calloc(P, 1, Polyhedron);
  P->index = -1;
  P->even = -1;
  P->K = 2;
  return P;
}

void *matalloc(int rows, int row_size)
{
  void **mat;
  int i = 0;
  if (!(mat = (void **)malloc(rows * sizeof(void *))))
    return 0;
  while ((mat[i] = malloc(row_size)) && ++i < rows)
    ;
  if (i == rows)
    return (void *)mat;
  while (--i >= 0)
    free(mat[i]);
  free(mat);
  return 0;
}

// RK - allow Kaleido, Maeder, Coxeter, or Wenninger indexes
// return #(meader number)
int process_index(char *sym)
{
  int ret = 1;
  char c;

  // remove leading spaces
  while ((c = *sym++) && isspace(c))
    ;

  // find index
  int idx = atoi(sym);

  // range is never checked. do it here.
  if (c == '#') {
    if (idx < 1 || idx > 80)
      ret = 0;
  }
  else if (c == 'k' || c == 'K') {
    // kaleido has prisms at the begining of the list
    if (idx < 1 || idx > 80)
      ret = 0;
    else {
      idx -= 5;
      if (idx < 1)
        idx += 80;
      sprintf(--sym, "#%d", idx);
    }
  }
  else if (c == 'u' || c == 'U') {
    // already in Maeder order
    if (idx < 1 || idx > 80)
      ret = 0;
    else
      sprintf(--sym, "#%d", idx);
  }
  else if (c == 'c' || c == 'C') {
    // Coxeter number
    Uniform *uniform;
    int last_uniform = get_uniform_list(&uniform);
    bool found = false;
    int i;
    for (i = 0; i < last_uniform; i++) {
      if (uniform[i].Coxeter == idx) {
        found = true;
        break;
      }
    }
    if (found)
      sprintf(--sym, "#%d", i + 1);
    else
      ret = 0;
  }
  else if (c == 'w' || c == 'W') {
    // Wenninger number
    Uniform *uniform;
    int last_uniform = get_uniform_list(&uniform);
    bool found = false;
    int i;
    for (i = 0; i < last_uniform; i++) {
      if (uniform[i].Wenninger == idx) {
        found = true;
        break;
      }
    }
    if (found)
      sprintf(--sym, "#%d", i + 1);
    else
      ret = 0;
  }
  // else return original string
  else {
    sym--;
  }

  return ret;
}

int unpacksym(char *sym, Polyhedron *P, Uniform *uniform, int last_uniform)
{
  int i = 0, n, d, bars = 0;
  char c;

  // RK - added functionality
  if (!process_index(sym))
    Err("index out of range");

  while ((c = *sym++) && isspace(c))
    ;
  if (!c)
    Err("no data");

  char wyth[MSG_SZ] = {0};
  if (c == '#') {
    while ((c = *sym++) && isspace(c))
      ;
    if (!c)
      Err("no digit after #");
    if (!isdigit(c))
      Err("not a digit");
    n = c - '0';
    while ((c = *sym++) && isdigit(c))
      n = n * 10 + c - '0';
    if (!n)
      Err("zero index");
    if (n > last_uniform)
      Err("index too big");
    sym--;
    while ((c = *sym++) && isspace(c))
      ;
    if (c)
      Err("data exceeded");
    strcpy_msg(wyth, uniform[P->index = n - 1].Wythoff);
    sym = wyth;
  }
  else
    sym--;
  for (;;) {
    while ((c = *sym++) && isspace(c))
      ;
    if (!c) {
      if (i == 4 && (bars || P->index == GR_DIRH_IDX))
        return 1;
      if (!bars)
        Err("no bars");
      Err("not enough fractions");
    }
    if (i == 4)
      Err("data exceeded");
    if (c == '|') {
      if (++bars > 1)
        Err("too many bars");
      P->p[i++] = 0;
      continue;
    }
    if (!isdigit(c))
      Err("not a digit");
    n = c - '0';
    while ((c = *sym++) && isdigit(c))
      n = n * 10 + c - '0';
    if (c && isspace(c))
      while ((c = *sym++) && isspace(c))
        ;
    if (c != '/') {
      sym--;
      if ((P->p[i++] = n) <= 1)
        Err("fraction<=1");
      continue;
    }
    while ((c = *sym++) && isspace(c))
      ;
    if (!c || !isdigit(c))
      return 0;
    d = c - '0';
    while ((c = *sym++) && isdigit(c))
      d = d * 10 + c - '0';
    if (!d)
      Err("zero denominator");
    sym--;
    if ((P->p[i++] = (double)n / d) <= 1)
      Err("fraction<=1");
  }
}

/*
 * Using Wythoff symbol (p|qr, pq|r, pqr| or |pqr), find the Moebius triangle
 * (2 3 K) (or (2 2 n)) of the Schwarz triangle (pqr), the order g of its
 * symmetry group, its Euler characteristic chi, and its covering density D.
 * g is the number of copies of (2 3 K) covering the sphere, i.e.,
 *
 *  g * pi * (1/2 + 1/3 + 1/K - 1) = 4 * pi
 *
 * D is the number of times g copies of (pqr) cover the sphere, i.e.
 *
 *  D * 4 * pi = g * pi * (1/p + 1/q + 1/r - 1)
 *
 * chi is V - E + F, where F = g is the number of triangles, E = 3*g/2 is the
 * number of triangle edges, and V = Vp+ Vq+ Vr, with Vp = g/(2*np) being the
 * number of vertices with angle pi/p (np is the numerator of p).
 */
int moebius(Polyhedron *P)
{
  int twos = 0, j, len = 1;
  /*
   * Arrange Wythoff symbol in a presentable form. In the same time check the
   * restrictions on the three fractions: They all have to be greater then one,
   * and the numerators 4 or 5 cannot occur together.  We count the ocurrences
   * of
   * 2 in `two', and save the largest numerator in `P->K', since they reflect on
   * the symmetry group.
   */
  P->K = 2;
  if (P->index == GR_DIRH_IDX) {
    Malloc(P->polyform, ++len, char);
    strcpy(P->polyform, "|");
  }
  else
    Calloc(P->polyform, len, char);
  for (j = 0; j < 4; j++) {
    if (P->p[j]) {
      char *s;
      Sprintfrac(s, P->p[j]);
      if (j && P->p[j - 1]) {
        Realloc(P->polyform, len += strlen(s) + 1, char);
        strcat(P->polyform, " ");
      }
      else
        Realloc(P->polyform, len += strlen(s), char);
      strcat(P->polyform, s);
      free(s);
      if (P->p[j] != 2) {
        int k;
        if ((k = numerator(P->p[j])) > P->K) {
          if (P->K == 4)
            break;
          P->K = k;
        }
        else if (k < P->K && k == 4)
          break;
      }
      else
        twos++;
    }
    else {
      Realloc(P->polyform, ++len, char);
      strcat(P->polyform, "|");
    }
  }
  /*
   * Find the symmetry group P->K (where 2, 3, 4, 5 represent the dihedral,
   * tetrahedral, octahedral and icosahedral groups, respectively), and its
   * order
   * P->g.
   */
  if (twos >= 2) { /* dihedral */
    P->g = 4 * P->K;
    P->K = 2;
  }
  else {
    if (P->K > 5)
      Err("numerator too large");
    P->g = 24 * P->K / (6 - P->K);
  }
  /*
   * Compute the nominal density P->D and Euler characteristic P->chi.
   * In few exceptional cases, these values will be modified later.
   */
  if (P->index != GR_DIRH_IDX) {
    int i;
    P->D = P->chi = -P->g;
    for (j = 0; j < 4; j++)
      if (P->p[j]) {
        P->chi += i = P->g / numerator(P->p[j]);
        P->D += i * denominator(P->p[j]);
      }
    P->chi /= 2;
    P->D /= 4;
    if (P->D <= 0)
      Err("nonpositive density");
  }
  return 1;
}

/*
 * Decompose Schwarz triangle into N right triangles and compute the vertex
 * count V and the vertex valency M.  V is computed from the number g of
 * Schwarz triangles in the cover, divided by the number of triangles which
 * share a vertex. It is halved for one-sided polyhedra, because the
 * kaleidoscopic construction really produces a double orientable covering of
 * such polyhedra. All q' q|r are of the "hemi" type, i.e. have equatorial {2r}
 * faces, and therefore are (except 3/2 3|3 and the dihedra 2 2|r) one-sided. A
 * well known example is 3/2 3|4, the "one-sided heptahedron". Also, all p q r|
 * with one even denominator have a crossed parallelogram as a vertex figure,
 * and thus are one-sided as well.
 */
int decompose(Polyhedron *P)
{
  int j, J, *s, *t;
  if (!P->p[1]) { /* p|q r */
    P->N = 2;
    P->M = 2 * numerator(P->p[0]);
    P->V = P->g / P->M;
    Malloc(P->n, P->N, double);
    Malloc(P->m, P->N, double);
    Malloc(P->rot, P->M, int);
    s = P->rot;
    for (j = 0; j < 2; j++) {
      P->n[j] = P->p[j + 2];
      P->m[j] = P->p[0];
    }
    for (j = P->M / 2; j--;) {
      *s++ = 0;
      *s++ = 1;
    }
  }
  else if (!P->p[2]) { /* p q|r */
    P->N = 3;
    P->M = 4;
    P->V = P->g / 2;
    Malloc(P->n, P->N, double);
    Malloc(P->m, P->N, double);
    Malloc(P->rot, P->M, int);
    s = P->rot;
    P->n[0] = 2 * P->p[3];
    P->m[0] = 2;
    for (j = 1; j < 3; j++) {
      P->n[j] = P->p[j - 1];
      P->m[j] = 1;
      *s++ = 0;
      *s++ = j;
    }
    if (fabs(P->p[0] - compli(P->p[1])) < DBL_EPSILON) { /* p = q' */
      /* P->p[0]==compli(P->p[1]) should work.  However, MSDOS
       * yeilds a 7e-17 difference! Reported by Jim Buddenhagen
       * <jb1556@daditz.sbc.com> */
      P->hemi = 1;
      P->D = 0;
      if (P->p[0] != 2 && !(P->p[3] == 3 && (P->p[0] == 3 || P->p[1] == 3))) {
        P->onesided = 1;
        P->V /= 2;
        P->chi /= 2;
      }
    }
  }
  else if (!P->p[3]) { /* p q r| */
    P->M = P->N = 3;
    P->V = P->g;
    Malloc(P->n, P->N, double);
    Malloc(P->m, P->N, double);
    Malloc(P->rot, P->M, int);
    s = P->rot;
    for (j = 0; j < 3; j++) {
      if (!(denominator(P->p[j]) % 2)) {
        /* what happens if there is more then one even denominator? */
        if (P->p[(j + 1) % 3] != P->p[(j + 2) % 3]) { /* needs postprocessing */
          P->even = j; /* memorize the removed face */
          P->chi -= P->g / numerator(P->p[j]) / 2;
          P->onesided = 1;
          P->D = 0;
        }
        else { /* for p = q we get a double 2 2r|p */
               /* noted by Roman Maeder <maeder@inf.ethz.ch> for 4 4 3/2| */
               /* Euler characteristic is still wrong */
          P->D /= 2;
        }
        P->V /= 2;
      }
      P->n[j] = 2 * P->p[j];
      P->m[j] = 1;
      *s++ = j;
    }
  }
  else { /* |p q r - snub polyhedron */
    P->N = 4;
    P->M = 6;
    P->V = P->g / 2; /* Only "white" triangles carry a vertex */
    Malloc(P->n, P->N, double);
    Malloc(P->m, P->N, double);
    Malloc(P->rot, P->M, int);
    Malloc(P->snub, P->M, int);
    s = P->rot;
    t = P->snub;
    P->m[0] = P->n[0] = 3;
    for (j = 1; j < 4; j++) {
      P->n[j] = P->p[j];
      P->m[j] = 1;
      *s++ = 0;
      *s++ = j;
      *t++ = 1;
      *t++ = 0;
    }
  }
  /*
   * Sort the fundamental triangles (using bubble sort) according to decreasing
   * n[i], while pushing the trivial triangles (n[i] = 2) to the end.
   */
  J = P->N - 1;
  while (J) {
    int last;
    last = J;
    J = 0;
    for (j = 0; j < last; j++) {
      if ((P->n[j] < P->n[j + 1] || P->n[j] == 2) && P->n[j + 1] != 2) {
        int i;
        double temp;
        temp = P->n[j];
        P->n[j] = P->n[j + 1];
        P->n[j + 1] = temp;
        temp = P->m[j];
        P->m[j] = P->m[j + 1];
        P->m[j + 1] = temp;
        for (i = 0; i < P->M; i++) {
          if (P->rot[i] == j)
            P->rot[i] = j + 1;
          else if (P->rot[i] == j + 1)
            P->rot[i] = j;
        }
        if (P->even != -1) {
          if (P->even == j)
            P->even = j + 1;
          else if (P->even == j + 1)
            P->even = j;
        }
        J = j;
      }
    }
  }
  /*
   *  Get rid of repeated triangles.
   */
  for (J = 0; J < P->N && P->n[J] != 2; J++) {
    int k, i;
    for (j = J + 1; j < P->N && P->n[j] == P->n[J]; j++)
      P->m[J] += P->m[j];
    k = j - J - 1;
    if (k) {
      for (i = j; i < P->N; i++) {
        P->n[i - k] = P->n[i];
        P->m[i - k] = P->m[i];
      }
      P->N -= k;
      for (i = 0; i < P->M; i++) {
        if (P->rot[i] >= j)
          P->rot[i] -= k;
        else if (P->rot[i] > J)
          P->rot[i] = J;
      }
      if (P->even >= j)
        P->even -= k;
    }
  }
  /*
   * Get rid of trivial triangles.
   */
  if (!J)
    J = 1; /* hosohedron */
  if (J < P->N) {
    int i;
    P->N = J;
    for (i = 0; i < P->M; i++) {
      if (P->rot[i] >= P->N) {
        for (j = i + 1; j < P->M; j++) {
          P->rot[j - 1] = P->rot[j];
          if (P->snub)
            P->snub[j - 1] = P->snub[j];
        }
        P->M--;
      }
    }
  }
  /*
   * Truncate arrays
   */
  Realloc(P->n, P->N, double);
  Realloc(P->m, P->N, double);
  Realloc(P->rot, P->M, int);
  if (P->snub)
    Realloc(P->snub, P->M, int) return 1;
}

/*
 * Solve the fundamental right spherical triangles.
 * If need_approx is set, print iterations on standard error.
 */
// RK: mods
// local function so we can have the -x option
// only write to *ofile from opts
int newton(Polyhedron *P, int need_approx, FILE *fp)
{
  /*
   * First, we find initial approximations.
   */
  int j;
  double cosa;
  Malloc(P->gamma, P->N, double);
  if (P->N == 1) {
    P->gamma[0] = M_PI / P->m[0];
    return 1;
  }
  for (j = 0; j < P->N; j++)
    P->gamma[j] = M_PI / 2 - M_PI / P->n[j];
  errno = 0; /* may be non-zero from some reason */
             /*
              * Next, iteratively find closer approximations for gamma[0] and compute
              * other gamma[j]'s from Napier's equations.
              */
  if (need_approx)
    fprintf(fp, "Solving %s\n", P->polyform);
  for (;;) {
    double delta = M_PI, sigma = 0;
    for (j = 0; j < P->N; j++) {
      if (need_approx)
        fprintf(fp, "%-20.15f", P->gamma[j]);
      delta -= P->m[j] * P->gamma[j];
    }
    if (need_approx)
      fprintf(fp, "(%g)\n", delta);
    if (fabs(delta) < 11 * DBL_EPSILON)
      return 1;
    /* On a RS/6000, fabs(delta)/DBL_EPSILON may occilate between 8 and
     * 10. Reported by David W. Sanderson <dws@ssec.wisc.edu> */
    for (j = 0; j < P->N; j++)
      sigma += P->m[j] * tan(P->gamma[j]);
    // sigma can be zero resulting in P->gamma[0] = infinity
    P->gamma[0] += delta * tan(P->gamma[0]) / sigma;
    // previously Err() had a return statement (return 0)
    if (P->gamma[0] < 0 || P->gamma[0] > M_PI)
      // Err("gamma out of bounds");
      return 0;
    cosa = cos(M_PI / P->n[0]) / sin(P->gamma[0]);
    for (j = 1; j < P->N; j++)
      P->gamma[j] = asin(cos(M_PI / P->n[j]) / cosa);
    // RK - output error message
    if (errno)
      // Err(0);
      printErrorMessage(errno, "newton");
  }
}

/*
 * Postprocess pqr| where r has an even denominator (cf. Coxeter &al. Sec.9).
 * Remove the {2r} and add a retrograde {2p} and retrograde {2q}.
 */
int exceptions(Polyhedron *P)
{
  int j;
  if (P->even != -1) {
    P->M = P->N = 4;
    Realloc(P->n, P->N, double);
    Realloc(P->m, P->N, double);
    Realloc(P->gamma, P->N, double);
    Realloc(P->rot, P->M, int);
    for (j = P->even + 1; j < 3; j++) {
      P->n[j - 1] = P->n[j];
      P->gamma[j - 1] = P->gamma[j];
    }
    P->n[2] = compli(P->n[1]);
    P->gamma[2] = -P->gamma[1];
    P->n[3] = compli(P->n[0]);
    P->m[3] = 1;
    P->gamma[3] = -P->gamma[0];
    P->rot[0] = 0;
    P->rot[1] = 1;
    P->rot[2] = 3;
    P->rot[3] = 2;
  }

  /*
   * Postprocess the last polyhedron |3/2 5/3 3 5/2 by taking a |5/3 3 5/2,
   * replacing the three snub triangles by four equatorial squares and adding
   * the missing {3/2} (retrograde triangle, cf. Coxeter &al. Sec. 11).
   */
  if (P->index == GR_DIRH_IDX) {
    P->N = 5;
    P->M = 8;
    Realloc(P->n, P->N, double);
    Realloc(P->m, P->N, double);
    Realloc(P->gamma, P->N, double);
    Realloc(P->rot, P->M, int);
    Realloc(P->snub, P->M, int);
    P->hemi = 1;
    P->D = 0;
    for (j = 3; j; j--) {
      P->m[j] = 1;
      P->n[j] = P->n[j - 1];
      P->gamma[j] = P->gamma[j - 1];
    }
    P->m[0] = P->n[0] = 4;
    P->gamma[0] = M_PI / 2;
    P->m[4] = 1;
    P->n[4] = compli(P->n[1]);
    P->gamma[4] = -P->gamma[1];
    for (j = 1; j < 6; j += 2)
      P->rot[j]++;
    P->rot[6] = 0;
    P->rot[7] = 4;
    P->snub[6] = 1;
    P->snub[7] = 0;
  }
  return 1;
}

/*
 * Compute edge and face counts, and update D and chi.  Update D in the few
 * cases the density of the polyhedron is meaningful but different than the
 * density of the corresponding Schwarz triangle (cf. Coxeter &al., p. 418 and
 * p. 425).
 * In these cases, spherical faces of one type are concave (bigger than a
 * hemisphere), and the actual density is the number of these faces less the
 * computed density.  Note that if j != 0, the assignment gamma[j] = asin(...)
 * implies gamma[j] cannot be obtuse.  Also, compute chi for the only
 * non-Wythoffian polyhedron.
 */
int count(Polyhedron *P)
{
  int j, temp;
  Malloc(P->Fi, P->N, int);
  for (j = 0; j < P->N; j++) {
    P->E += temp = P->V * numerator(P->m[j]);
    P->F += P->Fi[j] = temp / numerator(P->n[j]);
  }
  P->E /= 2;
  if (P->D && P->gamma[0] > M_PI / 2)
    P->D = P->Fi[0] - P->D;
  if (P->index == GR_DIRH_IDX)
    P->chi = P->V - P->E + P->F;
  return 1;
}

/*
 * Generate a printable vertex configuration symbol.
 */
int configuration(Polyhedron *P)
{
  int j, len = 2;
  for (j = 0; j < P->M; j++) {
    char *s;
    Sprintfrac(s, P->n[P->rot[j]]);
    len += strlen(s) + 1;
    if (!j) {
      Malloc(P->config, len, char);
      strcpy(P->config, "(");
    }
    else {
      Realloc(P->config, len, char);
      strcat(P->config, ".");
    }
    strcat(P->config, s);
    free(s);
  }
  strcat(P->config, ")");
  if ((j = denominator(P->m[0])) != 1) {
    char s[MAXDIGITS + 2];
    sprintf(s, "/%d", j);
    Realloc(P->config, len + strlen(s), char);
    strcat(P->config, s);
  }
  return 1;
}

/*
 * Compute polyhedron vertices and vertex adjecency lists.
 * The vertices adjacent to v[i] are v[adj[0][i], v[adj[1][i], ...
 * v[adj[M-1][i], ordered counterclockwise.  The algorith is a BFS on the
 * vertices, in such a way that the vetices adjacent to a givem vertex are
 * obtained from its BFS parent by a cyclic sequence of rotations. firstrot[i]
 * points to the first  rotaion in the sequence when applied to v[i]. Note that
 * for non-snub polyhedra, the rotations at a child are opposite in sense when
 * compared to the rotations at the parent. Thus, we fill adj[*][i] from the
 * end to signify clockwise rotations. The firstrot[] array is not needed for
 * display thus it is freed after being used for face computations below.
 */
int vertices(Polyhedron *P)
{
  int i, newV = 2;
  double cosa;
  Malloc(P->v, P->V, Vector);
  Matalloc(P->adj, P->M, P->V, int);
  Malloc(P->firstrot, P->V, int); /* temporary , put in Polyhedron
        structure so that may be freed on error */
  cosa = cos(M_PI / P->n[0]) / sin(P->gamma[0]);
  P->v[0].x = 0;
  P->v[0].y = 0;
  P->v[0].z = 1;
  P->firstrot[0] = 0;
  P->adj[0][0] = 1;
  P->v[1].x = 2 * cosa * sqrt(1 - cosa * cosa);
  P->v[1].y = 0;
  P->v[1].z = 2 * cosa * cosa - 1;
  if (!P->snub) {
    P->firstrot[1] = 0;
    P->adj[0][1] = -1; /* start the other side */
    P->adj[P->M - 1][1] = 0;
  }
  else {
    P->firstrot[1] = P->snub[P->M - 1] ? 0 : P->M - 1;
    P->adj[0][1] = 0;
  }
  for (i = 0; i < newV; i++) {
    int j, k;
    int last, one, start, limit;
    if (P->adj[0][i] == -1) {
      one = -1;
      start = P->M - 2;
      limit = -1;
    }
    else {
      one = 1;
      start = 1;
      limit = P->M;
    }
    k = P->firstrot[i];
    for (j = start; j != limit; j += one) {
      Vector temp;
      int J;
      temp = rotate(P->v[P->adj[j - one][i]], P->v[i],
                    one * 2 * P->gamma[P->rot[k]]);
      for (J = 0; J < newV && !same(P->v[J], temp, BIG_EPSILON); J++)
        ; /* noop */
      P->adj[j][i] = J;
      last = k;
      if (++k == P->M)
        k = 0;
      if (J == newV) { /* new vertex */
        if (newV == P->V)
          Err("too many vertices");
        P->v[newV++] = temp;
        if (!P->snub) {
          P->firstrot[J] = k;
          if (one > 0) {
            P->adj[0][J] = -1;
            P->adj[P->M - 1][J] = i;
          }
          else {
            P->adj[0][J] = i;
          }
        }
        else {
          P->firstrot[J] =
              !P->snub[last] ? last : !P->snub[k] ? (k + 1) % P->M : k;
          P->adj[0][J] = i;
        }
      }
    }
  }
  return 1;
}

/*
 * Compute polyhedron faces (dual vertices) and incidence matrices.
 * For orientable polyhedra, we can distinguish between the two faces meeting
 * at a given directed edge and identify the face on the left and the face on
 * the right, as seen from the outside.  For one-sided polyhedra, the vertex
 * figure is a papillon (in Coxeter &al.  terminology, a crossed parallelogram)
 * and the two faces meeting at an edge can be identified as the side face
 * (n[1] or n[2]) and the diagonal face (n[0] or n[3]).
 */
int faces(Polyhedron *P)
{
  int i, newF = 0;
  Malloc(P->f, P->F, Vector);
  Malloc(P->ftype, P->F, int);
  Matalloc(P->incid, P->M, P->V, int);
  P->minr = 1 / fabs(tan(M_PI / P->n[P->hemi]) * tan(P->gamma[P->hemi]));
  for (i = P->M; --i >= 0;) {
    int j;
    for (j = P->V; --j >= 0;)
      P->incid[i][j] = -1;
  }
  for (i = 0; i < P->V; i++) {
    int j;
    for (j = 0; j < P->M; j++) {
      int i0, J;
      int pap = 0; /* papillon edge type */
      if (P->incid[j][i] != -1)
        continue;
      P->incid[j][i] = newF;
      if (newF == P->F)
        Err("too many faces");
      P->f[newF] = pole(P->minr, P->v[i], P->v[P->adj[j][i]],
                        P->v[P->adj[mod(j + 1, P->M)][i]]);
      P->ftype[newF] = P->rot[mod(
          P->firstrot[i] + (P->adj[0][i] < P->adj[P->M - 1][i] ? j : -j - 2),
          P->M)];
      if (P->onesided)
        pap = (P->firstrot[i] + j) % 2;
      i0 = i;
      J = j;
      for (;;) {
        int k;
        k = i0;
        if ((i0 = P->adj[J][k]) == i)
          break;
        for (J = 0; J < P->M && P->adj[J][i0] != k; J++)
          ; /* noop */
        if (J == P->M)
          Err("too many faces");
        if (P->onesided && (J + P->firstrot[i0]) % 2 == pap) {
          P->incid[J][i0] = newF;
          if (++J >= P->M)
            J = 0;
        }
        else {
          if (--J < 0)
            J = P->M - 1;
          P->incid[J][i0] = newF;
        }
      }
      newF++;
    }
  }
  Free(P->firstrot) Free(P->rot) Free(P->snub);
  return 1;
}

Polyhedron *kaleido(char *sym, Uniform *uniform, int last_uniform)
{
  Polyhedron *P;
  /*
   * Allocate a Polyhedron structure P.
   */
  if (!(P = polyalloc()))
    return 0;
  /*
   * Unpack input symbol into P.
   */
  if (!unpacksym(sym, P, uniform, last_uniform))
    return 0;
  /*
   * Find Mebius triangle, its density and Euler characteristic.
   */
  if (!moebius(P))
    return 0;
  /*
   * Decompose Schwarz triangle.
   */
  if (!decompose(P))
    return 0;
  /*
   * Solve Fundamental triangles, optionally printing approximations.
   */
  if (!newton(P, 0, NULL))
    return 0;
  /*
   * Deal with exceptional polyhedra.
   */
  if (!exceptions(P))
    return 0;
  /*
   * Count edges and faces, update density and characteristic if needed.
   */
  if (!count(P))
    return 0;
  /*
   * Generate printable vertex configuration.
   */
  if (!configuration(P))
    return 0;
  /*
   * Compute coordinates.
   */
  if (!vertices(P))
    return 0;
  if (!faces(P))
    return 0;

  return P;
}

// RK: Ported Functions for Printing and VRML of Kaleido

/*
 * Print polyhedron data
 */
// RK: mods
// uniform list is referenced (must have been global) so pass it from main
// 'static char *' becomes 'static const char *' (deprication)
// only write to *ofile from opts
int printit(Polyhedron *P, int need_coordinates, int just_list, int digits,
            int more_lines, Uniform *uniform_list, FILE *fp)
{
  int j, i;
  double cosa;
  static const char *group[] = {"di", "tetra", "octa", "icosa"};
  static const char *alias[] = {"D", "A4", "S4", "A5"};
  /*
   * Print polyhedron name, Wythoff symbol, and reference figures.
   */
  if (P->index != -1)
    fprintf(fp, "%d) ", P->index + 1);
  fprintf(fp, "%s %s", P->name, P->polyform);
  if (P->index != -1 && uniform_list[P->index].Coxeter)
    fprintf(fp, " [%d,%d]", uniform_list[P->index].Coxeter,
            uniform_list[P->index].Wenninger);
  fprintf(fp, "\n");
  if (just_list)
    return 1;
  /*
   * Print combinatorial description.
   */
  fprintf(fp, "\tdual: %s\n\t%s, %shedral group %s", P->dual_name, P->config,
          group[P->K - 2], alias[P->K - 2]);
  if (P->K == 2)
    fprintf(fp, "%d", P->g / 4);
  fprintf(fp, ", chi=%d", P->chi);
  if (P->D)
    fprintf(fp, ", D=%d", P->D);
  else if (P->onesided)
    fprintf(fp, ", one-sided");
  else
    fprintf(fp, ", two-sided");
  fprintf(fp, "\n\tV=%d, E=%d, F=%d=", P->V, P->E, P->F);
  for (j = 0; j < P->N; j++) {
    if (j)
      fprintf(fp, "+");
    fprintf(fp, "%d{", P->Fi[j]);
    printfrac(P->n[j], fp);
    fprintf(fp, "}");
  }
  /*
   * Print solution.
   */
  if (P->index == -1 && P->K == 2) {
    char *s;
    Sprintfrac(s, P->gon);
    i = strlen(s) + 2;
    if (i < 6)
      i = 6;
    free(s);
  }
  else
    i = 6;
  fprintf(fp, "\n%*s%6s%*s%*s%*s%*s%*s%*s%*s%*s\n", i, "", "alpha", digits + 3,
          "gamma", digits + 1, "a", digits + 1, "b", digits + 1, "c",
          digits + 3, "rho/R", digits + 3, "r/rho", digits + 3, "l/rho",
          digits + 3, "h/r");
  cosa = cos(M_PI / P->n[0]) / sin(P->gamma[0]);
  for (j = 0; j < P->N; j++) {
    double cosc = cos(P->gamma[j]) / sin(M_PI / P->n[j]);
    char *s, *t;
    Sprintfrac(s, P->n[j]);
    Malloc(t, strlen(s) + 3, char);
    sprintf(t, "{%s}", s);
    free(s);
    fprintf(fp, "%*s%6.1f%*.*f%*.*f%*.*f%*.*f%*.*f%*.*f", i, t, 180. / P->n[j],
            digits + 3, digits - 2, DEG * P->gamma[j], digits + 1, digits - 4,
            DEG * acos(cosa), digits + 1, digits - 4, DEG * acos(cosa * cosc),
            digits + 1, digits - 4, DEG * acos(cosc), digits + 3, digits, cosa,
            digits + 3, digits, cosc);
    if (log10(fabs(cosa)) < -digits)
      fprintf(fp, "%*s", digits + 3, "infinity");
    else
      fprintf(fp, "%*.*f", digits + 3, digits, sqrt(1 - cosa * cosa) / cosa);
    if (log10(fabs(cosc)) < -digits)
      fprintf(fp, "%*s", digits + 3, "infinity");
    else
      fprintf(fp, "%*.*f", digits + 3, digits, sqrt(1 - cosc * cosc) / cosc);
    fprintf(fp, "\n");
    free(t);
  }
  fprintf(fp, "\n");
  more();
  if (!need_coordinates)
    return 1;
  /*
   * Print vertices
   */
  fprintf(fp, "vertices:\n");
  for (i = 0; i < P->V; i++) {
    fprintf(fp, "v%-3d (%*.*f,%*.*f,%*.*f)", i + 1, digits + 3, digits,
            P->v[i].x, digits + 3, digits, P->v[i].y, digits + 3, digits,
            P->v[i].z);
    for (j = 0; j < P->M; j++)
      fprintf(fp, " v%-3d", P->adj[j][i] + 1);
    fprintf(fp, "\n%*s", 3 * digits + 20, "");
    for (j = 0; j < P->M; j++)
      fprintf(fp, " f%-3d", P->incid[j][i] + 1);
    fprintf(fp, "\n");
    if (!((i + 1) % (more_lines / 2)))
      more();
  }
  if (P->V % (more_lines / 2))
    more();
  /*
   * Print faces.
   */
  fprintf(fp, "faces (RHS=%*.*f", digits + 2, digits, P->minr);
  fprintf(fp, "):\n");
  for (i = 0; i < P->F; i++) {
    fprintf(fp, "f%-3d (%*.*f,%*.*f,%*.*f) {", i + 1, digits + 3, digits,
            P->f[i].x, digits + 3, digits, P->f[i].y, digits + 3, digits,
            P->f[i].z);
    printfrac(P->n[P->ftype[i]], fp);
    fprintf(fp, "}");
    if (P->hemi && !P->ftype[i])
      fprintf(fp, "*");
    fprintf(fp, "\n");
    if (!((i + 1) % more_lines))
      more();
  }
  if (P->F % more_lines)
    more();
  fprintf(fp, "\n");
  return 1;
}

/*
 * rotate the standard frame
 */
void rotframe(double azimuth, double elevation, double angle)
{
  static Vector X = {1, 0, 0}, Y = {0, 1, 0}, Z = {0, 0, 1};
  Vector axis;

  axis = rotate(rotate(X, Y, elevation), Z, azimuth);
  x = rotate(X, axis, angle);
  y = rotate(Y, axis, angle);
  z = rotate(Z, axis, angle);
}

/*
 * rotate an array of n Vectors
 */
void rotarray(Vector *newvec, Vector *old, int n)
{
  while (n--) {
    *newvec++ = sum3(scale(old->x, x), scale(old->y, y), scale(old->z, z));
    old++;
  }
}

void printvec(FILE *fp, Vector v, int digits)
{
  fprintf(fp, "\t\t\t\t%*.*f %*.*f %*.*f,\n", digits + 3, digits, v.x,
          digits + 3, digits, v.y, digits + 3, digits, v.z);
}

/*
 * Choose a color for a given face valency.
 * 3-,4- and 5-sided polygons have the simple colors R, G, and B.
 * 6-,8- and 10-sided polygons, which also occur in the standard polyhedra,
 * have the blended colors RG (yellow), RB (magenta) and GB (cyan).
 * All othe polygons (which occur in kaleido only if a Whythoff formula is
 * entered) are colored pink (its RGB value is taken from X11's rgb.txt).
 */
void rgbcolor(FILE *fp, int n)
{
  double R, G, B;
  switch (n) {
  case 3: /* red */
    R = 1;
    G = 0;
    B = 0;
    break;
  case 4: /* green */
    R = 0;
    G = 1;
    B = 0;
    break;
  case 5: /* blue */
    R = 0;
    G = 0;
    B = 1;
    break;
  case 6: /* yellow */
    R = 1;
    G = 1;
    B = 0;
    break;
  case 8: /* magenta */
    R = 1;
    G = 0;
    B = 1;
    break;
  case 10: /* cyan */
    R = 0;
    G = 1;
    B = 1;
    break;
  default: /* pink */
    R = 1;
    G = 192. / 255.;
    B = 203. / 255.;
    break;
  }
  fprintf(fp, "%g %g %g,", R, G, B);
}

// RK: mods
// message lines 0 and 2 put in as strings
// k, *hit not initialized
// char *star add const per compiler warning
// prefix no longer used, removed
// pass last_uniform from main (added)
// only write to *ofile from opts
// fixed if statements per clang warning
int vrmodel(Polyhedron *P, Vector *v, int V, Vector *f, int F, char *name,
            const char *star, int digits, double azimuth, double elevation,
            double freeze, FILE *fp)
{
  int i, j, l, ll, ii, facelets;
  int *hit = 0;
  int k = 0;
  Vector *temp;

  // RK - Last uniform has changed to Maeder index
  int last_uniform = 75;

  // RK: last remaining part of original usage message
  string kaleido_message0 =
      "Kaleidoscopic Construction of Uniform Polyhedra, $Revision: 3.27 $";
  string kaleido_message2 =
      "Copyright \302\251 1991-2002 Dr. Zvi Har'El <rl@math.technion.ac.il>";

  /*
     FILE *fp;
     char *fn;
     Malloc(fn, strlen(prefix) + 8, char);
     for (i = 1; i < 1000; i++) {
        sprintf(fn, "%s%03d.wrl", prefix, i);
        if (access(fn, 0)) {       // file doesn't exist
           fp = fopen(fn, "w");
           break;
        }
     }
     if (i == 1000)
        Err("too many files");
     if (!fp)
        Err(0);
  */

  /*
   * Rotate polyhedron
   */
  rotframe(azimuth, elevation, freeze);

  Malloc(temp, V, Vector) rotarray(temp, v, V);
  v = temp;
  Malloc(temp, F, Vector) rotarray(temp, f, F);
  f = temp;
  /*
   * File header
   */
  fprintf(fp, "#VRML V2.0 utf8\n");
  fprintf(fp, "WorldInfo{\n");
  fprintf(fp, "\tinfo[\n");
  fprintf(fp, "\t\t\"%s\"\n", kaleido_message0.c_str());
  fprintf(fp, "\t\t\"%s\"\n", kaleido_message2.c_str());
  fprintf(fp, "\t]\n");
  fprintf(fp, "\ttitle \"%d%s) %s %s %s\"\n", P->index + 1, star, name,
          P->polyform, P->config);
  fprintf(fp, "}\n");
  fprintf(fp, "NavigationInfo {\n");
  fprintf(fp, "\ttype \"EXAMINE\"\n");
  /* fprintf(fp, "\theadlight TRUE\n"); */
  fprintf(fp, "}\n");
  fprintf(fp, "Shape{\n");
  fprintf(fp, "\tappearance Appearance{\n");
  fprintf(fp, "\t\tmaterial Material{\n");
  fprintf(fp, "\t\t\tshininess 1\n");
  fprintf(fp, "\t\t}\n");
  fprintf(fp, "\t}\n");
  fprintf(fp, "\tgeometry IndexedFaceSet{\n");
  if (P->D != 1) {
    fprintf(fp, "\t\tconvex FALSE\n");
    fprintf(fp, "\t\tsolid FALSE\n");
  }
  fprintf(fp, "\t\tcreaseAngle 0\n");
  fprintf(fp, "\t\tcolorPerVertex FALSE\n");
  /*
   * Color map
   * Face colors are assigned as a function of the face valency.
   * Thus, pentagons and pentagrams will be colored alike.
   */
  fprintf(fp, "\t\tcolor Color{\n");
  fprintf(fp, "\t\t\tcolor[");
  if (*star)
    rgbcolor(fp, P->M);
  else
    for (i = 0; i < P->N; i++)
      rgbcolor(fp, numerator(P->n[i]));
  fprintf(fp, "]\n");
  fprintf(fp, "\t\t}\n");
  /*
   * Vertex list
   */
  fprintf(fp, "\t\tcoord Coordinate{\n");
  fprintf(fp, "\t\t\tpoint[\n");
  for (i = 0; i < V; i++)
    printvec(fp, v[i], digits);
  /*
   * Auxiliary vertices (needed because current VRML browsers cannot handle
   * non-simple polygons, i.e., ploygons with self intersections):
   * Each non-simple face is assigned an auxiliary vertex. By connecting it to
   * the rest of the vertices the face is triangulated. The circum-center is
   * used
   * for the regular star faces of uniform polyhedra. The in-center is used for
   * the pentagram (#79) and hexagram (#77) of the high-density snub duals, and
   * for the pentagrams (#40, #58) and hexagram (#52) of the stellated duals
   * with
   * configuration (....)/2. Finally, the self-intersection of the crossed
   * parallelogram is used for duals with form p q r| with an even denominator.
   *
   * This method do not work for the hemi-duals, whose faces are not star-shaped
   * and have two self-intersections each.
   * Thus, for each face we need six auxiliary vertices: The self intersections
   * and the terminal points of the truncations of the infinite edges. The ideal
   * vertices are listed, but are not used by the face-list.
   * Note that the face of the last dual (#80) is octagonal, and constists of
   * two
   * quadrilaterals of the infinite type.
   */
  fprintf(fp, "\t\t\t\t# auxiliary vertices:\n");
  if (*star && P->even != -1)
    Malloc(hit, F, int);
  for (i = 0; i < F; i++)
    if ((!*star &&
         (frac(P->n[P->ftype[i]]), frax.d != 1 && frax.d != frax.n - 1)) ||
        (*star && ((P->K == 5 && P->D > 30) || denominator(P->m[0]) != 1))) {
      /* find the center of the face */
      double h;
      if (!*star && P->hemi && !P->ftype[i])
        h = 0;
      else
        h = P->minr / dot(f[i], f[i]);
      printvec(fp, scale(h, f[i]), digits);
    }
    else if (*star && P->even != -1) {
      /* find the self-intersection of a crossed parallelogram.
       * hit is set if v0v1 intersects v2v3*/
      Vector v0, v1, v2, v3, c0, c1, p;
      double d0, d1;
      v0 = v[P->incid[0][i]];
      v1 = v[P->incid[1][i]];
      v2 = v[P->incid[2][i]];
      v3 = v[P->incid[3][i]];
      d0 = sqrt(dot(diff(v0, v2), diff(v0, v2)));
      d1 = sqrt(dot(diff(v1, v3), diff(v1, v3)));
      c0 = scale(d1, sum(v0, v2));
      c1 = scale(d0, sum(v1, v3));
      p = scale(0.5 / (d0 + d1), sum(c0, c1));
      printvec(fp, p, digits);
      p = cross(diff(p, v2), diff(p, v3));
      hit[i] = (dot(p, p) < 1e-6);
    }
    else if (*star && P->hemi && P->index != last_uniform - 1) {
      /* find the terminal points of the truncation and the
       * self-intersections.
       *  v23       v0       v21
       *  |  \     /  \     /  |
       *  |   v0123    v0321   |
       *  |  /     \  /     \  |
       *  v01       v2       v03
       */
      Vector v0, v1, v2, v3, v01, v03, v21, v23, v0123, v0321;
      Vector u;
      double t = 1.5; /* truncation adjustment factor */
      j = !P->ftype[P->incid[0][i]];
      v0 = v[P->incid[j][i]];     /* real vertex */
      v1 = v[P->incid[j + 1][i]]; /* ideal vertex (unit vector) */
      v2 = v[P->incid[j + 2][i]]; /* real */
                                  /* ideal */
      v3 = v[P->incid[(j + 3) % 4][i]];
      /* compute intersections
       * this uses the following linear algebra:
       * v0123 = v0 + a v1 = v2 + b v3
       * v0 x v3 + a (v1 x v3) = v2 x v3
       * a (v1 x v3) = (v2 - v0) x v3
       * a (v1 x v3) . (v1 x v3) = (v2 - v0) x v3 . (v1 x v3)
       */
      u = cross(v1, v3);
      v0123 = sum(v0, scale(dot(cross(diff(v2, v0), v3), u) / dot(u, u), v1));
      v0321 = sum(v0, scale(dot(cross(diff(v0, v2), v1), u) / dot(u, u), v3));
      /* compute truncations */
      v01 = sum(v0, scale(t, diff(v0123, v0)));
      v23 = sum(v2, scale(t, diff(v0123, v2)));
      v03 = sum(v0, scale(t, diff(v0321, v0)));
      v21 = sum(v2, scale(t, diff(v0321, v2)));
      printvec(fp, v01, digits);
      printvec(fp, v23, digits);
      printvec(fp, v0123, digits);
      printvec(fp, v03, digits);
      printvec(fp, v21, digits);
      printvec(fp, v0321, digits);
    }
    else if (*star && P->index == last_uniform - 1) {
      /* find the terminal points of the truncation and the
       * self-intersections.
       *  v23       v0       v21
       *  |  \     /  \     /  |
       *  |   v0123    v0721   |
       *  |  /     \  /     \  |
       *  v01       v2       v07
       *
       *  v65       v4       v67
       *  |  \     /  \     /  |
       *  |   v4365    v4567   |
       *  |  /     \  /     \  |
       *  v43       v6       v45
       */
      Vector v0, v1, v2, v3, v4, v5, v6, v7, v01, v07, v21, v23;
      Vector v43, v45, v65, v67, v0123, v0721, v4365, v4567;
      double t = 1.5; /* truncation adjustment factor */
      Vector u;
      for (j = 0; j < 8; j++)
        if (P->ftype[P->incid[j][i]] == 3)
          break;
      v0 = v[P->incid[j][i]]; /* real {5/3} */
                              /* ideal */
      v1 = v[P->incid[(j + 1) % 8][i]];
      /* real {3} */
      v2 = v[P->incid[(j + 2) % 8][i]];
      /* ideal */
      v3 = v[P->incid[(j + 3) % 8][i]];
      /* real {5/2} */
      v4 = v[P->incid[(j + 4) % 8][i]];
      /* ideal */
      v5 = v[P->incid[(j + 5) % 8][i]];
      /* real {3/2} */
      v6 = v[P->incid[(j + 6) % 8][i]];
      /* ideal */
      v7 = v[P->incid[(j + 7) % 8][i]];
      /* compute intersections */
      u = cross(v1, v3);
      v0123 = sum(v0, scale(dot(cross(diff(v2, v0), v3), u) / dot(u, u), v1));
      u = cross(v7, v1);
      v0721 = sum(v0, scale(dot(cross(diff(v2, v0), v1), u) / dot(u, u), v7));
      u = cross(v5, v7);
      v4567 = sum(v4, scale(dot(cross(diff(v6, v4), v7), u) / dot(u, u), v5));
      u = cross(v3, v5);
      v4365 = sum(v4, scale(dot(cross(diff(v6, v4), v5), u) / dot(u, u), v3));
      /* compute truncations */
      v01 = sum(v0, scale(t, diff(v0123, v0)));
      v23 = sum(v2, scale(t, diff(v0123, v2)));
      v07 = sum(v0, scale(t, diff(v0721, v0)));
      v21 = sum(v2, scale(t, diff(v0721, v2)));
      v45 = sum(v4, scale(t, diff(v4567, v4)));
      v67 = sum(v6, scale(t, diff(v4567, v6)));
      v43 = sum(v4, scale(t, diff(v4365, v4)));
      v65 = sum(v6, scale(t, diff(v4365, v6)));
      printvec(fp, v01, digits);
      printvec(fp, v23, digits);
      printvec(fp, v0123, digits);
      printvec(fp, v07, digits);
      printvec(fp, v21, digits);
      printvec(fp, v0721, digits);
      printvec(fp, v45, digits);
      printvec(fp, v67, digits);
      printvec(fp, v4567, digits);
      printvec(fp, v43, digits);
      printvec(fp, v65, digits);
      printvec(fp, v4365, digits);
    }
  fprintf(fp, "\t\t\t]\n");
  fprintf(fp, "\t\t}\n");
  /*
   * Face list:
   * Each face is printed in a separate line, by listing the indices of its
   * vertices. In the non-simple case, the polygon is represented by the
   * triangulation, each triangle consists of two polyhedron vertices and one
   * auxiliary vertex.
   */
  fprintf(fp, "\t\tcoordIndex[\n");
  ii = V;
  facelets = 0;
  for (i = 0; i < F; i++) {
    fprintf(fp, "\t\t\t");
    if (*star) {
      if ((P->K == 5 && P->D > 30) || denominator(P->m[0]) != 1) {
        for (j = 0; j < P->M - 1; j++) {
          fprintf(fp, "%d,%d,%d,-1,", P->incid[j][i], P->incid[j + 1][i], ii);
          facelets++;
        }
        fprintf(fp, "%d,%d,%d,-1,", P->incid[j][i], P->incid[0][i], ii++);
        facelets++;
      }
      else if (P->even != -1) {
        if (hit[i]) {
          fprintf(fp, "%d,%d,%d,-1,%d,%d,%d,-1,", P->incid[3][i],
                  P->incid[0][i], ii, P->incid[1][i], P->incid[2][i], ii);
        }
        else {
          fprintf(fp, "%d,%d,%d,-1,%d,%d,%d,-1,", P->incid[0][i],
                  P->incid[1][i], ii, P->incid[2][i], P->incid[3][i], ii);
        }
        ii++;
        facelets += 2;
      }
      else if (P->hemi && P->index != last_uniform - 1) {
        j = !P->ftype[P->incid[0][i]];
        fprintf(fp, "%d,%d,%d,-1,%d,%d,%d,%d,-1,%d,%d,%d,-1,", ii, ii + 1,
                ii + 2, P->incid[j][i], ii + 2, P->incid[j + 2][i], ii + 5,
                ii + 3, ii + 4, ii + 5);
        ii += 6;
        facelets += 3;
      }
      else if (P->index == last_uniform - 1) {
        for (j = 0; j < 8; j++)
          if (P->ftype[P->incid[j][i]] == 3)
            break;
        fprintf(fp, "%d,%d,%d,-1,%d,%d,%d,%d,-1,%d,%d,%d,-1,", ii, ii + 1,
                ii + 2, P->incid[j][i], ii + 2, P->incid[(j + 2) % 8][i],
                ii + 5, ii + 3, ii + 4, ii + 5);
        fprintf(fp, "%d,%d,%d,-1,%d,%d,%d,%d,-1,%d,%d,%d,-1,", ii + 6, ii + 7,
                ii + 8, P->incid[(j + 4) % 8][i], ii + 8,
                P->incid[(j + 6) % 8][i], ii + 11, ii + 9, ii + 10, ii + 11);
        ii += 12;
        facelets += 6;
      }
      else {
        for (j = 0; j < P->M; j++)
          fprintf(fp, "%d,", P->incid[j][i]);
        fprintf(fp, "-1,");
        facelets++;
      }
    }
    else {
      int split =
          (frac(P->n[P->ftype[i]]), frax.d != 1 && frax.d != frax.n - 1);
      for (j = 0; j < V; j++) {
        for (k = 0; k < P->M; k++)
          if (P->incid[k][j] == i)
            break;
        if (k != P->M)
          break;
      }
      if (!split)
        fprintf(fp, "%d,", j);
      ll = j;
      for (l = P->adj[k][j]; l != j; l = P->adj[k][l]) {
        for (k = 0; k < P->M; k++)
          if (P->incid[k][l] == i)
            break;
        if (P->adj[k][l] == ll)
          k = mod(k + 1, P->M);
        if (!split)
          fprintf(fp, "%d,", l);
        else {
          fprintf(fp, "%d,%d,%d,-1,", ll, l, ii);
          facelets++;
        }
        ll = l;
      }
      if (!split) {
        fprintf(fp, "-1,");
        facelets++;
      }
      else {
        fprintf(fp, "%d,%d,%d,-1,", ll, j, ii++);
        facelets++;
      }
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "\t\t]\n");
  /*
   * Face color indices - for polyhedra with multiple face types
   * For non-simple faces, the index is repeated as many times as needed by the
   * triangulation.
   */
  fprintf(fp, "\t\tcolorIndex[");
  if (!*star && P->N != 1) {
    for (i = 0; i < F; i++)
      if (frac(P->n[P->ftype[i]]), frax.d == 1 || frax.d == frax.n - 1)
        fprintf(fp, "%d,", P->ftype[i]);
      else
        for (j = 0; j < frax.n; j++)
          fprintf(fp, "%d,", P->ftype[i]);
  }
  else {
    for (i = 0; i < facelets; i++)
      fprintf(fp, "0,");
  }
  fprintf(fp, "]\n");
  fprintf(fp, "\t}\n");
  if (*star && P->even != -1)
    free(hit);
  free(v);
  free(f);
  fprintf(fp, "}\n");
  if (ferror(fp))
    Err(0);
  fclose(fp);
  // fprintf(stderr, "%s: written on %s\n", name, fn);
  if (P->index != -1)
    fprintf(stderr, "%s\n", name);
  return 1;
}

// end original Kaleido code

// RK - from unipoly ported form antiprism 19.1 by Adrian Rossiter
int get_poly(Geometry &geom, Polyhedron *poly)
{
  vector<Vec3d> &verts = geom.raw_verts();
  // vector<vector<int> > &edges = raw_edges();
  vector<vector<int>> &faces = geom.raw_faces();

  faces.resize(poly->F);
  vector<int> edges;
  for (int i = 0; i < poly->V; i++) {
    verts.push_back(Vec3d(poly->v[i].x, poly->v[i].y, poly->v[i].z));

    for (int j = 0; j < poly->M; j++) {
      faces[poly->incid[j][i]].push_back(i);
      if (i < poly->adj[j][i]) {
        edges.push_back(i);
        edges.push_back(poly->adj[j][i]);
      }
    }
  }

  for (unsigned int i = 0; i < faces.size(); i++) {
    set<int> vs;
    for (unsigned int j = 0; j < faces[i].size(); j++)
      vs.insert(faces[i][j]);

    vector<int> lns;
    for (unsigned int j = 0; j < edges.size(); j += 2)
      if (vs.find(edges[j]) != vs.end() && vs.find(edges[j + 1]) != vs.end()) {
        lns.push_back(edges[j]);
        lns.push_back(edges[j + 1]);
      }

    unsigned int num_edges = faces[i].size();
    faces[i].clear();
    faces[i].push_back(lns[0]);
    int pt = lns[1];
    for (unsigned int j = 1; j < num_edges; j++) {
      faces[i].push_back(pt);
      for (unsigned int k = j; k < num_edges; k++) {
        if (lns[k * 2] == pt) {
          pt = lns[k * 2 + 1];
          swap(lns[j * 2], lns[k * 2 + 1]);
          swap(lns[j * 2 + 1], lns[k * 2]);
          break;
        }
        else if (lns[k * 2 + 1] == pt) {
          pt = lns[k * 2];
          swap(lns[j * 2], lns[k * 2]);
          swap(lns[j * 2 + 1], lns[k * 2 + 1]);
          break;
        }
      }
    }
  }
  geom.orient();

  return 1;
}

// color faces by N sides
void color_off(Geometry &geom)
{
  auto *col_map = new ColorMapMap;
  col_map->set_col(3, Color(255, 0, 0));    // red
  col_map->set_col(4, Color(0, 255, 0));    // green
  col_map->set_col(5, Color(0, 0, 255));    // blue
  col_map->set_col(6, Color(255, 255, 0));  // yellow
  col_map->set_col(8, Color(255, 0, 255));  // magenta
  col_map->set_col(10, Color(0, 255, 255)); // cyan

  const vector<vector<int>> &faces = geom.faces();
  for (unsigned int i = 0; i < faces.size(); i++) {
    Color col = col_map->get_col(faces[i].size());
    geom.colors(FACES).set(i, col);
  }
}

int main(int argc, char *argv[])
{
  kaleido_opts opts;
  opts.process_command_line(argc, argv);

  FILE *ofile = stdout; // write to stdout by default
  if (opts.just_list || opts.need_coordinates || opts.need_approx ||
      opts.model == 2) {
    if (opts.ofile != "") {
      ofile = fopen(opts.ofile.c_str(), "w");
      if (ofile == 0)
        opts.error("could not open output file \'" + opts.ofile + "\'");
    }
  }

  Uniform *uniform;
  int last_uniform = get_uniform_list(&uniform);

  char sym[MSG_SZ] = {0};
  Polyhedron *P = nullptr;

  int first = 1;
  int last = 1;
  if (opts.symbol == "0") {
    first = 1;
    last = last_uniform;
  }
  else {
    sprintf(sym, "%s", opts.symbol.c_str());
    P = kaleido(sym, uniform, last_uniform);
    if (P == NULL)
      Err("kaleido function failed");
    if (P->index != -1) {
      first = P->index + 1;
      last = first;
    }
  }

  for (int i = first; i <= last; i++) {
    // size of char strings used in kaleido
    char sym[MSG_SZ] = {0};
    if (opts.symbol == "0") {
      sprintf(sym, "#%d", i);
      P = kaleido(sym, uniform, last_uniform);
    }

    // P->name and P->dual_name are not set. do it here
    Malloc(P->name, strlen(uniform[i - 1].name) + 1, char);
    strcpy(P->name, uniform[i - 1].name);
    Malloc(P->dual_name, strlen(uniform[i - 1].dual) + 1, char);
    strcpy(P->dual_name, uniform[i - 1].dual);

    // if listing
    if (opts.just_list || opts.need_coordinates || opts.need_approx) {
      // high number calls more() less often, must be 2 or higher or will cause
      // crash
      int more_lines = 1000000;
      if (opts.need_approx) {
        // make a temporary copy of Polyhedron so P will not change
        Polyhedron Q = *P;
        newton(&Q, opts.need_approx, ofile);
      }
      printit(P, opts.need_coordinates, opts.just_list, opts.sig_digits,
              more_lines, uniform, ofile);
    }
    else {
      // off model
      if (opts.model == 1) {
        Geometry geom;
        get_poly(geom, P);

        if (opts.base == 2) {
          // center on origin
          geom.transform(Trans3d::translate(-centroid(geom.verts())));
          Vec3d cent = centroid(geom.verts());

          // from std_polys.cc
          GeometryInfo info(geom);
          info.set_center(cent);
          double rad = info.iedge_dist_lims().sum / info.num_iedges();
          Geometry dual;
          const double inf = 1200;
          get_dual(dual, geom, rad, cent, inf);
          add_extra_ideal_elems(dual, cent,
                                0.95 * inf); // limit closer than inf

          geom = dual;
          geom.orient();
        }

        // color faces by N sides
        color_off(geom);

        opts.write_or_error(geom, opts.ofile, opts.sig_digits);
      }
      // vrml model with original Kaleido code
      else if (opts.model == 2) {
        if (opts.base == 1)
          vrmodel(P, P->v, P->V, P->f, P->F, P->name, "", opts.sig_digits,
                  opts.azimuth, opts.elevation, opts.freeze, ofile);
        else if (opts.base == 2)
          vrmodel(P, P->f, P->F, P->v, P->V, P->dual_name, "*", opts.sig_digits,
                  opts.azimuth, opts.elevation, opts.freeze, ofile);
      }
    }
  }

  return 0;
}
