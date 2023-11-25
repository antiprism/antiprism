/*
   Copyright (c) 2003-2023, Adrian Rossiter, Roger Kaufman
   Includes ideas and algorithms by George W. Hart, http://www.georgehart.com

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
   Name: canonical.cc
   Description: canonicalize a polyhedron
                Uses George Hart's two canonicalization algorithm
                http://library.wolfram.com/infocenter/Articles/2012/
                http://www.georgehart.com/virtual-polyhedra/conway_notation.html
   Project: Antiprism - http://www.antiprism.com
*/

#include "../base/antiprism.h"
#include "canonical_common.h"
#include "color_common.h"

#include <cstdio>
#include <string>
#include <vector>

using std::string;
using std::vector;

using namespace anti;

// constants for color subscripts
namespace {
/*
  -O <args> output b - base, d - dual, i - intersection points (default: b)
               edge nearpoints, n - base, m - dual; C - base/dual convex hull
               edge nearpoints centroid, p - base, q - dual; o - origin point
               tangent sphere, u - minimum, U - maximum
               incircles, s - base, t - dual; as rings, S - base, T - dual
*/
unsigned const col_b = 0;
unsigned const col_d = 1;
unsigned const col_i = 2;
unsigned const col_n = 3;
unsigned const col_m = 4;
unsigned const col_c = 5;
unsigned const col_p = 6;
unsigned const col_q = 7;
unsigned const col_o = 8;
unsigned const col_u = 9;
unsigned const col_s = 10;
unsigned const col_t = 11;
unsigned const col_total = 12;
}; // namespace

int find_assignment(const string &assignment)
{
  string assigned_chars("bdinmcpqoust");
  return ((int)assigned_chars.find(assignment));
}

class cn_opts : public ProgramOpts {
public:
  IterationControl it_ctrl;

  string ifile;
  string ofile;

  char edge_distribution = '\0';    // can project onto sphere
  string shuffle_model_idxs;        // shuffle indexes can be vef
  char target_model = 'b';          // work on base model is default
  char planarize_method = '\0';     // no algorithm for planar is default
  int num_iters_planar = -1;        // unlimited iterations for planar
  char canonical_method = 'c';      // circle packings algorithm is default
  int num_iters_canonical = -1;     // unlimited iterations for canonical
  bool use_symmetry = false;        // don't use symmetry
  bool realign = false;             // realign output to full symmetry
  double edge_factor = NAN;         // mathematica algorithm variable
  double plane_factor = NAN;        // mathematica algorithm variable
  double radius_range_percent = -1; // percent expansion of model
  double factor = NAN;              // initial adjustment factor
  double factor_max = NAN;          // maximum adjustment factor
  char initial_point_type = 'c';    // c - edge centroids, n - edge near points
  string output_parts = "b";        // parts of output model
  double offset = 0;                // incircle offset from faces
  int roundness = 8;                // roundness of tangency sphere

  OffColor off_color[col_total] = {OffColor("")}; // for color/operator input
  int off_opacity[col_total][3]; // for transparency element's operator

  double eps = anti::epsilon;

  cn_opts() : ProgramOpts("canonical")
  {
    it_ctrl.set_max_iters(-1);
    it_ctrl.set_status_checks("1000,1");
    it_ctrl.set_sig_digits(int(-log(anti::epsilon) / log(10) + 0.5));
  }

  void process_command_line(int argc, char **argv);
  void usage();

  void parse_color_string(const char *optarg, const char c);
};

void extended_help()
{
  fprintf(stdout, R"(
Calculating Canonical Polyhedra (description of mathematica algorithm)

A relaxation algorithm is presented to determine a canonical form for an
arbitrary convex polyhedron.

by George W. Hart

One important role of high-level mathematical software such as Mathematica is
that it easily allows for the testing of experimental algorithms. Here we
explore a method of finding a canonical form of a polyhedron. The canonical
form is a polyhedron topologically equivalent to an input polyhedron, but with
all edges tangent to a unit sphere and with the center of gravity of the
tangent points being the origin.  A variety of examples will illustrate this
procedure. One application is that the algorithm constructs a geometrically
self-dual polyhedron given one which is only combinatorially self-dual.

A theorem of Schramm [Schramm 1992] states that for any given convex polyhedron
or 3-connected planar graph, there is a topologically equivalent polyhedron
with the following properties:

 1) each edge is tangent to the unit sphere,

 2) the center of gravity of the points of tangency is the origin.

The solution is unique, up to rotations and reflections of the sphere and so
provides a kind of canonical representation. Although we will not pursue it
here, the theorem is actually more general in that it allows the sphere to be
replaced with an arbitrary smooth convex body, e.g., an egg-shaped
quasi-ellipsoid. The theorem is closely related to the Koebe-Andreev-Thurston
Circle Packing theorem for planar graphs. See Ziegler [Ziegler 1995, p. 118]
for discussion and other references. 

For example, if our input polyhedron is a geometrically distorted (but
topologically unchanged) form of any Platonic or Archimedean solid, and we
calculate its canonical form, our algorithm should output its undistorted form
centered at the origin, with a midradius of 1, and with an arbitrary rotation.
E.g., given any parallelepiped as input, the canonical form output is a cube of
edge Sqrt[2]. For an arbitrary polyhedron, the canonical form is a way of
illustrating its combinatorial or topological structure, which often lets one
immediately see and understand its structure and symmetry.  

As an illustration, we will use the algorithm to find three geometrically
self-dual polyhedra, given starting points which are only combinatorially
self-dual. It follows from the above that the dual polyhedron also shares
properties (1) and (2). Thus, for these three self-dual examples, we can make
a compound of the polyhedron with its dual, i.e., with itself, in which both
polyhedra have the same points of tangencies and at these points their edges
cross each other's at right angles. 

The proofs cited above are of existence and not constructive, so we are
interested in an algorithm for determining the canonical form of an input
polyhedron. The algorithm proposed here operates by relaxation to iteratively
move the vertices of the given polyhedron along a trajectory which converges
at the canonical form. Although we begin and end with a polyhedron, during this
relaxation, the object defined by the vertices is likely not to be a polyhedron
geometrically.  Sets of points which belong to a face (combinatorially) are
likely not to be coplanar. So we add a third condition for a solution, to the
two above:

 3) the faces are planar.


Algorithm

Our algorithm inputs a polyhedron and iteratively adjusts its vertices to
slightly improve its conformance with the three conditions above. Within a
couple of hundred iterations, it typically finds all the conditions are
satisfied within a small tolerance, and stops. Three simple operations are all
that is needed, corresponding to the three conditions:

 A) For each edge, the closest point to the origin is determined; call it p.
    If p lies at unit distance from the origin, condition (1) is satisfied for
    that edge. If not, a small fraction of p is added to the two vertices
    which define the edge, (in proportion to the sign and amount of the error)
    so that at the next iteration the edge will be closer to tangency with the
    unit sphere. 

 B) The center of gravity of all the points p is determined. If it is zero,
    condition (2) is satisfied. If not, it is subtracted from all the vertices. 

 C) For each face, if the vertices lie in some plane in space, condition (3) is
    satisfied. If not, a plane which approximates it is computed. Each vertex
    of the face is then moved along a normal towards the plane. 

Iterating, it sometimes takes many steps for a correction at one part of the
polyhedron to percolate its way around and equalize the conditions everywhere.
On the examples below, between 50 and 100 iterations were sufficient to have
all conditions satisfied within a tolerance of 10^-8. Although the individual
steps are quite simple, this takes a minute or so in these examples and would
take longer on more complex polyhedra. Refinement and optimization are left as
future work.


Additional Work by Adrian Rossiter:

Edge near-point / circle-packing canonicalisation algorithm
===========================================================

Approach
--------

A polyhedron with a midsphere corresponds to two circle packings on the
same sphere: the incircles of the base faces and the incircles of the dual
faces. The circles from the two packs intersect at the edge tangency points
of the base (or dual) polyhedron in two orthogonal tangent pairs. By moving
a set of points until they satisfy the conditions for being the intersection
points of these circles, a model corresponding to a polyhedron with a
midsphere can be made.

If the midsphere is the unit sphere, and the edge tangency points have
their centroid at the sphere centre then the polyhedron is canonical.


Processing model
----------------

The processing model has a vertex for each edge tangency point, and a face
for each circle of the two packs. Each vertex is therefore surrounded by four
faces, and these faces correspond, in opposing pairs, to each of the two
circle packs.

The canonical model is solved when the processing model satisfies
the following conditions:

*  the vertices are at distance 1 from the origin

*  the vertices have their centroid at the origin

*  the faces are planar, hence each set of vertices corresponding to a
   base/dual face lies on a circle of the base/dual circle pack

*  each vertex lies on the two planes through the origin containing the
   normals of opposing faces, hence each pair of circles meeting at the
   vertex are tangent (and this is sufficient to ensure that they are
   also orthogonal)


Algorithm
---------

Initialisation

1. Translate the base model to carry the vertex centroid to the origin.

2. Converted to an 'ambo' form. The vertices are truncated to a single
   point on each edge. These points are initially set to either the
   centroid of the edge vertices (better for a general input), or the
   point on the edge line nearest to the origin (better for a
   near-canonical input).

3. Make a list of the cycle of four faces around each vertex. Alternate
   faces will be opposing faces.

4. Choose a small termination value, that if the minimum vertex
   movement is less than this then the iteration terminates.

5. The amount of vertex movement is controlled by an adjustment factor,
   which can change each iteration. Choose a starting value (e.g. .1),
   and a maximum value (e.g. .5).


Iteration

An offset will be calculated for each vertex, and will be applied
near the end of the iteration loop.

For each vertex:

1. Initialize the offset to zero.

2. Adjust for the tangency point centroid:
   Add
      -vertices_centroid
   to the offset.

3. Adjust for coplanar / circular points:
   Calculate the projection the vertex onto its four surrounding planes,
   and then calculate the centroid of these four projection points. Add
      (projection_point_centroid - vertex) * factor
   to the offset.

4. Adjust for mutual tangency / orthogonality:
   For each pair of opposing faces:
      find a normal to their normals (a base or dual edge direction)
      calculate the projection of the vertex onto a plane through the origin
      with this normal
   Add
      (projection_point_centroid - vertex) * factor
   to the offset.

5. Adjust for unscrambling:
   If a vertex does not lie inside the cycle of its four neighbouring
   vertices, make the following adjustment. Add
      (neighbour_vertices_centroid - vertex) * 0.5 * factor
   to the offset.

6. Add the offset to the vertex

7. Scale the vertex to have a length of exactly 1

8. Adjust the adjustment factor:
   If the maximum offset length is less than that of the last iteration
   then scale the factor by 1.01 (if this will not exceed the maximum
   specified), but if it is greater then scale the factor by 0.995.

9. Terminate iteration if the maximum offset length is less than the
   termination value.


Final model

The processing model has faces that correspond to base faces and those
that correspond to dual faces, which also correspond to base vertices as
they were the faces produced by vertex truncation.

The final base model retains the faces of the original base model, but
each vertex is set to the polar reciprocal of the corresponding processing
model face plane.

For each vertex in the base model

1. Determine the corresponding processing model face.

2. Calculate the point nearest to the origin on the face plane, and its
   distance from the origin, and the final vertex position is
     position = point / distance^2


Symmetry optimisation
---------------------

The algorithm is suitable for use with a symmetry optimisation. This
also forces the original symmetry to be maintained, as repeated processing
of the vertices may otherwise cause them to wander from the original
symmetry.

Use the symmetry group of the base model. In the processing model just one
vertex from each orbit is processed. The faces surrounding the vertex have
their other vertex positions calculated, once per iteration, by a symmetry
transformation of a processed vertex.

To avoid calculating all the vertex positions for the centroid, it can be
calculated as: the centroid of the processed vertices, each projected onto
the subspace left invariant by the subgroup that fixes it, and weighted by
the number of vertices in its orbit.

-------------------------------------------------------------------------------
From George Hart: https://www.georgehart.com/virtual-polyhedra/canonical.html

An interesting theorem states that there exists a "canonical form" of any given
convex polyhedron. This canonical form is a possibly distorted version of the
given polyhedron in which the vertices are positioned in space to satisfy the
following properties:

1) all the edges are tangent to the unit sphere,
2) the origin is the center of gravity of the points at which the
   edge touch the sphere,
3) the faces are flat (i.e. the vertices of each face lie in some plane),
   but are not necessarily regular.

It follows that a dual to the canonical polyhedron can be constructed which has
the above three properties as well, and the edges of the canonical polyhedron
and its dual cross at right angles. The representation is unique except for its
rotations and reflections.

Antiprism Note: These properties are measured for making a canonical model.
While the algorithms calculate vertex locations within a distance given by
-l n, meaning 10 to the -n power, the mathematical error of the measure of
intersections between the base model and its dual accumulates. The precision
of the intersections is measured to be within 10 to the -(n-2) power. The
mathematical error of the planar measure of a face is multiplied by the face's
number of edges. To achieve a measurably planar model may be beyond precision
limits of the processor.

Tip: If the dual of a model has faces with less sides, it may easier to achieve
a planar model by instead working on the dual using option -t c. Face side
counts can be found with off_report -C s

If a canonical base is produced as an off file, the reciprocal will need to be
generated using the midsphere with pol_recip -c M

)");
}

void cn_opts::usage()
{
  fprintf(
      stdout, R"(
Usage: %s [options] [input_file]

Read a polyhedron from a file in OFF format. Canonicalize or planarize it.
Uses algorithms by George W. Hart, http://www.georgehart.com/
http://www.georgehart.com/virtual-polyhedra/conway_notation.html
http://www.georgehart.com/virtual-polyhedra/canonical.html
If input_file is not given the program reads from standard input.

Options
%s
  -H        documention on algorithm
  -z <nums> number of iterations between status reports (implies termination
            check) (0 for final report only, -1 for no report), optionally
            followed by a comma and the number of iterations between
            termination checks (0 for report checks only) (default: %d,%d)
  -l <lim>  minimum distance change to terminate, as negative exponent
               (default: %d giving %.0e)
            WARNING: high values can cause non-terminal behaviour. Use -n
  -o <file> write output to file (default: write to standard output)

Canonical and Planarization Options
  -e <opt>  edge distribution (default : none, canonicalization only)
               s - project vertices onto a sphere
  -s <opt>  shuffle model indexes
               v - vertices, e - edges, f - faces, a - all (default: none)
  -t <opt>  target model
               b - work on base only (default)
               p - work on dual for planarization only
               c - work on dual for planarization and canonicalization
  -p <opt>  planarization only, or initial planarization if -c is set
              (default: no planarization)
               p - poly_form planarize (poly_form -a p)
               m - mathematica planarize
               a - moving edge planarize
               b - base/dual (reciprocate on face centroids magnitude squared)
  -i <itrs> maximum planarize iterations. -1 for unlimited (default: %d)
            WARNING: unstable models may not finish unless -i is set
  -c <opt>  canonicalization
               c - circle packings (default)
               m - mathematica version
               a - moving edge version
               b - base/dual version (reciprocate on face normals)
               x - none (default, if -p is set)
  -n <itrs> maximum canonical iterations. -1 for unlimited (default: %d)
            WARNING: unstable models may not finish unless -n is set
  -d <perc> radius test. percent difference between minimum and maximum radius
               checks if polyhedron is collapsing. 0 for no test 
               (default: 80 for canonicalizing, not used for planarizing)
  -y        maintain symmetry and alignment of the base model during
               processing (with -p p, -c c) 

Extra Options
  -Q <perc> percentage to scale edge tangency (default: 50) (-c m)
  -P <perc> percentage to scale face planarity (default: 20) (-c m, -p m, -p p)
  -f <adj>  initial percent adjustment factor, optionally followed by a comma
            and a maximum percent adjustment (default: 1,50) (-c c)
  -C        continue processing a near-canonical model (the initial
            intermediate processing model will preserves the geometry
            of the base model rather than avoid scrambling) (-c c)
            
Scene Options
  -O <args> output b - base, d - dual, i - intersection points (default: b)
               edge nearpoints, n - base, m - dual; c - base/dual convex hull
               edge nearpoints centroid, p - base, q - dual; o - origin point
               tangent sphere, u - minimum, U - maximum
               incircles, s - base, t - dual; as rings, S - base, T - dual
  -q <dist> incircles offset to avoid coplanarity e.g 0.0001 (default: 0)
  -g <opt>  roundness of tangent sphere, positive integer n (default: 8)
  -Y        align output model geometry to full symmetry

Coloring Options (run 'off_util -H color' for help on color formats)
keyword: none - sets no color, u - input base model, color of element unchanged
  -F <col>  face color/operation, assignment (required), transparency (optional)
              assignments to are: b - base, d - dual, i - intersection points
                edge nearpoints, n - base, m - dual; c - base/dual convex hull
                edge nearpoints centroid, p - base, q - dual; o - origin point
                u - tangent sphere, incircles (or rings), s - base, t - dual
                (one element, required; multiple -F as needed)
              -F d,b - face color by convexity, -F b,d - color from base verts
              -F s,f take color from faces, -F t,f take color from faces
              defaults: b - 'u' unchanged, d - 'b' take color from base verts,
                i - yellow, n,p - red, m,q - darkgreen, c - white,
                o - yellow, u - white, s,t - 'f' color of face
  -E <col>  edge color (same format as for faces) (defaults: unchanged)
            -E d,b - edge color by convexity
  -V <col>  vertex color (same format as for faces) (defaults: unchanged)
            transparency: valid range from 0 (invisible) to 255 (opaque)
  -m <maps> a comma separated list of color maps used to transform color
            indexes (default: colorful), a part consisting of letters from
            v, e, f, selects the element types to apply the map list to
            (default 'vef'). use map name of 'index' to output index numbers
               convexity:  white,gray50,gray25 (for -F d,b, -E d,b)

)",
      prog_name(), help_ver_text, it_ctrl.get_status_check_and_report_iters(),
      it_ctrl.get_status_check_only_iters(), it_ctrl.get_sig_digits(),
      it_ctrl.get_test_val(), it_ctrl.get_max_iters(), it_ctrl.get_max_iters());
}

void cn_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;
  int num;

  Split parts;
  Color col;

  bool p_set = false;
  bool c_set = false;

  off_color[col_d].set_f_col_op('b');
  off_color[col_s].set_f_col_op('f');               // take colors from faces
  off_color[col_t].set_f_col_op('f');               // take colors from faces
  off_color[col_i].set_v_col(Color(255, 255, 0));   // yellow
  off_color[col_n].set_v_col(Color(255, 0, 0));     // red
  off_color[col_m].set_v_col(Color(0, 100, 0));     // darkgreen
  off_color[col_p].set_v_col(Color(255, 0, 0));     // red
  off_color[col_q].set_v_col(Color(0, 100, 0));     // darkgreen
  off_color[col_o].set_v_col(Color(255, 255, 0));   // yellow
  off_color[col_c].set_f_col(Color(255, 255, 255)); // white
  off_color[col_u].set_f_col(Color(255, 255, 255)); // white
  off_color[col_u].set_e_col(Color::invisible);     // invisible
  off_color[col_u].set_v_col(Color::invisible);     // invisible

  // initialize opacity array
  for (int i = 0; i < 3; i++)
    std::fill(off_opacity[i], off_opacity[i] + col_total, -1);

  handle_long_opts(argc, argv);

  while ((c = getopt(argc, argv,
                     ":hHe:s:t:p:i:c:n:yO:q:g:Q:P:f:Cd:Yz:V:E:F:m:l:o:")) !=
         -1) {
    if (common_opts(c, optopt))
      continue;

    switch (c) {
    case 'H':
      extended_help();
      exit(0);

    case 'e':
      if (strlen(optarg) == 1 && strchr("s", int(*optarg)))
        edge_distribution = *optarg;
      else
        error("edge_distribution method type must be s", c);
      break;

    case 's':
      if (strspn(optarg, "vefa") != strlen(optarg))
        error(msg_str("shuffle parts are '%s' must be any or all from "
                      "v, e, f, a",
                      optarg),
              c);
      shuffle_model_idxs = optarg;
      break;

    case 't':
      if (strlen(optarg) == 1 && strchr("bpc", int(*optarg)))
        target_model = *optarg;
      else
        error("target model type must be b, p, c", c);
      break;

    case 'p':
      p_set = true;
      if (strlen(optarg) == 1 && strchr("pmab", int(*optarg)))
        planarize_method = *optarg;
      else
        error("planarize method type must be p, a, m, b", c);
      break;

    case 'i':
      print_status_or_exit(read_int(optarg, &num_iters_planar), c);
      if (num_iters_planar < -1)
        error("number of iterations for planarization must be -1 or greater",
              c);
      break;

    case 'c':
      c_set = true;
      if (strlen(optarg) == 1 && strchr("cmabx", int(*optarg)))
        canonical_method = *optarg;
      else
        error("canonical method type must be c, m, a, b, x", c);
      break;

    case 'n':
      print_status_or_exit(read_int(optarg, &num_iters_canonical), c);
      if (num_iters_canonical < -1)
        error("number of iterations for canonical must be -1 or greater", c);
      break;

    case 'y':
      use_symmetry = true;
      break;

    case 'Y':
      realign = true;
      break;

    case 'O':
      if (strspn(optarg, "bdinmopqsStTuUc") != strlen(optarg))
        error(msg_str("output parts are '%s' must be any or all from "
                      "b, d, i, n, m, o, p, q, s, S, t, T, u, U, c",
                      optarg),
              c);
      output_parts = optarg;
      break;

    case 'q':
      print_status_or_exit(read_double(optarg, &offset), c);
      break;

    case 'g':
      print_status_or_exit(read_int(optarg, &roundness), c);
      break;

    case 'Q':
      print_status_or_exit(read_double(optarg, &edge_factor), c);
      if (edge_factor <= 0 || edge_factor >= 100)
        warning("edge factor not inside range 0 to 100", c);
      break;

    case 'P':
      print_status_or_exit(read_double(optarg, &plane_factor), c);
      if (plane_factor <= 0 || plane_factor >= 100) {
        warning("plane factor not inside range 0 to 100", c);
      }
      break;

    case 'f': {
      vector<double> nums;
      print_status_or_exit(read_double_list(optarg, nums), c);
      if (nums.size() < 1 || nums.size() > 2)
        error(msg_str("must give one or two numbers (%lu were given)",
                      (unsigned long)nums.size()),
              c);
      factor = nums[0];
      if (factor < 0 || factor > 100)
        warning("initial factor not in range 0 to 100", c);
      factor_max = 110;
      if (nums.size() > 1) {
        factor_max = nums[1];
        if (factor_max < 0)
          warning("maximum factor negative", c);
      }
      if (factor_max < factor) {
        warning("maximum factor less than initial factor, setting to initial "
                "factor value",
                c);
        factor_max = factor;
      }
      break;
    }

    case 'C':
      initial_point_type = 'n';
      break;

    case 'd':
      print_status_or_exit(read_double(optarg, &radius_range_percent), c);
      if (radius_range_percent < 0)
        error("percentage must not be less than 0", c);
      break;

    case 'z':
      print_status_or_exit(it_ctrl.set_status_checks(optarg), c);
      break;

    case 'V':
      parse_color_string(optarg, c);
      break;

    case 'E':
      parse_color_string(optarg, c);
      break;

    case 'F':
      parse_color_string(optarg, c);
      break;

    case 'm':
      print_status_or_exit(read_colorings(off_color[col_b].clrngs, optarg), c);
      break;

    case 'l':
      print_status_or_exit(read_int(optarg, &num), c);
      print_status_or_exit(it_ctrl.set_sig_digits(num), c);
      break;

    case 'o':
      ofile = optarg;
      break;

    default:
      error("unknown command line error");
    }
  }

  // if planarizing only do not canonicalize
  if (p_set && !c_set) {
    canonical_method = 'x';
    warning("-c x in force, planarizing only", 'c');
  }

  if (target_model == 'p' && !planarize_method)
    error("target is for planarization but no method is selected", 't');

  if (target_model == 'c' && canonical_method == 'x')
    error("target is for canonicalization but no method is selected", 't');

  // set default variables
  if (canonical_method == 'm') {
    if (std::isnan(edge_factor))
      edge_factor = 50.0;
  }
  else if (!std::isnan(edge_factor))
    warning("set, but not used for this algorithm", 'E');

  if (canonical_method == 'm' || planarize_method == 'm' ||
      planarize_method == 'p') {
    if (std::isnan(plane_factor))
      plane_factor = 20.0;
  }
  else if (!std::isnan(plane_factor))
    warning("set, but not used for this algorithm", 'P');

  if (canonical_method == 'c') {
    if (std::isnan(factor))
      factor = 1.0;
    if (std::isnan(factor_max))
      factor_max = 50.0;
  }
  else if (!std::isnan(factor))
    warning("set, but not used for this algorithm", 'F');

  if ((canonical_method == 'x') && planarize_method &&
      (radius_range_percent > -1))
    warning("set, but not used for planarization", 'd');

  if (use_symmetry) {
    if (planarize_method && planarize_method != 'p')
      warning("set, but not used for this planarize algorithm", 'y');
    if (canonical_method != 'c')
      warning("set, but not used for this canonical algorithm", 'y');
  }

  if (argc - optind > 1)
    error("too many arguments");

  if (argc - optind == 1)
    ifile = argv[optind];

  // fill in missing maps
  string default_map_name = "convexity";
  for (unsigned int i = 0; i < 3; i++) {
    string map_name = default_map_name;
    // if map is already set, skip
    if (off_color[col_b].clrngs[i].get_cmaps().size())
      continue;
    off_color[col_b].clrngs[i].add_cmap(colormap_from_name(map_name.c_str()));
  }

  eps = it_ctrl.get_test_val();
}

// parse color entry. assingment is required
// 2 args - color name/index/operation, assignment
// 3 args - color name/index/operation, assignment, transparency
// 4 args - color r,g,b/index/operation, assignment
// 5 args - color r,g,b/index/operation, assignment, transparency
// assign color or operation to assingment character in optarg
// c is parameter character to display if an error is displayed
void cn_opts::parse_color_string(const char *optarg, const char c)
{
  Split parts(optarg, ",");
  unsigned int parts_sz = parts.size();

  int dummy;
  Status is_int = read_int(parts[0], &dummy);
  if (is_int.is_ok()) {
    if (parts_sz == 3)
      error("assigment is required, got only color", c);
    else if (parts_sz > 5)
      error("the argument is a color value and has more than 5 parts", c);
  }
  else if (!is_int.is_ok()) {
    if (parts_sz == 1)
      error("assigment is required, got only color", c);
    else if (parts_sz > 3)
      error("the argument is a color name and has more than 3 parts", c);
  }

  Color col;
  bool valid_entry = false;
  unsigned int next_parms_idx = 1;

  // see if entry is a valid r,g,b
  if (parts_sz >= 3) {
    string color_str_tmp = string(parts[0]) + "," + parts[1] + "," + parts[2];
    if (col.read(color_str_tmp.c_str())) {
      if (col.is_set())
        valid_entry = true;
      next_parms_idx = 3;
    }
  }

  // check to see if it is a valid color name
  if (!valid_entry) {
    if (col.read(parts[0])) {
      if (col.is_set())
        valid_entry = true;
      next_parms_idx = 1;
    }
  }

  // allow valid operators
  bool is_operator = false;
  if (!valid_entry) {
    // for faces:
    // u for b
    // d for b
    // b for d
    // f for s
    // f for t
    // borrow off_color[0]
    if (off_color[0].f_op_check((char *)parts[0], "udbf")) {
      is_operator = true;
      valid_entry = true; // don't check for color
    }
  }

  // if the color is not valid, output error
  if (!valid_entry)
    print_status_or_exit(col.read(parts[0]), c);

  // find what entity color or operator is assigned to
  string assignment;
  int assign_idx = -1;
  if (parts_sz > next_parms_idx) {
    // if the next part is an assignment
    assignment = parts[next_parms_idx];
    bool err = false;
    if (assignment.length() > 1)
      err = true;
    else {
      assign_idx = find_assignment(assignment);
      if (assign_idx < 0)
        err = true;
    }
    if (err)
      error(msg_str("assignment is '%s' and must be one character of 'b, c, d, "
                    "i, n, m, o, p, q, s, t, u'",
                    assignment.c_str()),
            c);
    else
      next_parms_idx++;
  }

  if (is_operator) {
    valid_entry = false; // reset
    // check for valid operators in context
    if (c == 'F') {
      // base colors unchanged
      if (!strcmp(parts[0], "u") && assignment == "b") {
        off_color[col_b].set_f_col_op('u');
        valid_entry = true;
      }
      // base color by convexity of faces
      else if (!strcmp(parts[0], "d") && assignment == "b") {
        off_color[col_b].set_f_col_op('d');
        valid_entry = true;
      }
      // if dual keeping colors from the base vertices, do nothing
      else if (!strcmp(parts[0], "b") && assignment == "d") {
        off_color[col_d].set_f_col_op('\0');
        valid_entry = true;
      }
      // if base incircles take colors from faces
      else if (!strcmp(parts[0], "f") && assignment == "s") {
        off_color[col_s].set_f_col_op('f');
        valid_entry = true;
      }
      // if dual incircles take colors from faces
      else if (!strcmp(parts[0], "f") && assignment == "t") {
        off_color[col_t].set_f_col_op('f');
        valid_entry = true;
      }
    }
    else if (c == 'E') {
      // base colors unchanged
      if (!strcmp(parts[0], "u") && assignment == "b") {
        off_color[col_b].set_e_col_op('u');
        valid_entry = true;
      }
      // base color by convexity of edges
      else if (!strcmp(parts[0], "d") && assignment == "b") {
        off_color[col_b].set_e_col_op('d');
        valid_entry = true;
      }
    }
    else if (c == 'V') {
      // base colors unchanged
      if (!strcmp(parts[0], "u") && assignment == "b") {
        off_color[col_b].set_v_col_op('u');
        valid_entry = true;
      }
    }

    if (!valid_entry)
      error(msg_str("invalid operator '%s' for assignment '%s'", parts[0],
                    assignment.c_str()),
            c);
  }

  // check for transparency
  int opacity = 255;
  if (parts_sz > next_parms_idx) {
    // if the next part is an integer
    if (read_int(parts[next_parms_idx], &opacity)) {
      if (opacity < 0 || opacity > 255)
        error(msg_str("transparency is '%d' and must be between 0 and 255",
                      opacity),
              c);
      if (!col.set_alpha(opacity))
        warning("transparency has no effect on map indexes or invisible", c);
      next_parms_idx++;
    }
  }

  if (c == 'V') {
    if (!is_operator)
      off_color[assign_idx].set_v_col(col);
    off_opacity[assign_idx][VERTS] = opacity;
  }
  else if (c == 'E') {
    if (!is_operator)
      off_color[assign_idx].set_e_col(col);
    off_opacity[assign_idx][EDGES] = opacity;
  }
  else if (c == 'F') {
    if (!is_operator)
      off_color[assign_idx].set_f_col(col);
    off_opacity[assign_idx][FACES] = opacity;
  }
  else
    // should never happend
    error("parse_color_string must be called from 'V', 'E' or 'F'", c);
}

bool is_nonoverlap_single_cover(const Geometry &geom)
{
  // Get vertices projected onto a sphere
  vector<Vec3d> verts;
  std::transform(geom.verts().cbegin(), geom.verts().cend(),
                 std::back_inserter(verts),
                 [](const Vec3d &v) { return v.unit(); });

  double area_sum = 0; // for sum of unsigned spherical areas
  for (const auto &face : geom.faces()) {
    const auto C = geom.face_cent(face).unit(); // face cent proj onto sphere
    const int f_sz = face.size();
    // Consider each edge of the face to make a triangle with the face
    // centroid
    for (int i = 0; i < f_sz; i++) {
      vector<Vec3d> tri = {verts[face[i]], verts[face[(i + 1) % f_sz]], C};
      // get the three angles and sum them
      double tri_angle_sum = 0.0;
      for (int tri_v = 0; tri_v < 3; tri_v++) {
        const auto x20 = vcross(tri[(2 + tri_v) % 3], tri[tri_v]);
        const auto x10 = vcross(tri[(1 + tri_v) % 3], tri[tri_v]);
        tri_angle_sum += acos(safe_for_trig(vdot(x20.unit(), x10.unit())));
      }
      area_sum += fabs(tri_angle_sum) - M_PI; // triangle area is angle excess
    }
  }
  return double_eq(area_sum, 4 * M_PI, 1e-5);
}

Geometry get_dual(const Geometry &base)
{
  Geometry dual;
  get_dual(dual, base, 1);
  return dual;
}

void check_model(const Geometry &geom, string s, const cn_opts &opts)
{
  double epsilon_local1 = 1e-4;     // coincident elements
  double epsilon_local2 = opts.eps; // face area (default is 1e-12)

  Geometry geom_merged = geom;
  merge_coincident_elements(geom_merged, "vef", 0, epsilon_local1);

  if (geom.verts().size() != geom_merged.verts().size())
    opts.warning(
        msg_str("possible coincident vertices found in %s", s.c_str()));

  if (geom.faces().size() != geom_merged.faces().size())
    opts.warning(msg_str("possible coincident faces found in %s", s.c_str()));

  GeometryInfo info(geom);
  vector<double> areas = info.get_f_areas();
  for (unsigned int i = 0; i < areas.size(); i++) {
    // fprintf(stderr, "%.17lf\n", areas[i]);
    if (areas[i] < epsilon_local2) {
      opts.warning(msg_str(
          "possible faces of near 0 area found in %s (face %d)", s.c_str(), i));
      break;
    }
  }
}

void check_coincidence(const Geometry &base, const Geometry &dual,
                       const cn_opts &opts)
{
  string s = "base";
  check_model(base, s, opts);

  s = "dual";
  check_model(dual, s, opts);
}

/*
// check convex hull for congruency to geom
bool check_hull(const Geometry &geom)
{
  Geometry hull = geom;
  hull.set_hull();
  //return (check_coincidence(geom, hull));
  return (geom.faces().size() == hull.faces().size() &&
          geom.verts().size() == hull.verts().size());
}

bool check_planarity(const Geometry &geom)
{
  bool convex = true;
  for (unsigned int i = 0; i < geom.faces().size(); i++) {
    vector<int> vec(1, i);
    Geometry polygon = faces_to_geom(geom, vec);
    polygon.set_hull();
    if (polygon.faces().size() > 1) {
// fprintf(stderr, "%d faces tested\n", i + 1);
      convex = false;
      break;
    }
  }
  return convex;
}
*/

void planarity_info(const Geometry &geom, const cn_opts &opts)
{
  GeometryInfo rep(geom);

  double max_nonplanar = std::numeric_limits<double>::min();
  // double sum_nonplanar = 0;
  unsigned int sz = geom.faces().size();
  for (unsigned int i = 0; i < sz; i++) {
    double nonplanar = rep.get_f_max_nonplanars()[i];
    // sum_nonplanar += nonplanar;
    if (nonplanar > max_nonplanar)
      max_nonplanar = nonplanar;
  }
  opts.message(msg_str("maximum nonplanarity = %.0e", max_nonplanar));
  // fprintf(stderr, "average_nonplanarity = %.17g\n", sum_nonplanar / sz);
  // fprintf(stderr, "isoperimetric quotient = %.17g\n",
  // rep.isoperimetric_quotient());

  return;
}

void generate_points(const Geometry &base, const Geometry &dual,
                     vector<Vec3d> &base_nearpts, vector<Vec3d> &dual_nearpts,
                     vector<Vec3d> &ips, const double &epsilon_local)
{
  vector<vector<int>> base_edges;
  vector<vector<int>> dual_edges;

  base.get_impl_edges(base_edges);
  dual.get_impl_edges(dual_edges);

  for (auto &base_edge : base_edges)
    base_nearpts.push_back(base.edge_nearpt(base_edge, Vec3d::zero));

  for (auto &dual_edge : dual_edges)
    dual_nearpts.push_back(dual.edge_nearpt(dual_edge, Vec3d::zero));

  for (unsigned int i = 0; i < base_edges.size(); i++) {
    Vec3d b0 = base.verts(base_edges[i][0]);
    Vec3d b1 = base.verts(base_edges[i][1]);
    for (unsigned int j = 0; j < dual_edges.size(); j++) {
      Vec3d d0 = dual.verts(dual_edges[j][0]);
      Vec3d d1 = dual.verts(dual_edges[j][1]);
      // does base edge intersect with dual edge?
      Vec3d intersection_point =
          segments_intersection(b0, b1, d0, d1, epsilon_local);
      if (intersection_point.is_set()) {
        ips.push_back(intersection_point);
        break;
      }
    }
  }
}

// find nearpoints error. center is unchanged
double edge_nearpoints_error(const Geometry &geom, double &min, double &max,
                             Vec3d &center)
{
  min = std::numeric_limits<double>::max();
  max = std::numeric_limits<double>::min();

  vector<vector<int>> edges;
  geom.get_impl_edges(edges);

  vector<Vec3d> near_pts;

  double nearpt_error = 0;
  for (auto &edge : edges) {
    Vec3d P = geom.edge_nearpt(edge, Vec3d::zero);
    near_pts.push_back(P);

    double l = fabs(P.len() - 1.0);
    nearpt_error += l;
    if (l < min)
      min = l;
    if (l > max)
      max = l;
  }

  center = centroid(near_pts);

  // adds computational error, not using average
  return nearpt_error / double(edges.size());
}

double nearpoint_report(const Geometry &geom, const vector<Vec3d> &nearpts,
                        string str, bool &perfect_score,
                        const double &epsilon_local, const cn_opts &opts)
{
  double min = 0;
  double max = 0;
  Vec3d center;
  edge_nearpoints_error(geom, min, max, center);
  opts.message(msg_str("%s range of edge nearpoint error is %.0e to %.0e",
                       str.c_str(), min, max));

  int radius_count = nearpts.size();
  int nearpts_size = nearpts.size();
  for (int i = 0; i < nearpts_size; i++) {
    double l = fabs(nearpts[i].len() - 1.0);
    if (double_ne(l, 0.0, epsilon_local))
      radius_count--;
  }
  if (radius_count != nearpts_size)
    perfect_score = false;

  double np_pct = radius_count / (double)nearpts_size * 100;
  opts.message(msg_str(
      "%d out of %d %s edge nearpoint radii lie within %.0e of length 1 "
      "(%g%%)",
      radius_count, nearpts_size, str.c_str(), epsilon_local, np_pct));

  double max_error = std::numeric_limits<int>::min();
  for (unsigned int i = 0; i < 3; i++) {
    double err = fabs(center[i]);
    if (err > max_error)
      max_error = err;
  }

  if (max_error > max)
    max = max_error;

  string lstr = "lies ";
  if (double_ne(max_error, 0.0, epsilon_local)) {
    lstr = "does not lie ";
    perfect_score = false;
  }
  opts.message(msg_str(
      "%s nearpoint centroid (%.0e, %.0e, %.0e) %swithin %.0e of the origin",
      str.c_str(), center[0], center[1], center[2], lstr.c_str(),
      epsilon_local));

  return max;
}

bool canonical_report(const Geometry &base, vector<Vec3d> &base_nearpts,
                      vector<Vec3d> &ips, const double &epsilon_local,
                      const cn_opts &opts)
{
  bool perfect_score = true;

  double max_error = std::numeric_limits<int>::min();

  // intersection point count is not used in score
  int ip_size = ips.size();
  int bn_size = base_nearpts.size();
  double pct = ip_size / (double)bn_size * 100;
  opts.message(
      msg_str("%d out of %d base/dual edge intersection points found (%g%%)",
              ip_size, bn_size, pct));
  if (ip_size != bn_size)
    perfect_score = false;

  max_error = nearpoint_report(base, base_nearpts, "base", perfect_score,
                               epsilon_local, opts);
  opts.message(
      msg_str("base canonical model maximum error is %.0e", max_error));

  // report planarity
  planarity_info(base, opts);
  // check if extra faces are produced on a convex hull
  string s = "no ";
  if (get_convexity(base) != Convexity::convex_strict) {
    s = "";
    perfect_score = false;
  }
  opts.message(msg_str("%sextra faces produced by convex hull", s.c_str()));

  string str;
  if (is_nonoverlap_single_cover(base))
    str = "no ";
  else
    perfect_score = false;
  opts.message(
      msg_str("base canonical model has %soverlap error", str.c_str()));

  return perfect_score;
}

// is model already canonical?
bool precheck(const Geometry &base, const double &epsilon_local)
{
  bool perfect_score = true;

  Geometry dual = get_dual(base);

  vector<Vec3d> base_nearpts;
  vector<Vec3d> dual_nearpts;
  vector<Vec3d> ips;
  generate_points(base, dual, base_nearpts, dual_nearpts, ips, epsilon_local);

  int ip_size = ips.size();
  int bn_size = base_nearpts.size();
  if (ip_size != bn_size) {
    perfect_score = false;
    return perfect_score;
  }

  int base_radius_count = base_nearpts.size();
  int base_nearpts_size = base_nearpts.size();
  for (int i = 0; i < base_nearpts_size; i++) {
    double l = fabs(base_nearpts[i].len() - 1.0);
    if (double_ne(l, 0.0, epsilon_local))
      base_radius_count--;
  }
  if (base_radius_count != base_nearpts_size) {
    perfect_score = false;
    return perfect_score;
  }

  // check if extra faces are produced on a convex hull
  if (get_convexity(base) != Convexity::convex_strict)
    perfect_score = false;

  return perfect_score;
}

Geometry unit_circle(int polygon_size, const Color &incircle_color, bool filled)
{
  Geometry incircle;
  double arc = deg2rad(360.0 / (double)polygon_size);

  double angle = 0.0;
  for (int i = 0; i < polygon_size; i++) {
    incircle.add_vert(Vec3d(cos(angle), sin(angle), 0.0), Color::invisible);
    incircle.add_edge(make_edge(i, (i + 1) % polygon_size),
                      (filled ? Color::invisible : incircle_color));
    angle += arc;
  }

  if (filled) {
    vector<int> face;
    for (int i = 0; i < polygon_size; i++)
      face.push_back(i);
    incircle.add_face(face, incircle_color);
  }

  return (incircle);
}

// get the minimum incircle radius. for canonicalized they are all equal
double incircle_radius(const Geometry &geom, Vec3d &center, const int face_no)
{
  vector<Vec3d> near_pts;

  // get the near points to measure length from center
  unsigned int fsz = geom.faces(face_no).size();
  for (unsigned int i = 0; i < fsz; i++) {
    int v1 = geom.faces(face_no)[i];
    int v2 = geom.faces(face_no)[(i + 1) % fsz];

    Vec3d P = geom.edge_nearpt(make_edge(v1, v2), Vec3d::zero);
    near_pts.push_back(P);
  }

  // get the minimum radius
  double radius = std::numeric_limits<double>::max();
  for (unsigned int i = 0; i < fsz; i++) {
    double l = (center - near_pts[i]).len();
    if (l < radius)
      radius = l;
  }

  return radius;
}

Geometry incircles(const Geometry &geom, const char incircle_color_method,
                   const Color &incircle_color, bool filled, double offset)
{
  Geometry incircles_geom;
  Geometry circle = unit_circle(60, incircle_color, filled);

  for (unsigned int i = 0; i < geom.faces().size(); i++) {
    // find incircle rotation place
    Vec3d face_centroid = anti::centroid(geom.verts(), geom.faces(i));
    Vec3d face_normal = face_norm(geom.verts(), geom.faces(i)).unit();
    Vec3d center = face_normal * vdot(face_centroid, face_normal);

    // find radius of incircle, and make incircle of radius
    Geometry incircle = circle;

    // color incircles here if based on face color
    if (incircle_color_method && strchr("fF", incircle_color_method)) {
      Color col = geom.colors(FACES).get(i);
      Coloring(&incircle).e_one_col(col);
      Coloring(&incircle).f_one_col(col);
    }

    // set the radius of the incircle
    double radius = incircle_radius(geom, center, i);
    incircle.transform(Trans3d::scale(radius));

    // set depth of incircle
    double depth = center.len();
    incircle.transform(Trans3d::translate(Vec3d(0, 0, depth)));

    // rotate incircle into place
    incircle.transform(Trans3d::rotate(Vec3d::Z, center));
    incircles_geom.append(incircle);
  }

  // scale whole geom instead of each indi
  if (offset)
    incircles_geom.transform(Trans3d::scale(offset + 1.0));

  return incircles_geom;
}

// average radius rather than maximum has more reliability than max
void unitize_vertex_radius(Geometry &geom)
{
  GeometryInfo info(geom);
  info.set_center(geom.centroid());
  // geom.transform(Trans3d::scale(1 / info.vert_dist_lims().max));
  double avg = info.vert_dist_lims().sum / info.num_verts();
  geom.transform(Trans3d::scale(1 / avg));
}

// models from incomplete processing and other situations be large or small
// set radius to 1 for get_dual call, if necessary
void reset_model_size(Geometry &geom, const double &epsilon_local,
                      cn_opts &opts)
{
  double radius = edge_nearpoints_radius(geom);
  if (double_ne(radius, 1.0, epsilon_local)) {
    opts.message(
        msg_str("Resetting base nearpoints average radius to within %.0e of 1",
                epsilon_local));
    unitize_nearpoints_radius(geom);
  }
}

// changed base is returned
void construct_model(Geometry &base, Geometry &dual,
                     vector<Vec3d> &base_nearpts, vector<Vec3d> &dual_nearpts,
                     vector<Vec3d> &ips, cn_opts &opts)
{
  // get statistics before model is changed, if needed
  double min = 0;
  double max = 0;
  Vec3d center;
  if (opts.output_parts.find_first_of("uU") != string::npos)
    edge_nearpoints_radius(base, min, max, center);

  // base elements
  char op = opts.off_color[col_b].get_f_col_op();
  if (op && strchr("uU", op)) {
    // apply transparency to operation
    Status stat = apply_transparency(base, opts.off_opacity[col_b][FACES]);
    if (stat.is_warning())
      opts.warning(stat.msg(), col_b);
    // clear operation
    opts.off_color[col_b].set_f_col_op('\0');
  }
  else if (op && strchr("dD", op)) {
    color_faces_by_convexity(base, opts.off_color[col_b].clrngs[FACES],
                             (op == 'D'), opts.eps);
    // apply transparency to operation
    Status stat = apply_transparency(base, opts.off_opacity[col_b][FACES]);
    if (stat.is_warning())
      opts.warning(stat.msg(), col_b);
    // clear operation
    opts.off_color[col_b].set_f_col_op('\0');
  }

  op = opts.off_color[col_b].get_e_col_op();
  if (op && strchr("uU", op)) {
    // apply transparency to operation
    Status stat =
        apply_transparency(base, opts.off_opacity[col_b][EDGES], EDGES);
    if (stat.is_warning())
      opts.warning(stat.msg(), col_b);
    // clear operation
    opts.off_color[col_b].set_e_col_op('\0');
  }
  else if (op && strchr("dD", op)) {
    color_edges_by_dihedral(base, opts.off_color[col_b].clrngs[EDGES],
                            (op == 'D'), opts.eps);
    // apply transparency to operation
    Status stat =
        apply_transparency(base, opts.off_opacity[col_b][EDGES], EDGES);
    if (stat.is_warning())
      opts.warning(stat.msg(), col_b);
    // clear operation
    opts.off_color[col_b].set_e_col_op('\0');
  }

  op = opts.off_color[col_b].get_v_col_op();
  if (op && strchr("uU", op)) {
    // apply transparency to operation
    Status stat =
        apply_transparency(base, opts.off_opacity[col_b][VERTS], VERTS);
    if (stat.is_warning())
      opts.warning(stat.msg(), col_b);
    // clear operation
    opts.off_color[col_b].set_v_col_op('\0');
  }

  opts.off_color[col_b].off_color_main(base);

  // dual elements
  opts.off_color[col_d].off_color_main(dual);
  op = opts.off_color[col_d].get_f_col_op();
  if (op && strchr("bB", op)) {
    // apply transparency to operation
    Status stat = apply_transparency(dual, opts.off_opacity[col_d][FACES]);
    if (stat.is_warning())
      opts.warning(stat.msg(), col_d);
    // clear operation
    opts.off_color[col_d].set_f_col_op('\0');
  }

  // base incircles
  op = opts.off_color[col_s].get_f_col_op();
  Color col = opts.off_color[col_s].get_f_col();
  Geometry base_incircles;
  if (opts.output_parts.find_first_of("sS") != string::npos) {
    bool filled = (opts.output_parts.find("s") != string::npos) ? true : false;
    base_incircles = incircles(base, op, col, filled, opts.offset);
  }

  // color base incircle elements, if op is f, clear operation
  if (op && strchr("fF", op)) {
    // apply transparency to operation
    Status stat =
        apply_transparency(base_incircles, opts.off_opacity[col_s][FACES]);
    if (stat.is_warning())
      opts.warning(stat.msg(), col_s);
    opts.off_color[col_s].set_f_col_op('\0');
  }
  opts.off_color[col_s].off_color_main(base_incircles);

  // dual incircles
  op = opts.off_color[col_t].get_f_col_op();
  col = opts.off_color[col_t].get_f_col();
  Geometry dual_incircles;
  if (opts.output_parts.find_first_of("tT") != string::npos) {
    bool filled = (opts.output_parts.find("t") != string::npos) ? true : false;
    dual_incircles = incircles(dual, op, col, filled, opts.offset);
  }

  // color dual incircle elements, if op is f, clear operation
  if (op && strchr("fF", op)) {
    // apply transparency to operation
    Status stat =
        apply_transparency(dual_incircles, opts.off_opacity[col_t][FACES]);
    if (stat.is_warning())
      opts.warning(stat.msg(), col_t);
    opts.off_color[col_t].set_f_col_op('\0');
  }

  // build convex hull of base and dual, color
  Geometry base_dual = base;
  if (opts.output_parts.find("c") != string::npos) {
    base_dual.append(dual);
    base_dual.set_hull();
    opts.off_color[col_c].off_color_main(base_dual);
  }

  // build the model...

  // clear base if not using
  if (opts.output_parts.find_first_of("b") == string::npos)
    base.clear_all();

  // append incircles, pre-colored
  if (base_incircles.verts().size())
    base.append(base_incircles);

  if (dual_incircles.verts().size())
    base.append(dual_incircles);

  // append dual
  if (opts.output_parts.find_first_of("d") != string::npos)
    base.append(dual);

  // add intersection points
  Geometry tmp;
  if (opts.output_parts.find("i") != string::npos) {
    for (auto &ip : ips)
      tmp.add_vert(ip);
    opts.off_color[col_i].off_color_main(tmp);
    base.append(tmp);
  }

  // add base near points
  tmp.clear_all();
  if (opts.output_parts.find("n") != string::npos) {
    for (auto &base_nearpt : base_nearpts)
      tmp.add_vert(base_nearpt);
    opts.off_color[col_n].off_color_main(tmp);
    base.append(tmp);
  }

  // add dual near points
  tmp.clear_all();
  if (opts.output_parts.find("m") != string::npos) {
    for (auto &dual_nearpt : dual_nearpts)
      tmp.add_vert(dual_nearpt);
    opts.off_color[col_m].off_color_main(tmp);
    base.append(tmp);
  }

  // add base near points centroid
  tmp.clear_all();
  if (opts.output_parts.find("p") != string::npos) {
    tmp.add_vert(centroid(base_nearpts));
    opts.off_color[col_p].off_color_main(tmp);
    base.append(tmp);
  }

  // add dual near points centroid
  tmp.clear_all();
  if (opts.output_parts.find("q") != string::npos) {
    tmp.add_vert(centroid(dual_nearpts));
    opts.off_color[col_q].off_color_main(tmp);
    base.append(tmp);
  }

  // add origin point
  tmp.clear_all();
  if (opts.output_parts.find("o") != string::npos) {
    tmp.add_vert(Vec3d::zero);
    opts.off_color[col_o].off_color_main(tmp);
    base.append(tmp);
  }

  // add unit sphere on origin
  tmp.clear_all();
  if (opts.output_parts.find_first_of("uU") != string::npos) {
    string geo_str = "geo_" + std::to_string(opts.roundness) + "_" +
                     std::to_string(opts.roundness);
    tmp.read_resource(geo_str);
    tmp.transform(Trans3d::translate(-centroid(tmp.verts())));
    unitize_vertex_radius(tmp);

    if (opts.output_parts.find("u") != string::npos)
      tmp.transform(Trans3d::scale(min));
    else
      tmp.transform(Trans3d::scale(max));

    opts.off_color[col_u].off_color_main(tmp);
    base.append(tmp);
  }

  // append convex hull
  if (opts.output_parts.find("c") != string::npos)
    base.append(base_dual);
}

vector<int> geom_deal(Geometry &geom, const int pack_size)
{
  Coloring clrng;
  clrng.set_geom(&geom);
  string map_name = "deal" + std::to_string(pack_size);
  ColorMap *pack = colormap_from_name(map_name.c_str());
  clrng.add_cmap(pack);

  vector<int> deal;
  for (int i = 0; i < pack->effective_size(); i++)
    deal.push_back(pack->get_col(i).get_index());

  return (deal);
}

void shuffle_model_indexes(Geometry &geom, const cn_opts &opts)
{
  // Geometry shuffled = geom;

  bool shuffle_verts =
      (opts.shuffle_model_idxs.find_first_of("va") != string::npos);
  bool shuffle_faces =
      (opts.shuffle_model_idxs.find_first_of("fa") != string::npos);
  bool shuffle_edges =
      (opts.shuffle_model_idxs.find_first_of("ea") != string::npos);

  map<int, int> new_verts;
  if (!shuffle_verts) {
    for (unsigned int i = 0; i < geom.verts().size(); i++)
      new_verts[i] = i;
  }
  else {
    opts.message("shuffle model indexes: vertices", "s");
    unsigned int sz = geom.verts().size();
    vector<int> deal = geom_deal(geom, sz);
    map<int, Color> new_cols;
    vector<Vec3d> shuffled_verts;
    for (unsigned int i = 0; i < sz; i++) {
      int v_new = deal[i];
      new_verts[v_new] = i;
      shuffled_verts.push_back(geom.verts(v_new));
      new_cols[i] = geom.colors(VERTS).get(v_new);
    }
    geom.clear(VERTS);
    for (unsigned int i = 0; i < sz; i++)
      geom.add_vert(shuffled_verts[i], new_cols[i]);
  }

  if (shuffle_faces || shuffle_verts) {
    unsigned int sz = geom.faces().size();
    map<int, int> face_order;
    if (!shuffle_faces) {
      for (unsigned int i = 0; i < sz; i++)
        face_order[i] = i;
    }
    else {
      opts.message("shuffle model indexes: faces", "s");
      vector<int> deal = geom_deal(geom, sz);
      for (unsigned int i = 0; i < sz; i++)
        face_order[deal[i]] = i;
    }
    map<int, Color> new_cols;
    vector<vector<int>> shuffled_faces;
    for (unsigned int i = 0; i < face_order.size(); i++) {
      vector<int> face;
      for (unsigned int j = 0; j < geom.faces(face_order[i]).size(); j++) {
        int v = geom.faces(face_order[i])[j];
        int v_new = new_verts[v];
        face.push_back(v_new);
      }
      shuffled_faces.push_back(face);
      new_cols[i] = geom.colors(FACES).get(face_order[i]);
    }
    geom.clear(FACES);
    for (unsigned int i = 0; i < sz; i++)
      geom.add_face(shuffled_faces[i], new_cols[i]);
  }

  if (!geom.edges().size()) {
    opts.warning("shuffle model indexes: no explicit edges", "s");
  }
  else if (shuffle_edges || shuffle_verts) {
    unsigned int sz = geom.edges().size();
    map<int, int> edge_order;
    if (!shuffle_edges)
      for (unsigned int i = 0; i < sz; i++) {
        edge_order[i] = i;
      }
    else {
      opts.message("shuffle model indexes: edges", "s");
      vector<int> deal = geom_deal(geom, sz);
      for (unsigned int i = 0; i < sz; i++)
        edge_order[deal[i]] = i;
    }
    map<int, Color> new_cols;
    vector<vector<int>> shuffled_edges;
    for (unsigned int i = 0; i < sz; i++) {
      vector<int> edge;
      for (unsigned int j = 0; j < geom.edges(edge_order[i]).size(); j++) {
        int v = geom.edges(edge_order[i])[j];
        int v_new = new_verts[v];
        edge.push_back(v_new);
      }
      shuffled_edges.push_back(edge);
      new_cols[i] = geom.colors(EDGES).get(edge_order[i]);
    }
    geom.clear(EDGES);
    for (unsigned int i = 0; i < sz; i++)
      geom.add_edge(shuffled_edges[i], new_cols[i]);
  }
}

// Implementation of George Hart's canonicalization algorithm
// http://library.wolfram.com/infocenter/Articles/2012/
// RK - the model will possibly become non-convex early in the loops.
// if it contorts too badly, the model will implode. Having the input
// model at a radius of near 1 minimizes this problem
bool canonicalize_mm(Geometry &geom, IterationControl it_ctrl,
                     const double edge_factor, const double plane_factor,
                     const double radius_range_percent,
                     const bool planarize_only)
{
  bool completed = false;
  it_ctrl.set_finished(false);

  vector<Vec3d> &verts = geom.raw_verts();

  vector<vector<int>> edges;
  geom.get_impl_edges(edges);

  double test_val = it_ctrl.get_test_val();
  double max_diff2 = 0;

  for (it_ctrl.start_iter(); !it_ctrl.is_done(); it_ctrl.next_iter()) {
    vector<Vec3d> verts_last = verts;

    if (!planarize_only) {
      vector<Vec3d> near_pts;
      for (auto &edge : edges) {
        Vec3d P = geom.edge_nearpt(edge, Vec3d::zero);
        near_pts.push_back(P);
        Vec3d offset = edge_factor * (P.len() - 1) * P;
        verts[edge[0]] -= offset;
        verts[edge[1]] -= offset;
      }

      // re-center for drift
      Vec3d cent_near_pts = centroid(near_pts);
      for (unsigned int i = 0; i < verts.size(); i++)
        verts[i] -= cent_near_pts;
    }

    // Make a copy of verts into vs and zero out
    // Accumulate vertex changes instead of altering vertices in place
    // This can help relieve when a vertex is pushed towards one plane
    // and away from another
    vector<Vec3d> vs = verts;
    for (auto &v : vs)
      v = Vec3d::zero;

    for (unsigned int f = 0; f < geom.faces().size(); f++) {
      if (geom.faces(f).size() == 3)
        continue;
      Vec3d face_normal = face_norm(geom.verts(), geom.faces(f)).unit();
      Vec3d face_centroid = geom.face_cent(f);
      // make sure face_normal points outward
      if (vdot(face_normal, face_centroid) < 0)
        face_normal *= -1.0;
      // place a planar vertex over or under verts[v]
      // adds or subtracts it to get to the planar verts[v]
      for (int v : geom.faces(f))
        vs[v] += vdot(plane_factor * face_normal, face_centroid - verts[v]) *
                 face_normal;
    }

    // adjust vertices post-loop
    for (unsigned int i = 0; i < vs.size(); i++)
      verts[i] += vs[i];

    string finish_msg;
    if (it_ctrl.is_status_check_iter()) {
      // len2() for difference value to minimize internal sqrt() calls
      max_diff2 = 0;
      for (unsigned int i = 0; i < verts.size(); i++) {
        double diff2 = (verts[i] - verts_last[i]).len2();
        if (diff2 > max_diff2)
          max_diff2 = diff2;
      }

      if (sqrt(max_diff2) < test_val) {
        completed = true;
        it_ctrl.set_finished();
        finish_msg = "solved, test value achieved";
      }
      else if (it_ctrl.is_last_iter()) {
        // reached last iteration without solving
        it_ctrl.set_finished();
        finish_msg = "not solved, test value not achieved";
      }

      // check if radius is expanding or contracting unreasonably,
      // but only for the purpose of finishing early
      // if minimum and maximum radius are differing, the polyhedron is
      // crumpling
      if (radius_range_percent &&
          canonical_radius_range_test(geom, radius_range_percent)) {
        if (!it_ctrl.is_finished())
          it_ctrl.set_finished();
        finish_msg = "breaking out: radius range detected. try increasing -d";
      }
    }

    if (it_ctrl.is_status_report_iter()) {
      if (it_ctrl.is_finished())
        it_ctrl.print("Final iteration (%s):\n", finish_msg.c_str());

      it_ctrl.print("%-12u max_diff:%17.15e\n", it_ctrl.get_current_iter(),
                    sqrt(max_diff2));
    }
  }

  return completed;
}

// P and Q are modified
void move_line_to_point(Vec3d &P, Vec3d &Q, const Vec3d &X)
{
  Vec3d Y = X + (Q - P);
  Vec3d V = P + X;
  Vec3d P2 = lines_intersection(P, V, X, Y, 0);
  if (P2.is_set()) {
    Q += P2 - P;
    P = P2;
  }
}

// RK - edge near points of base seek 1
bool canonicalize_unit(Geometry &geom, IterationControl it_ctrl,
                       const double radius_range_percent,
                       const bool planarize_only)
{
  bool completed = false;
  it_ctrl.set_finished(false);

  vector<vector<int>> edges;
  geom.get_impl_edges(edges);

  vector<Vec3d> &verts = geom.raw_verts();

  double test_val = it_ctrl.get_test_val();
  double max_diff2 = 0;

  for (it_ctrl.start_iter(); !it_ctrl.is_done(); it_ctrl.next_iter()) {
    vector<Vec3d> verts_last = verts;

    if (!planarize_only) {
      vector<Vec3d> near_pts;
      for (auto &edge : edges) {
        // unit near point
        Vec3d P = geom.edge_nearpt(edge, Vec3d::zero).unit();
        near_pts.push_back(P);
        move_line_to_point(verts[edge[0]], verts[edge[1]], P);
      }

      // re-center for drift
      Vec3d cent_near_pts = centroid(near_pts);
      for (unsigned int i = 0; i < verts.size(); i++)
        verts[i] -= cent_near_pts;
    }

    for (unsigned int f = 0; f < geom.faces().size(); f++) {
      /*
            // give polygon its own geom. face index needs to reside in vector
            vector<int> face_idxs(1);
            face_idxs[0] = f;
            Geometry polygon = faces_to_geom(geom, face_idxs);
            plane_face(polygon);

            // map vertices back into original geom
            // the numerical order of vertex list in polygon geom is preserved
            vector<int> v_idx;
            for (int v : geom.faces(f))
              v_idx.push_back(v);
            sort(v_idx.begin(), v_idx.end());
            int j = 0;
            for (int v : v_idx)
              verts[v] = polygon.verts(j++);
      */
      // RK - this does formulaically what the above does by brute force
      Vec3d face_normal = face_norm(geom.verts(), geom.faces(f)).unit();
      Vec3d face_centroid = geom.face_cent(f);
      // make sure face_normal points outward
      if (vdot(face_normal, face_centroid) < 0)
        face_normal *= -1.0;
      // place a planar vertex over or under verts[v]
      // adds or subtracts it to get to the planar verts[v]
      for (int v : geom.faces(f))
        verts[v] += vdot(face_normal, face_centroid - verts[v]) * face_normal;
    }

    string finish_msg;
    if (it_ctrl.is_status_check_iter()) {
      // len2() for difference value to minimize internal sqrt() calls
      max_diff2 = 0;
      for (unsigned int i = 0; i < verts.size(); i++) {
        double diff2 = (verts[i] - verts_last[i]).len2();
        if (diff2 > max_diff2)
          max_diff2 = diff2;
      }

      if (sqrt(max_diff2) < test_val) {
        completed = true;
        it_ctrl.set_finished();
        finish_msg = "solved, test value achieved";
      }
      else if (it_ctrl.is_last_iter()) {
        // reached last iteration without solving
        it_ctrl.set_finished();
        finish_msg = "not solved, test value not achieved";
      }

      // check if radius is expanding or contracting unreasonably,
      // but only for the purpose of finishing early
      // if minimum and maximum radius are differing, the polyhedron is
      // crumpling
      if (radius_range_percent &&
          canonical_radius_range_test(geom, radius_range_percent)) {
        if (!it_ctrl.is_finished())
          it_ctrl.set_finished();
        finish_msg = "breaking out: radius range detected. try increasing -d";
      }
    }

    if (it_ctrl.is_status_report_iter()) {
      if (it_ctrl.is_finished())
        it_ctrl.print("Final iteration (%s):\n", finish_msg.c_str());

      it_ctrl.print("%-12u max_diff:%17.15e\n", it_ctrl.get_current_iter(),
                    sqrt(max_diff2));
    }
  }

  return completed;
}

// RK - wrapper for basic planarization with unit algorithm
// meant to be called with finite num_iters (not -1)
bool planarize_unit(Geometry &geom, IterationControl it_ctrl)
{
  double radius_range_percent = 0;
  bool planarize_only = true;
  return canonicalize_unit(geom, it_ctrl, radius_range_percent, planarize_only);
}

void realign_output(Geometry &base, const cn_opts &opts)
{
  Symmetry sym(base);
  opts.message(msg_str("the symmetry of output model is %s (realigned)\n",
                       sym.get_symbol().c_str()));
  base.transform(sym.get_to_std());
}

int main(int argc, char *argv[])
{
  cn_opts opts;
  opts.process_command_line(argc, argv);

  Geometry base;
  opts.read_or_error(base, opts.ifile);

  // the program can take many non-convex input, so no need
  // if (!check_convexity(base))
  //  opts.warning("input model may not be convex");

  // epsilon setting for checking edge intersections
  double epsilon_local = opts.eps * 100;
  // don't let epsilon fall below 1e-11
  if (epsilon_local > anti::epsilon)
    epsilon_local = anti::epsilon * 10;
  if (opts.planarize_method || opts.canonical_method != 'x')
    opts.message(
        msg_str("analyzing edge intersections at a fixed epsilon of %.0e",
                epsilon_local));

  if (opts.shuffle_model_idxs.length())
    shuffle_model_indexes(base, opts);

  Symmetry sym;
  if (opts.use_symmetry) {
    opts.print_status_or_exit(sym.init(base), 'y');
    // find the nearest point to the origin in the subspace fixed by sym
    // the distance of this point from the origin should be within the
    // precision limit, or the canonical algorithm may not complete
    auto near_pt = sym.get_fixed_subspace().nearest_point(Vec3d::zero);
    if (near_pt.len() > opts.it_ctrl.get_test_val())
      opts.error("polyhedron not centred on origin (try aligning first with "
                 "off_trans -y full)",
                 'y');
  }

  if (opts.target_model != 'b') {
    opts.message("converting target to dual for processing", "t");
    base = get_dual(base);
  }

  bool perfect_score = true;
  bool completed = false;
  if (opts.planarize_method) {
    string planarize_str;
    if (opts.planarize_method == 'q')
      planarize_str = "face centroids magnitude squared";
    else if (opts.planarize_method == 'm')
      planarize_str = "mathematica";
    else if (opts.planarize_method == 'a')
      planarize_str = "sand and fill";
    else if (opts.planarize_method == 'p')
      planarize_str = "poly_form -a p";
    opts.message(msg_str("planarize: %s method\n", planarize_str.c_str()), "p");

    bool planarize_only = true;
    opts.it_ctrl.set_max_iters(opts.num_iters_planar);
    double radius_range_pct =
        0; // RK: not used when planarizing
           // (opts.radius_range_percent < 0) ? 0 : opts.radius_range_percent;

    if (opts.planarize_method == 'b') {
      completed = canonicalize_bd(base, opts.it_ctrl, radius_range_pct / 100,
                                  planarize_only);
    }
    else if (opts.planarize_method == 'm') {
      completed = canonicalize_mm(base, opts.it_ctrl, opts.edge_factor / 100,
                                  opts.plane_factor / 100,
                                  radius_range_pct / 100, planarize_only);
    }
    else if (opts.planarize_method == 'a') {
      completed = canonicalize_unit(base, opts.it_ctrl, radius_range_pct / 100,
                                    planarize_only);
    }
    else if (opts.planarize_method == 'p') {
      Status stat;
      stat = make_planar(base, opts.it_ctrl, opts.plane_factor / 100, sym);
      completed = stat.is_ok(); // true if completed;
    }

    fprintf(stderr, "\n");
    opts.message(msg_str("the plarization algorithm %s",
                         (completed ? "completed" : "did not complete")),
                 "p");

    if (opts.target_model == 'p') {
      opts.message("converting target back to base after planarization", "t");
      base = get_dual(base);
    }

    // report planarity, report convex hull test only if canonicalizing
    planarity_info(base, opts);
    string s = "no ";
    if (get_convexity(base) != Convexity::convex_strict) {
      s = "";
      perfect_score = false;
    }
    opts.message(msg_str("%sextra faces produced by convex hull", s.c_str()));
  }

  // check if model is already canonical
  bool input_is_canonical = false;
  if ((opts.canonical_method != 'x') && (opts.num_iters_canonical != 0)) {
    perfect_score = false;
    input_is_canonical = precheck(base, epsilon_local);
    if (input_is_canonical) {
      opts.warning(msg_str("INPUT MODEL IS ALREADY CANONICAL at the at -l %d\n",
                           opts.it_ctrl.get_sig_digits()),
                   'l');
      perfect_score = true;
      // cancel canonicalization
      opts.canonical_method = 'x';
      completed = true;
    }
  }

  if ((opts.canonical_method != 'x') && (opts.num_iters_canonical != 0)) {
    if (opts.edge_distribution) {
      opts.message("edge distribution: project onto sphere", "e");
      if (opts.edge_distribution == 's')
        project_onto_sphere(base);
    }

    // make it possible to just measure input
    if (opts.num_iters_canonical != 0) {
      opts.message("centering: edge near points centroid moved to origin");
      base.transform(
          Trans3d::translate(-edge_nearpoints_centroid(base, Vec3d::zero)));

      opts.message("starting radius: average edge near points to 1");
      unitize_nearpoints_radius(base);
    }

    string canonicalize_str;
    if (opts.canonical_method == 'm')
      canonicalize_str = "mathematica";
    else if (opts.canonical_method == 'c')
      canonicalize_str = "circle packings";
    else if (opts.canonical_method == 'b')
      canonicalize_str = "base/dual";
    else if (opts.canonical_method == 'a')
      canonicalize_str = "moving edge";
    opts.message(msg_str("canonicalize: %s method\n", canonicalize_str.c_str()),
                 "c");

    bool planarize_only = false;
    opts.it_ctrl.set_max_iters(opts.num_iters_canonical);
    double radius_range_pct =
        (opts.radius_range_percent < 0) ? 80 : opts.radius_range_percent;

    if (opts.canonical_method == 'm') {
      completed = canonicalize_mm(base, opts.it_ctrl, opts.edge_factor / 100,
                                  opts.plane_factor / 100,
                                  radius_range_pct / 100, planarize_only);
    }
    else if (opts.canonical_method == 'c') {
      Status stat =
          make_canonical(base, opts.it_ctrl, opts.factor / 100,
                         opts.factor_max / 100, opts.initial_point_type, sym);
      if (stat.is_error())
        opts.print_status_or_exit(stat);

      completed = stat.is_ok(); // ok - completed, warning - not completed
    }
    else if (opts.canonical_method == 'b') {
      completed = canonicalize_bd(base, opts.it_ctrl, radius_range_pct / 100,
                                  planarize_only);
    }
    else if (opts.canonical_method == 'a') {
      completed = canonicalize_unit(base, opts.it_ctrl, radius_range_pct / 100,
                                    planarize_only);
    }
  }

  if ((opts.canonical_method != 'x') && (opts.num_iters_canonical != 0)) {
    fprintf(stderr, "\n");
    opts.message(msg_str("the canonical algorithm %s",
                         (completed ? "completed" : "did not complete")),
                 "c");
  }

  if (opts.target_model == 'c') {
    opts.message("converting target back to base after canonicalization", "t");
    base = get_dual(base);
  }

  // standardize model radius if needed since dual is reciprocated on 1
  if ((opts.canonical_method != 'x') && (opts.num_iters_canonical != 0))
    reset_model_size(base, epsilon_local, opts);

  // generate dual once
  Geometry dual = get_dual(base);

  // coincidence check for simulataneous vertices and faces
  check_coincidence(base, dual, opts);

  vector<Vec3d> base_nearpts;
  vector<Vec3d> dual_nearpts;
  vector<Vec3d> ips;
  generate_points(base, dual, base_nearpts, dual_nearpts, ips, epsilon_local);

  if (opts.canonical_method != 'x' || input_is_canonical)
    perfect_score =
        canonical_report(base, base_nearpts, ips, epsilon_local, opts);

  if (opts.canonical_method != 'x' || input_is_canonical ||
      opts.planarize_method) {
    if (perfect_score)
      opts.message(msg_str("model PASSED ALL TESTS at -l %d",
                           opts.it_ctrl.get_sig_digits()));
    else {
      opts.message("model could be improved by raising -l or raising -n");
      if (opts.it_ctrl.get_sig_digits() >= 15)
        opts.warning("if -l is raised higher than 15, algorithms "
                     "may not halt");
    }
  }

  if (opts.realign)
    realign_output(base, opts);

  // parts to output
  construct_model(base, dual, base_nearpts, dual_nearpts, ips, opts);

  opts.write_or_error(base, opts.ofile);

  return 0;
}
