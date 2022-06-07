/*
   Copyright (c) 2003-2021, Adrian Rossiter, Roger Kaufman
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

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using std::string;
using std::vector;

using namespace anti;

class cn_opts : public ProgramOpts {
public:
  IterationControl it_ctrl;

  string ifile;
  string ofile;

  char edge_distribution = '\0';    // can project onto sphere
  string shuffle_model_indexes;     // shuffle indexes can be vef
  char target_model = 'b';          // work on base model is default
  char planarize_method = '\0';     // no algorithm for planar is default
  int num_iters_planar = -1;        // unlimited iterations for planar
  char canonical_method = 'c';      // circle packings algorithm is default
  int num_iters_canonical = -1;     // unlimited iterations for canonical
  bool use_symmetry = false;        // don't use symmetry
  double edge_factor = NAN;         // mathematica algorithm variable
  double plane_factor = NAN;        // mathematica algorithm variable
  double radius_range_percent = -1; // percent expansion of model
  double factor = NAN;              // initial adjustment factor
  double factor_max = NAN;          // maximum adjustment factor
  char initial_point_type = 'c';    // c - edge centroids, n - edge near points
  string output_parts = "b";        // parts of output model
  int face_opacity = -1;            // transparency
  double offset = 0;                // offset of tangency sphere
  int roundness = 8;                // roundness of tangency sphere

  double eps = anti::epsilon;

  Color ipoints_col = Color(255, 255, 0);    // yellow
  Color base_nearpts_col = Color(255, 0, 0); // red
  Color dual_nearpts_col = Color(0, 100, 0); // darkgreen
  char base_incircles_color_method = 'f';    // take face colors is default
  Color base_incircles_col = Color();
  char dual_incircles_color_method = 'f'; // take face colors is default
  Color dual_incircles_col = Color();
  char dual_face_color_method = 'b'; // take colors from base vertices
  Color dual_face_col = Color();
  Color base_edge_col = Color();
  Color dual_edge_col = Color();
  Color sphere_col = Color(255, 255, 255); // white

  cn_opts() : ProgramOpts("canonical")
  {
    it_ctrl.set_max_iters(-1);
    it_ctrl.set_status_checks("1000,1");
    it_ctrl.set_sig_digits(int(-log(anti::epsilon) / log(10) + 0.5));
  }

  void process_command_line(int argc, char **argv);
  void usage();
};

void extended_help()
{
  fprintf(stdout, R"(
Calculating Canonical Polyhedra

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

)");
}

void cn_opts::usage()
{
  fprintf(stdout, R"(
Usage: %s [options] [input_file]

Read a polyhedron from a file in OFF format. Canonicalize or planarize it.
Uses algorithms by George W. Hart, http://www.georgehart.com/
http://www.georgehart.com/virtual-polyhedra/conway_notation.html
http://www.georgehart.com/virtual-polyhedra/canonical.html
If input_file is not given the program reads from standard input.

Options
%s
  -H        documention on algorithm 
  -e <opt>  edge distribution (default : none)
               s - project vertices onto a sphere
  -s <opt>  shuffle model indexes
               v - vertices, e - edges, f - faces, a - all (default: none)
  -t <opt>  target model
               b - work on base only (default)
               p - work on dual for planarization only
               c - work on dual for planarization and canonicalization
  -p <opt>  planarization (done before canoncalization. default: none)
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
  -y        maintain symmetry of the base model (-p p, -c c)
  -O <args> output b - base, d - dual, i - intersection points (default: b)
               edge nearpoints, n - base, m - dual; C - base/dual convex hull
               edge nearpoints centroid, p - base, q - dual; o - origin point
               tangent sphere, u - minimum, U - maximum
               incircles, s - base, t - dual; as rings, S - base, T - dual
  -q <dist> offset for incircles to avoid coplanarity e.g 0.0001 (default: 0)
  -g <opt>  roundness of tangent sphere, positive integer n (default: 8)
  -d <perc> radius test. percent difference between minimum and maximum radius
               checks if polyhedron is collapsing. 0 for no test 
               (default: 80 for canonicalizing, not used for planarizing)
  -l <lim>  minimum distance change to terminate, as negative exponent
               (default: %d giving %.0e)
            WARNING: high values can cause non-terminal behaviour. Use -n
  -z <nums> number of iterations between status reports (implies termination
            check) (0 for final report only, -1 for no report), optionally
            followed by a comma and the number of iterations between
            termination checks (0 for report checks only) (default: %d,%d)
  -o <file> write output to file (default: write to standard output)

Extra Options
  -E <perc> percentage to scale edge tangency (default: 50) (-c m)
  -P <perc> percentage to scale face planarity (default: 20) (-c m, -p m, -p p)
  -f <adj>  initial percent adjustment factor, optionally followed by a comma
            and a maximum percent adjustment (default: 1,50) (-c c)
  -C        continue processing a near-canonical model (the initial
            intermediate processing model will preserves the geometry
            of the base model rather than avoid scrambling) (-c c)

Coloring Options (run 'off_util -H color' for help on color formats)
  -I <col>  intersection points and/or origin color (default: yellow)
  -N <col>  base near points, centroid color (default: red)
  -M <col>  dual near points, centroid color (default: darkgreen)
  -S <col>  base incircles color. keyword: f take color of face (default)
  -R <col>  dual incircles color. keyword: f take color of face (default)
  -D <col>  dual face color. keyword: b take color from base vertices (default)
  -J <col>  base edge color (default: unchanged)
  -K <col>  dual edge color (default: unchanged)
  -U <col>  unit sphere and/or convex hull color (default: white)
  -T <tran> base/dual transparency. range from 0 (invisible) to 255 (opaque)

)",
          prog_name(), help_ver_text, it_ctrl.get_max_iters(),
          it_ctrl.get_max_iters(), it_ctrl.get_sig_digits(),
          it_ctrl.get_test_val(), it_ctrl.get_status_check_and_report_iters(),
          it_ctrl.get_status_check_only_iters());
}

void cn_opts::process_command_line(int argc, char **argv)
{
  opterr = 0;
  int c;
  int num;

  bool p_set = false;
  bool c_set = false;

  handle_long_opts(argc, argv);

  while ((c = getopt(
              argc, argv,
              ":hHe:s:t:p:i:c:n:yO:q:g:E:P:F:Cd:z:I:N:M:S:R:D:J:K:U:T:l:o:")) !=
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
        error(msg_str("suffle parts are '%s' must be any or all from "
                      "v, e, f, a",
                      optarg),
              c);
      shuffle_model_indexes = optarg;
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

    case 'O':
      if (strspn(optarg, "bdinmopqsStTuUC") != strlen(optarg))
        error(msg_str("output parts are '%s' must be any or all from "
                      "b, d, i, n, m, o, p, q, s, S, t, T, u, U, C",
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

    case 'E':
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

    case 'F': {
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

    case 'I':
      print_status_or_exit(ipoints_col.read(optarg), c);
      break;

    case 'N':
      print_status_or_exit(base_nearpts_col.read(optarg), c);
      break;

    case 'M':
      print_status_or_exit(dual_nearpts_col.read(optarg), c);
      break;

    case 'S':
      base_incircles_color_method = '\0';
      if (strchr("f", *optarg))
        base_incircles_color_method = *optarg;
      else
        print_status_or_exit(base_incircles_col.read(optarg), c);
      break;

    case 'R':
      dual_incircles_color_method = '\0';
      if (strchr("f", *optarg))
        dual_incircles_color_method = *optarg;
      else
        print_status_or_exit(dual_incircles_col.read(optarg), c);
      break;

    case 'D':
      dual_face_color_method = '\0';
      if (strchr("b", *optarg))
        dual_face_color_method = *optarg;
      else
        print_status_or_exit(dual_face_col.read(optarg), c);
      break;

    case 'J':
      print_status_or_exit(base_edge_col.read(optarg), c);
      break;

    case 'K':
      print_status_or_exit(dual_edge_col.read(optarg), c);
      break;

    case 'U':
      print_status_or_exit(sphere_col.read(optarg), c);
      break;

    case 'T':
      print_status_or_exit(read_int(optarg, &face_opacity), c);
      if (face_opacity < 0 || face_opacity > 255) {
        error("face transparency must be between 0 and 255", c);
      }
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
    if (planarize_method &&
        (planarize_method != 'p' && planarize_method != 'e'))
      warning("set, but not used for this planarize algorithm", 'y');
    if (canonical_method != 'c')
      warning("set, but not used for this canonical algorithm", 'y');
  }

  if (argc - optind > 1)
    error("too many arguments");

  if (argc - optind == 1)
    ifile = argv[optind];

  eps = it_ctrl.get_test_val();
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
    // Consider each edge of the face to make a triangle with the face centroid
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
  //return (check_congruence(geom, hull));
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

void planarity_info(Geometry &geom)
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
  fprintf(stderr, "\n");
  fprintf(stderr, "maximum nonplanarity = %.0e\n", max_nonplanar);
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
    base_nearpts.push_back(base.edge_nearpt(base_edge, Vec3d(0, 0, 0)));

  for (auto &dual_edge : dual_edges)
    dual_nearpts.push_back(dual.edge_nearpt(dual_edge, Vec3d(0, 0, 0)));

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
    Vec3d P = geom.edge_nearpt(edge, Vec3d(0, 0, 0));
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
                        string str, const double &epsilon_local)
{
  double min = 0;
  double max = 0;
  Vec3d center;
  edge_nearpoints_error(geom, min, max, center);
  fprintf(stderr, "%s range of edge nearpoint error is %.0e to %.0e\n",
          str.c_str(), min, max);

  int radius_count = nearpts.size();
  int nearpts_size = nearpts.size();
  for (int i = 0; i < nearpts_size; i++) {
    double l = fabs(nearpts[i].len() - 1.0);
    if (double_ne(l, 0.0, epsilon_local))
      radius_count--;
  }

  double np_pct = radius_count / (double)nearpts_size * 100;
  fprintf(stderr,
          "%d out of %d %s edge nearpoint radii lie within %.0e of length 1 "
          "(%g%%)\n",
          radius_count, nearpts_size, str.c_str(), epsilon_local, np_pct);

  double max_error = std::numeric_limits<int>::min();
  for (unsigned int i = 0; i < 3; i++) {
    double err = fabs(center[i]);
    if (err > max_error)
      max_error = err;
  }

  if (max_error > max)
    max = max_error;

  string lstr = "lies ";
  if (double_ne(max_error, 0.0, epsilon_local))
    lstr = "does not lie ";
  fprintf(
      stderr,
      "%s nearpoint centroid (%.0e, %.0e, %.0e) %swithin %.0e of the origin\n",
      str.c_str(), center[0], center[1], center[2], lstr.c_str(),
      epsilon_local);

  return max;
}

void canonical_report(const Geometry &base, const Geometry &dual,
                      vector<Vec3d> &base_nearpts, vector<Vec3d> &dual_nearpts,
                      vector<Vec3d> &ips, const bool completed,
                      const double &epsilon_local)
{
  double max_error = std::numeric_limits<int>::min();
  double err = 0;

  fprintf(stderr, "the canonical algorithm %s\n",
          (completed ? "completed" : "did not complete"));

  // Symmetry sym(base);
  // fprintf(stderr, "the symmetry of base model is %s\n",
  //        sym.get_symbol().c_str());

  // intersection point count is not used in score
  int ip_size = ips.size();
  int bn_size = base_nearpts.size();
  double pct = ip_size / (double)bn_size * 100;
  fprintf(stderr,
          "%d out of %d base/dual edge intersection points found (%g%%)\n",
          ip_size, bn_size, pct);

  err = nearpoint_report(base, base_nearpts, "base", epsilon_local);
  if (err > max_error)
    max_error = err;

  nearpoint_report(dual, dual_nearpts, "dual", epsilon_local);
  // don't consider dual statistics in score
  // if (err > max_error)
  //  max_error = err;

  fprintf(stderr, "base canonical model maximum error is %.0e\n", max_error);

  string str;
  if (is_nonoverlap_single_cover(base))
    str = "no ";
  fprintf(stderr, "base canonical model has %soverlap error\n", str.c_str());

  // fprintf(stderr, "note the midcenter is the origin\n");
  fprintf(stderr, "\n");
}

/*
Vec3d face_edge_nearpoints_centroid(const Geometry &geom, double &radius, const
int face_no)
{
  vector<Vec3d> near_pts;

  unsigned int fsz = geom.faces(face_no).size();
  for (unsigned int i = 0; i < fsz; i++) {
    int v1 = geom.faces(face_no)[i];
    int v2 = geom.faces(face_no)[(i + 1) % fsz];

    Vec3d P = geom.edge_nearpt(make_edge(v1,v2), Vec3d(0, 0, 0));
    near_pts.push_back(P);
  }

  Vec3d face_edge_nearpt_centroid = centroid(near_pts);
//  face_edge_nearpt_centroid = geom.face_cent(face_no);

  // get the minimum radius
  double min = DBL_MAX;
  for (unsigned int i = 0; i < fsz; i++) {
    double l = (face_edge_nearpt_centroid - near_pts[i]).len();
    if (l < min)
      min = l;
  }
  radius = min;

  return face_edge_nearpt_centroid;
}

vector<Vec3d> face_edge_nearpoints_centroids(const Geometry &base)
{
  double rad;
  Geometry near_pts;
  unsigned int fsz = base.faces().size();
  for (unsigned int i = 0; i < fsz; i++)
    near_pts.add_vert(face_edge_nearpoints_centroid(base, rad, i));

  return(near_pts.verts());
}
*/

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

    Vec3d P = geom.edge_nearpt(make_edge(v1, v2), Vec3d(0, 0, 0));
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

Geometry incircles(const Geometry &geom, const char &incircle_color_method,
                   const Color &incircle_color, bool filled, double offset)
{
  Geometry incircles;
  Geometry circle = unit_circle(60, incircle_color, filled);

  for (unsigned int i = 0; i < geom.faces().size(); i++) {
    // find incircle rotation place
    Vec3d face_centroid = anti::centroid(geom.verts(), geom.faces(i));
    Vec3d face_normal = face_norm(geom.verts(), geom.faces(i)).unit();
    Vec3d center = face_normal * vdot(face_centroid, face_normal);
    Color col = geom.colors(FACES).get(i);

    // find radius of incircle, and make incircle of radius
    Geometry incircle = circle;
    if (incircle_color_method == 'f') {
      Coloring(&incircle).e_one_col(col);
      Coloring(&incircle).f_one_col(col);
    }

    double radius = incircle_radius(geom, center, i);
    incircle.transform(Trans3d::scale(radius));

    // set depth of incircle
    double depth = center.len() + offset;
    incircle.transform(Trans3d::translate(Vec3d(0, 0, depth)));

    // rotate incircle into place
    incircle.transform(Trans3d::rotate(Vec3d(0, 0, 1), center));
    incircles.append(incircle);
  }
  return incircles;
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
void reset_model_size(Geometry &geom, const double &epsilon_local)
{
  double radius = edge_nearpoints_radius(geom);
  if (double_ne(radius, 1.0, epsilon_local)) {
    fprintf(stderr,
            "Resetting base nearpoints average radius to within %.0e of 1\n",
            epsilon_local);
    unitize_nearpoints_radius(geom);
  }
}

void set_edge_colors(Geometry &geom, const Color col)
{
  // it is possible unset faces already exist
  if (col.is_set()) {
    geom.add_missing_impl_edges();
    Coloring(&geom).e_one_col(col);
  }
}

void construct_model(Geometry &base, Geometry &dual,
                     vector<Vec3d> &base_nearpts, vector<Vec3d> &dual_nearpts,
                     vector<Vec3d> &ips, const cn_opts &opts)
{
  // get statistics before model is changed, if needed
  double min = 0;
  double max = 0;
  Vec3d center;
  if (opts.output_parts.find_first_of("uU") != string::npos)
    edge_nearpoints_radius(base, min, max, center);

  if (opts.dual_face_color_method != 'b')
    Coloring(&dual).f_one_col(opts.dual_face_col);

  // base incircles
  Geometry base_incircles;
  if (opts.output_parts.find_first_of("sS") != string::npos) {
    bool filled = (opts.output_parts.find("s") != string::npos) ? true : false;
    base_incircles = incircles(base, opts.base_incircles_color_method,
                               opts.base_incircles_col, filled, opts.offset);
  }

  // dual incircles
  Geometry dual_incircles;
  if (opts.output_parts.find_first_of("tT") != string::npos) {
    bool filled = (opts.output_parts.find("t") != string::npos) ? true : false;
    dual_incircles = incircles(dual, opts.dual_incircles_color_method,
                               opts.dual_incircles_col, filled, opts.offset);
  }

  if (opts.output_parts.find_first_of("bC") != string::npos)
    // set edge colors here
    set_edge_colors(base, opts.base_edge_col);
  else
    // clear base if not using
    base.clear_all();

  // append dual
  if (opts.output_parts.find_first_of("dC") != string::npos) {
    // set edge colors here
    set_edge_colors(dual, opts.dual_edge_col);
    base.append(dual);
  }

  // add intersection points
  if (opts.output_parts.find("i") != string::npos) {
    for (auto &ip : ips)
      base.add_vert(ip, opts.ipoints_col);
  }

  // add base near points
  if (opts.output_parts.find("n") != string::npos) {
    for (auto &base_nearpt : base_nearpts)
      base.add_vert(base_nearpt, opts.base_nearpts_col);
  }

  // add dual near points
  if (opts.output_parts.find("m") != string::npos) {
    for (auto &dual_nearpt : dual_nearpts)
      base.add_vert(dual_nearpt, opts.dual_nearpts_col);
  }

  // add base near points centroid
  if (opts.output_parts.find("p") != string::npos) {
    base.add_vert(centroid(base_nearpts), opts.base_nearpts_col);
  }

  // add dual near points centroid
  if (opts.output_parts.find("q") != string::npos) {
    base.add_vert(centroid(dual_nearpts), opts.dual_nearpts_col);
  }

  // add origin point
  if (opts.output_parts.find("o") != string::npos) {
    base.add_vert(Vec3d(0, 0, 0), opts.ipoints_col);
  }

  // apply transparency
  Status stat = Coloring(&base).apply_transparency(opts.face_opacity);
  if (stat.is_warning())
    opts.warning(stat.msg(), 'T');

  // add unit sphere on origin
  if (opts.output_parts.find_first_of("uU") != string::npos) {
    Geometry sgeom;
    string geo_str = "geo_" + std::to_string(opts.roundness) + "_" +
                     std::to_string(opts.roundness);
    sgeom.read_resource(geo_str);
    sgeom.transform(Trans3d::translate(-centroid(sgeom.verts())));
    unitize_vertex_radius(sgeom);

    if (opts.output_parts.find("u") != string::npos)
      sgeom.transform(Trans3d::scale(min));
    else
      sgeom.transform(Trans3d::scale(max));

    Coloring(&sgeom).vef_one_col(Color::invisible, Color::invisible,
                                 opts.sphere_col);
    base.append(sgeom);
  }

  // add convex hull
  if (opts.output_parts.find("C") != string::npos) {
    Geometry base_dual = base;
    base_dual.append(dual);
    base_dual.set_hull();
    Coloring(&base_dual).f_one_col(opts.sphere_col);
    base.append(base_dual);
  }

  if (base_incircles.verts().size())
    base.append(base_incircles);

  if (dual_incircles.verts().size())
    base.append(dual_incircles);
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
  bool shuffle_verts =
      (opts.shuffle_model_indexes.find_first_of("va") != string::npos);
  bool shuffle_faces =
      (opts.shuffle_model_indexes.find_first_of("fa") != string::npos);
  bool shuffle_edges =
      (opts.shuffle_model_indexes.find_first_of("ea") != string::npos);

  map<int, int> new_verts;
  if (!shuffle_verts) {
    for (unsigned int i = 0; i < geom.verts().size(); i++)
      new_verts[i] = i;
  }
  else {
    fprintf(stderr, "shuffle model indexes: vertices\n");
    vector<int> deal = geom_deal(geom, geom.verts().size());
    map<int, Color> new_cols;
    vector<Vec3d> shuffled_verts;
    for (unsigned int i = 0; i < geom.verts().size(); i++) {
      int v_new = deal[i];
      new_verts[v_new] = i;
      shuffled_verts.push_back(geom.verts(v_new));
      new_cols[i] = geom.colors(VERTS).get(v_new);
    }
    geom.raw_verts() = shuffled_verts;
    for (unsigned int i = 0; i < geom.verts().size(); i++)
      geom.colors(VERTS).set(i, new_cols[i]);
  }

  if (shuffle_faces || shuffle_verts) {
    map<int, int> face_order;
    if (!shuffle_faces) {
      for (unsigned int i = 0; i < geom.faces().size(); i++)
        face_order[i] = i;
    }
    else {
      fprintf(stderr, "shuffle model indexes: faces\n");
      vector<int> deal = geom_deal(geom, geom.faces().size());
      for (unsigned int i = 0; i < geom.faces().size(); i++)
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
    geom.raw_faces() = shuffled_faces;
    for (unsigned int i = 0; i < geom.faces().size(); i++)
      geom.colors(FACES).set(i, new_cols[i]);
  }

  if (!geom.edges().size()) {
    fprintf(stderr, "shuffle model indexes: no explicit edges\n");
  }
  else if (shuffle_edges || shuffle_verts) {
    map<int, int> edge_order;
    if (!shuffle_edges)
      for (unsigned int i = 0; i < geom.edges().size(); i++) {
        edge_order[i] = i;
      }
    else {
      fprintf(stderr, "shuffle model indexes: edges\n");
      vector<int> deal = geom_deal(geom, geom.edges().size());
      for (unsigned int i = 0; i < geom.edges().size(); i++)
        edge_order[deal[i]] = i;
    }
    map<int, Color> new_cols;
    vector<vector<int>> shuffled_edges;
    for (unsigned int i = 0; i < geom.edges().size(); i++) {
      vector<int> edge;
      for (unsigned int j = 0; j < geom.edges(edge_order[i]).size(); j++) {
        int v = geom.edges(edge_order[i])[j];
        int v_new = new_verts[v];
        edge.push_back(v_new);
      }
      shuffled_edges.push_back(edge);
      new_cols[i] = geom.colors(EDGES).get(edge_order[i]);
    }
    geom.raw_edges() = shuffled_edges;
    for (unsigned int i = 0; i < geom.edges().size(); i++)
      geom.colors(EDGES).set(i, new_cols[i]);
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
        Vec3d P = geom.edge_nearpt(edge, Vec3d(0, 0, 0));
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
      v = Vec3d(0, 0, 0);

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

Geometry base_to_ambo(const Geometry &base, char point_type)
{
  Geometry ambo = base;
  truncate_verts(ambo, 0.5);
  if (point_type == 'n') {
    vector<vector<int>> impl_edges;
    base.get_impl_edges(impl_edges);
    int i = 0;
    for (const auto &edge : impl_edges)
      ambo.verts(i++) = nearest_point(Vec3d(0, 0, 0), base.verts(), edge);
  }
  ambo.transform(Trans3d::translate(-ambo.centroid()));
  return ambo;
}

void update_base_from_ambo(Geometry &base, const Geometry &ambo)
{
  Geometry ambo_d;
  get_pol_recip_verts(ambo_d, ambo, 1, Vec3d::zero);
  auto v_sz = base.verts().size();
  base.raw_verts().clear();
  std::copy(ambo_d.verts().begin(), ambo_d.verts().begin() + v_sz,
            std::back_inserter(base.raw_verts()));
  return;
}

Status make_planar_unit(Geometry &base_geom, IterationControl it_ctrl,
                        double factor, double factor_max, Symmetry sym)
{
  // chosen by experiment
  const double readjust_up = 1.01;    // to adjust adjustment factor up
  const double readjust_down = 0.995; // to adjust adjustment factor down
  const double unit_mult = 1;         // factor for unit adjustment
  const double orth_mult = 0.5;       // extra multiplier for orthogonality adj
  const double overlap_mult = 0.5;    // extra multiplier for overlap adj

  Status stat;

  base_geom.orient(1); // positive orientation

  // transform to centre on centroid, adjust to_std to match
  auto initial_centroid = base_geom.centroid();
  base_geom.transform(Trans3d::translate(base_geom.centroid()));
  sym.set_to_std(sym.get_to_std() * Trans3d::translate(initial_centroid));

  bool using_symmetry = (sym.get_sym_type() > Symmetry::C1);
  auto fixed_subspace = sym.get_fixed_subspace(); // used to find centroid
  SymmetricUpdater sym_updater((using_symmetry) ? base_geom : Geometry(), sym);
  const Geometry &geom =
      (using_symmetry) ? sym_updater.get_geom_working() : base_geom;

  // No further processing if no faces, but not an error
  if (geom.faces().size() == 0)
    return stat;

  const vector<Vec3d> &verts = geom.verts();
  const vector<vector<int>> &faces = geom.faces();

  // List of faces that a vertex is part of
  auto vert_faces = geom.get_info().get_dual().faces();
  vector<vector<int>> vert_figs;
  {
    auto vfigs = geom.get_info().get_vert_figs();
    for (const auto &vfig : vfigs) {
      if (vfig[0].size() != 4) {
        stat.set_error(msg_str(
            "intermediate (ambo) model has vertex with order %d instead of 4",
            (int)vfig[0].size()));
        return stat;
      }
      vert_figs.push_back(vfig[0]);
    }
  }

  // Get a list of the faces that contain a principal vertex of any type
  vector<int> principal_verts;  // first vertex in each vertex orbit
  vector<int> verts_to_update;  // vertices accessed during iteration
  vector<int> faces_to_process; // faces processed during iteration
  if (using_symmetry) {
    principal_verts = sym_updater.get_principal(VERTS);
    faces_to_process = sym_updater.get_associated_elems(vert_faces);
    verts_to_update =
        SymmetricUpdater::get_included_verts(faces_to_process, geom.faces());
    auto vfig_verts_to_update =
        SymmetricUpdater::get_included_verts(principal_verts, vert_figs);
    verts_to_update.insert(verts_to_update.end(), vfig_verts_to_update.begin(),
                           vfig_verts_to_update.end());
    SymmetricUpdater::to_unique_index_list(verts_to_update);
  }
  else { // not using_symmetry
    principal_verts =
        SymmetricUpdater::sequential_index_list(geom.verts().size());
    faces_to_process =
        SymmetricUpdater::sequential_index_list(geom.faces().size());
  }

  // Use oversized arrays to avoid mapping
  vector<Vec3d> offsets(verts.size()); // Vertex adjustments
  vector<Vec3d> norms(faces.size());   // Face normals
  vector<Vec3d> cents(faces.size());   // Face centroids

  double test_val = it_ctrl.get_test_val();
  double last_max_diff2 = 0.0;
  for (it_ctrl.start_iter(); !it_ctrl.is_done(); it_ctrl.next_iter()) {
    std::fill(offsets.begin(), offsets.end(), Vec3d::zero);

    if (using_symmetry) {
      // Ensure that the vertices used in adjustment are up to date
      for (auto v_idx : verts_to_update)
        sym_updater.update_from_principal_vertex(v_idx);
    }

    // Initialize face data for just the necessary faces
    for (auto f_idx : faces_to_process) {
      norms[f_idx] = geom.face_norm(f_idx).unit();
      cents[f_idx] = geom.face_cent(f_idx);
    }

    Vec3d centroid = Vec3d::zero;
    if (using_symmetry) {
      // For each orbit, project a vertex weighted by the orbit size onto
      // the fixed subspace
      const auto &vorbits = sym_updater.get_equiv_sets(VERTS);
      for (const auto &vorbit : vorbits)
        centroid += fixed_subspace.nearest_point(verts[*vorbit.begin()]) *
                    vorbit.size();
      centroid /= verts.size(); // centroid of weighted projected vertices
    }
    else {
      centroid = geom.centroid();
    }

    double max_diff2 = 0.0;
    for (auto v_idx : principal_verts) {
      const auto &vfaces = vert_faces[v_idx];
      const int vf_sz = vfaces.size();
      // target vertex is centroid of projection of vertex onto planes
      for (int f0 = 0; f0 < vf_sz; f0++) {
        int f0_idx = vfaces[f0];
        offsets[v_idx] +=
            nearpoint_on_plane(verts[v_idx], cents[f0_idx], norms[f0_idx]);
      }
      offsets[v_idx] = (offsets[v_idx] / vf_sz - verts[v_idx]) * factor;

      // adjust for centroid
      offsets[v_idx] -= centroid;

      // adjust for orthogonality
      for (int i = 0; i < 2; i++) {
        auto n = vcross(norms[vfaces[i + 2]], norms[vfaces[i]]).unit();
        const auto v_ideal = nearpoint_on_plane(verts[v_idx], Vec3d::zero, n);
        const auto offset = (v_ideal - verts[v_idx]) * factor * orth_mult;
        offsets[v_idx] += offset;
      }

      // adjust for non-overlap
      const auto &vfig = vert_figs[v_idx];
      for (int i = 0; i < 4; i++) {
        if (vtriple(verts[v_idx], verts[vfig[i]], verts[vfig[(i + 1) % 4]]) >
            0) {
          auto v_ideal = ::centroid(
              {verts[vfig[0]], verts[vfig[1]], verts[vfig[2]], verts[vfig[3]]});
          offsets[v_idx] += (v_ideal - verts[v_idx]) * overlap_mult;
          break;
        }
      }

      auto diff2 = offsets[v_idx].len2();
      if (diff2 > max_diff2)
        max_diff2 = diff2;
    }

    // adjust vertices post-loop
    if (using_symmetry) {
      // adjust principal vertices
      for (int v_idx : principal_verts) {
        auto new_v = verts[v_idx] + offsets[v_idx];
        double new_v_len = new_v.len();
        new_v *= 1 + (1 / new_v_len - 1) * unit_mult;
        sym_updater.update_principal_vertex(v_idx, new_v);
      }
    }
    else { // not using_symmetry
      // adjust all vertices
      for (unsigned int i = 0; i < verts.size(); i++) {
        auto new_v = verts[i] + offsets[i];
        double new_v_len = new_v.len();
        new_v *= 1 + (1 / new_v_len - 1) * unit_mult;
        base_geom.raw_verts()[i] = new_v;
      }
    }

    // adjust plane factor
    if (max_diff2 < last_max_diff2)
      factor *= readjust_up;
    else
      factor *= readjust_down;
    if (factor > factor_max)
      factor = factor_max;
    last_max_diff2 = max_diff2;

    string finish_msg;
    if (it_ctrl.is_status_check_iter()) {
      sym_updater.update_all();
      double width = BoundBox(verts).max_width();
      if (sqrt(max_diff2) / width < test_val) {
        it_ctrl.set_finished();
        finish_msg = "solved, test value achieved";
        stat.set_ok(finish_msg);
      }
      else if (it_ctrl.is_last_iter()) {
        // reached last iteration without solving
        it_ctrl.set_finished();
        finish_msg = "not solved, test value not achieved";
        stat.set_warning(finish_msg);
      }
    }

    if (it_ctrl.is_status_report_iter()) {
      if (it_ctrl.is_finished())
        it_ctrl.print("Final iteration (%s):\n", finish_msg.c_str());

      it_ctrl.print("%-12u max_diff:%17.15e  -f %-10.5f\n",
                    it_ctrl.get_current_iter(), sqrt(max_diff2), 100 * factor);
    }
  }

  if (using_symmetry)
    base_geom = sym_updater.get_geom_final();

  return stat;
}

bool make_canonical_enp(Geometry &geom, IterationControl it_ctrl, double factor,
                        double factor_max, char initial_point_type,
                        Symmetry sym)
{
  Geometry ambo = base_to_ambo(geom, initial_point_type);
  Status stat = make_planar_unit(ambo, it_ctrl, factor, factor_max, sym);
  update_base_from_ambo(geom, ambo);
  return stat.is_ok(); // true if completed;
}

int main(int argc, char *argv[])
{
  cn_opts opts;
  opts.process_command_line(argc, argv);

  Geometry base;
  opts.read_or_error(base, opts.ifile);

  if (!check_convexity(base))
    opts.warning("input model may not be convex");

  if (opts.edge_distribution) {
    fprintf(stderr, "edge distribution: project onto sphere\n");
    if (opts.edge_distribution == 's')
      project_onto_sphere(base);
  }

  if (opts.shuffle_model_indexes.length())
    shuffle_model_indexes(base, opts);

  fprintf(stderr, "\n");
  fprintf(stderr, "starting radius: average edge near points\n");
  unitize_nearpoints_radius(base);

  fprintf(stderr, "centering: edge near points centroid moved to origin\n");
  base.transform(
      Trans3d::translate(-edge_nearpoints_centroid(base, Vec3d(0, 0, 0))));

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
    fprintf(stderr, "converting target to dual for planarization\n");
    base = get_dual(base);
  }

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
    fprintf(stderr, "planarize: %s method\n", planarize_str.c_str());

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

    if (opts.target_model == 'p') {
      fprintf(stderr, "converting target back to base after planarization\n");
      base = get_dual(base);
    }

    // report planarity, report convex hull test if planarizing only
    planarity_info(base);
    if (opts.canonical_method == 'x')
      fprintf(stderr, "convex hull test: %s\n",
              (check_convexity(base)) ? "passed"
                                      : "triangulated. trying raising -l");
    fprintf(stderr, "\n");
  }

  if (opts.canonical_method != 'x') {
    string canonicalize_str;
    if (opts.canonical_method == 'm')
      canonicalize_str = "mathematica";
    else if (opts.canonical_method == 'c')
      canonicalize_str = "circle packings";
    else if (opts.canonical_method == 'b')
      canonicalize_str = "base/dual";
    else if (opts.canonical_method == 'a')
      canonicalize_str = "moving edge";
    fprintf(stderr, "canonicalize: %s method\n", canonicalize_str.c_str());

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
      completed = make_canonical_enp(base, opts.it_ctrl, opts.factor / 100,
                                     opts.factor_max / 100,
                                     opts.initial_point_type, sym);
    }
    else if (opts.canonical_method == 'b') {
      completed = canonicalize_bd(base, opts.it_ctrl, radius_range_pct / 100,
                                  planarize_only);
    }
    else if (opts.canonical_method == 'a') {
      completed = canonicalize_unit(base, opts.it_ctrl, radius_range_pct / 100,
                                    planarize_only);
    }

    if (opts.target_model == 'c') {
      fprintf(stderr,
              "converting target back to base after canonicalization\n");
      base = get_dual(base);
    }

    // report planarity
    planarity_info(base);
    fprintf(stderr, "convex hull test: %s\n",
            (check_convexity(base)) ? "passed"
                                    : "triangulated. trying raising -l");
    fprintf(stderr, "\n");
  }

  // use fixed epsilon for quality comparisons
  double epsilon_local = anti::epsilon * 10;
  if (opts.canonical_method != 'x')
    fprintf(stderr, "analyzing result at a fixed epsilon of %.0e\n",
            epsilon_local);

  // standardize model radius if needed since dual is reciprocated on 1
  reset_model_size(base, epsilon_local);

  // generate dual once
  Geometry dual = get_dual(base);

  // coincidence check for simulataneous vertices and faces
  check_coincidence(base, dual, opts);

  vector<Vec3d> base_nearpts;
  vector<Vec3d> dual_nearpts;
  vector<Vec3d> ips;
  generate_points(base, dual, base_nearpts, dual_nearpts, ips, epsilon_local);

  if (opts.canonical_method != 'x')
    canonical_report(base, dual, base_nearpts, dual_nearpts, ips, completed,
                     epsilon_local);

  // parts to output
  construct_model(base, dual, base_nearpts, dual_nearpts, ips, opts);

  opts.write_or_error(base, opts.ofile);

  return 0;
}
