/*
   Copyright (c) 2012-2020, Adrian Rossiter

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

#include "geometryutils.h"
#include "tiling.h"
#include "utils.h"

#include <cstring>
#include <map>
#include <string>
#include <vector>

using std::map;
using std::pair;
using std::string;
using std::to_string;
using std::vector;

using namespace anti;

namespace { // unnamed
// nomalize point, if integer coordinates remove any common factor
static Vec3d normalize_point(const Vec3d &Q)
{
  vector<long> coords(3);
  for (int i = 0; i < 3; i++) {
    coords[i] = std::lround(Q[i]);
    if (!double_eq(coords[i], Q[i]))
      return Vec3d(Q); // not integer coordinates
  }

  long div = gcd(coords[0], gcd(coords[1], coords[2]));
  if (div == 0) // will not happen for barycentric coordinates
    div = 1;

  return Vec3d(coords[0] / div, coords[1] / div, coords[2] / div);
}

static string coord_string(Vec3d v)
{
  string VEF = "VEF";
  string coords;
  for (int i = 0; i < 3; i++)
    if (v[i]) {
      if (v[i] != 1.0) {
        auto str = msg_str("%g", v[i]);
        if (str.find('e') != string::npos)
          str = msg_str("%f", v[i]);
        coords += str;
      }
      coords += VEF[i];
    }

  return coords;
}

static string M_pattern(int N)
{
  if (N < 1)
    return ""; // number out of range

  string pat = "[F";
  for (int i = 0; i < N + !is_even(N); i++) {
    double e_num = 2 * i;
    double v_num = N + 1 - e_num;
    pat += "," + coord_string(normalize_point(Vec3d(v_num, e_num, 0)));
  }
  int last_idx = (N + 3) / 2;
  pat += "]0_2_1e2e";

  for (int i = 2; i < last_idx; i++)
    pat += msg_str(",*0_%d_%d", i, i + 1);

  if (is_even(N)) {
    pat += msg_str(",%d_0v%dv", last_idx, last_idx);
    pat += msg_str(",%dE", last_idx);
  }

  return pat;
}

static string m_pattern(int N)
{
  if (N < 1)
    return ""; // number out of range

  string pat = "[F";
  for (int i = 0; i < N + 1; i += 2) {
    double e_num = i;
    double v_num = N - e_num;
    pat += "," + coord_string(normalize_point(Vec3d(v_num, e_num, 0)));
  }
  int last_idx = N / 2 + 1;

  pat += "]";

  for (int i = 1; i < last_idx; i++)
    pat += msg_str(",*0_%d_%d", i, i + 1);

  if (!is_even(N)) {
    pat += msg_str(",%d_0v%dv", last_idx, last_idx);
    pat += msg_str(",%dE", last_idx);
  }

  return pat;
}

static string b_pattern(int N)
{
  if (N < 1)
    return ""; // number out of range

  string pat = "[";
  for (int b = 1; b <= N + N % 2; b += 2)
    pat += coord_string(Vec3d((N - b), b, 1)) + ",";
  pat.back() = ']';

  pat += "0e0f,";

  for (int b = 0; b < N + N % 2; b += 2)
    pat += msg_str("%d_", b / 2);
  pat.back() = 'v';
  for (int b = 0; b < N - 1; b += 2)
    pat += msg_str("%d_", N / 2 - b / 2 - 1);
  if (pat.back() == '_')
    pat.pop_back();
  pat += 'e';

  for (int b = 0; b < N - 2 + N % 2; b += 2)
    pat += msg_str(",%d_%df%d_%df", b / 2, b / 2 + 1, b / 2 + 1, b / 2);
  if (N % 2)
    pat += msg_str(",%dE", N / 2);
  else
    pat += msg_str(",%dv%df", N / 2 - 1, N / 2 - 1);

  return pat;
}

static string g_pattern(int N)
{
  if (N < 1)
    return ""; // number out of range

  string pat = "[V";
  int divs = 2 * N + 1;
  for (int b = 0; b < N; b++) {
    const int e_coord = 2 * (b + 1);
    pat += "," + coord_string(Vec3d(divs - e_coord, e_coord, 0));
  }
  pat += ",F]";

  int F_idx = N + 1;

  if (N == 1)
    pat += msg_str("1_2F1_0V1E");
  else
    pat += msg_str("%d_1_0e1_2e", F_idx);

  for (int b = 0; b < N - 1; b++) {
    pat += msg_str(",%d", F_idx);
    int div_start = 2 * b + 1;
    bool past_center = false;
    for (int i = 0; i < 3; i++) {
      int div = div_start + i;
      char op = '_';
      if (div > N && !past_center) {
        op = 'v';
        past_center = true;
      }
      pat += msg_str("%c%d", op, (div <= N) ? div : 2 * N + 1 - div);
    }
    if (past_center)
      pat += 'v';
  }

  pat += msg_str(",%dE", N);

  return pat;
}

static string s_pattern(int N)
{
  if (N < 1)
    return ""; // number out of range

  string pat = "[";
  int divs = N - 1;
  for (int b = 0; b < divs / 2 + 1; b++) {
    const int e_coord = 2 * b;
    pat += coord_string(Vec3d(divs - e_coord, e_coord, 1)) + ',';
  }
  pat.back() = ']';

  auto div2idx = [](int idx, int div) {
    return (idx <= div / 2) ? idx : div - idx;
  };

  pat += msg_str("0V,%dE", divs / 2);

  if (divs) {
    pat += ",";
    bool past_center = false;
    for (int b = 0; b < divs; b++) {
      string op = (b) ? "_" : "";
      if (2 * b > divs && !past_center) {
        op = "v";
        past_center = true;
      }
      pat += msg_str("%s%d", op.c_str(), div2idx(b, divs));
    }
    if (past_center)
      pat += 'v';

    pat += 'F';
  }

  for (int b = 0; b < divs / 2; b++) {
    pat += msg_str(",%d_%df%df", div2idx(b, divs), div2idx(b + 1, divs),
                   div2idx(divs - (b + 1), divs));
    pat += msg_str(",%d_f%d_%df", div2idx(b, divs),
                   div2idx(divs - (b + 1), divs), div2idx(divs - b, divs));
  }
  if (is_even(N))
    pat += msg_str(",%df%dv%dvf", div2idx(divs / 2, divs),
                   div2idx(divs - divs / 2, divs),
                   div2idx(divs - (divs / 2 + 1), divs));

  return pat;
}

static string l_pattern(int N)
{
  if (N < 0)
    return ""; // number out of range

  string pat = "[";
  int divs = N + 1;
  for (int i = 0; i < divs; i++) {
    pat += coord_string(normalize_point(Vec3d(divs - i, 0, i))) + ',';
  }
  pat.back() = ']';

  pat += msg_str("%dF,", divs - 1);

  for (int i = 0; i < divs - 1; i++)
    pat += msg_str("%d_%dv%d_%dv,", i, i + 1, i + 1, i);

  pat += msg_str("0E");

  return pat;
}

static string L_pattern(int N)
{
  if (N < 0)
    return ""; // number out of range

  if (N == 0)
    return "[V,EF]1F,1e1_0e,1_0E"; // L_0 is its own operator

  string pat = "[";
  int divs = N + 1;
  for (int i = 0; i < divs; i++) {
    if (i == 0)
      pat += coord_string(Vec3d(1, 0, 0)) + ',';
    else if (i == 1)
      pat += coord_string(Vec3d(0, 5, 1)) + ',';
    else if (i % 2 == 0)
      pat += coord_string(Vec3d(1, 0, pow(5, i - 1))) + ',';
    else if (i % 2 == 1)
      pat += coord_string(Vec3d(0, 1, pow(2, i / 2) * pow(5, i - 2))) + ',';
  }
  pat.back() = ']';

  pat += msg_str("%dF,", divs - 1);

  for (int i = 0; i < divs - 1; i++) {
    const char *c = (i % 2) ? "ve" : "ev";
    pat += msg_str("%d%c%d_%d%c,", i + 1, c[0], i + 1, i, c[0]);
    pat += msg_str("%d%c%d_%d%c,", i, c[1], i, i + 1, c[1]);
  }

  pat += msg_str("0E");

  return pat;
}

class Tri {
public:
  vector<Vec3d> verts;
  string op;
  Tri(const vector<Vec3d> &verts, string op) : verts(verts), op(op) {}
  bool cart2bary(const Vec3d &P, Vec3d &Q);
};

bool Tri::cart2bary(const Vec3d &P, Vec3d &Q)
{
  auto &A = verts[0];
  auto &B = verts[1];
  auto &C = verts[2];
  double det = (B[1] - C[1]) * (A[0] - C[0]) + (C[0] - B[0]) * (A[1] - C[1]);
  Q[0] = ((B[1] - C[1]) * (P[0] - C[0]) + (C[0] - B[0]) * (P[1] - C[1])) / det;
  Q[1] = ((C[1] - A[1]) * (P[0] - C[0]) + (A[0] - C[0]) * (P[1] - C[1])) / det;
  Q[2] = 1 - Q[0] - Q[1];

  for (int i = 0; i < 3; i++) {
    if (double_eq(Q[i], 0))
      Q[i] = 0;
    if (double_lt(Q[i], 0))
      return false;
  }

  return true;
}

static string invert_op(char op)
{
  const char *c = strchr("VEF", op);
  if (c == nullptr) // mirror or no_op
    return string(1, op);
  else if (*c == 'V')
    return "fe";
  else if (*c == 'E')
    return "vf";
  else if (*c == 'F')
    return "ev";
  else
    return string(); // error
}

static string invert_op(const string &op)
{
  string inv;
  for (auto it = op.rbegin(); it < op.rend(); it++)
    inv += invert_op(*it);
  return inv;
}

static string u_pattern(int N, int M)
{
  // One of the parameters can be 0, but not both
  if (N < 1 && M < 1)
    return "";

  // Colour indexes for labelling tiles
  enum {
    tile_invalid = -1,   // invalid
    tile_unused = 0,     // outside the meta triangles of interest
    tile_center = 1,     // centred on the face centre
    tile_inner_edge = 2, // on an internal edge between meta triangles
    tile_normal = 3      // any other tile in the meta triangles
  };

  // Make a geodesic division of an equilateral triangle
  // Base model includes a triangle on each edge to capture the points
  // of the geodesic division outside the central triangle
  Geometry base;
  base.add_vert({-2, 0, 0});       // central triangle has vertex to the left
  base.add_vert({1, sqrt(3), 0});  // label V0
  base.add_vert({1, -sqrt(3), 0}); // label V1
  for (int i = 0; i < 3; i++)
    base.add_vert(-2 * base.verts(i));

  base.add_face({0, 1, 2}, Color(tile_unused)); // only central face has index
  base.add_face({3, 2, 1});
  base.add_face({4, 0, 2});
  base.add_face({5, 1, 0});
  Geometry geo;
  make_geodesic_planar(geo, base, M, N);

  // Find tiles that lie on a positive and negative meta triangle
  for (size_t i = 0; i < geo.faces().size(); i++) {
    auto cent = geo.face_cent(i); // tile centre, use for testing position
    bool on_central_face = (geo.colors(FACES).get(i).is_index());
    if (on_central_face && double_ge(cent[0], 0)) { // in sector to the right
      int tile_type = tile_invalid;
      if (double_eq(cent.len(), 0.0))
        tile_type = tile_center;
      else if (double_eq(sqrt(3) * cent[0], cent[1]))
        tile_type = tile_inner_edge; // top edge only (avoid duplication)
      else if (double_gt(sqrt(3) * fabs(cent[0]), fabs(cent[1])))
        tile_type = tile_normal; // inside sector

      // The right edge of tiles has no tiles shared with right original face
      // and so there is no possibility of duplication

      if (tile_type > tile_unused)
        geo.colors(FACES).set(i, tile_type);
    }
  }

  // Label some points of the meta triangle tiling of the base triangles
  const Vec3d V0 = base.verts(1);  // 60
  const Vec3d V1 = base.verts(2);  // -60
  const Vec3d V2 = base.verts(0);  // 180
  const Vec3d E01 = (V0 + V1) / 2; // 0
  const Vec3d E02 = (V0 + V2) / 2; // 120
  const Vec3d E12 = (V1 + V2) / 2; // -120
  const Vec3d F(0, 0, 0);          // centre (origin)
  const Vec3d F01 = 2 * E01;       // 0
  auto F_r = F - V2;               // Face centre to right

  // These are all the meta triangles which might include points
  vector<Tri> meta_triangles({
      Tri({V0, E01, F}, "_"),               // base
      Tri({V1, E01, F}, "v"),               // base negative
      Tri({V0, E02, F}, "e"),               // fan anticlockwise from base
      Tri({V0, E01, F01}, "f"),             // outside base
      Tri({V0, F_r + 0.5 * V0, F_r}, "ef"), // outside outside base
      Tri({V1, E01, F01}, "vf"),            // outside base negative
      Tri({V1, E12, F}, "ve"),              // fan clockwise from base negative
  });

  // base vertex index -> meta triangle index and barycentric coords
  map<int, pair<int, Vec3d>> vert_idx2tri_pt;

  // barycentric coords -> pattern point index
  auto vec_cmp = [](const Vec3d &v1, const Vec3d &v2) {
    return compare(v1, v2) == -1;
  };
  map<Vec3d, size_t, decltype(vec_cmp)> b_coords2pat_pt(vec_cmp);

  const int to_int_factor = M * M + M * N + N * N; // for integer coordinates
  string pattern_paths;
  string pattern_edge_paths;

  for (size_t i = 0; i < geo.faces().size(); i++) {
    auto col = geo.colors(FACES).get(i);
    int tile_type = col.is_index() ? col.get_index() : tile_invalid;
    int edge_start_idx = -1;       // face offset where an edge tile starts
    if (tile_type > tile_unused) { // a valid tile
      vector<int> tile_pts;
      vector<string> tile_ops;
      for (int v_idx = 0; v_idx < 3; v_idx++) {
        const auto P = geo.face_v(i, v_idx);

        Vec3d Q;    // point in barycentric coordinates for containing meta
        size_t tri; // will hold index of meta triangle where point found
        for (tri = 0; tri < meta_triangles.size(); tri++)
          if (meta_triangles[tri].cart2bary(P, Q)) // P is in meta triangle
            break;

        if (tri < meta_triangles.size()) { // point found in a meta triangle
          Q = normalize_point(Q * to_int_factor); // to low integer coordinates
          // If Q has already been seen, get the pattern index number,
          // otherwise set it
          auto it_bool = b_coords2pat_pt.insert({Q, b_coords2pat_pt.size()});
          auto pat_pt_idx = it_bool.first->second; // idx is iterator value

          string op = meta_triangles[tri].op;
          if (tile_type == tile_center) {
            tile_ops.push_back((op == "_" || is_even(op.size())) ? "_" : "v");
            tile_pts.push_back(pat_pt_idx);
            break; // only need one point for the central triangle
          }
          else { // tile_inner_edge or tile_normal
            tile_ops.push_back(op);
            tile_pts.push_back(pat_pt_idx);
          }
        }
        else {                          // point not found in a meta triangle
          if (tile_type != tile_center) // not an error for centre triangle
            return "";                  // INTERNAL ERROR
        }

        // test for edge on E, to convert to tile that winds E
        if (compare(geo.edge_cent(
                        {geo.faces(i, v_idx), geo.faces_mod(i, v_idx + 1)}),
                    E01) == 0)
          edge_start_idx = v_idx;
      }

      // Have found tiles as paths of points on meta triangles, need to
      // convert to paths of points flipping between meta triangles
      string tile_str;
      if (tile_type == tile_center) {
        // No extra handling required (handled when point was found)
        if (tile_ops[0] == "v")
          tile_str += '-';
        tile_str += to_string(tile_pts[0]) + "F";
      }
      else { // tile_inner_edge or tile_normal
        // handling will depend on whether path is positive or negative
        // get the index of the first _, otherwise the first v
        int first = -1; // index of first _ if found, or first v
        for (int i = 0; i < 3; i++) {
          if (tile_ops[i] == "_") { // if found then done
            first = i;
            break;
          }
          else if (tile_ops[i] == "v" && first == -1) // only record first v
            first = i;
        }

        if (first == -1) { // shouldn't happen
          return "";       // INTERNAL ERROR
        }
        else if (tile_ops[first] == "_") {
          for (int i = 0; i < 3; i++) {
            int idx = (i + first) % 3; // index starts at operator '_'
            string start_op = (i) ? tile_ops[idx] : ""; // skip leading op
            tile_str +=
                start_op + to_string(tile_pts[idx]) + invert_op(start_op);
          }
        }
        else if (tile_ops[first] == "v") {
          tile_str = "-";
          for (int i = 0; i < 3; i++) {
            int idx = (i + first) % 3; // index starts at operator 'v'
            // conjugate the operation with 'v' to reuse the base operators
            string start_op = (i) ? "v" + tile_ops[idx] : ""; // skip leading op
            tile_str +=
                start_op + to_string(tile_pts[idx]) + invert_op(start_op);
          }
        }
        else {
          return ""; // INTERNAL ERROR, should not occur
        }

        if (edge_start_idx != -1) { // edge tile found
          int e0 = edge_start_idx;
          int e1 = (e0 + 1) % 3;
          // choose e0 to prefer '_' over 'v'
          if ((tile_ops[e1] == "_") ||
              ((tile_ops[e1] == "v") && (tile_ops[e0] != "_")))
            std::swap(e0, e1);
          pattern_edge_paths = ((tile_ops[e0] == "v") ? "-" : "") +
                               to_string(tile_pts[e0]) + "E,";
        }
      }
      if (!tile_str.empty())
        pattern_paths += tile_str + ",";
    }
  }

  if (M + N == 1)
    pattern_edge_paths = "0E,"; // isn't set by the algorithm for 1,0

  // make an array of the pattern points in index order
  vector<Vec3d> pat_points(b_coords2pat_pt.size());
  for (const auto &kv : b_coords2pat_pt)
    pat_points[kv.second] = kv.first;

  // write the points in the pattern format
  string pattern = "[";
  for (const auto &coords : pat_points)
    pattern += coord_string(coords) + ",";
  pattern.back() = ']';

  pattern += pattern_paths + pattern_edge_paths;
  pattern.pop_back(); // remove the final ','
  pattern = Tile::normalize_pattern(pattern);

  return pattern;
}

static string oe_pattern(char oper, int N, int M)
{
  // One of the parameters can be 0, but not both
  if (N < 1 && M < 1)
    return "";

  // Colour indexes for labelling tiles
  enum {
    tile_invalid = -1, // invalid
    tile_face = 0,     // centred on a face centre
    tile_edge = 1,     // centred on an edge centre
    tile_vert = 2,     // centred on a vertex centre
    tile_normal = 3,   // any other tile in a meta triangle
  };

  // Label some points of the meta triangle tiling of the base triangles
  const Vec3d V0(1, 1, 0);         // 45
  const Vec3d V1(1, -1, 0);        // -45
  const Vec3d V2(-1, -1, 0);       // -135
  const Vec3d V3(-1, 1, 0);        // 135
  const Vec3d E01 = (V0 + V1) / 2; // 0
  const Vec3d E30 = (V3 + V0) / 2; // 90
  const Vec3d E12 = (V1 + V2) / 2; // -90
  const Vec3d F(0, 0, 0);          // centre (origin)
  const Vec3d F01 = 2 * E01;       // 0

  // These are all the meta triangles which might include points
  vector<Tri> meta_triangles({
      Tri({V0, E01, F}, "_"),    // base
      Tri({V1, E01, F}, "v"),    // base negative
      Tri({V0, E30, F}, "e"),    // fan anticlockwise from base
      Tri({V0, E01, F01}, "f"),  // outside base
      Tri({V1, E01, F01}, "vf"), // outside base negative
      Tri({V1, E12, F}, "ve"),   // fan clockwise from base negative
  });

  // barycentric coords -> pattern point index
  auto vec_cmp = [](const Vec3d &v1, const Vec3d &v2) {
    return compare(v1, v2) == -1;
  };
  map<Vec3d, size_t, decltype(vec_cmp)> b_coords2pat_pt(vec_cmp);

  const int to_int_factor = 2 * (M * M + N * N); // XXXfor integer coordinates
  string pattern_paths;
  string pattern_edge_paths;

  Geometry geom;
  geom.add_verts({V0, V1, V2, V3, F});
  geom.add_face({0, 1, 4});
  double tile_edge_len = 2 / sqrt(N * N + M * M);
  Vec3d vx = Vec3d(N, M, 0).with_len(tile_edge_len);
  Vec3d vy = Vec3d(-M, N, 0).with_len(tile_edge_len);
  for (int i = -(N + M + 1); i < N + M + 1; i++) {   // enough
    for (int j = -(N + M + 1); j < N + M + 1; j++) { // enough
      Vec3d off = -0.5 * (vx + vy);
      auto C = is_even(M + N + (oper == 'e')) * off + i * vx + j * vy;
      const auto x = C[0];
      const auto y = C[1];
      int tile_type = tile_invalid;
      if (double_ge(x, 0) && double_le(x, 1) && // 0 <= x <= 1
          double_le(std::abs(y), 1)) {          // -1 <= y <= 1
        if (compare(C, F) == 0)
          tile_type = tile_face;
        else if (compare(C, E01) == 0)
          tile_type = tile_edge;
        else if (compare(C, V0) == 0)
          tile_type = tile_vert;
        else if ((double_lt(x, 1) && double_gt(x, std::abs(y))) || // interior
                 double_eq(x, y) ||                                // e mirror
                 (double_eq(x, 1) && double_ge(y, 0)))             // f mirror
          tile_type = tile_normal;
      }
      if (tile_type == tile_invalid)
        continue;
      // The tile is valid

      vector<Vec3d> points = {C + (vx + vy) / 2, C + (vx - vy) / 2,
                              C + (-vx - vy) / 2, C + (-vx + vy) / 2};

      int vsz = geom.verts().size();
      geom.add_verts(points);
      geom.add_face({vsz, vsz + 1, vsz + 2, vsz + 3});

      int edge_start_idx = -1; // face offset where an edge tile starts
      vector<int> tile_pts;
      vector<string> tile_ops;
      for (int v_idx = 0; v_idx < 4; v_idx++) {
        const auto P = points[v_idx];

        Vec3d Q;    // point in barycentric coordinates for containing meta
        size_t tri; // will hold index of meta triangle where point found
        for (tri = 0; tri < meta_triangles.size(); tri++)
          if (meta_triangles[tri].cart2bary(P, Q)) // P is in meta triangle
            break;

        if (tri < meta_triangles.size()) { // point found in a meta triangle
          Q = normalize_point(Q * to_int_factor); // to low integer coords
          // If Q has already been seen, get the pattern index number,
          // otherwise set it
          auto it_bool = b_coords2pat_pt.insert({Q, b_coords2pat_pt.size()});
          auto pat_pt_idx = it_bool.first->second; // idx is iterator value

          string op = meta_triangles[tri].op;
          if (tile_type == tile_normal) {
            tile_ops.push_back(op);
            tile_pts.push_back(pat_pt_idx);
          }
          else { // tile cycles an element
            tile_ops.push_back((op == "_" || is_even(op.size())) ? "_" : "v");
            tile_pts.push_back(pat_pt_idx);
            if ((tile_type != tile_edge || tile_ops.size() > 1) && oper != 'X')
              break; // need one point for F and V, and two for E
          }
        }
        else {                            // point not found in a meta triangle
          if (tile_type == tile_normal) { // not an error for centre triangle
            return "";                    // INTERNAL ERROR
          }
        }

        if (edge_start_idx == -1) { // edge element not yet found on tile
          // test for edge on E, to convert to tile that winds E
          int next_v_idx = (v_idx + 1 + (oper == 'X')) % 4;
          if (edge_start_idx == -1 &&
              compare((points[v_idx] + points[next_v_idx]) / 2, E01) == 0)
            edge_start_idx = v_idx;
        }
      }

      // Have found tiles as paths of points on meta triangles, need to
      // convert to paths of points flipping between meta triangles
      string tile_str;
      if (tile_type != tile_normal && oper != 'X') {
        // No extra handling required (handled when point was found)
        if (tile_ops[0] == "v")
          tile_str += '-';
        string elems("FEV");
        tile_str += to_string(tile_pts[0]);
        if (tile_type == tile_edge)
          tile_str += "v" + to_string(tile_pts[1]) + "v";
        tile_str += elems[tile_type];
      }
      else { // tile_normal
        // handling will depend on whether path is positive or negative
        // get the index of the first _, otherwise the first v
        int first = -1; // index of first _ if found, or first v
        for (int i = 0; i < 4; i++) {
          if (tile_ops[i] == "_") { // if found then done
            first = i;
            break;
          }
          else if (tile_ops[i] == "v" && first == -1) // only record first v
            first = i;
        }

        // Is this tile a cross operator triangle
        bool is_cross_tri_tile =
            (oper == 'X') &&
            (tile_type == tile_normal || tile_type == tile_edge) &&
            (double_eq(x, 1) && double_ge(y, 0));

        if (first == -1) { // shouldn't happen
          return "";       // INTERNAL ERROR
        }
        else if (tile_ops[first] == "_") {
          if (is_cross_tri_tile) { // split quad in two
            for (int t = 0; t < 2 - (tile_type == tile_edge); t++) {
              if (t > 0)
                tile_str += ","; // more than one tile, so need separator
              for (int i = 0; i < 3; i++) {
                int idx = (t) ? i + first : 4 - i + first; // starts at op '_'
                idx %= 4;
                string start_op = (i) ? tile_ops[idx] : ""; // skip leading op
                tile_str +=
                    start_op + to_string(tile_pts[idx]) + invert_op(start_op);
              }
            }
          }
          else { // usual quad
            for (int i = 0; i < 4; i++) {
              int idx = (i + first) % 4; // index starts at operator '_'
              string start_op = (i) ? tile_ops[idx] : ""; // skip leading op
              tile_str +=
                  start_op + to_string(tile_pts[idx]) + invert_op(start_op);
            }
          }
        }
        else if (tile_ops[first] == "v") {
          tile_str = "-";
          for (int i = 0; i < 4; i++) {
            int idx = (i + first) % 4; // index starts at operator 'v'
            // conjugate the operation with 'v' to reuse the base
            // operators
            string start_op = (i) ? "v" + tile_ops[idx] : ""; // skip leading op
            tile_str +=
                start_op + to_string(tile_pts[idx]) + invert_op(start_op);
          }
        }
        else {
          return ""; // INTERNAL ERROR, should not occur
        }

        if (edge_start_idx != -1) { // edge tile found
          int e0 = edge_start_idx;
          int e1 = (e0 + 1 + (oper == 'X')) % 4;
          // choose e0 to prefer '_' over 'v'
          if ((tile_ops[e1] == "_") ||
              ((tile_ops[e1] == "v") && (tile_ops[e0] != "_")))
            std::swap(e0, e1);
          pattern_edge_paths = ((tile_ops[e0] == "v") ? "-" : "") +
                               to_string(tile_pts[e0]) + "E,";
        }
      }
      if (!tile_str.empty())
        pattern_paths += tile_str + ",";
    }
  }

  if (M + N == 1)
    pattern_edge_paths = "0E,"; // isn't set by the algorithm for 1,0

  // make an array of the pattern points in index order
  vector<Vec3d> pat_points(b_coords2pat_pt.size());
  for (const auto &kv : b_coords2pat_pt)
    pat_points[kv.second] = kv.first;

  // write the points in the pattern format
  string pattern = "[";
  for (const auto &coords : pat_points)
    pattern += coord_string(coords) + ",";
  pattern.back() = ']';

  pattern += pattern_paths + pattern_edge_paths;
  pattern.pop_back(); // remove the final ','
  pattern = Tile::normalize_pattern(pattern);

  return pattern;

  // geom.write("tmp.off");
}

struct ConwayOperator {
  std::string operator_short;
  std::string operator_name;
  std::string pattern;
  int num_params;
};

// clang-format off
ConwayOperator conway_operator_list[]{
    // Equivalent: d, a, s
    {"d",   "dual",           "[F]0V,0E", 0},
    {"a",   "ambo",           "[E]0F,0V", 0},
    {"S",   "seed",           "[V]0E,0F", 0},

    {"j",   "join",           "[F,V]0_1E", 0},

    // Equivalent: k, n, u
    {"k",   "kis",            "[F,V]0_1v1v,1E", 0},
    {"n",   "needle",         "[V,F]0_1f1f,1E", 0},
    {"u",   "subdivide",      "[V,E]1F,0_1e1e", 2},

    // Equivalent: t, z, e (tile order to match e0=z and e1=e
    {"t",   "truncate",       "[VE]0V0E,0V,0E", 0},
    {"z",   "zip",            "[EF]0E0F,0F,0E", 0},
    {"e",   "expand",         "[FV]0V,0F,0F0V", 2},

    // Symmetric: s, m, b
    {"s",   "snub",           "[VEF]0V,0E,0F,0V0E0F", 1},
    {"m",   "meta",           "[V,E,F]*0_1_2", 1},
    {"b",   "bevel",          "[VEF]0e0f,0v0e,0f0v", 1},

    {"o",   "ortho",          "[V,E,F]1_0e1_2e", 2},
    {"g",   "gyro",           "[F,VE,V]1_0F1_2V1E,1E", 1},
    {"c",   "chamfer",        "[V,VF]1F,0_1v1f", 0},
    {"l",   "loft",           "[V,VF]1F,0_1v1_0v,0E", 1},
    {"p",   "propellor",      "[V,VEF]1F,1_0V1E1F,1E", 0},
    {"q",   "quinto",         "[V,E,EF]2F,0_1_2e2_1e", 0},
    {"L_0", "joined-lace",    "[V,EF]1F,1e1_0e,1_0E", 0},
    {"G",   "opposite-lace",  "[V,EF]1F,1e1_0e,0_1f1f,1E", 0},
    {"L",   "lace",           "[V,EF]1F,1e1_0e,1_0v0v,0E", 1},
    {"K",   "stake",          "[V,EF,F]0_1_2e1e,1_0v0v,0E", 0},
    {"M",   "edge-medial",    "[F,3V,V2E]0_2_1e2e,2_0v2v,2E", 1},
    {"J",   "joined-medial",  "[F,V,EF]*0_1_2,1_2E", 0},
    {"X",   "cross",          "[V,E,F,VF]3_1v3_2v,*0_1_3", 1},
    {"w",   "whirl",          "[VF,VE,V]0F,0_1V2_1E1_0F,1E", 0},
    {"E",   "ethel",          "[V,VE,VF]0_1_2e1e,2F,1_2v2f", 0},
    {"W",   "waffle",         "[V,E,F,V2E,VF]0_4_3f4f,2_4_3v3_4v,3E", 0},
    {"B",   "bowtie",         "[V,E,F,VE,EF]1_3_4,0_3_4_2e4_1_3e", 0},
};
// clang-format on

}; // namespace

namespace anti {

Status Tiling::read_conway(const string &op)
{
  const auto unset = std::numeric_limits<int>::max(); // unset parameter value
  int param1 = unset;
  int param2 = unset;
  char op_char; // Conway operation letter
  char buff;
  bool valid_format =
      (sscanf(op.c_str(), "%c%c", &op_char, &buff) == 1 ||             // X
       sscanf(op.c_str(), "%c%d%c", &op_char, &param1, &buff) == 2 ||  // X1
       sscanf(op.c_str(), "%c_%d%c", &op_char, &param1, &buff) == 2 || // X_1
       sscanf(op.c_str(), "%c_%d_%d%c", &op_char, &param1, &param2,
              &buff) == 3 // X_1_2
      );
  if (!valid_format)
    return Status::error("Conway operator '" + op + "': invalid format");

  // Look up operator
  int op_idx = -1;             // index number of operator in list
  vector<string> param_ops(3); // op letters for each number of paramieters
  int last_op = sizeof(conway_operator_list) / sizeof(conway_operator_list[0]);
  for (int i = 0; i < last_op; i++) {
    const char c = conway_operator_list[i].operator_short[0];
    param_ops[conway_operator_list[i].num_params] += c;
    if (c == op_char)
      op_idx = i;
  }

  auto prefix = msg_str("Conway operator %c: ", op_char);
  if (op_idx < 0) // not found
    return Status::error(prefix + "operator not recognised");

  int number_of_params = (param1 != unset) + (param2 != unset);
  if (number_of_params > conway_operator_list[op_idx].num_params) {
    if (conway_operator_list[op_idx].num_params == 0)
      return Status::error(prefix + "does not take a parameter");
    else
      return Status::error(prefix + "too many parameters for operator");
  }
  if (param1 < 0 || param2 < 0)
    return Status::error(prefix + "parameter cannot be negative");

  string pat;            // final wythoff pattern
  if (param1 == unset) { // no parameter given, use base pattern
    pat = conway_operator_list[op_idx].pattern;
  }
  else { // parameter given, use sequence processing
    if (strchr(param_ops[1].c_str(), op_char)) {
      if (op_char == 'M')
        pat = M_pattern(param1);
      else if (op_char == 'm')
        pat = m_pattern(param1);
      else if (op_char == 'b')
        pat = b_pattern(param1);
      else if (op_char == 'g')
        pat = g_pattern(param1);
      else if (op_char == 's')
        pat = s_pattern(param1);
      else if (op_char == 'l')
        pat = l_pattern(param1);
      else if (op_char == 'L')
        pat = L_pattern(param1);
      else if (op_char == 'X')
        pat = oe_pattern('X', param1, param1);
      else
        return Status::error(
            prefix + "one-parameter operator not handled (internal error)");

      if (pat == "")
        return Status::error(prefix + msg_str("invalid parameter %d", param1));
    }
    else if (strchr(param_ops[2].c_str(), op_char)) {
      if (op_char == 'u')
        pat = u_pattern((param1 == unset) ? 3 : param1,
                        (param2 == unset) ? 0 : param2);
      if (op_char == 'o')
        pat = oe_pattern('o', (param1 == unset) ? 2 : param1,
                         (param2 == unset) ? 0 : param2);
      if (op_char == 'e')
        pat = oe_pattern('e', (param1 == unset) ? 2 : param1,
                         (param2 == unset) ? 0 : param2);
      if (pat == "") {
        if (param2 == unset)
          return Status::error(prefix +
                               msg_str("invalid parameter %d", param1));
        else
          return Status::error(
              prefix + msg_str("invalid parameters %d, %d", param1, param2));
      }
    }
  }

  return read_pattern(pat);
}

string Tiling::pattern_string()
{
  string pat = "[";
  for (auto v : points)
    pat += coord_string(v.first) + ",";
  pat.back() = ']';
  for (auto path : pat_paths)
    pat += path.tile_string() + ",";

  pat.pop_back();
  return pat;
}

void Tiling::print_conway_list(FILE *ofile)
{
  fprintf(ofile, "%-4s%-5s%-15s%s\n", "Op", "Args", "Description",
          "Tiling Pattern");
  fprintf(ofile, "%-4s%-5s%-15s%s\n", "--", "----", "-----------",
          "--------------");
  int last_op = sizeof(conway_operator_list) / sizeof(conway_operator_list[0]);
  for (int i = 0; i < last_op; i++) {
    ConwayOperator op = conway_operator_list[i];
    string params = op.num_params ? msg_str("[%d]", op.num_params) : "   ";
    fprintf(ofile, "%-4s%-5s%-15s%s\n", op.operator_short.c_str(),
            params.c_str(), op.operator_name.c_str(), op.pattern.c_str());
  }
  fprintf(ofile, R"(
Some operators are part of a sequence and take an optional integer parameter.
  g, l, L:               1 is the base, 0 is a lower level operator
                         (except L_0, which is a standalone operator)
  m, b, M, s, u, o, e:   2 is the base, 1 is a lower level operator and
                         0 is invalid

Operators g, o, e optionally take two integer parameters, not both zero.

Examples e, M, m_3, g_3, s_1, u_3, u_2_3, o_3_3, e_0_4

Some Conway operators, like t and k, take a number to filter the elements
that the pattern will be applied to, but this is not supported.
)");
}

} // namespace anti
